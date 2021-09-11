/**
 * \file HexahedronMesh.h
 * \author Huayi Wei/weihuayi@xtu.edu.cn
 * \date 2021.09.08
 * \brief 非结构四面体网格类
 * \details
 * \note
 *
 */
#ifndef HexahedronMesh_h
#define HexahedronMesh_h

#include <vector>
#include <array>
#include <map>

#include "MeshToplogy.h"


namespace OF {
namespace Mesh {

/**
 * \brief 非结构六面体网格类
 *
 *
 * \note 这里用 std::vector 做为容器, 用std::array 整数数组代表 edge, face, 和 cell
 *  实体.
 *
 */
template<typename GK, typename NODE, typename VECTOR>
class HexahedronMesh 
{
public:
  typedef NODE Node;
  typedef VECTOR Vector;

  typedef typename GK::Int I;
  typedef typename GK::Float F;

  typedef typename GK::Vector_2 Vector_2;


  // 实体类型
  typedef typename std::array<I, 2> Edge;
  typedef typename std::array<I, 4> Face;
  typedef typename std::array<I, 8> Cell;

  // 规则化的实体关系类型，这里的规则化是指一种实体有固定个数的另一种邻接实体
  // 每个单元有 6 个四边形面
  // 每个单元有 12 条边
  typedef typename std::array<I, 4> Face2cell; // 如果是边界面， 则左右单元的编号一样
  typedef typename std::array<I, 6> Cell2face;
  typedef typename std::array<I, 6> Cell2cell;
  typedef typename std::array<I, 12> Cell2edge;

  // 非规则化的拓扑关系， 如共享一个节点的单元个数是不固定的
  // 共享一条边的单元个数也是不固定的
  typedef MeshToplogy<I, std::vector<I> > Toplogy;

  // 迭代子类型
  typedef std::vector<Node> NodeArray;
  typedef std::vector<Edge> EdgeArray;
  typedef std::vector<Face> FaceArray;
  typedef std::vector<Cell> CellArray;

  typedef typename NodeArray::iterator NodeIterator;
  typedef typename NodeArray::const_iterator ConstNodeIterator;

  typedef typename EdgeArray::iterator EdgeIterator;
  typedef typename EdgeArray::const_iterator ConstEdgeIterator;

  typedef typename FaceArray::iterator FaceIterator;
  typedef typename FaceArray::const_iterator ConstFaceIterator;

  typedef typename CellArray::iterator CellIterator;
  typedef typename CellArray::const_iterator ConstCellIterator;

public:

  /**
   * \brief 构造函数
   *
   */
  HexahedronMesh()
  {
  }

  /**
   * \brief 插入节点到网格中 
   *
   */
  void insert(const Node & node)
  {
    m_node.push_back(node);
  }

  /**
   * \brief 插入单元到网格中 
   *
   */
  void insert(const Cell & cell)
  {
    m_cell.push_back(cell);
  }

  /**
   * \brief 获得网格节点个数
   *
   */
  I number_of_nodes()
  {
    return m_node.size();
  }

  /**
   * \brief 获得网格单元个数
   *
   */
  I number_of_cells()
  {
    return m_cell.size();
  }

  /**
   * \brief 获得网格边的个数
   *
   */
  I number_of_edges()
  {
    return m_edge.size();
  }

  /**
   * \brief 获得网格面的个数
   *
   */
  I number_of_faces()
  {
    return m_face.size();
  }


  /**
   * \brief 获得每个网格单元顶点的个数
   *
   */
  static I number_of_nodes_of_each_cell()
  {
      return 8;
  }

  /**
   * \brief 获得每个网格单元顶点的个数
   * 
   */
  static I number_of_vertices_of_each_cell()
  {
      return 8;
  }

  /**
   * \brief 获得网格几何维数
   * 
   */
  static I geo_dimension()
  {
      return 3;
  }

  /**
   * \brief 获得网格拓扑维数
   * 
   */
  static I top_dimension()
  {
      return 3;
  }

  /**
   * \brief 获得六面体单元顶点对应的 VTK 单元编号规则
   * 
   */
  int* vtk_cell_index()
  {
    return m_vtkidx;
  }

  /**
   * \brief 获得六面体单元对应的 VTK 单元类型
   * 
   */
  I vtk_cell_type(I TD=3)
  {
    return 12; // VTK_HEXAHEDRON = 10
  }
  
  /**
   * \brief 获得六面体单元顶点对应的 Abaqus 单元编号规则
   * 
   */
  int* abaqus_cell_index()
  {
    return m_vtkidx;
  }


  /**
   * \brief 初始化网格的拓扑邻接关系 
   * 
   * \note 实现上提取所有单元的六个面，并对面的顶点进行排序，并利用 map
   * 顶点排序后的面一个唯一连续的整数编号
   */
  void init_top()
  {
      m_face2cell.clear();
      m_cell2face.clear();

      auto NN = number_of_nodes();
      auto NC = number_of_cells();
      m_face2cell.reserve(2*NC); //TODO: 欧拉公式?
      m_cell2face.resize(NC);
      std::map<Face, I> fidxmap; // sorted face 到整体编号的映射

      I NF = 0;
      // 遍历所有单元
      for(I i = 0; i < NC; i++)
      {
          for(I j = 0; j < 6; j++)
          {
             auto f = local_sorted_face(i, j); // 第 i 个单元的第 j 个 sorted face
             auto it = fidxmap.find(f);
             if(it == fidxmap.end())
             {
                m_cell2face[i][j] = NF;
                fidxmap.insert(std::pair<Face, I>(f, NF));
                m_face2cell.push_back(Face2cell{i, i, j, j});
                NF++;
             }
             else
             {
                m_cell2face[i][j] = it->second;
                m_face2cell[it->second][1] = i;
                m_face2cell[it->second][3] = j;
             }
          }
      }
      fidxmap.clear();

      m_face.resize(NF);
      for(I i = 0; i < NF; i++)
      {
          auto & c = m_cell[m_face2cell[i][0]];
          auto j = m_face2cell[i][2];
          m_face[i][0] = c[m_localface[j][0]];
          m_face[i][1] = c[m_localface[j][1]];
          m_face[i][2] = c[m_localface[j][2]];
          m_face[i][3] = c[m_localface[j][3]];
      }

      std::map<Edge, I> eidxmap; // sorted edge 到整体唯一编号的映射
      m_cell2edge.resize(NC);
      I NE = 0;
      for(I i = 0; i < NC; i++)
      {
          for(I j = 0; j < 12; j++)
          { 
              auto & c = m_cell[i];
              auto e = local_sorted_edge(i, j); 
              auto it = eidxmap.find(e);
              if(it == eidxmap.end())
              {
                  m_cell2edge[i][j] = NE;
                  eidxmap.insert(std::pair<Edge, I>(e, NE));
                  m_edge.push_back(Edge{c[m_localedge[j][0]], c[m_localedge[j][1]]});
                  NE++; 
              }
             else
             {
                m_cell2edge[i][j] = it->second;
             }
          }
      }
  }


  /**
   * \defgroup 实体迭代子
   * \{
   */
  NodeIterator node_begin()
  {
      return m_node.begin();
  }

  NodeIterator node_end()
  {
      return m_node.end();
  }

  EdgeIterator edge_begin()
  {
      return m_edge.begin();
  }

  EdgeIterator edge_end()
  {
      return m_edge.end();
  }

  FaceIterator face_begin()
  {
      return m_face.begin();
  }

  FaceIterator face_end()
  {
      return m_face.end();
  }

  CellIterator cell_begin()
  {
      return m_cell.begin();
  }

  CellIterator cell_end()
  {
      return m_cell.end();
  }
  /**
   * \}
   */


  std::vector<Node> & nodes()
  {
      return m_node;
  }

  std::vector<Cell> & cells()
  {
      return m_cell;
  }

  Node & node(const I i)
  {
      return m_node[i];
  }

  Edge & edge(const I i)
  {
      return m_edge[i];
  }

  Face & face(const I i)
  {
      return m_face[i];
  }

  Cell & cell(const I i)
  {
      return m_cell[i];
  }


  Face2cell & face_to_cell(const I i)
  {
      return m_face2cell[i];
  }

  Cell2face & cell_to_face(const I i)
  {
      return m_cell2face[i];
  }

  Cell2edge & cell_to_edge(const I i)
  {
      return m_cell2edge[i];
  }

  /**
   * \brief 计算第 i 条边的长度
   */
  F edge_measure(const I i)
  {
      auto & e = m_edge[i];
      auto v = m_node[e[1]] - m_node[e[0]];
      return std::sqrt(v.squared_length());
  }

  /**
   * \brief 计算所有边的长度
   */
  void edge_measure(std::vector<F> & measure)
  {
      auto NE = number_of_edges();
      measure.resize(NE);
      for(I i = 0; i < NE; i++)
          measure[i] = edge_measure(i);
  }

  /*
   * \brief 计算第 i 个四边形面的面积
   *
   * \note 一个四边形面可以由如下的映射来定义
   *
   * \f[
   *    x = x_0\phi_0(u, v) + x_1\phi_1(u, v) + x_2\phi_2(u, v) 
   *      + x_3\phi_3(u, v)
   * \f]
   *
   * 其中 \f$(u, v)\f$ 是参考单元 \f$[0, 1]^2\f$ 上的自变量，\f$\phi_i\f$
   * 是四个顶点对应的形函数
   */
  F face_measure(const I i)
  {
      return 0.0; 
  }

  /*
   * \brief 计算所有四边形面的面积
   */
  void face_measure(std::vector<F> & measure)
  {
    auto NF = number_of_faces();
    measure.resize(NF);
    for(I i = 0; i < NF; i++)
        measure[i] = face_measure(i);
  }

  /*
   * \brief 计算六面体单元的体积
   *
   * \note
   *
   * \f[
   *    x = x_0\phi_0(u, v, w) + x_1\phi_1(u, v, w) + x_2\phi_2(u, v, w) 
   *      + x_3\phi_3(u, v, w) + x_4\phi_4(u, v, w) + x_5\phi_5(u, v, w)
   *      + x_6\phi_6(u, v, w) + x_7\phi_7(u, v, w)
   * \f]
   */
  F cell_measure(const I i)
  {
    return 0.0;
  }

  /*
   * 计算所有六面体单元的体积
   *
   * Notes
   * -----
   *
   */
  void cell_measure(std::vector<F> & measure)
  {
      auto NC = number_of_cells();
      measure.resize(NC);
      for(I i = 0; i < NC; i++)
          measure[i] = cell_measure(i);
  }

  Node edge_barycenter(const I i)
  {
      auto & e = m_edge[i];
      F x = (m_node[e[0]][0] + m_node[e[1]][0])/2.0;
      F y = (m_node[e[0]][1] + m_node[e[1]][1])/2.0;
      F z = (m_node[e[0]][2] + m_node[e[1]][2])/2.0;
      return Node(x, y, z);
  }

  void edge_barycenter(const I i, Node & node)
  {
      auto & e = m_edge[i];
      node[0] = (m_node[e[0]][0] + m_node[e[1]][0])/2.0;
      node[1] = (m_node[e[0]][1] + m_node[e[1]][1])/2.0;
      node[2] = (m_node[e[0]][2] + m_node[e[1]][2])/2.0;
  }

  Node face_barycenter(const I i)
  {
      auto & f = m_face[i];
      F x = (m_node[f[0]][0] + m_node[f[1]][0] + m_node[f[2]][0] + m_node[f[3]][0])/4.0;
      F y = (m_node[f[0]][1] + m_node[f[1]][1] + m_node[f[2]][1] + m_node[f[3]][1])/4.0;
      F z = (m_node[f[0]][2] + m_node[f[1]][2] + m_node[f[2]][2] + m_node[f[3]][2])/4.0;
      return Node(x, y, z);
  }

  void face_barycenter(const I i, Node & node)
  {
      auto & f = m_face[i];
      node[0] = (m_node[f[0]][0] + m_node[f[1]][0] + m_node[f[2]][0] + m_node[f[3]][0])/4.0;
      node[1] = (m_node[f[0]][1] + m_node[f[1]][1] + m_node[f[2]][1] + m_node[f[3]][1])/4.0;
      node[2] = (m_node[f[0]][2] + m_node[f[1]][2] + m_node[f[2]][2] + m_node[f[3]][2])/4.0;
  }

  Node cell_barycenter(const I i)
  {
    auto & c = m_cell[i];
    F x = m_node[c[0]][0];
    F y = m_node[c[0]][1];
    F z = m_node[c[0]][2];
    for(int j = 1; j < 8; j++)
    {
      x += m_node[c[j]][0];
      y += m_node[c[j]][1];
      z += m_node[c[j]][2];
    }
    return Node(x/8.0, y/8.0, z/8.0);
  }

  void cell_barycenter(const I i, Node & node)
  {
    auto & c = m_cell[i];
    auto GD = geo_dimension();
    for(int d = 0; d < GD; d++)
    {
      node[d] = m_node[c[0]][d];
      for(int j = 1; j < 8; j++)
        node[d] += m_node[c[j]][d];
      node[d] /= 8.0;
    }
  }

  /*
   * 计算第 i 个单元所有形函数在参考积分点处的值
   *
   *
   *  x = phi_0 x_0 + phi_1 x_1 + ... + phi_7 x_7
   *
   *  (u, v, w) 
   *  
   *
   */

  template<typename RNode>
  void cell_shape_function(I i, const RNode & rnode, std::vector<F> & phi)
  {
    assert(rnode.size() == 3);
    phi.resize(8);

    phi[0] = (1-rnode[0])*(1-rnode[1])*(1-rnode[2]);
    phi[1] = rnode[0]*(1 - rnode[1])*(1-rnode[2]);
    phi[2] = rnode[0]*rnode[1]*(1 - rnode[2]);
    phi[3] = (1-rnode[0])*rnode[1]*(1-rnode[2]);

    phi[4] = (1-rnode[0])*(1-rnode[1])*rnode[2];
    phi[5] = rnode[0]*(1-rnode[1])*rnode[2];
    phi[6] = rnode[0]*rnode[1]*rnode[2];
    phi[7] = (1-rnode[0])*rnode[1]*rnode[2];
  }


  template<typename RNode>
  void cell_grad_shape_function(I i, const RNode & rnode, std::vector<Vector> & gphi)
  {

    assert(rnode.size() == 3);

    gphi.resize(8);

    gphi[0][0] = -(1-rnode[1])*(1-rnode[2]);
    gphi[0][1] = -(1-rnode[0])*(1-rnode[2]);
    gphi[0][2] = -(1-rnode[0])*(1-rnode[1]);

    gphi[1][0] = (1 - rnode[1])*(1-rnode[2]);
    gphi[1][1] = -rnode[0]*(1-rnode[2]);
    gphi[1][2] = -rnode[0]*(1 - rnode[1]);

    gphi[2][0]= rnode[1]*(1 - rnode[2]);
    gphi[2][1]= rnode[0]*(1 - rnode[2]);
    gphi[2][2]= -rnode[0]*rnode[1];

    gphi[3][0]= -rnode[1]*(1-rnode[2]);
    gphi[3][1]= (1-rnode[0])*(1-rnode[2]);
    gphi[3][2]= -(1-rnode[0])*rnode[1];

    gphi[4][0]= -(1-rnode[1])*rnode[2];
    gphi[4][1]= -(1-rnode[0])*rnode[2];
    gphi[4][2]= (1-rnode[0])*(1-rnode[1]);

    gphi[5][0]= (1-rnode[1])*rnode[2];
    gphi[5][1]= -rnode[0]*rnode[2];
    gphi[5][2]= rnode[0]*(1-rnode[1]);

    gphi[6][0]= rnode[1]*rnode[2];
    gphi[6][1]= rnode[0]*rnode[2];
    gphi[6][2]= rnode[0]*rnode[1];

    gphi[7][0]= -rnode[1]*rnode[2];
    gphi[7][1]= (1-rnode[0])*rnode[2];
    gphi[7][2]= (1-rnode[0])*rnode[1];
  }

  /*
   * 计算第 i 个单元在参考积分点处的 Jacobi 矩阵 
   *
   *
   * Notes
   * -----
   *  这里的 Jacobi 是一个 3X3 二维数组
   *
   *  x = phi_0 x_0 + phi_1 x_1 + ... + phi_7 x_7
   *
   *  J 要默认初始化 0 矩阵
   *
   *  J[0] = x_u
   *  J[1] = x_v
   *  J[2] = x_w
   */
  template<typename RNode>
  void cell_jacobi_matrix(I i, RNode & rnode, std::vector<Vector> & J)
  {
    J.resize(3);
    std::vector<Vector> gphi;
    cell_grad_shape_function(i, rnode, gphi);
    for(I j = 0; j < 8; j++)
    {
      auto & node = m_node[m_cell[i][j]];
      auto & grad = gphi[j];
      for(I m=0; m < 3; m++) 
        for(I n=0; n < 3; n++)
          J[m][n] += grad[m]*node[n];
    }
  }

  /*
   * 计算第 i 个面所有形函数在参考积分点处的值
   *
   *
   */
  template<typename RNode>
  void face_shape_function(I i, RNode & rnode, std::vector<F> & phi)
  {
    phi.resize(4);
    phi[0] = (1-rnode[0])*(1-rnode[1]);
    phi[1] = rnode[0]*(1 - rnode[1]);
    phi[2] = rnode[0]*rnode[1];
    phi[3] = (1-rnode[0])*rnode[1];
  }

  /*
   * 计算第 i 个面所有形函数在参考积分点处的导数值
   *
   *
   * Notes
   * -----
   *
   */
  template<class RNode>
  void face_grad_shape_function(I i, RNode & rnode, std::vector<Vector_2> & gphi)
  {
    gphi.resize(4);
    gphi[0][0] = -(1-rnode[1]);
    gphi[0][1] = -(1-rnode[0]);

    gphi[1][0] = (1 - rnode[1]);
    gphi[1][1] = -rnode[0];

    gphi[2][0] = rnode[1];
    gphi[2][1] = rnode[0];

    gphi[3][0] = -rnode[1];
    gphi[3][1] = 1-rnode[0];
  }

  /*
   * 计算第 i 个面在参考积分点 rnode 处的 Jacobi 矩阵 
   *
   *
   * Notes
   * -----
   *  这里的 Jacobi 是一个 3X3 二维数组
   *
   *  x = phi_0 x_0 + phi_1 x_1 + phi_2 x_2 + phi_3 x_3
   *
   *  注意，J 要默认初始化 0 矩阵
   *
   *  J[0] = x_u
   *  J[1] = x_v
   */
  template<class RNode>
  void face_jacobi_matrix(I i, RNode & rnode, std::vector<Vector> & J)
  {
    J.resize(2);
    std::vector<Vector_2> gphi;
    face_grad_shape_function(i, rnode, gphi);
    for(I j = 0; j < 4; j++)
    {
      auto & node = m_node[m_face[i][j]];
      auto & grad = gphi[j];
      for(I m = 0; m < 2; m++) 
        for(I n = 0; n < 3; n++)
          J[m][n] += grad[m]*node[n];
    }

  }

  void uniform_refine(const I n=1)
  {
      for(I i=0; i < n; i++)
      {
          auto NN = number_of_nodes();
          auto NE = number_of_edges();
          auto NF = number_of_faces();
          auto NC = number_of_cells();
          m_node.resize(NN + NE + NF + NC);
          for(I j = 0; j < NE; j++)
          {
             edge_barycenter(j, m_node[NN+j]); 
          }
          for(I j = 0; j < NF; j++)
          {
             face_barycenter(j, m_node[NN+NE+j]); 
          }
          for(I j = 0; j < NC; j++)
          {
             cell_barycenter(j, m_node[NN+NE+NF+j]); 
          }

          //单元加密时, 子单元顶点的局部编号, 每个单元的顺序是{n, e, f, e, e, f, c, f}
          m_cell.resize(8*NC);
          for(I j = 0; j < NC; j++)
          { 
              auto c = m_cell[j]; 
              for(I k = 0; k < 8; k++)
              {
                m_cell[k*NC + j][0] = c[k];
                m_cell[k*NC + j][1] = NN + m_cell2edge[j][m_refine[k][1]];
                m_cell[k*NC + j][2] = NN + NE + m_cell2face[j][m_refine[k][2]];
                m_cell[k*NC + j][3] = NN + m_cell2edge[j][m_refine[k][3]];
                m_cell[k*NC + j][4] = NN + m_cell2edge[j][m_refine[k][4]];
                m_cell[k*NC + j][5] = NN + NE + m_cell2face[j][m_refine[k][5]];
                m_cell[k*NC + j][6] = NN + NE + NF + j ;
                m_cell[k*NC + j][7] = NN + NE + m_cell2face[j][m_refine[k][7]];
              }
          }

          m_edge.clear();
          m_face.clear();
          m_cell2edge.clear();
          m_face2cell.clear();
          m_cell2face.clear();
          init_top(); 
      }
  }

  void is_boundary_face(std::vector<bool> & isBdFace)
  {
      auto NF = number_of_faces();
      isBdFace.resize(NF);
      for(int i = 0; i < NF; i++)
      {
          isBdFace[i] = false;
          if(face_to_cell(i)[0]==face_to_cell(i)[1])
          {
              isBdFace[i] = true;
          }
      }
  }

  void is_boundary_node(std::vector<bool> & isBdNode)
  {
      auto NN = number_of_nodes();
      auto NF = number_of_faces();
      isBdNode.resize(NN);
      for(int i = 0; i < NF; i++)
      {
          if(face_to_cell(i)[0]==face_to_cell(i)[1])
          {
              isBdNode[face(i)[0]] = true;
              isBdNode[face(i)[1]] = true;
              isBdNode[face(i)[2]] = true;
              isBdNode[face(i)[3]] = true;
          }
      }
  }

  void cell_to_node(Toplogy & top)
  {
      auto NC = number_of_cells();
      auto NN = number_of_nodes();
      auto nn = number_of_vertices_of_each_cell();

      top.init(3, 0, NC, NN);

      auto & loc = top.locations();
      loc.resize(NC+1, 0);
      auto & nei = top.neighbors();
      nei.resize(NC*nn);

      for(I i = 0; i < NC; i++)
      {
          loc[i+1] = loc[i] + nn;
          nei[nn*i] = m_cell[i][0];
          nei[nn*i + 1] = m_cell[i][1];
          nei[nn*i + 2] = m_cell[i][2];
          nei[nn*i + 3] = m_cell[i][3];
          nei[nn*i + 4] = m_cell[i][4];
          nei[nn*i + 5] = m_cell[i][5];
          nei[nn*i + 6] = m_cell[i][6];
          nei[nn*i + 7] = m_cell[i][7];
      }
  }

  void cell_to_cell(Toplogy & top)
  {

      auto NC = number_of_cells();
      auto NF = number_of_faces();
      top.init(3, 3, NC, NC);
      auto & loc = top.locations();
      loc.resize(NC+1, 0);
      for(I i=0; i < NF; i++)
      {
          loc[m_face2cell[i][0]+1] += 1;
          if(m_face2cell[i][0] != m_face2cell[i][1])
          {
              loc[m_face2cell[i][1]+1] += 1;
          }
      }
      for(I i=0; i < NC; i++)
      {
          loc[i+1] += loc[i];
      }

      auto & nei = top.neighbors();
      nei.resize(loc[NC]);
      std::vector<I> start = loc;
      for(I i = 0; i < NF; i++)
      {
          nei[start[m_face2cell[i][0]]++] = m_face2cell[i][1];
          if(m_face2cell[i][0] != m_face2cell[i][1])
          {
              nei[start[m_face2cell[i][1]]++] = m_face2cell[i][0];
          }
      }
  }

  void node_to_node(Toplogy & top)
  {
      auto NN = number_of_nodes();
      auto NE = number_of_edges();
      top.init(3, 3, NN, NN);
      auto & loc = top.locations();
      loc.resize(NN+1, 0);
      for(I i=0; i < NE; i++)
      {
          loc[m_edge[i][0]+1] += 1;
          loc[m_edge[i][1]+1] += 1;
      }
      for(I i=0; i < NN; i++)
      {
          loc[i+1] += loc[i];
      }

      auto & nei = top.neighbors();
      nei.resize(loc[NN]);
      std::vector<I> start(loc);
      for(I i = 0; i < NE; i++)
      {
          nei[start[m_edge[i][0]]++] = m_edge[i][1];
          nei[start[m_edge[i][1]]++] = m_edge[i][0];
      }
  }

  void node_to_cell(Toplogy & top)
  {
    auto NN = number_of_nodes();
    auto NC = number_of_cells();

    auto & loc = top.locations();
    loc.resize(NN+1, 0);
    for(I i=0; i < NC; i++)
    {
      for(auto v : m_cell[i])
        loc[v+1] += 1;
    }
    for(I i=0; i < NN; i++)
    {
      loc[i+1] += loc[i];
    }

    auto & nei = top.neighbors();
    nei.resize(loc[NN]);
    std::vector<I> start(loc);
    for(I i = 0; i < NC; i++)
    {
      for(auto v : m_cell[i])
        nei[start[v]++] = i;
    }
  }

  void print()
  {
      std::cout << "Nodes:" << std::endl;
      print_entities(m_node);

      std::cout << "Edges:" << std::endl;
      print_entities(m_edge);

      std::cout << "Faces:" << std::endl;
      print_entities(m_face);

      std::cout << "Cells:" << std::endl;
      print_entities(m_cell);

      std::cout << "Face2cell:" << std::endl;
      print_entities(m_face2cell);

      std::cout << "Cell2face:" << std::endl;
      print_entities(m_cell2face);
  }

  template<typename Entities>
  void print_entities(Entities & entities)
  {
      auto N = entities.size();
      for(I i = 0; i < N; i++)
      {
          auto & e = entities[i];
          auto n = e.size();
          std::cout << i << ":";
          for(I j = 0; j < n; j++)
          {
              std::cout << " " << e[j];
          }
          std::cout << std::endl;
      }
  }

private:
    /**
     * \brief 计算第 i 个 cell 的第 j 个 face 的全局唯一整数编号 
     *
     * \param[in] i 单元编号
     * \param[in] j 局部面的编号
     *
     * \warning 返回值可能会发生溢出
     */
    I local_face_index(I i, I j)
    {
        I f[4] = 
        {
            m_cell[i][m_localface[j][0]], 
            m_cell[i][m_localface[j][1]], 
            m_cell[i][m_localface[j][2]],
            m_cell[i][m_localface[j][3]]
        };
        std::sort(f, f+3);
        return  f[0] + f[1]*(f[1]+1)/2 + f[2]*(f[2]+1)*(f[2]+2)/6 + f[2]*(f[2]+1)*(f[2]+2)*(f[3]+3)/24;
    }

    /*
     * \brief 计算第 i 个 cell 的 j 条 edge 全局唯一整数编号
     *
     * \param[in] i 单元编号
     * \param[in] j 局部边的编号
     *
     * \warning 返回值可能会发生溢出
     *
     */
    I local_edge_index(I i, I j)
    {
        I e[2] = {m_cell[i][m_localedge[j][0]], m_cell[i][m_localedge[j][1]]};
        std::sort(e, e+2);
        return  e[0] + e[1]*(e[1]+1)/2;
    }

    /**
     * \brief 计算第 i 个 cell 的 j 个 edge（顶点编号从大到小进行了排序）， 
     *
     * \param[in] i 单元编号
     * \param[in] j 局部边的编号
     *
     * \return 返回一个顶点编号从小到大排序的 Edge 实体
     *
     * \note
     *
     */
    Edge local_sorted_edge(I i, I j)
    {
      Edge e = {m_cell[i][m_localedge[j][0]], m_cell[i][m_localedge[j][1]]};
      std::sort(e.begin(), e.end());
      return e;
    }

    /**
     * \brief 获得第 i 个 cell 的第 j 个 face， 
     *
     * \param[in] i cell 编号
     * \param[in] j 局部 face 的编号
     *
     * \return 返回一个顶点编号从小到大排序的 Face 实体
     *
     * \note
     *
     */
    Face local_sorted_face(I i, I j)
    {
      Face f = 
      {
          m_cell[i][m_localface[j][0]], 
          m_cell[i][m_localface[j][1]], 
          m_cell[i][m_localface[j][2]],
          m_cell[i][m_localface[j][3]]
      };
      std::sort(f.begin(), f.end());
      return  f;
    }

private:
    static int m_localedge[12][2];
    static int m_localface[6][6];
    static int m_localface2edge[6][4];
    static int m_refine[8][8];
    static int m_vtkidx[8];

    std::vector<Node> m_node;
    std::vector<Cell> m_cell; 
    std::vector<Edge> m_edge;
    std::vector<Face> m_face;
    std::vector<Face2cell> m_face2cell;
    std::vector<Cell2face> m_cell2face;
    std::vector<Cell2edge> m_cell2edge;
};

template<typename GK, typename Node, typename Vector>
int HexahedronMesh<GK, Node, Vector>::m_localedge[12][2] = {
    {0, 1}, {1, 2}, {2, 3}, {0, 3},
    {0, 4}, {1, 5}, {2, 6}, {3, 7},
    {4, 5}, {5, 6}, {6, 7}, {4, 7}
};

template<typename GK, typename Node, typename Vector>
int HexahedronMesh<GK, Node, Vector>::m_localface[6][6] = {
    {0, 3, 2, 1}, {4, 5, 6, 7}, // bottom and top faces
    {0, 4, 7, 3}, {1, 2, 6, 5}, // left and right faces  
    {0, 1, 5, 4}, {2, 3, 7, 6}  // front and back faces
};

template<typename GK, typename Node, typename Vector>
int HexahedronMesh<GK, Node, Vector>::m_localface2edge[6][4] = {
    {3,  2, 1, 0}, {8, 9, 10, 11},
    {4, 11, 7, 3}, {1, 6,  9,  5},
    {0,  5, 8, 4}, {2, 7, 10,  6}
};

//单元加密时, 子单元顶点的局部编号, 每个单元的顺序是{n, e, f, e, e, f, c, f}
template<typename GK, typename Node, typename Vector>
int HexahedronMesh<GK, Node, Vector>::m_refine[8][8] = {
    {0, 0, 0, 3, 4, 4, 0, 2}, {1, 1, 0, 0, 5, 3, 0, 4},
    {2, 2, 0, 1, 6, 5, 0, 3}, {3, 3, 0, 2, 7, 2, 0, 5},
    {4, 11, 1, 8, 4, 2, 0, 4}, {5, 8, 1, 9, 5, 4, 0, 3},
    {6, 9, 1, 10, 6, 3, 0, 5}, {7, 10, 1, 11, 7, 5, 0, 2}
};


template<typename GK, typename Node, typename Vector>
int HexahedronMesh<GK, Node, Vector>::m_vtkidx[8] = {0, 1, 2, 3, 4, 5, 6, 7};

} // end of namespace Mesh 

} // end of namespace OF

#endif // end of HexahedronMesh_h
