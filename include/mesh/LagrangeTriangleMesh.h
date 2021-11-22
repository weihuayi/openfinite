
#ifndef LagrangeTriangleMesh_h
#define LagrangeTriangleMesh_h

#include <vector>
#include <array>
#include <map>
#include <math.h>
#include <algorithm>
#include <memory>
#include <iostream>

#include "MeshToplogy.h"
#include "NodeData.h"
#include "common/Multindex.h"
#include "geometry/Vector_2.h"

namespace OF {
namespace Mesh {

/*
 * 
 *
 * Notes
 * -----
 *  三角形网格类, 用 std::vector 做为容器, 用整数数组代表 edge, face, 和 cell
 *  实体, 这里 face 实体即为 edge 实体.
 *
 */
template<typename GK, typename NODE, typename VECTOR, int p>
class LagrangeTriangleMesh 
{
public:
  typedef NODE Node;
  typedef VECTOR Vector;

  typedef typename GK::Vector_2 Vector_2; 

  typedef typename GK::Int I;
  typedef typename GK::Float F;

  typedef typename std::array<I, (p+1)*(p+2)/2> Cell;
  typedef typename std::array<I, p+1> Edge;

  typedef typename std::array<I, 4> Edge2cell;
  typedef typename std::array<I, 3> Cell2edge;

  // 非规则化的拓扑关系， 如共享一个节点的单元个数是不固定的
  // 共享一条边的单元个数也是不固定的
  typedef MeshToplogy<I> Toplogy;

  typedef std::vector<Node> NodeArray;
  typedef std::vector<Edge> EdgeArray;
  typedef std::vector<Cell> CellArray;

  typedef typename NodeArray::iterator NodeIterator;
  typedef typename NodeArray::const_iterator ConstNodeIterator;

  typedef typename EdgeArray::iterator EdgeIterator;
  typedef typename EdgeArray::const_iterator ConstEdgeIterator;

  typedef typename CellArray::iterator CellIterator;
  typedef typename CellArray::const_iterator ConstCellIterator;

  typedef typename OF::Common::MultiIndex MultiIndex;

  static int m_localedgeindex[3];
  static int m_localedge[3][2];
  static int m_vtkidx[10];
  static int m_num[3][3];
  double runtime;

public:
  template<class Mesh>
  LagrangeTriangleMesh(std::shared_ptr<Mesh> linemesh): runtime(0)
  {
    m_holes = 1; // 默认一个外部无界区域的洞
    m_genus = 0;

    m_edge2cell = linemesh->edge_to_cell();
    m_cell2edge = linemesh->cell_to_edge();
    m_NN = linemesh->number_of_nodes();

    construct_mesh_entity(linemesh);
  }

  void insert(const Node & node)
  {
    m_node.push_back(node);
  }

  void insert(const Cell & cell)
  {
    m_cell.push_back(cell);
  }

  int & number_of_holes()
  {
      return m_holes;
  }

  int & number_of_genus()
  {
      return m_genus;
  }

  static int number_of_nodes_of_each_cell()
  {
      return (p+1)*(p+2)/2;
  }

  static int number_of_vertices_of_each_cell()
  {
      return 3;
  }

  static int number_of_edges_of_each_cell()
  {
      return 3;
  }

  static int number_of_nodes_of_each_edge()
  {
      return p+1;
  }

  I number_of_nodes()
  {
      return m_NN;
  }

  I number_of_cells()
  {
      return m_cell.size();
  }

  I number_of_edges()
  {
      return m_edge.size();
  }

  I number_of_faces()
  {
      return m_edge.size();
  }

  I vtk_cell_type(int TD=2)
  {
      if(TD == 2)
          return 69; // VTK_LAGRANGE_TRIANGLE
      else if(TD == 1)//TODO
          return 3; // VTK_LINE
      else
        return 0;// VTK_EMPLTY_CELL = 0
  }

  int* vtk_read_cell_index()
  {
    return m_vtkidx;
  }

  int* vtk_write_cell_index()
  {
    return m_vtkidx;
  }

  int geo_dimension()
  {
      return Node::dimension();
  }

  static int top_dimension()
  {
      return 2;
  }

  /*
   * 
   * Notes
   * -----
   *  TODO: 考虑如何在原来拓扑的基础上重建拓扑信息
   *        原来的边每个变成 2 条, 每个单元内部增加 3 条边
   *        每个单元变成 4 个单元, 这样会提高程序的效率吗?
   */
  void update_top()
  {
      return;
  }

  /**
   * \brief 生成网格的节点, 边, 单元
   */
  template<class Mesh>
  void construct_mesh_entity(std::shared_ptr<Mesh> linemesh)
  {
    const auto & node = linemesh->nodes();
    const auto & edge = linemesh->edges();
    const auto & cell = linemesh->cells();

    int NN = linemesh->number_of_nodes();
    int NE = linemesh->number_of_edges();
    int NC = linemesh->number_of_cells();
    int NNE = p-1;
    int NNC = (p-1)*(p-2)/2;

    int NN0 = NN + NE*NNE + NC*NNC; /**< 高阶网格的节点个数 */
    m_node.resize(NN0);
    m_edge.resize(NE);
    m_cell.resize(NC);

    // 构建 node
    std::copy(node.begin(), node.end(), m_node.begin());

    //构建edge
    for(int i = 0; i < NE; i++)
    {
      m_edge[i][0] = edge[i][0];
      for(int j = 1; j < p; j++)
      {
        m_edge[i][j] = NN+i*NNE+j-1;
        m_node[m_edge[i][j]] = (j*node[edge[i][1]] + (p-j)*node[edge[i][0]])/p;
      }
      m_edge[i][p] = edge[i][1];
    }

    //构建cell
    std::vector<MultiIndex> mindexs;
    MultiIndex::gen_all_index(p, 3, mindexs);

    std::vector<std::array<int, 3>> dpoint;
    get_domain_point(dpoint);
    for(int i = 0; i < NC; i++)
    {
      auto & idx = m_localedgeindex; /**< 网格局部边的编号与递增序列编号的索引 TODO */
      bool eflag[3]; /**< 三条边的方向是否正确 */
      eflag[0] = edge[m_cell2edge[i][idx[0]]][0] == cell[i][0];
      eflag[1] = edge[m_cell2edge[i][idx[1]]][0] == cell[i][0];
      eflag[2] = edge[m_cell2edge[i][idx[2]]][0] == cell[i][1];

      for(int j = 0; j < (int)mindexs.size(); j++)/**< 循环所有域点 */
      {
        if(dpoint[j][0]==0) /**< 0维 域点 */
        {
          m_cell[i][j] = cell[i][dpoint[j][1]];
        }
        else if(dpoint[j][0]==1) /**< 1维 域点 */
        {
          if(eflag[dpoint[j][1]])
            m_cell[i][j] = m_edge[m_cell2edge[i][idx[dpoint[j][1]]]][dpoint[j][2]+1];
          else
            m_cell[i][j] = m_edge[m_cell2edge[i][idx[dpoint[j][1]]]][p-dpoint[j][2]-1];
        }
        else /**< 2维 域点 */
        {
          m_cell[i][j] = NN + NE*NNE + i*NNC + dpoint[j][2];
          m_node[m_cell[i][j]] = (mindexs[j][0]*node[cell[i][0]] + 
              mindexs[j][1]*node[cell[i][1]] + mindexs[j][2]*node[cell[i][2]])/p;
        }
      }
    }
  }

  /**
   * \brief 获取单元的域点的信息 [维数, 子单形编号, 域点在子单形中的编号]
   */
  void get_domain_point(std::vector<std::array<int, 3>> & domain)
  {
    domain.resize(number_of_nodes_of_each_cell());
    domain[0] = {0, 0, 0};

    int num = 1, num0 = 0;
    for(int i = 0; i < p-1; i++)
    {
      domain[num++] = {1, 0, i};
      for(int j = 0; j < i; j++)
      {
        domain[num++] = {2, 0, num0++}; 
      }
      domain[num++] = {1, 1, i};
    }
    domain[num++] = {0, 1, 0};
    for(int i = 0; i < p-1; i++)
    {
      domain[num++] = {1, 2, i};
    }
    domain[num++] = {0, 2, 0};
  }

  template<class Container>
  void BMatrix(const Container & bc, std::vector<std::array<double, 3>> & B)
  {
    clock_t s = clock();
    /** 计算 B */
    B[0] = {1, 1, 1};
    for(int i = 0; i < p; i++)/**< 赋值 */
    {
      for(int j = 0; j < 3; j++)
        B[i+1][j] = p*bc[j] - i;
    }

    for(int i = 0; i < p; i++)/**< 连乘 */
    {
      for(int j = 0; j < 3; j++)
        B[i+1][j] *= B[i][j];
    }

    int P = 1;
    for(int i = 0; i < p; i++)/**< 乘 P */
    {
      P *= i+1;
      for(int j = 0; j < 3; j++)
        B[i+1][j] /= P;
    }
    clock_t e = clock();
    runtime += (double)(e-s)/CLOCKS_PER_SEC;
  }

  template<class Container>
  void nabla_BMatrix(const Container & bc, std::vector<std::array<Vector_2, 3>> & NB)
  {
    /** 计算 D^i */
    std::vector<std::array<double, p>> D0(p);
    std::vector<std::array<double, p>> D1(p);
    std::vector<std::array<double, p>> D2(p);
    for(int i = 0; i < p; i++)/**< 赋值 */
    {
      for(int j = 0; j < p; j++)
      {
        D0[i][j] = p*bc[0]-i;
        D1[i][j] = p*bc[1]-i;
        D2[i][j] = p*bc[2]-i;
      }
      D0[i][i] = p;
      D1[i][i] = p;
      D2[i][i] = p;
    }

    for(int i = 1; i < p; i++)/**< 连乘 */
    {
      for(int j = 0; j < p; j++)
      {
        D0[i][j] *= D0[i-1][j];
        D1[i][j] *= D1[i-1][j];
        D2[i][j] *= D2[i-1][j];
      }
    }

    /** 计算 D */
    std::vector<std::array<double, 3> > D(p+1, {0, 0, 0});
    for(int i = 0; i < p; i++)
    {
      for(int j = 0; j < i+1; j++)
      {
        D[i+1][0] += D0[i][j];
        D[i+1][1] += D1[i][j];
        D[i+1][2] += D2[i][j];
      }
    }

    Vector_2 GLambda[3];
    GLambda[0] = Vector_2{-1.0, -1.0};
    GLambda[1] = Vector_2{1.0, 0.0};
    GLambda[2] = Vector_2{0.0, 1.0};

    int P = 1;
    for(int i = 0; i < p; i++)/**< 乘 P */
    {
      P *= i+1;
      for(int j = 0; j < 3; j++)
      {
        NB[i+1][j] = D[i+1][j]*GLambda[j]/P;
      }
    }
  }

  /**
   * \brief 计算 Lagrange 形函数在重心坐标为 bc 的点处的值.
   * 前三个形函数是重心坐标函数其数值就是 bc[0], bc[1], bc[2].
   */
  template<class Container>
  void shape_function(const Container & bc, std::vector<double> & val)
  {
    /** 计算 B */
    std::vector<std::array<double, 3>> B(p+1, {1, 1, 1});
    BMatrix(bc, B);

    /** 计算形函数值 val */
    MultiIndex midx(p, 3);
    val.resize((p+1)*(p+2)/2);
    for(int i = 0; i < (int)val.size(); i++)
    {
      val[i] = B[midx[0]][0]*B[midx[1]][1]*B[midx[2]][2];
      midx.nex_multi_index();
    }
  }

  /**
   * \brief 计算 Lagrange 形函数在重心坐标为 bc 的点处的值.
   * 前三个形函数是重心坐标函数其数值就是 bc[0], bc[1], bc[2].
   */
  template<class Container>
  void grad_shape_function(const Container & bc, std::vector<Vector_2> & val)
  {
    /** 计算 B */
    std::vector<std::array<double, 3>> B(p+1, {1, 1, 1});
    BMatrix(bc, B);

    /** 计算 \nabla B */
    std::vector<std::array<Vector_2, 3>> NB(p+1);
    nabla_BMatrix(bc, NB);

    /** 计算形函数梯度值 val */
    MultiIndex midx(p, 3);
    val.resize((p+1)*(p+2)/2);
    for(int i = 0; i < (int)val.size(); i++)
    {
      val[i] = B[midx[1]][1]*B[midx[2]][2]*NB[midx[0]][0] + 
        B[midx[0]][0]*B[midx[2]][2]*NB[midx[1]][1] + 
        B[midx[0]][0]*B[midx[1]][1]*NB[midx[2]][2];
      midx.nex_multi_index();
    }
  }

  template<class Container0, class Container1>
  void jacobi_matrix(int c, const Container0 & bc,  Container1 & val)
  {
    std::vector<Vector_2> gsf;
    grad_shape_function(bc, gsf);
    val.resize(geo_dimension(), {0, 0});
    for(int i = 0; i < number_of_nodes_of_each_cell(); i++)
    {
      for(int j = 0; j < geo_dimension(); j++)
      {
        for(int k = 0; k < 2; k++)
        {
          val[j][k] += m_node[m_cell[c][i]][j] * gsf[i][k];
        }
      }
    }
  }

  Node bc_to_point(int i, const std::array<double, 3> & bc)
  {
    std::vector<double> sf;
    shape_function(bc, sf);
    auto & c = cell(i);
    Node point = sf[0]*m_node[c[0]];
    for(int i = 1; i < number_of_nodes_of_each_cell(); i++)
    {
      point += sf[i]*m_node[c[i]];
    }
    return point;
  }

  Node bc_to_edge_point(); //TODO
  
  void push_dual_vector(int c, const std::array<double, 3> & bc, std::vector<Vector_2> & v0, std::vector<Vector> & v1)
  {
    auto GD = geo_dimension();
    std::vector<std::array<double, 2>> J(GD, {0, 0});
    jacobi_matrix(c, bc, J);

    v1.resize(v0.size());
    if(GD == 2)/**< 推前矩阵为 (J^T)^{-1} */
    {
      double det_J = J[0][0]*J[1][1] - J[0][1]*J[1][0];
      for(int i = 0; i < (int)v0.size(); i++)
      {
        v1[i][0] = (J[1][1] * v0[i][0] - J[1][0] * v0[i][1])/det_J;
        v1[i][1] = (-J[0][1] * v0[i][0] + J[0][0] * v0[i][1])/det_J;
      }
    }
    else/**< 推前矩阵为 J(J^T J)^{-1} */
    {
      double G[3] = {0, 0, 0};
      for(int i = 0; i < GD; i++)
      {
        G[0] += J[i][0]*J[i][0];
        G[1] += J[i][1]*J[i][0];
        G[2] += J[i][1]*J[i][1];
      }
      double det_G = G[0]*G[2] - G[1]*G[1];
      for(int i = 0; i < GD; i++)
      {
        double a = (J[i][0]*G[2] - J[i][1]*G[1])/det_G;
        double b = (-J[i][0]*G[1] + J[i][1]*G[0])/det_G;
        J[i][0] = a;
        J[i][1] = b;
      }/**< 获得推前矩阵 */
      for(int i = 0; i < (int)v0.size(); i++)
      {
        for(int j = 0; j < GD; j++)
        {
          v1[i][j] = (J[j][0] * v0[i][0] + J[j][1] * v0[i][1]);
        }
      }
    }
  }

  void push_vector(int c, const std::array<double, 3> & bc, std::vector<Vector_2> & v0, std::vector<Vector> & v1)
  {
    auto J = jacobi_matrix(c, bc);
    int GD = geo_dimension();
    v1.resize(v0.size());
    for(int i = 0; i < (int)v0.size(); i++)
    {
      for(int j = 0; j < GD; j++)
      {
        v1[i][j] = J[j][0]*v0[i][0] + J[j][1]*v0[i][0];
      }
    }
  }

  template<class Container0>
  double det_jacobi(int c, const Container0 & bc)
  {
    std::vector<std::array<double, 2>> J; 
    jacobi_matrix(c, bc, J);
    double det[3] = {0.0, 0.0, 0.0};
    for(int k = 0; k < (int)J.size(); k++)
    {
      det[0] += J[k][0]*J[k][0];
      det[1] += J[k][0]*J[k][1];
      det[2] += J[k][1]*J[k][1];
    }
    return std::sqrt(det[0]*det[2] - det[1]*det[1]);
  }

  bool is_boundary_edge(int i)
  {
    return edge_to_cell(i)[0] == edge_to_cell(i)[1];
  }

  void is_boundary_edge(std::vector<bool> & isBdEdge)
  {
    auto NE = number_of_edges();
    isBdEdge.resize(NE, false);
    for(int i = 0; i < NE; i++)
    {
      isBdEdge[i] = edge_to_cell(i)[0]==edge_to_cell(i)[1];
    }
  }

  void is_boundary_node(std::vector<bool> & isBdNode)
  {
    auto NN = number_of_nodes();
    auto NE = number_of_edges();
    isBdNode.resize(NN);
    for(int i = 0; i < NE; i++)
    {
      if(edge_to_cell(i)[0]==edge_to_cell(i)[1])
      {
        isBdNode[edge(i)[0]] = true;
        isBdNode[edge(i)[1]] = true;
      }
    }
  }

  void cell_to_vertex(int i, int j)
  {
    int idx[3] = {0, (p+1)*(p+2)/2 - p -1, (p+1)*(p+2)/2 -1};
    return m_cell[i][idx[j]];
  }

  void edge_to_vertex(int i, int j)
  {
    int idx[2] = {0, p};
    return m_edge[i][idx[j]];
  }

  void cell_to_vertex(int i, std::array<int, 3> & c2v)
  {
    c2v[0] = m_cell[i][0];
    c2v[1] = m_cell[i][(p+1)*(p+2)/2 - p - 1];
    c2v[2] = m_cell[i][(p+1)*(p+2)/2 - 1];
  }

  void edge_to_vertex(int i, std::array<int, 2> & e2v)
  {
    e2v[0] = m_edge[i][0];
    e2v[1] = m_edge[i][p];
  }

  /**
   * \brief 判断单元 i 的边的局部方向是否与全局方向相同
   */
  void cell_to_edge_sign(int i, std::array<double, 3> & flag)
  {
    int NCN = number_of_nodes_of_each_cell();
    auto & c2e = cell_to_edge(i);
    flag[0] = m_edge[c2e[0]][0] == m_cell[NCN-p-1]; 
    flag[1] = m_edge[c2e[1]][0] == m_cell[NCN-1]; 
    flag[2] = m_edge[c2e[2]][0] == m_cell[0]; 
  }

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

  CellIterator cell_begin()
  {
    return m_cell.begin();
  }

  CellIterator cell_end()
  {
    return m_cell.end();
  }

  std::vector<Node> & nodes()
  {
    return m_node;
  }

  std::vector<Cell> & cells()
  {
    return m_cell;
  }

  std::vector<Edge> & edges()
  {
    return m_edge;
  }

  Node & node(const I i)
  {
    return m_node[i];
  }

  Cell & cell(const I i)
  {
    return m_cell[i];
  }

  Edge & edge(const I i)
  {
    return m_edge[i];
  }

  Edge2cell & edge_to_cell(const I i)
  {
    return m_edge2cell[i];
  }

  Cell2edge & cell_to_edge(const I i)
  {
    return m_cell2edge[i];
  }

  void print()
  {
      std::cout << "Nodes:" << std::endl;
      print_entities(m_node);

      std::cout << "Edges:" << std::endl;
      print_entities(m_edge);

      std::cout << "Cells:" << std::endl;
      print_entities(m_cell);

      std::cout << "Edge2cell:" << std::endl;
      print_entities(m_edge2cell);

      std::cout << "Cell2edge:" << std::endl;
      print_entities(m_cell2edge);
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
  /*
   *
   * Notes
   *  计算第 i 个单元的第 j 条边全局唯一的一个整数索引
   */
  I local_edge_index(I i, I j)
  {
      I e[2] = {m_cell[i][m_localedge[j][0]], m_cell[i][m_localedge[j][1]]};
      std::sort(e, e+2);
      return  e[0] + e[1]*(e[1]+1)/2;
  }

  Edge local_edge_index0(I i, I j)
  {
      I e[2] = {m_cell[i][m_localedge[j][0]], m_cell[i][m_localedge[j][1]]};
      std::sort(e, e+2);
      return Edge({e[0], e[1]});
  }

private:
  int m_NN;
  int m_holes; // 网格中洞的个数
  int m_genus; // 网格表示曲面的亏格数
  std::vector<Node> m_node;
  std::vector<Edge> m_edge;
  std::vector<Cell> m_cell; 
  std::vector<Edge2cell> m_edge2cell;
  std::vector<Cell2edge> m_cell2edge;
};


template<typename GK, typename NODE, typename VECTOR, int p>
int LagrangeTriangleMesh<GK, NODE, VECTOR, p>::m_localedgeindex[3] = {
    2, 1, 0
};

template<typename GK, typename NODE, typename VECTOR, int p>
int LagrangeTriangleMesh<GK, NODE, VECTOR, p>::m_localedge[3][2] = {
    {0, 1}, {2, 0}, {1, 2}
};

template<typename GK, typename NODE, typename VECTOR, int p>
int LagrangeTriangleMesh<GK, NODE, VECTOR, p>::m_vtkidx[10] = {
  9, 0, 6, 5, 2, 1, 3, 7, 8, 4
};

template<typename GK, typename NODE, typename VECTOR, int p>
int LagrangeTriangleMesh<GK, NODE, VECTOR, p>::m_num[3][3] = {
  {0, 1, 2}, {1, 2, 0}, {2, 0, 1}
};

} // end of namespace Mesh 

} // end of namespace OF

#endif // end of LagrangeTriangleMesh_h
