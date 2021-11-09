
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

  static int m_localedge[3][2];
  static int m_vtkidx[3];
  static int m_num[3][3];

public:
  template<class Mesh>
  LagrangeTriangleMesh(std::shared_ptr<Mesh> linemesh)
  {
    m_holes = 1; // 默认一个外部无界区域的洞
    m_genus = 0;

    m_edge2cell = linemesh->edge_to_cell();
    m_cell2edge = linemesh->cell_to_edge();

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
      return m_node.size();
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
          return 5; // VTK_TRIANGLE
      else if(TD == 1)
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
        m_edge[i][j] = NN+i*NE+j;
        m_node[m_edge[i][j]] = (j*node[edge[i][0]] + (p-j)*node[edge[i][1]])/p;
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
      bool eflag[3]; /**< 三条边的方向是否正确 */
      eflag[0] = m_edge2cell[m_cell2edge[i][0]][0]==i;
      eflag[1] = m_edge2cell[m_cell2edge[i][1]][0]==i;
      eflag[2] = m_edge2cell[m_cell2edge[i][2]][0]==i;

      for(int j = 0; j < (int)mindexs.size(); j++)/**< 循环所有域点 */
      {
        if(dpoint[j][0]==0) /**< 0维 域点 */
        {
          m_cell[i][j] = cell[i][dpoint[j][1]];
        }
        else if(dpoint[j][0]==1) /**< 1维 域点 */
        {
          if(eflag[dpoint[j][1]])
            m_cell[i][j] = m_edge[m_cell2edge[i][dpoint[j][1]]][dpoint[j][2]+1];
          else
            m_cell[i][j] = m_edge[m_cell2edge[i][dpoint[j][1]]][p-dpoint[j][2]-1];
        }
        else /**< 2维 域点 */
        {
          m_cell[i][j] = NN + NE*NNE + i*NNC + dpoint[j][2];
          m_node[m_cell[i][j]] = ((double)mindexs[j][0]*node[cell[i][0]] + 
              (double)mindexs[j][1]*node[cell[i][1]]+
              (double)mindexs[j][2]*node[cell[i][2]])/(double)p;
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

  void is_boundary_edge(std::vector<bool> & isBdEdge)
  {
      auto NE = number_of_edges();
      isBdEdge.resize(NE);
      for(int i = 0; i < NE; i++)
      {
          isBdEdge[i] = false;
          if(edge_to_cell(i)[0]==edge_to_cell(i)[1])
              isBdEdge[i] = true;
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


  void cell_to_cell(Toplogy & top)
  {
      auto NC = number_of_cells();
      auto NE = number_of_edges();
      auto & loc = top.locations();
      loc.resize(NC+1, 0);
      for(I i=0; i < NE; i++)
      {
          loc[m_edge2cell[i][0]+1] += 1;
          if(m_edge2cell[i][0] != m_edge2cell[i][1])
          {
              loc[m_edge2cell[i][1]+1] += 1;
          }
      }
      for(I i=0; i < NC; i++)
      {
          loc[i+1] += loc[i];
      }

      auto & nei = top.neighbors();
      nei.resize(loc[NC]);
      std::vector<I> start = loc;
      for(I i = 0; i < NE; i++)
      {
          nei[start[m_edge2cell[i][0]]++] = m_edge2cell[i][1];
          if(m_edge2cell[i][0] != m_edge2cell[i][1])
          {
              nei[start[m_edge2cell[i][1]]++] = m_edge2cell[i][0];
          }
      }
  }

  void node_to_node(Toplogy & top)
  {
    auto NN = number_of_nodes();
    auto NE = number_of_edges();
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
    auto & locid = top.local_indices();
    nei.resize(loc[NN]);
    locid.resize(loc[NN]);
    std::vector<I> start(loc);
    for(I i = 0; i < NC; i++)
    {
      for(int j = 0; j < 3; j++)
      {
        auto v = m_cell[i][j];
        locid[start[v]] = j;
        nei[start[v]++] = i;
      }
    }
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

  NodeIntData & get_node_int_data()
  {
    return m_nodeintdata;
  }

  NodeDoubleData & get_node_double_data()
  {
    return m_nodedoubledata;
  }
 
  F cell_measure(const I i)
  {//TODO: 需要考虑 inline 函数吗?
    auto & c = m_cell[i];
    auto v1 = m_node[c[1]] - m_node[c[0]];
    auto v2 = m_node[c[2]] - m_node[c[0]];
    return 0.5*cross(v1, v2);
  }

  void cell_measure(std::vector<F> & measure)
  {
    auto NC = number_of_cells();
    measure.resize(NC);
    for(I i = 0; i < NC; i++)
        measure[i] = cell_measure(i);
  }

  void cell_size(std::vector<F> & cellsize)
  {
    auto NC = number_of_cells();
    cellsize.resize(NC);
    for(I i = 0; i < NC; i++)
        cellsize[i] = std::sqrt(cell_measure(i));
  }

  F cell_size(I cidx)
  {
    return std::sqrt(cell_measure(cidx));
  }

  Node edge_barycenter(const I i)
  {
    auto & e = m_edge[i];
    F x = (m_node[e[0]][0] + m_node[e[1]][0])/2.0;
    F y = (m_node[e[0]][1] + m_node[e[1]][1])/2.0;
    return Node(x, y);
  }

  void edge_barycenter(const I i, Node & node)
  {
    auto & e = m_edge[i];
    for(int i = 0; i < geo_dimension(); i++)
      node[i] = (m_node[e[0]][i] + m_node[e[1]][i])/2.0;
  }

  void edge_barycenter(std::vector<Node> & edgebarycenter)
  {
    int NE = number_of_edges();
    edgebarycenter.resize(NE);
    for(int i = 0; i < NE; i++)
      edge_barycenter(i, edgebarycenter[i]);
  }


  void cell_barycenter(const I i, Node & node)
  {
    auto & c = m_cell[i];
    for(int i = 0; i < geo_dimension(); i++)
      node[i] = (m_node[c[0]][i] + m_node[c[1]][i] + m_node[c[2]][i])/3.0;
  }

  void cell_barycenter(std::vector<Node> & cellbarycenter)
  {
    int NC = number_of_cells();
    cellbarycenter.resize(NC);
    for(int i = 0; i < NC; i++)
      cell_barycenter(i, cellbarycenter[i]);
  }

  Vector edge_normal(const I & i)
  {
    auto v = edge_tangent(i);
    return Vector(v[1], -v[0]);
  }

  Vector edge_tangent(const I & i)
  {
    auto e = m_edge[i];
    auto v = m_node[e[1]] - m_node[e[0]];
    return v/std::sqrt(v.squared_length());
  }

  F edge_measure(const I i)
  {
    auto & e = m_edge[i];
    auto v = m_node[e[1]] - m_node[e[0]];
    return std::sqrt(v.squared_length());
  }

  void edge_measure(std::vector<F> & measure)
  {
    auto NE = number_of_edges();
    measure.resize(NE);
    for(I i = 0; i < NE; i++)
      measure[i] = edge_measure(i);
  }

  F face_measure(const I i)
  {
    auto & e = m_edge[i];
    auto v = m_node[e[1]] - m_node[e[0]];
    return std::sqrt(v.squared_length());
  }

  void face_measure(std::vector<F> & measure)
  {
    auto NE = number_of_edges();
    measure.resize(NE);
    for(I i = 0; i < NE; i++)
      measure[i] = edge_measure(i);
  }

  void cell_to_node(Toplogy & top)
  {
      auto NC = number_of_cells();
      auto NN = number_of_nodes();
      auto nn = number_of_vertices_of_each_cell();

      top.init(2, 0, NC, NN);

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
      }
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
  int m_holes; // 网格中洞的个数
  int m_genus; // 网格表示曲面的亏格数
  std::vector<Node> m_node;
  std::vector<Edge> m_edge;
  std::vector<Cell> m_cell; 
  std::vector<Edge2cell> m_edge2cell;
  std::vector<Cell2edge> m_cell2edge;
  NodeIntData m_nodeintdata;
  NodeDoubleData m_nodedoubledata;
};


template<typename GK, typename NODE, typename VECTOR, int p>
int LagrangeTriangleMesh<GK, NODE, VECTOR, p>::m_localedge[3][2] = {
    {0, 1}, {2, 0}, {1, 2}
};

template<typename GK, typename NODE, typename VECTOR, int p>
int LagrangeTriangleMesh<GK, NODE, VECTOR, p>::m_vtkidx[3] = {0, 1, 2};

template<typename GK, typename NODE, typename VECTOR, int p>
int LagrangeTriangleMesh<GK, NODE, VECTOR, p>::m_num[3][3] = {
  {0, 1, 2}, {1, 2, 0}, {2, 0, 1}
};

} // end of namespace Mesh 

} // end of namespace OF

#endif // end of LagrangeTriangleMesh_h
