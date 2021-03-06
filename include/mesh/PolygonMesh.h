#ifndef PolygonMesh_h
#define PolygonMesh_h

#include "MeshToplogy.h"

namespace OF {
namespace Mesh {

/*
 * 
 *
 * Notes
 * -----
 *  多边形网格类, 用 std::vector 做为容器, 用整数数组代表 edge, face, 和 cell
 *  实体, 这里 face 实体即为 edge 实体.
 *
 */
template<typename GK, typename NODE, typename VECTOR>
class PolygonMesh 
{
public:
  typedef NODE Node;
  typedef VECTOR Vector;

  typedef typename GK::Int I;
  typedef typename GK::Float F;

  // 非规则化的拓扑关系， 如共享一个节点的单元个数是不固定的
  // 共享一条边的单元个数也是不固定的
  typedef MeshToplogy<I> Toplogy;

  typedef Toplogy CellArray;
  typedef typename CellArray::AdjEntitiesIterator CellIterator;
  typedef CellIterator Cell;

  typedef typename std::array<I, 2> Edge;
  typedef Edge Face;
  typedef std::vector<Node> NodeArray;
  typedef std::vector<Edge> EdgeArray;
  typedef std::vector<Face> FaceArray;

  typedef typename NodeArray::iterator NodeIterator;
  typedef typename EdgeArray::iterator EdgeIterator;
  typedef typename FaceArray::iterator FaceIterator;

  typedef typename std::array<I, 4> Edge2cell;
  typedef typename std::array<I, 4> Face2cell;

public:
  PolygonMesh()
  {
      m_holes = 1; // 默认一个外部无界区域的洞
      m_genus = 0;
  }

  void insert(const Node & node)
  {
      m_node.push_back(node);
  }


  int & number_of_holes()
  {
      return m_holes;
  }

  int & number_of_genus()
  {
      return m_genus;
  }

  int number_of_nodes_of_each_cell(const I i)
  {
      return m_cell.number_of_adj_entities(i);
  }

  int number_of_vertices_of_each_cell(const I i)
  {
      return m_cell.number_of_adj_entities(i);
  }

  int number_of_edges_of_each_cell(const I i)
  {
      return m_cell.number_of_adj_entities(i);
  }

  int number_of_faces_of_each_cell(const I i)
  {
      return m_cell.number_of_adj_entities(i);
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
          return 7; // VTK_POLYGON
      else if(TD == 1)
          return 3; // VTK_LINE
      else
        return 0;// VTK_EMPLTY_CELL = 0
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
   */
  void update_top()
  {
      return;
  }

  void init_top()
  {
      auto NN = number_of_nodes();
      auto NC = number_of_cells();
      m_cell2edge.resize(m_cell.size());
      m_edge2cell.clear();
      // 在知道网格代表曲面的洞和亏格的个数后, 可以准确计算边的个数
      m_edge2cell.reserve(NN + NC + m_holes + 2*m_genus - 2);
      std::map<I, I> idxmap;

      I NE = 0;
      // 偏历所有单元
      for(I i = 0; i < m_cell.size(); i++)
      {
          for(I j = 0; j < 4; j++)
          {
             auto s = local_edge_index(i, j);
             auto it = idxmap.find(s);
             if(it == idxmap.end())
             {
                m_cell2edge[i][j] = NE;
                idxmap.insert(std::pair<I, I>(s, NE));
                m_edge2cell.push_back(Edge2cell{i, i, j, j});
                NE++;
             }
             else
             {
                m_cell2edge[i][j] = it->second;
                m_edge2cell[it->second][1] = i;
                m_edge2cell[it->second][3] = j;
             }
          }
      }

      m_edge.resize(NE);
      for(I i = 0; i < NE; i++)
      {
          auto & c = m_cell[m_edge2cell[i][0]];
          auto j = m_edge2cell[i][2];
          m_edge[i][0] = c[m_localedge[j][0]];
          m_edge[i][1] = c[m_localedge[j][1]];
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

  Node_iterator node_begin()
  {
    return m_node.begin();
  }

  Node_iterator node_end()
  {
    return m_node.end();
  }

  Edge_iterator edge_begin()
  {
    return m_edge.begin();
  }

  Edge_iterator edge_end()
  {
    return m_edge.end();
  }

  Face_iterator face_begin()
  {
    return m_edge.begin();
  }

  Face_iterator face_end()
  {
    return m_edge.end();
  }

  Cell_iterator cell_begin()
  {
    return m_cell.begin();
  }

  Cell_iterator cell_end()
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

  std::vector<Face> & faces()
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

  F cell_measure(const I i)
  {
    F a=0.0;
    for(I j=0; j < 4; j++)
    {
      auto & p0 = node(m_cell[i][m_localedge[j][0]]);
      auto & p1 = node(m_cell[i][m_localedge[j][1]]);
      auto v = p1 - p0;
      a += p0[0]*v[1];
      a -= p0[1]*v[0];
    }
    return 0.5*a;
  }

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
    return Node(x, y);
  }

  void edge_barycenter(const I i, Node & node)
  {
    auto & e = m_edge[i];
    node[0] = (m_node[e[0]][0] + m_node[e[1]][0])/2.0;
    node[1] = (m_node[e[0]][1] + m_node[e[1]][1])/2.0;
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

  void uniform_refine(const int n=1)
  {
      for(I i=0; i < n; i++)
      {
          auto NN = number_of_nodes();
          auto NE = number_of_edges();
          std::cout << NN << " " << NE << std::endl;
          m_node.resize(NN + NE);
          std::cout << m_node.size() << std::endl;
          for(I j = 0; j < NE; j++)
          {
             edge_barycenter(j, m_node[NN+j]); 
          }
          auto NC = number_of_cells();
          m_cell.resize(4*NC);
          for(I j = 0; j < NC; j++)
          { //TODO: 考虑不同的排列顺序是否影响程序的效率
              auto c = m_cell[j]; 
              m_cell[j][0] = c[0];
              m_cell[j][1] = m_cell2edge[j][2] + NN;
              m_cell[j][2] = m_cell2edge[j][1] + NN;

              m_cell[NC + j][0] = c[1];
              m_cell[NC + j][1] = m_cell2edge[j][0] + NN;
              m_cell[NC + j][2] = m_cell2edge[j][2] + NN;

              m_cell[2*NC + j][0] = c[2];
              m_cell[2*NC + j][1] = m_cell2edge[j][1] + NN;
              m_cell[2*NC + j][2] = m_cell2edge[j][0] + NN;

              m_cell[3*NC + j][0] = m_cell2edge[j][0] + NN; 
              m_cell[3*NC + j][1] = m_cell2edge[j][1] + NN;
              m_cell[3*NC + j][2] = m_cell2edge[j][2] + NN;
          }
          m_edge.clear();
          m_cell2edge.clear();
          m_edge2cell.clear();
          init_top(); 
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
    auto s = m_cell[i].size();
    auto j0 = j;
    auto j1 = (j+1)%s;
    I e[2] = {m_cell[i][j0], m_cell[i][j1]};
    std::sort(e, e+2);
    return  e[0] + e[1]*(e[1]+1)/2;
  }

private:
  int m_holes; // 网格中洞的个数
  int m_genus; // 网格表示曲面的亏格数
  CellArray m_cell; 
  NodeArray m_node;
  EdgeArray m_edge;
  std::vector<Edge2cell> m_edge2cell;
};

} // end of namespace Mesh

} // end of namespace OF
#endif // end of PolygonMesh_h
