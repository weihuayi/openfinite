
/**
 * \file TetrahedronMesh.h
 * \author Huayi Wei, Chunyu Chen
 *
 * \date 2022.02.15
 *
 * \brief 四面体网格
 *
 */
#ifndef TetrahedronMesh_h
#define TetrahedronMesh_h

#include <vector>
#include <array>
#include <map>
#include <math.h>
#include <iostream>
#include <algorithm>

#include "mesh/MeshToplogy.h"
#include "RArray2d.h"
#include "femcore.h"

namespace OF {
namespace Mesh {

/*
 * 
 *
 * Notes
 * -----
 *  四面体网格类, 用 std::vector 做为容器, 用整数数组代表 edge, face, 和 cell
 *  实体.
 *
 */
template<typename GK>
class TetrahedronMesh 
{
public:
  typedef typename GK::Int I;
  typedef typename GK::Float F;

  // 非规则化的拓扑关系， 如共享一个节点的单元个数是不固定的
  // 共享一条边的单元个数也是不固定的
  typedef MeshToplogy<I, std::vector<I> > Toplogy;

  typedef typename RArray2d<I>::iterator CellIterator;
  typedef typename RArray2d<I>::iterator FaceIterator;
  typedef typename RArray2d<I>::iterator EdgeIterator;
  typedef typename RArray2d<F>::iterator NodeIterator;

  static int m_localedge[6][2];
  static int m_localface[4][3];
  static int m_localface2edge[4][3];
  static int m_refine[3][6];
  static int m_index[12][4];
  static int m_vtkidx[4];
  static int m_num[4][4];

public:

  TetrahedronMesh(){}

  void insert(const std::initializer_list<F> & node)
  {
    m_node.push_back(node);
  }

  void insert(const std::initializer_list<I> & cell)
  {
    m_cell.push_back(cell);
  }

  I number_of_nodes()
  {
    return m_node.shape(0);
  }

  I number_of_cells()
  {
    return m_cell.shape(0);
  }

  I number_of_edges()
  {
    return m_edge.shape(0);
  }

  I number_of_faces()
  {
    return m_face.shape(0);
  }

  static I number_of_nodes_of_each_cell()
  {
      return 4;
  }

  static I number_of_vertices_of_each_cell()
  {
      return 4;
  }

  static I geo_dimension()
  {
      return 3;
  }

  static I top_dimension()
  {
      return 3;
  }

  int* vtk_read_cell_index()
  {
    return m_vtkidx;
  }

  int* vtk_write_cell_index()
  {
    return m_vtkidx;
  }

  I vtk_cell_type(I TD=3)
  {
      if(TD == 3)
          return 10; // VTK_TETRA = 10
      else if(TD == 2)
          return 5; // VTK_TRIANGLE = 5
      else if(TD == 1)
          return 3; // VTK_LINE = 1
      else
          return 0; // VTK_EMPLTY_CELL = 0
  }

  void init_top()
  {
    auto NN = number_of_nodes();
    auto NC = number_of_cells();
    RArray2d<I> allFace(NC*4, 3);
    RArray2d<I> allEdge(NC*6, 2); /**< 可以由 allFace resize 得到 */

    /**< 构造所有的面和边 */
    for(int i = 0; i < NC; i++)
    {
      for(int j = 0; j < 4; j++)
      {
        allFace(i*4+j, 0) = m_cell(i, m_localface[j][0]); 
        allFace(i*4+j, 1) = m_cell(i, m_localface[j][1]); 
        allFace(i*4+j, 2) = m_cell(i, m_localface[j][2]); 
      }
      for(int j = 0; j < 6; j++)
      {
        allEdge(i*6+j, 0) = m_cell(i, m_localedge[j][0]);
        allEdge(i*6+j, 1) = m_cell(i, m_localedge[j][1]);
      }
    }

    /** 构造 face */
    RArray2d<I> sortFace(allFace);
    for(int i = 0; i < NC*4; i++)
      std::sort(sortFace[i], sortFace[i]+3);

    std::vector<I> idx;
    OF::Core::sort_RArray2d(sortFace, idx);

    m_face.reserve(NC*2, 3);
    m_face2cell.reserve(NC*2, 4);
    m_cell2face.init(NC, 4);
    int NF = 0;
    for(int i = 0; i < NC*4; i++)
    {
      if(i==0 || !OF::Core::equal(sortFace[idx[i]], sortFace[idx[i-1]], 3))
      {
        m_face.push_back({allFace(idx[i], 0), allFace(idx[i], 1), allFace(idx[i], 2)});
        m_face2cell.push_back({idx[i]/4, idx[i]/4, idx[i]%4, idx[i]%4});
        m_cell2face(idx[i]/4, idx[i]%4) = NF;
        NF++;
      }
      else
      {
        m_face2cell(NF-1, 1) = idx[i]/4;
        m_face2cell(NF-1, 3) = idx[i]%4;
        m_cell2face(idx[i]/4, idx[i]%4) = NF-1;
      }
    }

    /** 构造 edge */
    RArray2d<I> sortEdge(allEdge);
    for(int i = 0; i < NC*6; i++)
      std::sort(sortEdge[i], sortEdge[i]+2);

    OF::Core::sort_RArray2d(sortEdge, idx);

    m_edge.reserve(NC*2, 2);
    m_cell2edge.init(NC, 6);
    int NE = 0;
    for(int i = 0; i < NC*6; i++)
    {
      if(i==0 || !OF::Core::equal(sortEdge[idx[i]], sortEdge[idx[i-1]], 2))
      {
        m_edge.push_back({allEdge(idx[i], 0), allEdge(idx[i], 1)});
        m_cell2edge(idx[i]/6, idx[i]%6) = NE;
        NE++;
      }
      else
      {
        m_cell2edge(idx[i]/6, idx[i]%6) = NE-1;
      }
    }
  }

  F cell_quality(const I i)
  {
    auto s = cell_surface_area(i);
    auto d = direction(i, 0);
    auto l = std::sqrt(d.squared_length());
    auto vol = cell_measure(i);
    auto R = l/vol/12.0;
    auto r = 3.0*vol/s;
    return r*3.0/R;
  }

  void cell_quality(std::vector<F> & q)
  {
    auto NC = number_of_cells();
    q.resize(NC);
    for(I i = 0; i < NC; i++)
      q[i] = cell_quality(i);
  }

  F cell_surface_area(const I i)
  {
    auto s = face_measure(m_cell2face[i][0]);
    s += face_measure(m_cell2face[i][1]);
    s += face_measure(m_cell2face[i][2]);
    s += face_measure(m_cell2face[i][3]);
    return s;
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

  RArray2d<F> & nodes()
  {
      return m_node;
  }

  RArray2d<I> & cells()
  {
      return m_cell;
  }

  RArray2d<I> & faces()
  {
      return m_face;
  }

  RArray2d<I> & edges()
  {
      return m_edge;
  }


  I * face_to_cell(const I i)
  {
      return m_face2cell[i];
  }

  I * cell_to_face(const I i)
  {
      return m_cell2face[i];
  }

  I * cell_to_edge(const I i)
  {
      return m_cell2edge[i];
  }

  // 实体测度 
  /**
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
      auto & f = m_face[i];
      auto v1 = m_node[f[1]] - m_node[f[0]];
      auto v2 = m_node[f[2]] - m_node[f[0]];
      return 0.5*std::sqrt(cross(v1, v2).squared_length());
  }

  void face_measure(std::vector<F> & measure)
  {
      auto NF = number_of_faces();
      measure.resize(NF);
      for(I i = 0; i < NF; i++)
          measure[i] = face_measure(i);
  }

  F cell_measure(const I i)
  {
      auto & c = m_cell[i];
      auto v01 = m_node[c[1]] - m_node[c[0]];
      auto v02 = m_node[c[2]] - m_node[c[0]];
      auto v03 = m_node[c[3]] - m_node[c[0]];
      return dot(cross(v01, v02), v03)/6.0;
  }

  void cell_measure(std::vector<F> & measure)
  {
      auto NC = number_of_cells();
      measure.resize(NC);
      for(I i = 0; i < NC; i++)
          measure[i] = cell_measure(i);
  }
  */

  // 实体重心
  void edge_barycenter(const I i, F * node)
  {
    auto e = m_edge[i];
    node[0] = (m_node(e[0], 0) + m_node(e[1], 0))/2.0;
    node[1] = (m_node(e[0], 1) + m_node(e[1], 1))/2.0;
    node[2] = (m_node(e[0], 2) + m_node(e[1], 2))/2.0;
  }

  /**
  Node face_barycenter(const I i)
  {
      auto & f = m_face[i];
      F x = (m_node[f[0]][0] + m_node[f[1]][0] + m_node[f[2]][0])/3.0;
      F y = (m_node[f[0]][1] + m_node[f[1]][1] + m_node[f[2]][1])/3.0;
      F z = (m_node[f[0]][2] + m_node[f[1]][2] + m_node[f[2]][2])/3.0;
      return Node(x, y, z);
  }

  void cell_barycenter(const I i, Node & node)
  {
    auto & c = m_cell[i];
    for(int i = 0; i < geo_dimension(); i++)
      node[i] = (m_node[c[0]][i] + m_node[c[1]][i] + m_node[c[2]][i] + m_node[c[3]][i])/4.0;
  }
  */

  void uniform_refine(const I n=1)
  {
    for(I i=0; i < n; i++)
    {
      auto NN = number_of_nodes();
      auto NE = number_of_edges();
      m_node.init(NN + NE, 3);
      for(I j = 0; j < NE; j++)
      {
        edge_barycenter(j, m_node[NN+j]); 
      }
      auto NC = number_of_cells();
      m_cell.init(8*NC, 4);
      for(I j = 0; j < NC; j++)
      { 
        I c[4] = {m_cell(j, 0), m_cell(j, 1), m_cell(j, 2), m_cell(j, 3)};
        m_cell(j, 0) = c[0];
        m_cell(j, 1) = m_cell2edge(j, 0) + NN;
        m_cell(j, 2) = m_cell2edge(j, 1) + NN;
        m_cell(j, 3) = m_cell2edge(j, 2) + NN;

        m_cell(NC + j, 0) = c[1];
        m_cell(NC + j, 1) = m_cell2edge(j, 3) + NN;
        m_cell(NC + j, 2) = m_cell2edge(j, 0) + NN;
        m_cell(NC + j, 3) = m_cell2edge(j, 4) + NN;

        m_cell(2*NC + j, 0) = c[2];
        m_cell(2*NC + j, 1) = m_cell2edge(j, 1) + NN;
        m_cell(2*NC + j, 2) = m_cell2edge(j, 3) + NN;
        m_cell(2*NC + j, 3) = m_cell2edge(j, 5) + NN;

        m_cell(3*NC + j, 0) = c[3];
        m_cell(3*NC + j, 1) = m_cell2edge(j, 4) + NN;
        m_cell(3*NC + j, 2) = m_cell2edge(j, 2) + NN;
        m_cell(3*NC + j, 3) = m_cell2edge(j, 5) + NN;

        m_cell(3*NC + j, 0) = c[3];
        m_cell(3*NC + j, 1) = m_cell2edge(j, 4) + NN;
        m_cell(3*NC + j, 2) = m_cell2edge(j, 2) + NN;
        m_cell(3*NC + j, 3) = m_cell2edge(j, 5) + NN;

        RArray2d<F> v(3, 3);
        for(int k = 0; k < 3; k++)
        {
          for(int l = 0; l < 3; l++)
            v(k, l) = m_node(m_cell2edge(j, k)+NN, l)-m_node(m_cell2edge(j, 5-k)+NN, l);
        }

        F tem[3] = {OF::Core::dot(v[0], v[0], 3), OF::Core::dot(v[1], v[1], 3), 
          OF::Core::dot(v[2], v[2], 3)};
        auto idx = std::distance(tem, std::min_element(tem, tem+3));

        m_cell(4*NC + j, 0) = m_cell2edge(j, m_refine[idx][0]) + NN; 
        m_cell(4*NC + j, 1) = m_cell2edge(j, m_refine[idx][1]) + NN;
        m_cell(4*NC + j, 2) = m_cell2edge(j, m_refine[idx][4]) + NN;
        m_cell(4*NC + j, 3) = m_cell2edge(j, m_refine[idx][5]) + NN;
                                                             
        m_cell(5*NC + j, 0) = m_cell2edge(j, m_refine[idx][1]) + NN; 
        m_cell(5*NC + j, 1) = m_cell2edge(j, m_refine[idx][2]) + NN;
        m_cell(5*NC + j, 2) = m_cell2edge(j, m_refine[idx][4]) + NN;
        m_cell(5*NC + j, 3) = m_cell2edge(j, m_refine[idx][5]) + NN;
                                                             
        m_cell(6*NC + j, 0) = m_cell2edge(j, m_refine[idx][2]) + NN; 
        m_cell(6*NC + j, 1) = m_cell2edge(j, m_refine[idx][3]) + NN;
        m_cell(6*NC + j, 2) = m_cell2edge(j, m_refine[idx][4]) + NN;
        m_cell(6*NC + j, 3) = m_cell2edge(j, m_refine[idx][5]) + NN;
                                                             
        m_cell(7*NC + j, 0) = m_cell2edge(j, m_refine[idx][3]) + NN; 
        m_cell(7*NC + j, 1) = m_cell2edge(j, m_refine[idx][0]) + NN;
        m_cell(7*NC + j, 2) = m_cell2edge(j, m_refine[idx][4]) + NN;
        m_cell(7*NC + j, 3) = m_cell2edge(j, m_refine[idx][5]) + NN;
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
              isBdNode[m_face(i, 0)] = true;
              isBdNode[m_face(i, 1)] = true;
              isBdNode[m_face(i, 2)] = true;
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
    auto & locid = top.local_indices();
    nei.resize(loc[NN]);
    locid.resize(loc[NN]);
    std::vector<I> start(loc);
    for(I i = 0; i < NC; i++)
    {
      for(auto j = 0; j < 4; j++)
      {
        auto v = m_cell[i][j];
        locid[start[v]] = j;
        nei[start[v]++] = i;
      }
    }
  }

  void print()
  {
      std::cout << "Nodes:" << std::endl;
      std::cout << m_node << std::endl;

      std::cout << "Edges:" << std::endl;
      std::cout << m_edge << std::endl;

      std::cout << "Faces:" << std::endl;
      std::cout << m_face << std::endl;

      std::cout << "Cells:" << std::endl;
      std::cout << m_cell << std::endl;

      std::cout << "Face2cell:" << std::endl;
      std::cout << m_face2cell << std::endl;

      std::cout << "Cell2face:" << std::endl;
      std::cout << m_cell2face << std::endl;
  }

private:
    RArray2d<F> m_node;
    RArray2d<I> m_cell; 
    RArray2d<I> m_edge;
    RArray2d<I> m_face;
    RArray2d<I> m_face2cell;
    RArray2d<I> m_cell2face;
    RArray2d<I> m_cell2edge;
};

template<typename GK>
int TetrahedronMesh<GK>::m_localedge[6][2] = {
    {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
};

template<typename GK>
int TetrahedronMesh<GK>::m_localface[4][3] = {
    {1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}
};

template<typename GK>
int TetrahedronMesh<GK>::m_localface2edge[4][3] = {
    {5, 4, 3}, {5, 1, 2}, {4, 2, 0}, {3, 0, 1}
};

template<typename GK>
int TetrahedronMesh<GK>::m_refine[3][6] = {
    {1, 3, 4, 2, 5, 0}, {0, 2, 5, 3, 4, 1}, {0, 4, 5, 1, 3, 2}
};

template<typename GK>
int TetrahedronMesh<GK>::m_index[12][4] = {
    {0, 1, 2, 3}, {0, 2, 3, 1}, {0, 3, 1, 2},
    {1, 2, 0, 3}, {1, 0, 3, 2}, {1, 3, 2, 0},
    {2, 0, 1, 3}, {2, 1, 3, 0}, {2, 3, 0, 1},
    {3, 0, 2, 1}, {3, 2, 1, 0}, {3, 1, 0, 2}
};

template<typename GK>
int TetrahedronMesh<GK>::m_vtkidx[4] = {0, 1, 2, 3};

template<typename GK>
int TetrahedronMesh<GK>::m_num[4][4] = {
    {0, 1, 2, 3}, {1, 0, 3, 2}, {2, 0, 1, 3}, {3, 0, 2, 1}
};

} // end of namespace Mesh 

} // end of namespace OF

#endif // end of TetrahedronMesh_h
