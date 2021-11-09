#ifndef Mesh2d_h
#define Mesh2d_h

#include <vector>
#include <map>

namespace OF {

namespace Mesh {

template<typename GK, typename VECTOR, typename NODE, typename CELL>
class Mesh2d
{
public:
  typedef NODE Node;
  typedef VECTOR Vector;

  typedef typename GK::Int I;
  typedef typename GK::Float F;

  typedef typename CELL::Edge Edge;
  typedef typename CELL::Cell Cell;

  typedef typename std::vector<Node>::iterator NodeIterator;
  typedef typename std::vector<Edge>::iterator EdgeIterator;
  typedef typename std::vector<Cell>::iterator CellIterator;

public:

  Mesh2d()
  {
  }

  I number_of_nodes()
  {
    return m_nodes.size();
  }

  I number_of_cells()
  {
    return m_cells.size();
  }

  I number_of_edges()
  {
    return m_edges.size();
  }

  I number_of_faces()
  {
    return m_edges.size();
  }

  int geo_dimension()
  {
      return NODE::dimension();
  }

  int top_dimension()
  {
      return CELL::dimension();
  }

  void insert(const Node & node)
  {
      m_nodes.push_back(node);
  }

  void insert(const Cell & cell)
  {
      m_cells.push_back(cell);
  }

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
          for(I j = 0; j < 4; j++)
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
          auto & c = m_cells[m_face2cell[i][0]];
          auto j = m_face2cell[i][2];
          m_face[i][0] = c[m_localface[j][0]];
          m_face[i][1] = c[m_localface[j][1]];
          m_face[i][2] = c[m_localface[j][2]];
      }

      std::map<Edge, I> eidxmap; // sorted edge 到整体唯一编号的映射
      m_cell2edge.resize(NC);
      I NE = 0;
      for(I i = 0; i < NC; i++)
      {
          for(I j = 0; j < 6; j++)
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

  void print()
  {
    for(I i = 0; i < m_NN; i++)
    {
      std::cout<< "point " << i << ":" <<  m_nodes[i] << std::endl;
    }

    for(I i = 0; i < m_NF; i++)
    {
      std::cout<< "face2node " << i << ":";
      for(int j = 0; j < Cell::NV[2]; j++)
      {
        std::cout<<  " " << face(i)[j];
      }
      std::cout << std::endl;
    }

    for(I i = 0; i < m_NE; i++)
    {
      std::cout<< "edge2node " << i << ":";
      for(int j = 0; j < Cell::NV[1]; j++)
      {
        std::cout<<  " " << edge(i)[j];
      }
      std::cout << std::endl;
    }

    auto TD = Cell::dimension();
    for(I i = 0; i < m_NC; i++)
    {
      std::cout <<"cell2node "<<  i << ":";
      for(I j = 0; j < Cell::ND[0]; j++)
        std::cout <<  " " << m_cells[i].node(j);
      std::cout << std::endl;
    }

    for(I i = 0; i < m_NC; i++)
    {
      std::cout <<"cell2edge "<<  i << ":";
      for(I j = 0; j < Cell::ND[1]; j++)
        std::cout << " " << m_cells[i].edge(j);
      std::cout << std::endl;
    }

    for(I i = 0; i < m_NC; i++)
    {
      std::cout <<"cell2face "<<  i << ":";
      for(I j = 0; j < Cell::ND[2]; j++)
        std::cout << " " << m_cells[i].face(j);
      std::cout << std::endl;
    }
  }


private:
    std::vector<Node> m_nodes;
    std::vector<Edge> m_edges;
    std::vector<Cell> m_cells;
};

} // end of namespace Mesh

} // end of namespace OF
#endif // end of Mesh_h
