#ifndef Triangle_h
#define Triangle_h

#include <array>

#include "CellType.h"

namespace OF {

namespace Mesh {

template<typename I>
class Triangle
{
public:
  typedef typename std::array<I, 3> Cell;
  typedef typename std::array<I, 2> Edge;

  typedef typename std::array<I, 4> Edge2cell;
  typedef typename std::array<I, 3> Cell2edge;

public:

  static CellType type()
  {
    return m_type;
  }

  static int dimension()
  {
    return 2;
  }

  static int order()
  {
    return 1;
  }

  static int number_of_nodes_of_cell()
  {
    return 3;
  }

  static int number_of_edges_of_cell()
  {
    return 3;
  }

  static int number_of_faces_of_cell()
  {
    return 3;
  }

  static int number_of_nodes_of_edge()
  {
    return 2;
  }

  static int number_of_nodes_of_face()
  {
    return 2;
  }

  static int number_of_edges_of_face()
  {
    return 1;
  }

  template<typename Cells>
  Edge local_sorted_edge(Cells & cells, I i, I j)
  {
    Edge e = {cells[i][edge[j][0]], cells[i][edge[j][1]]};
    std::sort(e.begin(), e.end());
    return e;
  }

public:
  static CellType type;
  static int order;
  static int localedge[3][2];
  static int vtkidx[3];
  static int index[3][3];
};

template<typename I>
int Triangle<I>::order = 1;

template<typename I>
CellType Triangle<I>::type = TRIANGLE; /**< 单元类型 */

template<typename I>
int Triangle<I>::localedge[3][2] = { /**< */
    {1, 2}, {2, 0}, {0, 1}
};

template<typename I>
int Triangle<I>::vtkidx[3] = {0, 1, 2};

template<typename I>
int Triangle<I>::index[3][3] = {
  {0, 1, 2}, {1, 2, 0}, {2, 0, 1}
};
} // end of namespace Mesh

} // end of namespace OF
#endif // end of Triangle_h
