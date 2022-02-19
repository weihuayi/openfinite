#ifndef ParallelMesh_h
#define ParallelMesh_h

#include <utility>
#include <vector>
#include <map>
#include <set>

namespace OF {
namespace Mesh {

struct ShardEntity
{
  int id; /**< 编号实体 */
  int master; /**< 主进程 */
  std::set<std::pair<int, int> > pset; /** (进程集, 在该进程中的编号) */

  ShardEntity(): id(0), master(0), pset(){}
};

template<typename GK, typename Mesh>
class ParallelMesh: public Mesh
{
public:
  typedef typename GK::Int I;
  typedef typename GK::Float F;
  typedef typename Mesh::Node Node;
  typedef typename Mesh::Edge Edge;
  typedef typename Mesh::Face Face;
  typedef typename Mesh::Cell Cell;

  typedef typename Mesh::Face2cell Face2cell;
  typedef typename Mesh::Cell2face Cell2face;
  typedef typename Mesh::Cell2edge Cell2edge;

  typedef typename Mesh::Toplogy Toplogy;

  typedef typename Mesh::NodeIterator NodeIterator;
  typedef typename Mesh::EdgeIterator EdgeIterator;
  typedef typename Mesh::FaceIterator FaceIterator;
  typedef typename Mesh::CellIterator CellIterator;

public:

  ParallelMesh(int id):
    m_id(id), m_sharedNode(), m_sharedEdge(), m_sharedFace()
  {}

  int id()
  {
    return m_id;
  }

  std::vector<ShardEntity> & shared_node()
  {
    return m_sharedNode;
  }

  std::vector<ShardEntity> & shared_edge()
  {
    return m_sharedEdge;
  }

  std::vector<ShardEntity> & shared_face()
  {
    return m_sharedFace;
  }

private:
  I m_id; // 网格块编号(进程的编号）
  std::vector<ShardEntity> m_sharedNode;
  std::vector<ShardEntity> m_sharedEdge;
  std::vector<ShardEntity> m_sharedFace;
};

} // end of namespace Mesh

} // end of namespace OF
#endif // end of ParallelMesh_h
