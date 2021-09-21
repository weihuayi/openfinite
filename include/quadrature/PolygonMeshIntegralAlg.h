#ifndef PolygonMeshIntegralAlg_h
#define PolygonMeshIntegralAlg_h

#include <math.h>
#include <memory>
#include <vector>
#include <array>

namespace OF {
namespace Quadrature{

template<typename Mesh>
class PolygonMeshIntegralAlg
{
public:

  typedef typename Mesh::F F;
  typedef typename Mesh::I I;
  typedef typename Mesh::Node Node;
  typedef typename Mesh::Vector Vector;

public:
  PolygonMeshIntegralAlg(std::shared_ptr<Mesh> mesh, int q): m_mesh(mesh), m_q(q) {}




private:
  int m_q;
  std::shared_ptr<Mesh> m_mesh;
};// end of PolygonMeshIntegralAlg

}//end of Quadrature
}//end of OF
#endif // end of PolygonMeshIntegralAlg_h
