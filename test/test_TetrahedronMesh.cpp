#include <math.h>
#include <algorithm>
#include "core/femcore.h"
#include "core/RArray2d.h"
#include "core/TetrahedronQuadrature.h"
#include "mesh/TetrahedronMesh.h"
#include "core/barycentric.h"

#include "geometry/Geometry_kernel.h"
#include <chrono>
using namespace std::chrono;

typedef OF::Geometry_kernel<double, int> GK;

typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef OF::Mesh::TetrahedronMesh<GK, Node, Vector> TetMesh;

template<typename Mesh>
void cube_tetrahedron_mesh_RArray(Mesh & mesh)
{
  auto & node = mesh.nodes();
  auto & cell = mesh.cells();

  node.push_back({0.0, 0.0, 0.0});
  node.push_back({1.0, 0.0, 0.0});
  node.push_back({1.0, 1.0, 0.0});
  node.push_back({0.0, 1.0, 0.0});
  node.push_back({0.0, 0.0, 1.0});
  node.push_back({1.0, 0.0, 1.0});
  node.push_back({1.0, 1.0, 1.0});
  node.push_back({0.0, 1.0, 1.0});

  cell.push_back({0, 1, 2, 6});
  cell.push_back({0, 5, 1, 6});
  cell.push_back({0, 4, 5, 6});
  cell.push_back({0, 7, 4, 6});
  cell.push_back({0, 3, 7, 6});
  cell.push_back({0, 2, 3, 6});
  mesh.init_top();
  mesh.print();
}

void one_tetrahedron_mesh(TetMesh & mesh)
{
  auto & node = mesh.nodes();
  auto & cell = mesh.cells();

  node.push_back({0.0, 0.0, 0.0});
  node.push_back({1.0, 0.0, 0.0});
  node.push_back({1.0, 1.0, 0.0});
  node.push_back({0.0, 1.0, 0.0});

  cell.push_back({0, 1, 2, 3});
  mesh.init_top();
  mesh.print();
}

int main()
{
  TetMesh mesh;
  one_tetrahedron_mesh(mesh);
  mesh.uniform_refine(7);
  mesh.clear();
  auto start = high_resolution_clock::now();
  mesh.init_top();
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  std::cout<< duration.count()/1000000.0 <<std::endl;
}

