
#include <math.h>

#include <iostream>
#include <memory>
#include <vector>
#include <cmath>
#include <iomanip>
#include <time.h>

#include "TestMacro.h"

#include "algebra/Algebra_kernel.h"
#include "geometry/Geometry_kernel.h"
#include "quadrature/TriangleMeshQuadrature.h"
#include "quadrature/QuadrilateralMeshQuadrature.h"

#include "mesh/QuadrilateralMesh.h"
#include "mesh/TriangleMesh.h"
#include "mesh/TetrahedronMesh.h"
#include "mesh/HexahedronMesh.h"

#include "functionspace/FirstKindNedecFiniteElementSpace2d.h"
#include "quadrature/TriangleQuadrature.h"

typedef OF::Algebra_kernel<double, int> AK;
typedef AK::Matrix Matrix;
typedef OF::Geometry_kernel<double, int> GK;
typedef GK::Point_2 Node;
typedef GK::Vector_2 Vector;
typedef OF::Mesh::TriangleMesh<GK, Node, Vector> TriMesh;
typedef OF::Quadrature::TriangleMeshQuadrature<TriMesh> TriQuadrature;
typedef OF::FunctionSpace::FirstKindNedecFiniteElementSpace2d<TriMesh, AK, TriQuadrature> Space;

void test_coeff(int p = 3)
{
  typedef TriMesh::Cell Cell;
  auto mesh = std::make_shared<TriMesh>();
  auto & nodes = mesh->nodes();
  auto & cells = mesh->cells();

  nodes.resize(4);
  cells.resize(2);
  nodes[0] = Node({0.0, 0.0});
  nodes[1] = Node({2.0, 0.0});
  nodes[2] = Node({2.0, 2.0});
  nodes[3] = Node({0.0, 2.0});
  cells[0] = Cell({0, 1, 2});
  cells[1] = Cell({0, 2, 3});
  mesh->init_top();

  auto space = Space(mesh, p);
  Matrix coeff; 
  space.basis_coefficients(0, coeff);
  std::cout<< coeff <<std::endl;
}

int main(int argc, char **argv)
{
  int p = std::stoi(argv[1]);
  test_coeff(p);
}
