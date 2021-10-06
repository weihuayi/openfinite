
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

void test_basis(int p = 3)
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
  std::vector<Vector> vb;

  Node point(4.0/3.0, 2.0/3.0);
  space.basis(0, point, vb);
  for(auto v:vb)
    std::cout<< v <<std::endl;
}

void test_curl_basis(int p = 3)
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
  std::vector<double> vb;

  Node point(10.0/6.0, 6.0/6.0);
  space.curl_basis(0, point, vb);
  for(auto v:vb)
    std::cout<< v <<std::endl;
}

void test_mass_matrix(int p = 3)
{
  typedef TriMesh::Cell Cell;
  typedef Space::BSMatrix BSMatrix;
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
  Matrix mat;
  space.mass_matrix(0, mat);
  std::cout<< mat <<std::endl;

  BSMatrix mata;
  space.mass_matrix(mata);
  std::cout<< " curl_matrix = "<<std::endl;
  for(unsigned int  i = 0; i < mata.size(); i++)
  {
    for(auto k : mata[i])
      std::cout<< "("  << i << ", " << k.first << "): " << k.second << std::endl;
  }
}

void test_curl_matrix(int p = 3)
{
  typedef TriMesh::Cell Cell;
  typedef Space::BSMatrix BSMatrix;
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

  for(int i = 0; i < mesh->number_of_edges(); i ++)
  {
    std::cout<< mesh->edge(i)[0] << " " << mesh->edge(i)[1] <<std::endl;
  }

  auto space = Space(mesh, p);
  Matrix mat;
  space.curl_matrix(1, mat);
  std::cout<< " curl_matrix = " << mat <<std::endl;

  BSMatrix mata;
  space.curl_matrix(mata);
  std::cout<< " curl_matrix = "<<std::endl;
  for(unsigned int  i = 0; i < mata.size(); i++)
  {
    for(auto k : mata[i])
      std::cout<< "("  << i << ", " << k.first << "): " << k.second << std::endl;
  }
}

void test_dof(int p = 3)
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
  auto dof = space.get_dof();
  for(int i = 0; i < 2; i++)
  {
    for(int j = 0; j < dof->number_of_local_dofs(); j++)
      std::cout<< "dof "<< i << " " << j << " = " << dof->cell_to_dof(i, j) << std::endl;
  }

  std::cout<< "\n" <<std::endl;
  for(int i = 0; i < 2; i++)
  {
    for(int j = 0; j < dof->number_of_local_dofs(); j++)
      std::cout<< "dof "<< i << " " << j << " = " << dof->cell_to_dof(i)[j] << std::endl;
  }
}

int main(int argc, char **argv)
{
  int p = std::stoi(argv[1]);
  //test_coeff(p);
  //test_basis(p);
  //test_curl_basis(p);
  test_mass_matrix(p);
  //test_curl_matrix(p);
  //test_dof(p);
}
