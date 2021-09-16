

/**
 * \file test_TriangleQuadrature.cpp
 * \author Chunyu Chen
 * \date 2021/09/08
 * \brief TriangleQuadrature.h 测试文件
 */
#include <array>
#include <vector>
#include <algorithm>
#include <math.h>
#include <functional>

#include "TestMacro.h"
#include "geometry/Point_2.h"
#include "geometry/Vector_2.h"
#include "mesh/TriangleMesh.h"
#include "quadrature/TriangleQuadrature.h"
#include "quadrature/TriangleMeshQuadrature.h"

#include "algebra/Algebra_kernel.h"
#include "geometry/Geometry_kernel.h"

#include "mesh/QuadrilateralMesh.h"
#include "mesh/TriangleMesh.h"
#include "mesh/TetrahedronMesh.h"
#include "mesh/HexahedronMesh.h"
#include "quadrature/TriangleQuadrature.h"

typedef OF::Algebra_kernel<double, int> AK;
typedef AK::Matrix Matrix;
typedef OF::Geometry_kernel<double, int> GK;
typedef GK::Point_2 Node;
typedef GK::Vector_2 Vector;

typedef OF::GeometryObject::Point_2<double> Point;
typedef OF::Quadrature::TriangleQuadrature TriangleQuadrature;

typedef OF::Mesh::TriangleMesh<GK, Node, Vector> TriMesh;
typedef OF::Quadrature::TriangleMeshQuadrature<TriMesh> TriQuadrature;

/** 
 * \brief 被积函数
 * \param p 函数的参数
 */
double f(const Point & p)
{
  return p[0];
}

void test_integral(int p = 3)
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

  TriQuadrature quadr(mesh, p);
  std::function<void(const Node&, Matrix&)> phi_jk = []
    (const Node & p, Matrix & m)->void
  {
    auto x = p[0];
    auto y = p[1];
    m.reshape(2, 2);
    m[0][0] = x*x;
    m[0][1] = x*y;
    m[1][0] = x;
    m[1][1] = 1;
  };
  Matrix mat;
  quadr.integral(0, phi_jk, mat);

  Matrix mat0;
  quadr.edge_integral(1, phi_jk, mat0);
  std::cout<< mat <<std::endl;
  std::cout<< mat0 <<std::endl;
}

int main(int args, char *argv[])
{
  int q = 4;
  if(args>1)
     q = std::stoi(argv[1]);
  test_integral(q);
}
