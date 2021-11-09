
#include <iostream>
#include <memory>
#include <vector>
#include <cmath>
#include <iomanip>
#include <time.h>

#include "TestMacro.h"

#include "geometry/Geometry_kernel.h"
#include "mesh/TriangleMesh.h"
#include "mesh/LagrangeTriangleMesh.h"

using namespace OF;

typedef Geometry_kernel<double, int> GK;
typedef GK::Point_2 Node;
typedef GK::Vector_2 Vector;
typedef Mesh::TriangleMesh<GK, Node, Vector> TriMesh;
typedef Mesh::LagrangeTriangleMesh<GK, Node, Vector, 1> LTMesh;
typedef TriMesh::Cell Cell;
typedef TriMesh::Edge Edge;

void test_insert_construction()
{
  std::shared_ptr<TriMesh> mesh = std::make_shared<TriMesh>();

  mesh->insert(Node{0.0, 0.0});
  mesh->insert(Node{1.0, 0.0});
  mesh->insert(Node{1.0, 1.0});
  mesh->insert(Node{0.0, 1.0});
  mesh->insert(Cell{1, 2, 0});
  mesh->insert(Cell{3, 0, 2});
  mesh->init_top();

  LTMesh lmesh(mesh);
  /**
  
  //测试节点, 边, 单元数量
  ASSERT_EQUAL(4, mesh->number_of_nodes());
  ASSERT_EQUAL(2, mesh->number_of_cells());
  ASSERT_EQUAL(5, mesh->number_of_edges());

  //测试 edge
  std::vector<std::array<int, 2> > edge_test(5);
  edge_test[0] = std::array<int, 2>({2, 0});
  edge_test[1] = std::array<int, 2>({0, 1});
  edge_test[2] = std::array<int, 2>({1, 2});
  edge_test[3] = std::array<int, 2>({2, 3});
  edge_test[4] = std::array<int, 2>({3, 0});
  for(int i = 0; i < 5; i++)
  {
    ASSERT_EQUAL(edge_test[i][0], mesh->edge(i)[0]);
    ASSERT_EQUAL(edge_test[i][1], mesh->edge(i)[1]);
  }
  */
}

int main(int argc, char **argv)
{
  test_insert_construction();
}
