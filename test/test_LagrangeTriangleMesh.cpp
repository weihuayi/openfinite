
#include <iostream>
#include <memory>
#include <vector>
#include <cmath>
#include <iomanip>
#include <time.h>

#include "TestMacro.h"

#include "geometry/Geometry_kernel.h"
#include "mesh/VTKMeshWriter.h"
#include "mesh/TriangleMesh.h"
#include "mesh/LagrangeTriangleMesh.h"

using namespace OF;

typedef Geometry_kernel<double, int> GK;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef Mesh::TriangleMesh<GK, Node, Vector> TriMesh;
typedef Mesh::LagrangeTriangleMesh<GK, Node, Vector, 3> LTMesh;
typedef Mesh::VTKMeshWriter Writer;
typedef TriMesh::Cell Cell;
typedef TriMesh::Edge Edge;

void test_construction()
{
  std::shared_ptr<TriMesh> mesh = std::make_shared<TriMesh>();

  mesh->insert(Node{0.0, 0.0, 0.0});
  mesh->insert(Node{1.0, 0.0, 0.0});
  mesh->insert(Node{1.0, 1.0, 0.0});
  mesh->insert(Node{0.0, 1.0, 0.0});
  mesh->insert(Node{0.0, 0.0, 1.0});
  mesh->insert(Node{1.0, 0.0, 1.0});
  mesh->insert(Node{1.0, 1.0, 1.0});
  mesh->insert(Node{0.0, 1.0, 1.0});
  mesh->insert(Cell{1, 2, 0});
  mesh->insert(Cell{3, 0, 2});
  mesh->insert(Cell{1, 2, 6});
  mesh->insert(Cell{6, 5, 1});
  mesh->insert(Cell{6, 2, 3});
  mesh->insert(Cell{7, 6, 3});
  mesh->insert(Cell{0, 4, 7});
  mesh->insert(Cell{3, 0, 7});
  mesh->insert(Cell{0, 1, 5});
  mesh->insert(Cell{0, 5, 4});
  mesh->insert(Cell{4, 5, 6});
  mesh->insert(Cell{4, 6, 7});
  mesh->init_top();
  mesh->uniform_refine(2);

  LTMesh lmesh(mesh);
  auto & node = lmesh.nodes();
  for(auto & n : node)
  {
    n -= Node{0.5, 0.5, 0.5};
    n /= std::sqrt((n-Node{0.0, 0.0, 0.0}).squared_length());
  }
  lmesh.print();

  Writer writer;                                                                
  writer.set_points(lmesh);                                                     
  writer.set_cells(lmesh);                                                      
  writer.write("lmesh.vtu");                        
}

void test_grad_shape_function()
{
  std::shared_ptr<TriMesh> mesh = std::make_shared<TriMesh>();

  mesh->insert(Node{0.0, 0.0, 1.0});
  mesh->insert(Node{1.0, 0.0, 1.0});
  mesh->insert(Node{1.0, 1.0, 1.0});
  mesh->insert(Node{0.0, 1.0, 1.0});
  mesh->insert(Cell{1, 2, 0});
  mesh->insert(Cell{3, 0, 2});
  mesh->init_top();

  LTMesh lmesh(mesh);
  std::vector<double> bc{0.3, 0.3, 0.4};
  std::vector<double> val;
  lmesh.shape_function(bc, val);
  for(int i = 0; i < (int)val.size(); i++)
    std::cout<< "val[" << i << "] = " << val[i] <<std::endl;

  std::cout<< "\n" <<std::endl;
  std::vector<GK::Vector_2> gval;
  lmesh.grad_shape_function(bc, gval);
  for(int i = 0; i < (int)gval.size(); i++)
    std::cout<< "gval[" << i << "] = " << gval[i] <<std::endl;

  std::cout<< "\n" <<std::endl;
  std::vector<GK::Vector_2> jval;
  lmesh.jacobi_matrix(0, bc, jval);
  for(int i = 0; i < (int)jval.size(); i++)
    std::cout<< "jval[" << i << "] = " << jval[i] <<std::endl;
}

int main()
{
  //test_construction();
  test_grad_shape_function();
}
