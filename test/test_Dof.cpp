#include <algorithm>
#include <memory>
#include <time.h>
#include <initializer_list>
#include <iostream>

#include "geometry/Geometry_kernel.h"
#include "mesh/TetrahedronMesh0.h"
#include "mesh/MeshFactory.h"
#include "TestMacro.h"
#include "functionspace/Dof.h"

typedef OF::Geometry_kernel<double, int> GK;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef OF::Mesh::TetrahedronMesh<GK, Node, Vector> TetMesh;
typedef OF::FunctionSpace::Dof<3, TetMesh> Dof;
typedef OF::Mesh::MeshFactory MF;
typedef TetMesh::Cell Cell;

int main(int argc, char ** argv)
{
  std::shared_ptr<TetMesh> mesh = std::make_shared<TetMesh>();
  auto & node = mesh->nodes();
  auto & cell = mesh->cells();
  node.resize(5);
  cell.resize(2);
  node[0] = Node{0.0, 0.0, 0.0};
  node[1] = Node{1.0, 0.0, 0.0};
  node[2] = Node{0.0, 1.0, 0.0};
  node[3] = Node{0.0, 0.0, 1.0};
  node[4] = Node{0.0, 0.0, -1.0};
  cell[0] = Cell{0, 1, 2, 3};
  cell[1] = Cell{0, 2, 1, 4};
  mesh->init_top();
  //mesh->uniform_refine(5);
  
  clock_t s0 = clock();
  Dof dofmang(std::stoi(argv[1]), mesh);
  clock_t e0 = clock();
  std::cout<< "time = " <<(double)(e0-s0)/CLOCKS_PER_SEC <<std::endl;

  //Dof dofmang(3, mesh);
  //dofmang.print();
  int NC = mesh->number_of_cells();
  std::cout<< "单元数 = " << NC <<std::endl;

  clock_t s = clock();
  for(int i = 0; i < NC; i++)
  {
    std::vector<int> c2d;
    dofmang.cell_to_dof(i, c2d);
    //std::cout<< "c2d[" << i << "] = " <<std::endl;
    //std::cout<< "[";
    //for(auto d:c2d)
    //{
      //std::cout<< d << " ";
    //}
    //std::cout<< "]" << std::endl;
  }
  clock_t e = clock();
  std::cout<< "time = " <<(double)(e-s)/CLOCKS_PER_SEC <<std::endl;
}








