
#include <algorithm>
#include <memory>
#include <time.h>
#include <initializer_list>
#include <iostream>

#include "core/femcore.h"
#include "geometry/Geometry_kernel.h"
#include "core/TetrahedronMesh.h"
#include "TestMacro.h"
#include "core/BDMdof.h"

typedef OF::Geometry_kernel<double, int> GK;
typedef OF::Mesh::TetrahedronMesh<GK> TetMesh;
typedef OF::Core::BDMDof3d<TetMesh> Dof;

void two_tetrahedron_mesh(std::shared_ptr<TetMesh> & mesh)
{
  auto & node = mesh->nodes();
  auto & cell = mesh->cells();
  node.reserve(5, 3);
  cell.reserve(2, 4);
  node.push_back({0.0, -0.323, 0.0});
  node.push_back({1.0, 0.0, 0.0});
  node.push_back({0.0, 1.0, -0.2});
  node.push_back({0.0, 0.0, 1.0});
  node.push_back({0.0, 0.0, -1.0});
  //node.push_back({0.0, 0.0, 0.0});
  //node.push_back({1.0, 0.0, 0.0});
  //node.push_back({0.0, 1.0, 0.0});
  //node.push_back({0.0, 0.0, 1.0});
  //node.push_back({0.0, 0.0, -1.0});
  cell.push_back({0, 1, 2, 3});
  cell.push_back({0, 2, 1, 4});
  mesh->init_top();
  mesh->print();
  return;
}

int main(int argc, char ** argv)
{
  argc++;
  std::shared_ptr<TetMesh> mesh = std::make_shared<TetMesh>();
  two_tetrahedron_mesh(mesh);
  mesh->uniform_refine(5);
  
  clock_t s0 = clock();
  std::cout<< "s" <<std::endl;
  Dof dof(std::stoi(argv[1]), mesh);
  std::cout<< "d" <<std::endl;
  clock_t e0 = clock();
  std::cout<< "time = " <<(double)(e0-s0)/CLOCKS_PER_SEC <<std::endl;
  dof.print();

  int NC = mesh->number_of_cells();
  std::cout<< "单元数 = " << NC <<std::endl;

  int ldof = dof.number_of_local_dofs()*3;
  clock_t s = clock();
  int c2d[ldof];
  for(int i = 0; i < NC; i++)
  {
    dof.cell_to_dof(i, c2d);
    //std::cout<< "c2d[" << i << "] = " <<std::endl;
    //std::cout<< "[";
    //for(int j = 0; j < ldof; j++)
    //{
    //  std::cout<< c2d[j] << " ";
    //}
    //std::cout<< "]" << std::endl;
  }
  clock_t e = clock();
  std::cout<< "time = " <<(double)(e-s)/CLOCKS_PER_SEC <<std::endl;


  /** dof_to_vector */
  RArray2d<double> d2v;
  RArray2d<double> b2v;
  RArray2d<double> metric(3, 3); 
  ldof /= 3;
  for(int i = 0; i < 1; i++)
  {
    dof.basis_vector(i, b2v);
    std::cout<< "b2v = " << b2v <<std::endl;
    dof.dof_vector(i, d2v);
    std::cout<< "d2v = " << d2v <<std::endl;
    for(int j = 0; j < ldof; j++)
    {
      for(int k = 0; k < 3; k++)
      {
        for(int m = 0; m < 3; m++)
          metric(k, m) = OF::Core::dot(b2v[j+ldof*k], d2v[j+ldof*m]);
      }
      std::cout<< "j = " << j << " " << metric << std::endl;
    }
  }
}






