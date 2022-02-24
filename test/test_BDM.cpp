
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
#include "quadrature/LagrangeTriangleMeshQuadrature.h"

#include "mesh/LagrangeTriangleMesh.h"
#include "mesh/TriangleMesh.h"

#include "functionspace/BDMFEMSpace2d.h"

typedef OF::Algebra_kernel<double, int> AK;
typedef AK::Matrix Matrix;
typedef OF::Geometry_kernel<double, int> GK;
typedef GK::Point_2 Node;
typedef GK::Vector_2 Vector;

void test_basis(int p=1)
{
  typedef OF::Mesh::TriangleMesh<GK, Node, Vector> TMesh;
  typedef TMesh::Cell Cell;
  typedef OF::Mesh::LagrangeTriangleMesh<GK, Node, Vector, 1> Mesh;
  typedef OF::Quadrature::LagrangeTriangleMeshQuadrature<Mesh> TQuadrature;
  typedef OF::FunctionSpace::BDMSpace2d<Mesh, AK, TQuadrature> Space;
  typedef Space::BSMatrix BSMatrix;

  auto tmesh = std::make_shared<TMesh>();
  auto & nodes = tmesh->nodes();
  auto & cells = tmesh->cells();

  nodes.resize(4);
  cells.resize(2);
  nodes[0] = Node({0.0, 0.0});
  nodes[1] = Node({1.0, 0.0});
  nodes[2] = Node({1.0, 1.0});
  nodes[3] = Node({0.0, 1.0});
  cells[0] = Cell({0, 1, 2});
  cells[1] = Cell({0, 2, 3});
  tmesh->init_top();

  tmesh->uniform_refine(1);
  int NE = tmesh->number_of_edges();
  for(int i = 0; i < NE; i++)
    std::cout<< i << " : [" << tmesh->edge(i)[0] << ", " << tmesh->edge(i)[1] << "]" <<std::endl;

  auto mesh = std::make_shared<Mesh>(tmesh);
  std::shared_ptr<Space> space = std::make_shared<Space>(mesh, p);

  std::cout<< "开始测试 space cell_to_dof" << std::endl;
  auto dof = space->get_dof();
  int ldof = dof->number_of_local_dofs();
  std::cout<< "ldof = " << ldof <<std::endl;
  for(int i = 0; i < ldof; i++)
    std::cout<< " cell[0][" << i << "] = " << dof->cell_to_dof(0)[i] <<std::endl;
  for(int i = 0; i < ldof; i++)
    std::cout<< " cell[1][" << i << "] = " << dof->cell_to_dof(1)[i] <<std::endl;

  int edof = dof->number_of_edge_local_dofs();
  std::cout<< "edof = " << edof <<std::endl;
  for(int i = 0; i < edof; i++)
    std::cout<< " edge[0][" << i << "] = " << dof->edge_to_dof(0)[i] <<std::endl;
  for(int i = 0; i < edof; i++)
    std::cout<< " edge[1][" << i << "] = " << dof->edge_to_dof(1)[i] <<std::endl;

  const auto & dir = dof->get_direction();
  for(int i = 0; i < (int)dir.size(); i++)
    std::cout<< i << " : " << dir[i] <<std::endl; 

  const auto & tan = dof->get_tan_flag();
  for(int i = 0; i < (int)tan.size(); i++)
    std::cout<< i << " : " << tan[i] <<std::endl; 

  std::array<double, 3> bc{0.3, 0.3, 0.4};
  //std::array<double, 3> bc{1, 0, 0};

  std::cout<< "开始测试 space basis" << std::endl;
  std::vector<Vector> val;
  space->basis(0, bc, val);
  for(unsigned int i = 0; i < val.size(); i++)
  {
    std::cout<< "0 val[" << i << "] = " << val[i] <<std::endl;
  }
  space->basis(1, bc, val);
  for(unsigned int i = 0; i < val.size(); i++)
  {
    std::cout<< "1 val[" << i << "] = " << val[i] <<std::endl;
  }

  /*
  std::cout<< "开始测试 插值" << std::endl;
  std::function<Vector(const Node &)> test_f = [] (const Node & p)->Vector
  {
    //double a = std::pow(p[0], 2)*p[1]+20*std::sin(p[0]+p[1]);
    //double b = std::pow(p[0], 2)*p[1]+20*std::sin(p[0]+p[1]);
    double b = -p[0]*p[0];
    double a = -p[1]*p[1];
    return Vector(a, b);
  };

  tmesh->print();
  for(int i = 0; i < 1; i++)
  {
    //tmesh->uniform_refine();
    mesh = std::make_shared<Mesh>(tmesh);
    space = std::make_shared<Space>(mesh, p);
    auto intep_f = space->inteprate(test_f);
    int gdof = space->get_dof()->number_of_dofs();
    for(int j = 0; j < gdof; j++)
      std::cout<< j << " : " << intep_f[j] <<std::endl;
    double val = space->L2_error(intep_f, test_f);
    std::cout<< "加密次数: " << i << " 误差为: " << val <<std::endl;
  }

  std::cout<< "开始测试 space mass mairix" << std::endl;
  BSMatrix mass;
  space->mass_matrix(mass);
  for(int i = 0; i < (int)mass.size(); i++)
  {
    for(const auto & pa : mass[i])
      std::cout<<  "(" << i << ", " << pa.first << "):  " << pa.second <<std::endl;
  }

  std::cout<< "开始测试 space source vector" << std::endl;
  std::vector<double> sou_vec;
  space->source_vector(test_f, sou_vec);
  for(int i = 0; i < (int)sou_vec.size(); i++)
  {
    std::cout<< "v[" << i << "] = " << sou_vec[i] <<std::endl;
  }

  std::cout<< "开始测试 space grad basis" << std::endl;
  std::vector<Vector> gval;
  space->grad_basis(0, bc, gval);
  for(unsigned int i = 0; i < gval.size(); i++)
  {
    std::cout<< "val[" << i << "] = " << gval[i] <<std::endl;
  }


  std::cout<< "开始测试 space stiff mairix" << std::endl;
  BSMatrix stiff;
  space->stiff_matrix(stiff);
  for(int i = 0; i < (int)stiff.size(); i++)
  {
    for(const auto & pa : stiff[i])
      std::cout<<  "(" << i << ", " << pa.first << "):  " << pa.second <<std::endl;
  }



  std::cout<< "开始测试 插值" << std::endl;
  std::function<double(const Node &)> test_f = [] (const Node & p)->double 
  {
    return std::pow(p[0], 2)*p[1]+20*std::sin(p[0]+p[1]);
  };

  tmesh->print();
  for(int i = 0; i < 4; i++)
  {
    tmesh->uniform_refine();
    mesh = std::make_shared<Mesh>(tmesh);
    space = std::make_shared<Space>(mesh, p);
    auto intep_f = space->inteprate(test_f);
    double val = space->L2_error(intep_f, test_f);
    std::cout<< "加密次数: " << i << " 误差为: " << val <<std::endl;
  }
  */
}


int main(int argc, char **argv)
{
  argc++;
  int p = std::stoi(argv[1]);
  test_basis(p);
}
