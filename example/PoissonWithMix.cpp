
#include <math.h>

#include <iostream>
#include <memory>
#include <vector>
#include <cmath>
#include <iomanip>
#include <time.h>

#include "algebra/Algebra_kernel.h"
#include "geometry/Geometry_kernel.h"
#include "quadrature/LagrangeTriangleMeshQuadrature.h"

#include "mesh/LagrangeTriangleMesh.h"
#include "mesh/TriangleMesh.h"
#include "Python.h"

#include "functionspace/BDMFEMSpace2d.h"

typedef OF::Algebra_kernel<double, int> AK;
typedef AK::Matrix Matrix;
typedef OF::Geometry_kernel<double, int> GK;
typedef GK::Point_2 Node;
typedef GK::Vector_2 Vector;

void solve(std::vector<int> & I, std::vector<int> & J, 
    std::vector<double> & data, std::vector<double> & b, std::vector<double> & x)
{
  int N0 = data.size();
  int N1 = b.size();

  Py_Initialize();
  PyRun_SimpleString("import sys");
  PyRun_SimpleString("sys.path.append('../example')");
  PyRun_SimpleString("import axb");

  PyObject* pModule = PyImport_ImportModule("axb");
  PyObject* pFunc = PyObject_GetAttrString(pModule, "Solve");//要运行的函数

  PyObject* plist0 = PyList_New(N0);//函数的参数是一个list
  PyObject* plist1 = PyList_New(N0);//函数的参数是一个list
  PyObject* plist2 = PyList_New(N0);//函数的参数是一个list
  PyObject* plist3 = PyList_New(N1);//函数的参数是一个list

  PyObject* ptuple0 = PyTuple_New(4);//把参数用 tuple 装起来

  for(int i = 0; i < N0; i++)
  {
    PyObject*  pra = Py_BuildValue("i", I[i]);
    PyList_SetItem(plist0, i, pra);
  }

  for(int i = 0; i < N0; i++)
  {
    PyObject*  pra = Py_BuildValue("i", J[i]);
    PyList_SetItem(plist1, i, pra);
  }

  for(int i = 0; i < N0; i++)
  {
    PyObject*  pra = Py_BuildValue("d", data[i]);
    PyList_SetItem(plist2, i, pra);
  }

  for(int i = 0; i < N1; i++)
  {
    PyObject*  pra = Py_BuildValue("d", b[i]);
    PyList_SetItem(plist3, i, pra);
  }

  PyTuple_SetItem(ptuple0, 0, plist0);
  PyTuple_SetItem(ptuple0, 1, plist1);
  PyTuple_SetItem(ptuple0, 2, plist2);
  PyTuple_SetItem(ptuple0, 3, plist3);

  PyObject *pRet = PyObject_CallObject(pFunc, ptuple0);//运行函数

  int count = (int)PyList_Size(pRet);
  for(int i = 0 ; i < count ; i++ )
  {
    PyObject *Item = PyList_GetItem(pRet, i);
    PyArg_Parse(Item, "d", &x[i]);
  }

	Py_Finalize();
}
void test_basis(int p=1, int iter = 0)
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

  tmesh->uniform_refine(iter);
  auto mesh = std::make_shared<Mesh>(tmesh);
  std::shared_ptr<Space> space = std::make_shared<Space>(mesh, p);

  std::cout<< "开始测试 space cell_to_dof" << std::endl;
  std::cout<< "开始测试 插值" << std::endl;
  std::function<Vector(const Node &)> test_f = [] (const Node & p)->Vector
  {
    double b = -p[0] - 1;
    double a = -p[1];
    return Vector(a, b);
  };

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

  int ND = mass.size();
  std::vector<int> I;
  std::vector<int> J;
  std::vector<double> data;
  for(int m = 0; m < ND; m++)
  {
    for(auto it : mass[m])
    {
      I.push_back(m);
      J.push_back(it.first);
      data.push_back(it.second);
    }
  }
  int gdof = space->get_dof()->number_of_dofs();
  std::vector<double> bb(gdof);
  for(int i = 0; i < gdof; i++)
  {
    bb[i] = sou_vec[i];
  }

  std::vector<double> val(gdof);
  solve(I, J, data, bb, val);

  auto uh = space->function(); 
  uh.set_data(val);

  double err = space->L2_error(uh, test_f);
  std::cout<< "err = " << err <<std::endl;

  /*
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
  int iter = std::stoi(argv[2]);
  test_basis(p, iter);
}
