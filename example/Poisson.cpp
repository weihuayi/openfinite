/**
 * \file Poisson.cpp
 * \brief Poisson 方程求解的例子
 * \author chenchunyu
 */
#include <math.h>

#include <iostream>
#include <memory>
#include <vector>
#include <cmath>
#include <iomanip>
#include <time.h>
#include <chrono>

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

#include "algebra/Algebra_kernel.h"
#include "geometry/Geometry_kernel.h"
#include "quadrature/LagrangeTriangleMeshQuadrature.h"

#include "mesh/LagrangeTriangleMesh.h"
#include "mesh/TriangleMesh.h"
#include "Python.h"

#include "functionspace/LagrangeFiniteElementSpace2d.h"

typedef OF::Algebra_kernel<double, int> AK;
typedef OF::Geometry_kernel<double, int> GK;
typedef GK::Point_2 Node;
typedef GK::Vector_2 Vector;

typedef OF::Mesh::TriangleMesh<GK, Node, Vector> TMesh;
typedef TMesh::Cell Cell;
typedef OF::Mesh::LagrangeTriangleMesh<GK, Node, Vector, 1> Mesh;
typedef OF::Quadrature::LagrangeTriangleMeshQuadrature<Mesh> TQuadrature;
typedef OF::FunctionSpace::LagrangeFiniteElementSpace2d<Mesh, AK, TQuadrature> Space;
typedef Space::BSMatrix BSMatrix;

typedef Eigen::Triplet<double> Triplet;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> CSRMatrix;
typedef Eigen::VectorXd VectorXd;
typedef Eigen::ConjugateGradient<CSRMatrix, Eigen::Lower|Eigen::Upper> Solver;

const double pi = M_PI;

using namespace std::chrono;

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

/**
 * \brief 右端项函数
 */
double source_function(const Node & p)
{
  return 2*pi*pi*std::sin(pi*p[0])*std::sin(pi*p[1]);
}

/**
 * \brief 真解
 */
double trust_solution(const Node & p)
{
  return std::sin(pi*p[0])*std::sin(pi*p[1]);
}

template<typename Triplet>
void to_coo(BSMatrix & K, std::vector<Triplet> & tlist)
{
  int NN = K.size();
  // 获得非零元的个数
  for(int m = 0; m < NN; m++)
  {
    for(auto it : K[m])
    {
      if(it.second>1e-100 || it.second<-1e-100)
        tlist.push_back(Triplet(m, it.first, it.second));
    }
  }
}

int main(int argc, char ** argv)
{
  argc++;
  int p = std::stoi(argv[1]);

  // 准备线性网格
  
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
  tmesh->uniform_refine(std::stoi(argv[2]));

  // 生成高阶网格
  auto mesh = std::make_shared<Mesh>(tmesh);

  // 生成有限元空间
  std::shared_ptr<Space> space = std::make_shared<Space>(mesh, p);
  int gdof = space->get_dof()->number_of_dofs();

  // 生成刚度矩阵的 BSMatrix
  std::cout<< "开始 stiff mairix" << std::endl;

  auto start = high_resolution_clock::now();
  BSMatrix stiff;
  space->stiff_matrix(stiff);
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  std::cout << "BSMatrix 刚度矩阵组装时间: "
       << duration.count()/1e+6 << " seconds" << std::endl;

  start = high_resolution_clock::now();
  BSMatrix mass;
  space->mass_matrix(mass);
  stop = high_resolution_clock::now();
  duration = duration_cast<microseconds>(stop - start);
  std::cout << "BSMatrix 刚度矩阵组装时间: "
       << duration.count()/1e+6 << " seconds" << std::endl;

  // 设置右端向量
  VectorXd x(gdof);
  VectorXd b(gdof);
  for(int i = 0; i < gdof; i++)
    b[i] = 0;
  std::function<double(const Node &)> F = [] (const Node & p)->double 
  {
    return source_function(p);
  };
  space->source_vector(F, b);

  // 处理边界条件
  std::vector<bool> isbddof;
  space->get_dof()->is_boundary_dof(isbddof);
  for(int i = 0; i < (int)isbddof.size(); i++)
  {
    if(isbddof[i])
    {
      b[i] = 0.0;
      stiff[i].clear();
      stiff[i][i] = 1.0;
    }
  }

  //准备数据调用 numpy 计算 Ax = b
  int ND = stiff.size();
  std::vector<int> I;
  std::vector<int> J;
  std::vector<double> data;
  for(int m = 0; m < ND; m++)
  {
    for(auto it : stiff[m])
    {
      I.push_back(m);
      J.push_back(it.first);
      data.push_back(it.second);
    }
  }
  std::vector<double> bb(gdof);
  for(int i = 0; i < gdof; i++)
  {
    bb[i] = b[i];
  }

  /*
  // 将 BSMatrix 转化为 CSRMatrix
  std::vector<Triplet> tlist;
  to_coo(stiff, tlist);

  CSRMatrix M(gdof, gdof);
  M.setFromTriplets(tlist.begin(), tlist.end());

  // 解方程
  Eigen::SparseLU<CSRMatrix, Eigen::COLAMDOrdering<int> > solver0;
  solver0.analyzePattern(M);
  solver0.factorize(M);
  x = solver0.solve(b)

  // 设置数值解函数
  */
  std::vector<double> val(gdof);
  solve(I, J, data, bb, val);

  auto uh = space->function(); 
  uh.set_data(val);

  std::function<double(const Node &)> Trust_s = [] (const Node & p)->double 
  {
    return trust_solution(p);
  };
  double err = space->L2_error(uh, Trust_s);
  std::cout<< "误差为: " << err <<std::endl;
  std::cout<< "组装BM: " << space->runtime <<std::endl;
  std::cout<< "组装BM: " << space->get_mesh()->runtime <<std::endl;
}


