#include <math.h>

#include <algorithm>
#include "core/femcore.h"
#include "core/RArray2d.h"
#include "core/TetrahedronQuadrature.h"
#include "core/TetrahedronMesh.h"
#include "core/barycentric.h"

#include "geometry/Geometry_kernel.h"

#include <chrono>
using namespace std::chrono;

typedef OF::Core::Barycentric<3> Barycentric;
typedef OF::Core::TetrahedronQuadrature<Barycentric> Quadrature;

typedef OF::Geometry_kernel<double, int> GK;
typedef OF::Mesh::TetrahedronMesh<GK> TetMesh;

void cube_tetrahedron_mesh_RArray(TetMesh & mesh)
{
  auto & node = mesh.nodes();
  auto & cell = mesh.cells();
  node.reserve(8, 3);
  cell.reserve(6, 4);

  node.push_back({0.0, 0.0, 0.0});
  node.push_back({1.0, 0.0, 0.0});
  node.push_back({1.0, 1.0, 0.0});
  node.push_back({0.0, 1.0, 0.0});
  node.push_back({0.0, 0.0, 1.0});
  node.push_back({1.0, 0.0, 1.0});
  node.push_back({1.0, 1.0, 1.0});
  node.push_back({0.0, 1.0, 1.0});

  cell.push_back({0, 1, 2, 6});
  cell.push_back({0, 5, 1, 6});
  cell.push_back({0, 4, 5, 6});
  cell.push_back({0, 7, 4, 6});
  cell.push_back({0, 3, 7, 6});
  cell.push_back({0, 2, 3, 6});
  mesh.init_top();
  mesh.print();
  return;
}

void one_tetrahedron_mesh(TetMesh & mesh)
{
  auto & nodes = mesh.nodes();
  auto & cells = mesh.cells();
  nodes.reserve(4, 3);
  cells.reserve(1, 4);

  nodes.push_back({0.0, 0.0, 0.0});
  nodes.push_back({1.0, 0.0, 0.0});
  nodes.push_back({0.5, std::sqrt(3.0)/2.0, 0.0});
  nodes.push_back({0.5, std::sqrt(3.0)/6.0, std::sqrt(2.0/3.0)}); 

  cells.push_back({0, 1, 2, 3});
  mesh.init_top();
}

void test_mass_and_stiff(int p)
{
  /** 多重指标 */
  RArray2d<int> multiindex;
  OF::Core::get_all_multiindex(p, 4, multiindex);

  /** 积分公式 */
  Quadrature integrator(p+3);
  auto & qpoints = integrator.quadrature_points();
  auto & qws = integrator.quadrature_weights();
  for(auto & point : qpoints)
  {
    point.compute_B(p);
    point.compute_F(p);
  }

  /** 网格 */
  TetMesh mesh;
  //cube_tetrahedron_mesh_RArray(mesh);
  one_tetrahedron_mesh(mesh);
  mesh.uniform_refine(5);

  RArray2d<double> M, A;
  OF::Core::mass_matrix(qpoints, qws, multiindex, 0.11785113, M);

  RArray2d<int> localFace(4, 3);
  for(int i = 0; i < 4; i++)
  {
    for(int j = 0; j < 3; j++)
      localFace(i, j) = mesh.m_localface[i][j];
  }

  RArray2d<double> glambda;
  int NC = mesh.number_of_cells();
  auto start = high_resolution_clock::now();
  for(int i = 0; i < NC; i++)
  {
    double vol = OF::Core::grad_lambda(mesh.nodes(), mesh.cells(), i, localFace, glambda);
    OF::Core::stiff_matrix(qpoints, qws, vol, multiindex, glambda, A);
  }
  std::cout<< "NC = " << NC <<std::endl;
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  std::cout<< duration.count()/1000000.0 <<std::endl;
  //std::cout<< "glambda = " << glambda <<std::endl;
  //std::cout<< "A = " << A <<std::endl;
}

void test_sort()
{
  std::vector<int> a({0, 1, 3, 4, 3, 2, 0, 2, 1, 4, 3, 3, 5, 6, 4, 4, 3, 2});
  RArray2d<int> R(6, 3);
  R.data() = a;

  int idx[6] = {0, 1, 2, 3, 4, 5};
  int idx_inv[6] = {0, 1, 2, 3, 4, 5};

  auto f0 = [&R](int a, int b) {return !OF::Core::greater(R[a], R[b], 3);};
  auto f1 = [&idx](int a, int b) {return idx[a]<idx[b];};
  std::sort(idx, idx+6, f0);
  std::sort(idx_inv, idx_inv+6, f1);

  std::cout<< R <<std::endl;
  for(int i = 0; i < 6; i++)
  {
    std::cout<< "i = "<< i << " " << idx[i] << " " << idx_inv[i] <<std::endl;
  }
  std::cout<< "\n" <<std::endl;

  std::vector<int> iidx;
  OF::Core::sort_RArray2d(R, iidx);
  std::cout<< R <<std::endl;
  for(int i = 0; i < 6; i++)
  {
    std::cout<< "i = "<< i << " " << iidx[i] << " " <<std::endl;
  }
}

int main(int argc, char ** argv)
{
  argc++;
  int p = std::stoi(argv[1]);
  test_mass_and_stiff(p);
  //test_sort();
}




















