#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <mpi.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/QuadrilateralMesh.h"
#include "mesh/ParallelMesh.h"
#include "mesh/ParallelMesher.h"
#include "mesh/GhostFillingAlg.h"
#include "mesh/ParallelMeshColoringAlg0.h"
#include "mesh/VTKMeshWriter.h"

typedef OF::Geometry_kernel<double, int> GK;
typedef GK::Point_2 Node;
typedef GK::Vector_2 Vector;
typedef OF::Mesh::QuadrilateralMesh<GK, Node, Vector> QuadMesh;
typedef OF::Mesh::ParallelMesh<GK, QuadMesh> PMesh;
typedef OF::Mesh::GhostFillingAlg<PMesh> SetGhostAlg;
typedef OF::Mesh::ParallelMesher<PMesh> PMesher;
typedef OF::Mesh::ParallelMeshColoringAlg<PMesh> PCA;
typedef OF::Mesh::VTKMeshWriter Writer;

int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);

  int rank, nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  PMesher pmesher(argv[1], ".vtu", MPI_COMM_WORLD);
  auto mesh = pmesher.get_mesh();

  auto set_ghost_alg = std::make_shared<SetGhostAlg>(mesh, MPI_COMM_WORLD);
  auto color_alg = std::make_shared<PCA>(mesh, set_ghost_alg, MPI_COMM_WORLD);

  auto NN = mesh->number_of_nodes();
  std::vector<int> color(NN);

  std::cout<< "a" <<std::endl;
  color_alg->coloring(true);

  std::stringstream ss;
  ss << "color_"<< mesh->id() << ".vtu";
  Writer writer;
  writer.set_points(*mesh);
  writer.set_cells(*mesh);
  writer.set_point_data(color_alg->m_color, 1, "color");
  writer.write(ss.str());
  
  MPI_Finalize();
  return 0;
}
