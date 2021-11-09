#include <sstream> 
#include <cmath> 

#include <fstream>
#include <metis.h>

#include <vtkDoubleArray.h>
#include <vtkIntArray.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/QuadrilateralMesh.h"
#include "mesh/ParallelMesh.h"
#include "mesh/MeshFactory.h"
#include "mesh/VTKMeshWriter.h"
#include "mesh/VTKMeshReader.h"

typedef OF::Geometry_kernel<double, int> GK;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef OF::Mesh::QuadrilateralMesh<GK, Node, Vector> QuadMesh;
typedef OF::Mesh::ParallelMesh<GK, QuadMesh> PMesh;

typedef OF::Mesh::VTKMeshReader<PMesh> Reader;
typedef OF::Mesh::VTKMeshWriter Writer;
typedef OF::Mesh::MeshFactory MF;
typedef QuadMesh::Cell Cell;

template<class Mesh>
void Mesher(std::shared_ptr<Mesh> mesh, double h1 = 0.1, int n0 = 10, int n1 = 10) 
{
  QuadMesh qmesh;
  auto & node = qmesh.nodes();
  auto & cell = qmesh.cells();
  node.resize(8);
  node[0] = Node{0.0, 0.0, 0.0};
  node[1] = Node{1.0, 0.0, 0.0};
  node[2] = Node{1.0, 1.0, 0.0};
  node[3] = Node{0.0, 1.0, 0.0};
  node[4] = Node{0.0, 0.0, 1.0};
  node[5] = Node{1.0, 0.0, 1.0};
  node[6] = Node{1.0, 1.0, 1.0};
  node[7] = Node{0.0, 1.0, 1.0};
  
  cell.resize(6);
  cell[0] = Cell{0, 1, 2, 3};
  cell[1] = Cell{4, 5, 6, 7};
  cell[2] = Cell{0, 1, 5, 4};
  cell[3] = Cell{2, 3, 7, 6};
  cell[4] = Cell{1, 2, 6, 5};
  cell[5] = Cell{3, 0, 4, 7};
  mesh.init_top();

  auto & edge = mesh.edges();
  for(auto e : edge)
  {}
  


}




















