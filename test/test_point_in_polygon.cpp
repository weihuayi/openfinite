#include <iostream> 
#include <cmath> 
#include <vector>
#include <fstream>
#include <metis.h>

#include "geometry/Geometry_kernel.h"
#include "geometry/Point_2.h"

typedef OF::Geometry_kernel<double, int> GK;
typedef GK::Point_2 Node;
//typedef GK::Vector_3 Vector;
//typedef OF::Mesh::QuadrilateralMesh<GK, Node, Vector> QuadMesh;
//typedef OF::Mesh::ParallelMesh<GK, QuadMesh> PMesh;

//typedef OF::Mesh::VTKMeshReader<PMesh> Reader;
//typedef OF::Mesh::VTKMeshWriter Writer;
//typedef OF::Mesh::MeshFactory MF;


int main(int argc, char * argv[])
{
  Node p1 = Node{0.5,0.5};
  Node p2 = Node{1.2,1.2};
  std::vector<Node> vertices(4);
  vertices[0] = Node{0.0,0.0};
  vertices[1] = Node{1.0,0.0};
  vertices[2] = Node{1.0,1.0};
  vertices[3] = Node{0.0,1.0};
  
  int n = 4;
  int wn1 = GK::wn_PnPoly(p1,vertices,n);
  int cn1 = GK::cn_PnPoly(p1,vertices,n);
  int wn2 = GK::wn_PnPoly(p2,vertices,n);
  int cn2 = GK::cn_PnPoly(p2,vertices,n);
  std::cout << wn1 << std::endl;
  std::cout << cn1 << std::endl;
  std::cout << wn2 << std::endl;
  std::cout << cn2 << std::endl;
  
}
