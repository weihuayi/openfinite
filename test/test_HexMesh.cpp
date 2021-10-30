/**
 * \file test_HexMesh
 * \author ccy
 * \date 09/29/2021
 * \brief 六面体网格测试文件
 */
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <time.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/HexahedronMesh.h"
#include "mesh/VTKMeshWriter.h"
#include "mesh/MeshFactory.h"

typedef OF::Geometry_kernel<double, int> GK;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;
typedef OF::Mesh::HexahedronMesh<GK, Node, Vector> HexMesh;
typedef OF::Mesh::VTKMeshWriter Writer;
typedef OF::Mesh::MeshFactory MF;
typedef HexMesh::Cell Cell;

int main()
{
    HexMesh mesh;

    mesh.insert(Node{0.0, 0.0, 0.0});
    mesh.insert(Node{1.0, 0.0, 0.0});
    mesh.insert(Node{1.0, 1.0, 0.0});
    mesh.insert(Node{0.0, 1.0, 0.0});
    mesh.insert(Node{0.0, 0.0, 1.0});
    mesh.insert(Node{1.0, 0.0, 1.0});
    mesh.insert(Node{1.0, 1.0, 1.0});
    mesh.insert(Node{0.0, 1.0, 1.0});

    auto & node = mesh.nodes();
    auto & cell = mesh.cells();

    node.resize(12);
    cell.resize(2);
    node[8] = Node({2.0, 0.0, 0.0});
    node[9] = Node({2.0, 1.0, 0.0});
    node[10] = Node({2.0, 1.0, 1.0});
    node[11] = Node({2.0, 0.0, 1.0});

    cell[0] = Cell({0, 1, 2, 3, 4, 5, 6, 7});
    cell[1] = Cell({1, 8, 9, 2, 5, 11, 10, 6});

    mesh.init_top();
    
    std::vector<double> q;
    mesh.uniform_refine(1);
    
    int NN = mesh.number_of_nodes();
    for(int i = 0; i < NN; i++)
      std::cout<< i << " : " << node[i] <<std::endl; 

    std::vector<int> gtag(NN, 0);
    std::vector<int> gdof(NN, 0);
    node[35] = Node{1, 0.4, 0.4};
    gdof[35] = 3;
    node[43] = Node{0.4, 0.4, 0.4};
    gdof[43] = 3;
    node[44] = Node{1.4, 0.6, 0.4};
    gdof[44] = 3;

    auto NC = mesh.number_of_cells();
    std::vector<int> zdata(NC);

    for(int i = 0; i < NC; i++)
    {
      zdata[i] = mesh.cell_barycenter(i)[2]*10000;
    }
    Writer writer;
    writer.set_points(mesh);
    writer.set_cells(mesh);
    writer.set_cell_data(zdata, 1, "z");
    writer.set_point_data(gtag, 1, "gtag");
    writer.set_point_data(gdof, 1, "gdof");
    writer.write("hex_test.vtu");
    return 0;
}
