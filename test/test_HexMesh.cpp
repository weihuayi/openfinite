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

int main()
{
    HexMesh mesh;

    MF::cube_hexahedron_mesh(mesh);
    
    std::vector<double> q;
    mesh.uniform_refine(5);

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
    writer.write("hex_test.vtu");
    return 0;
}
