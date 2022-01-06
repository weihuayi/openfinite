#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <chrono>



#include "geometry/Geometry_kernel.h"
#include "mesh/TetrahedronMesh.h"
#include "mesh/VTKMeshWriter.h"



typedef OF::Geometry_kernel<double, int> GK;
typedef GK::Point_3 Node;
typedef GK::Vector_3 Vector;

typedef OF::Mesh::TetrahedronMesh<GK, Node, Vector> TetMesh;
typedef TetMesh::Cell Cell;

typedef OF::Mesh::VTKMeshWriter VMWriter;



/**
 * \brief 生成一个长方体上四面体网格, 并指定固定的点
 *
 * \param[in] box [xmin, xmax, ymin, ymax, zmin, zmax]
 * \param[in] nx x 方向剖分的段数
 * \param[in] ny y 方向剖分的段数 
 * \param[in] nz z 方向剖分的段数
 * \param[out] nodes 四面体网格节点
 * \param[out] cells 四面体网格单元
 * \param[out] nidxs 网格中的固定点编号
 */
void boxmesh3d(
    std::array<double, 6> & box, 
    int nx, 
    int ny, 
    int nz, 
    std::vector<double> & nodes, 
    std::vector<int> & cells, 
    std::vector<int> & nidxs
    )
{

  // 1. 生成网格节点坐标
  std::cout << "开始构建网格节点......" << std::endl;
  auto start = std::chrono::steady_clock::now();

  double hx = (box[1] - box[0])/nx;
  double hy = (box[3] - box[2])/ny;
  double hz = (box[5] - box[4])/nz;

  auto NN = (nx+1)*(ny+1)*(nz+1);
  nodes.resize(3*NN);
  int idx=0;
  for(int i = 0; i < nx+1; i++)
  {
    for(int j = 0; j < ny+1; j++)
    {
      for(int k = 0; k < nz+1; k++)
      {
        nodes[3*idx + 0] = i*hx;
        nodes[3*idx + 1] = j*hy;
        nodes[3*idx + 2] = k*hz;
        if(i == 0 || i == nx)
          nidxs.push_back(idx);
        idx++;
      }
    }
  }

  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "节点构建完成，花费时间: "
       << elapsed_seconds.count() << " s\n" << std::endl;

  std::cout << "开始构建四面体网格单元......" << std::endl;
  start = std::chrono::steady_clock::now();
  idx=0; // 计数器置 0 
  std::vector<std::array<int, 8> > hexcell(nx*ny*nz);
  for(int i = 0; i < nx; i++)
  {
    for(int j = 0; j < ny; j++)
    {
      for(int k = 0; k < nz; k++)
      {
        hexcell[idx][0] = i*(ny+1)*(nz+1) + j*(nz+1) + k;
        hexcell[idx][1] = hexcell[idx][0] + (ny + 1)*(nz + 1);
        hexcell[idx][2] = hexcell[idx][1] + nz + 1;
        hexcell[idx][3] = hexcell[idx][0] + nz + 1;

        hexcell[idx][4] = hexcell[idx][0] + 1;
        hexcell[idx][5] = hexcell[idx][4] + (ny + 1)*(nz + 1);
        hexcell[idx][6] = hexcell[idx][5] + nz + 1;
        hexcell[idx][7] = hexcell[idx][4] + nz + 1;
        idx++;
      }
    }
  }

  std::vector<std::array<int, 4> > localCell;
  localCell.push_back({0, 1, 2, 6});
  localCell.push_back({0, 5, 1, 6});
  localCell.push_back({0, 4, 5, 6});
  localCell.push_back({0, 7, 4, 6});
  localCell.push_back({0, 3, 7, 6});
  localCell.push_back({0, 2, 3, 6});

  cells.resize(24*hexcell.size());
  idx = 0; // 计数器置 0 
  for(auto & c : hexcell)
  {
    for(int i = 0; i < 6; i++)
    {
      cells[4*idx + 0] = c[localCell[i][0]];
      cells[4*idx + 1] = c[localCell[i][1]];
      cells[4*idx + 2] = c[localCell[i][2]];
      cells[4*idx + 3] = c[localCell[i][3]];
      idx++;
    }
  }

  end = std::chrono::steady_clock::now();
  elapsed_seconds = end-start;
  std::cout << "单元构建完成，花费时间: "
       << elapsed_seconds.count() << " s\n" << std::endl;
}

/**
 * 这里人工构造了一个简单的模型, 一根梁长 10 m, 宽 1 m, 高 2 m,
 * 上表面承受 6700 Pa 的压力, 两端固定, 并考虑重力作用
 *
 */
int main(int argc, char **argv)
{
  std::string fname = "TetMesh.vtu";
  VMWriter writer;
  return 0;
}
