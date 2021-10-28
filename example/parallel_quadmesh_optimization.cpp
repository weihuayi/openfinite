#include <fstream>
#include <string>
#include <iostream>
#include <list>
#include <vector>
#include <set>

#include <mpi.h>
#include <time.h>

#include "geometry/Geometry_kernel.h"
#include "geometry/RectangleWithHole.h"
#include "geometry/RectangleWithTwoHoles.h"
#include "mesh/QuadrilateralMesh.h"
#include "mesh/ParallelMesh.h"
#include "mesh/ParallelMesher.h"
#include "mesh/VTKMeshReader.h"
#include "mesh/VTKMeshWriter.h"
#include "mesh/GhostFillingAlg.h"
#include "mesh/ParallelMeshOptAlg.h"
#include "mesh/QuadJacobiQuality.h"
#include "mesh/QuadPositiveJacobiQuality.h"
#include "mesh/SumNodePatchObjectFunction.h"
#include "mesh/MaxNodePatchObjectFunction.h"
#include "mesh/MeshFactory.h"
#include "Python.h"

typedef OF::Geometry_kernel<double, int> GK;
typedef OF::GeometryModel::RectangleWithHole<GK> Model;
//typedef OF::GeometryModel::RectangleWithTwoHoles<GK> Model;
typedef GK::Point_2 Node;
typedef GK::Vector_2 Vector;
typedef OF::Mesh::QuadrilateralMesh<GK, Node, Vector> QMesh;
typedef OF::Mesh::ParallelMesh<GK, QMesh> PMesh;
//typedef OF::Mesh::QuadJacobiQuality<PMesh> CellQuality;
typedef OF::Mesh::QuadPositiveJacobiQuality<PMesh> CellQuality;
typedef OF::Mesh::SumNodePatchObjectFunction<PMesh, CellQuality> ObjectFunction;
//typedef OF::Mesh::MaxNodePatchObjectFunction<PMesh, CellQuality> ObjectFunction;
typedef OF::Mesh::ParallelMesher<PMesh> PMesher;
typedef OF::Mesh::ParallelMeshOptAlg<PMesh, ObjectFunction, Model> PMeshOpt;
typedef OF::Mesh::VTKMeshWriter Writer;
typedef OF::Mesh::VTKMeshReader<PMesh> Reader;

template<typename I>
void plot(std::vector<I> & data0, std::vector<I> & data1)
{
  int N0 = data0.size();
  int N1 = data1.size();

  Py_Initialize();
  PyRun_SimpleString("import sys");
  PyRun_SimpleString("sys.path.append('../example')");

  PyObject* pModule = PyImport_ImportModule("plot");
  PyObject* pFunc = PyObject_GetAttrString(pModule, "Histogram_plot");//要运行的函数

  PyObject* plist0 = PyList_New(N0);//函数的参数是一个list
  PyObject* ptuple0 = PyTuple_New(1);//把参数用 tuple 装起来
  PyObject* plist1 = PyList_New(N1);//函数的参数是一个list
  PyObject* ptuple1 = PyTuple_New(1);//把参数用 tuple 装起来

  for(int i = 0; i < N0; i++)
  {
    PyObject*  pra = Py_BuildValue("d", data0[i]);
    PyList_SetItem(plist0, i, pra);
  }

  for(int i = 0; i < N1; i++)
  {
    PyObject*  pra = Py_BuildValue("d", data1[i]);
    PyList_SetItem(plist1, i, pra);
  }

  PyTuple_SetItem(ptuple0, 0, plist0);
  PyTuple_SetItem(ptuple1, 0, plist1);

	PyObject_CallObject(pFunc, ptuple0);//运行函数
	PyObject_CallObject(pFunc, ptuple1);//运行函数
	Py_Finalize();
}

int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);

  int rank, nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  PMesher pmesher(argv[1], ".vtu", MPI_COMM_WORLD);
  auto mesh = pmesher.get_mesh();

  auto quad = std::make_shared<Model>();

  auto NC = mesh->number_of_cells();
  std::vector<double> cellQualityInit(NC);
  std::vector<double> cellQualityOpt(NC);

  CellQuality mq(mesh);
  mq.quality_of_mesh(cellQualityInit);

  PMeshOpt optAlg(mesh, quad, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
  auto start = MPI_Wtime();
  optAlg.optimization(1e-5, 100);//优化
  MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
  auto end = MPI_Wtime();

  mq.quality_of_mesh(cellQualityOpt);

  // 数据展示 
  if(mesh->id() != 0)
  {
    MPI_Send(&NC, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(cellQualityInit.data(), NC, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    MPI_Send(cellQualityOpt.data(), NC, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
    MPI_Send(mesh->cell_global_id().data(), NC, MPI_INT, 0, 3, MPI_COMM_WORLD);
  }
  if(mesh->id() == 0)
  {
    int NAC = NC;
    int ND[nprocs];
    ND[0] = NC;
    for(int i = 1; i < nprocs; i++)
    {
      MPI_Recv(&ND[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      NAC += ND[i];
    }
    std::vector<double> initdata(NAC); 
    std::vector<double> optdata(NAC); 

    NAC = 0;
    for(int i = 0; i < nprocs; i++)
    {
      std::vector<int> cgid(ND[i]);
      std::vector<double> init_i(ND[i]);
      std::vector<double> opt_i(ND[i]);
      if(i == 0)
      {
        cgid = mesh->cell_global_id();
        init_i = cellQualityInit;
        opt_i = cellQualityOpt;
      }
      else
      {
        MPI_Recv(init_i.data(), ND[i], MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(opt_i.data(), ND[i], MPI_INT, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(cgid.data(), ND[i], MPI_INT, i, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      
      for(int j = 0; j < ND[i]; j++)
      {
        initdata[cgid[j]] = init_i[j];
        optdata[cgid[j]] = opt_i[j];
      }

      int N_i = *std::max_element(cgid.begin(), cgid.end());
      NAC = NAC>N_i?NAC:N_i;
    }
    initdata.resize(NAC);
    optdata.resize(NAC);

    double i_min_q = * std::min_element(initdata.begin(), initdata.end());
    double i_max_q = * std::max_element(initdata.begin(), initdata.end());
    double i_ave_q = std::accumulate(initdata.begin(), initdata.end(), 0)/NAC;
    double o_min_q = * std::min_element(optdata.begin(), optdata.end());
    double o_max_q = * std::max_element(optdata.begin(), optdata.end());
    double o_ave_q = std::accumulate(optdata.begin(), optdata.end(), 0)/NAC;

    std::cout<< "总运行时间:" << end-start << std::endl;
    std::cout<< "初始网格单元最差质量为: " << i_max_q <<std::endl; 
    std::cout<< "初始网格单元最优质量为: " << i_min_q <<std::endl; 
    std::cout<< "初始网格单元平均质量为: " << i_ave_q <<std::endl; 
    std::cout<< "优化网格单元最差质量为: " << o_max_q <<std::endl; 
    std::cout<< "优化网格单元最优质量为: " << o_min_q <<std::endl; 
    std::cout<< "优化网格单元平均质量为: " << o_ave_q <<std::endl; 

    std::ofstream wfile;
    wfile.open("quad_test.txt", ios::app);
    wfile<< argv[1] << "\n\n";
    wfile<< "总运行时间为: " << end-start << "\n";
    wfile<< "单元数目为: " << NAC << "\n";
    wfile<< "进程数为: " << nprocs << "\n";
    wfile<< "\n";
    wfile<< "初始网格单元最差质量为: " << i_max_q << "\n";
    wfile<< "初始网格单元最优质量为: " << i_min_q << "\n";
    wfile<< "初始网格单元平均质量为: " << i_ave_q << "\n";
    wfile<< "\n";
    wfile<< "优化网格单元最差质量为: " << o_max_q << "\n";
    wfile<< "优化网格单元最优质量为: " << o_min_q << "\n";
    wfile<< "优化网格单元平均质量为: " << o_ave_q << "\n";
    wfile<< "\n\n\n\n";
    wfile.close();
    //plot(cellQualityInit, cellQualityOpt);
  }
  std::stringstream ss;
  ss << "opt_"<< mesh->id() << ".vtu";
  std::vector<int> qInit(NC);                                                   
  std::vector<int> qOpt(NC);                                                    
  for(int i = 0; i < NC; i++)                                                   
  {                                                                             
    qInit[i] = (int)(cellQualityInit[i]*1000);                                  
    qOpt[i] = (int)(cellQualityOpt[i]*1000);                                    
  }                                                                             

  Writer writer;                                                                
  writer.set_points(*mesh);                                                     
  writer.set_cells(*mesh);                                                      
  writer.set_point_data(mesh->get_node_int_data()["nid"], 1, "nid");            
  writer.set_cell_data(qOpt, 1, "q_opt");                                       
  writer.set_cell_data(qInit, 1, "q_init");                                     
  writer.set_cell_data(cellQualityOpt, 1, "opt");                               
  writer.set_cell_data(cellQualityInit, 1, "init");                             
  writer.write(ss.str());                        
  MPI_Finalize();
  return 0;
}
