/**
 * \file Dof.h
 * \author ccy
 * \date 2021/10/04
 * \brief Lagrange 有限元空间自由度管理
 */
#ifndef Dof_h
#define Dof_h

#include <time.h>
#include <map>
#include <array>
#include <memory>
#include <vector>
#include <initializer_list>
#include <chrono>

#include "common/Multindex.h"
#include "common/Permutation.h"
#include "common/IncreaseSequence.h"

#include "algebra/Matrix.h"

namespace OF{
namespace FunctionSpace{

typedef OF::Common::MultiIndex MultiIndex;
typedef OF::Common::Permutation Permutation;
typedef OF::Common::IncreaseSequence IncreaseSequence;
typedef OF::AlgebraObject::Matrix<int, int> Matrix;

using namespace std::chrono;
/**
 * \brief 子单形上的自由度.
 */
class SubsimplxDof
{
public:
  SubsimplxDof() {}

  SubsimplxDof(int n, int p)
  {
    dof = new Matrix[n+1]; 
    for(int i = 0; i < n+1; i++)
    {
      int N0 = C(n+1, i+1); /**< i 维单形的个数 */
      int N1 = C(p-1, p-i-1); /**< 每个 i 维单形上的 dof 个数 */
      dof[i] = Matrix(N0, N1);
    }
  }

  ~SubsimplxDof()
  {
    delete [] dof;
  }

  Matrix & operator[](const int & i)
  {
    return dof[i];
  }

public:
  Matrix * dof;
};

/**
 * \brief Lagrange 有限元空间自由度管理 
 * \param n 维数
 * \note 要求 Mesh 有两个接口
 * 1. cell_to_subsimplex(int i, int j, int k), 即第 i 个单元的第 k 个 j 维单形的编号.
 */
template<int n, typename Mesh>
class Dof
{
public:
  Dof(int p, std::shared_ptr<Mesh> mesh):m_p(p), m_mesh(mesh)
  {
    int ** num = new int*[n+1]; /**< 每个维度, 每个子单形, 自由度的个数 */
    for(int i = 0; i < n+1; i++)
      num[i] = new int[C(n+1, i+1)]();

    int ldof = number_of_subsimplexdofs(n);
    m_dofs.resize(ldof, {-1, 0, 0});
    MultiIndex domainPoint(p, n+1);
    for(int k = 0; k < ldof; k++)
    {
      for(int i = 0; i < n+1; i++)
      {
        if(domainPoint[i]>0)
          m_dofs[k][0]++;
      }

      IncreaseSequence subsimplex(m_dofs[k][0]+1);
      int idx = 0;
      for(int i = 0; i < n+1; i++)
      {
        if(domainPoint[i]>0)
          subsimplex[idx++] = i;
      }
      m_dofs[k][1] = subsimplex.to_number(n+1);
      m_dofs[k][2] = num[m_dofs[k][0]][m_dofs[k][1]]++;
      domainPoint.nex_multi_index();
    }
    for(int i = 0; i < n+1; i++)
      delete [] num[i];
    delete [] num;

    for(int i = 0; i < (n>p?p:n); i++)
    {
      int N0 = A(i+1, i+1);
      int N1 = C(p-1, p-i-1);
      m_idx[i] = Matrix(N0, N1, 0); 
      Permutation idx(i+1);
      std::vector<MultiIndex> midxs;
      MultiIndex::gen_all_index(p-i-1, i+1, midxs);
      for(int j = 0; j < N0; j++)
      {
        for(int k = 0; k < N1; k++)
        {
          m_idx[i][j][k] = midxs[k].to_number(idx);
        }
        idx.nex_permutation();
      }
    }
  }

  void cell_to_dof(int c, std::vector<int> & c2d)
  {
    auto start = high_resolution_clock::now();

    const auto & cell = m_mesh->cell(c);
    SubsimplxDof ssdof(n, m_p);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<nanoseconds>(stop - start);
    std::cout << "1时间: "<< duration.count() << std::endl;

    start = high_resolution_clock::now();

    int N[4];
    N[0] = m_mesh->number_of_nodes();
    N[1] = m_mesh->number_of_edges(); 
    N[2] = m_mesh->number_of_faces(); 
    N[3] = m_mesh->number_of_cells(); 

    for(int j = 0; j < ssdof[0].shape[0]; j++)
      ssdof[0][j][0] = cell[j];

    stop = high_resolution_clock::now();
    duration = duration_cast<nanoseconds>(stop - start);
    std::cout << "2时间: "<< duration.count() << std::endl;

    start = high_resolution_clock::now();

    int ND = N[0];
    for(int i = 1; i < n+1; i++)
    {
      int N0 = ssdof[i].shape[0];
      int N1 = ssdof[i].shape[1];
      for(int j = 0; j < N0; j++)
      {
        int ii = m_mesh->cell_to_subsimplex(c, i, j); /**< 单形编号 */
        for(int k = 0; k < N1; k++)
        {
          ssdof[i][j][k] = ND+ii*N1+k;
        }
      }
      ND += N[i]*N1;
    }

    stop = high_resolution_clock::now();
    duration = duration_cast<nanoseconds>(stop - start);
    std::cout << "3时间: "<< duration.count() << std::endl;

    start = high_resolution_clock::now();

    /** 判断各单行旋转情况 */
    int ** pertype = new int*[n-1];
    for(int i = 0; i < n-1; i++)
    {
      int N0 = C(n+1, i+2);
      pertype[i] = new int[N0];
      std::vector<IncreaseSequence> incseq;
      IncreaseSequence::gen_all_sequence(n+1, i+2, incseq);
      for(int j = 0; j < N0; j++)
      {
        auto sid = m_mesh->cell_to_subsimplex(c, i+1, j);
        std::vector<int> glospx;
        m_mesh->subsimplex(i+1, sid, glospx);/**< 一个单形的节点表示 */
        std::vector<int> locspx(i+2);
        for(int k = 0; k < i+2; k++)
          locspx[k] = cell[incseq[j][k]];
        pertype[i][j] = Permutation::chack(glospx, locspx, i+2).to_number();
      }
    }

    stop = high_resolution_clock::now();
    duration = duration_cast<nanoseconds>(stop - start);
    std::cout << "4时间: "<< duration.count() << std::endl;

    start = high_resolution_clock::now();

    /** 按照规则填进数组 */
    int ldof = number_of_subsimplexdofs(n);
    c2d.resize(ldof);

    stop = high_resolution_clock::now();
    duration = duration_cast<nanoseconds>(stop - start);
    std::cout << "5时间: "<< duration.count() << std::endl;

    start = high_resolution_clock::now();

    int * sssss = new int[ldof]();

    stop = high_resolution_clock::now();
    duration = duration_cast<nanoseconds>(stop - start);
    std::cout << "6时间: "<< duration.count() << std::endl;

    start = high_resolution_clock::now();

    for(int i = 0; i < ldof; i++)
    {
      auto & dof = m_dofs[i];
      if(dof[0]==0 || dof[0]==n)
      {
        c2d[i] = ssdof[dof[0]][dof[1]][dof[2]];
      }
      else
      {
        int ptype = pertype[dof[0]-1][dof[1]];
        const auto & idx = m_idx[dof[0]][ptype];
        c2d[i] = ssdof[dof[0]][dof[1]][idx[dof[2]]]; 
      }
    }

    stop = high_resolution_clock::now();
    duration = duration_cast<nanoseconds>(stop - start);
    std::cout << "7时间: "<< duration.count() << std::endl;

    start = high_resolution_clock::now();

    for(int i = 0; i < n-1; i++)
      delete [] pertype[i];
    delete [] pertype;

    stop = high_resolution_clock::now();
    duration = duration_cast<nanoseconds>(stop - start);
    std::cout << "8时间: "<< duration.count() << std::endl;
  }

  int number_of_subsimplexdofs(const int & d)
  {
    return C(d+m_p, m_p);
  }

  void print()
  {
    int L = m_dofs.size();
    std::cout<< "dof = " <<std::endl;
    for(int i = 0; i < L; i++)
    {
      std::cout<< "[" << m_dofs[i][0] << ", " << m_dofs[i][1] << ", " << m_dofs[i][2] << "]" << std::endl;
    }
    std::cout<< "idx = " <<std::endl;
    for(int i = 0; i < n; i++)
    {
      std::cout<< m_idx[i] <<std::endl;
    }
  }

private:
  int m_p;
  std::vector<std::array<int, 3>> m_dofs;/**< [单形维数, 单形编号, dof 在单形中的编号] */
  std::shared_ptr<Mesh> m_mesh;
  Matrix m_idx[n]; /**< 每一维, 不同顶点排序的单形, 单行上的域点的对应关系 */
};

} // end of namespace FunctionSpace 

} // end of namespace OF

#endif // end of Dof_h
