
/**
 * \file BDMdof.h
 * \brief BDM 有限元自由度管理
 * \author 陈春雨
 * \data 2022/3/16
 */
#ifndef BDMDof_h
#define BDMDof_h

#include <math.h>
#include <memory>
#include <vector>
#include <array>
#include <map>
#include <numeric>
#include <time.h>

#include <functional>

#include "RArray2d.h"
#include "femcore.h"

namespace OF{
namespace Core{

template<int N>
struct DofMap
{
  RArray2d<int> index;

  DofMap(int p): index() 
  {
    RArray2d<int> multiindex;
    get_all_multiindex(p, N+1, multiindex);

    int ntype = A(N+1, N+1);
    int ldof = multiindex.shape(0);
    index.init(ntype, ldof);

    int type[N+1];
    for(int i = 0; i < N+1; i++)
      type[i] = i;

    for(int i = 0; i < ntype; i++)
    {
      for(int j = 0; j < ldof; j++)
        index(i, j) = multiindex_to_number(N+1, p, multiindex[j], type);

      std::next_permutation(type, type+N+1);
    }
  }

  int get_type(const int * E0, const int * E1) const 
  {
    int type[N+1] = {-1};
    for(int i = 0; i < N+1; i++)
    {
      for(int j = 0; j < N+1; j++)
      {
        if(E0[i] == E1[j])
        {
          type[i] = j;
          break;
        }
      }
    }
    return permutation_to_number(type, N+1);
  }
};

template<int N>
struct DofClassify
{
  RArray2d<int> boundary;
  RArray2d<int> bdof; // 每个边界上的自由度的编号
  std::vector<int> indof; // 内部自由度的编号
  RArray2d<int> fip; // 每个面内部插值点的编号
  RArray2d<int> eip; // 每个边内部插值点的编号
  std::vector<int> vip; // 每个顶点上插值点编号

  DofClassify(int p): boundary(), bdof(), fip(), eip(), vip(0), indof(0)
  {
    RArray2d<int> multiindex;
    get_all_multiindex(p, N+1, multiindex);

    int l0 = C(N-1+p, N-1), ldof = multiindex.shape(0), l1 = ldof*3 - l0*4;
    boundary.init(N+1, N);
    bdof.init(N+1, l0);
    fip.init(N+1, l0 - p*3);
    eip.init(6, p-1);
    vip.reserve(4);
    indof.reserve(l1);

    /** 构建 boundary */
    for(int i = 0; i < N+1; i++)
    {
      int k = 0;
      for(int j = 0; j < N+1; j++)
      {
        if(j != i)
          boundary(i, k++) = j;
      }
    }

    /** 构建 bdof, indof */
    int l[N+1] = {0}; // l[i] 是当前第 i 个面在已经编号的自由度的个数. 
    for(int i = 0; i < ldof; i++)
    {
      int f = 0; // 当前插值点周围面的个数
      for(int j = N; j > -1; j--)
      {
        if(multiindex(i, j)==0)
        {
          bdof(j, l[j]++) = i + ldof*(f++);
        }
      }
      for(int j = f; j < 3; j++)
        indof.push_back(i + ldof*j);
    }

    /** get ipoint in face, edge or vertex */
    int le[6] = {0}, lf[4] = {0};
    for(int i = 0; i < ldof; i++)
    {
      int d = 0, idx[4] = {0};
      for(int j = 0; j < 4; j++)
      {
        if(multiindex(i, j)==0)
          idx[d++] = j;
      }
      if(d == 2)
      {
        int e = 5 - idx[0]- idx[1];
        if(idx[0]==0)
          e++;
        eip(e, le[e]++) = i;
      }
      else if(d == 3)
        vip.push_back(i);
      else if(d==1)
      {
        int f = idx[0];
        fip(f, lf[f]++) = i;
      }
    }
  }
};

template<typename Mesh>
class BDMDof3d
{
public:
  BDMDof3d(int p, std::shared_ptr<Mesh> mesh): m_p(p), 
    m_mesh(mesh), m_faceMap(p), m_celldof(p) 
  {
    get_all_multiindex(p, 4, multiindex3);
    get_all_multiindex(p, 3, multiindex2);
  }
  
  void face_to_dof(int f, int * dof)
  {
    int nfdof = number_of_local_dofs(2);
    int N = f*nfdof;
    for(int i = 0; i < nfdof; i++)
    {
      dof[i] = N+i;
    }
  }

  void cell_to_dof(int c, int * dof)
  {
    int NF = m_mesh->number_of_faces();
    int nfldof = number_of_local_dofs(2);
    int ncldof = number_of_local_dofs(3)*3 - nfldof*4;

    const int * cell = m_mesh->cell(c);

    const auto & em = m_faceMap;
    const auto & idof = m_celldof.indof;
    const auto & bdof = m_celldof.bdof;
    const auto & bd = m_celldof.boundary; 

    int fdof[nfldof];
    for(int i = 0; i < 4; i++)
    {
      int f = m_mesh->cell_to_face(c)[i];
      face_to_dof(f, fdof);

      int face[3] = {cell[bd(i, 0)], cell[bd(i, 1)], cell[bd(i, 2)]};
      int t = em.get_type(face, m_mesh->face(f));
      const int * idx = em.index[t];

      for(int j = 0; j < nfldof; j++)
      {
        dof[bdof(i, j)] = fdof[idx[j]];
      }
    }
    for(int i = 0; i < (int)idof.size(); i++)
      dof[idof[i]] = NF*nfldof + c*ncldof + i;
  }

  int number_of_local_dofs(int d = 3)
  {
    return C(d+m_p, d);
  }

  int number_of_dofs()
  {
    int NF = m_mesh->number_of_faces();
    int NC = m_mesh->number_of_cells();
    int nfdof = C(2+m_p, 2), ncdof = C(m_p+3, 3)*4 - nfdof*4;
    return NF*nfdof + NC*ncdof;
  }

  void print()
  {
    auto & bou = m_celldof.boundary;
    auto & bdof = m_celldof.bdof;
    auto & indof = m_celldof.indof;
    auto & fdof = m_celldof.fip;
    auto & edof = m_celldof.eip;
    auto & vdof = m_celldof.vip;
    std::cout<< bou << std::endl;
    std::cout<< bdof << std::endl;
    std::cout<< "fdof: " <<  fdof << std::endl;
    std::cout<< edof << std::endl;
    for(auto & bb : indof)
      std::cout<< bb << std::endl;

    std::cout<< "vdof: " << std::endl;
    for(auto & bb : vdof)
      std::cout<< bb << std::endl;
  }

  void basis_vector(int c, RArray2d<double> & b2v)
  {
    RArray2d<double> c2n;
    RArray2d<double> c2t;
    m_mesh->cell_normal(c, c2n);
    m_mesh->cell_tangent(c, c2t);
    const auto & lf2e = m_mesh->localface2edge;
    const auto & bd = m_celldof.boundary;

    int ldof = multiindex3.shape(0);
    b2v.init(ldof*3, 3);

    /** 单元内部基函数的方向 */
    for(int i = 0; i < ldof; i++)
    {
      for(int j = 0; j < 3; j++)
        b2v(j*ldof+i, j) = 1;
    }

    /** 面内部的基函数的方向 */
    const auto & fip = m_celldof.fip;
    for(int i = 0; i < 4; i++)
    {
      double v[3] = {0};
      OF::Core::cross(c2n[i], c2t[lf2e[i][0]], v);
      for(int j = 0; j < fip.shape(1); j++)
      {
        std::copy_n(c2n[i], 3, b2v[fip(i, j)]);
        std::copy_n(c2t[lf2e[i][0]], 3, b2v[fip(i, j)+ldof]);
        std::copy_n(v, 3, b2v[fip(i, j)+ldof*2]);
      }
    }

    /** 边内部的基函数的方向 */
    const auto & le2f = m_mesh->localedge;
    const auto & eip = m_celldof.eip;
    for(int i = 0; i < 6; i++)
    {
      RArray2d<double> v(2, 3);
      OF::Core::cross(c2n[le2f[5-i][0]], c2t[i], v[0]);
      OF::Core::cross(c2n[le2f[5-i][1]], c2t[i], v[1]);
      double val = OF::Core::dot(v[0], c2n[le2f[5-i][1]]);
      for(int j = 0; j < 3; j++)
      {
        v(0, j) /= val; 
        v(1, j) /= -val; 
      }
      for(int j = 0; j < eip.shape(1); j++)
      {
        std::copy_n(v[0], 3, b2v[eip(i, j)]);
        std::copy_n(v[1], 3, b2v[eip(i, j)+ldof]);
        std::copy_n(c2t[i], 3, b2v[eip(i, j)+ldof*2]);
      }
    }

    /** 顶点上的基函数的方向 */
    int lv2e[4][3] = {{2, 1, 0}, {4, 3, 0}, {5, 3, 1}, {5, 4, 2}};
    const auto & vip = m_celldof.vip;
    for(int i = 0; i < 4; i++)
    {
      double val0 = OF::Core::dot(c2t[lv2e[i][0]], c2n[bd[i][2]]);
      double val1 = OF::Core::dot(c2t[lv2e[i][1]], c2n[bd[i][1]]);
      double val2 = OF::Core::dot(c2t[lv2e[i][2]], c2n[bd[i][0]]);
      for(int j = 0; j < 3; j++)
      {
        b2v(vip[i], j) = c2t(lv2e[i][0], j)/val0;
        b2v(vip[i]+ldof, j) = c2t(lv2e[i][1], j)/val1;
        b2v(vip[i]+ldof*2, j) =c2t(lv2e[i][2], j)/val2;
      }
    }
  }

  void dof_vector(int c, RArray2d<double> & d2v)
  {
    RArray2d<double> normal(4, 3);
    RArray2d<double> tangent(6, 3);
    m_mesh->cell_normal(c, normal);
    m_mesh->cell_tangent(c, tangent);
    const auto & lf2e = m_mesh->localface2edge;

    RArray2d<double> c2n(4, 3);
    RArray2d<double> c2t(6, 3);
    m_mesh->cell_normal(c, c2n);
    m_mesh->cell_tangent(c, c2t);

    /** 单元内部自由度的方向 */
    int ldof = multiindex3.shape(0);
    d2v.init(ldof*3, 3);
    for(int i = 0; i < ldof; i++)
    {
      for(int j = 0; j < 3; j++)
        d2v(j*ldof+i, j) = 1;
    }

    /** 面上的法相自由度 */
    const auto & bdof = m_celldof.bdof;
    for(int i = 0; i < 4; i++)
    {
      for(int j = 0; j < bdof.shape(1); j++)
        std::copy_n(c2n[i], 3, d2v[bdof(i, j)]);
    }

    /** 面内部的自由度的方向 */
    const auto & fip = m_celldof.fip;
    for(int i = 0; i < 4; i++)
    {
      double v[3] = {0};
      OF::Core::cross(c2n[i], c2t[lf2e[i][0]], v);
      for(int j = 0; j < fip.shape(1); j++)
      {
        std::copy_n(c2t[lf2e[i][0]], 3, d2v[fip(i, j)+ldof]);
        std::copy_n(v, 3, d2v[fip(i, j)+ldof*2]);
      }
    }

    /** 边内部的自由度的方向 */
    const auto & eip = m_celldof.eip;
    for(int i = 0; i < 6; i++)
    {
      for(int j = 0; j < eip.shape(1); j++)
        std::copy_n(c2t[i], 3, d2v[eip(i, j)+ldof*2]);
    }
  }

private:
  int m_p;
  std::shared_ptr<Mesh> m_mesh;
  RArray2d<int> multiindex2; 
  RArray2d<int> multiindex3; 
  DofMap<2> m_faceMap;
  DofClassify<3> m_celldof;
};

}//end of Core
}//end of OF
#endif // end of BDMDof_h
