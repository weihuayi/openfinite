
/**
 * \file NedelecDof.h
 * \brief Nedelec 有限元自由度管理
 * \author 陈春雨
 * \data 2022/3/19
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

struct FaceDofMap
{
  RArray2d<int> index;

  FaceDofMap(int p): index() 
  {
    RArray2d<int> multiindex;
    get_all_multiindex(p, 3, multiindex);

    int ntype = A(3, 3);
    int ldof = multiindex.shape(0);
    index.init(ntype, ldof*2);

    int type[3];
    for(int i = 0; i < 3; i++)
      type[i] = i;

    int idx[3] = {0, ldof-1-p, ldof-1};
    int edge[3][2] = {{1, 2}, {0, 2}, {0, 1}};
    for(int i = 0; i < ntype; i++)
    {
      for(int j = 0; j < ldof; j++)
      {
        index(i, j) = multiindex_to_number(3, p, multiindex[j], type);
        index(i, j+ldof) = index(i, j)+ldof;
      }
      for(int j = 0; j < 3; j++)
      {
        if(type[edge[j][0]] > type[edge[j][1]])
        {
          index(i, idx[j]) += ldof;
          index(i, idx[j]+ldof) -= ldof;
        }
      }
      std::next_permutation(type, type+3);
    }
  }

  int get_type(const int * E0, const int * E1) const 
  {
    int type[3] = {-1};
    for(int i = 0; i < 3; i++)
    {
      for(int j = 0; j < 3; j++)
      {
        if(E0[i] == E1[j])
        {
          type[i] = j;
          break;
        }
      }
    }
    return permutation_to_number(type, 3);
  }
};

struct EdgeDofMap
{
  RArray2d<int> index;

  EdgeDofMap(int p): index() 
  {
    RArray2d<int> multiindex;
    get_all_multiindex(p, 2, multiindex);

    int ldof = multiindex.shape(0);
    index.init(2, ldof);

    for(int j = 0; j < ldof; j++)
    {
      index(0, j) = j;
      index(1, j) = ldof-1-j;
    }
  }

  int get_type(const int * E0, const int * E1) const 
  {
    return E0[0]==E1[0]? 0 : 1;
  }
};

struct InterplationPointInCell
{
  RArray2d<int> fip; // 每个面内部插值点的编号
  RArray2d<int> eip; // 每个边内部插值点的编号
  std::vector<int> vip; // 每个顶点上插值点编号
  InterplationPointInCell(int p): fip(), eip(), vip(0)
  {
    RArray2d<int> multiindex;
    get_all_multiindex(p, 4, multiindex);
    int ldof = multiindex.shape(0);

    int l0 = C(2+p, 2);
    fip.init(4, l0 - p*3);
    eip.init(6, p-1);
    vip.reserve(4);

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
  };
};

struct CellDofClassify
{
  RArray2d<int> boundary;
  RArray2d<int> bdof; // 每个边界上的自由度的编号
  std::vector<int> indof; // 内部自由度的编号

  CellDofClassify(int p): boundary(), bdof(), indof(0)
  {
    RArray2d<int> multiindex;
    get_all_multiindex(p, 4, multiindex);

    int l0 = C(2+p, 2), ldof = multiindex.shape(0);
    boundary.init(4, 3);
    bdof.init(4, l0*2);
    indof.reserve(ldof*3);

    /** 构建 boundary */
    for(int i = 0; i < 4; i++)
    {
      int k = 0;
      for(int j = 0; j < 4; j++)
      {
        if(j != i)
          boundary(i, k++) = j;
      }
    }

    std::vector<bool> isbdof(ldof*3, 0);
    std::cout<< "aaa" <<std::endl;
    int m[3][2] = {{0, 1}, {0, 2}, {1, 2}};
    /** 构建 bdof, indof */
    int l[4] = {0}; // l[i] 是当前第 i 个面在已经编号的自由度的个数. 
    for(int i = 0; i < ldof; i++)
    {
      int d = 0, idx[4] = {0};/**< 3-d: 所在单形的维度; idx: 相邻面编号. */
      for(int j = 0; j < 4; j++)
      {
        if(multiindex(i, j)==0)
          idx[d++] = j;
      }
      for(int j = 0; j < d; j++)
      {
        bdof(idx[d-1-j], l[idx[d-1-j]]) = i + ldof*m[j][0];
        bdof(idx[d-1-j], (l[idx[d-1-j]]++)+l0) = i + ldof*m[j][1];
        isbdof[i + ldof*m[j][0]] = true;
        isbdof[i + ldof*m[j][1]] = true;
      }
    }
    for(int i = 0; i < ldof*3; i++)
    {
      if(!isbdof[i])
        indof.push_back(i);
    }
  }
};

struct FaceDofClassify
{
  RArray2d<int> boundary;
  RArray2d<int> bdof;
  std::vector<int> indof;

  FaceDofClassify(int p): boundary(), bdof(), indof(0)
  {
    RArray2d<int> multiindex;
    get_all_multiindex(p, 3, multiindex);

    int l0 = p+1, ldof = multiindex.shape(0);
    boundary.init(3, 2);
    bdof.init(3, l0);
    indof.reserve(ldof*2);

    /** 构建 boundary */
    for(int i = 0; i < 3; i++)
    {
      int k = 0;
      for(int j = 0; j < 3; j++)
      {
        if(j != i)
          boundary(i, k++) = j;
      }
    }

    /** 构建 bdof, indof */
    std::vector<bool> isbdof(ldof*2, false);
    int l[3] = {0}; // l[i] 是当前第 i 个面在已经编号的自由度的个数. 
    for(int i = 0; i < ldof; i++)
    {
      int d = 0, idx[3] = {0};/**< 2-d: 所在单形的维度; idx: 相邻边编号. */
      for(int j = 0; j < 3; j++)
      {
        if(multiindex(i, j)==0)
          idx[d++] = j;
      }
      for(int j = 0; j < d; j++)
      {
        isbdof[i + ldof*j] = true;
        bdof(idx[d-1-j], l[idx[d-1-j]]++) = i + ldof*j;
      }
    }
    for(int i = 0; i < ldof*2; i++)
    {
      if(!isbdof[i])
        indof.push_back(i);
    }
  }
};

template<typename Mesh>
class NedelecDof3d
{
public:
  NedelecDof3d(int p, std::shared_ptr<Mesh> mesh): m_p(p), 
    m_mesh(mesh), m_edgeMap(p), m_faceMap(p), m_facedof(p), m_celldof(p), m_cellip(p) 
  {
    get_all_multiindex(p, 4, multiindex3);
    get_all_multiindex(p, 3, multiindex2);
  }
  
  void edge_to_dof(int e, int * dof)
  {
    int nedof = number_of_local_dofs(1);
    int N = e*nedof;
    for(int i = 0; i < nedof; i++)
      dof[i] = N+i;
  }

  void face_to_dof(int f, int * dof)
  {
    int NE = m_mesh->number_of_edges();
    int nedof = m_p+1;
    int nfdof = number_of_local_dofs(2)*2-nedof*3;
    int N = f*nfdof + NE*nedof;

    const int * face = m_mesh->face(f);

    const auto & em = m_edgeMap;
    const auto & idof = m_facedof.indof;
    const auto & bdof = m_facedof.bdof;
    const auto & bd = m_facedof.boundary; 

    int edof[nedof];
    for(int i = 0; i < 3; i++)
    {
      int e = m_mesh->face_to_edge(f, i);
      edge_to_dof(e, edof);

      int edge[2] = {face[bd(i, 0)], face[bd(i, 1)]};
      const int t = em.get_type(edge, m_mesh->edge(e));
      const int * idx = em.index[t];

      for(int j = 0; j < nedof; j++)
        dof[bdof(i, j)] = edof[idx[j]];
    }
    for(int i = 0; i < (int)idof.size(); i++)
      dof[idof[i]] = N + i;
  }

  void cell_to_dof(int c, int * dof)
  {
    int NE = m_mesh->number_of_edges();
    int NF = m_mesh->number_of_faces();
    int N = NF*(number_of_local_dofs(2)*2-(m_p+1)*3) + NE*(m_p+1);

    int nfldof = number_of_local_dofs(2);
    int ncldof = number_of_local_dofs(3)*3 - nfldof*4;

    const int * cell = m_mesh->cell(c);

    const auto & em = m_faceMap;
    const auto & idof = m_celldof.indof;
    const auto & bdof = m_celldof.bdof;
    const auto & bd = m_celldof.boundary; 

    int fdof[nfldof*2];
    for(int i = 0; i < 4; i++)
    {
      int f = m_mesh->cell_to_face(c)[i];
      face_to_dof(f, fdof);

      int face[3] = {cell[bd(i, 0)], cell[bd(i, 1)], cell[bd(i, 2)]};
      int t = em.get_type(face, m_mesh->face(f));
      const int * idx = em.index[t];

      for(int j = 0; j < nfldof*2; j++)
        dof[bdof(i, j)] = fdof[idx[j]];
    }
    for(int i = 0; i < (int)idof.size(); i++)
      dof[idof[i]] = N + i;
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
    auto & fip = m_cellip.fip;
    auto & eip = m_cellip.eip;
    auto & vip = m_cellip.vip;
    std::cout<< "bou: " << bou << std::endl;
    std::cout<< "bdof: " << bdof << std::endl;
    std::cout<< "fip: " <<  fip << std::endl;
    std::cout<< "eip: " << eip << std::endl;
    std::cout<< "indof: " << std::endl;
    for(auto & bb : indof)
      std::cout<< bb << std::endl;

    std::cout<< "vdof: " << std::endl;
    for(auto & bb : vip)
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
    const auto & fip = m_cellip.fip;
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
    const auto & eip = m_cellip.eip;
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
    const auto & vip = m_cellip.vip;
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
    const auto & fip = m_cellip.fip;
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
    const auto & eip = m_cellip.eip;
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
  EdgeDofMap m_edgeMap;
  FaceDofMap m_faceMap;
  FaceDofClassify m_facedof;
  CellDofClassify m_celldof;
  InterplationPointInCell m_cellip;
};

}//end of Core
}//end of OF
#endif // end of BDMDof_h
