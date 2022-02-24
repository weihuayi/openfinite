
/**
 * \file dof.h
 * \brief 有限元自由度管理
 * \author 陈春雨
 * \data 2022/2/19
 */
#ifndef Dof_h
#define Dof_h

#include <math.h>
#include <memory>
#include <vector>
#include <array>
#include <map>
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

  int get_type(int * E0, int * E1)
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
  RArray2d<int> bdof;
  std::vector<int> indof;

  DofClassify(int p): boundary(), bdof(), indof(0)
  {
    RArray2d<int> multiindex;
    get_all_multiindex(p, N+1, multiindex);

    int l0 = C(N-1+p, N-1), l1 = C(N+p-2, N), ldof = multiindex.shape(0);
    boundary.init(N+1, N);
    bdof.init(N+1, l0);
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
    int l[N+1] = {0};
    for(int i = 0; i < ldof; i++)
    {
      bool flag = true;
      for(int j = 0; j < N+1; j++)
      {
        if(multiindex(i, j)==0)
        {
          flag = false;
          bdof(j, l[j]++) = i;
        }
      }
      if(flag)
        indof.push_back(i);
    }
  }
};

template<typename Mesh>
class LFEMDof3d
{
public:
  LFEMDof3d(int p, std::shared_ptr<Mesh> mesh): m_p(p), 
    m_mesh(mesh), m_faceMap(p), m_edgeMap(p), m_celldof(p), m_facedof(p) 
  {
    get_all_multiindex(p, 4, multiindex3);
    get_all_multiindex(p, 3, multiindex2);
  }
  
  void edge_to_dof(int e, int * dof)
  {
    int NN = m_mesh->number_of_nodes();
    int edof = number_of_local_dofs(1);
    dof[0] = m_mesh->edge(e, 0);

    for(int i = 1; i < edof-1; i++)
      dof[i] = NN + e*(m_p-1) + i;

    dof[edof-1] = m_mesh->edge(e, 1);
  }

  void face_to_dof(int f, int * dof)
  {
    int NN = m_mesh->number_of_nodes();
    int NE = m_mesh->number_of_edges();
    int eidof = m_p-1, fidof = (m_p-1)*(m_p-2)/2;

    const int * face = m_mesh->face(f);

    const auto & em = m_edgeMap;
    const auto & idof = m_facedof.indof;
    const auto & bdof = m_facedof.bdof;
    const auto & bd = m_facedof.boundary; 

    int edof[eidof];
    for(int i = 0; i < 3; i++)
    {
      int e = m_mesh->face_to_edge(f, i);
      edge_to_dof(e, edof);

      int edge[2] = {face[bd(i, 0)], face[bd(i, 1)]};
      int t = em.type(edge, m_mesh->edge(e));
      int * idx = em.index[t];

      for(int j = 0; j < eidof; j++)
        dof[bdof(i, j)] = edof[idx[j]];
    }
    for(int i = 0; i < (int)idof.size(); i++)
      dof[idof[i]] = NN + NE*eidof + f*fidof + i;
  }

  void cell_to_dof(int c, int * dof)
  {
    int NN = m_mesh->number_of_nodes();
    int NE = m_mesh->number_of_edges();
    int NF = m_mesh->number_of_faces();
    int eidof = m_p-1, fidof = eidof*(m_p-2)/2, cidof = fidof*(m_p-3)/3;

    const int * cell = m_mesh->cell(c);

    const auto & em = m_faceMap;
    const auto & idof = m_celldof.indof;
    const auto & bdof = m_celldof.bdof;
    const auto & bd = m_celldof.boundary; 

    int fdof[fidof];
    for(int i = 0; i < 4; i++)
    {
      int f = m_mesh->cell_to_face(c, i);
      face_to_dof(f, fdof);

      int face[3] = {cell[bd(i, 0)], cell[bd(i, 1)], cell[bd(i, 2)]};
      int t = em.type(face, m_mesh->face(f));
      int * idx = em.index[t];

      for(int j = 0; j < fidof; j++)
        dof[bdof(i, j)] = fdof[idx[j]];
    }
    for(int i = 0; i < (int)idof.size(); i++)
      dof[idof[i]] = NN + NE*eidof + NF*fidof + c*cidof + i;
  }

  int number_of_local_dofs(int d = 3)
  {
    return C(d+m_p, d);
  }

private:
  int m_p;
  std::shared_ptr<Mesh> m_mesh;
  RArray2d<int> multiindex2; 
  RArray2d<int> multiindex3; 
  DofMap<2> m_faceMap;
  DofMap<1> m_edgeMap;
  DofClassify<3> m_celldof;
  DofClassify<2> m_facedof;
};

}//end of Core
}//end of OF
#endif // end of Dof_h
