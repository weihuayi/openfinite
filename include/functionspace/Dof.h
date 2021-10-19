/**
 * \file Dof.h
 * \author ccy
 * \date 2021/10/04
 * \brief Lagrange 有限元空间自由度管理
 */
#ifndef Dof_h
#define Dof_h

#include <map>
#include <array>
#include <vector>
#include <initializer_list>

#include "common/Multindex.h"

namespace OF{
namespace FunctionSpace{

/**
 * \brief 子单形上的自由度.
 */
class SubsimplxDof
{
public:
  SubsimplxDof() {}

  SubsimplxDof(int n, int p)
  {
    dof.resize(C(n+p, n));
  }

  int & at(int i)
  {
    return dof[i];
  }

public:
  std::vector<int> dof;
};

/**
 * \brief Lagrange 有限元空间自由度管理 
 * \param n 维数
 */
template<int n>
class Dof
{
public:
  Dof(){}

  Dof(int p):m_p(p)
  {
    /*
    std::vector<std::array<int, n+1>> index;
    muindex(p, index);
    int ord = index.size();
    m_dof = new SubsimplxDof[ord];
    for(int i = 0; i < ord; i++)
    {
      int dim = -1;
      for(int j = 0; j <= n; j++)
      {
        if(index[i][j]>0)
        {
          dim++;
        }
      }
    }
    */
  }
  void construct_subsimplex_to_local_dof()
  {
    std::vector<std::array<int, n+1>> index;
    std::map<int, int> idx2num;
    muindex(index);
    int ord = index.size();
    for(int i = 0; i < ord; i++)
    {
      int key = array_to_unique_number(index[i]);
      idx2num[key] = i;
    }
  }

private:
  int m_p;
  SubsimplxDof * m_dof;
};

} // end of namespace FunctionSpace 

} // end of namespace OF

#endif // end of Dof_h
