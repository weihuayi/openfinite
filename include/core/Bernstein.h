#include <math.h>
#include <algorithm>
#include <numeric>
#include <iostream>

#include "femcore.h"
#include "RArray2d.h"

namespace OF{
namespace FEM{

class Bernstein
{
public:
  Bernstein(RArray2d<int> & multiindex): m_multiindex(multiindex), m_tree(), 
    m_inverse_tree()
  {
    int ldof = multiindex.shape(0);
    int d = multiindex.shape(1)-1;

    m_tree.init(ldof, d);
    m_inverse_tree.init(ldof, d, -1);

    std::vector<int> idx(ldof, 0);
    for(int i = 0; i < ldof; i++)
    {
      if(multiindex(i, 0)==0)
        break;
      int t = i;
      int p = multiindex(0, 0);
      for(int j = 1; j < d; j++)
      {
        p -= multiindex(i, j-1);
        t += OF::Core::C(p+d-j, d-j); 
        m_tree(i, j-1) = t; 
        m_inverse_tree(t, idx[t]++) = i;
      }
      m_tree(i, d-1) = m_tree(i, d-2)+1;
      m_inverse_tree(m_tree(i, d-1), idx[m_tree(i, d-1)]++) = i;
    }
  }

  double de_casteljau(double * coef, double * lambda)
  {
    int ldof = m_multiindex.shape(0);
    int d = m_multiindex.shape(1)-1;
    int p = m_multiindex(0, 0);
    for(int i = 0; i < p+1; i++)
    {
      for(int j = 0; j < ldof; j++)
      {
        if(m_multiindex(j, 0)==i)
          break;
        coef[j] *= lambda[0];
        for(int k = 0; k < d; k++)
          coef[j] += lambda[k+1]*coef[m_tree(j, k)];
      }
    }
    return coef[0];
  };

private:
  RArray2d<int> m_multiindex;
  RArray2d<int> m_tree;
  RArray2d<int> m_inverse_tree;
};

double de_casteljau(double * coef, double * lambda, int p)
{
  int N = p+1;
  for(int i = 0; i < p+1; i++)
  {
    N--; 
    for(int j = 0; j < N; j++)
    {
      coef[j] *= lambda[0];
      coef[j] += lambda[1]*coef[j+1];
    }
  }
  return coef[0];
};

class Bernstein2d
{
public:
  Bernstein2d(RArray2d<int> & multiindex): m_multiindex(multiindex) {}

  double de_casteljau(double * coef, double * lambda)
  {
    int N = m_multiindex.shape(0);
    int p = m_multiindex(0, 0);
    for(int i = 0; i < p+1; i++)
    {
      N -= p+1-i; 
      for(int j = 0; j < N; j++)
      {
        /**!
         * (a, b, c): j, 
         * (a-1, b+1, c): m = j + p-a+1, 
         * (a-1, b, c+1): n = m+1
         * */
        coef[j] *= lambda[0];
        int m = j + p - m_multiindex(j, 0) + 1;
        coef[j] += lambda[1]*coef[m];
        coef[j] += lambda[2]*coef[m+1];
      }
    }
    return coef[0];
  };

private:
  RArray2d<int> m_multiindex;
};

class Bernstein3d
{
public:
  Bernstein3d(RArray2d<int> & multiindex): m_multiindex(multiindex) {}

  double de_casteljau(double * coef, double * lambda)
  {
    int N = m_multiindex.shape(0);
    int p = m_multiindex(0, 0);
    for(int i = 0; i < p+1; i++)
    {
      N -= (p+1-i)*(p+2-i)/2; 
      for(int j = 0; j < N; j++)
      {
        /**!
         * (a, b, c, d): j, 
         * (a-1, b+1, c, d): l = j+(p-a+1)(p-a+2)/2, 
         * (a-1, b, c+1, d): m = l + p-a+1-b
         * (a-1, b, c+1, d): n = m+1
         * */
        coef[j] *= lambda[0];
        int nq = p - m_multiindex(j, 0) + 1;
        int l = j + (nq+1)*nq/2;
        coef[j] += lambda[1]*coef[l]; 
        int m = l + nq - m_multiindex(j, 1);
        coef[j] += lambda[2]*coef[m]; 
        coef[j] += lambda[3]*coef[m+1];
      }
    }
    return coef[0];
  };

private:
  RArray2d<int> m_multiindex;
};

}
}
