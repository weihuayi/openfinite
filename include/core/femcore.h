
/**
 * \file femcore.h
 * \brief 有限元计算中的核心函数
 * \author 陈春雨
 * \data 2021/11/13
 */
#ifndef FEMCore_h
#define FEMCore_h

#include <math.h>
#include <memory>
#include <vector>
#include <array>
#include <map>
#include <time.h>

#include <functional>

#include "RArray2d.h"

namespace OF{
namespace Core{

/**
 *\brief 计算组合数 A_n^m  = n*(n-1)*...*(n-m+1)
 */
int A(int n, int m)
{
  if(m>-1)
  {
    int a = 1;
    for(int i = 0; i < m; i++)
      a *= n-i;
    return a;
  }
  else
    return 0;
}

/**
 *\brief 计算组合数 C_n^m  = n*(n-1)*...*(n-m+1)/m!
 */
int C(int n, int m)
{
  if(m >= 0 && n >= m)
  {
    double a = 1;
    for(int i = 0; i < m; i++)
      a *= (double)(n-i)/(i+1);
    return (int)a;
  }
  else
    return 0;
}

void BMatrix(const int p, const RArray2d<double> & bc, RArray2d<double> & B)
{
  /** 计算 B */
  int l = bc.shape(1);
  B.init(p+1, l);
  std::fill(B[0], B[0]+l, 1.0);
  for(int i = 0; i < p; i++)/**< 赋值 */
  {
    for(int j = 0; j < l; j++)
      B(i+1, j) = p*bc(0, j) - i;
  }

  for(int i = 0; i < p; i++)/**< 连乘 */
  {
    for(int j = 0; j < l; j++)
      B(i+1, j) *= B(i, j);
  }

  int P = 1;
  for(int i = 0; i < p; i++)/**< 乘 P */
  {
    P *= i+1;
    for(int j = 0; j < l; j++)
      B(i+1, j)/= P;
  }
}

void FMatrix(int p, const RArray2d<double> & bc, RArray2d<double> & F)
{
  /** 计算 D^i */
  const int l = bc.shape(1);
  F.init(p+1, l);
  std::fill(F[0], F[0]+l, 0.0);
  RArray2d<double> D(p, p);
  for(int k = 0; k < l; k++)
  {
    for(int i = 0; i < p; i++)/**< 赋值 */
    {
      std::fill(D[i], D[i]+p, p*bc(0, k)-i);
      D(i, i) = p;
    }
    for(int i = 1; i < p; i++)/**< 连乘 */
    {
      for(int j = 0; j < p; j++)
        D(i, j) *= D(i-1, j);
    }
    for(int i = 0; i < p; i++)/**< 按行求和 */
    {
      for(int j = 0; j < i+1; j++)
        F(i+1, k) += D(i, j);
    }
  }

  int P = 1;
  for(int i = 0; i < p; i++)/**< 乘 P */
  {
    P *= i+1;
    for(int j = 0; j < l; j++)
      F(i+1, j)/= P;
  }
}

void LagrangeBasis(const RArray2d<double> & B, const RArray2d<int> & multiindex, 
    RArray2d<double> & val)
{
  int l = multiindex.shape(1);
  int NPhi = multiindex.shape(0);
  val.init(NPhi, 1);
  for(int i = 0; i < NPhi; i++)
  {
    val(i, 0) = B(multiindex(i, 0), 0);
    for(int j = 1; j < l; j++)
      val(i, 0) *= B(multiindex(i, j), j);
  }
}

/**
 * \brief 单元上基函数梯度在中心坐标 lambda 处的函数值, 其中 B, F 是 lambda
 * 有关的矩阵
 */
void LagrangeGradBasis(const RArray2d<double> & B, const RArray2d<double> & F, 
    const RArray2d<int> & multiindex, const RArray2d<double> glambda, 
    RArray2d<double> & val)
{
  int GD = glambda.shape(1);
  int l = multiindex.shape(1);
  int NPhi = multiindex.shape(0); /**< 基函数个数 */
  val.init(NPhi, GD);
  for(int i = 0; i < NPhi; i++) /**< 第 i 个基函数的梯度 */
  {
    double phi = B(multiindex(i, 0), 0);
    for(int j = 1; j < l; j++)
      phi *= B(multiindex(i, j), j); /**< 第 i 个基函数的值 */

    for(int j = 0; j < l; j++)
    {
      double m = F(multiindex(i, j), j)*phi/B(multiindex(i, j), j); /**< \nabla \lambda 的系数 */
      for(int k = 0; k < GD; k++)
        val(i, k) += m*glambda(j, k);
    }
  }
}

/**
 * \brief 单元上的质量矩阵
 */
template<typename Barycentric>
void mass_matrix(std::vector<Barycentric> & bcs, const std::vector<double> & ws, 
    RArray2d<int> & multiindex, double J, RArray2d<double> & M)
{
  int NQ = bcs.size();
  int ldof = multiindex.shape(0);
  M.init(ldof, ldof);
  for(int q = 0; q < NQ; q++)
  {
    RArray2d<double> val;
    LagrangeBasis(bcs[q].B, multiindex, val);
    for(int i = 0; i < ldof; i++)
    {
      for(int j = 0; j < ldof; j++)
        M(i, j) += J*ws[q]*val(i, 0)*val(j, 0);
    }
  }
}

/**
 * \brief 单元上的刚度矩阵
 */
template<typename Barycentric>
void stiff_matrix(const std::vector<Barycentric> & bcs, const std::vector<double> & ws, 
    const double & J, const RArray2d<int> & multiindex, 
    const RArray2d<double> & glambda, RArray2d<double> & A)
{
  int NQ = bcs.size();
  int ldof = multiindex.shape(0);
  A.init(ldof, ldof);
  for(int q = 0; q < NQ; q++)
  {
    RArray2d<double> val;
    LagrangeGradBasis(bcs[q].B, bcs[q].F, multiindex, glambda, val);
    for(int i = 0; i < ldof; i++)
    {
      for(int j = 0; j < ldof; j++)
      {
        double s = 0.0;
        for(int k = 0; k < (int)val.shape(1); k++)
          s += val(i, k)*val(j, k);
        A(i, j) += J*ws[q]*s;
      }
    }
  }
}

/**
 * \brief 求下一个多重指标
 * \param nex_index 若下一个多重指标存在, 将其存储在 nex_index 中
 * \return re 若存在下一个指标则返回 true, 否则返回 false
 */
bool nex_multi_index(int l, int * index, int * nexIndex)
{
  bool re = false;
  int non_zero_idx = l-1;
  for(int i = l-2; i > -1; i--)
  {
    if(index[i]>0 && !re)
    {
      non_zero_idx = i;
      re = true;
    }
    nexIndex[i] = index[i];
  }
  nexIndex[l-1] = 0;
  if(re)
  {
    nexIndex[non_zero_idx]--;
    nexIndex[non_zero_idx+1] = index[l-1]+1;
  }
  return re;
}

/**
 * \brief 生成所有阶数为 p 的多重指标
 * \param p 多重指标的阶数
 * \param midxs 将所有多重指标存在这里
 */
void get_all_multiindex(int p, int l, RArray2d<int> & midxs)
{
  int s = C(l-1+p, l-1);
  midxs.init(s, l);
  midxs(0, 0) = p;
  for(int i = 0; i < s-1; i++)
    nex_multi_index(l, midxs[i], midxs[i+1]);
}

/**
 * \brief 计算自身是第几个排列 如: 
 * [2, 0, 0]->0
 * [1, 1, 0]->1
 * [1, 0, 1]->2
 * [0, 2, 0]->3
 * [0, 1, 1]->4
 * [0, 0, 2]->5
 */
int multiindex_to_number(int l, int p, int * multiindex)
{
  int num = 0, ord = p;
  for(int i = 0; i < l-1; i++)
  {
      ord -= multiindex[i];
      num += C(l-i-2+ord, l-i-1);
  }
  return num;
}

/**
 * \brief 计算 m_data[index] 是第几个排列 如: 
 */
int multiindex_to_number(int l, int p, int * multiindex, int * index)
{
  int num = 0, ord = p;
  for(int i = 0; i < l-1; i++)
  {
    ord -= multiindex[index[i]];
    num += C(l-i-2+ord, l-i-1);
  }
  return num;
}

int permutation_to_number(int * per, int l)
{
  int num = 0, idx[l];
  for(int i = 0; i < l; i++)
    idx[i] = i;
  for(int i = 0; i < l; i++)
  {
    num += idx[per[i]]*A(l-1-i, l-1-i);
    int k = 0;
    for(int j = 0; j < l; j++)
    {
      if(j != per[i] && idx[j]>=0)
      {
        idx[j] = k++;
      }
      else if(j == per[i])
      {
        idx[j] = -1;
      }
    }
  }
  return num;
}

void cross(double * a, double * b, double * c)
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

template<typename T>
T dot(T * a, T * b, int l = 3)
{
  T s = 0;
  for(int i = 0; i < l; i++)
    s += a[i]*b[i];
  return s;
}

/**
 * \brief 判断两个数组大小 a > b 返回 true, 否则返回 false.
 */
template<typename T>
bool greater(T * a, T * b, int l)
{
  bool f = true;
  for(int i = 0; i < l; i++)
  {
    if(a[i]<b[i])
    {
      f = false;
      break;
    }
    else if(a[i]==b[i])
    {
      continue;
    }
    else
    {
      break;
    }
  }
  return f;
}

bool equal(int * a, int * b, int l)
{
  bool f = true;
  for(int i = 0; i < l; i++)
  {
    if(a[i]!=b[i])
    {
      f = false;
      break;
    }
  }
  return f;
}

double cell_volume(const RArray2d<double> & node, const RArray2d<int> & cell, int c)
{
  RArray2d<double> v(3, 3);
  double cro[3];
  for(int i = 0; i < 3; i++)
  {
    for(int j = 0; j < 3; j++)
      v(i, j) = node(cell(c, i+1), j) - node(cell(c, 0), j);
  }
  cross(v[0], v[1], cro);
  return dot(cro, v[2])/6.0;
}

double grad_lambda(const RArray2d<double> & node, const RArray2d<int> & cell, 
    int c, RArray2d<int> localFace, RArray2d<double> & Dlambda)
{
  RArray2d<double> v(2, 3); /**< 每个面上的两个向量*/
  Dlambda.init(4, 3);
  double vol = cell_volume(node, cell, c);
  for(int i = 0; i < 4; i++)
  {
    int j = localFace[i][0], k = localFace[i][1], m = localFace[i][2];
    for(int l = 0; l < 3; l++)
    {
      v(0, l) = node(cell(c, k), l) - node(cell(c, j), l);
      v(1, l) = node(cell(c, m), l) - node(cell(c, j), l);
    }
    cross(v[1], v[0], Dlambda[i]);
    for(int l = 0; l < 3; l++)
      Dlambda(i, l) /= 6*vol;
  }
  return vol;
}

template<typename T>
void sort_RArray2d(RArray2d<T> & Arr, std::vector<int> & idx)
{
  int N = Arr.shape(0);
  idx.resize(N);
  for(int i = 0; i < N; i++)
    idx[i] = i;

  auto f0 = [&Arr](int a, int b) {return !OF::Core::greater(Arr[a], Arr[b], 3);};
  std::sort(idx.begin(), idx.end(), f0);
}

}//end of Core
}//end of OF
#endif // end of FEMCore_h
