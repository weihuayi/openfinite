
/**
 * \file femcore.h
 * \brief 中心坐标类
 * \author 陈春雨
 * \date 2021/11/13
 */
#ifndef Barycentric_h
#define Barycentric_h

#include <array>
#include "RArray2d.h"

namespace OF{
namespace Core{

template<int N>
class Barycentric
{
public:
  std::array<double, N+1> data;
  RArray2d<double> B;
  RArray2d<double> F;

public:
  Barycentric(): data(), B(), F(){}

  template<typename Container>
  Barycentric(Container list): data()
  {
    for(int i = 0; i < N+1; i++)
      data[i] = list[i];
  }

  int dim()
  {
    return N;
  }

  int size()
  {
    return N+1;
  }

  double & operator[](int  i) 
  {
    return data[i];
  }

  const double & operator[](int i) const
  {
    return data[i];
  }

  void compute_B(const int p)
  {
    /** 计算 B */
    int l = N+1;
    B.init(p+1, l);
    std::fill(B[0], B[0]+l, 1.0);
    for(int i = 0; i < p; i++)/**< 赋值 */
    {
      for(int j = 0; j < l; j++)
        B(i+1, j) = p*data[j] - i;
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

  void compute_F(int p)
  {
    /** 计算 D^i */
    int l = N+1;
    F.init(p+1, l);
    std::fill(F[0], F[0]+l, 0.0);
    RArray2d<double> D(p, p);
    for(int k = 0; k < l; k++)
    {
      for(int i = 0; i < p; i++)/**< 赋值 */
      {
        std::fill(D[i], D[i]+p, p*data[k]-i);
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
};

}//end of Core
}//end of OF
#endif // end of Barycentric_h
