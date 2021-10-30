/**
 * \file Function.h
 * \author ccy
 * \date 2021/10/04
 * \brief 常用函数
 */
#ifndef Function_h
#define Function_h

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
  if(m>-1)
  {
    double a = 1;
    for(int i = 0; i < m; i++)
      a *= (double)(n-i)/(i+1);
    return (int)a;
  }
  else
    return 0;
}


#endif // end of Function_h
