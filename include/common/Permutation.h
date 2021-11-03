
/**
 * \file Permutation.h
 * \author ccy
 * \date 2021/10/04
 * \brief 排列组合类
 */
#ifndef Permutation_h
#define Permutation_h

#include <algorithm>
#include <initializer_list>
#include <map>
#include <array>
#include <vector>

#include "Function.h"

namespace OF{
namespace Common{

/** \brief 排列组合类 */
class Permutation
{
public:

  /** 构造函数 */
  Permutation(): m_data(NULL), m_len(0){}

  /** 
   * \brief 构造函数, 构造第一个阶数为 p 的排列组合
   * \param l 排列组合长度
   * \param p 排列组合的最大值
   */
  Permutation(int l): m_data(NULL), m_len(l)
  {
    m_data = new int[l]();
    for(int i = 0; i < l; i++)
      m_data[i] = i;
  }

  /** 
   * \brief 构造函数
   * \param list 初始化列表数据类型, 用来初始化排列组合
   */
  Permutation(std::initializer_list<int> list): m_data(NULL), m_len(list.size())
  {
    m_data = new int[m_len]();
    int i = 0;
    for(int val : list)
      m_data[i++] = val;
  }

  /** \brief 复制构造函数 */
  Permutation(const Permutation & mu):m_len(mu.size())
  {
    m_data = new int[m_len]();
    for(int i = 0; i < m_len; i++)
      m_data[i] = mu[i];
  }

  /** \brief 析构函数 */
  ~Permutation()
  {
    if(m_data != NULL)
    {
      delete [] m_data;
    }
  }

  /**
   * \brief 生成 0 - p-1 中 l 个数的全排列
   * \param p 元素个数
   * \param pers 将所有排列存在这里
   */
  static void gen_all_permutation(const int & l, std::vector<Permutation> & pers)
  {
    int s = A(l, l);
    pers.resize(s, l);
    for(int i = 0; i < s-1; i++)
      pers[i].nex_permutation(pers[i+1]);
  }

  /**
   * \brief 获取 A1 在 A0 中的编号
   */
  template<typename Arr0, typename Arr1>
  static Permutation chack(const Arr0 & A0, const Arr1 & A1, int n)
  {
    Permutation per(n);
    for(int i = 0; i < n; i++)
    {
      for(int j = 0; j < n; j++)
      {
        if(A0[i]==A1[j])
          per[i] = j;
      }
    }
    return per;
  }

  /**
   * \brief [] 操作符重载
   * \param i 索引
   * \return 排列组合的第 i 个整数
   */
  int & operator[](const int & i)const
  {
    return m_data[i];
  }

  /**
   * \brief 求下一个排列组合
   * \param nex_index 若下一个排列组合存在, 将其存储在 nex_index 中
   * \return re 若存在下一个指标则返回 true, 否则返回 false
   */
  bool nex_permutation(Permutation & nex_per)
  {
    int * nex_data = nex_per.get_data();
    std::copy(m_data, m_data+m_len, nex_data);
    return std::next_permutation(nex_data, nex_data+m_len);
  }

  bool nex_permutation()
  {
    return std::next_permutation(m_data, m_data+m_len);
  }

  /**
   * \brief 计算自身是第几个排列 如: 
   * [0, 1, 2]->0
   * [0, 2, 1]->1
   * [1, 0, 2]->2
   * [1, 2, 0]->3
   * [2, 0, 1]->4
   * [2, 1, 0]->5
   */
  int to_number()
  {
    int num = 0, l = m_len, idx[m_len];
    for(int i = 0; i < l; i++)
      idx[i] = i;
    for(int i = 0; i < l; i++)
    {
      num += idx[m_data[i]]*A(l-1-i, l-1-i);
      int k = 0;
      for(int j = 0; j < l; j++)
      {
        if(j != m_data[i] && idx[j]>=0)
        {
          idx[j] = k++;
        }
        else if(j == m_data[i])
        {
          idx[j] = -1;
        }
      }
    }
    return num;
  }

  int * get_data()
  {
    return m_data;
  }

  const int & size()const
  {
    return m_len;
  }

private:
  int * m_data; /**< 排列的数据 */
  int m_len; /**< 排列的长度 */ 
};

} // end of namespace Common 

} // end of namespace OF

#endif // end of Permutation_h
