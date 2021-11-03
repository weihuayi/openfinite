/**
 * \file IncreaseSequence.h
 * \author ccy
 * \date 2021/10/04
 * \brief 递增序列
 */
#ifndef IncreaseSequence_h
#define IncreaseSequence_h

#include <initializer_list>
#include <map>
#include <array>
#include <vector>

#include "Function.h"

namespace OF{
namespace Common{

/**
 * \brief 递增序列类
 */
class IncreaseSequence
{
public:

  /** 构造函数 */
  IncreaseSequence(): m_data(NULL), m_len(0){}

  /** 
   * \brief 构造函数, 构造一个阶数为 p 的递增序列
   * \param l 递增序列长度
   */
  IncreaseSequence(int l): m_data(NULL), m_len(l)
  {
    m_data = new int [l]();
    for(int i = 0; i < l; i++)
      m_data[i] = i;
  }

  /** 
   * \brief 构造函数
   * \param list 初始化列表数据类型, 用来初始化递增序列
   */
  IncreaseSequence(std::initializer_list<int> list):m_data(NULL), m_len(list.size())
  {
    m_data = new int [m_len]{0};
    int i = 0;
    for(int val : list)
      m_data[i++] = val;
  }

  /** \brief 复制构造函数 */
  IncreaseSequence(const IncreaseSequence & mu):m_len(mu.size())
  {
    m_data = new int [m_len]();
    for(int i = 0; i < m_len; i++)
      m_data[i] = mu[i];
  }

  /** \brief 析构函数 */
  ~IncreaseSequence()
  {
    if(m_data != NULL)
    {
      delete [] m_data;
    }
  }

  /**
   * \brief 生成所有阶数为 p 的递增序列
   * \param p 递增序列的阶数
   * \param seqs 将所有递增序列存在这里
   */
  static void gen_all_sequence(const int & n, int l, 
      std::vector<IncreaseSequence> & seqs)
  {
    int s = C(n, l);
    seqs.resize(s, l);
    for(int i = 0; i < s-1; i++)
      seqs[i].nex_increaseing_sequence(n, seqs[i+1]);
  }

  /**
   * \brief [] 操作符重载
   * \param i 索引
   * \return 递增序列的第 i 个整数
   */
  int & operator[](const int & i)const
  {
    return m_data[i];
  }

  /**
   * \brief 求下一个递增序列
   * \param nex_index 若下一个递增序列存在, 将其存储在 nex_index 中
   * \return re 若存在下一个序列则返回 true, 否则返回 false
   */
  bool nex_increaseing_sequence(const int & n, IncreaseSequence & nex_seq)
  {
    int l = m_len, idx = -1;
    for(int i = 0; i < l; i++)
    {
      if(m_data[l-1-i]<n-1-i)
      {
        idx = l-1-i;
        break;
      }
    }
    if(idx>-1)
    {
      for(int j = 0; j < idx; j++)
        nex_seq[j] = m_data[j];
      for(int j = idx; j < l; j++)
        nex_seq[j] = m_data[idx]+j-idx+1; 
    }
    return idx>0;
  }

  /**
   * \brief 计算自身是第几个序列 如: 
   * [0, 1, 2]->0
   * [0, 1, 3]->1
   * [0, 2, 3]->2
   * [1, 2, 3]->3
   */
  int to_number(int n)
  {
    int num = 0, l = m_len, k = 0;
    if(l>1)
    {
      for(int i = 0; i < l-1; i++)
      {
        for(int j = k; j < m_data[i]; j++)
        {
          num += C(n-j-1, l-i-1);
        }
        k = m_data[i]+1;
      }
      num += m_data[l-1]-m_data[l-2]-1;
      return num;
    }
    else
    {
      return m_data[0];
    }
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
  int * m_data; /**< 递增序列的数据 */
  int m_len; /**< 递增序列长度 */
};

} // end of namespace Common 

} // end of namespace OF

#endif // end of IncreaseSequence_h
