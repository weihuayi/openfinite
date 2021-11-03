
/**
 * \file Multindex.h
 * \author ccy
 * \date 2021/10/04
 * \brief 多重指标
 */
#ifndef Multindex_h
#define Multindex_h

#include <initializer_list>
#include <map>
#include <array>
#include <vector>

#include "Function.h"

namespace OF{
namespace Common{

/**
 * \brief 多重指标类
 */
class MultiIndex
{
public:

  /** 构造函数 */
  MultiIndex(): m_data(NULL), m_len(0) {}

  /** 
   * \brief 构造函数, 构造一个阶数为 p 的多重指标
   * \param l 多重指标长度
   */
  MultiIndex(int l): m_data(NULL), m_len(l)
  {
    m_data = new int [l]();
  }

  /** 
   * \brief 构造函数, 构造第一个阶数为 p 的多重指标
   * \param p 初始化列表数据类型, 用来初始化多重指标
   * \param l 多重指标长度
   */
  MultiIndex(int p, int l): m_data(NULL), m_len(l)
  {
    m_data = new int [l]();
    m_data[0] = p;
  }

  /** 
   * \brief 构造函数
   * \param list 初始化列表数据类型, 用来初始化多重指标
   */
  MultiIndex(std::initializer_list<int> list):m_data(NULL), m_len(list.size())
  {
    m_data = new int [m_len]();
    int i = 0;
    for(int val : list)
      m_data[i++] = val;
  }

  /** \brief 复制构造函数 */
  MultiIndex(const MultiIndex & mu):m_len(mu.size())
  {
    m_data = new int [m_len]();
    for(int i = 0; i < m_len; i++)
      m_data[i] = mu[i];
  }

  /** \brief 析构函数 */
  ~MultiIndex()
  {
    if(m_data != NULL)
    {
      delete [] m_data;
    }
  }

  /**
   * \brief 生成所有阶数为 p 的多重指标
   * \param p 多重指标的阶数
   * \param midxs 将所有多重指标存在这里
   */
  static void gen_all_index(const int & p, int l, std::vector<MultiIndex> & midxs)
  {
    int s = C(l-1+p, l-1);
    midxs.resize(s, l);
    midxs[0][0] = p;
    for(int i = 0; i < s-1; i++)
      midxs[i].nex_multi_index(midxs[i+1]);
  }

  /**
   * \brief [] 操作符重载
   * \param i 索引
   * \return 多重指标的第 i 个整数
   */
  int & operator[](const int & i)const
  {
    return m_data[i];
  }

  /**
   * \brief 求下一个多重指标
   * \param nex_index 若下一个多重指标存在, 将其存储在 nex_index 中
   * \return re 若存在下一个指标则返回 true, 否则返回 false
   */
  bool nex_multi_index(MultiIndex & nex_index)
  {
    bool re = false;
    int non_zero_idx = m_len-1, l = m_len;
    for(int i = l-2; i > -1; i--)
    {
      if(m_data[i]>0 && !re)
      {
        non_zero_idx = i;
        re = true;
      }
      nex_index[i] = m_data[i];
    }
    nex_index[l-1] = 0;
    if(re)
    {
      nex_index[non_zero_idx]--;
      nex_index[non_zero_idx+1] = m_data[l-1]+1;
    }
    return re;
  }

  bool nex_multi_index()
  {
    MultiIndex copy_idx(*this);
    return copy_idx.nex_multi_index(*this);
  }

  /** \brief 计算阶数 */
  int order()
  {
    int ord = 0;
    for(int i = 0; i < m_len; i++)
      ord += m_data[i];
    return ord;
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
  int to_number()
  {
    int num = 0, l = m_len;
    int ord = order();
    for(int i = 0; i < l-1; i++)
    {
        ord -= m_data[i];
        num += C(l-i-2+ord, l-i-1);
    }
    return num;
  }

  /**
   * \brief 计算 m_data[index] 是第几个排列 如: 
   */
  template<typename Arr>
  int to_number(const Arr & index)
  {
    int num = 0, l = m_len;
    int ord = order();
    for(int i = 0; i < l-1; i++)
    {
      ord -= m_data[index[i]];
      num += C(l-i-2+ord, l-i-1);
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
  int * m_data; /**< 多重指标的数据 */
  int m_len; /**< 多重指标长度 */
};

} // end of namespace Common 

} // end of namespace OF

#endif // end of Multindex_h
