
/**
 * \file Permutation.h
 * \author ccy
 * \date 2021/10/04
 * \brief Lagrange 有限元空间自由度管理
 */
#ifndef Permutation_h
#define Permutation_h

#include <algorithm>
#include <initializer_list>
#include <map>
#include <array>
#include <vector>

namespace OF{
namespace Common{

/**
 *\brief 计算组合数 C_n^m  = n*(n-1)*...*(n-m+1)/m!
 */
int C(int n, int m)
{
  double a = 1;
  for(int i = 0; i < m; i++)
    a *= (double)(n-i)/(i+1);
  return (int)a;
}

/**
 * \brief 多重指标类
 * \param l 多重指标的长度
 */
template<int l>
class Permutation
{
public:

  /** 
   * \brief 构造函数, 构造第一个排列
   */
  Permutation():data{0} 
  {
    for(int i = 1; i < l; i++)
      data[i] = i;
  }

  /** 
   * \brief 构造函数
   * \param list 初始化列表数据类型, 用来初始化多重指标
   */
  Permutation(std::initializer_list<int> list):data{0}
  {
    int i = 0;
    for(int val : list)
      data[i++] = val;
  }

  /**
   * \brief 生成 0 - p-1 中 l 个数的全排列
   * \param p 元素个数
   * \param pers 将所有排列存在这里
   */
  static void gen_all_index(const int & p, std::vector<Permutation<l> > & pers)
  {
    int size = C(l-1+p, l-1);
    pers.resize(size);
    pers[0][0] = p;
    for(int i = 0; i < size-1; i++)
      pers[i].nex_permutation(pers[i+1]);
  }

  /**
   * \brief [] 操作符重载
   * \param i 索引
   * \return 多重指标的第 i 个整数
   */
  int & operator[](const int & i)
  {
      return data[i];
  }

  /**
   * \brief 求下一个多重指标
   * \param nex_index 若下一个多重指标存在, 将其存储在 nex_index 中
   * \return re 若存在下一个指标则返回 true, 否则返回 false
   */
  bool nex_permutation(Permutation<l> & nex_index)
  {
    auto nex_data = nex_index.get_data();
    memcpy(nex_data, &data, &data+l);
    bool re = std::next_permutation(nex_data, nex_data+l);
    return re;
  }

  /** \brief 计算阶数 */
  int order()
  {
    int ord = 0;
    for(auto val : data)
      ord += val;
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
    int num = 0;
    int ord = order();
    for(int i = 0; i < l-1; i++)
    {
        ord -= data[i];
        num += C(l-i-2+ord, l-i-1);
    }
    return num;
  }

  /**
   * \brief 计算 data[index] 是第几个排列 如: 
   */
  int to_number(int * index)
  {
    int num = 0;
    int ord = order();
    for(int i = 0; i < l-1; i++)
    {
        ord -= data[index[i]];
        num += C(l-i-2+ord, l-i-1);
    }
    return num;
  }

  int * get_data()
  {
    return &data;
  }

private:
  int data[l]; /**< 多重指标的数据 */
};

} // end of namespace Common 

} // end of namespace OF

#endif // end of Permutation_h
