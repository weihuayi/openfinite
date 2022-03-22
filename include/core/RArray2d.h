#ifndef RArray2d_h
#define RArray2d_h

#include <vector>
#include <array>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <utility>
#include <cassert>


template<class T>
class RArray2d 
{
public:
  // 基本类型 
  using Data  = std::vector<T>;

  using value_type = typename Data::value_type;
  using size_type  = typename Data::size_type;

  // 引用 
  using reference       = typename Data::reference;
  using const_reference = typename Data::const_reference;

  // 迭代子 
  using iterator       = typename Data::iterator;
  using const_iterator = typename Data::const_iterator;

  // 反向迭代子 
  using reverse_iterator       = typename Data::reverse_iterator;
  using const_reverse_iterator = typename Data::const_reverse_iterator;

  using Shape = typename std::array<size_type, 2>;

  // 空构造函数 
  RArray2d() = default;

  // 默认插入 m*n 个值 
  RArray2d(size_type m, size_type n)
    : m_shape{m, n}, m_data(m*n)
  {}

  // 复制初始化二维数组 m*n
  RArray2d(size_type m, size_type n, const_reference val)
    : m_shape{m, n}, m_data(m*n, val)
  {}

  // 复制构造函数
  RArray2d(const RArray2d<T> & val)
    : m_shape{val.shape(0), val.shape(1)}, m_data(val.shape(0)*val.shape(1))
  {
    for(int i = 0; i < m_shape[0]; i++)
    {
      for(int j = 0; j < m_shape[1]; j++)
        at(i, j) = val(i, j);
    }
  }

  void init(const size_type m, const size_type n)
  {
    m_shape[0] = m;
    m_shape[1] = n;
    m_data.resize(m*n, 0.0);
  }

  // 1 维迭代子 
  iterator begin() { return m_data.begin(); }
  iterator end() { return m_data.end(); }

  const_iterator begin() const { return m_data.begin(); }
  const_iterator end() const { return m_data.end(); }

  const_iterator cbegin() const { return m_data.cbegin(); }
  const_iterator cend() const { return m_data.cend(); }

  reverse_iterator rbegin() { return m_data.rbegin(); }
  reverse_iterator rend() { return m_data.rend(); }

  const_reverse_iterator rbegin() const { return m_data.rbegin(); }
  const_reverse_iterator rend() const { return m_data.rend(); }

  const_reverse_iterator crbegin() const { return m_data.crbegin(); }
  const_reverse_iterator crend() const { return m_data.crend(); }

  // 元素访问 (行优先)
  reference operator() (size_type const i, size_type const j)
  {
    assert(i >=0 && i < m_shape[0] && j >= 0 && j < m_shape[1]);
    return m_data[m_shape[1]*i + j];
  }

  const_reference operator() (size_type const i, size_type const j) const
  {
    assert(i >=0 && i < m_shape[0] && j >= 0 && j < m_shape[1]);
    return m_data[m_shape[1]*i + j];
  }

  // 返回第 (i, j) 个位置上的元素引用
  reference at(size_type const i, size_type const j)
  {
    assert(i >=0 && i < m_shape[0] && j >= 0 && j < m_shape[1]);
    return m_data.at(m_shape[1]*i + j);
  }

  // 返回第 (i, j) 个位置上的元素常引用
  const_reference at(size_type const i, size_type const j) const
  {
    assert(i >=0 && i < m_shape[0] && j >= 0 && j < m_shape[1]);
    return m_data.at(m_shape[1]*i + j);
  }

  // 返回第 i 行的首指针
  T * operator[](size_type const  i) 
  {
      return &m_data[m_shape[1]*i];
  }

  // 返回第 i 行的首指针 const 版本
  const T * operator[](size_type const i) const
  {
      return &m_data[m_shape[1]*i];
  }

  // 增加一行
  void push_back(const std::initializer_list<T> & l)
  {
    assert(l.size() == m_shape[1]);
    for(auto val:l)
      m_data.push_back(val);
    m_shape[0] += 1;
  }

//  void resize(size_type m, size_type n)
//  {
//    // 新的二维数组 m*n 
//    RArray2d tmp(m, n);
//
//    auto mc = std::min(m_shape[0], m);
//    auto mr = std::min(m_shape[1], n);
//
//    for (size_type i(0U); i < mr; ++i)
//    {
//      auto row = begin() + i*m_shape[1];
//      auto tmp_row = tmp.begin() + i*nc;
//      std::move(row, row + mc, tmp_row);
//    }
//    *this = std::move(tmp);
//  }

  void clear()
  {
    m_shape = {0, 0};
    m_data.clear();
  }

  size_type data_size() const { return m_data.size(); }

  // 预留内存 
  void reserve(const size_type m, const size_type n)
  {
    m_shape[1] = n;
    m_data.reserve(m*n);
  }

  bool empty() const { return m_data.empty(); }

  size_type shape(unsigned i) const 
  {
    assert(i >=0 && i < 2);
    return m_shape[i];
  }

  Shape & shape()
  {
    return m_shape;
  }

  /*!
   * \brief 返回原始数组的指针
   */
  value_type * row_data() {return m_data.data();}

  /*!
   * \brief 返回内部数据的引用，这里是 std::vector 
   */
  Data & data() {return m_data;}

private:

  Shape m_shape{0u, 0u};
  Data m_data{};
};

template<class T>
void swap(RArray2d<T> & lhs, RArray2d<T> & rhs)
{
  lhs.swap(rhs);
}

template<class T>
bool operator == (RArray2d<T> const &a, RArray2d<T> const &b)
{
  if (a.shape(0) != b.shape(0) || a.shape(1) != b.shape(1))
  {
    return false;
  }
  return std::equal(a.begin(), a.end(), b.begin(), b.end());
}

template<class T>
bool operator != (RArray2d<T> const &a, RArray2d<T> const &b)
{
  return !(a == b);
}

template<class T>
std::ostream& operator << (std::ostream & os, const RArray2d<T> & a)
{
  std::cout << "RArray2d("<<  a.shape(0) << "," << a.shape(1) << "):" << std::endl;
  for(typename RArray2d<T>::size_type i = 0U; i < a.shape(0); i++)
  {
      os << i << " : "; 
      for(typename RArray2d<T>::size_type j = 0U; j < a.shape(1); j++)
      {
          os << a(i, j) << " "; 
      }
      os << "\n";
  }
  os << "\n";
  return os;
}

#endif // end of RArray2d_h 
