#ifndef RArray3d_h
#define RArray3d_h

#include <vector>
#include <array>
#include <algorithm>
#include <iterator>
#include <utility>
#include <cassert>


namespace TOPT {
namespace AlgebraObject {

template<class T>
class RArray3d 
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

  using Shape = typename std::array<size_type, 3>;

  // 空构造函数 
  RArray3d() = default;

  // 默认插入 m*n*p 个值 
  RArray3d(size_type m, size_type n, size_type p)
    : m_shape{m, n, p}, m_data(m*n*p)
  {}

  // 复制初始化三维数组 m*n*p
  RArray3d(size_type m, size_type n, size_type p, const_reference val)
    : m_shape{m, n, p}, m_data(m*n*p, val)
  {}

  void init(const size_type m, const size_type n, const size_type p)
  {
    m_shape[0] = m;
    m_shape[1] = n;
    m_shape[2] = p;
    m_data.resize(m*n*p);
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
  reference operator() (size_type const i, size_type const j, size_type const k)
  {
    return m_data[m_shape[2]*(m_shape[1]*i + j) + k];
  }

  const_reference operator() (size_type const i, size_type const j, size_type const k) const
  {
    return m_data[m_shape[2]*(m_shape[1]*i + j) + k];
  }

  // 返回第 (i, j, k) 个位置上的元素引用
  reference at(size_type const i, size_type const j, size_type const k)
  {
    return m_data.at(m_shape[2]*(m_shape[1]*i + j) + k);
  }

  // 返回第 (i, j, k) 个位置上的元素常引用
  const_reference at(size_type const i, size_type const j, size_type const k) const
  {
    return m_data.at(m_shape[2]*(m_shape[1]*i + j) + k);
  }

  size_type data_size() const { return m_data.size(); }


  // 预留内存 
  void reserve(const size_type m, const size_type n, const size_type p)
  {
    m_data.reserve(m*n*p);
  }

  bool empty() const { return m_data.empty(); }

  size_type shape(unsigned i) const 
  {
    assert(i >=0 && i < 3);
    return m_shape[i];
  }

  Shape & shape()
  {
    return m_shape;
  }

  // 返回原始数组
  value_type * row_data() {return m_data.data();}

  /*!
   * \brief 返回内部数据的引用，这里是 std::vector 数据的引用 
   */
  Data & data() {return m_data;}

private:

  Shape m_shape{0u, 0u, 0u};
  Data m_data{};
};

template<class T>
void swap(RArray3d<T> & lhs, RArray3d<T> & rhs)
{
  lhs.swap(rhs);
}

template<class T>
bool operator == (RArray3d<T> const &a, RArray3d<T> const &b)
{
  if (a.shape(0) != b.shape(0) || a.shape(1) != b.shape(1) || a.shape(2) != b.shape(2))
  {
    return false;
  }
  return std::equal(a.begin(), a.end(), b.begin(), b.end());
}

template<class T>
bool operator != (RArray3d<T> const &a, RArray3d<T> const &b)
{
  return !(a == b);
}

template<class T>
std::ostream& operator << (std::ostream & os, const RArray3d<T> & a)
{
  std::cout << "RArray3d("<<  a.shape(0) << "," << a.shape(1) << "," << a.shape(1) << "):" << std::endl;
  for(typename RArray3d<T>::size_type i = 0U; i < a.shape(0); i++)
  {
    for(typename RArray3d<T>::size_type j = 0U; j < a.shape(1); j++)
    {
      for(typename RArray3d<T>::size_type k = 0U; k < a.shape(2); k++)
      {
          os << a(i, j, k) << " "; 
      }
      os << "\n";
    }
    os << "\n";
  }
  return os;
}
} // end of namespace AlgebraObject

} // end of namespace TOPT 

#endif // end of RArray3d_h
