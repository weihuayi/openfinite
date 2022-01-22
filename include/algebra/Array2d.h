#ifndef RArray2d_h
#define RArray2d_h

#include <vector>

template<class T>
class RArray2d 
{
public:
  // 基本类型 
  using data_type  = std::vector<T>;
  using value_type = typename std::vector<T>::value_type;
  using size_type  = typename std::vector<T>::size_type;

  // 引用 
  using reference       = typename std::vector<T>::reference;
  using const_reference = typename std::vector<T>::const_reference;

  // 迭代子 
  using iterator       = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;

  // 反向迭代子 
  using reverse_iterator       = typename std::vector<T>::reverse_iterator;
  using const_reverse_iterator = typename std::vector<T>::const_reverse_iterator;

  // 空构造函数 
  RArray2d() = default;

  // 默认插入 rows*cols 个值 
  RArray2d(size_type rows, size_type cols)
    : m_rows(rows), m_cols(cols), m_data(rows*cols)
  {}

  // 复制初始化二维数组 rows*cols
  RArray2d(size_type rows, size_type cols, const_reference val)
    : m_rows(rows), m_cols(cols), m_data(rows*cols, val)
  {}

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

  // 行迭代子 

  // 元素访问 (行优先)
  reference operator() (size_type const i,
    size_type const j)
  {
    return m_data[m_cols*i + j];
  }

  const_reference operator() (size_type const i,
    size_type const j) const
  {
    return m_data[m_cols*i + j];
  }

  reference at() (size_type const i, size_type const j)
  {
    return m_data.at(m_cols*i + j);
  }
  const_reference at() (size_type const i, size_type const j) const
  {
    return m_data.at(m_cols*i + j);
  }

  // resizing
  void resize(size_type new_rows, size_type new_cols)
  {
    // 新的二维数组  new_rows*new_cols
    Array2d tmp(new_rows, new_cols);

    // select smaller row and col size
    auto mc = std::min(m_cols, new_cols);
    auto mr = std::min(m_rows, new_rows);

    for (size_type i(0U); i < mr; ++i)
    {
      // iterators to begin of rows
      auto row = begin() + i*m_cols;
      auto tmp_row = tmp.begin() + i*new_cols;
      // move mc elements to tmp
      std::move(row, row + mc, tmp_row);
    }
    // move assignment to this
    *this = std::move(tmp);
  }

  size_type size() const { return m_data.size(); }
  size_type max_size() const { return m_data.max_size(); }

  bool empty() const { return m_data.empty(); }

  // 行数 
  size_type rows() const { return m_rows; }

  // 列数
  size_type cols() const { return m_cols; }

  // 数据交换 
  void swap(RArray2d & rhs)
  {
    using std::swap;
    m_data.swap(rhs.m_data);
    swap(m_rows, rhs.m_rows);
    swap(m_cols, rhs.m_cols);
  }

private:

  size_type m_rows{0u}; // 无符号的整型
  size_type m_cols{0u}; 
  data_type m_data{};
};

template<class T>
void swap(RArray2d<T> & lhs, RArray2d<T> & rhs)
{
  lhs.swap(rhs);
}

template<class T>
bool operator == (RArray2d<T> const &a, RArray2d<T> const &b)
{
  if (a.rows() != b.rows() || a.cols() != b.cols())
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
  auto rows = a.rows();
  auto cols = a.cols();

  std::cout << "RArray2d("<<  rows << "," << cols << "):" << std::endl;
  for(RArray2d<T>::size_type i = 0U; i < rows; i++)
  {
      for(RArray2d<T>::size_type j = 0U; j < cols; j++)
      {
          os << m(i, j) << " "; 
      }
      os << "\n";
  }
  os << "\n";
  return os;
}

} // end of namespace AlgebraObject

} // end of namespace OF

#endif // end of Array2d_h 
