
#include <iostream>
#include <cmath>
#include <iomanip>

#include "algebra/Matrix.h"
#include "algebra/Vector.h"

typedef OF::AlgebraObject::Matrix<double, int> Matrix;
typedef OF::AlgebraObject::Vector<double, int> Vector;

void test_copy_op()
{
  Matrix m = {{1.0, 2.0, 4.0}, {3.0, 4.0, 5.0}, {6.0, 7.0, 8.0}};
  Matrix n;
  n = m;
  std::cout<< m <<std::endl;
  std::cout<< n <<std::endl;
}

void test_transpose_mutiply()
{
  Matrix m = {{1.0, 2.0, 4.0}, {3.0, 4.0, 5.0}, {6.0, 7.0, 8.0}};
  Matrix n;
  n = m;
  std::cout<< m <<std::endl;
  std::cout<< n <<std::endl;
  std::cout<< n.transpose_multiply(m) <<std::endl;
  std::cout<< " adsadsa" <<std::endl;
}

int main()
{
  test_copy_op();
  test_transpose_mutiply();
  return 0;
}
