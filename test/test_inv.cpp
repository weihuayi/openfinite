#include <math.h>
#include <iostream>
#include "algebra/Algebra_kernel.h"
#include "algebra/linalg.h"
#include "algebra/Matrix.h"

#include "TestMacro.h"

typedef OF::AlgebraObject::Matrix<double, int> Matrix;

void test_inv()
{

  Matrix mat({{0, 1, 2}, {2, 3, 4}, {1, 2, 12}});
  Matrix mat_inv({{-1.55555556,  0.44444444,  0.11111111},
                  { 1.11111111,  0.11111111, -0.22222222},
                  {-0.05555556, -0.05555556,  0.11111111}});
  OF::AlgebraAlgrithom::matInv(mat);
  ASSERT_THROW((mat-mat_inv).row_norm_l2(0)<1e-5);
  ASSERT_THROW((mat-mat_inv).row_norm_l2(1)<1e-5);
  ASSERT_THROW((mat-mat_inv).row_norm_l2(2)<1e-5);
}

int main()
{
  test_inv();
  return 0;
}

