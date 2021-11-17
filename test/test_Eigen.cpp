
#include <math.h>
#include <memory>
#include <vector>
#include <array>
#include <map>
#include <iostream>

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>


typedef Eigen::Triplet<double> Triplet;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> CSRMatrix;
typedef Eigen::VectorXd VectorXd;
typedef Eigen::ConjugateGradient<CSRMatrix, Eigen::Lower|Eigen::Upper> Solver;


int main()
{
  std::vector<Triplet> tlist(10);
  tlist[0] = Triplet(0, 0, 1.0);
  tlist[1] = Triplet(1, 1, 1.0);
  tlist[2] = Triplet(2, 2, 1.0);
  tlist[3] = Triplet(3, 3, 1.0);
  tlist[4] = Triplet(4, 4, 1.0);
  tlist[5] = Triplet(2, 4, -1.0);
  tlist[6] = Triplet(2, 3, -1.0);
  tlist[7] = Triplet(2, 1, -1.0);
  tlist[8] = Triplet(2, 0, -1.0);
  tlist[9] = Triplet(3, 2, 0);

  VectorXd x(5);
  VectorXd b(5);
  b[3] = -1.2;

  CSRMatrix M(5, 5);
  M.setFromTriplets(tlist.begin(), tlist.end());

  std::cout<< M <<std::endl;

  // 解方程
  Eigen::SparseLU<CSRMatrix, Eigen::COLAMDOrdering<int> > solver0;
  solver0.analyzePattern(M);
  solver0.factorize(M); /**< Compute the numerical factorization */
  x = solver0.solve(b); /**< Use the factors to solve the linear system */
}



