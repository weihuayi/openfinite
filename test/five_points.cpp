#include <iostream>
#include <math.h>
#include <algorithm>
#include <numeric>

#include "core/RArray2d.h"

struct PDE
{
  double a;
  double b;

  PDE(): a(0), b(M_PI){}

  double source(double x,double y)
  {
    double z = std::sin(x)*std::sin(y);
    return z;
  }

  double solution(double x, double y)
  {
    return 0.5*source(x, y);
  }
};

int main(int argc, char ** argv)
{
    argc++;
    int n = std::stoi(argv[0]);
    int N = (n+1)*(n+1);
    RArray2d<double> A(N, N);

    PDE pde;
    double h = (pde.b - pde.a)/n;

    for (int k = 0; k < N; k++) 
    {
      int i = k%n, j = k/n;
      if(i == 0 || i == n || j == 0 || j == n)
      {
        A(k, k) = 1;
      }
      else 
      {
        A(k, k+n) = 1/h/h;
        A(k, k-n) = A(k, k+n);
        A(k, k) = -4*A(k, k+n);
      }
    }
}
