#include <iostream>
#include "core/Bernstein.h"
#include "core/femcore.h"

typedef OF::FEM::Bernstein Bernstein;
typedef OF::FEM::Bernstein2d Bernstein2d;
typedef OF::FEM::Bernstein3d Bernstein3d;

void test_2d()
{
  RArray2d<int> multiindex;
  OF::Core::get_all_multiindex(100, 3, multiindex);

  Bernstein b(multiindex);
  Bernstein2d b2(multiindex);
  int ldof = multiindex.shape(0);
  double lambda[3] = {0.1, 0.2, 0.7};
  for(int i = 0; i < ldof; i++)
  {
    std::vector<double> v(ldof, 0);
    v[i] = 1;
    double a = b.de_casteljau(v.data(), lambda);
    std::vector<double> v1(ldof, 0);
    v1[i] = 1;
    double c = b2.de_casteljau(v1.data(), lambda);
    std::cout<< "i = " << i << " a = " << a << " c = " << c << "scale = " << (a-c)/a << std::endl;
  }
}

void test_3d()
{
  RArray2d<int> multiindex;
  OF::Core::get_all_multiindex(20, 4, multiindex);

  Bernstein b(multiindex);
  Bernstein3d b3(multiindex);
  int ldof = multiindex.shape(0);
  double lambda[4] = {0.1, 0.2, 0.3, 0.4};
  for(int i = 0; i < ldof; i++)
  {
    std::vector<double> v(ldof, 0);
    v[i] = 1;
    double a = b.de_casteljau(v.data(), lambda);
    std::vector<double> v1(ldof, 0);
    v1[i] = 1;
    double c = b3.de_casteljau(v1.data(), lambda);
    std::cout<< "i = " << i << " a = " << a << " c = " << c << "scale = " << (a-c)/a << std::endl;
  }
  clock_t s0 = clock();
  for(int i = 0; i < ldof; i++)
  {
    std::vector<double> v(ldof, 0);
    v[i] = 1;
    double a = b.de_casteljau(v.data(), lambda);
  }
  clock_t e0 = clock();
  std::cout<< "time = " <<(double)(e0-s0)/CLOCKS_PER_SEC <<std::endl;

  s0 = clock();
  for(int i = 0; i < ldof; i++)
  {
    std::vector<double> v(ldof, 0);
    v[i] = 1;
    double a = b3.de_casteljau(v.data(), lambda);
  }
  e0 = clock();
  std::cout<< "time = " <<(double)(e0-s0)/CLOCKS_PER_SEC <<std::endl;
}


int main()
{
  //test_2d();
  test_3d();
  return 0;
}



