
#include <time.h>
#include <map>
#include <vector>
#include <iostream>

void test_vector(int n, int iter)
{
  std::vector<int> a(n, 3);
  clock_t S = clock();
  for(int i = 0; i < iter; i++)
  {
    for(int j = 0; j < n; j++)
      auto b = a[j]*1*2*3*4*5*6*7*8*9*10; 
  }
  clock_t E = clock();
  double runtime1 = double (E-S)/CLOCKS_PER_SEC;
  std::cout<< "vector 用时: " << runtime1 << " 秒 " <<std::endl; 
}

void test_map(int n, int iter)
{
  std::map<int, int> a;
  for(int i = 0; i < n; i++)
    a[i] = i;

  clock_t S = clock();
  for(int i = 0; i < iter; i++)
  {
    for(int j = 0; j < n; j++)
      auto b = a[j]; 
  }
  clock_t E = clock();
  double runtime1 = double (E-S)/CLOCKS_PER_SEC;
  std::cout<< "map 用时: " << runtime1 << " 秒 " <<std::endl; 
}

template<typename Arr>
void test_for(const Arr & a, const Arr & b, int * index, int l, int iter)
{
  clock_t S = clock();
  for(int k = 0; k < iter; k++)
  {
    for(int i = 0; i < l; i++)
    {
      for(int j = 0; j < l; j++)
      {
        if(a[i]==b[j])
        {
          index[i] = j;
          break;
        }
      }
    }
  }
  clock_t E = clock();
  double runtime1 = double (E-S)/CLOCKS_PER_SEC;
  std::cout<< "for 用时: " << runtime1 << " 秒 " <<std::endl; 
  for(int i = 0; i < l; i++)
    std::cout<< "index[i] = " << index[i] <<std::endl;;
}

void test_map0(int * a, int * b, int * index, int l, int iter)
{
  clock_t S = clock();
  for(int k = 0; k < iter; k++)
  {
    std::map<int, int> m;
    for(int i = 0; i < l; i++)
      m[a[i]] = i;
    for(int i = 0; i < l; i++)
      index[i] = m[b[i]];
  }
  clock_t E = clock();
  double runtime1 = double (E-S)/CLOCKS_PER_SEC;
  std::cout<< "map 用时: " << runtime1 << " 秒 " <<std::endl; 
  for(int i = 0; i < l; i++)
    std::cout<< "index[" << i << "] = " << index[i] <<std::endl;;
}

void test_double_array()
{
  std::cout<< "12dw" <<std::endl;
  std::map<std::array<double, 3> *, int> M;
  std::cout<< "123123dwsdw" <<std::endl;
  std::array<double, 3> a{0.0012012010201, 0.012031023401203120312, 2.412412412312323};
  std::cout<< a[0] <<std::endl;
  M[&a] = 1;
  std::array<double, 3> b(a);
  b[0] += 1e-130;
  M[&b] = 2;
  std::cout<< "M = " << M[&a] <<std::endl;
  std::cout<< "M = " << M[&b] <<std::endl;
}

int main(int argc, char ** argv)
{
  int n = std::stoi(argv[1]);
  int iter = std::stoi(argv[2]);
  //test_vector(n, iter);
  //test_map(n, iter);
  int a[3] = {10, 23, 34};
  int b[3] = {34, 23, 10};
  int index[3] = {0};
  //test_for(a, b, index, 3, iter);
  //test_map0(a, b, index, 3, iter);
  test_double_array();
}
