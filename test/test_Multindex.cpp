
#include <time.h>
#include <initializer_list>
#include <iostream>

#include "TestMacro.h"
#include "common/Multindex.h"

void test_C()
{
  ASSERT_THROW(C(3, 2)==3);
  ASSERT_THROW(C(4, 2)==6);
  ASSERT_THROW(C(5, 2)==10);
  ASSERT_THROW(C(5, 3)==10);
  ASSERT_THROW(C(6, 2)==15);
}

void test_mutindex()
{
  OF::Common::MultiIndex idx({0, 1, 2, 3});
  idx[2] = 100000;
  for(int i = 0; i < 4; i++)
    std::cout<< idx[i] <<std::endl;

  std::cout<< "he" <<std::endl;;
  std::vector<OF::Common::MultiIndex> idx1;
  std::cout<< "he12312" <<std::endl;;
  clock_t S = clock();
  OF::Common::MultiIndex::gen_all_index(4, 4, idx1);
  clock_t E = clock();
  double runtime1 = double (E-S)/CLOCKS_PER_SEC;
  std::cout<< "生成多重指标用时: " << runtime1 << " 秒 " <<std::endl; 

  std::cout<< " idx1 = " <<std::endl;;
  for(int j = 0; j < (int)idx1.size(); j++)
  {
    std::cout<< " [";
    for(int i = 0; i < 4; i++)
    {
      std::cout<< idx1[j][i];
      if(i < 3)
        std::cout<<", ";
    }
    std::cout<< "] : "<< idx1[j].to_number() <<std::endl;
  }
}

int main()
{
  test_C();
  std::cout<< "he" <<std::endl;;
  test_mutindex();
}
