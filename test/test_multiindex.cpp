#include "core/femcore.h"
#include "core/RArray2d.h"

int main(int argc, char ** argv)
{
  int p = std::stoi(argv[1]);
  int l = std::stoi(argv[2]);
  RArray2d<int> idx;
  OF::Core::get_all_multiindex(p, l, idx);
  std::cout<< idx << std::endl;
  int N = idx.shape(0);
  for(int i = 0; i < N; i++)
  {
    std::cout<< i << " " << OF::Core::multiindex_to_number(l, p, idx[i]) <<std::endl;
  }

}
