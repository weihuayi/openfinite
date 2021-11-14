#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <time.h>

#include "TestMacro.h"
#include "geometry/Point_2.h"

using namespace OF;

typedef GeometryObject::Point_2<double> Node;

void test_copy()
{
  Node n({0.0, 0.0});
  std::vector<Node> nodes;
  nodes.resize(100, n);
}

int main(int argc, char **argv)
{
  test_copy();
}
