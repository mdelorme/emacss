#include <iostream>

#include "emacss.h"

int main() {
  dynamics d;
  node n;

  d.mynode = &n;
  n.galaxy.gamma = 1.0;
  n.galaxy.R = 10.0;
  n.galaxy.scale = 1.0;
  n.galaxy.M = 1.0;
  n.galaxy.type = 3;
  n.galaxy.f = 1.0;
  n.rh = 10.0;
  n.N = 4000;
  n.mm = 1.0;
  n.rj = n.r_jacobi(d);

  
  std::cout << d.dr2dt() << std::endl;
  
  return 0;
}
