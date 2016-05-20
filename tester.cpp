#include <iostream>

#include "emacss.h"

int main() {
  dynamics d;
  node n;

  d.mynode = &n;
  n.galaxy.gamma = 1;
  n.galaxy.R = 10.0;
  n.galaxy.scale = 1.0;
  n.galaxy.M = 1.0;
  n.galaxy.type = 3;
  n.galaxy.f = 1.0;
  n.rh = 10.0;
  n.N = 4000;
  n.mm = 1.0;
  n.s = 1;

  cout << "s = " << n.s << endl;
  
  n.galaxy.R = 1e-2;
  while (n.galaxy.R < 10.0) {
    n.galaxy.R += 1e-2;
    n.rj = n.r_jacobi(d);
    //cout << "dr2dt " << d.dr2dt() << endl;
    //cout << "rtide " << n.rj << endl;
  }
  //std::cout << d.dr2dt() << std::endl; d.fv(1e-12)
  
  return 0;
}
