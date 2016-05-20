#include <iostream>

#include "emacss.h"


int main() {
  dynamics d;
  node n;

  d.mynode = &n;
  n.galaxy.gamma = 0;
  n.galaxy.R = 1.0e3;         // pc
  n.galaxy.scale = 1.0e3;     // pc
  n.galaxy.M = 1.0e8;         // Msun
  n.galaxy.type = 3; 
  n.galaxy.f = 1.0;  
  n.rh = 10.0;                // pc
  n.N = 4000;                 
  n.mm = 1.0;                 // Msun
  n.s = 1;

  
  n.G_star = 0.00449857;                                     // Grav constant (pc^3M_sun^-1Myr^-2)


  // Cluster units
  
  n.M_star = n.mm*n.N;                                       // Initial Mass of cluster (M_sun)
  n.R_star = n.rh/0.78;                                      // Virial radii (parsec / N-body)
  n.T_star = sqrt(pow(n.R_star,3)/(n.M_star*n.G_star));        // N-body time (Myr)
  n.galaxy.R /= n.R_star;
  n.galaxy.M /= n.M_star;
  n.galaxy.scale /= n.R_star;
  n.rh /= n.R_star;
  n.mm = 1.0/n.N;

  // Galaxy units
  /*n.M_star = n.galaxy.M;                                       // Initial Mass of cluster (mass of the galaxy)
  n.R_star = n.galaxy.scale;                                   // Virial radii (parsec / N-body) (scale of the galaxy)
  n.T_star = sqrt(pow(n.R_star,3)/(n.M_star*n.G_star));        // N-body time (Myr)
  n.galaxy.R /= n.R_star;
  n.galaxy.M /= n.M_star;
  n.galaxy.scale /= n.R_star;
  n.rh /= n.R_star;
  n.mm /= n.M_star;*/
  
  d.dr2dt();
  
  return 0;
}
