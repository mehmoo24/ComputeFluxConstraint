#include "../include/NuECrossSection.h"

#include <cmath>
#include <stdexcept>
#include "TMath.h"

double GetY(double Ee, double Enu){
   // Returns the y value corresponding to particular electron energy and neutrino energy, total electron energy being used
   return (Ee - kElectronMass) / Enu;
}

Couplings GetCouplings(int pdg){
   Couplings c;
   if(std::abs(pdg) == 14){
     c.cLL = -0.2723;
   }
   else if(std::abs(pdg) == 12) c.cLL = 0.7277;
   else throw std::runtime_error("Unsupported PDG!");
   c.cLR = 0.2277;
   c.isAntinu = (pdg < 0);
  
   return c;
}

double GetMandelstamInvariant_s(double Enu){
   return (kElectronMass*kElectronMass + 2*kElectronMass*Enu);
} // end GetMandelstamInvariant_s

double DSigmaDY(int pdg, double Enu, double y){
   Couplings c = GetCouplings(pdg);
   double factor = (GetMandelstamInvariant_s(Enu) * (kG_F * kG_F)) / kPi;
   double dSigma_dy;
   if (y < 0.0 || y > 1.0) throw std::runtime_error("y must be between 0 and 1, physically y = T_e/Enu");
   if (Enu <= 0.0) throw std::runtime_error("Enu must be positive");

   if (!c.isAntinu)
     dSigma_dy = factor * (c.cLL*c.cLL + (c.cLR*c.cLR)*((1-y)*(1-y)) - c.cLL*c.cLR*kElectronMass*y/Enu);
   else // antineutrino case
     dSigma_dy = factor * (c.cLR*c.cLR + (c.cLL*c.cLL)*((1-y)*(1-y)) - c.cLR*c.cLL*kElectronMass*y/Enu);  

   return dSigma_dy;
} // end DSigmaDY

double IntegrateDSigmaDY(int pdg, double Enu, double ymin, double ymax, int N){
  // Doing numerical integration using Simpson's rule
  if (N % 2 != 0) throw std::runtime_error("N must be even to use Simpson's rule!");

  double h = (ymax - ymin)/N;
  double f_min = DSigmaDY(pdg, Enu, ymin);
  double f_max = DSigmaDY(pdg, Enu, ymax);

  double sum_odd = 0.0;
  double sum_even = 0.0;

  for(int i=1; i<N; i++){
     // odd i
     if (i % 2 != 0) 
        sum_odd = sum_odd + 4*DSigmaDY(pdg, Enu, ymin+i*h);
     else // even i
        sum_even = sum_even + 2*DSigmaDY(pdg, Enu, ymin+i*h);
  }
  double result = h/3.0 * (f_min + f_max + sum_odd + sum_even);
  return result;

} // end IntegrateDSigmaDY
