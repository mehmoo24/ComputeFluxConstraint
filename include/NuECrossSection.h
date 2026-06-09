// NuECrossSection.h

#ifndef NUECROSSSECTION_H
#define NUECROSSSECTION_H

// Constants
constexpr double kElectronMass = 0.000511; // GeV
constexpr double kG_F = 1.1663787E-5; // GeV^-2, fermi coupling constant
constexpr double kPi = 3.14159265358979323846;

struct Couplings{
  double cLL;
  double cLR;
  bool isAntinu; 
};

// Function Headers
double GetY(double Ee, double Enu); // Ee and Enu in GeV

Couplings GetCouplings(int pdg);

double GetMandelstamInvariant_s(double Enu); // Enu in GeV

double DSigmaDY(int pdg, double Enu, double y); // Enu in GeV

double IntegrateDSigmaDY(int pdg, double Enu, double ymin, double ymax, int N=1000); // for Simpson's rule, N must be even! Enu in GeV

#endif


