#include "NuECrossSection.h"
#include "MatrixBuilder.h"
#include "ElectronEnergyDistribution.h"
#include <iostream>
#include <TFile.h>


int main(int argc, char *argv[]){

  TString flux_file = argv[1];
  TString output_path = argv[2];
  TString output_name = argv[3];

 
  std::vector<double> bins = {0.8,2,3,5,7,9,20}; // Electron Energy Distribution

  EeDistributionCreator creator(flux_file,bins);    
  creator.LoadFluxHistos();
  creator.CreateXSecMatrices();
  creator.SetupEeDistributions();
  creator.PopulateEeDistributions();
 

  // Ensure calling this AFTER LoadFluxHistos(), CreateXSecMatrices(), SetupEeDistributions()

  creator.WriteEverythingToROOTFile(output_path, output_name);



/*
   std::cout << GetY(2,4) << std::endl;

   Couplings test = GetCouplings(14);
   std::cout << test.cLL << " " << test.cLR << " " << test.isAntinu << std::endl;
   std::cout << "Mandelstam s: "<< GetMandelstamInvariant_s(2) << std::endl;

   std::cout << "DSigmaDY, pdg=14: " << DSigmaDY(14,4,0.8) << std::endl;

   std::cout << "DSigmaDY, pdg=-14: " << DSigmaDY(-14,3,0.8) << std::endl;
   std::cout << "DSigmaDY, pdg=12: " << DSigmaDY(12,3,0.8) << std::endl;
   std::cout << "DSigmaDY, pdg=-12: " << DSigmaDY(-12,3,0.8) << std::endl;
   
   std::cout << "Simpson's rule for numerical integration: "<< IntegrateDSigmaDY(14, 5, 0, 0.5, 1000) << std::endl;

   std::cout << "Matrix Builder Testing Now!" << std::endl;
   MatrixK m = MatrixK("Test", 3, 0.0, 100.0, 4, 0., 100.);
   std::cout << "Ee Bins: " << m.GetNumEeBins() << std::endl;
//   m.SetNumEeBins(6);
//   std::cout << "New Ee Bins: " << m.GetNumEeBins() << std::endl;
//   std::cout << "New Enu Bins: " << m.GetNumEnuBins() << std::endl;

   MatrixK q = MatrixK("Test2", {2,3,4}, 6,0,10 );

   auto histo = m.GetMatrixHisto();
   auto histoq = q.GetMatrixHisto();
   TFile outfile("test.root", "RECREATE");
   outfile.cd();
   histo->Write();
   histoq->Write();
   q.PopulateMatrixHisto(14);
   histoq = q.GetMatrixHisto();
   histoq->Write("AfterPopulate");

   outfile.Close();
*/
   


   return 0.0;


}



