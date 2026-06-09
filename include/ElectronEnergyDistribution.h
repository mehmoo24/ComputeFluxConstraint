// ElectronEnergyDistribution.h

#ifndef ELECTRONENERGYDISTRIBUTION_H
#define ELECTRONENERGYDISTRIBUTION_H

#include <TFile.h>
#include <TH2D.h>
#include "PlotUtils/MnvH1D.h"
#include "MatrixBuilder.h"

constexpr double kNumElectrons = 1.0; // number of electrons

class EeDistributionCreator{
   public:
      EeDistributionCreator(TString fluxFilePath, std::vector<double> EeBinEdges);      

      ~EeDistributionCreator(); // destructor

      void CreateXSecMatrices();

      TH2D* GetXSecMatrix(int pdg);

      void SetupEeDistributions();

      std::vector<double> GetNeutrinoEnergyBinning(int pdg);
       
      void WriteEverythingToROOTFile(TString path, TString output_name);

      void LoadFluxHistos();

      void PopulateEeDistributions();
      PlotUtils::MnvH1D* GetTotalEeDistribution();

   private:
      TString kfluxFilePath;
      // electron energy distributions for each flavour
      PlotUtils::MnvH1D* kEeDistr_PDG14;
      PlotUtils::MnvH1D* kEeDistr_PDGminus14;
      PlotUtils::MnvH1D* kEeDistr_PDG12;
      PlotUtils::MnvH1D* kEeDistr_PDGminus12;
      // flux distributions for each flavour
      PlotUtils::MnvH1D* kFluxDistr_PDG14;
      PlotUtils::MnvH1D* kFluxDistr_PDGminus14;
      PlotUtils::MnvH1D* kFluxDistr_PDG12;
      PlotUtils::MnvH1D* kFluxDistr_PDGminus12;
      // xsec matrices for each flavour
//      TH2D* kXSecMatrix_PDG14;
//      TH2D* kXSecMatrix_PDGminus14;
//      TH2D* kXSecMatrix_PDG12;
//      TH2D* kXSecMatrix_PDGminus12;

      MatrixK* kXSecMatrix_PDG14;
      MatrixK* kXSecMatrix_PDGminus14;
      MatrixK* kXSecMatrix_PDG12;
      MatrixK* kXSecMatrix_PDGminus12;

      std::vector<double> kEeBinEdges;

}; // end class EeDistributionCreator



#endif 
