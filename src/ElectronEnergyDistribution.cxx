#include "../include/ElectronEnergyDistribution.h"
#include "../include/MatrixBuilder.h"

EeDistributionCreator::EeDistributionCreator(TString fluxFilePath,std::vector<double> EeBinEdges)
   : kfluxFilePath(fluxFilePath), kXSecMatrix_PDG14(nullptr),
     kXSecMatrix_PDGminus14(nullptr), 
     kXSecMatrix_PDG12(nullptr),
     kXSecMatrix_PDGminus12(nullptr),
     kEeBinEdges(EeBinEdges){

   kEeDistr_PDG14 = nullptr;
   kEeDistr_PDGminus14 = nullptr;
   kEeDistr_PDG12 = nullptr;
   kEeDistr_PDGminus12 = nullptr;

   kFluxDistr_PDG14 = nullptr;
   kFluxDistr_PDGminus14 = nullptr;
   kFluxDistr_PDG12  = nullptr;
   kFluxDistr_PDGminus12 = nullptr;


} // end constructor

EeDistributionCreator::~EeDistributionCreator(){
   delete kXSecMatrix_PDG14;
   delete kXSecMatrix_PDGminus14;
   delete kXSecMatrix_PDG12;
   delete kXSecMatrix_PDGminus12;

   delete kEeDistr_PDG14;
   delete kEeDistr_PDGminus14;
   delete kEeDistr_PDG12;
   delete kEeDistr_PDGminus12;
}

PlotUtils::MnvH1D* EeDistributionCreator::GetTotalEeDistribution(){

   auto totalEeDistr = (PlotUtils::MnvH1D*)(kEeDistr_PDG14->Clone("kEeDistr_TOTAL"));
   totalEeDistr->Add(kEeDistr_PDGminus14);
   totalEeDistr->Add(kEeDistr_PDG12);
   totalEeDistr->Add(kEeDistr_PDGminus12);
   return totalEeDistr;
}
   


void EeDistributionCreator::CreateXSecMatrices(){
//   kXSecMatrix_PDG14 = new MatrixK("MatrixK_PDG14", {5,6,7,8}, 6,0,10 );
   auto bins_PDG14 = GetNeutrinoEnergyBinning(14);
   auto bins_PDGminus14 = GetNeutrinoEnergyBinning(-14);
   auto bins_PDG12 = GetNeutrinoEnergyBinning(12);
   auto bins_PDGminus12 = GetNeutrinoEnergyBinning(-12);

   // for raw pointers, delete any old objects first
   delete kXSecMatrix_PDG14;
   delete kXSecMatrix_PDGminus14;
   delete kXSecMatrix_PDG12;
   delete kXSecMatrix_PDGminus12; 

   kXSecMatrix_PDG14 = new MatrixK("kXSecMatrix_PDG14",bins_PDG14, kEeBinEdges);
   kXSecMatrix_PDGminus14 = new MatrixK("kXSecMatrix_PDGminus14", bins_PDGminus14, kEeBinEdges);
   kXSecMatrix_PDG12 = new MatrixK("kXSecMatrix_PDG12", bins_PDG12, kEeBinEdges);
   kXSecMatrix_PDGminus12 = new MatrixK("kXSecMatrix_PDGminus12", bins_PDGminus12, kEeBinEdges);

   kXSecMatrix_PDG14->PopulateMatrixHisto(14);
   kXSecMatrix_PDGminus14->PopulateMatrixHisto(-14);
   kXSecMatrix_PDG12->PopulateMatrixHisto(12);
   kXSecMatrix_PDGminus12->PopulateMatrixHisto(-12);

} // end CreateXSecMatrix

void EeDistributionCreator::SetupEeDistributions(){
   // raw pointers
   delete kEeDistr_PDG14;
   delete kEeDistr_PDGminus14;
   delete kEeDistr_PDG12;
   delete kEeDistr_PDGminus12;

   // setup the distributions now
   kEeDistr_PDG14 = new PlotUtils::MnvH1D("kEeDistr_PDG14", "kEeDistr_PDG14", kEeBinEdges.size()-1, kEeBinEdges.data());
   kEeDistr_PDGminus14 = new PlotUtils::MnvH1D("kEeDistr_PDGminus14", "kEeDistr_PDGminus14", kEeBinEdges.size()-1, kEeBinEdges.data());
   kEeDistr_PDG12 = new PlotUtils::MnvH1D("kEeDistr_PDG12", "kEeDistr_PDG12", kEeBinEdges.size()-1, kEeBinEdges.data());
   kEeDistr_PDGminus12 = new PlotUtils::MnvH1D("kEeDistr_PDGminus12", "kEeDistr_PDGminus12", kEeBinEdges.size()-1, kEeBinEdges.data());

} // end SetupEeDistributions

void EeDistributionCreator::PopulateEeDistributions(){
  double val = 0.0;
  int globalBin = 0.0;
  for (int j=1; j<=kEeDistr_PDG14->GetNbinsX(); j++){
   // loop over neutrino energy bins for that flavour
   for (int i=1; i<=kXSecMatrix_PDG14->GetMatrixHisto()->GetNbinsX(); i++){
     globalBin = kXSecMatrix_PDG14->GetMatrixHisto()->GetBin(i,j);
     val = val + kXSecMatrix_PDG14->GetMatrixHisto()->GetBinContent(globalBin) * kFluxDistr_PDG14->GetBinContent(i);
   } // end neutrino energy bins for loop
   kEeDistr_PDG14->SetBinContent(j, val);
   val = 0.0;
  } // end electron energy bins for loop

} // end PopulateEeDistributions


void EeDistributionCreator::LoadFluxHistos(){
   auto fluxFile = new TFile(kfluxFilePath, "READ");
   std::cout << "NOTE: INPUT CORRECT FLUX HISTO NAMES" << std::endl;
   kFluxDistr_PDG14 = (PlotUtils::MnvH1D*) fluxFile->Get("LE_FHC_minerva1_numu");
   kFluxDistr_PDGminus14 = (PlotUtils::MnvH1D*)fluxFile->Get("LE_FHC_minerva1_numu");
   kFluxDistr_PDG12 = (PlotUtils::MnvH1D*)fluxFile->Get("LE_FHC_minerva1_numu");
   kFluxDistr_PDGminus12 = (PlotUtils::MnvH1D*)fluxFile->Get("LE_FHC_minerva1_numu");

} // end LoadFluxHistos




std::vector<double> EeDistributionCreator::GetNeutrinoEnergyBinning(int pdg){
   std::vector<double> bins;
   TH1D* flux = nullptr;
   if (pdg == 14){ flux = kFluxDistr_PDG14; }
   else if (pdg == -14) { flux = kFluxDistr_PDGminus14; }
   else if (pdg == 12) { flux = kFluxDistr_PDG12; }
   else if (pdg == -12) { flux = kFluxDistr_PDGminus12; }

   if (!flux) return bins; // don't want to call GetNbins if nullptr

   for(int i=1;i<=flux->GetNbinsX()+1;i++){
         bins.push_back(flux->GetXaxis()->GetBinLowEdge(i));
   }
   return bins;
}


TH2D* EeDistributionCreator::GetXSecMatrix(int pdg){
   if(pdg == 14){return kXSecMatrix_PDG14->GetMatrixHisto();}
//   else if (pdg == -14){return kXSecMatrix_PDGminus14.GetMatrixHisto();}
//   else if (pdg == 12) {return kXSecMatrix_PDG12.GetMatrixHisto();}
//   else if (pdg == -12) {return kXSecMatrix_PDGminus12.GetMatrixHisto();}
   else {throw std::runtime_error("Invalid PDG!");}
}

void EeDistributionCreator::WriteEverythingToROOTFile(TString path, TString output_name){
   auto filename = path+output_name+".root";
   TFile f(filename, "RECREATE");
   f.cd();
   std::cout << "yes! " << std::endl;
   kFluxDistr_PDG14->Write("kFluxDistr_PDG14");  
   kFluxDistr_PDGminus14->Write("kFluxDistr_PDGminus14");
   kFluxDistr_PDG12->Write("kFluxDistr_PDG12");
   kFluxDistr_PDGminus12->Write("kFluxDistr_PDGminus12");

   kXSecMatrix_PDG14->GetMatrixHisto()->Write("kXSecMatrix_PDG14");
   kXSecMatrix_PDGminus14->GetMatrixHisto()->Write("kXSecMatrix_PDGminus14");
   kXSecMatrix_PDG12->GetMatrixHisto()->Write("kXSecMatrix_PDG12");
   kXSecMatrix_PDGminus12->GetMatrixHisto()->Write("kXSecMatrix_PDGminus12");   
   kEeDistr_PDG14->Write("kEeDistr_PDG14");   
   kEeDistr_PDGminus14->Write("kEeDistr_PDGminus14");
   kEeDistr_PDG12->Write("kEeDistr_PDG12");
   kEeDistr_PDGminus12->Write("kEeDistr_PDGminus12");

   auto kEeDistr_TOTAL = GetTotalEeDistribution();
   kEeDistr_TOTAL->Write("kEeDistr_TOTAL");
   delete kEeDistr_TOTAL;

   f.Close();

}
