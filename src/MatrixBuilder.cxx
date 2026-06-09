#include "../include/MatrixBuilder.h"
#include "../include/NuECrossSection.h"
#include <algorithm>


MatrixK::MatrixK(const char* matrixName, int numEnuBins, double Enu_min, double Enu_max, int numEeBins, double Ee_min, double Ee_max)
   : fnumEeBins(numEeBins), fnumEnuBins(numEnuBins) {

   matrixHisto = new TH2D(matrixName, matrixName, fnumEnuBins, Enu_min, Enu_max, fnumEeBins, Ee_min, Ee_max);

} // end of constructor

MatrixK::MatrixK(const char* matrixName, std::vector<Double_t> EnuBinEdges, std::vector<Double_t> EeBinEdges)
   : fnumEnuBins(EnuBinEdges.size()-1), fnumEeBins(EeBinEdges.size()-1) {
//   fnumEnuBins = EnuBinEdges.size()-1;
//   fnumEeBins = EeBinEdges.size()-1;
   matrixHisto = new TH2D(matrixName, matrixName, fnumEnuBins, EnuBinEdges.data(), fnumEeBins, EeBinEdges.data());
   
} // end of constructor w/ bin edges array

MatrixK::MatrixK(const char* matrixName, std::vector<Double_t> EnuBinEdges, int numEeBins, double Ee_min, double Ee_max)
   : fnumEnuBins(EnuBinEdges.size()-1), fnumEeBins(numEeBins) {
  
   // determine the bin edges for electron energy w/ info given
   std::vector<Double_t> EeBinEdges(numEeBins+1);
   double width = (Ee_max-Ee_min)/numEeBins;
   for(int i=0; i<=numEeBins; i++) EeBinEdges[i]=Ee_min+i*width;

   matrixHisto = new TH2D(matrixName, matrixName, fnumEnuBins, EnuBinEdges.data(), fnumEeBins, EeBinEdges.data());
}


MatrixK::~MatrixK(){ delete matrixHisto; } 

double MatrixK::GetNumEeBins() {return fnumEeBins;}

double MatrixK::GetNumEnuBins() {return fnumEnuBins;}

void MatrixK::SetNumEeBins(int numEeBins) {fnumEeBins = numEeBins;}

void MatrixK::SetNumEnuBins(int numEnuBins) {fnumEnuBins = numEnuBins;}

TH2D* MatrixK::GetMatrixHisto(){ return matrixHisto;}   

void MatrixK::PopulateMatrixHisto(int pdg, int N){
  int globalBin; double Enu; double ymin; double ymax; double val;
  double Emin, Emax;
  for (int x=1; x <= fnumEnuBins; x++){
     // use the bin center for the Enu value
     Enu = matrixHisto->GetXaxis()->GetBinCenter(x);
     for (int y=1; y <= fnumEeBins; y++){
        globalBin = matrixHisto->GetBin(x,y);
        Emin = matrixHisto->GetYaxis()->GetBinLowEdge(y);
        Emax = matrixHisto->GetYaxis()->GetBinUpEdge(y);
        ymin = GetY(Emin, Enu);
        ymax = GetY(Emax, Enu);

        // only integrate for the y values within the physical bounds: between 0 and 1:
 
        if (ymin >=1 || ymax <= 0) {val = 0.0;}
        else {
           ymin = std::max(ymin, 0.0);
           ymax = std::min(1.0, ymax);
           val = IntegrateDSigmaDY(pdg, Enu, ymin, ymax, N);}
        matrixHisto->SetBinContent(globalBin, val);
        matrixHisto->SetBinError(globalBin, 0.0);
     } // end looping over EeBins (ybins)
  } // end looping over EnuBins (xbins)

} // end PopulateMatrixHisto



