//user code
#include "TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"   //Needed to load TRootMCParticle & TRootJet, which is used in TFCreation.h

//Specific code for anomalous couplings analysis:
#include "TopTreeAnalysis/AnomCouplings/interface/TFCreation.h"
#include "TopTreeAnalysis/AnomCouplings/src/TFCreation.cc"

using namespace std;

int main ()
{
    TFile *fout = new TFile ("TransferFunctionsFromFile.root", "RECREATE");
    clock_t start = clock();
  
    cout << "***********************************************" << endl;  
    cout << " Beginning of the program for fitting the TF ! " << endl;
    cout << "*********************************************** \n" << endl;
  
    //Used classes
    TFCreation tfCreation;

    //Define all variables!
    const int NrFitHistos = 2;
    bool useROOTClass = false;
    bool useStartValues = true;
    int histoNrForStartValues = 0;
    bool useStartArray = true;
    TFile* file = new TFile("PlotsForTransferFunctions_AllEvts_UpdatedElAndMu.root","READ");

    //Define all histograms which need to be fitted!
    std::cout << " Will look at " << NrFitHistos << " different histograms to fit! \n" << std::endl;
    string Histo[NrFitHistos] = {"Mu_DiffPtVsGenPt","Light_DiffPtVsGenPt"};

    float startValues[6] = {1,2,3,4,5,6};

    for(int ii = 0; ii < NrFitHistos; ii++){
        TH2F* histoForFit = (TH2F*) file->Get( ("2D_histograms_graphs/"+Histo[ii]).c_str() );    
        tfCreation.CalculateTFFromFile(histoForFit, useStartValues, histoNrForStartValues, useROOTClass, useStartArray, startValues, fout);
    }
    
    //Close the root file where all histograms are saved!
    fout->Close();

    //Delete the used pointers:
    delete fout, file;
  
    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "           hasn't crashed yet ;-)           " << endl;
    cout << "********************************************" << endl;
  
    return 0;
}
