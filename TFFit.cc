#include "TStyle.h"

//user code
//#include "PersonalClasses/Style.C"                                                 //CHECK if this works!
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "TopTreeAnalysisBase/MCInformation/interface/MCWeighter.h"
#include "TopTreeAnalysisBase/Selection/interface/ElectronPlotter.h"
#include "TopTreeAnalysisBase/Selection/interface/MuonPlotter.h"
#include "TopTreeAnalysisBase/Selection/interface/JetPlotter.h"
#include "TopTreeAnalysisBase/Selection/interface/VertexPlotter.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "TopTreeAnalysisBase/MCInformation/interface/ResolutionFit.h"
#include "TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "PersonalClasses/Style.C"                                                 //CHECK if this works!
#include "TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"
#include "TopTreeAnalysisBase/Tools/interface/LeptonTools.h"

//Specific code for anomalous couplings analysis:
#include "TopTreeAnalysis/AnomCouplings/interface/TFCreation.h"
#include "TopTreeAnalysis/AnomCouplings/src/TFCreation.cc"

using namespace std;

int main (int argc, char *argv[])
{
    TFile *fout = new TFile ("TransferFunctionsFromFile.root", "RECREATE");
    clock_t start = clock();
  
    cout << "***********************************************" << endl;  
    cout << " Beginning of the program for fitting the TF ! " << endl;
    cout << "*********************************************** \n" << endl;
  
    //SetStyle if needed
    //setMyStyle();

    /////////////////////
    //  Used classes   //
    /////////////////////  
    TFCreation tfCreation;

    //Define all variables!
    const int NrFitHistos = 2;
    bool useROOTClass = false;
    bool useStartValues = true;
    int histoNrForStartValues = 0;
    TFile* file = new TFile("PlotsForTransferFunctions_AllEvts_UpdatedElAndMu.root","READ");

    //Define all histograms which need to be fitted!
    std::cout << " \n Will look at " << NrFitHistos << " different histograms to fit! \n" << std::endl;
    string Histo[NrFitHistos] = {"Mu_DiffPtVsGenPt","Light_DiffPtVsGenPt"};

    for(int ii = 0; ii < NrFitHistos; ii++){
        TH2F* histoForFit = (TH2F*) file->Get( ("2D_histograms_graphs/"+Histo[ii]).c_str() );    
        tfCreation.CalculateTFFromFile(histoForFit, useStartValues, histoNrForStartValues, useROOTClass, fout);
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
