//user code
//#include "TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"   //Needed to load TRootMCParticle & TRootJet, which is used in TFCreation.h
#include "TFile.h"
#include "TH2.h"
#include <iostream>
#include <fstream>

#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TF2.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TBranch.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TDirectory.h"

//Specific code for anomalous couplings analysis:
#include "AnomalousCouplings/PersonalClasses/interface/TFCreation.h"
#include "AnomalousCouplings/PersonalClasses/interface/TFnTuple.h"

using namespace std;
//using namespace TopTree;

int main ()
{
    TFile *fout = new TFile ("TransferFunctionsFromFile.root", "RECREATE");
    clock_t start = clock();
  
    cout << "***********************************************" << endl;  
    cout << " Beginning of the program for fitting the TF ! " << endl;
    cout << "*********************************************** \n" << endl;
  
    //Used classes
    TFCreation tfCreation;

    ////////////////////////////////////////////////////////////////////
    //  Choose whether created plots are used or Tree information !!  //
    ////////////////////////////////////////////////////////////////////
    bool CreateTFFromTree = true;

    if(CreateTFFromTree){
        //Load the TFTree information
        vector<string> inputTFRoot;
        inputTFRoot.push_back("TFInformation/TransferFunctionTree.root");

        for(unsigned int iDataSet = 0; iDataSet <inputTFRoot.size(); iDataSet++){
            TFile* inputTFFile = new TFile(inputTFRoot[iDataSet].c_str(),"READ");

            TTree* inputTFTree = (TTree*) inputTFFile->Get("TFTree");
            TBranch* m_br = (TBranch*) inputTFTree->GetBranch("TheTFTree");

            int nEvent = inputTFTree->GetEntries(); 
            TFnTuple* tfNTuple = 0;
            m_br->SetAddress(&tfNTuple);

            //Get Pt of Quark1:
            for(unsigned int iEvt = 0; iEvt < nEvent; iEvt++){
                inputTFTree->GetEvent(iEvt);
                std::cout << " Pt of Quark1 (gen,reco) = " << (tfNTuple->genVectorLight1()).Pt() << " , " << (tfNTuple->recoVectorLight1()).Pt() << std::endl;
            }

            inputTFFile->Close();
            delete inputTFFile;
        }
    }
    else{
        //Define all needed variables:
        TFile* file = new TFile("PlotsForTransferFunctions_AllEvts_UpdatedElAndMu.root","READ");
        //Define all histograms which need to be fitted!
        const int NrFitHistos = 2;
        const int NrParamsDblGaus = 6;
        std::cout << " Will look at " << NrFitHistos << " different histograms to fit! " << std::endl;
        string Histo[NrFitHistos] = {"Mu_DiffPtVsGenPt","Light_DiffPtVsGenPt"};
        float StartValues[NrFitHistos][NrParamsDblGaus] = { {1,2,3,4,5,6},
                                                            {-10,0,500,0,10,10}};
        float FitRangeDblGaus[NrFitHistos][2] ={ {-10,10},{-5,5}};

        //Set the booleans!
        bool useROOTClass = false;
        bool useStartValues = true;
        int histoNrForStartValues = 0;
        bool useStartArray = true;
        bool changeFitRange = true;
 
        ofstream myTF, myTFCard;
        myTF.open("TFInformation/TransferFunctions_TABLE.txt");
        myTFCard.open("TFInformation/transfer_card_user.dat");
    
        float startValues[NrParamsDblGaus], fitRangeDblGaus[2];
        for(int iHisto = 0; iHisto < NrFitHistos; iHisto++){
            TH2F* histoForFit = (TH2F*) file->Get( ("2D_histograms_graphs/"+Histo[iHisto]).c_str() );    
        
            //Set the correct startValues and fit the distribution
            for(int jj = 0; jj < NrParamsDblGaus; jj++) startValues[jj] = StartValues[iHisto][jj];
            for(int jj = 0; jj < 2; jj++) fitRangeDblGaus[jj] = FitRangeDblGaus[iHisto][jj];
            tfCreation.CalculateTFFromFile(histoForFit, useStartValues, histoNrForStartValues, useROOTClass, useStartArray, startValues, changeFitRange, fitRangeDblGaus, fout);
    
            //Set the caption correct:
            string CaptionName, BlockName, PartName, KinVarName;
            // -- 1) which particle
            if(Histo[iHisto].find("Light_") == 0)    {PartName = "light jets"; BlockName = "TF_nonbjet_";}
            else if(Histo[iHisto].find("BJet") == 0) {PartName = "b-jets";     BlockName = "TF_bjet_";}
            else if(Histo[iHisto].find("Mu_") == 0)  {PartName = "muons";      BlockName = "TF_muon_";}
            else if(Histo[iHisto].find("El_") == 0)  {PartName = "electrons";  BlockName = "TF_electron_";}
            // -- 2) which kinematic variable
            if(Histo[iHisto].find("DiffPt") <= Histo[iHisto].size())         {CaptionName = PartName+" transverse momentum";            KinVarName += "PT";}
            else if(Histo[iHisto].find("DiffTheta") <= Histo[iHisto].size()) {CaptionName = PartName+" polar angle \\theta";            KinVarName += "THETA";}
            else if(Histo[iHisto].find("DiffPhi") <= Histo[iHisto].size())   {CaptionName = PartName+" azimuthal angle \\phi";          KinVarName += "PHI";}
            else if(Histo[iHisto].find("DiffInvPt") <= Histo[iHisto].size()) {CaptionName = PartName+" inverse of transverse momentum"; KinVarName += "InvPt";}
            BlockName = BlockName + KinVarName;
    
            //Write the TF's in a table and in a MadWeight card!:
            myTF<< endl;
            myTF<<" \n \\begin{table}[h!]" << endl;
            myTF<<"\\caption{Parameters of the transfer function for " << CaptionName << "}" << endl;
            myTF<<"\\label{tab::" << Histo[iHisto] << "}" << endl;
            myTF<<"\\centering" << endl;
            myTF<<"\\begin{tabular}{c|ccc}" << endl;
            myTF<<"\\hline" << endl;
            myTF << "Type      & $a_{i0}$ & $a_{i1}$ ($\\sqrt{E}$) & $a_{i2}$ ($E$)" << "\\\\" << endl;
            myTF<<"\\hline" << endl;
    
            //Only write this for the Pt- or 1/Pt-parameters (= start of the TF Card!)
            if(Histo[iHisto].find("DiffPt") <= Histo[iHisto].size() || Histo[iHisto].find("DiffInvPt") <= Histo[iHisto].size() ){
                myTFCard<<"#+-----------------------------------------------------------------------------------+" <<endl;
                myTFCard<<"#|     Parameter for particles: "<<PartName << endl; 
                myTFCard<<"#|      --> Used formula: Double Gaussian fit with parameters depending on momentum" << endl;
                myTFCard<<"#|      --> Dependency defined as: A + B*sqrt("<<KinVarName<<") + C*"<<KinVarName<< endl;
                myTFCard<<"#+-----------------------------------------------------------------------------------+" <<endl;
            }
            myTFCard<<"BLOCK "<<BlockName << endl;
    
            tfCreation.WriteTF(histoForFit, myTF, myTFCard, fout);
    
            myTF<<"\\hline" << endl;
            myTF<<"\\end{tabular}"<<endl;
            myTF<<"\\end{table} \n"<<endl;
        }
        
        //Close the root file where all histograms are saved together with the output files!
        fout->Close();
        myTF.close();
        myTFCard.close();
    
        //Delete the used pointers:
        delete fout, file;
    }//End of TF calculation when ROOT file is used!
  
    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "           hasn't crashed yet ;-)           " << endl;
    cout << "********************************************" << endl;
  
    return 0;
}
