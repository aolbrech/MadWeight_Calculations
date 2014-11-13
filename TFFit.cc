//user code
#include "TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"   //Needed to load TRootMCParticle & TRootJet, which is used in TFCreation.h
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
#include <TApplication.h>

//Specific code for anomalous couplings analysis:
#include "AnomalousCouplings/PersonalClasses/interface/TFCreation.h"
#include "AnomalousCouplings/PersonalClasses/interface/TFnTuple.h"

using namespace std;
//using namespace TopTree;

int main (int argc, char **argv)
{
    TApplication theApp("App", &argc, argv); //Needed to run on local Linux!
    clock_t start = clock();
  
    cout << "***********************************************" << endl;  
    cout << " Beginning of the program for fitting the TF ! " << endl;
    cout << "*********************************************** \n" << endl;
  
    ////////////////////////////////////////////////////////////////////
    //  Choose whether created plots are used or Tree information !!  //
    ////////////////////////////////////////////////////////////////////
    bool CreateTFFromTree = true;
    bool RunFitForTF = true;
    int nEtaBins = 4;
    //Used classes
    TFCreation tfCreation(nEtaBins);

    if(CreateTFFromTree){
        //Load the TFTree information
        vector<string> inputTFRoot;
        inputTFRoot.push_back("TFInformation/TransferFunctionTree_AllEvts.root");

        for(unsigned int iDataSet = 0; iDataSet <inputTFRoot.size(); iDataSet++){
            TFile* inputTFFile = new TFile(inputTFRoot[iDataSet].c_str(),"READ");

            TTree* inputTFTree = (TTree*) inputTFFile->Get("TFTree");
            TBranch* m_br = (TBranch*) inputTFTree->GetBranch("TheTFTree");
            TFnTuple* tfNTuple = 0;
            m_br->SetAddress(&tfNTuple);

            //Set the number of selected events (for loop on events):
            int nEvent = inputTFTree->GetEntries(); 
            //int nEvent = 1000;
            std::cout << " *** Looking at dataset " << iDataSet+1 << "/" << inputTFRoot.size() << " with " << nEvent << " selected events! \n " << std::endl;

            //Initialize the TFCreation class (create all histograms):
            tfCreation.InitializeVariables(nEtaBins); //Add option of nr eta bins here!

            //Read in the TLorenztVectors:
            TLorentzVector genPart[5], recoPart[5];
            enum DecayChannel_t {isSemiMu, isSemiEl};
            DecayChannel_t decayChannel;
            for(unsigned int iEvt = 0; iEvt < nEvent; iEvt++){
                if(iEvt%10000 == 0)
                    std::cout<<"Processing the "<<iEvt<<"th event (" << ((double)iEvt/(double)inputTFTree->GetEntries())*100  << "%)" << flush<<"\r";

                inputTFTree->GetEvent(iEvt);
                genPart[0] = tfNTuple->genVectorLight1();
                genPart[1] = tfNTuple->genVectorLight2();
                genPart[2] = tfNTuple->genVectorHadrB();
                genPart[3] = tfNTuple->genVectorLeptB();
                genPart[4] = tfNTuple->genVectorLepton();

                recoPart[0] = tfNTuple->recoVectorLight1();
                recoPart[1] = tfNTuple->recoVectorLight2();
                recoPart[2] = tfNTuple->recoVectorHadrB();
                recoPart[3] = tfNTuple->recoVectorLeptB();
                recoPart[4] = tfNTuple->recoVectorLepton();

                if(genPart[4].M() <= 0.05) decayChannel = isSemiEl; //Electron channel --> decayChannel == 1
                else                       decayChannel = isSemiMu; //Muon     channel --> decayChannel == 0

                //Fill the histograms of the TFCreation class!
                tfCreation.FillHistograms( &genPart[0], &genPart[1], &genPart[2], &genPart[3], &genPart[4], &recoPart[0], &recoPart[1], &recoPart[2], &recoPart[3], &recoPart[4], decayChannel, nEtaBins);
            }//Loop on events

            TFile* fillFile = new TFile("TFInformation/PlotsForTransferFunctions_FromTree.root","RECREATE");
            std::cout << "    ----> Information writen in file : " << fillFile->GetName() << std::endl << std::endl;
            tfCreation.WritePlots(fillFile);
            fillFile->Close();
            delete fillFile;

            inputTFFile->Close();
            delete inputTFFile;
        }
    }//CreateTFFromTree = true loop

    if(RunFitForTF == true){

        std::cout << " *** Starting to perform the double Gaussian fits  \n " << std::endl;

        //Set which TFFile should be used
        TFile *readFile, *writeFile;
        if(CreateTFFromTree == false){
		readFile = new TFile("TFInformation/PlotsForTransferFunctions_AllEvts_UpdatedElAndMu.root","READ");
		writeFile = new TFile("TFInformation/CreatedTFFromDistributions_AllEvts_UpdatedElAndMu.root","RECREATE");
	}
        else{
		readFile = new TFile("TFInformation/PlotsForTransferFunctions_FromTree.root","READ");
		writeFile = new TFile("TFInformation/CreatedTFFromDistributions_FromTree.root","RECREATE");
	}
	//Also draw the 2D histograms!	

        //Define all histograms which need to be fitted!
        const int NrFitHistos = 12;
        int ConsHisto = 1;
        const int NrParamsDblGaus = 6;
        std::cout << " Will look at " << NrFitHistos << " different histograms to fit! " << std::endl;
        string HistoInfo[12][1+NrParamsDblGaus+2] = { "BJet_DiffPhiVsGenPt",    "0.0002", "0.022", "8000",  "0.0002",   "0.06", "3000",  "-0.12",  "0.12",
                                                      "BJet_DiffPtVsGenPt",         "10",   "-12","20000",      "13",    "10", "-5000",    "-30",    "45",
				                      "BJet_DiffThetaVsGenPt",       "0",  "0.04", "2000",       "0",  "0.013", "6000",   "-0.1",   "0.1",
				                      "El_DiffPhiVsGenPt",           "0", "0.006",  "600",       "0", "0.0012", "1500", "-0.012", "0.012",
				                      "El_DiffPtVsGenPt",            "0",    "-2",  "600",       "0",    "0.9", "1500",     "-4",     "5",
				                      "El_DiffThetaVsGenPt",         "0", "0.007",  "700",       "0", "0.0013", "2500", "-0.018", "0.018",
				                      "Light_DiffPhiVsGenPt",        "0", "0.022", "8000",  "0.0004",  "0.002", "3000",  "-0.14",  "0.14",
				                      "Light_DiffPtVsGenPt",         "0",     "8", "4000",       "0",     "12", "4000",    "-28",    "35",
				                      "Light_DiffThetaVsGenPt",      "0", "-0.05", "2000",       "0", "-0.014", "6000",  "-0.12",  "0.12",
				                      "Mu_DiffPhiVsGenInvPt",        "0", "0.008",  "700",       "0", "0.0015",  "600", "-0.004", "0.004",
				                      "Mu_DiffInvPtVsGenInvPt",      "0","0.0003", "2000",       "0", "0.0006",  "500","-0.0015", "0.001",
				                      "Mu_DiffThetaVsGenInvPt",      "0", "0.002",  "500",       "0", "0.0004",  "500","-0.0035","0.0035"};

        //Set the booleans!
        bool useROOTClass = false;
        bool useStartValues = true;
        int histoNrForStartValues = NrFitHistos; //Not needed if useStartValues = false
        bool useStartArray = true;
        bool changeFitRange = true;
 
        ofstream myTF, myTFCard;
        myTF.open("TFInformation/TransferFunctions_TABLE.txt");
        myTFCard.open("TFInformation/transfer_card_user.dat");
    		
        float startValues[NrParamsDblGaus], fitRanges[2];
        for(int iHisto = 0; iHisto < NrFitHistos; iHisto++){
            if(NrFitHistos == 1) iHisto = ConsHisto;

            //Set the correct startValues and fit the distribution
            for(int jj = 0; jj < NrParamsDblGaus; jj++) startValues[jj] = std::stof(HistoInfo[iHisto][1+jj]);
            for(int jj = 0; jj < 2; jj++) fitRanges[jj] = std::stof(HistoInfo[iHisto][NrParamsDblGaus+1+jj]);
            std::cout << " Values from HistoInfo : " << HistoInfo[iHisto][NrParamsDblGaus+1] << " & " << HistoInfo[iHisto][NrParamsDblGaus+2] << std::endl;
            std::cout << " Fit range set to : " << fitRanges[0] << " & " << fitRanges[1] << " --> in analyzer ! " << std::endl;
    
            for(int iEtaBin = 0; iEtaBin <= nEtaBins; iEtaBin++)
                tfCreation.CalculateTFFromFile(HistoInfo[iHisto][0], useStartValues, histoNrForStartValues, useROOTClass, useStartArray, startValues, changeFitRange, fitRanges, writeFile, iEtaBin, readFile);
    
            //Set the caption correct:
            string CaptionName, BlockName, PartName, KinVarName;
            // -- 1) which particle
            if(HistoInfo[iHisto][0].find("Light_") == 0)    {PartName = "light jets"; BlockName = "TF_nonbjet_";}
            else if(HistoInfo[iHisto][0].find("BJet") == 0) {PartName = "b-jets";     BlockName = "TF_bjet_";}
            else if(HistoInfo[iHisto][0].find("Mu_") == 0)  {PartName = "muons";      BlockName = "TF_muon_";}
            else if(HistoInfo[iHisto][0].find("El_") == 0)  {PartName = "electrons";  BlockName = "TF_electron_";}
            // -- 2) which kinematic variable
            if(HistoInfo[iHisto][0].find("DiffPt") <= HistoInfo[iHisto][0].size())         {CaptionName = PartName+" transverse momentum";            KinVarName += "PT";}
            else if(HistoInfo[iHisto][0].find("DiffTheta") <= HistoInfo[iHisto][0].size()) {CaptionName = PartName+" polar angle \\theta";            KinVarName += "THETA";}
            else if(HistoInfo[iHisto][0].find("DiffPhi") <= HistoInfo[iHisto][0].size())   {CaptionName = PartName+" azimuthal angle \\phi";          KinVarName += "PHI";}
            else if(HistoInfo[iHisto][0].find("DiffInvPt") <= HistoInfo[iHisto][0].size()) {CaptionName = PartName+" inverse of transverse momentum"; KinVarName += "INVPT";}
            BlockName = BlockName + KinVarName;
    
            //Write the TF's in a table and in a MadWeight card!:
            myTF<< endl;
            myTF<<" \n \\begin{table}[h!]" << endl;
            myTF<<"\\caption{Parameters of the transfer function for " << CaptionName << "}" << endl;
            myTF<<"\\label{tab::" << HistoInfo[iHisto][0] << "}" << endl;
            myTF<<"\\centering" << endl;
            myTF<<"\\begin{tabular}{c|ccc}" << endl;
            myTF<<"\\hline" << endl;
            myTF << "Type      & $a_{i0}$ & $a_{i1}$ ($\\sqrt{E}$) & $a_{i2}$ ($E$)" << "\\\\" << endl;
            myTF<<"\\hline" << endl;
    
            //Only write this for the Pt- or 1/Pt-parameters (= start of the TF Card!)
            if(HistoInfo[iHisto][0].find("DiffPt") <= HistoInfo[iHisto][0].size() || HistoInfo[iHisto][0].find("DiffInvPt") <= HistoInfo[iHisto][0].size() ){
                myTFCard<<"#+-----------------------------------------------------------------------------------+" <<endl;
                myTFCard<<"#|     Parameter for particles: "<<PartName << endl; 
                myTFCard<<"#|      --> Used formula: Double Gaussian fit with parameters depending on momentum" << endl;
                myTFCard<<"#|      --> Dependency defined as: A + B*sqrt("<<KinVarName<<") + C*"<<KinVarName<< endl;
                myTFCard<<"#+-----------------------------------------------------------------------------------+" <<endl;
            }
            myTFCard<<"BLOCK "<<BlockName << endl;
    
            //tfCreation.WriteTF(HistoInfo[iHisto][0], myTF, myTFCard, whichEtaBin);
    
            myTF<<"\\hline" << endl;
            myTF<<"\\end{tabular}"<<endl;
            myTF<<"\\end{table} \n"<<endl;
        }             
        //Close the root file where all histograms are saved together with the output files!
        readFile->Close();
        writeFile->Close();
        myTF.close();
        myTFCard.close();
    
        //Delete the used pointers:
        delete readFile,writeFile;
    }//End of TF calculation when ROOT file is used!
  
    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "           hasn't crashed yet ;-)           " << endl;
    cout << "********************************************" << endl;
  
    return 0;
}
