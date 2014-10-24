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

    //TFile *fout = new TFile ("TFInformation/CreatedTFFromDistributions.root", "RECREATE");
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
    bool RunFitForTF = true;

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
            //int nEvent = 10000;
            std::cout << " *** Looking at dataset " << iDataSet+1 << "/" << inputTFRoot.size() << " with " << nEvent << " selected events! \n " << std::endl;

            //Initialize the TFCreation class (create all histograms):
            TFCreation tfCreation;
            tfCreation.InitializeVariables(); //Add option of nr eta bins here!

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
                tfCreation.FillHistograms( &genPart[0], &genPart[1], &genPart[2], &genPart[3], &genPart[4], &recoPart[0], &recoPart[1], &recoPart[2], &recoPart[3], &recoPart[4], decayChannel);
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
        TDirectory* th2dir = writeFile->mkdir("2D_histograms_graphs");

        //Define all histograms which need to be fitted!
        const int NrFitHistos = 12;
        const int NrParamsDblGaus = 6;
        std::cout << " Will look at " << NrFitHistos << " different histograms to fit! " << std::endl;
        string Histo[NrFitHistos] = { "BJet_DiffPhiVsGenPt",
				      "BJet_DiffPtVsGenPt",
				      "BJet_DiffThetaVsGenPt",
				      "El_DiffPhiVsGenPt",
				      "El_DiffPtVsGenPt",
				      "El_DiffThetaVsGenPt",
				      "Light_DiffPhiVsGenPt",
				      "Light_DiffPtVsGenPt",
				      "Light_DiffThetaVsGenPt",
				      "Mu_DiffPhiVsGenInvPt",
				      "Mu_DiffInvPtVsGenInvPt",
				      "Mu_DiffThetaVsGenInvPt"};

	//Meaning of the six different variables:                                                    
        float StartValues[NrFitHistos][NrParamsDblGaus] = { {      0, 0.038, 77,   0.004,  0.011, 6.5},
                                                            {     -8,    18, 63,       0,    8.6, 4.1},
							    {      0, 0.038, 77,   0.004,  0.011, 6.5},
                                                            {      0, 0.038, 77,   0.004,  0.011, 6.5},
                                                            {     -8,    18, 63,       0,    8.6, 4.1},
                                                            {      0, 0.038, 77,   0.004,  0.011, 6.5},
                                                            {      0, 0.038, 77,   0.004,  0.011, 6.5},
                                                            {     -8,    18, 63,       0,    8.6, 4.1},
                                                            {      0, 0.038, 77,   0.004,  0.011, 6.5},
							    {    0.0,  0.01, 24,       0,  0.001,   4},
							    {-0.0008, 0.001, 24, -0.0001, 0.0001,   4} };
        float FitRangeDblGaus[NrFitHistos][2] ={ {-.15,0.15},{-30,40},{-0.1,0.1},{-0.02,0.02},{-4,6},{-0.02,0.02},{-0.15,0.15},{-40,40},{-0.1,0.1},{-0.03,0.03},{-0.004,0.004},{-0.02,0.02}};

        //Set the booleans!
        bool useROOTClass = false;
        bool useStartValues = true;
        int histoNrForStartValues = NrFitHistos; //Not needed if useStartValues = false
        bool useStartArray = true;
        bool changeFitRange = true;
 
        ofstream myTF, myTFCard;
        myTF.open("TFInformation/TransferFunctions_TABLE.txt");
        myTFCard.open("TFInformation/transfer_card_user.dat");
    		
        float startValues[NrParamsDblGaus], fitRangeDblGaus[2];
        for(int iHisto = 0; iHisto < NrFitHistos; iHisto++){
            TH2F* histoForFit = (TH2F*) readFile->Get( ("2D_histograms_graphs/"+Histo[iHisto]).c_str() );

            std::cout << " Updating y-axis range to : " << FitRangeDblGaus[iHisto][0] << " & " << FitRangeDblGaus[iHisto][1] << std::endl;
            if(changeFitRange){
                histoForFit->GetYaxis()->SetRangeUser(FitRangeDblGaus[iHisto][0],FitRangeDblGaus[iHisto][1]);
            }

            //Save the 2D histogram used for the fit!            	
            th2dir->cd();
            histoForFit->Write();
            writeFile->cd();
  
            //Set the correct startValues and fit the distribution
            for(int jj = 0; jj < NrParamsDblGaus; jj++) startValues[jj] = StartValues[iHisto][jj];
            tfCreation.CalculateTFFromFile(histoForFit, useStartValues, histoNrForStartValues, useROOTClass, useStartArray, startValues, writeFile);
    
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
    
            tfCreation.WriteTF(histoForFit, myTF, myTFCard);
    
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

