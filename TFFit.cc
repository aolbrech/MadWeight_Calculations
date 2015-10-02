//user code
#include "TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"   //Needed to load TRootMCParticle & TRootJet, which is used in TFCreation.h
#include "TFile.h"
#include "TH2.h"
#include <iostream>
#include <fstream>
#include <cstdlib>

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
#include "AnomalousCouplings/PersonalClasses/interface/AnomCoupLight.h"

using namespace std;
//using namespace TopTree;

template <typename T> string tostr(const T& t) { ostringstream os; os<<t; return os.str(); }

int main (int argc, char **argv)
{
  TApplication theApp("App", &argc, argv); //Needed to run on local Linux!
  clock_t start = clock();

  gROOT->SetBatch();
  
  cout << "***********************************************" << endl;  
  cout << " Beginning of the program for fitting the TF ! " << endl;
  cout << "*********************************************** \n" << endl;
  
  ////////////////////////////////////////////////////////////////////
  //  Choose whether created plots are used or Tree information !!  //
  ////////////////////////////////////////////////////////////////////
  bool CreateTFFromTree = false;
  bool RunFitForTF = true; 
  int nEtaBins = 4;
  bool TFForPhi = false;
  bool TFForTheta = false;
  enum TFDependency_t {EDepTF, PtDepTF};
  TFDependency_t tfDependency = EDepTF;

  std::string EtaConsidered = "";
  if(nEtaBins == 4) EtaConsidered = "_etaBins";

  //Used classes
  TFCreation tfCreation(nEtaBins, EtaConsidered, RunFitForTF);

  if(CreateTFFromTree){
    //Load the TFTree information
    vector<string> inputTFRoot;
    inputTFRoot.push_back("LightTree/AnomCoupLight_TTbarJets_SemiLept_AllTTbarEvents_19Aug2015.root"); 

    for(unsigned int iDataSet = 0; iDataSet <inputTFRoot.size(); iDataSet++){
      TFile* inputTFFile = new TFile(inputTFRoot[iDataSet].c_str(),"READ");

      TTree* inputTFTree = (TTree*) inputTFFile->Get("LightTree");
      TBranch* m_br = (TBranch*) inputTFTree->GetBranch("TheAnomCoupLight");
      AnomCoupLight* light = 0;
      m_br->SetAddress(&light);

      //Set the number of selected events (for loop on events):
      int nEvent = inputTFTree->GetEntries(); 
      //int nEvent = 5000;
      std::cout << " *** Looking at dataset " << iDataSet+1 << "/" << inputTFRoot.size() << " with " << nEvent << " selected events! \n " << std::endl;

      //Initialize the TFCreation class (create all histograms):
      tfCreation.InitializeVariables(); 
      //Read in the TLorenztVectors:
      TLorentzVector genPart[5], recoPart[5];
      enum DecayChannel_t {isSemiMu, isSemiEl};
      DecayChannel_t decayChannel;
      for(unsigned int iEvt = 0; iEvt < nEvent; iEvt++){
	if(iEvt%10000 == 0)
	  std::cout<<"Processing the "<<iEvt<<"th event (" << ((double)iEvt/(double)inputTFTree->GetEntries())*100  << "%)" << flush<<"\r";

	inputTFTree->GetEvent(iEvt);
        vector<TLorentzVector> selectedJets = light->selectedJets();
        vector<int> correctJetCombi = light->correctJetCombi();
        int correctLeptBIndex = correctJetCombi[0];
        int correctHadrBIndex = correctJetCombi[1];
        int correctQuark1Index = correctJetCombi[2];
        int correctQuark2Index = correctJetCombi[3];
        
        if( correctJetCombi[0] != 9999 && correctJetCombi[1] != 9999 && correctJetCombi[2] != 9999 && correctJetCombi[3] != 9999){
    	  recoPart[0] = selectedJets[correctQuark1Index];
	  recoPart[1] = selectedJets[correctQuark2Index];
	  recoPart[2] = selectedJets[correctHadrBIndex];
	  recoPart[3] = selectedJets[correctLeptBIndex];
	  recoPart[4] = light->selectedLepton();

	  genPart[0] = light->genVectorLight1();
	  genPart[1] = light->genVectorLight2();
	  genPart[2] = light->genVectorHadrB();
	  genPart[3] = light->genVectorLeptB();
	  genPart[4] = light->genVectorLepton();
  
	  if(genPart[4].M() <= 0.05) decayChannel = isSemiEl; //Electron channel --> decayChannel == 1
	  else                       decayChannel = isSemiMu; //Muon     channel --> decayChannel == 0

	  //Fill the histograms of the TFCreation class!
	  tfCreation.FillHistograms( &genPart[0], &genPart[1], &genPart[2], &genPart[3], &genPart[4], &recoPart[0], &recoPart[1], &recoPart[2], &recoPart[3], &recoPart[4], decayChannel);

        }//Only for matched particles reconstructed!
      }//Loop on events

      TFile* fillFile = new TFile(("TFInformation/PlotsForTransferFunctions_FromLightTree"+EtaConsidered+".root").c_str(),"RECREATE");
      fillFile->cd();
      std::cout << "    ----> Information writen in file : " << fillFile->GetName() << std::endl << std::endl;
      tfCreation.WritePlots(fillFile);
      fillFile->Close();
      delete fillFile;

      inputTFFile->Close();
      delete inputTFFile;
    }
  }//CreateTFFromTree = true loop

  if(RunFitForTF == true){

    std::cout << " *** Starting to perform the double Gaussian fits  " << std::endl;

    //Set which TFFile should be used
    TFile *readFile, *writeFile;
    if(CreateTFFromTree == false){
      //readFile = new TFile("TFInformation/PlotsForTransferFunctions_AllEvts_LightTree.root","READ");
      writeFile = new TFile(("TFInformation/CreatedTFFromDistributions_AllEvts_LightTree"+EtaConsidered+".root").c_str(),"RECREATE");
      readFile = new TFile(("TFInformation/PlotsForTransferFunctions_FromLightTree"+EtaConsidered+".root").c_str(),"READ");
      //writeFile = new TFile("TFInformation/CreatedTFFromDistributions_FromLightTree.root","RECREATE");
    }
    else{
      readFile = new TFile(("TFInformation/PlotsForTransferFunctions_FromLightTree"+EtaConsidered+".root").c_str(),"READ");
      writeFile = new TFile(("TFInformation/CreatedTFFromDistributions_FromLightTree"+EtaConsidered+".root").c_str(),"RECREATE");
    }
    //Also draw the 2D histograms!	

    //Define all histograms which need to be fitted!
    const int NrFitHistos = 24;
    int ConsHisto = 1;
    const int NrParamsDblGaus = 5;
    std::cout << " --> Will look at " << NrFitHistos << " different histograms to fit! " << std::endl;
    string HistoInfo[NrFitHistos][1+NrParamsDblGaus] = { "BJet_DiffPhiVsGenE",        "0",      "0",    "1",     "20",     "25",
						         "BJet_DiffEVsGenE",          "0",     "10",    "1",     "60",     "65",
						         "BJet_DiffThetaVsGenE",      "0",     "10",    "1",      "0",   "0.02",
						         "El_DiffPhiVsGenE",          "0",  "0.006",   "0.",      "0", "0.0012",
						         "El_DiffEVsGenE",           "30",    "-20",   "0.",      "0",      "1",
						         "El_DiffThetaVsGenE",        "0",  "0.007",   "0.",      "0", "0.0013",
						         "Light_DiffPhiVsGenE",       "0",  "0.022",   "0.", "0.0004",  "0.002",
						         "Light_DiffEVsGenE",        "-1",     "12",  "0.2",    "-30",     "50",
						         "Light_DiffThetaVsGenE",     "0", " -0.05",   "0.",      "0", "-0.014",
						         "Mu_DiffPhiVsGenE",          "0", "0.0026",   "0.",      "0", "0.0004",
						         "Mu_DiffEVsGenE",           "30",     "18",    "1",      "1",      "1",
						         "Mu_DiffThetaVsGenE",        "0",  "0.002",   "0.",      "0", "0.0004",
                                                         /// ---         Transition from E-dependent to Pt-dependent      --- ///
						         "BJet_DiffPhiVsGenPt",   "0.002",  "0.022",   "0.", "0.0002",   "0.06",
						         "BJet_DiffPtVsGenPt",       "10",    "-12",   "0.",     "13",     "10",
						         "BJet_DiffThetaVsGenPt",     "0",   "0.04",   "0.",      "0",  "0.013",
						         "El_DiffPhiVsGenPt",         "0",  "0.006",   "0.",      "0", "0.0012",
						         "El_DiffPtVsGenPt",          "0",     "-2",   "0.",      "0",    "0.9",
						         "El_DiffThetaVsGenPt",       "0",  "0.007",   "0.",      "0", "0.0013",
						         "Light_DiffPhiVsGenPt",      "0",  "0.022",   "0.", "0.0004",  "0.002",
						         "Light_DiffPtVsGenPt",       "0",      "8",   "0.",      "0",     "12",
						         "Light_DiffThetaVsGenPt",    "0", " -0.05",   "0.",      "0", "-0.014",
						         "Mu_DiffPhiVsGenInvPt",      "0", "0.0026",   "0.",      "0", "0.0004",
						         "Mu_DiffInvPtVsGenInvPt",    "0", "0.0006",   "0.",      "0", "0.0003",
						         "Mu_DiffThetaVsGenInvPt",    "0",  "0.002",   "0.",      "0", "0.0004"};

    //Set the booleans!
    bool useROOTClass = false;
    bool useStartValues = true;
    int histoNrForStartValues = NrFitHistos; //Not needed if useStartValues = false
    bool useStartArray = true;
    bool changeFitRange = true;
    //Use a enum to distinguish between EDep and PtDep!! 
 
    ofstream myTFTable, myTransferCard, myTF, myLaTeX;
    myTFTable.open("TFInformation/TransferFunctions_TABLE.tex");
    if(nEtaBins == 1){
      myTF.open("TFInformation/TF_user.dat");
      myTransferCard.open("TFInformation/transfer_card_user.dat");
    }
    else{
      myTF.open("TFInformation/TF_user_etaBins.dat");
      myTransferCard.open("TFInformation/transfer_card_user_etaBins.dat");
    }
    myLaTeX.open(("TFInformation/TFSummary"+EtaConsidered+".tex").c_str());
    myLaTeX << "\\documentclass[a4paper,10pt]{article} " << endl;
    myLaTeX << "\\usepackage[utf8]{inputenc} " << endl;
    myLaTeX << "\\usepackage[margin=1in]{geometry} \n \\usepackage{graphicx} " << endl;
    myLaTeX << "\\title{Transfer Functions} \n\\begin{document} \n\\maketitle " << endl;
    myLaTeX << "\\begin{abstract} " << endl;
    myLaTeX << " Summary of the obtained Transfer Functions which will be implemented in Madweight as 'Double Gaussian TF' \\\\" << endl;
    myLaTeX << "\\end{abstract} " << endl;

    myTransferCard<<"#+-----------------------------------------------------------------------+"<<endl;
    myTransferCard<<"#|                         TRANSFER_CARD.DAT                             |"<<endl;
    myTransferCard<<"#|                                                                       |"<<endl;
    myTransferCard<<"#|     Author: Annik Olbrechts (VUB)                                     |"<<endl;
    myTransferCard<<"#|             02 October 2015                                           |"<<endl;
    myTransferCard<<"#+-----------------------------------------------------------------------+"<<endl;
    myTransferCard<<"#|     This file is generated automaticly by MADWEIGHT                   |"<<endl;
    myTransferCard<<"#|     card generation version: 1.0.0                                    |"<<endl;
    myTransferCard<<"#+-----------------------------------------------------------------------+"<<endl;
    myTransferCard<<"#|                                                                       |"<<endl;
    myTransferCard<<"#|    To change the transfer function run ./bin/change_tf.py             |"<<endl;
    if(nEtaBins == 1) myTransferCard<<"#|    Current parametrization : DblGaus_E                                |"<<endl;
    else              myTransferCard<<"#|    Current parametrization : DblGausEtaBins_E                        |"<<endl;
    myTransferCard<<"#|    Contains full double Gaussian for all kinematics and particles     |"<<endl;
    if(nEtaBins == 1) myTransferCard<<"#|    ** Information for all eta-bins **                                 |"<<endl;
    else              myTransferCard<<"#|    ** Information for the "<<nEtaBins<<" considered eta-bins separately **         |"<<endl;
    myTransferCard<<"#+-----------------------------------------------------------------------+"<<endl;	

    myTF<<"<file>## ##################################################################"<<endl;
    myTF<<"##                                          ##"<<endl;
    myTF<<"##                          Matrix Element                               ##"<<endl;
    myTF<<"##                          ==============                               ##"<<endl;
    myTF<<"##                                                                       ##"<<endl;
    myTF<<"##		    Generate the transfer functions                         ##"<<endl;
    myTF<<"##	             -------------------------------                        ##"<<endl;
    myTF<<"## 	     			                                            ##"<<endl;
    myTF<<"##				                                            ##"<<endl;
    myTF<<"##    Author: Annik Olbrechts (VUB)                                      ##"<<endl;
    myTF<<"##   			                                            ##"<<endl;
    myTF<<"##    Version:     1.0.0                         		            ##"<<endl;
    myTF<<"##    Last change: 02/02/15			                            ##"<<endl;
    myTF<<"##					                                    ##"<<endl;
    myTF<<"###########################################################################"<<endl;
    myTF<<"###########################################################################"<<endl;
    myTF<<"##				                                            ##"<<endl;
    myTF<<"##    Instructions:			                                    ##"<<endl;
    myTF<<"##								            ##"<<endl;
    myTF<<"##	- This program  creates transfer functions in THETA/PHI/E           ##"<<endl;
    myTF<<"##	- Those functions must be defined in f77 standard                   ##"<<endl;
    myTF<<"##	- In addition to each transfer function(START_TF), you MUST give    ##"<<endl;
    myTF<<"##	   the typical width associated to your function (START_WIDTH)	    ##"<<endl;
    myTF<<"##      - If a transfer functions is not defined here it will by default ##"<<endl;
    myTF<<"##          - equals to one for neutrino/neutralino                      ##"<<endl;
    myTF<<"##          - a delta functions for the rest                             ##"<<endl;
    myTF<<"###########################################################################"<<endl;
    myTF<<"##                                                                       ##"<<endl;
    myTF<<"##   Syntax/variables:                                                   ## "<<endl;
    myTF<<"##                                                                       ##"<<endl;
    myTF<<"##  - a definition for transfer functions should define the variable tf  ##"<<endl;
    myTF<<"##    while a definition for the width shoud define the variable width   ## "<<endl;                              
    myTF<<"##	- You can use all standard f77 functions. (All variables are	    ##"<<endl;
    myTF<<"##		in double precision format). 	                            ##"<<endl;
    myTF<<"##	- The experimental event is  defined by the variable pexp(i)        ##"<<endl;
    myTF<<"##		i=0->3 (0->E,1->Px,2->Py,3->Pz)	                            ##"<<endl;
    myTF<<"##	- The partonic event is defined by the variable p(i)   	            ##"<<endl;
    myTF<<"##		i=0->3 (0->E,1->Px,2->Py,3->Pz)			            ##"<<endl;
    myTF<<"##		sigma can not depend on those variables		            ##"<<endl;
    myTF<<"##	- You can use 10 local variables			            ##"<<endl;
    myTF<<"##		(double precision):  prov1,prov2,...,prov10	            ##"<<endl;
    myTF<<"##	- You can call specific functions on p or pexp:	                    ##"<<endl;
    myTF<<"##		-pt(p)   : transverse momenta 			            ##"<<endl;
    myTF<<"##		-eta(p)  : pseudo-rapidity			            ##"<<endl;
    myTF<<"##		-rap(p)  : rapidity					    ##"<<endl;
    myTF<<"##		-theta(p): polar angle				            ##"<<endl;
    myTF<<"##		-phi(p)	 : azimuthal angle		  	     	    ##"<<endl;
    myTF<<"##	- The whole LHCO information is available.                          ##"<<endl;
    myTF<<"##              -run_number,trigger                       	 	    ##"<<endl;
    myTF<<"##		-eta_init(N),phi_init(N),pt_init(N)                         ##"<<endl;
    myTF<<"##              -j_mass(N),ntrk(N),btag(N),had_em(N)     	            ##"<<endl;
    myTF<<"##		-dummy1(N),dummy2(N)		                            ##"<<endl;
    myTF<<"##	    N is the LHCO tag(first column)	                            ##"<<endl;
    myTF<<"##		- current tag is n_lhco  				    ##"<<endl;
    myTF<<"##		- tag for missing ET is met_lhco			    ##"<<endl;
    myTF<<"##				  					    ##"<<endl;
    myTF<<"##	- You can incorporate parameters that will be passed through        ##"<<endl;
    myTF<<"##	        the transfert_card.dat. Those ones must have the            ##"<<endl;
    myTF<<"##		following syntax: #1,#2,#3,.. You can restart 		    ##"<<endl;
    myTF<<"##		the assignement for each different transfer function	    ##"<<endl;
    myTF<<"##	- In addition to each transfer function(tf_), you MUST give	    ##"<<endl;
    myTF<<"##		the typical width associated to your function (sigma_)	    ##"<<endl;
    myTF<<"##		This is needed for the phase space generator      	    ##"<<endl;
    myTF<<"##									    ##"<<endl;
    myTF<<"###########################################################################"<<endl;
    myTF<<"###########################################################################"<<endl;

    float startValues[NrParamsDblGaus];
    for(int iHisto = 0; iHisto < NrFitHistos; iHisto++){
      int histoSize = HistoInfo[iHisto][0].size();
      std::string histoTitle = HistoInfo[iHisto][0];

      //Select the histograms which need to be considered!!
      if(!( ( (tfDependency == 1 && (histoTitle.find("VsGenPt") <= histoSize || histoTitle.find("VsGenInvPt") <= histoSize ) ) || (tfDependency == 0 && histoTitle.find("VsGenE") <= histoSize ) ) 
          && ( ( TFForPhi && histoTitle.find("DiffPhi") <= histoSize ) || (TFForTheta && histoTitle.find("DiffTheta") <= histoSize ) || histoTitle.find("DiffE") <= histoSize || histoTitle.find("DiffInvE") <= histoSize ) ) )
        continue;
  
      //Set the correct startValues and fit the distribution
      for(int jj = 0; jj < NrParamsDblGaus; jj++) startValues[jj] = atof((HistoInfo[iHisto][1+jj]).c_str());

      std::cout << "  ** Caclulating TF for histo : " <<  HistoInfo[iHisto][0] << std::endl;
      for(int iEtaBin = 1; iEtaBin <= nEtaBins; iEtaBin++)
	tfCreation.CalculateTFFromFile(HistoInfo[iHisto][0], useStartValues, histoNrForStartValues, useROOTClass, useStartArray, startValues, changeFitRange, writeFile, iEtaBin, readFile);
    
      //Set the caption correct:
      string CaptionName, BlockName, PartName, KinVarName, particles, widthType;
      // -- 1) which particle
      if(     histoTitle.find("Light_") == 0) {PartName = "nonbjet";  BlockName = "TF_nonbjet_";  particles = "u,d,s,c,g"; widthType = "thin";}
      else if(histoTitle.find("BJet")   == 0) {PartName = "bjet";     BlockName = "TF_bjet_";     particles = "b";         widthType = "thin";}
      else if(histoTitle.find("Mu_")    == 0) {PartName = "muon";     BlockName = "TF_muon_";     particles = "mu";        widthType = "thin";}
      else if(histoTitle.find("El_")    == 0) {PartName = "electron"; BlockName = "TF_electron_"; particles = "el";        widthType = "thin";}
      // -- 2) which kinematic variable
      if(     histoTitle.find("DiffE") <=     histoSize) {CaptionName = PartName+" transverse momentum";            KinVarName = "E";   }
      else if(histoTitle.find("DiffTheta") <= histoSize) {CaptionName = PartName+" polar angle \\theta";            KinVarName = "THETA";}
      else if(histoTitle.find("DiffPhi") <=   histoSize) {CaptionName = PartName+" azimuthal angle \\phi";          KinVarName = "PHI";  }
      else if(histoTitle.find("DiffInvE") <=  histoSize) {CaptionName = PartName+" inverse of transverse momentum"; KinVarName = "E";   } //MadWeight only knows 'PT'
      BlockName = BlockName + KinVarName;
    
      //Write the TF's in a table and in a MadWeight card!:
      myTFTable<< endl;
      myTFTable<<" \n \\begin{table}[h!]" << endl;
      myTFTable<<"\\caption{Parameters of the transfer function for " << CaptionName << "}" << endl;
      myTFTable<<"\\label{tab::" << HistoInfo[iHisto][0] << "}" << endl;
      myTFTable<<"\\centering" << endl;
      myTFTable<<"\\begin{tabular}{c|ccc}" << endl;
      myTFTable<<"\\hline" << endl;
      if(tfDependency == 0 || (tfDependency == 1 && PartName != "muon") ) myTFTable<< "Type      & $a_{i0}$ & $a_{i1}$ ($\\sqrt{E}$) & $a_{i2}$ ($E$)" << "\\\\" << endl;
      else                                                                myTFTable<< "Type      & $a_{i0}$ & $a_{i1}$ ($\\sqrt{\\frac{1}{E}}$) & $a_{i2}$ ($\\frac{1}{E}$)" << "\\\\" << endl;
      myTFTable<<"\\hline" << endl;
    
      if( (TFForPhi == false && TFForTheta == false) || 
          (TFForTheta == true && TFForPhi == false && histoTitle.find("DiffE") <= histoSize) || 
          (TFForTheta == true && TFForPhi == true  && histoTitle.find("DiffPhi") <= histoSize) ){   //Only write this for the first variable! 

	myTransferCard<<"#+--------------------------------------------------------------------------------------+" <<endl;
	myTransferCard<<"#|     Parameter for particles: "<<PartName << endl; 
	myTransferCard<<"#|      --> Used formula: Double Gaussian fit with parameters depending on momentum" << endl;
	if(tfDependency == 0 || (tfDependency == 1 && PartName != "muon") ) myTransferCard<<"#|      --> Dependency defined as: A + B*sqrt(E) + C*E  for width of gaussians "<< endl;
	else                                                                myTransferCard<<"#|      --> Dependency defined as: A + B*sqrt(1/E) + C*1/E  for width of gaussians "<< endl;
	if(tfDependency == 0 || (tfDependency == 1 && PartName != "muon") ) myTransferCard<<"#|      -->                        A + B*E + C*E² + D*E³ + F*sqrt(E)"<< endl;
	else                                                                myTransferCard<<"#|      -->                        A + B*1/E + C*1/E² + D*1/E³ + F*1/E^4"<< endl;
	myTransferCard<<"#+--------------------------------------------------------------------------------------+" <<endl;

	myTF<<"\n##**********************************************************************##"<<endl;
	myTF<<"##             TF for "<<PartName<<"                                      "<<endl;
	myTF<<"##**********************************************************************##"<<endl;
	myTF<<"<block name='"<<PartName<<"'>   #name can be anything"<<endl;
	if(tfDependency == 0 || (tfDependency == 1 && PartName != "muon") ) myTF<<"  <info> double gaussian with parameter depending on the energy </info>"<<endl;
	else                                                                myTF<<"  <info> double gaussian with parameter depending on the inverse of energy </info>"<<endl;
	myTF<<"  <particles> "<<particles<<" </particles>"<<endl;
	myTF<<"  <width_type> "<<widthType<<" </width_type>"<<endl;
	myTF<<"  #width type should be thin or large (thin is for energy accurate up to 5-10%)";
      }

      //Output for all variables!
      myTransferCard<<"BLOCK "<<BlockName << endl;
                
      myTF<<"\n  <variable name='"<<KinVarName<<"'>"<<endl;
      myTF<<"    <tf>";
  
      tfCreation.WriteTF(myTFTable, myTransferCard, myTF, myLaTeX, KinVarName, PartName, tfDependency);
      if(TFForTheta == false || (TFForTheta == true && histoTitle.find("DiffTheta") <= histoSize) ){
	if( (tfDependency == 0 && histoTitle.find("Mu_DiffE") <= histoSize) || (tfDependency == 1 && histoTitle.find("Mu_DiffInvE") <= histoSize) ){ myTF << "\n</block>\n</file>";}  //Muon is the last of the 4 particles!!
	else{                                                                                                                                        myTF << "\n</block>";         }
      }
    
      myTFTable<<"\\hline" << endl;
      myTFTable<<"\\end{tabular}"<<endl;
      myTFTable<<"\\end{table} \n"<<endl;

    }             
    //Close the root file where all histograms are saved together with the output files!
    readFile->Close();
    writeFile->Close();
    myTF.close();
    myTransferCard.close();
    myTFTable.close();

    myLaTeX << "\\end{document}"<<endl;
    myLaTeX.close();
    
    //Delete the used pointers:
    delete readFile;
    delete writeFile;
  }//End of TF calculation when ROOT file is used!

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "           hasn't crashed yet ;-)           " << endl;
  cout << "********************************************" << endl;
  
  return 0;
}
