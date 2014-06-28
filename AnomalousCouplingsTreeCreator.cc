///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// Analysis Skeleton.cc: This macro is intended to be an example analysis macro which works out of the box. /////
/////       It should serve as the first port of call for new users of the TopBrussels framework.              /////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

//user code
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
#include "Style.C"
#include "TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"

//Specific code for anomalous couplings analysis:
#include "TopTreeAnalysis/AnomCouplings/interface/LHCOOutput.h"
#include "TopTreeAnalysis/AnomCouplings/src/LHCOOutput.cc"
#include "TopTreeAnalysis/AnomCouplings/interface/BTagStudy.h"
#include "TopTreeAnalysis/AnomCouplings/src/BTagStudy.cc"       //--> Need to fix the makefile (don't include .cc files ...)
#include "TopTreeAnalysis/AnomCouplings/src/MlbStudy.cc"
#include "TopTreeAnalysis/AnomCouplings/interface/MlbStudy.h"
#include "TLorentzVector.h"

using namespace std;
using namespace reweight;
using namespace TopTree;

int main (int argc, char *argv[])
{
  string rootFileName = "AnomCouplings.root";
  
  clock_t start = clock();
  
  cout << "********************************************************" << endl;
  cout << " Beginning of the program for creating the AnomCoupl Trees ! " << endl;
  cout << "********************************************************" << endl;
  
  //SetStyle if needed
  //setTDRStyle(); 
  setMyStyle();

  /////////////////////
  // Configuration
  /////////////////////

  //xml file
  string xmlFileName ="config/myAnomCouplConfig.xml";

  if (argc > 3)
    xmlFileName = (string)argv[3];
  
  const char *xmlfile = xmlFileName.c_str();
  cout << "used config file: " << xmlfile << endl;

  //Configuration output format
  TTree *configTree = new TTree("configTree","configuration Tree");
  TClonesArray* tcdatasets = new TClonesArray("Dataset",1000);
  configTree->Branch("Datasets","TClonesArray",&tcdatasets);
  TClonesArray* tcAnaEnv = new TClonesArray("AnalysisEnvironment",1000);
  configTree->Branch("AnaEnv","TClonesArray",&tcAnaEnv);

  //////////////////////////
  // Verbosity for output //
  //////////////////////////
  int verbosity                 = 1;
  //0 muet
  //1 Main Info
  //2 mcParticlesMatchin Info
  //3 
  //4 Info for each event
  //5 Debug
 
  ///////////////////////////////
  //  Run specific parts only  //
  ///////////////////////////////
  bool GenLHCOOutput = false;
  bool RecoLHCOOutput = false;  

  bool FinalEventSelectionChoiceIsMade = false;
  int CorrectEventFound5Jets = 0, CorrectEventFound4Jets = 0;

  //Values needed for bTag study (select which of the 6 b-tag options is optimal!)
  const int NrConsideredBTagOptions = 1;   //Make sure this number is also the same in the bTagStudy class!!
  const int ChosenBTagOption = 3;  //2 T b-tags!!

  int CorrectEventMlbMqqb[NrConsideredBTagOptions];
  int WrongEventMlbMqqb[NrConsideredBTagOptions];
  for(int ii = 0; ii < NrConsideredBTagOptions; ii++){
     CorrectEventMlbMqqb[ii] = 0;
     WrongEventMlbMqqb[ii] = 0;
  }

  float BJetWP[6] = {0.244,0.679,0.679,0.898,0.898,0.898};
  float LightJetWP[6] = {0.244,0.679,0.244,0.898,0.679,0.244};
  
  std::string OptionName[6] = {"2 L b-tags,             ", //#0
    			       "2 M b-tags,             ", //#1
			       "2 M b-tags, light L-veto", //#2
			       "2 T b-tags,             ", //#3
			       "2 T b-tags, light M-veto", //#4
			       "2 T b-tags, light L-veto"};//#5

  //Values needed for comparison between 4-jet and 5-jet case (Hence having the choice between three light jets!)
  int NrConsideredJets = 5;  //Options are 4 or 5

  ////////////////////////////////////
  /// AnalysisEnvironment  
  ////////////////////////////////////

  AnalysisEnvironment anaEnv;
  cout<<"Loading environment ..."<<endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
  new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv);
  int verbose = anaEnv.Verbose;
  float oldLuminosity = anaEnv.Luminosity;	// in 1/pb
 
  cout << "analysis environment luminosity for rescaling "<< oldLuminosity << endl;
  
  /////////////////////
  // Load Datasets
  /////////////////////

  TTreeLoader treeLoader;
  cout << " - Load datasets ..." << endl;
  vector < Dataset* > datasets;
  vector < Dataset* > datasetsMu;
  vector < Dataset* > datasetsEl;
  vector < Dataset* > datasetsPlot;

  treeLoader.LoadDatasets (datasets, xmlfile);
  for(unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);
  
  float Luminosity = oldLuminosity;

  float LuminosityMu = oldLuminosity;
  float LuminosityEl = oldLuminosity;

  bool isSemiMu = false;
  bool isSemiE = false;

  bool foundMu = false;
  bool foundEl = false;

  for (unsigned int d = 0; d < datasets.size (); d++) {
    
    if(Luminosity > datasets[d]->EquivalentLumi() ) Luminosity = datasets[d]->EquivalentLumi();
    string dataSetName = datasets[d]->Name();
    
    if(dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) {
      LuminosityMu = datasets[d]->EquivalentLumi();
      foundMu=true;
    }  
    if(dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) {
      LuminosityEl = datasets[d]->EquivalentLumi();
      foundEl=true;  
    }  

    if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) { // when only using one channel
      Luminosity = datasets[d]->EquivalentLumi();
    }   

    if( dataSetName.find("QCD") == 0 ) datasets[d]->SetColor(kYellow);
    if( dataSetName.find("TT") == 0 ) datasets[d]->SetColor(kRed+1);
    if( dataSetName.find("TTbarJets_Other") == 0 ) datasets[d]->SetColor(kRed-7);
    if( dataSetName.find("WJets") == 0 )
    {
      datasets[d]->SetTitle("W#rightarrowl#nu");
      datasets[d]->SetColor(kGreen-3);
    }
    if( dataSetName.find("ZJets") == 0 )
    {
      datasets[d]->SetTitle("Z/#gamma*#rightarrowl^{+}l^{-}");
      datasets[d]->SetColor(kAzure-2);
    }
    if( dataSetName.find("ST") == 0 || dataSetName.find("SingleTop") ==0 )
      datasets[d]->SetColor(kMagenta);
    
  }
  
  if(!foundMu && !foundEl && Luminosity != oldLuminosity) cout << "changed analysis environment luminosity to "<< Luminosity << endl;
  else {
    if(LuminosityMu != oldLuminosity) cout << "Muon PD: changed analysis environment luminosity to "<< LuminosityMu << endl;
    if(LuminosityEl != oldLuminosity) cout << "Electron PD: changed analysis environment luminosity to "<< LuminosityEl << endl;
  }
  
  // make a datasets vector only for 
  if (foundMu) {
    for (unsigned int d = 0; d < datasets.size (); d++) {
      string dataSetName = datasets[d]->Name();
      if ( ! (dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) ) 
	datasetsMu.push_back(datasets[d]);
    }
  }
  
  if (foundEl) {
    for (unsigned int d = 0; d < datasets.size (); d++) {
      string dataSetName = datasets[d]->Name();
      if ( ! (dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) ) 
	datasetsEl.push_back(datasets[d]);
    }
  }
  
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");
  
  //Global variable
  //TRootEvent* event = 0;
  
  //nof selected events
  double NEvtsData = 0;
  Double_t *nEvents = new Double_t[datasets.size()];

  //Chi-Squared or KinFit:
  bool applyKinFit = false;

  ////////////////////////////////////
  /// Normal Plots (TH1F* and TH2F*)
  ////////////////////////////////////
  //All histograms can be defined as pre-programmed maps which makes definitions and looping easier
  map<string,TH1F*> histo1D;     
  map<string,TH2F*> histo2D;  
  
  // Histograms needed to calculate the sigma and Mc mass (from mean value) for W and top mass distribution
  //   --> Comment out after initializing most recent values ( also lines 1046 and 1356 )
  histo1D["WMass"]= new TH1F("WMass","WMass", 200,0,160);
  histo1D["TopMass"]= new TH1F("TopMass","TopMass", 200,0,350);
  histo1D["MlbMass"]= new TH1F("MlbMass","MlbMass",200,0,300);
  histo1D["MqqbMass"]= new TH1F("MqqbMass","MlbMass",400,0,500);
  
  ////////////////////////////////////
  /// MultiSamplePlot
  ////////////////////////////////////

  map<string,MultiSamplePlot*> MSPlot;

  ////////////////////////////////////
  /// Selection table
  ////////////////////////////////////
      
  vector<string> CutsSelecTableSemiMu;
  CutsSelecTableSemiMu.push_back(string("preselected"));
  CutsSelecTableSemiMu.push_back(string("trigged"));
  CutsSelecTableSemiMu.push_back(string("Good PV"));
  CutsSelecTableSemiMu.push_back(string("1 selected muon"));
  CutsSelecTableSemiMu.push_back(string("Veto 2nd muon"));
  CutsSelecTableSemiMu.push_back(string("Veto electron"));
  
  vector<string> CutsSelecTableSemiEl;
  CutsSelecTableSemiEl.push_back(string("preselected"));
  CutsSelecTableSemiEl.push_back(string("trigged"));
  CutsSelecTableSemiEl.push_back(string("Good PV"));
  CutsSelecTableSemiEl.push_back(string("1 selected electron"));
  CutsSelecTableSemiEl.push_back(string("Veto muon"));
  CutsSelecTableSemiEl.push_back(string("Veto 2nd electron from Z-decay"));
  CutsSelecTableSemiEl.push_back(string("Conversion veto"));
  
  char LabelNJets[100];
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-3);
  CutsSelecTableSemiMu.push_back(string(LabelNJets));
  CutsSelecTableSemiEl.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-2);
  CutsSelecTableSemiMu.push_back(string(LabelNJets));
  CutsSelecTableSemiEl.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-1);
  CutsSelecTableSemiMu.push_back(string(LabelNJets));
  CutsSelecTableSemiEl.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets);
  CutsSelecTableSemiMu.push_back(string(LabelNJets));
  CutsSelecTableSemiEl.push_back(string(LabelNJets));

  if (verbose > 0)
    cout << " - CutsSelectionTable instantiated ..." << endl;
  SelectionTable selecTableSemiMu(CutsSelecTableSemiMu, datasets);
  selecTableSemiMu.SetLuminosity(LuminosityMu);
  SelectionTable selecTableSemiEl(CutsSelecTableSemiEl, datasets);
  selecTableSemiEl.SetLuminosity(LuminosityEl);

  if (verbose > 0)
    cout << " - SelectionTable instantiated ..." << endl;

  ////////////////////////
  // PileUp Reweighting //
  ////////////////////////

  //cout << Luminosity << endl;
  
  LumiReWeighting LumiWeights, LumiWeightsUp, LumiWeightsDown;

  LumiWeights = LumiReWeighting("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_Summer12_S10.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2012Data53X_UpToRun208357/nominal.root", "pileup", "pileup");
  LumiWeightsUp = LumiReWeighting("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_Summer12_S10.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2012Data53X_UpToRun208357/sys_up.root", "pileup", "pileup");
  LumiWeightsDown = LumiReWeighting("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_Summer12_S10.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2012Data53X_UpToRun208357/sys_down.root", "pileup", "pileup");

  cout << " Initialized LumiReWeighting stuff" << endl;
  
  ////////////////////////////////////
  //	Loop on datasets
  ////////////////////////////////////

  if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
  for (unsigned int d = 0; d < datasets.size (); d++) {
    
    string previousFilename = "";
    int iFile = -1;
    string dataSetName = datasets[d]->Name();
    
    int nSelectedMu=0, nSelectedMuLCSV=0, nSelectedMuMCSV=0, nSelectedMuTCSV=0;
    int nSelectedEl=0, nSelectedElLCSV=0, nSelectedElMCSV=0, nSelectedElTCSV=0;
    int nLargeLCSVEvents = 0, nLargeMCSVEvents = 0, nLargeTCSVEvents = 0;
    
    if (verbose > 1)
      cout << "   Dataset " << d << ": " << datasets[d]->Name () << "/ title : " << datasets[d]->Title () << endl;
    if (verbose > 1)
      std::cout<<"      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
    
    //open files and load
    cout<<"LoadEvent"<<endl;
    treeLoader.LoadDataset (datasets[d], anaEnv);
    cout<<"LoadEvent"<<endl;
    
    /////////////////////////////////////
    /// Initialize JEC factors
    /////////////////////////////////////
   	    
    vector<JetCorrectorParameters> vCorrParam;

    /*JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer12_V3_MC_L3Absolute_AK5PFchs.txt");
    JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer12_V3_MC_L2Relative_AK5PFchs.txt");
    JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer12_V3_MC_L1FastJet_AK5PFchs.txt");
    
    //  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
    vCorrParam.push_back(*L1JetPar);
    vCorrParam.push_back(*L2JetPar);
    vCorrParam.push_back(*L3JetPar);

    if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) { // DATA!
      JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer12_V3_DATA_L2L3Residual_AK5PFchs.txt");
      vCorrParam.push_back(*ResJetCorPar);
      }*/
    
    JetCorrectionUncertainty *jecUnc =new JetCorrectionUncertainty(*(new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Fall12_V6_DATA_UncertaintySources_AK5PFchs.txt", "Total")));
    
    // true means redo also the L1
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);
    
    /////////////////////////////////////////
    //  LHCO Output files + GeneratorInfo  //
    /////////////////////////////////////////
    ofstream EventInfoFile;
    EventInfoFile.open("EventNumberInformation.lhco");
    EventInfoFile << " Event Number  MuPos  MuNeg  ElPos  ElNeg  ChannelNumber  selectedEvent  selectedChannelNumber " << endl;

    LHCOOutput lhcoOutput; //Initialize class
    ofstream outFile[4];
    outFile[0].open("TTbarLHCO_PositiveMuon.lhco");
    outFile[1].open("TTbarLHCO_NegativeMuon.lhco");
    outFile[2].open("TTbarLHCO_PositiveElectron.lhco");
    outFile[3].open("TTbarLHCO_NegativeElectron.lhco");

    ofstream outFileReco[16];
    unsigned int NumberPosRecoMu = 0, NumberNegRecoMu = 0, NumberPosRecoEl = 0, NumberNegRecoEl = 0;
    
    unsigned int NumberCorrectEvents = 0; //Counts the number of semi-leptonic events
    unsigned int NumberNegativeElectrons = 0, NumberNegativeMuons = 0, NumberPositiveElectrons = 0, NumberPositiveMuons = 0;
    int EventContent[5]; //0:top; 1:b; 2: u,c,d,s; 3:W; 4:mu + neutrino
    
    //Cos Theta information
    TLorentzVector *sTop, *WLeptTRF, *leptonWRF;
    float standardCosTheta = 0;
    TH1F h_StandardCosThetaNoEvtSel("StCosTheta_BeforeEvtSel","StCosTheta_BeforeEvtSel",200,-1,1);
    TH1F h_StandardCosThetaNoBTag("StCosThetaNoBTag","StCosThetaNoBTag",200,-1,1);
    TH1F h_StandardCosThetaLCSV("StCosThetaLCSV","StCosThetaLCSV",200,-1,1);
    TH1F h_StandardCosThetaAllLCSV("StCosThetaAllLCSV","StCosThetaAllLCSV",200,-1,1);
    TH1F h_StandardCosThetaMCSV("StCosThetaMCSV","StCosThetaMCSV",200,-1,1);
    TH1F h_StandardCosThetaTCSV("StCosThetaTCSV","StCosThetaTCSV",200,-1,1);
    TH1F h_JetTypeLargeLCSVEvents("JetTypeLargeLCSVEvents","JetTypeLargeLCSVEvents",51,-25.5,25.5);
    TH1F h_JetTypeLargeLCSVLeadingPtEvents("JetTypeLargeLCSVLeadingPtEvents","JetTypeLargeLCSVLeadingPtEvents",51,-25.5,25.5);
    TH1F h_JetTypeLCSVLightJetsLeadingPt("JetTypeLCSVLightJetsLeadingPt","JetTypeLCSVLightJetsLeadingPt",51,-25.5,25.5);
    TH1F h_JetTypeLargeMCSVEvents("JetTypeLargeMCSVEvents","JetTypeLargeMCSVEvents",51,-25.5,25.5);
    TH1F h_JetTypeLargeMCSVLeadingPtEvents("JetTypeLargeMCSVLeadingPtEvents","JetTypeLargeMCSVLeadingPtEvents",51,-25.5,25.5);
    TH1F h_JetTypeLargeTCSVEvents("JetTypeLargeTCSVEvents","JetTypeLargeTCSVEvents",51,-25.5,25.5);
    TH1F h_JetTypeLargeTCSVLeadingPtEvents("JetTypeLargeTCSVLeadingPtEvents","JetTypeLargeTCSVLeadingPtEvents",51,-25.5,25.5);
    TH1F h_JetTypeLCSV("JetTypeLCSV","JetTypeLCSV",51,-25.5,25.5);
    TH1F h_JetTypeMCSV("JetTypeMCSV","JetTypeMCSV",51,-25.5,25.5);
    TH1F h_JetTypeTCSV("JetTypeTCSV","JetTypeTCSV",51,-25.5,25.5);
    TH1F h_JetTypeLCSVLightJets("JetTypeLCSVLightJets","JetTypeLCSVLightJets",51,-25.5,25.5);
    TH1F h_CSVDiscrLCSVLightJets("CSVDiscrLCSVLightJets","CSVDiscrLCSVLightJets",400,-2.5,1.5);
    TH1F h_CSVDiscrLCSVLightJetsLeadingPt("CSVDiscrLCSVLightJetsLeadingPt","CSVDiscrLCSVLightJetsLeadingPt", 400, -2.5, 1.5);
    TH1F h_CorrectBLeptCSVDiscr("CorrectBLeptCSVDiscr","CorrectBLeptCSVDiscr",400,-2.5,1.5);
    TH1F h_CorrectBHadrCSVDiscr("CorrectBHadrCSVDiscr","CorrectBHadrCSVDiscr",400,-2.5,1.5);
    TH1F h_CorrectQuark1CSVDiscr("CorrectQuark1CSVDiscr","CorrectQuark1CSVDiscr",400,-2.5,1.5);
    TH1F h_CorrectQuark2CSVDiscr("CorrectQuark2CSVDiscr","CorrectQuark2CSVDiscr",400,-2.5,1.5);

    TH1F h_Quark1JetNumber("Quark1JetNumber","Quark1JetNumber",12,-1.5,10.5);
    TH1F h_Quark2JetNumber("Quark2JetNumber","Quark2JetNumber",12,-1.5,10.5);

    TH1F h_CosThetaReco("CosThetaReco","CosThetaReco",200,-1,1);
    TH1F h_NeutrinoEta("NeutrinoEta","NeutrinoEta",200,-8,8);

    //Mlb and Mqqb information:
    //
    TH1F h_WMass("WMass","WMass", 200,0,160);
    TH1F h_TopMass("TopMass","TopMass", 200,0,350);
    TH1F h_MlbMass("MlbMass","MlbMass",200,0,300);
    TH1F h_MqqbMass("MqqbMass","MlbMass",400,0,500);
    
    TH2F h_MlbMqqbCorrectChosen("MlbMqqbCorrectChosen","MlbMqqbCorrectChosen",200,0,500,200,0,300);
    TH2F h_MlbMqqbCorrectAll("MlbMqqbCorrectAll","MlbMqqbCorrectAll",200,0,500,200,0,300);
    TH2F h_MlbMqqbWrongOne("MlbMqqbWrongOne","MlbMqqbWrongOne",200,0,500,200,0,300);
    TH2F h_MlbMqqbWrongTwo("MlbMqqbWrongTwo","MlbMqqbWrongTwo",200,0,500,200,0,300);

    bool FalseEventContent = false;
    cout << " FalseEventContent : " << FalseEventContent << endl;
    TRootMCParticle *Top,*TopBar,*Bottom, *BottomBar,*Lepton,*NeutrinoMC,*WPlus,*WMinus,*Light,*LightBar;

    ////////////////////////////
    //  Class for bTag study  //
    ////////////////////////////  
    BTagStudy bTagStudy;  //--> Should only be called before the event loop (otherwise the counters will not give the correct result)
    bTagStudy.InitializeBegin();
    MlbStudy mlbStudy;
    mlbStudy.initializeBegin();

    ////////////////////////////////////
    //	loop on events
    ////////////////////////////////////
    
    nEvents[d] = 0;
    int itriggerSemiMu = -1,itriggerSemiEl = -1, previousRun = -1;
    if (verbose > 1)
      cout << "	Loop over events " << endl;

    //for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++){
    for (unsigned int ievt = 0; ievt < 5000; ievt++){

      //Initialize all values:
      bTagStudy.InitializePerEvent();
      mlbStudy.initializePerEvent();

      if(verbosity > 3) std::cout << " Looking at event : " << ievt << std::endl;    
      vector < TRootVertex* > vertex;
      vector < TRootMuon* > init_muons;
      vector < TRootElectron* > init_electrons;
      vector < TRootJet* > init_jets_corrected;
      vector < TRootJet* > init_jets;
      vector < TRootMET* > mets;
      vector < TRootGenJet* > genjets;
      
      nEvents[d]++;
      
      if(ievt%1000 == 0)
	std::cout<<"Processing the "<<ievt<<"th event (" << ((double)ievt/(double)datasets[d]->NofEvtsToRunOver())*100  << "%)" << " +> # selected: " << nSelectedMu << " (mu+jets) " << nSelectedEl << " (e+jets)" << flush<<"\r";
      
      ////////////////
      // LOAD EVENT //
      ////////////////
      
      TRootEvent* event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets_corrected, mets);  

      if(! (dataSetName.find("Data")==0 || dataSetName.find("DATA")==0  || dataSetName.find("data")==0 ) ) {
        genjets = treeLoader.LoadGenJet(ievt,false);
        sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
      
      vector<TRootMCParticle*> mcParticles;      
      if(dataSetName.find("TTbarJets") == 0){
	
        treeLoader.LoadMCEvent(ievt, 0, 0, mcParticles,false);  
        sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
      
      // check with genEvent which ttbar channel it is
      if(dataSetName.find("TTbarJets") == 0)  {
	//cout << "LOADING GenEvent" << endl;
	TRootGenEvent* genEvt = treeLoader.LoadGenEvent(ievt,false);
	//std::cout << "genEvt: " << genEvt << std::endl;
	if( genEvt->isSemiLeptonic(TRootGenEvent::kMuon) ) {
	  isSemiMu=true;
	  isSemiE=false;
	}
	else if( genEvt->isSemiLeptonic(TRootGenEvent::kElec) ) {
	  isSemiMu=false;
	  isSemiE=true;
	}
	else {
	  isSemiMu=false;
	  isSemiE=false;
	}
      }

      /////////////////////////////////
      // DETERMINE EVENT SCALEFACTOR //
      /////////////////////////////////
      
      // scale factor for the event
      float scaleFactor = 1.;
      
      // Load the GenEvent and calculate the branching ratio correction
      /*if(dataSetName.find("TTbarJets") == 0)
	{
	//cout << "LOADING GenEvent" << endl;
	  TRootGenEvent* genEvt = treeLoader.LoadGenEvent(ievt,false);
	  if( genEvt->isSemiLeptonic() )
	    scaleFactor *= (0.108*9.)*(0.676*1.5);
	  else if( genEvt->isFullHadronic() )
	    scaleFactor *= (0.676*1.5)*(0.676*1.5);
	  else if( genEvt->isFullLeptonic() )
	    scaleFactor *= (0.108*9.)*(0.108*9.);

      }*/
      
      //////////////////////////////////////
      // Apply Jet Corrections on-the-fly //   
      //////////////////////////////////////

      // not needed for now, GT contains good stuff
      /*if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) {
      	//jetTools->correctJets(init_jets_corrected,event->kt6PFJetsPF2PAT_rho(),true); //last boolean: isData (needed for L2L3Residual...)
      } else {
	jetTools->correctJets(init_jets_corrected,event->kt6PFJets_rho(),false); //last boolean: isData (needed for L2L3Residual...)
	}*/

      // PU reweighting

      // old method

      double lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() );

      if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0)
	lumiWeight=1;
      
      // up syst -> lumiWeight = LumiWeightsUp.ITweight( (int)event->nTruePU() );
      // down syst -> lumiWeight = LumiWeightsDown.ITweight( (int)event->nTruePU() );

      scaleFactor = scaleFactor*lumiWeight;

      ///////////////////
      // TRIGGER SETUP //
      ///////////////////

      string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
      if(previousFilename != currentFilename){
	previousFilename = currentFilename;
	iFile++;
	cout<<"File changed!!! => iFile = "<<iFile<<endl;
      }
      
      int currentRun = event->runId();
      
      if(previousRun != currentRun) {
        previousRun = currentRun;
	
	//semi-mu
	if(dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) {
	  
	  // 2.7/fb recalib 
	  if( event->runId() <= 190738 )
	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v11"), currentRun, iFile);

	  else if( event->runId() >= 191043 && event->runId() <= 193621 )
	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v12"), currentRun, iFile);
	  
	  else if( event->runId() >= 193834 && event->runId() <= 196531 )
	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);

	  else if( event->runId() >= 198049 && event->runId() <= 199608)
	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v14"), currentRun, iFile);

	  else if( event->runId() >= 199698 && event->runId() <= 208357)
	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v15"), currentRun, iFile);

	  else {
	    cout << "Unknown run for HLTpath selection: " << event->runId() << endl;
	    exit(1);
	  }
	  
	  if( itriggerSemiMu == 9999 ){	    
            cout << "itriggerSemiMu == 9999 for SemiMu HLTpath selection: " << event->runId() << endl;
            exit(-1);
          }
	  
	} else {
	  itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun); // Summer12 DR53X
	}
	
	// semi-el
	// semi-electron
        if(dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0 ) {
	  
	  // 2.7/fb recalib 
	  if( event->runId() <= 190738 )
	    itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v8"), currentRun, iFile);
	  
	  else if( event->runId() >= 191043 && event->runId() <= 191411 )
	    itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v9"), currentRun, iFile);

	  else if( event->runId() >= 191695 && event->runId() <= 196531)
	    itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v10"), currentRun, iFile);
	  
	  else if( event->runId() >= 198049 && event->runId() <= 208357)
	    itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v11"), currentRun, iFile);
  
	  else { 
            cout << "Unknown run for SemiEl HLTpath selection: " << event->runId() << endl;
	    exit(1);
	  }
	  if( itriggerSemiEl == 9999 ){	    
	    cout << "itriggerSemiEl == 9999 for SemiEl HLTpath selection: " << event->runId() << endl;
	    exit(-1);
	  }
        }
        else{	  
	  itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v10"), currentRun); // Summer12 DR53X 
	}
      }
      
      if (itriggerSemiMu == 9999 && itriggerSemiEl == 9999) {
	cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT IN RUN " << event->runId() << endl;
	exit(1);
      }
      
      ///////////////////////////////////////////////////////////////////////////////////////
      //   JES SYSTEMATICS && SMEAR JET RESOLUTION TO MIMIC THE RESOLUTION IN DATA         //
      //     - JES Corrections: The nominal corrections are already applied at PAT level   //
      //       so these tools should only be used for studies of the effect of systematics //
      ///////////////////////////////////////////////////////////////////////////////////////

      //JER Smearing:
      if( ! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) ) {
	
	/*    From file James (FourTop_EventSelection.cc)
	//JER
	doJERShift == 0;
	if(doJERShift == 1)
	jetTools->correctJetJER(init_jets, genjets, mets[0], "minus");
	else if(doJERShift == 2)
	jetTools->correctJetJER(init_jets, genjets, mets[0], "plus");
	else
	jetTools->correctJetJER(init_jets, genjets, mets[0], "nominal");
        
	//     coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JER correction:");
        
	// JES sysematic!
	if (doJESShift == 1)
	jetTools->correctJetJESUnc(init_jets, mets[0], "minus");
	else if (doJESShift == 2)
	jetTools->correctJetJESUnc(init_jets, mets[0], "plus");
        
	//            coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JES correction:");
	*/                                                         	                                                                                                	    
	
	jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "nominal",false);
	
	//jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "minus",false); //false means don't use old numbers but newer ones...
	//jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "plus",false);
	
	// Example how to apply JES systematics
	 
	//jetTools->correctJetJESUnc(init_jets_corrected, "minus",1);
	//jetTools->correctJetJESUnc(init_jets_corrected, "plus",1);
      }

      ////////////////////////////////////////////////////////
      // Access particle information before event selection //
      // Write this information to LHCO Output for MW       //
      ////////////////////////////////////////////////////////
      for(int ll = 0;ll<5;ll++){EventContent[ll]=0;}
      
      //Loop over all the mcParticles
      for(unsigned int i=0; i<mcParticles.size(); i++){
	if( mcParticles[i]->status() != 3) continue;
	
	int partType=mcParticles[i]->type(); if(verbosity>4)cout<<"-->Type of mcParticle : "<<partType<<endl;
	
	if(fabs(partType)<7 || fabs(partType)==24 || (fabs(partType)<=14 && fabs(partType)>=11) ){ //Considering only the semileptonic particles
	  int motherType=mcParticles[i]->motherType(); 
	  int grannyType=mcParticles[i]->grannyType();
	  if(verbosity > 5)cout<<"Mother type of particle : "<<motherType<<", and granny type : "<<grannyType<<endl;
	  
	  if(partType == 6){      
	    Top   =(TRootMCParticle*) mcParticles[i]; EventContent[0]++; if(verbosity>4) cout<<"*Particle found: Top"<<endl;
	  }
	  else if(partType == -6){
	    TopBar=(TRootMCParticle*) mcParticles[i]; EventContent[0]++; if(verbosity>4) cout<<"*Particle found: AntiTop"<<endl;
	  }
	  
	  else if(fabs(partType) == 5 && fabs(motherType) == 6){
	    EventContent[1]++;
	    if(partType == 5){      Bottom =    (TRootMCParticle*) mcParticles[i]; if(verbosity>4) cout<<"*Particle found: Bottom"<<endl;}
	    else if(partType == -5){BottomBar = (TRootMCParticle*) mcParticles[i]; if(verbosity>4) cout<<"*Particle found: AntiBottom"<<endl;}
	  }//End of bottom particle identification  
	  
	  else if(fabs(partType) == 24 && fabs(motherType) == 6){//Check correct definition!!!
	    EventContent[3]++;
	    if(partType == 24){      WPlus =  (TRootMCParticle*) mcParticles[i]; if(verbosity>4) cout<<"*Particle found: WPlus"<<endl;}        
	    else if(partType == -24){WMinus = (TRootMCParticle*) mcParticles[i]; if(verbosity>4) cout<<"*Particle found: WMinus"<<endl;}
	  }//End of WBoson identification
	  
	  else if(fabs(partType) <=4 && fabs(motherType) == 24 && fabs(grannyType) == 6){
	    EventContent[2]++;
	    if(partType > 0){     Light =    (TRootMCParticle*) mcParticles[i]; if(verbosity>4) cout<<"*Particle found: Light"<<endl;}
	    else if(partType < 0){LightBar = (TRootMCParticle*) mcParticles[i]; if(verbosity>4) cout<<"*Particle found: AntiLight"<<endl;}
	  }//End of light particle identification
	  else if((fabs(partType) == 13 || fabs(partType) == 11 ) && fabs(motherType) == 24 && fabs(grannyType) == 6){
	    EventContent[4]++;
	    string leptonType="";
	    if(fabs(partType) == 13){      if(verbosity>4) leptonType = "*Particle found: Muon";}
	    else if(fabs(partType) == 11){ if(verbosity>4) leptonType = "*Particle found: Electron";}
	    Lepton = (TRootMCParticle*) mcParticles[i]; if(verbosity > 4) cout<<leptonType<<endl;
	  }//End of lepton identification
	  
	  else if((fabs(partType) == 14 || fabs(partType) == 12 ) && fabs(motherType) == 24 && fabs(grannyType) == 6){
	    EventContent[4]++;
	    string neutrinoType="";
	    if(fabs(partType) == 14){      if(verbosity>4) neutrinoType = "*Particle found: Muon-neutrino";}
	    else if(fabs(partType) == 12){ if(verbosity>4) neutrinoType = "*Particle found: Electron-neutrino";}
	    NeutrinoMC = (TRootMCParticle*) mcParticles[i]; if(verbosity > 4) cout<<neutrinoType<<endl;
	  }//End of neutrino identification		
	  
	}//End of looking at semi-leptonic particles inside event ==> Semileptonic event is completely created now!	
      }//End of loop over mcParticles inside one particular event

      //////////////////////////////////////////////////////////////////////
      //  Consider only events with correct event content (b b q q l vl)  //
      //////////////////////////////////////////////////////////////////////
      EventInfoFile << "         " << ievt << "    ";
      if(ievt == 871811) cout << " Lepton type : " << Lepton->type() << " \n " << std::endl;
      if(EventContent[0]==2 && EventContent[1]==2 && EventContent[2]==2 && EventContent[3]==2 && EventContent[4]==2){
	FalseEventContent = false;
	vector<TRootMCParticle*> LHCOVector(6);
	vector<int> MadGraphId(6,4);
	
	NumberCorrectEvents++;
	if(verbosity>3){
	  cout << " Event with correct event content found " << endl;
	  cout << " Mass of bottom quark : " << Bottom->M() << endl;
	  cout << " Mass of light quark : " << Light->M() << endl;
	  cout << " Mass of LightBar quark : " << LightBar->M() << endl;
	  cout << " Mass of BottomBar quark : " << BottomBar->M() << endl;
	  cout << " Mass of lepton : " << Lepton->M() << endl;
	  cout << " Mass of neutrino : " << NeutrinoMC->M() << endl;
	}
	
	//Create the lhco file for pp > t t~:
	if(Lepton->type() == 13 || Lepton->type() == 11){ //Negative lepton, hence t~ > b~ W-, W- > e/mu- ve/vm
	  LHCOVector[0] = Bottom;
	  LHCOVector[1] = Light;
	  LHCOVector[2] = LightBar;
	  LHCOVector[3] = BottomBar;
	  LHCOVector[4] = Lepton;
	  LHCOVector[5] = NeutrinoMC;
	  if(Lepton->type() == 11){           //Looking at negative electron events (index 3 for LHCO file)
	    MadGraphId[4] = 1; //MadGraphId of e = 1
	    MadGraphId[5] = 6; 
	    NumberNegativeElectrons++;
	    if(GenLHCOOutput == true){
		lhcoOutput.LHCOEventOutput(3, outFile[3], NumberNegativeElectrons,LHCOVector,MadGraphId);
	    	EventInfoFile << "  0      0       0       1      " << NumberNegativeElectrons << "   ";
	    }
	  }//Negative electron
	  else if(Lepton->type() == 13){       //Looking at negative muon events (index 1 for LHCO file)
	    MadGraphId[4] = 2; //MadGraphId of mu = 2
	    MadGraphId[5] = 6; 
	    NumberNegativeMuons++;
	    if(GenLHCOOutput == true){
		lhcoOutput.LHCOEventOutput(1, outFile[1], NumberNegativeMuons,LHCOVector,MadGraphId);
		EventInfoFile << "  0      1       0       0      " << NumberNegativeMuons << "   ";
	    }
	  }//Negative muon

	  if(verbosity>3){
	    cout<<" WMinus information : "<<WMinus->Px()<< ", "<<WMinus->Py()<<", "<< WMinus->Pz()<<", "<<WMinus->E()<< endl;
	    cout<<" TopBar information : "<<TopBar->Px()<< ", "<<TopBar->Py()<<", "<< TopBar->Pz()<<", "<<TopBar->E()<< endl;
	  }
	  WLeptTRF = (TLorentzVector*) WMinus;
	  sTop = (TLorentzVector*) TopBar;				
	  if(verbosity>3){
	    cout<<" WLeptTRF information : "<<WLeptTRF->Px()<<", "<<WLeptTRF->Py()<<", "<<WLeptTRF->Pz()<<", "<<WLeptTRF->E()<<endl;
	    cout<<" sTop information : "<<sTop->Px()<<", "<<sTop->Py()<<", "<<sTop->Pz()<<", "<<sTop->E()<<endl;
	  }
	}//Negative lepton
	else if(Lepton->type() == -13 || Lepton->type() == -11){ //Positive lepton, hence t > b W+, W+ > e/mu+ ve/vm
	  LHCOVector[0] = Bottom; 
	  LHCOVector[1] = Lepton;
	  LHCOVector[2] = NeutrinoMC;
	  LHCOVector[3] = BottomBar;
	  LHCOVector[4] = Light;
	  LHCOVector[5] = LightBar;
	  if(Lepton->type() == -11){            //Looking at positive electron events (index 2 for LHCO file)
	    MadGraphId[1] = 1; //MadGraphId of electron = 1
	    MadGraphId[2] = 6; 
	    NumberPositiveElectrons++;
	    if(GenLHCOOutput == true){
		lhcoOutput.LHCOEventOutput(2, outFile[2], NumberPositiveElectrons,LHCOVector,MadGraphId);
		EventInfoFile << "  0      0       1       0      " << NumberPositiveElectrons << "   ";
	    }
	  }//Positive electron
	  else if(Lepton->type() == -13){             //Looking at positive muon events (index 0 for LHCO file)
	    MadGraphId[1] = 2; //MadGraphId of muon = 2
	    MadGraphId[2] = 6; 
	    NumberPositiveMuons++;
	    if(GenLHCOOutput == true){
		lhcoOutput.LHCOEventOutput(0, outFile[0], NumberPositiveMuons,LHCOVector,MadGraphId);
		EventInfoFile << "  1      0       0       0      " << NumberPositiveMuons << "   ";
	    }
	  }//Positive muon
	  
	  if(verbosity>3){
	    cout << " WPlus information : "<<WPlus->Px()<< ", "<<WPlus->Py()<<", "<< WPlus->Pz()<<", "<<WPlus->E()<< endl;
	    cout << " Top information : "<<Top->Px()<< ", "<<Top->Py()<<", "<< Top->Pz()<<", "<<Top->E()<< endl;
	  }
	  WLeptTRF = (TLorentzVector*) WPlus;
	  sTop = (TLorentzVector*) Top;			
	}//Positive lepton
	
	//////////////////////////////////////
	//  Look at cos theta distribution  //
	////////////////////////////////////// 
	
	//-----    Applying boost on muon and W    -----//
	leptonWRF = Lepton;
	leptonWRF->Boost(-WLeptTRF->BoostVector());
	WLeptTRF->Boost(-sTop->BoostVector());
	if(verbosity>3){
	  cout<<" leptonWRF information : "<<leptonWRF->Px()<<", "<<leptonWRF->Py()<<", "<<leptonWRF->Pz()<<", "<<leptonWRF->E()<<endl;
	}
	
	//-----   Calculating cos theta:   -----
	standardCosTheta = ((WLeptTRF->Vect()).Dot(leptonWRF->Vect()))/(((WLeptTRF->Vect()).Mag())*((leptonWRF->Vect()).Mag()));
	if(verbosity>4) cout << " cos theta (gen): " << standardCosTheta << endl << endl;
	h_StandardCosThetaNoEvtSel.Fill(standardCosTheta);
	
      }//Correct event content found
      else{
	FalseEventContent = true;
	EventInfoFile << endl;
	if(verbosity>4){
	  cout << " Number of top quarks      : " << EventContent[0] << endl;
	  cout << " Number of bottom quarks   : " << EventContent[1] << endl;
	  cout << " Number of light quarks    : " << EventContent[2] << endl;
	  cout << " Number of W-bosons        : " << EventContent[3] << endl;
	  cout << " Number of lepton/neutrino : " << EventContent[4] << endl;
	}
      }			    
 
      /////////////////////
      // EVENT SELECTION //
      /////////////////////
      
      //Declare selection instance    
      Selection selection(init_jets_corrected, init_muons, init_electrons, mets, event->kt6PFJets_rho());
      selection.setJetCuts(40,2.5,0.01,1.,0.98,0.3,0.1);     
      selection.setMuonCuts(25,2.1,0.12,0.2,0.3,1,0.5,5,0); 
      selection.setElectronCuts(32,2.5,0.1,0.02,0.5,0.3,0); 
      selection.setLooseMuonCuts(10,2.5,0.2);
      selection.setLooseElectronCuts(20,2.5,0.15,0.);
      
      bool triggedSemiMu = false;
      bool triggedSemiEl = false;

      if( ! (dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) )
        triggedSemiMu = treeLoader.EventTrigged (itriggerSemiMu);
      if( ! (dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) )
        triggedSemiEl = treeLoader.EventTrigged (itriggerSemiEl);

      bool isGoodPV = selection.isPVSelected(vertex, 4, 24, 2.);

      vector<TRootJet*> selectedJets, selectedJetsNoMu, selectedJetsNoEl;
      vector<TRootMuon*> selectedMuons;
      vector<TRootElectron*> selectedElectrons;
      vector<TRootMuon*> vetoMuons = selection.GetSelectedLooseMuons(10,2.5,0.2);
      vector<TRootElectron*> vetoElectronsSemiMu = selection.GetSelectedLooseElectrons(20,2.5,0.15);
      vector<TRootElectron*> vetoElectronsSemiEl = selection.GetSelectedLooseElectrons(20,2.5,0.15);

      if( dataSetName.find("InvIso") != string::npos )  { // event selection for special Data TopTrees for ARI QCD
        vector<TRootMuon*> overlapMuons = selection.GetSelectedMuonsInvIso(0.2, vertex[0]);
        vector<TRootElectron*> overlapElectrons = selection.GetSelectedElectronsInvIso(0.2);
        selectedJetsNoMu = selection.GetSelectedJets(overlapMuons,true);
        selectedJetsNoEl = selection.GetSelectedJets(overlapElectrons,true);

	/*if (selectedJetsNoMu.size() >= 4) {
	  //cout << "ol" << endl;
	  if (selectedJetsNoMu[0]->Pt() < 45) selectedJetsNoMu.clear();
	  if (selectedJetsNoMu[1]->Pt() < 45) selectedJetsNoMu.clear();
	  if (selectedJetsNoMu[2]->Pt() < 40) selectedJetsNoMu.clear();
	  if (selectedJetsNoMu[3]->Pt() < 40) selectedJetsNoMu.clear();
	}

	if (selectedJetsNoEl.size() >= 4) {
	  //cout << "ol" << endl;
	  if (selectedJetsNoEl[0]->Pt() < 45) selectedJetsNoEl.clear();
	  if (selectedJetsNoEl[1]->Pt() < 45) selectedJetsNoEl.clear();
	  if (selectedJetsNoEl[2]->Pt() < 40) selectedJetsNoEl.clear();
	  if (selectedJetsNoEl[3]->Pt() < 40) selectedJetsNoEl.clear();
	  }*/

	//selectedJetsNoMu = selection.GetSelectedJets(true);
        //selectedJetsNoEl = selection.GetSelectedJets(true);

	selectedMuons = selection.GetSelectedMuonsInvIso(0.2, vertex[0], selectedJetsNoMu);
        selectedElectrons = selection.GetSelectedElectronsInvIso(0.2,selectedJetsNoEl);
      }
      else { // Normal selection
	selectedJets = selection.GetSelectedJets(true);
	
	/*if (selectedJets.size() >= 4) {
	//cout << "ol" << endl;
	  if (selectedJets[0]->Pt() < 45) selectedJets.clear();
	  if (selectedJets[1]->Pt() < 45) selectedJets.clear();
	  if (selectedJets[2]->Pt() < 40) selectedJets.clear();
	  if (selectedJets[3]->Pt() < 40) selectedJets.clear();
	  }*/

	//selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);	
	selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);
	selectedElectrons = selection.GetSelectedElectrons(selectedJets);
      }

      //Variables containing the right jet combination for TTJets
      int CorrectQuark1=999, CorrectQuark2=999, CorrectBHadronic=999, CorrectBLeptonic=999;
      
      ///////////////////////////////////////
      //  Initialize variables ChiSquared  //
      //  Look for correct combination     //
      //  --> Only used for b-quarks       //
      ///////////////////////////////////////
      float ChiSquared[2];  //Needed for chi squared caclulation      
      int UsedCombination;
      int BHadronicIndex[2];
      int BLeptIndex, BHadrIndex, QOneIndex, QTwoIndex;
      float ChiSquaredValue;
      	    
      //float MassW=83.6103;
      //float MassTop = 172.956;
      //float SigmaW=11.1534;  //Obtained from gaussian fit on Top and W distribution with simulated information
      //float SigmaTop=18.232;

      float MassMlb = 103.2552;
      float MassMqqb = 178.7309;
      float SigmaMlb = 27.0706;
      float SigmaMqqb = 18.2675;
      
      /////////////////////////////
      // Neutrino Reconstruction //
      /////////////////////////////
      float NeutrinoPx;
      float NeutrinoPy;
      float NeutrinoPz=99.;//with this value it can be distinguished in plot!
      TLorentzVector Neutrino;
      
      //////////////////////
      // Event selection  //
      //////////////////////
      bool eventselectedSemiMu = false;
      bool eventselectedSemiEl = false;
      
      // semi-mu selection
      if (triggedSemiMu) {
	if (isGoodPV) {
	  if(dataSetName.find("InvIso") != string::npos) selectedJets = selectedJetsNoMu;
	  if (selectedMuons.size() == 1) {
	    if( vetoMuons.size() == 1 || ( dataSetName.find("InvIso") != string::npos && vetoMuons.size() == 0 ) ) { // if InvertedIso, selected muon not part of vetoMuons vector!
	      if (vetoElectronsSemiMu.size() == 0) {
		if (selectedJets.size() >= 4) {
		  eventselectedSemiMu = true;
		}
	      }
	    }
	  }
	}
      }
      selecTableSemiEl.Fill(d,0,scaleFactor*lumiWeight);
      
      if( triggedSemiEl) {
	if(dataSetName.find("InvIso") != string::npos) selectedJets = selectedJetsNoEl;
	selecTableSemiEl.Fill(d,1,scaleFactor*lumiWeight);
	if (isGoodPV ) {
	  selecTableSemiEl.Fill(d,2,scaleFactor*lumiWeight);
	  if( selectedElectrons.size() == 1 ) {
	    selecTableSemiEl.Fill(d,3,scaleFactor*lumiWeight);
	    if( vetoMuons.size() == 0 ) {
	      selecTableSemiEl.Fill(d,4,scaleFactor*lumiWeight);
	      if (vetoElectronsSemiEl.size() == 1) {
		selecTableSemiEl.Fill(d,5,scaleFactor*lumiWeight);
		selecTableSemiEl.Fill(d,6,scaleFactor*lumiWeight);
		if( selectedJets.size()>=1 ) {
		  selecTableSemiEl.Fill(d,7,scaleFactor*lumiWeight);
		  if( selectedJets.size()>=2 ) {
		    selecTableSemiEl.Fill(d,8,scaleFactor*lumiWeight);
		    if( selectedJets.size()>=3 ) {
		      selecTableSemiEl.Fill(d,9,scaleFactor*lumiWeight);
		      if( selectedJets.size()>=4 ) {
			selecTableSemiEl.Fill(d,10,scaleFactor*lumiWeight);
			eventselectedSemiEl=true;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
      
      if( !eventselectedSemiMu && !eventselectedSemiEl && !FalseEventContent) EventInfoFile << endl;
      if (!eventselectedSemiMu && !eventselectedSemiEl) continue;
      if(FalseEventContent == 0) EventInfoFile << "             1          ";  //To avoid tau's which are reconstructed as muons!
      
      //Counting the number of events passing through the 'basic' event selection requirements    
      if (eventselectedSemiMu) nSelectedMu++;
      if (eventselectedSemiEl) nSelectedEl++;
      
      //-----------------//
      // do some data-mc //
      //-----------------//
      
      // when running both electron and muon data, pick the right dataset vector and lumi for the MSPlots
      if (!foundMu && !foundEl) datasetsPlot = datasets;
      else if (eventselectedSemiMu) {
	datasetsPlot = datasetsMu;
	Luminosity = LuminosityMu;
      }
      else if (eventselectedSemiEl) {
	datasetsPlot = datasetsEl;
	Luminosity = LuminosityEl;
      }
      
      // Selecting correct lepton
      TLorentzVector* selectedLepton;
      float LeptonRecoCharge;
      if (eventselectedSemiMu){
	selectedLepton = (TLorentzVector*)selectedMuons[0];
	LeptonRecoCharge = selectedMuons[0]->charge();
      }
      else if (eventselectedSemiEl){
	selectedLepton = (TLorentzVector*)selectedElectrons[0];
	LeptonRecoCharge = selectedElectrons[0]->charge();
      }

      string leptonFlav="_other";
      
      if (eventselectedSemiMu) leptonFlav="_mu";
      else if (eventselectedSemiEl) leptonFlav="_el";
      
      // MSPlots after 'basic' event selection (no b-tag)
      if (MSPlot.find("Selected_Events_pT_jet1"+leptonFlav) == MSPlot.end()){
	MSPlot["Selected_Events_pT_jet1"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_jet1"+leptonFlav, 30, 0, 600, "p_{T} (GeV)");
	MSPlot["Selected_Events_pT_jet2"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_jet2"+leptonFlav, 30, 0, 600, "p_{T} (GeV)");
	MSPlot["Selected_Events_pT_jet3"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_jet3"+leptonFlav, 30, 0, 600, "p_{T} (GeV)");
	MSPlot["Selected_Events_pT_jet4"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_jet4"+leptonFlav, 30, 0, 600, "p_{T} (GeV)");
	MSPlot["Selected_Events_pT_4leadingjets"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_4leadingjets"+leptonFlav, 30, 0, 600, "p_{T} (GeV)");
	MSPlot["Selected_Events_pT_alljets"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_alljets"+leptonFlav, 30, 0, 600, "p_{T} (GeV)");
      }
      
      MSPlot["Selected_Events_pT_jet1"+leptonFlav]->Fill(selectedJets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["Selected_Events_pT_jet2"+leptonFlav]->Fill(selectedJets[1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["Selected_Events_pT_jet3"+leptonFlav]->Fill(selectedJets[2]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["Selected_Events_pT_jet4"+leptonFlav]->Fill(selectedJets[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      
      for (unsigned int q=0; q<selectedJets.size(); q++) {
	MSPlot["Selected_Events_pT_alljets"+leptonFlav]->Fill(selectedJets[q]->Pt(), datasets[d], true, Luminosity*scaleFactor);       
	if (q<4)
	  MSPlot["Selected_Events_pT_4leadingjets"+leptonFlav]->Fill(selectedJets[q]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      }
      
      ////////////////////////////////////////////////////////////////////
      //   Use genEvent information to get the correct event topology   //
      ////////////////////////////////////////////////////////////////////    
      float CorrectRecMassW=0;      //? Still needed if no chi-squared is done?
      float CorrectRecMassTop=0;      
      float CorrectRecMassMlb=0;
      float CorrectRecMassMqqb=0;
      vector<int> jetCombi;
      vector<TRootMCParticle> mcParticlesMatching;      	
      vector< pair<unsigned int, unsigned int> > JetPartonPair, ISRJetPartonPair; // First one is jet number, second one is mcParticle number

      if(dataSetName.find("TTbarJets") == 0){
	
	pair<unsigned int, unsigned int> leptonicBJet_, hadronicBJet_, hadronicWJet1_, hadronicWJet2_; //First index is the JET number, second one is the parton
	leptonicBJet_ = hadronicBJet_ = hadronicWJet1_ = hadronicWJet2_ = pair<unsigned int, unsigned int>(9999,9999);
	vector<TLorentzVector> mcParticlesTLV, selectedJetsTLV;
	//vector<TRootMCParticle> mcParticlesMatching;
	bool muPlusFromTop = false, muMinusFromTop = false;
	bool elPlusFromTop = false, elMinusFromTop = false;
	if(verbosity>1) cout << " Looking at mcParticlesMatching " << endl;
      	for(unsigned int i=0; i<mcParticles.size(); i++){
      	  if( mcParticles[i]->status() != 3) continue;
	  
	  //Muon identification:
      	  if( mcParticles[i]->type() == 13 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -6 ){
      	    if(muMinusFromTop) cerr<<"muMinusFromTop was already true"<<endl;
      	    muMinusFromTop = true;
      	  }
      	  if( mcParticles[i]->type() == -13 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == 6 ){
      	    if(muPlusFromTop) cerr<<"muPlusFromTop was already true"<<endl;
      	    muPlusFromTop = true;
      	  }
	  
	  //Electron identification:
     	  if( mcParticles[i]->type() == 11 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -6 ){
      	    if(elMinusFromTop) cerr<<"elMinusFromTop was already true"<<endl;
      	    elMinusFromTop = true;
      	  }
      	  if( mcParticles[i]->type() == -11 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == 6 ){
      	    if(elPlusFromTop) cerr<<"elPlusFromTop was already true"<<endl;
      	    elPlusFromTop = true;
      	  }
	  
      	  if( abs(mcParticles[i]->type()) < 6 || abs(mcParticles[i]->type()) == 21 ){
      	    mcParticlesTLV.push_back(*mcParticles[i]);
      	    mcParticlesMatching.push_back(*mcParticles[i]);
      	  }
      	}
      	if(muPlusFromTop && muMinusFromTop)
      	  cerr<<"muPlusFromTop and muMinusFromTop are both true ?!\nCheck if you are using the right sample..."<<endl;
      	if(elPlusFromTop && elMinusFromTop)
      	  cerr<<"elPlusFromTop and elMinusFromTop are both true ?!\nCheck if you are using the right sample..."<<endl;
	
	if(muPlusFromTop && verbosity>1) cout << " found  muPlus from Top " << endl;
	else if(muMinusFromTop && verbosity>1) cout << " found muMinus from Top " << endl;
	
      	// take all the selectedJets_ to study the radiation stuff, selectedJets are already ordened in decreasing Pt()
      	for(unsigned int i=0; i<selectedJets.size(); i++)
      	  selectedJetsTLV.push_back(*selectedJets[i]);
	
      	JetPartonMatching matching = JetPartonMatching(mcParticlesTLV, selectedJetsTLV, 2, true, true, 0.3);
	
      	if(matching.getNumberOfAvailableCombinations() != 1)
      	  cerr << "matching.getNumberOfAvailableCombinations() = "<<matching.getNumberOfAvailableCombinations()<<"  This should be equal to 1 !!!"<<endl;
	
      	//vector< pair<unsigned int, unsigned int> > JetPartonPair, ISRJetPartonPair; // First one is jet number, second one is mcParticle number
	
      	for(unsigned int i=0; i<mcParticlesTLV.size(); i++){
      	  int matchedJetNumber = matching.getMatchForParton(i, 0);
      	  if(matchedJetNumber != -1)
      	    JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
      	}

      	for(unsigned int i=0; i<JetPartonPair.size(); i++){
      	  unsigned int j = JetPartonPair[i].second;
	 
	  if(verbosity > 3){ 
		std::cout <<" Jet number " << JetPartonPair[i].first << " can be matched with the mcParticle number " << JetPartonPair[i].second << " which is of the type : " << mcParticlesMatching[j].type() << std::endl;
		std::cout << " Type of mcParticlesMatching[JetPartonPair[i].first] = " << mcParticlesMatching[JetPartonPair[i].first].type() << std::endl;
  }
 
      	  if( fabs(mcParticlesMatching[j].type()) < 6 ){
      	    if( ( (muPlusFromTop || elPlusFromTop) && mcParticlesMatching[j].motherType() == -24 && mcParticlesMatching[j].grannyType() == -6 ) || ( (muMinusFromTop || elMinusFromTop) && mcParticlesMatching[j].motherType() == 24 && mcParticlesMatching[j].grannyType() == 6 ) ){
      	      if(hadronicWJet1_.first == 9999) 
      		hadronicWJet1_ = JetPartonPair[i];
      	      else if(hadronicWJet2_.first == 9999) 
      		hadronicWJet2_ = JetPartonPair[i];
      	      else cerr<<"Found a third jet coming from a W boson which comes from a top quark..."<<endl;
      	    }
      	  }
      	  if( fabs(mcParticlesMatching[j].type()) == 5 ){
      	    if( ( (muPlusFromTop || elPlusFromTop) && mcParticlesMatching[j].motherType() == -6 ) || ( (muMinusFromTop || elMinusFromTop) && mcParticlesMatching[j].motherType() == 6 ) ){
      	      hadronicBJet_ = JetPartonPair[i];
	      if(verbosity > 3) std::cout << " Hadronic b jet is matched to JetPartonPair : " << JetPartonPair[i].first << " , " << JetPartonPair[i].second << std::endl;
	    }
      	    else if( ( (muPlusFromTop || elPlusFromTop) && mcParticlesMatching[j].motherType() == 6 ) || ( (muMinusFromTop || elMinusFromTop) && mcParticlesMatching[j].motherType() == -6 ) ){
      	      leptonicBJet_ = JetPartonPair[i];
	      if( verbosity > 3) std::cout << " Leptonic b jet is matched to JetPartonPair : " << JetPartonPair[i].first << " , " << JetPartonPair[i].second << std::endl;
	    }
      	  }
      	}

	//Plot the jet number of the two light jets:
	if(hadronicWJet1_.first != 9999) h_Quark1JetNumber.Fill(hadronicWJet1_.first);
	else 				 h_Quark1JetNumber.Fill(-1);
	if(hadronicWJet2_.first != 9999) h_Quark2JetNumber.Fill(hadronicWJet2_.first);
	else				 h_Quark2JetNumber.Fill(-1);

	if(verbose > 3){
	//First index is the jet number, the second one the mcParticle number:
	std::cout << " Indices of quark 1 : " << hadronicWJet1_.first << " & " << hadronicWJet1_.second << std::endl;
	if(hadronicWJet1_.first != 9999){
	   std::cout << "      --- Pt of jet : " << selectedJets[hadronicWJet1_.first]->Pt() << std::endl;
	   std::cout << "      --- CSV value : " << selectedJets[hadronicWJet1_.first]->btag_combinedSecondaryVertexBJetTags() << " (light cut = 0.244) " << std::endl;
	}
	if(hadronicWJet1_.second != 9999){
	   std::cout << "      --- Pt of mc  : " << mcParticlesMatching[hadronicWJet1_.second].Pt() << std::endl;
	   std::cout << "      --- Type      : " << mcParticlesMatching[hadronicWJet1_.second].type() << std::endl;
	}
	std::cout << " " << std::endl;
        std::cout << " Indices of quark 2 : " << hadronicWJet2_.first << " & " << hadronicWJet2_.second << std::endl;
        if(hadronicWJet2_.first != 9999){
	   std::cout << "      --- Pt of jet : " << selectedJets[hadronicWJet2_.first]->Pt() << std::endl;
	   std::cout << "      --- CSV value : " << selectedJets[hadronicWJet2_.first]->btag_combinedSecondaryVertexBJetTags() << " (light cut = 0.244) " << std::endl;
	}
        if(hadronicWJet2_.second != 9999){
           std::cout << "      --- Pt of mc  : " << mcParticlesMatching[hadronicWJet2_.second].Pt() << std::endl;
           std::cout << "      --- Type      : " << mcParticlesMatching[hadronicWJet2_.second].type() << std::endl;
        }
	std::cout << " ************************************" << std::endl;	
	}

      	jetCombi.push_back(hadronicWJet1_.first);
      	jetCombi.push_back(hadronicWJet2_.first);
      	jetCombi.push_back(hadronicBJet_.first);
      	jetCombi.push_back(leptonicBJet_.first);
	
      	CorrectQuark1=jetCombi[0];
      	CorrectQuark2=jetCombi[1];
      	CorrectBHadronic = jetCombi[2];
      	CorrectBLeptonic = jetCombi[3];

	for(int ii = 0; ii < JetPartonPair.size(); ii++){
	   if(CorrectBLeptonic == JetPartonPair[ii].first){
		if(verbosity>3) cout<<"Type of lept b : "<<mcParticlesMatching[JetPartonPair[ii].second].type()<<" - from JetPartonPair "<<JetPartonPair[ii].first<<" , "<<JetPartonPair[ii].second<<endl;
  	   }
	   else if(CorrectBHadronic == JetPartonPair[ii].first){
		if(verbosity>3) cout<<"Type of hadr b : "<<mcParticlesMatching[JetPartonPair[ii].second].type()<<" - from JetPartonPair "<<JetPartonPair[ii].first<<" , "<<JetPartonPair[ii].second<<endl;
	   }
	   else if(CorrectBHadronic == 9999){
		if(verbosity>3) std::cout << " CorrectBHadronic not matched to jet " << std::endl;
	   }
	   else if(CorrectBLeptonic == 9999){
		if(verbosity>3) std::cout << " CorrectBLeptonic not matched to jet " << std::endl;
	   }
	}

	if(CorrectBLeptonic != 9999) h_CorrectBLeptCSVDiscr.Fill(selectedJets[CorrectBLeptonic]->btag_combinedSecondaryVertexBJetTags());
	else h_CorrectBLeptCSVDiscr.Fill(-2);
        if(CorrectBHadronic != 9999) h_CorrectBHadrCSVDiscr.Fill(selectedJets[CorrectBHadronic]->btag_combinedSecondaryVertexBJetTags());
	else h_CorrectBHadrCSVDiscr.Fill(-2);
	if(CorrectQuark1 != 9999) h_CorrectQuark1CSVDiscr.Fill(selectedJets[CorrectQuark1]->btag_combinedSecondaryVertexBJetTags());
	else h_CorrectQuark1CSVDiscr.Fill(-2);
	if(CorrectQuark2 != 9999) h_CorrectQuark2CSVDiscr.Fill(selectedJets[CorrectQuark2]->btag_combinedSecondaryVertexBJetTags());
	else h_CorrectQuark2CSVDiscr.Fill(-2);
 	
	//Working on generator level (i.e. jets level):  
	if(jetCombi[0]!=9999 && jetCombi[1]!=9999 && jetCombi[2]!=9999 && jetCombi[3]!=9999){    
	  CorrectRecMassW=(*selectedJets[jetCombi[0]]+*selectedJets[jetCombi[1]]).M();
	  CorrectRecMassTop=(*selectedJets[jetCombi[0]]+*selectedJets[jetCombi[1]]+*selectedJets[jetCombi[2]]).M();
	  CorrectRecMassMlb=(*selectedJets[CorrectBLeptonic]+*selectedLepton).M();
	  CorrectRecMassMqqb=(*selectedJets[CorrectQuark1]+*selectedJets[CorrectQuark2]+*selectedJets[CorrectBHadronic]).M();

	  histo1D["WMass"]->Fill(CorrectRecMassW);
	  histo1D["TopMass"]->Fill(CorrectRecMassTop);
	  histo1D["MlbMass"]->Fill(CorrectRecMassMlb);
	  histo1D["MqqbMass"]->Fill(CorrectRecMassMqqb);
	  h_WMass.Fill(CorrectRecMassW);
	  h_TopMass.Fill(CorrectRecMassTop);
	  h_MlbMass.Fill(CorrectRecMassMlb);
	  h_MqqbMass.Fill(CorrectRecMassMqqb);
	}	      	      	      	       	      
      }//if dataset Semi mu ttbar

      //----------------------------------------------------------------------------------------------------------------------------------- Start of bTagStudy class stuff!!

      /////////////////////////////
      //   B-Tag requirements    //
      /////////////////////////////

      if(NrConsideredBTagOptions > 1){   //Go through the different options and compare the efficiencies!

         for(int option = 0; option < NrConsideredBTagOptions; option++){
	    bTagStudy.CalculateJets(selectedJets, BJetWP[option], LightJetWP[option], option);  //First float is the working point for the b-jets, the second is the one for the light jets!! 
            //Added integer points to which of the considered options is currently active. This should allow to obtain the final numbers at the end of the file!

            if(verbose > 3){
	       cout<<" Size of bTaggedJets: "<<(bTagStudy.getbTaggedJets(option)).size()<<" , of NonbTaggedJets: "<<(bTagStudy.getNonbTaggedJets(option)).size()<<" & of lightJets: "<<(bTagStudy.getLightJets(option)).size()<<endl;
	       cout<<" Index of BHadronic: "<<CorrectBHadronic<<" , Index of BLeptonic: "<<CorrectBLeptonic<<" , Index of quark1: "<<CorrectQuark1<<" & Index of quark2: "<<CorrectQuark2<<endl;
            }

            //Also check the correct jet combi:
            if((bTagStudy.getbTaggedJets(option)).size() >= 2 && (bTagStudy.getLightJets(option)).size() >=2){
	       bTagStudy.CorrectJetCombi(CorrectBHadronic, CorrectBLeptonic, CorrectQuark1, CorrectQuark2, option);
            }
            else{
	       if(verbose > 3) std::cout << " Event doesn't have two b-tagged jets and/or two light jets ! " << std::endl;
            }
         }

         ///////////////////////////////////////////////////////
         //  Kinematic information for different b-tag cases  //
         ///////////////////////////////////////////////////////
         //
         if(dataSetName.find("TTbarJets") == 0){

	   for(int jj = 0; jj < 2; jj++){      //Only look at the two leading jets!


	      //Kinematic information for the light (= Non LCSV jets)
	      if((bTagStudy.getLightJets(0)).size()>=2){
	        h_CSVDiscrLCSVLightJets.Fill(selectedJets[(bTagStudy.getLightJets(0))[jj]]->btag_combinedSecondaryVertexBJetTags());   
		//--> check whether the so-called light jets don't all have discr -1 ...	

	        for(int ii = 0; ii<JetPartonPair.size(); ii++){ //Look at all the matched jets
	          if((bTagStudy.getLightJets(0))[jj] == JetPartonPair[ii].first){  //Check whether the considered jet can be matched!
		    h_JetTypeLCSVLightJets.Fill(mcParticlesMatching[JetPartonPair[ii].second].type());
	          }
	          else{    //Unmatched jets!
		    h_JetTypeLCSVLightJets.Fill(25.);
	          }
	        }
	      }//End of light (= non LCSV) jets

	      //********** Try to combine this with a loop !!
	      // ==> Use vector of TH1F's!
	      // + Add the event selection counter for each option!
/*	      TH1F CosThetaHistos[3] = {h_StandardCosThetaLCSV, h_StandardCosThetaMCSV, h_StandardCosThetaTCSV};      --> Empty histograms for the moment ! 
	      TH1F JetTypeHistos[3]  = {h_JetTypeLCSV         , h_JetTypeMCSV        ,  h_JetTypeTCSV         };
	      int  bTagOption[3]     = {0                     , 1                    ,  3                     };
	      int nSelectedEvtsMu[3] = {nSelectedMuLCSV       , nSelectedMuMCSV      , nSelectedMuTCSV        };
              int nSelectedEvtsEl[3] = {nSelectedElLCSV       , nSelectedElMCSV      , nSelectedElTCSV        };

	      for(int options = 0; options < 3; options++){
		if( (bTagStudy.getbTaggedJets(bTagOption[options])).size() >= 2){
		   CosThetaHistos[options].Fill(standardCosTheta);
		   if(eventselectedSemiMu) nSelectedEvtsMu[options]++;    //Maybe not possible to have existing integers as elements!
		   if(eventselectedSemiEl) nSelectedEvtsEl[options]++;

		   for(int ii = 0; ii < JetPartonPair.size(); ii++){                            //Look at all the matched jets
		      if( (bTagStudy.getbTaggedJets(options))[jj] == JetPartonPair[ii].first){  //Check whether the considered jet can be matched!	
			JetTypeHistos[options].Fill( mcParticlesMatching[JetPartonPair[ii].second].type() );
		      }
		      else
			JetTypeHistos[options].Fill(25.);
		   }
		}
	      }
*/	
	      //Kinematic information for the Loose b-jets
              if((bTagStudy.getbTaggedJets(0)).size() >=2){
                h_StandardCosThetaLCSV.Fill(standardCosTheta);
                for(int ii = 0; ii<JetPartonPair.size(); ii++){ //Look at all the matched jets
                   if((bTagStudy.getbTaggedJets(0))[jj] == JetPartonPair[ii].first){  //Check whether the considered jet can be matched!
                      h_JetTypeLCSV.Fill(mcParticlesMatching[JetPartonPair[ii].second].type());
                   }
                   else{
                      h_JetTypeLCSV.Fill(25.);
                   }
                }
              }//End of Loose b-jets 

	      //Kinematic information for the Medium b-jets
              if((bTagStudy.getbTaggedJets(1)).size() >=2){
                 h_StandardCosThetaMCSV.Fill(standardCosTheta);
                 for(int ii = 0; ii<JetPartonPair.size(); ii++){ //Look at all the matched jets
                    if((bTagStudy.getbTaggedJets(1))[jj] == JetPartonPair[ii].first){  //Check whether the considered jet can be matched!
                        h_JetTypeMCSV.Fill(mcParticlesMatching[JetPartonPair[ii].second].type());
                    }
                    else{
                        h_JetTypeMCSV.Fill(25.);
                    }
                 }
              }//End of Medium b-jets 

	      //Kinematic information for the Tight b-jets
	      if((bTagStudy.getbTaggedJets(3)).size() >=2){
                 h_StandardCosThetaTCSV.Fill(standardCosTheta);
                 for(int ii = 0; ii<JetPartonPair.size(); ii++){ //Look at all the matched jets
		    if((bTagStudy.getbTaggedJets(3))[jj] == JetPartonPair[ii].first){  //Check whether the considered jet can be matched!
                        h_JetTypeTCSV.Fill(mcParticlesMatching[JetPartonPair[ii].second].type());
		    }
		    else{
                        h_JetTypeTCSV.Fill(25.);
		    }
	         }
	      }//End of Tight b-jets 
        
	   }//Only look at the two leading jets!
        }//End of TTbarJets!

        /////////////////////////////////////////////
        //   Count the selected number of events   //
        /////////////////////////////////////////////
        // 2 Loose CSV b-tags
        if((bTagStudy.getbTaggedJets(0)).size() >= 2 && eventselectedSemiMu) nSelectedMuLCSV++;
        if((bTagStudy.getbTaggedJets(0)).size() >= 2 && eventselectedSemiEl) nSelectedElLCSV++;
        if((bTagStudy.getbTaggedJets(0)).size() > 2) nLargeLCSVEvents++;

        // 2 Medium CSV b-tags
        if((bTagStudy.getbTaggedJets(1)).size() >= 2 && eventselectedSemiMu) nSelectedMuMCSV++;
        if((bTagStudy.getbTaggedJets(1)).size() >= 2 && eventselectedSemiEl) nSelectedElMCSV++;
        if((bTagStudy.getbTaggedJets(1)).size() > 2) nLargeMCSVEvents++;

        // 2 Tight CSV b-tags
        if((bTagStudy.getbTaggedJets(3)).size() >= 2 && eventselectedSemiMu) nSelectedMuTCSV++;
        if((bTagStudy.getbTaggedJets(3)).size() >= 2 && eventselectedSemiEl) nSelectedElTCSV++;
        if((bTagStudy.getbTaggedJets(3)).size() > 2) nLargeTCSVEvents++;

      }
      ///////////////////////////////////////////////////////////////
      //   End of loop for NrConsideredBTagOptions larger than 1   //
      ///////////////////////////////////////////////////////////////
      else if(NrConsideredBTagOptions == 1){
	bTagStudy.CalculateJets(selectedJets, BJetWP[ChosenBTagOption], LightJetWP[ChosenBTagOption], ChosenBTagOption);  
	//First float is the working point for the b-jets, the second is the one for the light jets!! 
        //Added integer points to which of the considered options is currently active. This should allow to obtain the final numbers at the end of the file!

        if(verbose > 3){
	  cout<<" Size of bTaggedJets: "<<(bTagStudy.getbTaggedJets(ChosenBTagOption)).size()<<" , of NonbTaggedJets: "<<(bTagStudy.getNonbTaggedJets(ChosenBTagOption)).size()<<" & of lightJets: "<<(bTagStudy.getLightJets(ChosenBTagOption)).size()<<endl;
	  cout<<" Index of BHadronic: "<<CorrectBHadronic<<" , Index of BLeptonic: "<<CorrectBLeptonic<<" , Index of quark1: "<<CorrectQuark1<<" & Index of quark2: "<<CorrectQuark2<<endl;
        }

        //Also check the correct jet combi:  --> Do this for both the 4- and 5-jet case (will make it possible to compare values!)
        if((bTagStudy.getbTaggedJets(ChosenBTagOption)).size() >= 2 && (bTagStudy.getLightJets(ChosenBTagOption)).size() >=2){ 

	    if( (bTagStudy.getLightJets(ChosenBTagOption)).size() >= 2)
	      bTagStudy.CorrectJetCombi(CorrectBHadronic, CorrectBLeptonic, CorrectQuark1, CorrectQuark2, ChosenBTagOption);      //4-jet case
	    if( (bTagStudy.getLightJets(ChosenBTagOption)).size() > 2)
	      bTagStudy.CorrectJetCombi5Jets(CorrectBHadronic, CorrectBLeptonic, CorrectQuark1, CorrectQuark2, ChosenBTagOption); //5-jet case
        }
        else{
	    if(verbose > 3) std::cout << " Event doesn't have two b-tagged jets and/or two light jets ! " << std::endl;
        }

      }//End of loop for NrConsideredBTagOptions equal to 1!
      h_StandardCosThetaNoBTag.Fill(standardCosTheta); 
      //---------------------------------------------------------------------------------------------------------------------------- End of bTagStudy class stuff

      //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
      //  Event selection choice (17/06/2014)  //
      //   --> Continue with 2 T b-tags        //
      //   --> No veto on light jets!          //
      //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
      //
      if((bTagStudy.getbTaggedJets(ChosenBTagOption)).size() < 2) continue;
      
      // Count the number of events:
      if(eventselectedSemiMu) nSelectedMu++;
      if(eventselectedSemiEl) nSelectedEl++;
     
      //--> Now check whether improvement is obtained when adding another light jet (so look at three jet possibilities!)
      if( (bTagStudy.getLightJets(3)).size() >2 && CorrectQuark1 != 9999 && CorrectQuark2 != 9999 && CorrectBLeptonic != 9999 && CorrectBHadronic != 9999){
	 if(verbose > 3){
	   std::cout << " \n Check whether the third light jet can be matched with the correct quarks: " << std::endl;
	   std::cout << " Correct quark indices         : " << CorrectQuark1 << " & " << CorrectQuark2 << std::endl;
	   std::cout << " Three first light jet indices : " << (bTagStudy.getLightJets(3))[0] << " , " << (bTagStudy.getLightJets(3))[1] << " & " << (bTagStudy.getLightJets(3))[2] << std::endl;
 	 }
      }
  
      ////////////////////////////////
      //  Mlb and Mqqb information  //
      ////////////////////////////////
      //
      vector<int> LeptonicB, HadronicB, Quark1, Quark2;
      float MlbCorrect = 0, MqqbCorrect = 0;
      if(CorrectBLeptonic != 9999) MlbCorrect = (*selectedLepton+*selectedJets[CorrectBLeptonic]).M();
      if(CorrectBHadronic != 9999 && CorrectQuark1 != 9999 && CorrectQuark2 != 9999) MqqbCorrect = (*selectedJets[CorrectBHadronic] + *selectedJets[CorrectQuark1] + *selectedJets[CorrectQuark2]).M();
      h_MlbMqqbCorrectAll.Fill(MqqbCorrect,MlbCorrect);

      vector<int> CorrectValues;  //Find better solution!
      CorrectValues.clear();
      CorrectValues.push_back(CorrectBLeptonic);
      CorrectValues.push_back(CorrectBHadronic);
      CorrectValues.push_back(CorrectQuark1);
      CorrectValues.push_back(CorrectQuark2);
     
      //****************************************//
      //  Loop can be kept, since if NrConsideredBTagOptions == 1 only 1 case will be considered  //
      //    --> If this is the case, force the option to be the one chosen above!!                //
      //******************************************************************************************//
      for(int Option = 0; Option < NrConsideredBTagOptions; Option++){
	if(NrConsideredBTagOptions == 1) Option = ChosenBTagOption;    //Force Option to be equal to the one chosen!

	mlbStudy.calculateChiSquared(CorrectValues, bTagStudy.getbTaggedJets(Option), bTagStudy.getLightJets(Option), selectedLepton, selectedJets, MassMlb, SigmaMlb, MassMqqb, SigmaMqqb);

	if(bTagStudy.getbTaggedJets(Option).size() >1 && bTagStudy.getLightJets(Option).size() > 1){

	//TODO --> Update!!
	  if( (CorrectBLeptonic == (bTagStudy.getbTaggedJets(Option))[0] || CorrectBLeptonic == (bTagStudy.getbTaggedJets(Option))[1]) &
	      (CorrectBHadronic == (bTagStudy.getbTaggedJets(Option))[1] || CorrectBHadronic == (bTagStudy.getbTaggedJets(Option))[0]) &&
	      (CorrectQuark1 == (bTagStudy.getLightJets(Option))[0]      || CorrectQuark1 == (bTagStudy.getLightJets(Option))[1]) &&
	      (CorrectQuark2 == (bTagStudy.getLightJets(Option))[0]      || CorrectQuark2 == (bTagStudy.getLightJets(Option))[1])){ 
		h_MlbMqqbCorrectChosen.Fill(MqqbCorrect,MlbCorrect); 
	  } 
	
	  mlbStudy.calculateEfficiency(Option, CorrectValues, bTagStudy.getbTaggedJets(Option), bTagStudy.getLightJets(Option), NrConsideredBTagOptions );

	  // Match the jet number to the particles
	  mlbStudy.getIndices(mlbStudy.getLowestChiSqIndex());   //Or use .getLowestChiSq4JetsIndex() when only considering the 4-jet case!
	  LeptonicB.push_back(bTagStudy.getbTaggedJets(Option)[mlbStudy.getChosenBLept()]);   //!! Probably not possible to use push_back ...
	  HadronicB.push_back(bTagStudy.getbTaggedJets(Option)[mlbStudy.getChosenBHadr()]);
          Quark1.push_back(bTagStudy.getLightJets(Option)[mlbStudy.getChosenQuark1()]);
	  Quark2.push_back(bTagStudy.getLightJets(Option)[mlbStudy.getChosenQuark2()]);

	  //Looking how often the correct event is found (in the 5-jet case):    --> Have to ask .getIndices() since otherwise last call from class will be used!
	  if( (CorrectBLeptonic == bTagStudy.getbTaggedJets(Option)[mlbStudy.getChosenBLept()] ) &&
	      (CorrectBHadronic == bTagStudy.getbTaggedJets(Option)[mlbStudy.getChosenBHadr()] ) &&
	      (CorrectQuark1    == bTagStudy.getLightJets(Option)[mlbStudy.getChosenQuark1()] || CorrectQuark1 == bTagStudy.getLightJets(Option)[mlbStudy.getChosenQuark2()] ) &&
	      (CorrectQuark2    == bTagStudy.getLightJets(Option)[mlbStudy.getChosenQuark1()] || CorrectQuark2 == bTagStudy.getLightJets(Option)[mlbStudy.getChosenQuark2()] ) ){
		CorrectEventFound5Jets++;
	  }

	  mlbStudy.getIndices(mlbStudy.getLowestChiSq4JetsIndex());
          //Looking how often the correct event is found (in the 5-jet case):    --> Have to ask .getIndices() since otherwise last call from class will be used!
          if( (CorrectBLeptonic == bTagStudy.getbTaggedJets(Option)[mlbStudy.getChosenBLept()] ) &&
              (CorrectBHadronic == bTagStudy.getbTaggedJets(Option)[mlbStudy.getChosenBHadr()] ) &&
              (CorrectQuark1    == bTagStudy.getLightJets(Option)[mlbStudy.getChosenQuark1()] || CorrectQuark1 == bTagStudy.getLightJets(Option)[mlbStudy.getChosenQuark2()] ) &&
              (CorrectQuark2    == bTagStudy.getLightJets(Option)[mlbStudy.getChosenQuark1()] || CorrectQuark2 == bTagStudy.getLightJets(Option)[mlbStudy.getChosenQuark2()] ) )
                CorrectEventFound4Jets++;
        }
      }

      //************************************************************************************//   
      //  Only go on to produce the .lhco output if final event selection choice is made !! //
      //************************************************************************************//
      //
      if( FinalEventSelectionChoiceIsMade == false || RecoLHCOOutput == false) continue;

      //Array of b-tagged jets and light jets for the possible configurations
      const int NrCombi = 2;  //--> Is only 2 in the case of a chi-squared applied!
      TLorentzVector* LeptBArray[NrCombi] = {selectedJets[(bTagStudy.getbTaggedJets(ChosenBTagOption))[mlbStudy.getChosenBLept()]],
					     selectedJets[(bTagStudy.getbTaggedJets(ChosenBTagOption))[mlbStudy.getChosenBLept()]]};
      TLorentzVector* HadrBArray[NrCombi] = {selectedJets[(bTagStudy.getbTaggedJets(ChosenBTagOption))[mlbStudy.getChosenBHadr()]],
					     selectedJets[(bTagStudy.getbTaggedJets(ChosenBTagOption))[mlbStudy.getChosenBHadr()]]};
      TLorentzVector* Quark1Array[NrCombi] ={selectedJets[(bTagStudy.getLightJets(ChosenBTagOption))[mlbStudy.getChosenQuark1()]], 
					     selectedJets[(bTagStudy.getLightJets(ChosenBTagOption))[mlbStudy.getChosenQuark2()]]};
      TLorentzVector* Quark2Array[NrCombi] ={selectedJets[(bTagStudy.getLightJets(ChosenBTagOption))[mlbStudy.getChosenQuark2()]], 
					     selectedJets[(bTagStudy.getLightJets(ChosenBTagOption))[mlbStudy.getChosenQuark1()]]}; 


      /////////////////////////////////////////////
      //  Filling of LHCO files for reco events  //
      /////////////////////////////////////////////
      vector<TLorentzVector*> LHCORecoVector(6);
      vector<int> MadGraphRecoId(6,4);
      //Need to distinguish between charge and lepton type
      for(int ConsideredCombi = 0; ConsideredCombi <NrCombi; ConsideredCombi++){  
	      string JetCombiString = static_cast<ostringstream*>( &(ostringstream() << ConsideredCombi) )->str();

	      /////////////////////////////////////////////////////////
	      //   Reconstructing Lepton & Neutrino (partially)      //
	      //   -> Neutrino Pt and M needs to be known for .lhco  //
	      /////////////////////////////////////////////////////////  
	      //  --> Only do this part for the combination which is actually chosen!
//	      Neutrino.clear();
	      NeutrinoPx = -(*selectedLepton+*LeptBArray[ConsideredCombi]+*HadrBArray[ConsideredCombi]+*Quark1Array[ConsideredCombi]+*Quark2Array[ConsideredCombi]).Px();
	      NeutrinoPy = -(*selectedLepton+*LeptBArray[ConsideredCombi]+*HadrBArray[ConsideredCombi]+*Quark1Array[ConsideredCombi]+*Quark2Array[ConsideredCombi]).Py();
	      Neutrino.SetPxPyPzE(NeutrinoPx, NeutrinoPy, 0.0, sqrt(NeutrinoPx*NeutrinoPx + NeutrinoPy*NeutrinoPy));
	      Neutrino.SetPxPyPzE(NeutrinoPx,NeutrinoPy,0.0, Neutrino.Pt());  //Reset the Neutrino Energy to get the correct precision
	      if(verbosity > 3) std::cout << " Mass value for the neutrino : " << Neutrino.M() << " \n" << std::endl;

	      if(LeptonRecoCharge < 0.0 ){ //Negative lepton events
		      if(verbosity>4) cout << " Looking at negative lepton events for Reco LHCO files " << endl;
		      LHCORecoVector[0] = HadrBArray[ConsideredCombi]; 
		      LHCORecoVector[1] = Quark1Array[ConsideredCombi];
		      LHCORecoVector[2] = Quark2Array[ConsideredCombi];
		      LHCORecoVector[3] = LeptBArray[ConsideredCombi];
		      LHCORecoVector[4] = selectedLepton;
		      LHCORecoVector[5] = &Neutrino;
		      if(eventselectedSemiEl){//Negative electron
			      MadGraphRecoId[1] = 1;
			      MadGraphRecoId[2] = 6;
			      if(ConsideredCombi == 0) NumberNegRecoEl++;  //Only need to raise the eventNumber for one combination of the 4!!
			      if(verbosity > 4) cout << " Event : " << ievt << " with Number " << NumberNegRecoEl << " sent to LHCO Reco output (Negative electron) " << endl;
			      if(NumberNegRecoEl == 1) outFileReco[3*4+ConsideredCombi].open(("TTbarSemiLepton_Reco_NegativeElectron_JetCombi"+JetCombiString+".lhco").c_str());
			      lhcoOutput.LHCOEventRecoOutput(3,outFileReco[3*4+ConsideredCombi], NumberNegRecoEl, LHCORecoVector, MadGraphRecoId);
			      EventInfoFile << "     " << NumberNegRecoEl << endl;
		      }
		      if(eventselectedSemiMu){//Negative muon
			      MadGraphRecoId[1] = 2;
			      MadGraphRecoId[2] = 6;
			      if(ConsideredCombi == 0) NumberNegRecoMu++;
			      if(verbosity > 4) cout << " Event : " << ievt << " with Number " << NumberNegRecoMu << " sent to LHCO Reco output (Negative muon) " << endl;
			      if(NumberNegRecoMu == 1) outFileReco[1*4+ConsideredCombi].open(("TTbarSemiLepton_Reco_NegativeMuon_JetCombi"+JetCombiString+".lhco").c_str());
			      lhcoOutput.LHCOEventRecoOutput(1, outFileReco[1*4+ConsideredCombi], NumberNegRecoMu, LHCORecoVector, MadGraphRecoId);
			      EventInfoFile << "     " << NumberNegRecoMu << endl;
		      }
	      }//End of negative lepton

	      if(LeptonRecoCharge > 0.0 ){ //Positive lepton events
		      if(verbosity>4) cout << " Looking at positive lepton events for Reco LHCO files " << endl;
		      LHCORecoVector[0] = LeptBArray[ConsideredCombi];
		      LHCORecoVector[1] = selectedLepton;
		      LHCORecoVector[2] = &Neutrino;
		      LHCORecoVector[3] = HadrBArray[ConsideredCombi];
		      LHCORecoVector[4] = Quark1Array[ConsideredCombi];
		      LHCORecoVector[5] = Quark2Array[ConsideredCombi];
		      if(eventselectedSemiEl){//Positive electron
			      MadGraphRecoId[1] = 1;
			      MadGraphRecoId[2] = 6;
			      if(ConsideredCombi == 0) NumberPosRecoEl++;
			      if(verbosity > 4) cout << " Event : " << ievt << " with Number " << NumberPosRecoEl << " sent to LHCO Reco output (Positive electron) " << endl;
			      if(NumberPosRecoEl == 1) outFileReco[2*4+ConsideredCombi].open(("TTbarSemiLepton_Reco_PositiveElectron_JetCombi"+JetCombiString+".lhco").c_str());
			      lhcoOutput.LHCOEventRecoOutput(2,outFileReco[2*4+ConsideredCombi], NumberPosRecoEl, LHCORecoVector, MadGraphRecoId);	 
			      EventInfoFile << "     " << NumberPosRecoEl << endl;
		      }
		      if(eventselectedSemiMu){//Positive muon
			      MadGraphRecoId[1] = 2;
			      MadGraphRecoId[2] = 6;
			      if(ConsideredCombi == 0) NumberPosRecoMu++;
			      if(verbosity > 4) cout << " Event : " << ievt << " with Number " << NumberPosRecoMu << " sent to LHCO Reco output (Positive muon) " << endl;
			      if(NumberPosRecoMu == 1) outFileReco[0*4+ConsideredCombi].open(("TTbarSemiLepton_Reco_PositiveMuon_JetCombi"+JetCombiString+".lhco").c_str());
			      lhcoOutput.LHCOEventRecoOutput(0, outFileReco[0*4+ConsideredCombi], NumberPosRecoMu, LHCORecoVector, MadGraphRecoId);
			      EventInfoFile << "     " << NumberPosRecoMu << endl;
		      }
	      }//End of positive lepton

      }//End of loop over the different jet combinations  

    }			//loop on events

    //--------------------  Sigma for W Mass and Top Mass  --------------------
    //histo1D["WMass"]->Fit("gaus","Q");     
    //histo1D["TopMass"]->Fit("gaus","Q");
    //std::cout << " sigma values : " << histo1D["WMass"]->GetFunction("gaus")->GetParameter(2) << " " << histo1D["TopMass"]->GetFunction("gaus")->GetParameter(2) << std::endl;
    //std::cout << " mass values : " << histo1D["WMass"]->GetFunction("gaus")->GetParameter(1) << " " << histo1D["TopMass"]->GetFunction("gaus")->GetParameter(1) << std::endl;
    histo1D["MlbMass"]->Fit("gaus","Q");
    histo1D["MqqbMass"]->Fit("gaus","Q");
    std::cout << " values for Mlb :" << histo1D["MlbMass"]->GetFunction("gaus")->GetParameter(1) << " +- " << histo1D["MlbMass"]->GetFunction("gaus")->GetParameter(2) << std::endl;
    std::cout << " values for Mqqb :" << histo1D["MqqbMass"]->GetFunction("gaus")->GetParameter(1) << " +- " << histo1D["MqqbMass"]->GetFunction("gaus")->GetParameter(2) << std::endl;

    cout<<endl; 
    cout << "-> " << nSelectedMu << " mu+jets events where selected from which " << nSelectedMuLCSV << " have two or more Light wp CSV b-tags, " << nSelectedMuMCSV << " have two or more Medium wp CSV b-tags and " << nSelectedMuTCSV << " have two or more Tight wp CSV b-tags " << endl;
    cout << "-> " << nSelectedEl << " e+jets events where selected from which " << nSelectedElLCSV << " have two or more Light wp CSV b-tags, " << nSelectedElMCSV << " have two or more Medium wp CSV b-tags and " << nSelectedElTCSV << " have two or more Tight wp CSV b-tags " << endl;
    //cout << "-> " << nLargeMCSVEvents << " events with more than 2 Medium CSV b-tags ( " << nLargeTCSVEvents << " with 2 Tight CSV b-tags) --> Reject these events! " << endl;

    cout << " " << endl;
    if(verbosity>0) cout << "---> Number of events with correct semileptonic event content on generator level: " << NumberCorrectEvents << " (semiMuon, semiElec) : ( " << NumberPositiveMuons+NumberNegativeMuons << " , " << NumberPositiveElectrons+NumberNegativeElectrons << " ) " << endl;
    cout << " " << endl;
    
    //Output for 4 of 5 jet-case!
    std::cout << " ********************************************* " << std::endl;
    std::cout << " **  Output for 4- or 5-jet case choice !!  ** " << std::endl;
    std::cout << " ********************************************* " << std::endl;
    std::cout << " BTagged events with at least two light jets : " << bTagStudy.getNrEventsWithTwoLightJetsAndBTagged(ChosenBTagOption)   << std::endl;
    std::cout << " BTagged events with more than two light jets: " << bTagStudy.getNrEventsWithThreeLightJetsAndBTagged(ChosenBTagOption) << std::endl;
    std::cout << " Number of times the third jet is one of the correct ones : " << bTagStudy.getNrTimesThirdJetIsActualQuark(ChosenBTagOption)  ;
    std::cout << "  ==> In " << (float)((float)(bTagStudy.getNrTimesThirdJetIsActualQuark(ChosenBTagOption)*100.0)/((float)bTagStudy.getNrEventsWithThreeLightJetsAndBTagged(ChosenBTagOption))) << " \% of the cases an improvement is found when adding a third light jet !! " << std::endl;
    std::cout << " Number of times the third light jet is correct (no b's)  : " << bTagStudy.getNrTimesThirdJetIsCorrectQuark(ChosenBTagOption) << std::endl;
    std::cout << " ------------------------ " << std::endl;
    std::cout << " Four jets case (all correct vs one wrong): " << bTagStudy.getNrCorrectMatchedEvts(ChosenBTagOption) << " vs " << bTagStudy.getNrWrongMatchedEvts(ChosenBTagOption) << std::endl;
    std::cout << " Signal over background in case of 4 jets : " << bTagStudy.getSignalOverBkg(ChosenBTagOption) << std::endl;
    std::cout << " Five jets case (all correct vs one wrong): " << bTagStudy.getNrCorrectMatchedEvts5Jets(ChosenBTagOption) << " vs " << bTagStudy.getNrWrongMatchedEvts5Jets(ChosenBTagOption) << std::endl;
    std::cout << " Signal over background in case of 5 jets : " << bTagStudy.getSignalOverBkg5Jets(ChosenBTagOption) << std::endl;
    std::cout << " Signal over background for 4 & 5 jets    : " << (float)(bTagStudy.getNrCorrectMatchedEvts(ChosenBTagOption)+bTagStudy.getNrCorrectMatchedEvts5Jets(ChosenBTagOption))/(bTagStudy.getNrWrongMatchedEvts(ChosenBTagOption)+bTagStudy.getNrWrongMatchedEvts5Jets(ChosenBTagOption)) << std::endl;
    std::cout << " ------------------------------------------------------------------------- \n" << std::endl;

    //////////////////////////////
    //  Jet combination output  //
    //////////////////////////////
    //--> Save directly to .tex output as a table
    ofstream eventSelOutput;
    int BTagOptionOfInterest;
    std::string OptionName5Jets[6];
    std::string OptionName4Jets[6];
    if(NrConsideredBTagOptions > 1){
      eventSelOutput.open("/user/aolbrech/GitTopTree_Feb2014/TopBrussels/AnomalousCouplings/eventSelectionChoiceTables.tex");
      BTagOptionOfInterest = 7;
      for(int ii = 0; ii < 6; ii++){OptionName4Jets[ii] = OptionName[ii]; OptionName5Jets[ii] = OptionName[ii];}
    }
    if(NrConsideredBTagOptions == 1){
      eventSelOutput.open("/user/aolbrech/GitTopTree_Feb2014/TopBrussels/AnomalousCouplings/eventSelectionTableForChosenCombination.tex");
      for(int ii = 0; ii < 6; ii++){OptionName4Jets[ii] = " 4 jet case with "+OptionName[ii]; OptionName5Jets[ii] = " 5 jet case with "+OptionName[ii];}
      BTagOptionOfInterest = ChosenBTagOption;
    }

    bTagStudy.ReturnTable(OptionName4Jets, OptionName5Jets, 0, NrConsideredBTagOptions, eventSelOutput, BTagOptionOfInterest);  //0 stands for all 4 particles, 1 for b-jets and 2 for light jets!
    bTagStudy.ReturnTable(OptionName4Jets, OptionName5Jets, 1, NrConsideredBTagOptions, eventSelOutput, BTagOptionOfInterest);
    bTagStudy.ReturnTable(OptionName4Jets, OptionName5Jets, 2, NrConsideredBTagOptions, eventSelOutput, BTagOptionOfInterest);
    eventSelOutput.close();
        
    //////////////////////////////
    //  Mlb combination output  //
    //////////////////////////////
    std::cout << " Number of correct reconstructed events for 5-jet case using Mlb method : " << CorrectEventFound5Jets << std::endl;
    std::cout << " Number of correct reconstructed events for 4-jet case using Mlb method : " << CorrectEventFound4Jets << std::endl;

    ofstream mlbOutput;
    int OptionOfInterest;
    if(NrConsideredBTagOptions > 1 ){ 
	mlbOutput.open("/user/aolbrech/GitTopTree_Feb2014/TopBrussels/AnomalousCouplings/mlbChoiceTables.tex");              
	OptionOfInterest = 7;
    } 
    //In above case the OptionOfInterest variable will not be used inside the class! Is only used when only one bTag should be considered!
    if(NrConsideredBTagOptions == 1){ 
	mlbOutput.open("/user/aolbrech/GitTopTree_Feb2014/TopBrussels/AnomalousCouplings/mlbTableForChosenCombination.tex"); 
	OptionOfInterest = ChosenBTagOption;
	OptionName[ChosenBTagOption+2] = " Pure 5 jet case with "+OptionName[ChosenBTagOption];
        OptionName[ChosenBTagOption+1] = " 4 jet case with      "+OptionName[ChosenBTagOption];
	OptionName[ChosenBTagOption] =   " 5 jet case with      "+OptionName[ChosenBTagOption];
    }
    mlbStudy.saveNumbers(OptionName, 0, NrConsideredBTagOptions, mlbOutput, OptionOfInterest );  //All 4 jets correctly matched
    mlbStudy.saveNumbers(OptionName, 1, NrConsideredBTagOptions, mlbOutput, OptionOfInterest );  //Only b-jets correctly matched
    mlbOutput.close();

   //vector<std::string> HistoTitle = {"ChiSqCorrect","ChiSqCorrectFound","ChiSqMinimum","ChiSqNotMinimum","ChiSqWrong","ChiSqCorrectWhenMatched","ChiSqMinimumWhenMatched","ChiSqNotMinimumWhenMatched","ChiSqAllWhenNotMatched"};
   /*vector<std::string> Histoname = {"#chi^{2} distribution for the correct combination (all events)",
				    "#chi^{2} distribution for the correct combination (when it is one of the 6 considered)",
				    "#chi^{2} distribution of the minimal combination considered (all events)",
				    "#chi^{2} distribution of the non-minimal combinations considered (all events)",
				    "#chi^{2} distribution of the non-correct combinations considered (all events)",
				    "#chi^{2} distribution for the correct combination (matched events only)",
				    "#chi^{2} distribution for the minimal combination (matched events only)",
				    "#chi^{2} distribution for the non-minimal combinations (matched events only)",
				    "#chi^{2} distribution for all the combinations when the event is not matched"};
*/
   std::string Title[3] = {"5Jets","4Jets","Pure5Jets"};
   std::string Name[3] =  {" - 5 jets case) "," - 4 jets case) "," - pure 5 jets case) "};
 
   TH1F *h_ChiSqCorrect[3], *h_ChiSqCorrectFound[3], *h_ChiSqMinimum[3],* h_ChiSqNotMinimum[3], *h_ChiSqWrong[3];
   TH1F *h_ChiSqCorrectWhenMatched[3], *h_ChiSqMinimumWhenMatched[3], *h_ChiSqNotMinimumWhenMatched[3], *h_ChiSqAllWhenNotMatched[3]; 
   for(int ii = 0; ii < 3; ii++){
     h_ChiSqCorrect[ii]     =new TH1F(("ChiSqCorrect"+Title[ii]).c_str(),     ("#chi^{2} distribution for the correct combination (all events"+Name[ii]).c_str() ,              500,0,500);
     h_ChiSqCorrectFound[ii]=new TH1F(("ChiSqCorrectFound"+Title[ii]).c_str(),("#chi^{2} distribution for the correct combination (when found"+Name[ii]).c_str(),               500,0,500);
     h_ChiSqMinimum[ii]     =new TH1F(("ChiSqMinimum"+Title[ii]).c_str(),     ("#chi^{2} distribution of the minimal combination considered (all events"+Name[ii]).c_str(),     500,0,500);
     h_ChiSqNotMinimum[ii]  =new TH1F(("ChiSqNotMinimum"+Title[ii]).c_str(),  ("#chi^{2} distribution of the non-minimal combinations considered (all events"+Name[ii]).c_str(),500,0,500);
     h_ChiSqWrong[ii]       =new TH1F(("ChiSqWrong"+Title[ii]).c_str(),       ("#chi^{2} distribution of the non-correct combinations considered (all events"+Name[ii]).c_str(),500,0,500);
     
     h_ChiSqCorrectWhenMatched[ii]   =new TH1F(("ChiSqCorrectWhenMatched"+Title[ii]).c_str(),("#chi^{2} distribution for the correct combination (matched events only"+Name[ii]).c_str(),        100,0,100);
     h_ChiSqMinimumWhenMatched[ii]   =new TH1F(("ChiSqMinimumWhenMatched"+Title[ii]).c_str(),("#chi^{2} distribution for the minimal combination (matched events only"+Name[ii]).c_str(),        500,0,500);
     h_ChiSqNotMinimumWhenMatched[ii]=new TH1F(("ChiSqNotMinimumWhenMatched"+Title[ii]).c_str(),("#chi^{2} distribution for the non-minimal combinations (matched events only"+Name[ii]).c_str(),500,0,500);
     h_ChiSqAllWhenNotMatched[ii]    =new TH1F(("ChiSqAllWhenNotMatched"+Title[ii]).c_str(),    ("#chi^{2} distribution for all the combinations (non-matched evens only"+Name[ii]).c_str(),     500,0,500);
   }

   for(int jetCase = ChosenBTagOption; jetCase <(ChosenBTagOption+3); jetCase++){
      std::cout << " Size of ChiSqCorect : " << (mlbStudy.getChiSqCorrect(jetCase)).size() << " for jetCase : " << jetCase << std::endl;
      for(int ii = 0; ii < (mlbStudy.getChiSqCorrect(jetCase)).size(); ii++){
	//In mlb Class the ChosenBTagOption, ChosenBTagOption+1 and ChosenBTagOption+2 indices are filled (from the 6 options).
	//However to reduce the number of TH1F's made, this is translated to 0,1 and 2 here in the analyzer!!
	h_ChiSqCorrect[jetCase-ChosenBTagOption]->Fill((mlbStudy.getChiSqCorrect(jetCase))[ii]);
        h_ChiSqCorrectFound[jetCase-ChosenBTagOption]->Fill((mlbStudy.getChiSqCorrectFound(jetCase))[ii]);
	h_ChiSqMinimum[jetCase-ChosenBTagOption]->Fill((mlbStudy.getChiSqMinimum(jetCase))[ii]);
	h_ChiSqNotMinimum[jetCase-ChosenBTagOption]->Fill((mlbStudy.getChiSqNotMinimum(jetCase))[ii]);
	h_ChiSqWrong[jetCase-ChosenBTagOption]->Fill((mlbStudy.getChiSqWrong(jetCase))[ii]);

	h_ChiSqCorrectWhenMatched[jetCase-ChosenBTagOption]->Fill((mlbStudy.getChiSqCorrectWhenMatched(jetCase))[ii]);
	h_ChiSqMinimumWhenMatched[jetCase-ChosenBTagOption]->Fill((mlbStudy.getChiSqMinimumWhenMatched(jetCase))[ii]);
	h_ChiSqNotMinimumWhenMatched[jetCase-ChosenBTagOption]->Fill((mlbStudy.getChiSqNotMinimumWhenMatched(jetCase))[ii]);
	h_ChiSqAllWhenNotMatched[jetCase-ChosenBTagOption]->Fill((mlbStudy.getChiSqAllWhenNotMatched(jetCase))[ii]);
      }
   }   

    //Close the LHCO Output files!
    for(int ii = 0; ii<16; ii++){
      if(ii < 4) outFile[ii].close();	
      outFileReco[ii].close();
    }
    EventInfoFile.close();
    
    //TFile* fout = new TFile("GeneratorOutput.root","RECREATE");
    fout->cd();

    //Write down the different Chi-Sq histograms obtained from the mlb class!
    for(int ii = 0; ii < 3; ii++){
       h_ChiSqCorrect[ii]->Write();
       h_ChiSqCorrectFound[ii]->Write();
       h_ChiSqMinimum[ii]->Write();
       h_ChiSqNotMinimum[ii]->Write();
       h_ChiSqWrong[ii]->Write();

       h_ChiSqCorrectWhenMatched[ii]->Write();
       h_ChiSqMinimumWhenMatched[ii]->Write();
       h_ChiSqNotMinimumWhenMatched[ii]->Write();
       h_ChiSqAllWhenNotMatched[ii]->Write();
    }

    h_StandardCosThetaNoEvtSel.Write();
    h_StandardCosThetaNoBTag.Write();
    h_StandardCosThetaLCSV.Write();
    h_StandardCosThetaAllLCSV.Write();
    h_StandardCosThetaMCSV.Write();
    h_StandardCosThetaTCSV.Write();
    h_JetTypeLargeLCSVEvents.Write();
    h_JetTypeLargeLCSVLeadingPtEvents.Write();
    h_JetTypeLCSVLightJetsLeadingPt.Write();
    h_JetTypeLargeMCSVEvents.Write();
    h_JetTypeLargeMCSVLeadingPtEvents.Write();
    h_JetTypeLargeTCSVEvents.Write();
    h_JetTypeLargeTCSVLeadingPtEvents.Write();
    h_JetTypeLCSV.Write();
    h_JetTypeMCSV.Write();
    h_JetTypeTCSV.Write();
    h_JetTypeLCSVLightJets.Write();
    h_CSVDiscrLCSVLightJets.Write();
    h_CSVDiscrLCSVLightJetsLeadingPt.Write();
    h_CorrectBLeptCSVDiscr.Write();
    h_CorrectBHadrCSVDiscr.Write();
    h_CorrectQuark1CSVDiscr.Write();
    h_CorrectQuark2CSVDiscr.Write();

    h_Quark1JetNumber.Write();
    h_Quark2JetNumber.Write();

    h_WMass.Write();
    h_TopMass.Write();
    h_MlbMass.Write();
    h_MqqbMass.Write();
    
    h_CosThetaReco.Write();
    h_NeutrinoEta.Write();

    h_MlbMqqbCorrectAll.Write();
    h_MlbMqqbCorrectChosen.Write();
    h_MlbMqqbWrongOne.Write();
    h_MlbMqqbWrongTwo.Write();
    //fout->Close();
    
    //////////////
    // CLEANING //
    //////////////
    
    //Delte TLorentzVector used for standardCosTheta calculation (Still need to check whether this works correctly when using multiple datasets)
    //delete sTop;
    //delete WLeptTRF;
    //delete leptonWRF;
    //delete Top, TopBar, Bottom, BottomBar, Lepton, NeutrinoMC, WPlus, WMinus, Light, LightBar;
 

    if (jecUnc) delete jecUnc;
    if (jetTools) delete jetTools;
  
    //important: free memory
    treeLoader.UnLoadDataset();
    
  }				//loop on datasets
  
  //Once everything is filled ...
  if (verbose > 0)
    cout << " We ran over all the data ;-)" << endl;
  
  /////////////////////////
  // Write out the plots //
  /////////////////////////
  
  mkdir("DemoPlots",0777);
  
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++){    
    MultiSamplePlot *temp = it->second;
    string name = it->first;
    temp->Draw(name); //, true, true, true, true, true, 1, false);
    temp->Write(fout, name, true, "DemoPlots/");
  }
  
  //Selection tables
  selecTableSemiEl.TableCalculator(false, true, true, true, true);
  string selectiontableEl = "SelectionTable_BTAG_SEMIEL.tex";
  selecTableSemiEl.Write(selectiontableEl.c_str());
  
  // Do some special things with certain plots (normalize, BayesDivide, ... )
  if (verbose > 0)
    cout << "Treating the special plots." << endl;
  
  delete fout;
  delete tcdatasets;
  delete tcAnaEnv;
  delete configTree;
  
  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
  
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "           hasn't crashed yet ;-)           " << endl;
  cout << "********************************************" << endl;
  
  return 0;
}
