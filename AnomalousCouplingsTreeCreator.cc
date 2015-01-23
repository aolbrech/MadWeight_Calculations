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
#include "PersonalClasses/Style.C"                                                 //CHECK if this works!
#include "TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"
#include "TopTreeAnalysisBase/Tools/interface/LeptonTools.h"
#include "TLorentzVector.h"

//Specific code for anomalous couplings analysis:
#include "AnomalousCouplings/PersonalClasses/interface/LHCOOutput.h"
#include "AnomalousCouplings/PersonalClasses/interface/BTagStudy.h"
#include "AnomalousCouplings/PersonalClasses/interface/MlbStudy.h"
#include "AnomalousCouplings/PersonalClasses/interface/TFCreation.h"
#include "AnomalousCouplings/PersonalClasses/interface/TFnTuple.h"

using namespace std;
using namespace reweight;
using namespace TopTree;

int main (int argc, char *argv[])
{
  string rootFileName = "AnomCouplings.root";
  
  clock_t start = clock();
  
  cout << "*************************************************************" << endl;
  cout << " Beginning of the program for creating the AnomCoupl Trees ! " << endl;
  cout << "*************************************************************" << endl;
  
  //SetStyle if needed
  //setTDRStyle(); 
  setMyStyle();

  /////////////////////
  // Configuration
  /////////////////////

  //xml file
  string xmlFileName ="PersonalClasses/config/myAnomCouplConfig.xml";

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
  bool GenLHCOOutput = true;
  bool RecoLHCOOutput = false;  
  bool FinalEventSelectionChoiceIsMade = false;

  bool CalculateResolutions = false; // If false, the resolutions will be loaded from a previous calculation
  bool CalculateTF = false;
  bool CalculateBTag = false;

  //Values needed for bTag study (select which of the 6 b-tag options is optimal!)
  const int NrConsideredBTagOptions = 1;   //Make sure this number is also the same in the bTagStudy class!!
  int ChosenBTagOption;
  if(NrConsideredBTagOptions > 1){ ChosenBTagOption = 7; GenLHCOOutput = false; RecoLHCOOutput = false;}
  else ChosenBTagOption = 3;    //3 = 2 T b-tags!  

  int ChiSqCutValue =51;  //The Chi-sq values in the mlb method has to be larger than this value! (Put on 51 to include all events, since the chi-sq is set manually to a maximum of 49.5)

  /////////////////////////////////////
  //   Initializiation of counters   //
  ///////////////////////////////////// 
  int CorrectEventMlbMqqb[NrConsideredBTagOptions];
  int WrongEventMlbMqqb[NrConsideredBTagOptions];
  for(int ii = 0; ii < NrConsideredBTagOptions; ii++){
    CorrectEventMlbMqqb[ii] = 0;
    WrongEventMlbMqqb[ii] = 0;
  }

  float BJetWP[6] = {0.244,0.679,0.679,0.898,0.898,0.898};
  float LightJetWP[6] = {0.244,0.679,0.244,0.898,0.679,0.244};
  
  std::string OptionName[6] = {"2 L b-tags             ", //#0
    			       "2 M b-tags             ", //#1
			       "2 M b-tags, light L-veto", //#2
			       "2 T b-tags             ", //#3
			       "2 T b-tags, light M-veto", //#4
			       "2 T b-tags, light L-veto"};//#5

  std::string ChiSqCutValueStr;
  ostringstream convert; convert << ChiSqCutValue;
  if(ChiSqCutValue != 51) ChiSqCutValueStr = "_ChiSqSmallerThan"+convert.str();
  else ChiSqCutValueStr = "";

  //Counters needed for N(4) and N(5) categories:
  int FourJetsCategory = 0, FiveJetsCategory = 0, FourJetsMatched = 0, FiveJetsMatched = 0, FourJetsGoodCombiChosen = 0, FiveJetsAsFourGoodCombiChosen = 0, FiveJetsAsFiveGoodCombiChosen = 0;

  /////////////////////////////
  //   Which systematics ?   //
  /////////////////////////////
  std::string doJESShift = "nominal"; // Other options are "minus" and "plus"
  std::string doJERShift = "nominal";
  std::string doLeptonSFShift = "Nominal";  //Other options are "Minus" and "Plus" (Capital letters are needed!
  std::string doLumiWeightShift = "Nominal";  //Other options are "Minus" and "Plus"
 
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
  histo1D["MqqbMass"]= new TH1F("MqqbMass","MqqbMass",400,0,500);

  histo1D["genPt_Muon"] = new TH1F("genPt_Muon","genPt_Muon",400,0,200);
  histo1D["recoPt_Muon"] = new TH1F("recoPt_Muon","recoPt_Muon",400,0,200);
  histo1D["genPt_Elec"] = new TH1F("genPt_Elec","genPt_Elec",400,0,200);
  histo1D["recoPt_Elec"] = new TH1F("recoPt_Elec","recoPt_Elec",400,0,200);

  histo1D["StCosTheta_BeforeEvtSel"] = new TH1F("StCosTheta_BeforeEvtSel","StCosTheta_BeforeEvtSel",200,-1,1);
  histo1D["StCosThetaNoBTag"] = new TH1F("StCosThetaNoBTag","StCosThetaNoBTag",200,-1,1);
  histo1D["StCosThetaLCSV"] = new TH1F("StCosThetaLCSV","StCosThetaLCSV",200,-1,1);
  histo1D["StCosThetaAllLCSV"] = new TH1F("StCosThetaAllLCSV","StCosThetaAllLCSV",200,-1,1);
  histo1D["StCosThetaMCSV"] = new TH1F("StCosThetaMCSV","StCosThetaMCSV",200,-1,1);
  histo1D["StCosThetaTCSV"] = new TH1F("StCosThetaTCSV","StCosThetaTCSV",200,-1,1);
  histo1D["JetTypeLargeLCSVEvents"] = new TH1F("JetTypeLargeLCSVEvents","JetTypeLargeLCSVEvents",51,-25.5,25.5);
  histo1D["JetTypeLargeLCSVLeadingPtEvents"] = new TH1F("JetTypeLargeLCSVLeadingPtEvents","JetTypeLargeLCSVLeadingPtEvents",51,-25.5,25.5);
  histo1D["JetTypeLCSVLightJetsLeadingPt"] = new TH1F("JetTypeLCSVLightJetsLeadingPt","JetTypeLCSVLightJetsLeadingPt",51,-25.5,25.5);
  histo1D["JetTypeLargeMCSVEvents"] = new TH1F("JetTypeLargeMCSVEvents","JetTypeLargeMCSVEvents",51,-25.5,25.5);
  histo1D["JetTypeLargeMCSVLeadingPtEvents"] = new TH1F("JetTypeLargeMCSVLeadingPtEvents","JetTypeLargeMCSVLeadingPtEvents",51,-25.5,25.5);
  histo1D["JetTypeLargeTCSVEvents"] = new TH1F("JetTypeLargeTCSVEvents","JetTypeLargeTCSVEvents",51,-25.5,25.5);
  histo1D["JetTypeLargeTCSVLeadingPtEvents"] = new TH1F("JetTypeLargeTCSVLeadingPtEvents","JetTypeLargeTCSVLeadingPtEvents",51,-25.5,25.5);
  histo1D["JetTypeLCSV"] = new TH1F("JetTypeLCSV","JetTypeLCSV",51,-25.5,25.5);
  histo1D["JetTypeMCSV"] = new TH1F("JetTypeMCSV","JetTypeMCSV",51,-25.5,25.5);
  histo1D["JetTypeTCSV"] = new TH1F("JetTypeTCSV","JetTypeTCSV",51,-25.5,25.5);
  histo1D["JetTypeLCSVLightJets"] = new TH1F("JetTypeLCSVLightJets","JetTypeLCSVLightJets",51,-25.5,25.5);
  histo1D["CSVDiscrLCSVLightJets"] = new TH1F("CSVDiscrLCSVLightJets","CSVDiscrLCSVLightJets",400,-2.5,1.5);
  histo1D["CSVDiscrLCSVLightJetsLeadingPt"] = new TH1F("CSVDiscrLCSVLightJetsLeadingPt","CSVDiscrLCSVLightJetsLeadingPt", 400, -2.5, 1.5);
  histo1D["CorrectBLeptCSVDiscr"] = new TH1F("CorrectBLeptCSVDiscr","CorrectBLeptCSVDiscr",400,-2.5,1.5);
  histo1D["CorrectBHadrCSVDiscr"] = new TH1F("CorrectBHadrCSVDiscr","CorrectBHadrCSVDiscr",400,-2.5,1.5);
  histo1D["CorrectQuark1CSVDiscr"] = new TH1F("CorrectQuark1CSVDiscr","CorrectQuark1CSVDiscr",400,-2.5,1.5);
  histo1D["CorrectQuark2CSVDiscr"] = new TH1F("CorrectQuark2CSVDiscr","CorrectQuark2CSVDiscr",400,-2.5,1.5);

  histo1D["Quark1JetNumber"] = new TH1F("Quark1JetNumber","Quark1JetNumber",12,-1.5,10.5);
  histo1D["Quark2JetNumber"] = new TH1F("Quark2JetNumber","Quark2JetNumber",12,-1.5,10.5);

  histo1D["CosThetaReco"] = new TH1F("CosThetaReco","CosThetaReco",200,-1,1);
  histo1D["NeutrinoEta"] = new TH1F("NeutrinoEta","NeutrinoEta",200,-8,8);

  histo1D["lumiWeights"] = new TH1F("lumiWeights","lumiWeights", 200,-100,100);

  //Mlb and Mqqb information:
  histo2D["MlbMqqbCorrectChosen"] = new TH2F("MlbMqqbCorrectChosen","MlbMqqbCorrectChosen",200,0,500,200,0,300);
  histo2D["MlbMqqbCorrectAll"] = new TH2F("MlbMqqbCorrectAll","MlbMqqbCorrectAll",200,0,500,200,0,300);
  histo2D["MlbMqqbWrongOne"] = new TH2F("MlbMqqbWrongOne","MlbMqqbWrongOne",200,0,500,200,0,300);
  histo2D["MlbMqqbWrongTwo"] = new TH2F("MlbMqqbWrongTwo","MlbMqqbWrongTwo",200,0,500,200,0,300);

  ////////////////////////////////////
  /// MultiSamplePlot
  ////////////////////////////////////
  map<string,MultiSamplePlot*> MSPlot;
  MSPlot["InitJets_pT_jet1_beforeEvtSel"] = new MultiSamplePlot(datasets, "InitJets_pT_jet1_beforeEvtSel", 60, -150, 650, "p_{T} (GeV)");
  MSPlot["InitJets_pT_jet2_beforeEvtSel"] = new MultiSamplePlot(datasets, "InitJets_pT_jet2_beforeEvtSel", 60, -150, 650, "p_{T} (GeV)");
  MSPlot["InitJets_pT_jet3_beforeEvtSel"] = new MultiSamplePlot(datasets, "InitJets_pT_jet3_beforeEvtSel", 60, -150, 650, "p_{T} (GeV)");
  MSPlot["InitJets_pT_jet4_beforeEvtSel"] = new MultiSamplePlot(datasets, "InitJets_pT_jet4_beforeEvtSel", 60, -150, 650, "p_{T} (GeV)");

  MSPlot["InitJets_METPt"] = new MultiSamplePlot(datasets, "InitJets_METPt", 60,0,300,"p_{T} (GeV) ");
  MSPlot["InitJets_Pt_jet1"] = new MultiSamplePlot(datasets, "InitJets_Pt_jet1",60,0,400,"p_{T} (GeV)");
  MSPlot["InitJets_Pt_jet2"] = new MultiSamplePlot(datasets, "InitJets_Pt_jet2",60,0,400,"p_{T} (GeV)");
  MSPlot["InitJets_Pt_jet3"] = new MultiSamplePlot(datasets, "InitJets_Pt_jet3",60,0,400,"p_{T} (GeV)");
  MSPlot["InitJets_Pt_jet4"] = new MultiSamplePlot(datasets, "InitJets_Pt_jet4",60,0,400,"p_{T} (GeV)");
  MSPlot["InitJets_METPt_METTypeOneCorrected"] = new MultiSamplePlot(datasets, "InitJets_METPt_METTypeOneCorrected", 60,0,300, "p_{T} (GeV)");
  MSPlot["InitJets_METPt_JerSmearingApplied"] = new MultiSamplePlot(datasets, "InitJets_METPt_JerSmearingApplied", 60,0,300, "p_{T} (GeV)");

  string leptFlavs[3]={"_other","_mu","_el"};
  for(int ii = 0; ii < 3; ii++){
      string leptFlav = leptFlavs[ii];
      MSPlot["Selected_Events_pT_jet1"+leptFlav] = new MultiSamplePlot(datasets, "Selected_Events_pT_jet1"+leptFlav, 60, -150, 650, "p_{T} (GeV)");
      MSPlot["Selected_Events_pT_jet2"+leptFlav] = new MultiSamplePlot(datasets, "Selected_Events_pT_jet2"+leptFlav, 60, -150, 650, "p_{T} (GeV)");
      MSPlot["Selected_Events_pT_jet3"+leptFlav] = new MultiSamplePlot(datasets, "Selected_Events_pT_jet3"+leptFlav, 60, -150, 650, "p_{T} (GeV)");
      MSPlot["Selected_Events_pT_jet4"+leptFlav] = new MultiSamplePlot(datasets, "Selected_Events_pT_jet4"+leptFlav, 60, -150, 650, "p_{T} (GeV)");
      MSPlot["Selected_Events_pT_4leadingjets"+leptFlav] = new MultiSamplePlot(datasets, "Selected_Events_pT_4leadingjets"+leptFlav,60, 0, 600, "p_{T} (GeV)");
      MSPlot["Selected_Events_pT_alljets"+leptFlav] = new MultiSamplePlot(datasets, "Selected_Events_pT_alljets"+leptFlav, 60, 0, 600, "p_{T} (GeV)");
      MSPlot["Selected_Events_pT_lepton"+leptFlav] = new MultiSamplePlot(datasets, "Selected_Events_pT_lepton"+leptFlav,150,-100,350,"p_{t} (GeV)");
          
      MSPlot["nSelectedJets_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "nSelectedJets_BeforeBTag"+leptFlav,14, -3.5, 10.5, "# selected jets");
      MSPlot["nSelectedJets_AfterBTag"+leptFlav] = new MultiSamplePlot(datasets, "nSelectedJets_AfterBTag"+leptFlav,14, -3.5, 10.5, "# selected jets");
      MSPlot["nBTaggedJets_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "nBTaggedJets_BeforeBTag"+leptFlav,14, -3.5, 10.5, "# b-tagged jets");
      MSPlot["nBTaggedJets_AfterBTag"+leptFlav] = new MultiSamplePlot(datasets, "nBTaggedJets_AfterBTag"+leptFlav,14, -3.5, 10.5, "# b-tagged jets");
      MSPlot["nLightJets_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "nLightJets_BeforeBTag"+leptFlav,14, -3.5, 10.5, "# light jets");
      MSPlot["nLightJets_AfterBTag"+leptFlav] = new MultiSamplePlot(datasets, "nLightJets_AfterBTag"+leptFlav,14, -3.5, 10.5, "# light jets");
  }

  /////////////////////////////
  /// ResolutionFit Stuff
  /////////////////////////////

  ResolutionFit *resFitLightJets = 0, *resFitBJets = 0, *resFitMuon = 0, *resFitElectron = 0, *resFitNeutrinoMu = 0, *resFitNeutrinoEl = 0;
    
  resFitLightJets = new ResolutionFit("LightJet");
  resFitBJets = new ResolutionFit("BJet");
  resFitMuon = new ResolutionFit("Muon");
  resFitElectron = new ResolutionFit("Electron");
  resFitNeutrinoMu = new ResolutionFit("NeutrinoMu");
  resFitNeutrinoEl = new ResolutionFit("NeutrinoEl");

  if(!CalculateResolutions && !CalculateTF){
    resFitLightJets->LoadResolutions("PersonalClasses/resolutions/lightJetReso.root");
    resFitBJets->LoadResolutions("PersonalClasses/resolutions/bJetReso.root");
    resFitMuon->LoadResolutions("PersonalClasses/resolutions/muonReso.root");
    resFitNeutrinoMu->LoadResolutions("PersonalClasses/resolutions/neutrinoSemiMuReso.root");  //Once resolutions are newly created they will be split up for SemiMu and SemiEl for Neutrino !!
    resFitNeutrinoEl->LoadResolutions("PersonalClasses/resolutions/neutrinoSemiElReso.root");
    resFitElectron->LoadResolutions("PersonalClasses/resolutions/electronReso.root");
    std::cout << " Resolutions loaded " << std::endl;
  }

  if (verbose > 0)
    cout << " - ResolutionFit instantiated ..." << endl;  

  ////////////////////////////////////
  /// Selection table
  ////////////////////////////////////
      
  vector<string> CutsSelecTableSemiMu;
  CutsSelecTableSemiMu.push_back(string("Preselected"));
  CutsSelecTableSemiMu.push_back(string("Trigged"));
  CutsSelecTableSemiMu.push_back(string("Good PV"));
  CutsSelecTableSemiMu.push_back(string("1 selected muon"));
  CutsSelecTableSemiMu.push_back(string("Veto 2nd muon"));
  CutsSelecTableSemiMu.push_back(string("Veto electron"));
  
  vector<string> CutsSelecTableSemiEl;
  CutsSelecTableSemiEl.push_back(string("Preselected"));
  CutsSelecTableSemiEl.push_back(string("Trigged"));
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

  LumiReWeighting LumiWeights, LumiWeightsUp, LumiWeightsDown;

  LumiWeights = LumiReWeighting("PersonalClasses/Calibrations/PUReweighting/pileup_MC_Summer12_S10.root", "PersonalClasses/Calibrations/PUReweighting/pileup_2012Data53X_UpToRun208357/nominal.root", "pileup", "pileup");
  LumiWeightsUp = LumiReWeighting("PersonalClasses/Calibrations/PUReweighting/pileup_MC_Summer12_S10.root", "PersonalClasses/Calibrations/PUReweighting/pileup_2012Data53X_UpToRun208357/sys_up.root", "pileup", "pileup");
  LumiWeightsDown = LumiReWeighting("PersonalClasses/Calibrations/PUReweighting/pileup_MC_Summer12_S10.root", "PersonalClasses/Calibrations/PUReweighting/pileup_2012Data53X_UpToRun208357/sys_down.root", "pileup", "pileup");
  cout << " Initialized LumiReWeighting stuff" << endl;

  // initialize lepton SF
  LeptonTools* leptonTools = new LeptonTools(false);
  leptonTools->readMuonSF("PersonalClasses/Calibrations/LeptonSF/Muon_ID_iso_Efficiencies_Run_2012ABCD_53X.root", "PersonalClasses/Calibrations/LeptonSF/MuonEfficiencies_Run_2012A_2012B_53X.root", "PersonalClasses/Calibrations/LeptonSF/MuonEfficiencies_Run_2012C_53X.root", "PersonalClasses/Calibrations/LeptonSF/TriggerMuonEfficiencies_Run_2012D_53X.root");
  leptonTools->readElectronSF();

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
    /// Initialize JEC factors            --> Updated on 5/08/2014
    /////////////////////////////////////
    vector<JetCorrectorParameters> vCorrParam;

    if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) // Data!
      {
	JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("PersonalClasses/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L1FastJet_AK5PFchs.txt");
	vCorrParam.push_back(*L1JetCorPar);
	JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("PersonalClasses/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L2Relative_AK5PFchs.txt");
	vCorrParam.push_back(*L2JetCorPar);
	JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("PersonalClasses/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L3Absolute_AK5PFchs.txt");
	vCorrParam.push_back(*L3JetCorPar);
	JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("PersonalClasses/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L2L3Residual_AK5PFchs.txt");
	vCorrParam.push_back(*L2L3ResJetCorPar);
      }
    else
      {
	JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("PersonalClasses/Calibrations/JECFiles/START53_V23_Summer13_L1FastJet_AK5PFchs.txt");
	vCorrParam.push_back(*L1JetCorPar);
	JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("PersonalClasses/Calibrations/JECFiles/START53_V23_Summer13_L2Relative_AK5PFchs.txt");
	vCorrParam.push_back(*L2JetCorPar);
	JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("PersonalClasses/Calibrations/JECFiles/START53_V23_Summer13_L3Absolute_AK5PFchs.txt");
	vCorrParam.push_back(*L3JetCorPar);
      }
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("PersonalClasses/Calibrations/JECFiles/START53_V23_Summer13_Uncertainty_AK5PFchs.txt");
    //JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(*(new JetCorrectorParameters("JECFiles/Fall12_V7_DATA_UncertaintySources_AK5PFchs.txt", "SubTotalMC")));
    //JetCorrectionUncertainty *jecUncTotal = new JetCorrectionUncertainty(*(new JetCorrectorParameters("JECFiles/Fall12_V7_DATA_UncertaintySources_AK5PFchs.txt", "Total")));

    //True means also redo L1        
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);
   
    /////////////////////////////////////////
    //  LHCO Output files + GeneratorInfo  //
    /////////////////////////////////////////
    ofstream EventInfoFile;
    ofstream outFileReco[16];
    LHCOOutput lhcoOutput; //Initialize class
    ofstream outFile[4];
    if(GenLHCOOutput == true){	
      EventInfoFile.open("MadWeightInput/AnalyzerOutput/EventNumberInformation.lhco");
      EventInfoFile << " Event Number  MuPos  MuNeg   ElPos   ElNeg  ChannelNumber  FalseEventContent       selectedEvent      selectedChannelNumber " << endl;

      outFile[0].open("MadWeightInput/AnalyzerOutput/TTbarLHCO_PositiveMuon.lhco");
      outFile[1].open("MadWeightInput/AnalyzerOutput/TTbarLHCO_NegativeMuon.lhco");
      outFile[2].open("MadWeightInput/AnalyzerOutput/TTbarLHCO_PositiveElectron.lhco");
      outFile[3].open("MadWeightInput/AnalyzerOutput/TTbarLHCO_NegativeElectron.lhco");
    }
    unsigned int NumberPosRecoMu = 0, NumberNegRecoMu = 0, NumberPosRecoEl = 0, NumberNegRecoEl = 0;
    
    unsigned int NumberCorrectEvents = 0; //Counts the number of semi-leptonic events
    unsigned int NumberNegativeElectrons = 0, NumberNegativeMuons = 0, NumberPositiveElectrons = 0, NumberPositiveMuons = 0;
    int EventContent[5]; //0:top; 1:b; 2: u,c,d,s; 3:W; 4:mu + neutrino
    
    //Cos Theta information
    TLorentzVector *sTop, *WLeptTRF;       //Pointers allowed .. ?
    TLorentzVector leptonWRF;              //Don't use pointers here to avoid overwriting the original TLorentzVectors!
    float standardCosTheta = 0;

    bool FalseEventContent = false;
    cout << " FalseEventContent : " << FalseEventContent << endl;
    TRootMCParticle *Top,*TopBar,*Bottom, *BottomBar,*Lepton,*NeutrinoMC,*WPlus,*WMinus,*Light,*LightBar;

    /////////////////////
    //  Used classes   //
    /////////////////////  
    BTagStudy bTagStudy;  //--> Should only be called before the event loop (otherwise the counters will not give the correct result)
    MlbStudy mlbStudy(NrConsideredBTagOptions);
    //TFCreation tfCreation;
    //tfCreation.InitializeVariables();       //Should be called since constructor should work for both analyzers!
    TFnTuple* tfNTuple = 0;

    //Initialize TFnTuple specific stuff:
    TTree* TFTree = new TTree("TFTree","Tree containing the Transfer Function information");
    TFTree->Branch("TheTFTree","TFnTuple",&tfNTuple);
    TFile* TFTreeFile = new TFile("TFInformation/TransferFunctionTree.root","RECREATE");

    ////////////////////////////////////
    //	loop on events
    ////////////////////////////////////
    
    nEvents[d] = 0;
    int itriggerSemiMu = -1,itriggerSemiEl = -1, previousRun = -1;
    if (verbose > 1)
      cout << "	Loop over events " << endl;

    for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++){
    //for (unsigned int ievt = 0; ievt < 500000; ievt++){

      //if(ievt > 200000) GenLHCOOutput = false;	

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
      
      //TRootEvent* event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets_corrected, mets);  //Use uncorrected jets ...
      TRootEvent* event = treeLoader.LoadEvent(ievt, vertex, init_muons, init_electrons, init_jets, mets);
      
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
	TRootGenEvent* genEvt = treeLoader.LoadGenEvent(ievt,false);
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
      if(dataSetName.find("TTbarJets") == 0){
	//cout << "LOADING GenEvent" << endl;
	TRootGenEvent* genEvt = treeLoader.LoadGenEvent(ievt,false);
	if( genEvt->isSemiLeptonic() )
	scaleFactor *= (0.108*9.)*(0.676*1.5);
	else if( genEvt->isFullHadronic() )
	scaleFactor *= (0.676*1.5)*(0.676*1.5);
	else if( genEvt->isFullLeptonic() )
	scaleFactor *= (0.108*9.)*(0.108*9.);
      }
      MSPlot["InitJets_METPt"]->Fill(mets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      if(init_jets.size() > 0) MSPlot["InitJets_Pt_jet1"]->Fill(init_jets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      if(init_jets.size() > 1) MSPlot["InitJets_Pt_jet2"]->Fill(init_jets[1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      if(init_jets.size() > 2) MSPlot["InitJets_Pt_jet3"]->Fill(init_jets[2]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      if(init_jets.size() > 3) MSPlot["InitJets_Pt_jet4"]->Fill(init_jets[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      
      //////////////////////////////////////
      // Apply Jet Corrections on-the-fly //   
      //////////////////////////////////////
      if( dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ){
	  jetTools->unCorrectMETTypeOne(init_jets, mets[0], true);
	  jetTools->correctJets(init_jets, event->kt6PFJets_rho(), true);
	  jetTools->correctMETTypeOne(init_jets, mets[0], true);
      }
      else{
	  jetTools->unCorrectMETTypeOne(init_jets, mets[0], false);
	  jetTools->correctJets(init_jets, event->kt6PFJets_rho(), false);
	  jetTools->correctMETTypeOne(init_jets, mets[0], false);
      }
      MSPlot["InitJets_METPt_METTypeOneCorrected"]->Fill(mets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);

      ////////////////////////////////////////
      //  Beam scraping and PU reweighting
      ////////////////////////////////////////
      double lumiWeight = 1;
      if(! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") ){   
	 if(doLumiWeightShift == "Nominal")    lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() );
         else if(doLumiWeightShift == "Plus")  lumiWeight = LumiWeightsUp.ITweight( (int)event->nTruePU() );
         else if(doLumiWeightShift == "Minus") lumiWeight = LumiWeightsDown.ITweight( (int)event->nTruePU() );
      }
      histo1D["lumiWeights"]->Fill(lumiWeight);	

      // up syst -> lumiWeight = LumiWeightsUp.ITweight( (int)event->nTruePU() );
      // down syst -> lumiWeight = LumiWeightsDown.ITweight( (int)event->nTruePU() );

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
	  if( event->runId() >= 190456 && event->runId() <= 190738 )
	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v11"), currentRun, iFile);
	  else if( event->runId() >= 190782 && event->runId() <= 193621)
	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v12"), currentRun, iFile);
	  else if(event->runId() >= 193834  && event->runId() <= 196531 )
	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);
	  else if( event->runId() >= 198022  && event->runId() <= 199608)
	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v14"), currentRun, iFile);
	  else if( event->runId() >= 199698 && event->runId() <= 209151)
	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v15"), currentRun, iFile);
	  else
	    cout << "Unknown run for SemiMu HLTpath selection: " << event->runId() << endl;
	  if( itriggerSemiMu == 9999 ){
	    cout << "itriggerSemiMu == 9999 for SemiMu HLTpath selection: " << event->runId() << endl;
	    exit(-1);
	  }
        }
        else{ 
	  itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun); // Summer12 DR53X
	}
	
	// semi-electron
        if(dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0 ) {
	  
	  // 2.7/fb recalib 
	  if( event->runId() <= 190738 )
	    itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v8"), currentRun, iFile);
	  else if( event->runId() >= 190782 && event->runId() <= 191411 )
	    itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v9"), currentRun, iFile);
	  else if( event->runId() >= 191695 && event->runId() <= 196531)
	    itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v10"), currentRun, iFile);
	  else if( event->runId() >= 198049 && event->runId() <= 208686)
	    itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v11"), currentRun, iFile);
	  //else if( event->runId() > 208686)
	  //     itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v11"), currentRun, iFile);
	  else { 
	    cout << "Unknown run for SemiEl HLTpath selection: " << event->runId() << endl;
	    exit(1);
	  }
	  if( itriggerSemiEl == 9999 ){	    
	    cout << "itriggerSemiEl == 9999 for SemiEl HLTpath selection: " << event->runId() << endl;
	    exit(-1);
	  }
        }
        else	  
	  itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v10"), currentRun); // Summer12 DR53X 
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
	
	jetTools->correctJetJER(init_jets, genjets, mets[0], doJERShift, false); //false means don't use old numbers but newer ones ...
	if (doJESShift != "nominal")
	  jetTools->correctJetJESUnc(init_jets, mets[0], doJESShift, 1);  //1 = nSigma
      }
      MSPlot["InitJets_METPt_JerSmearingApplied"]->Fill(mets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);

      ////////////////////////////////////////////////////////
      // Access particle information before event selection //
      // Write this information to LHCO Output for MW       //
      ////////////////////////////////////////////////////////
      for(int ll = 0;ll<5;ll++){EventContent[ll]=0;}
      
      //Loop over all the mcParticles
      for(unsigned int i=0; i<mcParticles.size(); i++){
	if( mcParticles[i]->status() != 3) continue;
	
	int partType=mcParticles[i]->type(); if(verbosity>4)cout<<"-->Type of mcParticle : "<<partType<<endl;
        //if( fabs(partType) == 4) std::cout << " Mass of c-quark : " << mcParticles[i]->M() << " with mother : " << mcParticles[i]->motherType() << std::endl; --> Almost lways 1.5
	
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
      if(GenLHCOOutput == true) EventInfoFile << "     " << ievt << "         ";
      if(EventContent[0]==2 && EventContent[1]==2 && EventContent[2]==2 && EventContent[3]==2 && EventContent[4]==2){
	FalseEventContent = false;
	vector<TRootMCParticle*> LHCOVector(6);
	vector<int> MadGraphId(6,4);
        vector<float> MGBtag(6,0.0);
	
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
          MGBtag[0] = 2.0; MGBtag[3] = 2.0;
	  if(Lepton->type() == 11){           //Looking at negative electron events (index 3 for LHCO file)
	    MadGraphId[4] = 1; //MadGraph Id of e = 1
	    MadGraphId[5] = 6; //MadGraph Id of MET = 6
	    NumberNegativeElectrons++;
	    if(GenLHCOOutput == true){
	      lhcoOutput.LHCOEventOutput(3, outFile[3], NumberNegativeElectrons,LHCOVector,MadGraphId, MGBtag);
	      EventInfoFile << "  0      0       0       1        " << NumberNegativeElectrons << "     ";
	    }
	  }//Negative electron
	  else if(Lepton->type() == 13){       //Looking at negative muon events (index 1 for LHCO file)
	    MadGraphId[4] = 2; //MadGraphId of mu = 2
	    MadGraphId[5] = 6; 
	    NumberNegativeMuons++;
	    if(GenLHCOOutput == true){
	      lhcoOutput.LHCOEventOutput(1, outFile[1], NumberNegativeMuons,LHCOVector,MadGraphId, MGBtag);
	      EventInfoFile << "  0      1       0       0        " << NumberNegativeMuons << "     ";
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
          MGBtag[0] = 2.0; MGBtag[3] = 2.0;
	  if(Lepton->type() == -11){            //Looking at positive electron events (index 2 for LHCO file)
	    MadGraphId[1] = 1; //MadGraphId of electron = 1
	    MadGraphId[2] = 6; 
	    NumberPositiveElectrons++;
	    if(GenLHCOOutput == true){
	      lhcoOutput.LHCOEventOutput(2, outFile[2], NumberPositiveElectrons,LHCOVector,MadGraphId, MGBtag);
	      EventInfoFile << "  0      0       1       0        " << NumberPositiveElectrons << "     ";
	    }
	  }//Positive electron
	  else if(Lepton->type() == -13){             //Looking at positive muon events (index 0 for LHCO file)
	    MadGraphId[1] = 2; //MadGraphId of muon = 2
	    MadGraphId[2] = 6; 
	    NumberPositiveMuons++;
	    if(GenLHCOOutput == true){
	      lhcoOutput.LHCOEventOutput(0, outFile[0], NumberPositiveMuons,LHCOVector,MadGraphId, MGBtag);
	      EventInfoFile << "  1      0       0       0        " << NumberPositiveMuons << "     ";
	    }
	  }//Positive muon
	  
	  if(verbosity>3){
	    cout << " WPlus information : "<<WPlus->Px()<< ", "<<WPlus->Py()<<", "<< WPlus->Pz()<<", "<<WPlus->E()<< endl;
	    cout << " Top information : "<<Top->Px()<< ", "<<Top->Py()<<", "<< Top->Pz()<<", "<<Top->E()<< endl;
	  }
	  WLeptTRF = (TLorentzVector*) WPlus;
	  sTop = (TLorentzVector*) Top;			
	}//Positive lepton
	if(GenLHCOOutput == true) EventInfoFile << "        1               ";
	
	//////////////////////////////////////
	//  Look at cos theta distribution  //
	////////////////////////////////////// 

	//-----    Applying boost on muon and W    -----//
	leptonWRF = *Lepton;
	leptonWRF.Boost(-WLeptTRF->BoostVector());
	WLeptTRF->Boost(-sTop->BoostVector());
	if(verbosity>3){
	  cout<<" leptonWRF information : "<<leptonWRF.Px()<<", "<<leptonWRF.Py()<<", "<<leptonWRF.Pz()<<", "<<leptonWRF.E()<<endl;
	}
	
	//-----   Calculating cos theta:   -----
	standardCosTheta = ((WLeptTRF->Vect()).Dot(leptonWRF.Vect()))/(((WLeptTRF->Vect()).Mag())*((leptonWRF.Vect()).Mag()));
	if(verbosity>4) cout << " cos theta (gen): " << standardCosTheta << endl << endl;
	histo1D["StCosTheta_BeforeEvtSel"]->Fill(standardCosTheta);
	
      }//Correct event content found
      else{
	FalseEventContent = true;
	if(GenLHCOOutput == true) EventInfoFile << "                                                 0                ";
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

      if(init_jets.size() >=4){ // MSPlots before 'basic' event selection (no b-tag)
	MSPlot["InitJets_pT_jet1_beforeEvtSel"]->Fill(init_jets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
	MSPlot["InitJets_pT_jet2_beforeEvtSel"]->Fill(init_jets[1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
	MSPlot["InitJets_pT_jet3_beforeEvtSel"]->Fill(init_jets[2]->Pt(), datasets[d], true, Luminosity*scaleFactor);
	MSPlot["InitJets_pT_jet4_beforeEvtSel"]->Fill(init_jets[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      }
      
      //Declare selection instance    
      Selection selection(init_jets, init_muons, init_electrons, mets, event->kt6PFJets_rho());
      selection.setJetCuts(30,2.5,0.01,1.,0.98,0.3,0.1);     
      selection.setMuonCuts(26,2.1,0.12,0.2,0.3,1,0.5,5,0); 
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

      vector<TRootJet*> selectedJets;
      vector<TRootMuon*> selectedMuons;
      vector<TRootElectron*> selectedElectrons;
      vector<TRootMuon*> vetoMuons = selection.GetSelectedLooseMuons(10,2.5,0.2);
      vector<TRootElectron*> vetoElectronsSemiMu = selection.GetSelectedLooseElectrons(20,2.5,0.15);
      vector<TRootElectron*> vetoElectronsSemiEl = selection.GetSelectedLooseElectrons(20,2.5,0.15);

      selectedJets = selection.GetSelectedJets(true);
      selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);
      selectedElectrons = selection.GetSelectedElectrons(selectedJets);
      
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
      	    
      float MassMlb = 103.286;
      float MassMqqb = 178.722;
      float SigmaMlb = 26.7764;
      float SigmaMqqb = 18.1385;
      
      /////////////////////////////
      // Neutrino Reconstruction //
      /////////////////////////////
      float NeutrinoPx, NeutrinoPy;
      float NeutrinoPz=99.;//with this value it can be distinguished in plot!
      TLorentzVector Neutrino;
      
      //////////////////////
      // Event selection  //
      //////////////////////
      bool eventselectedSemiMu = false;
      bool eventselectedSemiEl = false;
      enum DecayChannel_t {semiMu, semiEl};
      DecayChannel_t decayChannel;

      if (dataSetName != "Data"&&  selectedElectrons.size() ==1 ) {
	scaleFactor = scaleFactor*leptonTools->getElectronSF(selectedElectrons[0]->Eta(), selectedElectrons[0]->Pt(), doLeptonSFShift );
	//histo1D["leptonScales"]->Fill(leptonTools->getElectronSF(selectedElectrons[0]->Eta(), selectedElectrons[0]->Pt(), doLeptonSFShift));
      }
     
      // semi-mu selection
      selecTableSemiMu.Fill(d,0,scaleFactor*lumiWeight);              
      if (triggedSemiMu) {
	selecTableSemiMu.Fill(d,1,scaleFactor*lumiWeight);
	if (isGoodPV) {
	  selecTableSemiMu.Fill(d,2,scaleFactor*lumiWeight);
	  if (selectedMuons.size() == 1) {
	    selecTableSemiMu.Fill(d,3,scaleFactor*lumiWeight);
	    if( vetoMuons.size() == 1 ) {
	      selecTableSemiMu.Fill(d,4,scaleFactor*lumiWeight);
	      if (vetoElectronsSemiMu.size() == 0) {
		selecTableSemiMu.Fill(d,5,scaleFactor*lumiWeight);
		if (selectedJets.size() >= 1) {
		  selecTableSemiMu.Fill(d,6,scaleFactor*lumiWeight);
		  if (selectedJets.size() >= 2) {
		    selecTableSemiMu.Fill(d,7,scaleFactor*lumiWeight);
		    if (selectedJets.size() >= 3) {
		      selecTableSemiMu.Fill(d,8,scaleFactor*lumiWeight);
		      if (selectedJets.size() >= 4) {
			selecTableSemiMu.Fill(d,9,scaleFactor*lumiWeight);
		 	eventselectedSemiMu = true;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }

      // semi-elec selection
      selecTableSemiEl.Fill(d,0,scaleFactor*lumiWeight);
      if( triggedSemiEl) {
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
 
      if( !eventselectedSemiMu && !eventselectedSemiEl && RecoLHCOOutput == true ) EventInfoFile << "     Evt sel failed   "<<endl;
      if (!eventselectedSemiMu && !eventselectedSemiEl) continue;
      if(RecoLHCOOutput == true) EventInfoFile << "             1          ";  
      
      //Counting the number of events passing through the 'basic' event selection requirements    
      if (eventselectedSemiMu){ nSelectedMu++; decayChannel = semiMu;}
      if (eventselectedSemiEl){ nSelectedEl++; decayChannel = semiEl;}

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
      MSPlot["Selected_Events_pT_jet1"+leptonFlav]->Fill(selectedJets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["Selected_Events_pT_jet2"+leptonFlav]->Fill(selectedJets[1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["Selected_Events_pT_jet3"+leptonFlav]->Fill(selectedJets[2]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["Selected_Events_pT_jet4"+leptonFlav]->Fill(selectedJets[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["Selected_Events_pT_lepton"+leptonFlav]->Fill(selectedLepton->Pt(), datasets[d], true, Luminosity*scaleFactor);
      
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
	if(hadronicWJet1_.first != 9999) histo1D["Quark1JetNumber"]->Fill(hadronicWJet1_.first);
	else 				 histo1D["Quark1JetNumber"]->Fill(-1);
	if(hadronicWJet2_.first != 9999) histo1D["Quark2JetNumber"]->Fill(hadronicWJet2_.first);
	else				 histo1D["Quark2JetNumber"]->Fill(-1);

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

	//Change order to match with class! (0 = BLeptonic, 1 = BHadronic, 2 = Quark1 & 3 =  Quark2!!)
	jetCombi.push_back(leptonicBJet_.first);
	jetCombi.push_back(hadronicBJet_.first);
	jetCombi.push_back(hadronicWJet1_.first);
	jetCombi.push_back(hadronicWJet2_.first);

	CorrectBLeptonic = jetCombi[0];
	CorrectBHadronic = jetCombi[1];
	CorrectQuark1 = jetCombi[2];
	CorrectQuark2 = jetCombi[3];	

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

	if(CorrectBLeptonic != 9999) histo1D["CorrectBLeptCSVDiscr"]->Fill(selectedJets[CorrectBLeptonic]->btag_combinedSecondaryVertexBJetTags());
	else histo1D["CorrectBLeptCSVDiscr"]->Fill(-2);
        if(CorrectBHadronic != 9999) histo1D["CorrectBHadrCSVDiscr"]->Fill(selectedJets[CorrectBHadronic]->btag_combinedSecondaryVertexBJetTags());
	else histo1D["CorrectBHadrCSVDiscr"]->Fill(-2);
	if(CorrectQuark1 != 9999) histo1D["CorrectQuark1CSVDiscr"]->Fill(selectedJets[CorrectQuark1]->btag_combinedSecondaryVertexBJetTags());
	else histo1D["CorrectQuark1CSVDiscr"]->Fill(-2);
	if(CorrectQuark2 != 9999) histo1D["CorrectQuark2CSVDiscr"]->Fill(selectedJets[CorrectQuark2]->btag_combinedSecondaryVertexBJetTags());
	else histo1D["CorrectQuark2CSVDiscr"]->Fill(-2);
 	
	//Working on generator level (i.e. jets level):  
	if(jetCombi[0]!=9999 && jetCombi[1]!=9999 && jetCombi[2]!=9999 && jetCombi[3]!=9999){    
	  CorrectRecMassW=(*selectedJets[jetCombi[2]]+*selectedJets[jetCombi[3]]).M();
	  CorrectRecMassTop=(*selectedJets[jetCombi[2]]+*selectedJets[jetCombi[3]]+*selectedJets[jetCombi[1]]).M();
	  CorrectRecMassMlb=(*selectedJets[CorrectBLeptonic]+*selectedLepton).M();
	  CorrectRecMassMqqb=(*selectedJets[CorrectQuark1]+*selectedJets[CorrectQuark2]+*selectedJets[CorrectBHadronic]).M();

	  histo1D["WMass"]->Fill(CorrectRecMassW);
	  histo1D["TopMass"]->Fill(CorrectRecMassTop);
	  histo1D["MlbMass"]->Fill(CorrectRecMassMlb);
	  histo1D["MqqbMass"]->Fill(CorrectRecMassMqqb);
	}	      	      	      	       	      

	if(hadronicWJet1_.first < 9999 && hadronicWJet2_.first < 9999 && hadronicBJet_.first < 9999 && leptonicBJet_.first < 9999 ){
          if(CalculateResolutions){
	    // Fill the resolution-stuff for events where the 4 ttbar semi-lep partons are all matched to jets
	    resFitLightJets->Fill(selectedJets[hadronicWJet1_.first], mcParticles[hadronicWJet1_.second]);
	    resFitLightJets->Fill(selectedJets[hadronicWJet2_.first], mcParticles[hadronicWJet2_.second]);
	    resFitBJets->Fill(selectedJets[hadronicBJet_.first], mcParticles[hadronicBJet_.second]);
	    resFitBJets->Fill(selectedJets[leptonicBJet_.first], mcParticles[leptonicBJet_.second]);
	    if(eventselectedSemiMu == true){
	      resFitMuon->Fill(selectedMuons[0], Lepton);
	      resFitNeutrinoMu->Fill(mets[0], NeutrinoMC);
	    }
	    else if(eventselectedSemiEl == true){
	      resFitElectron->Fill(selectedElectrons[0], Lepton);
	      resFitNeutrinoEl->Fill(mets[0], NeutrinoMC);
	    }
	  }//End of calculate Resolutions
	  if(CalculateTF){
	    //tfCreation.FillHistograms( (TLorentzVector*) &mcParticlesMatching[hadronicWJet1_.second], (TLorentzVector*) &mcParticlesMatching[hadronicWJet2_.second], (TLorentzVector*) &mcParticlesMatching[hadronicBJet_.second], (TLorentzVector*) &mcParticlesMatching[leptonicBJet_.second], (TLorentzVector*) Lepton, (TLorentzVector*) selectedJets[hadronicWJet1_.first], (TLorentzVector*) selectedJets[hadronicWJet2_.first], (TLorentzVector*) selectedJets[hadronicBJet_.first], (TLorentzVector*) selectedJets[leptonicBJet_.first], (TLorentzVector*) selectedLepton, decayChannel);

	    //Check the DeltaR vlaue between the different partons and reconstructed particles!:
	    if(eventselectedSemiMu == true){
	      histo1D["genPt_Muon"]->Fill( Lepton->Pt() );
	      histo1D["recoPt_Muon"]->Fill(selectedLepton->Pt());
	    }
            if(eventselectedSemiEl == true){
	      histo1D["genPt_Elec"]->Fill( Lepton->Pt() );
	      histo1D["recoPt_Elec"]->Fill( selectedLepton->Pt());
	    }

            //Fill the Transfer Function Tree file!
            tfNTuple = new TFnTuple();
            tfNTuple->setEventID( event->eventId() );
            tfNTuple->setRecoVectorLight1( (TLorentzVector) *selectedJets[hadronicWJet1_.first] );
            tfNTuple->setRecoVectorLight2( (TLorentzVector) *selectedJets[hadronicWJet2_.first] );
            tfNTuple->setRecoVectorHadrB(  (TLorentzVector) *selectedJets[hadronicBJet_.first]  );
            tfNTuple->setRecoVectorLeptB(  (TLorentzVector) *selectedJets[leptonicBJet_.first]  );
            tfNTuple->setRecoVectorLepton( (TLorentzVector) *selectedLepton                     );
            tfNTuple->setGenVectorLight1(  (TLorentzVector) mcParticlesMatching[hadronicWJet1_.second] );
            tfNTuple->setGenVectorLight2(  (TLorentzVector) mcParticlesMatching[hadronicWJet2_.second] );
            tfNTuple->setGenVectorHadrB(   (TLorentzVector) mcParticlesMatching[hadronicBJet_.second]  );
            tfNTuple->setGenVectorLeptB(   (TLorentzVector) mcParticlesMatching[leptonicBJet_.second]  );
            tfNTuple->setGenVectorLepton(  (TLorentzVector) *Lepton                                     );

            TFTree->Fill();
            delete tfNTuple;

	  }//End of calculate Transfer Functions
        }//End of matched particles reconstructed
      }//if dataset Semi mu ttbar

      if(!CalculateBTag) continue;
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
	      histo1D["CSVDiscrLCSVLightJets"]->Fill(selectedJets[(bTagStudy.getLightJets(0))[jj]]->btag_combinedSecondaryVertexBJetTags());   
	      //--> check whether the so-called light jets don't all have discr -1 ...	

	      for(int ii = 0; ii<JetPartonPair.size(); ii++){ //Look at all the matched jets
		if((bTagStudy.getLightJets(0))[jj] == JetPartonPair[ii].first){  //Check whether the considered jet can be matched!
		  histo1D["JetTypeLCSVLightJets"]->Fill(mcParticlesMatching[JetPartonPair[ii].second].type());
		}
		else{    //Unmatched jets!
		  histo1D["JetTypeLCSVLightJets"]->Fill(25.);
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
	      histo1D["StCosThetaLCSV"]->Fill(standardCosTheta);
	      for(int ii = 0; ii<JetPartonPair.size(); ii++){ //Look at all the matched jets
		if((bTagStudy.getbTaggedJets(0))[jj] == JetPartonPair[ii].first){  //Check whether the considered jet can be matched!
		  histo1D["JetTypeLCSV"]->Fill(mcParticlesMatching[JetPartonPair[ii].second].type());
		}
		else{
		  histo1D["JetTypeLCSV"]->Fill(25.);
		}
	      }
	    }//End of Loose b-jets 

	    //Kinematic information for the Medium b-jets
	    if((bTagStudy.getbTaggedJets(1)).size() >=2){
	      histo1D["StCosThetaMCSV"]->Fill(standardCosTheta);
	      for(int ii = 0; ii<JetPartonPair.size(); ii++){ //Look at all the matched jets
		if((bTagStudy.getbTaggedJets(1))[jj] == JetPartonPair[ii].first){  //Check whether the considered jet can be matched!
		  histo1D["JetTypeMCSV"]->Fill(mcParticlesMatching[JetPartonPair[ii].second].type());
		}
		else{
		  histo1D["JetTypeMCSV"]->Fill(25.);
		}
	      }
	    }//End of Medium b-jets 

	    //Kinematic information for the Tight b-jets
	    if((bTagStudy.getbTaggedJets(3)).size() >=2){
	      histo1D["StCosThetaTCSV"]->Fill(standardCosTheta);
	      for(int ii = 0; ii<JetPartonPair.size(); ii++){ //Look at all the matched jets
		if((bTagStudy.getbTaggedJets(3))[jj] == JetPartonPair[ii].first){  //Check whether the considered jet can be matched!
		  histo1D["JetTypeTCSV"]->Fill(mcParticlesMatching[JetPartonPair[ii].second].type());
		}
		else{
		  histo1D["JetTypeTCSV"]->Fill(25.);
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
      histo1D["StCosThetaNoBTag"]->Fill(standardCosTheta); 
      //---------------------------------------------------------------------------------------------------------------------------- End of bTagStudy class stuff

      //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
      //  Event selection choice (17/06/2014)  //
      //   --> Continue with 2 T b-tags        //
      //   --> No veto on light jets!          //
      //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//

      //MSPlots before and after #b-tagged and #light jets constraints
      MSPlot["nSelectedJets_BeforeBTag"+leptonFlav]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["nBTaggedJets_BeforeBTag"+leptonFlav]->Fill( (bTagStudy.getbTaggedJets(ChosenBTagOption)).size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["nLightJets_BeforeBTag"+leptonFlav]->Fill( (bTagStudy.getLightJets(ChosenBTagOption)).size(), datasets[d], true, Luminosity*scaleFactor);
      if( NrConsideredBTagOptions == 1 && ( (bTagStudy.getbTaggedJets(ChosenBTagOption)).size() < 2 || (bTagStudy.getLightJets(ChosenBTagOption)).size() < 2 ) ){
	if(RecoLHCOOutput == true) EventInfoFile<<"    B-tag failed "<<endl;
	continue;
      }      

      MSPlot["nSelectedJets_AfterBTag"+leptonFlav]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["nBTaggedJets_AfterBTag"+leptonFlav]->Fill( (bTagStudy.getbTaggedJets(ChosenBTagOption)).size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["nLightJets_AfterBTag"+leptonFlav]->Fill( (bTagStudy.getLightJets(ChosenBTagOption)).size(), datasets[d], true, Luminosity*scaleFactor);

      // Count the number of events:
      if(eventselectedSemiMu) nSelectedMu++;
      if(eventselectedSemiEl) nSelectedEl++;
     
      ////////////////////////////////
      //  Mlb and Mqqb information  //
      ////////////////////////////////
      float MlbCorrect = 0, MqqbCorrect = 0;
      if(CorrectBLeptonic != 9999 && CorrectBHadronic != 9999 && CorrectQuark1 != 9999 && CorrectQuark2 != 9999){
	MlbCorrect = (*selectedLepton+*selectedJets[CorrectBLeptonic]).M();
	MqqbCorrect = (*selectedJets[CorrectBHadronic] + *selectedJets[CorrectQuark1] + *selectedJets[CorrectQuark2]).M();
	histo2D["MlbMqqbCorrectAll"]->Fill(MqqbCorrect,MlbCorrect);
      }

      //******************************************************************************************//
      //  Loop can be kept, since if NrConsideredBTagOptions == 1 only 1 case will be considered  //
      //    --> If this is the case, force the option to be the one chosen above!!                //
      //******************************************************************************************//
      for(int Option = 0; Option < NrConsideredBTagOptions; Option++){
	if(NrConsideredBTagOptions == 1) Option = ChosenBTagOption;    //Force Option to be equal to the one chosen!
       	
	mlbStudy.calculateChiSquared(jetCombi, bTagStudy.getbTaggedJets(Option), bTagStudy.getLightJets(Option), selectedLepton, selectedJets, MassMlb, SigmaMlb, MassMqqb, SigmaMqqb);
	mlbStudy.calculateEfficiency(Option, jetCombi, bTagStudy.getbTaggedJets(Option), bTagStudy.getLightJets(Option), NrConsideredBTagOptions, ChiSqCutValue);

	//	//TODO --> Update!!
	//	  if( (CorrectBLeptonic == (bTagStudy.getbTaggedJets(Option))[0] || CorrectBLeptonic == (bTagStudy.getbTaggedJets(Option))[1]) &
	//	      (CorrectBHadronic == (bTagStudy.getbTaggedJets(Option))[1] || CorrectBHadronic == (bTagStudy.getbTaggedJets(Option))[0]) &&
	//	      (CorrectQuark1 == (bTagStudy.getLightJets(Option))[0]      || CorrectQuark1 == (bTagStudy.getLightJets(Option))[1]) &&
	//	      (CorrectQuark2 == (bTagStudy.getLightJets(Option))[0]      || CorrectQuark2 == (bTagStudy.getLightJets(Option))[1])){ 
	//		h_MlbMqqbCorrectChosen.Fill(MqqbCorrect,MlbCorrect); 
	//	  } 
      }

      //*********************************************************//
      //   Divide events in categories N(4 jets) and N(5 jets)   //
      //   --> Used to make 2 or 3 light jets choice in evtsel   //
      //*********************************************************//
      //cout << " \n ChiSq Indices (4 jets and 4+ jets): " << mlbStudy.getLowestChiSq4JetsIndex() << " & " << mlbStudy.getLowestChiSqIndex() << std::endl;
      //cout << " Correct ChiSq index               : " << mlbStudy.getCorrectChiSq() << endl;
      //cout << "  --> Correct indices (bLept, bHadr, quark1 & quark2) : " << jetCombi[0] << " , " << jetCombi[1] << " , " << jetCombi[2] << " & " << jetCombi[3] << endl;
      //cout << "  --> Available indices (b1, b2, q1, q2 & q3) : " << (bTagStudy.getbTaggedJets(ChosenBTagOption))[0] << " , " << (bTagStudy.getbTaggedJets(ChosenBTagOption))[1] << ", " << (bTagStudy.getLightJets(ChosenBTagOption))[0] << " , " << (bTagStudy.getLightJets(ChosenBTagOption))[1];
      //if( (bTagStudy.getLightJets(ChosenBTagOption)).size() > 2) cout << " & " << (bTagStudy.getLightJets(ChosenBTagOption))[2] << endl;
      //else cout << " " << endl;

      if(NrConsideredBTagOptions == 1){  //Only look at these categories for the chosen 2T b-tag option
	if( (bTagStudy.getLightJets(ChosenBTagOption)).size() == 2){
	  //Looking at the N(2 light jets) category!
	  FourJetsCategory++;
	  if(  jetCombi[0] != 9999 && jetCombi[1] != 9999 && jetCombi[2] != 9999 && jetCombi[3] != 9999 ){  //Event has been matched!
	    FourJetsMatched++;
	    //Include the ChiSq Mlb-Mqqb information for selecting the b-jets (and the 2 light jets in the N(3 light) case)
	    if(mlbStudy.getLowestChiSq4JetsIndex() == mlbStudy.getCorrectChiSq()) FourJetsGoodCombiChosen++; //Represents 's' and 'b' is defined as FourJetsCategory (all evts) - 's'
	  }
   	}//End of N(2 light) case
	else{
	  FiveJetsCategory++;
	  if(  jetCombi[0] != 9999 && jetCombi[1] != 9999 && jetCombi[2] != 9999 && jetCombi[3] != 9999 ){  //Event has been matched!
	    FiveJetsMatched++;
	    //Include the ChiSq Mlb-Mqqb information for selecting the b-jets (and the 2 light jets in the N(3 light) case)
	    //Represents 's', and 'b' is defined as FourJetsCategory (all evts) - 's'
	    if( mlbStudy.getLowestChiSq4JetsIndex() == mlbStudy.getCorrectChiSq() ) FiveJetsAsFourGoodCombiChosen++; 
	    if( mlbStudy.getLowestChiSqIndex() == mlbStudy.getCorrectChiSq() ) FiveJetsAsFiveGoodCombiChosen++; 
	  }
	}//End of N(3 light) case 
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
      float LeptBCSV = selectedJets[(bTagStudy.getbTaggedJets(ChosenBTagOption))[mlbStudy.getChosenBLept()]]->btag_combinedSecondaryVertexBJetTags();
      float HadrBCSV = selectedJets[(bTagStudy.getbTaggedJets(ChosenBTagOption))[mlbStudy.getChosenBHadr()]]->btag_combinedSecondaryVertexBJetTags();

      /////////////////////////////////////////////
      //  Filling of LHCO files for reco events  //
      /////////////////////////////////////////////
      vector<TLorentzVector*> LHCORecoVector(6);
      vector<int> MadGraphRecoId(6,4);
      vector<float> MGRecoBtagId(6,0.0);
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

          //Add CSV b-tag information:
          MGRecoBtagId[0] = HadrBCSV; MGRecoBtagId[3] = LeptBCSV;

	  if(eventselectedSemiEl){//Negative electron
	    MadGraphRecoId[1] = 1;
	    MadGraphRecoId[2] = 6;
	    if(RecoLHCOOutput == true){
	      if(ConsideredCombi == 0) NumberNegRecoEl++;  //Only need to raise the eventNumber for one combination of the 4!!
	      if(verbosity > 4) cout << " Event : " << ievt << " with Number " << NumberNegRecoEl << " sent to LHCO Reco output (Negative electron) " << endl;
	      if(NumberNegRecoEl==1)outFileReco[3*4+ConsideredCombi].open(("MadWeightInput/AnalyzerOutput/TTbarSemiLepton_Reco_NegativeElectron_JetCombi"+JetCombiString+".lhco").c_str());
	      lhcoOutput.LHCOEventRecoOutput(3,outFileReco[3*4+ConsideredCombi], NumberNegRecoEl, LHCORecoVector, MadGraphRecoId, MGRecoBtagId);
	      if(ConsideredCombi == 0) EventInfoFile << "     " << NumberNegRecoEl << endl;  //Only need the output once!
	    }
	  }
	  if(eventselectedSemiMu){//Negative muon
	    MadGraphRecoId[1] = 2;
	    MadGraphRecoId[2] = 6;
	    if(RecoLHCOOutput == true){
	      if(ConsideredCombi == 0) NumberNegRecoMu++;
	      if(verbosity > 4) cout << " Event : " << ievt << " with Number " << NumberNegRecoMu << " sent to LHCO Reco output (Negative muon) " << endl;
	      if(NumberNegRecoMu==1) outFileReco[1*4+ConsideredCombi].open(("MadWeightInput/AnalyzerOutput/TTbarSemiLepton_Reco_NegativeMuon_JetCombi"+JetCombiString+".lhco").c_str());
	      lhcoOutput.LHCOEventRecoOutput(1, outFileReco[1*4+ConsideredCombi], NumberNegRecoMu, LHCORecoVector, MadGraphRecoId, MGRecoBtagId);
	      if(ConsideredCombi == 0) EventInfoFile << "     " << NumberNegRecoMu << endl;
	    }	
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

          //Add CSV b-tag information:
          MGRecoBtagId[0] = HadrBCSV; MGRecoBtagId[3] = LeptBCSV;

	  if(eventselectedSemiEl){//Positive electron
	    MadGraphRecoId[1] = 1;
	    MadGraphRecoId[2] = 6;
	    if(RecoLHCOOutput == true){ 
	      if(ConsideredCombi == 0) NumberPosRecoEl++;
	      if(verbosity > 4) cout << " Event : " << ievt << " with Number " << NumberPosRecoEl << " sent to LHCO Reco output (Positive electron) " << endl;
	      if(NumberPosRecoEl==1)outFileReco[2*4+ConsideredCombi].open(("MadWeightInput/AnalyzerOutput/TTbarSemiLepton_Reco_PositiveElectron_JetCombi"+JetCombiString+".lhco").c_str());
	      lhcoOutput.LHCOEventRecoOutput(2,outFileReco[2*4+ConsideredCombi], NumberPosRecoEl, LHCORecoVector, MadGraphRecoId, MGRecoBtagId);	 
	      if(ConsideredCombi == 0) EventInfoFile << "     " << NumberPosRecoEl << endl;
	    }
	  }
	  if(eventselectedSemiMu){//Positive muon
	    MadGraphRecoId[1] = 2;
	    MadGraphRecoId[2] = 6;
	    if(RecoLHCOOutput == true){
	      if(ConsideredCombi == 0) NumberPosRecoMu++;
	      if(verbosity > 4) cout << " Event : " << ievt << " with Number " << NumberPosRecoMu << " sent to LHCO Reco output (Positive muon) " << endl;
	      if(NumberPosRecoMu==1)outFileReco[0*4+ConsideredCombi].open(("MadWeightInput/AnalyzerOutput/TTbarSemiLepton_Reco_PositiveMuon_JetCombi"+JetCombiString+".lhco").c_str());
	      lhcoOutput.LHCOEventRecoOutput(0, outFileReco[0*4+ConsideredCombi], NumberPosRecoMu, LHCORecoVector, MadGraphRecoId, MGRecoBtagId);
	      if(ConsideredCombi == 0) EventInfoFile << "     " << NumberPosRecoMu << endl;
	    }	
	  }
	}//End of positive lepton
      }//End of loop over the different jet combinations  
    } //loop on events

    // -------- Calculate TF MadWeight  --------//
    if(CalculateTF){
        //tfCreation.CalculateTF(true, true, false, true); //bool drawHistos, bool doFits, bool useROOTClass, bool useStartValues   --> Writing TF to file only possible using TFFit analyzer!

        TFTreeFile->cd();
      
        TTree *configTreeTFFile = new TTree("configTreeTFFile","configuration Tree in Transfer Function Tree file");
        TClonesArray* tcdatasetTFFile = new TClonesArray("Dataset",1);
        configTreeTFFile->Branch("Dataset","TClonesArray",&tcdatasetTFFile);
        TClonesArray* tcAnaEnvTFFile = new TClonesArray("AnalysisEnvironment",1);
        configTreeTFFile->Branch("AnaEnv","TClonesArray",&tcAnaEnvTFFile);
        new ((*tcdatasetTFFile)[0]) Dataset(*datasets[d]);
        new ((*tcAnaEnvTFFile)[0]) AnalysisEnvironment(anaEnv);

        configTreeTFFile->Fill();
        configTreeTFFile->Write();
        TFTree->Write();
        TFTreeFile->Close();
        delete TFTreeFile;

    }
    if(CalculateBTag){	
      //--------------------  Sigma for W Mass and Top Mass  --------------------
      histo1D["MlbMass"]->Fit("gaus","Q");
      histo1D["MqqbMass"]->Fit("gaus","Q");
      cout <<" values for Mlb :"<< histo1D["MlbMass"]->GetFunction("gaus")->GetParameter(1) << " +- " << histo1D["MlbMass"]->GetFunction("gaus")->GetParameter(2) << endl;
      cout <<" values for Mqqb :" << histo1D["MqqbMass"]->GetFunction("gaus")->GetParameter(1) << " +- " << histo1D["MqqbMass"]->GetFunction("gaus")->GetParameter(2) << endl;

      cout << "\n -> " << nSelectedMu << " mu+jets events where selected from which " << nSelectedMuLCSV << " have two or more Light wp CSV b-tags, " << nSelectedMuMCSV << " have two or more Medium wp CSV b-tags and " << nSelectedMuTCSV << " have two or more Tight wp CSV b-tags " << endl;
      cout << "-> " << nSelectedEl << " e+jets events where selected from which " << nSelectedElLCSV << " have two or more Light wp CSV b-tags, " << nSelectedElMCSV << " have two or more Medium wp CSV b-tags and " << nSelectedElTCSV << " have two or more Tight wp CSV b-tags " << endl;
      //cout << "-> " << nLargeMCSVEvents << " events with more than 2 Medium CSV b-tags ( " << nLargeTCSVEvents << " with 2 Tight CSV b-tags) --> Reject these events! " << endl;

      if(verbosity>0) cout << "\n ---> Number of events with correct semileptonic event content on generator level: " << NumberCorrectEvents << " (semiMuon, semiElec) : ( " << NumberPositiveMuons+NumberNegativeMuons << " , " << NumberPositiveElectrons+NumberNegativeElectrons << " ) \n" << endl;

      cout << " \n Output for N(2 light jets) and N(3 light jets) categories: " << endl;
      cout << "   -- Number of events in each category    : " << FourJetsCategory << "   --   " << FiveJetsCategory << endl;
      cout << "   -- Number of matched events             : " << FourJetsMatched <<  "   --   " << FiveJetsMatched << endl;
      cout<<"   -- Number of times good combi is chosen : "<<FourJetsGoodCombiChosen<< "   --   "<<FiveJetsAsFourGoodCombiChosen<<" & "<<FiveJetsAsFiveGoodCombiChosen<<endl;

      //////////////////////////////
      //  Jet combination output  //
      //////////////////////////////
      //--> Save directly to .tex output as a table
      ofstream eventSelOutput;
      int BTagOptionOfInterest;
      std::string OptionName5Jets[6], OptionName4Jets[6];
      if(NrConsideredBTagOptions > 1){
	eventSelOutput.open("/user/aolbrech/GitTopTree_Feb2014/TopBrussels/AnomalousCouplings/eventSelectionChoiceTables.tex");
	BTagOptionOfInterest = 7;
	for(int ii = 0; ii < 6; ii++){OptionName4Jets[ii] = OptionName[ii]; OptionName5Jets[ii] = OptionName[ii];}
      }
      if(NrConsideredBTagOptions == 1){
	eventSelOutput.open("/user/aolbrech/GitTopTree_Feb2014/TopBrussels/AnomalousCouplings/eventSelectionTableForChosenCombination.tex");
	for(int ii = 0; ii < 6; ii++){OptionName4Jets[ii] = " 4 jet case, "+OptionName[ii]; OptionName5Jets[ii] = " 5 jet case, "+OptionName[ii];}
	BTagOptionOfInterest = ChosenBTagOption;
      }

      bTagStudy.ReturnTable(OptionName4Jets, OptionName5Jets, 0, NrConsideredBTagOptions, eventSelOutput, BTagOptionOfInterest);  //0 stands for all 4 particles, 1 for b-jets and 2 for light jets!
      bTagStudy.ReturnTable(OptionName4Jets, OptionName5Jets, 1, NrConsideredBTagOptions, eventSelOutput, BTagOptionOfInterest);
      bTagStudy.ReturnTable(OptionName4Jets, OptionName5Jets, 2, NrConsideredBTagOptions, eventSelOutput, BTagOptionOfInterest);
      eventSelOutput.close();
        
      //////////////////////////////
      //  Mlb combination output  //
      //////////////////////////////
      //mlbStudy.saveNumbers(OptionName, 0, NrConsideredBTagOptions, ChosenBTagOption, ChiSqCutValueStr );  //All 4 jets correctly matched
      mlbStudy.saveNumbers(OptionName, 1, NrConsideredBTagOptions, ChosenBTagOption, ChiSqCutValueStr );  //Also get table for "only b-jets correctly matched"
      mlbStudy.WritePlots(fout);

    } //Only go through all of this output if the TF are not being calculated!

    //Close the LHCO Output files!
    for(int ii = 0; ii<16; ii++){
      if(ii < 4) outFile[ii].close();	
      outFileReco[ii].close();
    }
    if(GenLHCOOutput == true) EventInfoFile.close();

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
  }  //loop on datasets
  
  //Once everything is filled ...
  if (verbose > 0)
    cout << " We ran over all the data ;-)" << endl;
  
  /////////////////////////
  // Write out the plots //
  /////////////////////////
  // Fill the resolution histograms and calculate the resolutions
  //
  string pathPNG = "PlotsMacro";
  if(CalculateResolutions){    
    mkdir((pathPNG+"/resFit_LightJet/").c_str(),0777);
    mkdir((pathPNG+"/resFit_BJet/").c_str(),0777);
    mkdir((pathPNG+"/resFit_Muon/").c_str(),0777);
    mkdir((pathPNG+"/resFit_NeutrinoSemiMu/").c_str(),0777);     
    mkdir((pathPNG+"/resFit_NeutrinoSemiEl/").c_str(),0777);
    mkdir((pathPNG+"/resFit_Electron/").c_str(),0777);

    resFitMuon->WritePlots(fout, true, pathPNG+"resFit_Muon/");
    resFitMuon->WriteResolutions("PersonalClasses/resolutions/Calculated/muonReso.root");

    resFitNeutrinoMu->WritePlots(fout, true, pathPNG+"/resFit_NeutrinoSemiMu/");
    resFitNeutrinoMu->WriteResolutions("PersonalClasses/resolutions/Calculated/neutrinoSemiMuReso.root");

    resFitNeutrinoEl->WritePlots(fout, true, pathPNG+"/resFit_Neutrino/");
    resFitNeutrinoEl->WriteResolutions("PersonalClasses/resolutions/Calculated/neutrinoSemiElReso.root");

    resFitElectron->WritePlots(fout, true, pathPNG+"/resFit_Electron/");
    resFitElectron->WriteResolutions("PersonalClasses/resolutions/Calculated/electronReso.root");      
      
    resFitLightJets->WritePlots(fout, true, pathPNG+"/resFit_LightJet/");
    resFitLightJets->WriteResolutions("PersonalClasses/resolutions/Calculated/lightJetReso.root");

    resFitBJets->WritePlots(fout, true, pathPNG+"/resFit_BJet/");
    resFitBJets->WriteResolutions("PersonalClasses/resolutions/Calculated/bJetReso.root");
  }

  //Write away the MSPlots, the histo1D's and the histo2D's!!
  fout -> cd();
  mkdir((pathPNG+"/MSPlots").c_str(),0777);
 
  TDirectory* msdir = fout->mkdir("MSPlots");
  msdir->cd(); 
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++){    
    MultiSamplePlot *temp = it->second;
    string name = it->first;
    temp->Draw(name, 0, false, false, false, 1);     //string label, unsigned int RatioType, bool addRatioErrorBand, bool addErrorBand, bool ErrorBandAroundTotalInput, int scaleNPSSignal 
    temp->Write(fout, name, true, "PlotsMacro/MSPlots/", "png");
  }

  TDirectory* th1dir = fout->mkdir("1D_histograms");
  th1dir->cd();
  for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++){
    TH1F *temp = it->second;
    int N = temp->GetNbinsX();
    temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
    temp->SetBinContent(N+1,0);
    temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
    temp->Write();
  }
  TDirectory* th2dir = fout->mkdir("2D_histograms_graphs");
  th2dir->cd();
  for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++){
    TH2F *temp = it->second;
    temp->Write();
  }

  //Selection tables
  selecTableSemiMu.TableCalculator(false,true,true,true,true);
  string selectiontableMu = "SelectionTable_BTAG_SEMIMU.tex";
  selecTableSemiMu.Write(selectiontableMu.c_str());
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
