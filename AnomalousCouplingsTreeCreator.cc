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
#include "AnomalousCouplings/PersonalClasses/interface/TFnTuple.h"
#include "AnomalousCouplings/PersonalClasses/interface/KinematicFunctions.h"

using namespace std;
using namespace reweight;
using namespace TopTree;

int main (int argc, char *argv[]){

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
  bool RecoLHCOOutput = true;  
  bool CalculateResolutions = false; // If false, the resolutions will be loaded from a previous calculation
  bool CalculateTF = false;

  //Values needed for bTag study (select which of the 6 b-tag options is optimal!)
  int ChosenBTag;

  int ChiSqCutValue =51;  //The Chi-sq values in the mlb method has to be larger than this value! (Put on 51 to include all events, since the chi-sq is set manually to a maximum of 49.5)
  std::string ChiSqCutValueStr;
  ostringstream convert; convert << ChiSqCutValue;
  if(ChiSqCutValue != 51) ChiSqCutValueStr = "_ChiSqSmallerThan"+convert.str();
  else ChiSqCutValueStr = "";

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

    if(dataSetName.find("QCD") == 0)                                      datasets[d]->SetColor(kYellow);
    if(dataSetName.find("TT") == 0)                                       datasets[d]->SetColor(kRed+1);
    if(dataSetName.find("TTbarJets_Other") == 0)                          datasets[d]->SetColor(kRed-7);
    if(dataSetName.find("ST") == 0 || dataSetName.find("SingleTop") == 0) datasets[d]->SetColor(kMagenta);
    if(dataSetName.find("WJets") == 0){                                   datasets[d]->SetColor(kGreen-3); datasets[d]->SetTitle("W#rightarrowl#nu");}      
    if(dataSetName.find("ZJets") == 0){                                   datasets[d]->SetColor(kAzure-2); datasets[d]->SetTitle("Z/#gamma*#rightarrowl^{+}l^{-}");}
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
    
    if (verbose > 1){
      cout << "   Dataset " << d << ": " << datasets[d]->Name () << "/ title : " << datasets[d]->Title () << endl;
      std::cout<<"      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
    }
    
    //open files and load
    treeLoader.LoadDataset (datasets[d], anaEnv);
    
    /////////////////////////////////////
    /// Initialize JEC factors            --> Updated on 5/08/2014
    /////////////////////////////////////
    vector<JetCorrectorParameters> vCorrParam;
    
    if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ){// Data!
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("PersonalClasses/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L1FastJet_AK5PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("PersonalClasses/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L2Relative_AK5PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("PersonalClasses/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L3Absolute_AK5PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
      JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("PersonalClasses/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L2L3Residual_AK5PFchs.txt");
      vCorrParam.push_back(*L2L3ResJetCorPar);
    }
    else{
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
   
    /////////////////////
    //  Used classes   //
    /////////////////////  
    BTagStudy bTagStudy(verbose);  
    //MlbStudy mlbStudy(6);//NrConsideredBTagOptions);
    //mlbStudy.InitializeMlb();
    LHCOOutput lhcoOutput(verbose, GenLHCOOutput, RecoLHCOOutput); 
    KinematicFunctions kinFunctions;  //Variable accessible in KinematicFunctions using kinFunctions.CosTheta(TLorentzVector *Top, TLorentzVector *WLept, TLorentzVector *lepton)
    TFnTuple* tfNTuple = 0;

    //Initialize TFnTuple specific stuff:
    TTree* TFTree = new TTree("TFTree","Tree containing the Transfer Function information");
    TFTree->Branch("TheTFTree","TFnTuple",&tfNTuple);
    TFile* TFTreeFile = new TFile("TFInformation/TransferFunctionTree.root","RECREATE");

    /////////////////////////////////////////
    //  LHCO Output files + GeneratorInfo  //
    /////////////////////////////////////////
    ofstream EventInfoFile;
    if(GenLHCOOutput == true){	
      EventInfoFile.open("MadWeightInput/AnalyzerOutput/EventNumberInformation.lhco");
      EventInfoFile << " Event Number      Lepton Type       Event selection       selectedChannelNumber " << endl;
    }
    unsigned int NumberCorrectEvents = 0; //Counts the number of semi-leptonic events
    
    ////////////////////////////////////
    //	loop on events
    ////////////////////////////////////
    
    nEvents[d] = 0;
    int itriggerSemiMu = -1,itriggerSemiEl = -1, previousRun = -1;
    if (verbose > 1) cout << "	Loop over events " << endl;

    for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++){
    //for (unsigned int ievt = 0; ievt < 50000; ievt++){
      
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
      ////////////////////////////////////////////////////////         -->This part only needs mcParticles information so can be moved to class!!
      lhcoOutput.StoreGenInfo(mcParticles);
      
      //Accessing information and store in EventInfoFile    
      std::string leptonTypeString[4] = {"muPlus","muMinus","elPlus","elMinus"};      //Not possible to define this in header file of LHCOOutput ...
      
      if( GenLHCOOutput == true){
        EventInfoFile << "     " << ievt << "         ";
        if( lhcoOutput.GenEventContentCorrect() ){
	  NumberCorrectEvents++;
          EventInfoFile << "     " << leptonTypeString[lhcoOutput.getLeptonType()] << "    ";
	}//End of Gen event with correct content    
        else    
          EventInfoFile << "        0      ";  //Output when event content is wrong!
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
      
      ///////////////////////////////////////
      //  Initialize variables ChiSquared  //
      //  Look for correct combination     //
      //  --> Only used for b-quarks       //
      ///////////////////////////////////////
      int UsedCombination, BHadronicIndex[2], BLeptIndex, BHadrIndex, QOneIndex, QTwoIndex;
      float ChiSquared[2], ChiSquaredValue;
      
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
      
      //--- Only continue with events passing one of the two event selections!  ---//
      if (!eventselectedSemiMu && !eventselectedSemiEl){
        if(RecoLHCOOutput == true) EventInfoFile << "        failed       "<<endl;
        continue;
      }
      
      if(RecoLHCOOutput == true) EventInfoFile << "        passed            ";  
      //--- do some data-mc ---//
      if (!foundMu && !foundEl) datasetsPlot = datasets;    //Data!
      string leptonFlav="_other";
      TLorentzVector* selectedLepton;
      float LeptonRecoCharge;
      if (eventselectedSemiMu){
        nSelectedMu++; decayChannel = semiMu;
        selectedLepton = (TLorentzVector*)selectedMuons[0]; LeptonRecoCharge = selectedMuons[0]->charge();
	datasetsPlot = datasetsMu; Luminosity = LuminosityMu;
        leptonFlav="_mu";
      }
      else if (eventselectedSemiEl){
        nSelectedEl++; decayChannel = semiEl;
        selectedLepton = (TLorentzVector*)selectedElectrons[0]; LeptonRecoCharge = selectedElectrons[0]->charge();
	datasetsPlot = datasetsEl; Luminosity = LuminosityEl;
        leptonFlav="_el";
      }
      
      // MSPlots after 'basic' event selection (no b-tag)
      MSPlot["Selected_Events_pT_jet1"+leptonFlav]->Fill(selectedJets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["Selected_Events_pT_jet2"+leptonFlav]->Fill(selectedJets[1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["Selected_Events_pT_jet3"+leptonFlav]->Fill(selectedJets[2]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["Selected_Events_pT_jet4"+leptonFlav]->Fill(selectedJets[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["Selected_Events_pT_lepton"+leptonFlav]->Fill(selectedLepton->Pt(), datasets[d], true, Luminosity*scaleFactor);
      
      for (unsigned int q=0; q<selectedJets.size(); q++) {
	MSPlot["Selected_Events_pT_alljets"+leptonFlav]->Fill(selectedJets[q]->Pt(), datasets[d], true, Luminosity*scaleFactor);       
	if(q<4)
	  MSPlot["Selected_Events_pT_4leadingjets"+leptonFlav]->Fill(selectedJets[q]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      }
      
      ////////////////////////////////////////////////////////////////////
      //   Use genEvent information to get the correct event topology   //
      ////////////////////////////////////////////////////////////////////    
      if(dataSetName.find("TTbarJets") == 0){
        vector<int> jetCombi;
        vector<TRootMCParticle> mcParticlesMatching;      	
        vector< pair<unsigned int, unsigned int> > JetPartonPair, ISRJetPartonPair; // First one is jet number, second one is mcParticle number
	
	//First index is the jet number, the second one the mcParticle number:
	pair<unsigned int, unsigned int> leptonicBJet_, hadronicBJet_, hadronicWJet1_, hadronicWJet2_;
	leptonicBJet_ = hadronicBJet_ = hadronicWJet1_ = hadronicWJet2_ = pair<unsigned int, unsigned int>(9999,9999);
	vector<TLorentzVector> mcParticlesTLV, selectedJetsTLV;
	bool muPlusFromTop = false, muMinusFromTop = false, elPlusFromTop = false, elMinusFromTop = false;
	if(verbosity>1) cout << " Looking at mcParticlesMatching " << endl;

	for(unsigned int i=0; i<mcParticles.size(); i++){
	  if(mcParticles[i]->status() != 3) continue;
	  
	  //Muon identification:
	  if(mcParticles[i]->type() == 13 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -6){if(muMinusFromTop) cerr<<"muMinusFromTop already true"<<endl; muMinusFromTop = true;}
	  if(mcParticles[i]->type() == -13 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == 6){if(muPlusFromTop) cerr<<"muPlusFromTop already true"<<endl; muPlusFromTop = true;}
	  
	  //Electron identification:
	  if(mcParticles[i]->type() == 11 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -6){if(elMinusFromTop) cerr<<"elMinusFromTop already true"<<endl; elMinusFromTop = true;}
	  if(mcParticles[i]->type() == -11 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == 6){if(elPlusFromTop) cerr<<"elPlusFromTop already true"<<endl; elPlusFromTop = true;}
	  
	  if( abs(mcParticles[i]->type()) < 6 || abs(mcParticles[i]->type()) == 21 ){
	    mcParticlesTLV.push_back(*mcParticles[i]); mcParticlesMatching.push_back(*mcParticles[i]);
	  }
	}

        //--- Control checks ---//
	if(muPlusFromTop && muMinusFromTop) cerr<<"muPlusFromTop and muMinusFromTop are both true ?!\nCheck if you are using the right sample..."<<endl;
	if(elPlusFromTop && elMinusFromTop) cerr<<"elPlusFromTop and elMinusFromTop are both true ?!\nCheck if you are using the right sample..."<<endl;

        //--- MATCHING ---//	
	for(unsigned int i=0; i<selectedJets.size(); i++) selectedJetsTLV.push_back(*selectedJets[i]);  //Need selectedJets for radiation stuff!
	
	JetPartonMatching matching = JetPartonMatching(mcParticlesTLV, selectedJetsTLV, 2, true, true, 0.3);
	if(matching.getNumberOfAvailableCombinations() != 1) cerr << "matching.getNumberOfAvailableCombinations() = "<<matching.getNumberOfAvailableCombinations()<<"  This should be equal to 1 !!!"<<endl;
	for(unsigned int i=0; i<mcParticlesTLV.size(); i++){
	  int matchedJetNumber = matching.getMatchForParton(i, 0);
	  if(matchedJetNumber != -1) JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
	}
	
	for(unsigned int itJPPair=0; itJPPair<JetPartonPair.size(); itJPPair++){
	  unsigned int mcPartNr = JetPartonPair[itJPPair].second;
	  
	  if(verbosity > 3){ 
	    std::cout <<" Jet number " << JetPartonPair[itJPPair].first << " can be matched with the mcParticle number " << JetPartonPair[itJPPair].second << " which is of the type : " << mcParticlesMatching[mcPartNr].type() << std::endl;
	  }

          //--- Recover the light jets ---//	  
	  if( fabs(mcParticlesMatching[mcPartNr].type()) < 6 ){
	    if( ( (muPlusFromTop  || elPlusFromTop)  && mcParticlesMatching[mcPartNr].motherType() == -24 && mcParticlesMatching[mcPartNr].grannyType() == -6 ) || 
                ( (muMinusFromTop || elMinusFromTop) && mcParticlesMatching[mcPartNr].motherType() == 24  && mcParticlesMatching[mcPartNr].grannyType() == 6 )  ){
	      if(hadronicWJet1_.first == 9999) 
		hadronicWJet1_ = JetPartonPair[itJPPair];
	      else if(hadronicWJet2_.first == 9999) 
		hadronicWJet2_ = JetPartonPair[itJPPair];
	      else cerr<<"Found a third jet coming from a W boson which comes from a top quark..."<<endl;
	    }
	  }

          //--- Recover the b-jets ---//
	  if( fabs(mcParticlesMatching[mcPartNr].type()) == 5 ){
	    if( ( (muPlusFromTop  || elPlusFromTop)  && mcParticlesMatching[mcPartNr].motherType() == -6 ) || 
                ( (muMinusFromTop || elMinusFromTop) && mcParticlesMatching[mcPartNr].motherType() == 6  ) )
	      hadronicBJet_ = JetPartonPair[itJPPair];
	    
	    else if( ( (muPlusFromTop  || elPlusFromTop)  && mcParticlesMatching[mcPartNr].motherType() == 6 ) || 
                     ( (muMinusFromTop || elMinusFromTop) && mcParticlesMatching[mcPartNr].motherType() == -6) )
	      leptonicBJet_ = JetPartonPair[itJPPair];
	  }
	}
	
        //--- Initialize jetCombi vector which has the jet number information ---//
	jetCombi.push_back(leptonicBJet_.first);
	jetCombi.push_back(hadronicBJet_.first);
	jetCombi.push_back(hadronicWJet1_.first);
	jetCombi.push_back(hadronicWJet2_.first);
	if(verbose > 3) cout<<" Index of BHadronic: "<<jetCombi[1]<<" , Index of BLeptonic: "<<jetCombi[0]<<" , Index of quark1: "<<jetCombi[2]<<" & Index of quark2: "<<jetCombi[3]<<endl;
	
        //--- Make CSV discriminant and jet number plots ---//	
	if(jetCombi[0] != 9999) histo1D["CorrectBLeptCSVDiscr"]->Fill(selectedJets[jetCombi[0]]->btag_combinedSecondaryVertexBJetTags());
	else                    histo1D["CorrectBLeptCSVDiscr"]->Fill(-2);
	if(jetCombi[1] != 9999) histo1D["CorrectBHadrCSVDiscr"]->Fill(selectedJets[jetCombi[1]]->btag_combinedSecondaryVertexBJetTags());
	else                    histo1D["CorrectBHadrCSVDiscr"]->Fill(-2);
	if(jetCombi[2] != 9999){histo1D["CorrectQuark1CSVDiscr"]->Fill(selectedJets[jetCombi[2]]->btag_combinedSecondaryVertexBJetTags()); histo1D["Quark1JetNumber"]->Fill(jetCombi[2]);}
	else                   {histo1D["CorrectQuark1CSVDiscr"]->Fill(-2);                                                                histo1D["Quark1JetNumber"]->Fill(-1);}
	if(jetCombi[3] != 9999){histo1D["CorrectQuark2CSVDiscr"]->Fill(selectedJets[jetCombi[3]]->btag_combinedSecondaryVertexBJetTags()); histo1D["Quark2JetNumber"]->Fill(jetCombi[3]);}
	else                   {histo1D["CorrectQuark2CSVDiscr"]->Fill(-2);                                                                histo1D["Quark2JetNumber"]->Fill(-1);}
 	
	//Working on generator level (i.e. parton level):  
	if(jetCombi[0]!=9999 && jetCombi[1]!=9999 && jetCombi[2]!=9999 && jetCombi[3]!=9999){    
	  histo1D["MlbMass"]->Fill((*selectedJets[jetCombi[0]]+*selectedLepton).M());
	  histo1D["MqqbMass"]->Fill((*selectedJets[jetCombi[2]]+*selectedJets[jetCombi[3]]+*selectedJets[jetCombi[1]]).M());
	
	  if(CalculateResolutions){
	    // Fill the resolution-stuff for events where the 4 ttbar semi-lep partons are all matched to jets
	    resFitLightJets->Fill(selectedJets[jetCombi[2]], mcParticles[hadronicWJet1_.second]);
	    resFitLightJets->Fill(selectedJets[jetCombi[3]], mcParticles[hadronicWJet2_.second]);
	    resFitBJets->Fill(selectedJets[jetCombi[1]], mcParticles[hadronicBJet_.second]);
	    resFitBJets->Fill(selectedJets[jetCombi[0]], mcParticles[leptonicBJet_.second]);
	    if(decayChannel == 0){//semiMu case
	      resFitMuon->Fill(selectedMuons[0], lhcoOutput.getGenLepton());         resFitNeutrinoMu->Fill(mets[0], lhcoOutput.getGenNeutrino());}
	    else if(decayChannel == 1){
	      resFitElectron->Fill(selectedElectrons[0], lhcoOutput.getGenLepton()); resFitNeutrinoEl->Fill(mets[0], lhcoOutput.getGenNeutrino());}
	  }//End of calculate Resolutions

	  if(CalculateTF){
	    if(decayChannel == 0){ histo1D["genPt_Muon"]->Fill( lhcoOutput.getGenLepton()->Pt() ); histo1D["recoPt_Muon"]->Fill(selectedLepton->Pt());}	    
	    if(decayChannel == 1){ histo1D["genPt_Elec"]->Fill( lhcoOutput.getGenLepton()->Pt() ); histo1D["recoPt_Elec"]->Fill(selectedLepton->Pt());}	    
	    
	    //Fill the Transfer Function Tree file!
	    tfNTuple = new TFnTuple();
	    tfNTuple->setEventID( event->eventId() );
	    tfNTuple->setRecoVectorLight1( (TLorentzVector) *selectedJets[jetCombi[2]] );
	    tfNTuple->setRecoVectorLight2( (TLorentzVector) *selectedJets[jetCombi[3]] );
	    tfNTuple->setRecoVectorHadrB(  (TLorentzVector) *selectedJets[jetCombi[1]] );
	    tfNTuple->setRecoVectorLeptB(  (TLorentzVector) *selectedJets[jetCombi[0]] );
	    tfNTuple->setRecoVectorLepton( (TLorentzVector) *selectedLepton            );
	    tfNTuple->setGenVectorLight1(  (TLorentzVector) mcParticlesMatching[hadronicWJet1_.second] );
	    tfNTuple->setGenVectorLight2(  (TLorentzVector) mcParticlesMatching[hadronicWJet2_.second] );
	    tfNTuple->setGenVectorHadrB(   (TLorentzVector) mcParticlesMatching[hadronicBJet_.second]  );
	    tfNTuple->setGenVectorLeptB(   (TLorentzVector) mcParticlesMatching[leptonicBJet_.second]  );
	    tfNTuple->setGenVectorLepton(  *lhcoOutput.getGenLepton()                                   );
	    TFTree->Fill();
	    delete tfNTuple;
	  }//End of calculate Transfer Functions
	}//End of matched particles reconstructed
      
        //----------------------------------------------------------------------------------------------------------------------------------- Start of bTagStudy class stuff!!      	
	if(lhcoOutput.GenEventContentCorrect()) histo1D["StCosThetaNoBTag"]->Fill(kinFunctions.CosTheta(lhcoOutput.getGenLeptTop(), lhcoOutput.getGenLeptW(), lhcoOutput.getGenLepton())); 	
	
	//---  Get the b-jet and light information for all six b-tag options!  ---//
	bTagStudy.CalculateJets(selectedJets, jetCombi, selectedLepton);
        float bTagLeptIndex = bTagStudy.getBLeptIndex(3);

	//////////////////////////////////////////////////////
	//  Kinematic information for specific b-tag cases   //
	///////////////////////////////////////////////////////
	//
	/*
        for(int jj = 0; jj < 2; jj++){      //Only look at the two leading jets!
	    
	  //Kinematic information for the light (= Non LCSV) jets
	  if((bTagStudy.getLightJets(0)).size()>=2){
	    histo1D["CSVDiscrLCSVLightJets"]->Fill(selectedJets[(bTagStudy.getLightJets(0))[jj]]->btag_combinedSecondaryVertexBJetTags());   
    	    //--> check whether the so-called light jets don't all have discr -1 ...		
	    for(int ii = 0; ii<JetPartonPair.size(); ii++){ //Look at all the matched jets
	      if((bTagStudy.getLightJets(0))[jj] == JetPartonPair[ii].first)  //Check whether the considered jet can be matched!
	        histo1D["JetTypeLCSVLightJets"]->Fill(mcParticlesMatching[JetPartonPair[ii].second].type());	
	      else    //Unmatched jets!
	        histo1D["JetTypeLCSVLightJets"]->Fill(25.);	      
	    }
	  }//End of light (= non LCSV) jets

	  //Kinematic information for the Loose b-jets
	  if((bTagStudy.getbTaggedJets(0)).size() >=2){
	    if(lhcoOutput.GenEventContentCorrect()){ histo1D["StCosThetaLCSV"]->Fill(kinFunctions.CosTheta(lhcoOutput.getGenLeptTop(), lhcoOutput.getGenLeptW(), lhcoOutput.getGenLepton())); std::cout << " LCSV event found with cos theta = " << kinFunctions.CosTheta(lhcoOutput.getGenLeptTop(), lhcoOutput.getGenLeptW(), lhcoOutput.getGenLepton()) << std::endl;}
	    for(int ii = 0; ii<JetPartonPair.size(); ii++){ //Look at all the matched jets
	      if((bTagStudy.getbTaggedJets(0))[jj] == JetPartonPair[ii].first)  //Check whether the considered jet can be matched!
	        histo1D["JetTypeLCSV"]->Fill(mcParticlesMatching[JetPartonPair[ii].second].type());	    
	      else
	        histo1D["JetTypeLCSV"]->Fill(25.);	      
	    }
	  }//End of Loose b-jets 

	  //Kinematic information for the Medium b-jets
	  if((bTagStudy.getbTaggedJets(1)).size() >=2){
	    if(lhcoOutput.GenEventContentCorrect()) histo1D["StCosThetaMCSV"]->Fill(kinFunctions.CosTheta(lhcoOutput.getGenLeptTop(), lhcoOutput.getGenLeptW(), lhcoOutput.getGenLepton()));
	    for(int ii = 0; ii<JetPartonPair.size(); ii++){ //Look at all the matched jets
  	      if((bTagStudy.getbTaggedJets(1))[jj] == JetPartonPair[ii].first)  //Check whether the considered jet can be matched!
	        histo1D["JetTypeMCSV"]->Fill(mcParticlesMatching[JetPartonPair[ii].second].type());	    
	      else
	        histo1D["JetTypeMCSV"]->Fill(25.);	    
	    }
	  }//End of Medium b-jets 

	  //Kinematic information for the Tight b-jets
	  if((bTagStudy.getbTaggedJets(3)).size() >=2){
	    if(lhcoOutput.GenEventContentCorrect()) histo1D["StCosThetaTCSV"]->Fill(kinFunctions.CosTheta(lhcoOutput.getGenLeptTop(), lhcoOutput.getGenLeptW(), lhcoOutput.getGenLepton()));
	    for(int ii = 0; ii<JetPartonPair.size(); ii++){ //Look at all the matched jets
	      if((bTagStudy.getbTaggedJets(3))[jj] == JetPartonPair[ii].first)  //Check whether the considered jet can be matched!
	        histo1D["JetTypeTCSV"]->Fill(mcParticlesMatching[JetPartonPair[ii].second].type());	    
	      else
	        histo1D["JetTypeTCSV"]->Fill(25.);	      
	    }
	  }//End of Tight b-jets 
        }//Only look at the two leading jets!

        /////////////////////////////////////////////
        //   Count the selected number of events   //
        /////////////////////////////////////////////
        // 2 Loose CSV b-tags
        if((bTagStudy.getbTaggedJets(0)).size() >= 2 && decayChannel == 0) nSelectedMuLCSV++;
        if((bTagStudy.getbTaggedJets(0)).size() >= 2 && decayChannel == 1) nSelectedElLCSV++;
        if((bTagStudy.getbTaggedJets(0)).size() > 2) nLargeLCSVEvents++;

        // 2 Medium CSV b-tags
        if((bTagStudy.getbTaggedJets(1)).size() >= 2 && decayChannel == 0) nSelectedMuMCSV++;
        if((bTagStudy.getbTaggedJets(1)).size() >= 2 && decayChannel == 1) nSelectedElMCSV++;
        if((bTagStudy.getbTaggedJets(1)).size() > 2) nLargeMCSVEvents++;

        // 2 Tight CSV b-tags
        if((bTagStudy.getbTaggedJets(3)).size() >= 2 && decayChannel == 0) nSelectedMuTCSV++;
        if((bTagStudy.getbTaggedJets(3)).size() >= 2 && decayChannel == 1) nSelectedElTCSV++;
        if((bTagStudy.getbTaggedJets(3)).size() > 2) nLargeTCSVEvents++;
	*/
      
        //---------------------------------------------------------------------------------------------------------------------------- End of bTagStudy class stuff

        ////////////////////////////////
        //  Mlb and Mqqb information  //
        ////////////////////////////////
        float MlbCorrect = 0, MqqbCorrect = 0;
        if(jetCombi[0] != 9999 && jetCombi[1] != 9999 && jetCombi[2] != 9999 && jetCombi[3] != 9999){
	  MlbCorrect = (*selectedLepton+*selectedJets[jetCombi[0]]).M();
	  MqqbCorrect = (*selectedJets[jetCombi[1]] + *selectedJets[jetCombi[2]] + *selectedJets[jetCombi[3]]).M();
	  histo2D["MlbMqqbCorrectAll"]->Fill(MqqbCorrect,MlbCorrect);
        }
      
      }//End of TTbarJets!
      
      //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
      //  Event selection choice (17/06/2014)  //
      //   --> Continue with 2 T b-tags        //
      //   --> No veto on light jets!          //
      //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
      ChosenBTag = 3;

      //MSPlots with number of jets information before requiring at least two b-jets and at least 2 light jets!
      MSPlot["nSelectedJets_BeforeBTag"+leptonFlav]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["nBTaggedJets_BeforeBTag"+leptonFlav]->Fill( (bTagStudy.getbTaggedJets(ChosenBTag)).size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["nLightJets_BeforeBTag"+leptonFlav]->Fill( (bTagStudy.getLightJets(ChosenBTag)).size(), datasets[d], true, Luminosity*scaleFactor);

      //---  Require two b-jets and two light jets!  ---//
      if( (bTagStudy.getbTaggedJets(ChosenBTag)).size() < 2 || (bTagStudy.getLightJets(ChosenBTag)).size() < 2 ){
	if(RecoLHCOOutput == true) EventInfoFile<<"    B-tag failed "<<endl;
	continue;
      }      
            
      //Identical MSPlots with number of jets information after requiring at least two b-jets and at least 2 light jets!
      MSPlot["nSelectedJets_AfterBTag"+leptonFlav]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["nBTaggedJets_AfterBTag"+leptonFlav]->Fill( (bTagStudy.getbTaggedJets(ChosenBTag)).size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["nLightJets_AfterBTag"+leptonFlav]->Fill( (bTagStudy.getLightJets(ChosenBTag)).size(), datasets[d], true, Luminosity*scaleFactor);
      
      // Count the number of events:
      if(decayChannel == 0) nSelectedMu++;
      if(decayChannel == 1) nSelectedEl++;
      
      //****************************//   
      //  Produce Reco LHCO Output  //
      //****************************//
      if( RecoLHCOOutput == true)
        lhcoOutput.StoreRecoInfo(selectedLepton, selectedJets, bTagStudy.getBLeptIndex(ChosenBTag), bTagStudy.getBHadrIndex(ChosenBTag), bTagStudy.getLight1Index5Jets(ChosenBTag), bTagStudy.getLight2Index5Jets(ChosenBTag), decayChannel, LeptonRecoCharge);

    } //loop on events

    // -------- Calculate TF MadWeight  --------//
    if(CalculateTF){
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

    //---  Mlb and Mqqb fit result ---//
    histo1D["MlbMass"]->Fit("gaus","Q");
    histo1D["MqqbMass"]->Fit("gaus","Q");
    cout <<"\n values for Mlb :"<< histo1D["MlbMass"]->GetFunction("gaus")->GetParameter(1) << " +- " << histo1D["MlbMass"]->GetFunction("gaus")->GetParameter(2) << endl;
    cout <<" values for Mqqb :" << histo1D["MqqbMass"]->GetFunction("gaus")->GetParameter(1) << " +- " << histo1D["MqqbMass"]->GetFunction("gaus")->GetParameter(2) << endl;

    cout << "\n -> " << nSelectedMu << " mu+jets and " << nSelectedEl << " e+jets events where selected" << endl;

    //--- Get output from bTagStudy class ---//
    bTagStudy.ReturnBTagTable();
        
    //---  Mlb combination output  ---//
    //mlbStudy.saveNumbers(OptionName, 0, NrConsideredBTagOptions, ChosenBTag, ChiSqCutValueStr );  //All 4 jets correctly matched
    //mlbStudy.saveNumbers(OptionName, 1, 6, ChosenBTag, ChiSqCutValueStr); //NrConsideredBTagOptions, ChosenBTag, ChiSqCutValueStr );  //Also get table for "only b-jets correctly matched"
    //mlbStudy.WritePlots(fout);

    if(GenLHCOOutput == true) EventInfoFile.close();

    //////////////
    // CLEANING //
    //////////////
    
    if (jecUnc) delete jecUnc;
    if (jetTools) delete jetTools;
  
    //important: free memory
    treeLoader.UnLoadDataset();    
  }  //loop on datasets
  
  //Once everything is filled ...
  if(verbose > 0)
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
