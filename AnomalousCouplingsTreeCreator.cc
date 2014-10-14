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

//Specific code for anomalous couplings analysis:
#include "TopTreeAnalysis/AnomCouplings/interface/LHCOOutput.h"
#include "TopTreeAnalysis/AnomCouplings/src/LHCOOutput.cc"
#include "TopTreeAnalysis/AnomCouplings/interface/BTagStudy.h"
#include "TopTreeAnalysis/AnomCouplings/src/BTagStudy.cc"       //--> Need to fix the makefile (don't include .cc files ...)
#include "TopTreeAnalysis/AnomCouplings/src/MlbStudy.cc"
#include "TopTreeAnalysis/AnomCouplings/interface/MlbStudy.h"
#include "TopTreeAnalysis/AnomCouplings/interface/TFCreation.h"
#include "TopTreeAnalysis/AnomCouplings/src/TFCreation.cc"
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
  histo1D["MqqbMass"]= new TH1F("MqqbMass","MqqbMass",400,0,500);

  histo1D["DeltaR_Ana_Muon"] = new TH1F("DeltaR_Ana_Muon","DeltaR_Ana_Muon",200,0,1);
  histo1D["DeltaR_Ana_Elec"] = new TH1F("DeltaR_Ana_Elec","DeltaR_Ana_Elec",200,0,1);
  histo1D["DeltaR_Ana_Light1_mcPartMatch"] = new TH1F("DeltaR_Ana_Light1_mcPartMatch","DeltaR_Ana_Light1_mcPartMatch",200,0,1);
  histo1D["DeltaR_Ana_Light1_mcPart"] = new TH1F("DeltaR_Ana_Light1_mcPart","DeltaR_Ana_Light1_mcPart",200,0,10);
  histo1D["DeltaR_Ana_Light2_mcPartMatch"] = new TH1F("DeltaR_Ana_Light2_mcPartMatch","DeltaR_Ana_Light2_mcMatching",200,0,1);
  histo1D["DeltaR_Ana_Light2_mcPart"] = new TH1F("DeltaR_Ana_Light2_mcPart","DeltaR_Ana_Light2_mcPart",200,0,10);
  histo1D["DeltaR_Ana_BHadr_mcPartMatch"] = new TH1F("DeltaR_Ana_BHadr_mcPartMatch","DeltaR_Ana_BHadr_mcPartMatch",200,0,1);
  histo1D["DeltaR_Ana_BHadr_mcPart"] = new TH1F("DeltaR_Ana_BHadr_mcPart","DeltaR_Ana_BHadr_mcPart",200,0,10);
  histo1D["DeltaR_Ana_BLept_mcPartMatch"] = new TH1F("DeltaR_Ana_BLept_mcPartMatch","DeltaR_Ana_BLept_mcPartMatch",200,0,1);
  histo1D["DeltaR_Ana_BLept_mcPart"] = new TH1F("DeltaR_Ana_BLept_mcPart","DeltaR_Ana_BLept_mcPart",200,0,10);

  histo1D["genPt_Muon"] = new TH1F("genPt_Muon","genPt_Muon",400,0,200);
  histo1D["recoPt_Muon"] = new TH1F("recoPt_Muon","recoPt_Muon",400,0,200);
  histo1D["genPt_Elec"] = new TH1F("genPt_Elec","genPt_Elec",400,0,200);
  histo1D["recoPt_Elec"] = new TH1F("recoPt_Elec","recoPt_Elec",400,0,200);

  ////////////////////////////////////
  /// MultiSamplePlot
  ////////////////////////////////////

  map<string,MultiSamplePlot*> MSPlot;

  /////////////////////////////
  /// ResolutionFit Stuff
  /////////////////////////////

  bool CalculateResolutions = false; // If false, the resolutions will be loaded from a previous calculation
  bool CalculateTF = true;

  std::cout << " CalculateResolutions = " << CalculateResolutions << endl;

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

  LumiWeights = LumiReWeighting("PersonalClasses/Calibrations/PUReweighting/pileup_MC_Summer12_S10.root", "PersonalClassesCalibrations/PUReweighting/pileup_2012Data53X_UpToRun208357/nominal.root", "pileup", "pileup");
  LumiWeightsUp = LumiReWeighting("PersonalClasses/Calibrations/PUReweighting/pileup_MC_Summer12_S10.root", "PersonalClasses/Calibrations/PUReweighting/pileup_2012Data53X_UpToRun208357/sys_up.root", "pileup", "pileup");
  LumiWeightsDown = LumiReWeighting("PersonalClassesCalibrations/PUReweighting/pileup_MC_Summer12_S10.root", "PersonalClasses/Calibrations/PUReweighting/pileup_2012Data53X_UpToRun208357/sys_down.root", "pileup", "pileup");
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
    TH1F h_MqqbMass("MqqbMass","MqqbMass",400,0,500);
    
    TH2F h_MlbMqqbCorrectChosen("MlbMqqbCorrectChosen","MlbMqqbCorrectChosen",200,0,500,200,0,300);
    TH2F h_MlbMqqbCorrectAll("MlbMqqbCorrectAll","MlbMqqbCorrectAll",200,0,500,200,0,300);
    TH2F h_MlbMqqbWrongOne("MlbMqqbWrongOne","MlbMqqbWrongOne",200,0,500,200,0,300);
    TH2F h_MlbMqqbWrongTwo("MlbMqqbWrongTwo","MlbMqqbWrongTwo",200,0,500,200,0,300);

    bool FalseEventContent = false;
    cout << " FalseEventContent : " << FalseEventContent << endl;
    TRootMCParticle *Top,*TopBar,*Bottom, *BottomBar,*Lepton,*NeutrinoMC,*WPlus,*WMinus,*Light,*LightBar;

    /////////////////////
    //  Used classes   //
    /////////////////////  
    BTagStudy bTagStudy;  //--> Should only be called before the event loop (otherwise the counters will not give the correct result)
    bTagStudy.InitializeBegin();
    MlbStudy mlbStudy;
    mlbStudy.initializeBegin();
    TFCreation tfCreation;
    tfCreation.InitializeVariables();

    ////////////////////////////////////
    //	loop on events
    ////////////////////////////////////
    
    nEvents[d] = 0;
    int itriggerSemiMu = -1,itriggerSemiEl = -1, previousRun = -1;
    if (verbose > 1)
      cout << "	Loop over events " << endl;

    for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++){
      //for (unsigned int ievt = 0; ievt < 2000000; ievt++){

      if(ievt > 20000) GenLHCOOutput = false;	

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
      if( dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 )
	{
	  jetTools->unCorrectMETTypeOne(init_jets, mets[0], true);
	  jetTools->correctJets(init_jets, event->kt6PFJets_rho(), true);
	  jetTools->correctMETTypeOne(init_jets, mets[0], true);
	}
      else
	{
	  jetTools->unCorrectMETTypeOne(init_jets, mets[0], false);
	  jetTools->correctJets(init_jets, event->kt6PFJets_rho(), false);
	  jetTools->correctMETTypeOne(init_jets, mets[0], false);
	}

      ////////////////////////////////////////
      //  Beam scraping and PU reweighting
      ////////////////////////////////////////
      double lumiWeight = 1;
      if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA"){    //Does this mean that it will be applied on data ? (why not with .find method ) ???
	// Apply the scraping veto. (Is it still needed?)
	bool isBeamBG = true;
	if(event->nTracks() > 10){
	  if( ( (float) event->nHighPurityTracks() ) / ( (float) event->nTracks() ) > 0.25 )
	    isBeamBG = false;
	}
	if(isBeamBG) continue;
      }
      else{
	lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() );
	if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0){   //Loop will never get here ...
	  lumiWeight=1;
	}
	//scaleFactor = scaleFactor*lumiWeight;  //Correct to do this immediately ?? (Is not the case in analyzer code Stijn ...)
      }
      //histo1D["lumiWeights"]->Fill(scaleFactor);	

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
      if(GenLHCOOutput == true) EventInfoFile << "     " << ievt << "         ";
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
	      EventInfoFile << "  0      0       0       1        " << NumberNegativeElectrons << "     ";
	    }
	  }//Negative electron
	  else if(Lepton->type() == 13){       //Looking at negative muon events (index 1 for LHCO file)
	    MadGraphId[4] = 2; //MadGraphId of mu = 2
	    MadGraphId[5] = 6; 
	    NumberNegativeMuons++;
	    if(GenLHCOOutput == true){
	      lhcoOutput.LHCOEventOutput(1, outFile[1], NumberNegativeMuons,LHCOVector,MadGraphId);
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
	  if(Lepton->type() == -11){            //Looking at positive electron events (index 2 for LHCO file)
	    MadGraphId[1] = 1; //MadGraphId of electron = 1
	    MadGraphId[2] = 6; 
	    NumberPositiveElectrons++;
	    if(GenLHCOOutput == true){
	      lhcoOutput.LHCOEventOutput(2, outFile[2], NumberPositiveElectrons,LHCOVector,MadGraphId);
	      EventInfoFile << "  0      0       1       0        " << NumberPositiveElectrons << "     ";
	    }
	  }//Positive electron
	  else if(Lepton->type() == -13){             //Looking at positive muon events (index 0 for LHCO file)
	    MadGraphId[1] = 2; //MadGraphId of muon = 2
	    MadGraphId[2] = 6; 
	    NumberPositiveMuons++;
	    if(GenLHCOOutput == true){
	      lhcoOutput.LHCOEventOutput(0, outFile[0], NumberPositiveMuons,LHCOVector,MadGraphId);
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
	h_StandardCosThetaNoEvtSel.Fill(standardCosTheta);
	
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

      // MSPlots before 'basic' event selection (no b-tag)
      if (MSPlot.find("Init_Events_pT_jet1_beforeEvtSel") == MSPlot.end()){
	MSPlot["Init_Events_pT_jet1_beforeEvtSel"] = new MultiSamplePlot(datasets, "Init_Events_pT_jet1_beforeEvtSel", 60, 0, 600, "p_{T} (GeV)");
	MSPlot["Init_Events_pT_jet2_beforeEvtSel"] = new MultiSamplePlot(datasets, "Init_Events_pT_jet2_beforeEvtSel", 60, 0, 600, "p_{T} (GeV)");
	MSPlot["Init_Events_pT_jet3_beforeEvtSel"] = new MultiSamplePlot(datasets, "Init_Events_pT_jet3_beforeEvtSel", 60, 0, 600, "p_{T} (GeV)");
	MSPlot["Init_Events_pT_jet4_beforeEvtSel"] = new MultiSamplePlot(datasets, "Init_Events_pT_jet4_beforeEvtSel", 60, 0, 600, "p_{T} (GeV)");
      }

      if(init_jets_corrected.size() >=4){
	MSPlot["Init_Events_pT_jet1_beforeEvtSel"]->Fill(init_jets_corrected[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
	MSPlot["Init_Events_pT_jet2_beforeEvtSel"]->Fill(init_jets_corrected[1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
	MSPlot["Init_Events_pT_jet3_beforeEvtSel"]->Fill(init_jets_corrected[2]->Pt(), datasets[d], true, Luminosity*scaleFactor);
	MSPlot["Init_Events_pT_jet4_beforeEvtSel"]->Fill(init_jets_corrected[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      }

      
      //Declare selection instance    
      Selection selection(init_jets_corrected, init_muons, init_electrons, mets, event->kt6PFJets_rho());
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
	MSPlot["Selected_Events_pT_jet1"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_jet1"+leptonFlav, 60, 0, 600, "p_{T} (GeV)");
	MSPlot["Selected_Events_pT_jet2"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_jet2"+leptonFlav, 60, 0, 600, "p_{T} (GeV)");
	MSPlot["Selected_Events_pT_jet3"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_jet3"+leptonFlav, 60, 0, 600, "p_{T} (GeV)");
	MSPlot["Selected_Events_pT_jet4"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_jet4"+leptonFlav, 60, 0, 600, "p_{T} (GeV)");
	MSPlot["Selected_Events_pT_4leadingjets"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_4leadingjets"+leptonFlav,60, 0, 600, "p_{T} (GeV)");
	MSPlot["Selected_Events_pT_alljets"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_alljets"+leptonFlav, 60, 0, 600, "p_{T} (GeV)");

	MSPlot["Selected_Events_pT_lepton"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_lepton"+leptonFlav,150,0,300,"p_{t} (GeV)");
      }

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
	  CorrectRecMassW=(*selectedJets[jetCombi[2]]+*selectedJets[jetCombi[3]]).M();
	  CorrectRecMassTop=(*selectedJets[jetCombi[2]]+*selectedJets[jetCombi[3]]+*selectedJets[jetCombi[1]]).M();
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
	    tfCreation.FillHistograms(&mcParticlesMatching[hadronicWJet1_.second], &mcParticlesMatching[hadronicWJet2_.second], &mcParticlesMatching[hadronicBJet_.second], &mcParticlesMatching[leptonicBJet_.second], Lepton, selectedJets[hadronicWJet1_.first], selectedJets[hadronicWJet2_.first], selectedJets[hadronicBJet_.first], selectedJets[leptonicBJet_.first], selectedLepton, eventselectedSemiMu, eventselectedSemiEl);

	    //Check the DeltaR vlaue between the different partons and reconstructed particles!:
	    if(eventselectedSemiMu == true){
	      histo1D["DeltaR_Ana_Muon"]->Fill( sqrt( pow((selectedLepton->Phi() - Lepton->Phi()),2) + pow((selectedLepton->Eta() - Lepton->Eta()),2) ) );
	      histo1D["genPt_Muon"]->Fill( Lepton->Pt() );
	      histo1D["recoPt_Muon"]->Fill(selectedLepton->Pt());
	    }
            if(eventselectedSemiEl == true){
	      histo1D["DeltaR_Ana_Elec"]->Fill(sqrt( pow((selectedLepton->Phi() - Lepton->Phi()),2) + pow((selectedLepton->Eta() - Lepton->Eta()),2) ) );
	      histo1D["genPt_Elec"]->Fill( Lepton->Pt() );
	      histo1D["recoPt_Elec"]->Fill( selectedLepton->Pt());
	    }
	    histo1D["DeltaR_Ana_Light1_mcPartMatch"]->Fill( sqrt( pow((selectedJets[hadronicWJet1_.first]->Phi() - mcParticlesMatching[hadronicWJet1_.second].Phi()),2) + pow((selectedJets[hadronicWJet1_.first]->Eta() - mcParticlesMatching[hadronicWJet1_.second].Eta()),2) ) );
            histo1D["DeltaR_Ana_Light1_mcPart"]->Fill( sqrt( pow((selectedJets[hadronicWJet1_.first]->Phi() - mcParticles[hadronicWJet1_.second]->Phi()),2) + pow((selectedJets[hadronicWJet1_.first]->Eta() - mcParticles[hadronicWJet1_.second]->Eta()),2) ) );
            histo1D["DeltaR_Ana_Light2_mcPartMatch"]->Fill( sqrt( pow((selectedJets[hadronicWJet2_.first]->Phi() - mcParticlesMatching[hadronicWJet2_.second].Phi()),2) + pow((selectedJets[hadronicWJet2_.first]->Eta() - mcParticlesMatching[hadronicWJet2_.second].Eta()),2) ) );
            histo1D["DeltaR_Ana_Light2_mcPart"]->Fill( sqrt( pow((selectedJets[hadronicWJet2_.first]->Phi() - mcParticles[hadronicWJet2_.second]->Phi()),2) + pow((selectedJets[hadronicWJet2_.first]->Eta() - mcParticles[hadronicWJet2_.second]->Eta()),2) ) );
	    histo1D["DeltaR_Ana_BHadr_mcPartMatch"]->Fill( sqrt( pow((selectedJets[hadronicBJet_.first]->Phi() - mcParticlesMatching[hadronicBJet_.second].Phi()),2) + pow((selectedJets[hadronicBJet_.first]->Eta() - mcParticlesMatching[hadronicBJet_.second].Eta()),2) ) );
            histo1D["DeltaR_Ana_BHadr_mcPart"]->Fill( sqrt( pow((selectedJets[hadronicBJet_.first]->Phi() - mcParticles[hadronicBJet_.second]->Phi()),2) + pow((selectedJets[hadronicBJet_.first]->Eta() - mcParticles[hadronicBJet_.second]->Eta()),2) ) );
            histo1D["DeltaR_Ana_BLept_mcPartMatch"]->Fill( sqrt( pow((selectedJets[leptonicBJet_.first]->Phi() - mcParticlesMatching[leptonicBJet_.second].Phi()),2) + pow((selectedJets[leptonicBJet_.first]->Eta() - mcParticlesMatching[leptonicBJet_.second].Eta()),2) ) );
            histo1D["DeltaR_Ana_BLept_mcPart"]->Fill( sqrt( pow((selectedJets[leptonicBJet_.first]->Phi() - mcParticles[leptonicBJet_.second]->Phi()),2) + pow((selectedJets[leptonicBJet_.first]->Eta() - mcParticles[leptonicBJet_.second]->Eta()),2) ) );


	  }//End of calculate Transfer Functions
        }//End of matched particles reconstructed
      }//if dataset Semi mu ttbar

      if(CalculateTF) continue;
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

      if (MSPlot.find("nSelectedJets_BeforeBTag"+leptonFlav) == MSPlot.end()){
	//MSPlots before and after #b-tagged and #light jets constraints
	MSPlot["nSelectedJets_BeforeBTag"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "nSelectedJets_BeforeBTag"+leptonFlav,11, -0.5, 10.5, "# selected jets");
	MSPlot["nSelectedJets_AfterBTag"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "nSelectedJets_AfterBTag"+leptonFlav,11, -0.5, 10.5, "# selected jets");
	MSPlot["nBTaggedJets_BeforeBTag"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "nBTaggedJets_BeforeBTag"+leptonFlav,11, -0.5, 10.5, "# b-tagged jets");
	MSPlot["nBTaggedJets_AfterBTag"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "nBTaggedJets_AfterBTag"+leptonFlav,11, -0.5, 10.5, "# b-tagged jets");
	MSPlot["nLightJets_BeforeBTag"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "nLightJets_BeforeBTag"+leptonFlav,11, -0.5, 10.5, "# light jets");
	MSPlot["nLightJets_AfterBTag"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "nLightJets_AfterBTag"+leptonFlav,11, -0.5, 10.5, "# light jets");
      }

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
	h_MlbMqqbCorrectAll.Fill(MqqbCorrect,MlbCorrect);
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
	    if(RecoLHCOOutput == true){
	      if(ConsideredCombi == 0) NumberNegRecoEl++;  //Only need to raise the eventNumber for one combination of the 4!!
	      if(verbosity > 4) cout << " Event : " << ievt << " with Number " << NumberNegRecoEl << " sent to LHCO Reco output (Negative electron) " << endl;
	      if(NumberNegRecoEl==1)outFileReco[3*4+ConsideredCombi].open(("MadWeightInput/AnalyzerOutput/TTbarSemiLepton_Reco_NegativeElectron_JetCombi"+JetCombiString+".lhco").c_str());
	      lhcoOutput.LHCOEventRecoOutput(3,outFileReco[3*4+ConsideredCombi], NumberNegRecoEl, LHCORecoVector, MadGraphRecoId);
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
	      lhcoOutput.LHCOEventRecoOutput(1, outFileReco[1*4+ConsideredCombi], NumberNegRecoMu, LHCORecoVector, MadGraphRecoId);
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
	  if(eventselectedSemiEl){//Positive electron
	    MadGraphRecoId[1] = 1;
	    MadGraphRecoId[2] = 6;
	    if(RecoLHCOOutput == true){ 
	      if(ConsideredCombi == 0) NumberPosRecoEl++;
	      if(verbosity > 4) cout << " Event : " << ievt << " with Number " << NumberPosRecoEl << " sent to LHCO Reco output (Positive electron) " << endl;
	      if(NumberPosRecoEl==1)outFileReco[2*4+ConsideredCombi].open(("MadWeightInput/AnalyzerOutput/TTbarSemiLepton_Reco_PositiveElectron_JetCombi"+JetCombiString+".lhco").c_str());
	      lhcoOutput.LHCOEventRecoOutput(2,outFileReco[2*4+ConsideredCombi], NumberPosRecoEl, LHCORecoVector, MadGraphRecoId);	 
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
	      lhcoOutput.LHCOEventRecoOutput(0, outFileReco[0*4+ConsideredCombi], NumberPosRecoMu, LHCORecoVector, MadGraphRecoId);
	      if(ConsideredCombi == 0) EventInfoFile << "     " << NumberPosRecoMu << endl;
	    }	
	  }
	}//End of positive lepton
      }//End of loop over the different jet combinations  
    } //loop on events

    // -------- Calculate TF MadWeight  --------//
    if(CalculateTF){
      tfCreation.CalculateTF(true, false, true, false, true); //bool drawHistos, bool writeTF, bool doFits, bool useROOTClass, bool useStartValues);
    }
    else{	
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
      std::string OptionName5Jets[6];
      std::string OptionName4Jets[6];
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

      //Save the histograms belonging to the mlb output information! 
      std::string Title[3] = {"5Jets","4Jets","Pure5Jets"};
      std::string Name[3] =  {" - 5 jets case) "," - 4 jets case) "," - pure 5 jets case) "};
 
      TH1F *h_ChiSqCorrect[3], *h_ChiSqCorrectFound[3], *h_ChiSqMinimum[3],* h_ChiSqNotMinimum[3], *h_ChiSqWrong[3];
      TH1F *h_ChiSqCorrectWhenMatched[3], *h_ChiSqMinimumWhenMatched[3], *h_ChiSqNotMinimumWhenMatched[3], *h_ChiSqAllWhenNotMatched[3], *h_ChiSqMinimumWhenCorrect[3], *h_ChiSqMinimumWhenWrong[3], *h_ChiSqDiffWhenWrong[3]; 
      for(int ii = 0; ii < 3; ii++){
	h_ChiSqCorrect[ii]     =new TH1F(("ChiSqCorrect"+Title[ii]).c_str(),     ("#chi^{2} distribution for the correct combination (all events"+Name[ii]).c_str() ,              150,0,50);
	h_ChiSqCorrectFound[ii]=new TH1F(("ChiSqCorrectFound"+Title[ii]).c_str(),("#chi^{2} distribution for the correct combination (when found"+Name[ii]).c_str(),               150,0,50);
	h_ChiSqMinimum[ii]     =new TH1F(("ChiSqMinimum"+Title[ii]).c_str(),     ("#chi^{2} distribution of the minimal combination considered (all events"+Name[ii]).c_str(),     150,0,50);
	h_ChiSqNotMinimum[ii]  =new TH1F(("ChiSqNotMinimum"+Title[ii]).c_str(),  ("#chi^{2} distribution of the non-minimal combinations considered (all events"+Name[ii]).c_str(),150,0,50);
	h_ChiSqWrong[ii]       =new TH1F(("ChiSqWrong"+Title[ii]).c_str(),       ("#chi^{2} distribution of the non-correct combinations considered (all events"+Name[ii]).c_str(),150,0,50);
     
	h_ChiSqCorrectWhenMatched[ii]   =new TH1F(("ChiSqCorrectWhenMatched"+Title[ii]).c_str(),   ("#chi^{2} distribution for the correct combination (matched events only"+Name[ii]).c_str(),     150,0,50);
	h_ChiSqMinimumWhenMatched[ii]   =new TH1F(("ChiSqMinimumWhenMatched"+Title[ii]).c_str(),   ("#chi^{2} distribution for the minimal combination (matched events only"+Name[ii]).c_str(),     150,0,50);
	h_ChiSqNotMinimumWhenMatched[ii]=new TH1F(("ChiSqNotMinimumWhenMatched"+Title[ii]).c_str(),("#chi^{2} distribution for the non-minimal combinations (matched events only"+Name[ii]).c_str(),150,0,50);
	h_ChiSqMinimumWhenCorrect[ii]   =new TH1F(("ChiSqMinimumWhenCorrect"+Title[ii]).c_str(),   ("#chi^{2} distribution for the minimal combination (correct choice only"+Name[ii]).c_str(),     150,0,50);
	h_ChiSqMinimumWhenWrong[ii]     =new TH1F(("ChiSqMinimumWhenWrong"+Title[ii]).c_str(),     ("#chi^{2} distribution for the minimal combination (wrong choice only"+Name[ii]).c_str(),       150,0,50);
	h_ChiSqAllWhenNotMatched[ii]    =new TH1F(("ChiSqAllWhenNotMatched"+Title[ii]).c_str(),    ("#chi^{2} distribution for all the combinations (non-matched evens only"+Name[ii]).c_str(),     150,0,50);

	h_ChiSqDiffWhenWrong[ii] = new TH1F(("ChiSqDiffWhenWrong"+Title[ii]).c_str(), ("#chi^{2}_{correct} - #chi^{2}_{minimum} for the wrong choice"+Name[ii]).c_str(),500,-5,15);
      }

      for(int jetCase = ChosenBTagOption; jetCase <(ChosenBTagOption+3); jetCase++){
	//In mlb Class the ChosenBTagOption, ChosenBTagOption+1 and ChosenBTagOption+2 indices are filled (from the 6 options).
	//However to reduce the number of TH1F's made, this is translated to 0,1 and 2 here in the analyzer!!
	for(int ii = 0; ii < (mlbStudy.getChiSqCorrect(jetCase)).size();     ii++) h_ChiSqCorrect[jetCase-ChosenBTagOption]->Fill((mlbStudy.getChiSqCorrect(jetCase))[ii]);
	for(int ii = 0; ii < (mlbStudy.getChiSqCorrectFound(jetCase)).size();ii++) h_ChiSqCorrectFound[jetCase-ChosenBTagOption]->Fill((mlbStudy.getChiSqCorrectFound(jetCase))[ii]);
	for(int ii = 0; ii < (mlbStudy.getChiSqMinimum(jetCase)).size();     ii++) h_ChiSqMinimum[jetCase-ChosenBTagOption]->Fill((mlbStudy.getChiSqMinimum(jetCase))[ii]);
	for(int ii = 0; ii < (mlbStudy.getChiSqNotMinimum(jetCase)).size();  ii++) h_ChiSqNotMinimum[jetCase-ChosenBTagOption]->Fill((mlbStudy.getChiSqNotMinimum(jetCase))[ii]);
	for(int ii = 0; ii < (mlbStudy.getChiSqWrong(jetCase)).size();       ii++) h_ChiSqWrong[jetCase-ChosenBTagOption]->Fill((mlbStudy.getChiSqWrong(jetCase))[ii]);

	for(int ii = 0; ii <(mlbStudy.getChiSqCorrectWhenMatched(jetCase)).size();   ii++) h_ChiSqCorrectWhenMatched[jetCase-ChosenBTagOption]->Fill((mlbStudy.getChiSqCorrectWhenMatched(jetCase))[ii]);
	for(int ii = 0; ii <(mlbStudy.getChiSqMinimumWhenMatched(jetCase)).size();   ii++) h_ChiSqMinimumWhenMatched[jetCase-ChosenBTagOption]->Fill((mlbStudy.getChiSqMinimumWhenMatched(jetCase))[ii]);
	for(int ii = 0; ii <(mlbStudy.getChiSqNotMinimumWhenMatched(jetCase)).size();ii++) h_ChiSqNotMinimumWhenMatched[jetCase-ChosenBTagOption]->Fill((mlbStudy.getChiSqNotMinimumWhenMatched(jetCase))[ii]);
	for(int ii = 0; ii <(mlbStudy.getChiSqMinimumWhenCorrect(jetCase)).size();   ii++) h_ChiSqMinimumWhenCorrect[jetCase-ChosenBTagOption]->Fill((mlbStudy.getChiSqMinimumWhenCorrect(jetCase))[ii]);
	for(int ii = 0; ii < (mlbStudy.getChiSqMinimumWhenWrong(jetCase)).size();    ii++) h_ChiSqMinimumWhenWrong[jetCase-ChosenBTagOption]->Fill((mlbStudy.getChiSqMinimumWhenWrong(jetCase))[ii]);
	for(int ii = 0; ii < (mlbStudy.getChiSqAllWhenNotMatched(jetCase)).size();   ii++) h_ChiSqAllWhenNotMatched[jetCase-ChosenBTagOption]->Fill((mlbStudy.getChiSqAllWhenNotMatched(jetCase))[ii]); 

	for(int ii = 0; ii < (mlbStudy.getChiSqDiffWhenWrong(jetCase)).size(); ii++) h_ChiSqDiffWhenWrong[jetCase-ChosenBTagOption]->Fill((mlbStudy.getChiSqDiffWhenWrong(jetCase))[ii]);
      }   

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
	h_ChiSqMinimumWhenCorrect[ii]->Write();
	h_ChiSqMinimumWhenWrong[ii]->Write();
	h_ChiSqAllWhenNotMatched[ii]->Write();

	h_ChiSqDiffWhenWrong[ii]->Write();
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

      h_MlbMass.Write();
      h_MqqbMass.Write();
    
      h_CosThetaReco.Write();
      h_NeutrinoEta.Write();

      h_MlbMqqbCorrectAll.Write();
      h_MlbMqqbCorrectChosen.Write();
      h_MlbMqqbWrongOne.Write();
      h_MlbMqqbWrongTwo.Write();
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
    temp->Draw(name); //, true, true, true, true, true, 1, false);
    temp->Write(fout, name, true, "PlotsMacro/MSPlots/");
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
