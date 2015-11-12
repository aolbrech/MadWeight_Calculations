//#include "TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"   //Needed to load TRootMCParticle & TRootJet, which is used in TFCreation.h
#include "TTree.h"
#include "TBranch.h"
#include "TNtuple.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include "TStyle.h"
#include "TDirectory.h"
#include "TClonesArray.h"
#include "TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"
#include "TopTreeAnalysisBase/Tools/interface/LeptonTools.h"
#include <sys/stat.h>  //Needed for mkdir option

#include "TLorentzVector.h"

#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "PersonalClasses/Style.C"                                                 //CHECK if this works!

#include "AnomalousCouplings/PersonalClasses/interface/AnomCoupLight.h"
#include "AnomalousCouplings/PersonalClasses/interface/BTagStudy_OLD.h"
#include "AnomalousCouplings/PersonalClasses/interface/BTagStudy.h"
#include "AnomalousCouplings/PersonalClasses/interface/LHCOOutput.h"
#include "AnomalousCouplings/PersonalClasses/interface/ExtraEvtSelCuts.h"

using namespace std;
using namespace reweight;   //Need this in order to use LumiReWeighting!

//Convert integer to string!
template <typename T>
string NumberToString ( T Number ){
    stringstream ss;
    ss << Number;
    return ss.str();
}

//Add a timestamp!
std::string timestamp(){
    time_t ltime;
    struct tm *Tm;
 
    ltime=time(NULL);
    Tm=localtime(&ltime);

    std::string TIME = (NumberToString(Tm->tm_mday)+NumberToString(Tm->tm_mon+1)+NumberToString(Tm->tm_year+1900)+"_"+NumberToString(Tm->tm_hour)+NumberToString(Tm->tm_min)+NumberToString(Tm->tm_sec)).c_str();
    return TIME;
}

int main (int argc, char *argv[])
{
  clock_t start = clock();
  cout << "************************************************************" << endl;
  cout << "  Beginning of the program for analyzing the Light Trees !  " << endl;
  cout << "************************************************************" << endl;

  //Get the timestamp to add this in the output root file!
  std::string time = timestamp();

  //SetStyle if needed
  //setTDRStyle(); 
  setMyStyle();

  //--------------------------------//
  //  Get NrEvts from command line  //
  //--------------------------------//
  int nrEvtsComLine = -1;
  if(argc >= 2 && string(argv[1]) != "-1")
    nrEvtsComLine = atoi(argv[1]);

  //------------------------//
  //  Verbosity for output  //
  //  -->Used in bTagStudy  //
  //------------------------//
  int verbosity                 = 1;
  //1 Main Info
  //2 mcParticlesMatchin Info
  //3 
  //4 Info for each event
  //5 Debug

  //---------------------------//
  //  Run specific parts only  //
  //---------------------------//
  bool getLHCOOutput = true; 
  bool splitLeptonChargeLHCO = false;
  bool getCorrectAndWrongLHCO = true; 
  bool savePDF = false;
  bool bTagChoiceMade = true;  
  bool getMassFits = false;

  //-- Specific b-tag stuff!
  int ChosenBTag = 3;  //-->Corresponds to two T b-tags and no light veto!
  int NrBTags = 6;
  string bTitle[6] = {"LooseTags","MediumTags","MediumTagsLVeto","TightTags","TightTagsMVeto","TightTagsLVeto"};
  string numberBTags = "AllBTags";
  if(bTagChoiceMade){ NrBTags = 1; bTitle[0] = bTitle[ChosenBTag]; numberBTags = bTitle[ChosenBTag];}

  //-------------------------//
  //   Which systematics ?   //
  //-------------------------//
  std::string doLeptonSFShift = "Nominal";  //Other options are "Minus" and "Plus" (Capital letters are needed!)
  std::string doLumiWeightShift = "Nominal";  //Other options are "Minus" and "Plus"

  //----------------------------//
  // Input & output information //
  //----------------------------//
  //ROOT file for storing plots
  string pathPNG = "PlotsMacro_Light";
  TFile* outputFile = new TFile((pathPNG+"/AnomCoup_Analysis_"+numberBTags+"_"+time+".root").c_str(),"RECREATE");

  //Which datasets should be considered
  vector<string> inputFiles;
  vector<Dataset*> datasets;
  //inputFiles.push_back("LightTree/AnomCoupLight_TTbarJets_SemiLept_AllTTbarEvents_19Aug2015.root");
  std::string inputFileDir = "/user/aolbrech/PBS_ScriptRunning/Results/RESULTS_AnomCoup_08112015_200540/";
  inputFiles.push_back((inputFileDir+"AnomCoupLight_Data_Mu_Merged_22Jan2013_Nominal.root").c_str());                                                    //Summer13_v4
  //inputFiles.push_back("/user/aolbrech/PBS_ScriptRunning/Results/RESULTS_AnomCoup_10112015_144641/AnomCoupLight_Data_Mu_Merged_22Jan2013_Nominal.root"); //Winter14_v5
  //inputFiles.push_back("/user/aolbrech/PBS_ScriptRunning/Results/RESULTS_AnomCoup_10112015_210750/AnomCoupLight_Data_Mu_Merged_22Jan2013_Nominal.root"); //Winter14_v8
  inputFiles.push_back((inputFileDir+"AnomCoupLight_WJets_1jets_Nominal.root").c_str());
  inputFiles.push_back((inputFileDir+"AnomCoupLight_WJets_2jets_Nominal.root").c_str());
  inputFiles.push_back((inputFileDir+"AnomCoupLight_WJets_3jets_Nominal.root").c_str());
  inputFiles.push_back((inputFileDir+"AnomCoupLight_WJets_4jets_Nominal.root").c_str());
  inputFiles.push_back((inputFileDir+"AnomCoupLight_ZJets_1jets_Nominal.root").c_str());
  inputFiles.push_back((inputFileDir+"AnomCoupLight_ZJets_2jets_Nominal.root").c_str());
  inputFiles.push_back((inputFileDir+"AnomCoupLight_ZJets_3jets_Nominal.root").c_str());
  inputFiles.push_back((inputFileDir+"AnomCoupLight_ZJets_4jets_Nominal.root").c_str());
  inputFiles.push_back((inputFileDir+"AnomCoupLight_TTbarJets_SemiLept_Nominal.root").c_str());
  inputFiles.push_back((inputFileDir+"AnomCoupLight_TTbarJets_FullLept_Nominal.root").c_str());
  inputFiles.push_back((inputFileDir+"AnomCoupLight_TTbarJets_FullHadr_Nominal.root").c_str());
  if(verbosity > 0) std::cout << " - All ROOT files loaded " << std::endl;
	
  //-------------------------//
  // Set analysis luminosity //
  //-------------------------//
  float Luminosity = 19646.840; //IsoMu24_eta2p1 trigger for semiMu case!  
  if(verbosity > 0) std::cout << " - Analysis luminosity is set to : " << Luminosity << std::endl;

  //-------------------//
  // Load all datasets //
  //-------------------//
  for(unsigned int iDataSet=0; iDataSet<inputFiles.size(); iDataSet++){
    TFile* inputFile = new TFile(inputFiles[iDataSet].c_str(),"READ");

    TTree* inConfigTree = (TTree*) inputFile->Get("configTreeLightFile");
    TBranch* d_br = (TBranch*) inConfigTree->GetBranch("Dataset");
    TClonesArray* tc_dataset = new TClonesArray("Dataset",0);
    d_br->SetAddress(&tc_dataset);
    inConfigTree->GetEvent(0);
    Dataset* dataSet = (Dataset*) tc_dataset->At(0);
    //XSection = dataSet->Xsection();

    int color = 0;
    if(dataSet->Name().find("QCD") == 0)             color = kYellow;
    if(dataSet->Name().find("TT") == 0){             color = kRed+1; dataSet->SetTitle("t#bar{t}");}        //Giving the different ttbar samples the same title will merge them in the MSPlot!
    if(dataSet->Name().find("TTbarJets_Other") == 0) color = kRed-7;
    if(dataSet->Name().find("ST") == 0 || 
       dataSet->Name().find("SingleTop") == 0)       color = kMagenta;
    if(dataSet->Name().find("WJets") == 0){          color = kGreen-3; dataSet->SetTitle("W#rightarrowl#nu");}      
    if(dataSet->Name().find("ZJets") == 0){          color = kAzure-2; dataSet->SetTitle("Z/#gamma*#rightarrowl^{+}l^{-}");}
    if(dataSet->Name().find("Data_Mu") == 0){        color = kBlack; dataSet->SetTitle("Data"); dataSet->SetEquivalentLuminosity(Luminosity);}

    Dataset* tmpDS = new Dataset(dataSet->Name(), dataSet->Title(), dataSet->DoIt(), color, dataSet->LineStyle(), dataSet->LineWidth(), dataSet->NormFactor(), dataSet->Xsection());
    tmpDS->SetEquivalentLuminosity( dataSet->EquivalentLumi() );
    datasets.push_back( tmpDS );
  }
  if(verbosity > 0) std::cout << " - Datasets initialized " << std::endl;

  //--------------------------//
  //  Initialize MSPlots      //
  //  --> Needs dataset info  //
  //--------------------------//
  map<string,MultiSamplePlot*> MSPlot;
  string leptFlavs[2]={"_mu","_el"};

  //-- Histograms which should be made separately for the two possible decay channels
  for(int ii = 0; ii < 2; ii++){
    string leptFlav = leptFlavs[ii];

    //-- Histograms which should contain the considered btag name!
    for(int ibTag = 0; ibTag < NrBTags; ibTag++){
      MSPlot["nBTaggedJets_"+ bTitle[ibTag]+"_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "nBTaggedJets_"+ bTitle[ibTag]+"_BeforeBTag"+leptFlav,10, -0.5, 9.5, "# b-tagged jets");
      MSPlot["nLightJets_"+   bTitle[ibTag]+"_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "nLightJets_"+   bTitle[ibTag]+"_BeforeBTag"+leptFlav,10, -0.5, 9.5, "# light jets");
      MSPlot["nSelectedJets_"+bTitle[ibTag]+"_AfterBTag"+ leptFlav] = new MultiSamplePlot(datasets, "nSelectedJets_"+bTitle[ibTag]+"_AfterBTag"+ leptFlav,10, -0.5, 9.5, "# selected jets");
      MSPlot["nBTaggedJets_"+ bTitle[ibTag]+"_AfterBTag"+ leptFlav] = new MultiSamplePlot(datasets, "nBTaggedJets_"+ bTitle[ibTag]+"_AfterBTag"+ leptFlav,10, -0.5, 9.5, "# b-tagged jets");
      MSPlot["nLightJets_"+   bTitle[ibTag]+"_AfterBTag"+ leptFlav] = new MultiSamplePlot(datasets, "nLightJets_"+   bTitle[ibTag]+"_AfterBTag"+ leptFlav,10, -0.5, 9.5, "# light jets");
    }
      
    MSPlot["nSelectedJets_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "nSelectedJets_BeforeBTag"+leptFlav,10, -0.5, 9.5, "# selected jets");

    MSPlot["LeptonPt"+leptFlav] = new MultiSamplePlot(datasets, "LeptonPt"+leptFlav, 150, 0, 200,"Lepton p_{T} (GeV)");
    MSPlot["LeptonPt_BeforePU"+leptFlav] = new MultiSamplePlot(datasets, "LeptonPt_BeforePU"+leptFlav, 150, 0, 200,"Lepton p_{T} (GeV)");
    MSPlot["LeptonPt_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "LeptonPt_BeforeBTag"+leptFlav, 150, 0, 200,"Lepton p_{T} (GeV)");
    MSPlot["LeptonEta"+leptFlav] = new MultiSamplePlot(datasets, "LeptonEta"+leptFlav, 100, -2.6, 2.6,"Lepton #eta");
    MSPlot["LeptonCharge"+leptFlav] = new MultiSamplePlot(datasets, "LeptonCharge"+leptFlav, 150, 0, 200,"Lepton charge");
    MSPlot["JetPt_LeadingJet"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_LeadingJet"+leptFlav, 150, 0, 250,"Jet p_{T} (GeV)");
    MSPlot["JetPt_LeadingJet_BeforePU"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_LeadingJet_BeforePU"+leptFlav, 150, 0, 250,"Jet p_{T} (GeV)");
    MSPlot["JetPt_LeadingJet_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_LeadingJet_BeforeBTag"+leptFlav, 150, 0, 250,"Jet p_{T} (GeV)");
    MSPlot["JetEta_LeadingJet"+leptFlav] = new MultiSamplePlot(datasets, "JetEta_LeadingJet"+leptFlav, 150, 0, 200,"Jet #eta");
    MSPlot["nPV"+leptFlav] = new MultiSamplePlot(datasets, "nPV"+leptFlav, 50, 0, 50,"Number of primary vertices");
    MSPlot["nPV_BeforePU"+leptFlav] = new MultiSamplePlot(datasets, "nPV_BeforePU"+leptFlav, 50, 0, 50,"Number of primary vertices");
    MSPlot["nPV_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "nPV_BeforeBTag"+leptFlav, 50, 0, 50,"Number of primary vertices");
  }

  //-----------------------//
  // Initialize histograms //
  //-----------------------//
  map<string,TH1F*> histo1D;
  map<string,TH2F*> histo2D;  

  histo1D["Mlb_CorrTT"] = new TH1F("Mlb_CorrTT", "mass_{l,b} for the actual semi-leptonic ttbar events", 200, 0, 250);
  histo1D["HadrMTop_CorrTT"] = new TH1F("HadrMTop_CorrTT", "Hadronic m_{t} for the actual semi-leptonic ttbar events", 200, 50, 350);
  histo1D["HadrMW_CorrTT"] = new TH1F("HadrMW_CorrTT", "Hadronic m_{W} for the actual semi-leptonic ttbar events", 200, 0, 250);
  histo2D["MW_vs_MTop_CorrTT"] = new TH2F("MW_vs_MTop_CorrTT","Hadronic m_{W} versus m_{t} for the actual semi-leptonic ttbar jet combinations",200,50,350,200,0,250);

  //histo1D["CosTheta_SelEvts"] = new TH1F("CosTheta_SelectedEvents","Cos #theta^{*} distribution for selected RECO events",200, -1, 1);
  histo1D["CosTheta_Gen_SelEvts_Mu"] = new TH1F("CosTheta_Gen_SelEvts_Mu","Cos #theta^{*}_{gen} distribution for the selected events (muon channel)",200,-1,1);
  histo1D["CosTheta_Gen_SelEvts_El"] = new TH1F("CosTheta_Gen_SelEvts_El","Cos #theta^{*}_{gen} distribution for the selected events (elec channel)",200,-1,1);
  histo1D["CosTheta_Gen_SelEvts"] = new TH1F("CosTheta_Gen_SelEvts","Cos #theta^{*}_{gen} distribution for the selected events",200,-1,1);

  for(int ibTag = 0; ibTag < NrBTags; ibTag++){
    histo1D["Mlb_CorrTT_"+bTitle[ibTag]]        = new TH1F(("Mlb_CorrTT_"+bTitle[ibTag]).c_str(),  ("mass_{l,b} for the actual semi-leptonic ttbar events -- "+bTitle[ibTag]).c_str(), 200, 0, 250);
    histo1D["HadrMTop_CorrTT_"+bTitle[ibTag]]   = new TH1F(("HadrMTop_CorrTT_"+bTitle[ibTag]).c_str(), ("Hadronic m_{t} for the actual semi-leptonic ttbar events -- "+bTitle[ibTag]).c_str(), 200, 50, 350);
    histo1D["HadrMW_CorrTT_"+bTitle[ibTag]]     = new TH1F(("HadrMW_CorrTT_"+bTitle[ibTag]).c_str(),   ("Hadronic m_{W} for the actual semi-leptonic ttbar events -- "+bTitle[ibTag]).c_str(), 200, 0, 250);
    histo2D["MW_vs_MTop_CorrTT_"+bTitle[ibTag]] = new TH2F(("MW_vs_MTop_CorrTT_"+bTitle[ibTag]).c_str(), ("Hadronic m_{W} versus m_{t} for actual semi-lept ttbar events -- "+bTitle[ibTag]).c_str(), 200, 50, 350, 200, 0,250);
  }

  //--------------------------------//
  // Lumi reweighting and lepton SF //
  //--------------------------------//
  LumiReWeighting LumiWeights;
  LumiReWeighting LumiWeightsUp;
  LumiReWeighting LumiWeightsDown;
  LumiWeights = LumiReWeighting("PersonalClasses/Calibrations/PUReweighting/pileup_MC_Summer12_S10.root", "PersonalClasses/Calibrations/PUReweighting/pileup_2012Data53X_UpToRun208686_Mu/nominal.root", "pileup", "pileup");
  LumiWeightsUp = LumiReWeighting("PersonalClasses/Calibrations/PUReweighting/pileup_MC_Summer12_S10.root", "PersonalClasses/Calibrations/PUReweighting/pileup_2012Data53X_UpToRun208686_Mu/sys_up.root", "pileup", "pileup");
  LumiWeightsDown = LumiReWeighting("PersonalClasses/Calibrations/PUReweighting/pileup_MC_Summer12_S10.root", "PersonalClasses/Calibrations/PUReweighting/pileup_2012Data53X_UpToRun208686_Mu/sys_down.root", "pileup", "pileup");
  cout << " - LumiReWeighting instantiated ... " << endl;
  
  // initialize lepton SF (ROOT files taken from: IsoMu24_eta2p1)
  LeptonTools* leptonTools = new LeptonTools(false);
  //leptonTools->readMuonSF("PersonalClasses/Calibrations/LeptonSF/Muon_ID_iso_Efficiencies_Run_2012ABCD_53X.root", "PersonalClasses/Calibrations/LeptonSF/MuonEfficiencies_Run_2012A_2012B_53X.root", "PersonalClasses/Calibrations/LeptonSF/MuonEfficiencies_Run_2012C_53X.root", "PersonalClasses/Calibrations/LeptonSF/TriggerMuonEfficiencies_Run_2012D_53X.root"); //Order: ID/Iso/Trig files!
  leptonTools->readMuonSF("PersonalClasses/Calibrations/LeptonSF/MuonEfficiencies_Run2012ReReco_53X.root","PersonalClasses/Calibrations/LeptonSF/MuonEfficiencies_ISO_Run_2012ReReco_53X.root","PersonalClasses/Calibrations/LeptonSF/SingleMuonTriggerEfficiencies_eta2p1_Run2012ABCD_v5trees.root");
  leptonTools->readElectronSF();

  //-----------------------//
  // Load personal classes //
  //-----------------------//
  float Mlb  = 108.1841, S_Mlb  = 31.4213;
  float Mqqb = 174.6736, S_Mqqb = 17.5757;
  float MW = 83.8037, S_MW = 10.2385;
//  ExtraEvtSelCuts extraEvtSelCuts(Mqqb, S_Mqqb, MW, S_MW, bTagChoiceMade, 3, 2);
  BTagStudy bTagStudy(verbosity, datasets, bTagChoiceMade, ChosenBTag, Mlb, S_Mlb, Mqqb, S_Mqqb);

  //--------------------------------------//
  // Loop on datasets for actual analysis //
  //--------------------------------------//
  for(unsigned int iDataSet = 0; iDataSet < inputFiles.size(); iDataSet++){

    TFile* inputFile = new TFile(inputFiles[iDataSet].c_str(),"READ");

    TTree* inLightTree = (TTree*) inputFile->Get("LightTree");
    TBranch* light_br = (TBranch*) inLightTree->GetBranch("TheAnomCoupLight");
    AnomCoupLight* light = 0;
    light_br->SetAddress(&light);
  
    //Number of events that will be used in the "loop on events"
    int nrEvts = 0;
    if(nrEvtsComLine == -1) nrEvts = inLightTree->GetEntries();
    else                    nrEvts = nrEvtsComLine;
  
    Dataset* dataSet = datasets[iDataSet];//(Dataset*) tc_dataset->At(0);
    string dataSetName = dataSet->Name();
    if(verbosity > 0) std::cout << "   *** Looking at dataset "<< dataSetName << " (" << iDataSet+1 << "/" << inputFiles.size() << ") with " << nrEvts << " events (total = " << inLightTree->GetEntries() << ") ! " << std::endl;

    //-----------------------//
    // Load personal classes //
    //-----------------------//
    bTagStudy.InitializeDataSet(dataSetName);
    LHCOOutput lhcoOutput(verbosity, getLHCOOutput, splitLeptonChargeLHCO, getCorrectAndWrongLHCO);
//    extraEvtSelCuts.Initialize(bTitle[0], dataSetName);
    lhcoOutput.Initialize("Reco", dataSetName);

    ofstream EvtNrMatching;
    if(bTagChoiceMade && getLHCOOutput){
      EvtNrMatching.open(("MadWeightInput/AnalyzerOutput/EventNrMatching_"+dataSetName+".txt").c_str());
      EvtNrMatching << "  Event Nr     Extra cuts survived      Gen cos theta*        Main lhco file      MW Number (main)      TTbar splitting lhco file      MW Number (splitting)   " << endl;
    }

    // --------------------------- //
    //  Start looping over events  //
    // --------------------------- //    
    int nSelectedMu = 0, nSelectedEl = 0;
    for(unsigned int iEvt = 0; iEvt < nrEvts; iEvt++){
      inLightTree->GetEvent(iEvt);
      if(iEvt%5000 == 0)
	std::cout<<"    Processing the "<<iEvt<<"th event ("<< ((double)iEvt/(double)nrEvts)*100<<"%)"<<" -> # selected: "<<nSelectedMu<<" (mu+jets) "<<nSelectedEl<<" (e+jets)"<< flush<<"\r";

      //*** Start with corrections ***//
      // Beam scraping and PU reweighting
      double lumiWeight = 1;  
      if(! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) ){   
        if(doLumiWeightShift == "Nominal")    lumiWeight = LumiWeights.ITweight( (int)light->nTruePU() );
        else if(doLumiWeightShift == "Plus")  lumiWeight = LumiWeightsUp.ITweight( (int)light->nTruePU() );
        else if(doLumiWeightShift == "Minus") lumiWeight = LumiWeightsDown.ITweight( (int)light->nTruePU() );
      }

      int eventId = light->eventID();
      int runId = light->runID();
      int nTruePU = light->nTruePU();
      int nPrimVertices = light->nPV();
      float scaleFactor = light->scaleFactor();

      vector<float> bTagCSV = light->CSVbTag();
      vector<TLorentzVector> selJets = light->selectedJets();
      TLorentzVector selLepton = light->selectedLepton();
      TLorentzVector MET = light->met();
      int decayCh = light->decayChannel();  //0 = semiMu and 1 = semiEl
      std::string leptChannel;
      if(decayCh == 0) leptChannel = "_mu";
      else if(decayCh == 1) leptChannel = "_el";
      float leptCharge = light->leptonCharge();
      vector<int> correctJetCombi = light->correctJetCombi();    //0 = LeptB, 1 = HadrB, 2 = Quark1 & 3 = Quark2
      float genCosTheta = light->genCosTh();

      MSPlot["nPV_BeforePU"+leptChannel]->Fill( nPrimVertices, datasets[iDataSet], true, Luminosity*scaleFactor);
      MSPlot["JetPt_LeadingJet_BeforePU"+leptChannel]->Fill( selJets[0].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
      MSPlot["LeptonPt_BeforePU"+leptChannel]->Fill( selLepton.Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);

      //ooooooooOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOooooooooooooOOOOO
      //ooOOooOOoo      Reading out nTuples done           ooOOooOOoo
      //ooOOooOOoo-----------------------------------------ooOOooOOoo
      //ooOOooOOoo      Start of actual analysis           ooOOooOOoo
      //ooooooooOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOooooooooooooOOOOO
      bTagStudy.CalculateJets(selJets, bTagCSV, correctJetCombi, selLepton, datasets[iDataSet], Luminosity*scaleFactor);
 
      //MSPlots with number of jets information before requiring at least two b-jets and at least 2 light jets!
      MSPlot["nSelectedJets_BeforeBTag"+leptChannel]->Fill( selJets.size(),                    datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight);

      MSPlot["nPV_BeforeBTag"+leptChannel]->Fill( nPrimVertices, datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight);
      MSPlot["JetPt_LeadingJet_BeforeBTag"+leptChannel]->Fill( selJets[0].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight);
      MSPlot["LeptonPt_BeforeBTag"+leptChannel]->Fill( selLepton.Pt(), datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight);

      //Store the expected Mlb, Mqqb, Mt and MW distributions
      if(dataSetName.find("TTbarJets_SemiLept") == 0){
        if(correctJetCombi[0] != 9999 && correctJetCombi[1] != 9999 && correctJetCombi[2] != 9999 && correctJetCombi[3] != 9999){
          histo1D["Mlb_CorrTT"]->Fill( (selLepton+selJets[correctJetCombi[0]]).M());
          histo1D["HadrMTop_CorrTT"]->Fill( (selJets[correctJetCombi[1]]+selJets[correctJetCombi[2]]+selJets[correctJetCombi[3]]).M() );
          histo1D["HadrMW_CorrTT"]->Fill( (selJets[correctJetCombi[2]]+selJets[correctJetCombi[3]]).M());
          histo2D["MW_vs_MTop_CorrTT"]->Fill( (selJets[correctJetCombi[1]]+selJets[correctJetCombi[2]]+selJets[correctJetCombi[3]]).M(), (selJets[correctJetCombi[2]]+selJets[correctJetCombi[3]]).M() );
        }
      }
 
      for(int ibTag = 0; ibTag < NrBTags; ibTag++){
        vector<int> selJetCombi = bTagStudy.getIndices(ibTag);

        //--- Check whether the event is correctly reconstructed  ---// 
        //---  (jetCombi is initialized to 9999 for all dataSets) ---//
        int CWUIndex = 999;
        if(correctJetCombi[0] != 9999 && correctJetCombi[1] != 9999 && correctJetCombi[2] != 9999 && correctJetCombi[3] != 9999){
          if( selJetCombi[0] == correctJetCombi[0] && selJetCombi[1] == correctJetCombi[1]  &&
             (selJetCombi[2] == correctJetCombi[2] || selJetCombi[2] == correctJetCombi[3]) &&
             (selJetCombi[3] == correctJetCombi[2] || selJetCombi[3] == correctJetCombi[3]) )
            CWUIndex = 0;
          else
            CWUIndex = 1; 
        }
        else
          CWUIndex = 2;

        MSPlot["nBTaggedJets_"+ bTitle[ibTag]+"_BeforeBTag"+leptChannel]->Fill( bTagStudy.getNrBTaggedJets(ibTag), datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight);
        MSPlot["nLightJets_"+   bTitle[ibTag]+"_BeforeBTag"+leptChannel]->Fill( bTagStudy.getNrLightJets(ibTag),   datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight);

        //Apply the event selection
        if( bTagStudy.getNrBTaggedJets(ibTag) < 2 || bTagStudy.getNrLightJets(ibTag) < 2 ) continue;
        if(decayCh == 0) nSelectedMu += 1;
        else if(decayCh == 1) nSelectedEl += 1;

        MSPlot["nPV"+leptChannel]->Fill( nPrimVertices, datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight);
        MSPlot["JetPt_LeadingJet"+leptChannel]->Fill( selJets[0].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight);
        MSPlot["LeptonPt"+leptChannel]->Fill( selLepton.Pt(), datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight);
        MSPlot["JetEta_LeadingJet"+leptChannel]->Fill( selJets[0].Eta(), datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight);
        MSPlot["LeptonEta"+leptChannel]->Fill( selLepton.Eta(), datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight);
        MSPlot["LeptonCharge"+leptChannel]->Fill( leptCharge, datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight);

        //Check whether the Mlb, Mqqb, Mt and MW depends a lot on the considered b-tag and on the fact whether it is applied!
        if(dataSetName.find("TTbarJets_SemiLept") == 0){
          if(correctJetCombi[0] != 9999 && correctJetCombi[1] != 9999 && correctJetCombi[2] != 9999 && correctJetCombi[3] != 9999){
            histo1D["Mlb_CorrTT_"+bTitle[ibTag]]->Fill( (selLepton+selJets[correctJetCombi[0]]).M());
            histo1D["HadrMTop_CorrTT_"+bTitle[ibTag]]->Fill( (selJets[correctJetCombi[1]]+selJets[correctJetCombi[2]]+selJets[correctJetCombi[3]]).M() );
            histo1D["HadrMW_CorrTT_"+bTitle[ibTag]]->Fill( (selJets[correctJetCombi[2]]+selJets[correctJetCombi[3]]).M());
            histo2D["MW_vs_MTop_CorrTT_"+bTitle[ibTag]]->Fill((selJets[correctJetCombi[1]]+selJets[correctJetCombi[2]]+selJets[correctJetCombi[3]]).M(),(selJets[correctJetCombi[2]]+selJets[correctJetCombi[3]]).M());
          }
        }

        //Identical MSPlots with number of jets information after requiring at least two b-jets and at least 2 light jets!
        MSPlot["nSelectedJets_"+bTitle[ibTag]+"_AfterBTag"+leptChannel]->Fill( selJets.size(),                    datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight);
        MSPlot["nBTaggedJets_"+ bTitle[ibTag]+"_AfterBTag"+leptChannel]->Fill( bTagStudy.getNrBTaggedJets(ibTag), datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight);
        MSPlot["nLightJets_"+   bTitle[ibTag]+"_AfterBTag"+leptChannel]->Fill( bTagStudy.getNrLightJets(ibTag),   datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight);      

//        bool CutsSurvived = extraEvtSelCuts.KeepEvent(correctJetCombi, selLepton, selJets, selJetCombi, bTagStudy.getMlbMqqbChiSq(ibTag), CWUIndex, decayCh);
        bool CutsSurvived = true;
 
        //Write out the LHCO output!
        if( getLHCOOutput && bTagChoiceMade){

          //Keep track of the original event number and the madweight numbers:
          EvtNrMatching << "\n    " << iEvt;
          if( iEvt < 10) EvtNrMatching << " "; if( iEvt < 100) EvtNrMatching << " "; if( iEvt < 1000) EvtNrMatching << " "; if( iEvt < 10000) EvtNrMatching << " "; if( iEvt < 100000) EvtNrMatching << " ";
          EvtNrMatching << "              " << CutsSurvived << "                " << fixed << setprecision(9)<< genCosTheta << "";
          if(genCosTheta > 0) EvtNrMatching << " ";

          //if(CutsSurvived) histo1D["CosTheta_SelEvts"]->Fill(kinFunctions.CosTheta(lhcoOutput.getGenLeptTop(), lhcoOutput.getGenLeptW(), lhcoOutput.getGenLepton()));
          // --> Cannot add this since the neutrino is not completely reconstructed and thus the leptonic top is not known ...
          if(CutsSurvived) lhcoOutput.StoreRecoInfo(selLepton, selJets, selJetCombi, decayCh, leptCharge, EvtNrMatching, CWUIndex);
          if(CutsSurvived && genCosTheta != 2.0 && decayCh == 0) histo1D["CosTheta_Gen_SelEvts_Mu"]->Fill(genCosTheta);
          if(CutsSurvived && genCosTheta != 2.0 && decayCh == 1) histo1D["CosTheta_Gen_SelEvts_El"]->Fill(genCosTheta);
          if(CutsSurvived && genCosTheta != 2.0) histo1D["CosTheta_Gen_SelEvts"]->Fill(genCosTheta);
        }
        else if(getLHCOOutput && !bTagChoiceMade){
          cout << " ERROR : Not possible to write the lhco file when all b-tag options are still being considered!!  ==> Setting boolean to FALSE !" << endl;
          getLHCOOutput = false;
        }
      }

    }//End of loop on events
    std::cout<<"    Processed all "<<nrEvts<<" events  --> # selected: "<<nSelectedMu<<" (mu+jets) "<<nSelectedEl<<" (e+jets)"<< flush<<"\n";

    //------------------------------//
    //  Calculate the ChiSq values  //
    //------------------------------//
    if(dataSetName.find("TTbarJets_SemiLept") == 0 && getMassFits){
      std::cout << " -----------------     Fitting the different mass distributions    --------------------- \n" << endl;   
      histo1D["Mlb_CorrectTTEvts"]->Fit("gaus","Q");     
      histo1D["TopMass_Hadr_CorrectTTEvts"]->Fit("gaus","Q");
      histo1D["WMass_Hadr_CorrectTTEvts"]->Fit("gaus","Q");
      histo1D["Mlb_CorrectTTEvts"]->Fit("gaus","Q","", histo1D["Mlb_CorrectTTEvts"]->GetFunction("gaus")->GetParameter(1) - histo1D["Mlb_CorrectTTEvts"]->GetFunction("gaus")->GetParameter(2),
                                                       histo1D["Mlb_CorrectTTEvts"]->GetFunction("gaus")->GetParameter(1) + histo1D["Mlb_CorrectTTEvts"]->GetFunction("gaus")->GetParameter(2));
      histo1D["TopMass_Hadr_CorrectTTEvts"]->Fit("gaus","Q","",
                                                 histo1D["TopMass_Hadr_CorrectTTEvts"]->GetFunction("gaus")->GetParameter(1) - histo1D["TopMass_Hadr_CorrectTTEvts"]->GetFunction("gaus")->GetParameter(2),
                                                 histo1D["TopMass_Hadr_CorrectTTEvts"]->GetFunction("gaus")->GetParameter(1) + histo1D["TopMass_Hadr_CorrectTTEvts"]->GetFunction("gaus")->GetParameter(2));
      histo1D["WMass_Hadr_CorrectTTEvts"]->Fit("gaus","Q","",
                                               histo1D["WMass_Hadr_CorrectTTEvts"]->GetFunction("gaus")->GetParameter(1) - histo1D["WMass_Hadr_CorrectTTEvts"]->GetFunction("gaus")->GetParameter(2),
                                               histo1D["WMass_Hadr_CorrectTTEvts"]->GetFunction("gaus")->GetParameter(1) + histo1D["WMass_Hadr_CorrectTTEvts"]->GetFunction("gaus")->GetParameter(2));
      
      //Write out the mass values!
      std::cout << "   ** Mlb  = " << histo1D["Mlb_CorrectTTEvts"]->GetFunction("gaus")->GetParameter(1)          << " +- " << histo1D["Mlb_CorrectTTEvts"]->GetFunction("gaus")->GetParameter(2) << endl;
      std::cout << "   ** Mqqb = " << histo1D["TopMass_Hadr_CorrectTTEvts"]->GetFunction("gaus")->GetParameter(1) << " +- " << histo1D["TopMass_Hadr_CorrectTTEvts"]->GetFunction("gaus")->GetParameter(2) << endl;
      std::cout << "   ** MW   = " << histo1D["WMass_Hadr_CorrectTTEvts"]->GetFunction("gaus")->GetParameter(1)   << " +- " << histo1D["WMass_Hadr_CorrectTTEvts"]->GetFunction("gaus")->GetParameter(2) << std::endl;
    
      //Output for the different b-tags
      for(int ibTag = 0; ibTag < NrBTags; ibTag++){
        histo1D["Mlb_CorrectTTEvts_"+bTitle[ibTag]]->Fit("gaus","Q"); histo1D["TopMass_Hadr_CorrectTTEvts_"+bTitle[ibTag]]->Fit("gaus","Q"); histo1D["WMass_Hadr_CorrectTTEvts_"+bTitle[ibTag]]->Fit("gaus","Q");
        histo1D["Mlb_CorrectTTEvts_"+bTitle[ibTag]]->Fit("gaus","Q","", 
                                                         histo1D["Mlb_CorrectTTEvts_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(1) - histo1D["Mlb_CorrectTTEvts_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(2),
                                                         histo1D["Mlb_CorrectTTEvts_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(1) + histo1D["Mlb_CorrectTTEvts_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(2));
        histo1D["TopMass_Hadr_CorrectTTEvts_"+bTitle[ibTag]]->Fit("gaus","Q","",
                                          histo1D["TopMass_Hadr_CorrectTTEvts_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(1) - histo1D["TopMass_Hadr_CorrectTTEvts_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(2),
                                          histo1D["TopMass_Hadr_CorrectTTEvts_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(1) + histo1D["TopMass_Hadr_CorrectTTEvts_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(2));
        histo1D["WMass_Hadr_CorrectTTEvts_"+bTitle[ibTag]]->Fit("gaus","Q","",
                                               histo1D["WMass_Hadr_CorrectTTEvts_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(1) - histo1D["WMass_Hadr_CorrectTTEvts_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(2),
                                               histo1D["WMass_Hadr_CorrectTTEvts_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(1) + histo1D["WMass_Hadr_CorrectTTEvts_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(2));
      
        //Write out the mass values!
        cout << "\n   ** Mlb  -- " << bTitle[ibTag] << " = " << histo1D["Mlb_CorrectTTEvts_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(1)          << " +- " << histo1D["Mlb_CorrectTTEvts_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(2) << endl;
        cout << "   ** Mqqb -- " << bTitle[ibTag] << " = " << histo1D["TopMass_Hadr_CorrectTTEvts_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(1) << " +- " << histo1D["TopMass_Hadr_CorrectTTEvts_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(2) << endl;
        cout << "   ** MW   -- " << bTitle[ibTag] << " = " << histo1D["WMass_Hadr_CorrectTTEvts_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(1)   << " +- " << histo1D["WMass_Hadr_CorrectTTEvts_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(2) << std::endl;
      }
    }
 
    //--- Get output from bTagStudy class ---//
    if(dataSetName.find("TTbarJets_SemiLept") == 0) bTagStudy.ReturnBTagTable();
    bTagStudy.CreateHistograms(outputFile, savePDF, pathPNG, iDataSet);  //Security is added inside class such that MSPlots are only written when all datasets are considered!

    //--- Get output from LHCOOutput class ---//
    if(getLHCOOutput && bTagChoiceMade) lhcoOutput.WriteLHCOPlots(outputFile);

    //--- Get output form the ExtraEvtSel class ---//
//    if(dataSetName.find("TTbarJets_SemiLept") == 0) extraEvtSelCuts.StoreCutInfluence(outputFile);

    //---- Close the EventNrMatching output file for the considered dataset
    if(bTagChoiceMade && getLHCOOutput) EvtNrMatching.close();

    inputFile->Close();
    delete inputFile;
  }//End of loop on datasets

  /////////////////////////
  // Write out the plots //
  /////////////////////////
  outputFile -> cd();
  mkdir((pathPNG+"/MSPlots").c_str(),0777);

  if(verbosity > 0) std::cout << "\n - Making the plots in the LightAnomCoupAnalyzer file " << std::endl;
  TDirectory* msdir = outputFile->mkdir("MSPlots");
  msdir->cd(); 
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++){    
    MultiSamplePlot *temp = it->second;
    string name = it->first;
    temp->Draw(name, 1, false, false, false, 1);
    temp->Write(outputFile, name, savePDF, (pathPNG+"/MSPlots/").c_str(), "pdf");
  }

  TDirectory* th1dir = outputFile->mkdir("1D_histograms");
  th1dir->cd();
  for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++){
    TH1F *temp = it->second;
    int N = temp->GetNbinsX();
    temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
    temp->SetBinContent(N+1,0);
    temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
    temp->Write();
  }
  TDirectory* th2dir = outputFile->mkdir("2D_histograms_graphs");
  th2dir->cd();
  for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++){
    TH2F *temp = it->second;
    temp->Write();
  }

  delete outputFile;

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "           hasn't crashed yet ;-)           " << endl;
  cout << "********************************************" << endl;
  
  return 0;
}

/*
	if(jetCombi[0]!=9999 && jetCombi[1]!=9999 && jetCombi[2]!=9999 && jetCombi[3]!=9999){    
	  histo1D["MlbMass"]->Fill((*selJets[jetCombi[0]]+*selLepton).M());
	  histo1D["MqqbMass"]->Fill((*selJets[jetCombi[2]]+*selJets[jetCombi[3]]+*selJets[jetCombi[1]]).M());
        }

        ////////////////////////////////
        //  Mlb and Mqqb information  //
        ////////////////////////////////
        float MlbCorrect = 0, MqqbCorrect = 0;
        if(jetCombi[0] != 9999 && jetCombi[1] != 9999 && jetCombi[2] != 9999 && jetCombi[3] != 9999){
	  MlbCorrect = (*selLepton+*selJets[jetCombi[0]]).M();
	  MqqbCorrect = (*selJets[jetCombi[1]] + *selJets[jetCombi[2]] + *selJets[jetCombi[3]]).M();
	  histo2D["MlbMqqbCorrectAll"]->Fill(MqqbCorrect,MlbCorrect);
        }

      // Count the number of events:
      if(decayCh == 0) nSelectedMu++;
      if(decayCh == 1) nSelectedEl++;
      
      ////////////////////////////////   
      //  Produce Reco LHCO Output  //
      //  --> Last integer = mode   //
      //       * 0 ~ all events     //
      //       * 1 ~ good combi's   //
      //       * 2 ~ bad combi's    //
      ////////////////////////////////
      if( RecoLHCOOutput == true)
        lhcoOutput.StoreRecoInfo(selLepton, selJets, bTagStudy.getBLeptIndex(ChosenBTag), bTagStudy.getBHadrIndex(ChosenBTag), bTagStudy.getLight1Index5Jets(ChosenBTag), bTagStudy.getLight2Index5Jets(ChosenBTag), decayCh, LeptonRecoCharge, jetCombi);

    //---  Mlb and Mqqb fit result ---//
    histo1D["MlbMass"]->Fit("gaus","Q");
    histo1D["MqqbMass"]->Fit("gaus","Q");
    cout <<"\n values for Mlb :"<< histo1D["MlbMass"]->GetFunction("gaus")->GetParameter(1) << " +- " << histo1D["MlbMass"]->GetFunction("gaus")->GetParameter(2) << endl;
    cout <<" values for Mqqb :" << histo1D["MqqbMass"]->GetFunction("gaus")->GetParameter(1) << " +- " << histo1D["MqqbMass"]->GetFunction("gaus")->GetParameter(2) << endl;

*/
