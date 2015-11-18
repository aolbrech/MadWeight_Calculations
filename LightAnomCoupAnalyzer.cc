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
#include "AnomalousCouplings/PersonalClasses/interface/TFLight.h"
#include "AnomalousCouplings/PersonalClasses/interface/BTagStudy_OLD.h"
#include "AnomalousCouplings/PersonalClasses/interface/BTagStudy.h"
#include "AnomalousCouplings/PersonalClasses/interface/LHCOOutput.h"
#include "AnomalousCouplings/PersonalClasses/interface/ExtraEvtSelCuts.h"
#include "TopTreeAnalysisBase/Tools/interface/BTagWeightTools.h"

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

//Check if a file exists or not!
bool file_exist(const char *fileName){
    std::ifstream infile(fileName);
    return infile.good();
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
  bool bTagPlotsMade = true;         //Decide whether the bTag distributions for the efficiencies should be created or read from file
  bool createTFTree = false;         //Decide whether the TFTree should be made or just skipped during this run! 
  bool fillTFAfterCuts = false;      //Decide whether the TFTree is filled before or after the additional event selection cuts (b-tag, chi-sq and Mw-Mt)

  std::string tfFill = "";
  if(fillTFAfterCuts) tfFill = "_AfterExtraCuts";

  std::string bTagPlotsOutput = "PersonalClasses/Calibrations/BTagSF/BTagWeightPlots_CSVT_noTTbar.root";
  int stop = 0;
  if( file_exist(bTagPlotsOutput.c_str()) ){
    if(!bTagPlotsMade){
      std::cout << " \n According to the bTagPlotsMade boolean the plots are not yet made but the output file has been found ... " << std::endl;
      std::cout << "  --> Maybe the output file needs to be updated: " << bTagPlotsOutput << std::endl;
      std::cout << "  ==> Will stop the script until this issue is fixed ! \n" << std::endl;
      stop = 1;
    }
    else
      std::cout << " - Will apply the b-tag efficiencies from file : " << bTagPlotsOutput << std::endl;
  }
  else{
    if(bTagPlotsMade){
      std::cout << "\n According to the bTagPlotsMade boolean the plots should be created but the output file cannot be found ... " << std::endl;
      std::cout << " --> Maybe the wrong output file has been given : " << bTagPlotsOutput << std::endl;
      std::cout << " ==> Will stop the script until this issue is fixed ! \n" << std::endl;
      stop = 1;
    }
    else
      std::cout << " - Will create the b-tag efficiencies and store them in : " << bTagPlotsOutput << std::endl;
  }
  if(stop == 1) return 0;

  //-- Specific b-tag stuff!
  int ChosenBTag = 3;  //-->Corresponds to two T b-tags and no light veto!
  int NrBTags = 6;
  string bTitle[6] = {"LooseTags","MediumTags","MediumTagsLVeto","TightTags","TightTagsMVeto","TightTagsLVeto"};
  string numberBTags = "AllBTags";
  if(bTagChoiceMade){ NrBTags = 1; bTitle[0] = bTitle[ChosenBTag]; numberBTags = bTitle[ChosenBTag];}

  //-------------------------//
  //   Which systematics ?   //
  //-------------------------//
  std::string doLeptonSFShift = "Nominal";   //Other options are "Minus" and "Plus" (Capital letters are needed!)
  std::string doLumiWeightShift = "Nominal"; //Other options are "Minus" and "Plus"
  std::string systematic = "Nominal";

  //b-tag systematics
  int syst_btag = 0; //+1 for SF up, -1 for SF down
  if(systematic == "bTagMinus") syst_btag = -1;
  else if(systematic == "bTagPlus") syst_btag = 1;
  int syst_mistag = 0; //+1 for SF up, -1 for SF down
  if(systematic == "misTagMinus") syst_mistag = -1;
  else if(systematic == "misTagPlus") syst_mistag = 1;
  //syst_btag = 1; syst_mistag = 1;
        
/*  //lepton SF systematics
  string syst_muonSF = "Nominal";
  if(systematic == "MuonSFMinus") syst_muonSF = "Minus";
  else if(systematic == "MuonSFPlus") syst_muonSF = "Plus";
  string syst_electronSF = "Nominal";
  if(systematic == "ElectronSFMinus") syst_electronSF = "Minus";
  else if(systematic == "ElectronSFPlus") syst_electronSF = "Plus";
*/
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
  std::string inputFileDir = "/user/aolbrech/PBS_ScriptRunning/Results/RESULTS_AnomCoup_13112015_100848/";
  inputFiles.push_back((inputFileDir+"AnomCoupLight_Data_Mu_Merged_22Jan2013_Nominal.root").c_str());                                                     //Winter14_v8
  //inputFiles.push_back("/user/aolbrech/PBS_ScriptRunning/Results/RESULTS_AnomCoup_13112015_121239/AnomCoupLight_Data_Mu_Merged_22Jan2013_Nominal.root");  //Winter14_v5
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
    if(dataSet->Name().find("Data_Mu") == 0){        color = kBlack; dataSet->SetTitle("Data"); dataSet->SetName("Data_Mu_Merged"); dataSet->SetEquivalentLuminosity(Luminosity);}

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

    MSPlot["nBJets_BeforeBTag"+leptFlav]   = new MultiSamplePlot(datasets, "nBJets_BeforeBTag"+leptFlav,  10, -0.5, 9.5, "# b-jets");
    MSPlot["partFlav_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets,"partFlav_BeforeBTag"+leptFlav, 5, -0.5, 4.5, "Parton flavour of selected jets");

    MSPlot["nBJets_BTag"+leptFlav]    = new MultiSamplePlot(datasets, "nBJets_BTag"+leptFlav,    10, -0.5, 9.5, "# b-jets");
    MSPlot["nBJets_BTagSF"+leptFlav]  = new MultiSamplePlot(datasets, "nBJets_BTagSF"+leptFlav,  10, -0.5, 9.5, "# b-jets");
    MSPlot["nBJets_AllCuts"+leptFlav] = new MultiSamplePlot(datasets, "nBJets_AllCuts"+leptFlav, 10, -0.5, 9.5, "# b-jets");

    MSPlot["JetsPt_partFlavGluon_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "JetsPt_partFlavGluon_BeforeBTag"+leptFlav, 120, 0, 350, "Jet p_T for all jets with gluon parton flavour");
    MSPlot["JetsPt_partFlavLight_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "JetsPt_partFlavLight_BeforeBTag"+leptFlav, 120, 0, 350, "Jet p_T for all jets with light parton flavour");
    MSPlot["JetsPt_partFlavB_BeforeBTag"+leptFlav]     = new MultiSamplePlot(datasets, "JetsPt_partFlavB_BeforeBTag"+leptFlav, 120, 0, 350, "Jet p_T for all jets with b parton flavour");
    MSPlot["JetsPt_partFlavC_BeforeBTag"+leptFlav]     = new MultiSamplePlot(datasets, "JetsPt_partFlavC_BeforeBTag"+leptFlav, 120, 0, 350, "Jet p_T for all jets with c parton flavour");
    MSPlot["JetsPt_partFlavUndef_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "JetsPt_partFlavUndef_BeforeBTag"+leptFlav, 120, 0, 350, "Jet p_T for all jets with undefined parton flavour");

    MSPlot["JetsPt_partFlavGluon_BTag"+leptFlav]       = new MultiSamplePlot(datasets, "JetsPt_partFlavGluon_BTag"+leptFlav,       120, 0, 350, "Jet p_T for all jets with gluon parton flavour");
    MSPlot["JetsPt_partFlavGluon_BTagSF"+leptFlav]     = new MultiSamplePlot(datasets, "JetsPt_partFlavGluon_BTagSF"+leptFlav,     120, 0, 350, "Jet p_T for all jets with gluon parton flavour");
    MSPlot["JetsPt_partFlavGluon_AllCuts"+leptFlav]    = new MultiSamplePlot(datasets, "JetsPt_partFlavGluon_AllCuts"+leptFlav,    120, 0, 350, "Jet p_T for all jets with gluon parton flavour");
    MSPlot["BJetPt_partFlavGluon_BTag"+leptFlav]    = new MultiSamplePlot(datasets, "BJetPt_partFlavGluon_BTag"+leptFlav,    120, 0, 350, "Jet p_T for b-tagged jets with gluon parton flavour");
    MSPlot["BJetPt_partFlavGluon_BTagSF"+leptFlav]  = new MultiSamplePlot(datasets, "BJetPt_partFlavGluon_BTagSF"+leptFlav,  120, 0, 350, "Jet p_T for b-tagged jets with gluon parton flavour");
    MSPlot["BJetPt_partFlavGluon_AllCuts"+leptFlav] = new MultiSamplePlot(datasets, "BJetPt_partFlavGluon_AllCuts"+leptFlav, 120, 0, 350, "Jet p_T for b-tagged jets with gluon parton flavour");

    MSPlot["LeptonPt_BeforePU"+leptFlav] = new MultiSamplePlot(datasets, "LeptonPt_BeforePU"+leptFlav, 100, 0, 200,"Lepton p_{T} (GeV)");
    MSPlot["LeptonPt_PU_NoLeptonSF"+leptFlav] = new MultiSamplePlot(datasets, "LeptonPt_PU_NoLeptonSF"+leptFlav, 100, 0, 200,"Lepton p_{T} (GeV)");
    MSPlot["LeptonPt_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "LeptonPt_BeforeBTag"+leptFlav, 100, 0, 200,"Lepton p_{T} (GeV)");
    MSPlot["LeptonPt"+leptFlav] = new MultiSamplePlot(datasets, "LeptonPt"+leptFlav, 100, 0, 200,"Lepton p_{T} (GeV)");
    MSPlot["LeptonPt_BTagSF"+leptFlav] = new MultiSamplePlot(datasets, "LeptonPt_BTagSF"+leptFlav, 100, 0, 200,"Lepton p_{T} (GeV)");
    MSPlot["LeptonPt_AllCuts"+leptFlav] = new MultiSamplePlot(datasets, "LeptonPt_AllCuts"+leptFlav, 100, 0, 200,"Lepton p_{T} (GeV)");
    MSPlot["LeptonMass_AllCuts"+leptFlav] = new MultiSamplePlot(datasets, "LeptonMass_AllCuts"+leptFlav, 100, 0, 2,"Lepton p_{T} (GeV)");
    MSPlot["LeptonEta"+leptFlav] = new MultiSamplePlot(datasets, "LeptonEta"+leptFlav, 50, -2.6, 2.6,"Lepton #eta");
    MSPlot["LeptonEta_BTagSF"+leptFlav] = new MultiSamplePlot(datasets, "LeptonEta_BTagSF"+leptFlav, 50, -2.6, 2.6,"Lepton #eta");
    MSPlot["LeptonEta_AllCuts"+leptFlav] = new MultiSamplePlot(datasets, "LeptonEta_AllCuts"+leptFlav, 50, -2.6, 2.6,"Lepton #eta");
    MSPlot["LeptonCharge"+leptFlav] = new MultiSamplePlot(datasets, "LeptonCharge"+leptFlav, 10, -2, 2,"Lepton charge");
    MSPlot["LeptonCharge_BTagSF"+leptFlav] = new MultiSamplePlot(datasets, "LeptonCharge_BTagSF"+leptFlav, 10, -2, 2,"Lepton charge");
    MSPlot["LeptonCharge_AllCuts"+leptFlav] = new MultiSamplePlot(datasets, "LeptonCharge_AllCuts"+leptFlav, 10, -2, 2,"Lepton charge");

    MSPlot["JetPt_LeadingJet_BeforePU"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_LeadingJet_BeforePU"+leptFlav, 120, 0, 350,"Jet p_{T} (GeV)");
    MSPlot["JetPt_LeadingJet_PU_NoLeptonSF"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_LeadingJet_PU_NoLeptonSF"+leptFlav, 120, 0, 350,"Jet p_{T} (GeV)");
    MSPlot["JetPt_LeadingJet_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_LeadingJet_BeforeBTag"+leptFlav, 120, 0, 350,"Jet p_{T} (GeV)");
    MSPlot["JetPt_LeadingJet"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_LeadingJet"+leptFlav, 120, 0, 350,"Jet p_{T} (GeV)");
    MSPlot["JetPt_LeadingJet_AllCuts"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_LeadingJet_AllCuts"+leptFlav, 120, 0, 350,"Jet p_{T} (GeV)");
    MSPlot["JetMass_LeadingJet_AllCuts"+leptFlav] = new MultiSamplePlot(datasets, "JetMass_LeadingJet_AllCuts"+leptFlav, 75, 0, 50,"Jet p_{T} (GeV)");
    MSPlot["JetPt_LeadingJet_BTagSF"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_LeadingJet_BTagSF"+leptFlav, 120, 0, 350,"Jet p_{T} (GeV)");
    MSPlot["JetEta_LeadingJet"+leptFlav] = new MultiSamplePlot(datasets, "JetEta_LeadingJet"+leptFlav, 50, -2.6, 2.6,"Jet #eta");
    MSPlot["JetEta_LeadingJet_BTagSF"+leptFlav] = new MultiSamplePlot(datasets, "JetEta_LeadingJet_BTagSF"+leptFlav, 50, -2.6, 2.6,"Jet #eta");
    MSPlot["JetEta_LeadingJet_AllCuts"+leptFlav] = new MultiSamplePlot(datasets, "JetEta_LeadingJet_AllCuts"+leptFlav, 50, -2.6, 2.6,"Jet #eta");

    MSPlot["nPV_BeforePU"+leptFlav] = new MultiSamplePlot(datasets, "nPV_BeforePU"+leptFlav, 50, 0, 50,"Number of primary vertices");
    MSPlot["nPV_PU_NoLeptonSF"+leptFlav] = new MultiSamplePlot(datasets, "nPV_PU_NoLeptonSF"+leptFlav, 50, 0, 50,"Number of primary vertices");
    MSPlot["nPV_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "nPV_BeforeBTag"+leptFlav, 50, 0, 50,"Number of primary vertices");
    MSPlot["nPV"+leptFlav] = new MultiSamplePlot(datasets, "nPV"+leptFlav, 50, 0, 50,"Number of primary vertices");
    MSPlot["nPV_BTagSF"+leptFlav] = new MultiSamplePlot(datasets, "nPV_BTagSF"+leptFlav, 50, 0, 50,"Number of primary vertices");
    MSPlot["nPV_AllCuts"+leptFlav] = new MultiSamplePlot(datasets, "nPV_AllCuts"+leptFlav, 50, 0, 50,"Number of primary vertices");

    MSPlot["CSVDistr_BJets"+leptFlav] = new MultiSamplePlot(datasets, "CSVDistr_BJets"+leptFlav,50,-1,2,"CSV discriminant");
    MSPlot["CSVDistr_LightJets"+leptFlav] = new MultiSamplePlot(datasets, "CSVDistr_LightJets"+leptFlav,50,-1,2,"CSV discriminant");
  }

  //-----------------------//
  // Initialize histograms //
  //-----------------------//
  map<string,TH1F*> histo1D;
  map<string,TH2F*> histo2D;  


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

  //------------------------//
  // B-tagging scale factor //
  //------------------------//
  //Method 1a of https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods will be used!

  BTagWeightTools *bTagTool = new BTagWeightTools("PersonalClasses/Calibrations/BTagSF/SFb-pt_NOttbar_payload_EPS13.txt","CSVT");   //Standard EPS2013 is used!
  //BTagWeightTools *bTagTool = new BTagWeightsTools("PersonalClasses/Calibrations/BTagSF/SFb-pt_WITHttbar_payload_EPS13.txt","CSVT");   //Standard EPS2013 is used!

  //During a first run the plots need to be created!
  if(!bTagPlotsMade){
    bTagTool->InitializeMCEfficiencyHistos(50,30.,340.,2);  //How to get these histo-values?
  }
  else{
    bTagTool->ReadMCEfficiencyHistos("PersonalClasses/Calibrations/BTagSF/BTagWeightPlots_CSVT_noTTbar.root");
  }

  //-----------------------//
  // Load personal classes //
  //-----------------------//
  float Mlb  = 108.1841, S_Mlb  = 31.4213;
  float Mqqb = 174.6736, S_Mqqb = 17.5757;
  float MW = 83.8037, S_MW = 10.2385;
  ExtraEvtSelCuts extraEvtSelCuts(Mqqb, S_Mqqb, MW, S_MW, bTagChoiceMade, 3, 2);
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
    string dsName = dataSet->Name();
    string dsTitle = dataSet->Title();
    if(verbosity > 0) std::cout << "\n   *** Looking at dataset "<< dsName << " (" << iDataSet+1 << "/" << inputFiles.size() << ") with " << nrEvts << " events (total = " << inLightTree->GetEntries() << ") !" << std::endl;

    //**** Create the 1D and 2D histograms for each dataset specific! ****//
    histo1D["BTagEfficiency_BJets"+dsName] = new TH1F(("BTagEfficiency_BJets"+dsName).c_str(),("b-tag efficiency for b-jets for "+dsTitle).c_str(),200, 0,2);
    histo1D["BTagEfficiency_LightJets"+dsName] = new TH1F(("BTagEfficiency_LightJets"+dsName).c_str(),("b-tag efficiency for light jets for "+dsTitle).c_str(),200, 0,2);

    histo1D["BTagWeight_"+dsName] = new TH1F(("BTagWeight_"+dsName).c_str(),("b-tag reweighting distribution for "+dsTitle).c_str(),200,0,2);
    histo1D["BTagWeight_AltMethod_"+dsName] = new TH1F(("BTagWeight_AltMethod_"+dsName).c_str(),("b-tag reweighting distribution for "+dsTitle).c_str(),200,0,2);
    histo1D["lumiWeight_"+dsName] = new TH1F(("lumiWeight_"+dsName).c_str(),("lumi reweighting distribution for "+dsTitle).c_str(),200,0,2);
    histo1D["leptonSF_"+dsName] = new TH1F(("leptonSF_"+dsName).c_str(), ("lepton scale factor for "+dsTitle).c_str(),200,0,2);

    histo2D["bTagWeight_vs_BJetPt_"+dsName] = new TH2F(("bTagWeight_vs_BJetPt_"+dsName).c_str(),("btag weight versus b-jet PT for "+dsTitle).c_str(),800,0,800,150,0,1);
    histo2D["bTagWeight_vs_BJetPt_partFlavB_"+dsName] = new TH2F(("bTagWeight_vs_BJetPt_partFlavB_"+dsName).c_str(),("btag weight versus b-jet PT (partFlav = b-jet) for "+dsTitle).c_str(),800,0,800,150,0,1);
    histo2D["bTagWeight_vs_BJetPt_partFlavC_"+dsName] = new TH2F(("bTagWeight_vs_BJetPt_partFlavC_"+dsName).c_str(),("btag weight versus b-jet PT (partFlav = c-jet) for "+dsTitle).c_str(),800,0,800,150,0,1);
    histo2D["bTagWeight_vs_BJetPt_partFlavLight_"+dsName] = new TH2F(("bTagWeight_vs_BJetPt_partFlavLight_"+dsName).c_str(),("btag weight versus b-jet PT (partFlav = light jet) for "+dsTitle).c_str(),800,0,800,150,0,1);
    histo2D["bTagWeight_vs_LightJetPt_"+dsName] = new TH2F(("bTagWeight_vs_LightJetPt_"+dsName).c_str(),("btag weight versus light jet PT for "+dsTitle).c_str(),800,0,800,150,0,1);
    histo2D["bTagWeight_vs_LightJetPt_partFlavB_"+dsName] = new TH2F(("bTagWeight_vs_LightJetPt_partFlavB_"+dsName).c_str(),("btag weight versus light jet PT (partFlav = b-jet) for "+dsTitle).c_str(),800,0,800,150,0,1);
    histo2D["bTagWeight_vs_LightJetPt_partFlavC_"+dsName] = new TH2F(("bTagWeight_vs_LightJetPt_partFlavC_"+dsName).c_str(),("btag weight versus light jet PT (partFlav = c-jet) for "+dsTitle).c_str(),800,0,800,150,0,1);
    histo2D["bTagWeight_vs_LightJetPt_partFlavLight_"+dsName] = new TH2F(("bTagWeight_vs_LightJetPt_partFlavLight_"+dsName).c_str(),("btag weight versus light jet PT (partFlav = light jet) for "+dsTitle).c_str(),800,0,800,150,0,1);
    histo2D["bTagWeight_vs_JetPt_partFlavLight_"+dsName] = new TH2F(("bTagWeight_vs_JetPt_partFlavLight_"+dsName).c_str(),("btag weight versus jet PT (partFlav = light jet) for "+dsTitle).c_str(),800,0,800,150,0,1);
    histo1D["partonFlav_BJets_"+dsName] = new TH1F(("partonFlav_BJets_"+dsName).c_str(),("parton flavour for b-jets ("+dsTitle+")").c_str(),101,-0.5,100.5);
    histo1D["partonFlav_LightJets_"+dsName] = new TH1F(("partonFlav_LightJets_"+dsName).c_str(),("parton flavour for light jets ("+dsTitle+")").c_str(),101,-0.5,100.5);
    

    //Only want these to be filled for TTbarJets_SemiMu!  (what about FullLept and FullHadr ??)
    if(dsName.find("TTbarJets_SemiLept") == 0){ 
      histo1D["Mlb_CorrTT_"+dsName]        = new TH1F(("Mlb_CorrTT_"+dsName).c_str(),       ("mass_{l,b} for the actual "+dsTitle+" events").c_str(),     200,  0, 250);
      histo1D["HadrMTop_CorrTT_"+dsName]   = new TH1F(("HadrMTop_CorrTT_"+dsName).c_str(),  ("Hadronic m_{t} for the actual "+dsTitle+" events").c_str(), 200, 50, 350);
      histo1D["HadrMW_CorrTT_"+dsName]     = new TH1F(("HadrMW_CorrTT_"+dsName).c_str(),    ("Hadronic m_{W} for the actual "+dsTitle+" events").c_str(), 200,  0, 250);
      histo2D["MW_vs_MTop_CorrTT_"+dsName] = new TH2F(("MW_vs_MTop_CorrTT_"+dsName).c_str(),("Hadronic m_{W} versus m_{t} for the actual "+dsTitle+" jet combinations").c_str(),200,50,350,200,0,250);

      for(int ibTag = 0; ibTag < NrBTags; ibTag++){
        histo1D["Mlb_CorrTT_"+dsName+"_"+bTitle[ibTag]]        = new TH1F(("Mlb_CorrTT_"+dsName+"_"+bTitle[ibTag]).c_str(),  ("mass_{l,b} for the actual "+dsTitle+" events -- "+bTitle[ibTag]).c_str(), 200, 0, 250);
        histo1D["HadrMTop_CorrTT_"+dsName+"_"+bTitle[ibTag]]   = new TH1F(("HadrMTop_CorrTT_"+dsName+"_"+bTitle[ibTag]).c_str(), ("Hadronic m_{t} for the actual "+dsTitle+" events -- "+bTitle[ibTag]).c_str(), 200, 50, 350);
        histo1D["HadrMW_CorrTT_"+dsName+"_"+bTitle[ibTag]]     = new TH1F(("HadrMW_CorrTT_"+dsName+"_"+bTitle[ibTag]).c_str(),   ("Hadronic m_{W} for the actual "+dsTitle+" events -- "+bTitle[ibTag]).c_str(), 200, 0, 250);
        histo2D["MW_vs_MTop_CorrTT_"+dsName+"_"+bTitle[ibTag]] = new TH2F(("MW_vs_MTop_CorrTT_"+dsName+"_"+bTitle[ibTag]).c_str(), ("Hadronic m_{W} versus m_{t} for actual "+dsTitle+" events -- "+bTitle[ibTag]).c_str(), 200, 50, 350, 200, 0,250);
      }
    }
    
    for(int ii = 0; ii < 2; ii++){
      string leptFlav = leptFlavs[ii];
      histo1D["CosTheta_Gen_SelEvts_"+dsName+""+leptFlav] = new TH1F(("CosTheta_Gen_SelEvts_"+dsName+""+leptFlav).c_str(),("Cos #theta^{*}_{gen} distribution for the selected "+dsTitle+" events").c_str(),200,-1,1);
    }
    histo1D["CosTheta_Gen_SelEvts_"+dsName] = new TH1F(("CosTheta_Gen_SelEvts_"+dsName).c_str(),("Cos #theta^{*}_{gen} distribution for the selected "+dsTitle+" events").c_str(),200,-1,1);

    TFLight* tfLight_mu = 0;
    //TFLight* tfLight_el = 0;
    TTree* TFLightTree;
    TFile* TFLightFile = 0;
    if(dsName.find("TTbarJets_SemiLept") == 0 && createTFTree){
      std::cout << " --> Going into this loop for dataset : " << dsName  << std::endl;
      //Initialize LightTuple (TFTree) specific stuff:
      TFLightFile = new TFile(("TFTree/TFLight_"+dsName+"_"+systematic+""+tfFill+".root").c_str(),"RECREATE");  
      TFLightTree = new TTree("TFLightTree",("Tree containing the TFLight information for "+dsName+" at "+timestamp()).c_str());
      TFLightTree->Branch("TheTFLight_muCh","TFLight",&tfLight_mu);
      //TFLightTree->Branch("TheTFLight_elCh","TFLight",&tfLight_el);
    }

    //-----------------------//
    // Load personal classes //
    //-----------------------//
    bTagStudy.InitializeDataSet(dsName);
    LHCOOutput lhcoOutput(verbosity, getLHCOOutput, splitLeptonChargeLHCO, getCorrectAndWrongLHCO);
    extraEvtSelCuts.Initialize(bTitle[0], dsName);
    lhcoOutput.Initialize("Reco", dsName);

    ofstream EvtNrMatching;
    if(bTagChoiceMade && getLHCOOutput){
      EvtNrMatching.open(("MadWeightInput/AnalyzerOutput/EventNrMatching_"+dsName+".txt").c_str());
      EvtNrMatching << "  Event Nr     Extra cuts survived      Gen cos theta*        Main lhco file      MW Number (main)      TTbar splitting lhco file      MW Number (splitting)    MC scale factor" << endl;
    }

    // --------------------------- //
    //  Start looping over events  //
    // --------------------------- //    
    int nSelectedMu = 0, nSelectedEl = 0;
    for(unsigned int iEvt = 0; iEvt < nrEvts; iEvt++){
      inLightTree->GetEvent(iEvt);
      if(iEvt%5000 == 0)
	std::cout<<"    Processing the "<<iEvt<<"th event ("<<((double)iEvt/(double)nrEvts)*100<<"%)"<<" -> # selected: "<<nSelectedMu<<" (mu+jets) "<<nSelectedEl<<" (e+jets)"<< flush<<"\r";

      //*** Start with corrections ***//
      // Beam scraping and PU reweighting
      double lumiWeight = 1;  
      if(! (dsName.find("Data") == 0 || dsName.find("data") == 0 || dsName.find("DATA") == 0) ){   
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
      vector<int> partFlavour = light->selectedJetsPartonFlavour();
      TLorentzVector selLepton = light->selectedLepton();
      TLorentzVector MET = light->met();
      int decayCh = light->decayChannel();  //0 = semiMu and 1 = semiEl
      std::string leptChannel;
      if(decayCh == 0) leptChannel = "_mu";
      else if(decayCh == 1) leptChannel = "_el";
      float leptCharge = light->leptonCharge();
      vector<int> correctJetCombi = light->correctJetCombi();    //0 = LeptB, 1 = HadrB, 2 = Quark1 & 3 = Quark2
      float genCosTheta = light->genCosTh();

      //Distributions before PU reweighting is taken into account!
      MSPlot["nPV_BeforePU"+leptChannel]->Fill( nPrimVertices, datasets[iDataSet], true, Luminosity*scaleFactor);
      MSPlot["JetPt_LeadingJet_BeforePU"+leptChannel]->Fill( selJets[0].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
      MSPlot["LeptonPt_BeforePU"+leptChannel]->Fill( selLepton.Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
      scaleFactor = scaleFactor * lumiWeight;
      histo1D["lumiWeight_"+dsName]->Fill(lumiWeight);

      MSPlot["nPV_PU_NoLeptonSF"+leptChannel]->Fill( nPrimVertices, datasets[iDataSet], true, Luminosity*scaleFactor);
      MSPlot["JetPt_LeadingJet_PU_NoLeptonSF"+leptChannel]->Fill( selJets[0].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
      MSPlot["LeptonPt_PU_NoLeptonSF"+leptChannel]->Fill( selLepton.Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
      //Get the correct lepton scalefactor!
      double leptonSF = 1;
      if(decayCh == 0 && !(dsName.find("Data") == 0 || dsName.find("DATA") == 0 || dsName.find("data") == 0) )
        leptonSF = leptonTools->getMuonSF(selLepton.Eta(), selLepton.Pt(), systematic);
      else if(decayCh == 1 && !(dsName.find("Data") == 0 || dsName.find("DATA") == 0 || dsName.find("data") == 0) )
        leptonSF = leptonTools->getElectronSF(selLepton.Eta(), selLepton.Pt(), systematic);
      scaleFactor = scaleFactor * leptonSF;
      histo1D["leptonSF_"+dsName]->Fill(leptonSF); 

      //Fill the b-tag histo's in case they do not yet exist
      if(!bTagPlotsMade)
        bTagTool->FillMCEfficiencyHistos(selJets, partFlavour, bTagCSV);

      //ooooooooOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOooooooooooooOOOOO
      //ooOOooOOoo      Reading out nTuples done           ooOOooOOoo
      //ooOOooOOoo-----------------------------------------ooOOooOOoo
      //ooOOooOOoo      Start of actual analysis           ooOOooOOoo
      //ooooooooOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOooooooooooooOOOOO
      bTagStudy.CalculateJets(selJets, bTagCSV, correctJetCombi, selLepton, datasets[iDataSet], Luminosity*scaleFactor);
 
      //MSPlots with number of jets information before requiring at least two b-jets and at least 2 light jets!
      MSPlot["nSelectedJets_BeforeBTag"+leptChannel]->Fill( selJets.size(),                    datasets[iDataSet], true, Luminosity*scaleFactor);

      MSPlot["nPV_BeforeBTag"+leptChannel]->Fill( nPrimVertices, datasets[iDataSet], true, Luminosity*scaleFactor);
      MSPlot["JetPt_LeadingJet_BeforeBTag"+leptChannel]->Fill( selJets[0].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
      MSPlot["LeptonPt_BeforeBTag"+leptChannel]->Fill( selLepton.Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);

      //Store the expected Mlb, Mqqb, Mt and MW distributions
      if(dsName.find("TTbarJets_SemiLept") == 0){
        if(correctJetCombi[0] != 9999 && correctJetCombi[1] != 9999 && correctJetCombi[2] != 9999 && correctJetCombi[3] != 9999){
          histo1D["Mlb_CorrTT_"+dsName]->Fill( (selLepton+selJets[correctJetCombi[0]]).M());
          histo1D["HadrMTop_CorrTT_"+dsName]->Fill( (selJets[correctJetCombi[1]]+selJets[correctJetCombi[2]]+selJets[correctJetCombi[3]]).M() );
          histo1D["HadrMW_CorrTT_"+dsName]->Fill( (selJets[correctJetCombi[2]]+selJets[correctJetCombi[3]]).M());
          histo2D["MW_vs_MTop_CorrTT_"+dsName]->Fill( (selJets[correctJetCombi[1]]+selJets[correctJetCombi[2]]+selJets[correctJetCombi[3]]).M(), (selJets[correctJetCombi[2]]+selJets[correctJetCombi[3]]).M() );
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

        MSPlot["nBTaggedJets_"+ bTitle[ibTag]+"_BeforeBTag"+leptChannel]->Fill( bTagStudy.getNrBTaggedJets(ibTag), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["nLightJets_"+   bTitle[ibTag]+"_BeforeBTag"+leptChannel]->Fill( bTagStudy.getNrLightJets(ibTag),   datasets[iDataSet], true, Luminosity*scaleFactor);

        MSPlot["nBJets_BeforeBTag"+leptChannel]->Fill(bTagStudy.getNrBTaggedJets(ibTag), datasets[iDataSet], true, Luminosity*scaleFactor);
        for(int ii = 0; ii < selJets.size(); ii++){
          if(dsName.find("Data") != 0 && (partFlavour[ii] == 1 || partFlavour[ii] == 2 || partFlavour[ii] == 3) ){
            MSPlot["partFlav_BeforeBTag"+leptChannel]->Fill(1, datasets[iDataSet], true, Luminosity*scaleFactor);
            MSPlot["JetsPt_partFlavLight_BeforeBTag"+leptChannel]->Fill(selJets[ii].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
          }
          else if(dsName.find("Data") != 0 && partFlavour[ii] == 4 ){
            MSPlot["partFlav_BeforeBTag"+leptChannel]->Fill(2, datasets[iDataSet], true, Luminosity*scaleFactor);
            MSPlot["JetsPt_partFlavC_BeforeBTag"+leptChannel]->Fill(selJets[ii].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
          }
          else if(dsName.find("Data") != 0 && partFlavour[ii] == 5 ){
            MSPlot["partFlav_BeforeBTag"+leptChannel]->Fill(3, datasets[iDataSet], true, Luminosity*scaleFactor);
            MSPlot["JetsPt_partFlavB_BeforeBTag"+leptChannel]->Fill(selJets[ii].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
          }
          else if(dsName.find("Data") != 0 && partFlavour[ii] == 21 ){
            MSPlot["partFlav_BeforeBTag"+leptChannel]->Fill(4, datasets[iDataSet], true, Luminosity*scaleFactor);
            MSPlot["JetsPt_partFlavGluon_BeforeBTag"+leptChannel]->Fill(selJets[ii].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
          }
          else{
            MSPlot["partFlav_BeforeBTag"+leptChannel]->Fill(0, datasets[iDataSet], true, Luminosity*scaleFactor);
            MSPlot["JetsPt_partFlavUndef_BeforeBTag"+leptChannel]->Fill(selJets[ii].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
          }
        }
      
        if(dsName.find("TTbarJets_SemiLept") == 0 && !fillTFAfterCuts && createTFTree){  //Combine both muon and electron events!
          tfLight_mu = new TFLight();
            
          double fullScaleFactor = scaleFactor;   
          tfLight_mu->setFullScaleFactor(fullScaleFactor);
          tfLight_mu->setSelectedJets(selJets);
          tfLight_mu->setSelectedLepton(selLepton);
          tfLight_mu->setCorrectJetCombi(correctJetCombi);
	
          //Store the information needed for the TF (but only has value when dataset is ttbar)
          tfLight_mu->setGenVectorLight1( light->genVectorLight1() );
          tfLight_mu->setGenVectorLight2( light->genVectorLight2() );
          tfLight_mu->setGenVectorHadrB(  light->genVectorHadrB()  );
          tfLight_mu->setGenVectorLeptB(  light->genVectorLeptB()  );
          tfLight_mu->setGenVectorLepton( light->genVectorLepton() );

          TFLightTree->Fill();
          delete tfLight_mu;
          //----  End of Tree file filling (for TF's after evtSel)  ----//
        }

        //*****************//  
        // Apply the b-tag //
        //*****************//  
        if( bTagStudy.getNrBTaggedJets(ibTag) < 2 || bTagStudy.getNrLightJets(ibTag) < 2 ) continue;

        if(partFlavour[bTagStudy.getBLeptIndex(ibTag)] == 21) MSPlot["BJetPt_partFlavGluon_BTag"+leptChannel]->Fill( selJets[bTagStudy.getBLeptIndex(ibTag)].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        if(partFlavour[bTagStudy.getBHadrIndex(ibTag)] == 21) MSPlot["BJetPt_partFlavGluon_BTag"+leptChannel]->Fill( selJets[bTagStudy.getBHadrIndex(ibTag)].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        for(int ii = 0; ii < selJets.size(); ii++){
          if(partFlavour[ii] == 21) MSPlot["JetsPt_partFlavGluon_BTag"+leptChannel]->Fill(selJets[ii].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        }

        MSPlot["nBJets_BTag"+leptChannel]->Fill(bTagStudy.getNrBTaggedJets(ibTag), datasets[iDataSet], true, Luminosity*scaleFactor);
        vector<TLorentzVector> selJetsAfterBTag;
        vector<int> partFlavAfterBTag;
        vector<float> bTagCSVAfterBTag;
        selJetsAfterBTag.push_back(selJets[bTagStudy.getBLeptIndex(ibTag)]); partFlavAfterBTag.push_back(partFlavour[bTagStudy.getBLeptIndex(ibTag)]); bTagCSVAfterBTag.push_back(bTagCSV[bTagStudy.getBLeptIndex(ibTag)]);
        selJetsAfterBTag.push_back(selJets[bTagStudy.getBHadrIndex(ibTag)]); partFlavAfterBTag.push_back(partFlavour[bTagStudy.getBHadrIndex(ibTag)]); bTagCSVAfterBTag.push_back(bTagCSV[bTagStudy.getBHadrIndex(ibTag)]);
        selJetsAfterBTag.push_back(selJets[bTagStudy.getLight1Index(ibTag)]); partFlavAfterBTag.push_back(partFlavour[bTagStudy.getLight1Index(ibTag)]); bTagCSVAfterBTag.push_back(bTagCSV[bTagStudy.getLight1Index(ibTag)]);
        selJetsAfterBTag.push_back(selJets[bTagStudy.getLight2Index(ibTag)]); partFlavAfterBTag.push_back(partFlavour[bTagStudy.getLight2Index(ibTag)]); bTagCSVAfterBTag.push_back(bTagCSV[bTagStudy.getLight2Index(ibTag)]);

        if(decayCh == 0) nSelectedMu += 1;
        else if(decayCh == 1) nSelectedEl += 1;

        MSPlot["nPV"+leptChannel]->Fill( nPrimVertices, datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["JetPt_LeadingJet"+leptChannel]->Fill( selJets[0].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["LeptonPt"+leptChannel]->Fill( selLepton.Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["JetEta_LeadingJet"+leptChannel]->Fill( selJets[0].Eta(), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["LeptonEta"+leptChannel]->Fill( selLepton.Eta(), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["LeptonCharge"+leptChannel]->Fill( leptCharge, datasets[iDataSet], true, Luminosity*scaleFactor);

        //Get the bTag scaleFactor!
        double BTagWeight = 1;
        double BTagWeight2b = 1;
        if(bTagPlotsMade && !(dsName.find("Data") == 0 || dsName.find("DATA") == 0)){
          BTagWeight = bTagTool->getMCEventWeight(selJetsAfterBTag, partFlavAfterBTag, bTagCSVAfterBTag, syst_btag, syst_mistag);

          //Take into account the two b-tags (following 1b this time ...)
          vector<float> bTagSFs, bTagEffs;
          bTagSFs.clear(); bTagEffs.clear();
          double w0TagsMC = 1;
          double w0TagsData = 1;
          //std::cout << " w0Tags is = " << endl;
          for(int iJet = 0; iJet < selJetsAfterBTag.size(); iJet++){
            bTagSFs.push_back(bTagTool->getSF(selJetsAfterBTag[iJet].Pt(), selJetsAfterBTag[iJet].Eta(), partFlavAfterBTag[iJet], "CSVT", syst_btag, syst_mistag));
            bTagEffs.push_back(bTagTool->getTagEff(selJetsAfterBTag[iJet].Pt(), selJetsAfterBTag[iJet].Eta(), partFlavAfterBTag[iJet]));
            w0TagsMC = w0TagsMC * (1-bTagEffs[iJet]);
            w0TagsData = w0TagsData * (1-bTagEffs[iJet]*bTagSFs[iJet]);
            //std::cout << "  " << iJet << ") Jet with flavour " << partFlavAfterBTag[iJet] << ": 1 - " << bTagEffs[iJet] << " * " << bTagSFs[iJet] << " = " <<  1-bTagEffs[iJet]*bTagSFs[iJet] << " ==> w0 = " << w0TagsData << std::endl;
          }

          //Now calculate the event weight w(>= 2 btags):
          double w1TagIndivMC[4] = {0, 0, 0, 0};
          double w1TagIndivData[4] = {0, 0, 0, 0};
          double w1TagMC = 0;
          double w1TagData = 0;
          //std::cout << " w1Tag = " << endl;
          for(int i = 0; i< selJetsAfterBTag.size(); i++){
            w1TagIndivMC[i] = bTagEffs[i];
            w1TagIndivData[i] = bTagSFs[i]*bTagEffs[i];
            //std::cout << "  " << i << ") Jet with flavour " << partFlavAfterBTag[i] << ": " << bTagSFs[i] << " * " << bTagEffs[i] << " = " << w1TagIndivData[i] << std::endl;
            for(int j = 0; j < selJetsAfterBTag.size(); j++){
              if(j != i) w1TagIndivMC[i] = w1TagIndivMC[i] * (1-bTagEffs[j]);
              if(j != i) w1TagIndivData[i] = w1TagIndivData[i] * (1-bTagSFs[j]*bTagEffs[j]);
              //if(j != i) std::cout << "   " << j << ") Jet with flavour " << partFlavAfterBTag[j] << ": 1- " << bTagSFs[j] << " * " << bTagEffs[j] << " = " << 1-bTagSFs[j]*bTagEffs[j] << " ==> w1 = " << w1TagIndivData[i] << std::endl;
            }
            w1TagMC = w1TagMC + w1TagIndivMC[i];
            w1TagData = w1TagData + w1TagIndivData[i];
            //std::cout << "  ==> w1 = " << w1TagData << std::endl;
          }
          BTagWeight2b = 1.0-w0TagsData-w1TagData;
          //BTagWeight2b = (1-w0TagsData-w1TagData)/(1-w0TagsMC-w1TagData);
          float DataBTag = 1-(1-bTagEffs[0]*bTagSFs[0])*(1-bTagEffs[1]*bTagSFs[1])*(1-bTagEffs[2]*bTagSFs[2])*(1-bTagEffs[3]*bTagSFs[3]) - ( bTagEffs[0]*bTagSFs[0]*((1-bTagEffs[1]*bTagSFs[1])*(1-bTagEffs[2]*bTagSFs[2])*(1-bTagEffs[3]*bTagSFs[3])) + bTagEffs[1]*bTagSFs[1]*((1-bTagEffs[0]*bTagSFs[0])*(1-bTagEffs[2]*bTagSFs[2])*(1-bTagEffs[3]*bTagSFs[3])) + bTagEffs[2]*bTagSFs[2]*((1-bTagEffs[0]*bTagSFs[0])*(1-bTagEffs[1]*bTagSFs[1])*(1-bTagEffs[2]*bTagSFs[2])) + bTagEffs[3]*bTagSFs[3]*((1-bTagEffs[0]*bTagSFs[0])*(1-bTagEffs[1]*bTagSFs[1])*(1-bTagEffs[2]*bTagSFs[2])));
          float MCBTag = 1-(1-bTagEffs[0])*(1-bTagEffs[1])*(1-bTagEffs[2])*(1-bTagEffs[3]) - ( bTagEffs[0]*((1-bTagEffs[1])*(1-bTagEffs[2])*(1-bTagEffs[3])) + bTagEffs[1]*((1-bTagEffs[0])*(1-bTagEffs[2])*(1-bTagEffs[3])) + bTagEffs[2]*((1-bTagEffs[0])*(1-bTagEffs[1])*(1-bTagEffs[2])) + bTagEffs[3]*((1-bTagEffs[0])*(1-bTagEffs[1])*(1-bTagEffs[2])));
          float DataBTag1 = 1-(1-bTagEffs[0]*bTagSFs[0])*(1-bTagEffs[1]*bTagSFs[1])*(1-bTagEffs[2]*bTagSFs[2])*(1-bTagEffs[3]*bTagSFs[3]);
          float MCBTag1 = 1-(1-bTagEffs[0])*(1-bTagEffs[1])*(1-bTagEffs[2])*(1-bTagEffs[3]);
          //std::cout << " Full formula (2btag) 'Data' is " << DataBTag << " and 'MC' is " << MCBTag << " ==> Ratio is : " << DataBTag/MCBTag << std::endl;
          //std::cout << " Full formula (1btag) 'Data' is " << DataBTag1 << " and 'MC' is " << MCBTag1 << " ==> Ratio is : " << DataBTag1/MCBTag1 << std::endl; 
          BTagWeight2b = DataBTag/MCBTag;
          
          histo1D["BTagEfficiency_BJets"+dsName]->Fill(bTagEffs[bTagStudy.getBLeptIndex(ibTag)]);
          histo1D["BTagEfficiency_BJets"+dsName]->Fill(bTagEffs[bTagStudy.getBHadrIndex(ibTag)]);
          histo1D["BTagEfficiency_LightJets"+dsName]->Fill(bTagEffs[bTagStudy.getLight1Index(ibTag)]);
          histo1D["BTagEfficiency_LightJets"+dsName]->Fill(bTagEffs[bTagStudy.getLight2Index(ibTag)]);
        }
        histo1D["BTagWeight_AltMethod_"+dsName]->Fill(BTagWeight2b);
        histo1D["BTagWeight_"+dsName]->Fill(BTagWeight);

        if(bTagStudy.getBHadrIndex(ibTag) == 0 || bTagStudy.getBLeptIndex(ibTag) == 0){
          histo2D["bTagWeight_vs_BJetPt_"+dsName]->Fill(selJets[0].Pt(),BTagWeight);
          if(abs(partFlavour[0]) == 5)      histo2D["bTagWeight_vs_BJetPt_partFlavB_"+dsName]->Fill(selJets[0].Pt(),BTagWeight);
          else if(abs(partFlavour[0]) == 4) histo2D["bTagWeight_vs_BJetPt_partFlavC_"+dsName]->Fill(selJets[0].Pt(),BTagWeight);
          else                         histo2D["bTagWeight_vs_BJetPt_partFlavLight_"+dsName]->Fill(selJets[0].Pt(),BTagWeight);
        }
        if(bTagStudy.getLight1Index(ibTag) == 0 || bTagStudy.getLight2Index(ibTag) == 0){
          histo2D["bTagWeight_vs_LightJetPt_"+dsName]->Fill(selJets[0].Pt(),BTagWeight);
          if(abs(partFlavour[0]) == 5)      histo2D["bTagWeight_vs_LightJetPt_partFlavB_"+dsName]->Fill(selJets[0].Pt(),BTagWeight);
          else if(abs(partFlavour[0]) == 4) histo2D["bTagWeight_vs_LightJetPt_partFlavC_"+dsName]->Fill(selJets[0].Pt(),BTagWeight);
          else                         histo2D["bTagWeight_vs_LightJetPt_partFlavLight_"+dsName]->Fill(selJets[0].Pt(),BTagWeight);
        }
        if(abs(partFlavour[0]) != 4 && abs(partFlavour[0]) != 5) histo2D["bTagWeight_vs_JetPt_partFlavLight_"+dsName]->Fill(selJets[0].Pt(),BTagWeight);

        histo1D["partonFlav_BJets_"+dsName]->Fill(abs(partFlavour[bTagStudy.getBLeptIndex(ibTag)]));
        histo1D["partonFlav_BJets_"+dsName]->Fill(abs(partFlavour[bTagStudy.getBHadrIndex(ibTag)]));
        histo1D["partonFlav_LightJets_"+dsName]->Fill(abs(partFlavour[bTagStudy.getLight1Index(ibTag)]));
        histo1D["partonFlav_LightJets_"+dsName]->Fill(abs(partFlavour[bTagStudy.getLight2Index(ibTag)]));
        scaleFactor = scaleFactor * BTagWeight;

        if(partFlavour[bTagStudy.getBLeptIndex(ibTag)] == 21) MSPlot["BJetPt_partFlavGluon_BTagSF"+leptChannel]->Fill( selJets[bTagStudy.getBLeptIndex(ibTag)].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        if(partFlavour[bTagStudy.getBHadrIndex(ibTag)] == 21) MSPlot["BJetPt_partFlavGluon_BTagSF"+leptChannel]->Fill( selJets[bTagStudy.getBHadrIndex(ibTag)].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        for(int ii = 0; ii < selJets.size(); ii++){
          if(partFlavour[ii] == 21) MSPlot["JetsPt_partFlavGluon_BTagSF"+leptChannel]->Fill(selJets[ii].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        }
        MSPlot["nBJets_BTagSF"+leptChannel]->Fill(bTagStudy.getNrBTaggedJets(ibTag), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["nPV_BTagSF"+leptChannel]->Fill( nPrimVertices, datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["JetPt_LeadingJet_BTagSF"+leptChannel]->Fill( selJets[0].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["LeptonPt_BTagSF"+leptChannel]->Fill( selLepton.Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["JetEta_LeadingJet_BTagSF"+leptChannel]->Fill( selJets[0].Eta(), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["LeptonEta_BTagSF"+leptChannel]->Fill( selLepton.Eta(), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["LeptonCharge_BTagSF"+leptChannel]->Fill( leptCharge, datasets[iDataSet], true, Luminosity*scaleFactor);

        MSPlot["CSVDistr_BJets"+leptChannel]->Fill(bTagCSV[bTagStudy.getBHadrIndex(ibTag)], datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["CSVDistr_BJets"+leptChannel]->Fill(bTagCSV[bTagStudy.getBLeptIndex(ibTag)], datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["CSVDistr_LightJets"+leptChannel]->Fill(bTagCSV[bTagStudy.getLight1Index(ibTag)], datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["CSVDistr_LightJets"+leptChannel]->Fill(bTagCSV[bTagStudy.getLight2Index(ibTag)], datasets[iDataSet], true, Luminosity*scaleFactor);

        //Check whether the Mlb, Mqqb, Mt and MW depends a lot on the considered b-tag and on the fact whether it is applied!
        if(dsName.find("TTbarJets_SemiLept") == 0){
          if(correctJetCombi[0] != 9999 && correctJetCombi[1] != 9999 && correctJetCombi[2] != 9999 && correctJetCombi[3] != 9999){
            histo1D["Mlb_CorrTT_"+dsName+"_"+bTitle[ibTag]]->Fill( (selLepton+selJets[correctJetCombi[0]]).M());
            histo1D["HadrMTop_CorrTT_"+dsName+"_"+bTitle[ibTag]]->Fill( (selJets[correctJetCombi[1]]+selJets[correctJetCombi[2]]+selJets[correctJetCombi[3]]).M() );
            histo1D["HadrMW_CorrTT_"+dsName+"_"+bTitle[ibTag]]->Fill( (selJets[correctJetCombi[2]]+selJets[correctJetCombi[3]]).M());
            histo2D["MW_vs_MTop_CorrTT_"+dsName+"_"+bTitle[ibTag]]->Fill((selJets[correctJetCombi[1]]+selJets[correctJetCombi[2]]+selJets[correctJetCombi[3]]).M(),(selJets[correctJetCombi[2]]+selJets[correctJetCombi[3]]).M());
          }
        }

        //Identical MSPlots with number of jets information after requiring at least two b-jets and at least 2 light jets!
        MSPlot["nSelectedJets_"+bTitle[ibTag]+"_AfterBTag"+leptChannel]->Fill( selJets.size(),                    datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["nBTaggedJets_"+ bTitle[ibTag]+"_AfterBTag"+leptChannel]->Fill( bTagStudy.getNrBTaggedJets(ibTag), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["nLightJets_"+   bTitle[ibTag]+"_AfterBTag"+leptChannel]->Fill( bTagStudy.getNrLightJets(ibTag),   datasets[iDataSet], true, Luminosity*scaleFactor);      

        bool CutsSurvived = extraEvtSelCuts.KeepEvent(selLepton, selJets, selJetCombi, bTagStudy.getMlbMqqbChiSq(ibTag), CWUIndex, decayCh);

        //Now get the MSPlots after this additional cuts
        if( CutsSurvived ){
          MSPlot["nBJets_AllCuts"+leptChannel]->Fill(bTagStudy.getNrBTaggedJets(ibTag), datasets[iDataSet], true, Luminosity*scaleFactor);
          MSPlot["nPV_AllCuts"+leptChannel]->Fill( nPrimVertices, datasets[iDataSet], true, Luminosity*scaleFactor);
          MSPlot["JetPt_LeadingJet_AllCuts"+leptChannel]->Fill( selJets[0].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
          MSPlot["JetMass_LeadingJet_AllCuts"+leptChannel]->Fill( selJets[0].M(), datasets[iDataSet], true, Luminosity*scaleFactor);
          MSPlot["LeptonPt_AllCuts"+leptChannel]->Fill( selLepton.Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
          MSPlot["LeptonMass_AllCuts"+leptChannel]->Fill( selLepton.M(), datasets[iDataSet], true, Luminosity*scaleFactor);
          MSPlot["JetEta_LeadingJet_AllCuts"+leptChannel]->Fill( selJets[0].Eta(), datasets[iDataSet], true, Luminosity*scaleFactor);
          MSPlot["LeptonEta_AllCuts"+leptChannel]->Fill( selLepton.Eta(), datasets[iDataSet], true, Luminosity*scaleFactor);
          MSPlot["LeptonCharge_AllCuts"+leptChannel]->Fill( leptCharge, datasets[iDataSet], true, Luminosity*scaleFactor);

          if(partFlavour[bTagStudy.getBLeptIndex(ibTag)] == 21) MSPlot["BJetPt_partFlavGluon_AllCuts"+leptChannel]->Fill( selJets[bTagStudy.getBLeptIndex(ibTag)].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
          if(partFlavour[bTagStudy.getBHadrIndex(ibTag)] == 21) MSPlot["BJetPt_partFlavGluon_AllCuts"+leptChannel]->Fill( selJets[bTagStudy.getBHadrIndex(ibTag)].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
          for(int ii = 0; ii < selJets.size(); ii++){
            if(partFlavour[ii] == 21) MSPlot["JetsPt_partFlavGluon_AllCuts"+leptChannel]->Fill(selJets[ii].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
          }

          if(dsName.find("TTbarJets_SemiLept") == 0 && fillTFAfterCuts && createTFTree){  //Combine muon and electron channel!
            tfLight_mu = new TFLight();
            
            double fullScaleFactor = scaleFactor;   
            tfLight_mu->setFullScaleFactor(fullScaleFactor);
            tfLight_mu->setSelectedJets(selJets);
            tfLight_mu->setSelectedLepton(selLepton);
            tfLight_mu->setCorrectJetCombi(correctJetCombi);
	
            //Store the information needed for the TF (but only has value when dataset is ttbar)
            tfLight_mu->setGenVectorLight1( light->genVectorLight1() );
            tfLight_mu->setGenVectorLight2( light->genVectorLight2() );
            tfLight_mu->setGenVectorHadrB(  light->genVectorHadrB()  );
            tfLight_mu->setGenVectorLeptB(  light->genVectorLeptB()  );
            tfLight_mu->setGenVectorLepton( light->genVectorLepton() );

            TFLightTree->Fill();
            delete tfLight_mu;
            //----  End of Tree file filling (for TF's after evtSel)  ----//
          }
          /*else if(dsName.find("TTbarJets_SemiLept") == 0 && decayCh == 1 && createTFTree){
           tfLight_el = new TFLight();
            
            double fullScaleFactor = scaleFactor;   
            tfLight_el->setFullScaleFactor(fullScaleFactor);
            tfLight_el->setSelectedJets(selJets);
            tfLight_el->setSelectedLepton(selLepton);
            tfLight_el->setDecayChannel(decayCh);
            tfLight_el->setLeptonCharge(leptCharge);
            tfLight_el->setCorrectJetCombi(correctJetCombi);
            tfLight_el->setMET(MET);
	
            //Store the information needed for the TF (but only has value when dataset is ttbar)
            tfLight_el->setGenVectorLight1( light->genVectorLight1() );
            tfLight_el->setGenVectorLight2( light->genVectorLight2() );
            tfLight_el->setGenVectorHadrB(  light->genVectorHadrB()  );
            tfLight_el->setGenVectorLeptB(  light->genVectorLeptB()  );
            tfLight_el->setGenVectorLepton( light->genVectorLepton() );

            TFLightTree->Fill();
            delete tfLight_el;
            //----  End of Tree file filling (for TF's after evtSel)  ----//
          }*/
        }
 
        //Write out the LHCO output!
        if( getLHCOOutput && bTagChoiceMade){

          //Keep track of the original event number and the madweight numbers:
          EvtNrMatching << "\n    " << iEvt;
          if( iEvt < 10) EvtNrMatching << " "; if( iEvt < 100) EvtNrMatching << " "; if( iEvt < 1000) EvtNrMatching << " "; if( iEvt < 10000) EvtNrMatching << " "; if( iEvt < 100000) EvtNrMatching << " ";
          EvtNrMatching << "              " << CutsSurvived << "                " << fixed << setprecision(9)<< genCosTheta << "";
          if(genCosTheta > 0) EvtNrMatching << " ";

          //if(CutsSurvived) histo1D["CosTheta_SelEvts"]->Fill(kinFunctions.CosTheta(lhcoOutput.getGenLeptTop(), lhcoOutput.getGenLeptW(), lhcoOutput.getGenLepton()));
          // --> Cannot add this since the neutrino is not completely reconstructed and thus the leptonic top is not known ...
          if(CutsSurvived){
            lhcoOutput.StoreRecoInfo(selLepton, selJets, selJetCombi, decayCh, leptCharge, EvtNrMatching, CWUIndex);
            EvtNrMatching << "                 " << scaleFactor << std::endl;
          }
          
          if(CutsSurvived && genCosTheta != 2.0 && decayCh == 0) histo1D["CosTheta_Gen_SelEvts_"+dsName+"_mu"]->Fill(genCosTheta);
          if(CutsSurvived && genCosTheta != 2.0 && decayCh == 1) histo1D["CosTheta_Gen_SelEvts_"+dsName+"_el"]->Fill(genCosTheta);
          if(CutsSurvived && genCosTheta != 2.0) histo1D["CosTheta_Gen_SelEvts_"+dsName]->Fill(genCosTheta);
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
    if(dsName.find("TTbarJets_SemiLept") == 0 && getMassFits){
      std::cout << " -----------------     Fitting the different mass distributions    --------------------- \n" << endl;   
      histo1D["Mlb_CorrTT_"+dsName]->Fit("gaus","Q");     
      histo1D["HadrMTop_CorrTT_"+dsName]->Fit("gaus","Q");
      histo1D["HadrMW_CorrTT_"+dsName]->Fit("gaus","Q");
      histo1D["Mlb_CorrTT_"+dsName]->Fit("gaus","Q","", histo1D["Mlb_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(1)-histo1D["Mlb_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(2),
                                                        histo1D["Mlb_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(1)+histo1D["Mlb_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(2));
      histo1D["HadrMTop_CorrTT_"+dsName]->Fit("gaus","Q","", histo1D["HadrMTop_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(1)-histo1D["HadrMTop_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(2),
                                                             histo1D["HadrMTop_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(1)+histo1D["HadrMTop_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(2));
      histo1D["HadrMW_CorrTT_"+dsName]->Fit("gaus","Q","",histo1D["HadrMW_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(1)-histo1D["HadrMW_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(2),
                                                          histo1D["HadrMW_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(1)+histo1D["HadrMW_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(2));
      
      //Write out the mass values!
      std::cout << "   ** Mlb  = " << histo1D["Mlb_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(1)      << " +- " << histo1D["Mlb_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(2) << endl;
      std::cout << "   ** Mqqb = " << histo1D["HadrMTop_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(1) << " +- " << histo1D["HadrMTop_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(2) << endl;
      std::cout << "   ** MW   = " << histo1D["HadrMW_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(1)   << " +- " << histo1D["HadrMW_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(2) << std::endl;
    
      //Output for the different b-tags
      for(int ibTag = 0; ibTag < NrBTags; ibTag++){
        histo1D["Mlb_CorrTT_"+dsName+"_"+bTitle[ibTag]]->Fit("gaus","Q"); histo1D["HadrMTop_CorrTT_"+dsName+"_"+bTitle[ibTag]]->Fit("gaus","Q"); histo1D["HadrMW_CorrTT_"+dsName+"_"+bTitle[ibTag]]->Fit("gaus","Q");
        histo1D["Mlb_CorrTT_"+dsName+"_"+bTitle[ibTag]]->Fit("gaus","Q","", 
                                                      histo1D["Mlb_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(1)-histo1D["Mlb_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(2),
                                                      histo1D["Mlb_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(1)+histo1D["Mlb_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(2));
        histo1D["HadrMTop_CorrTT_"+dsName+"_"+bTitle[ibTag]]->Fit("gaus","Q","",
                                          histo1D["HadrMTop_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(1)-histo1D["HadrMTop_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(2),
                                          histo1D["HadrMTop_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(1)+histo1D["HadrMTop_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(2));
        histo1D["HadrMW_CorrTT_"+bTitle[ibTag]]->Fit("gaus","Q","",
                                               histo1D["HadrMW_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(1)-histo1D["HadrMW_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(2),
                                               histo1D["HadrMW_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(1)+histo1D["HadrMW_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(2));
      
        //Write out the mass values!
        cout << "\n   ** Mlb  -- " << bTitle[ibTag] << " = " << histo1D["Mlb_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(1)    << " +- " << histo1D["Mlb_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(2) << endl;
        cout << "   ** Mqqb -- " << bTitle[ibTag] << " = " << histo1D["HadrMTop_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(1) << " +- " << histo1D["HadrMTop_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(2) << endl;
        cout << "   ** MW   -- " << bTitle[ibTag] << " = " << histo1D["HadrMW_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(1)   << " +- " << histo1D["HadrMW_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(2) << std::endl;
      }
    }

    //--- Get output from bTagStudy class ---//
    if(dsName.find("TTbarJets_SemiLept") == 0) bTagStudy.ReturnBTagTable();
    bTagStudy.CreateHistograms(outputFile, savePDF, pathPNG, iDataSet);  //Security is added inside class such that MSPlots are only written when all datasets are considered!

    //--- Get output from LHCOOutput class ---//
    if(getLHCOOutput && bTagChoiceMade) lhcoOutput.WriteLHCOPlots(outputFile);

    //--- Get output form the ExtraEvtSel class ---//
    extraEvtSelCuts.StoreCutInfluence(outputFile);

    //---- Close the EventNrMatching output file for the considered dataset  ----//
    if(bTagChoiceMade && getLHCOOutput) EvtNrMatching.close();

    //---- Store the dataset specific 1D and 2D histograms ----//
    if(histo1D.size() > 0){
      TDirectory* th1dir = outputFile->GetDirectory("1D_histograms");   //Check whether directory already exists ..
      if(!th1dir) th1dir = outputFile->mkdir("1D_histograms");          // .. and otherwise create it!
      th1dir->cd();
      TDirectory* ds1Dir = th1dir->mkdir(dsName.c_str());
      ds1Dir->cd();
      for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++){
        TH1F *temp = it->second;
        int N = temp->GetNbinsX();
        temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
        temp->SetBinContent(N+1,0);
        temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
        temp->Write();
      }
    }

    if(histo2D.size() > 0){
      TDirectory* th2dir = outputFile->GetDirectory("2D_histograms_graphs");
      if(!th2dir) th2dir = outputFile->mkdir("2D_histograms_graphs");
      th2dir->cd();
      TDirectory* ds2Dir = th2dir->mkdir(dsName.c_str());
      ds2Dir->cd();
      for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++){
        TH2F *temp = it->second;
        temp->Write();
      }
    }
    histo1D.clear();
    histo2D.clear();

    //------ Store the TF tree ------//
    if(dsName.find("TTbarJets_SemiLept") == 0 && createTFTree){
      TFLightFile->cd();
      TTree* configTreeTFLightFile = new TTree("configTreeTFLightFile","configuration Tree in TFLight File");
      TClonesArray* tcdatasettflightfile = new TClonesArray("Dataset",1);
      configTreeTFLightFile->Branch("Dataset","TClonesArray",&tcdatasettflightfile);
      new ((*tcdatasettflightfile)[0]) Dataset(*datasets[iDataSet]);
      //TClonesArray* tcanaenvlightfile = new TClonesArray("AnalysisEnvironment",1);   --> Needed?
      //configTreeLightFile->Branch("AnaEnv","TClonesArray",&tcanaenvlightfile);
      //new ((*tcanaenvlightfile)[0]) AnalysisEnvironment(anaEnv);

      configTreeTFLightFile->Fill();
      configTreeTFLightFile->Write();
      TFLightTree->Write();
      TFLightFile->Close();
      //----  End of storing Tree  ----//
    }
    delete TFLightFile;  //Will need to delete this for each dataset since it is initialized for each of them!

    inputFile->Close();
    delete inputFile;
  }//End of loop on datasets

  //Store the b-tag histograms:
  if(!bTagPlotsMade)
    bTagTool->WriteMCEfficiencyHistos("PersonalClasses/Calibrations/BTagSF/BTagWeightPlots_CSVT_noTTbar.root");

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

  delete leptonTools;
  delete outputFile;

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "           hasn't crashed yet ;-)           " << endl;
  cout << "********************************************" << endl;
  
  return 0;
}
