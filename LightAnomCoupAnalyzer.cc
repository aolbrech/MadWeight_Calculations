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

int main (int argc, char *argv[])
{
  clock_t start = clock();
  cout << "************************************************************" << endl;
  cout << "  Beginning of the program for analyzing the Light Trees !  " << endl;
  cout << "************************************************************" << endl;

  //SetStyle if needed
  //setTDRStyle(); 
  setMyStyle();

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

  //----------------------------//
  // Input & output information //
  //----------------------------//
  //ROOT file for storing plots
  string pathPNG = "PlotsMacro_Light";
  TFile* outputFile = new TFile((pathPNG+"/AnomCoup_Analysis_"+numberBTags+".root").c_str(),"RECREATE");

  //Which datasets should be considered
  vector<string> inputFiles;
  vector<Dataset*> datasets;
  inputFiles.push_back("LightTree/AnomCoupLight_TTbarJets_SemiLept_AllTTbarEvents_19Aug2015.root");
  //inputFiles.push_back("LightTree/AnomalousCouplingsLight_Data_Mu.root");
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
    if(dataSet->Name().find("TT") == 0)              color = kRed+1;
    if(dataSet->Name().find("TTbarJets_Other") == 0) color = kRed-7;
    if(dataSet->Name().find("ST") == 0 || 
       dataSet->Name().find("SingleTop") == 0)       color = kMagenta;
    if(dataSet->Name().find("WJets") == 0){          color = kGreen-3; dataSet->SetTitle("W#rightarrowl#nu");}      
    if(dataSet->Name().find("ZJets") == 0){          color = kAzure-2; dataSet->SetTitle("Z/#gamma*#rightarrowl^{+}l^{-}");}

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
 
      MSPlot["nSelectedJets_"+bTitle[ibTag]+"_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "nSelectedJets_"+bTitle[ibTag]+"_BeforeBTag"+leptFlav,10, -0.5, 9.5, "# selected jets");
      MSPlot["nBTaggedJets_"+ bTitle[ibTag]+"_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "nBTaggedJets_"+ bTitle[ibTag]+"_BeforeBTag"+leptFlav,10, -0.5, 9.5, "# b-tagged jets");
      MSPlot["nLightJets_"+   bTitle[ibTag]+"_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "nLightJets_"+   bTitle[ibTag]+"_BeforeBTag"+leptFlav,10, -0.5, 9.5, "# light jets");
      MSPlot["nSelectedJets_"+bTitle[ibTag]+"_AfterBTag"+ leptFlav] = new MultiSamplePlot(datasets, "nSelectedJets_"+bTitle[ibTag]+"_AfterBTag"+ leptFlav,10, -0.5, 9.5, "# selected jets");
      MSPlot["nBTaggedJets_"+ bTitle[ibTag]+"_AfterBTag"+ leptFlav] = new MultiSamplePlot(datasets, "nBTaggedJets_"+ bTitle[ibTag]+"_AfterBTag"+ leptFlav,10, -0.5, 9.5, "# b-tagged jets");
      MSPlot["nLightJets_"+   bTitle[ibTag]+"_AfterBTag"+ leptFlav] = new MultiSamplePlot(datasets, "nLightJets_"+   bTitle[ibTag]+"_AfterBTag"+ leptFlav,10, -0.5, 9.5, "# light jets");
    }
  }

  //-----------------------//
  // Initialize histograms //
  //-----------------------//
  map<string,TH1F*> histo1D;
  map<string,TH2F*> histo2D;  

  histo1D["Mlb_CorrectTTEvts"] = new TH1F("Mlb_CorrectTTEvts", "mass_{l,b} for the actual semi-leptonic ttbar events", 200, 0, 250);
  histo1D["TopMass_Hadr_CorrectTTEvts"] = new TH1F("TopMass_Hadr_CorrectTTEvts", "Hadronic m_{top} for the actual semi-leptonic ttbar events", 200, 50, 350);
  histo1D["WMass_Hadr_CorrectTTEvts"] = new TH1F("WMass_Hadr_CorrectTTEvts", "Hadronic m_{W} for the actual semi-leptonic ttbar events", 200, 0, 250);
  histo2D["WMass_vs_TopMass_CorrectTTEvts"] = new TH2F("WMass_vs_TopMass_CorrectTTEvts","Hadronic m_{W} versus m_{top} for the actual semi-leptonic ttbar jet combinations",200,50,350,200,0,250);

  for(int ibTag = 0; ibTag < NrBTags; ibTag++){
    histo1D["Mlb_CorrectTTEvts_"+bTitle[ibTag]]          = new TH1F(("Mlb_CorrectTTEvts_"+bTitle[ibTag]).c_str(),          ("mass_{l,b} for the actual semi-leptonic ttbar events -- "+bTitle[ibTag]).c_str(), 200, 0, 250);
    histo1D["TopMass_Hadr_CorrectTTEvts_"+bTitle[ibTag]] = new TH1F(("TopMass_Hadr_CorrectTTEvts_"+bTitle[ibTag]).c_str(), ("Hadronic m_{top} for the actual semi-leptonic ttbar events -- "+bTitle[ibTag]).c_str(), 200, 50, 350);
    histo1D["WMass_Hadr_CorrectTTEvts_"+bTitle[ibTag]]   = new TH1F(("WMass_Hadr_CorrectTTEvts_"+bTitle[ibTag]).c_str(),   ("Hadronic m_{W} for the actual semi-leptonic ttbar events -- "+bTitle[ibTag]).c_str(), 200, 0, 250);
    histo2D["WMass_vs_TopMass_CorrectTTEvts_"+bTitle[ibTag]] = new TH2F(("WMass_vs_TopMass_CorrectTTEvts_"+bTitle[ibTag]).c_str(), ("Hadronic m_{W} versus m_{top} for actual semi-lept ttbar events -- "+bTitle[ibTag]).c_str(), 200, 50, 350, 200, 0,250);
  }

  //--------------------------------//
  // Lumi reweighting and lepton SF //
  //--------------------------------//
  
  //-----------------------//
  // Load personal classes //
  //-----------------------//
  BTagStudy bTagStudy(verbosity, datasets, bTagChoiceMade, ChosenBTag);
  ExtraEvtSelCuts extraEvtSelCuts;
  float Mlb  = 108.1841, S_Mlb  = 31.4213;
  float Mqqb = 174.6736, S_Mqqb = 17.5757;
  float MW = 83.8037, S_MW = 10.2385;
  //float Mlb  = 103.286, S_Mlb  = 26.7764;    --> OLD values (but used for CorrectReco file!!)
  //float Mqqb = 178.722, S_Mqqb = 18.1385;

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
    int nEvent = inLightTree->GetEntries();
    //int nEvent = 50000;
  
    Dataset* dataSet = datasets[iDataSet];//(Dataset*) tc_dataset->At(0);
    string dataSetName = dataSet->Name();
    if(verbosity > 0) std::cout << "   *** Looking at dataset "<< dataSetName << " (" << iDataSet+1 << "/" << inputFiles.size() << ") with " << nEvent << " selected events! " << std::endl;

    //-----------------------//
    // Load personal classes //
    //-----------------------//
    bTagStudy.InitializeDataSet(dataSetName);
    LHCOOutput lhcoOutput(verbosity, getLHCOOutput, splitLeptonChargeLHCO, getCorrectAndWrongLHCO);
    if(dataSetName.find("TTbarJets") == 0){
      lhcoOutput.Initialize("Reco");
      extraEvtSelCuts.Initialize(Mqqb, S_Mqqb, MW, S_MW, bTagChoiceMade, bTitle[0], 3, 2);
    }

    ofstream EvtNrMatching;
    EvtNrMatching.open(("MadWeightInput/AnalyzerOutput/EventNrMatching_"+dataSetName+".txt").c_str());
    EvtNrMatching << "  Event Nr     Extra cuts survived     Main lhco file      MW Number (main)      TTbar splitting lhco file      MW Number (splitting)   " << endl;


    // --------------------------- //
    //  Start looping over events  //
    // --------------------------- //    
    int nSelectedMu = 0, nSelectedEl = 0;
    for(unsigned int iEvt = 0; iEvt < nEvent; iEvt++){
      inLightTree->GetEvent(iEvt);

      if(iEvt%5000 == 0)
	std::cout<<"    Processing the "<<iEvt<<"th event ("<< ((double)iEvt/(double)nEvent)*100<<"%)"<<" -> # selected: "<<nSelectedMu<<" (mu+jets) "<<nSelectedEl<<" (e+jets)"<< flush<<"\r";

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

      //ooooooooOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOooooooooooooOOOOO
      //ooOOooOOoo      Reading out nTuples done           ooOOooOOoo
      //ooOOooOOoo-----------------------------------------ooOOooOOoo
      //ooOOooOOoo      Start of actual analysis           ooOOooOOoo
      //ooooooooOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOooooooooooooOOOOO
      bTagStudy.CalculateJets(selJets, bTagCSV, correctJetCombi, selLepton, datasets[iDataSet], Luminosity*scaleFactor);
 
      //Store the expected Mlb, Mqqb, Mt and MW distributions
      if(dataSetName.find("TTbarJets_SemiLept") == 0){
        if(correctJetCombi[0] != 9999 && correctJetCombi[1] != 9999 && correctJetCombi[2] != 9999 && correctJetCombi[3] != 9999){
          histo1D["Mlb_CorrectTTEvts"]->Fill( (selLepton+selJets[correctJetCombi[0]]).M());
          histo1D["TopMass_Hadr_CorrectTTEvts"]->Fill( (selJets[correctJetCombi[1]]+selJets[correctJetCombi[2]]+selJets[correctJetCombi[3]]).M() );
          histo1D["WMass_Hadr_CorrectTTEvts"]->Fill( (selJets[correctJetCombi[2]]+selJets[correctJetCombi[3]]).M());
          histo2D["WMass_vs_TopMass_CorrectTTEvts"]->Fill( (selJets[correctJetCombi[1]]+selJets[correctJetCombi[2]]+selJets[correctJetCombi[3]]).M(), (selJets[correctJetCombi[2]]+selJets[correctJetCombi[3]]).M() );
        }
      }
 
      for(int ibTag = 0; ibTag < NrBTags; ibTag++){
        vector<int> selJetCombi = bTagStudy.getIndices(ibTag);

        //MSPlots with number of jets information before requiring at least two b-jets and at least 2 light jets!
        MSPlot["nSelectedJets_"+bTitle[ibTag]+"_BeforeBTag"+leptChannel]->Fill( selJets.size(),                    datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["nBTaggedJets_"+ bTitle[ibTag]+"_BeforeBTag"+leptChannel]->Fill( bTagStudy.getNrBTaggedJets(ibTag), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["nLightJets_"+   bTitle[ibTag]+"_BeforeBTag"+leptChannel]->Fill( bTagStudy.getNrLightJets(ibTag),   datasets[iDataSet], true, Luminosity*scaleFactor);

        //Apply the event selection
        if( bTagStudy.getNrBTaggedJets(ibTag) < 2 || bTagStudy.getNrLightJets(ibTag) < 2 ) continue;
        if(decayCh == 0) nSelectedMu += 1;
        else if(decayCh == 1) nSelectedEl += 1;

        //Check whether the Mlb, Mqqb, Mt and MW depends a lot on the considered b-tag and on the fact whether it is applied!
        if(dataSetName.find("TTbarJets_SemiLept") == 0){
          if(correctJetCombi[0] != 9999 && correctJetCombi[1] != 9999 && correctJetCombi[2] != 9999 && correctJetCombi[3] != 9999){
            histo1D["Mlb_CorrectTTEvts_"+bTitle[ibTag]]->Fill( (selLepton+selJets[correctJetCombi[0]]).M());
            histo1D["TopMass_Hadr_CorrectTTEvts_"+bTitle[ibTag]]->Fill( (selJets[correctJetCombi[1]]+selJets[correctJetCombi[2]]+selJets[correctJetCombi[3]]).M() );
            histo1D["WMass_Hadr_CorrectTTEvts_"+bTitle[ibTag]]->Fill( (selJets[correctJetCombi[2]]+selJets[correctJetCombi[3]]).M());
            histo2D["WMass_vs_TopMass_CorrectTTEvts_"+bTitle[ibTag]]->Fill((selJets[correctJetCombi[1]]+selJets[correctJetCombi[2]]+selJets[correctJetCombi[3]]).M(),(selJets[correctJetCombi[2]]+selJets[correctJetCombi[3]]).M());
          }
        }

        EvtNrMatching << "\n    " << iEvt;
        if( iEvt < 10) EvtNrMatching << " "; if( iEvt < 100) EvtNrMatching << " "; if( iEvt < 1000) EvtNrMatching << " "; if( iEvt < 10000) EvtNrMatching << " "; if( iEvt < 100000) EvtNrMatching << " ";

        //Identical MSPlots with number of jets information after requiring at least two b-jets and at least 2 light jets!
        MSPlot["nSelectedJets_"+bTitle[ibTag]+"_AfterBTag"+leptChannel]->Fill( selJets.size(),                    datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["nBTaggedJets_"+ bTitle[ibTag]+"_AfterBTag"+leptChannel]->Fill( bTagStudy.getNrBTaggedJets(ibTag), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["nLightJets_"+   bTitle[ibTag]+"_AfterBTag"+leptChannel]->Fill( bTagStudy.getNrLightJets(ibTag),   datasets[iDataSet], true, Luminosity*scaleFactor);      

        bool CutsSurvived = extraEvtSelCuts.KeepEvent(correctJetCombi, selLepton, selJets, selJetCombi, bTagStudy.getMlbMqqbChiSq(ibTag));
        EvtNrMatching << "              " << CutsSurvived << "      ";
 
        //Write out the LHCO output!
        if( getLHCOOutput && bTagChoiceMade)
          lhcoOutput.StoreRecoInfo(selLepton, selJets, selJetCombi, decayCh, leptCharge, correctJetCombi, EvtNrMatching);
        else if(getLHCOOutput && !bTagChoiceMade)
          cout << " ERROR : Not possible to write the lhco file when all b-tag options are still being considered!! " << endl;
      }

    }//End of loop on events
    std::cout<<"    Processed all "<<nEvent<<" events  --> # selected: "<<nSelectedMu<<" (mu+jets) "<<nSelectedEl<<" (e+jets)"<< flush<<"\n";

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
      
      //Write out the sigma values!
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
      
        //Write out the sigma values!
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
    if(dataSetName.find("TTbarJets_SemiLept") == 0) extraEvtSelCuts.StoreCutInfluence(outputFile);

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
    temp->Draw(name, 0, false, false, false, 1);     //string label, unsigned int RatioType, bool addRatioErrorBand, bool addErrorBand, bool ErrorBandAroundTotalInput, int scaleNPSSignal 
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
