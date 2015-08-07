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
#include "AnomalousCouplings/PersonalClasses/interface/BTagStudy.h"
#include "AnomalousCouplings/PersonalClasses/interface/LHCOOutput.h"

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
  //0 muet
  //1 Main Info
  //2 mcParticlesMatchin Info
  //3 
  //4 Info for each event
  //5 Debug

  //----------------------------//
  // Input & output information //
  //----------------------------//
  //ROOT file for storing plots
  string pathPNG = "PlotsMacro_Light";
  TFile* outputFile = new TFile((pathPNG+"/AnomCoup_Analysis.root").c_str(),"RECREATE");
  std::cout << " Does MSPlot directory already exists at the start ?? --> " << outputFile->GetDirectory("MSPlots") << std::endl;

  //Which datasets should be considered
  vector<string> inputFiles;
  vector<Dataset*> datasets;
  inputFiles.push_back("LightTree/AnomalousCouplingsLight_TTbarJets_SemiLept.root");
  inputFiles.push_back("LightTree/AnomalousCouplingsLight_Data_Mu.root");

  //Which parts of the analysis should be performed
  bool getLHCOOutput = true;
  bool savePDF = false;

  //-------------------------//
  // Set analysis luminosity //
  //-------------------------//
  float Luminosity = 19646.840; //IsoMu24_eta2p1 trigger for semiMu case!  
  std::cout << " Analysis luminosity is set to : " << Luminosity << std::endl;


  //-----------------------//
  // Initialize histograms //
  //-----------------------//
  map<string,TH1F*> histo1D;
  map<string,TH2F*> histo2D;  

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

  //--------------------------//
  //  Initialize MSPlots      //
  //  --> Needs dataset info  //
  //--------------------------//
  map<string,MultiSamplePlot*> MSPlot;
  string leptFlavs[2]={"_mu","_el"};
  for(int ii = 0; ii < 2; ii++){
    string leptFlav = leptFlavs[ii];
 
    MSPlot["nSelectedJets_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "nSelectedJets_BeforeBTag"+leptFlav,14, -3.5, 10.5, "# selected jets");
    MSPlot["nBTaggedJets_BeforeBTag"+leptFlav]  = new MultiSamplePlot(datasets, "nBTaggedJets_BeforeBTag"+leptFlav, 14, -3.5, 10.5, "# b-tagged jets");
    MSPlot["nLightJets_BeforeBTag"+leptFlav]    = new MultiSamplePlot(datasets, "nLightJets_BeforeBTag"+leptFlav,   14, -3.5, 10.5, "# light jets");
    MSPlot["nSelectedJets_AfterBTag"+leptFlav]  = new MultiSamplePlot(datasets, "nSelectedJets_AfterBTag"+leptFlav, 14, -3.5, 10.5, "# selected jets");
    MSPlot["nBTaggedJets_AfterBTag"+leptFlav]   = new MultiSamplePlot(datasets, "nBTaggedJets_AfterBTag"+leptFlav,  14, -3.5, 10.5, "# b-tagged jets");
    MSPlot["nLightJets_AfterBTag"+leptFlav]     = new MultiSamplePlot(datasets, "nLightJets_AfterBTag"+leptFlav,    14, -3.5, 10.5, "# light jets");
  }

  //--------------------------------//
  // Lumi reweighting and lepton SF //
  //--------------------------------//

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
    //int nEvent = 10;
  
    Dataset* dataSet = datasets[iDataSet];//(Dataset*) tc_dataset->At(0);
    string dataSetName = dataSet->Name();
    std::cout << " *** Looking at dataset "<< dataSetName << " (" << iDataSet+1 << "/" << inputFiles.size() << ") with " << nEvent << " selected events! \n " << std::endl;

    //-----------------------//
    // Load personal classes //
    //-----------------------//
    BTagStudy bTagStudy(verbosity, datasets);
    LHCOOutput lhcoOutput(verbosity, getLHCOOutput);
    if(dataSetName.find("TTbarJets") == 0) lhcoOutput.Initialize("Reco");

    for(unsigned int iEvt = 0; iEvt < nEvent; iEvt++){
      inLightTree->GetEvent(iEvt);	

      int eventId = light->eventID();
      int runId = light->runID();
      int nTruePU = light->nTruePU();
      int nPrimVertices = light->nPV();
      float scaleFactor = light->scaleFactor();

      vector<float> bTagCSV = light->CSVbTag();
      vector<TLorentzVector> selectedJets = light->selectedJets();
      TLorentzVector selectedLepton = light->selectedLepton();
      TLorentzVector MET = light->met();
      int decayChannel = light->decayChannel();  //0 = semiMu and 1 = semiEl
      std::string leptChannel;
      if(decayChannel == 0) leptChannel = "_mu";
      else if(decayChannel == 1) leptChannel = "_el";
      float leptonCharge = light->leptonCharge();
      vector<int> correctJetCombi = light->correctJetCombi();    //0 = LeptB, 1 = HadrB, 2 = Quark1 & 3 = Quark2

      //ooooooooOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOooooooooooooOOOOO
      //ooOOooOOoo      Reading out nTuples done           ooOOooOOoo
      //ooOOooOOoo-----------------------------------------ooOOooOOoo
      //ooOOooOOoo      Start of actual analysis           ooOOooOOoo
      //ooooooooOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOooooooooooooOOOOO
      bTagStudy.CalculateJets(selectedJets, bTagCSV, correctJetCombi, selectedLepton, datasets[iDataSet], Luminosity*scaleFactor);
      float bTagLeptIndex = bTagStudy.getBLeptIndex(3);

      //---------------------------------//
      //  Decide on the event selection  //
      //---------------------------------//
      int ChosenBTag = 3;  //-->Corresponds to two T b-tags and no light veto!

      //MSPlots with number of jets information before requiring at least two b-jets and at least 2 light jets!
      MSPlot["nSelectedJets_BeforeBTag"+leptChannel]->Fill( selectedJets.size(),                           datasets[iDataSet], true, Luminosity*scaleFactor);
      MSPlot["nBTaggedJets_BeforeBTag"+leptChannel]->Fill(  (bTagStudy.getbTaggedJets(ChosenBTag)).size(), datasets[iDataSet], true, Luminosity*scaleFactor);
      MSPlot["nLightJets_BeforeBTag"+leptChannel]->Fill(    (bTagStudy.getLightJets(ChosenBTag)).size(),   datasets[iDataSet], true, Luminosity*scaleFactor);

      //Apply the event selection
      if( (bTagStudy.getbTaggedJets(ChosenBTag)).size() < 2 || (bTagStudy.getLightJets(ChosenBTag)).size() < 2 ) continue;

      //Identical MSPlots with number of jets information after requiring at least two b-jets and at least 2 light jets!
      MSPlot["nSelectedJets_AfterBTag"+leptChannel]->Fill( selectedJets.size(),                           datasets[iDataSet], true, Luminosity*scaleFactor);
      MSPlot["nBTaggedJets_AfterBTag"+leptChannel]->Fill(  (bTagStudy.getbTaggedJets(ChosenBTag)).size(), datasets[iDataSet], true, Luminosity*scaleFactor);
      MSPlot["nLightJets_AfterBTag"+leptChannel]->Fill(    (bTagStudy.getLightJets(ChosenBTag)).size(),   datasets[iDataSet], true, Luminosity*scaleFactor);      
 
      //Write out the LHCO output!
      if( getLHCOOutput == true && dataSetName.find("TTbarJets") == 0)
        lhcoOutput.StoreRecoInfo(selectedLepton, selectedJets, bTagStudy.getBLeptIndex(ChosenBTag), bTagStudy.getBHadrIndex(ChosenBTag), bTagStudy.getLight1Index5Jets(ChosenBTag), bTagStudy.getLight2Index5Jets(ChosenBTag), decayChannel, leptonCharge, correctJetCombi); 
    }//End of loop on events

    //--- Get output from bTagStudy class ---//
    bTagStudy.ReturnBTagTable(dataSetName);
    bTagStudy.CreateHistograms(outputFile, savePDF, pathPNG, iDataSet);

    //--- Get output from LHCOOutput class ---//
    if(dataSetName.find("TTbarJets") == 0) lhcoOutput.WriteLHCOPlots(outputFile);

  }//End of loop on datasets

  std::cout << " Does MSPlots directory exist (in LightAnomCoupAnalyzer) --> " << outputFile->GetDirectory("MSPlots") << std::endl;
  /////////////////////////
  // Write out the plots //
  /////////////////////////
  outputFile -> cd();
  mkdir((pathPNG+"/MSPlots").c_str(),0777);

  std::cout << " Making the plots in the LightAnomCoupAnalyzer file " << std::endl;
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
	  histo1D["MlbMass"]->Fill((*selectedJets[jetCombi[0]]+*selectedLepton).M());
	  histo1D["MqqbMass"]->Fill((*selectedJets[jetCombi[2]]+*selectedJets[jetCombi[3]]+*selectedJets[jetCombi[1]]).M());
        }

        ////////////////////////////////
        //  Mlb and Mqqb information  //
        ////////////////////////////////
        float MlbCorrect = 0, MqqbCorrect = 0;
        if(jetCombi[0] != 9999 && jetCombi[1] != 9999 && jetCombi[2] != 9999 && jetCombi[3] != 9999){
	  MlbCorrect = (*selectedLepton+*selectedJets[jetCombi[0]]).M();
	  MqqbCorrect = (*selectedJets[jetCombi[1]] + *selectedJets[jetCombi[2]] + *selectedJets[jetCombi[3]]).M();
	  histo2D["MlbMqqbCorrectAll"]->Fill(MqqbCorrect,MlbCorrect);
        }

      // Count the number of events:
      if(decayChannel == 0) nSelectedMu++;
      if(decayChannel == 1) nSelectedEl++;
      
      ////////////////////////////////   
      //  Produce Reco LHCO Output  //
      //  --> Last integer = mode   //
      //       * 0 ~ all events     //
      //       * 1 ~ good combi's   //
      //       * 2 ~ bad combi's    //
      ////////////////////////////////
      if( RecoLHCOOutput == true)
        lhcoOutput.StoreRecoInfo(selectedLepton, selectedJets, bTagStudy.getBLeptIndex(ChosenBTag), bTagStudy.getBHadrIndex(ChosenBTag), bTagStudy.getLight1Index5Jets(ChosenBTag), bTagStudy.getLight2Index5Jets(ChosenBTag), decayChannel, LeptonRecoCharge, jetCombi);

    //---  Mlb and Mqqb fit result ---//
    histo1D["MlbMass"]->Fit("gaus","Q");
    histo1D["MqqbMass"]->Fit("gaus","Q");
    cout <<"\n values for Mlb :"<< histo1D["MlbMass"]->GetFunction("gaus")->GetParameter(1) << " +- " << histo1D["MlbMass"]->GetFunction("gaus")->GetParameter(2) << endl;
    cout <<" values for Mqqb :" << histo1D["MqqbMass"]->GetFunction("gaus")->GetParameter(1) << " +- " << histo1D["MqqbMass"]->GetFunction("gaus")->GetParameter(2) << endl;

*/
