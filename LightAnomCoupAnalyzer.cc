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

#include "TLorentzVector.h"

#include "TROOT.h"
#include "TH1F.h"
#include "PersonalClasses/Style.C"                                                 //CHECK if this works!

#include "AnomalousCouplings/PersonalClasses/interface/AnomCoupLight.h"
#include "AnomalousCouplings/PersonalClasses/interface/BTagStudy.h"

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

  int verbose = 1;

  TFile* outputFile = new TFile("LightAnomCoup_Analysis.root","RECREATE");

  vector<string> inputFiles;
  vector<Dataset*> datasets;
  inputFiles.push_back("LightTree/AnomalousCouplingsLight_TTbarJets_SemiLept.root");

  BTagStudy bTagStudy(verbose);
  //-----------------------//
  // Initialize histograms //
  //-----------------------//

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
    //EqLumi = dataSet->EquivalentLumi();        

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
    std::cout << " *** Looking at dataset " << iDataSet+1 << "/" << inputFiles.size() << " with " << nEvent << " selected events! \n " << std::endl;

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
      vector<int> correctJetCombi = light->correctJetCombi();    //0 = LeptB, 1 = HadrB, 2 = Quark1 & 3 = Quark2

      //ooooooooOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOooooooooooooOOOOO
      //ooOOooOOoo      Reading out nTuples done           ooOOooOOoo
      //ooOOooOOoo-----------------------------------------ooOOooOOoo
      //ooOOooOOoo      Start of actual analysis           ooOOooOOoo
      //ooooooooOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOooooooooooooOOOOO
      bTagStudy.CalculateJets(selectedJets, bTagCSV, correctJetCombi, selectedLepton);
      float bTagLeptIndex = bTagStudy.getBLeptIndex(3);          //-->Filled with which value in case only 4 jets are considered ... ?
      std::cout << " bTagLeptIndex obtained from bTagStudy is : " << bTagLeptIndex << std::endl;
  
    }//End of loop on events

  }//End of loop on datasets

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

      //----------------------------------------------------------------------------------------------------------------------------------- Start of bTagStudy class stuff!!      	
      //---  Get the b-jet and light information for all six b-tag options!  ---//
      bTagStudy.CalculateJets(selectedJets, jetCombi, selectedLepton);
      float bTagLeptIndex = bTagStudy.getBLeptIndex(3);
      //---------------------------------------------------------------------------------------------------------------------------- End of bTagStudy class stuff
      
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

    //--- Get output from bTagStudy class ---//
    bTagStudy.ReturnBTagTable();
    bTagStudy.CreateHistograms(fout);

    //--- Get output from LHCOOutput class ---//
    lhcoOutput.WriteLHCOPlots(fout);	
*/
