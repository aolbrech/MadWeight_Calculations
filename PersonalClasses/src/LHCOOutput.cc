#include <iomanip>
#include "../interface/LHCOOutput.h"

#include <fstream>
#include <sstream>

LHCOOutput::LHCOOutput(int verbosity, bool GenOutput, bool RecoOutput){
    
  //Constructor: Define variables which have to be initialized in the beginning!
  NumberNegativeElectrons = 0; NumberNegativeMuons = 0; NumberPositiveElectrons = 0; NumberPositiveMuons = 0;
  WrongEvtCounter = 0;
  CorrectGenEvtContent = false;
  NumberNegRecoEl = 0;         NumberNegRecoMu = 0;     NumberPosRecoEl = 0;         NumberPosRecoMu = 0; 
  NrPosRecoMuCorrect = 0; NrPosRecoMuWrong = 0; 
  genOutput_ = GenOutput; recoOutput_ = RecoOutput;
  verbose_ = verbosity;

  if(genOutput_){
    GenOutFile[0].open("MadWeightInput/AnalyzerOutput/TTbarLHCO_PositiveMuon.lhco");
    GenOutFile[1].open("MadWeightInput/AnalyzerOutput/TTbarLHCO_NegativeMuon.lhco");
    GenOutFile[2].open("MadWeightInput/AnalyzerOutput/TTbarLHCO_PositiveElectron.lhco");
    GenOutFile[3].open("MadWeightInput/AnalyzerOutput/TTbarLHCO_NegativeElectron.lhco");
    WrongGenFile.open("MadWeightInput/AnalyzerOutput/WrongGeneratorEvents.lhco");
  }

  if(recoOutput_){
    RecoOutFile[0].open("MadWeightInput/AnalyzerOutput/TTbarSemiLepton_Reco_PositiveMuon.lhco");
    RecoOutFile[1].open("MadWeightInput/AnalyzerOutput/TTbarSemiLepton_Reco_NegativeMuon.lhco");
    RecoOutFile[2].open("MadWeightInput/AnalyzerOutput/TTbarSemiLepton_Reco_PositiveElec.lhco");
    RecoOutFile[3].open("MadWeightInput/AnalyzerOutput/TTbarSemiLepton_Reco_NegativeElec.lhco");
    WrongRecoMuPosFile.open("MadWeightInput/AnalyzerOutput/WrongRecoEvents_PositiveMuon.lhco");
    CorrectRecoMuPosFile.open("MadWeightInput/AnalyzerOutput/CorrectRecoEvents_PositiveMuon.lhco");
  }

  histo1D["WrongGen_TopMassLept"] = new TH1F("WrongGen_TopMassLept","Leptonic top-mass distribution for wrong generator events",250,0,500);
  histo1D["WrongGen_TopMassHadr"] = new TH1F("WrongGen_TopMassHadr","Hadronic top-mass distribution for wrong generator events",250,0,500);
  histo1D["WrongReco_TopMassLept"] = new TH1F("WrongReco_TopMassLept","Leptonic top-mass distribution for wrong reco events",250,0,500);
  histo1D["WrongReco_TopMassHadr"] = new TH1F("WrongReco_TopMassHadr","Hadronic top-mass distribution for wrong reco events",250,0,500);
  histo1D["CorrectReco_TopMassLept"] = new TH1F("CorrectReco_TopMassLept","Leptonic top-mass distribution for correct reco events",250,0,500);
  histo1D["CorrectReco_TopMassHadr"] = new TH1F("CorrectReco_TopMassHadr","Hadronic top-mass distribution for correct reco events",250,0,500);
}

LHCOOutput::~LHCOOutput(){

  //Close the output files:
  for(int ii = 0; ii < 4; ii++){
    if(genOutput_) GenOutFile[ii].close();
    if(recoOutput_) RecoOutFile[ii].close();
  }
  if(genOutput_) WrongGenFile.close();
  if(recoOutput_){ WrongRecoMuPosFile.close(); CorrectRecoMuPosFile.close();}

}

void LHCOOutput::StoreGenInfo(vector<TRootMCParticle*> mcParticles){

  //Initialize EventContent counter used for selecting correct particle content!
  int EventContent[5]; //0:top; 1:b; 2: u,c,d,s; 3:W; 4:mu + neutrino
  for(int ll = 0;ll<5;ll++){EventContent[ll]=0;}
 
  //Loop over all the mcParticles
  for(unsigned int i=0; i<mcParticles.size(); i++){
    if( mcParticles[i]->status() != 3) continue;
	
    int partType=mcParticles[i]->type(); if(verbose_>4)cout<<"-->Type of mcParticle : "<<partType<<endl;
    //if( fabs(partType) == 4) std::cout << " Mass of c-quark : " << mcParticles[i]->M() << " with mother : " << mcParticles[i]->motherType() << std::endl; --> Almost lways 1.5
	
    if(fabs(partType)<7 || fabs(partType)==24 || (fabs(partType)<=14 && fabs(partType)>=11) ){ //Considering only the semileptonic particles
      int motherType=mcParticles[i]->motherType(); 
      int grannyType=mcParticles[i]->grannyType();
      if(verbose_ > 5)cout<<"Mother type of particle : "<<motherType<<", and granny type : "<<grannyType<<endl;
	  
      if(partType == 6){      Top    = (TRootMCParticle*) mcParticles[i]; EventContent[0]++; if(verbose_>4) cout<<"*Particle found: Top"<<endl;    }      
      else if(partType == -6){TopBar = (TRootMCParticle*) mcParticles[i]; EventContent[0]++; if(verbose_>4) cout<<"*Particle found: AntiTop"<<endl;}	  	  

      else if(fabs(partType) == 5 && fabs(motherType) == 6){
	EventContent[1]++;
	if(partType == 5){      Bottom =    (TRootMCParticle*) mcParticles[i]; if(verbose_>4) cout<<"*Particle found: Bottom"<<endl;}
	else if(partType == -5){BottomBar = (TRootMCParticle*) mcParticles[i]; if(verbose_>4) cout<<"*Particle found: AntiBottom"<<endl;}
      }//End of bottom particle identification  
	  
      else if(fabs(partType) == 24 && fabs(motherType) == 6){//Check correct definition!!!
	EventContent[3]++;
	if(partType == 24){      WPlus =  (TRootMCParticle*) mcParticles[i]; if(verbose_>4) cout<<"*Particle found: WPlus"<<endl;}        
	else if(partType == -24){WMinus = (TRootMCParticle*) mcParticles[i]; if(verbose_>4) cout<<"*Particle found: WMinus"<<endl;}
      }//End of WBoson identification
	  
      else if(fabs(partType) <=4 && fabs(motherType) == 24 && fabs(grannyType) == 6){
	EventContent[2]++;
	if(partType > 0){     Light =    (TRootMCParticle*) mcParticles[i]; if(verbose_>4) cout<<"*Particle found: Light"<<endl;}
	else if(partType < 0){LightBar = (TRootMCParticle*) mcParticles[i]; if(verbose_>4) cout<<"*Particle found: AntiLight"<<endl;}
      }//End of light particle identification

      else if((fabs(partType) == 13 || fabs(partType) == 11 ) && fabs(motherType) == 24 && fabs(grannyType) == 6){
	EventContent[4]++;
	string leptonType="";
	if(fabs(partType) == 13){      if(verbose_>4) leptonType = "*Particle found: Muon";}
	else if(fabs(partType) == 11){ if(verbose_>4) leptonType = "*Particle found: Electron";}
	Lepton = (TRootMCParticle*) mcParticles[i]; if(verbose_ > 4) cout<<leptonType<<endl;
      }//End of lepton identification
	  
      else if((fabs(partType) == 14 || fabs(partType) == 12 ) && fabs(motherType) == 24 && fabs(grannyType) == 6){
	EventContent[4]++;
	string neutrinoType="";
	if(fabs(partType) == 14){      if(verbose_>4) neutrinoType = "*Particle found: Muon-neutrino";}
	else if(fabs(partType) == 12){ if(verbose_>4) neutrinoType = "*Particle found: Electron-neutrino";}
	NeutrinoMC = (TRootMCParticle*) mcParticles[i]; if(verbose_ > 4) cout<<neutrinoType<<endl;
      }//End of neutrino identification		
	  
    }//End of looking at semi-leptonic particles inside event ==> Semileptonic event is completely created now!	
  }//End of loop over mcParticles inside one particular event

  //---  Consider only events with correct event content (b b q q l vl)  ---//
  if(EventContent[0]==2 && EventContent[1]==2 && EventContent[2]==2 && EventContent[3]==2 && EventContent[4]==2){
    CorrectGenEvtContent = true;
    vector<TLorentzVector*> LHCOVector(6);
    vector<int> MadGraphId(6,4);
    vector<float> MGBtag(6,0.0);
    MGBtag[0] = 2.0; MGBtag[3] = 2.0;     //b-jets always on position 0 and 3!
	
    if(verbose_>3){
      cout << " Event with correct event content found "       << endl;
      cout << " Mass of bottom quark    : " << Bottom->M()     << endl;
      cout << " Mass of light quark     : " << Light->M()      << endl;
      cout << " Mass of LightBar quark  : " << LightBar->M()   << endl;
      cout << " Mass of BottomBar quark : " << BottomBar->M()  << endl;
      cout << " Mass of lepton          : " << Lepton->M()     << endl;
      cout << " Mass of neutrino        : " << NeutrinoMC->M() << endl;
    }

    //Save the lepton and neutrino as TLorentzVectors
    GenLepton = (TLorentzVector*) Lepton;
    GenNeutrino = (TLorentzVector*) NeutrinoMC;
	
    //Create the lhco file for pp > t t~:
    if(Lepton->type() == 13 || Lepton->type() == 11){ //Negative lepton, hence t~ > b~ W-, W- > e/mu- ve/vm
      LHCOVector[0] = (TLorentzVector*) Bottom;
      LHCOVector[1] = (TLorentzVector*) Light;
      LHCOVector[2] = (TLorentzVector*) LightBar;
      LHCOVector[3] = (TLorentzVector*) BottomBar;
      LHCOVector[4] = (TLorentzVector*) Lepton;
      LHCOVector[5] = (TLorentzVector*) NeutrinoMC;
      MadGraphId[5] = 6;                  //MadGraph Id of MET = 6
      if(Lepton->type() == 11){           //Looking at negative electron events (index 3 for LHCO file)
  	MadGraphId[4] = 1;                //MadGraph Id of e = 1
	if(genOutput_){
	  NumberNegativeElectrons++;
          leptonType = elMinus;      //Enum information	
	  LHCOEventOutput(3, GenOutFile[3], NumberNegativeElectrons,LHCOVector,MadGraphId, MGBtag);
	}
      }//Negative electron
      else if(Lepton->type() == 13){       //Looking at negative muon events (index 1 for LHCO file)
	MadGraphId[4] = 2; //MadGraphId of mu = 2
	if(genOutput_){
	  NumberNegativeMuons++;
          leptonType = muMinus;      //Enum information
	  LHCOEventOutput(1, GenOutFile[1], NumberNegativeMuons,LHCOVector,MadGraphId, MGBtag);
	}
      }//Negative muon

      //--- Store information needed for cos theta* ---//
      GenLeptonicW = (TLorentzVector*) WMinus;
      GenLeptonicTop = (TLorentzVector*) TopBar;
      //GenHadronicW = (TLorentzVector*) WPlus;
      //GenHadronicTop = (TLorentzVector*) Top;

    }//Negative lepton
    else if(Lepton->type() == -13 || Lepton->type() == -11){ //Positive lepton, hence t > b W+, W+ > e/mu+ ve/vm
      LHCOVector[0] = (TLorentzVector*) Bottom; 
      LHCOVector[1] = (TLorentzVector*) Lepton;
      LHCOVector[2] = (TLorentzVector*) NeutrinoMC;
      MadGraphId[2] = 6;          //MET always on position 2
      LHCOVector[3] = (TLorentzVector*) BottomBar;
      LHCOVector[4] = (TLorentzVector*) Light;
      LHCOVector[5] = (TLorentzVector*) LightBar;
      if(Lepton->type() == -11){            //Looking at positive electron events (index 2 for LHCO file)
	MadGraphId[1] = 1;                  //MadGraphId of electron = 1
	if(genOutput_){
	  NumberPositiveElectrons++;
          leptonType = elPlus;      //Enum information
	  LHCOEventOutput(2, GenOutFile[2], NumberPositiveElectrons,LHCOVector,MadGraphId, MGBtag);
	}
      }//Positive electron
      else if(Lepton->type() == -13){             //Looking at positive muon events (index 0 for LHCO file)
	MadGraphId[1] = 2;                        //MadGraphId of muon = 2
	if(genOutput_){
	  NumberPositiveMuons++;
          leptonType = muPlus;      //Enum information
	  LHCOEventOutput(0, GenOutFile[0], NumberPositiveMuons,LHCOVector,MadGraphId, MGBtag);
	}
      }//Positive muon
	  
      //--- Store information needed for cos theta* ---//
      GenLeptonicW = (TLorentzVector*) WPlus;
      GenLeptonicTop = (TLorentzVector*) Top;
      //GenHadronicW = (TLorentzVector*) WMinus;
      //GenHadronicTop = (TLorentzVector*) TopBar;

    }//Positive lepton
  }//Correct event content found
  else{
    CorrectGenEvtContent = false;
    WrongEvtCounter++;

    //Output in case wrong event content needs to be double-checked
    if(verbose_>4){
      cout << " Number of top quarks      : " << EventContent[0] << endl;
      cout << " Number of bottom quarks   : " << EventContent[1] << endl;
      cout << " Number of light quarks    : " << EventContent[2] << endl;
      cout << " Number of W-bosons        : " << EventContent[3] << endl;
      cout << " Number of lepton/neutrino : " << EventContent[4] << endl;
    }

    //---  Also create a .lhco file for wrong Gen events! ---//
    vector<TLorentzVector*> LHCOVectorWrongGen(6);
    vector<int> MadGraphIdWrongGen(6,4);
    vector<float> MGBtagWrongGen(6,0.0);
    MadGraphIdWrongGen[1] = 2; MadGraphIdWrongGen[2] = 6;    //Consider the event as a semi-mu (+) event
    MGBtagWrongGen[0] = 2.0; MGBtagWrongGen[3] = 2.0;     

    LHCOVectorWrongGen[0] = (TLorentzVector*) mcParticles[0];
    LHCOVectorWrongGen[1] = (TLorentzVector*) mcParticles[1];
    LHCOVectorWrongGen[2] = (TLorentzVector*) mcParticles[2];
    LHCOVectorWrongGen[3] = (TLorentzVector*) mcParticles[3];
    LHCOVectorWrongGen[4] = (TLorentzVector*) mcParticles[4];
    LHCOVectorWrongGen[5] = (TLorentzVector*) mcParticles[5];

    if(WrongEvtCounter <= 10000){
      LHCOEventOutput(2, WrongGenFile, WrongEvtCounter, LHCOVectorWrongGen, MadGraphIdWrongGen, MGBtagWrongGen);
      histo1D["WrongGen_TopMassLept"]->Fill( (*mcParticles[0]+*mcParticles[1]+*mcParticles[2]).M());
      histo1D["WrongGen_TopMassHadr"]->Fill( (*mcParticles[3]+*mcParticles[4]+*mcParticles[5]).M()); 	
    }
  }			    
}//End of class StoreGenInfo

void LHCOOutput::StoreRecoInfo(TLorentzVector* lepton, vector<TRootJet*> Jets, int bLeptIndex, int bHadrIndex, int light1Index, int light2Index, int decayChannel, float leptonCharge, vector<int> jetCombi){

  //--- Reconstruct neutrino partially ---//
  //---   (Neutrino Pt and M needed)   ---//
  //float NeutrinoPz=99.;                            //with this value it can be distinguished in plots!
  TLorentzVector Neutrino;
  float NeutrinoPx = -(*lepton+*Jets[bLeptIndex]+*Jets[bHadrIndex]+*Jets[light1Index]+*Jets[light2Index]).Px();
  float NeutrinoPy = -(*lepton+*Jets[bLeptIndex]+*Jets[bHadrIndex]+*Jets[light1Index]+*Jets[light2Index]).Py();
  Neutrino.SetPxPyPzE(NeutrinoPx, NeutrinoPy, 0.0, sqrt(NeutrinoPx*NeutrinoPx + NeutrinoPy*NeutrinoPy));
  Neutrino.SetPxPyPzE(NeutrinoPx, NeutrinoPy, 0.0, Neutrino.Pt());                        //Reset the Neutrino Energy to get the correct precision
  if(verbose_ > 3) std::cout << " Mass value for the neutrino : " << Neutrino.M() << " \n" << std::endl;
      
  //--- Initialize needed vectors ---//
  vector<TLorentzVector*> LHCORecoVector(6);
  vector<int> MadGraphRecoId(6,4);
  vector<float> MGRecoBtagId(6,0.0);
  if(Jets[bLeptIndex]->btag_combinedSecondaryVertexBJetTags() > 0 && Jets[bHadrIndex]->btag_combinedSecondaryVertexBJetTags() > 0)
    MGRecoBtagId[0] = 1; MGRecoBtagId[3] = 1;

  //--- Check whether the event is correctly reconstructed  ---// 
  //---  (jetCombi is initialized to 9999 for all dataSets) ---//
  bool jetCombiFound = false;
  bool EventCorrectlyMatched = false;
  if(jetCombi[0] != 9999 && jetCombi[1] != 9999 && jetCombi[2] != 9999 && jetCombi[3] != 9999){
    jetCombiFound = true;
    if( bLeptIndex == jetCombi[0] && bHadrIndex == jetCombi[1]    &&
       (light1Index == jetCombi[2] || light1Index == jetCombi[3]) &&
       (light2Index == jetCombi[3] || light2Index == jetCombi[3]) )
      EventCorrectlyMatched = true;
  }
	
  //---  Filling of LHCO files for reco events  ---//
  if(leptonCharge < 0.0 ){ //Negative lepton events
    LHCORecoVector[0] = Jets[bHadrIndex]; 
    LHCORecoVector[1] = Jets[light1Index];
    LHCORecoVector[2] = Jets[light2Index];
    LHCORecoVector[3] = Jets[bLeptIndex];
    LHCORecoVector[4] = lepton;
    LHCORecoVector[5] = &Neutrino;

    MadGraphRecoId[5] = 6;
    if(decayChannel == 1){//Negative electron
      MadGraphRecoId[4] = 1;
      NumberNegRecoEl++;  
      LHCOEventOutput(3,RecoOutFile[3], NumberNegRecoEl, LHCORecoVector, MadGraphRecoId, MGRecoBtagId);
    }
    if(decayChannel == 0){//Negative muon
      MadGraphRecoId[4] = 2;
      NumberNegRecoMu++;
      LHCOEventOutput(1, RecoOutFile[1], NumberNegRecoMu, LHCORecoVector, MadGraphRecoId, MGRecoBtagId);
    }
  }//End of negative lepton
  else if(leptonCharge > 0.0 ){ //Positive lepton events
    LHCORecoVector[0] = Jets[bLeptIndex];
    LHCORecoVector[1] = lepton;
    LHCORecoVector[2] = &Neutrino;
    LHCORecoVector[3] = Jets[bHadrIndex];
    LHCORecoVector[4] = Jets[light1Index];
    LHCORecoVector[5] = Jets[light2Index];

    MadGraphRecoId[2] = 6;
    if(decayChannel == 1){//Positive electron
      MadGraphRecoId[1] = 1;
      NumberPosRecoEl++;
      LHCOEventOutput(2,RecoOutFile[2], NumberPosRecoEl, LHCORecoVector, MadGraphRecoId, MGRecoBtagId);	 
    }
    if(decayChannel == 0){//Positive muon
      MadGraphRecoId[1] = 2;
      NumberPosRecoMu++;
      LHCOEventOutput(0, RecoOutFile[0], NumberPosRecoMu, LHCORecoVector, MadGraphRecoId, MGRecoBtagId);
      if(EventCorrectlyMatched == true && jetCombiFound == true){
        NrPosRecoMuCorrect++; 
        LHCOEventOutput(0, CorrectRecoMuPosFile, NrPosRecoMuCorrect, LHCORecoVector, MadGraphRecoId, MGRecoBtagId);
        histo1D["CorrectReco_TopMassLept"]->Fill( (*LHCORecoVector[0]+*LHCORecoVector[1]+*LHCORecoVector[2]).M());
        histo1D["CorrectReco_TopMassHadr"]->Fill( (*LHCORecoVector[3]+*LHCORecoVector[4]+*LHCORecoVector[5]).M());
      }
      else if(EventCorrectlyMatched == false && jetCombiFound == true){ 
        NrPosRecoMuWrong++;   
        LHCOEventOutput(0, WrongRecoMuPosFile, NrPosRecoMuWrong, LHCORecoVector, MadGraphRecoId, MGRecoBtagId);
        histo1D["WrongReco_TopMassLept"]->Fill( (*LHCORecoVector[0]+*LHCORecoVector[1]+*LHCORecoVector[2]).M());
        histo1D["WrongReco_TopMassHadr"]->Fill( (*LHCORecoVector[3]+*LHCORecoVector[4]+*LHCORecoVector[5]).M());
      }
    }
  }//End of positive lepton

}//End of class StoreRecoInfo 

void LHCOOutput::LHCOEventOutput(int LHCOIndex, ostream &outputFile, unsigned int EventNumber, std::vector<TLorentzVector*> vector, std::vector<int> MGId, std::vector<float> MGBtag){

  if(LHCOIndex == 0 || LHCOIndex == 2)
    LeptonCharge =1;
  else if(LHCOIndex == 1 || LHCOIndex == 3)
    LeptonCharge = -1;

  if(EventNumber == 1){
    outputFile << "#</MGPGSCard> " << endl;
    outputFile << "  #  typ      eta      phi       pt   jmas  ntrk  btag   had/em  dummy  dummy " << endl;
  }

  outputFile << " 0             " << EventNumber << "        6 " << endl;  //Start of a new event

  for(int ii = 0; ii < 6; ii++){
    outputFile << "  " << setprecision(1) << ii+1;
    outputFile << "    " << setprecision(1) << MGId[ii];
    cout.setf(ios::fixed,ios::floatfield);  //Add zero to obtain the asked number of digits
    cout.precision(4);
    outputFile << "  " << fixed << showpoint << setprecision(4) << vector[ii]->Eta();
    if(vector[ii]->Phi() < -3.14) outputFile << "  " << setprecision(4) << -1*(vector[ii]->Phi());
    else                          outputFile << "  " << setprecision(4) << vector[ii]->Phi();
    outputFile << "  " << setprecision(4) << vector[ii]->Pt();
    cout.setf(ios::fixed,ios::floatfield);  //Add zero to obtain the asked number of digits
    cout.precision(4);
    if(vector[ii]->M() > 0.0) outputFile << "  " << fixed << showpoint << setprecision(3) << vector[ii]->M();
    else                      outputFile << "  " << fixed << showpoint << setprecision(3) << " 0.00";
    if(MGId[ii] == 2 || MGId[ii] == 1) outputFile << "   " << setprecision(3) << LeptonCharge << "   ";       //MGId equal to 1 or 2 corresponds to muon or electron
    else                               outputFile << "    0.00";
    outputFile << " " << MGBtag[ii] << "     0.00  0.00  0.00" << endl;
  }
}

void LHCOOutput::WriteLHCOPlots(TFile* outfile){
  //--- Use this function to create ChiSq histograms ---//
  outfile->cd();
  std::cout << " Inside WriteLHCOPlots function of LHCOOutput class ! " << std::endl;
  std::cout << " Histograms will be filled in file : " << outfile->GetName() << " ************************************" << std::endl;

  TDirectory* th1dir = outfile->mkdir("1D_histograms_LHCOOutput");
  th1dir->cd();
  for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++){
    TH1F *temp = it->second;
    int N = temp->GetNbinsX();
    temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
    temp->SetBinContent(N+1,0);
    temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
    temp->Write();
  }
  TDirectory* th2dir = outfile->mkdir("2D_histograms_LHCOOutput");
  th2dir->cd();
  for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++){    
    TH2F *temp = it->second;
    temp->Write();
  }
  outfile->cd(); 

}
