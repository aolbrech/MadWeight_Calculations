#include <iomanip>
#include "../interface/LHCOOutput.h"

#include<fstream>
#include <sstream>

LHCOOutput::LHCOOutput(){
    
    //Constructor: Define variables which have to be initialized in the beginning!
    NumberNegativeElectrons = 0;
    NumberNegativeMuons = 0;
    NumberPositiveElectrons = 0;
    NumberPositiveMuons = 0;
    CorrectGenEvtContent = false;
    //std::string leptonTypeString[4] = {"muPlus","muMinus","elPlus","elMinus"};

    GenOutFile[0].open("MadWeightInput/AnalyzerOutput/TTbarLHCO_PositiveMuon.lhco");
    GenOutFile[1].open("MadWeightInput/AnalyzerOutput/TTbarLHCO_NegativeMuon.lhco");
    GenOutFile[2].open("MadWeightInput/AnalyzerOutput/TTbarLHCO_PositiveElectron.lhco");
    GenOutFile[3].open("MadWeightInput/AnalyzerOutput/TTbarLHCO_NegativeElectron.lhco");

}

LHCOOutput::~LHCOOutput(){

    //Close the output files:
    for(int ii = 0; ii < 4; ii++) GenOutFile[ii].close();	

}

void LHCOOutput::StoreGenInfo(vector<TRootMCParticle*> mcParticles, bool GenLHCOOutput, int verbosity){

    //Initialize EventContent counter used for selecting correct particle content!
    int EventContent[5]; //0:top; 1:b; 2: u,c,d,s; 3:W; 4:mu + neutrino
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
      if(EventContent[0]==2 && EventContent[1]==2 && EventContent[2]==2 && EventContent[3]==2 && EventContent[4]==2){
	CorrectGenEvtContent = true;
	vector<TRootMCParticle*> LHCOVector(6);
	vector<int> MadGraphId(6,4);
        vector<float> MGBtag(6,0.0);
        MGBtag[0] = 2.0; MGBtag[3] = 2.0;     //b-jets always on position 0 and 3!
	
	if(verbosity>3){
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
	  LHCOVector[0] = Bottom;
	  LHCOVector[1] = Light;
	  LHCOVector[2] = LightBar;
	  LHCOVector[3] = BottomBar;
	  LHCOVector[4] = Lepton;
	  LHCOVector[5] = NeutrinoMC;
	  MadGraphId[5] = 6;                  //MadGraph Id of MET = 6
	  if(Lepton->type() == 11){           //Looking at negative electron events (index 3 for LHCO file)
	    MadGraphId[4] = 1;                //MadGraph Id of e = 1
	    if(GenLHCOOutput == true){
	      NumberNegativeElectrons++;
              leptonType = elMinus;      //Enum information	
	      LHCOEventOutput(3, GenOutFile[3], NumberNegativeElectrons,LHCOVector,MadGraphId, MGBtag);
	    }
	  }//Negative electron
	  else if(Lepton->type() == 13){       //Looking at negative muon events (index 1 for LHCO file)
	    MadGraphId[4] = 2; //MadGraphId of mu = 2
	    if(GenLHCOOutput == true){
	      NumberNegativeMuons++;
              leptonType = muMinus;      //Enum information
	      LHCOEventOutput(1, GenOutFile[1], NumberNegativeMuons,LHCOVector,MadGraphId, MGBtag);
	    }
	  }//Negative muon

          //Store information needed for cos theta*
	  GenLeptonicW = (TLorentzVector*) WMinus;
	  GenLeptonicTop = (TLorentzVector*) TopBar;
	  //GenHadronicW = (TLorentzVector*) WPlus;
	  //GenHadronicTop = (TLorentzVector*) Top;

			
	}//Negative lepton
	else if(Lepton->type() == -13 || Lepton->type() == -11){ //Positive lepton, hence t > b W+, W+ > e/mu+ ve/vm
	  LHCOVector[0] = Bottom; 
	  LHCOVector[1] = Lepton;
	  LHCOVector[2] = NeutrinoMC;
	  MadGraphId[2] = 6;          //MET always on position 2
	  LHCOVector[3] = BottomBar;
	  LHCOVector[4] = Light;
	  LHCOVector[5] = LightBar;
	  if(Lepton->type() == -11){            //Looking at positive electron events (index 2 for LHCO file)
	    MadGraphId[1] = 1;                  //MadGraphId of electron = 1
	    if(GenLHCOOutput == true){
	      NumberPositiveElectrons++;
              leptonType = elPlus;      //Enum information
	      LHCOEventOutput(2, GenOutFile[2], NumberPositiveElectrons,LHCOVector,MadGraphId, MGBtag);
	    }
	  }//Positive electron
	  else if(Lepton->type() == -13){             //Looking at positive muon events (index 0 for LHCO file)
	    MadGraphId[1] = 2;                        //MadGraphId of muon = 2
	    if(GenLHCOOutput == true){
	      NumberPositiveMuons++;
              leptonType = muPlus;      //Enum information
	      LHCOEventOutput(0, GenOutFile[0], NumberPositiveMuons,LHCOVector,MadGraphId, MGBtag);
	    }
	  }//Positive muon
	  
          //Store information needed for cos theta*
	  GenLeptonicW = (TLorentzVector*) WPlus;
	  GenLeptonicTop = (TLorentzVector*) Top;
	  //GenHadronicW = (TLorentzVector*) WMinus;
	  //GenHadronicTop = (TLorentzVector*) TopBar;

	}//Positive lepton

	
      }//Correct event content found
      else{
	CorrectGenEvtContent = false;

        //Output in case wrong event content needs to be double-checked
	if(verbosity>4){
	  cout << " Number of top quarks      : " << EventContent[0] << endl;
	  cout << " Number of bottom quarks   : " << EventContent[1] << endl;
	  cout << " Number of light quarks    : " << EventContent[2] << endl;
	  cout << " Number of W-bosons        : " << EventContent[3] << endl;
	  cout << " Number of lepton/neutrino : " << EventContent[4] << endl;
	}
      }			    
}//End of class StoreGenInfo

void LHCOOutput::LHCOEventOutput(int LHCOIndex, ostream &outputFile, unsigned int EventNumber, std::vector<TRootMCParticle*> vector, std::vector<int> MGId, std::vector<float> MGBtag){

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
    else outputFile << "  " << setprecision(4) << vector[ii]->Phi();
    outputFile << "  " << setprecision(4) << vector[ii]->Pt();
    cout.setf(ios::fixed,ios::floatfield);  //Add zero to obtain the asked number of digits
    cout.precision(4);
    if(vector[ii]->M() > 0.0) outputFile << "  " << fixed << showpoint << setprecision(3) << vector[ii]->M();
    //if(vector[ii]->M() == 1.5) std::cout << " Particle Id of particle is : " << vector[ii]->type() << std::endl; --> Always c-quark (=4)
    else outputFile << "  " << fixed << showpoint << setprecision(3) << " 0.00";
    if(MGId[ii] == 2 || MGId[ii] == 1) outputFile << "   " << setprecision(3) << LeptonCharge << "   ";
    else outputFile << "    0.00";
    outputFile << " " << MGBtag[ii] << "     0.00  0.00  0.00" << endl;
  }
}

void LHCOOutput::LHCOEventRecoOutput(int LHCOIndex, ostream &outputFile, unsigned int EventNumber, std::vector<TLorentzVector*> vector, std::vector<int> MGId, std::vector<float> MGBtag){

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
    else outputFile << "  " << setprecision(4) << vector[ii]->Phi();
    outputFile << "  " << setprecision(4) << vector[ii]->Pt();
    cout.setf(ios::fixed,ios::floatfield);  //Add zero to obtain the asked number of digits
    cout.precision(4);
    if(vector[ii]->M() > 0.0) outputFile << "  " << fixed << showpoint << setprecision(3) << vector[ii]->M();
    else outputFile << "  " << fixed << showpoint << setprecision(3) << " 0.00";
    if(ii == 1) outputFile << "   " << setprecision(3) << LeptonCharge << "   ";
    else outputFile << "    0.00";
    outputFile << " " << MGBtag[ii] << "     0.00  0.00  0.00" << endl;
  }
}

