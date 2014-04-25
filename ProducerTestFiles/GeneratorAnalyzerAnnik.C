#include <iomanip>
#include "../interface/TRootGenEvent.h"
#include "../interface/TRootEvent.h"
#include "../interface/TRootRun.h"
#include "../interface/TRootParticle.h"
#include "../interface/TRootMCParticle.h"
#include "../interface/TRootVertex.h"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TBranch.h>
#include <TTree.h>

#include "TClonesArray.h"
#include<fstream>
#include <sstream>

//Needed for TCanvasCreator
//#include "../../TopTreeAnalysisBase/Tools/interface/PlottingTools.h"

using namespace TopTree;

class EventOutput{
	
	int LeptonCharge;
	public:
		void LHCOOutput(int LHCOIndex, ofstream &outputFile, unsigned int EventNumber, vector<TRootMCParticle*> vector, vector<int> MGId);

};

void EventOutput::LHCOOutput(int LHCOIndex, ofstream &outputFile, unsigned int EventNumber, vector<TRootMCParticle*> vector, vector<int> MGId){

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
		if(ii == 1) outputFile << "   " << setprecision(3) << LeptonCharge;
		else outputFile << "    0.00";
		outputFile << "  0.00     0.00  0.00  0.00" << endl;		
	}
}

int main(){
	clock_t start = clock();

        int verbosity                 = 1;
	//0 muet
	//1 Main Info
	//2
	//3 
	//4 Info for each event
	//5 Debug

	bool doMC                     = true;
	bool doGenEvent               = true;

	TFile* f = TFile::Open("TopTree_Skimmed_100.root");
	//user/aolbrech/AnomalousCouplings/GitTopTree/TopBrussels/TopTreeProducer/test/22042013_211038_TOPTREE_100_1_YYM.root");
	//dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/echabert/TopTree.root");
	TTree* runTree = (TTree*) f->Get("runTree");
	TTree* eventTree = (TTree*) f->Get("eventTree");

	TBranch* run_br = (TBranch *) runTree->GetBranch("runInfos");
	TRootRun* runInfos = 0;
	run_br->SetAddress(&runInfos);
	
	TBranch* event_br = (TBranch *) eventTree->GetBranch("Event");
	TRootEvent* event = 0;
	event_br->SetAddress(&event);
	
	//Declartion of Branches and TClonesArray
	TBranch* mcParticles_br;
	TClonesArray* mcParticles;
	TBranch* genEvents_br;
	TClonesArray* genEvents;

	//Access class to write output to correct .lhco file
	EventOutput eventOutput;
	ofstream outFile[4];
	outFile[0].open("TTbarLHCO_PositiveMuon.lhco");
	outFile[1].open("TTbarLHCO_NegativeMuon.lhco");
	outFile[2].open("TTbarLHCO_PositiveElectron.lhco");
	outFile[3].open("TTbarLHCO_NegativeElectron.lhco");
	
	if(doMC)
	{
		mcParticles_br = (TBranch *) eventTree->GetBranch("MCParticles");
		mcParticles = new TClonesArray("TopTree::TRootMCParticle", 0);
		mcParticles_br->SetAddress(&mcParticles);
	}
        
	if(doGenEvent)
	{
		genEvents_br = (TBranch *) eventTree->GetBranch("GenEvent");
		genEvents = new TClonesArray("TopTree::TRootGenEvent", 0);
		genEvents_br->SetAddress(&genEvents);
	}

	//Counting number of events for LHCO output
	unsigned int NumberCorrectEvents = 0; //Counts the number of semi-leptonic events
	unsigned int NumberNegativeElectrons = 0;
	unsigned int NumberNegativeMuons = 0;
	unsigned int NumberPositiveElectrons = 0;
	unsigned int NumberPositiveMuons = 0;
	int EventContent[5]; //0:top; 1:b; 2: u,c,d,s; 3:W; 4:mu + neutrino

	//Cos Theta information
	TLorentzVector *sTop, *WLeptTRF, *leptonWRF;
	float standardCosTheta = 0;
	TH1F h_StandardCosTheta("StCosTheta","StCosTheta",200,-1,1);

	//Loop over all events and access particle information
        unsigned int nEvents = (int)eventTree->GetEntries();
	if(verbosity>0) cout << " Total number of events: " << nEvents << endl;
	for(unsigned int ievt=0; ievt<nEvents; ievt++){
	        if(ievt%5000 == 0)
		  std::cout<<"Processing the "<<ievt<<"th event" <<flush<<"\r";

		bool FalseEventContent = false;
                eventTree->GetEvent(ievt);
		if(verbosity > 3) cout<<"\nLooking at event : "<<ievt<<endl;

	        //Access mcParticles information:
		if(doMC){
		  if(verbosity >4) cout<<"Access to mcParticles"<<endl;
		  if(mcParticles->GetEntriesFast()>0){
			TRootMCParticle *Top, *TopBar, *Bottom, *BottomBar, *Lepton, *Neutrino, *WPlus, *WMinus,*Light, *LightBar;
			for(int ll = 0;ll<5;ll++){EventContent[ll]=0;}

		  	for(int i = 0; i < mcParticles->GetEntriesFast();i++){
			  if(((TRootMCParticle*) mcParticles->At(i))->status() != 3) continue; //Only look at particles with status = 3!

			  int partType = ((TRootMCParticle*) mcParticles->At(i))->type();if(verbosity > 4) cout<<"-->Type of mcParticle : "<<partType<<endl;

			  if(fabs(partType) < 7 || fabs(partType) == 24 || (fabs(partType) <= 14 && fabs(partType) >= 11) ){ //Considering only the semileptonic particles
				int motherType=((TRootMCParticle*) mcParticles->At(i))->motherType(); int grannyType=((TRootMCParticle*) mcParticles->At(i))->grannyType();
				if(verbosity > 5) cout<<"Mother type of particle : "<< motherType<<", and granny type : "<<grannyType<<endl;

				if(partType == 6){      Top   =(TRootMCParticle*) mcParticles->At(i); EventContent[0]++; if(verbosity>4) cout<<"*Particle found: Top"<<endl;}
				else if(partType == -6){TopBar=(TRootMCParticle*) mcParticles->At(i); EventContent[0]++; if(verbosity>4) cout<<"*Particle found: AntiTop"<<endl;}

				else if(fabs(partType) == 5 && fabs(motherType) == 6){
					EventContent[1]++;
					if(partType == 5){      Bottom =    (TRootMCParticle*) mcParticles->At(i); if(verbosity>4) cout<<"*Particle found: Bottom"<<endl;}
					else if(partType == -5){BottomBar = (TRootMCParticle*) mcParticles->At(i); if(verbosity>4) cout<<"*Particle found: AntiBottom"<<endl;}
				}//End of bottom particle identification

				else if(fabs(partType) == 24 && fabs(motherType) == 6){//Check correct definition!!!
					EventContent[3]++;
					if(partType == 24){      WPlus =  (TRootMCParticle*) mcParticles->At(i); if(verbosity>4) cout<<"*Particle found: WPlus"<<endl;}        
					else if(partType == -24){WMinus = (TRootMCParticle*) mcParticles->At(i); if(verbosity>4) cout<<"*Particle found: WMinus"<<endl;}
				}//End of WBoson identification

				else if(fabs(partType) <=4 && fabs(motherType) == 24 && fabs(grannyType) == 6){
					EventContent[2]++;
					if(partType > 0){     Light =    (TRootMCParticle*) mcParticles->At(i); if(verbosity>4) cout<<"*Particle found: Light"<<endl;}
					else if(partType < 0){LightBar = (TRootMCParticle*) mcParticles->At(i); if(verbosity>4) cout<<"*Particle found: AntiLight"<<endl;}
				}//End of light particle identification

				else if((fabs(partType) == 13 || fabs(partType) == 11 ) && fabs(motherType) == 24 && fabs(grannyType) == 6){
					EventContent[4]++;
					string leptonType="";
					if(fabs(partType) == 13){      if(verbosity>4) leptonType = "*Particle found: Muon";}
					else if(fabs(partType) == 11){ if(verbosity>4) leptonType = "*Particle found: Electron";}
					Lepton = (TRootMCParticle*) mcParticles->At(i); if(verbosity > 4) cout<<leptonType<<endl;
				}//End of lepton identification

				else if((fabs(partType) == 14 || fabs(partType) == 12 ) && fabs(motherType) == 24 && fabs(grannyType) == 6){
					EventContent[4]++;
					string neutrinoType="";
					if(fabs(partType) == 14){      if(verbosity>4) neutrinoType = "*Particle found: Muon-neutrino";}
					else if(fabs(partType) == 12){ if(verbosity>4) neutrinoType = "*Particle found: Electron-neutrino";}
					Neutrino = (TRootMCParticle*) mcParticles->At(i); if(verbosity > 4) cout<<neutrinoType<<endl;
				}//End of neutrino identification		
		
			  }//End of looking at semi-leptonic particles inside event ==> Semileptonic event is completely created now!

			}//End of loop over mcParticles inside one particular event

			//////////////////////////////////////////////////////////////////////
			//  Consider only events with correct event content (b b q q l vl)  //
			//////////////////////////////////////////////////////////////////////
			if(EventContent[0]==2 && EventContent[1]==2 && EventContent[2]==2 && EventContent[3]==2 && EventContent[4]==2){
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
					cout << " Mass of neutrino : " << Neutrino->M() << endl <<endl;
			 	}

				//Create the lhco file for pp > t t~:
				if(Lepton->type() == 13 || Lepton->type() == 11){ //Negative lepton, hence t~ > b~ W-, W- > e/mu- ve/vm
					LHCOVector[0] = Bottom;
					LHCOVector[1] = Light;
					LHCOVector[2] = LightBar;
					LHCOVector[3] = BottomBar;
					LHCOVector[4] = Lepton;
					LHCOVector[5] = Neutrino;
					if(Lepton->type() == 11){           //Looking at negative electron events (index 3 for LHCO file)
						MadGraphId[4] = 1; //MadGraphId of e = 1
						MadGraphId[5] = 6; 
						NumberNegativeElectrons++;
						eventOutput.LHCOOutput(3, outFile[3], NumberNegativeElectrons,LHCOVector,MadGraphId);
					}//Negative electron
					else if(Lepton->type() == 13){       //Looking at negative muon events (index 1 for LHCO file)
						MadGraphId[4] = 2; //MadGraphId of mu = 2
						MadGraphId[5] = 6; 
						NumberNegativeMuons++;
						eventOutput.LHCOOutput(1, outFile[1], NumberNegativeMuons,LHCOVector,MadGraphId);
					}//Negative muon

					if(verbosity>3){
						cout<<" WMinus information : "<<WMinus->Px()<< ", "<<WMinus->Py()<<", "<< WMinus->Pz()<<", "<<WMinus->E()<< endl;
						cout<<" TopBar information : "<<TopBar->Px()<< ", "<<TopBar->Py()<<", "<< TopBar->Pz()<<", "<<TopBar->E()<< endl;
					}
					WLeptTRF = (TLorentzVector*) WMinus;
					sTop = (TLorentzVector*) TopBar;				
					//WLeptTRF->SetPxPyPzE((double)WMinus->Px(), (double)WMinus->Py(), (double)WMinus->Pz(), (double)WMinus->E());
					//sTop->SetPxPyPzE((double)TopBar->Px(), (double)TopBar->Py(), (double)TopBar->Pz(), (double)TopBar->E());
					if(verbosity>3){
						cout<<" WLeptTRF information : "<<WLeptTRF->Px()<<", "<<WLeptTRF->Py()<<", "<<WLeptTRF->Pz()<<", "<<WLeptTRF->E()<<endl;
						cout<<" sTop information : "<<sTop->Px()<<", "<<sTop->Py()<<", "<<sTop->Pz()<<", "<<sTop->E()<<endl;
					}
				}//Negative lepton
				else if(Lepton->type() == -13 || Lepton->type() == -11){ //Positive lepton, hence t > b W+, W+ > e/mu+ ve/vm
					LHCOVector[0] = Bottom; 
					LHCOVector[1] = Lepton;
					LHCOVector[2] = Neutrino;
					LHCOVector[3] = BottomBar;
					LHCOVector[4] = Light;
					LHCOVector[5] = LightBar;
					if(Lepton->type() == -11){            //Looking at positive electron events (index 2 for LHCO file)
						MadGraphId[1] = 1; //MadGraphId of electron = 1
						MadGraphId[2] = 6; 
						NumberPositiveElectrons++;
						eventOutput.LHCOOutput(2, outFile[2], NumberPositiveElectrons,LHCOVector,MadGraphId);
					}//Positive electron
					else if(Lepton->type() == -13){             //Looking at positive muon events (index 0 for LHCO file)
						MadGraphId[1] = 2; //MadGraphId of muon = 2
						MadGraphId[2] = 6; 
						NumberPositiveMuons++;
						eventOutput.LHCOOutput(0, outFile[0], NumberPositiveMuons,LHCOVector,MadGraphId);
					}//Positive muon

					if(verbosity>3){
						cout << " WPlus information : "<<WPlus->Px()<< ", "<<WPlus->Py()<<", "<< WPlus->Pz()<<", "<<WPlus->E()<< endl;
						cout << " Top information : "<<Top->Px()<< ", "<<Top->Py()<<", "<< Top->Pz()<<", "<<Top->E()<< endl;
					}
					WLeptTRF = (TLorentzVector*) WPlus;
					sTop = (TLorentzVector*) Top;			
				}//Positive lepton
				
				//////////////////////////////////////
				//  Look at cos theta distribution  //
				//////////////////////////////////////

				//-----    Applying boost on muon and W    -----//
				leptonWRF = Lepton;
				leptonWRF->Boost(-WLeptTRF->BoostVector());
				WLeptTRF->Boost(-sTop->BoostVector());
				if(verbosity>3){
					cout<<" leptonWRF information : "<<leptonWRF->Px()<<", "<<leptonWRF->Py()<<", "<<leptonWRF->Pz()<<", "<<leptonWRF->E()<<endl;
				}

				//-----   Calculating cos theta:   -----
			  	standardCosTheta = ((WLeptTRF->Vect()).Dot(leptonWRF->Vect()))/(((WLeptTRF->Vect()).Mag())*((leptonWRF->Vect()).Mag()));
				if(verbosity>4) cout << " cos theta : " << standardCosTheta << endl;
				h_StandardCosTheta.Fill(standardCosTheta);
									
			}//Correct event content found
			else{
				FalseEventContent = true;
				if(verbosity>4){
					cout << " Number of top quarks      : " << EventContent[0] << endl;
					cout << " Number of bottom quarks   : " << EventContent[1] << endl;
					cout << " Number of light quarks    : " << EventContent[2] << endl;
					cout << " Number of W-bosons        : " << EventContent[3] << endl;
					cout << " Number of lepton/neutrino : " << EventContent[4] << endl;
				}
			}			
	          }//Only looking at events with mcParticles!
		  else if(verbosity>0) cout<<"No access to mcParticles in this entry"<<endl;

		}//End of doMC

		//access to GenEvent
		if(doGenEvent){
		  if(verbosity>4) cout<<"Access to GenEvent"<<endl;
		  if(genEvents->GetEntriesFast()>0){  
		    TRootGenEvent* genEvent = (TRootGenEvent*) genEvents->At(0);
		    int semiLeptonicEvent = genEvent->isSemiLeptonic();
		    if(FalseEventContent == true && genEvent->semiLeptonicChannel() != 3){
			cout << " Problem event found!!! Should be semi-leptonic but lepton != tau!! " << endl;
			cout << " Semileptonic channel : " << genEvent->semiLeptonicChannel() << endl;
		    }
		  }
		  else if(verbosity>0) cout<<" No access to GenEvent in this entry"<<endl;
		}

	}// end of loop over evts
	if(verbosity>0) cout << "---> Number of events with correct semileptonic event content on generator level: " << NumberCorrectEvents << " (semiMuon, semiElec) : ( " << NumberPositiveMuons+NumberNegativeMuons << " , " << NumberPositiveElectrons+NumberNegativeElectrons << " ) " << endl;

	for(int ii = 0; ii<4; ii++)
		outFile[ii].close();	


	//delete mcParticles;
        delete mcParticles_br;

	if(verbosity>1) cout<<"Writting histograms in the root-file ... "<<endl;      
        TFile* fout = new TFile("GeneratorOutput.root","RECREATE");
        fout->cd();

	h_StandardCosTheta.Write();

        fout->Close();

        if(verbosity>1) cout<<"End of the Macro"<<endl;
        if(verbosity>0) cout << "It took " << ((double)clock() - start) / CLOCKS_PER_SEC << "s to run the program" << endl;     

        return(0);
}
