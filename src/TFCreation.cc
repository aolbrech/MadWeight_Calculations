#include "../interface/TFCreation.h"

void TFCreation::InitializeVariables(){
	histo2D["Epart_vs_Enonbjet"]                     =new TH2F("Epart_vs_Enonbjet","Energy of first light quark ( gen level (x-axis) vs reco level (y-axis) )", 150, 0, 350, 150, 0, 350);
	histo2D["Epart_vs_DiffEpartEnonbjet"]            =new TH2F("Epart_vs_DiffEpartEnonbjet","Energy of first light quark versus #Delta(E_{gen} - E_{reco})", 10, 20, 250, 160, -80, 80);
	histo2D["Epart_vs_Ebjet"]                        =new TH2F("Epart_vs_Ebjet","Epart_vs_Ebjet", 150, 0, 350, 150, 0, 350);
	histo2D["Epart_vs_DiffEpartEbjet"]               =new TH2F("Epart_vs_DiffEpartEbjet","Epart_vs_Epart-Ebjet", 10, 20, 250, 160, -80, 80);
	histo2D["EgenEl_vs_ErecEl"]                      =new TH2F("EgenEl_vs_ErecEl","EgenEl_vs_ErecEl", 100, 38, 42, 150, 20, 400);
	histo2D["EgenEl_vs_DiffEgenElErecEl"]            =new TH2F("EgenEl_vs_DiffEgenElErecEl","EgenEl_vs_EgenEl-ErecEl", 10, 38, 42, 75, -140, 25);
	histo2D["InvPtgenMu_vs_InvPtrecMu"]              =new TH2F("InvPtgenMu_vs_InvPtrecMu","InvPtgenMu_vs_InvPtrecMu", 10, 0.015, 0.05, 80, 0, 0.04);
	histo2D["InvPtgenMu_vs_DiffInvPtgenMuInvPtrecMu"]=new TH2F("InvPtgenMu_vs_DiffInvPtgenMuInvPtrecMu","InvPtgenMu_vs_InvPtgenMu-InvPtrecMu", 10, 0.022, 0.04, 80, -0.015, 0.04);
	histo2D["PtgenMu_vs_PtrecMu"]                    =new TH2F("PtgenMu_vs_PtrecMu","PtgenMu_vs_PtrecMu", 50, 0, 60, 150, 0, 250);
	histo2D["PtgenMu_vs_DiffPtgenMuPtrecMu"]         =new TH2F("PtgenMu_vs_DiffPtgenMuPtrecMu","PtgenMu_vs_PtgenMu-PtrecMu",10, 0, 60, 100, -120, 30);
	
	histo2D["Thpart_vs_Thnonbjet"]          =new TH2F("Thpart_vs_Thnonbjet","Thpart_vs_Thnonbjet", 60, 0, 3.15, 60, 0, 3.15);
	histo2D["Thpart_vs_DiffThpartThnonbjet"]=new TH2F("Thpart_vs_DiffThpartThnonbjet","Thpart_vs_Thpart-Thnonbjet", 10, 0.1, 3.1, 100, -0.15, 0.15);
	histo2D["Thpart_vs_Thbjet"]             =new TH2F("Thpart_vs_Thbjet","Thpart_vs_Thbjet", 60, 0, 3.15, 60, 0, 3.15);
	histo2D["Thpart_vs_DiffThpartThbjet"]   =new TH2F("Thpart_vs_DiffThpartThbjet","Thpart_vs_Thpart-Thbjet", 10, 0.1, 3.1, 100, -0.15, 0.15);
	histo2D["ThgenEl_vs_ThrecEl"]           =new TH2F("ThgenEl_vs_ThrecEl","ThgenEl_vs_ThrecEl", 60, 0, 3.15, 60, 0, 3.15);
	histo2D["ThgenEl_vs_DiffThgenElThrecEl"]=new TH2F("ThgenEl_vs_DiffThgenElThrecEl","ThgenEl_vs_ThgenEl-ThrecEl", 10, 0, 3.15, 100, -2, 2);
	histo2D["ThgenMu_vs_ThrecMu"]           =new TH2F("ThgenMu_vs_ThrecMu","ThgenMu_vs_ThrecMu", 60, 0, 3.15, 60, 0, 3.15);
	histo2D["ThgenMu_vs_DiffThgenMuThrecMu"]=new TH2F("ThgenMu_vs_DiffThgenMuThrecMu","ThgenMu_vs_ThgenMu-ThrecMu", 10, 0, 3.15, 100, -2, 2);
	histo2D["Epart_vs_Thnonbjet"]           =new TH2F("Epart_vs_Thnonbjet","Epart_vs_Thnonbjet", 120, 0, 300, 60, 0, 3.15);
	histo2D["Epart_vs_DiffThpartThnonbjet"] =new TH2F("Epart_vs_DiffThpartThnonbjet","Epart_vs_Thpart-Thnonbjet", 10, 20, 250, 150, -0.15, 0.15);
	histo2D["Epart_vs_Thbjet"]              =new TH2F("Epart_vs_Thbjet","Epart_vs_Thbjet", 120, 0, 300, 60, 0, 3.15);
	histo2D["Epart_vs_DiffThpartThbjet"]    =new TH2F("Epart_vs_DiffThpartThbjet","Epart_vs_Thpart-Thbjet", 10, 20, 250, 150, -0.15, 0.15);
	histo2D["EgenEl_vs_ThrecEl"]            =new TH2F("EgenEl_vs_ThrecEl","EgenEl_vs_ThrecEl", 100, 38, 42, 60, 0, 3.15);
	histo2D["EgenEl_vs_DiffThgenElThrecEl"] =new TH2F("EgenEl_vs_DiffThgenElThrecEl","EgenEl_vs_ThgenEl-ThrecEl", 10, 38, 42, 100, -2, 2);
	histo2D["InvPtgenMu_vs_ThrecMu"]        =new TH2F("InvPtgenMu_vs_ThrecMu","InvPtgenMu_vs_ThrecMu", 50, 0.015, 0.045, 60, 0, 3.15);
	histo2D["InvPtgenMu_vs_DiffThgenMuThrecMu"]=new TH2F("InvPtgenMu_vs_DiffThgenMuThrecMu","InvPtgenMu_vs_ThgenMu-ThrecMu", 10, 0.022, 0.04, 100, -2, 2);
	
	histo2D["Phipart_vs_Phinonbjet"]             =new TH2F("Phipart_vs_Phinonbjet","Phipart_vs_Phinonbjet", 60, -3.2, 3.2, 60, -3.2, 3.2);
	histo2D["Phipart_vs_DiffPhipartPhinonbjet"]  =new TH2F("Phipart_vs_DiffPhipartPhinonbjet","Phipart_vs_Phipart-Phinonbjet", 10, -3.2, 3.2, 100, -0.2, 0.2);
        histo2D["Phipart_vs_DiffPhipartPhinonbjetFULL"]=new TH2F("Phipart_vs_DiffPhipartPhinonbjetFULL","Phipart_vs_Phipart-Phinonbjet (outliers as well)", 10, -3.2, 3.2, 120, -6.2, 6.2);
	histo2D["Phipart_vs_Phibjet"]                =new TH2F("Phipart_vs_Phibjet","Phipart_vs_Phibjet", 60, -3.2, 3.2, 60, -3.2, 3.2);
	histo2D["Phipart_vs_DiffPhipartPhibjet"]     =new TH2F("Phipart_vs_DiffPhipartPhibjet","Phipart_vs_Phipart-Phibjet", 10, -3.2, 3.2, 100, -0.2, 0.2);
        histo2D["Phipart_vs_DiffPhipartPhibjetFULL"] =new TH2F("Phipart_vs_DiffPhipartPhibjetFULL","Phipart_vs_Phipart-Phibjet (outliers as well)", 10, -3.2, 3.2, 120, -6.2, 6.2);
	histo2D["PhigenEl_vs_PhirecEl"]              =new TH2F("PhigenEl_vs_PhirecEl","PhigenEl_vs_PhirecEl", 60, -3.2, 3.2, 60, -3.2, 3.2);
	histo2D["PhigenEl_vs_DiffPhigenElPhirecEl"]  =new TH2F("PhigenEl_vs_DiffPhigenElPhirecEl","PhigenEl_vs_PhigenEl-PhirecEl", 10, -3.2, 3.2, 75, -2, 2);
        histo2D["PhigenEl_vs_DiffPhigenElPhirecElFULL"]  =new TH2F("PhigenEl_vs_DiffPhigenElPhirecElFULL","PhigenEl_vs_PhigenEl-PhirecEl (outliers as well)", 10, -3.2, 3.2, 80, -3.2, 3.2);
	histo2D["PhigenMu_vs_PhirecMu"]              =new TH2F("PhigenMu_vs_PhirecMu","PhigenMu_vs_PhirecMu", 60, -3.2, 3.2, 60, -3.2, 3.2);
	histo2D["PhigenMu_vs_DiffPhigenMuPhirecMu"]  =new TH2F("PhigenMu_vs_DiffPhigenMuPhirecMu","PhigenMu_vs_PhigenMu-PhirecMu", 10, -3.2, 3.2, 75, -2, 2);
        histo2D["PhigenMu_vs_DiffPhigenMuPhirecMuFULL"]  =new TH2F("PhigenMu_vs_DiffPhigenMuPhirecMuFULL","PhigenMu_vs_PhigenMu-PhirecMu (outliers as well)", 10, -3.2, 3.2, 80, -3.2, 3.2);
	histo2D["Epart_vs_Phinonbjet"]               =new TH2F("Epart_vs_Phinonbjet","Epart_vs_Phinonbjet", 75, 0, 250, 60, -3.2, 3.2);
	histo2D["Epart_vs_DiffPhipartPhinonbjet"]    =new TH2F("Epart_vs_DiffPhipartPhinonbjet","Epart_vs_Phipart-Phinonbjet", 10, 20, 250, 100, -0.15, 0.15);
        histo2D["Epart_vs_DiffPhipartPhinonbjetFULL"]    =new TH2F("Epart_vs_DiffPhipartPhinonbjetFULL","Epart_vs_Phipart-Phinonbjet (outliers as well)", 10, 20, 250, 120, -6.2, 6.2);
	histo2D["Epart_vs_Phibjet"]                  =new TH2F("Epart_vs_Phibjet","Epart_vs_Phibjet", 75, 0, 250, 60, -3.2, 3.2);
	histo2D["Epart_vs_DiffPhipartPhibjet"]       =new TH2F("Epart_vs_DiffPhipartPhibjet","Epart_vs_Phipart-Phibjet", 10, 20, 250, 100, -0.15, 0.15);
        histo2D["Epart_vs_DiffPhipartPhibjetFULL"]       =new TH2F("Epart_vs_DiffPhipartPhibjetFULL","Epart_vs_Phipart-Phibjet (outliers as well)", 10, 0, 300, 120, -6.2, 6.2);
	histo2D["EgenEl_vs_PhirecEl"]                =new TH2F("EgenEl_vs_PhirecEl","EgenEl_vs_PhirecEl", 150, 38, 42, 60, -3.2, 3.2);
	histo2D["EgenEl_vs_DiffPhigenElPhirecEl"]    =new TH2F("EgenEl_vs_DiffPhigenElPhirecEl","EgenEl_vs_PhigenEl-PhirecEl", 10, 38, 42, 150, -3, 3);
        histo2D["EgenEl_vs_DiffPhigenElPhirecElFULL"]    =new TH2F("EgenEl_vs_DiffPhigenElPhirecElFULL","EgenEl_vs_PhigenEl-PhirecEl (outliers as well)", 10, 0, 400, 120, -6.2, 6.2);
	histo2D["InvPtgenMu_vs_PhirecMu"]            =new TH2F("InvPtgenMu_vs_PhirecMu","InvPtgenMu_vs_PhirecMu", 80, 0.015, 0.045, 60, -3.2, 3.2);
	histo2D["InvPtgenMu_vs_DiffPhigenMuPhirecMu"]=new TH2F("InvPtgenMu_vs_DiffPhigenMuPhirecMu","InvPtgenMu_vs_PhigenMu-PhirecMu", 10, 0.022, 0.04, 100, -2, 2);
        histo2D["InvPtgenMu_vs_DiffPhigenMuPhirecMuFULL"]=new TH2F("InvPtgenMu_vs_DiffPhigenMuPhirecMuFULL","InvPtgenMu_vs_PhigenMu-PhirecMu (outliers as well)", 10, 0.022, 0.04, 120, -6.2, 6.2);
}

void TFCreation::FillHistograms(TRootMCParticle* hadrWJet1, TRootMCParticle* hadrWJet2, TRootMCParticle* hadrBJet, TRootMCParticle* leptBJet, TRootMCParticle* lepton, TRootJet* selHadrWJet1, TRootJet* selHadrWJet2, TRootJet* selHadrBJet, TRootJet* selLeptBJet, TLorentzVector* selLepton, bool isSemiMu, bool isSemiEl){

	string GeneratorInfo[] = {"EtQs","ThQs"};
	
	histo2D["Epart_vs_Enonbjet"]->Fill(hadrWJet1->E(),selHadrWJet1->E());
	histo2D["Epart_vs_DiffEpartEnonbjet"]->Fill(hadrWJet1->E(),hadrWJet1->E()-selHadrWJet1->E());
	histo2D["Epart_vs_Enonbjet"]->Fill(hadrWJet2->E(),selHadrWJet2->E());
	histo2D["Epart_vs_DiffEpartEnonbjet"]->Fill(hadrWJet2->E(),hadrWJet2->E()-selHadrWJet2->E());
	histo2D["Epart_vs_Ebjet"]->Fill(hadrBJet->E(),selHadrBJet->E());
	histo2D["Epart_vs_DiffEpartEbjet"]->Fill(hadrBJet->E(),hadrBJet->E()-selHadrBJet->E());
	histo2D["Epart_vs_Ebjet"]->Fill(leptBJet->E(),selLeptBJet->E());
	histo2D["Epart_vs_DiffEpartEbjet"]->Fill(leptBJet->E(),leptBJet->E()-selLeptBJet->E());
	if(isSemiEl){
		histo2D["EgenEl_vs_ErecEl"]->Fill(lepton->E(),selLepton->E());
	    	histo2D["EgenEl_vs_DiffEgenElErecEl"]->Fill(lepton->E(),lepton->E()-selLepton->E());
	}
	if(isSemiMu){
		float InvPtgenMu = 1./lepton->Pt();
	    	float InvPtrecMu = 1./selLepton->Pt();
	    	histo2D["InvPtgenMu_vs_InvPtrecMu"]->Fill(InvPtgenMu,InvPtrecMu);
	    	histo2D["InvPtgenMu_vs_DiffInvPtgenMuInvPtrecMu"]->Fill(InvPtgenMu,InvPtgenMu-InvPtrecMu);
        	histo2D["PtgenMu_vs_PtrecMu"]->Fill(lepton->Pt(), selLepton->Pt());
	        histo2D["PtgenMu_vs_DiffPtgenMuPtrecMu"]->Fill(lepton->Pt(), lepton->Pt()-selLepton->Pt());
	}

	//angles
	histo2D["Thpart_vs_Thnonbjet"]->Fill(hadrWJet1->Theta(),selHadrWJet1->Theta());
	histo2D["Thpart_vs_DiffThpartThnonbjet"]->Fill(hadrWJet1->Theta(),hadrWJet1->Theta()-selHadrWJet1->Theta());
	histo2D["Thpart_vs_Thnonbjet"]->Fill(hadrWJet2->Theta(),selHadrWJet2->Theta());
	histo2D["Thpart_vs_DiffThpartThnonbjet"]->Fill(hadrWJet2->Theta(),hadrWJet2->Theta()-selHadrWJet2->Theta());
	histo2D["Thpart_vs_Thbjet"]->Fill(hadrBJet->Theta(),selHadrBJet->Theta());
	histo2D["Thpart_vs_DiffThpartThbjet"]->Fill(hadrBJet->Theta(),hadrBJet->Theta()-selHadrBJet->Theta());
	histo2D["Thpart_vs_Thbjet"]->Fill(leptBJet->Theta(),selLeptBJet->Theta());
	histo2D["Thpart_vs_DiffThpartThbjet"]->Fill(leptBJet->Theta(),leptBJet->Theta()-selLeptBJet->Theta());
	histo2D["Epart_vs_Thnonbjet"]->Fill(hadrWJet1->E(),selHadrWJet1->Theta());
	histo2D["Epart_vs_DiffThpartThnonbjet"]->Fill(hadrWJet1->E(),hadrWJet1->Theta()-selHadrWJet1->Theta());
	histo2D["Epart_vs_Thnonbjet"]->Fill(hadrWJet2->E(),selHadrWJet2->Theta());
	histo2D["Epart_vs_DiffThpartThnonbjet"]->Fill(hadrWJet2->E(),hadrWJet2->Theta()-selHadrWJet2->Theta());
	histo2D["Epart_vs_Thbjet"]->Fill(hadrBJet->E(),selHadrBJet->Theta());
	histo2D["Epart_vs_DiffThpartThbjet"]->Fill(hadrBJet->E(),hadrBJet->Theta()-selHadrBJet->Theta());
	histo2D["Epart_vs_Thbjet"]->Fill(leptBJet->E(),selLeptBJet->Theta());
	histo2D["Epart_vs_DiffThpartThbjet"]->Fill(leptBJet->E(),leptBJet->Theta()-selLeptBJet->Theta());
	if(isSemiEl){
		histo2D["ThgenEl_vs_ThrecEl"]->Fill(lepton->Theta(),selLepton->Theta());
	    	histo2D["ThgenEl_vs_DiffThgenElThrecEl"]->Fill(lepton->Theta(),lepton->Theta()-selLepton->Theta());
	    	histo2D["EgenEl_vs_ThrecEl"]->Fill(lepton->E(),selLepton->Theta());
	    	histo2D["EgenEl_vs_DiffThgenElThrecEl"]->Fill(lepton->E(),lepton->Theta()-selLepton->Theta());
	}
	if(isSemiMu){
		histo2D["ThgenMu_vs_ThrecMu"]->Fill(lepton->Theta(),selLepton->Theta());
		histo2D["ThgenMu_vs_DiffThgenMuThrecMu"]->Fill(lepton->Theta(),lepton->Theta()-selLepton->Theta());
	    	histo2D["InvPtgenMu_vs_ThrecMu"]->Fill(1./lepton->Pt(),selLepton->Theta());
	    	histo2D["InvPtgenMu_vs_DiffThgenMuThrecMu"]->Fill(1./lepton->Pt(),lepton->Theta()-selLepton->Theta());
	}
	
	histo2D["Phipart_vs_Phinonbjet"]->Fill(hadrWJet1->Phi(),selHadrWJet1->Phi());
	histo2D["Phipart_vs_Phinonbjet"]->Fill(hadrWJet2->Phi(),selHadrWJet2->Phi());
	histo2D["Epart_vs_Phinonbjet"]->Fill(hadrWJet1->E(),selHadrWJet1->Phi());
	histo2D["Epart_vs_Phinonbjet"]->Fill(hadrWJet2->E(),selHadrWJet2->Phi());
	float DeltaPhi_nonbjet1 = ROOT::Math::VectorUtil::DeltaPhi((TLorentzVector)*hadrWJet1,(TLorentzVector)*selHadrWJet1);
	histo2D["Phipart_vs_DiffPhipartPhinonbjet"]->Fill(hadrWJet1->Phi(),DeltaPhi_nonbjet1);
	histo2D["Epart_vs_DiffPhipartPhinonbjet"]->Fill(hadrWJet1->E(),DeltaPhi_nonbjet1);
        histo2D["Phipart_vs_DiffPhipartPhinonbjetFULL"]->Fill(hadrWJet1->Phi(),DeltaPhi_nonbjet1);
        histo2D["Epart_vs_DiffPhipartPhinonbjetFULL"]->Fill(hadrWJet1->E(),DeltaPhi_nonbjet1);
	float DeltaPhi_nonbjet2 = ROOT::Math::VectorUtil::DeltaPhi((TLorentzVector)*hadrWJet2,(TLorentzVector)*selHadrWJet2);
	histo2D["Phipart_vs_DiffPhipartPhinonbjet"]->Fill(hadrWJet2->Phi(),DeltaPhi_nonbjet2);
	histo2D["Epart_vs_DiffPhipartPhinonbjet"]->Fill(hadrWJet2->E(),DeltaPhi_nonbjet2);
	histo2D["Phipart_vs_DiffPhipartPhinonbjetFULL"]->Fill(hadrWJet2->Phi(),DeltaPhi_nonbjet2);
	histo2D["Epart_vs_DiffPhipartPhinonbjetFULL"]->Fill(hadrWJet2->E(),DeltaPhi_nonbjet2);
	histo2D["Phipart_vs_Phibjet"]->Fill(hadrBJet->Phi(),selHadrBJet->Phi());
	histo2D["Phipart_vs_Phibjet"]->Fill(leptBJet->Phi(),selLeptBJet->Phi());
	histo2D["Epart_vs_Phibjet"]->Fill(hadrBJet->E(),selHadrBJet->Phi());
	histo2D["Epart_vs_Phibjet"]->Fill(leptBJet->E(),selLeptBJet->Phi());
	float DeltaPhi_bjet1 = ROOT::Math::VectorUtil::DeltaPhi((TLorentzVector)*hadrBJet,(TLorentzVector)*selLeptBJet);				
	histo2D["Phipart_vs_DiffPhipartPhibjet"]->Fill(hadrBJet->Phi(),DeltaPhi_bjet1);
	histo2D["Epart_vs_DiffPhipartPhibjet"]->Fill(hadrBJet->E(),DeltaPhi_bjet1);
	histo2D["Phipart_vs_DiffPhipartPhibjetFULL"]->Fill(hadrBJet->Phi(),DeltaPhi_bjet1);
	histo2D["Epart_vs_DiffPhipartPhibjetFULL"]->Fill(hadrBJet->E(),DeltaPhi_bjet1);
	float DeltaPhi_bjet2 = ROOT::Math::VectorUtil::DeltaPhi((TLorentzVector)*leptBJet,(TLorentzVector)*selLeptBJet);				
	histo2D["Phipart_vs_DiffPhipartPhibjet"]->Fill(leptBJet->Phi(),DeltaPhi_bjet2);
	histo2D["Epart_vs_DiffPhipartPhibjet"]->Fill(leptBJet->E(),DeltaPhi_bjet2);
	histo2D["Phipart_vs_DiffPhipartPhibjetFULL"]->Fill(leptBJet->Phi(),DeltaPhi_bjet2);
	histo2D["Epart_vs_DiffPhipartPhibjetFULL"]->Fill(leptBJet->E(),DeltaPhi_bjet2);
	if(isSemiEl){ 
		histo2D["PhigenEl_vs_PhirecEl"]->Fill(lepton->Phi(),selLepton->Phi());
	    	histo2D["EgenEl_vs_PhirecEl"]->Fill(lepton->E(),selLepton->Phi());
	    	float DeltaPhi = ROOT::Math::VectorUtil::DeltaPhi((TLorentzVector)*lepton,(TLorentzVector)*selLepton);
	    	histo2D["PhigenEl_vs_DiffPhigenElPhirecEl"]->Fill(lepton->Phi(),DeltaPhi);
	    	histo2D["EgenEl_vs_DiffPhigenElPhirecEl"]->Fill(lepton->E(),DeltaPhi);
	    	histo2D["PhigenEl_vs_DiffPhigenElPhirecElFULL"]->Fill(lepton->Phi(),DeltaPhi);
	    	histo2D["EgenEl_vs_DiffPhigenElPhirecElFULL"]->Fill(lepton->E(),DeltaPhi);
	}
	if(isSemiMu){
		histo2D["PhigenMu_vs_PhirecMu"]->Fill(lepton->Phi(),selLepton->Phi());
	    	histo2D["InvPtgenMu_vs_PhirecMu"]->Fill(1./lepton->Pt(),selLepton->Phi());
	    	float DeltaPhi = ROOT::Math::VectorUtil::DeltaPhi((TLorentzVector)*lepton,(TLorentzVector)*selLepton);
	    	histo2D["PhigenMu_vs_DiffPhigenMuPhirecMu"]->Fill(lepton->Phi(),DeltaPhi);
	    	histo2D["InvPtgenMu_vs_DiffPhigenMuPhirecMu"]->Fill(1./lepton->Pt(),DeltaPhi);
	    	histo2D["PhigenMu_vs_DiffPhigenMuPhirecMuFULL"]->Fill(lepton->Phi(),DeltaPhi);
	    	histo2D["InvPtgenMu_vs_DiffPhigenMuPhirecMuFULL"]->Fill(1./lepton->Pt(),DeltaPhi);
	}	
}

void TFCreation::CalculateTF(bool drawHistos, bool writeTF, bool doFits, bool useROOTClass, bool useStartValues){
	TFile* file = new TFile("PlotsForTransferFunctions.root","RECREATE");
	file->cd();

	if(drawHistos == true) WritePlots(file);
	if(writeTF == true) WriteTF(file);
	if(doFits == true){
  	  ///////////////////////////////////////
	  //  Declare the used fit functions!  //
	  ///////////////////////////////////////
	  //
	  // 1) Double Gaussian --> its range depends on the jet/lepton energy range (hence, the Y-axis)
	  doubleGaussianFit = new TF1("doubleGaussianFit","[2]*(TMath::Exp(-TMath::Power((x-[0]),2)/(2*TMath::Power([1],2)))+[5]*TMath::Exp(-TMath::Power((x-[3]),2)/(2*TMath::Power([4],2))))");	
	  //give names to the parameters
	  const char* parnames[6]={"a1","a2","a3","a4","a5","a6"};
	  const int npar = doubleGaussianFit->GetNpar();

	  for(int ii = 0; ii < npar; ii++)
	    doubleGaussianFit->SetParName(ii,parnames[ii]);

	  //2) Calorimeter Energy formula (ai = ai0 + ai1*Ep + ai2*sqrt(Ep)) --> its range depends on the part energy range (hence, the X-axis)
	  caloEnergyFit = new TF1("caloEnergyFit", "[0]+[1]*sqrt(x)+[2]*x");	

	  ///////////////////////////////////////////
	  //  Choose the correct histogram to fit  //
	  ///////////////////////////////////////////
          TH2F* histoForFit;
	  for (unsigned int f=0; f<12;f++) {				
	    //if(f==2 || f==6 || f==10) continue; //electron plots not filled at the moment...
				
	    switch(f){
	      case 0:
		histoForFit=histo2D["Epart_vs_DiffEpartEnonbjet"]; break;
	      case 1:
		histoForFit=histo2D["Epart_vs_DiffEpartEbjet"]; break;
	      case 2:
		histoForFit=histo2D["EgenEl_vs_DiffEgenElErecEl"]; break;
	      case 3:
		histoForFit=histo2D["InvPtgenMu_vs_DiffInvPtgenMuInvPtrecMu"]; break;
	      case 4:
		histoForFit=histo2D["Epart_vs_DiffThpartThnonbjet"]; break;
	      case 5:
		histoForFit=histo2D["Epart_vs_DiffThpartThbjet"]; break;
	      case 6:
		histoForFit=histo2D["EgenEl_vs_DiffThgenElThrecEl"]; break;
	      case 7:
		histoForFit=histo2D["InvPtgenMu_vs_DiffThgenMuThrecMu"]; break;
	      case 8:
		histoForFit=histo2D["Epart_vs_DiffPhipartPhinonbjet"]; break;
	      case 9:
		histoForFit=histo2D["Epart_vs_DiffPhipartPhibjet"]; break;
	      case 10:
		histoForFit=histo2D["EgenEl_vs_DiffPhigenElPhirecEl"]; break;
	      case 11:
		histoForFit=histo2D["InvPtgenMu_vs_DiffPhigenMuPhirecMu"]; break;
	    }					

            TDirectory* histoFitDir = file->mkdir(histoForFit->GetName());
            histoFitDir->cd();
 
	    hlist = new TH1D*[npar];
	    if(useStartValues)
		SetStartValuesDoubleGaussian(f);

	    doubleGaussianFit->SetRange( histoForFit->GetYaxis()->GetXmin(), histoForFit->GetYaxis()->GetXmin() );
	    TObjArray aSlices;
	    if(useROOTClass){
		//Fit using the FitSliceY function of TF1!
		histoForFit->FitSlicesY(doubleGaussianFit, 0, -1, 0, "", &aSlices);
	        for(int ipar = 0; ipar <= npar; ipar++)
		    hlist[ipar] = (TH1D*) aSlices[ipar];
	    }
	    else
		FitSliceClassCode(histoForFit, npar, parnames);

	    //////////////////////////////////////////////////////////////////////////////////////////////
	    //   Now histogram with all parameters needs to be fitted with Calorimeter Energy formula   //
	    //////////////////////////////////////////////////////////////////////////////////////////////
	    caloEnergyFit->SetRange( histoForFit->GetXaxis()->GetXmin(), histoForFit->GetXaxis()->GetXmax() );						
	    for( int ipar = 0; ipar < npar; ipar++ ){

	  	//give names to the parameters		
		caloEnergyFit->SetParName(0, ( string(parnames[ipar])+"0").c_str() );
		caloEnergyFit->SetParName(1, ( string(parnames[ipar])+"1").c_str());
		caloEnergyFit->SetParName(2, ( string(parnames[ipar])+"2").c_str());

		caloEnergyFit->SetName( (string(histoForFit->GetName())+"_"+parnames[ipar]+"_Fit").c_str() );
		hlist[ipar]->SetName( (string(histoForFit->GetName())+"_"+parnames[ipar]+"_PointsAndFit").c_str() );

		hlist[ipar]->Fit(caloEnergyFit);
		hlist[ipar]->Write();            
		//caloEnergyFit->Write();
	    }
	    hlist[npar]->Write();
						
	  }//Loop over f						
	  delete caloEnergyFit, doubleGaussianFit;
	  delete histoForFit;
	  delete [] hlist;
	}                               //Boolean doFits = true
	file->Close();
}

void TFCreation::FitSliceClassCode(TH2F* histoForFit, int npar, const char* parnames[]){
	//------------------------------------------------------------------------------------------//
	// Main difference with the Root class FitSlicesY() is the plotting of histograms !        
	// In the Root class the distribution of each hlist histogram is not given!
	// --> Useful to use the own code when needing control histograms
	//
	// Other difference between the two codes have been removed!
	// Originally the treatment of the overflow bin was different, but is now made similar!
	//Create one histogram for each function parameter -> 6 histograms for each 2D plot
	const TArrayD *bins = histoForFit->GetXaxis()->GetXbins();
	for(int ipar=0 ; ipar < npar; ipar++){

	    hlist[ipar] = new TH1D( (string(histoForFit->GetName())+"_"+parnames[ipar]).c_str(), (string(histoForFit->GetName())+" : Fitted value of _"+parnames[ipar]).c_str(), histoForFit->GetXaxis()->GetNbins(), histoForFit->GetXaxis()->GetXmin(), histoForFit->GetXaxis()->GetXmax());
	    hlist[ipar]->GetXaxis()->SetTitle(histoForFit->GetXaxis()->GetTitle());
	}
	hlist[npar] = new TH1D( (string(histoForFit->GetName())+"_chi2").c_str(), (string(histoForFit->GetName())+": #chi^{2} distribution for "+string(doubleGaussianFit->GetExpFormula())).c_str(), histoForFit->GetXaxis()->GetNbins(), histoForFit->GetXaxis()->GetXmin(), histoForFit->GetXaxis()->GetXmax() );

	//Loop on all bins in X, generate a projection along Y and fit each bin separately!
	int cut = 0; // require a minimum number of bins in the slice to be filled --> Should this ever be larger than 0 ??
	int nbins = histoForFit->GetXaxis()->GetNbins();
	cout << " ** Looking at histogram : " << histoForFit->GetName() << "                            ** " << endl;
	for(int bin=1;bin <= nbins+1;bin ++) {
	    string projection_title = string(histoForFit->GetName())+"_sliceXbin"+tostr(bin);

	    TH1D *hp = histoForFit->ProjectionY(projection_title.c_str(),bin,bin,"e");
	    //if(bin==nbins) hp = histoForFit->ProjectionY(projection_title.c_str(),bin,bin+1,"e"); //include overflow in last bin
	    //if(bin==1) hp = histoForFit->ProjectionY(projection_title.c_str(),bin-1,bin,"e"); //include underflow in first bin

	    //Histogram doesn't have any memory space ...
	    if(hp == 0) continue;
	    if( float(hp->GetEntries()) <= 0){ delete hp; continue;} //|| float(hp->GetEntries()) < cut) {delete hp; continue;}

	    doubleGaussianFit->SetName((projection_title+"Fitted").c_str());
	    hp->Fit(doubleGaussianFit);

	    int npfits = doubleGaussianFit->GetNumberFitPoints();              //WHAT IS THIS .... ???
	    if(npfits > npar && npfits >= cut) {

		cout << "   --  Looking at bin : " << bin << endl;
		//Fill the hlist histogram for each parameter with the obtained Fit parameter and its uncertainty
	        //--> Each bin in this histogram represents a bin range in x-axis of considered 2D histogram!
	        for(int ipar=0; ipar<npar; ipar++ ){
		    hlist[ipar]->Fill(histoForFit->GetXaxis()->GetBinCenter(bin+1/2),doubleGaussianFit->GetParameter(ipar)); 
                    hlist[ipar]->SetBinError( (int) (bin+1/2) ,doubleGaussianFit->GetParError(ipar)); //WHY +1/2 .... (Is bin size always equal to 1 .. )?

	        }
		//Save hchi2 histogram as extra hlist!
	        hlist[npar]->Fill(histoForFit->GetXaxis()->GetBinCenter(bin+1/2),doubleGaussianFit->GetChisquare()/(npfits-npar));
	    }
	    hp->Write();
	    //doubleGaussianFit->Write();     //--> Interesting to draw this as well ?? Should appear in the histogram ...
	    delete hp;
	}//loop over bins!
}	

void TFCreation::SetStartValuesDoubleGaussian(int whichHisto){

	if(whichHisto==0 || whichHisto==1){ // for E transfer function of JETS
	    float StartValues[] = {-8,18,63,0,8.6,4.1};      //First three values are for the first broad gaussian (central, sigma and constant value respectively)
							     //Second three values are the same for the second narrow gaussian
	    for(int ii = 0; ii < 6; ii++)
		doubleGaussianFit->SetParameter(ii,StartValues[ii]);
	}
	else if (whichHisto==4 || whichHisto==5 || whichHisto==8 || whichHisto==9) { //for theta and phi transfer functions of JETS
	    float StartValues[] = {0,0.038,77,0.004,0.011,6.5};
	    for(int ii = 0; ii < 6; ii++)
		doubleGaussianFit->SetParameter(ii, StartValues[ii]);
	}
	else if (whichHisto==3){ //for 1/pt transfer function of muons
            float StartValues[] = {-0.0008,0.001,24,-0.0001,0.0001,4};
            for(int ii = 0; ii < 6; ii++)
                doubleGaussianFit->SetParameter(ii, StartValues[ii]);
	}
	else if (whichHisto==7 || whichHisto==11) { //for theta, phi transfer function of muons
            float StartValues[] = {0.0,0.01,24,0,0.001,4};
            for(int ii = 0; ii < 6; ii++)
                doubleGaussianFit->SetParameter(ii, StartValues[ii]);
	}
} 

void TFCreation::WriteTF(TFile* plotsFile){

	const int NrConsideredPlots = 9;
	const int NrConsideredPars = 6;
	const int NrConsideredCaloPars = 3;

	string histonames[9] = {"Epart_vs_DiffEpartEbjet",                   //Why is order suddenly changed ... ?   + Still need to include electrons !!
				"Epart_vs_DiffThpartThbjet",
				"Epart_vs_DiffPhipartPhibjet",
				"Epart_vs_DiffEpartEnonbjet",
				"Epart_vs_DiffThpartThnonbjet",
				"Epart_vs_DiffPhipartPhinonbjet",
				"InvPtgenMu_vs_DiffInvPtgenMuInvPtrecMu",
				"InvPtgenMu_vs_DiffThgenMuThrecMu",
				"InvPtgenMu_vs_DiffPhigenMuPhirecMu"};
				
	string histodescription[9] = {"b-jet energy", "b-jet theta", "b-jet phi", "non-b jet energy", "non-b jet theta", "non-b jet phi", "muon inv. pt", "muon theta", "muon phi",};
				
	ofstream myTFs;
	myTFs.open("TransferFunctions_TABLE.txt");

	ofstream myTFsForMadWeight;
	myTFsForMadWeight.open("transfer_card_user.dat");

	plotsFile->cd();
	string ParamName[NrConsideredPars] = {"Mean broad gaussian", "Width broad gaussian","Constant broad gaussian","Mean narrow gaussian","Width narrow gaussian","Constant narrow gaussian"};
	for(int iplot = 0; iplot < NrConsideredPlots; iplot++){

	  myTFs<< endl;
	  myTFs<<"\\begin{table}" << endl;
	  myTFs<<"\\caption{Parameters of the transfer function for the " << histodescription[iplot]  << "}" << endl;
	  myTFs<<"\\label{tab:}" << endl;
	  myTFs<<"\\centering" << endl;
	  myTFs<<"\\begin{tabular}{c|ccc}" << endl;
	  myTFs<<"\\hline" << endl;
	  myTFs << "Type	& $a_{i0}$ & $a_{i1}$ ($\\sqrt{E}$) & $a_{i2}$ ($E$)" << "\\\\" << endl;
	  myTFs<<"\\hline" << endl;
	
	  TF1 *TF_Created[6];                                                     //,*TF_par2,*TF_par3,*TF_par4,*TF_par5,*TF_par6;
	  for(int ipar = 1; ipar <= NrConsideredPars; ipar++){
		TF_Created[ipar-1] = (TF1*)plotsFile->Get( (histonames[iplot]+"_a"+tostr(ipar)+"_Fitted").c_str() );
		if( !TF_Created[ipar-1]) continue;

		for(int icalopar = 0; icalopar < NrConsideredCaloPars; icalopar++){
		    if(icalopar == 0) myTFs <<ParamName[ipar-1] <<" & $a_{" <<ipar <<icalopar <<"}$ = " <<TF_Created[ipar-1]->GetParameter(icalopar) <<"$\\pm$" <<TF_Created[ipar-1]->GetParError(icalopar);
		    else              myTFs <<                    " & $a_{" <<ipar <<icalopar <<"}$ = " <<TF_Created[ipar-1]->GetParameter(icalopar) <<"$\\pm$" <<TF_Created[ipar-1]->GetParError(icalopar);
		}
		myTFs << "\\\\" << endl;
	  }
	  //string name2 = histonames[i]+"_a2_Fitted";
	  //string name3 = histonames[i]+"_a3_Fitted";
	  //string name4 = histonames[i]+"_a4_Fitted";
	  //string name5 = histonames[i]+"_a5_Fitted";
	  //string name6 = histonames[i]+"_a6_Fitted";
					
	  //TF_par1 = (TF1*)plotsFile->Get(name1.c_str());
	  //TF_par2 = (TF1*)plotsFile->Get(name2.c_str());
	  //TF_par3 = (TF1*)plotsFile->Get(name3.c_str());
	  //TF_par4 = (TF1*)plotsFile->Get(name4.c_str());
	  //TF_par5 = (TF1*)plotsFile->Get(name5.c_str());
	  //TF_par6 = (TF1*)plotsFile->Get(name6.c_str());

	  //if (TF_par1 && TF_par2 && TF_par3 && TF_par4 && TF_par5 && TF_par6) {

	    /*myTFs << "Mean broad gaussian & $a_{10}$ = " << TF_par1->GetParameter(0) << "$\\pm$" << TF_par1->GetParError(0);
	    myTFs <<                    " & $a_{11}$ = " << TF_par1->GetParameter(1) << "$\\pm$" << TF_par1->GetParError(1);
	    myTFs <<                    " & $a_{12}$ = " << TF_par1->GetParameter(2) << "$\\pm$" << TF_par1->GetParError(2) << "\\\\" << endl;
	    myTFs << "Width broad gaussian & $a_{20}$ = " << TF_par2->GetParameter(0) << "$\\pm$" << TF_par2->GetParError(0);
	    myTFs <<                     " & $a_{21}$ = " << TF_par2->GetParameter(1) << "$\\pm$" << TF_par2->GetParError(1);
	    myTFs <<                     " & $a_{22}$ = " << TF_par2->GetParameter(2) << "$\\pm$" << TF_par2->GetParError(2) << "\\\\" << endl;
	    myTFs << "Constant broad gaussian & $a_{30}$ = " << TF_par2->GetParameter(0) << "$\\pm$" << TF_par2->GetParError(0);
	    myTFs <<                        " & $a_{31}$ = " << TF_par2->GetParameter(1) << "$\\pm$" << TF_par2->GetParError(1);
	    myTFs <<                        " & $a_{32}$ = " << TF_par2->GetParameter(2) << "$\\pm$" << TF_par2->GetParError(2) << "\\\\" << endl;
	    myTFs << "Mean narrow gaussian & $a_{40}$ = " << TF_par4->GetParameter(0) << "$\\pm$" << TF_par4->GetParError(0);
	    myTFs <<                     " & $a_{41}$ = " << TF_par4->GetParameter(1) << "$\\pm$" << TF_par4->GetParError(1);
	    myTFs <<                     " & $a_{42}$ = " << TF_par4->GetParameter(2) << "$\\pm$" << TF_par4->GetParError(2) << "\\\\" << endl;
	    myTFs << "Width narrow gaussian & $a_{50}$ = " << TF_par5->GetParameter(0) << "$\\pm$" << TF_par5->GetParError(0);
	    myTFs <<                      " & $a_{51}$ = " << TF_par5->GetParameter(1) << "$\\pm$" << TF_par5->GetParError(1);
	    myTFs <<                      " & $a_{52}$ = " << TF_par5->GetParameter(2) << "$\\pm$" << TF_par5->GetParError(2) << "\\\\" << endl;
	    myTFs << "Constant narrow gaussian & $a_{60}$ = " << TF_par6->GetParameter(0) << "$\\pm$" << TF_par6->GetParError(0);
	    myTFs <<                         " & $a_{61}$ = " << TF_par6->GetParameter(1) << "$\\pm$" << TF_par6->GetParError(1);
	    myTFs <<                         " & $a_{62}$ = " << TF_par6->GetParameter(2) << "$\\pm$" << TF_par6->GetParError(2) << "\\\\" << endl;
	    */
	    myTFs<<"\\hline" << endl;
	    myTFs<<"\\end{tabular}"<<endl;
	    myTFs<<"\\end{table}"<<endl;
	    myTFs<< endl;
						
	    if(iplot==0){
		myTFsForMadWeight<<"#+-----------------------------------------------------------------------+" << endl;
		myTFsForMadWeight<<"#|    Parameter for particles: b                                         |" << endl;
		myTFsForMadWeight<<"#+-----------------------------------------------------------------------+" << endl;
		myTFsForMadWeight<<"BLOCK TF_bjet_E" << endl;
	    }
	    else if(iplot==3){
		myTFsForMadWeight<<"#+-----------------------------------------------------------------------+" << endl;
		myTFsForMadWeight<<"#|    Parameter for particles: nonb                                      |" << endl;
		myTFsForMadWeight<<"#+-----------------------------------------------------------------------+" << endl;
		myTFsForMadWeight<<"BLOCK TF_nonbjet_E" << endl;
	    }
	    else if(iplot==6){
		myTFsForMadWeight<<"#+-----------------------------------------------------------------------+" << endl;
		myTFsForMadWeight<<"#|    Parameter for particles: muon                                      |" << endl;
		myTFsForMadWeight<<"#+-----------------------------------------------------------------------+" << endl;
		myTFsForMadWeight<<"BLOCK TF_muon_InvPt" << endl;
	    }
	    else if(iplot==1) myTFsForMadWeight<<"BLOCK TF_bjet_THETA" << endl;
	    else if(iplot==2) myTFsForMadWeight<<"BLOCK TF_bjet_PHI" << endl;
	    else if(iplot==4) myTFsForMadWeight<<"BLOCK TF_nonbjet_THETA" << endl;
	    else if(iplot==5) myTFsForMadWeight<<"BLOCK TF_nonbjet_PHI" << endl;
	    else if(iplot==7) myTFsForMadWeight<<"BLOCK TF_muon_THETA" << endl;
	    else if(iplot==8) myTFsForMadWeight<<"BLOCK TF_muon_PHI" << endl;

	    /*
	    if(iplot<6) { //for jets, should also work for electrons if we have them
	      for(int j = 0; j<3; j++)
		myTFsForMadWeight<<j+1 << " " << TF_par1->GetParameter(j) << "  # bias broad gaussian b1=#1+#2*sqrt(E)*#3*E" << endl;
	      for(int j = 3; j<6; j++)
		myTFsForMadWeight<<j+1 << " " << TF_par2->GetParameter(j-3) << "  # sigma broad gaussian s1=#4+#5*sqrt(E)*#6*E" << endl;
	      for(int j = 6; j<9; j++)
		myTFsForMadWeight<<j+1 << " " << TF_par3->GetParameter(j-6) << "  # constant broad gaussian c1=#7+#8*sqrt(E)*#9*E" << endl;
	      for(int j = 9; j<12; j++)
		myTFsForMadWeight<<j+1 << " " << TF_par4->GetParameter(j-9) << "  # bias narrow gaussian b2=#10+#11*sqrt(E)*#12*E" << endl;
	      for(int j = 12; j<15; j++)
		myTFsForMadWeight<<j+1 << " " << TF_par5->GetParameter(j-12) << "  # sigma narrow gaussian s2=#13+#14*sqrt(E)*#15*E" << endl;
	      for(int j = 15; j<18; j++)
		myTFsForMadWeight<<j+1 << " " << TF_par6->GetParameter(j-15) << "  # constant narrow gaussian c2=#16+#17*sqrt(E)*#18*E" << endl;
	    }
	    else{ //for muon (only difference: paramtrization according to InvPt instead of E -> something textual...)
	      for(int j = 0; j<3; j++)
		myTFsForMadWeight<<j+1 << " " << TF_par1->GetParameter(j) << "  # bias broad gaussian b1=#1+#2*sqrt(InvPt)*#3*InvPt" << endl;
	      for(int j = 3; j<6; j++)
		myTFsForMadWeight<<j+1 << " " << TF_par2->GetParameter(j-3) << "  # sigma broad gaussian s1=#4+#5*sqrt(InvPt)*#6*InvPt" << endl;
	      for(int j = 6; j<9; j++)
		myTFsForMadWeight<<j+1 << " " << TF_par3->GetParameter(j-6) << "  # constant broad gaussian c1=#7+#8*sqrt(InvPt)*#9*InvPt" << endl;
	      for(int j = 9; j<12; j++)
		myTFsForMadWeight<<j+1 << " " << TF_par4->GetParameter(j-9) << "  # bias narrow gaussian b2=#10+#11*sqrt(InvPt)*#12*InvPt" << endl;
	      for(int j = 12; j<15; j++)
		myTFsForMadWeight<<j+1 << " " << TF_par5->GetParameter(j-12) << "  # sigma narrow gaussian s2=#13+#14*sqrt(InvPt)*#15*InvPt" << endl;
	      for(int j = 15; j<18; j++)
		myTFsForMadWeight<<j+1 << " " << TF_par6->GetParameter(j-15) << "  # constant narrow gaussian c2=#16+#17*sqrt(InvPt)*#18*InvPt" << endl;
	    }
	    */
	  //}
	  //plotsFile->Close();
	}//End of loop over the different considered plots for the TF!
	myTFs.close();
	myTFsForMadWeight.close();			
}

void TFCreation::WritePlots(TFile* outfile){
	outfile->cd();

	TDirectory* th1dir = outfile->mkdir("1D_histograms");
	th1dir->cd();
	for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++){
	      TH1F *temp = it->second;
	      int N = temp->GetNbinsX();
	      temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
	      temp->SetBinContent(N+1,0);
	      temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
	      temp->Write();
	      //TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
	      //tempCanvas->SaveAs( (it->first+".png").c_str() );
	}
	TDirectory* th2dir = outfile->mkdir("2D_histograms_graphs");
	th2dir->cd();
	for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++){    
	      TH2F *temp = it->second;
	      temp->Write();
	      //TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
	      //tempCanvas->SaveAs( (it->first+".png").c_str() );
	}
	outfile->cd(); //Step out from 2D_histograms_graphs directory!
}
