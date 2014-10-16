#include "../interface/TFCreation.h"

TFCreation::TFCreation(){
    ///////////////////////////////////////
    //  Declare the used fit functions!  //
    ///////////////////////////////////////
    //
    // 1) Double Gaussian --> its range depends on the jet/lepton energy range (hence, the Y-axis)
    doubleGaussianFit = new TF1("doubleGaussianFit","[2]*(TMath::Exp(-TMath::Power((x-[0]),2)/(2*TMath::Power([1],2)))+[5]*TMath::Exp(-TMath::Power((x-[3]),2)/(2*TMath::Power([4],2))))");       

    //2) Calorimeter Energy formula (ai = ai0 + ai1*Ep + ai2*sqrt(Ep)) --> its range depends on the part energy range (hence, the X-axis)
    caloEnergyFit = new TF1("caloEnergyFit", "[0]+[1]*sqrt(x)+[2]*x");
}

TFCreation::~TFCreation(){
    delete caloEnergyFit, doubleGaussianFit;
}

void TFCreation::InitializeVariables(){
        histo1D["DeltaR_TFClass_Light1"] = new TH1F("DeltaR_TFClass_Light1","DeltaR_TFClass_Light1",200,0,1);
        histo1D["DeltaR_TFClass_Light2"] = new TH1F("DeltaR_TFClass_Light2","DeltaR_TFClass_Light2",200,0,1);
        histo1D["DeltaR_TFClass_HadrB"]  = new TH1F("DeltaR_TFClass_HadrB","DeltaR_TFClass_HadrB",200,0,1);
        histo1D["DeltaR_TFClass_LeptB"]  = new TH1F("DeltaR_TFClass_LeptB","DeltaR_TFClass_LeptB",200,0,1);
	histo1D["DeltaR_TFClass_Mu"]     = new TH1F("DeltaR_TFClass_Mu","DeltaR_TFClass_Mu",200,0,1);
	histo1D["DeltaR_TFClass_El"]     = new TH1F("DeltaR_TFClass_El","DeltaR_TFClass_El",200,0,1);

	histo2D["Light_RecoPtVsGenPt"]    = new TH2F("Light_RecoPtVsGenPt",   "Transverse momentum of light quarks (reco vs gen)",                        150,     0,  300, 150,      0,  300);
	histo2D["Light_DiffPtVsGenPt"]    = new TH2F("Light_DiffPtVsGenPt",   "p_{T} difference (gen-reco) versus p_{T,gen} for light quarks",             10,    25,  160,  80,    -40,   40);
	histo2D["BJet_RecoPtVsGenPt"]     = new TH2F("BJet_RecoPtVsGenPt",    "Transverse momentum of b-jets (reco vs gen level)",                        150,     0,  300, 150,      0,  300);
	histo2D["BJet_DiffPtVsGenPt"]     = new TH2F("BJet_DiffPtVsGenPt",    "p_{T} difference (gen-reco) versus p_{T,gen} for b-jets",                   10,    25,  160, 100,    -50,   50);
	histo2D["El_RecoPtVsGenPt"]       = new TH2F("El_RecoPtVsGenPt",      "Transverse momentum of electron (reco vs gen)",                            100,     0,  200, 100,      0,  200);
	histo2D["El_DiffPtVsGenPt"]       = new TH2F("El_DiffPtVsGenPt",      "p_{T} difference (gen-reco) versus p_{T,gen} for electron",                 10,    30,  150,  75,    -10,   10);
	histo2D["Mu_RecoInvPtVsGenInvPt"] = new TH2F("Mu_RecoInvPtVsGenInvPt","Inverse of transverse momentum of muon (reco vs gen)",                     100,     0, 0.05, 100,      0, 0.05);
	histo2D["Mu_DiffInvPtVsGenInvPt"] = new TH2F("Mu_DiffInvPtVsGenInvPt","#frac{1}{p_{T}} difference (gen-reco) versus #frac{1}{p_{T,gen}} for muon", 10, 0.005, 0.04,  80,  -0.01, 0.01);
        histo2D["Mu_DiffInvPtVsGenInvPt_All"] =new TH2F("Mu_DiffInvPtVsGenInvPt_All","#frac{1}{p_{T}} difference (gen-reco) versus #frac{1}{p_{T,gen}} for muon", 20, 0.005, 0.04,  180,  -0.1, 0.1);
	histo2D["Mu_RecoPtVsGenPt"]       = new TH2F("Mu_RecoPtVsGenPt",      "Transverse momentum of muon (reco vs gen)",                                150,     0,  200, 150,      0,  200);
	histo2D["Mu_DiffPtVsGenPt"]       = new TH2F("Mu_DiffPtVsGenPt",      "p_{T} difference (gen-reco) versus p_{T,gen} for muon",                     10,    30,  150, 100,    -10,   10);
	
	histo2D["Light_RecoThetaVsGenTheta"] = new TH2F("Light_RecoThetaVsGenTheta","Polar angle distribution of light quarks (reco vs gen)",                          60,    0,  3.15,  60,     0, 3.15);
	histo2D["Light_DiffThetaVsGenTheta"] = new TH2F("Light_DiffThetaVsGenTheta","#theta difference (gen-reco) versus #theta_{gen} for light quarks",               10,  0.1,   3.1, 100, -0.15, 0.15);
	histo2D["BJet_RecoThetaVsGenTheta"]  = new TH2F("BJet_RecoThetaVsGenTheta", "Polar angle distribution of b-jets (reco vs gen)",                                60,    0,  3.15,  60,     0, 3.15);
	histo2D["BJet_DiffThetaVsGenTheta"]  = new TH2F("BJet_DiffThetaVsGenTheta", "#theta difference (gen-reco) versus #theta_{gen} for b-jets",                     10,  0.1,   3.1, 100, -0.15, 0.15);
	histo2D["El_RecoThetaVsGenTheta"]    = new TH2F("El_RecoThetaVsGenTheta",   "Polar angle distribution of electron (reco vs gen)",                              60,    0,  3.15,  60,     0, 3.15);
	histo2D["El_DiffThetaVsGenTheta"]    = new TH2F("El_DiffThetaVsGenTheta",   "#theta difference (gen-reco) versus #theta_{gen} for electron",                   10,    0,  3.15, 100, -0.15, 0.15);
	histo2D["Mu_RecoThetaVsGenTheta"]    = new TH2F("Mu_RecoThetaVsGenTheta",   "Polar angle distribution of muon (reco vs gen)",                                  60,    0,  3.15,  60,     0, 3.15);
	histo2D["Mu_DiffThetaVsGenTheta"]    = new TH2F("Mu_DiffThetaVsGenTheta",   "#theta difference (gen-reco) versus #theta_{gen} for muon",                       10,    0,  3.15, 100, -0.15, 0.15);
	histo2D["Light_RecoThetaVsGenPt"]    = new TH2F("Light_RecoThetaVsGenPt",   "Polar angle #theta_{rec} versus transverse momentum p_{T,gen} for light quarks", 120,    0,   300,  60,     0, 3.15);
	histo2D["Light_DiffThetaVsGenPt"]    = new TH2F("Light_DiffThetaVsGenPt",   "#theta difference (gen-reco) versus p_{T,gen} for light quarks",                  10,   25,   160, 150, -0.15, 0.15);
	histo2D["BJet_RecoThetaVsGenPt"]     = new TH2F("BJet_RecoThetaVsGenPt",    "Polar angle #theta_{rec} versus transverse momentum p_{T,gen} for b-jets",       120,    0,   300,  60,     0, 3.15);
	histo2D["BJet_DiffThetaVsGenPt"]     = new TH2F("BJet_DiffThetaVsGenPt",    "#theta difference (gen-reco) versus p_{T,gen} for b-jets",                        10,   25,   160, 150, -0.15, 0.15);
	histo2D["El_RecoThetaVsGenPt"]       = new TH2F("El_RecoThetaVsGenPt",      "Polar angle #theta_{rec} versus transverse momentum p_{T,gen} for electron",     100,    0,   200,  60,     0, 3.15);
	histo2D["El_DiffThetaVsGenPt"]       = new TH2F("El_DiffThetaVsGenPt",      "#theta difference (gen-reco) versus p_{T,gen} for electron",                      10,   30,   150, 100, -0.04, 0.04);
	histo2D["Mu_RecoThetaVsGenInvPt"]    = new TH2F("Mu_RecoThetaVsGenInvPt",   "Polar angle #theta_{rec} versus #frac{1}{p_{T,gen}} for muon",                   100,    0,  0.05,  60,     0, 3.15);
        histo2D["Mu_RecoThetaVsGenPt"]       = new TH2F("Mu_RecoThetaVsGenPt",      "Polar angle #theta_{rec} versus transverse momentum p_{T,gen} for muon",         150,    0,   200,  60,     0, 3.15);
	histo2D["Mu_DiffThetaVsGenInvPt"]    = new TH2F("Mu_DiffThetaVsGenInvPt",   "#theta difference (gen-reco) versus #frac{1}{p_{T,gen}} for muon",                10,0.005,  0.04, 100, -0.04, 0.04);
        histo2D["Mu_DiffThetaVsGenPt"]       = new TH2F("Mu_DiffThetaVsGenPt",      "#theta difference (gen-reco) versus p_{T,gen} for muon",                          10,   30,   150, 100,  -0.1,  0.1);
	
	histo2D["Light_RecoPhiVsGenPhi"]     = new TH2F("Light_RecoPhiVsGenPhi",    "Azimuthal angle distribution of light quarks (reco vs gen)",                       60,  -3.2,   3.2,  60,  -3.2,  3.2);
	histo2D["Light_DiffPhiVsGenPhi"]     = new TH2F("Light_DiffPhiVsGenPhi",    "#phi difference (gen-reco) versus #phi_{gen} for light quarks",                    10,  -3.2,   3.2, 100,  -0.2,  0.2);
        histo2D["Light_DiffPhiVsGenPhi_All"] = new TH2F("Light_DiffPhiVsGenPhi_All","#phi difference (gen-reco) versus #phi_{gen} for light quarks",                    10,  -3.2,   3.2, 120,  -6.2,  6.2);
	histo2D["BJet_RecoPhiVsGenPhi"]      = new TH2F("BJet_RecoPhiVsGenPhi",     "Azimuthal angle distribution of b-jets (reco vs gen)",                             60,  -3.2,   3.2,  60,  -3.2,  3.2);
	histo2D["BJet_DiffPhiVsGenPhi"]      = new TH2F("BJet_DiffPhiVsGenPhi",     "#phi difference (gen-reco) versus #phi_{gen} for b-jets",                          10,  -3.2,   3.2, 100,  -0.2,  0.2);
        histo2D["BJet_DiffPhiVsGenPhi_All"]  = new TH2F("BJet_DiffPhiVsGenPhi_All", "#phi difference (gen-reco) versus #phi_{gen} for b-jets",                          10,  -3.2,   3.2, 120,  -6.2,  6.2);
	histo2D["El_RecoPhiVsGenPhi"]        = new TH2F("El_RecoPhiVsGenPhi",       "Azimuthal angle distribution of electron (reco vs gen)", 			        60,  -3.2,   3.2,  60,  -3.2,  3.2);
	histo2D["El_DiffPhiVsGenPhi"]        = new TH2F("El_DiffPhiVsGenPhi",       "#phi difference (gen-reco) versus #phi_{gen} for electron", 		        10,  -3.2,   3.2,  75, -0.15, 0.15);
        histo2D["El_DiffPhiVsGenPhi_All"]    = new TH2F("El_DiffPhiVsGenPhi_All",   "#phi difference (gen-reco) versus #phi_{gen} for electron", 		        10,  -3.2,   3.2,  80,  -3.2,  3.2);
	histo2D["Mu_RecoPhiVsGenPhi"]        = new TH2F("Mu_RecoPhiVsGenPhi",       "Azimuthal angle distribution of muon (reco vs gen)", 			        60,  -3.2,   3.2,  60,  -3.2,  3.2);
	histo2D["Mu_DiffPhiVsGenPhi"]        = new TH2F("Mu_DiffPhiVsGenPhi",       "#phi difference (gen-reco) versus #phi_{gen} for muon", 			        10,  -3.2,   3.2,  75,  -0.2,  0.2);
        histo2D["Mu_DiffPhiVsGenPhi_All"]    = new TH2F("Mu_DiffPhiVsGenPhi_All",   "#phi difference (gen-reco) versus #phi_{gen} for muon", 			        10,  -3.2,   3.2,  80,  -3.2,  3.2);
	histo2D["Light_RecoPhiVsGenPt"]      = new TH2F("Light_RecoPhiVsGenPt",     "Azimuthal angle #phi_{rec} versus transverse momentum p_{T,gen} for light quarks",150,     0,   300,  60,  -3.2,  3.2);
	histo2D["Light_DiffPhiVsGenPt"]      = new TH2F("Light_DiffPhiVsGenPt",     "#phi difference (gen-reco) versus p_{T,gen} for light quarks",                     10,    25,   160, 100, -0.15, 0.15);
        histo2D["Light_DiffPhiVsGenPt_All"]  = new TH2F("Light_DiffPhiVsGenPt_All", "#phi difference (gen-reco) versus p_{T,gen} for light quarks",                     10,     0,   250, 120,  -6.2,  6.2);
	histo2D["BJet_RecoPhiVsGenPt"]       = new TH2F("BJet_RecoPhiVsGenPt",      "Azimuthal angle #phi_{rec} versus transverse momentum p_{T,gen} for b-jets",      150,     0,   300,  60,  -3.2,  3.2);
	histo2D["BJet_DiffPhiVsGenPt"]       = new TH2F("BJet_DiffPhiVsGenPt",      "#phi difference (gen-reco) versus p_{T,gen} for b-jets",                           10,    25,   160, 100, -0.15, 0.15);
        histo2D["BJet_DiffPhiVsGenPt_All"]   = new TH2F("BJet_DiffPhiVsGenPt_All",  "#phi difference (gen-reco) versus p_{T,gen} for b-jets",                           10,     0,   250, 120,  -6.2,  6.2);
	histo2D["El_RecoPhiVsGenPt"]         = new TH2F("El_RecoPhiVsGenPt",        "Azimuthal angle #phi_{rec} versus transverse momentum p_{T,gen} for electron",    150,     0,   250,  60,  -3.2,  3.2);
	histo2D["El_DiffPhiVsGenPt"]         = new TH2F("El_DiffPhiVsGenPt",        "#phi difference (gen-reco) versus p_{T,gen} for electron",                         10,    30,   150,  80, -0.04, 0.04);
        histo2D["El_DiffPhiVsGenPt_All"]     = new TH2F("El_DiffPhiVsGenPt_All",    "#phi difference (gen-reco) versus p_{T,gen} for electron",                         17,     0,   250, 120,  -6.2,  6.2);
	histo2D["Mu_RecoPhiVsGenInvPt"]      = new TH2F("Mu_RecoPhiVsGenInvPt",     "Azimuthal angle #phi_{rec} versus #frac{1}{p_{T,gen}} for muon",                  100,     0,  0.05,  60,  -3.2,  3.2);
        histo2D["Mu_RecoPhiVsGenPt"]         = new TH2F("Mu_RecoPhiVsGenPt",        "Azimuthal angle #phi_{rec} versus transverse momentum p_{T,gen} for muon",        150,     0,   200,  60,  -3.2,  3.2);
	histo2D["Mu_DiffPhiVsGenInvPt"]      = new TH2F("Mu_DiffPhiVsGenInvPt",     "#phi difference (gen-reco) versus #frac{1}{p_{T,gen}} for muon",                   10, 0.005,  0.04, 100, -0.08, 0.08);
        histo2D["Mu_DiffPhiVsGenInvPt_All"]  = new TH2F("Mu_DiffPhiVsGenInvPt_All", "#phi difference (gen-reco) versus #frac{1}{p_{T,gen}} for muon",                   10,     0,  0.05, 120,  -6.2,  6.2);
        histo2D["Mu_DiffPhiVsGenPt"]         = new TH2F("Mu_DiffPhiVsGenPt",        "#phi difference (gen-reco) versus p_{T,gen} for muon",                             10,    30,   150, 100,  -0.2,  0.2);
        histo2D["Mu_DiffPhiVsGenPt_All"]     = new TH2F("Mu_DiffPhiVsGenPt_All",    "#phi difference (gen-reco) versus p_{T,gen} for muon",                             14,     0,   200, 100,  -6.2,  6.2);
}

void TFCreation::FillHistograms(TRootMCParticle* hadrWJet1, TRootMCParticle* hadrWJet2, TRootMCParticle* hadrBJet, TRootMCParticle* leptBJet, TRootMCParticle* lepton, TRootJet* selHadrWJet1, TRootJet* selHadrWJet2, TRootJet* selHadrBJet, TRootJet* selLeptBJet, TLorentzVector* selLepton, bool isSemiMu, bool isSemiEl){

	//Should use Pt information in stead of E!
	// --> Both concepts are identical in the case of CaloJets, but not in the case of PF
	// --> PF uses massive objects to construct particles!

	histo1D["DeltaR_TFClass_Light1"]->Fill( sqrt( pow(hadrWJet1->Phi() - selHadrWJet1->Phi(),2) + pow(hadrWJet1->Eta() - selHadrWJet1->Eta(),2) ) );
	histo1D["DeltaR_TFClass_Light2"]->Fill( sqrt( pow(hadrWJet2->Phi() - selHadrWJet2->Phi(),2) + pow(hadrWJet2->Eta() - selHadrWJet2->Eta(),2) ) );
        histo1D["DeltaR_TFClass_HadrB"]->Fill( sqrt( pow(hadrBJet->Phi() - selHadrBJet->Phi(),2) + pow(hadrBJet->Eta() - selHadrBJet->Eta(),2) ) );
        histo1D["DeltaR_TFClass_LeptB"]->Fill( sqrt( pow(leptBJet->Phi() - selLeptBJet->Phi(),2) + pow(leptBJet->Eta() - selLeptBJet->Eta(),2) ) );
	if(isSemiMu) histo1D["DeltaR_TFClass_Mu"]->Fill( sqrt( pow(lepton->Phi() - selLepton->Phi(),2) + pow(lepton->Eta() - selLepton->Eta(),2) ) );
	if(isSemiEl) histo1D["DeltaR_TFClass_El"]->Fill( sqrt( pow(lepton->Phi() - selLepton->Phi(),2) + pow(lepton->Eta() - selLepton->Eta(),2) ) );
	
	histo2D["Light_RecoPtVsGenPt"]->Fill(      hadrWJet1->Pt(),   selHadrWJet1->Pt()    );
	histo2D["Light_RecoThetaVsGenTheta"]->Fill(hadrWJet1->Theta(),selHadrWJet1->Theta() );
	histo2D["Light_RecoThetaVsGenPt"]->Fill(   hadrWJet1->Pt(),   selHadrWJet1->Theta() );
	histo2D["Light_RecoPhiVsGenPhi"]->Fill(    hadrWJet1->Phi(),  selHadrWJet1->Phi()   );
	histo2D["Light_RecoPhiVsGenPt"]->Fill(     hadrWJet1->Pt(),   selHadrWJet1->Phi()   );

	histo2D["Light_DiffPtVsGenPt"]->Fill(      hadrWJet1->Pt(),   hadrWJet1->Pt()    - selHadrWJet1->Pt()    );
	histo2D["Light_DiffThetaVsGenTheta"]->Fill(hadrWJet1->Theta(),hadrWJet1->Theta() - selHadrWJet1->Theta() );
	histo2D["Light_DiffThetaVsGenPt"]->Fill(   hadrWJet1->Pt(),   hadrWJet1->Theta() - selHadrWJet1->Theta() );
	histo2D["Light_DiffPhiVsGenPhi"]->Fill(    hadrWJet1->Phi(),  hadrWJet1->Phi()   - selHadrWJet1->Phi()   );
        histo2D["Light_DiffPhiVsGenPhi_All"]->Fill(hadrWJet1->Phi(),  hadrWJet1->Phi()   - selHadrWJet1->Phi()   );
	histo2D["Light_DiffPhiVsGenPt"]->Fill(     hadrWJet1->Pt(),   hadrWJet1->Phi()   - selHadrWJet1->Phi()   );
        histo2D["Light_DiffPhiVsGenPt_All"]->Fill( hadrWJet1->Pt(),   hadrWJet1->Phi()   - selHadrWJet1->Phi()   );

	histo2D["Light_RecoPtVsGenPt"]->Fill(      hadrWJet2->Pt(),   selHadrWJet2->Pt()    );
	histo2D["Light_RecoThetaVsGenTheta"]->Fill(hadrWJet2->Theta(),selHadrWJet2->Theta() );
	histo2D["Light_RecoThetaVsGenPt"]->Fill(   hadrWJet2->Pt(),   selHadrWJet2->Theta() );
	histo2D["Light_RecoPhiVsGenPhi"]->Fill(    hadrWJet2->Phi(),  selHadrWJet2->Phi()   );
	histo2D["Light_RecoPhiVsGenPt"]->Fill(     hadrWJet2->Pt(),   selHadrWJet2->Phi()   );

	float DeltaPhi_nonbjet2 = ROOT::Math::VectorUtil::DeltaPhi((TLorentzVector)*hadrWJet2,(TLorentzVector)*selHadrWJet2);
	histo2D["Light_DiffPtVsGenPt"]->Fill(      hadrWJet2->Pt(),   hadrWJet2->Pt()    - selHadrWJet2->Pt()    );
	histo2D["Light_DiffThetaVsGenTheta"]->Fill(hadrWJet2->Theta(),hadrWJet2->Theta() - selHadrWJet2->Theta() );
	histo2D["Light_DiffThetaVsGenPt"]->Fill(   hadrWJet2->Pt(),   hadrWJet2->Theta() - selHadrWJet2->Theta() );
	histo2D["Light_DiffPhiVsGenPhi"]->Fill(    hadrWJet2->Phi(),  hadrWJet2->Phi()   - selHadrWJet2->Phi()   );
	histo2D["Light_DiffPhiVsGenPhi_All"]->Fill(hadrWJet2->Phi(),  hadrWJet2->Phi()   - selHadrWJet2->Phi()   );
	histo2D["Light_DiffPhiVsGenPt"]->Fill(     hadrWJet2->Pt(),   hadrWJet2->Phi()   - selHadrWJet2->Phi()   );
	histo2D["Light_DiffPhiVsGenPt_All"]->Fill( hadrWJet2->Pt(),   hadrWJet2->Phi()   - selHadrWJet2->Phi()   );

	histo2D["BJet_RecoPtVsGenPt"]->Fill(      hadrBJet->Pt(),    selHadrBJet->Pt()    );
	histo2D["BJet_RecoThetaVsGenTheta"]->Fill(hadrBJet->Theta(), selHadrBJet->Theta() );
	histo2D["BJet_RecoThetaVsGenPt"]->Fill(   hadrBJet->Pt(),    selHadrBJet->Theta() );
	histo2D["BJet_RecoPhiVsGenPhi"]->Fill(    hadrBJet->Phi(),   selHadrBJet->Phi());
	histo2D["BJet_RecoPhiVsGenPt"]->Fill(     hadrBJet->Pt(),    selHadrBJet->Phi());

	float DeltaPhi_bjet1 = ROOT::Math::VectorUtil::DeltaPhi((TLorentzVector)*hadrBJet,(TLorentzVector)*selLeptBJet);				
	histo2D["BJet_DiffPtVsGenPt"]->Fill(      hadrBJet->Pt(),   hadrBJet->Pt()    - selHadrBJet->Pt()    );
	histo2D["BJet_DiffThetaVsGenTheta"]->Fill(hadrBJet->Theta(),hadrBJet->Theta() - selHadrBJet->Theta() );
	histo2D["BJet_DiffThetaVsGenPt"]->Fill(   hadrBJet->Pt(),   hadrBJet->Theta() - selHadrBJet->Theta() );
	histo2D["BJet_DiffPhiVsGenPhi"]->Fill(    hadrBJet->Phi(),  hadrBJet->Phi()   - selHadrBJet->Phi()   );
	histo2D["BJet_DiffPhiVsGenPhi_All"]->Fill(hadrBJet->Phi(),  hadrBJet->Phi()   - selHadrBJet->Phi()   );
	histo2D["BJet_DiffPhiVsGenPt"]->Fill(     hadrBJet->Pt(),   hadrBJet->Phi()   - selHadrBJet->Phi()   );
	histo2D["BJet_DiffPhiVsGenPt_All"]->Fill( hadrBJet->Pt(),   hadrBJet->Phi()   - selHadrBJet->Phi()   );

	histo2D["BJet_RecoPtVsGenPt"]->Fill(      leptBJet->Pt(),    selLeptBJet->Pt()    );
	histo2D["BJet_RecoThetaVsGenTheta"]->Fill(leptBJet->Theta(), selLeptBJet->Theta() );
	histo2D["BJet_RecoThetaVsGenPt"]->Fill(   leptBJet->Pt(),    selLeptBJet->Theta() );
	histo2D["BJet_RecoPhiVsGenPhi"]->Fill(    leptBJet->Phi(),   selLeptBJet->Phi()   );
	histo2D["BJet_RecoPhiVsGenPt"]->Fill(     leptBJet->Pt(),    selLeptBJet->Phi()   );

	histo2D["BJet_DiffPtVsGenPt"]->Fill(      leptBJet->Pt(),   leptBJet->Pt()    - selLeptBJet->Pt()    );
	histo2D["BJet_DiffThetaVsGenTheta"]->Fill(leptBJet->Theta(),leptBJet->Theta() - selLeptBJet->Theta() );
	histo2D["BJet_DiffThetaVsGenPt"]->Fill(   leptBJet->Pt(),   leptBJet->Theta() - selLeptBJet->Theta() );
	histo2D["BJet_DiffPhiVsGenPhi"]->Fill(    leptBJet->Phi(),  leptBJet->Phi()   - selLeptBJet->Phi()   );
	histo2D["BJet_DiffPhiVsGenPhi_All"]->Fill(leptBJet->Phi(),  leptBJet->Phi()   - selLeptBJet->Phi()   );
	histo2D["BJet_DiffPhiVsGenPt"]->Fill(     leptBJet->Pt(),   leptBJet->Phi()   - selLeptBJet->Phi()   );
	histo2D["BJet_DiffPhiVsGenPt_All"]->Fill( leptBJet->Pt(),   leptBJet->Phi()   - selLeptBJet->Phi()   );

	if(isSemiEl){
		histo2D["El_RecoPtVsGenPt"]->Fill(      lepton->Pt(),   selLepton->Pt()    );
		histo2D["El_RecoThetaVsGenTheta"]->Fill(lepton->Theta(),selLepton->Theta() );
	    	histo2D["El_RecoThetaVsGenPt"]->Fill(   lepton->Pt(),   selLepton->Theta() );
	    	histo2D["El_RecoPhiVsGenPt"]->Fill(     lepton->Pt(),   selLepton->Phi()   );
		histo2D["El_RecoPhiVsGenPhi"]->Fill(    lepton->Phi(),  selLepton->Phi()   );

	    	histo2D["El_DiffPtVsGenPt"]->Fill(      lepton->Pt(),   lepton->Pt()    - selLepton->Pt()    );
	    	histo2D["El_DiffThetaVsGenTheta"]->Fill(lepton->Theta(),lepton->Theta() - selLepton->Theta() );
	    	histo2D["El_DiffThetaVsGenPt"]->Fill(   lepton->Pt(),   lepton->Theta() - selLepton->Theta() );
	    	histo2D["El_DiffPhiVsGenPhi"]->Fill(    lepton->Phi(),  lepton->Phi()   - selLepton->Phi()   );
	    	histo2D["El_DiffPhiVsGenPt"]->Fill(     lepton->Pt(),   lepton->Phi()   - selLepton->Phi()   );
	    	histo2D["El_DiffPhiVsGenPhi_All"]->Fill(lepton->Phi(),  lepton->Phi()   - selLepton->Phi()   );
	    	histo2D["El_DiffPhiVsGenPt_All"]->Fill( lepton->Pt(),   lepton->Phi()   - selLepton->Phi()   );
	}
	if(isSemiMu){
		float InvPtgenMu = 1./lepton->Pt();
	    	float InvPtrecMu = 1./selLepton->Pt();
	    	histo2D["Mu_RecoInvPtVsGenInvPt"]->Fill(InvPtgenMu,     InvPtrecMu         );
        	histo2D["Mu_RecoPtVsGenPt"]->Fill(      lepton->Pt(),   selLepton->Pt()    );
		histo2D["Mu_RecoThetaVsGenTheta"]->Fill(lepton->Theta(),selLepton->Theta() );
	    	histo2D["Mu_RecoThetaVsGenInvPt"]->Fill(InvPtgenMu,     selLepton->Theta() );
		histo2D["Mu_RecoThetaVsGenPt"]->Fill(   lepton->Pt(),   selLepton->Theta() );
		histo2D["Mu_RecoPhiVsGenPhi"]->Fill(    lepton->Phi(),  selLepton->Phi()   );
	    	histo2D["Mu_RecoPhiVsGenInvPt"]->Fill(  InvPtgenMu,     selLepton->Phi()   );
		histo2D["Mu_RecoPhiVsGenPt"]->Fill(     lepton->Pt(),   selLepton->Phi()   );

	    	histo2D["Mu_DiffInvPtVsGenInvPt"]->Fill(    InvPtgenMu,     InvPtgenMu      - InvPtrecMu         );
		histo2D["Mu_DiffInvPtVsGenInvPt_All"]->Fill(InvPtgenMu,     InvPtgenMu      - InvPtrecMu         );
		histo2D["Mu_DiffPtVsGenPt"]->Fill(          lepton->Pt(),   lepton->Pt()    - selLepton->Pt()    );
		histo2D["Mu_DiffThetaVsGenTheta"]->Fill(    lepton->Theta(),lepton->Theta() - selLepton->Theta() );
	    	histo2D["Mu_DiffThetaVsGenInvPt"]->Fill(    InvPtgenMu,     lepton->Theta() - selLepton->Theta() );
		histo2D["Mu_DiffThetaVsGenPt"]->Fill(       lepton->Pt(),   lepton->Theta() - selLepton->Theta() );
	    	histo2D["Mu_DiffPhiVsGenPhi"]->Fill(        lepton->Phi(),  lepton->Phi()   - selLepton->Phi()   );
	    	histo2D["Mu_DiffPhiVsGenPhi_All"]->Fill(    lepton->Phi(),  lepton->Phi()   - selLepton->Phi()   );
	    	histo2D["Mu_DiffPhiVsGenInvPt"]->Fill(      InvPtgenMu,     lepton->Phi()   - selLepton->Phi()   );
	    	histo2D["Mu_DiffPhiVsGenInvPt_All"]->Fill(  InvPtgenMu,     lepton->Phi()   - selLepton->Phi()   );
		histo2D["Mu_DiffPhiVsGenPt"]->Fill(         lepton->Pt(),   lepton->Phi()   - selLepton->Phi()   );
		histo2D["Mu_DiffPhiVsGenPt_All"]->Fill(     lepton->Pt(),   lepton->Phi()   - selLepton->Phi()   );
	}
}

void TFCreation::CalculateTFFromFile(TH2F* fitHisto, bool useStartValues, int histoNr, bool useROOTClass, bool useStartArray, float startValues[], TFile* file){

    TDirectory* histoFitDir = file->mkdir(fitHisto->GetName());
    histoFitDir->cd();

    //Set parameter names
    const char* parnames[6]={"a1","a2","a3","a4","a5","a6"};
    const int npar = doubleGaussianFit->GetNpar();
    for(int ii = 0; ii < npar; ii++) doubleGaussianFit->SetParName(ii,parnames[ii]);    

    //Set range of fits
    doubleGaussianFit->SetRange( fitHisto->GetYaxis()->GetXmin(), fitHisto->GetYaxis()->GetXmin() );
    caloEnergyFit->SetRange( fitHisto->GetXaxis()->GetXmin(), fitHisto->GetXaxis()->GetXmax() );						

    //Initialize the start values if asked
    startValuesArray = startValues;
    if(useStartValues) SetStartValuesDoubleGaussian(histoNr, useStartArray);         //Can only be done after that doubleGaussianFit is initialized!

    //Choose the correct fit method:
    hlist = new TH1D*[npar];
    TObjArray aSlices;
    if(useROOTClass){
	fitHisto->FitSlicesY(doubleGaussianFit, 0, -1, 0, "", &aSlices);
	for(int ipar = 0; ipar <= npar; ipar++) hlist[ipar] = (TH1D*) aSlices[ipar];
    }
    else
	FitSliceClassCode(fitHisto, npar, parnames);

    //////////////////////////////////////////////////////////////////////////////////////////////
    //   Now histogram with all parameters needs to be fitted with Calorimeter Energy formula   //
    //////////////////////////////////////////////////////////////////////////////////////////////
    for( int ipar = 0; ipar < npar; ipar++ ){

        for(int ii = 0; ii < 3; ii++) caloEnergyFit->SetParName(ii, ( string(parnames[ipar])+tostr(ii)).c_str() ); //Name here since different for each doubleGaussian parameter!
	caloEnergyFit->SetName( (string(fitHisto->GetName())+"_"+parnames[ipar]+"_Fit").c_str() );
	hlist[ipar]->SetName( (string(fitHisto->GetName())+"_"+parnames[ipar]+"_PointsAndFit").c_str() );
	hlist[ipar]->Fit(caloEnergyFit);
        AllCaloEnergyFits[ipar] = *caloEnergyFit;       //caloEnergyFit is a pointer, but each member of the array should point to the corresponding value of the TF1!
	hlist[ipar]->Write();            
    }
    hlist[npar]->Write();
    delete [] hlist;
}

void TFCreation::CalculateTF(bool drawHistos, bool doFits, bool useROOTClass, bool useStartValues){
	TFile* file = new TFile("PlotsForTransferFunctions.root","RECREATE");
	file->cd();

	if(drawHistos == true) WritePlots(file);
	if(doFits == true){

	  //Set parameter names
	  const char* parnames[6]={"a1","a2","a3","a4","a5","a6"};
	  const int npar = doubleGaussianFit->GetNpar();
	  for(int ii = 0; ii < npar; ii++) doubleGaussianFit->SetParName(ii,parnames[ii]);

	  ///////////////////////////////////////////
	  //  Choose the correct histogram to fit  //
	  ///////////////////////////////////////////
          TH2F* histoForFit;
	  for (unsigned int f=0; f<12;f++) {				
				
	    switch(f){
	      case 0:
		histoForFit=histo2D["Light_DiffPtVsGenPt"]; break;
	      case 1:
		histoForFit=histo2D["BJet_DiffPtVsGenPt"]; break;
	      case 2:
		histoForFit=histo2D["El_DiffPtVsGenPt"]; break;
	      case 3:
		histoForFit=histo2D["Mu_DiffInvPtVsGenInvPt"]; break;
	      case 4:
		histoForFit=histo2D["Light_DiffThetaVsGenPt"]; break;
	      case 5:
		histoForFit=histo2D["BJet_DiffThetaVsGenPt"]; break;
	      case 6:
		histoForFit=histo2D["El_DiffThetaVsGenPt"]; break;
	      case 7:
		histoForFit=histo2D["Mu_DiffThetaVsGenInvPt"]; break;
	      case 8:
		histoForFit=histo2D["Light_DiffPhiVsGenPt"]; break;
	      case 9:
		histoForFit=histo2D["BJet_DiffPhiVsGenPt"]; break;
	      case 10:
		histoForFit=histo2D["El_DiffPhiVsGenPt"]; break;
	      case 11:
		histoForFit=histo2D["Mu_DiffPhiVsGenInvPt"]; break;
	    }					

            TDirectory* histoFitDir = file->mkdir(histoForFit->GetName());
            histoFitDir->cd();
 
	    hlist = new TH1D*[npar];
	    if(useStartValues)
		SetStartValuesDoubleGaussian(f, false);   //false means that normal start values are being used!

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
                for(int jj = 0; jj < 3; jj++) caloEnergyFit->SetParName(jj, ( string(parnames[ipar])+tostr(jj)).c_str() );
		caloEnergyFit->SetName( (string(histoForFit->GetName())+"_"+parnames[ipar]+"_Fit").c_str() );
		hlist[ipar]->SetName( (string(histoForFit->GetName())+"_"+parnames[ipar]+"_PointsAndFit").c_str() );

		hlist[ipar]->Fit(caloEnergyFit);
		hlist[ipar]->Write();            
	    }
	    hlist[npar]->Write();
						
	  }//Loop over f						
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

	    hlist[ipar] = new TH1D( (string(histoForFit->GetName())+"_"+parnames[ipar]).c_str(), (string(histoForFit->GetName())+" : Fitted value of "+parnames[ipar]).c_str(), histoForFit->GetXaxis()->GetNbins(), histoForFit->GetXaxis()->GetXmin(), histoForFit->GetXaxis()->GetXmax());
	    hlist[ipar]->GetXaxis()->SetTitle(histoForFit->GetXaxis()->GetTitle());
	}
	hlist[npar] = new TH1D( (string(histoForFit->GetName())+"_chi2").c_str(), (string(histoForFit->GetName())+": #chi^{2} distribution for "+string(doubleGaussianFit->GetExpFormula())).c_str(), histoForFit->GetXaxis()->GetNbins(), histoForFit->GetXaxis()->GetXmin(), histoForFit->GetXaxis()->GetXmax() );

	//Loop on all bins in X, generate a projection along Y and fit each bin separately!
	int cut = 0; // require a minimum number of bins in the slice to be filled --> Should this ever be larger than 0 ??
	int nbins = histoForFit->GetXaxis()->GetNbins();
	cout << "\n ** Looking at histogram : " << histoForFit->GetName() << "                               ******************************* " << endl;
	for(int bin=1;bin <= nbins+1;bin ++) {
	    cout << "   --  Looking at bin : " << bin << endl;
	    string projection_title = string(histoForFit->GetName())+"_sliceYbin"+tostr(bin);

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
	    delete hp;
	}//loop over bins!
}	

void TFCreation::SetStartValuesDoubleGaussian(int whichHisto, bool useStartArray){

    if(useStartArray == true){
        for(int ii = 0; ii < 6; ii++)
            doubleGaussianFit->SetParameter(ii, startValuesArray[ii]);
    }
    else{
	if(whichHisto==0 || whichHisto==1 || whichHisto == 2){ // for E transfer function of JETS (and elec -- added as test ...)
	    float StartValues[] = {-8,18,63,0,8.6,4.1};        //First three values are for the first broad gaussian (central, sigma and constant value respectively)
							       //Second three values are the same for the second narrow gaussian
	    for(int ii = 0; ii < 6; ii++)
		doubleGaussianFit->SetParameter(ii,StartValues[ii]);
	}
	else if (whichHisto==4 || whichHisto==5 || whichHisto==8 || whichHisto==9 || whichHisto == 6 || whichHisto == 10) { //for theta and phi transfer functions of JETS (and elec)
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
} 

void TFCreation::WriteTF(TH2F* fitHisto, ofstream &myTFs, ofstream &myTFCard, TFile* file){

    const int NrConsideredPlots = 9, NrConsideredPars = 6, NrConsideredCaloPars = 3;
    string ParamName[NrConsideredPars] = {"Mean broad gaussian", "Width broad gaussian","Constant broad gaussian","Mean narrow gaussian","Width narrow gaussian","Constant narrow gaussian"};

    for(int ipar = 0; ipar < NrConsideredPars; ipar++){

	for(int icalopar = 0; icalopar < NrConsideredCaloPars; icalopar++){
	    if(icalopar == 0) myTFs<<ParamName[ipar]<<" & $a_{" <<ipar <<icalopar <<"}$ = "<<AllCaloEnergyFits[ipar].GetParameter(icalopar)<<"$\\pm$"<<AllCaloEnergyFits[ipar].GetParError(icalopar);
	    else              myTFs<<                 " & $a_{" <<ipar <<icalopar <<"}$ = "<<AllCaloEnergyFits[ipar].GetParameter(icalopar)<<"$\\pm$"<<AllCaloEnergyFits[ipar].GetParError(icalopar);

            myTFCard<< (ipar*NrConsideredCaloPars)+icalopar+1 << "     " << AllCaloEnergyFits[ipar].GetParameter(icalopar)<< "     # " << ParamName[ipar] << endl;
	}
	myTFs << "\\\\" << endl;
    }
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
	}
	TDirectory* th2dir = outfile->mkdir("2D_histograms_graphs");
	th2dir->cd();
	for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++){    
	      TH2F *temp = it->second;
	      temp->Write();
	}
	outfile->cd(); //Step out from 2D_histograms_graphs directory!
}
