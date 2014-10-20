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
    delete caloEnergyFit;
    delete doubleGaussianFit;
}

void TFCreation::InitializeVariables(){
    histo1D["DeltaR_TFClass_Light1"] = new TH1F("DeltaR_TFClass_Light1","DeltaR_TFClass_Light1",200,0,1);
    histo1D["DeltaR_TFClass_Light2"] = new TH1F("DeltaR_TFClass_Light2","DeltaR_TFClass_Light2",200,0,1);
    histo1D["DeltaR_TFClass_HadrB"]  = new TH1F("DeltaR_TFClass_HadrB","DeltaR_TFClass_HadrB",200,0,1);
    histo1D["DeltaR_TFClass_LeptB"]  = new TH1F("DeltaR_TFClass_LeptB","DeltaR_TFClass_LeptB",200,0,1);
    histo1D["DeltaR_TFClass_Mu"]     = new TH1F("DeltaR_TFClass_Mu","DeltaR_TFClass_Mu",200,0,1);
    histo1D["DeltaR_TFClass_El"]     = new TH1F("DeltaR_TFClass_El","DeltaR_TFClass_El",200,0,1);

    histo2D["Light_GenPtVsRecoPt"]    = new TH2F("Light_GenPtVsRecoPt",   "Transverse momentum of light quarks (gen vs reco)",                        150,     0,  300, 150,      0,  300);
    histo2D["Light_DiffPtVsRecoPt"]    = new TH2F("Light_DiffPtVsRecoPt",   "p_{T} difference (gen-reco) versus p_{T,reco} for light quarks",             10,    25,  160,  80,    -40,   40);
    histo2D["BJet_GenPtVsRecoPt"]     = new TH2F("BJet_GenPtVsRecoPt",    "Transverse momentum of b-jets (reco vs gen level)",                        150,     0,  300, 150,      0,  300);
    histo2D["BJet_DiffPtVsRecoPt"]     = new TH2F("BJet_DiffPtVsRecoPt",    "p_{T} difference (gen-reco) versus p_{T,reco} for b-jets",                   10,    25,  160, 100,    -50,   50);
    histo2D["El_GenPtVsRecoPt"]       = new TH2F("El_GenPtVsRecoPt",      "Transverse momentum of electron (gen vs reco)",                            100,     0,  200, 100,      0,  200);
    histo2D["El_DiffPtVsRecoPt"]       = new TH2F("El_DiffPtVsRecoPt",      "p_{T} difference (gen-reco) versus p_{T,reco} for electron",                 10,    30,  150,  75,    -10,   10);
    histo2D["Mu_GenInvPtVsRecoInvPt"] = new TH2F("Mu_GenInvPtVsRecoInvPt","Inverse of transverse momentum of muon (gen vs reco)",                     100,     0, 0.05, 100,      0, 0.05);
    histo2D["Mu_DiffInvPtVsRecoInvPt"] = new TH2F("Mu_DiffInvPtVsRecoInvPt","#frac{1}{p_{T}} difference (gen-reco) versus #frac{1}{p_{T,reco}} for muon", 10, 0.005, 0.04,  80,  -0.01, 0.01);
    histo2D["Mu_DiffInvPtVsRecoInvPt_All"] =new TH2F("Mu_DiffInvPtVsRecoInvPt_All","#frac{1}{p_{T}} difference (gen-reco) versus #frac{1}{p_{T,reco}} for muon", 20, 0.005, 0.04,  180,  -0.1, 0.1);
    histo2D["Mu_GenPtVsRecoPt"]       = new TH2F("Mu_GenPtVsRecoPt",      "Transverse momentum of muon (gen vs reco)",                                150,     0,  200, 150,      0,  200);
    histo2D["Mu_DiffPtVsRecoPt"]       = new TH2F("Mu_DiffPtVsRecoPt",      "p_{T} difference (gen-reco) versus p_{T,reco} for muon",                     10,    30,  150, 100,    -10,   10);
	
    histo2D["Light_GenThetaVsRecoTheta"] = new TH2F("Light_GenThetaVsRecoTheta","Polar angle distribution of light quarks (gen vs reco)",                          60,    0,  3.15,  60,     0, 3.15);
    histo2D["Light_DiffThetaVsRecoTheta"] = new TH2F("Light_DiffThetaVsRecoTheta","#theta difference (gen-reco) versus #theta_{reco} for light quarks",               10,  0.1,   3.1, 100, -0.15, 0.15);
    histo2D["BJet_GenThetaVsRecoTheta"]  = new TH2F("BJet_GenThetaVsRecoTheta", "Polar angle distribution of b-jets (gen vs reco)",                                60,    0,  3.15,  60,     0, 3.15);
    histo2D["BJet_DiffThetaVsRecoTheta"]  = new TH2F("BJet_DiffThetaVsRecoTheta", "#theta difference (gen-reco) versus #theta_{reco} for b-jets",                     10,  0.1,   3.1, 100, -0.15, 0.15);
    histo2D["El_GenThetaVsRecoTheta"]    = new TH2F("El_GenThetaVsRecoTheta",   "Polar angle distribution of electron (gen vs reco)",                              60,    0,  3.15,  60,     0, 3.15);
    histo2D["El_DiffThetaVsRecoTheta"]    = new TH2F("El_DiffThetaVsRecoTheta",   "#theta difference (gen-reco) versus #theta_{reco} for electron",                   10,    0,  3.15, 100, -0.15, 0.15);
    histo2D["Mu_GenThetaVsRecoTheta"]    = new TH2F("Mu_GenThetaVsRecoTheta",   "Polar angle distribution of muon (gen vs reco)",                                  60,    0,  3.15,  60,     0, 3.15);
    histo2D["Mu_DiffThetaVsRecoTheta"]    = new TH2F("Mu_DiffThetaVsRecoTheta",   "#theta difference (gen-reco) versus #theta_{reco} for muon",                       10,    0,  3.15, 100, -0.15, 0.15);
    histo2D["Light_GenThetaVsRecoPt"]    = new TH2F("Light_GenThetaVsRecoPt",   "Polar angle #theta_{rec} versus transverse momentum p_{T,reco} for light quarks", 120,    0,   300,  60,     0, 3.15);
    histo2D["Light_DiffThetaVsRecoPt"]    = new TH2F("Light_DiffThetaVsRecoPt",   "#theta difference (gen-reco) versus p_{T,reco} for light quarks",                  10,   25,   160, 150, -0.15, 0.15);
    histo2D["BJet_GenThetaVsRecoPt"]     = new TH2F("BJet_GenThetaVsRecoPt",    "Polar angle #theta_{rec} versus transverse momentum p_{T,reco} for b-jets",       120,    0,   300,  60,     0, 3.15);
    histo2D["BJet_DiffThetaVsRecoPt"]     = new TH2F("BJet_DiffThetaVsRecoPt",    "#theta difference (gen-reco) versus p_{T,reco} for b-jets",                        10,   25,   160, 150, -0.15, 0.15);
    histo2D["El_GenThetaVsRecoPt"]       = new TH2F("El_GenThetaVsRecoPt",      "Polar angle #theta_{rec} versus transverse momentum p_{T,reco} for electron",     100,    0,   200,  60,     0, 3.15);
    histo2D["El_DiffThetaVsRecoPt"]       = new TH2F("El_DiffThetaVsRecoPt",      "#theta difference (gen-reco) versus p_{T,reco} for electron",                      10,   30,   150, 100, -0.04, 0.04);
    histo2D["Mu_GenThetaVsRecoInvPt"]    = new TH2F("Mu_GenThetaVsRecoInvPt",   "Polar angle #theta_{rec} versus #frac{1}{p_{T,reco}} for muon",                   100,    0,  0.05,  60,     0, 3.15);
    histo2D["Mu_GenThetaVsRecoPt"]       = new TH2F("Mu_GenThetaVsRecoPt",      "Polar angle #theta_{rec} versus transverse momentum p_{T,reco} for muon",         150,    0,   200,  60,     0, 3.15);
    histo2D["Mu_DiffThetaVsRecoInvPt"]    = new TH2F("Mu_DiffThetaVsRecoInvPt",   "#theta difference (gen-reco) versus #frac{1}{p_{T,reco}} for muon",                10,0.005,  0.04, 100, -0.04, 0.04);
    histo2D["Mu_DiffThetaVsRecoPt"]       = new TH2F("Mu_DiffThetaVsRecoPt",      "#theta difference (gen-reco) versus p_{T,reco} for muon",                          10,   30,   150, 100,  -0.1,  0.1);
	
    histo2D["Light_GenPhiVsRecoPhi"]     = new TH2F("Light_GenPhiVsRecoPhi",    "Azimuthal angle distribution of light quarks (gen vs reco)",                       60,  -3.2,   3.2,  60,  -3.2,  3.2);
    histo2D["Light_DiffPhiVsRecoPhi"]     = new TH2F("Light_DiffPhiVsRecoPhi",    "#phi difference (gen-reco) versus #phi_{reco} for light quarks",                    10,  -3.2,   3.2, 100,  -0.2,  0.2);
    histo2D["Light_DiffPhiVsRecoPhi_All"] = new TH2F("Light_DiffPhiVsRecoPhi_All","#phi difference (gen-reco) versus #phi_{reco} for light quarks",                    10,  -3.2,   3.2, 120,  -6.2,  6.2);
    histo2D["BJet_GenPhiVsRecoPhi"]      = new TH2F("BJet_GenPhiVsRecoPhi",     "Azimuthal angle distribution of b-jets (gen vs reco)",                             60,  -3.2,   3.2,  60,  -3.2,  3.2);
    histo2D["BJet_DiffPhiVsRecoPhi"]      = new TH2F("BJet_DiffPhiVsRecoPhi",     "#phi difference (gen-reco) versus #phi_{reco} for b-jets",                          10,  -3.2,   3.2, 100,  -0.2,  0.2);
    histo2D["BJet_DiffPhiVsRecoPhi_All"]  = new TH2F("BJet_DiffPhiVsRecoPhi_All", "#phi difference (gen-reco) versus #phi_{reco} for b-jets",                          10,  -3.2,   3.2, 120,  -6.2,  6.2);
    histo2D["El_GenPhiVsRecoPhi"]        = new TH2F("El_GenPhiVsRecoPhi",       "Azimuthal angle distribution of electron (gen vs reco)", 			        60,  -3.2,   3.2,  60,  -3.2,  3.2);
    histo2D["El_DiffPhiVsRecoPhi"]        = new TH2F("El_DiffPhiVsRecoPhi",       "#phi difference (gen-reco) versus #phi_{reco} for electron", 		        10,  -3.2,   3.2,  75, -0.15, 0.15);
    histo2D["El_DiffPhiVsRecoPhi_All"]    = new TH2F("El_DiffPhiVsRecoPhi_All",   "#phi difference (gen-reco) versus #phi_{reco} for electron", 		        10,  -3.2,   3.2,  80,  -3.2,  3.2);
    histo2D["Mu_GenPhiVsRecoPhi"]        = new TH2F("Mu_GenPhiVsRecoPhi",       "Azimuthal angle distribution of muon (gen vs reco)", 			        60,  -3.2,   3.2,  60,  -3.2,  3.2);
    histo2D["Mu_DiffPhiVsRecoPhi"]        = new TH2F("Mu_DiffPhiVsRecoPhi",       "#phi difference (gen-reco) versus #phi_{reco} for muon", 			        10,  -3.2,   3.2,  75,  -0.2,  0.2);
    histo2D["Mu_DiffPhiVsRecoPhi_All"]    = new TH2F("Mu_DiffPhiVsRecoPhi_All",   "#phi difference (gen-reco) versus #phi_{reco} for muon", 			        10,  -3.2,   3.2,  80,  -3.2,  3.2);
    histo2D["Light_GenPhiVsRecoPt"]      = new TH2F("Light_GenPhiVsRecoPt",     "Azimuthal angle #phi_{rec} versus transverse momentum p_{T,reco} for light quarks",150,     0,   300,  60,  -3.2,  3.2);
    histo2D["Light_DiffPhiVsRecoPt"]      = new TH2F("Light_DiffPhiVsRecoPt",     "#phi difference (gen-reco) versus p_{T,reco} for light quarks",                     10,    25,   160, 100, -0.15, 0.15);
    histo2D["Light_DiffPhiVsRecoPt_All"]  = new TH2F("Light_DiffPhiVsRecoPt_All", "#phi difference (gen-reco) versus p_{T,reco} for light quarks",                     10,     0,   250, 120,  -6.2,  6.2);
    histo2D["BJet_GenPhiVsRecoPt"]       = new TH2F("BJet_GenPhiVsRecoPt",      "Azimuthal angle #phi_{rec} versus transverse momentum p_{T,reco} for b-jets",      150,     0,   300,  60,  -3.2,  3.2);
    histo2D["BJet_DiffPhiVsRecoPt"]       = new TH2F("BJet_DiffPhiVsRecoPt",      "#phi difference (gen-reco) versus p_{T,reco} for b-jets",                           10,    25,   160, 100, -0.15, 0.15);
    histo2D["BJet_DiffPhiVsRecoPt_All"]   = new TH2F("BJet_DiffPhiVsRecoPt_All",  "#phi difference (gen-reco) versus p_{T,reco} for b-jets",                           10,     0,   250, 120,  -6.2,  6.2);
    histo2D["El_GenPhiVsRecoPt"]         = new TH2F("El_GenPhiVsRecoPt",        "Azimuthal angle #phi_{rec} versus transverse momentum p_{T,reco} for electron",    150,     0,   250,  60,  -3.2,  3.2);
    histo2D["El_DiffPhiVsRecoPt"]         = new TH2F("El_DiffPhiVsRecoPt",        "#phi difference (gen-reco) versus p_{T,reco} for electron",                         10,    30,   150,  80, -0.04, 0.04);
    histo2D["El_DiffPhiVsRecoPt_All"]     = new TH2F("El_DiffPhiVsRecoPt_All",    "#phi difference (gen-reco) versus p_{T,reco} for electron",                         17,     0,   250, 120,  -6.2,  6.2);
    histo2D["Mu_GenPhiVsRecoInvPt"]      = new TH2F("Mu_GenPhiVsRecoInvPt",     "Azimuthal angle #phi_{rec} versus #frac{1}{p_{T,reco}} for muon",                  100,     0,  0.05,  60,  -3.2,  3.2);
    histo2D["Mu_GenPhiVsRecoPt"]         = new TH2F("Mu_GenPhiVsRecoPt",        "Azimuthal angle #phi_{rec} versus transverse momentum p_{T,reco} for muon",        150,     0,   200,  60,  -3.2,  3.2);
    histo2D["Mu_DiffPhiVsRecoInvPt"]      = new TH2F("Mu_DiffPhiVsRecoInvPt",     "#phi difference (gen-reco) versus #frac{1}{p_{T,reco}} for muon",                   10, 0.005,  0.04, 100, -0.08, 0.08);
    histo2D["Mu_DiffPhiVsRecoInvPt_All"]  = new TH2F("Mu_DiffPhiVsRecoInvPt_All", "#phi difference (gen-reco) versus #frac{1}{p_{T,reco}} for muon",                   10,     0,  0.05, 120,  -6.2,  6.2);
    histo2D["Mu_DiffPhiVsRecoPt"]         = new TH2F("Mu_DiffPhiVsRecoPt",        "#phi difference (gen-reco) versus p_{T,reco} for muon",                             10,    30,   150, 100,  -0.2,  0.2);
    histo2D["Mu_DiffPhiVsRecoPt_All"]     = new TH2F("Mu_DiffPhiVsRecoPt_All",    "#phi difference (gen-reco) versus p_{T,reco} for muon",                             14,     0,   200, 100,  -6.2,  6.2);
}

void TFCreation::FillHistograms(TLorentzVector* hadrWJet1, TLorentzVector* hadrWJet2, TLorentzVector* hadrBJet, TLorentzVector* leptBJet, TLorentzVector* lepton, TLorentzVector* selHadrWJet1, TLorentzVector* selHadrWJet2, TLorentzVector* selHadrBJet, TLorentzVector* selLeptBJet, TLorentzVector* selLepton, bool isSemiMu, bool isSemiEl){

    //Should use Pt information in stead of E!
    // --> Both concepts are identical in the case of CaloJets, but not in the case of PF
    // --> PF uses massive objects to construct particles!

    histo1D["DeltaR_TFClass_Light1"]->Fill( hadrWJet1->DeltaR(selHadrWJet1) );
    histo1D["DeltaR_TFClass_Light2"]->Fill( hadrWJet2->DeltaR(selHadrWJet2) );
    histo1D["DeltaR_TFClass_HadrB"]->Fill( hadrBJet->DeltaR(selHadrBJet) );
    histo1D["DeltaR_TFClass_LeptB"]->Fill( leptBJet->DeltaR(selLeptBJet) );
    if(isSemiMu) histo1D["DeltaR_TFClass_Mu"]->Fill( lepton->DeltaR(selLepton) );
    if(isSemiEl) histo1D["DeltaR_TFClass_El"]->Fill( lepton->DeltaR(selLepton) );
	
    histo2D["Light_GenPtVsRecoPt"]->Fill(      selHadrWJet1->Pt(),   hadrWJet1->Pt()    );
    histo2D["Light_GenThetaVsRecoTheta"]->Fill(selHadrWJet1->Theta(),hadrWJet1->Theta() );
    histo2D["Light_GenThetaVsRecoPt"]->Fill(   selHadrWJet1->Pt(),   hadrWJet1->Theta() );
    histo2D["Light_GenPhiVsRecoPhi"]->Fill(    selHadrWJet1->Phi(),  hadrWJet1->Phi()   );
    histo2D["Light_GenPhiVsRecoPt"]->Fill(     selHadrWJet1->Pt(),   hadrWJet1->Phi()   );

    histo2D["Light_DiffPtVsRecoPt"]->Fill(      selHadrWJet1->Pt(),   (hadrWJet1 - selHadrWJet1)->Pt()    );
    histo2D["Light_DiffThetaVsRecoTheta"]->Fill(selHadrWJet1->Theta(),(hadrWJet1 - selHadrWJet1)->Theta() );
    histo2D["Light_DiffThetaVsRecoPt"]->Fill(   selHadrWJet1->Pt(),   (hadrWJet1 - selHadrWJet1)->Theta() );
    histo2D["Light_DiffPhiVsRecoPhi"]->Fill(    selHadrWJet1->Phi(),  hadrWJet1->DeltaPhi(selHadrWJet1) );
    histo2D["Light_DiffPhiVsRecoPhi_All"]->Fill(selHadrWJet1->Phi(),  hadrWJet1->DeltaPhi(selHadrWJet1) );
    histo2D["Light_DiffPhiVsRecoPt"]->Fill(     selHadrWJet1->Pt(),   hadrWJet1->DeltaPhi(selHadrWJet1) );
    histo2D["Light_DiffPhiVsRecoPt_All"]->Fill( selHadrWJet1->Pt(),   hadrWJet1->DeltaPhi(selHadrWJet1) );

    histo2D["Light_GenPtVsRecoPt"]->Fill(      selHadrWJet2->Pt(),   hadrWJet2->Pt()    );
    histo2D["Light_GenThetaVsRecoTheta"]->Fill(selHadrWJet2->Theta(),hadrWJet2->Theta() );
    histo2D["Light_GenThetaVsRecoPt"]->Fill(   selHadrWJet2->Pt(),   hadrWJet2->Theta() );
    histo2D["Light_GenPhiVsRecoPhi"]->Fill(    selHadrWJet2->Phi(),  hadrWJet2->Phi()   );
    histo2D["Light_GenPhiVsRecoPt"]->Fill(     selHadrWJet2->Pt(),   hadrWJet2->Phi()   );

    histo2D["Light_DiffPtVsRecoPt"]->Fill(      selHadrWJet2->Pt(),   (hadrWJet2 - selHadrWJet2)->Pt()    );
    histo2D["Light_DiffThetaVsRecoTheta"]->Fill(selHadrWJet2->Theta(),(hadrWJet2 - selHadrWJet2)->Theta() );
    histo2D["Light_DiffThetaVsRecoPt"]->Fill(   selHadrWJet2->Pt(),   (hadrWJet2 - selHadrWJet2)->Theta() );
    histo2D["Light_DiffPhiVsRecoPhi"]->Fill(    selHadrWJet2->Phi(),  hadrWJet2->DeltaPhi(selHadrWJet2)   );
    histo2D["Light_DiffPhiVsRecoPhi_All"]->Fill(selHadrWJet2->Phi(),  hadrWJet2->DeltaPhi(selHadrWJet2)   );
    histo2D["Light_DiffPhiVsRecoPt"]->Fill(     selHadrWJet2->Pt(),   hadrWJet2-Delta>Phi(selHadrWJet2)   );
    histo2D["Light_DiffPhiVsRecoPt_All"]->Fill( selHadrWJet2->Pt(),   hadrWJet2-Delta>Phi(selHadrWJet2)   );

    histo2D["BJet_GenPtVsRecoPt"]->Fill(      selHadrBJet->Pt(),    hadrBJet->Pt()    );
    histo2D["BJet_GenThetaVsRecoTheta"]->Fill(selHadrBJet->Theta(), hadrBJet->Theta() );
    histo2D["BJet_GenThetaVsRecoPt"]->Fill(   selHadrBJet->Pt(),    hadrBJet->Theta() );
    histo2D["BJet_GenPhiVsRecoPhi"]->Fill(    selHadrBJet->Phi(),   hadrBJet->Phi());
    histo2D["BJet_GenPhiVsRecoPt"]->Fill(     selHadrBJet->Pt(),    hadrBJet->Phi());

    histo2D["BJet_DiffPtVsRecoPt"]->Fill(      selHadrBJet->Pt(),   (hadrBJet - selHadrBJet)->Pt()    );
    histo2D["BJet_DiffThetaVsRecoTheta"]->Fill(selHadrBJet->Theta(),(hadrBJet - selHadrBJet)->Theta() );
    histo2D["BJet_DiffThetaVsRecoPt"]->Fill(   selHadrBJet->Pt(),   (hadrBJet - selHadrBJet)->Theta() );
    histo2D["BJet_DiffPhiVsRecoPhi"]->Fill(    selHadrBJet->Phi(),  hadrBJet->DeltaPhi(selHadrBJet)   );
    histo2D["BJet_DiffPhiVsRecoPhi_All"]->Fill(selHadrBJet->Phi(),  hadrBJet->DeltaPhi(selHadrBJet)   );
    histo2D["BJet_DiffPhiVsRecoPt"]->Fill(     selHadrBJet->Pt(),   hadrBJet->DeltaPhi(selHadrBJet)   );
    histo2D["BJet_DiffPhiVsRecoPt_All"]->Fill( selHadrBJet->Pt(),   hadrBJet->DeltaPhi(selHadrBJet)   );

    histo2D["BJet_GenPtVsRecoPt"]->Fill(      selLeptBJet->Pt(),    selLeptBJet->Pt()    );
    histo2D["BJet_GenThetaVsRecoTheta"]->Fill(selLeptBJet->Theta(), selLeptBJet->Theta() );
    histo2D["BJet_GenThetaVsRecoPt"]->Fill(   selLeptBJet->Pt(),    selLeptBJet->Theta() );
    histo2D["BJet_GenPhiVsRecoPhi"]->Fill(    selLeptBJet->Phi(),   selLeptBJet->Phi()   );
    histo2D["BJet_GenPhiVsRecoPt"]->Fill(     selLeptBJet->Pt(),    selLeptBJet->Phi()   );

    histo2D["BJet_DiffPtVsRecoPt"]->Fill(      selLeptBJet->Pt(),   leptBJet->Pt()    - selLeptBJet->Pt()    );
    histo2D["BJet_DiffThetaVsRecoTheta"]->Fill(selLeptBJet->Theta(),leptBJet->Theta() - selLeptBJet->Theta() );
    histo2D["BJet_DiffThetaVsRecoPt"]->Fill(   selLeptBJet->Pt(),   leptBJet->Theta() - selLeptBJet->Theta() );
    histo2D["BJet_DiffPhiVsRecoPhi"]->Fill(    selLeptBJet->Phi(),  leptBJet->Phi()   - selLeptBJet->Phi()   );
    histo2D["BJet_DiffPhiVsRecoPhi_All"]->Fill(selLeptBJet->Phi(),  leptBJet->Phi()   - selLeptBJet->Phi()   );
    histo2D["BJet_DiffPhiVsRecoPt"]->Fill(     selLeptBJet->Pt(),   leptBJet->Phi()   - selLeptBJet->Phi()   );
    histo2D["BJet_DiffPhiVsRecoPt_All"]->Fill( selLeptBJet->Pt(),   leptBJet->Phi()   - selLeptBJet->Phi()   );

    if(isSemiEl){
        histo2D["El_GenPtVsRecoPt"]->Fill(      selLepton->Pt(),   selLepton->Pt()    );
	histo2D["El_GenThetaVsRecoTheta"]->Fill(selLepton->Theta(),selLepton->Theta() );
       	histo2D["El_GenThetaVsRecoPt"]->Fill(   selLepton->Pt(),   selLepton->Theta() );
    	histo2D["El_GenPhiVsRecoPt"]->Fill(     selLepton->Pt(),   selLepton->Phi()   );
	histo2D["El_GenPhiVsRecoPhi"]->Fill(    selLepton->Phi(),  selLepton->Phi()   );

    	histo2D["El_DiffPtVsRecoPt"]->Fill(      selLepton->Pt(),   lepton->Pt()    - selLepton->Pt()    );
    	histo2D["El_DiffThetaVsRecoTheta"]->Fill(selLepton->Theta(),lepton->Theta() - selLepton->Theta() );
    	histo2D["El_DiffThetaVsRecoPt"]->Fill(   selLepton->Pt(),   lepton->Theta() - selLepton->Theta() );
    	histo2D["El_DiffPhiVsRecoPhi"]->Fill(    selLepton->Phi(),  lepton->Phi()   - selLepton->Phi()   );
    	histo2D["El_DiffPhiVsRecoPt"]->Fill(     selLepton->Pt(),   lepton->Phi()   - selLepton->Phi()   );
    	histo2D["El_DiffPhiVsRecoPhi_All"]->Fill(selLepton->Phi(),  lepton->Phi()   - selLepton->Phi()   );
    	histo2D["El_DiffPhiVsRecoPt_All"]->Fill( selLepton->Pt(),   lepton->Phi()   - selLepton->Phi()   );
    }
    if(isSemiMu){
	float InvPtgenMu = 1./lepton->Pt();
    	float InvPtrecMu = 1./selLepton->Pt();
    	histo2D["Mu_GenInvPtVsRecoInvPt"]->Fill(InvPtrecMu,     InvPtgenMu         );
       	histo2D["Mu_GenPtVsRecoPt"]->Fill(      selLepton->Pt(),   selLepton->Pt()    );
    	histo2D["Mu_GenThetaVsRecoTheta"]->Fill(selLepton->Theta(),selLepton->Theta() );
    	histo2D["Mu_GenThetaVsRecoInvPt"]->Fill(InvPtrecMu,     selLepton->Theta() );
    	histo2D["Mu_GenThetaVsRecoPt"]->Fill(   selLepton->Pt(),   selLepton->Theta() );
    	histo2D["Mu_GenPhiVsRecoPhi"]->Fill(    selLepton->Phi(),  selLepton->Phi()   );
    	histo2D["Mu_GenPhiVsRecoInvPt"]->Fill(  InvPtrecMu,     selLepton->Phi()   );
    	histo2D["Mu_GenPhiVsRecoPt"]->Fill(     selLepton->Pt(),   selLepton->Phi()   );

    	histo2D["Mu_DiffInvPtVsRecoInvPt"]->Fill(    InvPtrecMu,     InvPtgenMu      - InvPtrecMu         );
    	histo2D["Mu_DiffInvPtVsRecoInvPt_All"]->Fill(InvPtrecMu,     InvPtgenMu      - InvPtrecMu         );
    	histo2D["Mu_DiffPtVsRecoPt"]->Fill(          selLepton->Pt(),   lepton->Pt()    - selLepton->Pt()    );
    	histo2D["Mu_DiffThetaVsRecoTheta"]->Fill(    selLepton->Theta(),lepton->Theta() - selLepton->Theta() );
    	histo2D["Mu_DiffThetaVsRecoInvPt"]->Fill(    InvPtrecMu,     lepton->Theta() - selLepton->Theta() );
    	histo2D["Mu_DiffThetaVsRecoPt"]->Fill(       selLepton->Pt(),   lepton->Theta() - selLepton->Theta() );
    	histo2D["Mu_DiffPhiVsRecoPhi"]->Fill(        selLepton->Phi(),  lepton->Phi()   - selLepton->Phi()   );
    	histo2D["Mu_DiffPhiVsRecoPhi_All"]->Fill(    selLepton->Phi(),  lepton->Phi()   - selLepton->Phi()   );
    	histo2D["Mu_DiffPhiVsRecoInvPt"]->Fill(      InvPtrecMu,     lepton->Phi()   - selLepton->Phi()   );
    	histo2D["Mu_DiffPhiVsRecoInvPt_All"]->Fill(  InvPtrecMu,     lepton->Phi()   - selLepton->Phi()   );
	histo2D["Mu_DiffPhiVsRecoPt"]->Fill(         selLepton->Pt(),   lepton->Phi()   - selLepton->Phi()   );
	histo2D["Mu_DiffPhiVsRecoPt_All"]->Fill(     selLepton->Pt(),   lepton->Phi()   - selLepton->Phi()   );
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
          TH2F* histoForFit = new TH2F();
	  for (unsigned int f=0; f<12;f++) {				
				
	    switch(f){
	      case 0:
		histoForFit=histo2D["Light_DiffPtVsRecoPt"]; break;
	      case 1:
		histoForFit=histo2D["BJet_DiffPtVsRecoPt"]; break;
	      case 2:
		histoForFit=histo2D["El_DiffPtVsRecoPt"]; break;
	      case 3:
		histoForFit=histo2D["Mu_DiffInvPtVsRecoInvPt"]; break;
	      case 4:
		histoForFit=histo2D["Light_DiffThetaVsRecoPt"]; break;
	      case 5:
		histoForFit=histo2D["BJet_DiffThetaVsRecoPt"]; break;
	      case 6:
		histoForFit=histo2D["El_DiffThetaVsRecoPt"]; break;
	      case 7:
		histoForFit=histo2D["Mu_DiffThetaVsRecoInvPt"]; break;
	      case 8:
		histoForFit=histo2D["Light_DiffPhiVsRecoPt"]; break;
	      case 9:
		histoForFit=histo2D["BJet_DiffPhiVsRecoPt"]; break;
	      case 10:
		histoForFit=histo2D["El_DiffPhiVsRecoPt"]; break;
	      case 11:
		histoForFit=histo2D["Mu_DiffPhiVsRecoInvPt"]; break;
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

void TFCreation::FitSliceClassCode(TH2F* histoFit, int npar, const char* parnames[]){
	//------------------------------------------------------------------------------------------//
	// Main difference with the Root class FitSlicesY() is the plotting of histograms !        
	// In the Root class the distribution of each hlist histogram is not given!
	// --> Useful to use the own code when needing control histograms
	//
	// Other difference between the two codes have been removed!
	// Originally the treatment of the overflow bin was different, but is now made similar!
	//Create one histogram for each function parameter -> 6 histograms for each 2D plot
	for(int ipar=0 ; ipar < npar; ipar++){

	    hlist[ipar] = new TH1D( (string(histoFit->GetName())+"_"+parnames[ipar]).c_str(), (string(histoFit->GetName())+" : Fitted value of "+parnames[ipar]).c_str(), histoFit->GetXaxis()->GetNbins(), histoFit->GetXaxis()->GetXmin(), histoFit->GetXaxis()->GetXmax());
	    hlist[ipar]->GetXaxis()->SetTitle(histoFit->GetXaxis()->GetTitle());
	}
	hlist[npar] = new TH1D( (string(histoFit->GetName())+"_chi2").c_str(), (string(histoFit->GetName())+": #chi^{2} distribution for "+string(doubleGaussianFit->GetExpFormula())).c_str(), histoFit->GetXaxis()->GetNbins(), histoFit->GetXaxis()->GetXmin(), histoFit->GetXaxis()->GetXmax() );

	//Loop on all bins in X, generate a projection along Y and fit each bin separately!
	int cut = 0; // require a minimum number of bins in the slice to be filled --> Should this ever be larger than 0 ??
	int nbins = histoFit->GetXaxis()->GetNbins();
	cout << "\n ** Looking at histogram : " << histoFit->GetName() << "                               ******************************* " << endl;
	for(int bin=1;bin <= nbins+1;bin ++) {
	    cout << "   --  Looking at bin : " << bin << endl;
	    string projection_title = string(histoFit->GetName())+"_sliceYbin"+tostr(bin);

	    TH1D *hp = histoFit->ProjectionY(projection_title.c_str(),bin,bin,"e");
	    //if(bin==nbins) hp = histoFit->ProjectionY(projection_title.c_str(),bin,bin+1,"e"); //include overflow in last bin
	    //if(bin==1) hp = histoFit->ProjectionY(projection_title.c_str(),bin-1,bin,"e"); //include underflow in first bin

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
		    hlist[ipar]->Fill(histoFit->GetXaxis()->GetBinCenter(bin+1/2),doubleGaussianFit->GetParameter(ipar)); 
                    hlist[ipar]->SetBinError( (int) (bin+1/2) ,doubleGaussianFit->GetParError(ipar)); //WHY +1/2 .... (Is bin size always equal to 1 .. )?

	        }
		//Save hchi2 histogram as extra hlist!
	        hlist[npar]->Fill(histoFit->GetXaxis()->GetBinCenter(bin+1/2),doubleGaussianFit->GetChisquare()/(npfits-npar));
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

void TFCreation::WriteTF(TH2F* fitHisto, ostream &myTFs, ostream &myTFCard, TFile* file){

    const int NrConsideredPars = 6, NrConsideredCaloPars = 3;
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
