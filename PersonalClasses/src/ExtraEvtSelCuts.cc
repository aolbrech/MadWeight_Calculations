#include "../interface/ExtraEvtSelCuts.h"

ExtraEvtSelCuts::ExtraEvtSelCuts(){

  //Do not do anything here since this class should only be used in the case of the TTBarSemiLepton dataset!
  //This because it uses the correctJetCombi information in order to decide on the optimal cut values!
}

ExtraEvtSelCuts::~ExtraEvtSelCuts(){

}

void ExtraEvtSelCuts::Initialize(float TopMassHadr, float sTopMassHadr, float WMassHadr, float sWMassHadr, std::string bTagTitle){

  //Initialize the counters for the chi-sq and mW/mTop cuts!
  const int nrEvtSelCuts = 5;
  float ChiSqCutValues[nrEvtSelCuts] = {50, 30, 20, 10, 5};
  float MassWindowSigmas[nrEvtSelCuts] = {5, 3, 2, 1.5, 1}; 
  CorrEvts_ = 0; WrongEvts_ = 0; UnmatchEvts_ = 0;
  bTitle_ = bTagTitle;

  nrCuts_ = nrEvtSelCuts;
  //Now transfer this to variables existing throughout the entire class:
  for(int iCut = 0; iCut < nrEvtSelCuts; iCut++){
    CorrEvts_ChiSq_[iCut] = 0;      WrongEvts_ChiSq_[iCut] = 0;      UnmatchEvts_ChiSq_[iCut] = 0;
    CorrEvts_MassWindow_[iCut] = 0; WrongEvts_MassWindow_[iCut] = 0; UnmatchEvts_MassWindow_[iCut] = 0;
    CorrEvts_Comb_[iCut] = 0;       WrongEvts_Comb_[iCut] = 0;       UnmatchEvts_Comb_[iCut] = 0;

    //MassWindow_[iCut] = MassWindowSigmas[iCut]; 
    ChiSqCutVal_[iCut] = ChiSqCutValues[iCut];
    MTop_Down_[iCut] = TopMassHadr - MassWindowSigmas[iCut]*sTopMassHadr; MTop_Up_[iCut] = TopMassHadr + MassWindowSigmas[iCut]*sTopMassHadr; 
  }

  //Get the Mlb-Mqqb chiSq for correct and wrong events to decide on a cut-value
  histo1D["LowestChiSq_CorrectEvents_"+bTitle_]   = new TH1F(("LowestChiSq_CorrectEvents_"+bTitle_).c_str(),   ("#chi^{2} distribution for chosen jet-combination (correct events -- "+bTitle_+")").c_str(),   50,0,80);
  histo1D["LowestChiSq_WrongEvents_"+bTitle_]     = new TH1F(("LowestChiSq_WrongEvents_"+bTitle_).c_str(),     ("#chi^{2} distribution for chosen jet-combination (wrong events -- "+bTitle_+")").c_str(),     50,0,80);
  histo1D["LowestChiSq_UnmatchedEvents_"+bTitle_] = new TH1F(("LowestChiSq_UnmatchedEvents_"+bTitle_).c_str(), ("#chi^{2} distribution for chosen jet-combination (unmatched events -- "+bTitle_+")").c_str(), 50,0,80);
}

void ExtraEvtSelCuts::KeepEvent(vector<int> jetCombi, TLorentzVector lepton, vector<TLorentzVector> Jets, vector<int> selJetCombi, float chiSq){
        
  if(jetCombi[0] != 9999 && jetCombi[1] != 9999 && jetCombi[2] != 9999 && jetCombi[3] != 9999){

    //Correct/Wrong events:
    if( jetCombi[0] == selJetCombi[0] && jetCombi[1] == selJetCombi[1] && (jetCombi[2] == selJetCombi[2] || jetCombi[2] == selJetCombi[3]) && (jetCombi[3] == selJetCombi[2] || jetCombi[3] == selJetCombi[3]) ){
      histo1D["LowestChiSq_CorrectEvents_"+bTitle_]->Fill(chiSq);
      CorrEvts_++;

      for(int iCut = 0; iCut < nrCuts_; iCut++){
        if( chiSq < ChiSqCutVal_[iCut] ) CorrEvts_ChiSq_[iCut]++;
        if( (Jets[selJetCombi[1]]+Jets[selJetCombi[2]]+Jets[selJetCombi[3]]).M() > MTop_Down_[iCut] && (Jets[selJetCombi[1]]+Jets[selJetCombi[2]]+Jets[selJetCombi[3]]).M() < MTop_Up_[iCut] )
          CorrEvts_MassWindow_[iCut] ++;
      }
    }
    else{
      histo1D["LowestChiSq_WrongEvents_"+bTitle_]->Fill(chiSq);
      WrongEvts_++;

      for(int iCut = 0; iCut < nrCuts_; iCut++){
        if( chiSq < ChiSqCutVal_[iCut] ) WrongEvts_ChiSq_[iCut]++;
        if( (Jets[selJetCombi[1]]+Jets[selJetCombi[2]]+Jets[selJetCombi[3]]).M() > MTop_Down_[iCut] && (Jets[selJetCombi[1]]+Jets[selJetCombi[2]]+Jets[selJetCombi[3]]).M() < MTop_Up_[iCut] )
          WrongEvts_MassWindow_[iCut] ++;
      }
    }
  }
  else{
    histo1D["LowestChiSq_UnmatchedEvents_"+bTitle_]->Fill(chiSq);
    UnmatchEvts_++;

    for(int iCut = 0; iCut < nrCuts_; iCut++){
      if( chiSq < ChiSqCutVal_[iCut] ) UnmatchEvts_ChiSq_[iCut]++;
      if( (Jets[selJetCombi[1]]+Jets[selJetCombi[2]]+Jets[selJetCombi[3]]).M() > MTop_Down_[iCut] && (Jets[selJetCombi[1]]+Jets[selJetCombi[2]]+Jets[selJetCombi[3]]).M() < MTop_Up_[iCut] )
        UnmatchEvts_MassWindow_[iCut] ++;
    }
  }
}

void ExtraEvtSelCuts::StoreCutInfluence(TFile* outfile){

  //ofstream CutInfl.open("EventSelectionResults/AnalyzerOutput/ExtraEventSelectionCuts.tex");

  std::cout << "\n   -------------------------  Studying influence of additional event selection cuts  ------------------------------ " << endl; 
      std::cout << "      All Evts                    " << CorrEvts_ << "                          " << WrongEvts_ << "           " << CorrEvts_*100.0/(CorrEvts_+WrongEvts_) << "\n " << endl;

      std::cout << "  Chi-sq cut value          # correct events              # wrong events           s/b (%)" << std::endl;
      for(int iCut = 0; iCut < nrCuts_; iCut++){
        if( ChiSqCutVal_[iCut] < 10) std::cout << " ";
        std::cout << "     " << ChiSqCutVal_[iCut] << "                      " << CorrEvts_ChiSq_[iCut] << "                          " << WrongEvts_ChiSq_[iCut] << "           " << CorrEvts_ChiSq_[iCut]*100.0/(CorrEvts_ChiSq_[iCut]+WrongEvts_ChiSq_[iCut]) << endl;
      }

      std::cout << "\n  Mass-window (mT)             # correct events              # wrong events           s/b (%)" << std::endl;
      for(int iCut = 0; iCut < nrCuts_; iCut++){
        std::cout << " " << MTop_Down_[iCut] << " - " << MTop_Up_[iCut] << "             " << CorrEvts_MassWindow_[iCut] << "                          " << WrongEvts_MassWindow_[iCut] << "           " << CorrEvts_MassWindow_[iCut]*100.0/(CorrEvts_MassWindow_[iCut]+WrongEvts_MassWindow_[iCut]) << endl;
      }


  //Also store the plots!!
  outfile->cd();
    TDirectory* th1dir = outfile->GetDirectory("1D_histograms_ExtraEvtSelCuts");   //Check whether directory already exists ..
    if(!th1dir) th1dir = outfile->mkdir("1D_histograms_ExtraEvtSelCuts");          // .. and otherwise create it!
    th1dir->cd();
    for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++){
      TH1F *temp = it->second;
      int N = temp->GetNbinsX();
      temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
      temp->SetBinContent(N+1,0);
      temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
      temp->Write();
    }
    
    TDirectory* th2dir = outfile->GetDirectory("2D_histograms_ExtraEvtSelCuts");
    if(!th2dir) th2dir = outfile->mkdir("2D_histograms_ExtraEvtSelCuts");
    th2dir->cd();
    for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++){    
      TH2F *temp = it->second;
      temp->Write();
    }
    outfile->cd(); 
}
