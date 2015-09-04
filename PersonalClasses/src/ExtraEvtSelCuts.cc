#include "../interface/ExtraEvtSelCuts.h"

ExtraEvtSelCuts::ExtraEvtSelCuts(){

  //Do not do anything here since this class should only be used in the case of the TTBarSemiLepton dataset!
  //This because it uses the correctJetCombi information in order to decide on the optimal cut values!
}

ExtraEvtSelCuts::~ExtraEvtSelCuts(){

}

void ExtraEvtSelCuts::Initialize(float TopMassHadr, float sTopMassHadr, float WMassHadr, float sWMassHadr, bool oneBTag, std::string bTagTitle, int chiSqCut, int massWindow){

  //Initialize the counters for the chi-sq and mW/mTop cuts!
  float ChiSqCutValues[5] = {50, 30, 20, 10, 5};
  float MassWindowSigmas[5] = {5, 3, 2, 1.5, 1}; 
  CorrEvts_ = 0; WrongEvts_ = 0; UnmatchEvts_ = 0;
  bTitle_ = bTagTitle;
  oneBTag_ = oneBTag;

  nrCuts_ = 5;
  //if(chiSqCut != -1 && massWindow != -1) nrCuts_ = 1;
  //else if( (chiSqCut != -1 && massWindow == -1) || (chiSqCut == -1 && massWindow != -1) ) cout << " ERROR : Only one cut-value given, so all options will be considered! " << endl;
  //Maybe best to use this only for the _All option ==> Then it is possible to try out the most optimal combination!
  chosenChiSq_ = chiSqCut; chosenMassWindow_ = massWindow;
  CorrEvts_AllChosenCuts_ = 0;   WrongEvts_AllChosenCuts_ = 0;   UnmatchEvts_AllChosenCuts_ = 0;

  //Now transfer this to variables existing throughout the entire class:
  for(int iCut = 0; iCut < nrCuts_; iCut++){
    //if(nrCuts_ == 1){
    //  ChiSqCutVal_[iCut] = ChiSqCutValues[chiSqCut];
    //  MassWindow_Sigmas_[iCut] = MassWindowSigmas[massWindow];
    //}
    //else{
      ChiSqCutVal_[iCut] = ChiSqCutValues[iCut];
      MassWindow_Sigmas_[iCut] = MassWindowSigmas[iCut];
    //}

    CorrEvts_ChiSq_[iCut] = 0; WrongEvts_ChiSq_[iCut] = 0; UnmatchEvts_ChiSq_[iCut] = 0;
    CorrEvts_MT_[iCut] = 0;    WrongEvts_MT_[iCut] = 0;    UnmatchEvts_MT_[iCut] = 0;
    CorrEvts_MW_[iCut] = 0;    WrongEvts_MW_[iCut] = 0;    UnmatchEvts_MW_[iCut] = 0;
    CorrEvts_MComb_[iCut] = 0; WrongEvts_MComb_[iCut] = 0; UnmatchEvts_MComb_[iCut] = 0;

    MTop_Down_[iCut] = TopMassHadr - MassWindow_Sigmas_[iCut]*sTopMassHadr; MTop_Up_[iCut] = TopMassHadr + MassWindow_Sigmas_[iCut]*sTopMassHadr; 
    MW_Down_[iCut] = WMassHadr - MassWindow_Sigmas_[iCut]*sWMassHadr; MW_Up_[iCut] = WMassHadr + MassWindow_Sigmas_[iCut]*sWMassHadr; 
  }

  //Get the Mlb-Mqqb chiSq for correct and wrong events to decide on a cut-value
  histo1D["LowestChiSq_CorrectEvents_"+bTitle_]   = new TH1F(("LowestChiSq_CorrectEvents_"+bTitle_).c_str(),   ("#chi^{2} distribution for chosen jet-combination (correct events -- "+bTitle_+")").c_str(),   50,0,80);
  histo1D["LowestChiSq_WrongEvents_"+bTitle_]     = new TH1F(("LowestChiSq_WrongEvents_"+bTitle_).c_str(),     ("#chi^{2} distribution for chosen jet-combination (wrong events -- "+bTitle_+")").c_str(),     50,0,80);
  histo1D["LowestChiSq_UnmatchedEvents_"+bTitle_] = new TH1F(("LowestChiSq_UnmatchedEvents_"+bTitle_).c_str(), ("#chi^{2} distribution for chosen jet-combination (unmatched events -- "+bTitle_+")").c_str(), 50,0,80);
}

bool ExtraEvtSelCuts::KeepEvent(vector<int> jetComb, TLorentzVector lepton, vector<TLorentzVector> Jets, vector<int> selJetComb, float chiSq){

  bool keepEvt = false;  

  //Should only be doing all this when one b-tag should be considered!    
  if(oneBTag_){

    //Booleans in order to combine the different cuts!
    bool ChiSqSurv[5] = {false}, MTSurv[5] = {false}, MWSurv[5] = {false};

    if(jetComb[0] != 9999 && jetComb[1] != 9999 && jetComb[2] != 9999 && jetComb[3] != 9999){

      //Correct/Wrong events:
      if( jetComb[0] == selJetComb[0] && jetComb[1] == selJetComb[1] && (jetComb[2] == selJetComb[2] || jetComb[2] == selJetComb[3]) && (jetComb[3] == selJetComb[2] || jetComb[3] == selJetComb[3]) ){
        histo1D["LowestChiSq_CorrectEvents_"+bTitle_]->Fill(chiSq);
        CorrEvts_++;

        for(int iCut = 0; iCut < nrCuts_; iCut++){
          if(chiSq < ChiSqCutVal_[iCut] ){CorrEvts_ChiSq_[iCut]++; ChiSqSurv[iCut] = true;}
          if((Jets[selJetComb[1]]+Jets[selJetComb[2]]+Jets[selJetComb[3]]).M() > MTop_Down_[iCut] && (Jets[selJetComb[1]]+Jets[selJetComb[2]]+Jets[selJetComb[3]]).M() < MTop_Up_[iCut]){CorrEvts_MT_[iCut]++; MTSurv[iCut] = true;}
          if((Jets[selJetComb[2]]+Jets[selJetComb[3]]).M() > MW_Down_[iCut] && (Jets[selJetComb[2]]+Jets[selJetComb[3]]).M() < MW_Up_[iCut] ){CorrEvts_MW_[iCut]++; MWSurv[iCut] = true;}

          if(MTSurv[iCut] && MWSurv[iCut]) CorrEvts_MComb_[iCut]++;
        }
        if(ChiSqSurv[chosenChiSq_] && MTSurv[chosenMassWindow_] && MWSurv[chosenMassWindow_]){ CorrEvts_AllChosenCuts_++; keepEvt = true;}
      }
      else{
        histo1D["LowestChiSq_WrongEvents_"+bTitle_]->Fill(chiSq);
        WrongEvts_++;

        for(int iCut = 0; iCut < nrCuts_; iCut++){
          if(chiSq < ChiSqCutVal_[iCut] ){ WrongEvts_ChiSq_[iCut]++; ChiSqSurv[iCut] = true;}
          if((Jets[selJetComb[1]]+Jets[selJetComb[2]]+Jets[selJetComb[3]]).M() > MTop_Down_[iCut] && (Jets[selJetComb[1]]+Jets[selJetComb[2]]+Jets[selJetComb[3]]).M() < MTop_Up_[iCut]){WrongEvts_MT_[iCut]++;MTSurv[iCut] = true;}
          if((Jets[selJetComb[2]]+Jets[selJetComb[3]]).M() > MW_Down_[iCut] && (Jets[selJetComb[2]]+Jets[selJetComb[3]]).M() < MW_Up_[iCut]){ WrongEvts_MW_[iCut]++; MWSurv[iCut] = true;}
  
          if(MTSurv[iCut] && MWSurv[iCut]) WrongEvts_MComb_[iCut]++;
        }
        if(ChiSqSurv[chosenChiSq_] && MTSurv[chosenMassWindow_] && MWSurv[chosenMassWindow_]){ WrongEvts_AllChosenCuts_++; keepEvt = true;}
      }
    }
    else{
      histo1D["LowestChiSq_UnmatchedEvents_"+bTitle_]->Fill(chiSq);
      UnmatchEvts_++;

      for(int iCut = 0; iCut < nrCuts_; iCut++){
        if( chiSq < ChiSqCutVal_[iCut] ){ UnmatchEvts_ChiSq_[iCut]++; ChiSqSurv[iCut] = true;}
        if((Jets[selJetComb[1]]+Jets[selJetComb[2]]+Jets[selJetComb[3]]).M() > MTop_Down_[iCut] && (Jets[selJetComb[1]]+Jets[selJetComb[2]]+Jets[selJetComb[3]]).M() < MTop_Up_[iCut]){UnmatchEvts_MT_[iCut]++;MTSurv[iCut] = true;}
        if( (Jets[selJetComb[2]]+Jets[selJetComb[3]]).M() > MW_Down_[iCut] && (Jets[selJetComb[2]]+Jets[selJetComb[3]]).M() < MW_Up_[iCut]){ UnmatchEvts_MW_[iCut]++; MWSurv[iCut] = true;}

        if(MTSurv[iCut] && MWSurv[iCut]) UnmatchEvts_MComb_[iCut]++;
      }
      if(ChiSqSurv[chosenChiSq_] && MTSurv[chosenMassWindow_] && MWSurv[chosenMassWindow_]){ UnmatchEvts_AllChosenCuts_++; keepEvt = true;}
    }
  }

  return keepEvt;
}

void ExtraEvtSelCuts::StoreCutInfluence(TFile* outfile){

  if(oneBTag_){

    ofstream CutInfl;
    CutInfl.open("EventSelectionResults/AnalyzerOutput/ExtraEventSelectionCuts.tex");

    //Initialize the tex document
    CutInfl << "\\documentclass{article} \n\\usepackage[margin=0.5in]{geometry} \n\\begin{document} \n" << endl;

    CutInfl << " \\begin{abstract} \n " << endl;
    CutInfl << "   The tables in this document represent the influence of the additional event selection cuts that were applied in order to reduce the number of selected events for CPU reasons. \\\\ " << endl;
    CutInfl << "   The considered cuts are rather basic and are merely developed to reduce the number of so-called wrong events \\\\ " << endl;
    CutInfl << "   \\begin{itemize} \n     \\item Cut on Mlb-Mqqb $\\chi^{2}$ distribution \\\\ \n     \\item Cut on top and W-mass window \\\\ \n   \\end{itemize} \n " << endl;
    CutInfl << "   \\textbf{Created on :} \\today " << endl;
    CutInfl << " \\end{abstract} \n " << endl;
 
    //Statistics for all events 
    CutInfl << " \\begin{table}[h!t] \n  \\caption{Original number of events, before any additional event selection cuts are applied} \n  \\centering \n   \\begin{tabular}{c|c|c|c|c} " << endl;
    CutInfl << "     Correct evts    & Wrong evts     & Unmatched evts      &  Total evts & s/b ($\\%$)     \\\\ \n     \\hline" << endl;
    int TotalEvts = CorrEvts_+WrongEvts_+UnmatchEvts_;
    CutInfl << "     " << CorrEvts_ << "   &   " << WrongEvts_ << "  & " << UnmatchEvts_ << " & " << TotalEvts << "  & " << CorrEvts_*100.0/(CorrEvts_+WrongEvts_+UnmatchEvts_) << " \n " << endl;
    CutInfl << "   \\end{tabular} \n \\end{table} \n " << endl;

    //Numbers after chi-sq cut
    CutInfl << " \\begin{table}[h!t] \n  \\caption{Remaining number of events, after applying cuts on the Mlb-Mqqb $\\chi^{2}$.} \n  \\centering \n   \\begin{tabular}{c|c|c|c|c|c} " << endl;
    CutInfl << "     $\\chi^{2}$ cut-value    & Correct evts    & Wrong evts     & Unmatched evts  & Evt reduction ($\\%$)    & s/b ($\\%$)     \\\\ \n     \\hline" << endl;
    for(int iCut = 0; iCut < nrCuts_; iCut++){
      CutInfl << "     " << ChiSqCutVal_[iCut] << "  &   " << CorrEvts_ChiSq_[iCut] << "  &  " << WrongEvts_ChiSq_[iCut] << "  &   " << UnmatchEvts_ChiSq_[iCut] << " &  " << 100.0-((CorrEvts_ChiSq_[iCut]+WrongEvts_ChiSq_[iCut]+UnmatchEvts_ChiSq_[iCut])*100.0/TotalEvts) << " & " << CorrEvts_ChiSq_[iCut]*100.0/(CorrEvts_ChiSq_[iCut]+WrongEvts_ChiSq_[iCut]+UnmatchEvts_ChiSq_[iCut]);
      if(iCut < nrCuts_-1) CutInfl << " \\\\ " << endl;
    }
    CutInfl << "\n   \\end{tabular} \n \\end{table} \n " << endl;

    //Numbers after mT mass-window
    CutInfl << " \\begin{table}[h!t] \n  \\caption{Remaining number of events, after applying cuts on the $M_{top}$ mass-window.} \n  \\centering \n   \\begin{tabular}{c|c|c|c|c|c|c} " << endl;
    CutInfl << "     Nr $\\sigma$'s & Mass-window ($m_{top}$)    & Correct evts    & Wrong evts     & Unmatched evts   & Evt reduction ($\\%$)    & s/b ($\\%$)     \\\\ \n     \\hline" << endl;
    for(int iCut = 0; iCut < nrCuts_; iCut++){
      CutInfl << "     " << MassWindow_Sigmas_[iCut] << " & " << MTop_Down_[iCut] << " - " << MTop_Up_[iCut] << "  &   " << CorrEvts_MT_[iCut] << "  &  " << WrongEvts_MT_[iCut] << "  &   " << UnmatchEvts_MT_[iCut] << " & " << 100.0-((CorrEvts_MT_[iCut]+WrongEvts_MT_[iCut]+UnmatchEvts_MT_[iCut])*100.0/TotalEvts) << " &  " << CorrEvts_MT_[iCut]*100.0/(CorrEvts_MT_[iCut]+WrongEvts_MT_[iCut]+UnmatchEvts_MT_[iCut]);
      if(iCut < nrCuts_-1) CutInfl << " \\\\ " << endl;
    }
    CutInfl << "\n   \\end{tabular} \n \\end{table} \n " << endl;

    //Numbers after mW mass-window
    CutInfl << " \\begin{table}[h!t] \n  \\caption{Remaining number of events, after applying cuts on the $M_{W}$ mass-window.} \n  \\centering \n   \\begin{tabular}{c|c|c|c|c|c|c} " << endl;
    CutInfl << "     Nr $\\sigma$'s & Mass-window ($m_{top}$)    & Correct evts    & Wrong evts     & Unmatched evts   & Evt reduction ($\\%$)    & s/b ($\\%$)     \\\\ \n     \\hline" << endl;
    for(int iCut = 0; iCut < nrCuts_; iCut++){
      CutInfl << "     " << MassWindow_Sigmas_[iCut] << " & " << MW_Down_[iCut] << " - " << MW_Up_[iCut] << "  &   " << CorrEvts_MW_[iCut] << "  &  " << WrongEvts_MW_[iCut] << "  &   " << UnmatchEvts_MW_[iCut]  << " & " << 100.0-((CorrEvts_MW_[iCut]+WrongEvts_MW_[iCut]+UnmatchEvts_MW_[iCut])*100.0/TotalEvts) << " &  " << CorrEvts_MW_[iCut]*100.0/(CorrEvts_MW_[iCut]+WrongEvts_MW_[iCut]+UnmatchEvts_MW_[iCut]);
      if(iCut < nrCuts_-1) CutInfl << " \\\\ " << endl;
    }
    CutInfl << "\n   \\end{tabular} \n \\end{table} \n " << endl;

    //Numbers after combined mT-mW mass-window
    CutInfl << " \\begin{table}[h!t] \n  \\caption{Remaining number of events, after applying cuts on the combined $M_{top}-M_{W}$ mass-window.} \n  \\centering \n   \\begin{tabular}{c|c|c|c|c|c} " << endl;
    CutInfl << "     Mass-window $\\sigma$'s & Correct evts    & Wrong evts     & Unmatched evts  & Evt reduction ($\\%$)   & s/b ($\\%$)     \\\\ \n     \\hline" << endl;
    for(int iCut = 0; iCut < nrCuts_; iCut++){
      CutInfl << "     " << MassWindow_Sigmas_[iCut] << " & " << CorrEvts_MComb_[iCut] << "  &  " << WrongEvts_MComb_[iCut] << "  &   " << UnmatchEvts_MComb_[iCut] << " & " << 100.0-((CorrEvts_MComb_[iCut]+WrongEvts_MComb_[iCut]+UnmatchEvts_MComb_[iCut])*100.0/TotalEvts) << " &  " << CorrEvts_MComb_[iCut]*100.0/(CorrEvts_MComb_[iCut]+WrongEvts_MComb_[iCut]+UnmatchEvts_MComb_[iCut]);
      if(iCut < nrCuts_-1) CutInfl << " \\\\ " << endl;
    }
    CutInfl << "\n   \\end{tabular} \n \\end{table} \n " << endl;

    //Numbers after all cuts combined
    if(chosenChiSq_ != -1 && chosenMassWindow_ != -1){
      CutInfl << " \\begin{table}[h!t] \n  \\caption{Remaining number of events, after applying cuts on both the Mlb-Mqqb $\\chi^{2}$ and the combined $M_{top}-M_{W}$ mass-window.} \n  \\centering \n   \\begin{tabular}{c|c|c|c|c|c|c} " << endl;
      CutInfl << "     Mass-window $\\sigma$'s & $\\chi^{2}$ cut-value    & Correct evts    & Wrong evts     & Unmatched evts  & Evt reduction ($\\%$)    & s/b ($\\%$)     \\\\ \n     \\hline" << endl;
      CutInfl << "     " << MassWindow_Sigmas_[chosenMassWindow_] << " & " << ChiSqCutVal_[chosenChiSq_] << " & " << CorrEvts_AllChosenCuts_ << " & " << WrongEvts_AllChosenCuts_ << " & " << UnmatchEvts_AllChosenCuts_ << " & " <<100.0-((CorrEvts_AllChosenCuts_+WrongEvts_AllChosenCuts_+UnmatchEvts_AllChosenCuts_)*100.0/TotalEvts)  << " & " << CorrEvts_AllChosenCuts_*100.0/(CorrEvts_AllChosenCuts_+WrongEvts_AllChosenCuts_+UnmatchEvts_AllChosenCuts_) << endl;
      CutInfl << "\n   \\end{tabular} \n \\end{table} \n " << endl;
    }
    else{
      cout << " No cut-values given for chi-sq and mass-window " << endl;
    }
  
    CutInfl << "\\end{document} " << endl;
    CutInfl.close();

    //Also store the plots!!
    outfile->cd();
    if(histo1D.size() > 0){
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
    }
   
    if(histo2D.size() > 0){ 
      TDirectory* th2dir = outfile->GetDirectory("2D_histograms_ExtraEvtSelCuts");
      if(!th2dir) th2dir = outfile->mkdir("2D_histograms_ExtraEvtSelCuts");
      th2dir->cd();
      for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++){    
        TH2F *temp = it->second;
        temp->Write();
      }
    }
    outfile->cd(); 
  }
}
