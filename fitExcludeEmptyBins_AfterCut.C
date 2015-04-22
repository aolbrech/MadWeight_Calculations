#include "TH1.h"
#include "TF1.h"
#include "TList.h"
#include "TFile.h"
#include "TCanvas.h"

Double_t Minimum(TH1F* histo)
{
  Double_t minimum = histo->GetBinContent(1);
  for(int ii = 0; ii < histo->GetNbinsX(); ii++){
    if( histo->GetBinContent(ii+1) < minimum && histo->GetBinContent(ii+1) > 0.0)
      minimum = histo->GetBinContent(ii+1);
  }
  return minimum;
}

void fitExcludeEmptyBins_AfterCut() {
  //Get the histogram from the ROOT file   --> Name gets changed by pythonScript!!
  TFile* InputFile = new TFile(" Events/TopMass_RecoUnmatched_SingleGausTF_10000Evts/FitDeviation_UnmatchedReco_MTop.root ","r");
  TFile* OutputFile = new TFile(" Events/TopMass_RecoUnmatched_SingleGausTF_10000Evts/LimitedFitResult_UnmatchedReco_MTop_ScdDerBothCutApplied.root ","RECREATE");

  bool GenLevel = false; 

  TH1F* LL    = (TH1F*) InputFile->Get("LikelihoodAfterCuts/SignSecondDerivative/LLPosScdDerBoth"); 
  TH1F* LLXS  = (TH1F*) InputFile->Get("LikelihoodAfterCuts/SignSecondDerivative/LLXSPosScdDerBoth"); 
  TH1F* LLAcc = (TH1F*) InputFile->Get("LikelihoodAfterCuts/SignSecondDerivative/LLAccPosScdDerBoth"); 

  TF1 *func = new TF1("func","pol2",171,175); 
  func->SetName("fit_LL");
  LL->Fit(func);
  TCanvas* LLCanvas = new TCanvas("LLCanvas","Distribution of -ln(likelihood) (no normalisation)");
  LLCanvas->cd();
  LL->SetMinimum(Minimum(LL)*0.999);
  func->Draw();
  func->SetTitle(LL->GetTitle());
  LL->SetMarkerSize(1.2);
  LL->SetMarkerStyle(20);
  LL->Draw("SAMEP");
  LLCanvas->Write();
  func->Write();

  TF1 *funcXS = new TF1("funcXS","pol2",171,175); 
  funcXS->SetName("fit_LLXS");
  LLXS->Fit(funcXS);
  TCanvas* LLXSCanvas = new TCanvas("LLXSCanvas","Distribution of -ln(likelihood) (XS normalisation)");
  LLXSCanvas->cd();
  LLXS->SetMarkerStyle(21);
  LLXS->SetMarkerColor(2);
  LLXS->SetMinimum(Minimum(LLXS)*0.999);
  funcXS->Draw();
  funcXS->SetTitle(LLXS->GetTitle());
  LLXS->Draw("PSAME");
  LLXSCanvas->Write();
  funcXS->Write();

  if(GenLevel == false){
  TF1 *funcAcc = new TF1("funcAcc","pol2",171,175); 
    funcAcc->SetName("fit_LLAcc");
    LLAcc->Fit(funcAcc);
    TCanvas* LLAccCanvas = new TCanvas("LLAccCanvas","Distribution of -ln(likelihood) (Acc normalisation)");
    LLAccCanvas->cd();
    LLAcc->SetMinimum(Minimum(LLAcc)*0.999);
    LLAcc->SetMarkerColor(3);
    LLAcc->SetMarkerStyle(22);
    funcAcc->Draw();
    funcAcc->SetTitle(LLAcc->GetTitle());
    LLAcc->Draw("PSAME");
    LLAccCanvas->Write();
    funcAcc->Write();
  }

  OutputFile->Close();
}
