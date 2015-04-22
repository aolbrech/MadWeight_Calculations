#include "TH1.h"
#include "TF1.h"
#include "TList.h"
#include "TFile.h"
#include "TCanvas.h"

Bool_t reject;
Double_t fpol(Double_t *x, Double_t *par)
{
  if (reject && ( x[0] != 153 && x[0] != 163 && x[0] != 171 && x[0] != 172 && x[0] != 173 && x[0] != 174 && x[0] != 175 && x[0] != 183 && x[0] != 193) ){ 
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

Double_t Minimum(TH1F* histo)
{
  Double_t minimum = histo->GetBinContent(1);
  for(int ii = 0; ii < histo->GetNbinsX(); ii++){
    if( histo->GetBinContent(ii+1) < minimum && histo->GetBinContent(ii+1) > 0.0)
      minimum = histo->GetBinContent(ii+1);
  }
  return minimum;
}

void fitExcludeEmptyBins() {
  //Get the histogram from the ROOT file   --> Name gets changed by pythonScript!!
  TFile* InputFile = new TFile(" Events_ContainingMTopScansForRecLevelEvts/TopMass_Gen_ManyEvts/Likelihood_TopMass_Gen_UsingWeights.root ","r");
  TFile* OutputFile = new TFile(" Events_ContainingMTopScansForRecLevelEvts/TopMass_Gen_ManyEvts/LimitedFitResult_Gen.root ","RECREATE");

  bool GenLevel = true; 

  TH1F* LL    = (TH1F*) InputFile->Get("LL");
  TH1F* LLXS  = (TH1F*) InputFile->Get("LL_XS");
  TH1F* LLAcc = (TH1F*) InputFile->Get("LL_Acc");

  TF1 *func = new TF1("func",fpol,153,193,3); 
  func->SetName("fit_LL");
  reject = kTRUE;
  LL->Fit(func);
  reject = kFALSE;
  TCanvas* LLCanvas = new TCanvas("LLCanvas","Distribution of -ln(likelihood) (no normalisation)");
  LLCanvas->cd();
  LL->SetMinimum(Minimum(LL)*0.999);
  func->Draw();
  func->SetTitle(LL->GetTitle());
  LL->Draw("SAMEP");
  LLCanvas->Write();
  func->Write();

  TF1 *funcXS = new TF1("funcXS",fpol,153,193,3); 
  funcXS->SetName("fit_LLXS");
  reject = kTRUE;
  LLXS->Fit(funcXS);
  reject = kFALSE;
  TCanvas* LLXSCanvas = new TCanvas("LLXSCanvas","Distribution of -ln(likelihood) (XS normalisation)");
  LLXSCanvas->cd();
  LLXS->SetMinimum(Minimum(LLXS)*0.99);
  funcXS->Draw();
  funcXS->SetTitle(LLXS->GetTitle());
  LLXS->Draw("PSAME");
  LLXSCanvas->Write();
  funcXS->Write();

  if(GenLevel == false){
  TF1 *funcAcc = new TF1("funcAcc",fpol,153,193,3); 
    funcAcc->SetName("fit_LLAcc");
    reject = kTRUE;
    LLAcc->Fit(funcAcc);
    reject = kFALSE;
    TCanvas* LLAccCanvas = new TCanvas("LLAccCanvas","Distribution of -ln(likelihood) (Acc normalisation)");
    LLAccCanvas->cd();
    LLAcc->SetMinimum(Minimum(LLAcc)*0.99);
    funcAcc->Draw();
    funcAcc->SetTitle(LLAcc->GetTitle());
    LLAcc->Draw("PSAME");
    LLAccCanvas->Write();
    funcAcc->Write();
  }

  OutputFile->Close();
}
