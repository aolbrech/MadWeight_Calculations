//-----------------------------------------------------------------------------------------------------------------------------//
//                                                                                                                             //
//  Create a script which will compare the measured minimum position with the variable used for creating the MadGraph sample   //
//                                                                                                                             //
//-----------------------------------------------------------------------------------------------------------------------------//

#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <iostream>
#include "TCanvas.h"

TH2F* h_MinCurve_FstPol = new TH2F("MinimumCurveFirstPol","Measured minimum versus used value of MG sample (using pol2 on all points)",100,-0.5,0.5,100,-0.5,0.5);
TH2F* h_MinCurve_ScdPol = new TH2F("MinimumCurveSecondPol","Measured minimum versus used value of MG sample (using pol2 on reduced points)",100,-0.5,0.5,100,-0.5,0.5);
TH2F* h_CalCurve_FstPol = new TH2F("CalibrationCurveFirstPol","Deviation from expected minimum versus used value of MG sample (using pol2 on all points)",100,-0.5,0.5,100,-0.5,0.5);
TH2F* h_CalCurve_ScdPol = new TH2F("CalibrationCurveSecondPol","Deviation from expected minimum versus used value of MG sample (using pol2 on reduced points)",100,-0.5,0.5,100,-0.5,0.5);

//Make a straight line which indicates where you expect the points!
double x[9] = {-0.3, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2, 0.3};
double yMin[9] = {-0.3, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2, 0.3};
double yCal[9] = {0.0};
TGraph* gr_LineMin = new TGraph(9, x, yMin); 
TGraph* gr_LineCal = new TGraph(9, x, yCal);

void getMinimum(TH1F* h_Fst, TH1F* h_Scd, float rangeVal, int range, double MinFst[], double MinScd[], double eMinFst[], double eMinScd[]){
  
  TF1* polFit_Fst = new TF1("polFit_Fst","pol2",-0.5,0.5); 
  h_Fst->Fit(polFit_Fst,"Q");

  MinFst[range] = polFit_Fst->GetMinimumX();
  double Neg1Sigma = polFit_Fst->GetX(h_Fst->GetBinContent(h_Fst->FindBin(polFit_Fst->GetMinimumX()))+0.5,-0.5,polFit_Fst->GetMinimumX());
  double Pos1Sigma = polFit_Fst->GetX(h_Fst->GetBinContent(h_Fst->FindBin(polFit_Fst->GetMinimumX()))+0.5,polFit_Fst->GetMinimumX(), 0.5);
  eMinFst[range] = Pos1Sigma - Neg1Sigma;
  
  //h_MinCurve_FstPol->Fill(rangeVal, polFit_Fst->GetMinimumX());
  //Find 1sigma interval around minimum (~ y-axis difference of 0.5!)
  //h_CalCurve_FstPol->Fill(rangeVal, polFit_Fst->GetMinimumX() - rangeVal);

  //TF1* polFit_Scd = new TF1("polFit_Scd","pol2",-0.5,0.5); 
  //h_Scd->Fit(polFit_Scd,"Q");

  //h_MinCurve_ScdPol->Fill(rangeVal, polFit_Scd->GetMinimumX());
  //h_CalCurve_ScdPol->Fill(rangeVal, polFit_Scd->GetMinimumX() - rangeVal);
}

void CalibrationCurve(){

  double Range[] = {-0.3, -0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.3};//, 0.0}; //-0.2}; //-0.3, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2, 0.3};
  std::string RangeName[] = {"Neg03","Neg02", "Neg01","Neg005", "SM", "Pos005", "Pos01", "Pos02", "Pos03"};//,"SM"}; //"Neg03","Neg02","Neg01","Neg005","SM","Pos005","Pos01","Pos02","Pos03"};
  int sizeRange = sizeof(Range)/sizeof(Range[0]);

  //Initialize the arrays for the TGraphErrors:
  double min_Fst[sizeRange], min_Scd[sizeRange], cal_Fst[sizeRange], cal_Scd[sizeRange], calRel_Fst[sizeRange], calRel_Scd[sizeRange];
  double Emin_Fst[sizeRange], Emin_Scd[sizeRange], Ecal_Fst[sizeRange], Ecal_Scd[sizeRange], EcalRel_Fst[sizeRange], EcalRel_Scd[sizeRange];

  //Set the correct directory where all these files are saved:
  std::string whichDir = "Events_RgRScan_AllDeltaTF_EDep";
  std::string cuts = "AllCuts";

  //Now loop over all the configurations that need to be considered and get the correct distributions
  for( int iRange = 0; iRange < sizeRange; iRange++){
    TFile* file = new TFile((whichDir+"/AllDeltaTF_RgRMGSample"+RangeName[iRange]+"_2500Evts_NarrowWideRange_"+cuts+"/FitOptimizations_AllDeltaTF_RgRMGSample"+RangeName[iRange]+"_2500Evts_NarrowWideRange_"+cuts+"_2500Evts.root").c_str(),"READ");
  
    //Now get the relevant histograms from this file:
    TH1F* h_FirstPol = (TH1F*) file->Get("OriginalDistributions/FirstPolAcc_Summed");
    TH1F* h_SecondPol = (TH1F*) file->Get("OriginalDistributions/SecondPolAcc_Summed");

    //Now send this to the class which gets the minimum
    getMinimum(h_FirstPol, h_SecondPol, Range[iRange], iRange, min_Fst, min_Scd, Emin_Fst, Emin_Scd);

    //Now fill the other arrays:
    cal_Fst[iRange] = min_Fst[iRange] - Range[iRange];  Ecal_Fst[iRange] = Emin_Fst[iRange];
    cal_Scd[iRange] = min_Scd[iRange] - Range[iRange];  Ecal_Scd[iRange] = Emin_Scd[iRange];
    calRel_Fst[iRange] = cal_Fst[iRange]/Range[iRange]; EcalRel_Fst[iRange] = Emin_Fst[iRange]/Range[iRange];
    calRel_Scd[iRange] = cal_Scd[iRange]/Range[iRange]; EcalRel_Scd[iRange] = Emin_Scd[iRange]/Range[iRange];
  }

  //Make the TGraphError!
  TGraphErrors *gr_MinFst = new TGraphErrors(sizeRange, Range, min_Fst, 0, Emin_Fst);
  TGraphErrors *gr_CalFst = new TGraphErrors(sizeRange, Range, cal_Fst, 0, Ecal_Fst);

  TCanvas* canv_MinFst = new TCanvas("MinimumCurve_FirstFit","MinimumCurve_FirstFit"); canv_MinFst->cd();
  h_MinCurve_FstPol->GetXaxis()->SetTitle("gR value used for MG sample generation");
  h_MinCurve_FstPol->GetYaxis()->SetTitle("Minimum obtained from fit (using all points)");
  //h_MinCurve_FstPol->SetMarkerStyle(33); h_MinCurve_FstPol->Draw();
  gr_MinFst->Draw("A*");
  gr_LineMin->SetLineColor(3); gr_LineMin->Draw("C"); 
  canv_MinFst->SaveAs((whichDir+"/MinimumCurveFirstPol_"+cuts+".pdf").c_str());	

  TCanvas* canv_CalFst = new TCanvas("CalibrationCurve_FirstFit","CalibrationCurve_FirstFit"); canv_CalFst->cd();
  h_CalCurve_FstPol->GetXaxis()->SetTitle("gR value used for MG sample generation");
  h_CalCurve_FstPol->GetYaxis()->SetTitle("Fit Minimum - expected minimum (using all points)");
  //h_CalCurve_FstPol->SetMarkerStyle(33); h_CalCurve_FstPol->Draw();
  gr_CalFst->Draw("A*");
  gr_LineCal->SetLineColor(3); gr_LineCal->Draw("C"); 
  canv_CalFst->SaveAs((whichDir+"/CalibrationCurveFirstPol_"+cuts+".pdf").c_str());	

/*  TCanvas* canv_MinScd = new TCanvas("MinimumCurve_SecondFit","MinimumCurve_SecondFit"); canv_MinScd->cd();
  h_MinCurve_ScdPol->GetXaxis()->SetTitle("gR value used for MG sample generation");
  h_MinCurve_ScdPol->GetYaxis()->SetTitle("Minimum obtained from fit (using reduced points)");
  h_MinCurve_ScdPol->SetMarkerStyle(33); h_MinCurve_ScdPol->Draw();
  gr_LineMin->SetLineColor(3); gr_LineMin->Draw("C"); 
  canv_MinScd->SaveAs((whichDir+"/MinimumCurveSecondPol_"+cuts+".pdf").c_str());	

  TCanvas* canv_CalScd = new TCanvas("CalibrationCurve_SecondFit","CalibrationCurve_SecondFit"); canv_CalScd->cd();
  h_CalCurve_ScdPol->GetXaxis()->SetTitle("gR value used for MG sample generation");
  h_CalCurve_ScdPol->GetYaxis()->SetTitle("Fit Minimum - expected minimum (using reduced points)");
  h_CalCurve_ScdPol->SetMarkerStyle(33); h_CalCurve_ScdPol->Draw();
  gr_LineCal->SetLineColor(3); gr_LineCal->Draw("C"); 
  canv_CalScd->SaveAs((whichDir+"/CalibrationCurveSecondPol_"+cuts+".pdf").c_str());	
*/
}
