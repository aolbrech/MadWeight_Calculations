#include <iostream>
#include <fstream>
#include <sstream>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TDirectory.h"
#include <string>
#include <map>
#include <cmath>
#include <algorithm>

std::string VarValues[] = {"Re(V_{R}) = -0.3","Re(V_{R}) = -0.2","Re(V_{R}) = -0.1","Re(V_{R}) = -0.05","Re(V_{R}) = 0.0","Re(V_{R}) = 0.05","Re(V_{R}) = 0.1","Re(V_{R}) = 0.2","Re(V_{R}) = 0.3"}; 
float Var[] = {-0.3,-0.2,-0.1,-0.05,0.0,0.05,0.1,0.2,0.3}; 
float MGXS[] = {13.3944,12.06555,11.25909,11.02784,10.90059,10.88228,10.97767,11.49883,12.49056}; 
float MGXSCut[] = {13.3944,12.06555,11.25909,11.02784,10.90059,10.88228,10.97767,11.49883,12.49056}; 
int xBin = 13; 
float xLow = -0.325; 
float xHigh = 0.325; 
int xMinValue[] = {4,4,10}; 
std::string KinVar = "Re(V_{R})"; 
int VarWindow = 2; 
int xPos[] = {5,6}; 
int xNeg[] = {3,2}; 
std::string title = "Gen_RVR"; 
TFile* Tfile = new TFile("Test.root","RECREATE");

const int NrConfigs = 9;
const int nEvts = 10; 
const int TotalNrEvts = 10000;

void ReadTest(){
  
  int xMin = xMinValue[VarWindow-1];
  float xStep[] = {Var[xNeg[0]]-Var[xNeg[1]], Var[xMin]-Var[xNeg[0]], Var[xPos[0]]-Var[xMin], Var[xPos[1]]-Var[xPos[0]] };
  std::cout << " Value of xStep is : " << xStep[0] << std::endl;

  //--- Initialize histograms ---//
  TH1F* LnLikAll    = new TH1F("LnLikAll",   ("Distribution of -ln(L) using all events (no norm -- "+title+" evts)").c_str(), xBin,xLow,xHigh);
  TH1F* LnLikXSAll  = new TH1F("LnLikXSAll", ("Distribution of -ln(L) using all events (XS norm -- "+title+" evts)").c_str(), xBin,xLow,xHigh);
  TH1F* LnLikAccAll = new TH1F("LnLikAccAll",("Distribution of -ln(L) using all events (Acc norm -- "+title+" evts)").c_str(),xBin,xLow,xHigh);

  TH1F* YPlusGausTest   =new TH1F("YPlusGausTest",   ("Comparison of fit deviation in "+KinVar+" = "+Var[xPos[1]]+" and "+KinVar+" = "+Var[xPos[0]]+" (no norm -- "+title+" evts)").c_str(), 250,-0.15,0.15);
  TH1F* YPlusGausTestXS =new TH1F("YPlusGausTestXS", ("Comparison of fit deviation in "+KinVar+" = "+Var[xPos[1]]+" and "+KinVar+" = "+Var[xPos[0]]+" (XS norm -- "+title+" evts)").c_str(), 250,-0.15,0.15);
  TH1F* YPlusGausTestAcc=new TH1F("YPlusGausTestAcc",("Comparison of fit deviation in "+KinVar+" = "+Var[xPos[1]]+" and "+KinVar+" = "+Var[xPos[0]]+" (Acc norm -- "+title+" evts)").c_str(),250,-0.15,0.15);
  TH1F* YPlusGausTestPosScdDer    = new TH1F("YPlusGausTestPosScdDer",   ("Gaussian fit deviation using "+KinVar+" = "+Var[xPos[1]]+" & "+Var[xPos[0]]+" (no norm -- outer 2nd der > 0 -- "+title+" evts)").c_str(), 250,-0.1,0.1);
  TH1F* YPlusGausTestXSPosScdDer  = new TH1F("YPlusGausTestXSPosScdDer", ("Gaussian fit deviation using "+KinVar+" = "+Var[xPos[1]]+" & "+Var[xPos[0]]+" (XS norm -- outer 2nd der > 0 -- "+title+" evts)").c_str(), 250,-0.1,0.1);
  TH1F* YPlusGausTestAccPosScdDer = new TH1F("YPlusGausTestAccPosScdDer",("Gaussian fit deviation using "+KinVar+" = "+Var[xPos[1]]+" & "+Var[xPos[0]]+" (Acc norm -- outer 2nd der > 0 -- "+title+" evts)").c_str(),250,-0.1,0.1);
  TH1F* YPlusGausTestNegScdDer    = new TH1F("YPlusGausTestNegScdDer",   ("Gaussian fit deviation using "+KinVar+" = "+Var[xPos[1]]+" & "+Var[xPos[0]]+" (no norm -- outer 2nd der < 0 -- "+title+" evts)").c_str(), 250,-0.1,0.1);
  TH1F* YPlusGausTestXSNegScdDer  = new TH1F("YPlusGausTestXSNegScdDer", ("Gaussing fit deviation using "+KinVar+" = "+Var[xPos[1]]+" & "+Var[xPos[0]]+" (XS norm -- outer 2nd der < 0 -- "+title+" evts)").c_str(), 250,-0.1,0.1);
  TH1F* YPlusGausTestAccNegScdDer = new TH1F("YPlusGausTestAccNegScdDer",("Gaussian fit deviation using "+KinVar+" = "+Var[xPos[1]]+" & "+Var[xPos[0]]+" (Acc norm -- outer 2nd der < 0 -- "+title+" evts)").c_str(),250,-0.1,0.1);

  TH1F* YPlus    = new TH1F("YPlus",   ("Deviation from parabolic fit for "+KinVar+" = "+Var[xPos[0]]+" (no norm -- "+title+" evts)").c_str() ,150,-0.25,0.25);
  TH1F* YPlusXS  = new TH1F("YPlusXS", ("Deviation from parabolic fit for "+KinVar+" = "+Var[xPos[0]]+" (XS norm -- "+title+" evts)").c_str() ,150,-0.25,0.25);
  TH1F* YPlusAcc = new TH1F("YPlusAcc",("Deviation from parabolic fit for "+KinVar+" = "+Var[xPos[0]]+" (Acc norm -- "+title+" evts)").c_str(),150,-0.25,0.25);
  TH1F* YPlusPosScdDer    = new TH1F("YPlusPosScdDer",   ("Deviation from parabolic fit for "+KinVar+" = "+Var[xPos[0]]+" (no norm -- outer 2nd der > 0 -- "+title+" evts)").c_str() ,150,-0.2,0.2);
  TH1F* YPlusXSPosScdDer  = new TH1F("YPlusXSPosScdDer", ("Deviation from parabolic fit for "+KinVar+" = "+Var[xPos[0]]+" (XS norm -- outer 2nd der > 0 -- "+title+" evts)").c_str() ,150,-0.2,0.2);
  TH1F* YPlusAccPosScdDer = new TH1F("YPlusAccPosScdDer",("Deviation from parabolic fit for "+KinVar+" = "+Var[xPos[0]]+" (Acc norm -- outer 2nd der > 0 -- "+title+" evts)").c_str(),150,-0.2,0.2);
  TH1F* YPlusNegScdDer    = new TH1F("YPlusNegScdDer",   ("Deviation from parabolic fit for "+KinVar+" = "+Var[xPos[0]]+" (no norm -- outer 2nd der < 0 -- "+title+" evts)").c_str() ,150,-0.2,0.2);
  TH1F* YPlusXSNegScdDer  = new TH1F("YPlusXSNegScdDer", ("Deviation from parabolic fit for "+KinVar+" = "+Var[xPos[0]]+" (XS norm -- outer 2nd der < 0 -- "+title+" evts)").c_str() ,150,-0.2,0.2);
  TH1F* YPlusAccNegScdDer = new TH1F("YPlusAccNegScdDer",("Deviation from parabolic fit for "+KinVar+" = "+Var[xPos[0]]+" (Acc norm -- outer 2nd der < 0 -- "+title+" evts)").c_str(),150,-0.2,0.2);

  TH1F* YPlusPlus    = new TH1F("YPlusPlus",   ("Deviation from parabolic fit for "+KinVar+" = "+Var[xPos[1]]+" (no norm -- "+title+" evts)").c_str() ,150,-0.5,0.5);
  TH1F*  YPlusPlusXS  = new TH1F("YPlusPlusXS", ("Deviation from parabolic fit for "+KinVar+" = "+Var[xPos[1]]+" (XS norm -- "+title+" evts)").c_str() ,150,-0.5,0.5);
  TH1F* YPlusPlusAcc = new TH1F("YPlusPlusAcc",("Deviation from parabolic fit for "+KinVar+" = "+Var[xPos[1]]+" (Acc norm -- "+title+" evts)").c_str(),150,-0.5,0.5);
  TH1F*  YPlusPlusPosScdDer    = new TH1F("YPlusPlusPosScdDer",   ("Deviation from parabolic fit for "+KinVar+" = "+Var[xPos[1]]+" (no norm -- outer 2nd der > 0 -- "+title+" evts)").c_str() ,150,-0.3,0.3);
  TH1F*  YPlusPlusXSPosScdDer  = new TH1F("YPlusPlusXSPosScdDer", ("Deviation from parabolic fit for "+KinVar+" = "+Var[xPos[1]]+" (XS norm -- outer 2nd der > 0 -- "+title+" evts)").c_str() ,150,-0.3,0.3);
  TH1F*  YPlusPlusAccPosScdDer = new TH1F("YPlusPlusAccPosScdDer",("Deviation from parabolic fit for "+KinVar+" = "+Var[xPos[1]]+" (Acc norm -- outer 2nd der > 0 -- "+title+" evts)").c_str(),150,-0.3,0.3);
  TH1F*  YPlusPlusNegScdDer    = new TH1F("YPlusPlusNegScdDer",   ("Deviation from parabolic fit for "+KinVar+" = "+Var[xPos[1]]+" (no norm -- outer 2nd der < 0 -- "+title+" evts)").c_str() ,150,-0.3,0.3);
  TH1F*  YPlusPlusXSNegScdDer  = new TH1F("YPlusPlusXSNegScdDer", ("Deviation from parabolic fit for "+KinVar+" = "+Var[xPos[1]]+" (XS norm -- outer 2nd der < 0 -- "+title+" evts)").c_str() ,150,-0.3,0.3);
  TH1F*  YPlusPlusAccNegScdDer = new TH1F("YPlusPlusAccNegScdDer",("Deviation from parabolic fit for "+KinVar+" = "+Var[xPos[1]]+" (Acc norm -- outer 2nd der < 0 -- "+title+" evts)").c_str(),150,-0.3,0.3);

  TH1F* YMin    = new TH1F("YMin",   ("Deviation from parabolic fit for "+KinVar+" = "+Var[xNeg[0]]+" (no norm -- "+title+" evts)").c_str() ,150,-0.25,0.25);
  TH1F* YMinXS  = new TH1F("YMinXS", ("Deviation from parabolic fit for "+KinVar+" = "+Var[xNeg[0]]+" (XS norm -- "+title+" evts)").c_str() ,150,-0.25,0.25);
  TH1F* YMinAcc = new TH1F("YMinAcc",("Deviation from parabolic fit for "+KinVar+" = "+Var[xNeg[0]]+" (Acc norm -- "+title+" evts)").c_str(),150,-0.25,0.25);
  TH1F* YMinPosScdDer    = new TH1F("YMinPosScdDer",   ("Deviation from parabolic fit for "+KinVar+" = "+Var[xNeg[0]]+" (no norm -- outer 2nd der > 0 -- "+title+" evts)").c_str() ,150,-0.2,0.2);
  TH1F* YMinXSPosScdDer  = new TH1F("YMinXSPosScdDer", ("Deviation from parabolic fit for "+KinVar+" = "+Var[xNeg[0]]+" (XS norm -- outer 2nd der > 0 -- "+title+" evts)").c_str() ,150,-0.2,0.2);
  TH1F* YMinAccPosScdDer = new TH1F("YMinAccPosScdDer",("Deviation from parabolic fit for "+KinVar+" = "+Var[xNeg[0]]+" (Acc norm -- outer 2nd der > 0 -- "+title+" evts)").c_str(),150,-0.2,0.2);
  TH1F* YMinNegScdDer    = new TH1F("YMinNegScdDer",   ("Deviation from parabolic fit for "+KinVar+" = "+Var[xNeg[0]]+" (no norm -- outer 2nd der < 0 -- "+title+" evts)").c_str() ,150,-0.2,0.2);
  TH1F* YMinXSNegScdDer  = new TH1F("YMinXSNegScdDer", ("Deviation from parabolic fit for "+KinVar+" = "+Var[xNeg[0]]+" (XS norm -- outer 2nd der < 0 -- "+title+" evts)").c_str() ,150,-0.2,0.2);
  TH1F* YMinAccNegScdDer = new TH1F("YMinAccNegScdDer",("Deviation from parabolic fit for "+KinVar+" = "+Var[xNeg[0]]+" (Acc norm -- outer 2nd der < 0 -- "+title+" evts)").c_str(),150,-0.2,0.2);
  
  TH1F* YMinMin    = new TH1F("YMinMin",   ("Deviation from parabolic fit for "+KinVar+" = "+Var[xNeg[1]]+" (no norm -- "+title+" evts)").c_str() ,150,-0.5,0.5);
  TH1F* YMinMinXS  = new TH1F("YMinMinXS", ("Deviation from parabolic fit for "+KinVar+" = "+Var[xNeg[1]]+" (XS norm -- "+title+" evts)").c_str() ,150,-0.5,0.5);
  TH1F* YMinMinAcc = new TH1F("YMinMinAcc",("Deviation from parabolic fit for "+KinVar+" = "+Var[xNeg[1]]+" (Acc norm -- "+title+" evts)").c_str(),150,-0.5,0.5);
  TH1F* YMinMinPosScdDer    = new TH1F("YMinMinPosScdDer",   ("Deviation from parabolic fit for "+KinVar+" = "+Var[xNeg[1]]+" (no norm -- outer 2nd der > 0 -- "+title+" evts)").c_str() ,150,-0.3,0.3);
  TH1F* YMinMinXSPosScdDer  = new TH1F("YMinMinXSPosScdDer", ("Deviation from parabolic fit for "+KinVar+" = "+Var[xNeg[1]]+" (XS norm -- outer 2nd der > 0 -- "+title+" evts)").c_str() ,150,-0.3,0.3);
  TH1F* YMinMinAccPosScdDer = new TH1F("YMinMinAccPosScdDer",("Deviation from parabolic fit for "+KinVar+" = "+Var[xNeg[1]]+" (Acc norm -- outer 2nd der > 0 -- "+title+" evts)").c_str(),150,-0.3,0.3);
  TH1F* YMinMinNegScdDer    = new TH1F("YMinMinNegScdDer",   ("Deviation from parabolic fit for "+KinVar+" = "+Var[xNeg[1]]+" (no norm -- outer 2nd der < 0 -- "+title+" evts)").c_str() ,150,-0.3,0.3);
  TH1F* YMinMinXSNegScdDer  = new TH1F("YMinMinXSNegScdDer", ("Deviation from parabolic fit for "+KinVar+" = "+Var[xNeg[1]]+" (XS norm -- outer 2nd der < 0 -- "+title+" evts)").c_str() ,150,-0.3,0.3);
  TH1F* YMinMinAccNegScdDer = new TH1F("YMinMinAccNegScdDer",("Deviation from parabolic fit for "+KinVar+" = "+Var[xNeg[1]]+" (Acc norm -- outer 2nd der < 0 -- "+title+" evts)").c_str(),150,-0.3,0.3);
  
  TH1F* YRelPlus    = new TH1F("YRelPlus",   ("Relative deviation from parabolic fit for "+KinVar+" = "+Var[xNeg[0]]+" (no norm -- "+title+" evts)").c_str() ,250,-0.1,0.1);
  TH1F* YRelPlusXS  = new TH1F("YRelPlusXS", ("Relative deviation from parabolic fit for "+KinVar+" = "+Var[xNeg[0]]+" (XS norm -- "+title+" evts)").c_str() ,250,-0.1,0.1);
  TH1F* YRelPlusAcc = new TH1F("YRelPlusAcc",("Relative deviation from parabolic fit for "+KinVar+" = "+Var[xNeg[0]]+" (Acc norm -- "+title+" evts)").c_str(),250,-0.1,0.1);
  
  TH1F* YRelMin    = new TH1F("YRelMin",   ("Relative deviation from parabolic fit for "+KinVar+" = "+Var[xNeg[0]]+" (no norm -- "+title+" evts)").c_str() ,150,-0.1,0.1);
  TH1F* YRelMinXS  = new TH1F("YRelMinXS", ("Relative deviation from parabolic fit for "+KinVar+" = "+Var[xNeg[0]]+" (XS norm -- "+title+" evts)").c_str() ,150,-0.1,0.1);
  TH1F* YRelMinAcc = new TH1F("YRelMinAcc",("Relative deviation from parabolic fit for "+KinVar+" = "+Var[xNeg[0]]+" (Acc norm -- "+title+" evts)").c_str(),150,-0.1,0.1);
  
  TH1F* FstDer    = new TH1F("FirstDer",   ("First derivative of -ln(likelihood) distribution -- "+title+" evts").c_str(), 4,0,4);
  //std::string LabelOne = ("(yDATA(x="+Var[xNeg[1]]+") - yDATA(x="+Var[xNeg[0]]+"))/"+xStep[0]).c_str();
  //(FstDer->GetXaxis())->SetBinLabel(1,("(yDATA(x="+Var[xNeg[1]]+") - yDATA(x="+Var[xNeg[0]]+"))/"+xStep[0]));
  //FstDer->GetXaxis()->SetBinLabel(2,"bin2 test"); //("(y_{DATA}(x="+Var[xNeg[0]]+") - y_{DATA}(x="+Var[xMin]+"))/"+xStep[1]).c_str());
  //FstDer->GetXaxis()->SetBinLabel(3,"bin3 test"); //("(y_{DATA}(x="+Var[xMin]+") - y_{DATA}(x="+Var[xPos[0]]+"))/"+xStep[2]).c_str());
  //FstDer->GetXaxis()->SetBinLabel(4,"bin4 test"); //("(y_{DATA}(x="+Var[xPos[0]]+") - y_{DATA}(x="+Var[xPos[1]]+"))/"+xStep[3]).c_str());
  //TH1F* FstDerXS  = new TH1F("FirstDerivativeXS", ("First derivative of -ln(likelihood) distribution (XS norm -- "+title+" evts)").c_str(), 5,-0.25,0.25);
  //TH1F* FstDerAcc = new TH1F("FirstDerivativeAcc",("First derivative of -ln(likelihood) distribution (Acc norm -- "+title+" evts)").c_str(),5,-0.25,0.25);
  TH1F* ScdDerInner    = new TH1F("SecondDerivativeInner",   ("Second derivative of -ln(likelihood) distribution (no norm -- using inner points -- "+title+" evts)").c_str(), 250,-5,5);
  TH1F* ScdDerXSInner  = new TH1F("SecondDerivativeXSInner", ("Second derivative of -ln(likelihood) distribution (XS norm -- using inner points -- "+title+" evts)").c_str(), 250,-5,5);
  TH1F* ScdDerAccInner = new TH1F("SecondDerivativeAccInner",("Second derivative of -ln(likelihood) distribution (Acc norm -- using inner points -- "+title+" evts)").c_str(),250,-5,5);
  TH1F* ScdDerOuter    = new TH1F("SecondDerivativeOuter",   ("Second derivative of -ln(likelihood) distribution (no norm -- using outer points -- "+title+" evts)").c_str(), 250,-5,5);
  TH1F* ScdDerXSOuter  = new TH1F("SecondDerivativeXSOuter", ("Second derivative of -ln(likelihood) distribution (XS norm -- using outer points -- "+title+" evts)").c_str(), 250,-5,5);
  TH1F* ScdDerAccOuter = new TH1F("SecondDerivativeAccOuter",("Second derivative of -ln(likelihood) distribution (Acc norm -- using outer points -- "+title+" evts)").c_str(),250,-5,5);
  TH2F* ScdDerScatter    = new TH2F("ScdDerScatterPlot",   ("Second derivative of -ln(L) using inner points versus using outer points (no norm -- "+title+" evts)").c_str(), 250,-5,5,250,-5,5);
  TH2F* ScdDerXSScatter  = new TH2F("ScdDerXSScatterPlot", ("Second derivative of -ln(L) using inner points versus using outer points (XS norm -- "+title+" evts)").c_str(), 250,-5,5,250,-5,5);
  TH2F* ScdDerAccScatter = new TH2F("ScdDerAccScatterPlot",("Second derivative of -ln(L) using inner points versus using outer points (Acc norm -- "+title+" evts)").c_str(),250,-5,5,250,-5,5);
  
  TH1F* ProcContrDist = new TH1F("ProcentualContribution","Procentual contribution of each of the events to the total likelihood",250,0,1);
  
  TH1F* TotalFitDevDist    = new TH1F("TotalFitDeviation",   ("Sum of difference between likelihood and fit value for each point in range (no norm -- "+title+" events)").c_str(), 250,0,5);
  TH1F* TotalFitDevXSDist  = new TH1F("TotalFitXSDeviation", ("Sum of difference between likelihood and fit value for each point in range (XS norm -- "+title+" events)").c_str(), 250,0,5);
  TH1F* TotalFitDevAccDist = new TH1F("TotalFitAccDeviation",("Sum of difference between likelihood and fit value for each point in range (Acc norm -- "+title+" events)").c_str(),250,0,5);
  TH1F* TotalFctDevDist    = new TH1F("TotalFctDeviation",   ("Sum of difference between likelihood and function value for each point in range (no norm -- "+title+" events)").c_str(), 250,0,5);
  TH1F* TotalFctDevXSDist  = new TH1F("TotalFctXSDeviation", ("Sum of difference between likelihood and function value for each point in range (XS norm -- "+title+" events)").c_str(), 250,0,5);
  TH1F* TotalFctDevAccDist = new TH1F("TotalFctAccDeviation",("Sum of difference between likelihood and function value for each point in range (Acc norm -- "+title+" events)").c_str(),250,0,5);

  TH1F* LnLikDist = new TH1F("LnLik","title",xBin,xLow,xHigh);
  LnLikDist->SetMarkerStyle(20); LnLikDist->SetLineColor(1); LnLikDist->SetMarkerColor(1); LnLikDist->SetMarkerSize(1.2);
  TH1F* LnLikXSDist = new TH1F("LnLikXS","title",xBin,xLow,xHigh);
  LnLikXSDist->SetMarkerStyle(21); LnLikXSDist->SetLineColor(3); LnLikXSDist->SetMarkerColor(3); LnLikXSDist->SetMarkerSize(1.2);
  TH1F* LnLikAccDist = new TH1F("LnLikAcc","title",xBin,xLow,xHigh);
  LnLikAccDist->SetMarkerStyle(22); LnLikAccDist->SetLineColor(4); LnLikAccDist->SetMarkerColor(4); LnLikAccDist->SetMarkerSize(1.2);
  
  TDirectory* FitComp = Tfile->mkdir("FitComparison");
  TDirectory* LnLikDir = Tfile->mkdir("LnLikDist");
  TDirectory* LnLikStackDir = LnLikDir->mkdir("LnLikStackDist");
  TDirectory* LnLikXSDir = Tfile->mkdir("LnLikXSDist");
  TDirectory* LnLikAccDir = Tfile->mkdir("LnLikAccDist");
  TDirectory* LnLikAccDirVarVsUnc = Tfile->mkdir("LnLikAccDist_VarLargerThanAvgUnc");
  TDirectory* LnLikAccDirVarVsDUnc = Tfile->mkdir("LnLikAccDist_VarLargerThanTwiceAvgUnc");
  TDirectory* FstDerDir = Tfile->mkdir("FirstDerivativeDist");

  //-----------------------------------//
  //----  Specific cut histograms  ----//
  //-----------------------------------//
  TH1F* LLSmallFctDevDist = new TH1F("LLSmallFctDeviation","-ln(L) when sum of deviations between -ln(L)- and function-value < 0.5 (acc norm -- outer points used)",xBin,xLow,xHigh);

  TH1F* LLPosScdDerDistBoth     = new TH1F("LLPosScdDerBoth",    "-ln(L) when both 2nd derivatives > 0 (no norm -- ", xBin,xLow,xHigh);
  TH1F* LLXSPosScdDerDistBoth   = new TH1F("LLXSPosScdDerBoth",  "-ln(L) when both 2nd derivatives > 0 (XS norm -- ", xBin,xLow,xHigh);
  TH1F* LLAccPosScdDerDistBoth  = new TH1F("LLAccPosScdDerBoth", "-ln(L) when both 2nd derivatives > 0 (Acc norm -- ",xBin,xLow,xHigh);
  TH1F* LLPosScdDerDistInner    = new TH1F("LLPosScdDerInner",   "-ln(L) when inner 2nd derivative > 0 (no norm -- ", xBin,xLow,xHigh);
  TH1F* LLXSPosScdDerDistInner  = new TH1F("LLXSPosScdDerInner", "-ln(L) when inner 2nd derivative > 0 (XS norm -- ", xBin,xLow,xHigh);
  TH1F* LLAccPosScdDerDistInner = new TH1F("LLAccPosScdDerInner","-ln(L) when inner 2nd derivative > 0 (Acc norm -- ",xBin,xLow,xHigh);
  TH1F* LLPosScdDerDistOuter    = new TH1F("LLPosScdDerOuter",   "-ln(L) when outer 2nd derivative > 0 (no norm -- ", xBin,xLow,xHigh);
  TH1F* LLXSPosScdDerDistOuter  = new TH1F("LLXSPosScdDerOuter", "-ln(L) when outer 2nd derivative > 0 (XS norm -- ", xBin,xLow,xHigh);
  TH1F* LLAccPosScdDerDistOuter = new TH1F("LLAccPosScdDerOuter","-ln(L) when outer 2nd derivative > 0 (Acc norm -- ",xBin,xLow,xHigh);

  //Initialize the event counters:
  vector<int> EvtsWithYPlusGausSmall, EvtsWithYPlusGausSmallXS, EvtsWithYPlusGausSmallAcc;
  int EvtsWithPosScdDerInner = 0, EvtsWithPosScdDerXSInner = 0, EvtsWithPosScdDerAccInner = 0;
  int EvtsWithPosScdDerOuter = 0, EvtsWithPosScdDerXSOuter = 0, EvtsWithPosScdDerAccOuter = 0;
  int EvtsWithPosScdDerBoth  = 0, EvtsWithPosScdDerXSBoth  = 0, EvtsWithPosScdDerAccBoth  = 0;
  //vector<int> EvtsPosScdDerInner, EvtsPosScdDerXSInner, EvtsPosScdDerAccInner;
  //vector<int> EvtsPosScdDerOuter, EvtsPosScdDerXSOuter, EvtsPosScdDerAccOuter;
  vector<int> EvtsWithSmallFctDev;
  //double weightsValue[nEvts][NrConfigs];
  int consEvts = 0;
  int NrCanvas = 0, xDivide = 3, yDivide = 3;
  TCanvas* StackCanvas = new TCanvas("StackCanvas_Nr0","StackedCanvas");
  StackCanvas->Divide(xDivide,yDivide);
  TGraph* LnLikFctOuter[TotalNrEvts];

  //--- Read all likelihood values ! ---//
  std::ifstream ifs ("Events/RVR_Gen_SingleGausTFHalfWidth_10000Evts/weights.out", std::ifstream::in); 
  std::cout << " Value of ifs : " << ifs.eof() << std::endl;
  std::string line;
  int evt,config,tf;
  double weight, weightUnc;
  while( std::getline(ifs,line) && consEvts < nEvts){
    std::istringstream iss(line);
    if( iss >> evt >> config >> tf >> weight >> weightUnc){
    
      if(config == 1){ std::cout << " Looking at event : " << evt << std::endl; consEvts++;}

      stringstream ssEvt; ssEvt << evt; string sEvt = ssEvt.str();
      //--- Initialize the event-per-event variables! ---//
      if( config == 1){
        double aHat[2] = {0.0, 0.0}, aHatXS[2] = {0.0, 0.0}, aHatAcc[2] = {0.0, 0.0};
        double bHat[2] = {0.0, 0.0}, bHatXS[2] = {0.0, 0.0}, bHatAcc[2] = {0.0, 0.0};
        double cHat[2] = {0.0, 0.0}, cHatXS[2] = {0.0, 0.0}, cHatAcc[2] = {0.0, 0.0};
        double LnLik[NrConfigs] = {0.0}, LnLikXS[NrConfigs] = {0.0}, LnLikAcc[NrConfigs] = {0.0};
        //double Lik[NrConfigs] = {0.0},   LikXS[NrConfigs] = {0.0},   LikAcc[NrConfigs] = {0.0};

        LnLikDist->SetName(("LnLik_Evt"+sEvt).c_str());       LnLikDist->SetTitle(("LnLik distribution for event "+sEvt+" -- "+title+" evts").c_str());
        LnLikXSDist->SetName(("LnLikXS_Evt"+sEvt).c_str());   LnLikXSDist->SetTitle(("LnLikXS distribution for event "+sEvt+" -- "+title+" evts").c_str());
        LnLikAccDist->SetName(("LnLikAcc_Evt"+sEvt).c_str()); LnLikAccDist->SetTitle(("LnLikAcc distribution for event "+sEvt+" -- "+title+" evts").c_str());
      }
      //Lik[config-1] = weight;         LikXS[config-1] = weight/MGXS[config-1];              LikAcc[config-1] = weight/MGXSCut[config-1];
      LnLik[config-1] = -log(weight); LnLikXS[config-1] = -log(weight)+log(MGXS[config-1]); LnLikAcc[config-1] = -log(weight)+log(MGXSCut[config-1]);
      //weightsValue[evt-1][config-1] = weight;

      //---  Fill the LnLik histograms for each event and for all events together  ---//
      LnLikDist->SetBinContent(   LnLikDist->FindBin(Var[config-1]),    LnLik[config-1]);
      LnLikXSDist->SetBinContent( LnLikXSDist->FindBin(Var[config-1]),  LnLikXS[config-1]);
      LnLikAccDist->SetBinContent(LnLikAccDist->FindBin(Var[config-1]), LnLikAcc[config-1]);
      LnLikAll->SetBinContent(   LnLikAll->FindBin(Var[config-1]),    LnLikAll->GetBinContent(LnLikAll->FindBin(Var[config-1]))       + LnLik[config-1]   );
      LnLikXSAll->SetBinContent( LnLikXSAll->FindBin(Var[config-1]),  LnLikXSAll->GetBinContent(LnLikXSAll->FindBin(Var[config-1]))   + LnLikXS[config-1] );
      LnLikAccAll->SetBinContent(LnLikAccAll->FindBin(Var[config-1]), LnLikAccAll->GetBinContent(LnLikAccAll->FindBin(Var[config-1])) + LnLikAcc[config-1]);

      //---  Only perform the fit after all configurations are considered!  ---//
      if( config == NrConfigs){
        double LnLikFunction[2][NrConfigs], LnLikXSFunction[2][NrConfigs], LnLikAccFunction[2][NrConfigs];
        vector<double> FctDevOuter, RelFctDevOuter,FctDevXSOuter, RelFctDevXSOuter , FctDevAccOuter, RelFctDevAccOuter;
        double TotalFctDevOuter = 0, TotalRelFctDevOuter = 0, TotalFctDevXSOuter = 0, TotalRelFctDevXSOuter = 0, TotalFctDevAccOuter = 0, TotalRelFctDevAccOuter = 0;

        for(int ii=0; ii < 2; ii++){
          cHat[ii] = LnLik[xNeg[ii]]*Var[xPos[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]])) - LnLik[xMin]*Var[xNeg[ii]]*Var[xPos[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLik[xPos[ii]]*Var[xNeg[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]));
          cHatXS[ii] = LnLikXS[xNeg[ii]]*Var[xPos[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]])) - LnLikXS[xMin]*Var[xNeg[ii]]*Var[xPos[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLikXS[xPos[ii]]*Var[xNeg[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]));
          cHatAcc[ii] = LnLikAcc[xNeg[ii]]*Var[xPos[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]])) - LnLikAcc[xMin]*Var[xNeg[ii]]*Var[xPos[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLikAcc[xPos[ii]]*Var[xNeg[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]));

          bHat[ii] = LnLik[xMin]*(Var[xPos[ii]]+Var[xNeg[ii]])/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLik[xPos[ii]]*(Var[xNeg[ii]]+Var[xMin])/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]))-LnLik[xNeg[ii]]*(Var[xMin]+Var[xPos[ii]])/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]]));
          bHatXS[ii] = LnLikXS[xMin]*(Var[xPos[ii]]+Var[xNeg[ii]])/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLikXS[xPos[ii]]*(Var[xNeg[ii]]+Var[xMin])/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]))-LnLikXS[xNeg[ii]]*(Var[xMin]+Var[xPos[ii]])/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]]));
          bHatAcc[ii] = LnLikAcc[xMin]*(Var[xPos[ii]]+Var[xNeg[ii]])/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLikAcc[xPos[ii]]*(Var[xNeg[ii]]+Var[xMin])/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]))-LnLikAcc[xNeg[ii]]*(Var[xMin]+Var[xPos[ii]])/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]]));

          aHat[ii] = LnLik[xPos[ii]]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLik[xMin]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLik[xNeg[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xNeg[ii]]));
          aHatXS[ii] = LnLikXS[xPos[ii]]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLikXS[xMin]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLikXS[xNeg[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xNeg[ii]]));
          aHatAcc[ii] = LnLikAcc[xPos[ii]]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLikAcc[xMin]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLikAcc[xNeg[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xNeg[ii]]));

          for( int iConfig = 0; iConfig < NrConfigs; iConfig++){
            LnLikFunction[ii][iConfig] = aHat[ii]*Var[iConfig]*Var[iConfig]+bHat[ii]*Var[iConfig]+cHat[ii];
            LnLikXSFunction[ii][iConfig] = aHatXS[ii]*Var[iConfig]*Var[iConfig]+bHatXS[ii]*Var[iConfig]+cHatXS[ii];
            LnLikAccFunction[ii][iConfig] = aHatAcc[ii]*Var[iConfig]*Var[iConfig]+bHatAcc[ii]*Var[iConfig]+cHatAcc[ii];

            if(ii == 1){
              FctDevOuter.push_back(abs(LnLikFunction[1][iConfig]-LnLik[iConfig]));       RelFctDevOuter.push_back(FctDevOuter[iConfig]/LnLik[iConfig]);
              FctDevXSOuter.push_back(abs(LnLikXSFunction[1][iConfig]-LnLikXS[iConfig]));    RelFctDevXSOuter.push_back(FctDevXSOuter[iConfig]/LnLikXS[iConfig]);
              FctDevAccOuter.push_back(abs(LnLikAccFunction[1][iConfig]-LnLikAcc[iConfig])); RelFctDevAccOuter.push_back(FctDevAccOuter[iConfig]/LnLikAcc[iConfig]);
              TotalFctDevOuter    = TotalFctDevOuter+FctDevOuter[iConfig];       TotalRelFctDevOuter    = TotalRelFctDevOuter+RelFctDevOuter[iConfig];
              TotalFctDevXSOuter  = TotalFctDevXSOuter+FctDevXSOuter[iConfig];   TotalRelFctDevXSOuter  = TotalRelFctDevXSOuter+RelFctDevXSOuter[iConfig];
              TotalFctDevAccOuter = TotalFctDevAccOuter+FctDevAccOuter[iConfig]; TotalRelFctDevAccOuter = TotalRelFctDevAccOuter+RelFctDevAccOuter[iConfig];
            }
          }
        }//Looping over inner and outer parabolic function!

        //---  Save the LnLik distributions (event-per-event) together with the function "fit" ---//
        //-- Transformation needed in order to work with TGraphs --//
        double xVar[NrConfigs], yLnLik[NrConfigs], yLnLikXS[NrConfigs], yLnLikAcc[NrConfigs];
        for( int i = 0; i < NrConfigs; i++){ xVar[i] = Var[i]; yLnLik[i] = LnLikFunction[1][i]; yLnLikXS[i] = LnLikXSFunction[1][i]; yLnLikAcc[i] = LnLikAccFunction[1][i]; }

        //std::cout << " TotalFctDevOuter value is : " << TotalFctDevOuter << std::endl;
        stringstream ssTotalFctDevOuter;       ssTotalFctDevOuter << TotalFctDevOuter;             string sTotalFctDevOuter       = ssTotalFctDevOuter.str();
        stringstream ssTotalRelFctDevOuter;    ssTotalRelFctDevOuter << TotalRelFctDevOuter;       string sTotalRelFctDevOuter    = ssTotalRelFctDevOuter.str();
        stringstream ssTotalFctDevXSOuter;     ssTotalFctDevXSOuter << TotalFctDevXSOuter;         string sTotalFctDevXSOuter     = ssTotalFctDevXSOuter.str();
        stringstream ssTotalRelFctDevXSOuter;  ssTotalRelFctDevXSOuter << TotalRelFctDevXSOuter;   string sTotalRelFctDevXSOuter  = ssTotalRelFctDevXSOuter.str();
        stringstream ssTotalFctDevAccOuter;    ssTotalFctDevAccOuter << TotalFctDevAccOuter;       string sTotalFctDevAccOuter    = ssTotalFctDevAccOuter.str();
        stringstream ssTotalRelFctDevAccOuter; ssTotalRelFctDevAccOuter << TotalRelFctDevAccOuter; string sTotalRelFctDevAccOuter = ssTotalRelFctDevAccOuter.str();

        //TGraph* LnLikFctOuter[nEvts];
        LnLikFctOuter[evt]    = new TGraph(NrConfigs,xVar,yLnLik);    LnLikFctOuter[evt]->SetMarkerColor(2);    LnLikFctOuter[evt]->SetLineColor(2);
        TGraph* LnLikXSFctOuter  = new TGraph(NrConfigs,xVar,yLnLikXS);  LnLikXSFctOuter->SetMarkerColor(2);  LnLikXSFctOuter->SetLineColor(2);
        TGraph* LnLikAccFctOuter = new TGraph(NrConfigs,xVar,yLnLikAcc); LnLikAccFctOuter->SetMarkerColor(2); LnLikAccFctOuter->SetLineColor(2);
        LnLikFctOuter[evt]->SetName(("LnLikFctOuter_Evt"+sEvt).c_str());

        LnLikFctOuter[evt]->SetTitle(("LnLik for event "+sEvt+" -- Outer points used (Fct deviation is "+sTotalFctDevOuter+" -- "+sTotalRelFctDevOuter+")").c_str());
        LnLikXSFctOuter->SetTitle(("LnLikXS for event "+sEvt+" -- Outer points used (Fct deviation is "+sTotalFctDevXSOuter+" -- "+sTotalRelFctDevXSOuter+")").c_str());
        LnLikAccFctOuter->SetTitle(("LnLikAcc for event "+sEvt+" -- Outer points used (Fct deviation is "+sTotalFctDevAccOuter+" -- "+sTotalRelFctDevAccOuter+")").c_str());
        TCanvas* LnLikCanv =    new TCanvas(("LnLikCanv_"+sEvt).c_str(),"LnLik");      LnLikCanv->cd();   LnLikFctOuter[evt]->Draw("AC*");   LnLikDist->Draw("samep");   LnLikDir->cd();   LnLikCanv->Write();
        TCanvas* LnLikXSCanv =  new TCanvas(("LnLikXSCanv_"+sEvt).c_str(),"LnLikXS");  LnLikXSCanv->cd(); LnLikXSFctOuter->Draw("AC*"); LnLikXSDist->Draw("samep"); LnLikXSDir->cd(); LnLikXSCanv->Write();
        TCanvas* LnLikAccCanv = new TCanvas(("LnLikAccCanv_"+sEvt).c_str(),"LnLikAcc");LnLikAccCanv->cd();LnLikAccFctOuter->Draw("AC*");LnLikAccDist->Draw("samep");LnLikAccDir->cd();LnLikAccCanv->Write();

        //Save 20 of these histograms in one TCanvas!
        stringstream ssNrCanvas; ssNrCanvas << NrCanvas; string sNrCanvas = ssNrCanvas.str();
        std::cout << " Comparing " << consEvts << " with " << xDivide*yDivide*(NrCanvas+1) << std::endl;
        if( consEvts == (xDivide*yDivide*(NrCanvas+1))){
          std::cout << " ------------------ NrCanvas value = " << NrCanvas << std::endl;
          LnLikStackDir->cd(); 
          StackCanvas->Write(); 
          StackCanvas->SetName( ("StackCanvas_Nr"+sNrCanvas).c_str() );
          StackCanvas->Divide(xDivide,yDivide);
          NrCanvas++;
        }
        StackCanvas->cd(consEvts - (xDivide*yDivide*NrCanvas) );
        std::cout << " -- Storing on canvas nr : " << consEvts - (xDivide*yDivide*(NrCanvas)) << std::endl;
        LnLikFctOuter[evt]->Draw("AC*");   LnLikDist->Draw("samep");

        //---  Calculate the first derivative distribution and save them in event-by-event plot  ---//
        FstDer->SetName(("FirstDerivative_Evt"+sEvt).c_str());
        FstDer->SetTitle(("First derivative of -ln(likelihood) distribution for event"+sEvt+" (no normalisation) -- "+title+" evts").c_str());
        for(int ii = 0; ii < FstDer->GetNbinsX(); ii++)
          FstDer->SetBinContent(ii+1, (LnLik[ii+2] - LnLik[ii+3])/xStep[ii]);
        FstDerDir->cd(); FstDer->Write();

        //---  Calculate the positive and negative deviation  ---//
        double yPlusPlus[3] = {LnLik[xPos[1]] - LnLikFunction[0][xPos[1]], LnLikXS[xPos[1]] - LnLikXSFunction[0][xPos[1]], LnLikAcc[xPos[1]] - LnLikAccFunction[0][xPos[1]]};
        double yMinMin[3]   = {LnLik[xNeg[1]] - LnLikFunction[0][xNeg[1]], LnLikXS[xNeg[1]] - LnLikXSFunction[0][xNeg[1]], LnLikAcc[xNeg[1]] - LnLikAccFunction[0][xNeg[1]]};
        double yPlus[3] = {LnLik[xPos[0]] - LnLikFunction[1][xPos[0]], LnLikXS[xPos[0]] - LnLikXSFunction[1][xPos[0]], LnLikAcc[xPos[0]] - LnLikAccFunction[1][xPos[0]]};
        double yMin[3]  = {LnLik[xNeg[0]] - LnLikFunction[1][xNeg[0]], LnLikXS[xNeg[0]] - LnLikXSFunction[1][xNeg[0]], LnLikAcc[xNeg[0]] - LnLikAccFunction[1][xNeg[0]]};
        double scdDerInner[3] = {LnLik[xNeg[0]]-2*LnLik[xMin]+LnLik[xPos[0]], LnLikXS[xNeg[0]]-2*LnLikXS[xMin]+LnLikXS[xPos[0]], LnLikAcc[xNeg[0]]-2*LnLikAcc[xMin]+LnLikAcc[xPos[0]]};
        double scdDerOuter[3]   = {LnLik[xNeg[1]]-2*LnLik[xMin]+LnLik[xPos[1]], LnLikXS[xNeg[1]]-2*LnLikXS[xMin]+LnLikXS[xPos[1]], LnLikAcc[xNeg[1]]-2*LnLikAcc[xMin]+LnLikAcc[xPos[1]]};

        //--- Fill the histograms  ---//
        YPlus->Fill(yPlus[0]);                                     YPlusXS->Fill(yPlus[1]);                                     YPlusAcc->Fill(yPlus[2]);
        YPlusPlus->Fill(yPlusPlus[0]);                             YPlusPlusXS->Fill(yPlusPlus[1]);                             YPlusPlusAcc->Fill(yPlusPlus[2]);
        YPlusGausTest->Fill(yPlus[0] + yPlusPlus[0]/4);            YPlusGausTestXS->Fill(yPlus[1] + yPlusPlus[1]/4);            YPlusGausTestAcc->Fill(yPlus[2] + yPlusPlus[2]/4);
        //YRelPlus->Fill(yPlus[0]/LnLik[6]);                         YRelPlusXS->Fill(yPlus[1]/LnLikXS[6]);                       YRelPlusAcc->Fill(yPlus[2]/LnLikAcc[6]);
        YMin->Fill(yMin[0]);                                       YMinXS->Fill(yMin[1]);                                       YMinAcc->Fill(yMin[2]);
        YMinMin->Fill(yMinMin[0]);                                 YMinMinXS->Fill(yMinMin[1]);                                 YMinMinAcc->Fill(yMinMin[2]);
        //YRelMin->Fill(yMin[0]/LnLik[2]);                           YRelMinXS->Fill(yMin[1]/LnLikXS[2]);                         YRelMinAcc->Fill(yMin[2]/LnLikAcc[2]);
        ScdDerInner->Fill(scdDerInner[0]);                         ScdDerXSInner->Fill(scdDerInner[1]);                         ScdDerAccInner->Fill(scdDerInner[2]);
        ScdDerOuter->Fill(scdDerOuter[0]);                         ScdDerXSOuter->Fill(scdDerOuter[1]);                         ScdDerAccOuter->Fill(scdDerOuter[2]);
        ScdDerScatter->Fill(scdDerOuter[0], scdDerInner[0]);       ScdDerXSScatter->Fill(scdDerOuter[1], scdDerInner[1]);       ScdDerAccScatter->Fill(scdDerOuter[2], scdDerInner[2]);
        if( TotalFctDevAccOuter > 5) std::cout << "Overflow found for TotalFctDevDist : " << TotalFctDevAccOuter << std::endl;

        //-- Apply cut on YPlusGausTest --//
        if( (yPlus[0] + yPlusPlus[0]/4) <= 0.025 && (yPlus[0] + yPlusPlus[0]/4) >= -0.025) EvtsWithYPlusGausSmall.push_back(evt);
        if( (yPlus[1] + yPlusPlus[1]/4) <= 0.025 && (yPlus[1] + yPlusPlus[1]/4) >= -0.025) EvtsWithYPlusGausSmallXS.push_back(evt);
        if( (yPlus[2] + yPlusPlus[2]/4) <= 0.025 && (yPlus[2] + yPlusPlus[2]/4) >= -0.025) EvtsWithYPlusGausSmallAcc.push_back(evt);

        //------------------------------------------------------------//
        //---  Now apply cuts and save the interesting histograms  ---//
        //------------------------------------------------------------//
        for(int iConf = 0; iConf < NrConfigs; iConf++){

          //-- Apply cut on TotalFctDeviation --//
          if( TotalFctDevAccOuter <= 0.1){
            if(iConf == 0) EvtsWithSmallFctDev.push_back(evt);  //Still an interesting variable??
            LLSmallFctDevDist->SetBinContent(LLSmallFctDevDist->FindBin(Var[iConf]), LLSmallFctDevDist->GetBinContent(LLSmallFctDevDist->FindBin(Var[iConf])) + LnLikAcc[iConf]);
          }
          //-- Apply cut on ScdDer (using inner Var points) --//
          if( scdDerInner[0] > 0.0){
            if(iConf == 0) EvtsWithPosScdDerInner++; 
            LLPosScdDerDistInner->SetBinContent(LLPosScdDerDistInner->FindBin(Var[iConf]), LLPosScdDerDistInner->GetBinContent(LLPosScdDerDistInner->FindBin(Var[iConf])) + LnLik[iConf] );          
          }
          if( scdDerInner[1] > 0.0){
            if(iConf == 0) EvtsWithPosScdDerXSInner++;
            LLXSPosScdDerDistInner->SetBinContent(LLXSPosScdDerDistInner->FindBin(Var[iConf]), LLXSPosScdDerDistInner->GetBinContent(LLXSPosScdDerDistInner->FindBin(Var[iConf])) + LnLikXS[iConf] );
          }
          if( scdDerInner[2] > 0.0){
            if(iConf == 0) EvtsWithPosScdDerAccInner++;
            LLAccPosScdDerDistInner->SetBinContent(LLAccPosScdDerDistInner->FindBin(Var[iConf]), LLAccPosScdDerDistInner->GetBinContent(LLAccPosScdDerDistInner->FindBin(Var[iConf])) + LnLikAcc[iConf] );
          }

          //-- Apply cut on both ScdDer --//        
          if( scdDerInner[0] > 0.0 && scdDerOuter[0] > 0.0){
            if(iConf == 0) EvtsWithPosScdDerBoth++;
            LLPosScdDerDistBoth->SetBinContent(LLPosScdDerDistBoth->FindBin(Var[iConf]), LLPosScdDerDistBoth->GetBinContent(LLPosScdDerDistBoth->FindBin(Var[iConf])) + LnLik[iConf] );
          }
          if( scdDerInner[1] > 0.0 && scdDerOuter[1] > 0.0){
            if(iConf == 0) EvtsWithPosScdDerXSBoth++;
            LLXSPosScdDerDistBoth->SetBinContent(LLXSPosScdDerDistBoth->FindBin(Var[iConf]), LLXSPosScdDerDistBoth->GetBinContent(LLXSPosScdDerDistBoth->FindBin(Var[iConf])) + LnLikXS[iConf] );
          }
          if( scdDerInner[2] > 0.0 && scdDerOuter[2] > 0.0){
            if(iConf == 0) EvtsWithPosScdDerAccBoth++;
            LLAccPosScdDerDistBoth->SetBinContent(LLAccPosScdDerDistBoth->FindBin(Var[iConf]), LLAccPosScdDerDistBoth->GetBinContent(LLAccPosScdDerDistBoth->FindBin(Var[iConf])) + LnLikAcc[iConf] );
          }

          //-- Apply cut on ScdDer (using outer Var points) --//
          if( scdDerOuter[0] > 0.0){
            if(iConf == 0){
              EvtsWithPosScdDerOuter++;
              YPlusGausTestPosScdDer->Fill(yPlus[0] + yPlusPlus[0]/4);
              YPlusPosScdDer->Fill(yPlus[0]); YPlusPlusPosScdDer->Fill(yPlusPlus[0]); YMinPosScdDer->Fill(yMin[0]); YMinMinPosScdDer->Fill(yMinMin[0]);
            }
            LLPosScdDerDistOuter->SetBinContent(LLPosScdDerDistOuter->FindBin(Var[iConf]), LLPosScdDerDistOuter->GetBinContent(LLPosScdDerDistOuter->FindBin(Var[iConf])) + LnLik[iConf] );
          }
          else{  //Are these distributions for negative scdDer still interesting?
            if(iConf == 0){
              YPlusGausTestNegScdDer->Fill(yPlus[0] + yPlusPlus[0]/4);
              YPlusNegScdDer->Fill(yPlus[0]); YPlusPlusNegScdDer->Fill(yPlusPlus[0]); YMinNegScdDer->Fill(yMin[0]); YMinMinNegScdDer->Fill(yMinMin[0]);
            }
          }

          if( scdDerOuter[1] > 0.0){
            if(iConf == 0){
              EvtsWithPosScdDerXSOuter++;
              YPlusGausTestXSPosScdDer->Fill(yPlus[1] + yPlusPlus[1]/4);
              YPlusXSPosScdDer->Fill(yPlus[1]); YPlusPlusXSPosScdDer->Fill(yPlusPlus[1]); YMinXSPosScdDer->Fill(yMin[1]); YMinMinXSPosScdDer->Fill(yMinMin[1]);
            }
            LLXSPosScdDerDistOuter->SetBinContent(LLXSPosScdDerDistOuter->FindBin(Var[iConf]), LLXSPosScdDerDistOuter->GetBinContent(LLXSPosScdDerDistOuter->FindBin(Var[iConf])) + LnLikXS[iConf] );
          }
          else{
            if(iConf == 0){
              YPlusGausTestXSNegScdDer->Fill(yPlus[1] + yPlusPlus[1]/4);
              YPlusXSNegScdDer->Fill(yPlus[1]); YPlusPlusXSNegScdDer->Fill(yPlusPlus[1]); YMinXSNegScdDer->Fill(yMin[1]); YMinMinXSNegScdDer->Fill(yMinMin[1]);
            }
          }

          if( scdDerOuter[2] > 0.0){
            if(iConf == 0){
              EvtsWithPosScdDerAccOuter++;
              YPlusGausTestAccPosScdDer->Fill(yPlus[2] + yPlusPlus[2]/4);
              YPlusAccPosScdDer->Fill(yPlus[2]); YPlusPlusAccPosScdDer->Fill(yPlusPlus[2]); YMinAccPosScdDer->Fill(yMin[2]); YMinMinAccPosScdDer->Fill(yMinMin[2]);
            }
            LLAccPosScdDerDistOuter->SetBinContent(LLAccPosScdDerDistOuter->FindBin(Var[iConf]), LLAccPosScdDerDistOuter->GetBinContent(LLAccPosScdDerDistOuter->FindBin(Var[iConf])) + LnLikAcc[iConf] );
          }
          else{
            if(iConf == 0){  
              YPlusGausTestAccNegScdDer->Fill(yPlus[2] + yPlusPlus[2]/4);
              YPlusAccNegScdDer->Fill(yPlus[2]); YPlusPlusAccNegScdDer->Fill(yPlusPlus[2]); YMinAccNegScdDer->Fill(yMin[2]); YMinMinAccNegScdDer->Fill(yMinMin[2]);
            }
          }
        }

      }
    }
  }
  ifs.close();

  //--- Save all the histograms containing information about all the events! ---//
  Tfile->cd();
  LnLikAll->Write();                LnLikXSAll->Write();                LnLikAccAll->Write();
  YPlusGausTest->Write();           YPlusGausTestXS->Write();           YPlusGausTestAcc->Write();
  YPlus->Write();                   YPlusXS->Write();                   YPlusAcc->Write();
  YPlusPlus->Write();               YPlusPlusXS->Write();               YPlusPlusAcc->Write();
  YRelPlus->Write();                YRelPlusXS->Write();                YRelPlusAcc->Write();
  YMin->Write();                    YMinXS->Write();                    YMinAcc->Write();
  YMinMin->Write();                 YMinMinXS->Write();                 YMinMinAcc->Write();
  YRelMin->Write();                 YRelMinXS->Write();                 YRelMinAcc->Write();
  ScdDerInner->Write();             ScdDerXSInner->Write();             ScdDerAccInner->Write();
  ScdDerOuter->Write();             ScdDerXSOuter->Write();             ScdDerAccOuter->Write();
  ScdDerScatter->Write();           ScdDerXSScatter->Write();           ScdDerAccScatter->Write();

  //---  Draw the likelihood distribution separately for events surviving and passing the cuts!  ---//
  std::cout << "Nr of events with 2nd derivative > 0 (LnLik, LnLikXS & LnLikAcc -- using x = " << Var[xNeg[0]] << "/" << Var[xMin] << "/" << Var[xPos[0]] << ") : " << EvtsWithPosScdDerInner << ", " << EvtsWithPosScdDerXSInner << " & " << EvtsWithPosScdDerAccInner << std::endl;
  std::cout << "Nr of events with 2nd derivative > 0 (LnLik, LnLikXS & LnLikAcc -- using x = " << Var[xNeg[1]] << "/" << Var[xMin] << "/" << Var[xPos[1]] << ") : " << EvtsWithPosScdDerOuter << ", " << EvtsWithPosScdDerXSOuter << " & " << EvtsWithPosScdDerAccOuter << std::endl;
  std::cout << "Nr of events with both 2nd derivatives > 0 (LnLik, LnLikXS & LnLikAcc) : " << EvtsWithPosScdDerOuter << ", " << EvtsWithPosScdDerXSOuter << " & " << EvtsWithPosScdDerAccOuter << std::endl;
  std::cout << "Nr of events with small Gaussian-test deviation for plus-case (LnLik, LnLikXS & LnLikAcc) " << EvtsWithYPlusGausSmall.size() << ", " << EvtsWithYPlusGausSmallXS.size() << " & " << EvtsWithYPlusGausSmallAcc.size() << std::endl;
  std::cout << "Nr of events with total function deviation < 0.5 : " << EvtsWithSmallFctDev.size() << std::endl;

  TDirectory* AppliedCutsDir = Tfile->mkdir("LikelihoodAfterCuts");
  TDirectory* SmallFctDevDir = AppliedCutsDir->mkdir("SmallFunctionDeviation");
  SmallFctDevDir->cd(); LLSmallFctDevDist->Write();

  stringstream ssEvtsWithPosScdDerOuter;    ssEvtsWithPosScdDerOuter << EvtsWithPosScdDerOuter;    string sEvtsWithPosScdDerOuter    = ssEvtsWithPosScdDerOuter.str();
  stringstream ssEvtsWithPosScdDerXSOuter;  ssEvtsWithPosScdDerXSOuter << EvtsWithPosScdDerOuter;  string sEvtsWithPosScdDerXSOuter  = ssEvtsWithPosScdDerXSOuter.str();
  stringstream ssEvtsWithPosScdDerAccOuter; ssEvtsWithPosScdDerAccOuter << EvtsWithPosScdDerOuter; string sEvtsWithPosScdDerAccOuter = ssEvtsWithPosScdDerAccOuter.str();
  stringstream ssEvtsWithPosScdDerInner;    ssEvtsWithPosScdDerInner << EvtsWithPosScdDerInner;    string sEvtsWithPosScdDerInner    = ssEvtsWithPosScdDerInner.str();
  stringstream ssEvtsWithPosScdDerXSInner;  ssEvtsWithPosScdDerXSInner << EvtsWithPosScdDerInner;  string sEvtsWithPosScdDerXSInner  = ssEvtsWithPosScdDerXSInner.str();
  stringstream ssEvtsWithPosScdDerAccInner; ssEvtsWithPosScdDerAccInner << EvtsWithPosScdDerInner; string sEvtsWithPosScdDerAccInner = ssEvtsWithPosScdDerAccInner.str();
  stringstream ssEvtsWithPosScdDerBoth;     ssEvtsWithPosScdDerBoth << EvtsWithPosScdDerBoth;      string sEvtsWithPosScdDerBoth     = ssEvtsWithPosScdDerBoth.str();
  stringstream ssEvtsWithPosScdDerXSBoth;   ssEvtsWithPosScdDerXSBoth << EvtsWithPosScdDerBoth;    string sEvtsWithPosScdDerXSBoth   = ssEvtsWithPosScdDerXSBoth.str();
  stringstream ssEvtsWithPosScdDerAccBoth;  ssEvtsWithPosScdDerAccBoth << EvtsWithPosScdDerBoth;   string sEvtsWithPosScdDerAccBoth  = ssEvtsWithPosScdDerAccBoth.str();
  stringstream ssNEvts; ssNEvts << nEvts; string sNEvts = ssNEvts.str();

  LLPosScdDerDistBoth->SetTitle( ("-ln(L) when both 2nd derivatives > 0 (no norm -- "+sEvtsWithPosScdDerBoth+"/"+sNEvts+" evts -- "+title+")").c_str());
  LLXSPosScdDerDistBoth->SetTitle( ("-ln(L) when both 2nd derivatives > 0 (XS norm -- "+sEvtsWithPosScdDerXSBoth+"/"+sNEvts+" evts -- "+title+")").c_str());
  LLAccPosScdDerDistBoth->SetTitle( ("-ln(L) when both 2nd derivatives > 0 (Acc norm -- "+sEvtsWithPosScdDerAccBoth+"/"+sNEvts+" evts -- "+title+")").c_str());
  LLPosScdDerDistInner->SetTitle( ("-ln(L) when inner 2nd derivative > 0 (no norm -- "+sEvtsWithPosScdDerInner+"/"+sNEvts+" evts -- "+title+")").c_str()    );
  LLXSPosScdDerDistInner->SetTitle( ("-ln(L) when inner 2nd derivative > 0 (XS norm -- "+sEvtsWithPosScdDerXSInner+"/"+sNEvts+" evts -- "+title+")").c_str() );
  LLAccPosScdDerDistInner->SetTitle( ("-ln(L) when inner 2nd derivative > 0 (Acc norm -- "+sEvtsWithPosScdDerAccInner+"/"+sNEvts+" evts -- "+title+")").c_str() );
  LLPosScdDerDistOuter->SetTitle( ("-ln(L) when outer 2nd derivative > 0 (no norm -- "+sEvtsWithPosScdDerOuter+"/"+sNEvts+" evts -- "+title+")").c_str() );
  LLXSPosScdDerDistOuter->SetTitle( ("-ln(L) when outer 2nd derivative > 0 (XS norm -- "+sEvtsWithPosScdDerXSOuter+"/"+sNEvts+" evts -- "+title+")").c_str() );
  LLAccPosScdDerDistOuter->SetTitle( ("-ln(L) when outer 2nd derivative > 0 (Acc norm -- "+sEvtsWithPosScdDerAccOuter+"/"+sNEvts+" evts -- "+title+")").c_str() );

  TDirectory* SignScdDerDir = AppliedCutsDir->mkdir("SignSecondDerivative"); SignScdDerDir->cd();
  LLPosScdDerDistInner->Write(); LLXSPosScdDerDistInner->Write(); LLAccPosScdDerDistInner->Write();
  LLPosScdDerDistOuter->Write(); LLXSPosScdDerDistOuter->Write(); LLAccPosScdDerDistOuter->Write();
  LLPosScdDerDistBoth->Write();  LLXSPosScdDerDistBoth->Write();  LLAccPosScdDerDistBoth->Write();

  /*for(int iEvt = 0; iEvt < nEvts; iEvt++){
    if(std::find(EvtsWithSmallFctDev.begin(), EvtsWithSmallFctDev.end(), iEvt) != EvtsWithSmallFctDev.end() ){
      std::cout << iEvt << ") Event number is found : " << EvtsWithSmallFctDev[iEvt]  << std::endl;
    }
  }*/


}


