#include <iostream>
#include <fstream>
#include <sstream>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include <string>
#include <map>
#include <stdlib.h>

std::string VarValues[] = {"Re(V_{R}) = -0.1","Re(V_{R}) = -0.09","Re(VR) = -0.08","Re(V_{R}) = -0.07","Re(VR) = -0.06","Re(V_{R}) = -0.05","Re(V_{R}) = -0.04","Re(V_{R}) = -0.03","Re(V_{R}) = -0.02","Re(V_{R}) = -0.01","Re(VR) = 0.0","Re(V_{R}) = 0.01","Re(VR) = 0.02","Re(V_{R}) = 0.03","Re(V_{R}) = 0.04","Re(V_{R}) = 0.05","Re(V_{R}) = 0.06","Re(V_{R}) = 0.07","Re(V_{R}) = 0.08","Re(V_{R}) = 0.09","Re(V_{R}) = 0.1"};
double Var[] = {-0.1,     -0.09,    -0.08,    -0.07,    -0.06,    -0.05,    -0.04,    -0.03,   -0.02,    -0.01,    0.0,      0.01,     0.02,     0.03,     0.04,     0.05,     0.06,     0.07,     0.08,     0.09,    0.1};
double MGXS[] = {11.25909, 11.19599, 11.14075, 11.10667, 11.06145, 11.02784, 10.99474, 10.9534, 10.93846, 10.91787, 10.90059, 10.89643, 10.86829, 10.87792, 10.87266, 10.88228, 10.87684, 10.89534, 10.91641, 10.9317, 10.97767};
double MGXSCut[] = {2.43546,  2.42182,  2.41059,  2.40057,  2.38839,  2.38187,  2.36976,  2.36513, 2.35512,  2.35666,  2.35415,  2.35694,  2.35174,  2.34909,  2.34392,  2.35108,  2.34767,  2.35477,  2.36148,  2.3643,  2.37424};
int xBin = 21;
float xLow = -0.105;
float xHigh = 0.105;
int xMinValue[] = {4, 4, 10};
std::string KinVar = "Re(V_{R})";
int VarWindow = 3;
float xPos[] = {xMinValue[VarWindow-1]+1,xMinValue[VarWindow-1]+2};
float xNeg[] = {xMinValue[VarWindow-1]-1,xMinValue[VarWindow-1]-2};
std::string title = "title";
TFile* Tfile = new TFile("Test.root","RECREATE");

int NrConfigs = 21;
int nEvts = 10; 

void ReadTest(){
  
  float xMin = xMinValue[VarWindow-1];
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
  (FstDer->GetXaxis())->SetBinLabel(1,"test"); //("(y_{DATA}(x="+Var[xNeg[1]]+") - y_{DATA}(x="+Var[xNeg[0]]+"))/"+xStep[0]).c_str());
  FstDer->GetXaxis()->SetBinLabel(2,"bin2 test"); //("(y_{DATA}(x="+Var[xNeg[0]]+") - y_{DATA}(x="+Var[xMin]+"))/"+xStep[1]).c_str());
  FstDer->GetXaxis()->SetBinLabel(3,"bin3 test"); //("(y_{DATA}(x="+Var[xMin]+") - y_{DATA}(x="+Var[xPos[0]]+"))/"+xStep[2]).c_str());
  FstDer->GetXaxis()->SetBinLabel(4,"bin4 test"); //("(y_{DATA}(x="+Var[xPos[0]]+") - y_{DATA}(x="+Var[xPos[1]]+"))/"+xStep[3]).c_str());
  TH1F* FstDerXS  = new TH1F("FirstDerivativeXS", ("First derivative of -ln(likelihood) distribution (XS norm -- "+title+" evts)").c_str(), 5,-0.25,0.25);
  TH1F* FstDerAcc = new TH1F("FirstDerivativeAcc",("First derivative of -ln(likelihood) distribution (Acc norm -- "+title+" evts)").c_str(),5,-0.25,0.25);
  TH1F* ScdDerInner    = new TH1F("SecondDerivativeInner",   ("Second derivative of -ln(likelihood) distribution (no norm -- using inner points -- "+title+" evts)").c_str(), 250,-5,5);
  TH1F* ScdDerXSInner  = new TH1F("SecondDerivativeXSInner", ("Second derivative of -ln(likelihood) distribution (XS norm -- using inner points -- "+title+" evts)").c_str(), 250,-5,5);
  TH1F* ScdDerAccInner = new TH1F("SecondDerivativeAccInner",("Second derivative of -ln(likelihood) distribution (Acc norm -- using inner points -- "+title+" evts)").c_str(),250,-5,5);
  TH1F* ScdDerOuter    = new TH1F("SecondDerivativeOuter",   ("Second derivative of -ln(likelihood) distribution (no norm -- using outer points -- "+title+" evts)").c_str(), 250,-5,5);
  TH1F* ScdDerXSOuter  = new TH1F("SecondDerivativeXSOuter", ("Second derivative of -ln(likelihood) distribution (XS norm -- using outer points -- "+title+" evts)").c_str(), 250,-5,5);
  TH1F* ScdDerAccOuter = new TH1F("SecondDerivativeAccOuter",("Second derivative of -ln(likelihood) distribution (Acc norm -- using outer points -- "+title+" evts)").c_str(),250,-5,5);
  TH2F* ScdDerScatter    = new TH2F("ScdDerScatterPlot",   ("Second derivative of -ln(L) using inner points versus using outer points (no norm -- "+title+" evts)").c_str(), 250,-5,5,250,-5,5);
  TH2F* ScdDerXSScatter  = new TH2F("ScdDerXSScatterPlot", ("Second derivative of -ln(L) using inner points versus using outer points (XS norm -- "+title+" evts)").c_str(), 250,-5,5,250,-5,5);
  TH2F* ScdDerAccScatter = new TH2F("ScdDerAccScatterPlot",("Second derivative of -ln(L) using inner points versus using outer points (Acc norm -- "+title+" evts)").c_str(),250,-5,5,250,-5,5);
  
  TH1F* FstDerInnerPlusRelToUnc = new TH1F("FirstDer:Inner_Plus_RelativeToUnc",("First derivative of -ln(L) distribution wrt unc of "+KinVar+" = "+Var[xMin]+" (no norm -- pos inner point)").c_str(),250,0,25);
  TH1F* FstDerXSInnerPlusRelToUnc = new TH1F("FirstDerXSInner_Plus_RelativeToUnc",("First derivative of -ln(L) distribution wrt unc of "+KinVar+" = "+Var[xMin]+" (XS norm -- pos inner point)").c_str(),250,0,25);
  TH1F* FstDerAccInnerPlusRelToUnc = new TH1F("FirstDerAccInner_Plus_RelativeToUnc",("First derivative of -ln(L) distribution wrt unc of "+KinVar+" = "+Var[xMin]+" (Acc norm -- pos inner point)").c_str(),250,0,25);
  TH1F* FstDerInnerMinRelToUnc = new TH1F("FirstDerInner_Min_RelativeToUnc",("First derivative of -ln(L) distribution wrt unc of "+KinVar+" = "+Var[xMin]+" (no norm -- neg inner point)").c_str(),250,0,25);
  TH1F* FstDerXSInnerMinRelToUnc = new TH1F("FirstDerXSInner_Min_RelativeToUnc",("First derivative of -ln(L) distribution wrt unc of "+KinVar+" = "+Var[xMin]+" (XS norm -- neg inner point)").c_str(),250,0,25);
  TH1F* FstDerAccInnerMinRelToUnc = new TH1F("FirstDerAccInner_Min_RelativeToUnc",("First derivative of -ln(L) distribution wrt unc of "+KinVar+" = "+Var[xMin]+" (Acc norm -- neg inner point)").c_str(),250,0,25);
  TH1F* FstDerOuterPlusRelToUnc = new TH1F("FirstDerOuter_Plus_RelativeToUnc",("First derivative of -ln(L) distribution wrt unc of "+KinVar+" = "+Var[xMin]+" (no norm -- pos outer point)").c_str(),250,0,25);
  TH1F* FstDerXSOuterPlusRelToUnc = new TH1F("FirstDerXSOuter_Plus_RelativeToUnc",("First derivative of -ln(L) distribution wrt unc of "+KinVar+" = "+Var[xMin]+" (XS norm -- pos outer point)").c_str(),250,0,25);
  TH1F* FstDerAccOuterPlusRelToUnc = new TH1F("FirstDerAccOuter_Plus_RelativeToUnc",("First derivative of -ln(L) distribution wrt unc of "+KinVar+" = "+Var[xMin]+" (Acc norm -- pos outer point)").c_str(),250,0,25);
  TH1F* FstDerOuterMinRelToUnc = new TH1F("FirstDerOuter_Min_RelativeToUnc",("First derivative of -ln(L) distribution wrt unc of "+KinVar+" = "+Var[xMin]+" (no norm -- neg outer point)").c_str(),250,0,25);
  TH1F* FstDerXSOuterMinRelToUnc = new TH1F("FirstDerXSOuter_Min_RelativeToUnc",("First derivative of -ln(L) distribution wrt unc of "+KinVar+" = "+Var[xMin]+" (XS norm -- neg outer point)").c_str(),250,0,25);
  TH1F* FstDerAccOuterMinRelToUnc = new TH1F("FirstDerAccouter_Min_RelativeToUnc",("First derivative of -ln(L) distribution wrt unc of "+KinVar+" = "+Var[xMin]+" (Acc norm -- neg outer point)").c_str(),250,0,25);
  
  TH2F* WeightVsUnc = new TH2F("WeightVsUnc",("Distribution of weight vs uncertainty (ln applied) for "+KinVar+" = "+Var[xMin]).c_str(),250,20,100,250,0,0.2);
  TH2F* WeightVsUncPlus = new TH2F("WeightVsUnc_Plus",("Distribution of weight vs uncertainty (ln applied) for "+KinVar+" = "+Var[xPos[0]]).c_str(),250,20,100,250,0,0.2);
  TH2F* WeightVsUncMin = new TH2F("WeightVsUnc_Min",("Distribution of weight vs uncertainty (ln applied) for "+KinVar+" = "+Var[xNeg[0]]).c_str(),250,20,100,250,0,0.2);
  
  TH1F* LnLikUncDist = new TH1F("LnLikUncDist_RVR0",("Distribution of uncertainty on -ln(L) of "+KinVar+" = "+Var[xMin]+" (all norm -- approx -- "+title+" events)").c_str(),250,0,0.2);
  TH1F* LnLikUncDistPlus = new TH1F("LnLikUncDist_RVRPlus",("Distribution of uncertainty on -ln(L) of "+KinVar+" = "+Var[xPos[0]]+" (all norm -- approx -- "+title+" events)").c_str(),250,0,0.2);
  TH1F* LnLikUncDistMin = new TH1F("LnLikUncDist_RVRMin",("Distribution of uncertainty on -ln(L) of "+KinVar+" = "+Var[xNeg[0]]+" (all norm -- approx -- "+title+" events)").c_str(),250,0,0.2);
  TH1F* LikUncDist    = new TH1F("LikUncDist_RVR0",("Distribution of uncertainty on Likelihood of "+KinVar+" = "+Var[xMin]+" (no norm -- approx -- "+title+" events)").c_str(),250,0,0.00000000000000000002);
  TH1F* LikXSUncDist  = new TH1F("LikXSUncDist_RVR0",("Distribution of uncertainty on Likelihood of "+KinVar+" = "+Var[xMin]+" (XS norm -- approx -- "+title+" events)").c_str(),250,0,0.00000000000000000002);
  TH1F* LikAccUncDist = new TH1F("LikAccUncDist_RVR0",("Distribution of uncertainty on Likelihood of "+KinVar+" = "+Var[xMin]+" (Acc norm -- approx -- "+title+" events)").c_str(),250,0,0.00000000000000000002);
  TH1F* RelLnLikUncDist    = new TH1F("RelLnLikUncDist_RVR0",("Distribution of #sigma(-ln(L))/-ln(L) for "+KinVar+" = "+Var[xMin]+" (no norm -- "+title+" events)").c_str(),250,0,0.002);
  TH1F* RelLnLikXSUncDist  = new TH1F("RelLnLikXSUncDist_RVR0",("Distribution of #sigma(-ln(L))/-ln(L) for "+KinVar+" = "+Var[xMin]+" (XS norm -- "+title+" events)").c_str(),250,0,0.002);
  TH1F* RelLnLikAccUncDist = new TH1F("RelLnLikAccUncDist_RVR0",("Distribution of #sigma(-ln(L))/-ln(L) for "+KinVar+" = "+Var[xMin]+" (Acc norm -- "+title+" events)").c_str(),250,0,0.002);
  TH1F* RelLikUncDist    = new TH1F("RelLikUncDist_RVR0",("Distribution of #sigma(L)/L for "+KinVar+" = "+Var[xMin]+" (no norm -- "+title+" events)").c_str(),250,0,0.2);
  TH1F* RelLikXSUncDist  = new TH1F("RelLikXSUncDist_RVR0",("Distribution of #sigma(L)/L for "+KinVar+" = "+Var[xMin]+" (XS norm -- "+title+" events)").c_str(),250,0,0.2);
  TH1F* RelLikAccUncDist = new TH1F("RelLikAccUncDist_RVR0",("Distribution of #sigma(L)/L for "+KinVar+" = "+Var[xMin]+" (Acc norm -- "+title+" events)").c_str(),250,0,0.2);
  //TH1F* LnLikXSUncDist = new TH1F("LnLikUncDist_RVR0",("Distribution of uncertainty on -ln(L) of "+KinVar+" = "+Var[xMin]+" (no norm -- "+title+" events)").c_str(),250,0,50);
  //TH1F* LnLikAccUncDist = new TH1F("LnLikUncDist_RVR0",("Distribution of uncertainty on -ln(L) of "+KinVar+" = "+Var[xMin]+" (no norm -- "+title+" events)").c_str(),250,0,50);
  
  TH1F* SigmaVarianceDist = new TH1F("SigmaVariance","Spread of the uncertainty on the -ln(L) distribution",250,0,0.0002);
  TH1F* AverageSigmaDist = new TH1F("AverageSigma","Average value of the uncertainties for each of the considered configurations",250,0,0.2);
  TH1F* LnLikVariationDist = new TH1F("LnLikVariation","Difference between the maximum and the minimum value of the -ln(L) distribution",250,0,1);
  
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
  TDirectory* LnLikXSDir = Tfile->mkdir("LnLikXSDist");
  TDirectory* LnLikAccDir = Tfile->mkdir("LnLikAccDist");
  TDirectory* LnLikAccDirVarVsUnc = Tfile->mkdir("LnLikAccDist_VarLargerThanAvgUnc");
  TDirectory* LnLikAccDirVarVsDUnc = Tfile->mkdir("LnLikAccDist_VarLargerThanTwiceAvgUnc");
  TDirectory* FstDerDir = Tfile->mkdir("FirstDerivativeDist");
  TCanvas* LnLikFitCanvas = new TCanvas("name","title");
  std::string LnLikFitCanvasName = "LnLikFitCanvas_Evt", LnLikFitCanvasTitle = "Comparing fit function from ROOT fit and algebraic function for event -- '+title+' evts ";

  //--- Read all likelihood values ! ---//
  std::ifstream ifs ("weights.txt", std::ifstream::in);
  //ifs.open("weights.txt");
  std::cout << " Value of ifs : " << ifs.eof() << std::endl;
  std::string line;
  int evt,config,tf;
  double weight, weightUnc;
  bool weightFound = 0;
  while( std::getline(ifs,line) ){
    std::istringstream iss(line);
    if( iss >> evt >> config >> tf >> weight >> weightUnc){
      std::cout << " Looking at event : " << evt << std::endl;
      stringstream ssEvt; ssEvt << evt; string sEvt = ssEvt.str();
      //char *cEvt = itoa(evt); string sEvt = string(cEvt);
      //--- Initialize the event-per-event variables! ---//
      if( config == 1){
        vector<double> aHat(2,0.0), aHatXS(2,0.0), aHatAcc(2,0.0);
        vector<double> bHat(2,0.0), bHatXS(2,0.0), bHatAcc(2,0.0);
        vector<double> cHat(2,0.0), cHatXS(2,0.0), cHatAcc(2,0.0);
        vector<double> LnLik(NrConfigs,0.0), LnLikXS(NrConfigs,0.0), LnLikAcc(NrConfigs,0.0);
        vector<double> Lik(NrConfigs,0.0),   LikXS(NrConfigs,0.0),   LikAcc(NrConfigs,0.0);

        LnLikDist->SetName(("LnLik_Evt"+sEvt).c_str());       LnLikDist->SetTitle(("LnLik distribution for event "+sEvt+" -- "+title+" evts").c_str());
        LnLikXSDist->SetName(("LnLikXS_Evt"+sEvt).c_str());   LnLikXSDist->SetTitle(("LnLikXS distribution for event "+sEvt+" -- "+title+" evts").c_str());
        LnLikAccDist->SetName(("LnLikAcc_Evt"+sEvt).c_str()); LnLikAccDist->SetTitle(("LnLikAcc distribution for event "+sEvt+" -- "+title+" evts").c_str());
      }
      Lik[config-1] = weight;         LikXS[config-1] = weight/MGXS[config-1];              LikAcc[config-1] = weight/MGXSCut[config-1];
      LnLik[config-1] = -log(weight); LnLikXS[config-1] = -log(weight)+log(MGXS[config-1]); LnLikAcc[config-1] = -log(weight)+log(MGXSCut[config-1]);

      //---  Fill the LnLik histograms for each event and for all events together  ---//
      LnLikDist->SetBinContent(   LnLikDist->FindBin(Var[config-1]),    LnLik[config-1]);
      LnLikXSDist->SetBinContent( LnLikXSDist->FindBin(Var[config-1]),  LnLikXS[config-1]);
      LnLikAccDist->SetBinContent(LnLikAccDist->FindBin(Var[config-1]), LnLikAcc[config-1]);
      LnLikAll->SetBinContent(   LnLikAll->FindBin(Var[config-1]),    LnLikAll->GetBinContent(LnLikAll->FindBin(Var[config-1]))       + LnLik[config-1]   );
      LnLikXSAll->SetBinContent( LnLikXSAll->FindBin(Var[config-1]),  LnLikXSAll->GetBinContent(LnLikXSAll->FindBin(Var[config-1]))   + LnLikXS[config-1] );
      LnLikAccAll->SetBinContent(LnLikAccAll->FindBin(Var[config-1]), LnLikAccAll->GetBinContent(LnLikAccAll->FindBin(Var[config-1])) + LnLikAcc[config-1]);

      //---  Only perform the fit after all configurations are considered!  ---//
      if( config == NrConfigs){
        //vector<double> LnLikFunction[NrConfigs]
          //LnLikFunction = [array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d')]
          //LnLikXSFunction = [array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d')]
          //LnLikAccFunction = [array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d')]

        for(int ii=0; ii < 2; ii++){
          std::cout << " Looking at xPos and xNeg values of : " << xPos[ii] << " & " << xNeg[ii] << std::endl;
          std::cout << " LnLik values used : " << LnLik[xNeg[ii]] << " , " << LnLik[xMin] << " & " << LnLik[xPos[ii]] << std::endl;
          cHat[ii] = LnLik[xNeg[ii]]*Var[xPos[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]])) - LnLik[xMin]*Var[xNeg[ii]]*Var[xPos[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLik[xPos[ii]]*Var[xNeg[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]));
          cHatXS[ii] = LnLikXS[xNeg[ii]]*Var[xPos[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]])) - LnLikXS[xMin]*Var[xNeg[ii]]*Var[xPos[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLikXS[xPos[ii]]*Var[xNeg[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]));
          cHatAcc[ii] = LnLikAcc[xNeg[ii]]*Var[xPos[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]])) - LnLikAcc[xMin]*Var[xNeg[ii]]*Var[xPos[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLikAcc[xPos[ii]]*Var[xNeg[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]));

          bHat[ii] = LnLik[xMin]*(Var[xPos[ii]]+Var[xNeg[ii]])/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLik[xPos[ii]]*(Var[xNeg[ii]]+Var[xMin])/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]))-LnLik[xNeg[ii]]*(Var[xMin]+Var[xPos[ii]])/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]]));
          bHatXS[ii] = LnLikXS[xMin]*(Var[xPos[ii]]+Var[xNeg[ii]])/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLikXS[xPos[ii]]*(Var[xNeg[ii]]+Var[xMin])/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]))-LnLikXS[xNeg[ii]]*(Var[xMin]+Var[xPos[ii]])/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]]));
          bHatAcc[ii] = LnLikAcc[xMin]*(Var[xPos[ii]]+Var[xNeg[ii]])/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLikAcc[xPos[ii]]*(Var[xNeg[ii]]+Var[xMin])/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]))-LnLikAcc[xNeg[ii]]*(Var[xMin]+Var[xPos[ii]])/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]]));

          aHat[ii] = LnLik[xPos[ii]]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLik[xMin]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLik[xNeg[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xNeg[ii]]));
          aHatXS[ii] = LnLikXS[xPos[ii]]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLikXS[xMin]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLikXS[xNeg[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xNeg[ii]]));
          aHatAcc[ii] = LnLikAcc[xPos[ii]]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLikAcc[xMin]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLikAcc[xNeg[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xNeg[ii]]));

          //for( int iConfig = 0; iConfig < NrConfigs; iConfig++){
          //  LnLikFunction[ii].append( aHat[ii]*Var[var]*Var[var]+bHat[ii]*Var[var]+cHat[ii] );
          //}

        }

      }
    }
  }
  ifs.close();


  /*for(int iEvt = 0; iEvt < nEvts; iEvt++){
    vector<float> aHat = {0,0}, aHatXS = {0,0}, aHatAcc = {0,0};
    float bHat[2] = {0,0}, bHatXS[2] = {0,0}, bHatAcc[2]={0,0};
    float cHat[2], cHatXS[2], cHatAcc[2];
    vector<double> LnLik(NrConfigs,0.0), LnLikXS(NrConfigs,0.0), LnLikAcc(NrConfigs,0.0);
    vector<double> Lik(NrConfigs,0.0),   LikXS(NrConfigs,0.0),   LikAcc(NrConfigs,0.0);
    for(int iConfig = 0; iConfig < NrConfigs; iConfig++){
      //LnLik[iConfig] = readFromFile(iEvt+1, iConfig+1);
    }
    std::cout << iEvt+1 << ") Stored LnLik values are : " << LnLik[0] << ", " << LnLik[1] << " , " << LnLik[2] << " , " << LnLik[3] << " & " << LnLik[4] << std::endl;
  }*/

}


