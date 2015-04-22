#! python
import os
import sys
from math import log
from array import array
from ROOT import TH1F,TH2F,TFile,TCanvas,TLegend,gStyle,gROOT,TGraph
gROOT.SetBatch(True)                       #Don't print the histograms (necessary for fit)

#Get all the input from the command line:
if len(sys.argv) <= 1:
  print "Need to give the directory of interest, the considered kinematic variable and the number of events in the command line !"
  print " Correct syntax is : python FitDeviation.py Events/blabla/ MTop 100"
  sys.exit()
elif len(sys.argv) == 2:
  print "Need to specify the considered kinematic variable (MTop or RVR)"
  print " Correct syntax is : python FitDeviation.py Events/blabla/ MTop 100"
  sys.exit()
elif len(sys.argv) == 3:
  print "Need to give the number of considered events !"
  print " Correct syntax is : python FitDeviation.py Events/blabla/ MTop 100"
  sys.exit()

whichDir = sys.argv[1]
KinVariable = sys.argv[2]
nEvts = sys.argv[3]

if KinVariable != "MTop" and KinVariable != "RVR":
  print "Need to specify which kinematic variable should be considered (MTop or RVR are the only options!!)"
  KinVariable = raw_input('--> Choose one of the two : ')
print "Interested in directory :",whichDir," (using ",nEvts," events)"

if KinVariable == "RVR":
  #Information about the scanned RVR values and the corresponding cross-section
  VarValues =  ["Re(V_{R}) = -0.5","Re(V_{R}) = -0.3","Re(V_{R}) = -0.2","Re(V_{R}) = -0.1","Re(V_{R}) = -0.05","Re(V_{R}) = 0.0","Re(V_{R}) = 0.05","Re(V_{R}) = 0.1","Re(V_{R}) = 0.2","Re(V_{R}) = 0.3", "Re(V_{R}) = 0.5"]
  Var =        array('d',[-0.5,            -0.3,             -0.2,            -0.1,             -0.05,    0.0,              0.05,     0.1,              0.2,              0.3,              0.5            ])
  MGXS =       array('d',[17.9275,         13.3944,          12.06555,        11.25909,         11.02784, 10.90059,         10.88228, 10.97767,         11.49883,         12.49056,         16.1508        ])
  MGXSe =      array('d',[0.0123100357311, 0.00995808028337, 0.0093464076837, 0.00836607833038, 0.0,      0.00822214433527, 0.0,      0.00847293509122, 0.00901976602967, 0.00874682264197, 0.0113081652137])
  Acceptance = array('d',[0.22164,         0.21742,          0.21672,         0.21737,          0.21614,  0.21670,          0.21531,  0.21677,          0.21437,          0.21793,          0.22205        ])

  #Select which window of RVR values was considered!
  VarWindow = raw_input('** Choose the correct RVR-window corresponding to the studied file ** \n** Different options are : \n  1) Wide   : [-0.5, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.5] \n  2) Narrow : [-0.3, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2, 0.3] \n --> Choose the correct number : ')

  xMinValue = [4, 4]
  KinVar = "Re(V_{R})"

  if VarWindow == "1":
    VarValues.pop(6), Var.pop(6), MGXS.pop(6), MGXSe.pop(6), Acceptance.pop(6)
    VarValues.pop(4), Var.pop(4), MGXS.pop(4), MGXSe.pop(4), Acceptance.pop(4)
    xBin, xLow, xHigh = 11, -0.55, 0.55
    xBinZoom, xLowZoom, xHighZoom = 7, -0.35, 0.35 
  elif VarWindow == "2":
    VarValues.pop(10), Var.pop(10), MGXS.pop(10), MGXSe.pop(10), Acceptance.pop(10)
    VarValues.pop(0), Var.pop(0), MGXS.pop(0), MGXSe.pop(0), Acceptance.pop(0)
    xBin, xLow, xHigh = 13, -0.325, 0.325
    xBinZoom, xLowZoom, xHighZoom = xBin, xLow, xHigh

elif KinVariable == "MTop":
  #Information about the scanned MTop values and the corresponding cross-section
  VarValues = ["m_{top} = 153",  "m_{top} = 163",  "m_{top} = 170",  "m_{top} = 171",  "m_{top} = 172",  "m_{top} = 173",  "m_{top} = 174",  "m_{top} = 175","m_{top} = 183",  "m_{top} = 193"]
  Var =        array('d',[153,              163,              170,              171,              172,              173,              174,              175,            183,              193            ])
  MGXS =       array('d',[8.20916,          9.6299,           10.57123,         10.70485,         10.8257,          10.96469,         11.08428,         11.22448,       12.18068,         13.3046        ])
  MGXSe =      array('d',[0.00641266950107, 0.00775899974932, 0.00857143063963, 0.00814303595657, 0.00878899028501, 0.00816801717126, 0.00904797742371, 0.008653800078, 0.00931290317946, 0.0103310001752])
  Acceptance = array('d',[0.16203,          0.19152,          0.21008,          0.21460,          0.21735,          0.21290,          0.21752,          0.22185,        0.23941,          0.26413        ])

  #Select which window of masses was considered!
  VarWindow = raw_input('** Choose the correct mass-window corresponding to the studied file ** \n** Different options are : \n  1) Wide   : [153, 163, 170, 171, 172, 173, 174, 175, 183, 193] \n  2) Normal : [153, 163, 171, 172, 173, 174, 175, 183, 193] \n  3) Narrow : [171, 172, 173, 174, 175] \n --> Choose the correct number : ')

  xMinValue = [5, 4, 2]
  KinVar = "m_{top}"

  if VarWindow == "1":
    xBin, xLow, xHigh = 41, 152.5, 193.5
  elif VarWindow == "2":
    VarValues.pop(2), Var.pop(2), MGXS.pop(2), MGXSe.pop(2), Acceptance.pop(2)
    xBin, xLow, xHigh = 41, 152.5, 193.5
  elif VarWindow == "3":
    VarValues.pop(9), Var.pop(9), MGXS.pop(9), MGXSe.pop(9), Acceptance.pop(9)
    VarValues.pop(8), Var.pop(8), MGXS.pop(8), MGXSe.pop(8), Acceptance.pop(8)
    VarValues.pop(2), Var.pop(2), MGXS.pop(2), MGXSe.pop(2), Acceptance.pop(2)
    VarValues.pop(1), Var.pop(1), MGXS.pop(1), MGXSe.pop(1), Acceptance.pop(1)
    VarValues.pop(0), Var.pop(0), MGXS.pop(0), MGXSe.pop(0), Acceptance.pop(0)
    xBin, xLow, xHigh = 5, 170.5, 175.5
  xBinZoom, xLowZoom, xHighZoom = 5, 170.5, 175.5

print " List of considered Var values is : ",Var    
NrConfigs = len(Var)
xPos = [xMinValue[int(VarWindow)-1]+1,xMinValue[int(VarWindow)-1]+2]
xNeg = [xMinValue[int(VarWindow)-1]-1,xMinValue[int(VarWindow)-1]-2]
xMin = xMinValue[int(VarWindow)-1]
xStep = [Var[xNeg[0]]-Var[xNeg[1]], Var[xMin]-Var[xNeg[0]], Var[xPos[0]]-Var[xMin], Var[xPos[1]]-Var[xPos[0]] ]
print "Step size = ",xStep

#File of interest:
list_dir = os.listdir(whichDir)
WeightsFileArray = []
weightsFileCounter = 0
for file in list_dir:
  if file.endswith(".out") and file.startswith("weights"): # eg: '.txt'
    weightsFileCounter += 1
    WeightsFileArray.append(file)

if int(weightsFileCounter) == 1:
  WeightsFile = open(os.path.join(whichDir+''+WeightsFileArray[0]),'r')
  LikelihoodFile = open(os.path.join(whichDir+''+WeightsFileArray[0]),'r')
elif int(weightsFileCounter) == 0:
  print "No weights file found in this directory !"
  sys.exit()
elif int(weightsFileCounter) > 1:
  for ii in range(len(WeightsFileArray)):
    print " ",ii," ) ",WeightsFileArray[ii]
  fileNr = raw_input('Choose the number of the file of interest! : ')
  WeightsFile = open(os.path.join(whichDir+''+WeightsFileArray[int(fileNr)]),'r')
  LikelihoodFile = open(os.path.join(whichDir+''+WeightsFileArray[int(fileNr)]),'r')
print "Will be using file : ",WeightsFile

#Set title of root file!
title = ""
if whichDir.find("Correct") <= len(whichDir)     and whichDir.find("Correct") > 0:   title = "Correct"
elif whichDir.find("Wrong") <= len(whichDir)     and whichDir.find("Wrong") > 0:     title = "Wrong"
elif whichDir.find("Unmatched") <= len(whichDir) and whichDir.find("Unmatched") > 0: title = "Unmatched"

if whichDir.find("Reco") <= len(whichDir)  and whichDir.find("Reco") > 0: title += "Reco"
elif whichDir.find("Gen") <= len(whichDir) and whichDir.find("Gen") > 0:  title += "Gen"
title = title+"_"+KinVariable

#ROOT file where all the information will be stored:
Tfile = TFile(os.path.join(whichDir+"FitDeviation_"+title+".root"),'recreate')
gStyle.SetOptStat(0)

YPlusGausTest    = TH1F('YPlusGausTest',   'Comparison of fit deviation in '+KinVar+' = '+str(Var[xPos[1]])+' and '+KinVar+' = '+str(Var[xPos[0]])+' (no normalisation -- '+title+' evts)', 250,-0.15,0.15)
YPlusGausTestXS  = TH1F('YPlusGausTestXS', 'Comparison of fit deviation in '+KinVar+' = '+str(Var[xPos[1]])+' and '+KinVar+' = '+str(Var[xPos[0]])+' (XS normalisation -- '+title+' evts)', 250,-0.15,0.15)
YPlusGausTestAcc = TH1F('YPlusGausTestAcc','Comparison of fit deviation in '+KinVar+' = '+str(Var[xPos[1]])+' and '+KinVar+' = '+str(Var[xPos[0]])+' (Acc normalisation -- '+title+' evts)',250,-0.15,0.15)
YPlusGausTestPosScdDer    = TH1F('YPlusGausTestPosScdDer',   'Gaussian fit deviation using '+KinVar+' = '+str(Var[xPos[1]])+' & '+str(Var[xPos[0]])+' (no norm -- outer 2nd der > 0 -- '+title+' evts)', 250,-0.1,0.1)
YPlusGausTestXSPosScdDer  = TH1F('YPlusGausTestXSPosScdDer', 'Gaussian fit deviation using '+KinVar+' = '+str(Var[xPos[1]])+' & '+str(Var[xPos[0]])+' (XS norm -- outer 2nd der > 0 -- '+title+' evts)', 250,-0.1,0.1)
YPlusGausTestAccPosScdDer = TH1F('YPlusGausTestAccPosScdDer','Gaussian fit deviation using '+KinVar+' = '+str(Var[xPos[1]])+' & '+str(Var[xPos[0]])+' (Acc norm -- outer 2nd der > 0 -- '+title+' evts)',250,-0.1,0.1)
YPlusGausTestNegScdDer    = TH1F('YPlusGausTestNegScdDer',   'Gaussian fit deviation using '+KinVar+' = '+str(Var[xPos[1]])+' & '+str(Var[xPos[0]])+' (no norm -- outer 2nd der < 0 -- '+title+' evts)', 250,-0.1,0.1)
YPlusGausTestXSNegScdDer  = TH1F('YPlusGausTestXSNegScdDer', 'Gaussing fit deviation using '+KinVar+' = '+str(Var[xPos[1]])+' & '+str(Var[xPos[0]])+' (XS norm -- outer 2nd der < 0 -- '+title+' evts)', 250,-0.1,0.1)
YPlusGausTestAccNegScdDer = TH1F('YPlusGausTestAccNegScdDer','Gaussian fit deviation using '+KinVar+' = '+str(Var[xPos[1]])+' & '+str(Var[xPos[0]])+' (Acc norm -- outer 2nd der < 0 -- '+title+' evts)',250,-0.1,0.1)

YPlus    = TH1F('YPlus',   'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[0]])+' (no normalisation -- '+title+' evts)' ,150,-0.25,0.25)
YPlusXS  = TH1F('YPlusXS', 'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[0]])+' (XS normalisation -- '+title+' evts)' ,150,-0.25,0.25)
YPlusAcc = TH1F('YPlusAcc','Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[0]])+' (Acc normalisation -- '+title+' evts)',150,-0.25,0.25)
YPlusPosScdDer    = TH1F('YPlusPosScdDer',   'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[0]])+' (no normalisation -- outer 2nd der > 0 -- '+title+' evts)' ,150,-0.2,0.2)
YPlusXSPosScdDer  = TH1F('YPlusXSPosScdDer', 'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[0]])+' (XS normalisation -- outer 2nd der > 0 -- '+title+' evts)' ,150,-0.2,0.2)
YPlusAccPosScdDer = TH1F('YPlusAccPosScdDer','Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[0]])+' (Acc normalisation -- outer 2nd der > 0 -- '+title+' evts)',150,-0.2,0.2)
YPlusNegScdDer    = TH1F('YPlusNegScdDer',   'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[0]])+' (no normalisation -- outer 2nd der < 0 -- '+title+' evts)' ,150,-0.2,0.2)
YPlusXSNegScdDer  = TH1F('YPlusXSNegScdDer', 'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[0]])+' (XS normalisation -- outer 2nd der < 0 -- '+title+' evts)' ,150,-0.2,0.2)
YPlusAccNegScdDer = TH1F('YPlusAccNegScdDer','Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[0]])+' (Acc normalisation -- outer 2nd der < 0 -- '+title+' evts)',150,-0.2,0.2)

YPlusPlus    = TH1F('YPlusPlus',   'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[1]])+' (no normalisation -- '+title+' evts)' ,150,-0.5,0.5)
YPlusPlusXS  = TH1F('YPlusPlusXS', 'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[1]])+' (XS normalisation -- '+title+' evts)' ,150,-0.5,0.5)
YPlusPlusAcc = TH1F('YPlusPlusAcc','Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[1]])+' (Acc normalisation -- '+title+' evts)',150,-0.5,0.5)
YPlusPlusPosScdDer    = TH1F('YPlusPlusPosScdDer',   'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[1]])+' (no normalisation -- outer 2nd der > 0 -- '+title+' evts)' ,150,-0.3,0.3)
YPlusPlusXSPosScdDer  = TH1F('YPlusPlusXSPosScdDer', 'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[1]])+' (XS normalisation -- outer 2nd der > 0 -- '+title+' evts)' ,150,-0.3,0.3)
YPlusPlusAccPosScdDer = TH1F('YPlusPlusAccPosScdDer','Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[1]])+' (Acc normalisation -- outer 2nd der > 0 -- '+title+' evts)',150,-0.3,0.3)
YPlusPlusNegScdDer    = TH1F('YPlusPlusNegScdDer',   'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[1]])+' (no normalisation -- outer 2nd der < 0 -- '+title+' evts)' ,150,-0.3,0.3)
YPlusPlusXSNegScdDer  = TH1F('YPlusPlusXSNegScdDer', 'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[1]])+' (XS normalisation -- outer 2nd der < 0 -- '+title+' evts)' ,150,-0.3,0.3)
YPlusPlusAccNegScdDer = TH1F('YPlusPlusAccNegScdDer','Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[1]])+' (Acc normalisation -- outer 2nd der < 0 -- '+title+' evts)',150,-0.3,0.3)

YMin    = TH1F('YMin',   'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (no normalisation -- '+title+' evts)' ,150,-0.25,0.25)
YMinXS  = TH1F('YMinXS', 'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (XS normalisation -- '+title+' evts)' ,150,-0.25,0.25)
YMinAcc = TH1F('YMinAcc','Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (Acc normalisation -- '+title+' evts)',150,-0.25,0.25)
YMinPosScdDer    = TH1F('YMinPosScdDer',   'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (no normalisation -- outer 2nd der > 0 -- '+title+' evts)' ,150,-0.2,0.2)
YMinXSPosScdDer  = TH1F('YMinXSPosScdDer', 'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (XS normalisation -- outer 2nd der > 0 -- '+title+' evts)' ,150,-0.2,0.2)
YMinAccPosScdDer = TH1F('YMinAccPosScdDer','Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (Acc normalisation -- outer 2nd der > 0 -- '+title+' evts)',150,-0.2,0.2)
YMinNegScdDer    = TH1F('YMinNegScdDer',   'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (no normalisation -- outer 2nd der < 0 -- '+title+' evts)' ,150,-0.2,0.2)
YMinXSNegScdDer  = TH1F('YMinXSNegScdDer', 'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (XS normalisation -- outer 2nd der < 0 -- '+title+' evts)' ,150,-0.2,0.2)
YMinAccNegScdDer = TH1F('YMinAccNegScdDer','Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (Acc normalisation -- outer 2nd der < 0 -- '+title+' evts)',150,-0.2,0.2)

YMinMin    = TH1F('YMinMin',   'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[1]])+' (no normalisation -- '+title+' evts)' ,150,-0.5,0.5)
YMinMinXS  = TH1F('YMinMinXS', 'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[1]])+' (XS normalisation -- '+title+' evts)' ,150,-0.5,0.5)
YMinMinAcc = TH1F('YMinMinAcc','Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[1]])+' (Acc normalisation -- '+title+' evts)',150,-0.5,0.5)
YMinMinPosScdDer    = TH1F('YMinMinPosScdDer',   'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[1]])+' (no normalisation -- outer 2nd der > 0 -- '+title+' evts)' ,150,-0.3,0.3)
YMinMinXSPosScdDer  = TH1F('YMinMinXSPosScdDer', 'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[1]])+' (XS normalisation -- outer 2nd der > 0 -- '+title+' evts)' ,150,-0.3,0.3)
YMinMinAccPosScdDer = TH1F('YMinMinAccPosScdDer','Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[1]])+' (Acc normalisation -- outer 2nd der > 0 -- '+title+' evts)',150,-0.3,0.3)
YMinMinNegScdDer    = TH1F('YMinMinNegScdDer',   'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[1]])+' (no normalisation -- outer 2nd der < 0 -- '+title+' evts)' ,150,-0.3,0.3)
YMinMinXSNegScdDer  = TH1F('YMinMinXSNegScdDer', 'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[1]])+' (XS normalisation -- outer 2nd der < 0 -- '+title+' evts)' ,150,-0.3,0.3)
YMinMinAccNegScdDer = TH1F('YMinMinAccNegScdDer','Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[1]])+' (Acc normalisation -- outer 2nd der < 0 -- '+title+' evts)',150,-0.3,0.3)

YRelPlus    = TH1F('YRelPlus',   'Relative deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (no normalisation -- '+title+' evts)' ,250,-0.1,0.1)
YRelPlusXS  = TH1F('YRelPlusXS', 'Relative deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (XS normalisation -- '+title+' evts)' ,250,-0.1,0.1)
YRelPlusAcc = TH1F('YRelPlusAcc','Relative deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (Acc normalisation -- '+title+' evts)',250,-0.1,0.1)

YRelMin    = TH1F('YRelMin',   'Relative deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (no normalisation -- '+title+' evts)' ,150,-0.1,0.1)
YRelMinXS  = TH1F('YRelMinXS', 'Relative deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (XS normalisation -- '+title+' evts)' ,150,-0.1,0.1)
YRelMinAcc = TH1F('YRelMinAcc','Relative deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (Acc normalisation -- '+title+' evts)',150,-0.1,0.1)

FstDer    = TH1F('FirstDer',   'First derivative of -ln(likelihood) distribution -- '+title+' evts', 4,0,4)
FstDer.GetXaxis().SetBinLabel(1,'(y_{DATA}(x='+str(Var[xNeg[1]])+') - y_{DATA}(x='+str(Var[xNeg[0]])+'))/'+str(xStep[0]))
FstDer.GetXaxis().SetBinLabel(2,'(y_{DATA}(x='+str(Var[xNeg[0]])+') - y_{DATA}(x='+str(Var[xMin])+'))/'+str(xStep[1]))
FstDer.GetXaxis().SetBinLabel(3,'(y_{DATA}(x='+str(Var[xMin])+') - y_{DATA}(x='+str(Var[xPos[0]])+'))/'+str(xStep[2]))
FstDer.GetXaxis().SetBinLabel(4,'(y_{DATA}(x='+str(Var[xPos[0]])+') - y_{DATA}(x='+str(Var[xPos[1]])+'))/'+str(xStep[3]))
FstDerXS  = TH1F('FirstDerivativeXS', 'First derivative of -ln(likelihood) distribution (XS normalisation -- '+title+' evts)', 5,-0.25,0.25)
FstDerAcc = TH1F('FirstDerivativeAcc','First derivative of -ln(likelihood) distribution (Acc normalisation -- '+title+' evts)',5,-0.25,0.25)
ScdDerInner    = TH1F('SecondDerivativeInner',   'Second derivative of -ln(likelihood) distribution (no normalisation -- using inner points -- '+title+' evts)', 250,-5,5)
ScdDerXSInner  = TH1F('SecondDerivativeXSInner', 'Second derivative of -ln(likelihood) distribution (XS normalisation -- using inner points -- '+title+' evts)', 250,-5,5)
ScdDerAccInner = TH1F('SecondDerivativeAccInner','Second derivative of -ln(likelihood) distribution (Acc normalisation -- using inner points -- '+title+' evts)',250,-5,5)
ScdDerOuter    = TH1F('SecondDerivativeOuter',   'Second derivative of -ln(likelihood) distribution (no normalisation -- using outer points -- '+title+' evts)', 250,-5,5)
ScdDerXSOuter  = TH1F('SecondDerivativeXSOuter', 'Second derivative of -ln(likelihood) distribution (XS normalisation -- using outer points -- '+title+' evts)', 250,-5,5)
ScdDerAccOuter = TH1F('SecondDerivativeAccOuter','Second derivative of -ln(likelihood) distribution (Acc normalisation -- using outer points -- '+title+' evts)',250,-5,5)
ScdDerScatter    = TH2F('ScdDerScatterPlot',   'Second derivative of -ln(L) using inner points versus using outer points (no normalisation -- '+title+' evts)', 250,-5,5,250,-5,5)
ScdDerXSScatter  = TH2F('ScdDerXSScatterPlot', 'Second derivative of -ln(L) using inner points versus using outer points (XS normalisation -- '+title+' evts)', 250,-5,5,250,-5,5)
ScdDerAccScatter = TH2F('ScdDerAccScatterPlot','Second derivative of -ln(L) using inner points versus using outer points (Acc normalisation -- '+title+' evts)',250,-5,5,250,-5,5)

LnLogDist = TH1F("LnLog","title",xBin,xLow,xHigh)
LnLogDist.SetMarkerStyle(20), LnLogDist.SetLineColor(1), LnLogDist.SetMarkerColor(1), LnLogDist.SetMarkerSize(1.2)
LnLogXSDist = TH1F("LnLogXS","title",xBin,xLow,xHigh)
LnLogXSDist.SetMarkerStyle(21), LnLogXSDist.SetLineColor(3), LnLogXSDist.SetMarkerColor(3), LnLogXSDist.SetMarkerSize(1.2)
LnLogAccDist = TH1F("LnLogAcc","title",xBin,xLow,xHigh)
LnLogAccDist.SetMarkerStyle(22), LnLogAccDist.SetLineColor(4), LnLogAccDist.SetMarkerColor(4), LnLogAccDist.SetMarkerSize(1.2)

FitComp = Tfile.mkdir("FitComparison")
LnLogDir = Tfile.mkdir("LnLogDistr")
LnLogXSDir = Tfile.mkdir("LnLogXSDistr")
LnLogAccDir = Tfile.mkdir("LnLogAccDistr")
FstDerDir = Tfile.mkdir("FirstDerivativeDistr")
LnLogFitCanvas = TCanvas('name','title')
LnLogFitCanvasName = 'LnLogFitCanvas_Evt'
LnLogFitCanvasTitle = 'Comparing fit function from ROOT fit and algebraic function for event -- '+title+' evts '

#Create the arrays where the likelihood values will be stored
LnLog, LnLogXS, LnLogAcc = [], [], []
for ii in range(NrConfigs):
  LnLog.append(0), LnLogXS.append(0), LnLogAcc.append(0)

#---  Create arrays where the events passing specific cuts will be stored  ---#
EvtsWithPosScdDerInner, EvtsWithPosScdDerXSInner, EvtsWithPosScdDerAccInner = [], [], []
EvtsWithPosScdDerOuter, EvtsWithPosScdDerXSOuter, EvtsWithPosScdDerAccOuter = [], [], []
EvtsPosScdDerInner, EvtsPosScdDerXSInner, EvtsPosScdDerAccInner = [], [], []
EvtsPosScdDerOuter, EvtsPosScdDerXSOuter, EvtsPosScdDerAccOuter = [], [], []
EvtsWithYPlusGausSmall, EvtsWithYPlusGausSmallXS, EvtsWithYPlusGausSmallAcc = [], [], []

#Loop over all lines in weights file:
for WeightLine in WeightsFile:
  WeightWord = WeightLine.split()
  #Only interested in files starting with a number
  if str(WeightWord[0]) != "#" and str(WeightWord[3]) != "0.0" :
    for iEvt in range(int(nEvts)):
      if str(WeightWord[0]) == str(iEvt+1):                                                               #Look at one single event!
        if str(WeightWord[1]) == "1":
          cHat, cHatXS, cHatAcc = [0,0],[0,0],[0,0]
          bHat, bHatXS, bHatAcc = [0,0],[0,0],[0,0]
          aHat, aHatXS, aHatAcc = [0,0],[0,0],[0,0]
          LnLogDist.SetName("LnLog_Evt"+str(int(iEvt)+1)),       LnLogDist.SetTitle('LnLog distribution for event '+str(int(iEvt)+1)+' -- '+title+' evts')
          LnLogXSDist.SetName("LnLogXS_Evt"+str(int(iEvt)+1)),   LnLogXSDist.SetTitle('LnLogXS distribution for event '+str(int(iEvt)+1)+' -- '+title+' evts')
          LnLogAccDist.SetName("LnLogAcc_Evt"+str(int(iEvt)+1)), LnLogAccDist.SetTitle('LnLogAcc distribution for event '+str(int(iEvt)+1)+' -- '+title+' evts')
        LnLog[int(WeightWord[1])-1] = -log(float(WeightWord[3]))
        LnLogXS[int(WeightWord[1])-1] = -log(float(WeightWord[3])) + log(MGXS[int(WeightWord[1])-1])
        LnLogAcc[int(WeightWord[1])-1] = -log(float(WeightWord[3])) + log(MGXS[int(WeightWord[1])-1]) + log(Acceptance[int(WeightWord[1])-1])
        #---  Fill the LnLog histograms for each event  ---#
        LnLogDist.SetBinContent(LnLogDist.FindBin(Var[int(WeightWord[1])-1]), LnLog[int(WeightWord[1])-1])
        LnLogXSDist.SetBinContent(LnLogXSDist.FindBin(Var[int(WeightWord[1])-1]), LnLogXS[int(WeightWord[1])-1])
        LnLogAccDist.SetBinContent(LnLogAccDist.FindBin(Var[int(WeightWord[1])-1]), LnLogAcc[int(WeightWord[1])-1])
        #---  Only perform the fit after all configurations are considered!  ---#
        if str(WeightWord[1]) == str(NrConfigs):
          #---  Calculate the fit parameters (polynomial)  ---#
          LnLogFunction = [array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d')]
          LnLogXSFunction = [array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d')]
          LnLogAccFunction = [array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d')]
          for ii in range(2):
            cHat[ii] = LnLog[xNeg[ii]]*Var[xPos[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]])) - LnLog[xMin]*Var[xNeg[ii]]*Var[xPos[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLog[xPos[ii]]*Var[xNeg[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]))
            cHatXS[ii] = LnLogXS[xNeg[ii]]*Var[xPos[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]])) - LnLogXS[xMin]*Var[xNeg[ii]]*Var[xPos[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLogXS[xPos[ii]]*Var[xNeg[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]))
            cHatAcc[ii] = LnLogAcc[xNeg[ii]]*Var[xPos[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]])) - LnLogAcc[xMin]*Var[xNeg[ii]]*Var[xPos[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLogAcc[xPos[ii]]*Var[xNeg[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]))

            bHat[ii] = LnLog[xMin]*(Var[xPos[ii]]+Var[xNeg[ii]])/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLog[xPos[ii]]*(Var[xNeg[ii]]+Var[xMin])/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]))-LnLog[xNeg[ii]]*(Var[xMin]+Var[xPos[ii]])/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]]))
            bHatXS[ii] = LnLogXS[xMin]*(Var[xPos[ii]]+Var[xNeg[ii]])/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLogXS[xPos[ii]]*(Var[xNeg[ii]]+Var[xMin])/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]))-LnLogXS[xNeg[ii]]*(Var[xMin]+Var[xPos[ii]])/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]]))
            bHatAcc[ii] = LnLogAcc[xMin]*(Var[xPos[ii]]+Var[xNeg[ii]])/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLogAcc[xPos[ii]]*(Var[xNeg[ii]]+Var[xMin])/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]))-LnLogAcc[xNeg[ii]]*(Var[xMin]+Var[xPos[ii]])/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]]))

            aHat[ii] = LnLog[xPos[ii]]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLog[xMin]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLog[xNeg[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xNeg[ii]]))
            aHatXS[ii] = LnLogXS[xPos[ii]]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLogXS[xMin]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLogXS[xNeg[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xNeg[ii]]))
            aHatAcc[ii] = LnLogAcc[xPos[ii]]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLogAcc[xMin]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLogAcc[xNeg[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xNeg[ii]]))

            #-- Now calculate each of the 9 Var configurations using this parabola --#
            for rvr in range(NrConfigs):
              LnLogFunction[ii].append(aHat[ii]*Var[rvr]*Var[rvr]+bHat[ii]*Var[rvr]+cHat[ii])
              LnLogXSFunction[ii].append( aHatXS[ii]*Var[rvr]*Var[rvr]+bHatXS[ii]*Var[rvr]+cHatXS[ii])
              LnLogAccFunction[ii].append( aHatAcc[ii]*Var[rvr]*Var[rvr]+bHatAcc[ii]*Var[rvr]+cHatAcc[ii])
          LnLogGraphInner = TGraph(9, Var, LnLogFunction[0])
          LnLogGraphInner.SetTitle('Comparison between ROOT fit and algebraic method -- '+title+' evts')
          LnLogGraphInner.GetYaxis().SetTitle('#Delta ln(Likelihood)')
          LnLogGraphInner.GetXaxis().SetTitle(KinVar)
          LnLogGraphInner.SetMarkerColor(6), LnLogGraphInner.SetLineColor(6), LnLogGraphInner.SetMarkerSize(1.2)
          LnLogGraphOuter = TGraph(9,Var, LnLogFunction[1])
          LnLogGraphOuter.SetMarkerColor(7), LnLogGraphOuter.SetLineColor(7), LnLogGraphOuter.SetMarkerSize(1.2)	
          #---  Compare with TF1 fit from ROOT  ---#
          LnLogDist.Fit("pol2","Q")
          LnLogFitCanvas.SetName(LnLogFitCanvasName+str(int(iEvt)+1)), LnLogFitCanvas.SetTitle(LnLogFitCanvasTitle+str(int(iEvt)+1))
          legendLnLog = TLegend(0.75,0.7,0.95,0.95)
          LnLogFitCanvas.cd()
          legendLnLog.AddEntry(LnLogGraphInner,'Fit result using algebraic function (x = '+str(Var[xNeg[0]])+'/'+str(Var[xMin])+'/'+str(Var[xPos[0]])+')',"p")
          legendLnLog.AddEntry(LnLogGraphOuter,'Fit result using algebraic function (x = '+str(Var[xNeg[1]])+'/'+str(Var[xMin])+'/'+str(Var[xPos[1]])+')',"p")
          legendLnLog.AddEntry(LnLogDist,'Obtained LnLog values with ROOT fit',"p")
          LnLogGraphInner.Draw("AC*"), LnLogGraphOuter.Draw("C*"), LnLogDist.Draw("samep"), legendLnLog.Draw()
          FitComp.cd(), LnLogFitCanvas.Write()

          #---  Calculate the first derivative distribution and save them in event-by-event plot  ---#
          FstDer.SetName('FirstDerivative_Evt'+str(int(iEvt)+1))
          FstDer.SetTitle('First derivative of -ln(likelihood) distribution for event'+str(int(iEvt)+1)+' (no normalisation) -- '+title+' evts') 
          for ii in range(FstDer.GetNbinsX()):
            FstDer.SetBinContent(ii+1, (LnLog[ii+2] - LnLog[ii+3])/xStep[ii])
          FstDerDir.cd(), FstDer.Write()
          #---  Calculate the positive and negative deviation  ---#
          yPlusPlus = [LnLog[xPos[1]] - LnLogFunction[0][xPos[1]], LnLogXS[xPos[1]] - LnLogXSFunction[0][xPos[1]], LnLogAcc[xPos[1]] - LnLogAccFunction[0][xPos[1]]]
          yMinMin   = [LnLog[xNeg[1]] - LnLogFunction[0][xNeg[1]], LnLogXS[xNeg[1]] - LnLogXSFunction[0][xNeg[1]], LnLogAcc[xNeg[1]] - LnLogAccFunction[0][xNeg[1]]]
          yPlus = [LnLog[xPos[0]] - LnLogFunction[1][xPos[0]], LnLogXS[xPos[0]] - LnLogXSFunction[1][xPos[0]], LnLogAcc[xPos[0]] - LnLogAccFunction[1][xPos[0]]]
          yMin  = [LnLog[xNeg[0]] - LnLogFunction[1][xNeg[0]], LnLogXS[xNeg[0]] - LnLogXSFunction[1][xNeg[0]], LnLogAcc[xNeg[0]] - LnLogAccFunction[1][xNeg[0]]]
          scdDerInner = [LnLog[xNeg[0]]-2*LnLog[xMin]+LnLog[xPos[0]], LnLogXS[xNeg[0]]-2*LnLogXS[xMin]+LnLogXS[xPos[0]], LnLogAcc[xNeg[0]]-2*LnLogAcc[xMin]+LnLogAcc[xPos[0]]]	
          scdDerOuter   = [LnLog[xNeg[1]]-2*LnLog[xMin]+LnLog[xPos[1]], LnLogXS[xNeg[1]]-2*LnLogXS[xMin]+LnLogXS[xPos[1]], LnLogAcc[xNeg[1]]-2*LnLogAcc[xMin]+LnLogAcc[xPos[1]]]	

          #--- Fill the histograms  ---#
          YPlus.Fill(yPlus[0]),                               YPlusXS.Fill(yPlus[1]),                               YPlusAcc.Fill(yPlus[2])
          YPlusPlus.Fill(yPlusPlus[0]),                       YPlusPlusXS.Fill(yPlusPlus[1]),                       YPlusPlusAcc.Fill(yPlusPlus[2])
          YPlusGausTest.Fill(yPlus[0] + yPlusPlus[0]/4),      YPlusGausTestXS.Fill(yPlus[1] + yPlusPlus[1]/4),      YPlusGausTestAcc.Fill(yPlus[2] + yPlusPlus[2]/4)
          #YRelPlus.Fill(yPlus[0]/LnLog[6]),                   YRelPlusXS.Fill(yPlus[1]/LnLogXS[6]),                 YRelPlusAcc.Fill(yPlus[2]/LnLogAcc[6])
          YMin.Fill(yMin[0]),                                 YMinXS.Fill(yMin[1]),                                 YMinAcc.Fill(yMin[2])
          YMinMin.Fill(yMinMin[0]),                           YMinMinXS.Fill(yMinMin[1]),                           YMinMinAcc.Fill(yMinMin[2])
          #YRelMin.Fill(yMin[0]/LnLog[2]),                     YRelMinXS.Fill(yMin[1]/LnLogXS[2]),                   YRelMinAcc.Fill(yMin[2]/LnLogAcc[2])
          ScdDerInner.Fill(scdDerInner[0]),                   ScdDerXSInner.Fill(scdDerInner[1]),                   ScdDerAccInner.Fill(scdDerInner[2])
          ScdDerOuter.Fill(scdDerOuter[0]),                   ScdDerXSOuter.Fill(scdDerOuter[1]),                   ScdDerAccOuter.Fill(scdDerOuter[2])
          ScdDerScatter.Fill(scdDerOuter[0], scdDerInner[0]), ScdDerXSScatter.Fill(scdDerOuter[1], scdDerInner[1]), ScdDerAccScatter.Fill(scdDerOuter[2], scdDerInner[2])
          #---  Save the LnLog distributions (event-per-event) and the Y-deviations in histograms  ---#
          LnLogDir.cd(),    LnLogDist.Write()
          LnLogXSDir.cd(),  LnLogXSDist.Write()
          LnLogAccDir.cd(), LnLogAccDist.Write()
          #-- Apply cut on YPlusGausTest --#
          if (yPlus[0] + yPlusPlus[0]/4) <= 0.025 and (yPlus[0] + yPlusPlus[0]/4) >= -0.025: EvtsWithYPlusGausSmall.append(iEvt+1)
          if (yPlus[1] + yPlusPlus[1]/4) <= 0.025 and (yPlus[1] + yPlusPlus[1]/4) >= -0.025: EvtsWithYPlusGausSmallXS.append(iEvt+1)
          if (yPlus[2] + yPlusPlus[2]/4) <= 0.025 and (yPlus[2] + yPlusPlus[2]/4) >= -0.025: EvtsWithYPlusGausSmallAcc.append(iEvt+1)
          #-- Apply cut on ScdDer (using inner Var points) --#
          if scdDerInner[0] > 0.0: EvtsPosScdDerInner.append(iEvt+1)
          if scdDerInner[1] > 0.0: EvtsPosScdDerXSInner.append(iEvt+1)
          if scdDerInner[2] > 0.0: EvtsPosScdDerAccInner.append(iEvt+1)
          #-- Apply cut on ScdDer (using outer Var points) --#
          if scdDerOuter[0] > 0.0:
            EvtsPosScdDerOuter.append(iEvt+1)
            YPlusGausTestPosScdDer.Fill(yPlus[0] + yPlusPlus[0]/4)
            YPlusPosScdDer.Fill(yPlus[0]), YPlusPlusPosScdDer.Fill(yPlusPlus[0])
            YMinPosScdDer.Fill(yMin[0]),   YMinMinPosScdDer.Fill(yMinMin[0])
          else:
            YPlusGausTestNegScdDer.Fill(yPlus[0] + yPlusPlus[0]/4)
            YPlusNegScdDer.Fill(yPlus[0]), YPlusPlusNegScdDer.Fill(yPlusPlus[0])
            YMinNegScdDer.Fill(yMin[0]),   YMinMinNegScdDer.Fill(yMinMin[0])

          if scdDerOuter[1] > 0.0:
            EvtsPosScdDerXSOuter.append(iEvt+1)
            YPlusGausTestXSPosScdDer.Fill(yPlus[1] + yPlusPlus[1]/4)
            YPlusXSPosScdDer.Fill(yPlus[1]), YPlusPlusXSPosScdDer.Fill(yPlusPlus[1])
            YMinXSPosScdDer.Fill(yMin[1]),   YMinMinXSPosScdDer.Fill(yMinMin[1])
          else:
            YPlusGausTestXSNegScdDer.Fill(yPlus[1] + yPlusPlus[1]/4)
            YPlusXSNegScdDer.Fill(yPlus[1]), YPlusPlusXSNegScdDer.Fill(yPlusPlus[1])
            YMinXSNegScdDer.Fill(yMin[1]),   YMinMinXSNegScdDer.Fill(yMinMin[1])

          if scdDerOuter[2] > 0.0:
            EvtsPosScdDerAccOuter.append(iEvt+1)
            YPlusGausTestAccPosScdDer.Fill(yPlus[2] + yPlusPlus[2]/4)
            YPlusAccPosScdDer.Fill(yPlus[2]), YPlusPlusAccPosScdDer.Fill(yPlusPlus[2])
            YMinAccPosScdDer.Fill(yMin[2]),   YMinMinAccPosScdDer.Fill(yMinMin[2])
          else:
            YPlusGausTestAccNegScdDer.Fill(yPlus[2] + yPlusPlus[2]/4)
            YPlusAccNegScdDer.Fill(yPlus[2]), YPlusPlusAccNegScdDer.Fill(yPlusPlus[2])
            YMinAccNegScdDer.Fill(yMin[2]),   YMinMinAccNegScdDer.Fill(yMinMin[2])

#--- Save all the histograms containing information about all the events! ---#
Tfile.cd()
YPlusGausTest.Write(), YPlusGausTestXS.Write(), YPlusGausTestAcc.Write()
YPlus.Write(),         YPlusXS.Write(),         YPlusAcc.Write()
YPlusPlus.Write(),     YPlusPlusXS.Write(),     YPlusPlusAcc.Write()
YRelPlus.Write(),      YRelPlusXS.Write(),      YRelPlusAcc.Write()
YMin.Write(),          YMinXS.Write(),          YMinAcc.Write()
YMinMin.Write(),       YMinMinXS.Write(),       YMinMinAcc.Write()
YRelMin.Write(),       YRelMinXS.Write(),       YRelMinAcc.Write()
ScdDerInner.Write(),   ScdDerXSInner.Write(),   ScdDerAccInner.Write()
ScdDerOuter.Write(),   ScdDerXSOuter.Write(),   ScdDerAccOuter.Write()
ScdDerScatter.Write(), ScdDerXSScatter.Write(), ScdDerAccScatter.Write()

#---  Draw the likelihood distribution separately for events surviving and passing the cuts!  ---#
print "Nr of events with 2nd derivative > 0 (LnLog, LnLogXS & LnLogAcc -- using x = ",str(Var[xNeg[0]]),"/",str(Var[xMin]),"/",str(Var[xPos[0]]),") :",len(EvtsPosScdDerInner),", ",len(EvtsPosScdDerXSInner)," & ",len(EvtsPosScdDerAccInner)
print "Nr of events with 2nd derivative > 0 (LnLog, LnLogXS & LnLogAcc -- using x = ",str(Var[xNeg[1]]),"/",str(Var[xMin]),"/",str(Var[xPos[1]]),") :",len(EvtsPosScdDerOuter),", ",len(EvtsPosScdDerXSOuter)," & ",len(EvtsPosScdDerAccOuter)
print "Nr of events with Gaussiaanse vergelijking voor + (LnLog, LnLogXS & LnLogAcc) ", len(EvtsWithYPlusGausSmall),", ",len(EvtsWithYPlusGausSmallXS)," & ",len(EvtsWithYPlusGausSmallAcc)

LLPosScdDerInner, LLNegScdDerInner, LLXSPosScdDerInner, LLXSNegScdDerInner, LLAccPosScdDerInner, LLAccNegScdDerInner = [],[],[],[],[],[]
LLPosScdDerOuter, LLNegScdDerOuter, LLXSPosScdDerOuter, LLXSNegScdDerOuter, LLAccPosScdDerOuter, LLAccNegScdDerOuter = [],[],[],[],[],[]
LLPosScdDerBoth,  LLNegScdDerBoth,  LLXSPosScdDerBoth,  LLXSNegScdDerBoth,  LLAccPosScdDerBoth,  LLAccNegScdDerBoth  = [],[],[],[],[],[]
for ii in range(NrConfigs):
  LLPosScdDerInner.append(0), LLNegScdDerInner.append(0), LLXSPosScdDerInner.append(0), LLXSNegScdDerInner.append(0), LLAccPosScdDerInner.append(0), LLAccNegScdDerInner.append(0)
  LLPosScdDerOuter.append(0), LLNegScdDerOuter.append(0), LLXSPosScdDerOuter.append(0), LLXSNegScdDerOuter.append(0), LLAccPosScdDerOuter.append(0), LLAccNegScdDerOuter.append(0)
  LLPosScdDerBoth.append(0),  LLNegScdDerBoth.append(0),  LLXSPosScdDerBoth.append(0),  LLXSNegScdDerBoth.append(0),  LLAccPosScdDerBoth.append(0),  LLAccNegScdDerBoth.append(0)

EvtsPosScdDerBoth, EvtsPosScdDerXSBoth, EvtsPosScdDerAccBoth = 0, 0, 0
#Loop over all lines in weights file:
for LikelihoodLine in LikelihoodFile:
  LWord = LikelihoodLine.split()
  #Only interested in files starting with a number
  if str(LWord[0]) != "#" and str(LWord[3]) != "0.0" :
    #---  Separate the events with positive and negative second derivative (using both inner and outer Var points) ---#
    if int(LWord[0]) in EvtsPosScdDerInner and int(LWord[0]) in EvtsPosScdDerOuter:
      LLPosScdDerBoth[int(LWord[1])-1] = LLPosScdDerBoth[int(LWord[1])-1]-log(float(LWord[3]))
      if LWord[1] == "1": EvtsPosScdDerBoth += 1
    if int(LWord[0]) in EvtsPosScdDerXSInner and int(LWord[0]) in EvtsPosScdDerXSOuter:
      LLXSPosScdDerBoth[int(LWord[1])-1] = LLXSPosScdDerBoth[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])
      if LWord[1] == "1": EvtsPosScdDerXSBoth += 1
    if int(LWord[0]) in EvtsPosScdDerAccInner and int(LWord[0]) in EvtsPosScdDerAccOuter:
      LLAccPosScdDerBoth[int(LWord[1])-1] = LLAccPosScdDerBoth[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])+log(Acceptance[int(LWord[1])-1])
      if LWord[1] == "1": EvtsPosScdDerAccBoth += 1
    #---  Separate the events with positive and negative second derivative (using inner Var points) ---#
    if int(LWord[0]) in EvtsPosScdDerInner:    LLPosScdDerInner[int(LWord[1])-1] = LLPosScdDerInner[int(LWord[1])-1]-log(float(LWord[3]))
    else:                                      LLNegScdDerInner[int(LWord[1])-1] = LLNegScdDerInner[int(LWord[1])-1]-log(float(LWord[3]))
    if int(LWord[0]) in EvtsPosScdDerXSInner:  LLXSPosScdDerInner[int(LWord[1])-1] = LLXSPosScdDerInner[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])
    else:                                      LLXSNegScdDerInner[int(LWord[1])-1] = LLXSNegScdDerInner[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])
    if int(LWord[0]) in EvtsPosScdDerAccInner: LLAccPosScdDerInner[int(LWord[1])-1] = LLAccPosScdDerInner[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])+log(Acceptance[int(LWord[1])-1])
    else:                                      LLAccNegScdDerInner[int(LWord[1])-1] = LLAccNegScdDerInner[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])+log(Acceptance[int(LWord[1])-1])
    #---  Separate the events with positive and negative second derivative (using outer Var points) ---#
    if int(LWord[0]) in EvtsPosScdDerOuter:    LLPosScdDerOuter[int(LWord[1])-1] = LLPosScdDerOuter[int(LWord[1])-1]-log(float(LWord[3]))
    else:                                      LLNegScdDerOuter[int(LWord[1])-1] = LLNegScdDerOuter[int(LWord[1])-1]-log(float(LWord[3]))
    if int(LWord[0]) in EvtsPosScdDerXSOuter:  LLXSPosScdDerOuter[int(LWord[1])-1] = LLXSPosScdDerOuter[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])
    else:                                      LLXSNegScdDerOuter[int(LWord[1])-1] = LLXSNegScdDerOuter[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])
    if int(LWord[0]) in EvtsPosScdDerAccOuter: LLAccPosScdDerOuter[int(LWord[1])-1] = LLAccPosScdDerOuter[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])+log(Acceptance[int(LWord[1])-1])
    else:                                      LLAccNegScdDerOuter[int(LWord[1])-1] = LLAccNegScdDerOuter[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])+log(Acceptance[int(LWord[1])-1])

NrEvtsNegScdDerInner, NrEvtsNegScdDerXSInner, NrEvtsNegScdDerAccInner = int(nEvts)-len(EvtsPosScdDerInner), int(nEvts) - len(EvtsPosScdDerXSInner), int(nEvts) - len(EvtsPosScdDerAccInner)
NrEvtsNegScdDerOuter, NrEvtsNegScdDerXSOuter, NrEvtsNegScdDerAccOuter = int(nEvts)-len(EvtsPosScdDerOuter), int(nEvts) - len(EvtsPosScdDerXSOuter), int(nEvts) - len(EvtsPosScdDerAccOuter)
LLPosScdDerDistBoth     = TH1F('LLPosScdDerBoth',    '-ln(L) when both 2nd derivatives > 0 (no norm -- '+str(EvtsPosScdDerBoth)+'/'+nEvts+' evts -- '+title+')',    xBinZoom,xLowZoom,xHighZoom)
LLXSPosScdDerDistBoth   = TH1F('LLXSPosScdDerBoth',  '-ln(L) when both 2nd derivatives > 0 (XS norm -- '+str(EvtsPosScdDerXSBoth)+'/'+nEvts+' evts -- '+title+')',  xBinZoom,xLowZoom,xHighZoom)
LLAccPosScdDerDistBoth  = TH1F('LLAccPosScdDerBoth', '-ln(L) when both 2nd derivatives > 0 (Acc norm -- '+str(EvtsPosScdDerAccBoth)+'/'+nEvts+' evts -- '+title+')',xBinZoom,xLowZoom,xHighZoom)
LLPosScdDerDistInner    = TH1F('LLPosScdDerInner',   '-ln(L) when inner 2nd derivative > 0 (no norm -- '+str(len(EvtsPosScdDerInner))+'/'+nEvts+' evts -- '+title+')',    xBinZoom,xLowZoom,xHighZoom)
LLXSPosScdDerDistInner  = TH1F('LLXSPosScdDerInner', '-ln(L) when inner 2nd derivative > 0 (XS norm -- '+str(len(EvtsPosScdDerXSInner))+'/'+nEvts+' evts -- '+title+')',  xBinZoom,xLowZoom,xHighZoom)
LLAccPosScdDerDistInner = TH1F('LLAccPosScdDerInner','-ln(L) when inner 2nd derivative > 0 (Acc norm -- '+str(len(EvtsPosScdDerAccInner))+'/'+nEvts+' evts -- '+title+')',xBinZoom,xLowZoom,xHighZoom)
LLPosScdDerDistOuter    = TH1F('LLPosScdDerOuter',   '-ln(L) when outer 2nd derivative > 0 (no norm -- '+str(len(EvtsPosScdDerOuter))+'/'+nEvts+' evts -- '+title+')',    xBinZoom,xLowZoom,xHighZoom)
LLXSPosScdDerDistOuter  = TH1F('LLXSPosScdDerOuter', '-ln(L) when outer 2nd derivative > 0 (XS norm -- '+str(len(EvtsPosScdDerXSOuter))+'/'+nEvts+' evts -- '+title+')',  xBinZoom,xLowZoom,xHighZoom)
LLAccPosScdDerDistOuter = TH1F('LLAccPosScdDerOuter','-ln(L) when outer 2nd derivative > 0 (Acc norm -- '+str(len(EvtsPosScdDerAccOuter))+'/'+nEvts+' evts -- '+title+')',xBinZoom,xLowZoom,xHighZoom)
LLNegScdDerDistInner    = TH1F('LLNegScdDerInner',   '-ln(L) when inner 2nd derivative < 0 (no norm -- '+str(NrEvtsNegScdDerInner)+'/'+nEvts+' evts -- '+title+')',    xBinZoom,xLowZoom,xHighZoom)
LLXSNegScdDerDistInner  = TH1F('LLXSNegScdDerInner', '-ln(L) when inner 2nd derivative < 0 (XS norm -- '+str(NrEvtsNegScdDerXSInner)+'/'+nEvts+' evts -- '+title+')',  xBinZoom,xLowZoom,xHighZoom)
LLAccNegScdDerDistInner = TH1F('LLAccNegScdDerInner','-ln(L) when inner 2nd derivative < 0 (Acc norm -- '+str(NrEvtsNegScdDerAccInner)+'/'+nEvts+' evts -- '+title+')',xBinZoom,xLowZoom,xHighZoom)
LLNegScdDerDistOuter    = TH1F('LLNegScdDerOuter',   '-ln(L) when outer 2nd derivative < 0 (no norm -- '+str(NrEvtsNegScdDerOuter)+'/'+nEvts+' evts -- '+title+')',    xBinZoom,xLowZoom,xHighZoom)
LLXSNegScdDerDistOuter  = TH1F('LLXSNegScdDerOuter', '-ln(L) when outer 2nd derivative < 0 (XS norm -- '+str(NrEvtsNegScdDerXSOuter)+'/'+nEvts+' evts -- '+title+')',  xBinZoom,xLowZoom,xHighZoom)
LLAccNegScdDerDistOuter = TH1F('LLAccNegScdDerOuter','-ln(L) when outer 2nd derivative < 0 (Acc norm -- '+str(NrEvtsNegScdDerAccOuter)+'/'+nEvts+' evts -- '+title+')',xBinZoom,xLowZoom,xHighZoom)

for ii in range(NrConfigs):
  LLPosScdDerDistInner.SetBinContent(LLPosScdDerDistInner.FindBin(Var[ii]), float(LLPosScdDerInner[ii]))
  LLXSPosScdDerDistInner.SetBinContent(LLXSPosScdDerDistInner.FindBin(Var[ii]),float(LLXSPosScdDerInner[ii]))
  LLAccPosScdDerDistInner.SetBinContent(LLAccPosScdDerDistInner.FindBin(Var[ii]),float(LLAccPosScdDerInner[ii]))
  LLNegScdDerDistInner.SetBinContent(LLNegScdDerDistInner.FindBin(Var[ii]), float(LLNegScdDerInner[ii]))
  LLXSNegScdDerDistInner.SetBinContent(LLXSNegScdDerDistInner.FindBin(Var[ii]),float(LLXSNegScdDerInner[ii]))
  LLAccNegScdDerDistInner.SetBinContent(LLAccNegScdDerDistInner.FindBin(Var[ii]),float(LLAccNegScdDerInner[ii]))

  LLPosScdDerDistOuter.SetBinContent(LLPosScdDerDistOuter.FindBin(Var[ii]), float(LLPosScdDerOuter[ii]))
  LLXSPosScdDerDistOuter.SetBinContent(LLXSPosScdDerDistOuter.FindBin(Var[ii]),float(LLXSPosScdDerOuter[ii]))
  LLAccPosScdDerDistOuter.SetBinContent(LLAccPosScdDerDistOuter.FindBin(Var[ii]),float(LLAccPosScdDerOuter[ii]))
  LLNegScdDerDistOuter.SetBinContent(LLNegScdDerDistOuter.FindBin(Var[ii]), float(LLNegScdDerOuter[ii]))
  LLXSNegScdDerDistOuter.SetBinContent(LLXSNegScdDerDistOuter.FindBin(Var[ii]),float(LLXSNegScdDerOuter[ii]))
  LLAccNegScdDerDistOuter.SetBinContent(LLAccNegScdDerDistOuter.FindBin(Var[ii]),float(LLAccNegScdDerOuter[ii]))

  LLPosScdDerDistBoth.SetBinContent(LLPosScdDerDistBoth.FindBin(Var[ii]), float(LLPosScdDerBoth[ii]))
  LLXSPosScdDerDistBoth.SetBinContent(LLXSPosScdDerDistBoth.FindBin(Var[ii]),float(LLXSPosScdDerBoth[ii]))
  LLAccPosScdDerDistBoth.SetBinContent(LLAccPosScdDerDistBoth.FindBin(Var[ii]),float(LLAccPosScdDerBoth[ii]))

AppliedCutsDir = Tfile.mkdir("LikelihoodAfterCuts")
SignScdDerDir = AppliedCutsDir.mkdir("SignSecondDerivative")
SignScdDerDir.cd()
LLPosScdDerDistInner.Write(),   LLXSPosScdDerDistInner.Write(),   LLAccPosScdDerDistInner.Write()
LLNegScdDerDistInner.Write(),   LLXSNegScdDerDistInner.Write(),   LLAccNegScdDerDistInner.Write()
LLPosScdDerDistOuter.Write(),   LLXSPosScdDerDistOuter.Write(),   LLAccPosScdDerDistOuter.Write()
LLNegScdDerDistOuter.Write(),   LLXSNegScdDerDistOuter.Write(),   LLAccNegScdDerDistOuter.Write()
LLPosScdDerDistBoth.Write(), LLXSPosScdDerDistBoth.Write(), LLAccPosScdDerDistBoth.Write()
#-- Save the variables separate of scd Der sign  --#
#--   --> Can give hint for good cut!            --#
#- Save the YPlusGausTest variables -#
YPlusGausTestPosScdDer.SetLineColor(3), YPlusGausTestXSPosScdDer.SetLineColor(3), YPlusGausTestAccPosScdDer.SetLineColor(3)
YPlusGausTestNegScdDer.SetLineColor(2), YPlusGausTestXSNegScdDer.SetLineColor(2), YPlusGausTestAccNegScdDer.SetLineColor(2)
YPlusGausTestCanvas = TCanvas('YPlusGausTestCanvas','Combined distribution of YPlusGausTest variable for events with positive and negative scd derivative (no normalisation -- '+title+')')
YPlusGausTestCanvas.cd(), YPlusGausTestPosScdDer.Draw(), YPlusGausTestNegScdDer.Draw("same"), YPlusGausTestCanvas.Write()
YPlusGausTestXSCanvas = TCanvas('YPlusGausTestXSCanvas','Combined distribution of YPlusGausTest variable for events with positive and negative scd derivative (XS normalisation -- '+title+')')
YPlusGausTestXSCanvas.cd(), YPlusGausTestXSNegScdDer.Draw(), YPlusGausTestXSPosScdDer.Draw("same"), YPlusGausTestXSCanvas.Write()
YPlusGausTestAccCanvas = TCanvas('YPlusGausTestAccCanvas','Combined distribution of YPlusGausTest variable for events with positive and negative scd derivative (Acc normalisation -- '+title+')')
YPlusGausTestAccCanvas.cd(), YPlusGausTestAccNegScdDer.Draw(), YPlusGausTestAccPosScdDer.Draw("same"), YPlusGausTestAccCanvas.Write()
#- Save the YPlus variables -#
YPlusPosScdDer.SetLineColor(3), YPlusXSPosScdDer.SetLineColor(3), YPlusAccPosScdDer.SetLineColor(3)
YPlusNegScdDer.SetLineColor(2), YPlusXSNegScdDer.SetLineColor(2), YPlusAccNegScdDer.SetLineColor(2)
YPlusCanvas = TCanvas('YPlusCanvas','Combined distribution of YPlus variable for events with positive and negative scd derivative (no normalisation -- '+title+')')
YPlusCanvas.cd(), YPlusPosScdDer.Draw(), YPlusNegScdDer.Draw("same"), YPlusCanvas.Write()
YPlusXSCanvas = TCanvas('YPlusXSCanvas','Combined distribution of YPlus variable for events with positive and negative scd derivative (XS normalisation -- '+title+')')
YPlusXSCanvas.cd(), YPlusXSNegScdDer.Draw(), YPlusXSPosScdDer.Draw("same"), YPlusXSCanvas.Write()
YPlusAccCanvas = TCanvas('YPlusAccCanvas','Combined distribution of YPlus variable for events with positive and negative scd derivative (Acc normalisation -- '+title+')')
YPlusAccCanvas.cd(), YPlusAccNegScdDer.Draw(), YPlusAccPosScdDer.Draw("same"), YPlusAccCanvas.Write()
#- Save the YPlusPlus variables -#
YPlusPlusPosScdDer.SetLineColor(3), YPlusPlusXSPosScdDer.SetLineColor(3), YPlusPlusAccPosScdDer.SetLineColor(3)
YPlusPlusNegScdDer.SetLineColor(2), YPlusPlusXSNegScdDer.SetLineColor(2), YPlusPlusAccNegScdDer.SetLineColor(2)
YPlusPlusCanvas = TCanvas('YPlusPlusCanvas','Combined distribution of YPlusPlus variable for events with positive and negative scd derivative (no normalisation -- '+title+')')
YPlusPlusCanvas.cd(), YPlusPlusPosScdDer.Draw(), YPlusPlusNegScdDer.Draw("same"), YPlusPlusCanvas.Write()
YPlusPlusXSCanvas = TCanvas('YPlusPlusXSCanvas','Combined distribution of YPlusPlus variable for events with positive and negative scd derivative (XS normalisation -- '+title+')')
YPlusPlusXSCanvas.cd(), YPlusPlusXSNegScdDer.Draw(), YPlusPlusXSPosScdDer.Draw("same"), YPlusPlusXSCanvas.Write()
YPlusPlusAccCanvas = TCanvas('YPlusPlusAccCanvas','Combined distribution of YPlusPlus variable for events with positive and negative scd derivative (Acc normalisation -- '+title+')')
YPlusPlusAccCanvas.cd(), YPlusPlusAccNegScdDer.Draw(), YPlusPlusAccPosScdDer.Draw("same"), YPlusPlusAccCanvas.Write()
#- Save the YMin variables -#
YMinPosScdDer.SetLineColor(3), YMinXSPosScdDer.SetLineColor(3), YMinAccPosScdDer.SetLineColor(3)
YMinNegScdDer.SetLineColor(2), YMinXSNegScdDer.SetLineColor(2), YMinAccNegScdDer.SetLineColor(2)
YMinCanvas = TCanvas('YMinCanvas','Combined distribution of YMin variable for events with positive and negative scd derivative (no normalisation -- '+title+')')
YMinCanvas.cd(), YMinPosScdDer.Draw(), YMinNegScdDer.Draw("same"), YMinCanvas.Write()
YMinXSCanvas = TCanvas('YMinXSCanvas','Combined distribution of YMin variable for events with positive and negative scd derivative (XS normalisation -- '+title+')')
YMinXSCanvas.cd(), YMinXSNegScdDer.Draw(), YMinXSPosScdDer.Draw("same"), YMinXSCanvas.Write()
YMinAccCanvas = TCanvas('YMinAccCanvas','Combined distribution of YMin variable for events with positive and negative scd derivative (Acc normalisation -- '+title+')')
YMinAccCanvas.cd(), YMinAccNegScdDer.Draw(), YMinAccPosScdDer.Draw("same"), YMinAccCanvas.Write()
#- Save the YMinMin variables -#
YMinMinPosScdDer.SetLineColor(3), YMinMinXSPosScdDer.SetLineColor(3), YMinMinAccPosScdDer.SetLineColor(3)
YMinMinNegScdDer.SetLineColor(2), YMinMinXSNegScdDer.SetLineColor(2), YMinMinAccNegScdDer.SetLineColor(2)
YMinMinCanvas = TCanvas('YMinMinCanvas','Combined distribution of YMinMin variable for events with positive and negative scd derivative (no normalisation -- '+title+')')
YMinMinCanvas.cd(), YMinMinPosScdDer.Draw(), YMinMinNegScdDer.Draw("same"), YMinMinCanvas.Write()
YMinMinXSCanvas = TCanvas('YMinMinXSCanvas','Combined distribution of YMinMin variable for events with positive and negative scd derivative (XS normalisation -- '+title+')')
YMinMinXSCanvas.cd(), YMinMinXSNegScdDer.Draw(), YMinMinXSPosScdDer.Draw("same"), YMinMinXSCanvas.Write()
YMinMinAccCanvas = TCanvas('YMinMinAccCanvas','Combined distribution of YMinMin variable for events with positive and negative scd derivative (Acc normalisation -- '+title+')')
YMinMinAccCanvas.cd(), YMinMinAccNegScdDer.Draw(), YMinMinAccPosScdDer.Draw("same"), YMinMinAccCanvas.Write()

