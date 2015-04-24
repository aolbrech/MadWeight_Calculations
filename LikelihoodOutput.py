#! python
import os
import sys
from re import search           #used for search!
from array import array
from math import log
from ROOT import TH1F,TFile,TCanvas,TLegend,gStyle

#Get all the input from the command line:
if len(sys.argv) <= 1:
  print "Need to give the directory of interest, the considered kinematic variable and whether scdDer cuts should be applied in the command line !"
  print " Correct syntax is : python FitLikelihood.py Events/blabla/ MTop (WithCut)"
  sys.exit()
elif len(sys.argv) == 2:
  print "Need to specify the considered kinematic variable (MTop or RVR)"
  print " Correct syntax is : python FitLikelihood.py Events/blabla/ MTop (WithCut)"
  sys.exit()
whichDir = sys.argv[1]
KinVariable = sys.argv[2]

if KinVariable != "MTop" and KinVariable != "RVR":
  print "Need to specify which kinematic variable should be considered (MTop or RVR are the only options!!)"
  KinVariable = raw_input('--> Choose one of the two : ')
RVRValues =  ["Re(VR) = -0.25","Re(VR) = -0.2","Re(VR) = -0.15","Re(VR) = -0.1","Re(VR) = -0.05","Re(VR) = 0.0","Re(VR) = 0.05","Re(VR) = 0.1","Re(VR) = 0.15","Re(VR) = 0.2","Re(VR) = 0.25"]
RVR =        [-0.25,           -0.2,           -0.15,           -0.1,           -0.05,            0.0,          0.05,           0.1,           0.15,           0.2,           0.25           ]
MGXS =       [12.66,           12.06555,       11.61,           11.25909,       11.03,            10.90059,     10.88,          10.97767,      11.17,          11.49883,      11.91          ]
Acceptance = [0.21843,         0.21672,        0.21367,         0.21737,        0.21614,          0.21670,      0.21531,        0.21677,       0.21408,        0.21437,       0.21778        ]
#MGXSe =      [X                0.0093464076837, X               0.00836607833038, X               0.00822214433527, X           0.00847293509122, X            0.00901976602967, X           ]

if KinVariable == "RVR":
  #Information about the scanned RVR values and the corresponding cross-section
  VarValues =  ["Re(V_{R}) = -0.5","Re(V_{R}) = -0.3","Re(VR) = -0.25","Re(V_{R}) = -0.2","Re(VR) = -0.15","Re(V_{R}) = -0.1","Re(V_{R}) = -0.05","Re(V_{R}) = 0.0","Re(V_{R}) = 0.05","Re(V_{R}) = 0.1","Re(VR) = 0.15","Re(V_{R}) = 0.2","Re(VR) = 0.25","Re(V_{R}) = 0.3", "Re(V_{R}) = 0.5"]
  Var =        array('d',[-0.5,    -0.3,    -0.25,    -0.2,     -0.15,    -0.1,     -0.05,    0.0,      0.05,     0.1,      0.15,     0.2,      0.25,     0.3,       0.5    ])
  MGXS =       array('d',[17.9275, 13.3944, 12.66011, 12.06555, 11.60956, 11.25909, 11.02784, 10.90059, 10.88228, 10.97767, 11.17366, 11.49883, 11.90668, 12.49056,  16.1508])
  MGXSe =      array('d',[0.01231, 0.00995, 0.008608, 0.009347, 0.009240, 0.008366, 0.0,      0.008222, 0.0,      0.008472, 0.008225, 0.009017, 0.008763, 0.0087467, 0.01130])
  Acceptance = array('d',[0.22164, 0.21742, 0.21843,  0.21672,  0.21367,  0.21737,  0.21614,  0.21670,  0.21531,  0.21677,  0.21408,  0.21437,  0.21778,  0.21793,   0.22205])

  #Select which window of RVR values was considered!
  VarWindow = raw_input('** Choose the correct RVR-window corresponding to the studied file ** \n** Different options are : \n  1) Wide   : [-0.5, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.5] \n  2) Narrow : [-0.3, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2, 0.3] \n  3) Zoomed : [-0.25, -0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.25] \n --> Choose the correct number : ')
  KinVar = "Re(V_{R})"

  if VarWindow == "1":
    VarValues.pop(12), Var.pop(12), MGXS.pop(12), MGXSe.pop(12), Acceptance.pop(12)
    VarValues.pop(10), Var.pop(10), MGXS.pop(10), MGXSe.pop(10), Acceptance.pop(10)
    VarValues.pop(8),  Var.pop(8),  MGXS.pop(8),  MGXSe.pop(8),  Acceptance.pop(8)
    VarValues.pop(6),  Var.pop(6),  MGXS.pop(6),  MGXSe.pop(6),  Acceptance.pop(6)
    VarValues.pop(4),  Var.pop(4),  MGXS.pop(4),  MGXSe.pop(4),  Acceptance.pop(4)
    VarValues.pop(2),  Var.pop(2),  MGXS.pop(2),  MGXSe.pop(2),  Acceptance.pop(2)
    xBin, xLow, xHigh = 11, -0.55, 0.55
  elif VarWindow == "2":
    VarValues.pop(14), Var.pop(14), MGXS.pop(14), MGXSe.pop(14), Acceptance.pop(14)
    VarValues.pop(12), Var.pop(12), MGXS.pop(12), MGXSe.pop(12), Acceptance.pop(12)
    VarValues.pop(10), Var.pop(10), MGXS.pop(10), MGXSe.pop(10), Acceptance.pop(10)
    VarValues.pop(4),  Var.pop(4),  MGXS.pop(4),  MGXSe.pop(4),  Acceptance.pop(4)
    VarValues.pop(2),  Var.pop(2),  MGXS.pop(2),  MGXSe.pop(2),  Acceptance.pop(2)
    VarValues.pop(0),  Var.pop(0),  MGXS.pop(0),  MGXSe.pop(0),  Acceptance.pop(0)
    xBin, xLow, xHigh = 13, -0.325, 0.325
  elif VarWindow == "3":
    VarValues.pop(14), Var.pop(14), MGXS.pop(14), MGXSe.pop(14), Acceptance.pop(14)
    VarValues.pop(13), Var.pop(13), MGXS.pop(13), MGXSe.pop(13), Acceptance.pop(13)
    VarValues.pop(1),  Var.pop(1),  MGXS.pop(1),  MGXSe.pop(1),  Acceptance.pop(1)
    VarValues.pop(0),  Var.pop(0),  MGXS.pop(0),  MGXSe.pop(0),  Acceptance.pop(0)
    xBin, xLow, xHigh = 11, -0.275, 0.275

elif KinVariable == "MTop":
  #Information about the scanned MTop values and the corresponding cross-section
  VarValues = ["m_{top} = 153",  "m_{top} = 163",  "m_{top} = 170",  "m_{top} = 171",  "m_{top} = 172",  "m_{top} = 173",  "m_{top} = 174",  "m_{top} = 175","m_{top} = 183",  "m_{top} = 193"]
  Var =        array('d',[153,              163,              170,              171,              172,              173,              174,              175,            183,              193            ])
  MGXS =       array('d',[8.20916,          9.6299,           10.57123,         10.70485,         10.8257,          10.96469,         11.08428,         11.22448,       12.18068,         13.3046        ])
  MGXSe =      array('d',[0.00641266950107, 0.00775899974932, 0.00857143063963, 0.00814303595657, 0.00878899028501, 0.00816801717126, 0.00904797742371, 0.008653800078, 0.00931290317946, 0.0103310001752])
  Acceptance = array('d',[0.16203,          0.19152,          0.21008,          0.21460,          0.21735,          0.21290,          0.21752,          0.22185,        0.23941,          0.26413        ])

  #Select which window of masses was considered!
  VarWindow = raw_input('** Choose the correct mass-window corresponding to the studied file ** \n** Different options are : \n  1) Wide   : [153, 163, 170, 171, 172, 173, 174, 175, 183, 193] \n  2) Normal : [153, 163, 171, 172, 173, 174, 175, 183, 193] \n  3) Narrow : [171, 172, 173, 174, 175] \n --> Choose the correct number : ')

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

print " List of considered Var values is : ",Var    

#File of interest:
list_dir = []
list_dir = os.listdir(whichDir)
WeightsFileArray = []
weightsFileCounter = 0
for file in list_dir:
  if file.endswith(".out") and file.startswith("weights"): # eg: '.txt'
    weightsFileCounter += 1
    WeightsFileArray.append(file)

if int(weightsFileCounter) == 1:
  WeightsFile = open(os.path.join(whichDir+''+WeightsFileArray[0]),'r')
elif int(weightsFileCounter) == 0:
  print "No weights file found in this directory !"
  sys.exit()
elif int(weightsFileCounter) > 1:
  for ii in range(len(WeightsFileArray)):
    print " ",ii," ) ",WeightsFileArray[ii]
  fileNr = raw_input('Choose the number of the file of interest! : ')
  WeightsFile = open(os.path.join(whichDir+''+WeightsFileArray[int(fileNr)]),'r')
print "Will be using file : ",WeightsFile

#Set title of root file!
title = ""
if whichDir.find("Correct") <= len(whichDir)     and whichDir.find("Correct") > 0:   title = "Correct"
elif whichDir.find("Wrong") <= len(whichDir)     and whichDir.find("Wrong") > 0:     title = "Wrong"
elif whichDir.find("Unmatched") <= len(whichDir) and whichDir.find("Unmatched") > 0: title = "Unmatched"

GenLevel = False
if whichDir.find("Reco") <= len(whichDir)  and whichDir.find("Reco") > 0: title += "Reco"
elif whichDir.find("Gen") <= len(whichDir) and whichDir.find("Gen") > 0:
  title += "Gen"
  GenLevel = True
title = title+"_"+KinVariable

#ROOT file where all the information will be stored:
Tfile = TFile(os.path.join(whichDir+"Likelihood_"+title+".root"),'recreate')
LLDist   = TH1F('LL',    '-ln(L) for Single Gaus TF (no normalisation -- '+title+' evts)',xBin,xLow,xHigh)
LLXSDist = TH1F('LL_XS', '-ln(L) for Single Gaus TF (XS normalisation -- '+title+' evts)',xBin,xLow,xHigh)
if GenLevel == False: LLAccDist = TH1F('LL_Acc','-ln(L) for Single Gaus TF (Acc normalisation -- '+title+' evts)',xBin,xLow,xHigh)
gStyle.SetOptStat(0)

#Set the bin content originally very low such that these points don't show up in the histograms
#--> Should still try whether TGraph could be used?
for ii in range(LLDist.GetNbinsX()):
  LLDist.SetBinContent(ii+1,-1000)
  LLXSDist.SetBinContent(ii+1,-1000)
  if GenLevel == False: LLAccDist.SetBinContent(ii+1,-1000)  

#Initialize the LL and LLXS which will be created from the weights to 0
LL, LLError = [], []
LLXS, LLXSError = [], []
LLAcc, LLAccError = [], []
NrConsEvents = []
for ii in range(len(Var)):
  LL.append(0), LLXS.append(0) 
  LLError.append(0), LLXSError.append(0)
  if GenLevel == False: LLAcc.append(0), LLAccError.append(0)
  NrConsEvents.append(0)

#Loop over all lines in weights file:
for WeightLine in WeightsFile:
  WeightWord = WeightLine.split()
  #Only interested in files starting with a number
  if str(WeightWord[0]) != "#" and str(WeightWord[3]) != "0.0" :
    LL[int(WeightWord[1])-1] = LL[int(WeightWord[1])-1] -log(float(WeightWord[3]))
    LLXS[int(WeightWord[1])-1] = LLXS[int(WeightWord[1])-1] - log(float(WeightWord[3])) + log(MGXS[int(WeightWord[1])-1])
    if GenLevel == False: LLAcc[int(WeightWord[1])-1] = LLAcc[int(WeightWord[1])-1] - log(float(WeightWord[3])) + log(MGXS[int(WeightWord[1])-1]) + log(Acceptance[int(WeightWord[1])-1])
    NrConsEvents[int(WeightWord[1])-1] +=1 # NrConsEvents[int(WeightWord[1])-1] +1

print "Obtained LL output = ",LL
print "Obtained LLXS output =",LLXS
if GenLevel == False: print "Obtained LLAcc output =",LLAcc
print "Nr of events considered :",NrConsEvents

#Fill the histograms and calculate the minimums in order to draw the Delta ChiSq
LLMin = float(LL[0])
LLXSMin = float(LLXS[0])
if GenLevel == False: LLAccMin = float(LLAcc[0])
for ii in range(len(Var)):
  LLDist.SetBinContent(LLDist.FindBin(Var[ii]), float(LL[ii]))
  LLXSDist.SetBinContent(LLXSDist.FindBin(Var[ii]),float(LLXS[ii]))
  if GenLevel == False: LLAccDist.SetBinContent(LLAccDist.FindBin(Var[ii]),float(LLAcc[ii]))
  if float(LL[ii]) < LLMin:
    LLMin = float(LL[ii])
  if float(LLXS[ii]) < LLXSMin:
    LLXSMin = float(LLXS[ii])  
  if float(LLAcc[ii]) < LLAccMin and GenLevel == False:
    LLAccMin = float(LLAcc[ii])  

LLDist.SetMinimum(0),   LLDist.SetMarkerStyle(20),   LLDist.SetLineColor(1),   LLDist.SetMarkerColor(1),   LLDist.Write()
LLXSDist.SetMinimum(0), LLXSDist.SetMarkerStyle(21), LLXSDist.SetLineColor(2), LLXSDist.SetMarkerColor(2), LLXSDist.Write()
if GenLevel == False: LLAccDist.SetMinimum(0), LLAccDist.SetMarkerStyle(22), LLAccDist.SetLineColor(3), LLAccDist.SetMarkerColor(3), LLAccDist.Write()
print " Found minimums (LL, LLXS & LLAcc): ",LLMin,", ",LLXSMin," and ",LLAccMin
if GenLevel == False: print " Found minimums (LL & LLXS): ",LLMin," & ",LLXSMin

#Change bin content for combined canvas!
for ii in range(len(Var)):
  LLDist.SetBinContent(LLDist.FindBin(Var[ii]), 2*(float(LLDist.GetBinContent(LLDist.FindBin(Var[ii]))) - LLMin))
  LLXSDist.SetBinContent(LLXSDist.FindBin(Var[ii]), 2*(float(LLXSDist.GetBinContent(LLXSDist.FindBin(Var[ii]))) - LLXSMin))
  if GenLevel == False: LLAccDist.SetBinContent(LLAccDist.FindBin(Var[ii]), 2*(float(LLAccDist.GetBinContent(LLAccDist.FindBin(Var[ii]))) - LLAccMin))

#Prepare everything for the Delta ChiSq canvas!
LLDist.GetYaxis().SetTitle("#Delta #chi^{2} = -2ln(L/L_{max})")
LLDist.GetYaxis().SetTitleOffset(1.3)

LLCanvas = TCanvas('LLDistributions','Comparing influence of XS and acceptance normalisation on -ln(L) -- '+title+' evts')
legend = TLegend(0.75,0.7,0.95,0.95)
LLCanvas.cd()
legend.AddEntry(LLDist,'Non-normalized -ln(likelihood)',"p")
legend.AddEntry(LLXSDist,'Cross-section normalized -ln(likelihood)',"p")
if GenLevel == False: legend.AddEntry(LLAccDist,'Cross-section and acceptance normalized -ln(likelihood)',"p")
LLDist.SetTitle(LLCanvas.GetTitle())
LLDist.Draw("P")
LLXSDist.Draw("SAMEP")
if GenLevel == False: LLAccDist.Draw("SAMEP")
legend.Draw()
LLCanvas.Write()                 
