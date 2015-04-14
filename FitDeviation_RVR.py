#! python
import os
import sys
from math import log
#import ROOT
from array import array
from ROOT import TH1F,TH2F,TFile,TCanvas,TLegend,gStyle,TDirectory,gROOT,TGraph
gROOT.SetBatch(True)                       #Don't print the histograms (necessary for fit)

#Directory of interest:
if len(sys.argv) <= 1:
  print "Need to give the directory of interest in command line !"
  sys.exit()
whichDir = sys.argv[1] 
print " Interested in directory : ",whichDir

#Information about the scanned RVR values and the corresponding cross-section
RVRValues =  ["Re(V_{R}) = -0.5","Re(V_{R}) = -0.3","Re(V_{R}) = -0.2","Re(V_{R}) = -0.1", "Re(V_{R}) = 0.0","Re(V_{R}) = 0.1","Re(V_{R}) = 0.2","Re(V_{R}) = 0.3", "Re(V_{R}) = 0.5"]
RVR =        array('d',[-0.5,              -0.3,              -0.2,              -0.1,               0.0,              0.1,              0.2,              0.3,               0.5              ])
MGXS =       array('d',[17.9275,           13.3944,           12.06555,          11.25909,           10.90059,         10.97767,         11.49883,         12.49056,          16.1508          ])
MGXSe =      array('d',[0.0123100357311,   0.00995808028337,  0.0093464076837,   0.00836607833038,   0.00822214433527, 0.00847293509122, 0.00901976602967, 0.00874682264197,  0.0113081652137  ])
Acceptance = array('d',[0.22164,           0.21742,           0.21672,           0.21737,            0.21670,          0.21677,          0.21437,          0.21793,           0.22205          ])

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
  LikelihoodFile = open(os.path.join(whichDir+''+WeightsFileArray[int(fileNr)]),'r')

print "Will be using file : ",WeightsFile

#ROOT file where all the information will be stored:
Tfile = TFile(os.path.join(whichDir+"FitDeviation.root"),'recreate')
gStyle.SetOptStat(0)
YPlusGausTest    = TH1F('YPlusGausTest',   'Comparison of fit deviation in Re(VR) = 0.2 and Re(VR) = 0.1 (no normalisation)', 250,-0.5,0.5)
YPlusGausTestXS  = TH1F('YPlusGausTestXS', 'Comparison of fit deviation in Re(VR) = 0.2 and Re(VR) = 0.1 (XS normalisation)', 250,-0.5,0.5)
YPlusGausTestAcc = TH1F('YPlusGausTestAcc','Comparison of fit deviation in Re(VR) = 0.2 and Re(VR) = 0.1 (Acc normalisation)',250,-0.5,0.5)
YPlusGausTestPosScdDer    = TH1F('YPlusGausTestPosScdDer',   'Comparison of fit deviation in Re(VR) = 0.2 and Re(VR) = 0.1 (no normalisation -- second der > 0)', 250,-0.2,0.2)
YPlusGausTestXSPosScdDer  = TH1F('YPlusGausTestXSPosScdDer', 'Comparison of fit deviation in Re(VR) = 0.2 and Re(VR) = 0.1 (XS normalisation -- second der > 0)', 250,-0.2,0.2)
YPlusGausTestAccPosScdDer = TH1F('YPlusGausTestAccPosScdDer','Comparison of fit deviation in Re(VR) = 0.2 and Re(VR) = 0.1 (Acc normalisation -- second der > 0)',250,-0.2,0.2)
YPlusGausTestNegScdDer    = TH1F('YPlusGausTestNegScdDer',   'Comparison of fit deviation in Re(VR) = 0.2 and Re(VR) = 0.1 (no normalisation -- second der < 0)', 250,-0.2,0.2)
YPlusGausTestXSNegScdDer  = TH1F('YPlusGausTestXSNegScdDer', 'Comparison of fit deviation in Re(VR) = 0.2 and Re(VR) = 0.1 (XS normalisation -- second der < 0)', 250,-0.2,0.2)
YPlusGausTestAccNegScdDer = TH1F('YPlusGausTestAccNegScdDer','Comparison of fit deviation in Re(VR) = 0.2 and Re(VR) = 0.1 (Acc normalisation -- second der < 0)',250,-0.2,0.2)

YPlus    = TH1F('YPlus',   'Deviation from parabolic fit for Re(VR) = 0.1 (no normalisation)' ,150,-0.5,0.5)
YPlusXS  = TH1F('YPlusXS', 'Deviation from parabolic fit for Re(VR) = 0.1 (XS normalisation)' ,150,-0.5,0.5)
YPlusAcc = TH1F('YPlusAcc','Deviation from parabolic fit for Re(VR) = 0.1 (Acc normalisation)',150,-0.5,0.5)
YPlusPosScdDer    = TH1F('YPlusPosScdDer',   'Deviation from parabolic fit for Re(VR) = 0.1 (no normalisation -- second der > 0)' ,150,-0.2,0.2)
YPlusXSPosScdDer  = TH1F('YPlusXSPosScdDer', 'Deviation from parabolic fit for Re(VR) = 0.1 (XS normalisation -- second der > 0)' ,150,-0.2,0.2)
YPlusAccPosScdDer = TH1F('YPlusAccPosScdDer','Deviation from parabolic fit for Re(VR) = 0.1 (Acc normalisation -- second der > 0)',150,-0.2,0.2)
YPlusNegScdDer    = TH1F('YPlusNegScdDer',   'Deviation from parabolic fit for Re(VR) = 0.1 (no normalisation -- second der < 0)' ,150,-0.2,0.2)
YPlusXSNegScdDer  = TH1F('YPlusXSNegScdDer', 'Deviation from parabolic fit for Re(VR) = 0.1 (XS normalisation -- second der < 0)' ,150,-0.2,0.2)
YPlusAccNegScdDer = TH1F('YPlusAccNegScdDer','Deviation from parabolic fit for Re(VR) = 0.1 (Acc normalisation -- second der < 0)',150,-0.2,0.2)
YPlusPosScdDerPosYPlusPlus    = TH1F('YPlusPosScdDerPosYPlusPlus',   'Deviation from parabolic fit for Re(VR) = 0.1 (no normalisation -- second der > 0 & y++ > 0)' ,150,-0.2,0.2)
YPlusXSPosScdDerPosYPlusPlus  = TH1F('YPlusXSPosScdDerPosYPlusPlus', 'Deviation from parabolic fit for Re(VR) = 0.1 (XS normalisation -- second der > 0 & y++ > 0)' ,150,-0.2,0.2)
YPlusAccPosScdDerPosYPlusPlus = TH1F('YPlusAccPosScdDerPosYPlusPlus','Deviation from parabolic fit for Re(VR) = 0.1 (Acc normalisation -- second der > 0 & y++ > 0)',150,-0.2,0.2)

YPlusPlus    = TH1F('YPlusPlus',   'Deviation from parabolic fit for Re(VR) = 0.2 (no normalisation)' ,150,-0.5,0.5)
YPlusPlusXS  = TH1F('YPlusPlusXS', 'Deviation from parabolic fit for Re(VR) = 0.2 (XS normalisation)' ,150,-0.5,0.5)
YPlusPlusAcc = TH1F('YPlusPlusAcc','Deviation from parabolic fit for Re(VR) = 0.2 (Acc normalisation)',150,-0.5,0.5)
YPlusPlusPosScdDer    = TH1F('YPlusPlusPosScdDer',   'Deviation from parabolic fit for Re(VR) = 0.2 (no normalisation -- second der > 0)' ,150,-0.2,0.2)
YPlusPlusXSPosScdDer  = TH1F('YPlusPlusXSPosScdDer', 'Deviation from parabolic fit for Re(VR) = 0.2 (XS normalisation -- second der > 0)' ,150,-0.2,0.2)
YPlusPlusAccPosScdDer = TH1F('YPlusPlusAccPosScdDer','Deviation from parabolic fit for Re(VR) = 0.2 (Acc normalisation -- second der > 0)',150,-0.2,0.2)
YPlusPlusNegScdDer    = TH1F('YPlusPlusNegScdDer',   'Deviation from parabolic fit for Re(VR) = 0.2 (no normalisation -- second der < 0)' ,150,-0.2,0.2)
YPlusPlusXSNegScdDer  = TH1F('YPlusPlusXSNegScdDer', 'Deviation from parabolic fit for Re(VR) = 0.2 (XS normalisation -- second der < 0)' ,150,-0.2,0.2)
YPlusPlusAccNegScdDer = TH1F('YPlusPlusAccNegScdDer','Deviation from parabolic fit for Re(VR) = 0.2 (Acc normalisation -- second der < 0)',150,-0.2,0.2)

YMin    = TH1F('YMin',   'Deviation from parabolic fit for Re(VR) = -0.1 (no normalisation)' ,150,-0.5,0.5)
YMinXS  = TH1F('YMinXS', 'Deviation from parabolic fit for Re(VR) = -0.1 (XS normalisation)' ,150,-0.5,0.5)
YMinAcc = TH1F('YMinAcc','Deviation from parabolic fit for Re(VR) = -0.1 (Acc normalisation)',150,-0.5,0.5)
YMinPosScdDer    = TH1F('YMinPosScdDer',   'Deviation from parabolic fit for Re(VR) = -0.1 (no normalisation -- second der > 0)' ,150,-0.2,0.2)
YMinXSPosScdDer  = TH1F('YMinXSPosScdDer', 'Deviation from parabolic fit for Re(VR) = -0.1 (XS normalisation -- second der > 0)' ,150,-0.2,0.2)
YMinAccPosScdDer = TH1F('YMinAccPosScdDer','Deviation from parabolic fit for Re(VR) = -0.1 (Acc normalisation -- second der > 0)',150,-0.2,0.2)
YMinNegScdDer    = TH1F('YMinNegScdDer',   'Deviation from parabolic fit for Re(VR) = -0.1 (no normalisation -- second der < 0)' ,150,-0.2,0.2)
YMinXSNegScdDer  = TH1F('YMinXSNegScdDer', 'Deviation from parabolic fit for Re(VR) = -0.1 (XS normalisation -- second der < 0)' ,150,-0.2,0.2)
YMinAccNegScdDer = TH1F('YMinAccNegScdDer','Deviation from parabolic fit for Re(VR) = -0.1 (Acc normalisation -- second der < 0)',150,-0.2,0.2)

YMinMin    = TH1F('YMinMin',   'Deviation from parabolic fit for Re(VR) = -0.2 (no normalisation)' ,150,-0.5,0.5)
YMinMinXS  = TH1F('YMinMinXS', 'Deviation from parabolic fit for Re(VR) = -0.2 (XS normalisation)' ,150,-0.5,0.5)
YMinMinAcc = TH1F('YMinMinAcc','Deviation from parabolic fit for Re(VR) = -0.2 (Acc normalisation)',150,-0.5,0.5)
YMinMinPosScdDer    = TH1F('YMinMinPosScdDer',   'Deviation from parabolic fit for Re(VR) = -0.2 (no normalisation -- second der > 0)' ,150,-0.2,0.2)
YMinMinXSPosScdDer  = TH1F('YMinMinXSPosScdDer', 'Deviation from parabolic fit for Re(VR) = -0.2 (XS normalisation -- second der > 0)' ,150,-0.2,0.2)
YMinMinAccPosScdDer = TH1F('YMinMinAccPosScdDer','Deviation from parabolic fit for Re(VR) = -0.2 (Acc normalisation -- second der > 0)',150,-0.2,0.2)
YMinMinNegScdDer    = TH1F('YMinMinNegScdDer',   'Deviation from parabolic fit for Re(VR) = -0.2 (no normalisation -- second der < 0)' ,150,-0.2,0.2)
YMinMinXSNegScdDer  = TH1F('YMinMinXSNegScdDer', 'Deviation from parabolic fit for Re(VR) = -0.2 (XS normalisation -- second der < 0)' ,150,-0.2,0.2)
YMinMinAccNegScdDer = TH1F('YMinMinAccNegScdDer','Deviation from parabolic fit for Re(VR) = -0.2 (Acc normalisation -- second der < 0)',150,-0.2,0.2)

YRelPlus    = TH1F('YRelPlus',   'Relative deviation from parabolic fit for Re(VR) = 0.2 (no normalisation)' ,250,-0.1,0.1)
YRelPlusXS  = TH1F('YRelPlusXS', 'Relative deviation from parabolic fit for Re(VR) = 0.2 (XS normalisation)' ,250,-0.1,0.1)
YRelPlusAcc = TH1F('YRelPlusAcc','Relative deviation from parabolic fit for Re(VR) = 0.2 (Acc normalisation)',250,-0.1,0.1)

YRelMin    = TH1F('YRelMin',   'Relative deviation from parabolic fit for Re(VR) = -0.2 (no normalisation)' ,150,-0.1,0.1)
YRelMinXS  = TH1F('YRelMinXS', 'Relative deviation from parabolic fit for Re(VR) = -0.2 (XS normalisation)' ,150,-0.1,0.1)
YRelMinAcc = TH1F('YRelMinAcc','Relative deviation from parabolic fit for Re(VR) = -0.2 (Acc normalisation)',150,-0.1,0.1)

FstDer    = TH1F('FirstDer',   'First derivative of -ln(likelihood) distribution', 4,0,4)
FstDer.GetXaxis().SetBinLabel(1,"(y_{DATA}(x=-0.2) - y_{DATA}(x=-0.1))/0.1")
FstDer.GetXaxis().SetBinLabel(1,"(y_{DATA}(x=-0.1) - y_{DATA}(x=0.0))/0.1")
FstDer.GetXaxis().SetBinLabel(1,"(y_{DATA}(x=0.0) - y_{DATA}(x=0.1))/0.1")
FstDer.GetXaxis().SetBinLabel(1,"(y_{DATA}(x=0.1) - y_{DATA}(x=0.2))/0.1")
FstDerXS  = TH1F('FirstDerivativeXS', 'First derivative of -ln(likelihood) distribution (XS normalisation)', 5,-0.25,0.25)
FstDerAcc = TH1F('FirstDerivativeAcc','First derivative of -ln(likelihood) distribution (Acc normalisation)',5,-0.25,0.25)
ScdDer01    = TH1F('SecondDerivative01',   'Second derivative of -ln(likelihood) distribution (no normalisation -- using x = -0.1/0.0/0.1)', 250,-20,20)
ScdDerXS01  = TH1F('SecondDerivativeXS01', 'Second derivative of -ln(likelihood) distribution (XS normalisation -- using x = -0.1/0.0/0.1)', 250,-20,20)
ScdDerAcc01 = TH1F('SecondDerivativeAcc01','Second derivative of -ln(likelihood) distribution (Acc normalisation -- using x = -0.1/0.0/0.1)',250,-20,20)
ScdDer02    = TH1F('SecondDerivative02',   'Second derivative of -ln(likelihood) distribution (no normalisation -- using x = -0.2/0.0/0.2)', 250,-20,20)
ScdDerXS02  = TH1F('SecondDerivativeXS02', 'Second derivative of -ln(likelihood) distribution (XS normalisation -- using x = -0.2/0.0/0.2)', 250,-20,20)
ScdDerAcc02 = TH1F('SecondDerivativeAcc02','Second derivative of -ln(likelihood) distribution (Acc normalisation -- using x = -0.2/0.0/0.2)',250,-20,20)
ScdDerScatter    = TH2F('ScdDerScatterPlot',   'Second derivative of -ln(L) using Re(VR) = -0.1/0.0/0.1 versus using Re(VR) = -0.2/0.0/0.2 (no normalisation)', 250,-15,15,250,-15,15)
ScdDerXSScatter  = TH2F('ScdDerXSScatterPlot', 'Second derivative of -ln(L) using Re(VR) = -0.1/0.0/0.1 versus using Re(VR) = -0.2/0.0/0.2 (XS normalisation)', 250,-15,15,250,-15,15)
ScdDerAccScatter = TH2F('ScdDerAccScatterPlot','Second derivative of -ln(L) using Re(VR) = -0.1/0.0/0.1 versus using Re(VR) = -0.2/0.0/0.2 (Acc normalisation)',250,-15,15,250,-15,15)

LnLogDist = TH1F("LnLog","title",7,-0.35,0.35)
LnLogDist.SetMarkerStyle(20), LnLogDist.SetLineColor(1), LnLogDist.SetMarkerColor(1), LnLogDist.SetMarkerSize(1.2)
LnLogXSDist = TH1F("LnLogXS","title",7,-0.35,0.35)
LnLogXSDist.SetMarkerStyle(21), LnLogXSDist.SetLineColor(3), LnLogXSDist.SetMarkerColor(3), LnLogXSDist.SetMarkerSize(1.2)
LnLogAccDist = TH1F("LnLogAcc","title",7,-0.35,0.35)
LnLogAccDist.SetMarkerStyle(22), LnLogAccDist.SetLineColor(4), LnLogAccDist.SetMarkerColor(4), LnLogAccDist.SetMarkerSize(1.2)

FitComp = Tfile.mkdir("FitComparison")
LnLogDir = Tfile.mkdir("LnLogDistr")
LnLogXSDir = Tfile.mkdir("LnLogXSDistr")
LnLogAccDir = Tfile.mkdir("LnLogAccDistr")
FstDerDir = Tfile.mkdir("FirstDerivativeDistr")
LnLogFitCanvas = TCanvas('name','title')
LnLogFitCanvasName = "LnLogFitCanvas_Evt"
LnLogFitCanvasTitle = "Comparing fit function from ROOT fit and algebraic function for event "

#Create the arrays where the likelihood values will be stored
NrConfigs = 9
LnLog, LnLogXS, LnLogAcc = [], [], []
for ii in range(NrConfigs):
  LnLog.append(0), LnLogXS.append(0), LnLogAcc.append(0)

#---  Create arrays where the events passing specific cuts will be stored  ---#
EvtsWithPosScdDer01, EvtsWithPosScdDerXS01, EvtsWithPosScdDerAcc01 = [], [], []
EvtsWithPosScdDer02, EvtsWithPosScdDerXS02, EvtsWithPosScdDerAcc02 = [], [], []
EvtsWithYPlusGausSmall, EvtsWithYPlusGausSmallXS, EvtsWithYPlusGausSmallAcc = [], [], []
EvtsWithBothYPos, EvtsWithBothYPosXS, EvtsWithBothYPosAcc = [],[],[]

nEvts = 1000
print " **** Will be using ",nEvts," events!!"
#Loop over all lines in weights file:
for WeightLine in WeightsFile:
  WeightWord = WeightLine.split()
  #Only interested in files starting with a number
  if str(WeightWord[0]) != "#" and str(WeightWord[3]) != "0.0" :
    for iEvt in range(nEvts):
      if str(WeightWord[0]) == str(iEvt+1):                                                               #Look at one single event!
        if str(WeightWord[1]) == "1":
          #LogMin = [-log(float(WeightWord[3])),-log(float(WeightWord[3]))+log(MGXS[int(WeightWord[1])-1]), -log(float(WeightWord[3]))+log(MGXS[int(WeightWord[1])-1])+log(Acceptance[int(WeightWord[1])-1])]
          cHat, bHat01, bHat02, aHat01, aHat02 = [0,0,0], [0,0,0], [0,0,0], [0,0,0], [0,0,0]
          LnLogDistFit = TH1F("LnLogFit_Evt"+str(int(iEvt)+1),"LnLog distribution for fit (event "+str(int(iEvt)+1)+")",11,-0.55,0.55)
          LnLogDistFit.SetMarkerStyle(20), LnLogDistFit.SetLineColor(1), LnLogDistFit.SetMarkerColor(1), LnLogDistFit.SetMarkerSize(1.2)
          LnLogDist.SetName("LnLog_Evt"+str(int(iEvt)+1)),       LnLogDist.SetTitle("LnLog distribution for event "+str(int(iEvt)+1))
          LnLogXSDist.SetName("LnLogXS_Evt"+str(int(iEvt)+1)),   LnLogXSDist.SetTitle("LnLogXS distribution for event "+str(int(iEvt)+1))
          LnLogAccDist.SetName("LnLogAcc_Evt"+str(int(iEvt)+1)), LnLogAccDist.SetTitle("LnLogAcc distribution for event "+str(int(iEvt)+1))
        LnLog[int(WeightWord[1])-1] = -log(float(WeightWord[3]))
        LnLogXS[int(WeightWord[1])-1] = -log(float(WeightWord[3])) + log(MGXS[int(WeightWord[1])-1])
        LnLogAcc[int(WeightWord[1])-1] = -log(float(WeightWord[3])) + log(MGXS[int(WeightWord[1])-1]) + log(Acceptance[int(WeightWord[1])-1])
        #---  Fill the LnLog histograms for each event  ---#
        LnLogDistFit.SetBinContent(LnLogDistFit.FindBin(RVR[int(WeightWord[1])-1]), LnLog[int(WeightWord[1])-1])
        if str(WeightWord[1]) != "1" or str(WeightWord[1]) != "9":                                             #Don't fill the bins coresponding to -0.5 & 0.5 (= too wide range)
          LnLogDist.SetBinContent(LnLogDist.FindBin(RVR[int(WeightWord[1])-1]), LnLog[int(WeightWord[1])-1])
          LnLogXSDist.SetBinContent(LnLogXSDist.FindBin(RVR[int(WeightWord[1])-1]), LnLogXS[int(WeightWord[1])-1])
          LnLogAccDist.SetBinContent(LnLogAccDist.FindBin(RVR[int(WeightWord[1])-1]), LnLogAcc[int(WeightWord[1])-1])
        #---  Only perform the fit after all 9 configurations are considered!  ---#
        if str(WeightWord[1]) == str(NrConfigs):
          #---  Calculate the fit parameters (polynomial)  ---#
          cHat   = [LnLog[4],                            LnLogXS[4],                                LnLogAcc[4]                                 ]
          bHat01 = [5*(LnLog[5]-LnLog[3]),               5*(LnLogXS[5]-LnLogXS[3]),                 5*(LnLogAcc[5]-LnLogAcc[3])                 ]
          aHat01 = [50*(LnLog[3]+LnLog[5]-2*LnLog[4]),   50*(LnLogXS[3]+LnLogXS[5]-2*LnLogXS[4]),   50*(LnLogAcc[3]+LnLogAcc[5]-2*LnLogAcc[4])  ] 
          bHat02 = [5*(LnLog[6]-LnLog[2])/2,             5*(LnLogXS[6]-LnLogXS[2])/2,               5*(LnLogAcc[6]-LnLogAcc[2])/2               ]
          aHat02 = [25*(LnLog[2]+LnLog[6]-2*LnLog[4])/2, 25*(LnLogXS[2]+LnLogXS[6]-2*LnLogXS[4])/2, 50*(LnLogAcc[2]+LnLogAcc[6]-2*LnLogAcc[4])/2] 
          LnLogFunction01, LnLogFunction02 = array('d'), array('d')
          for ii in range(9):
            LnLogFunction01.append( aHat01[0]*float(RVR[ii])*float(RVR[ii])+bHat01[0]*float(RVR[ii])+cHat[0])
            LnLogFunction02.append( aHat02[0]*float(RVR[ii])*float(RVR[ii])+bHat02[0]*float(RVR[ii])+cHat[0])
            #LnLogXSFunction.append( aHat01[1]*float(RVR[ii])*float(RVR[ii])+bHat01[1]*float(RVR[ii])+cHat[1])
            #LnLogAccFunction.append( aHat01[2]*float(RVR[ii])*float(RVR[ii])+bHat01[2]*float(RVR[ii])+cHat[2])
          LnLogGraph01 = TGraph(9, RVR, LnLogFunction01)
          LnLogGraph01.SetTitle('Comparison between ROOT fit and algebraic method')
          LnLogGraph01.GetXaxis().SetTitle('#Delta ln(Likelihood)')
          LnLogGraph01.GetYaxis().SetTitle('Number of events')
          LnLogGraph01.SetMarkerColor(6), LnLogGraph01.SetLineColor(6), LnLogGraph01.SetMarkerSize(1.2)
          LnLogGraph02 = TGraph(9,RVR, LnLogFunction02)
          LnLogGraph02.SetMarkerColor(7), LnLogGraph02.SetLineColor(7), LnLogGraph02.SetMarkerSize(1.2)	
          #---  Compare with TF1 fit from ROOT  ---#
          LnLogDistFit.Fit("pol2","Q","",-0.35, 0.35)
          LnLogFitCanvas.SetName(LnLogFitCanvasName+str(int(iEvt)+1))
          LnLogFitCanvas.SetTitle(LnLogFitCanvasTitle+str(int(iEvt)+1))
          legendLnLog = TLegend(0.75,0.7,0.95,0.95)
          LnLogFitCanvas.cd()
          legendLnLog.AddEntry(LnLogGraph01,'Fit result using algebraic function (x = -0.1/0.0/0.1)',"p")
          legendLnLog.AddEntry(LnLogGraph02,'Fit result using algebraic function (x = -0.2/0.0/0.2)',"p")
          legendLnLog.AddEntry(LnLogDistFit,'Obtained LnLog values with ROOT fit',"p")
          LnLogGraph01.Draw("AC*"), LnLogGraph02.Draw("C*"), LnLogDistFit.Draw("samep"), legendLnLog.Draw()
          FitComp.cd(), LnLogFitCanvas.Write()
          #---  Calculate the first derivative distribution and save them in event-by-event plot  ---#
          FstDer.SetName("FirstDerivative_Evt"+str(int(iEvt)+1))
          FstDer.SetTitle("First derivative of -ln(likelihood) distribution for event"+str(int(iEvt)+1)+" (no normalisation)") 
          for ii in range(FstDer.GetNbinsX()):
            FstDer.SetBinContent(ii+1, (LnLog[ii+2] - LnLog[ii+3])/0.1)
          FstDerDir.cd(), FstDer.Write()
          #---  Calculate the positive and negative deviation  ---#
          FitOutputPlus, FitOutputMin = [0,0,0],[0,0,0]
          FitOutputPlusPlus, FitOutputMinMin = [0,0,0],[0,0,0]
          for ii in range(3):
            FitOutputPlusPlus[ii] = (aHat01[ii]*float(RVR[6])*float(RVR[6])+bHat01[ii]*float(RVR[6])+cHat[ii])
            FitOutputMinMin[ii]   = (aHat01[ii]*float(RVR[2])*float(RVR[2])+bHat01[ii]*float(RVR[2])+cHat[ii])
            FitOutputPlus[ii] = (aHat02[ii]*float(RVR[5])*float(RVR[5])+bHat02[ii]*float(RVR[5])+cHat[ii])
            FitOutputMin[ii]  = (aHat02[ii]*float(RVR[3])*float(RVR[3])+bHat02[ii]*float(RVR[3])+cHat[ii])
          yPlusPlus = [LnLog[6] - FitOutputPlusPlus[0], LnLogXS[6] - FitOutputPlusPlus[1], LnLogAcc[6] - FitOutputPlusPlus[2]]
          yMinMin   = [LnLog[2] - FitOutputMinMin[0],  LnLogXS[2] - FitOutputMinMin[1],  LnLogAcc[2] - FitOutputMinMin[2] ]
          yPlus = [LnLog[5] - FitOutputPlus[0], LnLogXS[5] - FitOutputPlus[1], LnLogAcc[5] - FitOutputPlus[2]]
          yMin  = [LnLog[3] - FitOutputMin[0],  LnLogXS[3] - FitOutputMin[1],  LnLogAcc[3] - FitOutputMin[2] ]
          #--- Fill the histograms  ---#
          YPlus.Fill(yPlus[0]),                             YPlusXS.Fill(yPlus[1]),                                   YPlusAcc.Fill(yPlus[2])
          YPlusPlus.Fill(yPlusPlus[0]),                     YPlusPlusXS.Fill(yPlusPlus[1]),                           YPlusPlusAcc.Fill(yPlusPlus[2])
          YPlusGausTest.Fill(yPlus[0] + yPlusPlus[0]/4),    YPlusGausTestXS.Fill(yPlus[1] + yPlusPlus[1]/4),          YPlusGausTestAcc.Fill(yPlus[2] + yPlusPlus[2]/4)
          #YRelPlus.Fill(yPlus[0]/LnLog[6]),                 YRelPlusXS.Fill(yPlus[1]/LnLogXS[6]),                     YRelPlusAcc.Fill(yPlus[2]/LnLogAcc[6])
          YMin.Fill(yMin[0]),                               YMinXS.Fill(yMin[1]),                                     YMinAcc.Fill(yMin[2])
          YMinMin.Fill(yMinMin[0]),                         YMinMinXS.Fill(yMinMin[1]),                               YMinMinAcc.Fill(yMinMin[2])
          #YRelMin.Fill(yMin[0]/LnLog[2]),                   YRelMinXS.Fill(yMin[1]/LnLogXS[2]),                       YRelMinAcc.Fill(yMin[2]/LnLogAcc[2])
          ScdDer01.Fill(10*(LnLog[3]-2*LnLog[4]+LnLog[5])), ScdDerXS01.Fill(10*(LnLogXS[3]-2*LnLogXS[4]+LnLogXS[5])), ScdDerAcc01.Fill(10*(LnLogAcc[3]-2*LnLogAcc[4]+LnLogAcc[5]))
          ScdDer02.Fill(5*(LnLog[2]-2*LnLog[4]+LnLog[6])),  ScdDerXS02.Fill(5*(LnLogXS[2]-2*LnLogXS[4]+LnLogXS[6])),  ScdDerAcc02.Fill(5*(LnLogAcc[2]-2*LnLogAcc[4]+LnLogAcc[6]))
          ScdDerScatter.Fill(5*(LnLog[2]-2*LnLog[4]+LnLog[6]) , 10*(LnLog[3]-2*LnLog[4]+LnLog[5]) )
          ScdDerXSScatter.Fill(5*(LnLogXS[2]-2*LnLogXS[4]+LnLogXS[6]) , 10*(LnLogXS[3]-2*LnLogXS[4]+LnLogXS[5]) )
          ScdDerAccScatter.Fill(5*(LnLogAcc[2]-2*LnLogAcc[4]+LnLogAcc[6]) , 10*(LnLogAcc[3]-2*LnLogAcc[4]+LnLogAcc[5]) )
          #---  Save the LnLog distributions (event-per-event) and the Y-deviations in histograms  ---#
          LnLogDir.cd(),    LnLogDist.Write()
          LnLogXSDir.cd(),  LnLogXSDist.Write()
          LnLogAccDir.cd(), LnLogAccDist.Write()
          #-- Apply cut on YPlusGausTest --#
          if (yPlus[0] + yPlusPlus[0]/4) <= 0.025 and (yPlus[0] + yPlusPlus[0]/4) >= -0.025: EvtsWithYPlusGausSmall.append(iEvt+1)
          if (yPlus[1] + yPlusPlus[1]/4) <= 0.025 and (yPlus[1] + yPlusPlus[1]/4) >= -0.025: EvtsWithYPlusGausSmallXS.append(iEvt+1)
          if (yPlus[2] + yPlusPlus[2]/4) <= 0.025 and (yPlus[2] + yPlusPlus[2]/4) >= -0.025: EvtsWithYPlusGausSmallAcc.append(iEvt+1)
          #-- Apply cut on ScdDer (using x = -0.1/0.0/0.1) --#
          if (10*(LnLog[3]-2*LnLog[4]+LnLog[5])) > 0.0: EvtsWithPosScdDer01.append(iEvt+1)
          if (10*(LnLogXS[3]-2*LnLogXS[4]+LnLogXS[5])) > 0.0: EvtsWithPosScdDerXS01.append(iEvt+1)
          if (10*(LnLogAcc[3]-2*LnLogAcc[4]+LnLogAcc[5])) > 0.0: EvtsWithPosScdDerAcc01.append(iEvt+1)
          #-- Apply cut on ScdDer (using x = -0.2/0.0/0.2) --#
          if (5*(LnLog[2]-2*LnLog[4]+LnLog[6])) > 0.0:
            EvtsWithPosScdDer02.append(iEvt+1)
            YPlusGausTestPosScdDer.Fill(yPlus[0] + yPlusPlus[0]/4)
            YPlusPosScdDer.Fill(yPlus[0]), YPlusPlusPosScdDer.Fill(yPlusPlus[0])
            YMinPosScdDer.Fill(yMin[0]),   YMinMinPosScdDer.Fill(yMinMin[0])
            if yPlusPlus[0] > 0.0:
              YPlusPosScdDerPosYPlusPlus.Fill(yPlus[0])
              if yPlus[0] > 0.0: EvtsWithBothYPos.append(iEvt+1)
          else:
            YPlusGausTestNegScdDer.Fill(yPlus[0] + yPlusPlus[0]/4)
            YPlusNegScdDer.Fill(yPlus[0]), YPlusPlusNegScdDer.Fill(yPlusPlus[0])
            YMinNegScdDer.Fill(yMin[0]),   YMinMinNegScdDer.Fill(yMinMin[0])
          if (5*(LnLogXS[2]-2*LnLogXS[4]+LnLogXS[6])) > 0.0:
            EvtsWithPosScdDerXS02.append(iEvt+1)
            YPlusGausTestXSPosScdDer.Fill(yPlus[1] + yPlusPlus[1]/4)
            YPlusXSPosScdDer.Fill(yPlus[1]), YPlusPlusXSPosScdDer.Fill(yPlusPlus[1])
            YMinXSPosScdDer.Fill(yMin[1]),   YMinMinXSPosScdDer.Fill(yMinMin[1])
            if yPlusPlus[1] > 0.0:
              YPlusXSPosScdDerPosYPlusPlus.Fill(yPlus[1])
              if yPlus[1] > 0.0: EvtsWithBothYPosXS.append(iEvt+1)
          else:
            YPlusGausTestXSNegScdDer.Fill(yPlus[1] + yPlusPlus[1]/4)
            YPlusXSNegScdDer.Fill(yPlus[1]), YPlusPlusXSNegScdDer.Fill(yPlusPlus[1])
            YMinXSNegScdDer.Fill(yMin[1]),   YMinMinXSNegScdDer.Fill(yMinMin[1])
          if (5*(LnLogAcc[2]-2*LnLogAcc[4]+LnLogAcc[6])) > 0.0:
            EvtsWithPosScdDerAcc02.append(iEvt+1)
            YPlusGausTestAccPosScdDer.Fill(yPlus[2] + yPlusPlus[2]/4)
            YPlusAccPosScdDer.Fill(yPlus[2]), YPlusPlusAccPosScdDer.Fill(yPlusPlus[2])
            YMinAccPosScdDer.Fill(yMin[2]),   YMinMinAccPosScdDer.Fill(yMinMin[2])
            if yPlusPlus[2] > 0.0:
              YPlusAccPosScdDerPosYPlusPlus.Fill(yPlus[2])
              if yPlus[2] > 0.0: EvtsWithBothYPosAcc.append(iEvt+1)
          else:
            YPlusGausTestAccNegScdDer.Fill(yPlus[2] + yPlusPlus[2]/4)
            YPlusAccNegScdDer.Fill(yPlus[2]), YPlusPlusAccNegScdDer.Fill(yPlusPlus[2])
            YMinAccNegScdDer.Fill(yMin[2]),   YMinMinAccNegScdDer.Fill(yMinMin[2])


#--- Save all the histograms containing information about all the events! ---#
Tfile.cd()
YPlusPosScdDerPosYPlusPlus.Write(), YPlusXSPosScdDerPosYPlusPlus.Write(), YPlusAccPosScdDerPosYPlusPlus.Write()  #How can they be both positive ... (y+ & y++)
print "List of events with both y+ and y++ positive (no norm) = ",EvtsWithBothYPos
print "\n List of events with both y+ and y++ positive (XS norm) = ",EvtsWithBothYPosXS
print "\n List of events with both y+ and y++ positive (Acc norm) = ",EvtsWithBothYPosAcc

YPlusGausTest.Write(), YPlusGausTestXS.Write(), YPlusGausTestAcc.Write()
YPlus.Write(),         YPlusXS.Write(),         YPlusAcc.Write()
YPlusPlus.Write(),     YPlusPlusXS.Write(),     YPlusPlusAcc.Write()
YRelPlus.Write(),      YRelPlusXS.Write(),      YRelPlusAcc.Write()
YMin.Write(),          YMinXS.Write(),          YMinAcc.Write()
YMinMin.Write(),       YMinMinXS.Write(),       YMinMinAcc.Write()
YRelMin.Write(),       YRelMinXS.Write(),       YRelMinAcc.Write()
ScdDer01.Write(),      ScdDerXS01.Write(),      ScdDerAcc01.Write()
ScdDer02.Write(),      ScdDerXS02.Write(),      ScdDerAcc02.Write()
ScdDerScatter.Write(), ScdDerXSScatter.Write(), ScdDerAccScatter.Write()

#---  Draw the likelihood distribution separately for events surviving and passing the cuts!  ---#
print "Nr of events with positive second derivative (LnLog, LnLogXS & LnLogAcc -- using x = -0.1/0.0/0.1) ",len(EvtsWithPosScdDer01),", ",len(EvtsWithPosScdDerXS01)," & ",len(EvtsWithPosScdDerAcc01)
print "Nr of events with positive second derivative (LnLog, LnLogXS & LnLogAcc -- using x = -0.2/0.0/0.2) ",len(EvtsWithPosScdDer02),", ",len(EvtsWithPosScdDerXS02)," & ",len(EvtsWithPosScdDerAcc02)
print "Nr of events with Gaussianse vergelijking voor + (LnLog, LnLogXS & LnLogAcc) ", len(EvtsWithYPlusGausSmall),", ",len(EvtsWithYPlusGausSmallXS)," & ",len(EvtsWithYPlusGausSmallAcc)
LLPosScdDer01,   LLNegScdDer01,   LLXSPosScdDer01,   LLXSNegScdDer01,   LLAccPosScdDer01,   LLAccNegScdDer01 = [],[],[],[],[],[]
LLPosScdDer02,   LLNegScdDer02,   LLXSPosScdDer02,   LLXSNegScdDer02,   LLAccPosScdDer02,   LLAccNegScdDer02 = [],[],[],[],[],[]
LLPosScdDerBoth, LLNegScdDerBoth, LLXSPosScdDerBoth, LLXSNegScdDerBoth, LLAccPosScdDerBoth, LLAccNegScdDerBoth = [],[],[],[],[],[]
LLSmallGausDiff, LLBigGausDiff, LLXSSmallGausDiff, LLXSBigGausDiff, LLAccSmallGausDiff, LLAccBigGausDiff = [],[],[],[],[],[]
for ii in range(len(RVR)):
  LLPosScdDer01.append(0),   LLNegScdDer01.append(0),   LLXSPosScdDer01.append(0),   LLXSNegScdDer01.append(0),   LLAccPosScdDer01.append(0),   LLAccNegScdDer01.append(0)
  LLPosScdDer02.append(0),   LLNegScdDer02.append(0),   LLXSPosScdDer02.append(0),   LLXSNegScdDer02.append(0),   LLAccPosScdDer02.append(0),   LLAccNegScdDer02.append(0)
  LLPosScdDerBoth.append(0), LLNegScdDerBoth.append(0), LLXSPosScdDerBoth.append(0), LLXSNegScdDerBoth.append(0), LLAccPosScdDerBoth.append(0), LLAccNegScdDerBoth.append(0)
  LLSmallGausDiff.append(0), LLBigGausDiff.append(0), LLXSSmallGausDiff.append(0), LLXSBigGausDiff.append(0), LLAccSmallGausDiff.append(0), LLAccBigGausDiff.append(0)

EvtsWithPosScdDerBoth = 0
EvtsWithPosScdDerXSBoth = 0
EvtsWithPosScdDerAccBoth = 0
#Loop over all lines in weights file:
for LikelihoodLine in LikelihoodFile:
  LWord = LikelihoodLine.split()
  #Only interested in files starting with a number
  if str(LWord[0]) != "#" and str(LWord[3]) != "0.0" :
    #---  Separate the events with positive and negative second derivative (using both x = -0.1/0.0/0.1 and x = -0.2/0.0/0.2) ---#
    if int(LWord[0]) in EvtsWithPosScdDer01 and int(LWord[0]) in EvtsWithPosScdDer02:
      LLPosScdDerBoth[int(LWord[1])-1] = LLPosScdDerBoth[int(LWord[1])-1]-log(float(LWord[3]))
      if LWord[1] == "1": EvtsWithPosScdDerBoth += 1
    if int(LWord[0]) in EvtsWithPosScdDerXS01 and int(LWord[0]) in EvtsWithPosScdDerXS02:
      LLXSPosScdDerBoth[int(LWord[1])-1] = LLXSPosScdDerBoth[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])
      if LWord[1] == "1": EvtsWithPosScdDerXSBoth += 1
    if int(LWord[0]) in EvtsWithPosScdDerAcc01 and int(LWord[0]) in EvtsWithPosScdDerAcc02:
      LLAccPosScdDerBoth[int(LWord[1])-1] = LLAccPosScdDerBoth[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])+log(Acceptance[int(LWord[1])-1])
      if LWord[1] == "1": EvtsWithPosScdDerAccBoth += 1
      ##print EvtsWithPosScdDerAccBoth,") Looking at event : ",int(LWord[0])
    #---  Separate the events with positive and negative second derivative (using x = -0.1/0.0/0.1) ---#
    if int(LWord[0]) in EvtsWithPosScdDer01:    LLPosScdDer01[int(LWord[1])-1] = LLPosScdDer01[int(LWord[1])-1]-log(float(LWord[3]))
    else:                                       LLNegScdDer01[int(LWord[1])-1] = LLNegScdDer01[int(LWord[1])-1]-log(float(LWord[3]))
    if int(LWord[0]) in EvtsWithPosScdDerXS01:  LLXSPosScdDer01[int(LWord[1])-1] = LLXSPosScdDer01[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])
    else:                                       LLXSNegScdDer01[int(LWord[1])-1] = LLXSNegScdDer01[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])
    if int(LWord[0]) in EvtsWithPosScdDerAcc01: LLAccPosScdDer01[int(LWord[1])-1] = LLAccPosScdDer01[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])+log(Acceptance[int(LWord[1])-1])
    else:                                       LLAccNegScdDer01[int(LWord[1])-1] = LLAccNegScdDer01[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])+log(Acceptance[int(LWord[1])-1])
    #---  Separate the events with positive and negative second derivative (using x = -0.2/0.0/0.2) ---#
    if int(LWord[0]) in EvtsWithPosScdDer02:    LLPosScdDer02[int(LWord[1])-1] = LLPosScdDer02[int(LWord[1])-1]-log(float(LWord[3]))
    else:                                       LLNegScdDer02[int(LWord[1])-1] = LLNegScdDer02[int(LWord[1])-1]-log(float(LWord[3]))
    if int(LWord[0]) in EvtsWithPosScdDerXS02:  LLXSPosScdDer02[int(LWord[1])-1] = LLXSPosScdDer02[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])
    else:                                       LLXSNegScdDer02[int(LWord[1])-1] = LLXSNegScdDer02[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])
    if int(LWord[0]) in EvtsWithPosScdDerAcc02: LLAccPosScdDer02[int(LWord[1])-1] = LLAccPosScdDer02[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])+log(Acceptance[int(LWord[1])-1])
    else:                                       LLAccNegScdDer02[int(LWord[1])-1] = LLAccNegScdDer02[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])+log(Acceptance[int(LWord[1])-1])
    #---  Separate the events with small and big Gaussian difference!  ---#
    if int(LWord[0]) in EvtsWithYPlusGausSmall:    LLSmallGausDiff[int(LWord[1])-1] = LLSmallGausDiff[int(LWord[1])-1] -log(float(LWord[3]))
    else:                                          LLBigGausDiff[int(LWord[1])-1] = LLBigGausDiff[int(LWord[1])-1] -log(float(LWord[3]))
    if int(LWord[0]) in EvtsWithYPlusGausSmallXS:  LLXSSmallGausDiff[int(LWord[1])-1] = LLXSSmallGausDiff[int(LWord[1])-1] -log(float(LWord[3])) + log(MGXS[int(LWord[1])-1])
    else:                                          LLXSBigGausDiff[int(LWord[1])-1] = LLXSBigGausDiff[int(LWord[1])-1] -log(float(LWord[3])) + log(MGXS[int(LWord[1])-1])
    if int(LWord[0]) in EvtsWithYPlusGausSmallAcc: LLAccSmallGausDiff[int(LWord[1])-1] = LLAccSmallGausDiff[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])+log(Acceptance[int(LWord[1])-1])
    else:                                          LLAccBigGausDiff[int(LWord[1])-1] = LLAccBigGausDiff[int(LWord[1])-1]-log(float(LWord[3]))+ log(MGXS[int(LWord[1])-1]) + log(Acceptance[int(LWord[1])-1])

LLPosScdDerDistBoth   =TH1F('LLPosScdDerBoth',   '-ln(L) for events with 2nd derivative > 0 (no norm -- using x = 0.1 & 0.2 -- '+str(EvtsWithPosScdDerBoth)+'/'+str(nEvts)+' evts)',    7,-0.35,0.35)
LLXSPosScdDerDistBoth =TH1F('LLXSPosScdDerBoth', '-ln(L) for events with 2nd derivative > 0 (XS norm -- using x = 0.1 & 0.2 -- '+str(EvtsWithPosScdDerXSBoth)+'/'+str(nEvts)+' evts)',  7,-0.35,0.35)
LLAccPosScdDerDistBoth=TH1F('LLAccPosScdDerBoth','-ln(L) for events with 2nd derivative > 0 (Acc norm -- using x = 0.1 & 0.2 -- '+str(EvtsWithPosScdDerAccBoth)+'/'+str(nEvts)+' evts)',7,-0.35,0.35)
LLPosScdDerDist01   =TH1F('LLPosScdDer01',   '-ln(L) for events with 2nd derivative > 0 (no norm -- x = -0.1/0.0/0.1 -- '+str(len(EvtsWithPosScdDer01))+'/'+str(nEvts)+' evts)',    7,-0.35,0.35)
LLXSPosScdDerDist01 =TH1F('LLXSPosScdDer01', '-ln(L) for events with 2nd derivative > 0 (XS norm -- x = -0.1/0.0/0.1 -- '+str(len(EvtsWithPosScdDerXS01))+'/'+str(nEvts)+' evts)',  7,-0.35,0.35)
LLAccPosScdDerDist01=TH1F('LLAccPosScdDer01','-ln(L) for events with 2nd derivative > 0 (Acc norm -- x = -0.1/0.0/0.1 -- '+str(len(EvtsWithPosScdDerAcc01))+'/'+str(nEvts)+' evts)',7,-0.35,0.35)
LLPosScdDerDist02   =TH1F('LLPosScdDer02',   '-ln(L) for events with 2nd derivative > 0 (no norm -- x = -0.2/0.0/0.2 -- '+str(len(EvtsWithPosScdDer02))+'/'+str(nEvts)+' evts)',    7,-0.35,0.35)
LLXSPosScdDerDist02 =TH1F('LLXSPosScdDer02', '-ln(L) for events with 2nd derivative > 0 (XS norm -- x = -0.2/0.0/0.2 -- '+str(len(EvtsWithPosScdDerXS02))+'/'+str(nEvts)+' evts)',  7,-0.35,0.35)
LLAccPosScdDerDist02=TH1F('LLAccPosScdDer02','-ln(L) for events with 2nd derivative > 0 (Acc norm -- x = -0.2/0.0/0.2 -- '+str(len(EvtsWithPosScdDerAcc02))+'/'+str(nEvts)+' evts)',7,-0.35,0.35)
LLNegScdDerDist01   =TH1F('LLNegScdDer01',   '-ln(L) for events with 2nd derivative < 0 (no norm -- x = -0.1/0.0/0.1 -- '+str(int(nEvts)-len(EvtsWithPosScdDer01))+'/'+str(nEvts)+' evts)',    7,-0.35,0.35)
LLXSNegScdDerDist01 =TH1F('LLXSNegScdDer01', '-ln(L) for events with 2nd derivative < 0 (XS norm -- x = -0.1/0.0/0.1 -- '+str(int(nEvts)-len(EvtsWithPosScdDerXS01))+'/'+str(nEvts)+' evts)',  7,-0.35,0.35)
LLAccNegScdDerDist01=TH1F('LLAccNegScdDer01','-ln(L) for events with 2nd derivative < 0 (Acc norm -- x = -0.1/0.0/0.1 -- '+str(int(nEvts)-len(EvtsWithPosScdDerAcc01))+'/'+str(nEvts)+' evts)',7,-0.35,0.35)
LLNegScdDerDist02   =TH1F('LLNegScdDer02',   '-ln(L) for events with 2nd derivative < 0 (no norm -- x = -0.2/0.0/0.2 -- '+str(int(nEvts)-len(EvtsWithPosScdDer01))+'/'+str(nEvts)+' evts)',    7,-0.35,0.35)
LLXSNegScdDerDist02 =TH1F('LLXSNegScdDer02', '-ln(L) for events with 2nd derivative < 0 (XS norm -- x = -0.2/0.0/0.2 -- '+str(int(nEvts)-len(EvtsWithPosScdDerXS01))+'/'+str(nEvts)+' evts)',  7,-0.35,0.35)
LLAccNegScdDerDist02=TH1F('LLAccNegScdDer02','-ln(L) for events with 2nd derivative < 0 (Acc norm -- x = -0.2/0.0/0.2 -- '+str(int(nEvts)-len(EvtsWithPosScdDerAcc01))+'/'+str(nEvts)+' evts)',7,-0.35,0.35)

LLSmallGausDiffDist    = TH1F('LLSmallGausDiffDist',   '-ln(L) for events with small (-0.025 < x < 0.025) Gaussian difference (no normalisation -- '+str(len(EvtsWithYPlusGausSmall))+'/'+str(nEvts)+' events)',    7,-0.35,0.35)
LLXSSmallGausDiffDist  = TH1F('LLXSSmallGausDiffDist', '-ln(L) for events with small (-0.025 < x < 0.025) Gaussian difference (XS normalisation -- '+str(len(EvtsWithYPlusGausSmallXS))+'/'+str(nEvts)+' events)',  7,-0.35,0.35)
LLAccSmallGausDiffDist = TH1F('LLAccSmallGausDiffDist','-ln(L) for events with small (-0.025 < x < 0.025) Gaussian difference (Acc normalisation -- '+str(len(EvtsWithYPlusGausSmallAcc))+'/'+str(nEvts)+' events)',7,-0.35,0.35)
LLBigGausDiffDist    = TH1F('LLBigGausDiffDist',   '-ln(L) for events with big (x > +-0.025) Gaussian difference (no normalisation -- '+str(int(nEvts)-len(EvtsWithYPlusGausSmall))+'/'+str(nEvts)+' events)',    7,-0.35,0.35)
LLXSBigGausDiffDist  = TH1F('LLXSBigGausDiffDist', '-ln(L) for events with big (x > +-0.025) Gaussian difference (XS normalisation -- '+str(int(nEvts)-len(EvtsWithYPlusGausSmallXS))+'/'+str(nEvts)+' events)',  7,-0.35,0.35)
LLAccBigGausDiffDist = TH1F('LLAccBigGausDiffDist','-ln(L) for events with big (x > +-0.025) Gaussian difference (Acc normalisation -- '+str(int(nEvts)-len(EvtsWithYPlusGausSmallAcc))+'/'+str(nEvts)+' events)',7,-0.35,0.35)

for ii in range(len(RVR)):
  LLPosScdDerDist01.SetBinContent(LLPosScdDerDist01.FindBin(RVR[ii]), float(LLPosScdDer01[ii]))
  LLXSPosScdDerDist01.SetBinContent(LLXSPosScdDerDist01.FindBin(RVR[ii]),float(LLXSPosScdDer01[ii]))
  LLAccPosScdDerDist01.SetBinContent(LLAccPosScdDerDist01.FindBin(RVR[ii]),float(LLAccPosScdDer01[ii]))
  LLNegScdDerDist01.SetBinContent(LLNegScdDerDist01.FindBin(RVR[ii]), float(LLNegScdDer01[ii]))
  LLXSNegScdDerDist01.SetBinContent(LLXSNegScdDerDist01.FindBin(RVR[ii]),float(LLXSNegScdDer01[ii]))
  LLAccNegScdDerDist01.SetBinContent(LLAccNegScdDerDist01.FindBin(RVR[ii]),float(LLAccNegScdDer01[ii]))

  LLPosScdDerDist02.SetBinContent(LLPosScdDerDist02.FindBin(RVR[ii]), float(LLPosScdDer02[ii]))
  LLXSPosScdDerDist02.SetBinContent(LLXSPosScdDerDist02.FindBin(RVR[ii]),float(LLXSPosScdDer02[ii]))
  LLAccPosScdDerDist02.SetBinContent(LLAccPosScdDerDist02.FindBin(RVR[ii]),float(LLAccPosScdDer02[ii]))
  LLNegScdDerDist02.SetBinContent(LLNegScdDerDist02.FindBin(RVR[ii]), float(LLNegScdDer02[ii]))
  LLXSNegScdDerDist02.SetBinContent(LLXSNegScdDerDist02.FindBin(RVR[ii]),float(LLXSNegScdDer02[ii]))
  LLAccNegScdDerDist02.SetBinContent(LLAccNegScdDerDist02.FindBin(RVR[ii]),float(LLAccNegScdDer02[ii]))

  LLPosScdDerDistBoth.SetBinContent(LLPosScdDerDistBoth.FindBin(RVR[ii]), float(LLPosScdDerBoth[ii]))
  LLXSPosScdDerDistBoth.SetBinContent(LLXSPosScdDerDistBoth.FindBin(RVR[ii]),float(LLXSPosScdDerBoth[ii]))
  LLAccPosScdDerDistBoth.SetBinContent(LLAccPosScdDerDistBoth.FindBin(RVR[ii]),float(LLAccPosScdDerBoth[ii]))

  LLSmallGausDiffDist.SetBinContent(LLSmallGausDiffDist.FindBin(RVR[ii]),float(LLSmallGausDiff[ii]))
  LLBigGausDiffDist.SetBinContent(LLBigGausDiffDist.FindBin(RVR[ii]),float(LLBigGausDiff[ii]))
  LLXSSmallGausDiffDist.SetBinContent(LLXSSmallGausDiffDist.FindBin(RVR[ii]),float(LLXSSmallGausDiff[ii]))
  LLXSBigGausDiffDist.SetBinContent(LLXSBigGausDiffDist.FindBin(RVR[ii]),float(LLXSBigGausDiff[ii]))
  LLAccSmallGausDiffDist.SetBinContent(LLAccSmallGausDiffDist.FindBin(RVR[ii]),float(LLAccSmallGausDiff[ii]))
  LLAccBigGausDiffDist.SetBinContent(LLAccBigGausDiffDist.FindBin(RVR[ii]),float(LLAccBigGausDiff[ii]))

AppliedCutsDir = Tfile.mkdir("LikelihoodAfterCuts")
SignScdDerDir = AppliedCutsDir.mkdir("SignSecondDerivative")
SignScdDerDir.cd()
LLPosScdDerDist01.Write(),   LLXSPosScdDerDist01.Write(),   LLAccPosScdDerDist01.Write()
LLNegScdDerDist01.Write(),   LLXSNegScdDerDist01.Write(),   LLAccNegScdDerDist01.Write()
LLPosScdDerDist02.Write(),   LLXSPosScdDerDist02.Write(),   LLAccPosScdDerDist02.Write()
LLNegScdDerDist02.Write(),   LLXSNegScdDerDist02.Write(),   LLAccNegScdDerDist02.Write()
LLPosScdDerDistBoth.Write(), LLXSPosScdDerDistBoth.Write(), LLAccPosScdDerDistBoth.Write()
#-- Save the variables separate of scd Der sign  --#
#--   --> Can give hint for good cut!            --#
#- Save the YPlusGausTest variables -#
YPlusGausTestPosScdDer.SetLineColor(3), YPlusGausTestXSPosScdDer.SetLineColor(3), YPlusGausTestAccPosScdDer.SetLineColor(3)
YPlusGausTestNegScdDer.SetLineColor(2), YPlusGausTestXSNegScdDer.SetLineColor(2), YPlusGausTestAccNegScdDer.SetLineColor(2)
YPlusGausTestCanvas = TCanvas('YPlusGausTestCanvas','Combined distribution of YPlusGausTest variable for events with positive and negative scd derivative (no normalisation)')
YPlusGausTestCanvas.cd(), YPlusGausTestPosScdDer.Draw(), YPlusGausTestNegScdDer.Draw("same"), YPlusGausTestCanvas.Write()
YPlusGausTestXSCanvas = TCanvas('YPlusGausTestXSCanvas','Combined distribution of YPlusGausTest variable for events with positive and negative scd derivative (XS normalisation)')
YPlusGausTestXSCanvas.cd(), YPlusGausTestXSNegScdDer.Draw(), YPlusGausTestXSPosScdDer.Draw("same"), YPlusGausTestXSCanvas.Write()
YPlusGausTestAccCanvas = TCanvas('YPlusGausTestAccCanvas','Combined distribution of YPlusGausTest variable for events with positive and negative scd derivative (Acc normalisation)')
YPlusGausTestAccCanvas.cd(), YPlusGausTestAccNegScdDer.Draw(), YPlusGausTestAccPosScdDer.Draw("same"), YPlusGausTestAccCanvas.Write()
#- Save the YPlus variables -#
YPlusPosScdDer.SetLineColor(3), YPlusXSPosScdDer.SetLineColor(3), YPlusAccPosScdDer.SetLineColor(3)
YPlusNegScdDer.SetLineColor(2), YPlusXSNegScdDer.SetLineColor(2), YPlusAccNegScdDer.SetLineColor(2)
YPlusCanvas = TCanvas('YPlusCanvas','Combined distribution of YPlus variable for events with positive and negative scd derivative (no normalisation)')
YPlusCanvas.cd(), YPlusPosScdDer.Draw(), YPlusNegScdDer.Draw("same"), YPlusCanvas.Write()
YPlusXSCanvas = TCanvas('YPlusXSCanvas','Combined distribution of YPlus variable for events with positive and negative scd derivative (XS normalisation)')
YPlusXSCanvas.cd(), YPlusXSNegScdDer.Draw(), YPlusXSPosScdDer.Draw("same"), YPlusXSCanvas.Write()
YPlusAccCanvas = TCanvas('YPlusAccCanvas','Combined distribution of YPlus variable for events with positive and negative scd derivative (Acc normalisation)')
YPlusAccCanvas.cd(), YPlusAccNegScdDer.Draw(), YPlusAccPosScdDer.Draw("same"), YPlusAccCanvas.Write()
#- Save the YPlusPlus variables -#
YPlusPlusPosScdDer.SetLineColor(3), YPlusPlusXSPosScdDer.SetLineColor(3), YPlusPlusAccPosScdDer.SetLineColor(3)
YPlusPlusNegScdDer.SetLineColor(2), YPlusPlusXSNegScdDer.SetLineColor(2), YPlusPlusAccNegScdDer.SetLineColor(2)
YPlusPlusCanvas = TCanvas('YPlusPlusCanvas','Combined distribution of YPlusPlus variable for events with positive and negative scd derivative (no normalisation)')
YPlusPlusCanvas.cd(), YPlusPlusPosScdDer.Draw(), YPlusPlusNegScdDer.Draw("same"), YPlusPlusCanvas.Write()
YPlusPlusXSCanvas = TCanvas('YPlusPlusXSCanvas','Combined distribution of YPlusPlus variable for events with positive and negative scd derivative (XS normalisation)')
YPlusPlusXSCanvas.cd(), YPlusPlusXSNegScdDer.Draw(), YPlusPlusXSPosScdDer.Draw("same"), YPlusPlusXSCanvas.Write()
YPlusPlusAccCanvas = TCanvas('YPlusPlusAccCanvas','Combined distribution of YPlusPlus variable for events with positive and negative scd derivative (Acc normalisation)')
YPlusPlusAccCanvas.cd(), YPlusPlusAccNegScdDer.Draw(), YPlusPlusAccPosScdDer.Draw("same"), YPlusPlusAccCanvas.Write()
#- Save the YMin variables -#
YMinPosScdDer.SetLineColor(3), YMinXSPosScdDer.SetLineColor(3), YMinAccPosScdDer.SetLineColor(3)
YMinNegScdDer.SetLineColor(2), YMinXSNegScdDer.SetLineColor(2), YMinAccNegScdDer.SetLineColor(2)
YMinCanvas = TCanvas('YMinCanvas','Combined distribution of YMin variable for events with positive and negative scd derivative (no normalisation)')
YMinCanvas.cd(), YMinPosScdDer.Draw(), YMinNegScdDer.Draw("same"), YMinCanvas.Write()
YMinXSCanvas = TCanvas('YMinXSCanvas','Combined distribution of YMin variable for events with positive and negative scd derivative (XS normalisation)')
YMinXSCanvas.cd(), YMinXSNegScdDer.Draw(), YMinXSPosScdDer.Draw("same"), YMinXSCanvas.Write()
YMinAccCanvas = TCanvas('YMinAccCanvas','Combined distribution of YMin variable for events with positive and negative scd derivative (Acc normalisation)')
YMinAccCanvas.cd(), YMinAccNegScdDer.Draw(), YMinAccPosScdDer.Draw("same"), YMinAccCanvas.Write()
#- Save the YMinMin variables -#
YMinMinPosScdDer.SetLineColor(3), YMinMinXSPosScdDer.SetLineColor(3), YMinMinAccPosScdDer.SetLineColor(3)
YMinMinNegScdDer.SetLineColor(2), YMinMinXSNegScdDer.SetLineColor(2), YMinMinAccNegScdDer.SetLineColor(2)
YMinMinCanvas = TCanvas('YMinMinCanvas','Combined distribution of YMinMin variable for events with positive and negative scd derivative (no normalisation)')
YMinMinCanvas.cd(), YMinMinPosScdDer.Draw(), YMinMinNegScdDer.Draw("same"), YMinMinCanvas.Write()
YMinMinXSCanvas = TCanvas('YMinMinXSCanvas','Combined distribution of YMinMin variable for events with positive and negative scd derivative (XS normalisation)')
YMinMinXSCanvas.cd(), YMinMinXSNegScdDer.Draw(), YMinMinXSPosScdDer.Draw("same"), YMinMinXSCanvas.Write()
YMinMinAccCanvas = TCanvas('YMinMinAccCanvas','Combined distribution of YMinMin variable for events with positive and negative scd derivative (Acc normalisation)')
YMinMinAccCanvas.cd(), YMinMinAccNegScdDer.Draw(), YMinMinAccPosScdDer.Draw("same"), YMinMinAccCanvas.Write()

GausDiffYPlusDir = AppliedCutsDir.mkdir("CounteringDiffFromGausYPlus")
GausDiffYPlusDir.cd()
LLSmallGausDiffDist.Write(), LLXSSmallGausDiffDist.Write(), LLAccSmallGausDiffDist.Write()
LLBigGausDiffDist.Write(),   LLXSBigGausDiffDist.Write(),   LLAccBigGausDiffDist.Write()

