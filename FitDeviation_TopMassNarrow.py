#! python
import os
import sys
from math import log
#import ROOT
from array import array
from ROOT import TH1F,TFile,TCanvas,TLegend,gStyle,TDirectory,gROOT,TGraph
gROOT.SetBatch(True)                       #Don't print the histograms (necessary for fit)

#Directory of interest:
if len(sys.argv) <= 1:
  print "Need to give the directory of interest in command line !"
  sys.exit()
whichDir = sys.argv[1] 
print " Interested in directory : ",whichDir

#Information about the scanned RVR values and the corresponding cross-section
MTopValues = ["m_{top} = 171",  "m_{top} = 172",  "m_{top} = 173",  "m_{top} = 174",  "m_{top} = 175"]
MTop =       array('d',[171,              172,              173,              174,              175            ])
MGXS =       array('d',[10.70485,         10.8257,          10.96469,         11.08428,         11.22448       ])
MGXSe =      array('d',[0.00814303595657, 0.00878899028501, 0.00816801717126, 0.00904797742371, 0.008653800078 ])
Acceptance = array('d',[0.21460,          0.21735,          0.21290,          0.21752,          0.22185        ])

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

#ROOT file where all the information will be stored:
Tfile = TFile(os.path.join(whichDir+"FitDeviation_RVR.root"),'recreate')
gStyle.SetOptStat(0)
YPlus    = TH1F('YPlus',   'Deviation from parabolic fit for m(top) = 175 (no normalisation)' ,150,-0.8,0.8)
YPlusXS  = TH1F('YPlusXS', 'Deviation from parabolic fit for m(top) = 175 (XS normalisation)' ,150,-0.8,0.8)
YPlusAcc = TH1F('YPlusAcc','Deviation from parabolic fit for m(top) = 175 (Acc normalisation)',150,-0.8,0.8)
YMin    = TH1F('YMin',   'Deviation from parabolic fit for m(top) = 171 (no normalisation)' ,150,-0.8,0.8)
YMinXS  = TH1F('YMinXS', 'Deviation from parabolic fit for m(top) = 171 (XS normalisation)' ,150,-0.8,0.8)
YMinAcc = TH1F('YMinAcc','Deviation from parabolic fit for m(top) = 171 (Acc normalisation)',150,-0.8,0.8)
LnVariation    = TH1F('LnVariation',   'Variation of obtained likelihood values (no normalisation)' ,150,0,2)
LnVariationXS  = TH1F('LnVariatonXS',  'Variation of obtained likelihood values (XS normalisation)' ,150,0,2)
LnVariationAcc = TH1F('LnVariationAcc','Variation of obtained likelihood values (Acc normalisation)',150,0,2)

FitComp = Tfile.mkdir("FitComparison")
LnLogDir = Tfile.mkdir("LnLogDistr")
LnLogXSDir = Tfile.mkdir("LnLogXSDistr")
LnLogAccDir = Tfile.mkdir("LnLogAccDistr")
LnLogFitCanvas = TCanvas('name','title')
LnLogFitCanvasName = "LnLogFitCanvas_Evt"
LnLogFitCanvasTitle = "Comparing fit function from ROOT fit and algebraic function for event "

#Create the arrays where the likelihood values will be stored
NrConfigs = 5
LnLog = []
LnLogXS = []
LnLogAcc = []
NrConsEvents = []
for ii in range(NrConfigs):
  LnLog.append(0)
  LnLogXS.append(0)
  LnLogAcc.append(0)
  NrConsEvents.append(0)

Nevents = 1000
#Loop over all lines in weights file:
for WeightLine in WeightsFile:
  WeightWord = WeightLine.split()
  #Only interested in files starting with a number
  if str(WeightWord[0]) != "#" and str(WeightWord[3]) != "0.0" :
    for iEvt in range(Nevents):
      if str(WeightWord[0]) == str(iEvt+1):                                                               #Look at one single event!
        if str(WeightWord[1]) == "1":
          print "Looking at event : ",iEvt+1
          LogMin = [-log(float(WeightWord[3])),-log(float(WeightWord[3]))+log(MGXS[int(WeightWord[1])-1]), -log(float(WeightWord[3]))+log(MGXS[int(WeightWord[1])-1])+log(Acceptance[int(WeightWord[1])-1])]
          LogMax = [0,0,0]
          cHat, bHat, aHat = [0,0,0], [0,0,0], [0,0,0]
          LnLogDist = TH1F("LnLog_Evt"+str(int(iEvt)+1),"LnLog distribution for event "+str(int(iEvt)+1),11,-0.55,0.55)
          LnLogDist.SetMarkerStyle(20)
          LnLogDist.SetLineColor(1)
          LnLogDist.SetMarkerColor(1)
          LnLogDist.SetMarkerSize(1.2)
          LnLogXSDist = TH1F("LnLogXS_Evt"+str(int(iEvt)+1),"LnLogXS distribution for event "+str(int(iEvt)+1),11,-0.55,0.55)
          LnLogXSDist.SetMarkerStyle(21)
          LnLogXSDist.SetLineColor(3)
          LnLogXSDist.SetMarkerColor(3)
          LnLogXSDist.SetMarkerSize(1.2)
          LnLogAccDist = TH1F("LnLogAcc_Evt"+str(int(iEvt)+1),"LnLogAcc distribution for event "+str(int(iEvt)+1),11,-0.55,0.55)
          LnLogAccDist.SetMarkerStyle(22)
          LnLogAccDist.SetLineColor(4)
          LnLogAccDist.SetMarkerColor(4)
          LnLogAccDist.SetMarkerSize(1.2)
        LnLog[int(WeightWord[1])-1] = -log(float(WeightWord[3]))
        LnLogXS[int(WeightWord[1])-1] = -log(float(WeightWord[3])) + log(MGXS[int(WeightWord[1])-1])
        LnLogAcc[int(WeightWord[1])-1] = -log(float(WeightWord[3])) + log(MGXS[int(WeightWord[1])-1]) + log(Acceptance[int(WeightWord[1])-1])
        NrConsEvents[int(WeightWord[1])-1] +=1
        #Determine max and min for LnLog:
        if LnLog[int(WeightWord[1])-1] >= LogMax[0]:
          LogMax[0] = LnLog[int(WeightWord[1])-1]
        if LnLog[int(WeightWord[1])-1] < LogMin[0]:
          LogMin[0] = LnLog[int(WeightWord[1])-1]
        #Determine max and min for LnLogXS:
        if LnLogXS[int(WeightWord[1])-1] >= LogMax[1]:
          LogMax[1] = LnLogXS[int(WeightWord[1])-1]
        if LnLogXS[int(WeightWord[1])-1] < LogMin[1]:
          LogMin[1] = LnLogXS[int(WeightWord[1])-1]
        #Determine max and min for LnLogAcc:
        if LnLogAcc[int(WeightWord[1])-1] >= LogMax[2]:
          LogMax[2] = LnLogAcc[int(WeightWord[1])-1]
        if LnLogAcc[int(WeightWord[1])-1] < LogMin[2]:
          LogMin[2] = LnLogAcc[int(WeightWord[1])-1]
        #---  Fill the LnLog histograms for each event  ---#
        #if str(WeightWord[1]) != "1" or str(WeightWord[1]) != "9":
        LnLogDist.SetBinContent(LnLogDist.FindBin(RVR[int(WeightWord[1])-1]), LnLog[int(WeightWord[1])-1])
        LnLogXSDist.SetBinContent(LnLogXSDist.FindBin(RVR[int(WeightWord[1])-1]), LnLogXS[int(WeightWord[1])-1])
        LnLogAccDist.SetBinContent(LnLogAccDist.FindBin(RVR[int(WeightWord[1])-1]), LnLogAcc[int(WeightWord[1])-1])
        #---  Only perform the fit after all 9 configurations are considered!  ---#
        if str(WeightWord[1]) == str(NrConfigs):
          #---  Calculate the fit parameters (polynomial)  ---#
          cHat = [LnLog[2],                          LnLogXS[4],                              LnLogAcc[4]                               ]
          bHat = [5*(LnLog[5]-LnLog[3]),             5*(LnLogXS[5]-LnLogXS[3]),               5*(LnLogAcc[5]-LnLogAcc[3])               ]
          aHat = [50*(LnLog[3]+LnLog[5]-2*LnLog[4]), 50*(LnLogXS[3]+LnLogXS[5]-2*LnLogXS[4]), 50*(LnLogAcc[3]+LnLogAcc[5]-2*LnLogAcc[4])] 
          LnLogFunction, LnLogXSFunction, LnLogAccFunction = array('d'), array('d'), array('d')
          for ii in range(9):
            LnLogFunction.append( aHat[0]*float(RVR[ii])*float(RVR[ii])+bHat[0]*float(RVR[ii])+cHat[0])
            LnLogXSFunction.append( aHat[1]*float(RVR[ii])*float(RVR[ii])+bHat[1]*float(RVR[ii])+cHat[1])
            LnLogAccFunction.append( aHat[2]*float(RVR[ii])*float(RVR[ii])+bHat[2]*float(RVR[ii])+cHat[2])
          LnLogGraph = TGraph(9, RVR, LnLogFunction)
          LnLogGraph.SetTitle('Comparison between ROOT fit and algebraic method')
          LnLogGraph.GetXaxis().SetTitle('#Delta ln(Likelihood)')
          LnLogGraph.GetYaxis().SetTitle('Number of events')
          LnLogGraph.SetMarkerColor(6)
          LnLogGraph.SetLineColor(6)
          LnLogGraph.SetMarkerSize(1.2)
          #LnLogXSGraph = TGraph(9, RVR, LnLogXSFunction)
          #LnLogXSGraph.SetMarkerColor(6)
          #LnLogXSGraph.SetLineColor(6)
          #LnLogXSGraph.SetMarkerSize(1.2)
          #LnLogAccGraph = TGraph(9, RVR, LnLogAccFunction)
          #LnLogAccGraph.SetMarkerColor(6)
          #LnLogAccGraph.SetLineColor(6)
          #LnLogAccGraph.SetMarkerSize(1.2)
          #---  Compare with TF1 fit from ROOT  ---#
          LnLogDist.Fit("pol2","Q","",-0.35, 0.35)
          LnLogFitResult = LnLogDist.GetFunction("pol2")
          #LnLogXSDist.Fit("pol2","Q")
          #LnLogXSFitResult = LnLogXSDist.GetFunction("pol2")
          #LnLogAccDist.Fit("pol2","Q")
          #LnLogAccFitResult = LnLogAccDist.GetFunction("pol2")
          #cHatFit = [LnLogFitResult.GetParameter(0), LnLogXSFitResult.GetParameter(0), LnLogAccFitResult.GetParameter(0)]
          #bHatFit = [LnLogFitResult.GetParameter(1), LnLogXSFitResult.GetParameter(1), LnLogAccFitResult.GetParameter(1)]
          #aHatFit = [LnLogFitResult.GetParameter(2), LnLogXSFitResult.GetParameter(2), LnLogAccFitResult.GetParameter(2)]
          #-- Plot comparision --#
          LnLogFitCanvas.SetName(LnLogFitCanvasName+str(int(iEvt)+1))
          LnLogFitCanvas.SetTitle(LnLogFitCanvasTitle+str(int(iEvt)+1))
          legendLnLog = TLegend(0.75,0.7,0.95,0.95)
          LnLogFitCanvas.cd()
          legendLnLog.AddEntry(LnLogFitResult,'Fit result using ROOT Fit',"l")
          legendLnLog.AddEntry(LnLogGraph,'Fit result using algebraic function',"p")
          legendLnLog.AddEntry(LnLogDist,'Obtained LnLog values',"p")
          LnLogGraph.Draw("AC*")
          LnLogFitResult.Draw("same")
          LnLogDist.Draw("samep")
          FitComp.cd()
          LnLogFitCanvas.Write()	
          #---  Calculate the positive and negative deviation  ---#
          FitOutputPlus = [0,0,0]
          FitOutputMin  = [0,0,0]
          for ii in range(3):
            FitOutputPlus[ii] = (aHat[ii]*int(RVR[6])*int(RVR[6])+bHat[ii]*int(RVR[6])+cHat[ii])
            FitOutputMin[ii]  = (aHat[ii]*int(RVR[2])*int(RVR[2])+bHat[ii]*int(RVR[2])+cHat[ii])
          yPlus = [LnLog[6] - FitOutputPlus[0], LnLogXS[6] - FitOutputPlus[1], LnLogAcc[6] - FitOutputPlus[2]]
          yMin  = [LnLog[2] - FitOutputMin[0],  LnLogXS[2] - FitOutputMin[1],  LnLogAcc[2] - FitOutputMin[2] ]
          YPlus.Fill(yPlus[0])
          YPlusXS.Fill(yPlus[1])
          YPlusAcc.Fill(yPlus[2])
          YMin.Fill(yMin[0])
          YMinXS.Fill(yMin[1])
          YMinAcc.Fill(yMin[2])
          LnVariation.Fill(LogMax[0] - LogMin[0])
          LnVariationXS.Fill(LogMax[1] - LogMin[1])
          LnVariationAcc.Fill(LogMax[2] - LogMin[2])
          #---  Save the LnLog distributions (event-per-event) and the Y-deviations in histograms  ---#
          LnLogDir.cd()
          LnLogDist.Write()
          LnLogXSDir.cd()
          LnLogXSDist.Write()
          LnLogAccDir.cd()
          LnLogAccDist.Write()
Tfile.cd()
YPlus.Write()
YPlusXS.Write()
YPlusAcc.Write()
YMin.Write()
YMinXS.Write()
YMinAcc.Write()
LnVariation.Write()
LnVariationXS.Write()
LnVariationAcc.Write()
