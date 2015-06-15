#! python
import os
import sys
from math import log

#Directory of interest:
if len(sys.argv) <= 2:
  print "Need to give the directory of interest and the desired pT-cut in command line !"
  sys.exit()
if len(sys.argv) == 1:
  print "Need to give the directory of interest AND the desired pT-cut in command line !"
  sys.exit()
whichDir = sys.argv[1] 
ptCut = sys.argv[2]
print " Interested in directory : ",whichDir

#Output file of interest:
outfile = file(os.path.join(whichDir+'weights_NoLowPt_Cut'+str(ptCut)+'.out'),'w')
outfileCorr = file(os.path.join(whichDir+'weights_NoLowPt_Cut'+str(ptCut)+'_CosThetaCorrection.out'),'w')
outfileCorrLn = file(os.path.join(whichDir+'weights_NoLowPt_Cut'+str(ptCut)+'_CosThetaCorr_LnWeight.out'),'w')

VarConfig = ""
if   whichDir.find("Pos05") < len(whichDir) and whichDir.find("Pos05") >= 0: VarConfig = "Pos05"
elif whichDir.find("Pos03") < len(whichDir) and whichDir.find("Pos03") >= 0: VarConfig = "Pos03"
elif whichDir.find("Pos02") < len(whichDir) and whichDir.find("Pos02") >= 0: VarConfig = "Pos02"
elif whichDir.find("Pos01") < len(whichDir) and whichDir.find("Pos01") >= 0: VarConfig = "Pos01"
elif whichDir.find("SampleSM") < len(whichDir) and whichDir.find("SampleSM") >= 0: VarConfig = "SM" 
elif whichDir.find("Neg01") < len(whichDir) and whichDir.find("Neg01") >= 0: VarConfig = "Neg01"
elif whichDir.find("Neg02") < len(whichDir) and whichDir.find("Neg02") >= 0: VarConfig = "Neg02"
elif whichDir.find("Neg03") < len(whichDir) and whichDir.find("Neg03") >= 0: VarConfig = "Neg03"
elif whichDir.find("Neg05") < len(whichDir) and whichDir.find("Neg05") >= 0: VarConfig = "Neg05"
CosThetaCorrFile = file(os.path.join('/user/aolbrech/AnomalousCouplings/MadAnalysis_v112/RVRAcceptance/SampleAnalyzer/CosThetaReweighting/'+str(VarConfig)+'/CosThetaWeight_'+str(VarConfig)+'_PtCut'+str(ptCut)+'.txt'),'r')
outputTest = file(os.path.join(whichDir+'outputTest.txt'),'w')

#Input file of interest:
list_dir = []
list_dir = os.listdir(whichDir)
WeightsFileArray = []
weightsFileCounter = 0
for file in list_dir:
  if file.endswith("weights.out") and file.startswith("weights.out"): # eg: '.txt'
    weightsFileCounter += 1
    WeightsFileArray.append(file)
  if file.endswith(".lhco"):
    LHCOFile = open(os.path.join(whichDir+''+file),'r')

if int(weightsFileCounter) == 1:
  WeightsFile = open(os.path.join(whichDir+''+WeightsFileArray[0]),'r')
  secondWeightsFile = open(os.path.join(whichDir+''+WeightsFileArray[0]),'r')
elif int(weightsFileCounter) == 0:
  print "No weights file found in this directory !"
  sys.exit()
elif int(weightsFileCounter) > 1:
  for ii in range(len(WeightsFileArray)):
    print " ",ii," ) ",WeightsFileArray[ii]
  fileNr = raw_input('Choose the number of the file of interest! : ')
  WeightsFile = open(os.path.join(whichDir+''+WeightsFileArray[int(fileNr)]),'r')
#  secondWeightsFile = open(os.path.join(whichDir+''+WeightsFileArray[int(fileNr)]),'r')

print "Will be using file : ",WeightsFile," to delete events from !"
print "Will be using lhco file : ",LHCOFile

EventsToDelete = []
#Loop over all lines in weights file:
for lhcoLine in LHCOFile:
  lhcoWord = lhcoLine.split()
  #Only interested in files starting with a number
  if str(lhcoWord[0]) != "#":
    if len(lhcoWord) == 3:
      if str(lhcoWord[1]) != "0":
        if SurvivedCut == True:
          outputTest.write(EvtNr+' 1 \n')
        elif SurvivedCut == False:
          outputTest.write(EvtNr+' 0 \n')
      EvtNr = lhcoWord[1]
      SurvivedCut = True
    if len(lhcoWord) == 11:
      if SurvivedCut == True and str(lhcoWord[1]) != "6" and float(lhcoWord[4]) < float(ptCut):
        SurvivedCut = False
        EventsToDelete.append(int(EvtNr))

outputTest.close()
print "Events with a particle with pT < ",str(ptCut)," GeV (excluding MET) :",len(EventsToDelete)

CosThetaCorr = []
for corrLine in CosThetaCorrFile:
  corrWord = corrLine.split()
  CosThetaCorr.append(float(corrWord[1]))


for line in WeightsFile:
  word = line.split()
  if str(word[0]) != "#":
    if not int(word[0]) in EventsToDelete:
      outfile.write(line)
      #Get the correct correction factor for the considered event!
      outfileCorr.write(line.replace(word[3],str(float(word[3])*CosThetaCorr[int(word[0])])))
      outfileCorrLn.write(line.replace(word[3],str(-log(float(word[3]))*CosThetaCorr[int(word[0])])))
  else:
    outfile.write(line)   #Write the first two lines containing text!
    outfileCorr.write(line)
    outfileCorrLn.write(line)

outfile.close()
outfileCorr.close()
