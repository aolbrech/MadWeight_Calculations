############################################
#                                          #
#  RemoveLowPtEvents.py script             #
#  --> Has a two-fold goal:                #
#      1) Removing events below pt-value   #
#      2) Include cos theta weight in file #
#                                          #
############################################
#! python
import os
import sys
from math import log,pow

#--------------------------------#
#   Get input from command line  #
#--------------------------------#
if len(sys.argv) <= 3:
  print "Need to give the directory of interest, the anomalous couplings coefficient, the desired pT-cut and the value of the configuration in command line !"
  print "--> Syntax example: python RemoveLowPtEvents.py Events/blabla RVR 15/0/(Also)MET Pos05"
  sys.exit()
whichDir, AnomCoef, ptCut, VarConfig = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

#----------------------------------#
#   Set the correct output files   #
#----------------------------------#
#Set the correct names!
PTCut = ""
if ptCut == "0": PTCut, OutFileTitle = "NoPtCut", "NoPtCut"
elif ptCut == "AlsoMET" or ptCut == "MET": PTCut, OutFileTitle = "PtCutsApplied_AlsoOnMET","PtCutsApplied_AlsoOnMET"
else:  PTCut, OutFileTitle = "PtCut"+str(ptCut), 'NoLowPt_PtCut'+str(ptCut)

outfile = file(os.path.join(whichDir+'weights_'+OutFileTitle+'.out'),'w')
outfileCorrLn = file(os.path.join(whichDir+'weights_'+OutFileTitle+'_ApplyCosThetaReweighting.out'),'w')
CosThetaCorrFile = file(os.path.join('/user/aolbrech/AnomalousCouplings/MadAnalysis_v112/RVRAcceptance/SampleAnalyzer/CosThetaReweighting_'+str(AnomCoef)+'/'+str(VarConfig)+'_'+str(AnomCoef)+'/CosThetaWeight_'+str(VarConfig)+'_'+str(AnomCoef)+'_'+PTCut+'.txt'),'r')

print "Getting cos theta info from : ",CosThetaCorrFile.name

#----------------------------------#
#   Get the correct output files   #
#----------------------------------#
list_dir = os.listdir(whichDir)
WeightsFileArray, weightsFileCounter = [], 0
for file in list_dir:
  if file.endswith("weights.out") and file.startswith("weights.out"): # eg: '.txt'
    weightsFileCounter += 1
    WeightsFileArray.append(file)
  if file.endswith(".lhco"): LHCOFile = open(os.path.join(whichDir+''+file),'r')  #Only 1 lhco file should exist!

if int(weightsFileCounter) == 0:
  print "No weights file found in this directory !"
  sys.exit()
elif int(weightsFileCounter) == 1: WeightsFile = open(os.path.join(whichDir+''+WeightsFileArray[0]),'r')
elif int(weightsFileCounter) > 1:
  for ii in range(len(WeightsFileArray)):
    print " ",ii," ) ",WeightsFileArray[ii]
  fileNr = raw_input('Choose the number of the file of interest! : ')
  WeightsFile = open(os.path.join(whichDir+''+WeightsFileArray[int(fileNr)]),'r')

print "Will be using file : ",WeightsFile.name," to delete events from !"
print "Will be using lhco file : ",LHCOFile.name

#-------------------------------#
#   Perform the event deletion  #
#-------------------------------#
EventsToDelete = []
for lhcoLine in LHCOFile:
  lhcoWord = lhcoLine.split()
  if str(lhcoWord[0]) != "#":                                        #Only interested in lines starting with a number
    if len(lhcoWord) == 3: EvtNr, SurvivedCut = lhcoWord[1], True    #Initialize for each event!
    if len(lhcoWord) == 11:
      if SurvivedCut == True and str(lhcoWord[1]) != "6" and PTCut != "PtCutsApplied_AlsoOnMET" and float(lhcoWord[4]) < float(ptCut):   #Currently pT-cut is not applied to MET ...
        SurvivedCut = False
        EventsToDelete.append(int(EvtNr))
print "Events with a particle with pT < ",str(ptCut)," GeV (excluding MET) :",len(EventsToDelete)

#------------------------------#
#   Save the cos theta weight  #
#------------------------------#
#First get the info from the output file of MadAnalysis calculation!
CosThetaCorr = []
for corrLine in CosThetaCorrFile: CosThetaCorr.append(float(corrLine.split()[1]))

for line in WeightsFile:
  word = line.split()
  if str(word[0]) != "#":
    if not int(word[0]) in EventsToDelete:
      outfile.write(line)
      outfileCorrLn.write(line.replace(word[4],str(CosThetaCorr[int(word[0])])+'  '+word[4]))
  else: outfile.write(line), outfileCorrLn.write(line)   #Write the first two lines containing text!

outfile.close()
outfileCorrLn.close()

#Delete the outfile (~ output file with events removed) when no events have actually been deleted!
if len(EventsToDelete) == 0: 
  print "Removing file",outfile.name," since no events have been deleted so it is identical to weights.out!",
  os.remove(outfile.name)
