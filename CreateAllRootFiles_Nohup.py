#! python
import os
import sys

#Get the Events directory, the considered range and the correction from the input line
whichDir       = sys.argv[1]
whichVar       = sys.argv[2]
nrEvts         = sys.argv[3]
whichAnalysis  = sys.argv[4]
applyAccNorm   = sys.argv[5]
applyCosThCorr = sys.argv[6]
cosThetaInLn   = sys.argv[7]
whichRange     = sys.argv[8]
whichPtCut     = sys.argv[9]
if len(sys.argv) != 10:
  print " Need to give the Events directory, the considered range and the correction in the input line"
  #                                   0                          1         2          3         4                     5                         6                          7          8     9
  print " Correct input : python CreateAllRootFiles_Nohup.py Events_xx RVR/RgR/MTop #evts doublePolFitMacro.C ApplyAccNormForCuts(y/n) ApplyCosThetaReweighting(y/n) CosThInLn?(y/n) Wide 0/15/30 "
  sys.exit()

if str(applyCosThCorr) == "y" or str(applyCosThCorr) == "Y" or str(applyCosThCorr) == "yes":
  whichCorr = "ApplyCosThetaReweighting"
  if str(cosThetaInLn) == "y" or str(cosThetaInLn) == "Y" or str(cosThetaInLn) == "yes":
    whichCorr = "CosThetaMultipliedWithWeight"    
elif str(applyCosThCorr) == "n" or str(applyCosThCorr) == "N" or str(applyCosThCorr) == "no":
  whichCorr = "NoCosTheta"
else:
  print "Need to specify whether the costheta correction needs to be applied (y/n)"
  sys.exit()

list_dir = []
list_dir = os.listdir(whichDir)
for VarDir in list_dir:
  if float(whichPtCut) == 0:
    os.system('python FitDeviationScript.py '+str(whichDir)+''+str(VarDir)+'/ '+str(whichVar)+' '+str(nrEvts)+' y n '+str(whichAnalysis)+' '+str(whichRange)+' '+str(whichDir)+''+str(VarDir)+'/weights.out')
  else:
    #First need to make sure the correct file exist, and otherwise create it
    #--> Used file depends on whichCorr is given!
    if str(whichCorr) == "ApplyCosThetaReweighting" or str(whichCorr) == "CosThetaMultipliedWithWeight":
      if not (os.path.exists(os.path.join(whichDir+""+VarDir+"/weights_NoLowPt_Cut"+str(whichPtCut)+"_"+str(whichCorr)+".out"))):
        os.system('python RemoveLowPtEvents.py '+str(whichDir)+''+str(VarDir)+'/ '+str(whichPtCut))
      os.system('python FitDeviationScript.py '+str(whichDir)+''+str(VarDir)+'/ '+str(whichVar)+' '+str(nrEvts)+' y n '+str(whichAnalysis)+' '+str(whichRange)+' '+str(whichDir)+''+str(VarDir)+'/weights_NoLowPt_Cut'+str(whichPtCut)+'_'+str(whichCorr)+'.out'+' '+str(applyAccNorm))
    elif str(whichCorr) == "NoCosTheta":
      if not (os.path.exists(os.path.join(whichDir+""+VarDir+"/weights_NoLowPt_Cut"+str(whichPtCut)+".out"))):
        os.system('python RemoveLowPtEvents.py '+str(whichDir)+''+str(VarDir)+'/ '+str(whichPtCut))
      os.system('python FitDeviationScript.py '+str(whichDir)+''+str(VarDir)+'/ '+str(whichVar)+' '+str(nrEvts)+' y n '+str(whichAnalysis)+' '+str(whichRange)+' '+str(whichDir)+''+str(VarDir)+'/weights_NoLowPt_Cut'+str(whichPtCut)+'.out'+' '+str(applyAccNorm))
