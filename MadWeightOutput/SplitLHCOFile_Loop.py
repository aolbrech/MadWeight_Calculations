#!python
import linecache
import sys

EventMatchingFile = file('EventNumberInformation.txt','r')
MatchedEvents = [file('MatchedEventInfo_PositiveMuon.lhco','w'), file('MatchedEventInfo_NegativeMuon.lhco','w'), file('MatchedEventInfo_PositiveElectron.lhco','w'), file('MatchedEventInfo_NegativeElectron.lhco','w')]

GeneratorFileName = ['TTbarLHCO_PositiveMuon.lhco', 'TTbarLHCO_NegativeMuon.lhco','TTbarLHCO_PositiveElectron.lhco','TTbarLHCO_NegativeElectron.lhco']
GeneratorFile = [file('TTbarLHCO_PositiveMuon.lhco','r'), file('TTbarLHCO_NegativeMuon.lhco','r'), file('TTbarLHCO_PositiveElectron.lhco','r'), file('TTbarLHCO_NegativeElectron.lhco','r')]

MatchedGeneratorName = ['MatchedTTbarLHCO_PositiveMuon.lhco', 'MatchedTTbarLHCO_NegativeMuon.lhco', 'MatchedTTbarLHCO_PositiveElectron.lhco', 'MatchedTTbarLHCO_NegativeElectron.lhco']
MatchedGeneratorFile = [file('MatchedTTbarLHCO_PositiveMuon.lhco','w'), file('MatchedTTbarLHCO_NegativeMuon.lhco','w'), file('MatchedTTbarLHCO_PositiveElectron.lhco','w'), file('MatchedTTbarLHCO_NegativeElectron.lhco','w')]

ReconstructedName = ['TTbarSemiLepton_Reco_PositiveMuon.lhco', 'TTbarSemiLepton_Reco_NegativeMuon.lhco', 'TTbarSemiLepton_Reco_PositiveElectron.lhco', 'TTbarSemiLepton_Reco_NegativeElectron.lhco']
ReconstructedFile = [file('TTbarSemiLepton_Reco_PositiveMuon.lhco','r'), file('TTbarSemiLepton_Reco_NegativeMuon.lhco','r'), file('TTbarSemiLepton_Reco_PositiveElectron.lhco','r'), file('TTbarSemiLepton_Reco_NegativeElectron.lhco','r')]

MadWeightStartLine = "#</MGPGSCard> \n "
MadWeightSecondStartLine = "  #  typ      eta      phi       pt   jmas  ntrk  btag   had/em  dummy  dummy \n"


LeptonTypeString = ["Positive Muon","Negative Muon","Positive Electron","Negative Electron"]
if len(sys.argv) != 3:
	print ' Used syntax is : python SplitLHCOFile.py "LeptonType" "NumberEvents" "PerformMatching" "PerformSplitting'
	print ' Need to specifiy which lepton Type has to be considered for the splitting of the .lhco file and the number of events which should be used'
	print ' The used numbering is the following:'
	print '		-> 0 - Positive Muon '
	print '		-> 1 - Negative Muon '
	print '		-> 2 - Positive Electron '
	print '		-> 3 - Negative Electron '
	exit
elif len(sys.argv) == 3:
	LeptonType = int(sys.argv[1])
	print ' Splitting will be done for the ', LeptonTypeString[LeptonType], ' case '
	NumberEvents = int(sys.argv[2])
	print ' LHCO files with ', NumberEvents, ' events will be created'

	MatchedGeneratorFile[LeptonType].write(linecache.getline(GeneratorFileName[LeptonType],1))             #First line containing specific madweight output
	MatchedGeneratorFile[LeptonType].write(linecache.getline(GeneratorFileName[LeptonType],2))             #Second line containing specific madweight output

	#Obtain the reconstructed and generator event number which need to be matched!
	for MatchingLine in EventMatchingFile:
	    MatchingWords = MatchingLine.split()
	    NumberWords = len(MatchingWords)    #Only need to perform the matching if the length is equal to 8! (this means that reconstruction is complete)
	    				#--> Will need to be updated when neutrinos are not reconstructed separately!!!
	    if NumberWords == 8:
		selectedEvent = 0
		if LeptonType == 0 and MatchingWords[1] == "1" and MatchingWords[2] == "0" and MatchingWords[3] == "0" and MatchingWords[4] == "0":   #Positive muon
			GeneratorEventNumber = MatchingWords[5]
			ReconstructedEventNumber = MatchingWords[7]
                        MatchedEvents[LeptonType].write(MatchingLine)

                        #Create new generator file
                        MatchedLine = linecache.getline(GeneratorFileName[LeptonType],(int(GeneratorEventNumber)-1)*7+3)
                        MatchedWords = MatchedLine.split()
                        MatchedGeneratorFile[LeptonType].write(MatchedLine.replace(MatchedWords[1],ReconstructedEventNumber)) #Line with the event number
                        for lineNr in range(6):
                                MatchedGeneratorFile[LeptonType].write(linecache.getline(GeneratorFileName[LeptonType],(int(GeneratorEventNumber)-1)*7+3+1+lineNr)) #Lines with the event kinematics
                                #Additional +1 since the iterator starts at 0 and goes up to 5!

		elif LeptonType == 1 and MatchingWords[1] == "0" and MatchingWords[2] == "1" and MatchingWords[3] == "0" and MatchingWords[4] == "0": #Negative muon
			GeneratorEventNumber = MatchingWords[5]
			ReconstructedEventNumber = MatchingWords[7]
                        MatchedEvents[LeptonType].write(MatchingLine)

                        #Create new generator file
                        MatchedLine = linecache.getline(GeneratorFileName[LeptonType],(int(GeneratorEventNumber)-1)*7+3)
                        MatchedWords = MatchedLine.split()
                        MatchedGeneratorFile[LeptonType].write(MatchedLine.replace(MatchedWords[1],ReconstructedEventNumber)) #Line with the event number
                        for lineNr in range(6):
                                MatchedGeneratorFile[LeptonType].write(linecache.getline(GeneratorFileName[LeptonType],(int(GeneratorEventNumber)-1)*7+3+1+lineNr)) #Lines with the event kinematics
                                #Additional +1 since the iterator starts at 0 and goes up to 5!

		elif LeptonType == 2 and MatchingWords[1] == "0" and MatchingWords[2] == "0" and MatchingWords[3] == "1" and MatchingWords[4] == "0": #Positive electron
			GeneratorEventNumber = MatchingWords[5]
			ReconstructedEventNumber = MatchingWords[7]
			MatchedEvents[LeptonType].write(MatchingLine)
	
	                #Create new generator file
			MatchedLine = linecache.getline(GeneratorFileName[LeptonType],(int(GeneratorEventNumber)-1)*7+3)
			MatchedWords = MatchedLine.split()
			MatchedGeneratorFile[LeptonType].write(MatchedLine.replace(MatchedWords[1],ReconstructedEventNumber)) #Line with the event number
			for lineNr in range(6):
	                        MatchedGeneratorFile[LeptonType].write(linecache.getline(GeneratorFileName[LeptonType],(int(GeneratorEventNumber)-1)*7+3+1+lineNr)) #Lines with the event kinematics
				#Additional +1 since the iterator starts at 0 and goes up to 5!

		elif LeptonType == 3 and MatchingWords[1] == "0" and MatchingWords[2] == "0" and MatchingWords[3] == "0" and MatchingWords[4] == "0": #Negative electron
			GeneratorEventNumber = MatchingWords[5]
			ReconstructedEventNumber = MatchingWords[7]
                        MatchedEvents[LeptonType].write(MatchingLine)

                        #Create new generator file
                        MatchedLine = linecache.getline(GeneratorFileName[LeptonType],(int(GeneratorEventNumber)-1)*7+3)
                        MatchedWords = MatchedLine.split()
                        MatchedGeneratorFile[LeptonType].write(MatchedLine.replace(MatchedWords[1],ReconstructedEventNumber)) #Line with the event number
                        for lineNr in range(6):
                                MatchedGeneratorFile[LeptonType].write(linecache.getline(GeneratorFileName[LeptonType],(int(GeneratorEventNumber)-1)*7+3+1+lineNr)) #Lines with the event kinematics
                                #Additional +1 since the iterator starts at 0 and goes up to 5!

	#Split the .lhco files to the desired number of events
	print '\n Which files do you want to split? '
	print '      - GeneratorFile        ', (GeneratorFileName)		
	print '      - MatchedGeneratorFile ', (MatchedGeneratorName)
	print '      - ReconstructedFile    ', (ReconstructedName)
	filename = input()
	
	OutputLeptonType = ["PositiveMuon","NegativeMuon","PositiveElectron","NegativeElectron"]
	OutputEvent = ["Generator","MatchedGenerator","Reconstructed"]
	if filename == GeneratorFile:
		EventType = 0
	elif filename == MatchedGeneratorFile:
		EventType = 1
	elif filename == ReconstructedFile:
		EventType = 2

	for fileIterator in range(len(filename)):
		LineCounter = 0
		splitCounter = 0
		splitOutputFileName = str('SplitFiles/Splitted_'+str(NumberEvents)+'Evts_'+OutputEvent[EventType]+'_'+OutputLeptonType[fileIterator]+'_'+str(splitCounter)+'.lhco')
		splitOutputFile = file(str(splitOutputFileName),'w')
		for splitLine in filename[fileIterator]:
			LineCounter = LineCounter+1
			#Write to output file	
			splitOutputFile.write(splitLine)

			if LineCounter == NumberEvents*7*(splitCounter+1)+2:    #lineCounter starts at 1 (as linecache) and splitCounter starts at 0
				splitCounter = splitCounter+1
 		                splitOutputFileName = str("SplitFiles/Splitted_"+str(NumberEvents)+"Evts_"+OutputEvent[EventType]+"_"+OutputLeptonType[fileIterator]+"_"+str(splitCounter)+".lhco")
				splitOutputFile = file(str(splitOutputFileName),'w')
				splitOutputFile.write(MadWeightStartLine)
				splitOutputFile.write(MadWeightSecondStartLine)

else:
        print ' Used syntax is : python SplitLHCOFile.py "LeptonType" "NumberEvents" '

