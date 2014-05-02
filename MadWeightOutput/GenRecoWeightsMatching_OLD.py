#!python
import linecache

NumberOfParams = 5
Gen_MadWeight = ["0.0","0.0","0.0","0.0","0.0"]
Reco_MadWeight =["0.0","0.0","0.0","0.0","0.0"] 
if len(Reco_MadWeight) != NumberOfParams:
	print 'Still need to change the Number of Parameters of the MadWeight arrays !! '

GenMuPlusFile = file('Gen_PositiveMuon_MadWeight_8000Evts_5Params.out','r') #RVL parameter changed between 1.2 and 0.8 in steps of 0.1
#GenMuMinusFile = file('Gen_NegativeMuon_MadWeight.out','r')
#GenElPlusFile = file('Gen_PositiveElec_MadWeight.out','r')
#GenElMinusFile = file('Gen_NegativeElec_MadWeight.out','r')

#RecoMuPlusFile = file('Reco_PositiveMuon_MadWeight.out','r')
#RecoMuMinusFile = file('Reco_NegativeMuon_MadWeight.out','r')
#RecoElPlusFile = file('Reco_PositiveElectron_MadWeight.out','r')
#RecoElMinusFile = file('Reco_NegativeElectron_MadWeight.out','r')

#File containing the matching between the generator and reco events numbers!
EventInfo = file('EventNumberInformation.txt','r')

LineCount = 0

for line in EventInfo:
    words = line.split()
    NumberWords = len(words) #Only need to perform the matching if the length is equal to 8! (this means that reconstruction is complete)
    #--> Will need to be updated when neutrinos are not reconstructed separately!!!
    if NumberWords == 8:
        if int(float(words[1])) == 1: #Positive muon case 
            print 'Matching positive muon event %s (event id %s) with reconstructed event %s' %(words[5],words[0],words[7])
	    #Get the MadWeight values from the .out file
	    #Last +3 is to take into account the two lines of text at the top!
	    for iterator in range(0,NumberOfParams):
		#Need the third word (hence 2), since this is the MadWeight value (3 will give the error on the weight)
    	        Gen_MadWeight[iterator] = linecache.getline('Gen_PositiveMuon_MadWeight_8000Evts_5Params.out',(int(float(words[5]))-1)*NumberOfParams+iterator+3).split()[2]
		Reco_MadWeight[iterator] = linecache.getline('Gen_PositiveMuon_MadWeight_8000Evts_5Params.out',(int(float(words[7]))-1)*NumberOfParams+iterator+3).split()[2]
		#print 'Gen_MadWeight : ',Gen_MadWeight[iterator],'for iterator = ',iterator,' (Event number = ',words[5],' ) -- Looking at line : ',(int(float(words[5]))-1)*NumberOfParams+iterator+3
		#print 'Reco_MadWeight : ',Reco_MadWeight[iterator],'for iterator = ',iterator,' (Event number = ',words[7],' ) -- Looking at line : ',(int(float(words[7]))-1)*NumberOfParams+iterator+3
            
