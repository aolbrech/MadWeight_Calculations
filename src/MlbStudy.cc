#include "../interface/MlbStudy.h"

//MlbStudy::MlbStudy(){
//}

//MlbStudy::~MlbStudy(){

//}

void MlbStudy::initializePerEvent(){
 
  for(int ii = 0; ii < 6; ii++){
     ChiSquared[ii] = 999.;
     MlbValues[ii] = 999.;
     MqqbValues[ii] = 999.;
  }
  chosenBLept = 999; //This means that the hightest Pt jet is considered as the leptonic b-jet!
  chosenBHadr = 999;
  chosenQuark1 = 999;
  chosenQuark2 = 999;
  LowestChiSq = 999;
  LowestChiSq4Jets = 999;
  CorrectChiSq = 999;
}

void MlbStudy::initializeBegin(){

  for(int ii = 0; ii < 6; ii++){
   NumberMatchedEvents[ii] = 0;
   NumberNotMatchedEvents[ii] = 0;
   CorrectOptionAvailable[ii] = 0;
   CorrectOptionChosen[ii] = 0;
   WrongOptionChosen[ii] = 0;
   CorrectEventMlbMqqb[ii] = 0;
   WrongEventMlbMqqb[ii] = 0;
  }
}

void MlbStudy::calculateChiSquared(vector<int> CorrectValues, vector<int> bTaggedJets, vector<int> lightJets, TLorentzVector* lepton, vector<TRootJet*> Jets, float MassMlb, float SigmaMlb, float MassMqqb, float SigmaMqqb){

   if(bTaggedJets.size() > 1 && lightJets.size() > 1){                                            //Event needs to have two b-tagged jets and two light jets!

       float MlbValuesLocal[6]  = {(*lepton+*Jets[bTaggedJets[0]]).M(), (*lepton+*Jets[bTaggedJets[1]]).M(),
			           (*lepton+*Jets[bTaggedJets[0]]).M(), (*lepton+*Jets[bTaggedJets[1]]).M(),
			           (*lepton+*Jets[bTaggedJets[0]]).M(), (*lepton+*Jets[bTaggedJets[1]]).M()};

       float MqqbValuesLocal[6] = {(*Jets[bTaggedJets[1]] + *Jets[lightJets[0]] + *Jets[lightJets[1]]).M(), (*Jets[bTaggedJets[0]] + *Jets[lightJets[0]] + *Jets[lightJets[1]]).M(),
	         	           10000., 10000.,10000., 10000.};
       if(lightJets.size() > 2){  //5-jet case!  --> In case there are only two light jets available it automatically becomes the 4-jet case!
         MqqbValuesLocal[2] = (*Jets[bTaggedJets[1]] + *Jets[lightJets[0]] + *Jets[lightJets[2]]).M();
	 MqqbValuesLocal[3] = (*Jets[bTaggedJets[0]] + *Jets[lightJets[0]] + *Jets[lightJets[2]]).M();
	 MqqbValuesLocal[4] = (*Jets[bTaggedJets[1]] + *Jets[lightJets[1]] + *Jets[lightJets[2]]).M();
	 MqqbValuesLocal[5] = (*Jets[bTaggedJets[0]] + *Jets[lightJets[1]] + *Jets[lightJets[2]]).M();
       }

       //Copy events to general MlbValues and MqqbValues array!
       for(int ii = 0; ii < 6; ii++){
	  MlbValues[ii] = MlbValuesLocal[ii];
	  MqqbValues[ii]= MqqbValuesLocal[ii];
       }

       LowestChiSq = 0;
       LowestChiSq4Jets = 0;
       for(int ii = 0; ii<6; ii++){
           ChiSquared[ii] = ((MassMlb-MlbValues[ii])/SigmaMlb)*((MassMlb-MlbValues[ii])/SigmaMlb) + ((MassMqqb - MqqbValues[ii])/SigmaMqqb)*((MassMqqb - MqqbValues[ii])/SigmaMqqb);
	   if( ii > 0 && ChiSquared[ii] < ChiSquared[LowestChiSq]) LowestChiSq = ii;
	   if(ii == 1 && ChiSquared[0] > ChiSquared[1]) LowestChiSq4Jets = 1;
       }

       //Check which Chi-Sq value corresponds to the correct event topology!
       if(CorrectValues[0] == bTaggedJets[0] && CorrectValues[1] == bTaggedJets[1] && (CorrectValues[2] == lightJets[0] || CorrectValues[2] == lightJets[1]) && (CorrectValues[3] == lightJets[0] || CorrectValues[3] == lightJets[1])) CorrectChiSq = 0;
       if(CorrectValues[0] == bTaggedJets[1] && CorrectValues[0] == bTaggedJets[1] && (CorrectValues[2] == lightJets[0] || CorrectValues[2] == lightJets[1]) && (CorrectValues[3] == lightJets[0] || CorrectValues[3] == lightJets[1])) CorrectChiSq = 1;
       if(CorrectValues[0]==bTaggedJets[0]&&CorrectValues[1]==bTaggedJets[1]&&(CorrectValues[2]==lightJets[0]||CorrectValues[2]==lightJets[2])&&(CorrectValues[3]==lightJets[0]||CorrectValues[3]==lightJets[2]))CorrectChiSq=2;
       if(CorrectValues[0]==bTaggedJets[1]&&CorrectValues[0]==bTaggedJets[1]&&(CorrectValues[2]==lightJets[0]||CorrectValues[2]==lightJets[2])&&(CorrectValues[3]==lightJets[0]||CorrectValues[3]==lightJets[2]))CorrectChiSq=3;
       if(CorrectValues[0]==bTaggedJets[0]&&CorrectValues[1]==bTaggedJets[1]&&(CorrectValues[2]==lightJets[1]||CorrectValues[2]==lightJets[2])&&(CorrectValues[3]==lightJets[1]||CorrectValues[3]==lightJets[2]))CorrectChiSq=4;
       if(CorrectValues[0]==bTaggedJets[1]&&CorrectValues[0]==bTaggedJets[1]&&(CorrectValues[2]==lightJets[1]||CorrectValues[2]==lightJets[2])&&(CorrectValues[3]==lightJets[1]||CorrectValues[3]==lightJets[2]))CorrectChiSq=5;

   }
}

void MlbStudy::getIndices(int LowestChiSqIndex){                  //Making this a separate class allows to choice in the analyzer code whether the 4- or 5-jet case will be considered!

   //Match the correct indices to the quarks!
   // --> b-jets
   if(      LowestChiSqIndex == 0 || LowestChiSqIndex == 2 || LowestChiSqIndex == 4){ chosenBLept = 0; chosenBHadr = 1;}
   else if( LowestChiSqIndex == 1 || LowestChiSqIndex == 3 || LowestChiSqIndex == 5){ chosenBLept = 1; chosenBHadr = 0;}
   // --> light jets
   if      ( LowestChiSqIndex == 0 || LowestChiSqIndex == 1){ chosenQuark1 = 0; chosenQuark2 = 1;}
   else if ( LowestChiSqIndex == 2 || LowestChiSqIndex == 3){ chosenQuark1 = 0; chosenQuark2 = 2;}
   else if ( LowestChiSqIndex == 4 || LowestChiSqIndex == 5){ chosenQuark1 = 1; chosenQuark2 = 2;}

}

void MlbStudy::calculateEfficiency(int option, vector<int> CorrectValues, vector<int> bTaggedJets, vector<int> lightJets, int NrConsideredBTagOptions){
  
   if(bTaggedJets.size() > 1 && lightJets.size() > 1){                                                                  //Event has the correct amount of b-tagged and light jets!

     int NumberTimesLoopShouldBeRepeated = 1;
     if(NrConsideredBTagOptions == 1) NumberTimesLoopShouldBeRepeated = 3;

     for(int ii = 0; ii < NumberTimesLoopShouldBeRepeated; ii++){

	//For the 'pure' 5-jet case the number of light jets should be larger than 2!
	if(NrConsideredBTagOptions == 1 && ii == 2 && lightJets.size() < 3) continue;

	//Get the correct indices, plot the minimum chi-sq values and initialize the chi-sq counter!
	int NrChiSqs = 6;
	if(NrConsideredBTagOptions > 1){            getIndices(LowestChiSq4Jets);NrChiSqs = 0;} 
	if(NrConsideredBTagOptions == 1 && ii == 0){getIndices(LowestChiSq);                      h_ChiSqMinimum[option].push_back(ChiSquared[LowestChiSq]);}      //5-jet case
	if(NrConsideredBTagOptions == 1 && ii == 1){getIndices(LowestChiSq4Jets);option=option+1; h_ChiSqMinimum[option].push_back(ChiSquared[LowestChiSq4Jets]); NrChiSqs = 2;} //4-jet case
	if(NrConsideredBTagOptions == 1 && ii == 2){getIndices(LowestChiSq);     option=option+1; h_ChiSqMinimum[option].push_back(ChiSquared[LowestChiSq]);} //pure 5-jet case! (+1 since +1 done for 4jets)

	//Plot the distribution for the correct and wrong chi-sq values
        for(int ll = 0; ll < 6; ll++){ if(ChiSquared[ll] >=50) ChiSquared[ii] = 49.9;}  //Force the ChiSquared values to be in the ROOT histograms to see the overflow!
	if(CorrectChiSq != 999){ 
	   h_ChiSqCorrectFound[option].push_back(ChiSquared[CorrectChiSq]); 
	   h_ChiSqCorrect[option].push_back(ChiSquared[CorrectChiSq]);
	}
	else h_ChiSqCorrect[option].push_back(49.9);  //Otherwise the ChiSquared[999] is set to 0 .... which is not the result desired to compare the distributions!
	//Wrong case
	for(int jj = 0; jj < NrChiSqs; jj++){
	   if(jj != CorrectChiSq) h_ChiSqWrong[option].push_back(ChiSquared[jj]);
	}

	//Plot the maximum chi-sq values:
	for(int jj = 0; jj < NrChiSqs; jj++){
	  if(NrConsideredBTagOptions == 1 && ( (ii == 1 && jj != LowestChiSq4Jets) || ( (ii == 2 || ii == 0) && jj != LowestChiSq ) ) ) h_ChiSqNotMinimum[option].push_back(ChiSquared[jj]);
	}

	//----------------------------//
	//  Get the matched events!!  //
	//----------------------------//    
	if(  CorrectValues[0] != 9999 && CorrectValues[1] != 9999 && CorrectValues[2] != 9999 && CorrectValues[3] != 9999 ){   //Event has been matched!
           NumberMatchedEvents[option]++;

	   ///////////////////////////////////////////////
           //  How often is the correct option chosen?  //
	   ///////////////////////////////////////////////
   	   //--> 4-jet case
           if( ( (NrConsideredBTagOptions == 1 && ( ii == 1 || (ii == 0 && lightJets.size() < 3) ) ) || NrConsideredBTagOptions > 1 ) &&   //Will never look at 5 jets AND all the b-tag options ...
	       (CorrectValues[0] == bTaggedJets[0] || CorrectValues[0] == bTaggedJets[1])   &&
               (CorrectValues[1] == bTaggedJets[1] || CorrectValues[1] == bTaggedJets[0])   ){  //Correct option is available!
		  CorrectOptionAvailable[option]++;

		  h_ChiSqCorrectWhenMatched[option].push_back(ChiSquared[CorrectChiSq]);
		  h_ChiSqMinimumWhenMatched[option].push_back(ChiSquared[LowestChiSq4Jets]);
		  for(int jj = 0; jj < NrChiSqs; jj++){
                    if(NrConsideredBTagOptions == 1 && ( (ii == 1 && jj != LowestChiSq4Jets) || ( (ii == 2 || ii == 0) && jj != LowestChiSq ) ) )
			h_ChiSqNotMinimumWhenMatched[option].push_back(ChiSquared[jj]); 
		  }

		  //Correct option chosen or not?
	   	  if( CorrectValues[0] == bTaggedJets[chosenBLept] && CorrectValues[1] == bTaggedJets[chosenBHadr] && 
		     (CorrectValues[2] == lightJets[chosenQuark1] || CorrectValues[2] == lightJets[chosenQuark2] ) &&
		     (CorrectValues[3] == lightJets[chosenQuark2] || CorrectValues[3] == lightJets[chosenQuark2] ) ){
		        CorrectOptionChosen[option]++;
			h_ChiSqMinimumWhenCorrect[option].push_back(ChiSquared[LowestChiSq4Jets]);
		  }
		  else{
			WrongOptionChosen[option]++;		 
			h_ChiSqMinimumWhenWrong[option].push_back(ChiSquared[LowestChiSq4Jets]);
		  }
	   }

	   //--> 5-jet case
           if(  (NrConsideredBTagOptions == 1 && (ii == 0 || ii == 2) && lightJets.size() > 2 )  &&   
	        (CorrectValues[0] == bTaggedJets[0] || CorrectValues[0] == bTaggedJets[1])   &&
                (CorrectValues[1] == bTaggedJets[1] || CorrectValues[1] == bTaggedJets[0])   &&
		(CorrectValues[2] == lightJets[0] || CorrectValues[2] == lightJets[1] || CorrectValues[2] == lightJets[2]) &&
		(CorrectValues[3] == lightJets[0] || CorrectValues[3] == lightJets[1] || CorrectValues[3] == lightJets[2]) ){  //Correct option is available!
		  CorrectOptionAvailable[option]++;

		  h_ChiSqCorrectWhenMatched[option].push_back(ChiSquared[CorrectChiSq]);
		  h_ChiSqMinimumWhenMatched[option].push_back(ChiSquared[LowestChiSq4Jets]);
		  for(int jj = 0; jj < NrChiSqs; jj++){
                    if(NrConsideredBTagOptions == 1 && ( (ii == 1 && jj != LowestChiSq4Jets) || ( (ii == 2 || ii == 0) && jj != LowestChiSq ) ) ) 
			h_ChiSqNotMinimumWhenMatched[option].push_back(ChiSquared[jj]); 
		  }

		  //Correct option chosen or not?
	   	  if( CorrectValues[0] == bTaggedJets[chosenBLept] && CorrectValues[1] == bTaggedJets[chosenBHadr] && 
		     (CorrectValues[2] == lightJets[chosenQuark1]  || CorrectValues[2] == lightJets[chosenQuark2] ) &&
		     (CorrectValues[3] == lightJets[chosenQuark2]  || CorrectValues[3] == lightJets[chosenQuark2] ) ){
		        CorrectOptionChosen[option]++;
			h_ChiSqMinimumWhenCorrect[option].push_back(ChiSquared[LowestChiSq4Jets]);
                  }
		  else{
		 	WrongOptionChosen[option]++;		  
			h_ChiSqMinimumWhenWrong[option].push_back(ChiSquared[LowestChiSq4Jets]);
                  }
	   }	   
	
           //How often is the full event found and correctly reconstructed?
           if( CorrectValues[0] == bTaggedJets[chosenBLept] && CorrectValues[1] == bTaggedJets[chosenBHadr]  &&
             ( CorrectValues[2] == lightJets[chosenQuark1]  || CorrectValues[2] == lightJets[chosenQuark2] ) &&
             ( CorrectValues[3] == lightJets[chosenQuark1]  || CorrectValues[3] == lightJets[chosenQuark2] ) )
                CorrectEventMlbMqqb[option]++;
           else WrongEventMlbMqqb[option]++;
        }//End of loop over matched events !
	else{
	   NumberNotMatchedEvents[option]++;
	   for(int jj = 0; jj < NrChiSqs; jj++) h_ChiSqAllWhenNotMatched[option].push_back(ChiSquared[jj]);
	}

     }//End of loop over ii!
  }//End of loop over events with at least two b-tagged jets and two light jets
}

void MlbStudy::saveNumbers(std::string NameOfOption[6], int WhichJets, int NrOptionsConsidered, ofstream &output, int OptionOfInterest){

   std::string Title[2] = {" \\textbf{Option} (with $\\chi^{2}$ $m_{lb}$) & all 4 correct & $\\geq$ 1 wrong  & $\\frac{s}{\\sqrt{b}}$ & $\\frac{s}{b}$ & non-matched \\\\", 
			   " \\textbf{Option} (with $\\chi^{2}$ $m_{lb}$) & Correct b's   & Wrong b's & \\% good chosen b's & $\\frac{s}{b}$ & Correct option exists \\\\"};

   std::string Caption[2] = {" \\caption{Overview of correct and wrong reconstructed events for the different b-tags when a $\\chi^{2}$ $m_{lb}$ - $m_{qqb}$ method is applied} ", 
			     " \\caption{Overview of the number of times the correct b-jet combination is chosen when using a $\\chi^{2}$ $m_{lb}$ - $m_{qqb}$ method} "};

   output << " \\begin{table}[!h] \n \\begin{tabular}{c|c|c|c|c|c} " << endl;
   output << Title[WhichJets] << " \\hline " << endl;

   bool lookingAtOneBTagOption = false;
   if(NrOptionsConsidered == 1){ NrOptionsConsidered = 3; lookingAtOneBTagOption = true;}
   for(int ii = 0; ii < NrOptionsConsidered; ii++){
	if(ii == 0 && lookingAtOneBTagOption == true) ii = OptionOfInterest;
	if(ii == 1 && lookingAtOneBTagOption == true) ii = OptionOfInterest+1;
	if(ii == 2 && lookingAtOneBTagOption == true) ii = OptionOfInterest+2;

       int CorrectOnes[2] = {CorrectEventMlbMqqb[ii],   CorrectOptionChosen[ii] };
       int WrongOnes[2]   = {WrongEventMlbMqqb[ii],     WrongOptionChosen[ii]};
       int LastOnes[2]    = {NumberNotMatchedEvents[ii],CorrectOptionAvailable[ii]}; 
       float sOverSqrtBORPercentage[2] = {(float)(CorrectOnes[WhichJets])/(float)(sqrt(WrongOnes[WhichJets])), (float)(((float)(CorrectOptionChosen[ii])*100.0)/(float)(CorrectOptionAvailable[ii]))};

       output << NameOfOption[ii]                                << " & " <<
       CorrectOnes[WhichJets]                                    << " & " <<
       WrongOnes[WhichJets]                                      << " & " <<
       sOverSqrtBORPercentage[WhichJets]                         << " & " <<
       float(CorrectOnes[WhichJets])/float(WrongOnes[WhichJets]) << " & " <<
       LastOnes[WhichJets]                                       << " \\\\ " << endl;

       if( ii == OptionOfInterest   && lookingAtOneBTagOption == true) ii = 0;  //Otherwise loop will not continue!
       if( ii == OptionOfInterest+1 && lookingAtOneBTagOption == true) ii = 1;  //Otherwise loop will not continue!
     }
    output << " \\end{tabular} " << endl;
    output << Caption[WhichJets] << endl;
    output << " \\end{table} \n " << endl;

}

