#include "../interface/BTagStudy.h"

BTagStudy::BTagStudy(){
    BTagStudy::InitializeBegin();
}

BTagStudy::~BTagStudy(){

}

void BTagStudy::InitializeBegin(){

   for(int ii = 0; ii < 6; ii++){
      NotReconstructedEvent[ii]=0;
      allFourJetsCorrectlyMatched[ii]=0;
      atLeastOneWronglyMatched[ii] = 0;

      twoBTagsCorrectlyMatched[ii] = 0;
      atLeastOneBTagWronglyMatched[ii] = 0;
      twoLightJetsCorrectlyMatched[ii] = 0;
      atLeastOneLightJetWronglyMatched[ii] = 0;

      EventWithTwoLightJets[ii] = 0;
      EventWithThreeLightJets[ii] = 0;
      EventWithTwoLightJetsAndBTagged[ii] = 0;
      EventWithThreeLightJetsAndBTagged[ii] = 0;
      thirdJetIsActualQuark[ii] = 0; secondJetIsActualQuark[ii] = 0; firstJetIsActualQuark[ii] = 0;
      thirdJetIsCorrectQuark[ii] = 0;

      //selectedNumberEvents[ii] = 0;    //--> Possible to also count the number of selected events!! (then additional booleans are needed for semiMu of semiEl!)

      bTaggedJetNr[ii].clear();
      NonbTaggedJetNr[ii].clear();
      LightJetNr[ii].clear();
   }
}

void BTagStudy::InitializePerEvent(){

   for(int ii = 0; ii < 6; ii++){
      bTaggedJetNr[ii].clear();
      NonbTaggedJetNr[ii].clear();
      LightJetNr[ii].clear();
   }
}

void BTagStudy::CalculateJets(vector<TRootJet*> Jets, float BTagWorkingPoint, float LightWorkingPoint, int OptionNr){

     for(unsigned int ii = 0; ii<Jets.size();ii++){
        if(Jets[ii]->btag_combinedSecondaryVertexBJetTags() >= BTagWorkingPoint){
          bTaggedJetNr[OptionNr].push_back(ii);
        }
        else{
          NonbTaggedJetNr[OptionNr].push_back(ii);
        }

	//Calculate the light jets when an additional working point is required for these:
	if(BTagWorkingPoint != LightWorkingPoint && Jets[ii]->btag_combinedSecondaryVertexBJetTags() < LightWorkingPoint){
	  LightJetNr[OptionNr].push_back(ii);
	}
     }

     //Copy the Nonbtagged vector into the light one in case the two b-tags are the same!	
     if(BTagWorkingPoint == LightWorkingPoint){
	LightJetNr[OptionNr] = NonbTaggedJetNr[OptionNr];
     }

     //Count how often an event has three or two light jets (for comparison of 4- and 5-jet case)
     if( LightJetNr[OptionNr].size() >1 ){
	EventWithTwoLightJets[OptionNr]++;
        if(bTaggedJetNr[OptionNr].size() >1)
	   EventWithTwoLightJetsAndBTagged[OptionNr]++;
     }
     if( LightJetNr[OptionNr].size() >2 ){
	EventWithThreeLightJets[OptionNr]++;
	if(bTaggedJetNr[OptionNr].size() >1)
	   EventWithThreeLightJetsAndBTagged[OptionNr]++;
     } 
}

void BTagStudy::CorrectJetCombi(int BHadrIndex, int BLeptIndex, int Quark1Index, int Quark2Index, int OptionNr){

     if(BHadrIndex == 9999 || BLeptIndex == 9999 || Quark1Index == 9999 || Quark2Index == 9999){
	NotReconstructedEvent[OptionNr]++;
     }
     else{
        if( (BLeptIndex == bTaggedJetNr[OptionNr][0] || BLeptIndex == bTaggedJetNr[OptionNr][1]) &&
            (BHadrIndex == bTaggedJetNr[OptionNr][0] || BHadrIndex == bTaggedJetNr[OptionNr][1]) &&
            (Quark1Index == LightJetNr[OptionNr][0] || Quark1Index == LightJetNr[OptionNr][1]) &&
            (Quark2Index == LightJetNr[OptionNr][0] || Quark2Index == LightJetNr[OptionNr][1])){
                allFourJetsCorrectlyMatched[OptionNr]++;
        }
        if( (BLeptIndex != bTaggedJetNr[OptionNr][0] && BLeptIndex != bTaggedJetNr[OptionNr][1]) ||
            (BHadrIndex != bTaggedJetNr[OptionNr][0] && BHadrIndex != bTaggedJetNr[OptionNr][1]) ||
            (Quark1Index != LightJetNr[OptionNr][0] && Quark1Index != LightJetNr[OptionNr][1]) ||
            (Quark2Index != LightJetNr[OptionNr][0] && Quark2Index != LightJetNr[OptionNr][1])){
	 	atLeastOneWronglyMatched[OptionNr]++;
        }
	if( (BLeptIndex == bTaggedJetNr[OptionNr][0] || BLeptIndex == bTaggedJetNr[OptionNr][1]) &&
            (BHadrIndex == bTaggedJetNr[OptionNr][0] || BHadrIndex == bTaggedJetNr[OptionNr][1])){
		twoBTagsCorrectlyMatched[OptionNr]++;
	}
	if( (BLeptIndex != bTaggedJetNr[OptionNr][0] && BLeptIndex != bTaggedJetNr[OptionNr][1]) ||
            (BHadrIndex != bTaggedJetNr[OptionNr][0] && BHadrIndex != bTaggedJetNr[OptionNr][1])){
		atLeastOneBTagWronglyMatched[OptionNr]++;
	}
	if( (Quark1Index == LightJetNr[OptionNr][0] || Quark1Index == LightJetNr[OptionNr][1]) &&
            (Quark2Index == LightJetNr[OptionNr][0] || Quark2Index == LightJetNr[OptionNr][1])){
		twoLightJetsCorrectlyMatched[OptionNr]++;
	}
	if( (Quark1Index != LightJetNr[OptionNr][0] && Quark1Index != LightJetNr[OptionNr][1]) ||
            (Quark2Index != LightJetNr[OptionNr][0] && Quark2Index != LightJetNr[OptionNr][1])){
		atLeastOneLightJetWronglyMatched[OptionNr]++;
	}
     }
}

void BTagStudy::CorrectJetCombi5Jets(int BHadrIndex, int BLeptIndex, int Quark1Index, int Quark2Index, int OptionNr){

     if(BHadrIndex == 9999 || BLeptIndex == 9999 || Quark1Index == 9999 || Quark2Index == 9999){
	NotReconstructedEvent5Jets[OptionNr]++;
     }
     else{
        if( (BLeptIndex == bTaggedJetNr[OptionNr][0] || BLeptIndex == bTaggedJetNr[OptionNr][1]) &&
            (BHadrIndex == bTaggedJetNr[OptionNr][0] || BHadrIndex == bTaggedJetNr[OptionNr][1]) &&
            (Quark1Index == LightJetNr[OptionNr][0] || Quark1Index == LightJetNr[OptionNr][1] || Quark1Index == LightJetNr[OptionNr][2] ) &&
            (Quark2Index == LightJetNr[OptionNr][0] || Quark2Index == LightJetNr[OptionNr][1] || Quark2Index == LightJetNr[OptionNr][2] ) ){
                allFourJetsCorrectlyMatched5Jets[OptionNr]++;
 		if(Quark1Index == LightJetNr[OptionNr][2] || Quark2Index == LightJetNr[OptionNr][2]) thirdJetIsActualQuark[OptionNr]++;
		if(Quark1Index == LightJetNr[OptionNr][1] || Quark2Index == LightJetNr[OptionNr][1]) secondJetIsActualQuark[OptionNr]++;
		if(Quark1Index == LightJetNr[OptionNr][0] || Quark2Index == LightJetNr[OptionNr][0]) firstJetIsActualQuark[OptionNr]++;
        }
        if( (BLeptIndex != bTaggedJetNr[OptionNr][0] && BLeptIndex != bTaggedJetNr[OptionNr][1]) ||
            (BHadrIndex != bTaggedJetNr[OptionNr][0] && BHadrIndex != bTaggedJetNr[OptionNr][1]) ||
            (Quark1Index != LightJetNr[OptionNr][0] && Quark1Index != LightJetNr[OptionNr][1] && Quark1Index != LightJetNr[OptionNr][2] ) ||
            (Quark2Index != LightJetNr[OptionNr][0] && Quark2Index != LightJetNr[OptionNr][1] && Quark2Index != LightJetNr[OptionNr][2] ) ){
	 	atLeastOneWronglyMatched5Jets[OptionNr]++;
        }
	if( (BLeptIndex == bTaggedJetNr[OptionNr][0] || BLeptIndex == bTaggedJetNr[OptionNr][1]) &&
            (BHadrIndex == bTaggedJetNr[OptionNr][0] || BHadrIndex == bTaggedJetNr[OptionNr][1])){
		twoBTagsCorrectlyMatched5Jets[OptionNr]++;
		if(Quark1Index == LightJetNr[OptionNr][2] || Quark2Index == LightJetNr[OptionNr][2]) thirdJetIsGoodQuark[OptionNr]++;
	}
	if( (BLeptIndex != bTaggedJetNr[OptionNr][0] && BLeptIndex != bTaggedJetNr[OptionNr][1]) ||
            (BHadrIndex != bTaggedJetNr[OptionNr][0] && BHadrIndex != bTaggedJetNr[OptionNr][1])){
		atLeastOneBTagWronglyMatched5Jets[OptionNr]++;
	}
	if( (Quark1Index == LightJetNr[OptionNr][0] || Quark1Index == LightJetNr[OptionNr][1] || Quark1Index == LightJetNr[OptionNr][2] ) &&
            (Quark2Index == LightJetNr[OptionNr][0] || Quark2Index == LightJetNr[OptionNr][1] || Quark2Index == LightJetNr[OptionNr][2] ) ){
		twoLightJetsCorrectlyMatched5Jets[OptionNr]++;
 		if(Quark1Index == LightJetNr[OptionNr][2] || Quark2Index == LightJetNr[OptionNr][2]) thirdJetIsCorrectQuark[OptionNr]++;
	}
	if( (Quark1Index != LightJetNr[OptionNr][0] && Quark1Index != LightJetNr[OptionNr][1] && Quark1Index != LightJetNr[OptionNr][2] ) ||
            (Quark2Index != LightJetNr[OptionNr][0] && Quark2Index != LightJetNr[OptionNr][1] && Quark2Index != LightJetNr[OptionNr][2] ) ){
		atLeastOneLightJetWronglyMatched5Jets[OptionNr]++;
	}
     }
}

void BTagStudy::ReturnTable(std::string NameOfOption4Jets[6], std::string NameOfOption5Jets[6], int WhichJets, int NrOptionsConsidered, ostream &output, int OptionOfInterest){
	
    //Values of WhichJets are the following:
	// 0 = all 4 jets matched and compared
	// 1 = only the b-jets matched and compared
	// 2 = only the light jets matched and compared

//    std::string Title[3]= {"\\textbf{Option} (no $\\chi^{2}$ $m_{lb}$) & chosen jets are correct ($\\%$)       & $\\frac{s}{b}$ & 3rd jet is correct ($\\%$) \\\\", 
//			   "\\textbf{Option} (no $\\chi^{2}$ $m_{lb}$) & 2 b's chosen correctly ($\\%$)        & $\\frac{s}{b}$ & 3rd jet is correct ($\\%$) \\\\", 
//			   "\\textbf{Option} (no $\\chi^{2}$ $m_{lb}$) & chosen light jets are correct ($\\%$) & $\\frac{s}{b}$ & 3rd jet is correct ($\\%$) \\\\"};
    std::string Title[3]= {"\\textbf{Option} (no $\\chi^{2}$ $m_{lb}$) & all 4 correct & $\\geq$ 1 wrong & correct ($\\%$)       & $\\frac{s}{b}$ & non-matched \\\\", 
			   "\\textbf{Option} (no $\\chi^{2}$ $m_{lb}$) & 2 b's correct & $\\geq$ 1 b wrong & b's correct ($\\%$) & $\\frac{s}{b}$ & non-matched \\\\", 
			   "\\textbf{Option} (no $\\chi^{2}$ $m_{lb}$) & 2 light good  & $\\geq$ 1 light wrong & light correct ($\\%$) & $\\frac{s}{b}$ & non-matched \\\\"};

    std::string Caption[3] = {"\\caption{Overview of correct and wrong reconstructed events for the different b-tags without the use of a $\\chi^{2}$ $m_{lb}$ - $m_{qqb}$ method} ", 
			      "\\caption{Overview of correct and wrong reconstructed b-jets for the different b-tags without the use of a $\\chi^{2}$ $m_{lb}$ - $m_{qqb}$ method} ", 
			      "\\caption{Overview of correct and wrong reconstructed light jets for the different b-tags without the use of a $\\chi^{2}$ $m_{lb}$ - $m_{qqb}$ method} "};

    output << " \\begin{table}[!h] \n \\begin{tabular}{c|c|c|c|c} " << endl;
    output << Title[WhichJets] << " \\hline " << endl;
    for(int ii = 0; ii < NrOptionsConsidered; ii++){

	if(NrOptionsConsidered == 1){
	  ii = OptionOfInterest;
          int CorrectOnes5Jets[3] = {allFourJetsCorrectlyMatched5Jets[OptionOfInterest], twoBTagsCorrectlyMatched5Jets[OptionOfInterest] , twoLightJetsCorrectlyMatched5Jets[OptionOfInterest]};
          int WrongOnes5Jets[3] = {atLeastOneWronglyMatched5Jets[OptionOfInterest], atLeastOneBTagWronglyMatched5Jets[OptionOfInterest], atLeastOneLightJetWronglyMatched5Jets[OptionOfInterest]};
    	  int Correct3rdJet[3] = {thirdJetIsActualQuark[OptionOfInterest],thirdJetIsGoodQuark[OptionOfInterest], thirdJetIsCorrectQuark[OptionOfInterest]};
	  float sOverSqrtB5Jets[3], sOverB5Jets[3], CorrectPercentage5Jets[3], ThirdJetPercentage[3];
	  for(int jj = 0; jj < 3; jj++){
	     sOverSqrtB5Jets[jj] = (float)(CorrectOnes5Jets[jj])/(float)(sqrt(WrongOnes5Jets[jj]));
	     sOverB5Jets[jj] = (float)(CorrectOnes5Jets[jj])/(float)(WrongOnes5Jets[jj]);
	     CorrectPercentage5Jets[jj] = (float)(CorrectOnes5Jets[jj]*100.0)/(float)(CorrectOnes5Jets[jj]+WrongOnes5Jets[jj]);
	     ThirdJetPercentage[jj] = (float)(Correct3rdJet[jj]*100.0)/(float)(CorrectOnes5Jets[jj]);
	  }
	
	  //Print some additional output:	 
	  //cout << "\n Values for first, second and third light jet (Compared to total number of correctly matched events of " << allFourJetsCorrectlyMatched5Jets[OptionOfInterest] << " ) " << endl;
	  //cout << "    * First jet : " << firstJetIsActualQuark[OptionOfInterest] << endl;
	  //cout << "    * Second jet: " << secondJetIsActualQuark[OptionOfInterest] << endl;
	  //cout << "    * Third jet : " << thirdJetIsActualQuark[OptionOfInterest] << endl;

          // 5-jet output
          output << NameOfOption5Jets[ii]            << 
	  //" & " << CorrectOnes5Jets[WhichJets]       << 
	  //" & " << WrongOnes5Jets[WhichJets]         << 
	  " & " << CorrectPercentage5Jets[WhichJets] << 
	  " & " << sOverB5Jets[WhichJets]            << 
	  //" & " << NotReconstructedEvent5Jets[ii]    << 
	  //" & " << Correct3rdJet[WhichJets]          << 
	  " & " << ThirdJetPercentage[WhichJets]     << 
	  "\\\\ " << endl;
	} 
	int CorrectOnes[3] = {allFourJetsCorrectlyMatched[ii], twoBTagsCorrectlyMatched[ii]    , twoLightJetsCorrectlyMatched[ii]};
	int WrongOnes[3]   = {atLeastOneWronglyMatched[ii],    atLeastOneBTagWronglyMatched[ii], atLeastOneLightJetWronglyMatched[ii]};
        float sOverSqrtB[3], sOverB[3], CorrectPercentage[3];
        for(int jj = 0; jj < 3; jj++){
           sOverSqrtB[jj] = (float)(CorrectOnes[jj])/(float)(sqrt(WrongOnes[jj]));
           sOverB[jj] = (float)(CorrectOnes[jj])/(float)(WrongOnes[jj]);
           CorrectPercentage[jj] = (float)(CorrectOnes[jj]*100.0)/(float)(CorrectOnes[jj]+WrongOnes[jj]);
        }

        output << NameOfOption4Jets[ii]       << 
	" & " << CorrectOnes[WhichJets]       << 
	" & " << WrongOnes[WhichJets]         << 
	" & " << CorrectPercentage[WhichJets] << 
	" & " << sOverB[WhichJets]            << 
	" & " << NotReconstructedEvent[ii]    << 
	" & X \\\\ " << endl;
    }
    output << " \\end{tabular} " << endl;
    output << Caption[WhichJets] << endl;
    output << " \\end{table} \n " << endl;

}
