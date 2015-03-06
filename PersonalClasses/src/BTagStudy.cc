#include "../interface/BTagStudy.h"

BTagStudy::BTagStudy(int outputVerbose){
    BTagStudy::InitializeBegin();
    verbose = outputVerbose;
    eventSelectionOutput.open("/user/aolbrech/GitTopTree_Feb2014/TopBrussels/AnomalousCouplings/eventSelectionChoiceTables.tex");
}

BTagStudy::~BTagStudy(){
    eventSelectionOutput.close();
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

  //Try to include this in a more cleaner way ....
  float bjetWP[6] = {0.244,0.679,0.679,0.898,0.898,0.898};
  float lightjetWP[6] = {0.244,0.679,0.244,0.898,0.679,0.244};
  
  std::string optionName[6] = {"2 L b-tags             ",  //#0
  	                       "2 M b-tags             ",  //#1
	                       "2 M b-tags, light L-veto", //#2
	                       "2 T b-tags             ",  //#3
                               "2 T b-tags, light M-veto", //#4
	                       "2 T b-tags, light L-veto"};//#5
    for(int ii = 0; ii < 6; ii++){
        BJetWP[ii] = bjetWP[ii];
        LightJetWP[ii] = lightjetWP[ii];
        OptionName[ii] = optionName[ii];
    }
}

void BTagStudy::InitializePerEvent(){

   for(int ii = 0; ii < 6; ii++){
      bTaggedJetNr[ii].clear();
      NonbTaggedJetNr[ii].clear();
      LightJetNr[ii].clear();
   }
}

void BTagStudy::CalculateJets(vector<TRootJet*> Jets, vector<int> jetCombi){      //, float BTagWorkingPoint, float LightWorkingPoint, int OptionNr){

    for (unsigned int bTagOption = 0; bTagOption < 6; bTagOption++){

        for(unsigned int ii = 0; ii<Jets.size();ii++){
            if(Jets[ii]->btag_combinedSecondaryVertexBJetTags() >= BJetWP[bTagOption])
                bTaggedJetNr[bTagOption].push_back(ii);            
            else
                NonbTaggedJetNr[bTagOption].push_back(ii);            

	    //Calculate the light jets when an additional working point is required for these:
	    if(BJetWP[bTagOption] != LightJetWP[bTagOption] && Jets[ii]->btag_combinedSecondaryVertexBJetTags() < LightJetWP[bTagOption])
	        LightJetNr[bTagOption].push_back(ii);	
        }

        //Copy the Nonbtagged vector into the light one in case the two b-tags are the same!	
        if(BJetWP[bTagOption] == LightJetWP[bTagOption])
	    LightJetNr[bTagOption] = NonbTaggedJetNr[bTagOption];     

        //Count how often an event has three or two light jets (for comparison of 4- and 5-jet case)
        if( LightJetNr[bTagOption].size() >1 ){
	    EventWithTwoLightJets[bTagOption]++;
            if(bTaggedJetNr[bTagOption].size() >1)
	        EventWithTwoLightJetsAndBTagged[bTagOption]++;
        }
        if( LightJetNr[bTagOption].size() >2 ){
	    EventWithThreeLightJets[bTagOption]++;
	    if(bTaggedJetNr[bTagOption].size() >1)
	        EventWithThreeLightJetsAndBTagged[bTagOption]++;
        } 

        //--- Go to CorrectJetCombi class  ---//
        if( bTaggedJetNr[bTagOption].size() >= 2 && LightJetNr[bTagOption].size() >=2 ){
            if( LightJetNr[bTagOption].size() >= 2)
                CorrectJetCombi(jetCombi, bTagOption);      //4-jet case
            if( LightJetNr[bTagOption].size() > 2)
                CorrectJetCombi5Jets(jetCombi, bTagOption); //5-jet case
        }
        else
            if(verbose > 3) std::cout << " Event doesn't have two b-tagged jets and/or two light jets ! " << std::endl;

        //--- Additional output for debugging --//
        if(verbose > 3)
	    cout<<"(BTagStudy class) -- Size of bTaggedJets: "<<bTaggedJetNr[bTagOption].size()<<", of NonbTaggedJets: "<<NonbTaggedJetNr[bTagOption].size()<<" & of lightJets: "<<LightJetNr[bTagOption].size()<<endl;
    }//Loop over all btag options!
}

void BTagStudy::CorrectJetCombi(vector<int> jetCombi, int OptionNr){
    //jet combi order is : (0 = BLeptonic, 1 = BHadronic, 2 = Quark1 & 3 =  Quark2)

     if(jetCombi[1] == 9999 || jetCombi[0] == 9999 || jetCombi[2] == 9999 || jetCombi[3] == 9999){
	NotReconstructedEvent[OptionNr]++;
     }
     else{
        if( (jetCombi[0] == bTaggedJetNr[OptionNr][0] || jetCombi[0] == bTaggedJetNr[OptionNr][1]) &&
            (jetCombi[1] == bTaggedJetNr[OptionNr][0] || jetCombi[1] == bTaggedJetNr[OptionNr][1]) &&
            (jetCombi[2] == LightJetNr[OptionNr][0] || jetCombi[2] == LightJetNr[OptionNr][1]) &&
            (jetCombi[3] == LightJetNr[OptionNr][0] || jetCombi[3] == LightJetNr[OptionNr][1])){
                allFourJetsCorrectlyMatched[OptionNr]++;
        }
        if( (jetCombi[0] != bTaggedJetNr[OptionNr][0] && jetCombi[0] != bTaggedJetNr[OptionNr][1]) ||
            (jetCombi[1] != bTaggedJetNr[OptionNr][0] && jetCombi[1] != bTaggedJetNr[OptionNr][1]) ||
            (jetCombi[2] != LightJetNr[OptionNr][0] && jetCombi[2] != LightJetNr[OptionNr][1]) ||
            (jetCombi[3] != LightJetNr[OptionNr][0] && jetCombi[3] != LightJetNr[OptionNr][1])){
	 	atLeastOneWronglyMatched[OptionNr]++;
        }
	if( (jetCombi[0] == bTaggedJetNr[OptionNr][0] || jetCombi[0] == bTaggedJetNr[OptionNr][1]) &&
            (jetCombi[1] == bTaggedJetNr[OptionNr][0] || jetCombi[1] == bTaggedJetNr[OptionNr][1])){
		twoBTagsCorrectlyMatched[OptionNr]++;
	}
	if( (jetCombi[0] != bTaggedJetNr[OptionNr][0] && jetCombi[0] != bTaggedJetNr[OptionNr][1]) ||
            (jetCombi[1] != bTaggedJetNr[OptionNr][0] && jetCombi[1] != bTaggedJetNr[OptionNr][1])){
		atLeastOneBTagWronglyMatched[OptionNr]++;
	}
	if( (jetCombi[2] == LightJetNr[OptionNr][0] || jetCombi[2] == LightJetNr[OptionNr][1]) &&
            (jetCombi[3] == LightJetNr[OptionNr][0] || jetCombi[3] == LightJetNr[OptionNr][1])){
		twoLightJetsCorrectlyMatched[OptionNr]++;
	}
	if( (jetCombi[2] != LightJetNr[OptionNr][0] && jetCombi[2] != LightJetNr[OptionNr][1]) ||
            (jetCombi[3] != LightJetNr[OptionNr][0] && jetCombi[3] != LightJetNr[OptionNr][1])){
		atLeastOneLightJetWronglyMatched[OptionNr]++;
	}
     }
}

void BTagStudy::CorrectJetCombi5Jets(vector<int> jetCombi, int OptionNr){

    if(jetCombi[1] == 9999 || jetCombi[0] == 9999 || jetCombi[2] == 9999 || jetCombi[3] == 9999){
        NotReconstructedEvent5Jets[OptionNr]++;
    }
    else{
        if( (jetCombi[0] == bTaggedJetNr[OptionNr][0] || jetCombi[0] == bTaggedJetNr[OptionNr][1]) &&
            (jetCombi[1] == bTaggedJetNr[OptionNr][0] || jetCombi[1] == bTaggedJetNr[OptionNr][1]) &&
            (jetCombi[2] == LightJetNr[OptionNr][0] || jetCombi[2] == LightJetNr[OptionNr][1] || jetCombi[2] == LightJetNr[OptionNr][2] ) &&
            (jetCombi[3] == LightJetNr[OptionNr][0] || jetCombi[3] == LightJetNr[OptionNr][1] || jetCombi[3] == LightJetNr[OptionNr][2] ) ){
                allFourJetsCorrectlyMatched5Jets[OptionNr]++;
 		if(jetCombi[2] == LightJetNr[OptionNr][2] || jetCombi[3] == LightJetNr[OptionNr][2]) thirdJetIsActualQuark[OptionNr]++;
		if(jetCombi[2] == LightJetNr[OptionNr][1] || jetCombi[3] == LightJetNr[OptionNr][1]) secondJetIsActualQuark[OptionNr]++;
		if(jetCombi[2] == LightJetNr[OptionNr][0] || jetCombi[3] == LightJetNr[OptionNr][0]) firstJetIsActualQuark[OptionNr]++;
        }
        if( (jetCombi[0] != bTaggedJetNr[OptionNr][0] && jetCombi[0] != bTaggedJetNr[OptionNr][1]) ||
            (jetCombi[1] != bTaggedJetNr[OptionNr][0] && jetCombi[1] != bTaggedJetNr[OptionNr][1]) ||
            (jetCombi[2] != LightJetNr[OptionNr][0] && jetCombi[2] != LightJetNr[OptionNr][1] && jetCombi[2] != LightJetNr[OptionNr][2] ) ||
            (jetCombi[3] != LightJetNr[OptionNr][0] && jetCombi[3] != LightJetNr[OptionNr][1] && jetCombi[3] != LightJetNr[OptionNr][2] ) ){
	 	atLeastOneWronglyMatched5Jets[OptionNr]++;
        }
	if( (jetCombi[0] == bTaggedJetNr[OptionNr][0] || jetCombi[0] == bTaggedJetNr[OptionNr][1]) &&
            (jetCombi[1] == bTaggedJetNr[OptionNr][0] || jetCombi[1] == bTaggedJetNr[OptionNr][1])){
		twoBTagsCorrectlyMatched5Jets[OptionNr]++;
		//if(jetCombi[2] == LightJetNr[OptionNr][2] || jetCombi[3] == LightJetNr[OptionNr][2]) thirdJetIsGoodQuark[OptionNr]++;
	}
	if( (jetCombi[0] != bTaggedJetNr[OptionNr][0] && jetCombi[0] != bTaggedJetNr[OptionNr][1]) ||
            (jetCombi[1] != bTaggedJetNr[OptionNr][0] && jetCombi[1] != bTaggedJetNr[OptionNr][1])){
		atLeastOneBTagWronglyMatched5Jets[OptionNr]++;
	}
	if( (jetCombi[2] == LightJetNr[OptionNr][0] || jetCombi[2] == LightJetNr[OptionNr][1] || jetCombi[2] == LightJetNr[OptionNr][2] ) &&
            (jetCombi[3] == LightJetNr[OptionNr][0] || jetCombi[3] == LightJetNr[OptionNr][1] || jetCombi[3] == LightJetNr[OptionNr][2] ) ){
		twoLightJetsCorrectlyMatched5Jets[OptionNr]++;
 		//if(jetCombi[2] == LightJetNr[OptionNr][2] || jetCombi[3] == LightJetNr[OptionNr][2]) thirdJetIsCorrectQuark[OptionNr]++;
	}
	if( (jetCombi[2] != LightJetNr[OptionNr][0] && jetCombi[2] != LightJetNr[OptionNr][1] && jetCombi[2] != LightJetNr[OptionNr][2] ) ||
            (jetCombi[3] != LightJetNr[OptionNr][0] && jetCombi[3] != LightJetNr[OptionNr][1] && jetCombi[3] != LightJetNr[OptionNr][2] ) ){
		atLeastOneLightJetWronglyMatched5Jets[OptionNr]++;
	}
     }
}

void BTagStudy::ReturnTable(){ //std::string NameOfOption4Jets[6], std::string NameOfOption5Jets[6], int WhichJets){

    std::string OptionName4Jets[6], OptionName5Jets[6];
    for(int ii = 0; ii < 6; ii++){OptionName4Jets[ii] = " 4 jet case, "+OptionName[ii]; OptionName5Jets[ii] = " 5 jet case, "+OptionName[ii];}
	
    std::string Title[3]= {"   \\textbf{Option} (no $\\chi^{2}$ $m_{lb}$) & all 4 correct & $\\geq$ 1 wrong & correct ($\\%$)       & $\\frac{s}{b}$ & non-matched \\\\", 
			   "   \\textbf{Option} (no $\\chi^{2}$ $m_{lb}$) & 2 b's correct & $\\geq$ 1 b wrong & b's correct ($\\%$) & $\\frac{s}{b}$ & non-matched \\\\", 
			   "   \\textbf{Option} (no $\\chi^{2}$ $m_{lb}$) & 2 light good  & $\\geq$ 1 light wrong & light correct ($\\%$) & $\\frac{s}{b}$ & non-matched \\\\"};

    std::string Caption[3] = {"   \\caption{Overview of correct and wrong reconstructed events for the different b-tags (no $\\chi^{2}$ $m_{lb}$ - $m_{qqb}$ applied)} ", 
			      "   \\caption{Overview of correct and wrong reconstructed b-jets for the different b-tags (no $\\chi^{2}$ $m_{lb}$ - $m_{qqb}$ applied)} ", 
			      "   \\caption{Overview of correct and wrong reconstructed light jets for the different b-tags (no $\\chi^{2}$ $m_{lb}$ - $m_{qqb}$ applied)} "};

    std::string Title5Jets[3]= {"   \\textbf{Option} (no $\\chi^{2}$ $m_{lb}$, 5 jets) & all 4 correct & $\\geq$ 1 wrong & correct ($\\%$)       & $\\frac{s}{b}$ & non-matched \\\\", 
			        "   \\textbf{Option} (no $\\chi^{2}$ $m_{lb}$, 5 jets) & 2 b's correct & $\\geq$ 1 b wrong & b's correct ($\\%$) & $\\frac{s}{b}$ & non-matched \\\\", 
			        "   \\textbf{Option} (no $\\chi^{2}$ $m_{lb}$, 5 jets) & 2 light good  & $\\geq$ 1 light wrong & light correct ($\\%$) & $\\frac{s}{b}$ & non-matched \\\\"};

    std::string Caption5Jets[3] = {"   \\caption{Overview of correct and wrong reconstructed events for the different b-tags (no $\\chi^{2}$ $m_{lb}$ - $m_{qqb}$ applied, 5 jets considered)} ", 
			           "   \\caption{Overview of correct and wrong reconstructed b-jets for the different b-tags (no $\\chi^{2}$ $m_{lb}$ - $m_{qqb}$ applied, 5 jets considered)} ", 
			           "   \\caption{Overview of correct and wrong reconstructed light jets for the different b-tags (no $\\chi^{2}$ $m_{lb}$ - $m_{qqb}$ applied, 5 jets considered)} "};


    //Maybe create a separate output file for 4 and 5 jet case?
    //--> This avoids creating a double bTag loop!

    for(int itNrJets = 0; itNrJets < 3; itNrJets++){      //Replacing loop over whichJets! (for 5jet case)

        for(int itBTag = 0; itBTag < 6; itBTag++){

            //---  Store 5 jet information  ---//
            if(itBTag == 0) eventSelectionOutput << " \\begin{table}[!h] \n  \\begin{tabular}{c|c|c|c|c} " << endl;
            eventSelectionOutput << Title5Jets[itNrJets] << " \\hline " << endl;

            int CorrectOnes5Jets[3] = {allFourJetsCorrectlyMatched5Jets[itBTag], twoBTagsCorrectlyMatched5Jets[itBTag] ,    twoLightJetsCorrectlyMatched5Jets[itBTag]};
            int WrongOnes5Jets[3]   = {atLeastOneWronglyMatched5Jets[itBTag],    atLeastOneBTagWronglyMatched5Jets[itBTag], atLeastOneLightJetWronglyMatched5Jets[itBTag]};
            //int Correct3rdJet[3]    = {thirdJetIsActualQuark[itBTag],            thirdJetIsGoodQuark[OptionOfInterest], thirdJetIsCorrectQuark[OptionOfInterest]};
	    //float sOverSqrtB5Jets[3], sOverB5Jets[3], CorrectPercentage5Jets[3], ThirdJetPercentage[3];
	    /*for(int jj = 0; jj < 3; jj++){
	        sOverSqrtB5Jets[jj] = (float)(CorrectOnes5Jets[jj])/(float)(sqrt(WrongOnes5Jets[jj]));
	        sOverB5Jets[jj] = (float)(CorrectOnes5Jets[jj])/(float)(WrongOnes5Jets[jj]);
	        CorrectPercentage5Jets[jj] = (float)(CorrectOnes5Jets[jj]*100.0)/(float)(CorrectOnes5Jets[jj]+WrongOnes5Jets[jj]);
	        ThirdJetPercentage[jj] = (float)(Correct3rdJet[jj]*100.0)/(float)(CorrectOnes5Jets[jj]);
	    }*/
	
            eventSelectionOutput << OptionName5Jets[itNrJets]     << 
	    " & " << CorrectOnes5Jets[itNrJets]       << 
	    " & " << WrongOnes5Jets[itNrJets]         << 
	    //" & " << CorrectPercentage5Jets[WhichJets] << 
	    //" & " << sOverB5Jets[WhichJets]            << 
	    //" & " << NotReconstructedEvent5Jets[ii]    << 
	    //" & " << Correct3rdJet[WhichJets]          << 
	    //" & " << ThirdJetPercentage[WhichJets]     << 
	    "\\\\ " << endl;

            if(itBTag == 5){
                eventSelectionOutput << "  \\end{tabular} " << endl;
                eventSelectionOutput << Caption5Jets[itNrJets] << endl;
                eventSelectionOutput << " \\end{table} \n " << endl;
            }
	}
        

        //---  Store 4 jet information  ---//         
	/*int CorrectOnes[3] = {allFourJetsCorrectlyMatched[ii], twoBTagsCorrectlyMatched[ii]    , twoLightJetsCorrectlyMatched[ii]};
	int WrongOnes[3]   = {atLeastOneWronglyMatched[ii],    atLeastOneBTagWronglyMatched[ii], atLeastOneLightJetWronglyMatched[ii]};
        float sOverSqrtB[3], sOverB[3], CorrectPercentage[3];
        for(int jj = 0; jj < 3; jj++){
           sOverSqrtB[jj] = (float)(CorrectOnes[jj])/(float)(sqrt(WrongOnes[jj]));
           sOverB[jj] = (float)(CorrectOnes[jj])/(float)(WrongOnes[jj]);
           CorrectPercentage[jj] = (float)(CorrectOnes[jj]*100.0)/(float)(CorrectOnes[jj]+WrongOnes[jj]);
        }

        eventSelectionOutput << OptionName4Jets[ii]       << 
	" & " << CorrectOnes[WhichJets]       << 
	" & " << WrongOnes[WhichJets]         << 
	" & " << CorrectPercentage[WhichJets] << 
	" & " << sOverB[WhichJets]            << 
	" & " << NotReconstructedEvent[ii]    << 
	" & X \\\\ " << endl;*/
    }

}
