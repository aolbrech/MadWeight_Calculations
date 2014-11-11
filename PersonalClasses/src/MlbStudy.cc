#include "../interface/MlbStudy.h"

MlbStudy::MlbStudy(int NrBTags){
    MlbStudy::initializeBegin(NrBTags);
}

MlbStudy::~MlbStudy(){

}

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

void MlbStudy::initializeBegin(int NrBTags){

  for(int ii = 0; ii < 6; ii++){
   NumberMatchedEvents[ii] = 0;
   NumberNotMatchedEvents[ii] = 0;
   CorrectOptionAvailable[ii] = 0;
   CorrectEventChosen[ii] = 0;
   WrongEventChosen[ii] = 0;
   CorrectBOptionChosen[ii] = 0;
   WrongBOptionChosen[ii] = 0;

   ThirdQuarkShouldBeChosen[ii] = 0;
   ThirdQuarkChosen[ii] = 0; SecondQuarkChosen[ii] = 0; FirstQuarkChosen[ii] = 0;
   ThirdQuarkCorrectChosen[ii] = 0; SecondQuarkCorrectChosen[ii] = 0; FirstQuarkCorrectChosen[ii] = 0;
   ThirdQuarkCorrect[ii] = 0; SecondQuarkCorrect[ii] = 0; FirstQuarkCorrect[ii] = 0;
   CorrectLightJetsChosen[ii] = 0;
   CorrectLightJetsWithThirdChosen[ii] = 0;
  }

  //Save the histograms belonging to the mlb output information! 
  if(NrBTags == 1){ 
    Title[0] = "5Jets"; Title[1] = "4Jets"; Title[2] = "Pure5Jets";
    Name[0] =  " - 5 jets case) "; Name[1] = " - 4 jets case) "; Name[2] = " - pure 5 jets case) ";
 
    for(int ii = 0; ii < 3; ii++){
        histo1D["ChiSqCorrect"+Title[ii]] = new TH1F(("ChiSqCorrect"+Title[ii]).c_str(), ("#chi^{2} distribution for correct combi (all events"+Name[ii]).c_str(), 150, 0, 50);
        histo1D["ChiSqCorrectFound"+Title[ii]] = new TH1F(("ChiSqCorrectFound"+Title[ii]).c_str(), ("#chi^{2} distribution for correct combi (when found"+Name[ii]).c_str(), 150, 0, 50);
        histo1D["ChiSqMin"+Title[ii]] = new TH1F(("ChiSqMin"+Title[ii]).c_str(), ("#chi^{2} distribution of the minimal combi (all events"+Name[ii]).c_str(), 150, 0, 50);
        histo1D["ChiSqNotMin"+Title[ii]] = new TH1F(("ChiSqNotMin"+Title[ii]).c_str(),("#chi^{2} distribution of the non-minimal combis (all events"+Name[ii]).c_str(), 150, 0, 50);
        histo1D["ChiSqWrong"+Title[ii]] = new TH1F(("ChiSqWrong"+Title[ii]).c_str(), ("#chi^{2} distribution of the non-correct combis (all events"+Name[ii]).c_str(), 150, 0, 50);
         
        histo1D["ChiSqCorrWhenMatched"+Title[ii]] = new TH1F(("ChiSqCorrWhenMatched"+Title[ii]).c_str(), ("#chi^{2} distribution for correct combi (matched events only"+Name[ii]).c_str(), 150, 0, 50);
        histo1D["ChiSqMinWhenMatched"+Title[ii]] = new TH1F(("ChiSqMinWhenMatched"+Title[ii]).c_str(), ("#chi^{2} distribution for minimal combi (matched events only"+Name[ii]).c_str(), 150, 0, 50);
        histo1D["ChiSqNotMinWhenMatched"+Title[ii]] = new TH1F(("ChiSqNotMinWhenMatched"+Title[ii]).c_str(), ("#chi^{2} distribution for non-minimal combis (matched events only"+Name[ii]).c_str(), 150, 0, 50);
        histo1D["ChiSqMinWhenCorrect"+Title[ii]] = new TH1F(("ChiSqMinWhenCorrect"+Title[ii]).c_str(), ("#chi^{2} distribution for minimal combi (correct choice only"+Name[ii]).c_str(), 150, 0, 50);
        histo1D["ChiSqMinWhenWrong"+Title[ii]] = new TH1F(("ChiSqMinWhenWrong"+Title[ii]).c_str(), ("#chi^{2} distribution for minimal combi (wrong choice only"+Name[ii]).c_str(), 150, 0, 50);
        histo1D["ChiSqAllWhenNotMatched"+Title[ii]] = new TH1F(("ChiSqAllWhenNotMatched"+Title[ii]).c_str(), ("#chi^{2} distribution for all combis (non-matched evens only"+Name[ii]).c_str(), 150, 0, 50);
    
        histo1D["ChiSqDiffWhenWrong"+Title[ii]] = new TH1F(("ChiSqDiffWhenWrong"+Title[ii]).c_str(), ("#chi^{2}_{correct} - #chi^{2}_{minimum} for the wrong choice"+Name[ii]).c_str(), 500, -5, 15);
    }
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
       if( CorrectValues[0] == bTaggedJets[0] && CorrectValues[1] == bTaggedJets[1] && 
           (CorrectValues[2] == lightJets[0] || CorrectValues[2] == lightJets[1]) && (CorrectValues[3] == lightJets[0] || CorrectValues[3] == lightJets[1]) ) CorrectChiSq = 0;
       if( CorrectValues[0] == bTaggedJets[1] && CorrectValues[0] == bTaggedJets[1] && 
	   (CorrectValues[2] == lightJets[0] || CorrectValues[2] == lightJets[1]) && (CorrectValues[3] == lightJets[0] || CorrectValues[3] == lightJets[1]) ) CorrectChiSq = 1;
       if( CorrectValues[0] == bTaggedJets[0] && CorrectValues[1] == bTaggedJets[1] && lightJets.size() > 2 &&
	   (CorrectValues[2] == lightJets[0] || CorrectValues[2] == lightJets[2]) && (CorrectValues[3] == lightJets[0] || CorrectValues[3] == lightJets[2]) ) CorrectChiSq = 2;
       if( CorrectValues[0] == bTaggedJets[1] && CorrectValues[0] == bTaggedJets[1] && lightJets.size() > 2 &&
	   (CorrectValues[2] == lightJets[0] || CorrectValues[2] == lightJets[2]) && (CorrectValues[3] == lightJets[0] || CorrectValues[3] == lightJets[2]) ) CorrectChiSq = 3;
       if( CorrectValues[0] == bTaggedJets[0] && CorrectValues[1] == bTaggedJets[1] && lightJets.size() > 2 &&
	   (CorrectValues[2] == lightJets[1] || CorrectValues[2] == lightJets[2]) && (CorrectValues[3] == lightJets[1] || CorrectValues[3] == lightJets[2]) ) CorrectChiSq = 4;
       if( CorrectValues[0] == bTaggedJets[1] && CorrectValues[0] == bTaggedJets[1] && lightJets.size() > 2 &&
	   (CorrectValues[2] == lightJets[1] || CorrectValues[2] == lightJets[2]) && (CorrectValues[3] == lightJets[1] || CorrectValues[3] == lightJets[2]) ) CorrectChiSq = 5;

	//cout << "   --  End of mlbStudy CalculateChiSq function " << endl;  //Output corresponds to numbers obtained in analyzer code!
	//cout << "     --  ChiSq Indices (4 jets and 4+ jets case) : " << LowestChiSq4Jets << " & " << LowestChiSq << endl;
	//cout << "     --  ChiSq index for correct combination     : " << CorrectChiSq << endl;
	//cout << "     --  ChiSq values : " << ChiSquared[0] << " , " << ChiSquared[1] << " , " << ChiSquared[2] << " , " << ChiSquared[3] << " , " << ChiSquared[4] << " & " << ChiSquared[5] << endl;
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

void MlbStudy::calculateEfficiency(int option, vector<int> CorrectValues, vector<int> bTaggedJets, vector<int> lightJets, int NrConsideredBTagOptions, int ChiSqCutValue){
  
   if(bTaggedJets.size() > 1 && lightJets.size() > 1){                                                                  //Event has the correct amount of b-tagged and light jets!

     int NumberTimesLoopShouldBeRepeated = 1;
     if(NrConsideredBTagOptions == 1) NumberTimesLoopShouldBeRepeated = 3;

     for(int ii = 0; ii < NumberTimesLoopShouldBeRepeated; ii++){

	//For the 'pure' 5-jet case the number of light jets should be larger than 2!
	if(NrConsideredBTagOptions == 1 && ii == 2 && lightJets.size() < 3) continue;

	//Get the correct indices, plot the minimum chi-sq values and initialize the chi-sq counter!
	int NrChiSqs = 6;
	int UsedLowestChiSq = LowestChiSq;
	int UsedCorrectChiSq = CorrectChiSq;
	if(NrConsideredBTagOptions > 1){ UsedLowestChiSq = LowestChiSq4Jets; getIndices(LowestChiSq4Jets); NrChiSqs = 0; if(CorrectChiSq > 1) UsedCorrectChiSq = 999;}
	else if(NrConsideredBTagOptions == 1){
	   if(ii == 1){ option=option+1; NrChiSqs = 2; UsedLowestChiSq = LowestChiSq4Jets; if(CorrectChiSq > 1) UsedCorrectChiSq = 999;} //4-jet case
	   if(ii == 2){ option=option+1;                                                                                               } //pure 5-jet case! (+1 since +1 done for 4jets)
	   getIndices(UsedLowestChiSq);
	   histo1D["ChiSqMin"+Title[ii]]->Fill(ChiSquared[UsedLowestChiSq]);
	}
        for(int ll = 0; ll < 6; ll++){ if(ChiSquared[ll] >=50) ChiSquared[ll] = 49.9;}  //Force the ChiSquared values to be in the ROOT histograms to see the overflow!
	if(ChiSquared[UsedLowestChiSq] > ChiSqCutValue) continue;

	//Plot the distribution for the correct and wrong chi-sq values
        if(NrConsideredBTagOptions == 1){     //Only plot for this case!
	    if(UsedCorrectChiSq != 999){
	        histo1D["ChiSqCorrectFound"+Title[ii]]->Fill(ChiSquared[UsedCorrectChiSq]); 
	        histo1D["ChiSqCorrect"+Title[ii]]->Fill(ChiSquared[UsedCorrectChiSq]);
	    }
	    else histo1D["ChiSqCorrect"+Title[ii]]->Fill(49.9);  //Otherwise the ChiSquared[999] is set to 0 .... which is not the result desired to compare the distributions!
	    //Wrong case
	    for(int jj = 0; jj < NrChiSqs; jj++){
	        if(jj != UsedCorrectChiSq) histo1D["ChiSqWrong"+Title[ii]]->Fill(ChiSquared[jj]);
	    }

	    //Plot the maximum chi-sq values:
	    for(int jj = 0; jj < NrChiSqs; jj++){
	        if(jj != UsedLowestChiSq ) histo1D["ChiSqNotMin"+Title[ii]]->Fill(ChiSquared[jj]);
	    }
        }

	//----------------------------//
	//  Get the matched events!!  //
	//----------------------------//    
	if(  CorrectValues[0] != 9999 && CorrectValues[1] != 9999 && CorrectValues[2] != 9999 && CorrectValues[3] != 9999 ){   //Event has been matched!
           NumberMatchedEvents[option]++;

	   // First check whether the correct option is available
           if( ( (NrConsideredBTagOptions==1 && (ii==0 || ii==2) && lightJets.size() > 2) && (CorrectValues[0] == bTaggedJets[0] || CorrectValues[0] == bTaggedJets[1])   &&
		  									     (CorrectValues[1] == bTaggedJets[1] || CorrectValues[1] == bTaggedJets[0])   &&
      											     (CorrectValues[2] == lightJets[0] || CorrectValues[2] == lightJets[1] || CorrectValues[2] == lightJets[2]) &&
									  		     (CorrectValues[3] == lightJets[0] || CorrectValues[3] == lightJets[1] || CorrectValues[3] == lightJets[2]) )
	       || ( (NrConsideredBTagOptions > 1 || (NrConsideredBTagOptions == 1 && ((ii == 0 && lightJets.size() == 2) || ii == 1) )) &&
											     (CorrectValues[0] == bTaggedJets[0] || CorrectValues[0] == bTaggedJets[1]) &&
                                                                                             (CorrectValues[1] == bTaggedJets[1] || CorrectValues[1] == bTaggedJets[0]) &&
                                                                                             (CorrectValues[2] == lightJets[0]   || CorrectValues[2] == lightJets[1]  ) &&
                                                                                             (CorrectValues[3] == lightJets[0]   || CorrectValues[3] == lightJets[1]  ) )      
	       ){
		CorrectOptionAvailable[option]++;
		
		//See how often each of the light jets are chosen, and how often it is the correct choice!
		if(CorrectValues[2] == lightJets[2] || CorrectValues[3] == lightJets[2]) ThirdQuarkCorrect[option]++;
		if(CorrectValues[2] == lightJets[1] || CorrectValues[3] == lightJets[1]) SecondQuarkCorrect[option]++;
		if(CorrectValues[2] == lightJets[0] || CorrectValues[3] == lightJets[0]) FirstQuarkCorrect[option]++;

		if(chosenQuark1 == 2 || chosenQuark2 == 2){ ThirdQuarkChosen[option]++;
		    if(CorrectValues[2] == lightJets[2] || CorrectValues[3] == lightJets[2]) ThirdQuarkCorrectChosen[option]++;  //Not the number of correct ones, but the number of correct ones when this quark is actually chosen!!
		}
		if(chosenQuark1 == 1 || chosenQuark2 == 1){ SecondQuarkChosen[option]++;
		    if(CorrectValues[2] == lightJets[1] || CorrectValues[3] == lightJets[1]) SecondQuarkCorrectChosen[option]++;
		}
		if(chosenQuark1 == 0 || chosenQuark2 == 0){FirstQuarkChosen[option]++;
		    if(CorrectValues[2] == lightJets[0] || CorrectValues[3] == lightJets[0]) FirstQuarkCorrectChosen[option]++;
		}

		if(lightJets.size() > 2 && (CorrectValues[2] == lightJets[2] || CorrectValues[3] == lightJets[2]) && (chosenQuark1 != 2 && chosenQuark2 != 2)) ThirdQuarkShouldBeChosen[option]++;
		if((CorrectValues[2]==lightJets[chosenQuark1]||CorrectValues[2]==lightJets[chosenQuark2]) && (CorrectValues[3]==lightJets[chosenQuark2]||CorrectValues[3]==lightJets[chosenQuark2])){
		    CorrectLightJetsChosen[option]++;
		    if(chosenQuark1 == 2 || chosenQuark2 ==2) CorrectLightJetsWithThirdChosen[option]++;
		}
		
		if(NrConsideredBTagOptions == 1){ //Fill the histograms!
        	    histo1D["ChiSqCorrWhenMatched"+Title[ii]]->Fill(ChiSquared[UsedCorrectChiSq]);
	    	    histo1D["ChiSqMinWhenMatched"+Title[ii]]->Fill(ChiSquared[UsedLowestChiSq]);
		    for(int jj = 0; jj < NrChiSqs; jj++){
		        if( jj != UsedLowestChiSq){ histo1D["ChiSqNotMinWhenMatched"+Title[ii]]->Fill(ChiSquared[jj]); }
		    }
                }
		
		//Correct b-jets chosen or not?
	   	if(CorrectValues[0] == bTaggedJets[chosenBLept] && CorrectValues[1] == bTaggedJets[chosenBHadr] ){
		   CorrectBOptionChosen[option]++;
		  
		   // .. also correct light jets chosen (so complete event is correctly reconstructed)?
		   if((CorrectValues[2]==lightJets[chosenQuark1]||CorrectValues[2]==lightJets[chosenQuark2]) && (CorrectValues[3]==lightJets[chosenQuark2]||CorrectValues[3]==lightJets[chosenQuark2])){
		        CorrectEventChosen[option]++;
		        if(NrConsideredBTagOptions == 1) histo1D["ChiSqMinWhenCorrect"+Title[ii]]->Fill(ChiSquared[UsedLowestChiSq]);
                   }
		   else{
		        WrongEventChosen[option]++;
		        if(NrConsideredBTagOptions == 1) histo1D["ChiSqMinWhenWrong"+Title[ii]]->Fill(ChiSquared[UsedLowestChiSq]);
			if(NrConsideredBTagOptions == 1) histo1D["ChiSqDiffWhenWrong"+Title[ii]]->Fill(ChiSquared[UsedCorrectChiSq]-ChiSquared[UsedLowestChiSq]);
                    }
	   	}
		else{
		    WrongBOptionChosen[option]++;
		    WrongEventChosen[option]++;
		    if(NrConsideredBTagOptions == 1) histo1D["ChiSqMinWhenWrong"+Title[ii]]->Fill(ChiSquared[UsedLowestChiSq]);
		    if(NrConsideredBTagOptions == 1) histo1D["ChiSqDiffWhenWrong"+Title[ii]]->Fill(ChiSquared[UsedCorrectChiSq]-ChiSquared[UsedLowestChiSq]);
		}	   
	   }//Correct option is available
	   else{
		CorrectMatchingNotExisting[option]++;
		WrongEventChosen[option]++;
		//if(NrConsideredBTagOptions == 1){ for(int jj = 0; jj < NrChiSqs;jj++) histo1D["ChiSqAllMatchingNotExisting"+Title[ii]]->Fill(ChiSquared[jj]);}
	   }
 
        }//End of loop over matched events !
	else{
	   NumberNotMatchedEvents[option]++;
	   if(NrConsideredBTagOptions == 1){ for(int jj = 0; jj < NrChiSqs; jj++) histo1D["ChiSqAllWhenNotMatched"+Title[ii]]->Fill(ChiSquared[jj]);}
	}

     }//End of loop over ii!
  }//End of loop over events with at least two b-tagged jets and two light jets
}

typedef basic_ofstream<char> ofstream;

void MlbStudy::saveNumbers(std::string NameOfOption[6], int LookAtBJets, int NrOptionsConsidered, int ChosenOption, std::string StrForChiSqCutValue){

   //Initialize and define the output files:
   ofstream mlbOutput;
   ofstream mlbMatchingOutput;
   if(NrOptionsConsidered > 1){
	mlbOutput.open("/user/aolbrech/GitTopTree_Feb2014/TopBrussels/AnomalousCouplings/mlbChoiceTables.tex");
	mlbMatchingOutput.open("/user/aolbrech/GitTopTree_Feb2014/TopBrussels/AnomalousCouplings/mlbMatchingChoiceTables.tex");
   }
   else if(NrOptionsConsidered == 1){
	mlbOutput.open(("/user/aolbrech/GitTopTree_Feb2014/TopBrussels/AnomalousCouplings/mlbTableForChosenCombination"+StrForChiSqCutValue+".tex").c_str()); 
	mlbMatchingOutput.open(("/user/aolbrech/GitTopTree_Feb2014/TopBrussels/AnomalousCouplings/mlbMatchingTableForChosenCombination"+StrForChiSqCutValue+".tex").c_str());
	NameOfOption[ChosenOption+2] = " Pure 5 jet case, "+NameOfOption[ChosenOption];
        NameOfOption[ChosenOption+1] = " 4 jet case,      "+NameOfOption[ChosenOption];
	NameOfOption[ChosenOption] =   " 5 jet case,      "+NameOfOption[ChosenOption];
   }

   //Set the titles and captions of the tables!
   //std::string Title[2] = {"\\multirow{2}{*}{\\textbf{Option} (with $\\chi^{2}$ $m_{lb}$)} & 4 chosen jets & $\\frac{s}{b}$ & 3rd jet is one of the & 3rd jet is chosen \\\\ & are correct ($\\%$) & & 2 correct light jets ($\\%$) & and correct ($\\%$) \\\\",
   //			   " \\textbf{Option} (with $\\chi^{2}$ $m_{lb}$) & \\% b's correct & $\\frac{s}{b}$ & & \\\\"};
  
   std::string Title[2] = {"\\textbf{Option} (with $\\chi^{2}$ $m_{lb}$) &all 4 good &$\\geq$ 1 wrong &4 chosen jets good ($\\%$) &$\\frac{s}{b}$ & non-matched\\\\",
			   "\\textbf{Option} (with $\\chi^{2}$ $m_{lb}$) & Correct b's & Wrong b's & & \\% b's correct   & $\\frac{s}{b}$ & Correct option exists \\\\"};
   std::string Caption[2] = {" \\caption{Overview of correct and wrong reconstructed events for the different b-tags when a $\\chi^{2}$ $m_{lb}$ - $m_{qqb}$ method is applied} ", 
			     " \\caption{Overview of the number of times the correct b-jet combination is chosen when using a $\\chi^{2}$ $m_{lb}$ - $m_{qqb}$ method} "};

   std::string TitleMatching[2] = {"\\textbf{Option} (with $\\chi^{2}$ $m_{lb}$) & \\frac}{s}{b} matching found &$ \\frac{s}{b} matching good & \\frac{1}{\\sqrt{s}} matching found &$\\frac{}{\\sqrt{s}}$ matching good & non-matched\\\\",
                                   "\\textbf{Option} (with $\\chi^{2}$ $m_{lb}$) & & & & & &  \\\\"};
   std::string CaptionMatching[2] = {" \\caption{Matching information when a $\\chi^{2}$ $m_{lb}$ - $m_{qqb}$ method is applied} ", 
                                     " \\caption{Matching information for the light jets specific when using a $\\chi^{2}$ $m_{lb}$ - $m_{qqb}$ method} "};

   //Run b-jet info separately:
   bool lookingAtOneBTagOption = false;
   if(NrOptionsConsidered == 1){ NrOptionsConsidered = 3; lookingAtOneBTagOption = true;}
   for(int BJetAlso = 0; BJetAlso <= LookAtBJets; BJetAlso++){
     //Filling of the tables	
     mlbOutput << " \\begin{table}[!h] \n \\begin{tabular}{c|c|c|c|c} \n " << Title[BJetAlso] << "\n \\hline " << endl;
     mlbMatchingOutput << " \\begin{table}[!h] \n \\begin{tabular}{c|c|c|c|c} \n " << TitleMatching[BJetAlso] << "\n \\hline " << endl;
     for(int ii = 0; ii < NrOptionsConsidered; ii++){
         if(ii == 0 && lookingAtOneBTagOption == true) ii = ChosenOption;
         if(ii == 1 && lookingAtOneBTagOption == true) ii = ChosenOption+1;
         if(ii == 2 && lookingAtOneBTagOption == true) ii = ChosenOption+2;
 
         //Get all the relevant information
         int CorrectOnes[2] = {CorrectEventChosen[ii],    CorrectBOptionChosen[ii] };
         int WrongOnes[2]   = {WrongEventChosen[ii],      WrongBOptionChosen[ii]};
         int LastOnes[2]    = {NumberNotMatchedEvents[ii],CorrectOptionAvailable[ii]}; 
         //int ThirdJetChosen[2] = {ThirdQuarkChosen[ii], 0};
         //int ThirdJetCorrectChosen[2] = {ThirdQuarkCorrectChosen[ii], 0};
         //int ThirdJetShouldBeChosen[2] = {ThirdQuarkShouldBeChosen[ii], 0};
         //float TimesThirdJetIsACorrectJet[2] = {((float)CorrectLightJetsWithThirdChosen[ii]*100.0/(float)CorrectLightJetsChosen[ii]), 0};
         //float TimesThirdJetIsChosenANDCorrect[2] = {((float)ThirdQuarkCorrectChosen[ii]*100.0/(float)ThirdQuarkChosen[ii]), 0};
         float MatchingSOverB[2]      = {(float)(CorrectOptionAvailable[ii]*100.0)/(float)(CorrectMatchingNotExisting[ii]+NumberNotMatchedEvents[ii]),0};
         float MatchingExactSOverB[2] = {(float)(CorrectEventChosen[ii]*100.0)/(float)(WrongEventChosen[ii]+NumberNotMatchedEvents[ii])              ,0};
         float SqrtSInverse[2]        = {1/(float)(sqrt(CorrectOptionAvailable[ii]))                                                                 ,0};
         float SrtSExactInverse[2]    = {1/(float)(sqrt(CorrectEventChosen[ii]))                                                                     ,0};
         //float TimesEventIsMatched[2] = {((

         float sOverSqrtB[2], CorrectPercentage[2], sOverB[2];
         for(int jj = 0; jj < 2; jj++){
  	   sOverSqrtB[jj] = (float)(CorrectOnes[jj])/(float)(sqrt(WrongOnes[jj]));
	   CorrectPercentage[jj] = (float)(CorrectOnes[jj]*100.0)/(float)(CorrectOptionAvailable[ii]);
	   sOverB[jj] = (float)(CorrectOnes[jj])/(float)(WrongOnes[BJetAlso]);
         }

	 //Additional output about light jet percentages:
	 //cout << " Values when first, second and third jet are correct (Compared to the number of correctly matched events of " << CorrectOptionAvailable[ii] << " ) " << endl;
	 //cout << "     * First jet : " << FirstQuarkCorrect[ii] << endl;
	 //cout << "     * Second jet: " << SecondQuarkCorrect[ii] << endl;
	 //cout << "     * Third jet : " << ThirdQuarkCorrect[ii] << endl;

         //First table: Information about correct and wrong reconstructed events!
         mlbOutput << NameOfOption[ii]   << 
  	 " & " << CorrectOnes[BJetAlso]       << 
	 " & " << WrongOnes[BJetAlso]         << 
	 " & " << CorrectPercentage[BJetAlso] << 
	 " & " << sOverB[BJetAlso]             << 
	 " & " << LastOnes[BJetAlso]          << 
	 " \\\\ " << endl;

	 //Second table: Information about matching!
         mlbMatchingOutput << NameOfOption[ii] <<
         " & " << MatchingSOverB[BJetAlso] <<
         " & " << MatchingExactSOverB[BJetAlso] <<
         " & " << SqrtSInverse[BJetAlso] <<
         " & " << SrtSExactInverse[BJetAlso] <<
         " \\\\ " << endl;

         if( ii == ChosenOption   && lookingAtOneBTagOption == true) ii = 0;  //Otherwise loop will not continue!
         if( ii == ChosenOption+1 && lookingAtOneBTagOption == true) ii = 1;  //Otherwise loop will not continue!
     }
     mlbOutput << " \\end{tabular} \n " << Caption[BJetAlso] << "\n \\end{table} \n " << endl;
     mlbMatchingOutput << " \\end{tabular} \n " << CaptionMatching[BJetAlso] << "\n \\end{table} \n " << endl;
  }
  mlbOutput.close();
  mlbMatchingOutput.close();

}

void MlbStudy::WritePlots(TFile* outfile){
	outfile->cd();
	std::cout << " Inside WritePlots class ! " << std::endl;
	std::cout << " Histograms will be filled in file : " << outfile->GetName() << " ************************************" << std::endl;

	TDirectory* th1dir = outfile->mkdir("1D_histograms_Mlb");
	th1dir->cd();
	for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++){
	      TH1F *temp = it->second;
	      int N = temp->GetNbinsX();
	      temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
	      temp->SetBinContent(N+1,0);
	      temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
	      temp->Write();
	}
	TDirectory* th2dir = outfile->mkdir("2D_histograms_Mlb");
	th2dir->cd();
	for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++){    
	      TH2F *temp = it->second;
	      temp->Write();
	}
	outfile->cd(); //Step out from 2D_histograms_graphs directory!
}

