#include "../interface/BTagStudy.h"

BTagStudy::BTagStudy(int outputVerbose){
    BTagStudy::InitializeBegin();
    verbose = outputVerbose;
    evtSelOutput.open("/user/aolbrech/GitTopTree_Feb2014/TopBrussels/AnomalousCouplings/eventSelectionChoiceTables4JetCase.tex");
    evtSelOutput5Jets.open("/user/aolbrech/GitTopTree_Feb2014/TopBrussels/AnomalousCouplings/eventSelectionChoiceTables5JetCase.tex");
}

BTagStudy::~BTagStudy(){
    evtSelOutput.close();
    evtSelOutput5Jets.close();
}

void BTagStudy::InitializeBegin(){

  for(int ii = 0; ii < 6; ii++){
    for(int jj = 0; jj < 2; jj++){
      NotReconstructedEvent[ii][jj]=0;
      allFourJetsCorrectlyMatched[ii][jj]=0;
      atLeastOneWronglyMatched[ii][jj] = 0;

      twoBTagsCorrectlyMatched[ii][jj] = 0;
      atLeastOneBTagWronglyMatched[ii][jj] = 0;
      twoLightJetsCorrectlyMatched[ii][jj] = 0;
      atLeastOneLightJetWronglyMatched[ii][jj] = 0;

      /*EventWithTwoLightJets[ii] = 0;
      EventWithThreeLightJets[ii] = 0;
      EventWithTwoLightJetsAndBTagged[ii] = 0;
      EventWithThreeLightJetsAndBTagged[ii] = 0;*/
      thirdJetIsActualQuark[ii][jj] = 0; secondJetIsActualQuark[ii][jj] = 0; firstJetIsActualQuark[ii][jj] = 0;
      //thirdJetIsCorrectQuark[ii] = 0;
    }
  }

  //Try to include this in a more cleaner way ....
  float bjetWP[6]      = { 0.244,                    0.679,                    0.679,                     0.898,                    0.898,                     0.898                    };
  float lightjetWP[6]  = { 0.244,                    0.679,                                0.244,         0.898,                                0.679,                     0.244        };
  string optionName[6] = {"2 L b-tags             ","2 M b-tags             ","2 M b-tags, light L-veto","2 T b-tags             ","2 T b-tags, light M-veto","2 T b-tags, light L-veto"};
  for(int ii = 0; ii < 6; ii++){
    BJetWP[ii] = bjetWP[ii];
    LightJetWP[ii] = lightjetWP[ii];
    OptionName[ii] = optionName[ii];
  }

  //ChiSq Mlb and Mqqb information!
  Mlb = 103.286;
  SigmaMlb = 26.7764;
  Mqqb = 178.722;
  SigmaMqqb = 18.1385;
}

void BTagStudy::ResetEventArrays(){

  for(int ii = 0; ii < 6; ii++){
    bTaggedJetNr[ii].clear();
    NonbTaggedJetNr[ii].clear();
    LightJetNr[ii].clear();
    ChiSquaredMlb[ii].clear();
    ChiSquaredMqqb[ii].clear();
  }
}

void BTagStudy::CalculateJets(vector<TRootJet*> Jets, vector<int> jetCombi, TLorentzVector* lepton){      //, float BTagWorkingPoint, float LightWorkingPoint, int OptionNr){

  ResetEventArrays();

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

    /*//Count how often an event has three or two light jets (for comparison of 4- and 5-jet case)    --> Maybe move this 4/5 jet comparison to a separate function!
    if(LightJetNr[bTagOption].size() >1 ){
      EventWithTwoLightJets[bTagOption]++;
      if(bTaggedJetNr[bTagOption].size() >1)
	EventWithTwoLightJetsAndBTagged[bTagOption]++;
    }
    if(LightJetNr[bTagOption].size() >2 ){
      EventWithThreeLightJets[bTagOption]++;
      if(bTaggedJetNr[bTagOption].size() >1)
	EventWithThreeLightJetsAndBTagged[bTagOption]++;
    } 
    */

    if(bTaggedJetNr[bTagOption].size() >= 2 && LightJetNr[bTagOption].size() >=2 ){

      //--- Use Mlb chi-sq to select the b-jets ---//
      CalculateMlbChiSq(bTagOption, lepton, Jets);

      //--- For 4-jet case jet combi numbers can be calculated directly ---//
      light1Index4Jets[bTagOption] = (LightJetNr[bTagOption])[0];
      light2Index4Jets[bTagOption] = (LightJetNr[bTagOption])[1];
      CompareJetCombi(jetCombi, bTagOption, 0, light1Index4Jets[bTagOption], light2Index4Jets[bTagOption]);      //4-jet case

      //--- In case of 5 jets, select the two correct light ones ---//
      light1Index5Jets[bTagOption] = (LightJetNr[bTagOption])[0];
      light2Index5Jets[bTagOption] = (LightJetNr[bTagOption])[1];
      if(LightJetNr[bTagOption].size() > 2){ 
        vector<int> lightMqqbIndex;
        lightMqqbIndex = CalculateMqqbChiSq(bTagOption, Jets);
        light1Index5Jets[bTagOption] = light1IndexMqqb[bTagOption];
        light2Index5Jets[bTagOption] = light2IndexMqqb[bTagOption];
        std::cout << " lightMqqbIndex[0] = " << lightMqqbIndex[0] << " =?= light1IndexMqqb[bTagOption] = " << light1IndexMqqb[bTagOption] << std::endl;
        std::cout << " lightMqqbIndex[1] = " << lightMqqbIndex[1] << " =?= light2IndexMqqb[bTagOption] = " << light2IndexMqqb[bTagOption] << std::endl;
      }

      //--- Now use the two chosen light jets and go to CorrectJetCombi calculation! ---//
      CompareJetCombi(jetCombi, bTagOption, 1, light1Index5Jets[bTagOption], light2Index5Jets[bTagOption]); //4+5-jet case
    }
    else{
      if(verbose > 3) std::cout << " Event doesn't have two b-tagged jets and/or two light jets ! " << std::endl;
      bHadrIndex[bTagOption] = 999;
      bLeptIndex[bTagOption] = 999;
      LowestChiSqMlb[bTagOption] = 999;
      for(int ii = 0; ii < 2; ii++)
        ChiSquaredMlb[bTagOption].push_back(999);
    }
      
    //--- Additional output for debugging --//
    if(verbose > 3)
      cout<<"(BTagStudy class) -- Size of bTaggedJets: "<<bTaggedJetNr[bTagOption].size()<<", of NonbTaggedJets: "<<NonbTaggedJetNr[bTagOption].size()<<" & of lightJets: "<<LightJetNr[bTagOption].size()<<endl;
  }//Loop over all btag options!
}

void BTagStudy::CalculateMlbChiSq(int bTagNr, TLorentzVector* lepton, vector<TRootJet*> Jets){

  //---  This will distinguish the leptonic and hadronic b-jets  ---//
  //---     ==> Identical for 4- and 5-jet case                  ---//
  float MlbValues[2]  = {(*lepton+*Jets[(bTaggedJetNr[bTagNr])[0]]).M(), (*lepton+*Jets[(bTaggedJetNr[bTagNr])[1]]).M()};
 
  LowestChiSqMlb[bTagNr] = 0;
  bHadrIndex[bTagNr] = (bTaggedJetNr[bTagNr])[1];
  for(int ii = 0; ii<2; ii++){
    ChiSquaredMlb[bTagNr].push_back((float)(((Mlb-(float)(MlbValues[ii]))/SigmaMlb)*((Mlb-(float)(MlbValues[ii]))/SigmaMlb))); 
    if((ChiSquaredMlb[bTagNr])[ii] < (ChiSquaredMlb[bTagNr])[LowestChiSqMlb[bTagNr]]){
      LowestChiSqMlb[bTagNr] = ii;
      bHadrIndex[bTagNr] = (bTaggedJetNr[bTagNr])[0];
    }
  }
  bLeptIndex[bTagNr] = (bTaggedJetNr[bTagNr])[LowestChiSqMlb[bTagNr]];
}

vector<int> BTagStudy::CalculateMqqbChiSq(int bTagNr, vector<TRootJet*> Jets){
  //---  This option chooses 2 out of three light jets  ---//
  //---     ==> Only considered in 5-jet case           ---//

  vector<int> lightIndices;
  lightIndices.clear();
  float MqqbValues[3] = {(*Jets[bHadrIndex[bTagNr]] + *Jets[(LightJetNr[bTagNr])[0]] + *Jets[(LightJetNr[bTagNr])[1]]).M(), 
                         (*Jets[bHadrIndex[bTagNr]] + *Jets[(LightJetNr[bTagNr])[0]] + *Jets[(LightJetNr[bTagNr])[2]]).M(),
                         (*Jets[bHadrIndex[bTagNr]] + *Jets[(LightJetNr[bTagNr])[1]] + *Jets[(LightJetNr[bTagNr])[2]]).M()};
  
  LowestChiSqMqqb[bTagNr] = 0;
  for(int ii = 0; ii<3; ii++){
    ChiSquaredMqqb[bTagNr].push_back((float)(((Mqqb-(float)(MqqbValues[ii]))/SigmaMqqb)*((Mqqb-(float)(MqqbValues[ii]))/SigmaMqqb))); 
    if((ChiSquaredMqqb[bTagNr])[ii] < (ChiSquaredMqqb[bTagNr])[LowestChiSqMqqb[bTagNr]]){
      LowestChiSqMqqb[bTagNr] = ii;
    }
  }
  
  //Set the light jet indices:
  if(LowestChiSqMqqb[bTagNr] == 0)     { lightIndices.push_back( (LightJetNr[bTagNr])[0]); lightIndices.push_back( (LightJetNr[bTagNr])[1]); light1IndexMqqb[bTagNr] = (LightJetNr[bTagNr])[0]; light2IndexMqqb[bTagNr] = (LightJetNr[bTagNr])[1]; }
  else if(LowestChiSqMqqb[bTagNr] == 1){ lightIndices.push_back( (LightJetNr[bTagNr])[0]); lightIndices.push_back((LightJetNr[bTagNr])[2]); light1IndexMqqb[bTagNr] = (LightJetNr[bTagNr])[0]; light2IndexMqqb[bTagNr] = (LightJetNr[bTagNr])[2];}
  else if(LowestChiSqMqqb[bTagNr] == 2){ lightIndices.push_back( (LightJetNr[bTagNr])[1]); lightIndices.push_back((LightJetNr[bTagNr])[2]); light1IndexMqqb[bTagNr] = (LightJetNr[bTagNr])[1]; light2IndexMqqb[bTagNr] = (LightJetNr[bTagNr])[2];} 

  return lightIndices;
}

void BTagStudy::CompareJetCombi(vector<int> jetCombi, int OptionNr, int NrJets, int lightOne, int lightTwo){

  //jet combi order is : (0 = BLeptonic, 1 = BHadronic, 2 = Quark1 & 3 =  Quark2)
  if(jetCombi[1] == 9999 || jetCombi[0] == 9999 || jetCombi[2] == 9999 || jetCombi[3] == 9999){
    NotReconstructedEvent[OptionNr][NrJets]++;
  }
  else{
    if( jetCombi[0] == bLeptIndex[OptionNr]  && jetCombi[1] == bHadrIndex[OptionNr]   &&
       (jetCombi[2] == lightOne || jetCombi[2] == lightTwo) &&
       (jetCombi[3] == lightOne || jetCombi[3] == lightTwo) ){
        allFourJetsCorrectlyMatched[OptionNr][NrJets]++;
    }
    if( jetCombi[0] != bLeptIndex[OptionNr]  || jetCombi[1] != bHadrIndex[OptionNr]   ||
       (jetCombi[2] != lightOne && jetCombi[2] != lightTwo) ||
       (jetCombi[3] != lightOne && jetCombi[3] != lightTwo) ){
	atLeastOneWronglyMatched[OptionNr][NrJets]++;
    }
    if(jetCombi[0] == bLeptIndex[OptionNr] && jetCombi[1] == bHadrIndex[OptionNr])
      twoBTagsCorrectlyMatched[OptionNr][NrJets]++;

    if(jetCombi[0] != bLeptIndex[OptionNr] || jetCombi[1] != bHadrIndex[OptionNr] )
      atLeastOneBTagWronglyMatched[OptionNr][NrJets]++;

    if((jetCombi[2] == lightOne || jetCombi[2] == lightTwo) &&
       (jetCombi[3] == lightOne || jetCombi[3] == lightTwo) )
	twoLightJetsCorrectlyMatched[OptionNr][NrJets]++;
	
    if( (jetCombi[2] != lightOne && jetCombi[2] != lightTwo) ||
        (jetCombi[3] != lightOne && jetCombi[3] != lightTwo) )
	atLeastOneLightJetWronglyMatched[OptionNr][NrJets]++;	
  }
}


/*void BTagStudy::CorrectJetCombi(vector<int> jetCombi, int OptionNr){                           //Maybe add an additional integer to differentiate between 4 and 5 jet case!
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
*/
void BTagStudy::ReturnBTagTable(){ 

    string Title[3]= {"   \\textbf{Option} & all 4 correct & $\\geq$ 1 wrong       & correct ($\\%$)       & $\\frac{s}{b}$ & non-matched \\\\", 
		      "   \\textbf{Option} & 2 b's correct & $\\geq$ 1 b wrong     & b's correct ($\\%$)   & $\\frac{s}{b}$ & non-matched \\\\", 
		      "   \\textbf{Option} & 2 light good  & $\\geq$ 1 light wrong & light correct ($\\%$) & $\\frac{s}{b}$ & non-matched \\\\"};

    string Caption[3] = {"   \\caption{Overview of correct and wrong reconstructed events for the different b-tags (no $\\chi^{2}$ $m_{lb}$ - $m_{qqb}$ applied)} ", 
			 "   \\caption{Overview of correct and wrong reconstructed b-jets for the different b-tags (no $\\chi^{2}$ $m_{lb}$ - $m_{qqb}$ applied)} ", 
			 "   \\caption{Overview of correct and wrong reconstructed light jets for the different b-tags (no $\\chi^{2}$ $m_{lb}$ - $m_{qqb}$ applied)} "};

    string Title5Jets[3]= {"   \\textbf{Option} & all 4 correct & $\\geq$ 1 wrong       & correct ($\\%$)       & $\\frac{s}{b}$ & $3^{rd}$ jet = light & $3^{rd}$ jet gain($\\%s$) \\\\",
			   "   \\textbf{Option} & 2 b's correct & $\\geq$ 1 b wrong     & b's correct ($\\%$)   & $\\frac{s}{b}$ &                      &                           \\\\", 
			   "   \\textbf{Option} & 2 light good  & $\\geq$ 1 light wrong & light correct ($\\%$) & $\\frac{s}{b}$ &                      &                           \\\\"};

    string Caption5Jets[3] = {"   \\caption{Overview of correct and wrong reconstructed events for the different b-tags (no $\\chi^{2}$ $m_{lb}$ - $m_{qqb}$ applied, 5 jets considered)} ", 
			      "   \\caption{Overview of correct and wrong reconstructed b-jets for the different b-tags (no $\\chi^{2}$ $m_{lb}$ - $m_{qqb}$ applied, 5 jets considered)} ", 
			      "   \\caption{Overview of correct and wrong reconstructed light jets for the different b-tags (no $\\chi^{2}$ $m_{lb}$ - $m_{qqb}$ applied, 5 jets considered)} "};

    for(int itNrJets = 0; itNrJets < 3; itNrJets++){      //Specifying whether all jets, only b-jets or only light jets are considered!

        for(int itBTag = 0; itBTag < 6; itBTag++){

            /*//---  Store 5 jet information  ---//
            if(itBTag == 0) evtSelOutput5Jets << " \\begin{table}[!h] \n  \\begin{tabular}{c|c|c|c|c|c|c} \n " << Title5Jets[itNrJets] << " \\hline " << endl;

            int CorrectOnes5Jets[3] = {allFourJetsCorrectlyMatched5Jets[itBTag], twoBTagsCorrectlyMatched5Jets[itBTag] ,    twoLightJetsCorrectlyMatched5Jets[itBTag]};
            int WrongOnes5Jets[3]   = {atLeastOneWronglyMatched5Jets[itBTag],    atLeastOneBTagWronglyMatched5Jets[itBTag], atLeastOneLightJetWronglyMatched5Jets[itBTag]};
            int Correct3rdJet[3]    = {thirdJetIsActualQuark[itBTag],            0,                                         0 };
	    float sOverSqrtB5Jets[3], sOverB5Jets[3], CorrectPercentage5Jets[3], ThirdJetPercentage[3];
	    for(int jj = 0; jj < 3; jj++){
	        sOverSqrtB5Jets[jj] = (float)(CorrectOnes5Jets[jj])/(float)(sqrt(WrongOnes5Jets[jj]));
	        sOverB5Jets[jj] = (float)(CorrectOnes5Jets[jj])/(float)(WrongOnes5Jets[jj]);
	        CorrectPercentage5Jets[jj] = (float)(CorrectOnes5Jets[jj]*100.0)/(float)(CorrectOnes5Jets[jj]+WrongOnes5Jets[jj]);
	        ThirdJetPercentage[jj] = (float)(Correct3rdJet[jj]*100.0)/(float)(CorrectOnes5Jets[jj]);
	    }
	
            evtSelOutput5Jets << OptionName[itBTag] << 
	    " & " << CorrectOnes5Jets[itNrJets]          << 
	    " & " << WrongOnes5Jets[itNrJets]            << 
	    " & " << CorrectPercentage5Jets[itNrJets]    << 
	    " & " << sOverB5Jets[itNrJets]               << 
	    //" & " << NotReconstructedEvent5Jets[ii]      << 
	    " & " << Correct3rdJet[itNrJets]             << 
	    " & " << ThirdJetPercentage[itNrJets]        << 
	    "\\\\ " << endl;

            if(itNrJets == 0) std::cout << " Value of not-reconstructed events (for 5 jets) : " << NotReconstructedEvent5Jets[itBTag] << " -- for option " << OptionName[itBTag] << std::endl;
            if(itBTag == 5){
                evtSelOutput5Jets << "  \\end{tabular} " << endl;
                evtSelOutput5Jets << Caption5Jets[itNrJets] << endl;
                evtSelOutput5Jets << " \\end{table} \n " << endl;
            }*/
            
            //---  Store 4 jet information  ---//         
            if(itBTag == 0) evtSelOutput << " \\begin{table}[!h] \n  \\begin{tabular}{c|c|c|c|c|c} \n " << Title[itNrJets] << " \\hline " << endl;

    	    int CorrectOnes[3] = {allFourJetsCorrectlyMatched[itBTag][0], twoBTagsCorrectlyMatched[itBTag][0]    , twoLightJetsCorrectlyMatched[itBTag][0]};
	    int WrongOnes[3]   = {atLeastOneWronglyMatched[itBTag][0],    atLeastOneBTagWronglyMatched[itBTag][0], atLeastOneLightJetWronglyMatched[itBTag][0]};
            float sOverSqrtB[3], sOverB[3], CorrectPercentage[3];
            for(int jj = 0; jj < 3; jj++){
                sOverSqrtB[jj] = (float)(CorrectOnes[jj])/(float)(sqrt(WrongOnes[jj]));
                sOverB[jj] = (float)(CorrectOnes[jj])/(float)(WrongOnes[jj]);
                CorrectPercentage[jj] = (float)(CorrectOnes[jj]*100.0)/(float)(CorrectOnes[jj]+WrongOnes[jj]);
            }

            evtSelOutput << OptionName[itBTag] << 
	    " & " << CorrectOnes[itNrJets]          << 
	    " & " << WrongOnes[itNrJets]            << 
	    " & " << CorrectPercentage[itNrJets]    << 
	    " & " << sOverB[itNrJets]               << 
	    " & " << NotReconstructedEvent[itBTag][0] << 
	    " \\\\ " << endl;

            if(itBTag == 5){
                evtSelOutput << "  \\end{tabular} " << endl;
                evtSelOutput << Caption[itNrJets] << endl;
                evtSelOutput << " \\end{table} \n " << endl;
            }
        }            	
    }        
}
