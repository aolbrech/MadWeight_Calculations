#ifndef BTagStudy_h
#define BTagStudy_h

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

#include "TopTreeProducer/interface/TRootJet.h"

using namespace std;
using namespace TopTree;

class BTagStudy{

  public:
    BTagStudy(int);
    ~BTagStudy();

    void CalculateJets(vector<TRootJet*>, vector<int> jetCombi, TLorentzVector* lepton);
    void ReturnBTagTable();

    int getBHadrIndex(int bTagNr)  {return bHadrIndex[bTagNr];};
    int getBLeptIndex(int bTagNr)  {return bLeptIndex[bTagNr];};
    int getLight1Index4Jets(int bTagNr) {return light1Index4Jets[bTagNr];};   //Need to think how to return the light jets ... (not possible to differentiate between 4 and 5 jet case)
    int getLight2Index4Jets(int bTagNr) {return light2Index4Jets[bTagNr];};
    int getLight1Index5Jets(int bTagNr) {return light1Index5Jets[bTagNr];};  
    int getLight2Index5Jets(int bTagNr) {return light2Index5Jets[bTagNr];};

    float Mlb, Mqqb, SigmaMlb, SigmaMqqb;

    //Functions to ask everything inside the analyzer code!                          -->Is this getter still necessary?
    vector<int> getbTaggedJets(int OptionNr)   {return bTaggedJetNr[OptionNr];};
    vector<int> getNonbTaggedJets(int OptionNr){return NonbTaggedJetNr[OptionNr];};
    vector<int> getLightJets(int OptionNr)     {return LightJetNr[OptionNr];};

  private:
    void CompareJetCombi(vector<int> jetCombi, int bTagNr, int jetCase, int lightJetOne, int lightJetTwo);
    //void CorrectJetCombi(vector<int> jetCombi, int);
    //void CorrectJetCombi5Jets(vector<int> jetCombi, int);
    void InitializeBegin();
    void ResetEventArrays();
    void CalculateMlbChiSq( int bTagNr, TLorentzVector* lepton, vector<TRootJet*> Jets); 
    vector<int> CalculateMqqbChiSq( int bTagNr, vector<TRootJet*> Jets); 

    int LowestChiSqMlb[6], LowestChiSqMqqb[6];
    int bHadrIndex[6], bLeptIndex[6], light1Index4Jets[6], light2Index4Jets[6], light1Index5Jets[6], light2Index5Jets[6], light1IndexMqqb[6], light2IndexMqqb[6];
    vector<float> ChiSquaredMlb[6], ChiSquaredMqqb[6];

    int verbose;
    float BJetWP[6], LightJetWP[6];  
    std::string OptionName[6];
    ofstream evtSelOutput, evtSelOutput5Jets;
    vector<int> bTaggedJetNr[6], NonbTaggedJetNr[6], LightJetNr[6];

    //4- and 5-jet case comparison
    //int EventWithTwoLightJets[6];
    //int EventWithThreeLightJets[6];
    //int EventWithTwoLightJetsAndBTagged[6];
    //int EventWithThreeLightJetsAndBTagged[6];
    int thirdJetIsActualQuark[6][2], secondJetIsActualQuark[6][2], firstJetIsActualQuark[6][2];
    //int thirdJetIsCorrectQuark[6];
    //int thirdJetIsGoodQuark[6];    //3rd quark is one of the quarks when the two b-jets are correctly matched!

    int NotReconstructedEvent[6][2];
    int allFourJetsCorrectlyMatched[6][2];
    int atLeastOneWronglyMatched[6][2];
    int twoBTagsCorrectlyMatched[6][2];
    int atLeastOneBTagWronglyMatched[6][2];
    int twoLightJetsCorrectlyMatched[6][2];
    int atLeastOneLightJetWronglyMatched[6][2];
  
    /*int getNrEventsWithTwoLightJets(int OptionNr)             {return EventWithTwoLightJets[OptionNr];};
    int getNrEventsWithThreeLightJets(int OptionNr)           {return EventWithThreeLightJets[OptionNr];};
    int getNrEventsWithTwoLightJetsAndBTagged(int OptionNr)   {return EventWithTwoLightJetsAndBTagged[OptionNr];};
    int getNrEventsWithThreeLightJetsAndBTagged(int OptionNr) {return EventWithThreeLightJetsAndBTagged[OptionNr];};
    int getNrTimesThirdJetIsActualQuark(int OptionNr)         {return thirdJetIsActualQuark[OptionNr];};
    int getNrTimesThirdJetIsCorrectQuark(int OptionNr)         {return thirdJetIsCorrectQuark[OptionNr];};

    int getNrNotFoundEvts(int OptionNr)       {return NotReconstructedEvent[OptionNr];};
    int getNrCorrectMatchedEvts(int OptionNr) {return allFourJetsCorrectlyMatched[OptionNr];};
    int getNrWrongMatchedEvts(int OptionNr)   {return atLeastOneWronglyMatched[OptionNr];};
    int getNrCorrectBJetMatchedEvts(int OptionNr) {return twoBTagsCorrectlyMatched[OptionNr];};
    int getNrWrongBJetMatchedEvts(int OptionNr)   {return atLeastOneBTagWronglyMatched[OptionNr];};
    int getNrCorrectLightJetMatchedEvts(int OptionNr) {return twoLightJetsCorrectlyMatched[OptionNr];};
    int getNrWrongLightJetMatchedEvts(int OptionNr)   {return atLeastOneLightJetWronglyMatched[OptionNr];};

    float getSignalOverSqrtBkg(int OptionNr) {return allFourJetsCorrectlyMatched[OptionNr]/sqrt(atLeastOneWronglyMatched[OptionNr]);};
    float getSignalOverBkg(int OptionNr)     {return (float)allFourJetsCorrectlyMatched[OptionNr]/(float)atLeastOneWronglyMatched[OptionNr];};
    */

    //
    // 5-jet case:
    int NotReconstructedEvent5Jets[6];
    int allFourJetsCorrectlyMatched5Jets[6];
    int atLeastOneWronglyMatched5Jets[6];
    int twoBTagsCorrectlyMatched5Jets[6];
    int atLeastOneBTagWronglyMatched5Jets[6];
    int twoLightJetsCorrectlyMatched5Jets[6];
    int atLeastOneLightJetWronglyMatched5Jets[6];
  
    /*int getNrNotFoundEvts5Jets(int OptionNr)       {return NotReconstructedEvent5Jets[OptionNr];};
    int getNrCorrectMatchedEvts5Jets(int OptionNr) {return allFourJetsCorrectlyMatched5Jets[OptionNr];};
    int getNrWrongMatchedEvts5Jets(int OptionNr)   {return atLeastOneWronglyMatched5Jets[OptionNr];};
    int getNrCorrectBJetMatchedEvts5Jets(int OptionNr) {return twoBTagsCorrectlyMatched5Jets[OptionNr];};
    int getNrWrongBJetMatchedEvts5Jets(int OptionNr)   {return atLeastOneBTagWronglyMatched5Jets[OptionNr];};
    int getNrCorrectLightJetMatchedEvts5Jets(int OptionNr) {return twoLightJetsCorrectlyMatched5Jets[OptionNr];};
    int getNrWrongLightJetMatchedEvts5Jets(int OptionNr)   {return atLeastOneLightJetWronglyMatched5Jets[OptionNr];};

    float getSignalOverSqrtBkg5Jets(int OptionNr) {return allFourJetsCorrectlyMatched5Jets[OptionNr]/sqrt(atLeastOneWronglyMatched5Jets[OptionNr]);};
    float getSignalOverBkg5Jets(int OptionNr)     {return (float)allFourJetsCorrectlyMatched5Jets[OptionNr]/(float)atLeastOneWronglyMatched5Jets[OptionNr];}; 
    */
};

#endif
