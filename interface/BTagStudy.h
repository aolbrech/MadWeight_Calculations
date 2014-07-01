#ifndef BTagStudy_h
#define BTagStudy_h

#include <string>
#include <vector>
#include <iostream>
#include <math.h>

#include <fstream>
#include <sstream>

//using namespace std;

class BTagStudy{

  public:
   //BTagStudy();
   //~BTagStudy();

   void InitializePerEvent();
   void InitializeBegin();
   void CalculateJets(vector<TRootJet*>, float BTagWorkingPoint, float LightWorkingPoint, int OptionNr);
   void CorrectJetCombi(int, int, int, int, int);
   void CorrectJetCombi5Jets(int, int, int, int, int);
   void ReturnTable(std::string NameOfOption4Jets[6], std::string NameOfOption5Jets[6], int WhichJets, int NrOptionsConsidered, ofstream &output, int OptionOfInterest);//Automatically writes everything out into a .tex file!

   vector<int> bTaggedJetNr[6];
   vector<int> NonbTaggedJetNr[6];
   vector<int> LightJetNr[6];

   //4- and 5-jet case comparison
   int EventWithTwoLightJets[6];
   int EventWithThreeLightJets[6];
   int EventWithTwoLightJetsAndBTagged[6];
   int EventWithThreeLightJetsAndBTagged[6];
   int thirdJetIsActualQuark[6];
   int thirdJetIsCorrectQuark[6];
   int thirdJetIsGoodQuark[6];    //3rd quark is one of the quarks when the two b-jets are correctly matched!

   int NotReconstructedEvent[6];
   int allFourJetsCorrectlyMatched[6];
   int atLeastOneWronglyMatched[6];
   int twoBTagsCorrectlyMatched[6];
   int atLeastOneBTagWronglyMatched[6];
   int twoLightJetsCorrectlyMatched[6];
   int atLeastOneLightJetWronglyMatched[6];

   //Functions to ask everything inside the analyzer code!
   vector<int> getbTaggedJets(int OptionNr)   {return bTaggedJetNr[OptionNr];};
   vector<int> getNonbTaggedJets(int OptionNr){return NonbTaggedJetNr[OptionNr];};
   vector<int> getLightJets(int OptionNr)     {return LightJetNr[OptionNr];};

   int getNrEventsWithTwoLightJets(int OptionNr)             {return EventWithTwoLightJets[OptionNr];};
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

   //
   // 5-jet case:
   int NotReconstructedEvent5Jets[6];
   int allFourJetsCorrectlyMatched5Jets[6];
   int atLeastOneWronglyMatched5Jets[6];
   int twoBTagsCorrectlyMatched5Jets[6];
   int atLeastOneBTagWronglyMatched5Jets[6];
   int twoLightJetsCorrectlyMatched5Jets[6];
   int atLeastOneLightJetWronglyMatched5Jets[6];

   int getNrNotFoundEvts5Jets(int OptionNr)       {return NotReconstructedEvent5Jets[OptionNr];};
   int getNrCorrectMatchedEvts5Jets(int OptionNr) {return allFourJetsCorrectlyMatched5Jets[OptionNr];};
   int getNrWrongMatchedEvts5Jets(int OptionNr)   {return atLeastOneWronglyMatched5Jets[OptionNr];};
   int getNrCorrectBJetMatchedEvts5Jets(int OptionNr) {return twoBTagsCorrectlyMatched5Jets[OptionNr];};
   int getNrWrongBJetMatchedEvts5Jets(int OptionNr)   {return atLeastOneBTagWronglyMatched5Jets[OptionNr];};
   int getNrCorrectLightJetMatchedEvts5Jets(int OptionNr) {return twoLightJetsCorrectlyMatched5Jets[OptionNr];};
   int getNrWrongLightJetMatchedEvts5Jets(int OptionNr)   {return atLeastOneLightJetWronglyMatched5Jets[OptionNr];};

   float getSignalOverSqrtBkg5Jets(int OptionNr) {return allFourJetsCorrectlyMatched5Jets[OptionNr]/sqrt(atLeastOneWronglyMatched5Jets[OptionNr]);};
   float getSignalOverBkg5Jets(int OptionNr)     {return (float)allFourJetsCorrectlyMatched5Jets[OptionNr]/(float)atLeastOneWronglyMatched5Jets[OptionNr];};


};

#endif
