#ifndef MlbStudy_h
#define MlbStudy_h

#include <string>
#include <vector>
#include <iostream>
#include <math.h>

#include <fstream>
#include <sstream>
#include "TH1.h"
#include "TFile.h"

//using namespace std;

class MlbStudy{

  public:
   //MLbStudy();
   //~MlbStudy();

   void calculateChiSquared(vector<int> CorrectValues, vector<int> bTaggedJets, vector<int> lightJets, TLorentzVector* lepton, vector<TRootJet*> Jets, float MassMlb, float SigmaMlb, float MassMqqb, float SigmaMqqb);
   void calculateEfficiency(int option, vector<int> CorrectValues, vector<int> bTaggedJets, vector<int> lightJets);
   void saveNumbers(std::string NameOfOption[6], int WhichJets, int NrOptionsConsidered, ofstream &output);
   void initializePerEvent();
   void initializeBegin();

   int chosenBHadr, chosenBLept, chosenQuark1, chosenQuark2;
   int LowestChiSq;                //Since Mqqb Values are set artificially high, only the first two chi-sq will be the lowest one in the case of only 4 jets!

   float MlbValues[6], MqqbValues[6], ChiSquared[6];

   vector<float> h_CorrectOptionChiSq;
   vector<float> h_WrongOptionChiSq;

   int NumberMatchedEvents[6];
   int NumberNotMatchedEvents[6];
   int CorrectOptionAvailable[6];
   int CorrectOptionChosen[6];
   int WrongOptionChosen[6];
   int CorrectEventMlbMqqb[6];
   int WrongEventMlbMqqb[6];

   int getChosenBLept() {return chosenBLept;};
   int getChosenBHadr() {return chosenBHadr;};
   int getChosenQuark1(){return chosenQuark1;};
   int getChosenQuark2(){return chosenQuark2;};
   vector<float> getHistoCorrectOptionChiSq() {return h_CorrectOptionChiSq;};
   vector<float> getHistoWrongOptionChiSq()   {return h_WrongOptionChiSq;};

   int getMlbValue(int order)     {return MlbValues[order];};
   int getMqqbValue(int order)    {return MqqbValues[order];};
   float getChiSquared(int order) {return ChiSquared[order];};

   int getNrMatchedEvents( int OptionNr) {return NumberMatchedEvents[OptionNr];};
   int getNrCorrectOptionsAvailable(int OptionNr) {return CorrectOptionAvailable[OptionNr];};
   int getNrCorrectOptionChosen(int OptionNr) {return CorrectOptionChosen[OptionNr];};
   int getNrWrongOptionChosen(int OptionNr) {return WrongOptionChosen[OptionNr];};
};

#endif
