#ifndef BTagStudy_h
#define BTagStudy_h

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

#include "TopTreeProducer/interface/TRootJet.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Content/interface/Dataset.h"

using namespace std;
using namespace TopTree;

struct sort_pred {
  bool operator()(const std::pair<int,double> &left, const std::pair<int,double> &right) {
    return left.second < right.second;
  }
};

class BTagStudy{

  public:
    BTagStudy(int, vector<Dataset*>, bool, int);
    ~BTagStudy();

    void CalculateJets(vector<TLorentzVector>, vector<float> bTagValues, vector<int> jetCombi, TLorentzVector lepton, Dataset*, float);
    void ReturnBTagTable(std::string);
    void ReturnThirdJetTable();    //Still to fill
    void CreateHistograms(TFile*, bool, std::string, int); 
 
    int getBHadrIndex(int bTagNr)  {return bHadrIndex[bTagNr];};
    int getBLeptIndex(int bTagNr)  {return bLeptIndex[bTagNr];};
    int getLight1Index(int bTagNr) {return light1Index[bTagNr];};
    int getLight2Index(int bTagNr) {return light2Index[bTagNr];};

    int getNrBTaggedJets(int OptionNr)   {return bTagJetNr[OptionNr].size();};
    int getNrNonbTaggedJets(int OptionNr){return NonbTagJetNr[OptionNr].size();};
    int getNrLightJets(int OptionNr)     {return LightJetNr[OptionNr].size();};

    int getMlbMqqbChiSquared(int bTagNr) {return LowestChiSq[bTagNr];};

  private:
    void CompareJetCombi(vector<int> jetCombi, int bTagNr, bool jetCase, int lightJetOne, int lightJetTwo);
    void InitializeBegin(vector<Dataset*>);
    void ResetEventArrays();
    int getLowestMlbMqqbChiSquared( int bTagNr, vector<TLorentzVector> Jets, TLorentzVector lepton);

    int bHadrIndex[6], bLeptIndex[6], light1Index[6], light2Index[6];
    int NotReconstructedEvent[2][6], CorrectlyMatched[2][3][6], atLeastOneWrongMatch[2][3][6];
    int LowestChiSq[6];

    int verbose;
    float BJetWP[6], LightJetWP[6];  
    std::string OptionName[6], BTitle[6], BName[6];
    ofstream evtSelOutput[2];
    vector<int> bTagJetNr[6], NonbTagJetNr[6], LightJetNr[6];

    float Mlb, Mqqb, S_Mlb, S_Mqqb;
    map<std::string,TH1F*> histo1D;
    map<string,TH2F*> histo2D;
    map<string,MultiSamplePlot*> MSPlot;
    int nrDatasets_, chosenBTag_, nrBTags_;
    bool singleWP_, use5Jets_;  

    //4- and 5-jet case comparison
    //int EventWithTwoLightJets[6];
    //int EventWithThreeLightJets[6];
    //int EventWithTwoLightJetsAndBTagged[6];
    //int EventWithThreeLightJetsAndBTagged[6];
    //int thirdJetIsActualQuark[6][2], secondJetIsActualQuark[6][2], firstJetIsActualQuark[6][2];
    //int thirdJetIsCorrectQuark[6];
    //int thirdJetIsGoodQuark[6];    //3rd quark is one of the quarks when the two b-jets are correctly matched!
};

#endif
