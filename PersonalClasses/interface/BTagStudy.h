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
    BTagStudy(int, vector<Dataset*>, bool, int, float, float, float, float);
    ~BTagStudy();

    void InitializeDataSet(std::string);
    void CalculateJets(vector<TLorentzVector>, vector<float> bTagValues, vector<int> jetCombi, TLorentzVector lepton, Dataset*, float);
    void ReturnBTagTable();
    //void ReturnThirdJetTable();    //Still to fill
    void CreateHistograms(TFile*, bool, std::string, int iDS); 

    vector<int> getIndices(int bTagNr) {return jetIndices[bTagNr];}; 
    int getBHadrIndex(int bTagNr)  {return bHadrIndex[bTagNr];};
    int getBLeptIndex(int bTagNr)  {return bLeptIndex[bTagNr];};
    int getLight1Index(int bTagNr) {return light1Index[bTagNr];};
    int getLight2Index(int bTagNr) {return light2Index[bTagNr];};

    int getNrBTaggedJets(int OptionNr)   {return bTagJetNr[OptionNr].size();};
    int getNrNonbTaggedJets(int OptionNr){return NonbTagJetNr[OptionNr].size();};
    int getNrLightJets(int OptionNr)     {return LightJetNr[OptionNr].size();};

    vector<int> getLightJets(int OptionNr) {return LightJetNr[OptionNr];};

    float getMlbMqqbChiSq(int bTagNr) {return LowestChiSq[bTagNr];};

  private:
    void CompareJetCombi(vector<int> jetCombi, int bTagNr, bool jetCase, int lightJetOne, int lightJetTwo);
    void ResetEventArrays();
    int getLowestMlbMqqbChiSquared( int bTagNr, vector<TLorentzVector> Jets, TLorentzVector lepton);

    int bHadrIndex[9], bLeptIndex[9], light1Index[9], light2Index[9];
    vector<int> jetIndices[9];
    int NotReconstructedEvent[2][9], CorrectlyMatched[2][3][9], atLeastOneWrongMatch[2][3][9];
    float LowestChiSq[9];

    int verbose;
    float BJetWP[9], LightJetWP[9];  
    std::string OptionName[9], BTitle[9], BName[9];
    ofstream evtSelOutput[2];
    vector<int> bTagJetNr[9], NonbTagJetNr[9], LightJetNr[9];

    float Mlb, Mqqb, S_Mlb, S_Mqqb;
    int chosenBTag_, nrBTags_, nrDatasets_;
    bool singleWP_, use5Jets_;  
    std::string dataSetName_;

    map<std::string,TH1F*> histo1D;
    map<string,TH2F*> histo2D;
    map<string,MultiSamplePlot*> MSPlot;

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
