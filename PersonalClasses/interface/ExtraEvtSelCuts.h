#ifndef ExtraEvtSelCuts_h
#define ExtraEvtSelCuts_h

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include "TopTreeProducer/interface/TRootJet.h"
//#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
//#include "TopTreeAnalysisBase/Content/interface/Dataset.h"

using namespace std;
using namespace TopTree;

class ExtraEvtSelCuts{

  public:
    ExtraEvtSelCuts();
    ~ExtraEvtSelCuts();

    void Initialize(float TopMassHadr, float sTopMassHadr, float WMassHadr, float sWMassHadr, std::string bTagTitle);
    void KeepEvent(vector<int>, TLorentzVector, vector<TLorentzVector>, vector<int>, float);
    void StoreCutInfluence(TFile*);

  private:
    float MTopHadr_, MWHadr_, sig_MTopH_, sig_MWH_;
    float ChiSqCutVal_[5], MTop_Down_[5], MTop_Up_[5]; //MassWindow_[5];
    int nrCuts_;
    std::string bTitle_;

    int CorrEvts_, WrongEvts_, UnmatchEvts_;
    int CorrEvts_ChiSq_[5], WrongEvts_ChiSq_[5], UnmatchEvts_ChiSq_[5];
    int CorrEvts_MassWindow_[5], WrongEvts_MassWindow_[5], UnmatchEvts_MassWindow_[5];
    int CorrEvts_Comb_[5], WrongEvts_Comb_[5], UnmatchEvts_Comb_[5];

    map<string,TH1F*> histo1D;
    map<string,TH2F*> histo2D;
};

#endif
