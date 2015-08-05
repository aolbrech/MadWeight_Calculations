#ifndef AnomCoupLight_h
#define AnomCoupLight_h

#include "TObject.h"
#include "Rtypes.h"
#include "TLorentzVector.h"

using namespace std;

class AnomCoupLight : public TObject
{

  public:
    AnomCoupLight():
      TObject()
      ,eventID_(0)
      ,runID_(0)
      ,lumiBlockID_(0)
      ,nPV_(0)
      ,nTruePU_(0)
      ,scaleFactor_(0)
      ,selectedJets_()
      ,CSVbTag_()
      ,selectedLepton_()
      ,decayChannel_(0)
      ,leptonCharge_(0)
      ,hadrBJet_(0)
      ,leptBJet_(0)
      ,quark1_(0)
      ,quark2_(0)
      {;}
    ~AnomCoupLight() {;}

    unsigned int eventID()     const {return eventID_;}
    unsigned int runID()       const {return runID_;}
    unsigned int lumiBlockID() const {return lumiBlockID_;}
    unsigned int nPV()         const {return nPV_;}
    unsigned int nTruePU()     const {return nTruePU_;}
    unsigned int scaleFactor() const {return scaleFactor_;}

    vector<TLorentzVector> selectedJets() const {return selectedJets_;}
    vector<float> CSVbTag() const {return CSVbTag_;}
    TLorentzVector selectedLepton() const {return selectedLepton_;}
    unsigned int decayChannel() const {return decayChannel_;}
    float leptonCharge() const {return leptonCharge_;}
    int hadrBJet() const {return hadrBJet_;}
    int leptBJet() const {return leptBJet_;}
    int quark1() const {return quark1_;}
    int quark2() const {return quark2_;} 

    void setEventId(unsigned int eventID) {eventID_ = eventID;}
    void setRunId(unsigned int runID) {runID_ = runID;}
    void setLumiBlockId(unsigned int lumiBlockID) {lumiBlockID_ = lumiBlockID;}
    void setNPV(unsigned int nPV) {nPV_ = nPV;}
    void setNTruePU(unsigned int nTruePU) {nTruePU_ = nTruePU;}
    void setScaleFactor(unsigned int scaleFactor) {scaleFactor_ = scaleFactor;}

    void setSelectedJets(vector<TLorentzVector> selectedJets) {selectedJets_ = selectedJets;}
    void setBTagCSV(vector<float> CSVbTag) {CSVbTag_ = CSVbTag;}
    void setSelectedLepton( TLorentzVector selectedLepton) {selectedLepton_ = selectedLepton;}
    void setDecayChannel(unsigned int decayChannel) {decayChannel_ = decayChannel;}
    void setLeptonCharge( float leptonCharge) {leptonCharge_ = leptonCharge;}
    void setHadrBJet( int hadrBJet) {hadrBJet_ = hadrBJet;}
    void setLeptBJet( int leptBJet) {leptBJet_ = leptBJet;}
    void setQuark1(int quark1) {quark1_ = quark1;}
    void setQuark2(int quark2) {quark2_ = quark2;}

  protected:
    unsigned int eventID_;
    unsigned int runID_;
    unsigned int lumiBlockID_;
    unsigned int nPV_;
    unsigned int nTruePU_;
    unsigned int scaleFactor_;

    vector<TLorentzVector> selectedJets_;
    vector<float> CSVbTag_;
    TLorentzVector selectedLepton_;
    unsigned int decayChannel_;
    float leptonCharge_;
    int hadrBJet_;
    int leptBJet_;
    int quark1_;
    int quark2_;

  ClassDef(AnomCoupLight,2); 
};

#endif
