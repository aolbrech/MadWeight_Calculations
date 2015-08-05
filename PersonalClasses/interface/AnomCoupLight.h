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
      {;}

    ~AnomCoupLight() {;}

    unsigned int eventID()     const {return eventID_;}
    unsigned int runID()       const {return runID_;}
    unsigned int lumiBlockID() const {return lumiBlockID_;}
    unsigned int nPV()         const {return nPV_;}
    unsigned int nTruePU       const {return nTruePU_;}
    unsigned int scaleFactor   const {return scaleFactor_;}

    void setEventId(unsigned int eventID) {eventID_ = eventID;}
    void setRunId(unsigned int runID) {runID_ = runID;}
    void setLumiBlockId(unsigned int lumiBlockID_ = lumiBlockID;}
    void setNPV(unsigned int nPV) {nPV_ = nPV;}
    void setNTruePU(unsigned int nTruePU) {nTruePU_ = nTruePU;}
    void setScaleFactor(unsigned int scaleFactor) {scaleFactor_ = scaleFactor;}

  protected:
    unsigned int eventID_;
    unsigned int runID_;
    unsigned int lumiBlockID_;
    unsigned int nPV_;
    unsigned int nTruePU_;
    unsigned int scaleFactor_;

  ClassDef(AnomCoupLight,2); 
};

#endif
