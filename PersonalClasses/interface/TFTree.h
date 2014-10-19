#ifndef TFTree_h
#define TFTree_h

#include "TObject.h"
#include "Rtypes.h"

using namespace std;

class TFTree : public TObject
{

public:

TFTree():
    TObject()
    ,eventID_(0)
    ,EtaReco_()
    {;}

    ~TFTree() {;}

    unsigned int eventID() const { return eventID_; }
    float EtaReco() {return EtaReco_;}

    void setEventID(unsigned int eventID) { eventID_ = eventID; }
    void setEtaReco(float EtaReco) {EtaReco_ = EtaReco;}

protected:
    unsigned int eventID_;
    float EtaReco_;

ClassDef (TFTree,2);

};

#endif
