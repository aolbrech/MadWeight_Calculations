#ifndef TFnTuple_h
#define TFnTuple_h

#include "TObject.h"
#include "Rtypes.h"
#include "TLorentzVector.h"

using namespace std;

class TFnTuple : public TObject
{

public:

TFnTuple():
    TObject()
    ,eventID_(0)
    ,recoVectorLight1_()
    ,recoVectorLight2_()
    ,recoVectorHadrB_()
    ,recoVectorLeptB_()
    ,recoVectorLepton_()
    ,genVectorLight1_()
    ,genVectorLight2_()
    ,genVectorHadrB_()
    ,genVectorLeptB_()
    ,genVectorLepton_()
    {;}

    ~TFnTuple() {;}

    unsigned int eventID() const { return eventID_; }
    TLorentzVector recoVectorLight1() {return recoVectorLight1_;}
    TLorentzVector recoVectorLight2() {return recoVectorLight2_;}
    TLorentzVector recoVectorHadrB()  {return recoVectorHadrB_;}
    TLorentzVector recoVectorLeptB()  {return recoVectorLeptB_;}
    TLorentzVector recoVectorLepton() {return recoVectorLepton_;}
    TLorentzVector genVectorLight1()  {return genVectorLight1_;}
    TLorentzVector genVectorLight2()  {return genVectorLight2_;}
    TLorentzVector genVectorHadrB()   {return genVectorHadrB_;}
    TLorentzVector genVectorLeptB()   {return genVectorLeptB_;}
    TLorentzVector genVectorLepton()  {return genVectorLepton_;}

    void setEventID(unsigned int eventID) { eventID_ = eventID; }
    void setRecoVectorLight1( TLorentzVector recoVectorLight1) {recoVectorLight1_ = recoVectorLight1;}
    void setRecoVectorLight2( TLorentzVector recoVectorLight2) {recoVectorLight2_ = recoVectorLight2;}
    void setRecoVectorHadrB(  TLorentzVector recoVectorHadrB)  {recoVectorHadrB_  = recoVectorHadrB;}
    void setRecoVectorLeptB(  TLorentzVector recoVectorLeptB)  {recoVectorLeptB_  = recoVectorLeptB;}
    void setRecoVectorLepton( TLorentzVector recoVectorLepton) {recoVectorLepton_ = recoVectorLepton;}
    void setGenVectorLight1(  TLorentzVector genVectorLight1)  {genVectorLight1_  = genVectorLight1;}
    void setGenVectorLight2(  TLorentzVector genVectorLight2)  {genVectorLight2_  = genVectorLight2;}
    void setGenVectorHadrB(   TLorentzVector genVectorHadrB)   {genVectorHadrB_   = genVectorHadrB;}
    void setGenVectorLeptB(   TLorentzVector genVectorLeptB)   {genVectorLeptB_   = genVectorLeptB;}
    void setGenVectorLepton(  TLorentzVector genVectorLepton)  {genVectorLepton_  = genVectorLepton;}

protected:
    unsigned int eventID_;
    TLorentzVector recoVectorLight1_;
    TLorentzVector recoVectorLight2_;
    TLorentzVector recoVectorHadrB_;
    TLorentzVector recoVectorLeptB_;
    TLorentzVector recoVectorLepton_;
    TLorentzVector genVectorLight1_;
    TLorentzVector genVectorLight2_;
    TLorentzVector genVectorHadrB_;
    TLorentzVector genVectorLeptB_;
    TLorentzVector genVectorLepton_;

ClassDef (TFnTuple,2);

};

#endif
