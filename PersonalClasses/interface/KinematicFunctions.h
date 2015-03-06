#ifndef KinematicFunctions_h
#define KinematicFunctions_h

#include "TLorentzVector.h"

using namespace std;
//using namespace TopTree;

class KinematicFunctions{

 public:
    float CosTheta(TLorentzVector*, TLorentzVector*, TLorentzVector*);

};

#endif
