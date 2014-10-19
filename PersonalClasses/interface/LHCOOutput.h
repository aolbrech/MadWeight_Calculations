#ifndef LHCOOutput_h
#define LHCOOutput_h

#include <iostream>
#include "/user/aolbrech/GitTopTree_Feb2014/TopBrussels/TopTreeProducer/interface/TRootMCParticle.h"
#include "TLorentzVector.h"

using namespace std;
using namespace TopTree;

class LHCOOutput{

  int LeptonCharge;
public:
  void LHCOEventOutput(int LHCOIndex, ostream &outputFile, unsigned int EventNumber, std::vector<TRootMCParticle*> vector, std::vector<int> MGId); //TRootMCParticle info needed?
  void LHCOEventRecoOutput(int LHCOIndex, ostream &outputFile, unsigned int EventNumber, std::vector<TLorentzVector*> vector, std::vector<int> MGId); //TRootMCParticle info needed?

};

#endif
