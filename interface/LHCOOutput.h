#ifndef LHCOOutput_h
#define LHCOOutput_h

#include <fstream>
#include <sstream>
#include "TopTreeProducer/interface/TRootMCParticle.h"
#include "TLorentzVector.h"

using namespace std;

class LHCOOutput{

  int LeptonCharge;
public:
  void LHCOEventOutput(int LHCOIndex, ofstream &outputFile, unsigned int EventNumber, std::vector<TRootMCParticle*> vector, std::vector<int> MGId); //TRootMCParticle info needed?
  void LHCOEventRecoOutput(int LHCOIndex, ofstream &outputFile, unsigned int EventNumber, std::vector<TLorentzVector*> vector, std::vector<int> MGId); //TRootMCParticle info needed?

};

#endif
