#ifndef LHCOOutput_h
#define LHCOOutput_h

#include <iostream>
#include <fstream>
#include <sstream>
#include "TopTreeProducer/interface/TRootMCParticle.h"
#include "TLorentzVector.h"

using namespace std;
using namespace TopTree;

class LHCOOutput{

  int LeptonCharge;
public:
  LHCOOutput();
  ~LHCOOutput();

  void LHCOEventOutput(int LHCOIndex, ostream &outputFile, unsigned int EventNumber, std::vector<TRootMCParticle*> vector, std::vector<int> MGId, std::vector<float> MGBtag);
  void LHCOEventRecoOutput(int LHCOIndex, ostream &outputFile, unsigned int EventNumber, std::vector<TLorentzVector*> vector, std::vector<int> MGId, std::vector<float> MGBtag);
  void StoreGenInfo(vector<TRootMCParticle*> mcParticles, bool GenLHCOOutput, int verbosity);

  bool GenEventContentCorrect()    {return CorrectGenEvtContent;};
  int getLeptonType()              {return leptonType;};
  TLorentzVector* getGenLeptTop()  {return GenLeptonicTop;}; 
  TLorentzVector* getGenLeptW()    {return GenLeptonicW;};
  TLorentzVector* getGenLepton()   {return GenLepton;};
  TLorentzVector* getGenNeutrino() {return GenNeutrino;}; 
  //TLorentzVector* getGenHadrTop()  {return GenHadronicTop;}; 
  //TLorentzVector* getGenHadrW()    {return GenHadronicW;}; 

private:
  TRootMCParticle *Top,*TopBar,*Bottom, *BottomBar,*Lepton,*NeutrinoMC,*WPlus,*WMinus,*Light,*LightBar;
  TLorentzVector *GenLeptonicTop, *GenLeptonicW, *GenLepton, *GenNeutrino;
  //TLorentzVector *GenHadronicTop, *GenHadronicW;
  unsigned int NumberNegativeElectrons, NumberNegativeMuons, NumberPositiveElectrons, NumberPositiveMuons;
  bool CorrectGenEvtContent;
  ofstream GenOutFile[4];

  enum LeptonType_t {muPlus, muMinus, elPlus, elMinus};
  LeptonType_t leptonType;
  //static std::string leptonTypeString[4];// = {"muPlus","muMinus","elPlus","elMinus"};
};

#endif
