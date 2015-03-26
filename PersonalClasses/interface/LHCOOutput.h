#ifndef LHCOOutput_h
#define LHCOOutput_h

#include <iostream>
#include <fstream>
#include <sstream>
#include "TopTreeProducer/interface/TRootMCParticle.h"
#include "TopTreeProducer/interface/TRootJet.h"
#include "TLorentzVector.h"

using namespace std;
using namespace TopTree;

class LHCOOutput{

public:
  LHCOOutput(int, bool, bool);
  ~LHCOOutput();

  void StoreGenInfo(vector<TRootMCParticle*> mcParticles);
  void StoreRecoInfo(TLorentzVector* lepton, vector<TRootJet*> Jets,int bLept, int bHadr, int light1, int light2, int decayChannelEnum, float leptonCharge, vector<int> jetCombi); 

  bool GenEventContentCorrect()    {return CorrectGenEvtContent;};
  int getLeptonType()              {return leptonType;};
  TLorentzVector* getGenLeptTop()  {return GenLeptonicTop;}; 
  TLorentzVector* getGenLeptW()    {return GenLeptonicW;};
  TLorentzVector* getGenLepton()   {return GenLepton;};
  TLorentzVector* getGenNeutrino() {return GenNeutrino;}; 
  //TLorentzVector* getGenHadrTop()  {return GenHadronicTop;}; 
  //TLorentzVector* getGenHadrW()    {return GenHadronicW;}; 

private:
  void LHCOEventOutput(int LHCOIndex, ostream &outputFile, unsigned int EventNumber, std::vector<TLorentzVector*> vector, std::vector<int> MGId, std::vector<float> MGBtag);

  TRootMCParticle *Top,*TopBar,*Bottom, *BottomBar,*Lepton,*NeutrinoMC,*WPlus,*WMinus,*Light,*LightBar;
  TLorentzVector *GenLeptonicTop, *GenLeptonicW, *GenLepton, *GenNeutrino;
  //TLorentzVector *GenHadronicTop, *GenHadronicW;
  unsigned int NumberNegativeElectrons, NumberNegativeMuons, NumberPositiveElectrons, NumberPositiveMuons, WrongEvtCounter;
  unsigned int NumberNegRecoEl, NumberNegRecoMu, NumberPosRecoEl, NumberPosRecoMu, NrPosRecoMuCorrect, NrPosRecoMuWrong;
  bool CorrectGenEvtContent;
  ofstream GenOutFile[4], RecoOutFile[4], WrongGenFile, CorrectRecoMuPosFile, WrongRecoMuPosFile;
  int verbose_, LeptonCharge;
  bool genOutput_, recoOutput_;

  enum LeptonType_t {muPlus, muMinus, elPlus, elMinus};
  LeptonType_t leptonType;
};

#endif
