#ifndef TFCreation_h
#define TFCreation_h

#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include <map>
#include "TLorentzVector.h"
#include "TFile.h"
#include "TObjArray.h"

#include <fstream>
#include <sstream>

using namespace std;

class TFCreation{

    public:
        TFCreation();
        ~TFCreation();
	void InitializeVariables();
	void FillHistograms(TRootMCParticle* hadrWJet1, TRootMCParticle* hadrWJet2, TRootMCParticle* hadrBJet, TRootMCParticle* leptBJet, TRootMCParticle* lepton, TRootJet* selHadrWJet1, TRootJet* selHadrWJet2, TRootJet* selHadrBJet, TRootJet* selLeptBJet, TLorentzVector* selLepton, bool isSemiMu, bool isSemiEl);
	void CalculateTF(bool, bool, bool, bool, bool);
        void CalculateTFFromFile(TH2F*, bool, int, bool, bool, float[], TFile*);
	void FitSliceClassCode(TH2F*, int, const char*[]);
	void SetStartValuesDoubleGaussian(int, bool);
	void WriteTF(TFile*);
	void WritePlots(TFile*);

    private:
	map<string,TH1F*> histo1D;
	map<string,TH2F*> histo2D;

	template <typename T> string tostr(const T& t) { ostringstream os; os<<t; return os.str(); }

	TF1 *doubleGaussianFit, *caloEnergyFit;
	TH1D **hlist;
        float* startValuesArray;
};
#endif
//
