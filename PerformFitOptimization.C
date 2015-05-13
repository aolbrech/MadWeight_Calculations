#include <iostream>
#include <fstream>
#include <sstream>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TDirectory.h"
#include <string>
#include <map>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <cstdlib>

TFile *inFile = new TFile("FitDistributions_10000RecoEvts_ManySteps.root","READ");       //Should be filled by the python script!
//TDirectory *dir_FstFit = inFile->GetDirectory("FirstPolynomialFit"), *dir_FstFitXS = inFile->GetDirectory("FirstPolynomialFit_XS"), *dir_FstFitAcc = inFile->GetDirectory("FirstPolynomialFit_Acc");
//TDirectory *dir_ScdFit = inFile->GetDirectory("SecondPolynomialFit"),*dir_ScdFitXS = inFile->GetDirectory("SecondPolynomialFit_XS"),*dir_ScdFitAcc = inFile->GetDirectory("SecondPolynomialFit_Acc");
//TDirectory *dir_OrigLL = inFile->GetDirectory("OriginalLL"),         *dir_OrigLLXS = inFile->GetDirectory("OriginalLL_XS"),         *dir_OrigLLAcc = inFile->GetDirectory("OriginalLL_Acc");

TFile *outputFile = new TFile("FitOptimizations.root","RECREATE"); //Should be filled by the python script!

int NrEvts = 10; 
int xBin = 21;        //Copy this from the python file!
float xLow = -0.105;  //Copy this from the python file!
float xHigh = 0.105;  //Copy this from the python file!

void getIndividualDirObjects(TDirectory *dir){
  outputFile->cd();

  std::stringstream ssNkeys; ssNkeys << dir->GetNkeys(); std::string sNkeys = ssNkeys.str();
  //if(!(dir->GetNkeys() != NrEvts) ){
  if( dir->GetNkeys() > 0 ){
    
    TList *list_dir = dir->GetListOfKeys();
    TIter next(list_dir); 
    TObject* object_dir = 0;  

    //Convert name of directory (const char*) to a string:
    std::string dirName(dir->GetName(), 0, 100);

    TH1F *h_LLSum = 0, *h_FitSum = 0;
    //Identify which directory is being studied --> Containing TF1 or TH1F & type of normalisation!
    std::string NormType = "", FitNr = "", NrUsedPoints = "";
    bool lookAtFits = false;
    if(dirName.find("Acc") != std::string::npos)     NormType = "Acc";
    else if(dirName.find("XS") != std::string::npos) NormType = "XS";
    else                                             NormType = "";

    if(dirName.find("FirstPolynomial") != std::string::npos){       FitNr = "FirstPol";  NrUsedPoints = "all";    }
    else if(dirName.find("SecondPolynomial") != std::string::npos){ FitNr = "SecondPol"; NrUsedPoints = "reduced";}

    if (dirName.find("LL") != std::string::npos)
      h_LLSum = new TH1F(("LL"+NormType+"_Summed").c_str(), ("Distribution of "+dirName+" after summing over all "+sNkeys.c_str()+" events").c_str(), xBin, xLow, xHigh);
    else{
      lookAtFits = true;
      h_FitSum = new TH1F((FitNr+""+NormType+"_Summed").c_str(), ("Distribution of "+FitNr+" after summing over "+NrUsedPoints+" points").c_str(), xBin, xLow, xHigh);
    }
    
    std::string fitName = "";
    double fitXMin = 0, fitXMax = 0;
    while ((object_dir = next())) {         //Keep it general to TObject because both TF1 and TH1F are stored in this ROOT file!
      if(lookAtFits == false){  //Looking at TH1F objects!
        TH1F *h_iLL = (TH1F*) dir->Get(object_dir->GetName());
        h_LLSum->Add(h_iLL);
      }
      else if(lookAtFits == true){ //Looking at TF1 objects!
        TF1 *fitFunc = (TF1*) dir->Get(object_dir->GetName());
        h_FitSum->Add(fitFunc);
      }
    }

    if(lookAtFits == true)
      h_FitSum->Write();
    else 
      h_LLSum->Write();
  }
  else{
    std::cout << "        ### Skipped the empty directory " << dir->GetName() << " ### " << std::endl;
  }
}

void PerformFitOptimization(){

  std::cout << " Looking at TFile " << inFile->GetName() << " (which has " << inFile->GetNkeys() << " directories) " << std::endl;
  TList *list_file = inFile->GetListOfKeys();
  TIter next(list_file);
  TObject* object_file = 0;
  TDirectory* dir_file = 0;
  //while( (object_file = next() ) ){
  while( (object_file = next() ) ){
    std::cout << "  ** Name of directory : " << object_file->GetName() << std::endl;
    dir_file = inFile->GetDirectory(object_file->GetName());
    getIndividualDirObjects(dir_file);
  }

  inFile->Close();
  outputFile->Close();
}
