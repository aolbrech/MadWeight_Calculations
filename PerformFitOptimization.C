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

TFile *inFile = new TFile("Events/MTop_MGSampleCreatedWith174_SingleGausTF_10000Evts_Narrow/FitDistributions_MGSample_MTop_1000Evts.root","READ"); 
TFile *outputFile = new TFile("Events/MTop_MGSampleCreatedWith174_SingleGausTF_10000Evts_Narrow/FitOptimizations_MGSample_MTop_1000Evts.root","RECREATE"); 

int NrEvts = 10; 
int xBin = 5; 
float xLow = 170.5; 
float xHigh = 175.5; 

void PaintOverflow(TH1F *h, TFile *FileToWrite, std::string dirName){     // This function draws the histogram h with an extra bin for overflows
  Int_t nx    = h->GetNbinsX()+1;
  Double_t x1 = h->GetBinLowEdge(1), bw = h->GetBinWidth(nx), x2 = h->GetBinLowEdge(nx)+bw;

  //Define a temporary histogram having an extra bin for overflows
  char newTitle[100], newName[100];
  strcpy(newTitle,h->GetTitle()); strcat(newTitle," (under- and overflow added)" );
  strcpy(newName,h->GetName());  strcat(newName, "_Flow");
  TH1F *h_tmp = new TH1F(newName, newTitle, nx, x1, x2);

  // Fill the new histogram including the extra bin for overflows
  for (Int_t i=1; i<=nx; i++)
    h_tmp->Fill(h_tmp->GetBinCenter(i), h->GetBinContent(i));
  // Fill the underflows
  h_tmp->Fill(x1-1, h->GetBinContent(0));

  // Restore the number of entries
  h_tmp->SetEntries(h->GetEntries());

  //Set the correct path to save the file
  const char* path = 0;
  if(dirName != ""){
    TDirectory *dir = FileToWrite->GetDirectory(dirName.c_str());
    if (!dir)
      dir = FileToWrite->mkdir(dirName.c_str());
    path = dir->GetName();
  }

  FileToWrite->cd(path);
  h_tmp->Write();
  FileToWrite->cd();  //Reset to general directory! 
}                  

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

    TH1F *h_LLSum = 0;
    TH1F *h_FitSum = 0, *h_ChiSq = 0, *h_FitSumSmallChiSq = 0, *h_SlopeFit = 0, *h_FitSumSmallChiSqPosSlope = 0;
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
      h_SlopeFit = new TH1F((FitNr+""+NormType+"_SlopeVar").c_str(), ("Distribution of slope of "+FitNr).c_str(), 200, -5.0, 5.0);
      h_ChiSq = new TH1F( ("ChiSq"+NormType+"_"+FitNr).c_str(), ("ChiSquared distribution of "+FitNr+" (norm = "+NormType+")").c_str(),200, 0, 0.005);
      h_FitSumSmallChiSq = new TH1F((FitNr+""+NormType+"_SmallChiSq").c_str(), ("Distribution of "+FitNr+"after summing over "+NrUsedPoints+" points (ChiSq-cut applied)").c_str(),xBin,xLow,xHigh);
      h_FitSumSmallChiSqPosSlope = new TH1F((FitNr+""+NormType+"_SmallChiSqPosSlope").c_str(), ("Distribution of "+FitNr+"after summing over "+NrUsedPoints+" points (ChiSq-cut applied, positive slope required)").c_str(),xBin,xLow,xHigh);
    }
   
    //Initialize counters:
    int allEvts = 0, evtsSmallChiSq = 0, evtsSmallChiSqPosSlope = 0;
 
    std::string fitName = "";
    while ((object_dir = next())) {         //Keep it general to TObject because both TF1 and TH1F are stored in this ROOT file!
      if(lookAtFits == false){  //Looking at TH1F objects!
        TH1F *h_iLL = (TH1F*) dir->Get(object_dir->GetName());
        h_LLSum->Add(h_iLL);
      }
      else if(lookAtFits == true){ //Looking at TF1 objects!
        TF1 *fitFunc = (TF1*) dir->Get(object_dir->GetName());
        h_FitSum->Add(fitFunc);
        h_SlopeFit->Fill(fitFunc->GetParameter(2));
        allEvts++;
        h_ChiSq->Fill(fitFunc->GetChisquare());
        if(fitFunc->GetChisquare() < 0.0002){
          h_FitSumSmallChiSq->Add(fitFunc);
          evtsSmallChiSq++;
          if(fitFunc->GetParameter(2) > 0){
            h_FitSumSmallChiSqPosSlope->Add(fitFunc);
            evtsSmallChiSqPosSlope++;
          }
        }
      }
    }

    if(lookAtFits == true){
      h_FitSum->Write();
      PaintOverflow(h_SlopeFit,outputFile,"");
      PaintOverflow(h_ChiSq,outputFile,"ChiSq");
      h_FitSumSmallChiSq->Write();
      std::cout << " Histogram " << h_FitSumSmallChiSq->GetName() << " has " << evtsSmallChiSq << "/" << allEvts << " events selected !" << std::endl;
      h_FitSumSmallChiSqPosSlope->Write();
      std::cout << " Histogram " << h_FitSumSmallChiSqPosSlope->GetName() << " has " << evtsSmallChiSqPosSlope << "/" << allEvts << " events selected !" << std::endl;
    }
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
