#include <vector>
#include <string>
#include <map>
#include <iostream>
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"

const static int BUFFER = 50;
double Ebeam = 2.211;//2193;

int unfoldingExtractor(const char* inFileData, int run, const char* field_config ){ 

 
  TFile *fData; 
  TFile *fMC; 
  Char_t tmpstr[80];
  Double_t fraction;

  fData = new TFile(inFileData,"");   // Input File
  if(fData->IsZombie() ){
    cout<<"Input file doesn't exist!" << endl;
    cout<<"Exit program" << endl;
    return 0;
  }

  cout << "Reading from File: " << inFileData << endl;

  return 0;

}
