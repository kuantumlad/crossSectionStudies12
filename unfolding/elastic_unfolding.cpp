#include <TMath.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TFitter.h>
#include <TF1.h>
#include <TStyle.h>
#include <TVector.h>
#include <TGraph.h>

#include "TUnfoldDensity.h"

Double_t E_ele; Double_t px_ele; Double_t py_ele; Double_t pz_ele; 
Double_t E_ele_mc; Double_t px_ele_mc; Double_t py_ele_mc; Double_t pz_ele_mc;

int elastic_unfolding(const char *inSimFile, const char *outFile, int run ){

  TFile *fSim;
  TFile *fData;

  fSim = new TFile(inSimFile,"");
  //fData = new TFile(inDataFile,"");
  
  Char_t *inSimTreeName="events_epX";
  TTree *simTree=(TTree *) fSim->Get(inSimTreeName);
  /*simTree->SetBranchAddress("E_ele",&E_ele);
  simTree->SetBranchAddress("px_ele",&px_ele);
  simTree->SetBranchAddress("py_ele",&py_ele);
  simTree->SetBranchAddress("pz_ele",&pz_ele);

  simTree->SetBranchAddress("E_ele_mc",&E_ele_mc);
  simTree->SetBranchAddress("px_ele_mc",&px_ele_mc);
  simTree->SetBranchAddress("py_ele_mc",&py_ele_mc);
  simTree->SetBranchAddress("pz_ele_mc",&pz_ele_mc);
  */
  /*
  Char_t *inDataTreeName="events_epX";
  TTree *dataTree=(TTree *) fData->Get(inSimTreeName);
  dataTree->SetBranchAddress("E_ele",&E_ele);
  dataTree->SetBranchAddress("px_ele",&px_ele);
  dataTree->SetBranchAddress("py_ele",&py_ele);
  dataTree->SetBranchAddress("pz_ele",&pz_ele);
  */



  int nGen = 20;
  double xminGen = 0.0;
  double xmaxGen = 20.0;
  int nDet = 20;
  double xminDet = 0.0;
  double xmaxDet = 20.0;
  
  TH2D *h_unfolding_matrix = new TH2D("theta_unfolding_matrix",";thetagen;thetarec", nGen, xminGen, xmaxGen, nDet, xminDet, xmaxDet );

  for( int i = 0; i < simTree->GetEntriesFast(); i++ ){
    simTree->GetEntry(i);


    double el_gen_px = px_ele;
    double el_gen_py = py_ele;
    double el_gen_pz = pz_ele;

    double el_rec_px = px_ele_mc;
    double el_rec_py = py_ele_mc;
    double el_rec_pz = pz_ele_mc;

    TLorentzVector lv_el_gen;
    TLorentzVector lv_el_rec;
    
    std::cout << " gen px " << el_gen_px << std::endl;
    
    //h_unfolding_matrix->Fill(

  }
  
  return 0;
}
