#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h> 
#include "Riostream.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TMath.h"
#include <vector>
#include <TAxis.h>
#include <TLorentzRotation.h>
#include<vector>
#include <algorithm>
#include <functional>
#include <TCutG.h>
#include <TStyle.h>


// used to compare the electron hit position in data and simulation for a given run and magnetic field configuration
// TO USE :

//selectorComparisonPlotter("elastic_out_clas12_002587.10.99.root","mc_selector_clas12_2GeV.0.100_tm06sm06.root","sim_elastic_clas12_2GeV.0.100_tm06sm06.root",2587,"tm06sm06")
double Ebeam = 7.546;

int genPlotter( const char* inFileMC, const char* genName ){

  TFile *fMC;
  Char_t tmpstr[80];
  Double_t fraction;

  fMC = new TFile(inFileMC,"");
  if( fMC->IsZombie() ){   // Check if TFile exists!
    cout<<"Input file doesn't exist!" << endl;
    cout<<"Exit program" << endl;
    return 0;
  }

  cout << "Reading from File: " << inFileMC << " "  << endl;
  
  bool printAll = false;
  bool setMinZero = true;
  double pmin = Ebeam-2.0;
  if( !setMinZero ) pmin = 0.0;

  TCanvas *c0 = new TCanvas("c0","c0",900,900);
  c0->cd(1);
  TH1D *h_gen_p = (TH1D*)fMC->Get("mc/hist_mc_all_electron_p");
  TH1D *h_gen_theta = (TH1D*)fMC->Get("mc/hist_mc_all_electron_theta");
  TH1D *h_gen_phi = (TH1D*)fMC->Get("mc/hist_mc_all_electron_phi");
  TH2F *h_gen_thetap = (TH2F*)fMC->Get("mc/hist_mc_all_electron_p_vs_theta");
  TH2F *h_gen_thetaphi = (TH2F*)fMC->Get("mc/hist_mc_all_electron_theta_vs_phi");
  TH2F *h_gen_wtheta = (TH2F*)fMC->Get("mc/h_mc_el_wtheta");

  TCanvas *c1 = new TCanvas("c1","c1",1200,600);

  c1->Divide(2,1);
  c1->cd(1);
  gPad->SetLogz();
  h_gen_thetap->Draw("colz");
  c1->cd(2);
  gPad->SetLogz();
  h_gen_thetaphi->Draw("colz");
  c1->Update();
  c1->Print(Form("h_generator_kin_%s.pdf",genName),"pdf");
  c1->Clear();

  c1->SetWindowSize(900,900);
  c1->Divide(1,1);
  c1->cd(1);
  gPad->SetLogz();
  h_gen_wtheta->SetTitle("Generated W vs Theta");
  h_gen_wtheta->GetXaxis()->CenterTitle();
  h_gen_wtheta->GetYaxis()->CenterTitle();  
  h_gen_wtheta->Draw("colz");
  h_gen_wtheta->GetXaxis()->SetRangeUser(0.8, 1.4);
  c1->Update();
  c1->Print(Form("h_generator_wtheta_%s.pdf",genName),"pdf");
  c1->Clear();



  return 0;
}
