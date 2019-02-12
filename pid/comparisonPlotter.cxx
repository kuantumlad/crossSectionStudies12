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

int comparisonPlotter( const char* inFileData, const char* inFileMC, int run , const char* field_config ){

  TFile *fData; 
  TFile *fMC;
  Char_t tmpstr[80];
  Double_t fraction;

  fData = new TFile(inFileData,"");   // Input File
  fMC = new TFile(inFileMC,"");
  if(fData->IsZombie() || fMC->IsZombie()){   // Check if TFile exists!
    cout<<"Input file doesn't exist!" << endl;
    cout<<"Exit program" << endl;
    return 0;
  }

  cout << "Reading from File: " << inFileData << " and " << inFileMC << endl;
  
  TH2D *h_data_ec_pcal_hit_position_final = (TH2D*)fData->Get("FD_PID_electron_EC_plots/EC_PCAL_hit_position_cut_09");
  TH2D *h_sim_ec_pcal_hit_position_final = (TH2D*)fMC->Get("FD_PID_electron_EC_plots/EC_PCAL_hit_position_cut_09");

  TCanvas *c1 = new TCanvas("C1","c1",920, 470);
  gStyle->SetOptStat(0000);
  c1->Divide(2,1);
  c1->cd(1);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.12);
  h_data_ec_pcal_hit_position_final->SetTitle("Data: PCAL Hit Position");
  h_data_ec_pcal_hit_position_final->GetXaxis()->SetTitle("PCAL hit position X [cm]");
  h_data_ec_pcal_hit_position_final->GetYaxis()->SetTitle("PCAL hit position Y [cm]");
  h_data_ec_pcal_hit_position_final->GetXaxis()->CenterTitle();
  h_data_ec_pcal_hit_position_final->GetYaxis()->CenterTitle();
  h_data_ec_pcal_hit_position_final->Draw("colz");
  c1->cd(2);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.12);
  h_sim_ec_pcal_hit_position_final->SetTitle("SIM: PCAL Hit Position");
  h_sim_ec_pcal_hit_position_final->GetXaxis()->SetTitle("PCAL hit position X [cm]");
  h_sim_ec_pcal_hit_position_final->GetYaxis()->SetTitle("PCAL hit position Y [cm]");
  h_sim_ec_pcal_hit_position_final->GetXaxis()->CenterTitle();
  h_sim_ec_pcal_hit_position_final->GetYaxis()->CenterTitle();
  h_sim_ec_pcal_hit_position_final->Draw("colz");

  c1->SaveAs(Form("h2_r%d_f%s_comparison_pid_electron_ec_pcal_hit_position_cut09.pdf",run,field_config));



  return 0;
}

