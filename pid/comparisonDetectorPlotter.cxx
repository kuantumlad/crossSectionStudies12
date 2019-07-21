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

int comparisonDetectorPlotter( const char* inFileData, const char* inFileMC, int run , const char* field_config ){

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
  
  TH2D *h_data_ec_pcal_hit_position_final = (TH2D*)fData->Get("FD_PID_electron_EC_plots/EC_PCAL_hit_position_cut_09"); //09
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
  gPad->SetLogz();
  h_data_ec_pcal_hit_position_final->Draw("colz");
  c1->cd(2);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.12);
  gPad->SetLogz();
  h_sim_ec_pcal_hit_position_final->SetTitle("SIM: PCAL Hit Position");
  h_sim_ec_pcal_hit_position_final->GetXaxis()->SetTitle("PCAL hit position X [cm]");
  h_sim_ec_pcal_hit_position_final->GetYaxis()->SetTitle("PCAL hit position Y [cm]");
  h_sim_ec_pcal_hit_position_final->GetXaxis()->CenterTitle();
  h_sim_ec_pcal_hit_position_final->GetYaxis()->CenterTitle();
  h_sim_ec_pcal_hit_position_final->Draw("colz");

  c1->SaveAs(Form("h2_r%d_f%s_comparison_pid_electron_ec_pcal_hit_position_cut09.pdf",run,field_config));


  TCanvas *c_sf_eb = new TCanvas("c_sf_eb","c_sf_eb",900,900);
  c_sf_eb->Divide(2,3);
  gStyle->SetOptStat(0);
  for( int ss = 0; ss < 6; ss++ ){
    c_sf_eb->cd(ss+1);
    TH2F *h_sf_neg_sim = (TH2F*)fMC->Get(Form("FD_PID_electron_EC_plots/EC_total_sampling_fraction_sec%d_cut_09",ss+1));
    TH2F *h_sf_neg_data = (TH2F*)fData->Get(Form("FD_PID_electron_EC_plots/EC_total_sampling_fraction_sec%d_cut_09",ss+1));

    TH2F *h_sim_clone = (TH2F*)h_sf_neg_sim->Clone();
    TH2F *h_data_clone = (TH2F*)h_sf_neg_data->Clone();
    
    h_data_clone->SetTitle(Form("Sampling Fraction Sector %d",ss+1));
    h_data_clone->GetXaxis()->SetTitle("p (GeV)");
    h_data_clone->GetXaxis()->CenterTitle();
    h_data_clone->GetYaxis()->SetTitle("SF");
    h_data_clone->GetYaxis()->CenterTitle();
    
    h_data_clone->Draw("colz");
    h_data_clone->GetXaxis()->SetRange(0, 2000);

    h_sim_clone->Draw("same");
    h_sim_clone->GetXaxis()->SetRange(0, 300);

  }
  c_sf_eb->SaveAs(Form("h_r%d_f%s_comparison_pid_electron_sf_cut09.pdf",run,field_config)); //looks at all negative tracks

  TCanvas *c_htcc_eb = new TCanvas("c_htcc_eb","c_htcc_eb",900,900);
  //c_htcc_eb->Divide(2,3);
  gStyle->SetOptStat(0);
  std::vector<double> poisson_k_sim;
  std::vector<double> poisson_k_data;

  //for( int ss  = 0; ss < 6 ; ss++ ){
  c_htcc_eb->cd(1);//ss+1);
  TH1F *h_nphe_neg_data = (TH1F*)fData->Get("FD_PID_electron_HTCC_plots/HTCC_nphe_cut_00");
  TH1F *h_nphe_neg_sim = (TH1F*)fMC->Get("FD_PID_electron_HTCC_plots/HTCC_nphe_cut_00");
  h_nphe_neg_data->GetXaxis()->SetTitle("HTCC nphe");
  h_nphe_neg_data->GetXaxis()->CenterTitle();
  double scaler_htcc = h_nphe_neg_sim->GetMaximum()/h_nphe_neg_data->GetMaximum();
  h_nphe_neg_data->Scale(scaler_htcc);
  h_nphe_neg_data->SetLineColor(kRed);
  h_nphe_neg_data->Draw("hist");
  h_nphe_neg_sim->Draw("hist+same");
  TF1 *f1_data = new TF1("f1_data","[0]*TMath::Power(([1]/[2]),(x/[2]))*(TMath::Exp(-([1]/[2])))/TMath::Gamma((x/[2])+1.)", 4, 15); 
  TF1 *f1_sim = new TF1("f1_sim","[0]*TMath::Power(([1]/[2]),(x/[2]))*(TMath::Exp(-([1]/[2])))/TMath::Gamma((x/[2])+1.)", 20, 39); 
  // "xmin" = 0, "xmax" = 10
  f1_data->SetParameters(h_nphe_neg_data->GetMaximum(), 30, 1); // you MUST set non-zero initial values for parameters
  f1_sim->SetParameters(h_nphe_neg_sim->GetMaximum(), 10, 1); // you MUST set non-zero initial values for parameters
  std::cout << " fit data " << std::endl;
  h_nphe_neg_data->Fit("f1_data", "R"); // "R" = fit between "xmin" and "xmax" of the "f1"
  std::cout << " fit sim " << std::endl;
  h_nphe_neg_sim->Fit("f1_sim", "R"); // "R" = fit between "xmin" and "xmax" of the "f1"
  f1_data->SetLineColor(kRed);
  f1_sim->SetLineColor(kBlue);
  f1_data->Draw("same");
  f1_sim->Draw("same");
    
  poisson_k_sim.push_back(f1_sim->GetParameter(1)/f1_sim->GetParameter(2));
  poisson_k_data.push_back(f1_data->GetParameter(1)/f1_data->GetParameter(2));
      
  c_htcc_eb->SaveAs(Form("h_r%d_f%s_comparison_pid_electron_nphe_cut09.pdf",run,field_config));

  //print out fit parameters 
   std::cout << " simulation fit parameter k " << std::endl;
    std::cout << poisson_k_sim[0] << std::endl;
  std::cout << " data fit parameter k " << std::endl;
    std::cout << poisson_k_data[0] << std::endl;
  
  return 0;
}

