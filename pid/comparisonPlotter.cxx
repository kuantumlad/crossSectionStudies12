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

  TCanvas *c2 = new TCanvas("c2","c2",600,900);
  c2->Divide(2,3);
  for( int s = 0; s < 6; s++ ){
    TH1D *h_data_w = (TH1D*)fData->Get(Form("kinematics/W_sec_0%d",s+1));
    TH1D *h_sim_w = (TH1D*)fMC->Get(Form("kinematics/W_sec_0%d",s+1));
    c2->cd(s+1);
    h_sim_w->SetTitle("Data vs Sim: W (GeV)");
    h_sim_w->GetXaxis()->SetTitle("W (GeV)");
    h_sim_w->GetXaxis()->CenterTitle();
    //h_sim_w->Scale( 1.0/h_data_w->Integral() );
    //h_data_w->Scale( 1.0/h_sim_w->Integral() );
    h_sim_w->Draw();
    h_data_w->SetLineColor(kRed);
    h_data_w->Draw("same");    
  }
  c2->SaveAs(Form("h_r%d_f%s_comparison_pid_electron_W_per_sector.pdf",run,field_config));


  TH1D *h_data_p = (TH1D*)fData->Get("particles_identified_histograms_selected/hist_electron_p");
  TH1D *h_sim_p = (TH1D*)fMC->Get("particles_identified_histograms_selected/hist_electron_p");
  TH1D *h_data_theta = (TH1D*)fData->Get("particles_identified_histograms_selected/hist_electron_theta");
  TH1D *h_sim_theta = (TH1D*)fMC->Get("particles_identified_histograms_selected/hist_electron_theta");
  TH1D *h_data_phi = (TH1D*)fData->Get("particles_identified_histograms_selected/hist_electron_phi");
  TH1D *h_sim_phi = (TH1D*)fMC->Get("particles_identified_histograms_selected/hist_electron_phi");
  
  TCanvas *c3 = new TCanvas("c3","c3",900,900);
  c3->cd(1);
  h_sim_p->SetTitle("Data vs Sim: P (GeV)");
  h_sim_p->GetXaxis()->SetTitle("p (GeV)");
  h_sim_p->GetXaxis()->CenterTitle();
  double scaler_p = h_sim_p->GetMaximum()/h_data_p->GetMaximum();
  h_data_p->Scale( scaler_p );
  h_sim_p->Draw("HIST+same");
  h_data_p->SetLineColor(kRed);
  h_data_p->Draw("HIST+same");
  c3->SaveAs(Form("h_r%d_f%s_comparison_pid_electron_p.pdf",run,field_config));      

  TCanvas *c4 = new TCanvas("c4","c4",900,900);
  c4->cd(1);
  h_sim_theta->SetTitle("Data vs Sim: #theta (deg)");
  h_sim_theta->GetXaxis()->SetTitle("#theta (deg)");
  h_sim_theta->GetXaxis()->CenterTitle();
  double scaler_theta = h_sim_theta->GetMaximum()/h_data_theta->GetMaximum();
  h_data_theta->Scale( scaler_theta );
  h_sim_theta->Draw();
  h_data_theta->SetLineColor(kRed);
  h_data_theta->Draw("HIST+same");
  c4->SaveAs(Form("h_r%d_f%s_comparison_pid_electron_theta.pdf",run,field_config));      

  TCanvas *c5 = new TCanvas("c5","c5",900,900);
  c5->cd(1);
  h_sim_phi->SetTitle("Data vs Sim: #phia (deg)");
  h_sim_phi->GetXaxis()->SetTitle("#phi (deg)");
  h_sim_phi->GetXaxis()->CenterTitle();
  double scaler_phi = h_sim_phi->GetMaximum()/h_data_phi->GetMaximum();
  h_data_phi->Scale( scaler_phi );
  h_sim_phi->Draw();
  h_data_phi->SetLineColor(kRed);
  h_data_phi->Draw("HIST+same");
  c5->SaveAs(Form("h_r%d_f%s_comparison_pid_electron_phi.pdf",run,field_config));      
  

  TCanvas *c6 = new TCanvas("c6","c6",600,900);
  c6->Divide(2,3);
  for( int s = 0; s < 6; s++ ){
    c6->cd(s+1);
    TH1D* h_data_vz = (TH1D*)fData->Get(Form("FD_PID_electron_DC_plots/DC_z_vertex_sec%d_cut_09",s+1));
    TH1D* h_sim_vz = (TH1D*)fMC->Get(Form("FD_PID_electron_DC_plots/DC_z_vertex_sec%d_cut_09",s+1));
    
    h_sim_vz->SetTitle("Data vs Sim: Vz (cm)");
    h_sim_vz->GetXaxis()->SetTitle("Vz (cm)");
    h_sim_vz->GetXaxis()->CenterTitle();
    double scaler_vz = h_sim_vz->GetMaximum()/h_data_vz->GetMaximum();
    h_data_vz->Scale( scaler_vz );
    
    h_sim_vz->Draw();
    h_data_vz->SetLineColor(kRed);
    h_data_vz->Draw("HIST+same");       
  }
  c6->SaveAs(Form("h_r%d_f%s_comparison_pid_electron_vz_per_sector.pdf",run,field_config));    

  return 0;
}

