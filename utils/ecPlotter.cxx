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
using namespace std;

int ecPlotter(const char* input, int run){

  TFile *fIn = new TFile(input,"");


  double p_min = 1.5;
  double Ebeam = 2.221;

  if(Ebeam > 10) p_min = 1.5;
  if(Ebeam < 10) p_min = 1.0;
  if(Ebeam < 3)  p_min = 0.5;

  double sigma_range = 3;

  /// //////////////////////////////////////////////////////////////////////////////////////////////////////
  /// a) cut based on sampling fraction versus drift chamber momentum

  // 10.6 GeV (5b61)
  double p0mean[] = {0.106333, 0.113711, 0.107714, 0.113276, 0.115548, 0.11108};
  double p1mean[] = {-0.374003, 0.164037, -0.101566, 0.22524, 0.272903, 0.0370852};
  double p2mean[] = {0.00816235, 0.00390166, 0.00832663, 0.00324039, 0.00376747, 0.00899919};
  double p3mean[] = {-0.000789648, -0.000432195, -0.00091734, -0.00013304, -0.000173935, -0.000932962};
  double p0sigma[] = {0.0162041, 0.0256151, 0.00996036, 0.0174414, 0.0195056, 0.0115662};
  double p1sigma[] = {0.00472716, -0.00465669, 0.0140508, 0.00455405, 0.000429308, 0.0119683};

  // mean = p0mean[k] *( 1 + part_p[j]/sqrt(pow(part_p[j],2) + p1mean[k])) + p2mean[k] * part_p[j] + p3mean[k] * pow(part_p[j],2);
  //sigma = p0sigma[k] + p1sigma[k] / sqrt(part_p[j]);
  //upper_lim_total = mean + sigma_range * sigma;
  //lower_lim_total = mean - sigma_range * sigma;
  
  TCanvas *c_el_sf = new TCanvas("c_el_sf", "c_el_sf", 1800, 900);
  c_el_sf->Divide(3,2);

  for( int s = 1; s <= 6; s++ ){
    c_el_sf->cd(s);
    
    TCanvas *c_temp_sf = new TCanvas(Form("c_temp_sf_s%d",s),Form("c_temp_sf_s%d",s),900,600);


    TH2F *h_el_sf = (TH2F*)fIn->Get(Form("FD_PID_electron_EC_plots/EC_total_sampling_fraction_sec%d_cut_00",s));

    TF1 *fit_top  = new TF1("fit_top","[a]*( 1 + x/sqrt(x*x + [b])) + [c]*x + [d]*x*x + 3*([e] + [f]/x)", 0.50, Ebeam);
    fit_top->SetParameter(0, p0mean[s-1] );
    fit_top->SetParameter(1, p1mean[s-1] );
    fit_top->SetParameter(2, p2mean[s-1] );
    fit_top->SetParameter(3, p3mean[s-1] );
    fit_top->SetParameter(4, p0sigma[s-1] );
    fit_top->SetParameter(5, p1sigma[s-1] );

    TF1 *fit_bot  = new TF1("fit_bot","[a]*( 1 + x/sqrt(x*x + [b])) + [c]*x + [d]*x*x - 3*([e] + [f]/x)", 0.50, Ebeam);
    fit_bot->SetParameter(0, p0mean[s-1] );
    fit_bot->SetParameter(1, p1mean[s-1] );
    fit_bot->SetParameter(2, p2mean[s-1] );
    fit_bot->SetParameter(3, p3mean[s-1] );
    fit_bot->SetParameter(4, p0sigma[s-1] );
    fit_bot->SetParameter(5, p1sigma[s-1] );

    fit_top->SetLineStyle(0);
    fit_top->SetLineWidth(6);
    fit_top->SetLineColor(2);

    fit_bot->SetLineStyle(0);
    fit_bot->SetLineWidth(6);
    fit_bot->SetLineColor(2);

    h_el_sf->SetTitle(Form("Sampling Fraction vs p Sector %d",s));
    h_el_sf->GetXaxis()->SetTitle("p [GeV]");
    h_el_sf->GetXaxis()->CenterTitle();
    h_el_sf->GetYaxis()->SetTitle("SF");
    h_el_sf->GetYaxis()->CenterTitle();
    gStyle->SetOptStat(0);
    h_el_sf->Draw("colz");
    fit_top->Draw("same");
    fit_bot->Draw("same");

    c_temp_sf->SaveAs(Form("h_el_sf_sector_%d_run%d.pdf",s,run));

  }

  //c_el_sf->SaveAs(Form("h_el_sf_sector_%d.pdf",run));

  TCanvas *c_ec_hit = new TCanvas("c_ec_hit","c_ec_hit",900,900);
  gStyle->SetOptStat(0);
  gPad->SetLogz();
  TH2F *h_pcal_fid_before = (TH2F*)fIn->Get("FD_PID_electron_EC_plots/EC_PCAL_hit_position_cut_00");
  TH2F *h_pcal_fid_after = (TH2F*)fIn->Get("FD_PID_electron_EC_plots/EC_PCAL_hit_position_cut_04");
  h_pcal_fid_after->SetTitle("PCAL Hit Position");
  h_pcal_fid_after->GetXaxis()->SetTitle("PCAL X Hit Position [cm]");
  h_pcal_fid_after->GetYaxis()->SetTitle("PCAL Y Hit Position [cm]");
  h_pcal_fid_after->GetXaxis()->CenterTitle();
  h_pcal_fid_after->GetYaxis()->CenterTitle();

  
  h_pcal_fid_before->Draw();
  h_pcal_fid_after->Draw("colz");

  c_ec_hit->SaveAs(Form("h_el_echit_%d.pdf",run));


  return 0;
}

