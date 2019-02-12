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
#include "loadECParameters.C"
using namespace std;

int ecPlotter(const char* input, int run){

  TFile *fIn = new TFile(input,"");

  double p_min = 1.5;
  double Ebeam = 10.5;

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


  std::vector<double> p0_mean;
  std::vector<double> p1_mean;
  std::vector<double> p2_mean;
  std::vector<double> p3_mean;

  std::vector<double> p0_sigma;
  std::vector<double> p1_sigma;

  std::map<int, std::vector<double> > ec_mean;
  std::map<int, std::vector<double> > ec_sig;
  //sf_cut_mean_run2
  std::string parm_path = "/w/hallb-scifs17exp/clas12/bclary/CLAS12/electron_studies/pid/parameters/sf_cut_";
  std::cout << " loading ec sampling fraction parameters from " << parm_path << " for run " << run << std::endl;

  ec_mean = loadECParameters(parm_path + "mean_run" + std::to_string(run) + ".txt");
  ec_sig = loadECParameters(parm_path + "sigma_run" + std::to_string(run) + ".txt");
  
  for( int i = 0; i < 6; i++ ){   
    p0_mean.push_back( ec_mean[0][i] );
    p1_mean.push_back( ec_mean[1][i] );
    p2_mean.push_back( ec_mean[2][i] );
    p3_mean.push_back( ec_mean[3][i] );

    p0_sigma.push_back( ec_sig[0][i] );
    p1_sigma.push_back( ec_sig[1][i] );
  }


  // mean = p0mean[k] *( 1 + part_p[j]/sqrt(pow(part_p[j],2) + p1mean[k])) + p2mean[k] * part_p[j] + p3mean[k] * pow(part_p[j],2);
  //sigma = p0sigma[k] + p1sigma[k] / sqrt(part_p[j]);
  //upper_lim_total = mean + sigma_range * sigma;
  //lower_lim_total = mean - sigma_range * sigma;
  
  TCanvas *c_el_sf = new TCanvas("c_el_sf", "c_el_sf", 1800, 900);
  c_el_sf->Divide(3,2);

  for( int s = 1; s <= 6; s++ ){
    c_el_sf->cd(s);
    
    TCanvas *c_temp_sf = new TCanvas(Form("c_temp_sf_s%d",s),Form("c_temp_sf_s%d",s),900,600);


    TH2F *h_el_sf = (TH2F*)fIn->Get(Form("FD_PID_electron_EC_plots/EC_total_sampling_fraction_sec%d_cut_04",s));

    TF1 *fit_mid  = new TF1("fit_mid","[a]*( 1 + x/sqrt(x*x + [b])) + [c]*x + [d]*x*x + 3.0*([e] + [f]/x)", 0.8, 2.20);
    fit_mid->SetParameter(0, p0_mean[s-1] );
    fit_mid->SetParameter(1, p1_mean[s-1] );
    fit_mid->SetParameter(2, p2_mean[s-1] );
    fit_mid->SetParameter(3, p3_mean[s-1] );
    

    TF1 *fit_top  = new TF1("fit_top","[a]*( 1 + x/sqrt(x*x + [b])) + [c]*x + [d]*x*x + 3.0*([e] + [f]/x)", 0.8, 2.20);
    fit_top->SetParameter(0, p0_mean[s-1] );
    fit_top->SetParameter(1, p1_mean[s-1] );
    fit_top->SetParameter(2, p2_mean[s-1] );
    fit_top->SetParameter(3, p3_mean[s-1] );
    fit_top->SetParameter(4, p0_sigma[s-1] );
    fit_top->SetParameter(5, p1_sigma[s-1] );

    TF1 *fit_bot  = new TF1("fit_bot","[a]*( 1 + x/sqrt(x*x + [b])) + [c]*x + [d]*x*x - 3.0*([e] + [f]/x)", 0.8, 2.2);
    fit_bot->SetParameter(0, p0_mean[s-1] );
    fit_bot->SetParameter(1, p1_mean[s-1] );
    fit_bot->SetParameter(2, p2_mean[s-1] );
    fit_bot->SetParameter(3, p3_mean[s-1] );
    fit_bot->SetParameter(4, p0_sigma[s-1] );
    fit_bot->SetParameter(5, p1_sigma[s-1] );

    fit_mid->SetLineStyle(2);
    fit_mid->SetLineWidth(2);
    fit_mid->SetLineColor(2);

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
    fit_mid->Draw("same");
    fit_top->Draw("same");
    fit_bot->Draw("same");

    c_temp_sf->SaveAs(Form("h_el_sf_sector_%d_run%d.pdf",s,run));

  }

  c_el_sf->SaveAs(Form("h_el_sf_sector_%d.pdf",run));

  TCanvas *c_ec_hit = new TCanvas("c_ec_hit","c_ec_hit",900,900);
  gStyle->SetOptStat(0);
  gPad->SetLogz();
  c_ec_hit->SetLeftMargin(0.125);
  TH2F *h_pcal_fid_before = (TH2F*)fIn->Get("FD_PID_electron_EC_plots/EC_PCAL_hit_position_cut_00");
  TH2F *h_pcal_fid_after = (TH2F*)fIn->Get("FD_PID_electron_EC_plots/EC_PCAL_hit_position_cut_04");
  h_pcal_fid_before->Draw("scat");
  h_pcal_fid_after->Draw("same+colz");

  c_ec_hit->Update();
  h_pcal_fid_before->SetTitle("PCAL Hit Position");
  h_pcal_fid_before->GetXaxis()->SetTitle("PCAL X Hit Position [cm]");
  h_pcal_fid_before->GetYaxis()->SetTitle("PCAL Y Hit Position [cm]");
  h_pcal_fid_before->GetXaxis()->CenterTitle();
  h_pcal_fid_before->GetYaxis()->CenterTitle();
  c_ec_hit->Update();

  c_ec_hit->SaveAs(Form("h_el_echit_%d.pdf",run));


  return 0;
}

