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

int fiducialPlotter(const char* input, int run){

  TFile *fIn = new TFile(input,"");

  int windowX = 1350;
  bool plotHadrons = false;

  // Histograms with all negative tracks in DC R1 -> R3
  TH2F *h_el_dc_hit_position_r1 = (TH2F*)fIn->Get("FD_PID_electron_DC_plots/DC_hit_position_region1_cut_00"); 
  TH2F *h_el_dc_hit_position_r2 = (TH2F*)fIn->Get("FD_PID_electron_DC_plots/DC_hit_position_region2_cut_00");
  TH2F *h_el_dc_hit_position_r3 = (TH2F*)fIn->Get("FD_PID_electron_DC_plots/DC_hit_position_region3_cut_00");
  // Histograms with cuts on DC R1 and R3 and negative charge and EB PID
  TH2F *h_el_dc_hit_position_r1_cut = (TH2F*)fIn->Get("FD_PID_electron_DC_plots/DC_hit_position_region1_cut_05"); 
  TH2F *h_el_dc_hit_position_r2_cut = (TH2F*)fIn->Get("FD_PID_electron_DC_plots/DC_hit_position_region2_cut_05"); 
  TH2F *h_el_dc_hit_position_r3_cut = (TH2F*)fIn->Get("FD_PID_electron_DC_plots/DC_hit_position_region3_cut_07"); 

  gStyle->SetOptStat(0);

  TCanvas *cel = new TCanvas("cel","cel",800,800);
  h_el_dc_hit_position_r1->SetTitle("Electron DCR1 Fiducial Cut");
  h_el_dc_hit_position_r1->Draw();
  h_el_dc_hit_position_r1_cut->Draw("same+colz");
  h_el_dc_hit_position_r1->GetXaxis()->SetTitle("x (cm)");
  h_el_dc_hit_position_r1->GetXaxis()->CenterTitle();
  h_el_dc_hit_position_r1->GetYaxis()->SetTitle("y (cm)");
  h_el_dc_hit_position_r1->GetYaxis()->CenterTitle();
  cel->SaveAs(Form("h_el_dc_hit_position_r1_%d.pdf",run));

  // Histogram for hadron with 2212 from EB
  TH2F *h_pr_dc_pos_r1 = (TH2F*)fIn->Get("FD_PID_hadron_DC_fiducial_plot/DC_hit_position_region1_hadron_cut_00");
  TH2F *h_pr_dc_pos_r3 = (TH2F*)fIn->Get("FD_PID_hadron_DC_fiducial_plot/DC_hit_position_region3_hadron_cut_00");

  // Histogram for hadron with 2212 from EB, pos charge, and DCR1 Fiducial and DCR3 respectively
  TH2F *h_pr_dc_pos_r1_cut = (TH2F*)fIn->Get("FD_PID_hadron_DC_fiducial_plot/DC_hit_position_region1_hadron_cut_02");
  TH2F *h_pr_dc_pos_r3_cut = (TH2F*)fIn->Get("FD_PID_hadron_DC_fiducial_plot/DC_hit_position_region1_hadron_cut_03");

  // Histogram for hadron with 321 from EB
  TH2F *h_kp_dc_pos_r1 = (TH2F*)fIn->Get("FD_PID_hadron_DC_fiducial_plot/DC_hit_position_region1_hadron_cut_40");
  TH2F *h_kp_dc_pos_r3 = (TH2F*)fIn->Get("FD_PID_hadron_DC_fiducial_plot/DC_hit_position_region3_hadron_cut_40");

  // Histogram for hadron with 321 from EB, pos charge, and DCR1 Fiducial and DCR3 respectively
  TH2F *h_kp_dc_pos_r1_cut = (TH2F*)fIn->Get("FD_PID_hadron_DC_fiducial_plot/DC_hit_position_region1_hadron_cut_42");
  TH2F *h_kp_dc_pos_r3_cut = (TH2F*)fIn->Get("FD_PID_hadron_DC_fiducial_plot/DC_hit_position_region1_hadron_cut_43");

  // Histogram for hadron with -321 from EB
  TH2F *h_km_dc_pos_r1 = (TH2F*)fIn->Get("FD_PID_hadron_DC_fiducial_plot/DC_hit_position_region1_hadron_cut_50");
  TH2F *h_km_dc_pos_r3 = (TH2F*)fIn->Get("FD_PID_hadron_DC_fiducial_plot/DC_hit_position_region3_hadron_cut_50");

  // Histogram for hadron with -321 from EB, pos charge, and DCR1 Fiducial and DCR3 respectively
  TH2F *h_km_dc_pos_r1_cut = (TH2F*)fIn->Get("FD_PID_hadron_DC_fiducial_plot/DC_hit_position_region1_hadron_cut_52");
  TH2F *h_km_dc_pos_r3_cut = (TH2F*)fIn->Get("FD_PID_hadron_DC_fiducial_plot/DC_hit_position_region1_hadron_cut_53");


  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // PCAL Fiducial hit before cuts with all negative tracks
  TH2F *h_el_pcal_hit_position = (TH2F*)fIn->Get("FD_PID_electron_EC_plots/EC_PCAL_hit_position_cut_00"); 
  // PCAL fiducial hits after fid cut and negative charge and EB PID
  TH2F *h_el_pcal_hit_position_cut = (TH2F*)fIn->Get("FD_PID_electron_EC_plots/EC_PCAL_hit_position_cut_04");


  // Plot information now
  TCanvas *c_dcr1 = new TCanvas("c_dcr1","c_dcr1",800,800);
  h_el_dc_hit_position_r1->SetTitle("Electron DCR1 Fiducial Cut"); 
  h_el_dc_hit_position_r1->Draw();
  h_el_dc_hit_position_r1_cut->Draw("same+colz");
  h_el_dc_hit_position_r1->GetXaxis()->SetTitle("x (cm)");
  h_el_dc_hit_position_r1->GetXaxis()->CenterTitle();
  h_el_dc_hit_position_r1->GetYaxis()->SetTitle("y (cm)");
  h_el_dc_hit_position_r1->GetYaxis()->CenterTitle();
  c_dcr1->SaveAs(Form("h_el_dc_hit_position_r1_%d.pdf",run));

  TCanvas *c_dcr3 = new TCanvas("c_dcr3","c_dcr3",800,800);
  h_el_dc_hit_position_r3->SetTitle("Electron DCR3 Fiducial Cut"); 
  h_el_dc_hit_position_r3->Draw();
  h_el_dc_hit_position_r3_cut->Draw("same+colz");
  h_el_dc_hit_position_r3->GetXaxis()->SetTitle("x (cm)");
  h_el_dc_hit_position_r3->GetXaxis()->CenterTitle();
  h_el_dc_hit_position_r3->GetYaxis()->SetTitle("y (cm)");
  h_el_dc_hit_position_r3->GetYaxis()->CenterTitle();

  //defined by x,y and rad1,rad2
  double d = 400;
  double h = 200;
  double ang = 15.0 * (3.1415/180.0);
  double tilt = 5.0 * (3.1415/180.0);
  double p = fabs((d*TMath::Tan(ang) - h) / (TMath::Cos(tilt) - TMath::Tan(ang)*TMath::Sin(tilt)));
  std::cout << " d " << d << " h " << h << " ang " << ang << " tilt " << tilt << " p " << p << std::endl;
  TEllipse *el1 = new TEllipse(0.0, 0.0 , p, p );
  //el1->SetFillColorAlpha(kRed,0.1);
  //el1->Draw("same");
  

  c_dcr3->SaveAs(Form("h_el_dc_hit_position_r3_%d.pdf",run));

  if ( plotHadrons ){
    // Kaon Plus
    TCanvas *c_kpdc = new TCanvas("c_kpdc","c_kpdc", 800, 800);
    h_kp_dc_pos_r1->Draw();
    h_kp_dc_pos_r1_cut->Draw("same+colz");
    h_kp_dc_pos_r1->GetXaxis()->SetTitle("x (cm)");
    h_kp_dc_pos_r1->GetXaxis()->CenterTitle();
    h_kp_dc_pos_r1->GetYaxis()->SetTitle("y (cm)");
    h_kp_dc_pos_r1->GetYaxis()->CenterTitle();
    c_kpdc->SaveAs(Form("h_kp_dc_pos_r1_%d.pdf",run));

    TCanvas *c_kpdc3 = new TCanvas("c_kpdc3","c_kpdc3", 800, 800);
    h_kp_dc_pos_r3->Draw();
    h_kp_dc_pos_r3_cut->Draw("same+colz");
    h_kp_dc_pos_r3->GetXaxis()->SetTitle("x (cm)");
    h_kp_dc_pos_r3->GetXaxis()->CenterTitle();
    h_kp_dc_pos_r3->GetYaxis()->SetTitle("y (cm)");
    h_kp_dc_pos_r3->GetXaxis()->CenterTitle();

    c_kpdc3->SaveAs(Form("h_kp_dc_pos_r3_%d.pdf",run));

    // Kaon Minus
    TCanvas *c_kmdc = new TCanvas("c_kmdc","c_kmdc", 800, 800);
    h_km_dc_pos_r1->Draw();
    h_km_dc_pos_r1_cut->Draw("same+colz");
    h_km_dc_pos_r1->GetXaxis()->SetTitle("x (cm)");
    h_km_dc_pos_r1->GetXaxis()->CenterTitle();
    h_km_dc_pos_r1->GetYaxis()->SetTitle("y (cm)");
    h_km_dc_pos_r1->GetXaxis()->CenterTitle();

    c_kmdc->SaveAs(Form("h_km_dc_pos_r1_%d.pdf",run));

    TCanvas *c_kmdc3 = new TCanvas("c_kmdc3","c_kmdc3", 800, 800);
    h_km_dc_pos_r3->Draw();
    h_km_dc_pos_r3_cut->Draw("same+colz");
    h_km_dc_pos_r3->GetXaxis()->SetTitle("x (cm)");
    h_km_dc_pos_r3->GetXaxis()->CenterTitle();
    h_km_dc_pos_r3->GetYaxis()->SetTitle("y (cm)");
    h_km_dc_pos_r3->GetYaxis()->CenterTitle();
    c_kmdc3->SaveAs(Form("h_km_dc_pos_r3_%d.pdf",run));




    // Proton
    TCanvas *c_prdc = new TCanvas("c_prdc","c_prdc", 800, 800);
    h_pr_dc_pos_r1->Draw();
    h_pr_dc_pos_r1_cut->Draw("same+colz");
    h_pr_dc_pos_r1->GetXaxis()->SetTitle("x (cm)");
    h_pr_dc_pos_r1->GetXaxis()->CenterTitle();
    h_pr_dc_pos_r1->GetYaxis()->SetTitle("y (cm)");
    h_pr_dc_pos_r1->GetYaxis()->CenterTitle();
    c_prdc->SaveAs(Form("h_pr_dc_pos_r1_%d.pdf",run));

    TCanvas *c_prdc3 = new TCanvas("c_prdc3","c_prdc3", 800, 800);
    h_pr_dc_pos_r3->Draw();
    h_pr_dc_pos_r3_cut->Draw("same+colz");
    h_pr_dc_pos_r3->GetXaxis()->SetTitle("x (cm)");
    h_pr_dc_pos_r3->GetXaxis()->CenterTitle();
    h_pr_dc_pos_r3->GetYaxis()->SetTitle("y (cm)");
    h_pr_dc_pos_r3->GetYaxis()->CenterTitle();
    c_prdc3->SaveAs(Form("h_pr_dc_pos_r3_%d.pdf",run));
  }

  // Electron PCAL hit position
  TCanvas *c_pcal = new TCanvas("c_pcal","c_pcal", 800, 800);
  h_el_pcal_hit_position->SetTitle("Electron PCAL Fiducial Cut"); 
  h_el_pcal_hit_position->Draw();
  h_el_pcal_hit_position_cut->Draw("same+colz");
  h_el_pcal_hit_position->GetXaxis()->SetTitle("x (cm)");
  h_el_pcal_hit_position->GetXaxis()->CenterTitle();
  h_el_pcal_hit_position->GetYaxis()->SetTitle("y (cm)");
  h_el_pcal_hit_position->GetYaxis()->CenterTitle();
  c_pcal->SaveAs(Form("h_el_pcal_hit_position_%d.pdf",run));

  return 0;
}
