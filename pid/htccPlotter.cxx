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
#include "TLine.h"
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

int htccPlotter(const char* input, const char* config, int run){

  TFile *fIn = new TFile(input,"");

  /*
  TCanvas *c_htcc_1 = new TCanvas("c_htcc_1","c_htcc_1",900,900);

  TLine *nphe_cut = new TLine(2.0,0.0,2.0,200);
  nphe_cut->SetLineColor(kRed);
  nphe_cut->SetLineWidth(5);

  gStyle->SetOptStat(0);
  TH2F *h_htcc_nphe_p = (TH2F*)fIn->Get("HTCC_Nphe/hist_HTCC_Nphe_vs_momentum");
  h_htcc_nphe_p->SetTitle("HTCC Nphe vs Momentum");
  gStyle->SetTitleFontSize(0.1);
  c_htcc_1->SetTopMargin(0.14);
  c_htcc_1->SetRightMargin(0.12);

  //h_htcc_nphe_p->RebinY(2);

  h_htcc_nphe_p->GetXaxis()->SetTitle("Momentum [GeV]");
  h_htcc_nphe_p->GetXaxis()->CenterTitle();
  h_htcc_nphe_p->GetXaxis()->SetTitleOffset(0.5);
  h_htcc_nphe_p->GetXaxis()->SetTitleSize(0.06);

  h_htcc_nphe_p->GetYaxis()->SetTitle("Nphe");
  h_htcc_nphe_p->GetYaxis()->CenterTitle();
  h_htcc_nphe_p->GetYaxis()->SetTitleOffset(0.61);
  h_htcc_nphe_p->GetYaxis()->SetTitleSize(0.06);

  h_htcc_nphe_p->Draw("colz");
  nphe_cut->Draw("same");
  c_htcc_1->SaveAs(Form("hist_nphe_vs_p_%d.pdf",run));
*/

  TCanvas *c_htcc_2 = new TCanvas("c_htcc_2","c_htcc_2",900,900);
  TH1D *h_htcc_nphe = ((TH2F*)fIn->Get("HTCC_Nphe/hist_HTCC_Nphe_vs_momentum"))->ProjectionY();
  h_htcc_nphe->SetTitle("HTCC Nphe");
  h_htcc_nphe->GetXaxis()->SetRangeUser(0.0,16);
  h_htcc_nphe->GetXaxis()->SetTitle("Nphe");
  h_htcc_nphe->GetXaxis()->CenterTitle();
  h_htcc_nphe->GetXaxis()->SetTitleOffset(0.5);
  h_htcc_nphe->GetXaxis()->SetTitleSize(0.06);

  h_htcc_nphe->Draw();
  c_htcc_2->Update();
  double y_max = gPad->GetUymax();
  TLine *nphe_cut_2 = new TLine(2.0,0.0,2.0,y_max);
  nphe_cut_2->SetLineColor(kRed);
  nphe_cut_2->SetLineWidth(5);
  nphe_cut_2->Draw("same");
  c_htcc_2->SaveAs(Form("hist_nphe_%s_%d.pdf",config,run));
  

  /*  TCanvas *c_htcc_eb = new TCanvas("c_htcc_eb","c_htcc_eb",900,900);
  c_htcc_eb->Divide(2,3);
  gStyle->SetOptStat(0);
  for( int ss  = 0; ss < 6 ; ss++ ){
    c_htcc_eb->cd(ss+1);
    TH1F *h_nphe_neg = (TH1F*)fIn->Get("FD_PID_electron_HTCC_plots/HTCC_nphe_cut_00");
    h_nphe_neg->GetXaxis()->SetTitle("HTCC nphe");
    h_nphe_neg->GetXaxis()->CenterTitle();
    h_nphe_neg->Draw();

    c_htcc_eb->Update();
    h_nphe_neg->SetTitle("");
    double y_nphe_max = gPad->GetUymax();    
    TLine *nphe_cut = new TLine(2.0, 0.0, 2.0, y_nphe_max);
    nphe_cut->SetLineColor(kRed);
    nphe_cut->SetLineWidth(2);
    nphe_cut->Draw("same");
  }
  c_htcc_eb->SaveAs(Form("hist_nphe_per_sector_r%d.pdf",run));
  */

  // htcc fiducial plots
  TCanvas *c_htcc_3 = new TCanvas("c_htcc_3","c_htcc_3",900,900);
  c_htcc_3->cd();
  TH2F *h_neg_htcc_hit = (TH2F*)fIn->Get("FD_PID_electron_HTCC_plots/HTCC_x_vs_y_cut_00");
  TH2F *h_neg_htcc_hit_wnphe = (TH2F*)fIn->Get("FD_PID_electron_HTCC_plots/HTCC_x_vs_y_wnphe_cut_00");
  TH2F *h_neg_htcc_nphe_hit = new TH2F("h_neg_htcc_nphe_hit","h_neg_htcc_nphe_hit",1000,-150,150, 1000, -150, 150);   
  h_neg_htcc_nphe_hit->Divide(h_neg_htcc_hit_wnphe,h_neg_htcc_hit,1.0,1.0);
  h_neg_htcc_nphe_hit->GetXaxis()->SetTitle("x_{HTCC} (cm)");
  h_neg_htcc_nphe_hit->GetYaxis()->SetTitle("y_{HTCC} (cm)");
  h_neg_htcc_nphe_hit->GetXaxis()->CenterTitle();
  h_neg_htcc_nphe_hit->GetYaxis()->CenterTitle();
  h_neg_htcc_nphe_hit->SetMaximum(28.0);
  h_neg_htcc_nphe_hit->SetMinimum(0.0);
  h_neg_htcc_nphe_hit->Draw("colz");
  c_htcc_3->SaveAs(Form("hist_all_neg_hit_pos_wnphe_%s_r%d.pdf",config,run)); 

  // htcc fiducial plots all cuts
  TCanvas *c_htcc_4 = new TCanvas("c_htcc_4","c_htcc_4",900,900);
  c_htcc_4->cd();
  TH2F *h_final_htcc_hit = (TH2F*)fIn->Get("FD_PID_electron_HTCC_plots/HTCC_x_vs_y_cut_04");
  TH2F *h_final_htcc_hit_wnphe = (TH2F*)fIn->Get("FD_PID_electron_HTCC_plots/HTCC_x_vs_y_wnphe_cut_04");
  TH2F *h_final_htcc_nphe_hit = new TH2F("h_final_htcc_nphe_hit","h_final_htcc_nphe_hit",1000,-150,150, 1000, -150, 150);   
  h_final_htcc_nphe_hit->Divide(h_final_htcc_hit_wnphe,h_final_htcc_hit,1.0,1.0);
  h_final_htcc_nphe_hit->SetTitle("Electron HTCC Hit NPHE PCAL Cuts");
  h_final_htcc_nphe_hit->GetXaxis()->SetTitle("x_{HTCC} (cm)");
  h_final_htcc_nphe_hit->GetYaxis()->SetTitle("y_{HTCC} (cm)");
  h_final_htcc_nphe_hit->GetXaxis()->CenterTitle();
  h_final_htcc_nphe_hit->GetYaxis()->CenterTitle();
  h_final_htcc_nphe_hit->SetMaximum(28.0);
  h_final_htcc_nphe_hit->SetMinimum(0.0);
  h_final_htcc_nphe_hit->Draw("colz");
  c_htcc_4->SaveAs(Form("hist_pcal_el_hit_pos_wnphe_%s_r%d.pdf",config,run)); 

  TCanvas *c_htcc_4a = new TCanvas("c_htcc_4a","c_htcc_4a",900,900);
  c_htcc_4a->cd();
  TH2F *h_final_htcc_hit_dcr1 = (TH2F*)fIn->Get("FD_PID_electron_HTCC_plots/HTCC_x_vs_y_cut_05");
  TH2F *h_final_htcc_hit_wnphe_dcr1 = (TH2F*)fIn->Get("FD_PID_electron_HTCC_plots/HTCC_x_vs_y_wnphe_cut_05");
  TH2F *h_final_htcc_nphe_hit_dcr1 = new TH2F("h_final_htcc_nphe_hit_dcr3","h_final_htcc_nphe_hit_dcr3",1000,-150,150, 1000, -150, 150);   
  h_final_htcc_nphe_hit_dcr1->Divide(h_final_htcc_hit_wnphe_dcr1,h_final_htcc_hit_dcr1,1.0,1.0);
  h_final_htcc_nphe_hit_dcr1->SetTitle("Electron HTCC Hit Nphe DCR1 Cuts");
  h_final_htcc_nphe_hit_dcr1->GetXaxis()->SetTitle("x_{HTCC} (cm)");
  h_final_htcc_nphe_hit_dcr1->GetYaxis()->SetTitle("y_{HTCC} (cm)");
  h_final_htcc_nphe_hit_dcr1->GetXaxis()->CenterTitle();
  h_final_htcc_nphe_hit_dcr1->GetYaxis()->CenterTitle();
  h_final_htcc_nphe_hit_dcr1->SetMaximum(28.0);
  h_final_htcc_nphe_hit_dcr1->SetMinimum(0.0);
  h_final_htcc_nphe_hit_dcr1->Draw("colz");
  c_htcc_4a->SaveAs(Form("hist_dcr1_el_hit_pos_wnphe_%s_r%d.pdf",config,run)); 

  TCanvas *c_htcc_4b = new TCanvas("c_htcc_4b","c_htcc_4b",900,900);
  c_htcc_4b->cd();
  TH2F *h_final_htcc_hit_dcr3 = (TH2F*)fIn->Get("FD_PID_electron_HTCC_plots/HTCC_x_vs_y_cut_06");
  TH2F *h_final_htcc_hit_wnphe_dcr3 = (TH2F*)fIn->Get("FD_PID_electron_HTCC_plots/HTCC_x_vs_y_wnphe_cut_06");
  TH2F *h_final_htcc_nphe_hit_dcr3 = new TH2F("h_final_htcc_nphe_hit_dcr3","h_final_htcc_nphe_hit_dcr3",1000,-150,150, 1000, -150, 150);   
  h_final_htcc_nphe_hit_dcr3->Divide(h_final_htcc_hit_wnphe_dcr3,h_final_htcc_hit_dcr3,1.0,1.0);
  h_final_htcc_nphe_hit_dcr3->SetTitle("Electron HTCC Hit Nphe DCR3 Cuts");
  h_final_htcc_nphe_hit_dcr3->GetXaxis()->SetTitle("x_{HTCC} (cm)");
  h_final_htcc_nphe_hit_dcr3->GetYaxis()->SetTitle("y_{HTCC} (cm)");
  h_final_htcc_nphe_hit_dcr3->GetXaxis()->CenterTitle();
  h_final_htcc_nphe_hit_dcr3->GetYaxis()->CenterTitle();
  h_final_htcc_nphe_hit_dcr3->SetMaximum(28.0);
  h_final_htcc_nphe_hit_dcr3->SetMinimum(0.0);
  h_final_htcc_nphe_hit_dcr3->Draw("colz");
  c_htcc_4b->SaveAs(Form("hist_dcr3_el_hit_pos_wnphe_%s_r%d.pdf",config,run)); 


  // htcc fiducial plots
  TCanvas *c_htcc_5 = new TCanvas("c_htcc_5","c_htcc_5",900,900);
  c_htcc_5->cd();
  TH2F * h_neg_htcc_nphe_hit_zoom = (TH2F*)h_neg_htcc_nphe_hit->Clone();
  h_neg_htcc_nphe_hit_zoom->GetXaxis()->SetRangeUser(-37.0, 37.0);
  h_neg_htcc_nphe_hit_zoom->GetYaxis()->SetRangeUser(-37.0, 37.0);
  h_neg_htcc_nphe_hit_zoom->SetMaximum(28.0);
  h_neg_htcc_nphe_hit_zoom->SetMinimum(0.0);
  h_neg_htcc_nphe_hit_zoom->Draw("colz");
  c_htcc_5->SaveAs(Form("hist_all_neg_hit_pos_wnphe_zoomed_%s_r%d.pdf",config,run)); 
  

  return 0;


}
