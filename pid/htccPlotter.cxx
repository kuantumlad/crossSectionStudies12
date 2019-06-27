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

int htccPlotter(const char* input, int run){

  TFile *fIn = new TFile(input,"");


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
  c_htcc_2->SaveAs(Form("hist_nphe_%d.pdf",run));

  TCanvas *c_htcc_eb = new TCanvas("c_htcc_eb","c_htcc_eb",900,900);
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
  
  return 0;


}
