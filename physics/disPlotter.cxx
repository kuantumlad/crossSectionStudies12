#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <vector>
#include <iostream>
#include <TCanvas.h>
#include <TStyle.h>

void disPlotter(const char* infile, const char* config, int run){

  TFile *fin = new TFile(infile);

  gStyle->SetOptStat(0000000);

  TCanvas *c1 = new TCanvas("c1","c1",1800,1200);
  c1->Divide(3,2);

  for( int s = 1; s <= 6; s++ ){
    
    c1->cd(s);
    //c1->SetGrid();
    TH2F* h_el_q2w = (TH2F*)fin->Get(Form("kinematics/h_el_q2w_s%d",s));
    TH2F* h_el_q2w_final = (TH2F*)fin->Get(Form("kinematics/h_el_q2w_final_s%d",s));

    h_el_q2w_final->SetTitle(Form("Q2 vs W Sector %d",s));
    h_el_q2w_final->GetXaxis()->SetTitle("W (GeV)");
    h_el_q2w_final->GetYaxis()->SetTitle("Q^{2} (GeV)");    
    h_el_q2w_final->Draw("colz");
    
    double w_bins = h_el_q2w->GetXaxis()->GetNbins();
    double q_bins = h_el_q2w->GetYaxis()->GetNbins();

    double w_min = h_el_q2w->GetXaxis()->GetBinLowEdge(1);
    double w_max = h_el_q2w->GetXaxis()->GetBinUpEdge(w_bins);
    double w_width = (w_max-w_min)/w_bins;

    double q_min = h_el_q2w->GetYaxis()->GetBinLowEdge(1);
    double q_max = h_el_q2w->GetYaxis()->GetBinUpEdge(q_bins);
    double q_width = (q_max-q_min)/q_bins;
    
    for( int ww = 1; ww<w_bins; ww++ ){
      double w_x = ww*w_width + w_min;
      TLine *w_l  = new TLine(w_x, q_min, w_x, q_max);
      w_l->SetLineColorAlpha(17,0.35);
      w_l->SetLineStyle(1);
      w_l->Draw("same");
    }

    for( int qq = 1; qq<q_bins; qq++ ){
      double q_x = qq*q_width + q_min;
      TLine *q_l  = new TLine(w_min, q_x, w_max, q_x);
      q_l->SetLineColorAlpha(17, 0.35);
      q_l->SetLineStyle(1);
      q_l->Draw("same");
    }

  }

  c1->SaveAs(Form("h_dis_q2w_final_sectors_%s.pdf",config));


  //// plot w spectrum 
  TCanvas *c2 = new TCanvas("c2","c2",1800,1200);
  c2->Divide(3,2);

  for( int s = 1; s <= 6; s++ ){
    
    c2->cd(s);
    TH1F* h_el_w = (TH1F*)fin->Get(Form("kinematics/h_el_w_s%d",s));
    h_el_w->SetTitle(Form("W Sector %d",s));
    h_el_w->GetXaxis()->SetTitle("W (GeV)");
    h_el_w->Draw();

    c2->Update();
    TVirtualPad *p1 = c2->GetPad(1);
    double ymax = gPad->GetUymax();
    double ymin = gPad->GetUymin();
    std::cout << ymax << std::endl;
    TLine *l_deltaP = new TLine(1.232,ymin, 1.232, ymax);
    TLine *l_nP = new TLine(1.520,ymin, 1.520, ymax);
    TLine *l_nPP = new TLine(1.680,ymin, 1.680, ymax);
    l_deltaP->SetLineColor(kRed);
    l_deltaP->SetLineWidth(1);
    l_deltaP->SetLineStyle(2);

    l_nP->SetLineColor(kRed);
    l_nP->SetLineWidth(1);
    l_nP->SetLineStyle(2);

    l_nPP->SetLineColor(kRed);
    l_nPP->SetLineWidth(1);
    l_nPP->SetLineStyle(2);

    l_deltaP->Draw("same");
    l_nP->Draw("same");
    l_nPP->Draw("same");
  }

  c2->SaveAs(Form("h_el_w_dis_%s.pdf",config));


  return 0;


}
