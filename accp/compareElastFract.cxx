#include <TCanvas.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <string>
#include <vector>
#include <map>
#include <sstream>

int compareElastFract(){

  double el_mass = 0.000511;
  double pr_mass = 0.937282;
  double kp_mass = 0.493;
  double km_mass = 0.493;
  double phi_mass = 1.01946;
  TLorentzVector lv_ebeam(0, 0, 2.221, 2.221 );
  TLorentzVector lv_target(0,0,0,pr_mass);
  
  int el_fract = 8;

  TFile *fIn2 = new TFile("outCheckLund_fract0p2.root","");
  TFile *fIn4 = new TFile("outCheckLund_fract0p4.root","");
  TFile *fIn6 = new TFile("outCheckLund_fract0p6.root","");
  TFile *fIn8 = new TFile("outCheckLund_fract0p8.root","");
  TFile *fIn9 = new TFile("outCheckLund_fract0p9.root","");


  TH1D *h_el_w_2 = (TH1D*)fIn2->Get("h_el_w_f0p2");
  TH1D *h_el_w_4 = (TH1D*)fIn4->Get("h_el_w_f0p4");
  TH1D *h_el_w_6 = (TH1D*)fIn6->Get("h_el_w_f0p6");
  TH1D *h_el_w_8 = (TH1D*)fIn8->Get("h_el_w_f0p8");
  TH1D *h_el_w_9 = (TH1D*)fIn9->Get("h_el_w_f0p9");


  TCanvas *c1 = new TCanvas("c1","c1",900,900);
  gStyle->SetOptStat(00000);
  c1->cd(0);
  h_el_w_2->SetLineColor(kRed);
  h_el_w_4->SetLineColor(kBlue);
  h_el_w_6->SetLineColor(kBlack);
  h_el_w_8->SetLineColor(kGreen-3);
  h_el_w_9->SetLineColor(kPink);

  gPad->SetLogy();
  h_el_w_2->Draw();
  h_el_w_4->Draw("same");
  h_el_w_6->Draw("same");
  h_el_w_8->Draw("same");
  h_el_w_9->Draw("same");
  
  
  TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->AddEntry( h_el_w_2," fract 0.2" );
  leg->AddEntry( h_el_w_4," fract 0.4" );
  leg->AddEntry( h_el_w_6," fract 0.6" );
  leg->AddEntry( h_el_w_8," fract 0.8" );
  leg->AddEntry( h_el_w_9," fract 0.9" );
  leg->Draw("same");

  c1->SaveAs("h_el_w_compare_fract.pdf");


  return 0;
}
