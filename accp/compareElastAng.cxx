#include <TCanvas.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <string>
#include <vector>
#include <map>
#include <sstream>

int compareElastAng(){

  double el_mass = 0.000511;
  double pr_mass = 0.937282;
  double kp_mass = 0.493;
  double km_mass = 0.493;
  double phi_mass = 1.01946;
  TLorentzVector lv_ebeam(0, 0, 2.221, 2.221 );
  TLorentzVector lv_target(0,0,0,pr_mass);
  
  int el_fract = 9;

  TFile *fIn0 = new TFile("outCheckLund_fract0p1112.root","");
  TFile *fIn1 = new TFile("outCheckLund_fract0p2021.root","");

  TH1D *h_el_w_0 = (TH1D*)fIn0->Get("h_el_w_f0p1112");
  TH1D *h_el_w_1 = (TH1D*)fIn1->Get("h_el_w_f0p2021");

  TCanvas *c1 = new TCanvas("c1","c1",900,900);
  gStyle->SetOptStat(00000);
  c1->cd(0);
  h_el_w_0->SetLineColor(kRed);
  h_el_w_1->SetLineColor(kBlue);

  gPad->SetLogy();
  h_el_w_0->Draw();
  h_el_w_1->Draw("same");

  TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->AddEntry(h_el_w_0,"ang 11-12" );
  leg->AddEntry(h_el_w_1,"ang 20-21" );
  leg->Draw("same");

  c1->SaveAs("h_el_w_compare_ang.pdf");

  return 0;
}
