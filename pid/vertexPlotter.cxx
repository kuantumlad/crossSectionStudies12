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

int vertexPlotter(const char* input, int run){

  TFile *fIn = new TFile(input,"");

  double p_min = 1.5;
  double Ebeam = 2.221;

  if(Ebeam > 10) p_min = 1.5;
  if(Ebeam < 10) p_min = 1.0;
  if(Ebeam < 3)  p_min = 0.5;

  double sigma_range = 3;

  /// //////////////////////////////////////////////////////////////////////////////////////////////////////
  /// a) cut based on sampling fraction versus drift chamber momentum

  
  TCanvas *c_el_vz = new TCanvas("c_el_vz", "c_el_vz", 1800, 900);
  c_el_vz->Divide(3,2);

  for( int s = 1; s <= 6; s++ ){
    c_el_vz->cd(s);
    
    //    TCanvas *c_temp_vz = new TCanvas(Form("c_temp_vz_s%d",s),Form("c_temp_vz_s%d",s),900,600);
    TH1F *h_el_vz = (TH1F*)fIn->Get(Form("FD_PID_electron_DC_plots/DC_z_vertex_sec%d_cut_01",s));

    // get hist parameters to set fit limits

    double xlow,xhigh,histmax;
    int binlow,binhigh,binmax;
    binmax = h_el_vz->GetMaximumBin();
    histmax = h_el_vz->GetMaximum();
    binlow=binmax;
    binhigh=binmax;

    // The 0.65 parameter can be changed, this basically means start at the peak and work your way left and right
    // until you've gotten to 65% of the max amplitude.
    while(h_el_vz->GetBinContent(binhigh++) >= .65*histmax&&binhigh<=h_el_vz->GetNbinsX()){};
    while(h_el_vz->GetBinContent(binlow--) >= .65*histmax&&binlow>=1){};
    
    xlow = h_el_vz->GetBinLowEdge(binlow);
    xhigh = h_el_vz->GetBinLowEdge(binhigh+1);

    h_el_vz->Fit("gaus","","",xlow,xhigh);
    //get the fitted gaus to extract parameters from
    TF1 *fit_vz = (TF1*)h_el_vz->GetListOfFunctions()->FindObject("gaus");    
    double fit_mean = fit_vz->GetParameter(1);
    double fit_sig = fit_vz->GetParameter(2);
    
    std::cout << " sector " << s << "fitted vertex mean " << fit_mean << " fitted sigma " << fit_sig << std::endl;
    
    h_el_vz->SetTitle(Form("Sector %d",s));
    h_el_vz->GetXaxis()->CenterTitle();
    gStyle->SetOptStat(0);
    h_el_vz->SetLineColor(kBlack);
    h_el_vz->Draw();
    
    c_el_vz->Update();
    TVirtualPad *p1 = c_el_vz->GetPad(1);
    double ymax = gPad->GetUymax();
    double ymin = gPad->GetUymin();
    std::cout << ymax << std::endl;
    //plot the 0 vertex line 
    TLine *l = new TLine(0.0,ymin, 0.0, ymax);
    l->SetLineColor(kRed);
    l->SetLineWidth(1);
    l->SetLineStyle(2);
    
    // use 2.5 from fitted center to each side of the target (5cm long target) 
    double vz_min = fit_mean - 2.5;//sigma_range*fit_sig;
    double vz_max = fit_mean + 2.5;// sigma_range*fit_sig;
    
    TBox *box_vz = new TBox(vz_min, ymin, vz_max, ymax );
    box_vz->SetFillColorAlpha(kBlue-4,0.15);
    box_vz->Draw("same");


    std::cout << " title " <<  h_el_vz->GetTitle() << std::endl;
    h_el_vz->Draw("same");
    l->Draw("same");
    
  }

  c_el_vz->SaveAs(Form("vzCutRegionRun%d.pdf",run));

  return 0;
}
