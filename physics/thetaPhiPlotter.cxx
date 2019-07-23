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


int thetaPhiPlotter( const char *inFileData, const char* field_config ){

  TFile *fData;
  fData = new TFile(inFileData,"");   // Input File


  if(fData->IsZombie() ){return 0; }


  TCanvas *c1 = new TCanvas("c1","c1",9000,9000);
  c1->Divide(4,6);
  
  TH1F *h_p_temp = (TH1F*)fData->Get("acceptance/h_p_bins");
  int n_p_bins = h_p_temp->GetXaxis()->GetNbins();
  for( int bb = 10; bb < n_p_bins; bb++ ){
    gPad->SetLogz();

    TH2D *h_theta_phi_per_p = (TH2D*)fData->Get(Form("acceptance/h_el_theta_phi_per_mntm_b%d",bb));

    ////////////////////////////
    //get a phi theta per sector
    //then projection of the phi axis for each theta bin
    
    int theta_bins = h_theta_phi_per_p->GetYaxis()->GetNbins();
    int phi_bins = h_theta_phi_per_p->GetXaxis()->GetNbins();
    TCanvas *c_tb_temp = new TCanvas(Form("c_tb_temp%d",bb),Form("c_tb_temp%d",bb),900,900);
    c_tb_temp->Divide(5,6,0,0);

    for( int tb = 1; tb < theta_bins; tb++ ){
      TH1D *h_proj_phi = new TH1D(Form("tb%d",tb),Form("tb%d",tb), phi_bins, -180.0, 180.0);
      c_tb_temp->cd(tb);
      for( int pb = 1; pb < phi_bins; pb++){		
	h_proj_phi->SetBinContent(pb, h_theta_phi_per_p->GetBinContent(pb,tb));
      }      
      h_proj_phi->Draw("hist");
    }
    c_tb_temp->SaveAs(Form("h_el_phi_per_theta_p%d_%s",bb, field_config));

    c1->cd(bb+1);
    h_theta_phi_per_p->Draw("colz");
  }
  c1->SaveAs(Form("h_el_theta_phi_per_p_%s.pdf",field_config));

  TH1F* h_vz_temp =  (TH1F*)fData->Get("acceptance/h_vz_bins");
  int n_vz_bins = h_vz_temp->GetXaxis()->GetNbins();

  TCanvas *c2 = new TCanvas("c2","c2",900,900);
  c2->Divide(5,10);
  for( int bb = 0; bb < n_vz_bins; bb++ ){
    c2->cd(bb+1);
    gPad->SetLogz();
    TH2D *h_theta_phi_per_vz = (TH2D*)fData->Get(Form("acceptance/h_el_theta_phi_per_vz_b%d",bb));
    h_theta_phi_per_vz->Draw("colz");
  }
  c2->SaveAs(Form("h_el_theta_phi_per_vz_%s.pdf",field_config));
  



  return 0;


}
