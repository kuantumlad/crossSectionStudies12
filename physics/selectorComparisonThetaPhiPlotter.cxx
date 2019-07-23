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


// used to compare the electron hit position in data and simulation for a given run and magnetic field configuration
// TO USE :

//selectorComparisonPlotter("elastic_out_clas12_002587.10.99.root","mc_selector_clas12_2GeV.0.100_tm06sm06.root","sim_elastic_clas12_2GeV.0.100_tm06sm06.root",2587,"tm06sm06")


int selectorComparisonThetaPhiPlotter( const char* inFileData, const char* inFileSim, int run , const char* field_config ){

  TFile *fData; 
  TFile *fSim;
  Char_t tmpstr[80];
  Double_t fraction;

  fData = new TFile(inFileData,"");   // Input File
  fSim = new TFile(inFileSim,"");
  if(fData->IsZombie() || fSim->IsZombie()){   // Check if TFile exists!
    cout<<"Input file doesn't exist!" << endl;
    cout<<"Exit program" << endl;
    return 0;
  }

  cout << "Reading from File: " << inFileData << " and " << " "  << inFileSim << endl;
  TCanvas *c1 = new TCanvas("c1","c1",900,900);
  c1->Divide(4,5);


  TH1F *h_p_temp = (TH1F*)fData->Get("acceptance/h_p_bins");
  int n_p_bins = h_p_temp->GetXaxis()->GetNbins();
  for( int bb = 10; bb < n_p_bins; bb++ ){
    gPad->SetLogz();

    TH2D *h_theta_phi_per_p = (TH2D*)fData->Get(Form("acceptance/h_el_theta_phi_per_mntm_b%d",bb));
    TH2D *h_theta_phi_per_p_sim = (TH2D*)fSim->Get(Form("acceptance/h_el_theta_phi_per_mntm_b%d",bb));

    ////////////////////////////
    //get a phi theta per sector
    //then projection of the phi axis for each theta bin
    
    int theta_bins = h_theta_phi_per_p->GetYaxis()->GetNbins();
    int phi_bins = h_theta_phi_per_p->GetXaxis()->GetNbins();
    TCanvas *c_tb_temp = new TCanvas(Form("c_tb_temp%d",bb),Form("c_tb_temp%d",bb),900,900);
    c_tb_temp->Divide(5,6,0,0);

    for( int tb = 1; tb < theta_bins; tb++ ){
      TH1D *h_proj_phi = new TH1D(Form("p%d_tb%d",bb,tb),Form("p%d_tb%d",bb,tb), phi_bins, -180.0, 180.0);
      TH1D *h_proj_phi_sim = new TH1D(Form("p%d_tb%d_sim",bb,tb),Form("p%d_tb%d_sim",bb,tb), phi_bins, -180.0, 180.0);
      c_tb_temp->cd(tb);
      for( int pb = 1; pb < phi_bins; pb++){		
	h_proj_phi->SetBinContent(pb, h_theta_phi_per_p->GetBinContent(pb,tb));
	h_proj_phi_sim->SetBinContent(pb, h_theta_phi_per_p_sim->GetBinContent(pb,tb));
      }      

      if( !(h_proj_phi->GetEntries() > 0 ) && !(h_proj_phi_sim->GetEntries() > 0 ) ) continue;
      std::cout << "  h_proj_phi_sim->GetMaximum() " <<  h_proj_phi_sim->GetMaximum() << " data " << h_proj_phi->GetMaximum() <<std::endl;
      double phi_scaler = h_proj_phi_sim->GetMaximum()/h_proj_phi->GetMaximum();
      std::cout << " scale value " << phi_scaler << std::endl;
      h_proj_phi_sim->Scale(1.0/phi_scaler);
      h_proj_phi->SetLineColor(kRed);      
      h_proj_phi->Draw("hist");
      h_proj_phi_sim->Draw("hist+same");
    }
    c_tb_temp->SaveAs(Form("comparison_h_el_phi_per_theta_p%d_%s.pdf",bb, field_config));

    c1->cd(bb+1);
    h_theta_phi_per_p->Draw("colz");
    h_theta_phi_per_p_sim->Draw("same");
  }
  c1->SaveAs(Form("comparison_h_el_theta_phi_per_p_%s.pdf",field_config));


  
  
  

  return 0;

}
