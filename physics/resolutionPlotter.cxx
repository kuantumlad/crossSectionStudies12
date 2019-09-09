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
double Ebeam = 7.546;

int resolutionPlotter( const char* inFileMC, const char* config, std::string data_type ){

  TFile *fMC;
  Char_t tmpstr[80];
  Double_t fraction;

  fMC = new TFile(inFileMC,"");
  if( fMC->IsZombie() ){   // Check if TFile exists!
    cout<<"Input file doesn't exist!" << endl;
    cout<<"Exit program" << endl;
    return 0;
  }

  cout << "Reading from File: " << inFileMC << " "  << endl;
  
  
  bool printAll = false;
  bool setMinZero = true;
  double pmin = Ebeam-2.0;
  if( !setMinZero ) pmin = 0.0;

  std::vector< TH1F* > v_el_res_p;
  std::vector< TH1F* > v_el_res_theta;
  std::vector< TH1F* > v_el_res_phi;
  std::vector< TH2F* > v_el_res_phi_vs_p;


  if( data_type == "SIM" ){
    for( int ss = 0; ss < 6; ss++ ){
      v_el_res_theta.push_back( (TH1F*)fMC->Get(Form("kin_resolution/h_res_el_theta_s%d",ss)));
      v_el_res_phi.push_back( (TH1F*)fMC->Get(Form("kin_resolution/h_res_el_phi_s%d",ss)));
      v_el_res_phi_vs_p.push_back( (TH2F*)fMC->Get(Form("kin_resolution/h2_res_el_phi_vs_p_s%d",ss)));   
    }

    TCanvas *c0 = new TCanvas("c0","c0",650,975);
    c0->Divide(2,3);
    for( int ss = 0; ss<6; ss++ ){
      c0->cd(ss+1);
      v_el_res_theta[ss]->SetTitle(Form("#theta Res. S%d",ss+1));
      v_el_res_theta[ss]->GetXaxis()->SetTitle("#theta_{rec} - #theta_{gen} (deg)");
      v_el_res_theta[ss]->GetXaxis()->CenterTitle();   
      v_el_res_theta[ss]->Draw();
      v_el_res_theta[ss]->GetXaxis()->SetRangeUser(-0.2,0.2);
      v_el_res_theta[ss]->Draw();

      c0->Update();
      TVirtualPad *p1 = c0->GetPad(1);
      double ymax = gPad->GetUymax();
      double ymin = gPad->GetUymin();
      std::cout << ymax << std::endl;
      TLine *l_z = new TLine(0.0,ymin, 0.0, ymax);
      l_z->SetLineColor(kRed);
      l_z->SetLineWidth(1);
      l_z->SetLineStyle(2);
      l_z->Draw("same");


    }
    c0->SaveAs(Form("h_el_resolution_theta_sect_%s.pdf",config));

    TCanvas *c0a = new TCanvas("c0a","c0a",650,975);
    c0a->Divide(2,3);
    for( int ss = 0; ss<6; ss++ ){
      c0a->cd(ss+1);
      v_el_res_phi[ss]->SetTitle(Form("#phi Res. S%d",ss+1));
      v_el_res_phi[ss]->GetXaxis()->SetTitle("#phi_{rec} - #phi_{gen} (deg)");
      v_el_res_phi[ss]->GetXaxis()->CenterTitle();   
      v_el_res_phi[ss]->Draw();

      c0a->Update();
      TVirtualPad *p1 = c0a->GetPad(1);
      double ymax = gPad->GetUymax();
      double ymin = gPad->GetUymin();
      std::cout << ymax << std::endl;
      TLine *l_z = new TLine(0.0,ymin, 0.0, ymax);
      l_z->SetLineColor(kRed);
      l_z->SetLineWidth(1);
      l_z->SetLineStyle(2);
      l_z->Draw("same");


    }
    c0a->SaveAs(Form("h_el_resolution_phi_sect_%s.pdf",config));

    TCanvas *c0b = new TCanvas("c0b","c0b",650,975);
    c0b->Divide(2,3);
    for( int ss = 0; ss<6; ss++ ){
      c0b->cd(ss+1);
      gPad->SetLogz();
      v_el_res_phi_vs_p[ss]->SetTitle(Form("#phi vs p Res. S%d",ss+1));
      v_el_res_phi_vs_p[ss]->GetYaxis()->SetTitle("#phi_{rec} - #phi_{gen} (deg)");
      v_el_res_phi_vs_p[ss]->GetYaxis()->CenterTitle();   
      v_el_res_phi_vs_p[ss]->GetXaxis()->SetTitle("p (GeV)");
      v_el_res_phi_vs_p[ss]->GetXaxis()->CenterTitle();   
      v_el_res_phi_vs_p[ss]->Draw("colz");
      v_el_res_phi_vs_p[ss]->GetXaxis()->SetRangeUser(Ebeam-2.0,Ebeam+0.2);
      v_el_res_phi_vs_p[ss]->Draw("colz");
    }
    c0b->SaveAs(Form("h_el_resolution_phi_vs_p_sect_%s.pdf",config));
  }

  if( data_type == "DATA"){
    // extract resolution from data
    TH1F* h_el_res_data_theta = (TH1F*)fMC->Get("hist_all_el_delta_theta");
    TH1F* h_el_res_data_energy = (TH1F*)fMC->Get("hist_all_el_energy");
    TH1F* h_el_res_data_theta_p = (TH1F*)fMC->Get("hist_all_el_delta_theta_p");
    //TH1F* h_el_res_data_p_theta = (TH1F*)fMC->Get("hist_all_el_delta_p_theta");

    std::vector< TH2F* > v_el_res_data_delta_theta;
    std::vector< TH2F* > v_el_res_data_delta_p;
    std::vector< TH2F* > v_el_res_data_delta_theta_final;
    std::vector< TH2F* > v_el_res_data_delta_p_final;
    std::vector< TH1F* > v_el_meas_beam_sect;
    
    TH2F *h_el_delta_theta_theta = (TH2F*)fMC->Get("hist_all_el_delta_theta_theta");
    TH1F *h_el_meas_beam = (TH1F*)fMC->Get("hist_all_el_meas_beam");

    for( int ss=0; ss < 6; ss++ ){
      v_el_res_data_delta_theta.push_back( (TH2F*)fMC->Get(Form("h_el_delta_theta_s%d",ss)));
      v_el_res_data_delta_p.push_back( (TH2F*)fMC->Get(Form("h_el_delta_p_s%d",ss)));
      v_el_res_data_delta_theta_final.push_back( (TH2F*)fMC->Get(Form("h_el_delta_theta_s%d_final",ss)));
      v_el_res_data_delta_p_final.push_back( (TH2F*)fMC->Get(Form("h_el_delta_p_s%d_final",ss)));
      v_el_meas_beam_sect.push_back( (TH1F*)fMC->Get(Form("h_el_meas_beam_s%d_final",ss)) );
    }


    TCanvas *c1 = new TCanvas("c1","c1",900,900);
    c1->cd(0);
    h_el_res_data_theta->SetTitle("Resolution: #theta_{calc} - #theta_{meas}");
    h_el_res_data_theta->GetXaxis()->SetTitle("#theta_{calc} - #theta_{meas} (deg)");
    TF1 *fit_g = new TF1("my_guas","gaus",-0.9, 0.5);
    fit_g->SetParameter(0,h_el_res_data_theta->GetMaximum());
    fit_g->SetParameter(1,0.0);
    fit_g->SetParameter(2,0.5);
    h_el_res_data_theta->Fit("my_guas","R");
    double delta_theta_mean = fit_g->GetParameter(1);
    double delta_theta_sig = fit_g->GetParameter(2);
    std::cout << " delta theta mean " << delta_theta_mean << " delta theta sig " << delta_theta_sig << std::endl;
    h_el_res_data_theta->Draw();
    gPad->Modified();
    h_el_res_data_theta->GetXaxis()->SetRangeUser(-5.0, 2.0);
    c1->SaveAs(Form("h_res_data_theta_%s.pdf",config));

    TCanvas *c1a = new TCanvas("c1a","c1a",900,900);
    c1a->cd(0);
    h_el_res_data_energy->SetTitle("Resolution: E_{calc} - E_{meas}");
    h_el_res_data_energy->GetXaxis()->SetTitle("E_{calc} - E_{meas} (GeV)");
    TF1 *fit_g2 = new TF1("my_guas2","gaus",-0.1, 0.07);
    fit_g2->SetParameter(0,h_el_res_data_energy->GetMaximum());
    fit_g2->SetParameter(1,0.0);
    fit_g2->SetParameter(2,0.1);
    h_el_res_data_energy->Fit("my_guas2","R");
    double delta_e_mean = fit_g->GetParameter(1);
    double delta_e_sig = fit_g->GetParameter(2);
    h_el_res_data_energy->Draw();
    c1a->SaveAs(Form("h_res_data_energy_%s.pdf",config));
    std::cout << " delta E mean " << delta_e_mean << " delta E sig " << delta_e_sig << std::endl;

    TCanvas *c1b = new TCanvas("c1b","c1b",900,900);
    c1b->cd(0);
    h_el_res_data_theta_p->SetTitle("Resolution: #theta_{calc} - #theta_{meas} vs p");
    h_el_res_data_theta_p->GetYaxis()->SetTitle("#theta_{calc} - #theta_{meas} (deg)");
    h_el_res_data_theta_p->GetXaxis()->SetTitle("p (GeV)");
    gPad->SetLogz();
    h_el_res_data_theta_p->Draw("colz");
    c1b->SaveAs(Form("h_res_data_theta_vs_p_%s.pdf",config));

    TCanvas *c1c = new TCanvas("c1c","c1c",900,900);
    c1c->Divide(2,3);
    for( int ss = 0; ss < 6; ss++ ){
      c1c->cd(ss+1);
      v_el_res_data_delta_theta[ss]->SetTitle(Form("#Delta #theta vs #theta S%d",ss));
      v_el_res_data_delta_theta[ss]->GetXaxis()->SetTitle("#theta (deg)");
      v_el_res_data_delta_theta[ss]->GetYaxis()->SetTitle("#Delta #theta (deg)");
      v_el_res_data_delta_theta[ss]->Draw("colz");
 
      TLine *l0 =  new TLine(0.0, 0.0, v_el_res_data_delta_theta[ss]->GetXaxis()->GetXmax(), 0.0);
      l0->SetLineColor(kRed);
      l0->SetLineWidth(3);
      l0->SetLineStyle(1);
      l0->Draw("same");
    }
    c1c->SaveAs(Form("h_res_data_p_vs_theta_%s.pdf",config));

    TCanvas *c1d = new TCanvas("c1d","c1d",900,900);
    c1d->Divide(2,3);
    for( int ss = 0; ss < 6; ss++ ){
      c1d->cd(ss+1);
      v_el_res_data_delta_p[ss]->SetTitle(Form("#Delta p vs p S%d",ss));
      v_el_res_data_delta_p[ss]->GetXaxis()->SetTitle("p (GeV)");
      v_el_res_data_delta_p[ss]->GetYaxis()->SetTitle("#Delta p (GeV)");
      v_el_res_data_delta_p[ss]->Draw("colz");
 
      TLine *l0 =  new TLine(0.0, 0.0, v_el_res_data_delta_p[ss]->GetXaxis()->GetXmax(), 0.0);
      l0->SetLineColor(kRed);
      l0->SetLineWidth(3);
      l0->SetLineStyle(1);
      l0->Draw("same");
    }
    c1d->SaveAs(Form("h_res_data_p_vs_p_sectors_%s.pdf",config));

    /////////////////////////////////
    //plot the final resolution plots after applying the W cut
    TCanvas *c1e = new TCanvas("c1e","c1e",900,900);
    c1e->Divide(2,3);
    for( int ss = 0; ss < 6; ss++ ){
      c1e->cd(ss+1);
      v_el_res_data_delta_theta_final[ss]->SetTitle(Form("Final Sample #Delta #theta vs #theta S%d",ss));
      v_el_res_data_delta_theta_final[ss]->GetXaxis()->SetTitle("#theta (deg)");
      v_el_res_data_delta_theta_final[ss]->GetYaxis()->SetTitle("#Delta #theta (deg)");
      v_el_res_data_delta_theta_final[ss]->Draw("colz");
 
      TLine *l0 =  new TLine(0.0, 0.0, v_el_res_data_delta_theta_final[ss]->GetXaxis()->GetXmax(), 0.0);
      l0->SetLineColor(kRed);
      l0->SetLineWidth(3);
      l0->SetLineStyle(1);
      l0->Draw("same");
    }
    c1e->SaveAs(Form("h_res_data_theta_vs_theta__final_%s.pdf",config));

    TCanvas *c1f = new TCanvas("c1f","c1f",900,900);
    c1f->Divide(2,3);
    for( int ss = 0; ss < 6; ss++ ){
      c1f->cd(ss+1);
      v_el_res_data_delta_p_final[ss]->SetTitle(Form("#Delta p vs p S%d",ss));
      v_el_res_data_delta_p_final[ss]->GetXaxis()->SetTitle("p (GeV)");
      v_el_res_data_delta_p_final[ss]->GetYaxis()->SetTitle("#Delta p (GeV)");
      v_el_res_data_delta_p_final[ss]->Draw("colz");
 
      TLine *l0 =  new TLine(0.0, 0.0, v_el_res_data_delta_p_final[ss]->GetXaxis()->GetXmax(), 0.0);
      l0->SetLineColor(kRed);
      l0->SetLineWidth(3);
      l0->SetLineStyle(1);
      l0->Draw("same");
    }
    c1f->SaveAs(Form("h_res_data_p_vs_p_sectors_%s_final.pdf",config));

    TCanvas *c1g = new TCanvas("c1g","c1g",900,900);
    c1g->cd();
    h_el_delta_theta_theta->SetTitle("#Delta #theta vs #theta All Sect.");
    h_el_delta_theta_theta->GetXaxis()->SetTitle("#theta (deg)");
    h_el_delta_theta_theta->GetYaxis()->SetTitle("#Delta #theta (deg)");
    h_el_delta_theta_theta->Draw("colz");
    TLine *l0t =  new TLine(0.0, 0.0, h_el_delta_theta_theta->GetXaxis()->GetXmax(), 0.0);
    l0t->SetLineColor(kRed);
    l0t->SetLineWidth(3);
    l0t->SetLineStyle(1);
    l0t->Draw("same");
    
    TCanvas *c2 = new TCanvas("c2","c2",900,900);
    c2->Divide(2,3);
    for( int ss = 0; ss < 6; ss++ ){
      v_el_meas_beam_sect[ss]->SetTitle(Form("Meas. Beam Energy S%d",ss));
      v_el_meas_beam_sect[ss]->GetXaxis()->SetTitle("Meas. Beam Energy (GeV)");
      v_el_meas_beam_sect[ss]->Draw();
      c2->Update();
      double ymax = gPad-GetUymax();
      TLine *lb = new TLine(Ebmeam, 0.0, Ebeam,ymax);
      lb->SetLineColor(kRed);
      lb->SetLineWidth(3);
      lb->Draw("same");
    }
    c2->SaveAs(Form("h_el_meas_beam_energy_per_sector_%s",config));
    

    TCanvas *c2a = new TCanvas("c2a","c2a",900,900);
    c2a->Divide(2,3);
    h_el_meas_beam->SetTitle("Meas. Beam All Sectors");
    h_el_meas_beam->GetXaxis()->SetTitle("Measured Beam Energy (GeV)");
    h_el_meas_beam->Draw();
    c2a->Update();
    TLine *lb = new TLine(Ebeam, 0.0, Ebeam, gPad->GetUymax() );
    lb->SetLineColor(kRed);
    lb->SetLineWidth(3);
    lb->Draw("same");
    c2a->SaveAs(Form("h_el_meas_beam_energy_%s",config));

    //h_el_res_data_p_theta->SetTitle("Resolution: p_{calc} - p_{meas} vs theta");
    //h_el_res_data_p_theta->GetYaxis()->SetTitle("#p_{calc} - #p_{meas} (GeV)");
    //h_el_res_data_p_theta->GetXaxis()->SetTitle("#theta (GeV)");
    //gPad->SetLogz();
    //h_el_res_data_thetap->Draw("colz");

    


  }

  
  return 0;
}
