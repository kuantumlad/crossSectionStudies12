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


int selectorComparisonPlotter( const char* inFileData, const char* inFileMC, const char* inFileSim, int run , const char* field_config ){

  TFile *fData; 
  TFile *fMC;
  TFile *fSim;
  Char_t tmpstr[80];
  Double_t fraction;

  fData = new TFile(inFileData,"");   // Input File
  fMC = new TFile(inFileMC,"");
  fSim = new TFile(inFileSim,"");
  if(fData->IsZombie() || fMC->IsZombie() || fSim->IsZombie()){   // Check if TFile exists!
    cout<<"Input file doesn't exist!" << endl;
    cout<<"Exit program" << endl;
    return 0;
  }

  cout << "Reading from File: " << inFileData << " and " << inFileMC << " "  << inFileSim << endl;
  
  bool printAll = false;

  TCanvas *c0 = new TCanvas("c0","c0",900,900);
  c0->cd(1);
  TH1D *h_data_final_p = (TH1D*)fData->Get("particle_histograms_selected/hist_electron_p");
  TH1D *h_sim_final_p = (TH1D*)fSim->Get("particle_histograms_selected/hist_electron_p");
  double scale_p_factor = h_data_final_p->GetMaximum()/h_sim_final_p->GetMaximum();
  h_sim_final_p->Scale(scale_p_factor);
  
  h_data_final_p->SetLineColor(kRed); 
  h_sim_final_p->SetLineColor(kBlue-2); 
  h_data_final_p->GetXaxis()->SetTitle("Momentum (GeV)");
  h_data_final_p->GetXaxis()->CenterTitle();

  h_sim_final_p->Draw("HIST");
  h_sim_final_p->GetXaxis()->SetRangeUser(1.25,2.5);
  h_sim_final_p->Draw("HIST");


  h_data_final_p->Draw("HIST SAME");
  h_data_final_p->GetXaxis()->SetRangeUser(1.25,2.5);
  h_data_final_p->Draw("HIST SAME");



  c0->SaveAs(Form("comparison_hist_p_final_r%d_f%s.pdf",run,field_config));

  

  TCanvas *c0a = new TCanvas("c0a","c0a",900,900);
  c0a->cd(1);
  TH1D *h_data_final_theta = (TH1D*)fData->Get("particle_histograms_selected/hist_electron_theta");
  TH1D *h_sim_final_theta = (TH1D*)fSim->Get("particle_histograms_selected/hist_electron_theta");
  double scale_theta_factor = h_data_final_theta->GetMaximum()/h_sim_final_theta->GetMaximum();
  h_sim_final_theta->Scale(scale_theta_factor);
  
  h_data_final_theta->SetLineColor(kRed); 
  h_sim_final_theta->SetLineColor(kBlue-2); 
  h_data_final_theta->GetXaxis()->SetTitle("#theta (deg)");
  h_data_final_theta->GetXaxis()->CenterTitle();
  h_sim_final_theta->Draw("HIST");
  h_data_final_theta->Draw("HIST SAME");
  c0a->SaveAs(Form("comparison_hist_theta_final_r%d_f%s.pdf",run,field_config));


  TCanvas *c0b = new TCanvas("c0b","c0b", 900 ,900);
  c0b->cd(1);
  TH1D *h_data_final_phi = (TH1D*)fData->Get("particle_histograms_selected/hist_electron_phi");
  TH1D *h_sim_final_phi = (TH1D*)fSim->Get("particle_histograms_selected/hist_electron_phi");
  double scale_phi_factor = h_data_final_phi->GetMaximum()/h_sim_final_phi->GetMaximum();
  h_sim_final_phi->Scale(scale_phi_factor);
  
  h_data_final_phi->SetLineColor(kRed); 
  h_sim_final_phi->SetLineColor(kBlue-2); 
  h_data_final_phi->GetXaxis()->SetTitle("#phi (deg)");
  h_data_final_phi->GetXaxis()->CenterTitle();
  h_sim_final_phi->Draw("HIST");
  h_data_final_phi->Draw("HIST SAME");
  c0b->SaveAs(Form("comparison_hist_phi_final_r%d_f%s.pdf",run,field_config));

  

  TCanvas *c1 = new TCanvas("c1","c1",800,1200);
  c1->Divide(2,3);
  for( int s = 0; s < 6; s++ ){
    c1->cd(s+1);
    gStyle->SetOptStat(0000);

    TH1D *h_data_p = (TH1D*)fData->Get(Form("kinematics/h_el_p_s%d_final",s+1));
    TH1D *h_sim_p = (TH1D*)fSim->Get(Form("kinematics/h_el_p_s%d_final",s+1));

    double scale_factor = h_data_p->GetMaximum()/h_sim_p->GetMaximum();
    std::cout << " Total number of entries in p data histogram " << scale_factor << std::endl;
    //h_data_p->Scale(1.0/scale_factor);
    h_sim_p->Scale(scale_factor);

    h_data_p->SetLineColor(kRed);
    h_sim_p->SetLineColor(kBlue-2);

    h_data_p->SetTitle(Form("Electron Momentum S%d",s+1));
    h_data_p->Draw("HIST SAME");
    h_data_p->GetXaxis()->SetRangeUser(1.5,2.5);
    h_data_p->GetXaxis()->SetTitle("Momentum [GeV]");
    h_data_p->GetXaxis()->CenterTitle();
    h_data_p->Draw("HIST SAME");

    h_sim_p->Draw("HIST SAME");
    h_sim_p->GetXaxis()->SetRangeUser(1.5, 2.5); 
    h_sim_p->Draw("HIST SAME");

    TLegend *l1 = new TLegend(0.11,0.7,0.48,0.89);
    l1->SetBorderSize(0);
    l1->AddEntry(h_data_p,"DATA");
    l1->AddEntry(h_sim_p,"SIM");
    l1->Draw();
   
  }
  
  //TCanvas *c2 = new TCanvas("c2","c2",1200,800);
  c1->Update();
  std::string kin_name = "data_sim_mntm";
  if( printAll ){
    c1->Print(Form("data_sim_gen_comparison_r%d_f%s.pdf(",run,field_config),"pdf"); //"h1.pdf(","pdf");
  }
  else{  
    c1->Print(Form("comparison_%s_r%d_f%s.pdf",kin_name.c_str(),run,field_config),"pdf"); //"h1.pdf(","pdf");
  }
  c1->Clear(); 
  c1->Divide(2,3);


  for( int s = 0; s < 6; s++ ){
    c1->cd(s+1);
    gStyle->SetOptStat(0000);

    TH1D *h_data_theta = (TH1D*)fData->Get(Form("kinematics/h_el_theta_s%d_final",s+1));
    TH1D *h_sim_theta = (TH1D*)fSim->Get(Form("kinematics/h_el_theta_s%d_final",s+1));

    double scale_factor = h_data_theta->GetMaximum()/h_sim_theta->GetMaximum();
    std::cout << " Total number of entries in theta data histogram " << scale_factor << std::endl;
    //h_data_theta->Scale(1.0/scale_factor);
    h_sim_theta->Scale(scale_factor);

    h_data_theta->SetLineColor(kRed);
    h_sim_theta->SetLineColor(kBlue-2);

    h_data_theta->SetTitle(Form("Electron Scattering Angle #theta S%d",s+1));
    h_data_theta->Draw("HIST SAME");
    h_data_theta->GetXaxis()->SetTitle("#theta [deg]");
    h_data_theta->GetXaxis()->CenterTitle();
    h_sim_theta->Draw("HIST SAME");

    TLegend *l2 = new TLegend(0.7,0.7,0.9,0.89);
    l2->SetBorderSize(0);
    l2->AddEntry(h_data_theta,"DATA");
    l2->AddEntry(h_sim_theta,"SIM");
    l2->Draw();
    
  }
  c1->Update();

  kin_name="data_sim_theta";
  if( printAll ){
    c1->Print(Form("data_sim_gen_comparison_r%d_f%s.pdf",run,field_config),"pdf"); //"h1.pdf(","pdf");
  }
  else{  
    c1->Print(Form("comparison_%s_r%d_f%s.pdf",kin_name.c_str(),run,field_config),"pdf"); //"h1.pdf(","pdf");
  }

  c1->Clear(); 
  //c1->Divide(3,2);
  //c1->SaveAs(Form("data_sim_gen_comparison_r%d_f%s.pdfa",run,field_config));


  //TCanvas *c3 = new TCanvas("c3","c3",900,900);
  TCanvas *c2 = new TCanvas("c2","c2",900,900);
  c2->Divide(1,1);
  c2->cd(1);
  TH1D *h_data_phi = (TH1D*)fData->Get("particle_histograms_selected/hist_electron_phi");
  TH1D *h_sim_phi = (TH1D*)fSim->Get("particle_histograms_selected/hist_electron_phi");
  TH1D *h_mc_phi = (TH1D*)fMC->Get("hist_mc_all_electron_phi");
  //c3->cd();
  gStyle->SetOptStat(0000);

  double phi_scale_factor = h_data_phi->GetMaximum()/h_sim_phi->GetMaximum(); 
  double mc_phi_scale_factor = h_data_phi->GetMaximum()/h_mc_phi->GetMaximum();
  //h_data_phi->Scale(1.0/phi_scale_factor);
  h_sim_phi->Scale(phi_scale_factor);
  h_mc_phi->Scale(mc_phi_scale_factor);
  
  h_data_phi->SetLineColor(kRed);
  h_sim_phi->SetLineColor(kBlue-2);
  h_mc_phi->SetLineColor(kViolet-4);

  h_data_phi->SetTitle("Reconstructed vs Generated Electron #phi");
  h_data_phi->Draw("HIST SAME");
  h_data_phi->GetYaxis()->SetRangeUser(0.0, 40000.0);
  h_data_phi->GetXaxis()->SetTitle("#phi [deg]");
  h_data_phi->GetXaxis()->CenterTitle();
  h_data_phi->Draw("HIST SAME");

  h_sim_phi->Draw("HIST SAME");
  h_sim_phi->GetYaxis()->SetRangeUser(0.0, 40000.0);
  h_sim_phi->Draw("HIST SAME");

  h_mc_phi->Draw("HIST SAME");
  h_mc_phi->GetYaxis()->SetRangeUser(0.0, 40000.0);
  h_mc_phi->Draw("HIST SAME");

  TLegend *l3 = new TLegend(0.11,0.7,0.48,0.89);
  l3->SetBorderSize(0);
  l3->AddEntry(h_data_phi,"DATA");
  l3->AddEntry(h_sim_phi,"SIM");
  l3->AddEntry(h_mc_phi,"GEN");
  l3->Draw();

  c2->Update();

  kin_name="data_sim_gen_phi";
  if( printAll ){
    c2->Print(Form("data_sim_gen_comparison_r%d_f%s.pdf",run,field_config),"pdf"); //"h1.pdf(","pdf");
  }
  else{  
    c2->Print(Form("comparison_%s_r%d_f%s.pdf",kin_name.c_str(),run,field_config),"pdf"); //"h1.pdf(","pdf");
  }

  c2->Clear(); 

  //  TCanvas *c4 = new TCanvas("c4","c4",900,900);
  c2->Divide(1,1);
  c2->cd(1);
  TH1D *h_data_p = (TH1D*)fData->Get("particle_histograms_selected/hist_electron_p");
  TH1D *h_sim_p = (TH1D*)fSim->Get("particle_histograms_selected/hist_electron_p");
  TH1D *h_mc_p = (TH1D*)fMC->Get("hist_mc_all_electron_p");
  //c4->cd();
  gStyle->SetOptStat(0000);

  double p_scale_factor = h_data_p->GetMaximum()/h_sim_p->GetMaximum(); 
  double mc_p_scale_factor = h_data_p->GetMaximum()/h_mc_p->GetMaximum();
  //h_data_p->Scale(1.0/p_scale_factor);
  h_sim_p->Scale(p_scale_factor);
  //h_mc_p->Scale(mc_p_scale_factor);
  
  h_data_p->SetLineColor(kRed);
  h_sim_p->SetLineColor(kBlue-2);
  h_mc_p->SetLineColor(kViolet-4);

  //h_mc_p->Draw("h");  
  //h_mc_p->GetXaxis()->SetRangeUser(1.5,2.5);
  //h_mc_p->Draw("h");  

  h_data_p->SetTitle("Data vs Sim. Electron Momentum");
  h_data_p->Draw("HIST SAME");
  h_data_p->GetXaxis()->SetRangeUser(1.5,2.5);
  h_data_p->GetXaxis()->SetTitle("Momentum [GeV]");
  h_data_p->GetXaxis()->CenterTitle();  
  h_data_p->Draw("HIST SAME");

  h_sim_p->Draw("HIST SAME");
  h_sim_p->GetXaxis()->SetRangeUser(1.5,2.5);
  h_sim_p->Draw("HIST SAME");

  TLegend *l4 = new TLegend(0.11,0.7,0.48,0.89);
  l4->SetBorderSize(0);
  l4->AddEntry(h_data_p,"DATA");
  l4->AddEntry(h_sim_p,"SIM");
  //l4->AddEntry(h_mc_p,"GEN");
  l4->Draw();

  c2->Update();  

  kin_name="data_sim_all_mntm";
  if( printAll ){
    c2->Print(Form("data_sim_gen_comparison_r%d_f%s.pdf",run,field_config),"pdf"); //"h1.pdf(","pdf");
  }
  else{  
    c2->Print(Form("comparison_%s_r%d_f%s.pdf",kin_name.c_str(),run,field_config),"pdf"); //"h1.pdf(","pdf");
  }

  c2->Clear(); 

  //TCanvas *c5 = new TCanvas("c5","c5",900,900);
  c2->Divide(1,1);
  c2->cd(1);
  TH1D *h_data_theta = (TH1D*)fData->Get("particle_histograms_selected/hist_electron_theta");
  TH1D *h_sim_theta = (TH1D*)fSim->Get("particle_histograms_selected/hist_electron_theta");
  TH1D *h_mc_theta = (TH1D*)fMC->Get("hist_mc_all_electron_theta");
  //c5->cd();
  gStyle->SetOptStat(0000);

  double theta_scale_factor = h_data_theta->GetMaximum()/h_sim_theta->GetMaximum(); 
  //double mc_theta_scale_factor = h_data_theta->GetMaximum()/h_mc_theta->GetMaximum();
  //h_data_theta->Scale(1.0/theta_scale_factor);
  h_sim_theta->Scale(theta_scale_factor);
  //h_mc_theta->Scale(mc_theta_scale_factor);
  
  h_data_theta->SetLineColor(kRed);
  h_sim_theta->SetLineColor(kBlue-2);

  //h_mc_theta->Draw("h+same");  
  h_data_theta->SetTitle("Data vs Sim. Electron #theta");
  h_data_theta->Draw("HIST SAME");
  h_data_theta->GetXaxis()->SetTitle("#theta [deg]");
  h_data_theta->GetXaxis()->CenterTitle();
  h_sim_theta->Draw("HIST SAME");
  
  TLegend *l5 = new TLegend(0.7,0.7,0.89,0.89);
  l5->SetBorderSize(0);
  l5->AddEntry(h_data_theta,"DATA");
  l5->AddEntry(h_sim_theta,"SIM");
  //l5->AddEntry(h_mc_theta,"GEN");
  l5->Draw();
  
  //draw generated vs reconstructed simulated events
  c2->Update();
  kin_name="data_sim_all_theta";
  if( printAll ){
    c2->Print(Form("data_sim_gen_comparison_r%d_f%s.pdf",run,field_config),"pdf"); //"h1.pdf(","pdf");
  }
  else{  
    c2->Print(Form("comparison_%s_r%d_f%s.pdf",kin_name.c_str(),run,field_config),"pdf"); //"h1.pdf(","pdf");
  }

  c2->Clear(); 

  c2->Divide(1,1);
  c2->cd(1);
  
  gStyle->SetOptStat(0000);
  
  //double mc_theta_scale_factor = h_sim_theta->GetMaximum()/h_mc_theta->GetMaximum();
  //h_sim_theta->Scale(mc_theta_scale_factor);
  h_mc_theta->SetLineColor(kViolet-4);

  //h_mc_theta->GetYaxis()->SetRangeUser(0.0,h_sim_theta->GetMaximum());
  //h_mc_theta->Draw();
  gPad->SetLogy();
  h_mc_theta->SetTitle("Simulated vs Generated Electron Scattering Angle #theta");
  h_mc_theta->Draw("HIST SAME");
  h_mc_theta->GetXaxis()->SetTitle("#theta [deg]");
  h_mc_theta->GetXaxis()->CenterTitle();
  h_sim_theta->Scale(1.0/theta_scale_factor);
  h_sim_theta->Draw("HIST SAME");

  TLegend *l6 = new TLegend(0.7,0.7,0.89,0.89);
  l6->SetBorderSize(0);
  l6->AddEntry(h_mc_theta,"GEN");
  l6->AddEntry(h_sim_theta,"SIM");
  l6->Draw();  

  c2->Update();
  kin_name="gen_sim_theta";
  if( printAll ){
    c2->Print(Form("data_sim_gen_comparison_r%d_f%s.pdf",run,field_config),"pdf"); //"h1.pdf(","pdf");
  }
  else{  
    c2->Print(Form("comparison_%s_r%d_f%s.pdf",kin_name.c_str(),run,field_config),"pdf"); //"h1.pdf(","pdf");
  }
  c2->Clear(); 

  //draw generated momentum vs reconstructed simulated momentum
  c2->Divide(1,1);
  c2->cd(1);
  
  gStyle->SetOptStat(0000);
  
  //double mc_theta_scale_factor = h_sim_theta->GetMaximum()/h_mc_theta->GetMaximum();
  //h_sim_theta->Scale(mc_theta_scale_factor);
  h_mc_p->SetLineColor(kViolet-4);

  //h_mc_theta->GetYaxis()->SetRangeUser(0.0,h_sim_theta->GetMaximum());
  //h_mc_theta->Draw();
  gPad->SetLogy();
  h_mc_p->SetTitle("Simulated vs Generated Electron Momentum");
  h_mc_p->Draw("HIST SAME");  
  h_mc_p->GetXaxis()->SetRangeUser(1.5,2.5);
  h_mc_p->GetXaxis()->SetTitle("Momentum [GeV]");
  h_mc_p->GetXaxis()->CenterTitle();
  h_mc_p->Draw("HIST SAME");  

  h_sim_p->Scale(1.0/p_scale_factor);
  h_sim_p->Draw("HIST SAME");
  h_sim_p->GetXaxis()->SetRangeUser(1.5,2.5);
  h_sim_p->SetLineStyle(1);
  h_sim_p->Draw("HIST SAME");


  TLegend *l7 = new TLegend(0.7,0.7,0.89,0.89);
  l7->SetBorderSize(0);
  l7->AddEntry(h_mc_p,"GEN");
  l7->AddEntry(h_sim_p,"SIM");
  l7->Draw();
  
  c2->Update();
  kin_name="gen_sim_p";
  if( printAll ){
    c2->Print(Form("data_sim_gen_comparison_r%d_f%s.pdf",run,field_config),"pdf"); //"h1.pdf(","pdf");
  }
  else{  
    c2->Print(Form("comparison_%s_r%d_f%s.pdf",kin_name.c_str(),run,field_config),"pdf"); //"h1.pdf(","pdf");
  }


  c2->Clear(); 
  c2->SetCanvasSize(300,900);
  c2->SetWindowSize(340,910);
  //draw generated momentum vs reconstructed simulated and data theta
  c2->Divide(1,3);
  gStyle->SetOptStat(0000);
  c2->cd(1);
  h_mc_theta->SetTitle("Generated Electron Scattering Angle #theta");
  h_mc_theta->GetXaxis()->SetTitle("#theta [deg]");
  h_mc_theta->GetXaxis()->CenterTitle();
  h_mc_theta->Draw("HIST");
  
  c2->cd(2);
  h_sim_theta->SetTitle("Simulated Electron Scattering Angle #theta");
  h_sim_theta->GetXaxis()->SetTitle("#theta [deg]");
  h_sim_theta->GetXaxis()->CenterTitle();
  h_sim_theta->Draw("HIST");

  c2->cd(3);
  h_data_theta->SetTitle("Data Electron #theta");
  h_data_theta->Draw("HIST");
  h_data_theta->GetXaxis()->SetTitle("#theta [deg]");
  h_data_theta->GetXaxis()->CenterTitle();

  c2->Update();
  kin_name="gen_sim_data_theta3";
  if( printAll ){
    c2->Print(Form("data_sim_gen_comparison_r%d_f%s.pdf)",run,field_config),"pdf"); //"h1.pdf(","pdf");
  }
  else{  
    c2->Print(Form("comparison_%s_r%d_f%s.pdf",kin_name.c_str(),run,field_config),"pdf"); //"h1.pdf(","pdf");
  }


  //draw theta distribution in max bins
  TCanvas *c_data_sim_comp_theta = new TCanvas("c_data_sim_comp_theta","c_data_sim_comp_theta",800,1200);
  c_data_sim_comp_theta->Divide(2,3);
  gStyle->SetOptStat(0000);
  std::vector<int> max_phi_bins = { 39, 51, 63, 2, 14, 27 };
  TH2F *h_rc_all_el_theta_vs_phi = (TH2F*)fData->Get("particle_histograms_selected/hist_electron_theta_vs_phi");
  TH2F *h_sim_all_el_theta_vs_phi = (TH2F*)fSim->Get("particle_histograms_selected/hist_electron_theta_vs_phi");

  for( int s = 0; s < 6; s++ ){

    c_data_sim_comp_theta->cd(s+1);
    gPad->SetLogy();
    int bin_to_check=max_phi_bins[s];

    TH1D *h_data = h_rc_all_el_theta_vs_phi->ProjectionY(Form("data_proj_%s_s%d",h_rc_all_el_theta_vs_phi->GetTitle(),s),bin_to_check,bin_to_check);
    TH1D *h_sim = h_sim_all_el_theta_vs_phi->ProjectionY(Form("sim_proj_%s_s%d",h_sim_all_el_theta_vs_phi->GetTitle(),s),bin_to_check,bin_to_check);

    h_data->SetTitle(Form("Data #phi bin %d",bin_to_check));
    h_data->GetXaxis()->SetTitle("#theta [deg]");
    h_data->GetXaxis()->CenterTitle();

    h_data->SetLineColor(kRed);
    h_sim->SetLineColor(kBlue);

    h_data->Draw();
    //h_sim->Draw("same");

    TLegend *l1 = new TLegend(0.7, 0.7, 0.89, 0.89);
    l1->SetBorderSize(0);
    l1->AddEntry(h_data,"DATA");
    //l1->AddEntry(h_sim,"SIM");
    l1->Draw();
    
  }
  c_data_sim_comp_theta->SaveAs("data_theta_max_phi_bin.pdf");
 

  // Draw the simulated reconstructed vs data theta vs p and theta vs phi distributions  
  c2->Clear(); 
  c2->SetCanvasSize(600,900);
  c2->SetWindowSize(610,910);
  //draw generated momentum vs reconstructed simulated and data theta
  c2->Divide(2,3);
  gStyle->SetOptStat(0000);
  for( int s = 0; s < 6; s++ ){
    c2->cd(s+1);
    TH2D *h_ptheta_temp = (TH2D*)fData->Get(Form("kinematics/h_el_ptheta_s%d",s+1));
    h_ptheta_temp->SetTitle(Form("Data Mntm vs #theta Sector %d",s+1));
    h_ptheta_temp->GetYaxis()->SetTitle("#theta [deg]");
    h_ptheta_temp->GetXaxis()->SetTitle("momentum [GeV]");
    h_ptheta_temp->GetXaxis()->CenterTitle();
    h_ptheta_temp->GetYaxis()->CenterTitle();
    h_ptheta_temp->Draw("colz");     
  }

  c2->Update();
  kin_name="data_ptheta";
  if( printAll ){
    c2->Print(Form("data_sim_gen_comparison_r%d_f%s.pdf)",run,field_config),"pdf"); //"h1.pdf(","pdf");
  }
  else{  
    c2->Print(Form("comparison_%s_r%d_f%s.pdf",kin_name.c_str(),run,field_config),"pdf"); //"h1.pdf(","pdf");
  }

  // Draw SIMULATED P VS THETA
  c2->Clear(); 
  c2->SetCanvasSize(300,900);
  c2->SetWindowSize(340,910);

  c2->Divide(2,3);
  gStyle->SetOptStat(0000);
  for( int s = 0; s < 6; s++ ){
    c2->cd(s+1);
    TH2D *h_ptheta_temp = (TH2D*)fSim->Get(Form("kinematics/h_el_ptheta_s%d",s+1));
    h_ptheta_temp->SetTitle(Form("Sim. Mntm vs #theta Sector %d",s+1));
    h_ptheta_temp->GetYaxis()->SetTitle("#theta [deg]");
    h_ptheta_temp->GetXaxis()->SetTitle("momentum [GeV]");
    h_ptheta_temp->GetXaxis()->CenterTitle();
    h_ptheta_temp->GetYaxis()->CenterTitle();
    h_ptheta_temp->Draw("colz");     
  }

  c2->Update();
  kin_name="sim_ptheta";
  if( printAll ){
    c2->Print(Form("data_sim_gen_comparison_r%d_f%s.pdf)",run,field_config),"pdf"); //"h1.pdf(","pdf");
  }
  else{  
    c2->Print(Form("comparison_%s_r%d_f%s.pdf",kin_name.c_str(),run,field_config),"pdf"); //"h1.pdf(","pdf");
  }

  // Draw DATA PHI VS THETA
  c2->Clear(); 
  c2->SetCanvasSize(900,900);
  c2->SetWindowSize(910,910);
  c2->Divide(1,1);
  gStyle->SetOptStat(0000);
  h_rc_all_el_theta_vs_phi->SetTitle("Data Elastic Events #theta vs #phi");
  h_rc_all_el_theta_vs_phi->GetXaxis()->SetTitle("#phi [deg]");
  h_rc_all_el_theta_vs_phi->GetYaxis()->SetTitle("#theta [deg]");
  h_rc_all_el_theta_vs_phi->GetXaxis()->CenterTitle();
  h_rc_all_el_theta_vs_phi->GetYaxis()->CenterTitle();
  h_rc_all_el_theta_vs_phi->Draw("colz");

  c2->Update();
  kin_name="data_thetaphi";
  if( printAll ){
    c2->Print(Form("data_sim_gen_comparison_r%d_f%s.pdf)",run,field_config),"pdf"); //"h1.pdf(","pdf");
  }
  else{  
    c2->Print(Form("comparison_%s_r%d_f%s.pdf",kin_name.c_str(),run,field_config),"pdf"); //"h1.pdf(","pdf");
  }

  // Draw SIMULATION PHI VS THETA
  c2->Clear(); 
  c2->SetCanvasSize(900,900);
  c2->SetWindowSize(910,910);
  c2->Divide(1,1);
  gStyle->SetOptStat(0000);
  h_sim_all_el_theta_vs_phi->SetTitle("Simulated Elastic Events #theta vs #phi");
  h_sim_all_el_theta_vs_phi->GetXaxis()->SetTitle("#phi [deg]");
  h_sim_all_el_theta_vs_phi->GetYaxis()->SetTitle("#theta [deg]");
  h_sim_all_el_theta_vs_phi->GetXaxis()->CenterTitle();
  h_sim_all_el_theta_vs_phi->GetYaxis()->CenterTitle();
  h_sim_all_el_theta_vs_phi->Draw("colz");
  
  c2->Update();
  kin_name="sim_thetaphi";
  if( printAll ){
    c2->Print(Form("data_sim_gen_comparison_r%d_f%s.pdf)",run,field_config),"pdf"); //"h1.pdf(","pdf");
  }
  else{  
    c2->Print(Form("comparison_%s_r%d_f%s.pdf",kin_name.c_str(),run,field_config),"pdf"); //"h1.pdf(","pdf");
  }
   
  TCanvas *c3 = new TCanvas("c3","c3",900,900);
  c3->Divide(1,1);
  c3->cd(1);

  TH1D *h_el_vz = (TH1D*)fData->Get(Form("dc_detector/h_el_vz")); 
  TH1D *h_el_vz_sim = (TH1D*)fSim->Get(Form("dc_detector/h_el_vz"));

  double scale_vz_factor = h_el_vz->GetMaximum()/h_el_vz_sim->GetMaximum();
  h_el_vz->Scale(scale_vz_factor);

  gStyle->SetOptStat(0000);                                                                                                                                                                              
  h_el_vz->SetLineColor(kRed);
  h_el_vz_sim->SetLineColor(kBlue);
  h_el_vz->GetXaxis()->SetTitle("vz (cm)");
  h_el_vz->GetYaxis()->SetTitle("counts");
  

  h_el_vz->SetTitle("Selected Elastic Vz");
  h_el_vz_sim->Draw("HIST");
  h_el_vz->Draw("HIST SAME");

  TLegend *l8 = new TLegend(0.11,0.7,0.48,0.89);
  l8->SetBorderSize(0);
  l8->AddEntry(h_el_vz,"DATA");
  l8->AddEntry(h_el_vz_sim,"SIM");
  l8->Draw();

  c3->Update();

  kin_name="data_sim_gen_vz";
  if( printAll ){
    c3->Print(Form("data_sim_gen_comparison_r%d_f%s.pdf",run,field_config),"pdf"); //"h1.pdf(","pdf");
  }
  else{  
    c3->Print(Form("comparison_%s_r%d_f%s.pdf",kin_name.c_str(),run,field_config),"pdf"); //"h1.pdf(","pdf");
  }

  c3->Clear(); 
  c3->SetCanvasSize(900,900);
  c3->SetWindowSize(910,910);
  c3->Divide(1,1);
  gStyle->SetOptStat(0000);
 
  TH2D *h_el_vz_phi_sim = (TH2D*)fSim->Get("dc_detector/h_el_vz_phi");
  
  h_el_vz_phi_sim->SetTitle("Sim Vz vs #phi");
  h_el_vz_phi_sim->GetXaxis()->SetTitle("#phi (deg)");
  h_el_vz_phi_sim->GetYaxis()->SetTitle("vz (cm)");
  h_el_vz_phi_sim->Draw("colz");

  c3->Update();
  kin_name="sim_vz_phi";
  if( printAll ){
    c3->Print(Form("data_sim_gen_comparison_r%d_f%s.pdf)",run,field_config),"pdf"); //"h1.pdf(","pdf");
  }
  else{  
    c3->Print(Form("comparison_%s_r%d_f%s.pdf",kin_name.c_str(),run,field_config),"pdf"); //"h1.pdf(","pdf");
  }

  c3->Clear(); 
  c3->SetCanvasSize(900,900);
  c3->SetWindowSize(910,910);
  c3->Divide(1,1);
  gStyle->SetOptStat(0000);
 
  TH2D *h_el_vz_phi = (TH2D*)fData->Get("dc_detector/h_el_vz_phi");
  
  h_el_vz_phi->SetTitle("Data Vz vs #phi");
  h_el_vz_phi->GetXaxis()->SetTitle("#phi (deg)");
  h_el_vz_phi->GetYaxis()->SetTitle("vz (cm)");
  h_el_vz_phi->Draw("colz");

  c3->Update();
  kin_name="data_vz_phi";
  if( printAll ){
    c3->Print(Form("data_sim_gen_comparison_r%d_f%s.pdf)",run,field_config),"pdf"); //"h1.pdf(","pdf");
  }
  else{  
    c3->Print(Form("comparison_%s_r%d_f%s.pdf",kin_name.c_str(),run,field_config),"pdf"); //"h1.pdf(","pdf");
  }



  //DRAW VERTEX FOR EACH SECTOR
  c3->Clear(); 
  c3->SetCanvasSize(900,1200);
  c3->SetWindowSize(910,1200);
  c3->Divide(2,3);
  gStyle->SetOptStat(0000);

  for( int ss = 0; ss <= 5; ss++ ){
    c3->cd(ss+1);
    TH1D *h_el_vz_sect = (TH1D*)fData->Get(Form("dc_detector/h_el_vz_s%d",ss));
    TH1D *h_el_vz_sect_sim = (TH1D*)fSim->Get(Form("dc_detector/h_el_vz_s%d",ss));

    double scale_vz_factor_sect = h_el_vz_sect->GetMaximum()/h_el_vz_sect_sim->GetMaximum();
    h_el_vz_sect->Scale(1.0/scale_vz_factor);

    h_el_vz_sect->SetLineColor(kRed);
    h_el_vz_sect_sim->SetLineColor(kBlue-2);
    
    h_el_vz_sect->GetXaxis()->SetTitle("vz (cm)");
    h_el_vz_sect->GetYaxis()->SetTitle("Counts");

    h_el_vz_sect_sim->SetTitle(Form("Selected Electron Vz Sector %d",ss+1));
    h_el_vz_sect_sim->Draw("HIST");
    h_el_vz_sect->Draw("HIST SAME");

  }

  c3->Update();
  kin_name="data_vz_sector";
  if( printAll ){
    c3->Print(Form("data_sim_gen_comparison_r%d_f%s.pdf)",run,field_config),"pdf"); //"h1.pdf(","pdf");
  }
  else{  
    c3->Print(Form("comparison_%s_r%d_f%s.pdf",kin_name.c_str(),run,field_config),"pdf"); //"h1.pdf(","pdf");
  }
    

  //

  return 0;
}
