#include <vector>
#include <string>
#include <map>
#include <iostream>
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"

const static int BUFFER = 50;
double Ebeam = 2.22193;

int ele_sector[BUFFER];
TLorentzVector ele[BUFFER]; 
TLorentzVector prot[BUFFER];
TLorentzVector mc_ele[BUFFER];
TLorentzVector mc_prot[BUFFER];

int select_ele;

TH1F* AcceptanceRatio(TH1F*,TH1F*);
TH2F* AcceptanceRatio2D(TH2F*,TH2F*);
TH2F* WriteAcceptanceRatio2D(TH2F*,TH2F*);

TH1F* AcceptanceRatio(TH1F *h_rec, TH1F *h_gen ){


  int rec_bins = h_rec->GetNbinsX();
  int mc_bins = h_gen->GetNbinsX();
  if( mc_bins != rec_bins ) std::cout << " ERROR WITH BIN NUMBERS CAN NO PROCEED " << std::endl;

  std::cout << " rec bin " << rec_bins << " mc bins " << mc_bins << std::endl;
  
  std::cout << " min bin center " << h_rec->GetBinCenter(1) << std::endl;
  std::cout << " >> bin width/2 " << (h_rec->GetBinCenter(1) + h_rec->GetBinCenter(2))/2.0  << std::endl;
  std::cout << " root bin width " << h_rec->GetBinWidth(1) << std::endl;
  double bin_width =  h_rec->GetBinWidth(1);
  double max_range = h_rec->GetBinCenter(h_rec->GetNbinsX()) + bin_width/2.0;
  double min_range = h_rec->GetBinCenter(1) - bin_width/2.0;

  //create acceptance histrogram over interval min and max range with rec bins (assuming same as mc bins)
  std::cout << " creating acceptance histogram for " << h_rec->GetTitle() << std::endl;
  std::cout << " bin range from " << min_range << " to " << max_range << std::endl;
  TH1F *h_accp = new TH1F(Form("h_acc_%s",h_rec->GetTitle()),Form("h_acc_%s",h_rec->GetTitle()), rec_bins, min_range, max_range);


  for( int b = 1; b < rec_bins; b++ ){
    
    
    
    double accp_ratio = h_rec->GetBinContent(b)/h_gen->GetBinContent(b);    
    double accp_ratio_err = accp_ratio * sqrt( h_rec->GetBinContent(b)/pow(h_gen->GetBinContent(b),2)  + pow(h_rec->GetBinContent(b),2)/pow(h_gen->GetBinContent(b),3));
    if(h_gen->GetBinContent(b) == 0 ) accp_ratio = 0.0;
    if(h_gen->GetBinContent(b) == 0 || h_rec->GetBinContent(b) == 0 ) accp_ratio_err = 0.0; 
    

    std::cout << " bin " << b << " bin center " << h_rec->GetBinCenter(b) << " rec " << h_rec->GetBinContent(b) << " gen " << h_gen->GetBinContent(b) << std::endl;
    std::cout << " bin " << b << " bin error " <<  accp_ratio_err << std::endl;
    h_accp->SetBinContent(b,accp_ratio);
    h_accp->SetBinError(b,accp_ratio_err);
			  
  }

  return h_accp;

}

TH2F* AcceptanceRatio2D(TH2F* h_rec, TH2F* h_gen){


  if( h_rec->GetNbinsX() != h_gen->GetNbinsX() || h_rec->GetNbinsY() != h_gen->GetNbinsY() ){
    std::cout << " bins are not equal in REC and GEN histograms " << std::endl;
  }

  double max_rangeX = h_rec->ProjectionX()->GetBinCenter(h_rec->ProjectionX()->GetNbinsX()) + (h_rec->ProjectionX()->GetBinCenter(1) + h_rec->ProjectionX()->GetBinCenter(2))/2.0;
  double min_rangeX = h_rec->ProjectionX()->GetBinCenter(1) - (h_rec->ProjectionX()->GetBinCenter(1) + h_rec->ProjectionX()->GetBinCenter(2))/2.0;

  double max_rangeY = h_rec->ProjectionY()->GetBinCenter(h_rec->ProjectionY()->GetNbinsY()) + (h_rec->ProjectionY()->GetBinCenter(1) + h_rec->ProjectionY()->GetBinCenter(2))/2.0;
  double min_rangeY = h_rec->ProjectionY()->GetBinCenter(1) - (h_rec->ProjectionY()->GetBinCenter(1) + h_rec->ProjectionY()->GetBinCenter(2))/2.0;


  std::cout << " Rec Number of X bins " << h_rec->GetNbinsX() << " Number of Y bins " << h_rec->GetNbinsY() << std::endl;
  std::cout << " Gen Number of X bins " << h_gen->GetNbinsX() << " Number of Y bins " << h_gen->GetNbinsY() << std::endl;
  std::cout <<  " X range " << min_rangeX << " " << max_rangeX << std::endl;
  std::cout <<  " Y range " << min_rangeY << " " << max_rangeY << std::endl;

  TH2F *h2_accp = new TH2F(Form("h2_accp_%s",h_rec->GetTitle()),Form("h2_accp_%s",h_rec->GetTitle()), 73, -180.0, 180.0, 30, 0.0, 30.0 ); //h_rec->GetNbinsY(), min_rangeX, max_rangeX, min_rangeY, max_rangeY);
  
  for( int bx = 1; bx < h_rec->GetNbinsX(); bx++ ){
    for( int by = 1; by < h_rec->GetNbinsY(); by++ ){
      double accp_ratio = h_rec->GetBinContent(bx,by)/h_gen->GetBinContent(bx,by);
      if( h_gen->GetBinContent(bx,by) == 0 ) accp_ratio = 0.0;

      std::cout << " bx " << bx << " by " << by << " " <<  h_rec->GetBinContent(bx,by) << " " << h_gen->GetBinContent(bx,by) << " "  << accp_ratio << std::endl;
      h2_accp->SetBinContent(bx,by,accp_ratio);
    }
  }
  
			   
  return h2_accp;

}

//// write acceptance values to file 
void WriteAcceptanceRatio2D(TH2F* h_rec, TH2F* h_gen, int run, const char* field){


  if( h_rec->GetNbinsX() != h_gen->GetNbinsX() || h_rec->GetNbinsY() != h_gen->GetNbinsY() ){
    std::cout << " bins are not equal in REC and GEN histograms " << std::endl;
  }

  double max_rangeX = h_rec->ProjectionX()->GetBinCenter(h_rec->ProjectionX()->GetNbinsX()) + (h_rec->ProjectionX()->GetBinCenter(1) + h_rec->ProjectionX()->GetBinCenter(2))/2.0;
  double min_rangeX = h_rec->ProjectionX()->GetBinCenter(1) - (h_rec->ProjectionX()->GetBinCenter(1) + h_rec->ProjectionX()->GetBinCenter(2))/2.0;

  double max_rangeY = h_rec->ProjectionY()->GetBinCenter(h_rec->ProjectionY()->GetNbinsY()) + (h_rec->ProjectionY()->GetBinCenter(1) + h_rec->ProjectionY()->GetBinCenter(2))/2.0;
  double min_rangeY = h_rec->ProjectionY()->GetBinCenter(1) - (h_rec->ProjectionY()->GetBinCenter(1) + h_rec->ProjectionY()->GetBinCenter(2))/2.0;

  std::cout << " Rec Number of X bins " << h_rec->GetNbinsX() << " Number of Y bins " << h_rec->GetNbinsY() << std::endl;
  std::cout << " Gen Number of X bins " << h_gen->GetNbinsX() << " Number of Y bins " << h_gen->GetNbinsY() << std::endl;
  std::cout <<  " X range " << min_rangeX << " " << max_rangeX << std::endl;
  std::cout <<  " Y range " << min_rangeY << " " << max_rangeY << std::endl;

  TH2F *h2_accp = new TH2F(Form("h2_accp_%s",h_rec->GetTitle()),Form("h2_accp_%s",h_rec->GetTitle()), 73, -180.0, 180.0, 30, 0.0, 30.0 ); //h_rec->GetNbinsY(), min_rangeX, max_rangeX, min_rangeY, max_rangeY);

  ofstream outputAcceptance;
  ofstream outputAcceptanceError;
  std::string parentDirectory = "/w/hallb-scifs17exp/clas12/bclary/CLAS12/electron_studies/physics/parameters/";
  std::string f_out_name = parentDirectory+"elastic_theta_phi_acceptance_"+std::to_string(run)+"_f"+field+".txt";
  std::string f_out_name_err = parentDirectory+"elastic_theta_phi_acceptance_error_"+std::to_string(run)+"_f"+field+".txt";
  std::cout << " Creating Acceptance output file " << f_out_name_err << std::endl;
  outputAcceptance.open(f_out_name);
  outputAcceptanceError.open(f_out_name_err);
  
  for( int bx = 1; bx < h_rec->GetNbinsX(); bx++ ){
    for( int by = 1; by < h_rec->GetNbinsY(); by++ ){
      
      double accp_ratio = h_rec->GetBinContent(bx,by)/h_gen->GetBinContent(bx,by);
      double accp_ratio_err = accp_ratio * sqrt( h_rec->GetBinContent(bx,by)/pow(h_gen->GetBinContent(bx,by),2)  + pow(h_rec->GetBinContent(bx,by),2)/pow(h_gen->GetBinContent(bx,by),3));
      if(h_gen->GetBinContent(bx,by) == 0 || h_rec->GetBinContent(bx,by) == 0 ) accp_ratio_err = 0.0; 
      if( h_gen->GetBinContent(bx,by) == 0 ) accp_ratio = 0.0;
      
      std::cout <<  accp_ratio << " ";
      outputAcceptance << accp_ratio << " ";
      outputAcceptanceError << accp_ratio_err << " ";
    }
    std::cout << "\n";
    outputAcceptance << "\n";
    outputAcceptanceError << "\n";

  }  			   
}


int acceptanceExtractor(const char* inFileData, const char* inFileMC, int run, const char* field_config ){ // const char* outfile, int run ){

 
  TFile *fData; 
  TFile *fMC; 
  Char_t tmpstr[80];
  Double_t fraction;

  fData = new TFile(inFileData,"");   // Input File
  fMC = new TFile(inFileMC,"");
  if(fData->IsZombie() || fMC->IsZombie()){   // Check if TFile exists!
    cout<<"Input file doesn't exist!" << endl;
    cout<<"Exit program" << endl;
    return 0;
  }

  cout << "Reading from File: " << inFileData << " and " << inFileMC << endl;

  // get relevant histograms to determine efficiency or geometrical acceptance
  TH1F *h_rc_all_el_p = (TH1F*)fData->Get("particle_histograms_selected/hist_electron_p");
  TH1F *h_mc_all_el_p = (TH1F*)fMC->Get("hist_mc_all_electron_p");

  TH1F *h_rc_all_el_theta = (TH1F*)fData->Get("particle_histograms_selected/hist_electron_theta");
  TH1F *h_mc_all_el_theta = (TH1F*)fMC->Get("hist_mc_all_electron_theta");

  TH1F *h_rc_all_el_phi = (TH1F*)fData->Get("particle_histograms_selected/hist_electron_phi");
  TH1F *h_mc_all_el_phi = (TH1F*)fMC->Get("hist_mc_all_electron_phi");

  TH2F *h_rc_all_el_p_theta = (TH2F*)fData->Get("particle_histograms_selected/hist_electron_p_theta");
  TH2F *h_mc_all_el_p_theta = (TH2F*)fMC->Get("hist_mc_all_electron_p_theta");

  TH2F *h_rc_all_el_theta_vs_phi = (TH2F*)fData->Get("particle_histograms_selected/hist_electron_theta_vs_phi");
  TH2F *h_mc_all_el_theta_vs_phi = (TH2F*)fMC->Get("hist_mc_all_electron_theta_vs_phi");

  // change theta range starting from 5 
  TH1F *h_el_theta_rebin = new TH1F("h_el_theta_rebin","h_el_theta_rebin", h_rc_all_el_theta->GetNbinsX(), 0.0, 30.0);
  for( int b = 0; b< h_rc_all_el_theta->GetNbinsX(); b++ ){
    double bin_content = h_rc_all_el_theta->GetBinContent(b);
    if( h_rc_all_el_theta->GetBinCenter(b) < 6.0 ) bin_content = 0.0;
    std::cout << " > rebin value " << b << " " << bin_content << std::endl;
    h_el_theta_rebin->SetBinContent(b, bin_content);
  }


  TCanvas *c1 = new TCanvas("c1","c1",900,900);
  c1->Divide(1,1);
  c1->cd(1);
  //TCanvas *c_accP = new TCanvas("c_accP","c_accP",900,900);
  //c_accP->Divide(1,1);
  //c_accP->cd(1);
  
  TH1F *h_accp_p = AcceptanceRatio(h_rc_all_el_p, h_mc_all_el_p );
  h_accp_p->SetMarkerStyle(20);
  h_accp_p->SetMarkerColor(kBlack);
  h_accp_p->Draw();

  c1->Update();
  c1->Print(Form("h_acceptance_r%d_f%s.pdf(",run,field_config),"pdf"); //"h1.pdf(","pdf");
  c1->Clear(); 
  c1->Divide(2,1);

  //TCanvas *c_accTheta = new TCanvas("c_accTheta","c_accTheta",900,900);
  //c_accTheta->Divide(2,1);
  c1->cd(1);
  //c_accTheta->cd(1);
  h_el_theta_rebin->SetMarkerStyle(20);
  h_el_theta_rebin->SetMarkerColor(kRed);
  h_el_theta_rebin->Draw();
  h_mc_all_el_theta->SetMarkerStyle(20);
  h_mc_all_el_theta->SetMarkerColor(kBlue);
  h_mc_all_el_theta->Draw("same");

  //c_accTheta->cd(2);  
  c1->cd(2);
  TH1F *h_accp_theta = AcceptanceRatio(h_rc_all_el_theta, h_mc_all_el_theta);
  h_accp_theta->SetMarkerStyle(20);
  h_accp_theta->SetMarkerColor(kBlack);
  h_accp_theta->Draw("HIST");

  c1->Update();
  c1->Print(Form("h_acceptance_r%d_f%s.pdf",run,field_config),"pdf"); //"h1.pdf(","pdf");
  c1->Clear(); 
  c1->Divide(1,1);
  c1->cd(1);


  //TCanvas *c_accPhi = new TCanvas("c_accPhi","c_accPhi",900,900);
  //c_accPhi->Divide(1,1);
  //c_accPhi->cd(1);
  
  TH1F *h_accp_phi = AcceptanceRatio(h_rc_all_el_phi, h_mc_all_el_phi);
  h_accp_phi->SetMarkerStyle(20);
  h_accp_phi->SetMarkerColor(kBlack);
  h_accp_phi->Draw("HIST");

  c1->Update();
  c1->Print(Form("h_acceptance_r%d_f%s.pdf(",run,field_config),"pdf"); //"h1.pdf(","pdf");
  c1->Clear(); 
  c1->Divide(1,1);

  //TCanvas *c_accThetaPhi = new TCanvas("c_accThetaPhi","c_accThetaPhi",900,900);
  //c_accThetaPhi->Divide(1,1);
  //c_accThetaPhi->cd(1);

  c1->cd(1);
  
  TH2F *h_accp_theta_phi = AcceptanceRatio2D(h_rc_all_el_theta_vs_phi, h_mc_all_el_theta_vs_phi );
  gStyle->SetOptStat(0000);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.12);
  h_accp_theta_phi->SetMarkerStyle(20);
  h_accp_theta_phi->SetMarkerColor(kBlack);
  h_accp_theta_phi->GetXaxis()->SetTitle("#phi [deg]");
  h_accp_theta_phi->GetYaxis()->SetTitle("#theta [deg]");
  h_accp_theta_phi->Draw("colz");

  c1->Update();
  c1->Print(Form("h_acceptance_r%d_f%s.pdf)",run,field_config),"pdf"); //"h1.pdf(","pdf");

  //Get the theta acceptance in the corresponding phi bin used in the cross section (bin 39)
  TCanvas *c_accp_theta_proj = new TCanvas("c_accp_theta_proj ","c_accp_theta_proj ",910,910);
  c_accp_theta_proj->cd();
  gStyle->SetOptStat(0000);
  gPad->SetLeftMargin(0.12);
  TH1D *h_accep_theta_proj = h_accp_theta_phi->ProjectionY(Form("proj_%s",h_accp_theta_phi->GetTitle()),39,39);
  h_accep_theta_proj->SetTitle("#theta Acceptance");
  h_accep_theta_proj->GetXaxis()->SetTitle("#theta [deg]");
  h_accep_theta_proj->GetYaxis()->SetTitle("Acceptance");
  h_accep_theta_proj->GetXaxis()->CenterTitle();
  h_accep_theta_proj->GetYaxis()->CenterTitle();
  h_accep_theta_proj->Draw("HIST");
  c_accp_theta_proj->SaveAs("h_accep_theta_proj_phibin39.pdf");

  TCanvas *c_accp_err = new TCanvas("c_accp_err","c_accp_err",900,900);
  TH1F *h_accp_theta_err = AcceptanceRatio((TH1F*)h_rc_all_el_theta_vs_phi->ProjectionY("proj_rec",39,39),(TH1F*)h_mc_all_el_theta_vs_phi->ProjectionY("proj_gen",39,39));
  c_accp_err->cd();
  h_accp_theta_err->SetTitle("#theta Acceptance");
  h_accp_theta_err->GetXaxis()->SetTitle("#theta [deg]");
  h_accp_theta_err->GetYaxis()->SetTitle("Acceptance"); 
  h_accp_theta_err->GetXaxis()->CenterTitle(); 
  h_accp_theta_err->GetYaxis()->CenterTitle(); 
  h_accp_theta_err->Draw("HIST E");
  c_accp_err->Update();
  c_accp_err->SaveAs("h_accp_theta_proj_phibin39_err.pdf");

  //show the sim vs gen theta in each bin for each sector
  TCanvas *c_sim_gen_comp_theta = new TCanvas("c_sim_gen_comp_theta","c_sim_gen_comp_theta",620,920);
  c_sim_gen_comp_theta->Divide(2,3);
  gStyle->SetOptStat(0000);
  std::vector<int> max_phi_bins = { 39, 51, 63, 2, 14, 27 };

  for( int s = 0; s < 6; s++ ){

    c_sim_gen_comp_theta->cd(s+1);
    gPad->SetLogy();
    int bin_to_check=max_phi_bins[s];

    TH1D *h_sim = h_rc_all_el_theta_vs_phi->ProjectionY(Form("sim_proj_%s_s%d",h_rc_all_el_theta_vs_phi->GetTitle(),s),bin_to_check,bin_to_check);
    TH1D *h_gen = h_mc_all_el_theta_vs_phi->ProjectionY(Form("mc_proj_%s_s%d",h_mc_all_el_theta_vs_phi->GetTitle(),s),bin_to_check,bin_to_check);
    
    h_gen->SetTitle(Form("SIM vs GEN #phi bin %d",bin_to_check));
    h_gen->GetXaxis()->SetTitle("#theta [deg]");
    h_gen->GetXaxis()->CenterTitle();

    h_gen->SetLineColor(kViolet);
    h_sim->SetLineColor(kBlue);

    h_gen->Draw("HIST");
    h_sim->Draw("HIST SAME");

    TLegend *l1 = new TLegend(0.7,0.7,0.89,0.89);
    l1->SetBorderSize(0);
    l1->AddEntry(h_sim,"SIM");
    l1->AddEntry(h_gen,"GEN");
    l1->Draw();
    
  }
  c_sim_gen_comp_theta->SaveAs("sim_gen_comp_theta_max_phi_bin.pdf");

  //Write acceptance to file for use in cross section calculations
  WriteAcceptanceRatio2D(h_rc_all_el_theta_vs_phi, h_mc_all_el_theta_vs_phi, run, field_config);

  return 0;  



}


