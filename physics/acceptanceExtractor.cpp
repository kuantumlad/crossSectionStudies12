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
double Ebeam = 2.211;//2193;

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
        
    if(  h_rec->GetBinCenter(b) > 25.0 ) continue;

    double accp_ratio = h_rec->GetBinContent(b)/h_gen->GetBinContent(b);    
    //    double accp_ratio_err = accp_ratio * sqrt( h_rec->GetBinContent(b)/pow(h_gen->GetBinContent(b),2)  + pow(h_rec->GetBinContent(b),2)/pow(h_gen->GetBinContent(b),3));
    
    double accp_ratio_err =1.0/h_gen->GetBinContent(b)  * sqrt( h_rec->GetBinContent(b)*h_rec->GetBinContent(b) / h_gen->GetBinContent(b) + h_rec->GetBinContent(b));
    if(h_gen->GetBinContent(b) == 0 ) accp_ratio = 0.0;
    if(h_gen->GetBinContent(b) == 0 || h_rec->GetBinContent(b) == 0 ) accp_ratio_err = 0.0; 
    
    std::cout << " bin " << b << " bin center " << h_rec->GetBinCenter(b) << " rec " << h_rec->GetBinContent(b) << " gen " << h_gen->GetBinContent(b) << std::endl;
    std::cout << " bin " << b << " bin error " <<  accp_ratio_err << std::endl;
    h_accp->SetBinContent(b,accp_ratio);
    h_accp->SetBinError(b,accp_ratio_err);
			  
  }

  return h_accp;

}

TGraphErrors* AcceptancePerSector(TH2F *h_rec_in, TH2F *h_gen_in){

  TH1F *h_rec = (TH1F*)h_rec_in->ProjectionY();
  
  int n_bins_filled = -1;
  for( int bb = 0; bb < h_rec_in->GetNbinsX(); bb++ ){    
    int current_content = h_rec_in->GetBinContent(bb,20); // scan over x bins at y_bin 20
    if( current_content > 0 ){
      n_bins_filled++;
    }
  }

  std::cout<< " projecting bins from 0 to " << n_bins_filled << std::endl;
  // get same number of phi bins as that in the reconstructed histogram
  TH1F *h_gen = (TH1F*)h_gen_in->ProjectionY(h_gen_in->GetTitle(), 0, n_bins_filled);
  

  // calculate the acceptance integrated over a sector
  std::vector<double> v_theta_center;
  std::vector<double> v_accp;
  std::vector<double> v_theta_center_err;
  std::vector<double> v_accp_err;

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
        
    double theta_center = h_rec->GetBinCenter(b);
    if( theta_center < 9.0 || theta_center > 25.0 ) continue; 
    double accp_ratio = h_rec->GetBinContent(b)/h_gen->GetBinContent(b);    
    //double accp_ratio_err = accp_ratio * sqrt( h_rec->GetBinContent(b)/pow(h_gen->GetBinContent(b),2)  + pow(h_rec->GetBinContent(b),2)/pow(h_gen->GetBinContent(b),3));
    double accp_ratio_err =1.0/h_gen->GetBinContent(b)  * sqrt( h_rec->GetBinContent(b)*h_rec->GetBinContent(b) / h_gen->GetBinContent(b) + h_rec->GetBinContent(b));
    if(h_gen->GetBinContent(b) == 0 ) accp_ratio = 0.0;
    if(h_gen->GetBinContent(b) == 0 || h_rec->GetBinContent(b) == 0 ) accp_ratio_err = 0.0; 
    
    std::cout << " bin " << b << " bin center " << h_rec->GetBinCenter(b) << " rec " << h_rec->GetBinContent(b) << " gen " << h_gen->GetBinContent(b) << std::endl;
    std::cout << " bin " << b << " bin error " <<  accp_ratio_err << std::endl;
    h_accp->SetBinContent(b,accp_ratio);
    h_accp->SetBinError(b,accp_ratio_err);

    v_theta_center.push_back(theta_center);
    v_accp.push_back(accp_ratio);
    v_accp_err.push_back(accp_ratio_err);
    v_theta_center_err.push_back(0.0);
			  
  }

  TGraphErrors *graph = new TGraphErrors(v_theta_center.size(), &(v_theta_center[0]), &(v_accp[0]), &(v_theta_center_err[0]), &(v_accp_err[0]) );

  return graph;
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
  
  //TH2F *h_rec_clone = (TH2F*)h_rec->Clone();
  //h_rec_clone->Divide(h_gen);

  for( int bx = 1; bx < h_rec->GetNbinsX(); bx++ ){
    for( int by = 1; by < h_rec->GetNbinsY(); by++ ){
      
      double accp_ratio = h_rec->GetBinContent(bx,by)/h_gen->GetBinContent(bx,by);
      double accp_ratio_err =1.0/h_gen->GetBinContent(bx,by)  * sqrt( h_rec->GetBinContent(bx,by)*h_rec->GetBinContent(bx,by) / h_gen->GetBinContent(bx,by) + h_rec->GetBinContent(bx,by));  // accp_ratio * sqrt( h_rec->GetBinContent(bx,by)/pow(h_gen->GetBinContent(bx,by),2)  + pow(h_rec->GetBinContent(bx,by),2)/pow(h_gen->GetBinContent(bx,by),3));

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


int acceptanceExtractor(const char* inFileData, int run, const char* field_config ){ // const char* outfile, int run ){

 
  TFile *fData; 
  TFile *fMC; 
  Char_t tmpstr[80];
  Double_t fraction;

  fData = new TFile(inFileData,"");   // Input File
  if(fData->IsZombie() ){
    cout<<"Input file doesn't exist!" << endl;
    cout<<"Exit program" << endl;
    return 0;
  }

  cout << "Reading from File: " << inFileData << endl;

  // get relevant histograms to determine efficiency or geometrical acceptance
  TH1F *h_rc_all_el_p = (TH1F*)fData->Get("particle_histograms_selected/hist_electron_p");
  TH1F *h_mc_all_el_p = (TH1F*)fData->Get("mc/hist_mc_all_electron_p");

  TH1F *h_rc_all_el_theta = (TH1F*)fData->Get("particle_histograms_selected/hist_electron_theta");
  TH1F *h_mc_all_el_theta = (TH1F*)fData->Get("mc/hist_mc_all_electron_theta");

  TH1F *h_rc_all_el_phi = (TH1F*)fData->Get("particle_histograms_selected/hist_electron_phi");
  TH1F *h_mc_all_el_phi = (TH1F*)fData->Get("mc/hist_mc_all_electron_phi");

  TH2F *h_rc_all_el_p_theta = (TH2F*)fData->Get("particle_histograms_selected/hist_electron_p_theta");
  TH2F *h_mc_all_el_p_theta = (TH2F*)fData->Get("mc/hist_mc_all_electron_p_theta");

  TH2F *h_rc_all_el_theta_vs_phi = (TH2F*)fData->Get("particle_histograms_selected/hist_electron_theta_vs_phi");
  TH2F *h_mc_all_el_theta_vs_phi = (TH2F*)fData->Get("mc/hist_mc_all_electron_theta_vs_phi");

  // change theta range starting from 5 
  TH1F *h_el_theta_rebin = new TH1F("h_el_theta_rebin","h_el_theta_rebin", h_rc_all_el_theta->GetNbinsX(), 0.0, 30.0);
  for( int b = 0; b< h_rc_all_el_theta->GetNbinsX(); b++ ){
    double bin_content = h_rc_all_el_theta->GetBinContent(b);
    if( h_rc_all_el_theta->GetBinCenter(b) <= 6.0 ) bin_content = 0.0;
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

  TGraphErrors *g_accp_sector;

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



  TCanvas *c_g_accp_s = new TCanvas("c_g_accp_s","c_g_accp_s",1200,800);
  c_g_accp_s->Divide(2,3);  
  for ( int s = 0; s < 6; s++ ){
    c_g_accp_s->cd(s+1);
    TH2F *h_rec_temp  = (TH2F*)fData->Get(Form("kinematics/h_el_phitheta_s%d_final",s+1));
    TH2F *h_gen_temp = (TH2F*)fData->Get("mc/hist_mc_all_electron_theta_vs_phi");
    std::cout << " getting acceptance for sector " << s << std::endl;
    TGraphErrors *g_accp_sector = AcceptancePerSector(h_rec_temp, h_gen_temp);   
    g_accp_sector->SetTitle(Form(" Acceptance Sector %d",s+1));
    g_accp_sector->Draw("AP");
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // acceptance per phi bin
  // using the historgrams from the selector_elastic 
  TCanvas *c_pb = new TCanvas("c_pb","c_pb",2000,2000);
  c_pb->Divide(10,10);
  for( int pb = 0; pb < 73; pb++ ){
    c_pb->cd(pb+1);
    
    TH1F *h_rec_temp = (TH1F*)fData->Get(Form("acceptance/h_el_theta_rec_phi%d",pb));
    TH1F *h_gen_temp = (TH1F*)fData->Get(Form("acceptance/h_el_theta_gen_phi%d",pb));
    TH1F* h_accp_ratio_temp = AcceptanceRatio(h_rec_temp, h_gen_temp );
    h_accp_ratio_temp->Draw("ERR");
    
    std::vector< double > v_theta_bin;
    std::vector< double > v_accp;
    std::vector< double > v_theta_bin_err;
    std::vector< double > v_accp_err;
    
    std::cout << "creating accp across theta for phibin " << pb << std::endl;
    for( int tb = 0; tb < h_rec_temp->GetNbinsX(); tb++ ){
      
      double accp_val = h_accp_ratio_temp->GetBinContent(tb);
      double theta_bin_center = h_accp_ratio_temp->GetBinCenter(tb);
      double accp_val_err = h_accp_ratio_temp->GetBinError(tb);
      double theta_bin_center_err = 0.0;
      std::cout << "bin center " << theta_bin_center << std::endl;
      if( theta_bin_center < 8 || theta_bin_center > 21 ) continue; // there is an issue with tracking above 20deg 
      std::cout << " accp " << accp_val << " bin center " << theta_bin_center << " accp err " <<accp_val_err << " theta bin err " << theta_bin_center_err << std::endl;
      v_theta_bin.push_back(theta_bin_center);
      v_accp.push_back(accp_val);
      v_theta_bin_err.push_back(theta_bin_center_err);
      v_accp_err.push_back(accp_val_err);
    }

    std::cout <<" creating graph for accetpance now " << std::endl;
    TGraphErrors *g_accp_temp = new TGraphErrors(v_theta_bin.size(), &(v_theta_bin[0]), &(v_accp[0]), &(v_theta_bin_err[0]), &(v_accp_err[0]) );
    g_accp_temp->SetTitle(Form("Acceptance #phi Bin %d",pb));  
    g_accp_temp->GetXaxis()->SetTitle("#theta (deg)");
    g_accp_temp->GetYaxis()->SetTitle("Accp");
    g_accp_temp->GetXaxis()->CenterTitle();
    g_accp_temp->GetYaxis()->CenterTitle();
    g_accp_temp->SetMarkerStyle(10);
    g_accp_temp->SetMarkerSize(0.75);
    g_accp_temp->SetMarkerColor(kBlack);
    g_accp_temp->Draw("AP");       

    if( pb == 28 ){
      TCanvas *c_pb_focus = new TCanvas("c_pb_focus","c_pb_focus",900,900);
      c_pb_focus->cd(1);
      g_accp_temp->Draw("AP"); 
      c_pb_focus->SaveAs("g_acceptance_phi_bin_28.pdf");
    }
   
  }
  c_pb->SaveAs(Form("g_acceptance_per_phi_bin_%s.pdf",field_config));
  
  TMultiGraph *mg_accp = new TMultiGraph();

  TCanvas *c_accp_sector = new TCanvas("c_accp_sector","c_accp_sector",900,1200);
  c_accp_sector->Divide(2,3);
  for(int ss = 0; ss < 6; ss++){
    c_accp_sector->cd(ss+1);
    TH1F *h_rec_temp = (TH1F*)fData->Get(Form("acceptance/h_el_theta_rec_s%d",ss));
    TH1F *h_gen_temp = (TH1F*)fData->Get(Form("acceptance/h_el_theta_gen_s%d",0));
    
    TH1F* h_accp_ratio_temp = AcceptanceRatio(h_rec_temp, h_gen_temp );
    //    h_accp_ratio_temp->Draw("ERR");
    std::vector< double > v_theta_bin;
    std::vector< double > v_accp;
    std::vector< double > v_theta_bin_err;
    std::vector< double > v_accp_err;
    
    for( int tb = 0; tb < h_rec_temp->GetNbinsX(); tb++ ){
      double accp_val = h_accp_ratio_temp->GetBinContent(tb)*6.0; // times 6 to account for 6 sectors -> will need to do on a bin by bin basis when looking at entire sector
      double theta_bin_center = h_accp_ratio_temp->GetBinCenter(tb);
      double accp_val_err = h_accp_ratio_temp->GetBinError(tb);
      double theta_bin_center_err = 0.0;
      if( theta_bin_center < 8 || theta_bin_center > 21 ) continue; // there is an issue with tracking above 20deg 
      v_theta_bin.push_back(theta_bin_center);
      v_accp.push_back(accp_val);
      v_theta_bin_err.push_back(theta_bin_center_err);
      v_accp_err.push_back(accp_val_err);
    }
    
    TGraphErrors *g_accp_temp = new TGraphErrors(v_theta_bin.size(), &(v_theta_bin[0]), &(v_accp[0]), &(v_theta_bin_err[0]), &(v_accp_err[0]) );
    g_accp_temp->SetTitle(Form("Acceptance Sector %d",ss+1));  
    g_accp_temp->GetXaxis()->SetTitle("#theta (deg)");
    g_accp_temp->GetYaxis()->SetTitle("Accp");
    g_accp_temp->GetXaxis()->CenterTitle();
    g_accp_temp->GetYaxis()->CenterTitle();
    g_accp_temp->SetMarkerStyle(10);
    g_accp_temp->SetMarkerSize(0.75);
    g_accp_temp->SetMarkerColor(kBlack);
    g_accp_temp->Draw("AP");
  
    g_accp_temp->SetMarkerColor(kSpring+ss);
    //g_accp_temp->SetMarkerSize(1.2);
    mg_accp->Add(g_accp_temp);

  }
  c_accp_sector->SaveAs("c_accp_sector.pdf");

  TCanvas *c_mg_accp = new TCanvas("c_mg_accp","c_mg_accp",900,900);  
  mg_accp->Draw("AP");
  mg_accp->GetXaxis()->SetTitle("#theta (deg)");
  mg_accp->GetYaxis()->SetTitle("Accp");
  mg_accp->GetXaxis()->CenterTitle();
  mg_accp->GetYaxis()->CenterTitle();
  c_mg_accp->SaveAs(Form("c_mg_accp_per_sector_%s.pdf",field_config)); 
  //mg_accp->GetXaxis()->SetLimits(6.0,21.0);
  //mg_accp->Draw("AP");                                                                                                                                                                                   
  //c_mg_accp->Update();
  
 
  return 0;  



}


