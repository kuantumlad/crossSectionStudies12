#include <iostream>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TMath.h>
#include <TH2F.h>
#include <vector>
#include <map>


TH1D* getSectorCrossSection(TH1F *h_in, TH2F *h2_in){

  int phi_bins = h_in->GetNbinsX();
  std::cout << " Number of phi bins: " << phi_bins << std::endl;

  int temp_entries=-1;
  int temp_bin = -1;
  for( int b = 0; b <= phi_bins; b++ ){
    int b_entries = h_in->GetBinContent(b);
    std::cout << " bin entries " << b_entries << std::endl;

    if( b_entries > temp_entries ){
      temp_entries=b_entries;
      temp_bin = b;
    }   
  }

  std::cout << " Bin with max entries: " << temp_bin << " with entries: " << temp_entries << std::endl;
  int bin_min = temp_bin - 2;
  int bin_max = temp_bin + 2;

  std::cout << " Extending bins by 2 bins on each side; bin min: " << bin_min << " bin max: " << bin_max << std::endl;
  
  TH1D *h_theta_proj = h2_in->ProjectionY(Form("proj_%s",h_in->GetTitle()),bin_min, bin_max);
  h_theta_proj->Rebin(5);
 
  double lum_run2391 = 6091588.274652874; //609158.862236702; //6.0915886E6;
  double bin_phi_size = ( 2.0 * 3.141592658 ) / (double)h_in->GetNbinsX();
  double bin_theta_size = (h_theta_proj->GetBinCenter(2) - h_theta_proj->GetBinCenter(1)) * (3.1415/180.0);
  double theta_max = h_theta_proj->GetBinCenter(h_theta_proj->GetNbinsX()) + (h_theta_proj->GetBinCenter(2) - h_theta_proj->GetBinCenter(1))/2.0;

  std::cout << " >> theta max " << theta_max << std::endl;
  std::cout << " >> phi bin size " <<bin_phi_size << " theta bin size " << bin_theta_size << std::endl;
  std::cout << " new theta bin count " << h_theta_proj->GetNbinsX() << std::endl;
  TH1D *h_temp = new TH1D(Form("h_cs_%s",h_in->GetTitle()), Form("h_cs_%s",h_in->GetTitle()), h_theta_proj->GetNbinsX(), 0.0, 30.0);

  for( int b = 0; b < h_theta_proj->GetNbinsX(); b++ ){
    double sin_angle = TMath::Sin( h_theta_proj->GetBinCenter(b) * 3.1415/180.0 );
    double newbincontent= h_theta_proj->GetBinContent(b) / (lum_run2391 * bin_phi_size * bin_theta_size * sin_angle );

    h_temp->SetBinContent(b,newbincontent);
  }

  return h_temp;
}


int crossSectionExtractor(const char* infile, int run){

  TFile *fIn = new TFile(infile,"");
  
  if( fIn->IsZombie() ){
    std::cout << " Input file "<< fIn << " doesn't exist" << std::endl;
    std::cout << " bye " << std::endl;
    return 0;
  }

  std::vector< TH1F* > h_el_p_sect_final;
  std::vector< TH1F* > h_el_theta_sect_final;
  std::vector< TH1F* > h_el_phi_sect_final;

  std::vector< TH2F* > h_el_ptheta_sect_final;
  std::vector< TH2F* > h_el_phitheta_sect_final;

  TH2F *h_el_phitheta_final = (TH2F*)fIn->Get(Form("/kinematics/h_el_phitheta_final_run%d",run));
  
  for( int s = 1; s <= 6; s++ ){
    h_el_p_sect_final.push_back( (TH1F*)fIn->Get(Form("/kinematics/h_el_p_s%d_final",s) ) );
    h_el_theta_sect_final.push_back( (TH1F*)fIn->Get(Form("/kinematics/h_el_theta_s%d_final",s) ) );
    h_el_phi_sect_final.push_back( (TH1F*)fIn->Get(Form("/kinematics/h_el_phi_s%d_final",s) ) );

    h_el_ptheta_sect_final.push_back( (TH2F*)fIn->Get(Form("/kinematics/h_el_ptheta_s%d_final",s) ) );
    h_el_phitheta_sect_final.push_back( (TH2F*)fIn->Get(Form("/kinematics/h_el_phitheta_s%d_final",s) ) );
  }
  

  int phi_bins = h_el_phi_sect_final[0]->GetNbinsX();
  std::cout << " Number of phi bins: " << phi_bins << std::endl;

  int temp_entries=-1;
  int temp_bin = -1;
  for( int b = 0; b <= phi_bins; b++ ){
    int b_entries = h_el_phi_sect_final[0]->GetBinContent(b);
    std::cout << " bin entries " << b_entries << std::endl;

    if( b_entries > temp_entries ){
      temp_entries=b_entries;
      temp_bin = b;
    }   
  }

  std::cout << " Bin with max entries: " << temp_bin << " with entries: " << temp_entries << std::endl;
  int bin_min = temp_bin - 2;
  int bin_max = temp_bin + 2;

  std::cout << " Extending bins by 2 bins on each side; bin min: " << bin_min << " bin max: " << bin_max << std::endl;
  
  TH1D *h_theta_proj = h_el_phitheta_sect_final[0]->ProjectionY("proj_s1_2391",bin_min, bin_max);
  TCanvas *c1 = new TCanvas("c1","c1",1600,800);
  c1->Divide(2,1);
  c1->cd(1);
  h_el_phitheta_sect_final[0]->Draw("colz");
  c1->cd(2);
  h_theta_proj->Rebin(10);

  double lum_run2391 = 609158.862236702; //6.0915886E6;
  double bin_phi_size = (h_el_phi_sect_final[0]->GetBinCenter(bin_max) - h_el_phi_sect_final[0]->GetBinCenter(bin_min)) * 3.1415/180.0;
  double bin_theta_size = (h_theta_proj->GetBinCenter(2) - h_theta_proj->GetBinCenter(1)) * (3.1415/180.0);
  double theta_max = h_theta_proj->GetBinCenter(h_theta_proj->GetNbinsX()) +  (h_theta_proj->GetBinCenter(2) - h_theta_proj->GetBinCenter(1))/2.0;


  std::cout << " >> theta max " << theta_max << std::endl;
  std::cout << " >> phi bin size " <<bin_phi_size << " theta bin size " << bin_theta_size << std::endl;
  std::cout << " new theta bin count " << h_theta_proj->GetNbinsX() << std::endl;
  TH1D *h_temp = new TH1D("h_Temp","h_Temp", h_theta_proj->GetNbinsX(), 0.0, 30.0);

  for( int b = 0; b < h_theta_proj->GetNbinsX(); b++ ){
    double sin_angle = TMath::Sin( h_theta_proj->GetBinCenter(b) * 3.1415/180.0 );
    double newbincontent= h_theta_proj->GetBinContent(b) / (lum_run2391 * bin_phi_size * bin_theta_size * sin_angle * 0.1);

    h_temp->SetBinContent(b,newbincontent);
  }

  //h_theta_proj->Draw();
  h_temp->Draw();
  
  std::vector<TH1D*> cs_results;
  for( int s=0; s<6; s++ ){
    TH1D *h_temp_cs = getSectorCrossSection( h_el_phi_sect_final[s], h_el_phitheta_sect_final[s] );
    cs_results.push_back(h_temp_cs);

  }


  TCanvas *cs_out = new TCanvas("cs_out","cs_out",1200,800);
  cs_out->Divide(3,2);
  for( int s = 0; s < 6; s++ ){
    cs_out->cd(s+1);
    cs_results[s]->Draw();
  }


  return 0;
}
