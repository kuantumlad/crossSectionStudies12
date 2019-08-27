#include <vector>
#include <string>
#include <map>
#include <iostream>
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"


double getBinSize(TH1F*);

double getBinSize(TH1F *h_temp ){
  double p_bins = (double)h_temp->GetXaxis()->GetNbins();
  double p_min = h_temp->GetXaxis()->GetBinLowEdge(1);
  double p_max = h_temp->GetXaxis()->GetBinUpEdge(p_bins);
  std::cout << " nbins " << p_bins << std::endl;
  std::cout << " min " << p_min << std::endl;
  std::cout << " max " << p_max << std::endl;
  double bs  = (p_max-p_min)/p_bins;  
  std::cout << " bs " << bs << std::endl;
  return bs;
}


int acceptancePlotter(const char* inFileMC, int run, const char* field_config ){ 

 
  TFile *fData; 
  TFile *fMC; 
  Char_t tmpstr[80];
  Double_t fraction;

  fMC = new TFile(inFileMC,""); 
  if(fMC->IsZombie() ){
    cout<<"Input file doesn't exist!" << endl;
    cout<<"Exit program" << endl;
    return 0;
  }

  cout << "Reading from File: " << inFileMC << endl;

  // get relevant histograms to determine efficiency or geometrical acceptance
  //TH1F *h_rc_all_el_p = (TH1F*)fMC->Get("particle_histograms_selected/hist_electron_p");

  std::vector< TH1F* > v_h_el_purity_p;
  std::vector< TH1F* > v_h_el_purity_theta;
  std::vector< TH1F* > v_h_el_purity_phi;

  for( int bb =0; bb < 20; bb++) {
    TH1F *h_el_purity_p =(TH1F*)fMC->Get(Form("acceptance/h_el_purity_p_bb%d",bb));
    TH1F *h_el_purity_theta =(TH1F*)fMC->Get(Form("acceptance/h_el_purity_theta_bb%d",bb));
    TH1F *h_el_purity_phi =(TH1F*)fMC->Get(Form("acceptance/h_el_purity_phi_bb%d",bb));

    v_h_el_purity_p.push_back( h_el_purity_p );
    v_h_el_purity_theta.push_back( h_el_purity_theta );
    v_h_el_purity_phi.push_back( h_el_purity_phi );

  }


  TCanvas *c_p = new TCanvas("c_p","c_p",900,900);
  c_p->Divide(4,5);
  for( int bb = 0; bb < v_h_el_purity_p.size() ; bb++ ){
    c_p->cd(bb+1);
    double p_bin_size = getBinSize(v_h_el_purity_p[bb]);
    std::cout << " >. title " << v_h_el_purity_p[bb]->GetTitle() << std::endl;
    std::cout << " Plotting momentum purity plot with bin size " << p_bin_size << std::endl;
    v_h_el_purity_p[bb]->SetTitle(Form("P Purity #Delta %f GeV",p_bin_size));
    v_h_el_purity_p[bb]->GetXaxis()->SetTitle("Mntm (GeV)");
    v_h_el_purity_p[bb]->Draw();    
  }
  c_p->SaveAs(Form("h_el_purity_p_%s_r%d.pdf",field_config,run));

  TCanvas *c_th = new TCanvas("c_th","c_th",900,900);
  c_th->Divide(4,5);
  for( int bb = 0; bb < v_h_el_purity_theta.size() ; bb++ ){
    c_th->cd(bb+1);
    double bin_size = getBinSize(v_h_el_purity_theta[bb]);
    std::cout << " >. title " << v_h_el_purity_theta[bb]->GetTitle() << std::endl;
    std::cout << " Plotting momentum purity plot with bin size " << bin_size << std::endl;
    v_h_el_purity_theta[bb]->SetTitle(Form("Theta Purity #Delta %f deg",bin_size));
    v_h_el_purity_theta[bb]->GetXaxis()->SetTitle("#theta (deg)");
    v_h_el_purity_theta[bb]->Draw();    
  }
  c_th->SaveAs(Form("h_el_purity_theta_%s_r%d.pdf",field_config,run));

  TCanvas *c_phi = new TCanvas("c_phi","c_phi",900,900);
  c_phi->Divide(4,5);
  for( int bb = 0; bb < v_h_el_purity_phi.size() ; bb++ ){
    c_phi->cd(bb+1);
    double bin_size = getBinSize(v_h_el_purity_phi[bb]);
    std::cout << " >. title " << v_h_el_purity_phi[bb]->GetTitle() << std::endl;
    std::cout << " Plotting momentum purity plot with bin size " << bin_size << std::endl;
    v_h_el_purity_phi[bb]->SetTitle(Form("Phi Purity #Delta %f deg",bin_size));
    v_h_el_purity_phi[bb]->GetXaxis()->SetTitle("phi (deg)");
    v_h_el_purity_phi[bb]->Draw();    
  }
  c_phi->SaveAs(Form("h_el_purity_phi_%s_r%d.pdf",field_config,run));

  return 0;
}
