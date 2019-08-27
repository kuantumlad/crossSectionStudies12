#include <iostream>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TMath.h>
#include <TH2F.h>
#include <vector>
#include <map>
#include <TLine.h>



int wPlotterHE(const char* infile, int run){

  TFile *fIn = new TFile(infile,"");
  TFile *fOut = new TFile(Form("cs_test_run%dV2.root",run),"RECREATE");
  
  if( fIn->IsZombie() ){
    std::cout << " Input file "<< fIn << " doesn't exist" << std::endl;
    std::cout << " bye " << std::endl;
    return 0;
  }


  std::vector<TH2F*> h_el_wtheta_sect;
  for( int s = 1; s <= 6; s++ ){
    h_el_wtheta_sect.push_back( (TH2F*)fIn->Get(Form("/kinematics/h_el_wtheta_s%d_final",s) ) );
  }

  
  TLine *l_mp = new TLine(0.938,5.0, 0.938, 20.0);

  TCanvas *c0 = new TCanvas("c0","c0",900,900);
  c0->Divide(3,2);
  for( int s  = 0; s < 6; s++ ){ 
    c0->cd(s+1);

    h_el_wtheta_sect[s]->SetTitle(Form(" #theta vs W S%d",s+1));
    h_el_wtheta_sect[s]->GetXaxis()->SetTitle("W (GeV)");
    h_el_wtheta_sect[s]->GetYaxis()->SetTitle("#theta (deg)");
    h_el_wtheta_sect[s]->Draw("colz");

    l_mp->SetLineColor(kRed);
    l_mp->SetLineWidth(2);
    l_mp->Draw("same");
  }

  c0->SaveAs(Form("h_el_wtheta_per_sec_run%d.pdf",run));


  for( int s  = 0; s < 6; s++ ){
    int n_bins_w = h_el_wtheta_sect[s]->GetNbinsX();    
    int n_bins_theta = h_el_wtheta_sect[s]->GetNbinsY();    

    std::cout << " w bins " << n_bins_w << " theta bins " << n_bins_theta <<std::endl;
    
    TCanvas *c_tb = new TCanvas(Form("c_tb_s%d",s),Form("c_tb_s%d",s),900,900);

    c_tb->Divide(5,8);
    for( int tb = 0; tb < n_bins_theta; tb++ ){
      
      if ( tb < 40 ){
	c_tb->cd(tb+1);

	
	TH1D *h_w_per_theta = new TH1D(Form("h_w_per_theta_s%d_b%d",s,tb),Form("h_w_per_theta_s%d_b%d",s,tb), n_bins_w ,h_el_wtheta_sect[s]->GetXaxis()->GetBinLowEdge(1), h_el_wtheta_sect[s]->GetXaxis()->GetBinUpEdge(n_bins_w) );
	double bcy = ((TAxis*)h_el_wtheta_sect[s]->GetYaxis())->GetBinCenter(tb);
      
	for( int wb = 0; wb < n_bins_w ; wb++ ){
	  double w_bin_content = h_el_wtheta_sect[s]->GetBinContent(wb,tb);
	  double bcx = ((TAxis*)h_el_wtheta_sect[s]->GetXaxis())->GetBinCenter(wb);
	
	  h_w_per_theta->SetBinContent(wb, w_bin_content);		
	}


	h_w_per_theta->Draw();
	TF1 *fitter = new TF1("fitter","gaus + [3] + [4]*x + [5]*x*x ", 0.5, 1.3);  
	fitter->SetParName(0, "Norm");
	fitter->SetParName(1, "Factor");
	fitter->SetParName(2, "#sigma");
	fitter->SetParName(3, "const term");
	fitter->SetParName(4, "linear term");
	fitter->SetParName(5, "quad term");
	fitter->SetParameters(25.0, 0.9, 0.4, 1.0, 0.1, 0.1);
	fitter->SetLineColor(kRed);
	h_w_per_theta->Fit("fitter","R");

	TF1 *fit_mean = new TF1("fit_mean","gaus", 0.5, 1.3);
	fit_mean->SetParameter(0,fitter->GetParameter(0));
	fit_mean->SetParameter(1,fitter->GetParameter(1));
	fit_mean->SetParameter(2,fitter->GetParameter(2));
	fit_mean->SetLineColor(kGreen);
	fit_mean->Draw("same");
	
	TF1*fit_poly = new TF1("fit_poly","[0] + [1]*x + [2]*x*x",0.5, 1.3);
	fit_poly->SetParameter(0,fitter->GetParameter(3));
	fit_poly->SetParameter(1,fitter->GetParameter(4));
	fit_poly->SetParameter(2,fitter->GetParameter(5));
	fit_poly->SetLineColor(kBlue);
	fit_poly->Draw("same");	  
      }
    }
    c_tb->SaveAs(Form("h_el_w_per_theta_bin_sect%d_r%d.pdf",s,run));
  }

  return 0;
}



