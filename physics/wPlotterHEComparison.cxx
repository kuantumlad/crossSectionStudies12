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



int wPlotterHEComparison(const char* infile, const char* infileMC, int run){

  TFile *fIn = new TFile(infile,"");
  TFile *fMC = new TFile(infileMC,"");
  TFile *fOut = new TFile(Form("cs_test_run%dV2.root",run),"RECREATE");
  
  if( fIn->IsZombie() || fMC->IsZombie()  ){
    std::cout << " Input file "<< fIn << " doesn't exist" << std::endl;
    std::cout << " bye " << std::endl;
    return 0;
  }


  std::vector<TH2F*> h_el_wtheta_sect;
  std::vector<TH2F*> mc_h_el_wtheta_sect;
  for( int s = 1; s <= 6; s++ ){
    h_el_wtheta_sect.push_back( (TH2F*)fIn->Get(Form("/kinematics/h_el_wtheta_s%d_final",s) ) );
    mc_h_el_wtheta_sect.push_back( (TH2F*)fMC->Get(Form("/kinematics/h_el_wtheta_s%d_final",s) ) );
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
      
      if ( tb < 40 && tb > 4 ){
	c_tb->cd(tb+1);

	//assuming the nbins for MC and data are the same (should be to make any comparison)
	TH1D *h_w_per_theta = new TH1D(Form("h_w_per_theta_s%d_b%d",s,tb),Form("h_w_per_theta_s%d_b%d",s,tb), n_bins_w ,h_el_wtheta_sect[s]->GetXaxis()->GetBinLowEdge(1), h_el_wtheta_sect[s]->GetXaxis()->GetBinUpEdge(n_bins_w) );
	TH1D *h_w_per_theta_mc = new TH1D(Form("h_w_per_theta_mc_s%d_b%d",s,tb),Form("h_w_per_theta_mc_s%d_b%d",s,tb), n_bins_w, mc_h_el_wtheta_sect[s]->GetXaxis()->GetBinLowEdge(1), mc_h_el_wtheta_sect[s]->GetXaxis()->GetBinUpEdge(n_bins_w) );
	double bcy = ((TAxis*)h_el_wtheta_sect[s]->GetYaxis())->GetBinCenter(tb);
      
	for( int wb = 0; wb < n_bins_w ; wb++ ){
	  double w_bin_content = h_el_wtheta_sect[s]->GetBinContent(wb,tb);
	  double bcx = ((TAxis*)h_el_wtheta_sect[s]->GetXaxis())->GetBinCenter(wb);	

	  double w_bin_content_mc = mc_h_el_wtheta_sect[s]->GetBinContent(wb,tb);	 
	  h_w_per_theta->SetBinContent(wb, w_bin_content);		
	  h_w_per_theta_mc->SetBinContent(wb,w_bin_content_mc);

	}
	double low_edge = h_el_wtheta_sect[s]->GetXaxis()->GetBinLowEdge(tb);
	double up_edge = h_el_wtheta_sect[s]->GetXaxis()->GetBinUpEdge(tb);
	
	double scale_w = (double)h_w_per_theta->GetMaximum()/(double)h_w_per_theta_mc->GetMaximum();
	std::cout << " dat max " << h_w_per_theta->GetMaximum() << " mc " << h_w_per_theta_mc->GetMaximum()<< "scale " << scale_w << std::endl;
	h_w_per_theta->Scale(1.0/scale_w);


	h_w_per_theta->SetTitle(Form("S%d: Theta Bin %f - %f ", s, low_edge, up_edge));
	h_w_per_theta->GetXaxis()->SetTitle("W (GeV)");
	h_w_per_theta->SetLineColor(kRed);
	h_w_per_theta_mc->SetLineColor(kBlue);
	h_w_per_theta->Draw("hist");
	h_w_per_theta_mc->Draw("hist+same");

	c_tb->Update();
	TVirtualPad *p1 = c_tb->GetPad(1);
	double ymax = gPad->GetUymax();
	double ymin = gPad->GetUymin();
	std::cout << ymax << std::endl;
	TLine *l = new TLine(0.938,ymin, 0.938, ymax);
	l->SetLineColor(kRed);
	l->SetLineWidth(1);
	l->SetLineStyle(2);
	l->Draw("same");
	
	TF1 *fit_mean = new TF1("fit_mean","gaus", 0.5, 1.3);
	//fit_mean->SetParameter(0,fitter->GetParameter(0));
	//fit_mean->SetParameter(1,fitter->GetParameter(1));
	//fit_mean->SetParameter(2,fitter->GetParameter(2));
	//fit_mean->SetLineColor(kGreen);
	//fit_mean->Draw("same");	

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// get hist parameters to set fit limits    
	double xlow,xhigh,histmax;
	int binlow,binhigh,binmax;
	binmax = -1;//h_w_per_theta->GetMaximumBin();
	histmax = h_w_per_theta->GetMaximum();   
	//rather than get the maximum value in the entire histogram
	//we get the max over a specific range
	int target_bin = h_w_per_theta->FindBin(0.938);
	int target_bin_min = target_bin - 10;
	int target_bin_max = target_bin + 6;
	for( int ii = target_bin_min; ii < target_bin_max; ii++ ){
	  //std::cout<< " ii " << ii << "bin max " << binmax << " bin cont" <<  h_w_per_theta->GetBinContent(ii) << std::endl;
	  if ( h_w_per_theta->GetBinContent(ii) > h_w_per_theta->GetBinContent(binmax) ) binmax=ii;
	}
	binlow=binmax;
	binhigh=binmax;

	// The 0.65 parameter can be changed, this basically means start at the peak and work your way left and right
	// until you've gotten to 65% of the max amplitude.
	while(h_w_per_theta->GetBinContent(binhigh++) >= .65*histmax&&binhigh<=h_w_per_theta->GetNbinsX()){};
	while(h_w_per_theta->GetBinContent(binlow--) >= .65*histmax&&binlow>=1){};
    
	xlow = h_w_per_theta->GetBinLowEdge(binlow);
	xhigh = h_w_per_theta->GetBinLowEdge(binhigh+1);

	std::cout << " xlow " << xlow << " xhigh " << xhigh << std::endl;
	TF1 *fit_mean = new TF1("fit_mean_data","gaus", xlow, xhigh);
	h_w_per_theta->Fit("fit_mean_data","R");
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	/*TF1 *fitter = new TF1("fitter","gaus + [3] + [4]*x + [5]*x*x ", 0.5, 1.3);  
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
	*/
      }
    }
    c_tb->SaveAs(Form("h_el_w_per_theta_bin_sect%d_r%d.pdf",s,run));
  }

  return 0;
}



