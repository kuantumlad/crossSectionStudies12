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

TF1 *fitHistogramRange(TH1D*, int, int );

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

  int n_sectors =6;

  ////////////////
  /// containers to plot the mean and sigma per theta per sector
  ////////////////
  std::vector<TGraphErrors* > s_sim_mean;
  std::vector<TGraphErrors* > s_sim_sig;
  std::vector<TGraphErrors* > s_data_mean;
  std::vector<TGraphErrors* > s_data_sig;
  

  TLine *l_mp = new TLine(0.938,5.0, 0.938, 20.0);

  TCanvas *c0 = new TCanvas("c0","c0",900,900);
  c0->Divide(3,2);
  for( int s  = 0; s < n_sectors; s++ ){ 
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

  ofstream outputWCuts;
  ofstream outputWCuts2;

  for( int s  = 0; s < n_sectors; s++ ){
    int n_bins_w = h_el_wtheta_sect[s]->GetNbinsX();    
    int n_bins_theta = h_el_wtheta_sect[s]->GetNbinsY();    

    std::string parentDirectory = "/w/hallb-scifs17exp/clas12/bclary/CLAS12/electron_studies/physics/parameters/";
    std::string f_out_w_name = parentDirectory+"whe_cut_limits_run"+std::to_string(run)+"_s"+std::to_string(s)+".txt";
    std::cout << " Creating W cuts output file " << f_out_w_name << std::endl;
    outputWCuts.open(f_out_w_name);

    std::string f_out_w_name2 = parentDirectory+"whe_sim_cut_limits_run"+std::to_string(run)+"_s"+std::to_string(s)+".txt";
    std::cout << " Creating W cuts output file " << f_out_w_name2 << std::endl;
    outputWCuts2.open(f_out_w_name2);


    std::cout << " w bins " << n_bins_w << " theta bins " << n_bins_theta <<std::endl;
    
    TCanvas *c_tb = new TCanvas(Form("c_tb_s%d",s),Form("c_tb_s%d",s),900,900);

    c_tb->Divide(5,8);
    




    //containers for grapherrors 
    TH1F *h_temp_theta_wmean = new TH1F(Form("h_temp_theta_wmean_s%d",s),Form("h_temp_theta_wmean_s%d",s),n_bins_theta,h_el_wtheta_sect[s]->GetYaxis()->GetBinLowEdge(1),h_el_wtheta_sect[s]->GetYaxis()->GetBinUpEdge(n_bins_theta) );
    TH1F *h_temp_theta_wsig = new TH1F(Form("h_temp_theta_wsig_s%d",s),Form("h_temp_theta_wsig_s%d",s),n_bins_theta,h_el_wtheta_sect[s]->GetYaxis()->GetBinLowEdge(1),h_el_wtheta_sect[s]->GetYaxis()->GetBinUpEdge(n_bins_theta) );
    std::vector<double> s_tc;
    std::vector<double> s_tc_er;
    std::vector<double> s_me;
    std::vector<double> s_si;
    std::vector<double> s_me_er;
    std::vector<double> s_si_er;

    std::vector<double> s_mc_tc;
    std::vector<double> s_mc_tc_er;
    std::vector<double> s_mc_me;
    std::vector<double> s_mc_si;
    std::vector<double> s_mc_me_er;
    std::vector<double> s_mc_si_er;
    

    for( int tb = 0; tb < n_bins_theta; tb++ ){
      
      if ( tb < 40 && tb > 4 ){
	c_tb->cd(tb+1);

	//assuming the nbins for MC and data are the same (should be to make any comparison)
	TH1D *h_w_per_theta = new TH1D(Form("h_w_per_theta_s%d_b%d",s,tb),Form("h_w_per_theta_s%d_b%d",s,tb), n_bins_w ,h_el_wtheta_sect[s]->GetXaxis()->GetBinLowEdge(1), h_el_wtheta_sect[s]->GetXaxis()->GetBinUpEdge(n_bins_w) );
	TH1D *h_w_per_theta_mc = new TH1D(Form("h_w_per_theta_mc_s%d_b%d",s,tb),Form("h_w_per_theta_mc_s%d_b%d",s,tb), n_bins_w, mc_h_el_wtheta_sect[s]->GetXaxis()->GetBinLowEdge(1), mc_h_el_wtheta_sect[s]->GetXaxis()->GetBinUpEdge(n_bins_w) );
	double bcy = ((TAxis*)h_el_wtheta_sect[s]->GetYaxis())->GetBinCenter(tb);
	s_tc.push_back(bcy);
	s_tc_er.push_back(0);
	s_mc_tc.push_back(bcy);
	s_mc_tc_er.push_back(0);

	
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
	l->SetLineColor(kGreen);
	l->SetLineWidth(1);
	//l->SetLineStyle(1);
	l->Draw("same");
	
	TF1 *fit_mean = new TF1("fit_mean","gaus", 0.5, 1.3);
	//fit_mean->SetParameter(0,fitter->GetParameter(0));
	//fit_mean->SetParameter(1,fitter->GetParameter(1));
	//fit_mean->SetParameter(2,fitter->GetParameter(2));
	//fit_mean->SetLineColor(kGreen);
	//fit_mean->Draw("same");	

	TF1* fit_data = fitHistogramRange(h_w_per_theta, 10, 6);
	h_w_per_theta->Fit(fit_data,"R");
	TF1* fit_sim = fitHistogramRange(h_w_per_theta_mc, 5, 5);
	h_w_per_theta_mc->Fit(fit_sim,"R");
	
	fit_data->SetLineColor(kRed);
	fit_sim->SetLineColor(kBlue);
	fit_data->Draw("same");
	fit_sim->Draw("same");
	//write cut mean and sigma for each theta bin
	std::cout << " bin " <<tb<< " bcy " << bcy << " mean " << fit_data->GetParameter(1) << std::endl;
	std::cout << " bin bcy " << bcy << " sig " << fit_data->GetParameter(2) << std::endl;
	s_me.push_back(fit_data->GetParameter(1));
	s_si.push_back(fit_data->GetParameter(2));
	s_me_er.push_back(fit_data->GetParError(1));
	s_si_er.push_back(fit_data->GetParError(2));

	s_mc_me.push_back(fit_sim->GetParameter(1));
	s_mc_si.push_back(fit_sim->GetParameter(2));
	s_mc_me_er.push_back(fit_sim->GetParError(1));
	s_mc_si_er.push_back(fit_sim->GetParError(2));

	h_temp_theta_wmean->SetBinContent(tb, fit_data->GetParameter(1) );
	h_temp_theta_wsig->SetBinContent(tb, fit_data->GetParameter(2) );
	
	outputWCuts << tb << " " << fit_data->GetParameter(1) << " " << fit_data->GetParameter(2) << std::endl;
	outputWCuts2 << tb << " " << fit_sim->GetParameter(1) << " " << fit_sim->GetParameter(2) << std::endl;

      }
    }
    s_data_mean.push_back( new TGraphErrors(s_tc.size(), &s_tc[0], &s_me[0], &s_tc_er[0], &s_me_er[0]) );
    s_data_sig.push_back( new TGraphErrors(s_tc.size(), &s_tc[0], &s_si[0], &s_tc_er[0], &s_si_er[0]) );
    s_sim_mean.push_back( new TGraphErrors(s_mc_tc.size(), &s_mc_tc[0], &s_mc_me[0], &s_mc_tc_er[0], &s_mc_me_er[0]) );
    s_sim_sig.push_back( new TGraphErrors(s_mc_tc.size(), &s_mc_tc[0], &s_mc_si[0], &s_mc_tc_er[0], &s_mc_si_er[0]) );
    outputWCuts.close();
    outputWCuts2.close();

    c_tb->SaveAs(Form("h_el_w_per_theta_bin_sect%d_r%d.pdf",s,run));
  }

  TCanvas *c_mean = new TCanvas("c_mean","c_mean",900,900);
  c_mean->Divide(2,3);
  std::vector<TMultiGraph*> v_mg_mean;
  std::vector<TMultiGraph*> v_mg_sig;
  for(int ss = 0; ss<s_data_mean.size(); ss++ ){
    c_mean->cd(ss+1);
    v_mg_mean.push_back(new TMultiGraph() );
    s_data_mean[ss]->SetTitle(Form("S%d Mean",ss));
    s_data_mean[ss]->SetMarkerStyle(20);
    s_data_mean[ss]->SetMarkerSize(0.8);    
    s_data_mean[ss]->SetMarkerColor(kRed);
    s_sim_mean[ss]->SetMarkerStyle(20);
    s_sim_mean[ss]->SetMarkerSize(0.8);
    s_sim_mean[ss]->SetMarkerColor(kBlue);
    v_mg_mean[ss]->SetTitle(Form("Data vs Sim Mean S%d",ss));
    v_mg_mean[ss]->GetXaxis()->SetTitle("Theta (deg)");
    v_mg_mean[ss]->GetYaxis()->SetTitle("Mean (GeV)");
    v_mg_mean[ss]->Add(s_data_mean[ss]);
    v_mg_mean[ss]->Add(s_sim_mean[ss]);
    v_mg_mean[ss]->Draw("APE");
  }
  c_mean->SaveAs(Form("g_w_fit_mean_r%d.pdf",run));

  TCanvas *c_sigma = new TCanvas("c_sigma","c_sigma",900,900);
  c_sigma->Divide(2,3);
  for(int ss = 0; ss<s_data_sig.size(); ss++ ){
    c_sigma->cd(ss+1);
    v_mg_sig.push_back(new TMultiGraph() );
    s_data_sig[ss]->SetTitle(Form("S%d Sig",ss));
    s_data_sig[ss]->SetMarkerStyle(20);
    s_data_sig[ss]->SetMarkerSize(0.8);    
    s_data_sig[ss]->SetMarkerColor(kRed);
    s_sim_sig[ss]->SetMarkerStyle(20);
    s_sim_sig[ss]->SetMarkerSize(0.8);
    s_sim_sig[ss]->SetMarkerColor(kBlue);
    v_mg_sig[ss]->SetTitle(Form("Data vs Sim Mean S%d",ss));
    v_mg_sig[ss]->GetXaxis()->SetTitle("Theta (deg)");
    v_mg_sig[ss]->GetYaxis()->SetTitle("Sigma (GeV)");
    v_mg_sig[ss]->Add(s_data_sig[ss]);
    v_mg_sig[ss]->Add(s_sim_sig[ss]);
    v_mg_sig[ss]->Draw("APE"); 
  }
  c_sigma->SaveAs(Form("g_w_fit_sigma_r%d.pdf",run));

  return 0;
}

TF1 *fitHistogramRange(TH1D *h_temp, int s_min, int s_max){

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // get hist parameters to set fit limits    
  double xlow,xhigh,histmax;
  int binlow,binhigh,binmax;
  binmax = -1;//h_w_per_theta->GetMaximumBin();
  histmax = h_temp->GetMaximum();   
  //rather than get the maximum value in the entire histogram
  //we get the max over a specific range
  int target_bin = h_temp->FindBin(0.938);
  int target_bin_min = target_bin - s_min;
  int target_bin_max = target_bin + s_max;
  for( int ii = target_bin_min; ii < target_bin_max; ii++ ){
    //std::cout<< " ii " << ii << "bin max " << binmax << " bin cont" <<  h_w_per_theta->GetBinContent(ii) << std::endl;
    if ( h_temp->GetBinContent(ii) > h_temp->GetBinContent(binmax) ) binmax=ii;
  }
  binlow=binmax;
  binhigh=binmax;

  // The 0.65 parameter can be changed, this basically means start at the peak and work your way left and right
  // until you've gotten to 65% of the max amplitude.
  while(h_temp->GetBinContent(binhigh++) >= .65*histmax&&binhigh<=h_temp->GetNbinsX()){};
  while(h_temp->GetBinContent(binlow--) >= .65*histmax&&binlow>=1){};
    
  xlow = h_temp->GetBinLowEdge(binlow);
  xhigh = h_temp->GetBinLowEdge(binhigh+1);

  std::cout << " xlow " << xlow << " xhigh " << xhigh << std::endl;
  TF1 *fit_mean = new TF1("fit_mean_data","gaus", xlow, xhigh);
  return fit_mean;
}

