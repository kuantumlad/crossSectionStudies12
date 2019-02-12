#include <iostream>
#include <TF1.h>
#include <vector>
#include <TH2D.h>
#include <TH2F.h>
#include <TH1D.h>
#include <map>
#include <string>
#include <sstream>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TCanvas.h>

TF1* fitHisto(TH1D* htemp){
  //////////////////////////////////////////////////
  // Start Andrew's fitting method:
  double xlow,xhigh,histmax;
  int binlow,binhigh,binmax;
  binmax = htemp->GetMaximumBin();
  histmax = htemp->GetMaximum();
  binlow=binmax;
  binhigh=binmax;

  // The 0.65 parameter can be changed, this basically means start at the peak and work your way left and right
  // until you've gotten to 65% of the max amplitude.
  while(htemp->GetBinContent(binhigh++) >= .65*histmax&&binhigh<=htemp->GetNbinsX()){};
  while(htemp->GetBinContent(binlow--) >= .65*histmax&&binlow>=1){};
    
  xlow = htemp->GetBinLowEdge(binlow);
  xhigh = htemp->GetBinLowEdge(binhigh+1);
    
  htemp->Fit("gaus","","",xlow,xhigh);

  TF1 *ftemp = (TF1*) htemp->GetListOfFunctions()->FindObject("gaus");

  return ftemp;

  // End
  /////////////////////////////////////////////

}


TF1* fitHistogram( TH1D *h_temp ){
      
  //System.out.println(" >> FITTING HISTOGRAM " + h_temp.getName() );
  double xlow, xhigh, histmax;
  int binlow, binhigh, binmax;

  double percentofmax = 0.45; // was 45%
    
  TF1 *fit = NULL;

  //if( h_temp.getEntries() > 0 ){
  binmax = h_temp->GetMaximumBin();
  histmax = h_temp->GetMaximum();
  binlow = binmax;
  binhigh = binmax;
  
  while( h_temp->GetBinContent(binhigh++) >= percentofmax*histmax && binhigh <= h_temp->GetNbinsX() ){}
  while( h_temp->GetBinContent(binlow--) >= percentofmax*histmax && binlow > 1 ){}
  
  xlow = h_temp->GetXaxis()->GetBinCenter(binlow) - ( h_temp->GetXaxis()->GetBinCenter(binlow) - h_temp->GetXaxis()->GetBinCenter(binlow+1))/2.0; // needs to be low edge, only center now
  xhigh = h_temp->GetXaxis()->GetBinCenter(binhigh+1) - (h_temp->GetXaxis()->GetBinCenter(binhigh+1) - h_temp->GetXaxis()->GetBinCenter(binhigh))/2.0; // needs to be low edge, only center now
												  
  std::cout << " >> values used " << xlow <<  " " << xhigh << " " << histmax << std::endl;
  
  TF1 *fit_temp = new TF1("fit_temp","gaus", xlow, xhigh );
  fit_temp->SetParameter(0, histmax);
  fit_temp->SetParameter(1, h_temp->GetMean() );
  fit_temp->SetParameter(2, h_temp->GetRMS() );
  
  h_temp->Fit(fit_temp, "REQ"); //was only R at first
  fit = fit_temp;  
  
  std::cout << " >> PARAMETER SET " << fit_temp->GetParameter(0) << " " << fit_temp->GetParameter(1) << " " << fit_temp->GetParameter(2) << std::endl;


  return fit;	    
}



int fitSamplingFraction(const char* inFile, int run){

  
  TFile *fIn = new TFile(inFile,"");
  
  TCanvas *c_el_sf = new TCanvas("c_el_sf", "c_el_sf", 1800, 900);
  c_el_sf->Divide(3,2);

  std::vector<TGraphErrors*> g_el_sf_mean;
  std::vector<TGraphErrors*> g_el_sf_sigma;

  std::vector<TH2F*> h_el_sfp;
  std::vector<TCanvas*> v_can;


  //get SF from histograms
  for( int s = 1; s <= 6; s++ ){
    c_el_sf->cd(s);
   
    TH2F *h_el_sf = (TH2F*)fIn->Get(Form("FD_PID_electron_EC_plots/EC_total_sampling_fraction_sec%d_cut_00",s));
    std::cout << " number of entries for sector " << s << "  " << h_el_sf->GetEntries() << std::endl;
    h_el_sfp.push_back(h_el_sf);
    h_el_sf->Draw("colz");  
  
    v_can.push_back(new TCanvas(Form("c_s%d",s), Form("c_s%d",s), 800, 800 ));
    v_can[s-1]->Divide(5,5);

  }


  for( int s = 0; s < h_el_sfp.size(); s++ ){
    int nbins_x = h_el_sfp[s]->GetNbinsX();
    int ndiv = 20;
    int nfits = nbins_x/ndiv;
    
    int start_bin = 0;
    int end_bin = ndiv;

    std::vector<double> temp_bin_center;
    std::vector<double> temp_bin_center_err;
    std::vector<double> temp_fit_mean;
    std::vector<double> temp_fit_mean_err;
    std::vector<double> temp_fit_sigma;
    std::vector<double> temp_fit_sigma_err;


    std::cout << " performing fits for sector " << s << std::endl;
    std::cout << " N bins in X " << nbins_x << std::endl;

    
    for( int f = 0; f < nfits; f++ ){


      TH1D *h_sf = h_el_sfp[s]->ProjectionY(Form("h_el_sfp_s%d_f%d",s,f), start_bin, end_bin);
      h_sf->Rebin(2);
      std::cout << " get n bins again " << h_sf->GetNbinsX() << std::endl;
      std::cout << " numberof entries in slices " << h_sf->GetEntries() << std::endl;
      std::cout << " projection of ec sf sector " << s << " for fit range " << f << " " << start_bin << " end bin " << end_bin << std::endl;
      //fit histogram
      
      double min_p = h_el_sfp[s]->GetXaxis()->GetBinCenter(start_bin) - h_el_sfp[s]->GetXaxis()->GetBinWidth(start_bin)/2.0;
      double max_p = h_el_sfp[s]->GetXaxis()->GetBinCenter(end_bin) + h_el_sfp[s]->GetXaxis()->GetBinWidth(end_bin)/2.0;
      double slice_center = ( min_p + max_p )/2.0;
      std::cout << " bin center " << slice_center << std::endl;

      start_bin+=ndiv;
      end_bin+=ndiv;

      if ( f  < 15 ) continue;

      
      TF1 *fit_sf = fitHistogram(h_sf);

      //get bin center of slices
      double mean = fit_sf->GetParameter(1);
      double sig = fit_sf->GetParameter(2);
      double mean_err  = 0.0;
      double sig_err = 0.0;

      if( f > 15 && f < 45 ){

      	v_can[s]->cd(f);
	h_sf->Draw();
	fit_sf->Draw("same");
      }



      temp_bin_center.push_back(slice_center);
      temp_bin_center_err.push_back(0.0);
      temp_fit_mean.push_back(mean);
      temp_fit_sigma.push_back(sig);
      temp_fit_mean_err.push_back(mean_err);
      temp_fit_sigma_err.push_back(sig_err);
     
    }
    
    g_el_sf_mean.push_back( new TGraphErrors( temp_bin_center.size(), &(temp_bin_center[0]), &(temp_fit_mean[0]), &(temp_bin_center_err[0]), &(temp_fit_mean_err[0]) ) );
    g_el_sf_sigma.push_back( new TGraphErrors( temp_bin_center.size(), &(temp_bin_center[0]), &(temp_fit_sigma[0]), &(temp_bin_center_err[0]), &(temp_fit_sigma_err[0]) ) );
      
    

  }

  TCanvas *c_temp_sfmean = new TCanvas("c_temp_sfmean", "c_temp_sfmean",900,900);
  c_temp_sfmean->Divide(3,2);

  TCanvas *c_temp_sfsig = new TCanvas("c_temp_sfsig", "c_temp_sfsig",900,900);
  c_temp_sfsig->Divide(3,2);
  
  // now fit the mean sf and sigma sf distributions to get constants
  // after fitting write the constants to a text file for the run
  std::vector<std::string> sf_mean_parm;
  std::vector<std::string> sf_sigma_parm;
  std::vector<std::string> par_names_mean = {"A", "B", "C", "D"};
  std::vector<std::string> par_names_sig = {"A", "B"};
  
  for( int i = 0; i < par_names_mean.size(); i++){
    sf_mean_parm.push_back(par_names_mean[i]+" \t" + "6" + "\t" + "v" + "\t");   
  }

  for( int i = 0; i < par_names_sig.size(); i++){
    sf_sigma_parm.push_back(par_names_sig[i]+" \t" + "6" + "\t" + "v" + "\t");   
  }

  for(int s = 0; s < g_el_sf_mean.size(); s++ ){
    c_temp_sfmean->cd(s+1);

    TF1 *fit_mean  = new TF1(Form("fit_mean_s%d",s),"[0]*( 1 + x/sqrt(x*x + [1])) + [2]*x + [3]*x*x", 0.90, 2.16);    
    fit_mean->SetParameter(0,0.25);
    fit_mean->SetParameter(1,0.01);
    fit_mean->SetParameter(2,0.01);
    fit_mean->SetParameter(3,0.01);
    g_el_sf_mean[s]->Fit(Form("fit_mean_s%d",s), "R");
    g_el_sf_mean[s]->SetTitle(Form("EL SF MEAN FOR SECTOR %d",s));
    g_el_sf_mean[s]->SetMarkerColor(kBlack);
    g_el_sf_mean[s]->SetMarkerStyle(8);
    g_el_sf_mean[s]->Draw("AP");
    fit_mean->SetLineColor(kRed);
    fit_mean->Draw("same");

    c_temp_sfsig->cd(s+1);
    TF1 *fit_sigma = new TF1(Form("fit_sigma_s%d",s), "[0] + [1]/x", 0.9, 2.16);
    fit_sigma->SetParameter(0,0.02);
    fit_sigma->SetParameter(1,0.0);
    g_el_sf_sigma[s]->Fit(Form("fit_sigma_s%d",s), "R");
    g_el_sf_sigma[s]->SetTitle(Form("EL SF SIGMA FOR SECTOR %d",s));
    g_el_sf_sigma[s]->SetMarkerColor(kBlack);
    g_el_sf_sigma[s]->SetMarkerStyle(8);
    g_el_sf_sigma[s]->Draw("AP");
    fit_sigma->SetLineColor(kRed);
    fit_sigma->Draw("same");

    for( int p = 0; p < fit_mean->GetNpar(); p++ ){
      sf_mean_parm[p] = sf_mean_parm[p] + " " + std::to_string(fit_mean->GetParameter(p));      
    }
    
    for( int p = 0; p < fit_sigma->GetNpar(); p++ ){
      sf_sigma_parm[p] = sf_sigma_parm[p] + " " + std::to_string(fit_sigma->GetParameter(p));      
    }
    

  }
   


  /////////////////////////////////////////////////////
  
    
  // write out sampling fraction variables to txt file
  ofstream outputSFCuts;
  std::string parentDirectory = "/w/hallb-scifs17exp/clas12/bclary/CLAS12/electron_studies/pid/parameters/";
  std::string f_out_sf_name = parentDirectory+"sf_cut_mean_run"+std::to_string(run)+".txt";
  std::cout << " Creating SF cuts output file " << f_out_sf_name << std::endl;
  outputSFCuts.open(f_out_sf_name);
  for( int l = 0; l < sf_mean_parm.size(); l++ ){
    std::cout << " >> " << sf_mean_parm[l] << std::endl;
    outputSFCuts << sf_mean_parm[l] << std::endl;
  }

  ofstream outputSFSigmaCuts;
  std::string f_out_sf_sig_name = parentDirectory+"sf_cut_sigma_run"+std::to_string(run)+".txt";
  std::cout << " Creating SF SIGMA cuts output file " << f_out_sf_name << std::endl;
  outputSFSigmaCuts.open(f_out_sf_sig_name);
  for( int l = 0; l < sf_sigma_parm.size(); l++ ){
    outputSFSigmaCuts << sf_sigma_parm[l] << std::endl;
  }

  outputSFCuts.close();
  outputSFSigmaCuts.close();

  return 0;
}

