#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <string>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TF1.h>
#include <TGraphErrors.h>




double testFunction1(double x){
  return 1.0/pow(x,4);
}

int exBinCenter(){

  TFile *f_test = new TFile("test_bin_center.root","RECREATE");
  

  double bin_min = 8.1;
  double bin_max = 28.1;
  double n_bin_size = 1.0;
  double n_bin_size_coarse = 0.010 * 10;
  int n_bins = (bin_max - bin_min )/n_bin_size;
  int n_bins_coarse = (bin_max - bin_min )/n_bin_size_coarse;
  std::cout << " BIN MIN " << bin_min << " BIN MAX " << bin_max << " N BINS " << n_bins << std::endl;
  
  TF1 *f_cs = new TF1("f_cs","1.0/(x^4)",bin_min,bin_max);
  TH1D *h_test_cs = new TH1D("h_test_cs","h_test_cs",(int)n_bins, bin_min, bin_max );
   
  for( int b = 1; b < n_bins; b++ ){
    double test_cs = testFunction1(h_test_cs->GetBinCenter(b));
    h_test_cs->SetBinContent(b,test_cs);       
  }

  std::vector< double> v_nom_cs;
  std::vector< double> v_corr_cs;
  std::vector<double> v_nom_corr_diff;
  std::vector< double > v_x;
  std::vector< double > v_nom_cs_err;
  std::vector< double > v_corr_cs_err;
  std::vector< double > v_nom_corr_cs_err;
  std::vector< double > v_x_err;

  for( int b=1; b < n_bins; b++ ){
    
    for( int bf = 1; bf < n_bins_fine ;bf++ ){


    }
    
  
    double bin_left = h_test_cs->GetBinLowEdge(b);
    double bin_half = h_test_cs->GetBinWidth(b)/2.0;
    double bin_right = h_test_cs->GetBinCenter(b) + bin_half;
    double bin_center = h_test_cs->GetBinCenter(b);

    double funct_left = testFunction1(bin_left);
    double funct_right = testFunction1(bin_right);
    double funct_center = testFunction1(bin_center);

    double funct_avg = (funct_left + funct_right)/2.0;
    double funct_diff = fabs(funct_avg - funct_center);
    double funct_ratio = fabs(funct_avg/funct_center);
    double funct_ratio_inv = fabs(funct_center/funct_avg);
    
    std::cout << " LOWER BIN X " << bin_left << " FUNCTION " << funct_left << std::endl;
    std::cout << " UPPER BIN " << bin_right << " FUNCTION " << funct_right << std::endl;
    
    std::cout << " FUNCT CENTER " << h_test_cs->GetBinCenter(b) << " FUNCT AVG CENTER " << funct_avg << std::endl;

    v_nom_cs.push_back(funct_center);
    v_corr_cs.push_back(funct_avg);
    v_nom_corr_diff.push_back(funct_ratio_inv);
    v_x.push_back(h_test_cs->GetBinCenter(b));
    v_nom_cs_err.push_back(0.0);
    v_corr_cs_err.push_back(0.0);
    v_nom_corr_cs_err.push_back(0.0);
    v_x_err.push_back(0.0);
 
  }

  TGraphErrors *g_cs_nom = new TGraphErrors(v_x.size(),&(v_x[0]), &(v_nom_cs[0]),&(v_nom_cs_err[0]), &(v_x_err[0]) );
  g_cs_nom->SetTitle("graph_cs_nom");

  TGraphErrors *g_cs_corr = new TGraphErrors(v_x.size(),&(v_x[0]), &(v_corr_cs[0]),&(v_corr_cs_err[0]), &(v_x_err[0]) );
  g_cs_corr->SetTitle("graph_cs_corr");

  TGraphErrors *g_cs_diff = new TGraphErrors(v_x.size(),&(v_x[0]), &(v_nom_corr_diff[0]),&(v_nom_corr_cs_err[0]), &(v_x_err[0]) );
  g_cs_diff->SetTitle("graph_cs_diff");

  TCanvas *c1 = new TCanvas("c1","c1",900,900);
  c1->cd(0);
  g_cs_nom->SetMarkerStyle(22);
  g_cs_nom->SetMarkerSize(2); 
  g_cs_nom->Draw("AP");

  TCanvas *c2 = new TCanvas("c2","c2",900,900);
  c2->cd(0);
  g_cs_corr->SetMarkerStyle(22);     
  g_cs_corr->SetMarkerSize(2); 
  g_cs_corr->Draw("AP");

  TCanvas *c3 = new TCanvas("c3","c3",900,900);
  c3->cd(0);
  g_cs_diff->GetHistogram()->SetMaximum(1.05);
  g_cs_diff->GetHistogram()->SetMinimum(0.93);
  g_cs_diff->SetMarkerStyle(22);     
  g_cs_diff->SetMarkerSize(2); 
  g_cs_diff->Draw("AP");


  f_test->Write();
  f_test->Save();
  f_test->Close();



  return 0;
}
