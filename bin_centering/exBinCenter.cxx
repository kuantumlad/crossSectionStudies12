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




double testFunction1(double x){
  return 1.0/pow(x,4);
}

int exBinCenter(){

  TFile *f_test = new TFile("test_bin_center.root","RECREATE");
  

  double bin_min = 8.1;
  double bin_max = 28.1;
  double n_bin_size = 0.010;
  double n_bin_size_coarse = 0.010 * 10;
  int n_bins = (bin_max - bin_min )/n_bin_size;
  int n_bins_coarse = (bin_max - bin_min )/n_bin_size_coarse;
  std::cout << " BIN MIN " << bin_min << " BIN MAX " << bin_max << " N BINS " << n_bins << std::endl;
  
  TF1 *f_cs = new TF1("f_cs","1.0/(x^4)",bin_min,bin_max);
  TH1D *h_test_cs = new TH1D("h_test_cs","h_test_cs",(int)n_bins, bin_min, bin_max );
   
  for( int b = 1; b < n_bins; b++ ){
    double test_cs = 1.0/testFunction1(h_fine_cs->GetBinCenter(b));
    h_fine_cs->SetBinContent(b,test_cs);   
    
  }

  std::vector< double> v_nom_
  for( int b=1; b < n_bins; b++ ){
  
    double bin_left = h_fine_cs->GetBinLowEdge(b);
    double bin_half = h_fine_cs->GetBinWidth(b)/2.0;
    double bin_right = h_fine_cs->GetBinCenter(b) + bin_half;
    
    double funct_left = testFunction1(bin_left);
    double funct_right = testFunction1(bin_right);
        
    double funct_avg = (funct_left + funct_right)/2.0;
    
    std::cout << " LOWER BIN X " << bin_left << " FUNCTION " << funct_left << std::endl;
    std::cout << " UPPER BIN " << bin_right << " FUNCTION " << funct_right << std::endl;
    
    std::cout << " FUNCT CENTER " << h_fine_cs->GetBinCenter(b) << " FUNCT AVG CENTER " << funct_avg << std::endl;   
  
  }
  

  f_test->Write();
  f_test->Save();
  f_test->Close();



  return 0;
}
