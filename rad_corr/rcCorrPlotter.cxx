#include <TCanvas.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <TGraphErrors.h>

int rcCorrPlotter(){


  ifstream infile;
  std::string line;
  infile.open("/work/clas12/bclary/CLAS12/electron_studies/rad_corr/elastic_gen_7GeV_CORRECTION_FINAL.txt");
  
  std::vector<double> bin_center;
  std::vector<double> rc;

  std::vector<double> bin_center_err;
  std::vector<double> rc_err;
  
  
  while( std::getline(infile,line) ){
    double n0;
    double n1;
    std::istringstream ss(line);

    ss >> n0 >> n1 ;
    std::cout << n0 << " " << n1  << std::endl;

    //if ( n0 > 9 ){
      bin_center.push_back(n0);
      rc.push_back(n1);

      std::cout << " bin center " << n0 << " radCorr " << n1 << std::endl;
      
      bin_center_err.push_back(0.0);
      rc_err.push_back(0.0);
      //}
  }

  TGraphErrors *g_rc = new TGraphErrors(bin_center.size(), &(bin_center[0]), &(rc[0]), &(bin_center_err[0]), &(rc_err[0]));

  TCanvas *c_rc = new TCanvas("c_rc","c_rc",900,900);
  c_rc->cd(1);
  g_rc->SetTitle("Radiative Correction");
  g_rc->GetXaxis()->SetTitle("#theta (deg)");
  g_rc->GetXaxis()->CenterTitle();
  g_rc->GetYaxis()->SetTitle("RC");
  g_rc->GetYaxis()->CenterTitle();  
  g_rc->SetMarkerStyle(8);
  g_rc->SetMarkerSize(0.75);
  g_rc->GetHistogram()->SetMaximum(5.5);
  g_rc->GetHistogram()->SetMinimum(0.0);//2);
  g_rc->Draw("AP");

  c_rc->SaveAs("g_rc7GeV_values.pdf");

  return 0;
}
