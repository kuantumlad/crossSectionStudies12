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
#include <TRandom.h>

#include "BostedElasticWrapper.h"
#include m"../odel/Bosted


double testFunction1(double x){
  return 1.0/pow(x,4);
}

TH1D* GetHistoSim(int n_bins, double min, double max){

  TRandom3 *gRandom = new TRandom3();
  

  TH1D *h_sim = new TH1D(Form("h_sim_%d",n_bins),Form("h_sim_%d",n_bins),n_bins,min,max);
  for( int b = 1; b < n_bins; b++ ){
    h_sim->SetBinContent(b, testFunction1(h_sim->GetBinCenter(b)) + gRandom->Gaus(0.0,0.000001) );
  }

  return h_sim;

}

TH1D* GetHistoGen(int n_bins, double min, double max){

  TH1D *h_gen = new TH1D(Form("h_gen_%d",n_bins),Form("h_gen_%d",n_bins),n_bins,min,max);
  //bin 0 is underflow
  for( int b = 1; b < n_bins; b++ ){
    h_gen->SetBinContent(b, testFunction1(h_gen->GetBinCenter(b)));
  }

  return h_gen;

}

TGraphErrors *ConvertHisto(TH1D* h_temp){

  std::vector<double> v_x;
  std::vector<double> v_y;
  std::vector<double> v_x_err;
  std::vector<double> v_y_err;
  
  int n_bins = h_temp->GetNbinsX();
  for(int b = 0; b < n_bins; b++ ){
    v_y.push_back( h_temp->GetBinContent(b) );
    v_x.push_back( h_temp->GetBinCenter(b) );
    v_x_err.push_back( 0.0 ); 
    v_y_err.push_back( 0.0 ); 
  }

  TGraphErrors *g_temp = new TGraphErrors(v_x.size(), &(v_x[0]), &(v_y[0]), &(v_x_err[0]), &(v_y_err[0]) );  
  return g_temp;

}

TGraphErrors* GetBinCorrection( TH1D* h_temp, TH1D* h_temp_coarse, int n_bin_group, int n_bins_per_group, double min, double max){

  // h_temp is finely made histograms
  std::vector<double> v_x;
  std::vector<double> v_y;
  std::vector<double> v_x_err;
  std::vector<double> v_y_err;

  int start=1;
  int end=n_bins_per_group+1;
    
  for( int bg = 1; bg < n_bin_group; bg++ ){
    std::cout << " start " << start << " start bin center  " << h_temp->GetBinCenter(start) << " end " << end << " end bin center " <<h_temp->GetBinCenter(end) << std::endl;
    double sum=0;
    double avg=0;
    for( int bi = start; bi < end; bi++ ){
      double bin_center_content=h_temp->GetBinContent(bi);
      sum=sum+bin_center_content;
      std::cout <<" bin " << bi << " bin cnter " << h_temp->GetBinCenter(bi) << " bin center value " <<  bin_center_content << std::endl;    
    }
    double bin_x = h_temp_coarse->GetBinCenter(bg);
    start+=n_bins_per_group;
    end+=n_bins_per_group;

    double bin_corr = sum/n_bins_per_group;
    std::cout << " bin group index " << bg << " x " << bin_x << " bin_corr " << bin_corr << std::endl;
    
    v_x.push_back( bin_x );
    v_y.push_back( bin_corr);
    v_x_err.push_back(0.0);
    v_y_err.push_back(0.0);

    double thetaStep=1.0;
    double model_value=getRadiatedValue(5.498,(bg-1)*thetaStep+min, 5.0/865.0, 1.1);
    std::cout << " elastic model value: " << model_value << std::endl;
    
  }
  TGraphErrors *g_bin_corr = new TGraphErrors(v_x.size(), &(v_x[0]), &(v_y[0]), &(v_x_err[0]), &(v_y_err[0]) );
  g_bin_corr->SetTitle(Form("g_bin_corr_nbins%d",(int)v_x.size()));
  return g_bin_corr;


}
  

int binCenter(){

  
  double bin_min = 8.0;
  double bin_max = 28.0;
  double n_bin_size = 1.0;
  double n_bin_size_fine = 0.10;// * 10;
  int n_bins = (bin_max - bin_min )/n_bin_size;
  int n_bins_fine = (bin_max - bin_min )/n_bin_size_fine;

  TH1D *h_sim = GetHistoSim(n_bin_size,bin_min, bin_max);
  TH1D *h_gen_coarse = GetHistoGen(n_bins, bin_min, bin_max);     
  
  TH1D *h_gen_fine = GetHistoGen(n_bins_fine, bin_min, bin_max);
  

  //GetBinCorrection( TH1D* h_temp, TH1D* h_temp_coarse, int n_bin_groups, int n_bins_per_group, double min, double max){   
  TGraphErrors *g_b_corr = GetBinCorrection(h_gen_fine, h_gen_coarse, n_bins, 10, bin_min, bin_max);
  g_b_corr->SetMarkerStyle(21);
  g_b_corr->SetMarkerSize(0.5);

  TGraphErrors *g_gen = ConvertHisto(h_gen_coarse);
  g_gen->SetMarkerStyle(20);
  g_gen->SetMarkerSize(0.50);
  g_gen->SetMarkerColor(kBlue);
  
  TCanvas *c1 = new TCanvas("C1","C1",900,900);
  c1->cd(0);
  gPad->SetLogy();
  g_b_corr->SetMarkerColor(kRed);
  g_b_corr->Draw("AP");
  h_gen_coarse->SetMarkerStyle(20);
  h_gen_coarse->SetMarkerSize(0.5);
  h_gen_coarse->Draw("AP+SAME");

  TH1D *h_simulated = GetHistoSim(n_bins, bin_min, bin_max);
  std::vector<double> cs_bin_corr;  
  std::vector<double> cs_bin_corr_err;  
  // get the ratio of the corr to nom value
  std::vector<double> v_x;
  std::vector<double> v_y;
  std::vector<double> v_x_err;
  std::vector<double> v_y_err;
  TH1D *h_b_corr = (TH1D*)g_b_corr->GetHistogram();
  
  for(int b = 0; b < n_bins-1; b++ ){
    double y_corr;
    double x_corr;
    g_b_corr->GetPoint(b,x_corr,y_corr);

    std::cout << " x_crr " << x_corr << "y " << y_corr <<std::endl;

    double y_nom = h_gen_coarse->GetBinContent(b+1);
    double cross_section = h_simulated->GetBinContent(b+1);
    double ratio = y_nom/y_corr;
    std::cout << "bin " << b << " ratio " << ratio << " center " << h_gen_coarse->GetBinCenter(b+1) <<std::endl;

    v_x.push_back(h_gen_coarse->GetBinCenter(b+1));
    v_y.push_back(ratio);
    cs_bin_corr.push_back( (1.0/ratio) * cross_section);
    v_x_err.push_back(0.0);
    v_y_err.push_back(0.0);  
    cs_bin_corr_err.push_back(0.0);
  }


  TGraphErrors *g_ratio = new TGraphErrors(v_x.size(),&(v_x[0]), &(v_y[0]), &(v_x_err[0]), &(v_y_err[0]) ); 
  g_ratio->SetTitle("Ratio of NOM to CORRECTED");
  
  TCanvas *c2 = new TCanvas("c2","c2",900,900);
  c2->cd(0);
  g_ratio->SetMarkerStyle(21);
  g_ratio->SetMarkerSize(0.5);
  g_ratio->SetMarkerColor(kBlue);
  g_ratio->GetHistogram()->SetMinimum(0.95);
  g_ratio->GetHistogram()->SetMaximum(1.03);
  g_ratio->Draw("AP");

  TGraphErrors *g_bin_corr_cross_section = new TGraphErrors(v_x.size(),&(v_x[0]), &(cs_bin_corr[0]), &(v_x_err[0]), &(cs_bin_corr_err[0]) ); 
  TCanvas *c3 = new TCanvas("c3","c3",900,900);
  c3->cd(0);
  gPad->SetLogy();
  g_bin_corr_cross_section->SetTitle("Bin Corrected Cross Section");
  g_bin_corr_cross_section->SetMarkerStyle(21);
  g_bin_corr_cross_section->SetMarkerSize(0.5);
  g_bin_corr_cross_section->SetMarkerColor(kBlack);
  g_bin_corr_cross_section->Draw("AP");
 
  
  return 0;
}
