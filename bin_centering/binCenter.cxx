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
#include <TRandom3.h>
#include <TCanvas.h>

TH1D* GetHistoModel(std::string model, int n_bins, double min, double max, int n_bin_start, double min_start, double temp_delta_theta){

  int model_theta_bin=n_bin_start;
  double model_theta=min_start;
  double delta_theta=temp_delta_theta;
  double in_model_cs;
  double in_model_theta;
  string line;
  ifstream readFromEModel(model);//parentDirectory+"w_cut_limits_run"+std::to_string(run)+".txt");

  TH1D *h_model = new TH1D(Form("h_model_%d",n_bins),Form("h_model_%d",n_bins),n_bins,min,max);

  if( readFromEModel.is_open() ){
    std::cout << " OPENED FILES " << std::endl;	
    while(readFromEModel >> in_model_theta >> in_model_cs ) {
      if ( in_model_theta < 9 && in_model_theta < 40 ) {
	//std::cout << " dont plot model theta " << std::endl;
      }
      else{
	h_model->SetBinContent(model_theta_bin+1,in_model_cs);
      }
      //std::cout << " MODEL BIN " << model_theta_bin << " CALC THETA " << model_theta << " MODEL THETA " << in_model_theta << " HISTOGRAM THETA CENTER " << h_model->GetBinCenter(model_theta_bin+1) <<  " >> MODEL VALUE " << in_model_cs << std::endl;
      model_theta+=delta_theta;
      model_theta_bin+=1;
    
    }
  }
  return h_model;
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

TGraphErrors* GetBinCorrection( TH1D *h_model,  TH1D* h_temp, int n_bin_group, int n_bins_per_group, double min, double max, double min_coarse, double max_coarse){

  // h_temp is finely made histograms
  std::vector<double> v_x;
  std::vector<double> v_y;
  std::vector<double> v_x_err;
  std::vector<double> v_y_err;

  int start=min;
  int end=start+n_bins_per_group;
  ///std::cout << " GET BIN CORR FOR " << h_model->GetTitle() << " coarse " << h_temp->GetTitle() << " n bin group " << n_bin_group << " n bins per group " << n_bins_per_group << " min " <<  min << " max " << max << std::endl;
    
  int global_sub_bin_counter=0; //count the number of small bins
  for( int bg = min_coarse; bg < n_bin_group-10; bg++ ){
    // std::cout << " start " << start << " start bin center  " << h_temp->GetBinCenter(start) << " end " << end << " end bin center " <<h_temp->GetBinCenter(end) << std::endl;
    double sum=0;
    double avg=0;
    for( int bi = start; bi < end; bi++ ){
      double bin_center_content=h_temp->GetBinContent(bi);
      sum=sum+bin_center_content;
      
      //std::cout <<" bin " << bi << " bin cnter " << h_temp->GetBinCenter(bi) << " bin center value " <<  bin_center_content << std::endl;    
      //std::cout << global_sub_bin_counter << " & " << h_temp->GetBinCenter(bi) << " & " <<  bin_center_content << " \\\\ "  << std::endl;    
      global_sub_bin_counter++;
    }

    int to_center_shift = TMath::Floor(((double)n_bins_per_group - 1.0)/2.0); 
    int center_bin = bg;
    double bin_center_cs = h_model->GetBinContent(bg+1);
    //std::cout << " start_bin " << start << " center bin  " <<  center_bin << " coarse model value " << bin_center_cs <<std::endl;
    avg=sum/(double)n_bins_per_group;
    double bin_corr = bin_center_cs/avg;
    double bin_x = h_model->GetBinCenter(bg+1);

    //std::cout << " bin group index " << bg << " x " << bin_x << " bin_corr " << bin_corr << std::endl;
    std::cout << bin_x << " & " << bin_corr << " \\\\" << std::endl;

    start+=n_bins_per_group;
    end+=n_bins_per_group;
    
    v_x.push_back( bin_x );
    v_y.push_back( bin_corr);
    v_x_err.push_back(0.0);
    v_y_err.push_back(0.0);

    double thetaStep=1.0;
    //double model_value=elas_(5.498,(bg-1)*thetaStep+min);//, 5.0/865.0, 1.1);
    //std::cout << " elastic model value: " << model_value << std::endl;
    
  }
  
  TGraphErrors *g_bin_corr = new TGraphErrors(v_x.size(), &(v_x[0]), &(v_y[0]), &(v_x_err[0]), &(v_y_err[0]) );
  g_bin_corr->SetTitle(Form("g_bin_corr_nbins%d",(int)v_x.size()));
  return g_bin_corr;
}

int binCenter(const char* inModel, const char* inModelCoarse, const char* fOutName){

  TFile *fOut = new TFile(fOutName,"RECREATE");

  std::string model = inModel;
  std::string modelCoarse = inModelCoarse;

  TH1D *h_model_coarse = GetHistoModel(inModelCoarse, 47, 0.0, 47.0, 5, 5.5, 1.0 );
  TH1D *h_model = GetHistoModel(model, 402, 0.0, 40.2, 50, 5.1, 0.1);
  
  int n_bin_group= h_model_coarse->GetNbinsX();
  int n_bins_per_group = 10;
  int min = 91;
  int max = 500;
  int min_coarse = 9;
  int max_coarse = 50;
  TGraphErrors *g_corr = GetBinCorrection( h_model_coarse, h_model, n_bin_group, n_bins_per_group, min, max, min_coarse, max_coarse);

  TCanvas *c1 = new TCanvas("c1","c1",900,900);
  c1->cd(0);
  h_model->Draw();
  
  TCanvas *c2 = new TCanvas("c2","c2",900,900);
  c2->cd(0);
  g_corr->SetMarkerStyle(21);
  g_corr->SetMarkerSize(0.6);
  g_corr->SetTitle("Bin Correction");
  g_corr->GetXaxis()->SetTitle("#theta (deg)");
  g_corr->GetYaxis()->SetTitle("Bin Corr.");
  g_corr->Draw("AP");


  TCanvas *c3 = new TCanvas("c3","c3",900,900);
  c3->cd(0);
  gStyle->SetOptStat(0);
  h_model_coarse->SetLineColor(kBlue+2);
  h_model->SetLineColor(kRed);
  //h_model->SetLineStyle(10);
  h_model->SetLineWidth(2);
  h_model_coarse->SetLineWidth(2);
  h_model->SetAxisRange(5.0,9.0,"Y");
  h_model_coarse->SetAxisRange(5.0,9.0,"Y");
  h_model->SetAxisRange(9.9,11.1,"X");
  h_model_coarse->SetAxisRange(9.9,11.1,"X");
  h_model->GetXaxis()->SetTitle("#theta (deg)");
  h_model->GetYaxis()->SetTitle("#sigma (mb)");
  h_model->SetTitleOffset(1.4);
  h_model->SetTitle("");
  h_model->Draw("same");
  h_model_coarse->Draw("same");

  TPad *grid = new TPad("grid","",0,0,1,1); 
  grid->Draw();
  grid->cd();
  grid->SetGrid();
  grid->SetFillStyle(4000); 
  grid->SetFrameFillStyle(0);

  
  TH2 *hgrid = new TH2C("hgrid","",12, 9.9, 11.1, 4, 5.0, 9.0);   
  hgrid->Draw();
  hgrid->GetXaxis()->SetNdivisions(12);
  hgrid->GetYaxis()->SetNdivisions(16);
  hgrid->GetYaxis()->SetLabelOffset(999.); 
  hgrid->GetXaxis()->SetLabelOffset(999.); 


  c3->SaveAs("coarse_fine_bin_corr_elastic_10bptheta.pdf");



  fOut->Write();

  return 0;
}
