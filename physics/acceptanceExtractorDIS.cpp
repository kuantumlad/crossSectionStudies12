#include <vector>
#include <string>
#include <map>
#include <iostream>
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"

const static int BUFFER = 50;
double Ebeam = 7.546;//2193;

int ele_sector[BUFFER];
TLorentzVector ele[BUFFER]; 
TLorentzVector prot[BUFFER];
TLorentzVector mc_ele[BUFFER];
TLorentzVector mc_prot[BUFFER];

int select_ele;

TH1F* AcceptanceRatio(TH1F*,TH1F*);
TH2F* AcceptanceRatio2D(TH2F*,TH2F*);
TH2F* WriteAcceptanceRatio2D(TH2F*,TH2F*);

TH1F* AcceptanceRatio(TH1F *h_rec, TH1F *h_gen ){


  int rec_bins = h_rec->GetNbinsX();
  int mc_bins = h_gen->GetNbinsX();
  if( mc_bins != rec_bins ) std::cout << " ERROR WITH BIN NUMBERS CAN NO PROCEED " << std::endl;

  std::cout << " rec bin " << rec_bins << " mc bins " << mc_bins << std::endl;
  
  std::cout << " min bin center " << h_rec->GetBinCenter(1) << std::endl;
  std::cout << " >> bin width/2 " << (h_rec->GetBinCenter(1) + h_rec->GetBinCenter(2))/2.0  << std::endl;
  std::cout << " root bin width " << h_rec->GetBinWidth(1) << std::endl;
  double bin_width =  h_rec->GetBinWidth(1);
  double max_range = h_rec->GetBinCenter(h_rec->GetNbinsX()) + bin_width/2.0;
  double min_range = h_rec->GetBinCenter(1) - bin_width/2.0;

  //create acceptance histrogram over interval min and max range with rec bins (assuming same as mc bins)
  std::cout << " creating acceptance histogram for " << h_rec->GetTitle() << std::endl;
  std::cout << " bin range from " << min_range << " to " << max_range << std::endl;
  TH1F *h_accp = new TH1F(Form("h_acc_%s",h_rec->GetTitle()),Form("h_acc_%s",h_rec->GetTitle()), rec_bins, min_range, max_range);


  for( int b = 1; b < rec_bins; b++ ){
        
    if(  h_rec->GetBinCenter(b) > 25.0 ) continue;

    double accp_ratio = h_rec->GetBinContent(b)/h_gen->GetBinContent(b);    
    //    double accp_ratio_err = accp_ratio * sqrt( h_rec->GetBinContent(b)/pow(h_gen->GetBinContent(b),2)  + pow(h_rec->GetBinContent(b),2)/pow(h_gen->GetBinContent(b),3));
    
    double accp_ratio_err =1.0/h_gen->GetBinContent(b)  * sqrt( h_rec->GetBinContent(b)*h_rec->GetBinContent(b) / h_gen->GetBinContent(b) + h_rec->GetBinContent(b));
    if(h_gen->GetBinContent(b) == 0 ) accp_ratio = 0.0;
    if(h_gen->GetBinContent(b) == 0 || h_rec->GetBinContent(b) == 0 ) accp_ratio_err = 0.0; 
    
    std::cout << " bin " << b << " bin center " << h_rec->GetBinCenter(b) << " rec " << h_rec->GetBinContent(b) << " gen " << h_gen->GetBinContent(b) << std::endl;
    std::cout << " bin " << b << " bin error " <<  accp_ratio_err << std::endl;
    h_accp->SetBinContent(b,accp_ratio);
    h_accp->SetBinError(b,accp_ratio_err);
			  
  }

  return h_accp;

}

TGraphErrors* AcceptancePerSector(TH2F *h_rec_in, TH2F *h_gen_in){

  TH1F *h_rec = (TH1F*)h_rec_in->ProjectionY();
  
  int n_bins_filled = -1;
  for( int bb = 0; bb < h_rec_in->GetNbinsX(); bb++ ){    
    int current_content = h_rec_in->GetBinContent(bb,20); // scan over x bins at y_bin 20
    if( current_content > 0 ){
      n_bins_filled++;
    }
  }

  std::cout<< " projecting bins from 0 to " << n_bins_filled << std::endl;
  // get same number of phi bins as that in the reconstructed histogram
  TH1F *h_gen = (TH1F*)h_gen_in->ProjectionY(h_gen_in->GetTitle(), 0, n_bins_filled);
  

  // calculate the acceptance integrated over a sector
  std::vector<double> v_theta_center;
  std::vector<double> v_accp;
  std::vector<double> v_theta_center_err;
  std::vector<double> v_accp_err;

  int rec_bins = h_rec->GetNbinsX();
  int mc_bins = h_gen->GetNbinsX();
  if( mc_bins != rec_bins ) std::cout << " ERROR WITH BIN NUMBERS CAN NO PROCEED " << std::endl;

  std::cout << " rec bin " << rec_bins << " mc bins " << mc_bins << std::endl;
  
  std::cout << " min bin center " << h_rec->GetBinCenter(1) << std::endl;
  std::cout << " >> bin width/2 " << (h_rec->GetBinCenter(1) + h_rec->GetBinCenter(2))/2.0  << std::endl;
  std::cout << " root bin width " << h_rec->GetBinWidth(1) << std::endl;
  double bin_width =  h_rec->GetBinWidth(1);
  double max_range = h_rec->GetBinCenter(h_rec->GetNbinsX()) + bin_width/2.0;
  double min_range = h_rec->GetBinCenter(1) - bin_width/2.0;

  //create acceptance histrogram over interval min and max range with rec bins (assuming same as mc bins)
  std::cout << " creating acceptance histogram for " << h_rec->GetTitle() << std::endl;
  std::cout << " bin range from " << min_range << " to " << max_range << std::endl;
  TH1F *h_accp = new TH1F(Form("h_acc_%s",h_rec->GetTitle()),Form("h_acc_%s",h_rec->GetTitle()), rec_bins, min_range, max_range);


  for( int b = 1; b < rec_bins; b++ ){
        
    double theta_center = h_rec->GetBinCenter(b);
    if( theta_center < 9.0 || theta_center > 25.0 ) continue; 
    double accp_ratio = h_rec->GetBinContent(b)/h_gen->GetBinContent(b);    
    //double accp_ratio_err = accp_ratio * sqrt( h_rec->GetBinContent(b)/pow(h_gen->GetBinContent(b),2)  + pow(h_rec->GetBinContent(b),2)/pow(h_gen->GetBinContent(b),3));
    double accp_ratio_err =1.0/h_gen->GetBinContent(b)  * sqrt( h_rec->GetBinContent(b)*h_rec->GetBinContent(b) / h_gen->GetBinContent(b) + h_rec->GetBinContent(b));
    if(h_gen->GetBinContent(b) == 0 ) accp_ratio = 0.0;
    if(h_gen->GetBinContent(b) == 0 || h_rec->GetBinContent(b) == 0 ) accp_ratio_err = 0.0; 
    
    std::cout << " bin " << b << " bin center " << h_rec->GetBinCenter(b) << " rec " << h_rec->GetBinContent(b) << " gen " << h_gen->GetBinContent(b) << std::endl;
    std::cout << " bin " << b << " bin error " <<  accp_ratio_err << std::endl;
    h_accp->SetBinContent(b,accp_ratio);
    h_accp->SetBinError(b,accp_ratio_err);

    v_theta_center.push_back(theta_center);
    v_accp.push_back(accp_ratio);
    v_accp_err.push_back(accp_ratio_err);
    v_theta_center_err.push_back(0.0);
			  
  }

  TGraphErrors *graph = new TGraphErrors(v_theta_center.size(), &(v_theta_center[0]), &(v_accp[0]), &(v_theta_center_err[0]), &(v_accp_err[0]) );

  return graph;
}



TH2F* AcceptanceRatio2D(TH2F* h_rec, TH2F* h_gen){


  if( h_rec->GetNbinsX() != h_gen->GetNbinsX() || h_rec->GetNbinsY() != h_gen->GetNbinsY() ){
    std::cout << " bins are not equal in REC and GEN histograms " << std::endl;
  }

  double max_rangeX = h_rec->ProjectionX()->GetBinCenter(h_rec->ProjectionX()->GetNbinsX()) + (h_rec->ProjectionX()->GetBinCenter(1) + h_rec->ProjectionX()->GetBinCenter(2))/2.0;
  double min_rangeX = h_rec->ProjectionX()->GetBinCenter(1) - (h_rec->ProjectionX()->GetBinCenter(1) + h_rec->ProjectionX()->GetBinCenter(2))/2.0;

  double max_rangeY = h_rec->ProjectionY()->GetBinCenter(h_rec->ProjectionY()->GetNbinsY()) + (h_rec->ProjectionY()->GetBinCenter(1) + h_rec->ProjectionY()->GetBinCenter(2))/2.0;
  double min_rangeY = h_rec->ProjectionY()->GetBinCenter(1) - (h_rec->ProjectionY()->GetBinCenter(1) + h_rec->ProjectionY()->GetBinCenter(2))/2.0;


  std::cout << " Rec Number of X bins " << h_rec->GetNbinsX() << " Number of Y bins " << h_rec->GetNbinsY() << std::endl;
  std::cout << " Gen Number of X bins " << h_gen->GetNbinsX() << " Number of Y bins " << h_gen->GetNbinsY() << std::endl;
  std::cout <<  " X range " << min_rangeX << " " << max_rangeX << std::endl;
  std::cout <<  " Y range " << min_rangeY << " " << max_rangeY << std::endl;

  TH2F *h2_accp = new TH2F(Form("h2_accp_%s",h_rec->GetTitle()),Form("h2_accp_%s",h_rec->GetTitle()), h_rec->GetNbinsX(), min_rangeX, max_rangeX, h_rec->GetNbinsY(), min_rangeY, max_rangeY);

  
  for( int bx = 1; bx < h_rec->GetNbinsX(); bx++ ){
    for( int by = 1; by < h_rec->GetNbinsY(); by++ ){
      double accp_ratio = h_rec->GetBinContent(bx,by)/h_gen->GetBinContent(bx,by);
      if( h_gen->GetBinContent(bx,by) == 0 ) accp_ratio = 0.0;
      
      std::cout << " bx " << bx << " by " << by << " " <<  h_rec->GetBinContent(bx,by) << " " << h_gen->GetBinContent(bx,by) << " "  << accp_ratio << std::endl;
      h2_accp->SetBinContent(bx,by,accp_ratio);
    }
  }
  
			   
  return h2_accp;

}

//// write acceptance values to file 
void WriteAcceptanceRatio2D(TH2F* h_rec, TH2F* h_gen, int run, const char* field){


  if( h_rec->GetNbinsX() != h_gen->GetNbinsX() || h_rec->GetNbinsY() != h_gen->GetNbinsY() ){
    std::cout << " bins are not equal in REC and GEN histograms " << std::endl;
  }

  double max_rangeX = h_rec->ProjectionX()->GetBinCenter(h_rec->ProjectionX()->GetNbinsX()) + (h_rec->ProjectionX()->GetBinCenter(1) + h_rec->ProjectionX()->GetBinCenter(2))/2.0;
  double min_rangeX = h_rec->ProjectionX()->GetBinCenter(1) - (h_rec->ProjectionX()->GetBinCenter(1) + h_rec->ProjectionX()->GetBinCenter(2))/2.0;

  double max_rangeY = h_rec->ProjectionY()->GetBinCenter(h_rec->ProjectionY()->GetNbinsY()) + (h_rec->ProjectionY()->GetBinCenter(1) + h_rec->ProjectionY()->GetBinCenter(2))/2.0;
  double min_rangeY = h_rec->ProjectionY()->GetBinCenter(1) - (h_rec->ProjectionY()->GetBinCenter(1) + h_rec->ProjectionY()->GetBinCenter(2))/2.0;

  std::cout << " Rec Number of X bins " << h_rec->GetNbinsX() << " Number of Y bins " << h_rec->GetNbinsY() << std::endl;
  std::cout << " Gen Number of X bins " << h_gen->GetNbinsX() << " Number of Y bins " << h_gen->GetNbinsY() << std::endl;
  std::cout <<  " X range " << min_rangeX << " " << max_rangeX << std::endl;
  std::cout <<  " Y range " << min_rangeY << " " << max_rangeY << std::endl;

  TH2F *h2_accp = new TH2F(Form("h2_accp_%s",h_rec->GetTitle()),Form("h2_accp_%s",h_rec->GetTitle()), h_rec->GetNbinsX(), min_rangeX, max_rangeX, h_rec->GetNbinsY(), min_rangeY, max_rangeY);

  

 //h_rec->GetNbinsY(), min_rangeX, max_rangeX, min_rangeY, max_rangeY);

  ofstream outputAcceptance;
  ofstream outputAcceptanceError;
  std::string parentDirectory = "/w/hallb-scifs17exp/clas12/bclary/CLAS12/electron_studies/physics/parameters/";
  std::string f_out_name = parentDirectory+"dis_theta_phi_acceptance_"+std::to_string(run)+"_f"+field+".txt";
  std::string f_out_name_err = parentDirectory+"dis_theta_phi_acceptance_error_"+std::to_string(run)+"_f"+field+".txt";
  std::cout << " Creating Acceptance output file " << f_out_name_err << std::endl;
  outputAcceptance.open(f_out_name);
  outputAcceptanceError.open(f_out_name_err);
  
  //TH2F *h_rec_clone = (TH2F*)h_rec->Clone();
  //h_rec_clone->Divide(h_gen);

  for( int bx = 1; bx < h_rec->GetNbinsX(); bx++ ){
    for( int by = 1; by < h_rec->GetNbinsY(); by++ ){
      
      double accp_ratio = h_rec->GetBinContent(bx,by)/h_gen->GetBinContent(bx,by);
      double accp_ratio_err =1.0/h_gen->GetBinContent(bx,by)  * sqrt( h_rec->GetBinContent(bx,by)*h_rec->GetBinContent(bx,by) / h_gen->GetBinContent(bx,by) + h_rec->GetBinContent(bx,by));  // accp_ratio * sqrt( h_rec->GetBinContent(bx,by)/pow(h_gen->GetBinContent(bx,by),2)  + pow(h_rec->GetBinContent(bx,by),2)/pow(h_gen->GetBinContent(bx,by),3));

      if(h_gen->GetBinContent(bx,by) == 0 || h_rec->GetBinContent(bx,by) == 0 ) accp_ratio_err = 0.0; 
      if( h_gen->GetBinContent(bx,by) == 0 ) accp_ratio = 0.0;
      
      std::cout <<  accp_ratio << " ";
      outputAcceptance << accp_ratio << " ";
      outputAcceptanceError << accp_ratio_err << " ";
    }
    std::cout << "\n";
    outputAcceptance << "\n";
    outputAcceptanceError << "\n";

  }  			   
}


int acceptanceExtractorDIS(const char* inFileData, int run, const char* field_config ){ // const char* outfile, int run ){

 
  TFile *fData; 
  TFile *fMC; 
  Char_t tmpstr[80];
  Double_t fraction;

  fData = new TFile(inFileData,"");   // Input File
  if(fData->IsZombie() ){
    cout<<"Input file doesn't exist!" << endl;
    cout<<"Exit program" << endl;
    return 0;
  }

  cout << "Reading from File: " << inFileData << endl;

  // get relevant histograms to determine efficiency or geometrical acceptance

  TH2F *h_rc_all_el_w_vs_q2 = (TH2F*)fData->Get("particle_histograms_selected/hist_electron_q2_w");
  TH2F *h_mc_all_el_w_vs_q2 = (TH2F*)fData->Get("mc/hist_mc_all_electron_q2_vs_w");

  
  int q2_bins = h_rc_all_el_w_vs_q2->GetYaxis()->GetNbins();
  std::cout << " q2 bins " << q2_bins <<std::endl;
  int boxes = -1;
  int rows = 7;
  int col = 5;
  std::cout << " q2_bins % rows " << q2_bins % rows << std::endl;
  if( q2_bins % col == 0 ) rows = q2_bins/col;
  if( q2_bins % col != 0 ) rows = (q2_bins - (q2_bins%col))/col + 1;   
  std::cout << " col " << col << " rows " << rows << std::endl;   

  TCanvas *c1 = new TCanvas("c1","c1",900,900);
  c1->Divide(col,rows);
  
  for(int q2b = 1; q2b < q2_bins ; q2b++ ){
    c1->cd(q2b);
    
    //for each q2 bin make a W plot for accpetance
    std::vector<double > w_center;
    std::vector<double> q2_center;
    std::vector<double> q2w_acc;
    std::vector<double> w_err;
    std::vector<double > accp_err;
    

    int w_bins = h_rc_all_el_w_vs_q2->GetNbinsX();
    for( int wb = 1; wb < w_bins; wb++ ){
      double q2w_accp = (double)h_rc_all_el_w_vs_q2->GetBinContent(q2b,wb)/(double)h_mc_all_el_w_vs_q2->GetBinContent(q2b,wb);
      std::cout << " q2 b " << q2b << " content " << h_rc_all_el_w_vs_q2->GetBinContent(q2b,wb) << " w bin " << wb << " content " << h_mc_all_el_w_vs_q2->GetBinContent(q2b,wb) << " ratio " << q2w_accp << std::endl;
      double w_cent = h_rc_all_el_w_vs_q2->ProjectionX(Form("proj_q2b_%d",q2b),q2b,q2b)->GetBinCenter(wb);
      std::cout << " w center " << w_cent << std::endl;
      

      if((double)h_mc_all_el_w_vs_q2->GetBinContent(q2b,wb) == 0 ) q2w_accp = 0.0;

      double werr = 0.0;
      double accperr = 0.0;
      
      w_center.push_back(w_cent);
      q2w_acc.push_back(q2w_accp);
      w_err.push_back(werr);
      accp_err.push_back(accperr);
    }
    TGraphErrors *g_q2w_accp = new TGraphErrors(w_center.size(), &(w_center[0]),&(q2w_acc[0]), &(w_err[0]), &(accp_err[0]) );    
    g_q2w_accp->SetTitle(Form("W Acceptance Q^{2} bin %d",q2b));
    g_q2w_accp->SetMarkerStyle(8);
    g_q2w_accp->SetMarkerSize(0.75);    
    g_q2w_accp->GetXaxis()->SetTitle("W (GeV)");
    g_q2w_accp->GetYaxis()->SetTitle("Accp");
    g_q2w_accp->Draw("AP");
  }

  
  ///  WriteAcceptanceRatio2D(h_rc_all_el_theta_vs_phi, h_mc_all_el_theta_vs_phi, run, field_config);



  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // acceptance per sector
  // using the histograms from the selector_dis

  TH2F *h_mc_all_el_w_vs_q2_per_sect = (TH2F*)fData->Get("mc/hist_mc_all_electron_q2_vs_w");

  for( int ss =0; ss < 6; ss ++ ){
    TCanvas *c_s = new TCanvas(Form("cs%d",ss),Form("cs%d",ss),900,900);
    c_s->Divide(col,rows);
    

    TH2F *h_rc_all_el_w_vs_q2_per_sect = (TH2F*)fData->Get(Form("particle_histograms_selected/hist_electron_q2_w_s%d",ss+1));
    
    for(int q2b = 1; q2b < q2_bins ; q2b++ ){
      c_s->cd(q2b);
    
      //for each q2 bin make a W plot for accpetance
      std::vector<double > w_center;
      std::vector<double> q2_center;
      std::vector<double> q2w_acc;
      std::vector<double> w_err;
      std::vector<double > accp_err;
    

      int w_bins = h_rc_all_el_w_vs_q2_per_sect->GetNbinsX();
      for( int wb = 1; wb < w_bins; wb++ ){
	double q2w_accp = (double)h_rc_all_el_w_vs_q2_per_sect->GetBinContent(q2b,wb)/(double)h_mc_all_el_w_vs_q2_per_sect->GetBinContent(q2b,wb); // divide by 6 for sectors
	//std::cout << " q2 b " << q2b << " content " << h_rc_all_el_w_vs_q2_per_sect->GetBinContent(q2b,wb) << " w bin " << wb << " content " << h_mc_all_el_w_vs_q2_per_sect->GetBinContent(q2b,wb) << " ratio " << q2w_accp << std::endl;
	double w_cent = h_rc_all_el_w_vs_q2_per_sect->ProjectionX(Form("proj_q2b_%d_s%d",q2b,ss),q2b,q2b)->GetBinCenter(wb);
	//std::cout << " w center " << w_cent << std::endl;
      

	if((double)h_mc_all_el_w_vs_q2_per_sect->GetBinContent(q2b,wb) == 0 ) q2w_accp = 0.0;

	double werr = 0.0;
	double accperr = 0.0;
      
	w_center.push_back(w_cent);
	q2w_acc.push_back(q2w_accp * 6.0);
	w_err.push_back(werr);
	accp_err.push_back(accperr);
      }
      TGraphErrors *g_q2w_accp = new TGraphErrors(w_center.size(), &(w_center[0]),&(q2w_acc[0]), &(w_err[0]), &(accp_err[0]) );    
      g_q2w_accp->SetTitle(Form("W Acceptance Q^{2} bin %d Sector %d",q2b,ss));
      g_q2w_accp->SetMarkerStyle(8);
      g_q2w_accp->SetMarkerSize(0.75);    
      g_q2w_accp->GetXaxis()->SetTitle("W (GeV)");
      g_q2w_accp->GetYaxis()->SetTitle("Accp");
      g_q2w_accp->Draw("AP");
    }
  }
 
  
 
  return 0;  



}


