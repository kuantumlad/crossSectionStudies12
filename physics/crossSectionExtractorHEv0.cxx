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

std::map< int , std::vector<double> > GetAcceptancePerSector(const char *fName ){

  std::string parentDirectory = "/w/hallb-scifs17exp/clas12/bclary/CLAS12/electron_studies/physics/parameters/";


  double accp;
  int tb;
  int theta_bin = 0;
  std::map< int, std::vector<double> > s_accp;

  string line;
  ifstream readFromAcceptance(parentDirectory+fName);
  if( readFromAcceptance.is_open() ){
    std::cout << " OPENED FILES " << std::endl;	
    for( std::string line; getline( readFromAcceptance, line); ){
      //std::cout << line << std::endl;
      
      std::stringstream stream(line);
      std::string accp_val;      
      std::vector<double> theta_accp_val;
      while( stream >> tb >>  accp_val ){       
	theta_accp_val.push_back(std::atof(accp_val.c_str()));
	std::cout << " tb "<< tb << " accp_val " << std::atof(accp_val.c_str());	
	theta_bin = tb;
      }
      s_accp[theta_bin] = theta_accp_val;
    }
  }

  return s_accp;


}

std::map< int , std::vector<double> > GetWRangePerSector( const char* fName ){

  std::map<int , std::vector<double> > s_cutrange;

  std::string parentDirectory = "/w/hallb-scifs17exp/clas12/bclary/CLAS12/electron_studies/physics/parameters/";

  int tb;
  double mean;
  double sig;
  int theta_bin = 0;

  string line;
  ifstream readFromFile(parentDirectory+fName);
  if( readFromFile.is_open() ){
    std::cout << " OPENED FILES " << std::endl;	
    for( std::string line; getline( readFromFile, line); ){
      std::cout << line << std::endl;
      
      std::stringstream stream(line);
      std::string accp_val;
      std::string mean_val;
      std::string sig_val;
      std::string tb_val;
      //std::cout << " GETTING VALUES FOR TB " << theta_ bi<< std::endl;
      std::vector<double> cut_val;
      while( stream >> tb >> mean_val >> sig_val ){       
	cut_val.push_back(std::atof(mean_val.c_str()));
	cut_val.push_back(std::atof(sig_val.c_str()));

	std::cout << " tb "<< tb << " mean  " << std::atof(mean_val.c_str()) << std::endl;
	std::cout << " tb "<< tb << " sig  " << std::atof(sig_val.c_str()) << std::endl;
	theta_bin = tb;
      }
      s_cutrange[theta_bin] = cut_val;
    }
  }

  return s_cutrange;


}


void GetCrossSectionPerSector( TH2F *hin, std::map< int , std::vector<double> > m_accp, std::map< int , std::vector<double> > m_range, double lum_factor_ug ){

  int n_bins_w = hin->GetNbinsX();    
  int n_bins_theta = hin->GetNbinsY();    
  std::cout << " w bins " << n_bins_w << " theta bins " << n_bins_theta <<std::endl;
  
  for( int tb = 0; tb < n_bins_theta; tb++ ){
      
    if ( tb < 40 && tb > 4 ){

	
      TH1D *h_w_per_theta = new TH1D(Form("%s_b%d",hin->GetTitle(),tb),Form("%s_b%d",hin->GetTitle(),tb), n_bins_w ,hin->GetXaxis()->GetBinLowEdge(1), hin->GetXaxis()->GetBinUpEdge(n_bins_w) );
	double bcy = ((TAxis*)hin->GetYaxis())->GetBinCenter(tb);
      
	for( int wb = 0; wb < n_bins_w ; wb++ ){
	  double w_bin_content = hin->GetBinContent(wb,tb);
	  double bcx = ((TAxis*)hin->GetXaxis())->GetBinCenter(wb);
	  	 	 	
	  h_w_per_theta->SetBinContent(wb, w_bin_content);		
	}

	double mean = m_range[tb][0];
	double sig = m_range[tb][1];
	double upper_cut = mean + 2.0*sig;
	double lower_cut = mean - 2.0*sig;
	std::cout << " theta center " << bcy << " mean " << mean << " sig " << sig << std::endl;
	std::cout << " upper cut " << upper_cut << " lower cut " << lower_cut << std::endl;
	int b_up = h_w_per_theta->FindBin(upper_cut);
	int b_low = h_w_per_theta->FindBin(lower_cut);
	int w_integral = h_w_per_theta->Integral(b_up, b_low);       
	std::cout << " bin up " << b_up << " b_low " << b_low << " integral " << w_integral << std::endl;

	double accp = m_accp[tb][0];
	double w_bin_width = h_w_per_theta->GetBinWidth(1);
	double phi_bin_width = 60.0; // assume getting the entire sector that covers 60 deg
	double rad_corr = 1.0;
	double sin_theta = sin( h_w_per_theta->GetXaxis()->GetBinCenter(tb) * 3.14159/180.0 );

	double cs = w_integral / (w_bin_width * phi_bin_width * sin_theta * lum_factor_ug * accp * rad_corr );
	double cs_err = 0.0;
	std::cout << " >> -------------------------------------------------------- << " << std::endl;
	std::cout << " THETA BIN CENTER      " << bcy << std::endl;
	std::cout << " CROSS SECTION         " << cs << std::endl;
	std::cout << " CROSS SECTION ERR     " << cs_err << std::endl;
	std::cout << " W INTEGRAL RANGE LOW  " << lower_cut << std::endl;
	std::cout << " W INTEGRAL RANGE UP   " << upper_cut << std::endl;
	std::cout << " W INTEGRAL BIN LOW    " << b_low  << std::endl;
	std::cout << " W INTEGRAL BIN UP     " << b_up  <<  std::endl;
	std::cout << " W INTEGRAL SUM        " << w_integral << std::endl;
	std::cout << " THETA BIN WIDTH       " <<  w_bin_width << std::endl; 
	std::cout << " SIN THETA             " << sin_theta << std::endl;
	std::cout << " PHI BIN WIDTH         " << phi_bin_width << std::endl;
	std::cout << " ACCP CORR             " << accp << std::endl; 
	std::cout << " RAD CORR              " << rad_corr << std::endl;
	std::cout << " >> -------------------------------------------------------- << " << std::endl;

	
	
	
    }
  }
	



}


int crossSectionExtractorHEv0(const char* infile, int run){

  TFile *fIn = new TFile(infile,"");
  TFile *fOut = new TFile(Form("cs_test_run%dV2.root",run),"RECREATE");
  
  if( fIn->IsZombie() ){
    std::cout << " Input file "<< fIn << " doesn't exist" << std::endl;
    std::cout << " bye " << std::endl;
    return 0;
  }

  ///////////////////////////////////////////////////////////////
  // load acceptance correction and w bin range
  ///////////////////////////////////////////////////////////////
  std::map< int , std::vector<double> > m_s_accp_corr = GetAcceptancePerSector("elastic_projY_theta_acceptance_5700_fcj6p3p1_tp1sm1_elastrad.txt"); // added sector here
  std::map< int , std::vector<double> > m_s_wrange = GetWRangePerSector("whe_cut_limits_run570011_s0.txt");

  double target_density = 0.0701; // g/cm^3
  double atomic_mass_hydrogen = 1.00794; // g/mol
  double avogado_const = 6.0221367E23; // Number/mol	
  double target_length = 5.0; //cm
  double cm_to_microbarn = 1E30;
  double el_charge = 1.602177E-19; // Coulomb
  double nC_to_C = 1e-9;
  double ns_to_s = 1e-9;

  double n_el_g =  1.2782884631348472E15;// (tot_beam_charge * nC_to_C)/el_charge;
  double n_el_ug = 1.1344064006348862E15;
  double n_el = n_el_g + n_el_ug;
  double n_pr = ( ( target_length * target_density * avogado_const ) / atomic_mass_hydrogen ) ;       	
  double lum_factor_g = (n_el_g*n_pr)/cm_to_microbarn;       	
  double lum_factor_ug = (n_el_ug*n_pr)/cm_to_microbarn;       	
  double lum_factor_tot = (n_el*n_pr)/cm_to_microbarn;       	
  std::cout << " gated lumi " << lum_factor_g << std::endl;
  std::cout << " ungated lumi " << lum_factor_ug << std::endl;

  

  std::vector<TH2F*> h_el_wtheta_sect;
  for( int s = 1; s <= 6; s++ ){
    h_el_wtheta_sect.push_back( (TH2F*)fIn->Get(Form("/kinematics/h_el_wtheta_s%d_final",s) ) );
  }

  //double lumi_run = ;
  
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
    

    GetCrossSectionPerSector(h_el_wtheta_sect[s], m_s_accp_corr, m_s_wrange, lum_factor_ug);




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



