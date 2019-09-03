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

TH1D* GetModel(std::string);
TGraphErrors* GetCSRatio( TH2F*, std::map< int , std::vector<double> >, std::map< int , std::vector<double> >, std::vector<double>, double, TH1D*);

std::vector<double> GetRadiativeCorr(const char* fileName){
  //theta range is from 5 to 30
  std::cout << " getting radiative corrections " << std::endl;
  std::vector< double > v_rad_corr;
  double model_theta = 5.5;
  int model_theta_bin = 6;
  double rad_corr;
  double bin_center;
  
  std::string parentDir = "/work/clas12/bclary/CLAS12/electron_studies/rad_corr/";
  string line;
  ifstream readFromRadCorr(parentDir+fileName);
  if( readFromRadCorr.is_open() ){
    std::cout << " OPENED FILES " << std::endl;	
    while(readFromRadCorr >> bin_center >> rad_corr ) {
      //if ( model_theta < 9 ) {/
      //std::cout << " dont plot model theta " << std::endl;
      //}
      //else{
      std::cout << " bin center " << bin_center << " rad corr " << rad_corr << std::endl;
      v_rad_corr.push_back(rad_corr);
      //}
      //std::cout << " MODEL THETA " << model_theta <<  " >> BIN CORRECTION " << rad_corr << std::endl;
      model_theta+=1.0;
      model_theta_bin+=1;    
    }
  }
  return v_rad_corr;
}



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


TGraphErrors* GetCrossSectionPerSector( TH2F *hin, std::map< int , std::vector<double> > m_accp, std::map< int , std::vector<double> > m_range, std::vector<double> v_rc, double lum_factor_ug ){

  int n_bins_w = hin->GetNbinsX();    
  int n_bins_theta = hin->GetNbinsY();    
  //std::cout << " w bins " << n_bins_w << " theta bins " << n_bins_theta <<std::endl;
  
  std::vector<double> v_bc;
  std::vector<double> v_cs;
  std::vector<double> v_bc_err;
  std::vector<double> v_cs_err;

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
	//std::cout << " theta center " << bcy << " mean " << mean << " sig " << sig << std::endl;
	//std::cout << " upper cut " << upper_cut << " lower cut " << lower_cut << std::endl;
	int b_up = h_w_per_theta->FindBin(upper_cut);
	int b_low = h_w_per_theta->FindBin(lower_cut);
	int w_integral = h_w_per_theta->Integral();//b_up, b_low);       
	//std::cout << " bin up " << b_up << " b_low " << b_low << " integral " << w_integral << std::endl;

	double accp = m_accp[tb][0];
	double w_bin_width = h_w_per_theta->GetBinWidth(1);
	double phi_bin_width = 2.0*3.1415; // assume getting the entire sector that covers 60 deg
	double rad_corr = v_rc[tb-1];
	double sin_theta = sin( h_w_per_theta->GetXaxis()->GetBinCenter(tb) * 3.14159/180.0 );

	double cs = 1 * w_integral / (w_bin_width * phi_bin_width * sin_theta * lum_factor_ug * accp * rad_corr );
	double csErr = sqrt(w_integral)/(w_bin_width * phi_bin_width * sin_theta * lum_factor_ug * accp * rad_corr );
	double accp_theta_bin_err = 0;
	double cs_err = sqrt( (csErr/accp)*(csErr/accp) );// + (cs*accp_theta_bin_err)*(cs*accp_theta_bin_err) );

	/*
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
	*/
	v_bc.push_back(bcy);
	v_cs.push_back(cs);
	v_bc_err.push_back(0.0);
	v_cs_err.push_back(cs_err);	
	
    }
  }
	
  TGraphErrors *g_cs = new TGraphErrors(v_bc.size(), &v_bc[0], &v_cs[0], &v_bc_err[0], &v_cs_err[0]);
  return g_cs;

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
  // load model for given beam energy
  //////////////////////////////////////////////////////////////
  TH1D *h_model_cs = GetModel("/w/hallb-scifs17exp/clas12/bclary/CLAS12/electron_studies/models/Bosted/py_elast/cs_model_elastic7.546.txt");
  TGraph *g_model = new TGraph(h_model_cs);
  g_model->SetMarkerStyle(20);
  g_model->SetMarkerColor(kRed);

  ///////////////////////////////////////////////////////////////
  // load radiative corrections for each theta bin
  ///////////////////////////////////////////////////////////////
  std::vector<double> v_rad_corr = GetRadiativeCorr("elastic_gen_7GeV_CORRECTION_FINAL.txt");

  double target_density = 0.0701; // g/cm^3
  double atomic_mass_hydrogen = 1.00794; // g/mol
  double avogado_const = 6.0221367E23; // Number/mol	
  double target_length = 5.0; //cm
  double cm_to_microbarn = 1E30;
  double el_charge = 1.602177E-19; // Coulomb
  double nC_to_C = 1e-9;
  double ns_to_s = 1e-9;

  double n_el_g =  1.2782884631348472E15;// (tot_beam_charge * nC_to_C)/el_charge;
  double n_el_ug = (14800.72308 * nC_to_C ) / el_charge; //1.1344064006348862E15;
  std::cout << " n_el_ug " << n_el_ug << std::endl;
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
 
  std::vector<TGraphErrors*> v_cs_sectors;
  std::vector<TGraphErrors*> v_csratio_sectors;
  
  for( int s  = 0; s < 6; s++ ){
    int n_bins_w = h_el_wtheta_sect[s]->GetNbinsX();    
    int n_bins_theta = h_el_wtheta_sect[s]->GetNbinsY();    

    std::cout << " w bins " << n_bins_w << " theta bins " << n_bins_theta <<std::endl;
    

    ///////////////////////////////////////////////////////////////
    // load acceptance correction and w bin range
    ///////////////////////////////////////////////////////////////
    std::string s_sector = std::to_string(s);
    std::string file_accp = "elastic_projY_theta_acceptance_s"+s_sector+"_570011_fskim4inclusive_005700_v13_yloose_gencut.txt";
    std::string file_wrange = "whe_cut_limits_run570011_s"+s_sector+".txt";
    std::map< int , std::vector<double> > m_s_accp_corr = GetAcceptancePerSector(file_accp.c_str());
    std::map< int , std::vector<double> > m_s_wrange = GetWRangePerSector(file_wrange.c_str());

    TGraphErrors *g_cs_sector = GetCrossSectionPerSector(h_el_wtheta_sect[s], m_s_accp_corr, m_s_wrange, v_rad_corr, lum_factor_ug);
    TGraphErrors *g_ratio_sector = GetCSRatio(h_el_wtheta_sect[s], m_s_accp_corr, m_s_wrange, v_rad_corr, lum_factor_ug, h_model_cs);
    
    v_cs_sectors.push_back(g_cs_sector);
    v_csratio_sectors.push_back(g_ratio_sector);

    /*
    //TCanvas *c_tb = new TCanvas(Form("c_tb_s%d",s),Form("c_tb_s%d",s),900,900);

    //c_tb->Divide(5,8);
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


	//h_w_per_theta->Draw();
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
	
	TF1*fit_poly = new TF1("fit_poly","[0] + [1]*x + [2]*cx*x",0.5, 1.3);
	fit_poly->SetParameter(0,fitter->GetParameter(3));
	fit_poly->SetParameter(1,fitter->GetParameter(4));
	fit_poly->SetParameter(2,fitter->GetParameter(5));
	fit_poly->SetLineColor(kBlue);
	fit_poly->Draw("same");	  

      }
    }
	*/  
    //c_tb->SaveAs(Form("h_el_w_per_theta_bin_sect%d_r%d.pdf",s,run));
  }

  std::vector<TMultiGraph*> v_mg_cs;
  TCanvas *c_cs = new TCanvas("cs","cs",900,900);
  c_cs->Divide(2,3);
  for( int ss = 0; ss < 6; ss++){

    c_cs->cd(ss+1);
    v_cs_sectors[ss]->SetTitle(Form("Elastic CS S%d",ss));
    v_cs_sectors[ss]->GetXaxis()->SetTitle("#theta (deg)");
    v_cs_sectors[ss]->SetMarkerStyle(20);

    v_mg_cs.push_back( new TMultiGraph() );
    v_cs_sectors[ss]->SetMarkerSize(0.5);
    v_mg_cs[ss]->Add(v_cs_sectors[ss]);
    v_mg_cs[ss]->Add(g_model);
    v_mg_cs[ss]->Draw("APE");
    v_mg_cs[ss]->GetXaxis()->SetLimits(2.0, 20.5);//SetRangeUser(0.0,26.5);
    v_mg_cs[ss]->Draw("APE");
    c_cs->Update();
  }    
  c_cs->SaveAs("cs_wmodel_7GeV_r5700.pdf");

  std::vector<TMultiGraph*> v_mg_cs_logy;
  TCanvas *c_cs_logy = new TCanvas("cs_logy","cs_logy",900,900);
  c_cs_logy->Divide(2,3);
  for( int ss = 0; ss < 6; ss++){

    c_cs_logy->cd(ss+1);
    gPad->SetLogy();
    v_cs_sectors[ss]->SetTitle(Form("Elastic CS S%d",ss));
    v_cs_sectors[ss]->GetXaxis()->SetTitle("#theta (deg)");
    v_cs_sectors[ss]->SetMarkerStyle(20);

    v_mg_cs_logy.push_back( new TMultiGraph() );
    v_cs_sectors[ss]->SetMarkerSize(0.5);
    v_mg_cs_logy[ss]->Add(v_cs_sectors[ss]);
    v_mg_cs_logy[ss]->Add(g_model);
    v_mg_cs_logy[ss]->Draw("APE");
    v_mg_cs_logy[ss]->GetXaxis()->SetLimits(2.0, 20.5);//SetRangeUser(0.0,26.5);
    v_mg_cs_logy[ss]->Draw("APE");
    c_cs_logy->Update();
  }    
  c_cs_logy->SaveAs("cs_llogy_wmodel_7GeV_r5700.pdf");
  

  TCanvas *c_csr = new TCanvas("csr","csr",900,900);
  c_csr->Divide(2,3);
  for( int ss = 0; ss < 6; ss++){

    c_csr->cd(ss+1);
    v_csratio_sectors[ss]->SetTitle(Form("Elastic CS S%d RATIO",ss));
    v_csratio_sectors[ss]->GetXaxis()->SetTitle("#theta (deg)");
    v_csratio_sectors[ss]->SetMarkerStyle(20);
    v_csratio_sectors[ss]->SetMarkerSize(0.5);
    v_csratio_sectors[ss]->Draw("AP");
    v_csratio_sectors[ss]->GetXaxis()->SetLimits(2.0, 26.5);//SetRangeUser(0.0,26.5);
    v_csratio_sectors[ss]->GetHistogram()->SetMaximum(1.90);   // along          
    v_csratio_sectors[ss]->GetHistogram()->SetMinimum(0.0);  //   Y     
    v_csratio_sectors[ss]->Draw("AP");
    c_csr->Update();
    
    TBox *b10 = new TBox(4.0, 0.9, 26.5, 1.1);
    TBox *b15 = new TBox(4.0, 0.85, 26.5, 1.15);
    TBox *b20 = new TBox(4.0, 0.80, 26.5, 1.20);

    b10->SetFillColorAlpha(kBlue-4, 0.15);
    b15->SetFillColorAlpha(kGreen-4, 0.15);
    b20->SetFillColorAlpha(kRed-4, 0.15);

    b20->Draw("same");
    b15->Draw("same");
    b10->Draw("same");

  }    
  c_csr->SaveAs("cs_ratio_7GeV_r5700.pdf");



  return 0;
}


TH1D* GetModel( std::string in_file ){
  TH1D *h_model = new TH1D("h_elastic_model","h_elastic_model", 60, 5.0, 20.0);
  double model_value;
  double model_val;
  double model_theta;//=5.5;
  int model_theta_bin=6;
  std::string model = in_file; //"model_cross_section.txt";
  string line;
  ifstream readFromEModel(model);//parentDirectory+"w_cut_limits_run"+std::to_string(run)+".txt");
  if( readFromEModel.is_open() ){
    std::cout << " OPENED FILES " << std::endl;	
    while(readFromEModel >> model_theta >> model_value ) {//std::getline (readFromWCut, line) ){      
      ////readFromEModel >> model_val;
      //if ( model_theta < 9 ) {
      //std::cout << " dont plot model theta " << std::endl;
      //}
      //else{
      
      model_theta_bin = h_model->FindBin(model_theta );
      h_model->SetBinContent(model_theta_bin,model_value);
      //}
      std::cout << " MODEL THETA BIN " << model_theta_bin <<  " MODEL THETA " << model_theta <<  " >> MODEL VALUE " << model_value << std::endl;
      model_theta+=1.0;
      model_theta_bin+=1;   
    }
  }
  return h_model;

}

TGraphErrors* GetCSRatio( TH2F *hin, std::map< int , std::vector<double> > m_accp, std::map< int , std::vector<double> > m_range, std::vector<double> v_rc, double lum_factor_ug, TH1D *h_model ){

  int n_bins_w = hin->GetNbinsX();    
  int n_bins_theta = hin->GetNbinsY();    
  std::cout << " w bins " << n_bins_w << " theta bins " << n_bins_theta <<std::endl;
  
  std::vector<double> v_bc;
  std::vector<double> v_cs;
  std::vector<double> v_bc_err;
  std::vector<double> v_cs_err;

  for( int tb = 0; tb < n_bins_theta; tb++ ){
      
    if ( tb < 40 && tb > 4 ){

	
      TH1D *h_w_per_theta = new TH1D(Form("%s_b%d_ratio",hin->GetTitle(),tb),Form("%s_b%d_ratio",hin->GetTitle(),tb), n_bins_w ,hin->GetXaxis()->GetBinLowEdge(1), hin->GetXaxis()->GetBinUpEdge(n_bins_w) );
	double bcy = ((TAxis*)hin->GetYaxis())->GetBinCenter(tb);
      
	for( int wb = 0; wb < n_bins_w ; wb++ ){
	  double w_bin_content = hin->GetBinContent(wb,tb);
	  double bcx = ((TAxis*)hin->GetXaxis())->GetBinCenter(wb);
	  	 	 	
	  h_w_per_theta->SetBinContent(wb, w_bin_content);		
	}

	double cs_model=h_model->GetBinContent(tb);

	double mean = m_range[tb][0];
	double sig = m_range[tb][1];
	double upper_cut = mean + 2.0*sig;
	double lower_cut = mean - 2.0*sig;
	//std::cout << " theta center " << bcy << " mean " << mean << " sig " << sig << std::endl;
	//std::cout << " upper cut " << upper_cut << " lower cut " << lower_cut << std::endl;
	int b_up = h_w_per_theta->FindBin(upper_cut);
	int b_low = h_w_per_theta->FindBin(lower_cut);
	int w_integral = h_w_per_theta->Integral();//b_up, b_low);       
	std::cout << " bin up " << b_up << " b_low " << b_low << " integral " << w_integral << std::endl;

	double accp = m_accp[tb][0];
	double w_bin_width = h_w_per_theta->GetBinWidth(1);
	double phi_bin_width = 2.0*3.1415; //60.0; // assume getting the entire sector that covers 60 deg
	double rad_corr = v_rc[tb-1];
	double sin_theta = sin( h_w_per_theta->GetXaxis()->GetBinCenter(tb) * 3.14159/180.0 );

	double cs = 1 * w_integral / (w_bin_width * phi_bin_width * sin_theta * lum_factor_ug * accp * rad_corr );
	double csErr = sqrt(w_integral)/(w_bin_width * phi_bin_width * sin_theta * lum_factor_ug * accp * rad_corr );
	double accp_theta_bin_err = 0;
	double cs_err = sqrt( (csErr/accp)*(csErr/accp) );// + (cs*accp_theta_bin_err)*(cs*accp_theta_bin_err) );

	std::cout << " >> -------------------------------------------------------- << " << std::endl;
	std::cout << " THETA BIN CENTER      " << bcy << std::endl;
	std::cout << " CROSS SECTION         " << cs << std::endl;
	std::cout << " CROSS SECTION ERR     " << cs_err << std::endl;
	std::cout << " CROSS SECTION MODEL   " << cs_model << std::endl;
	std::cout << " CROSS SECTION RATIO   " << cs/cs_model << std::endl;
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

	v_bc.push_back(bcy);
	v_cs.push_back(cs/cs_model);
	v_bc_err.push_back(0.0);
	v_cs_err.push_back(cs_err);	
	
    }
  }
	
  TGraphErrors *g_cs = new TGraphErrors(v_bc.size(), &v_bc[0], &v_cs[0], &v_bc_err[0], &v_cs_err[0]);
  return g_cs;

}
