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

double alpha = 0.007299270;
double proton_mass = 0.93827203;
double beam_energy = 7.546;
double pi = 3.1415926;
int n_sectors = 6;

std::vector< std::vector<double> > GetModelValues(const char* fileName){
  
  std::map< int,  std::vector< std::vector< double > > > m_model_cs;
  double model_theta = 9.5;
  int model_theta_bin = 0;
  double q2;
  double w;
  double model_cs;
  double model_err;
  
  std::vector<double> temp_q2;
  std::vector<double> temp_w;
  std::vector<double> temp_model;
  std::vector<double> temp_err;

  std::vector< std::vector<double> > q2w_cs_model;

  int q2_bin = 1;
  int w_bin = 1;
  double old_q2 = 1;
  std::string parentDir = "/work/clas12/bclary/CLAS12/electron_studies/physics/";
  string line;
  std::cout << " opening file " << fileName << std::endl;
  ifstream readFromModel(parentDir+fileName);
  if( readFromModel.is_open() ){
    std::cout << " OPENED FILES " << std::endl;	
    while(readFromModel >> q2 >> w >> model_cs >> model_err ) {

      if( q2 > old_q2 ){
	old_q2 = q2;

	//std::vector<double> reverse_model;
	//for( int ww = temp_w.size(); ww >= 0; ww-- ){	  
	//  std::cout <<  temp_w[ww] << std::endl;
	//  reverse_model.push_back(temp_model(w));
	//	}
	q2w_cs_model.push_back(temp_model);

	q2_bin++;

	temp_q2.clear();
	temp_w.clear();
	temp_model.clear();
	temp_err.clear();

      }

      temp_q2.push_back(q2);
      temp_w.push_back(w);
      temp_model.push_back(model_cs);
      temp_err.push_back(model_err);


      std::cout << " MODEL Q2 " << q2 <<  " W " << w << " CS " << model_cs << " ERR " << model_err << std::endl;
    }
  }

  double t_q2_min=0;
  int temp_q2_bin = 0;
  int temp_w_bin = 0;
  
  for( int qq = 0; qq < q2w_cs_model.size(); qq++ ){
    for( int ww = 0; ww < q2w_cs_model[qq].size(); ww++ ){
      std::cout << q2w_cs_model[qq][ww] << std::endl;
    }
  }

  return q2w_cs_model;
}

std::vector< std::vector<double> > GetModelWCenterValues(const char* fileName){
  
  std::map< int,  std::vector< std::vector< double > > > m_model_cs;
  double model_theta = 9.5;
  int model_theta_bin = 0;
  double q2;
  double w;
  double model_cs;
  double model_err;
  
  std::vector<double> temp_q2;
  std::vector<double> temp_w;
  std::vector<double> temp_model;
  std::vector<double> temp_err;

  std::vector< std::vector<double> > q2w_wcenter_model;

  int q2_bin = 1;
  int w_bin = 1;
  double old_q2 = 1;
  std::string parentDir = "/work/clas12/bclary/CLAS12/electron_studies/physics/";
  string line;
  std::cout << " opening file " << fileName << std::endl;
  ifstream readFromModel(parentDir+fileName);
  if( readFromModel.is_open() ){
    std::cout << " OPENED FILES " << std::endl;	
    while(readFromModel >> q2 >> w >> model_cs >> model_err ) {

      if( q2 > old_q2 ){
	old_q2 = q2;

	//std::vector<double> reverse_model;
	//for( int ww = temp_w.size(); ww >= 0; ww-- ){	  
	//  std::cout <<  temp_w[ww] << std::endl;
	//  reverse_model.push_back(temp_model(w));
	//	}
	q2w_wcenter_model.push_back(temp_w);

	q2_bin++;

	temp_q2.clear();
	temp_w.clear();
	temp_model.clear();
	temp_err.clear();

      }
      
      temp_q2.push_back(q2);
      temp_w.push_back(w);
      temp_model.push_back(model_cs);
      temp_err.push_back(model_err);

      std::cout << " MODEL Q2 " << q2 <<  " W " << w << " CS " << model_cs << " ERR " << model_err << std::endl;
    }
  }

  double t_q2_min=0;
  int temp_q2_bin = 0;
  int temp_w_bin = 0;
  
  //for( int qq = 0; qq < q2w_cs_model.size(); qq++ ){
  // for( int ww = 0; ww < q2w_cs_model[qq].size(); ww++ ){
  //   std::cout << q2w_cs_model[qq][ww] << std::endl;
  // }
  //}
  return q2w_wcenter_model;
}



int crossSectionExtractorDIS(const char* infile, int run, const char* config){

  TFile *fIn = new TFile(infile,"");
  TFile *fOut = new TFile(Form("cs_test_run%dV2.root",run),"RECREATE");
  
  if( fIn->IsZombie() ){
    std::cout << " Input file "<< fIn << " doesn't exist" << std::endl;
    std::cout << " bye " << std::endl;
    return 0;
  }


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

  std::vector<TH2F*> h_el_q2w_sect;
  for( int s = 1; s <= 6 ; s++ ){
    // binning  40, 0.9, 2.1, 20, 1.0, 6.0     
    h_el_q2w_sect.push_back( (TH2F*)fIn->Get(Form("/kinematics/h_el_q2w_final_s%d",s) ));
  }

  
  std::map< int, std::vector< TH1F* > > h_el_q2w_cs;
  std::map< int, std::vector< std::vector< double > > > m_wb_center;
  std::map< int, std::vector< std::vector< double > > > m_cs_value;
  std::map< int, std::vector< std::vector< double > > > m_wb_err;
  std::map< int, std::vector< std::vector< double > > > m_cs_err; 
  std::map< int , std::vector<TGraphErrors* > > m_g_cs;

  std::map< int , std::vector<TGraphErrors* > > m_g_ln;
  std::map< int , std::vector<TGraphErrors* > > m_g_neln;

  TCanvas *c0 = new TCanvas("c0","c0",900,900);
  c0->Divide(3,2);
  for( int ss = 0; ss < n_sectors; ss++ ){
    c0->cd(ss+1);
    h_el_q2w_sect[ss]->Draw("colz");

    TCanvas *c1 = new TCanvas(Form("c1_s%d",ss),Form("c1_s%d",ss),900,900);
    c1->Divide(4,5);

    std::vector<TH1F*> v_h_w_per_q2_per_sect_temp;
    std::vector<std::vector<double> > v_wb_temp;
    std::vector<std::vector<double> > v_cs_value_temp;
    std::vector<std::vector<double> > v_wb_err_temp;
    std::vector<std::vector<double> > v_cs_err_temp;

    std::vector<TGraphErrors*> v_g_cs;
    std::vector<TGraphErrors*> v_g_rate_lum_norm;
    std::vector<TGraphErrors*> v_g_rate_nel_norm;

    for( int q2b = 0; q2b <= 20; q2b++){
      c1->cd(q2b+1);
      double w_bins = (double)h_el_q2w_sect[ss]->GetXaxis()->GetNbins();
      double w_min = h_el_q2w_sect[ss]->GetXaxis()->GetBinLowEdge(1);
      double w_max = h_el_q2w_sect[ss]->GetXaxis()->GetBinUpEdge(w_bins);
      double w_width = (w_max-w_min)/w_bins;
     

      double q2_bins = (double)h_el_q2w_sect[ss]->GetYaxis()->GetNbins();
      double q2_min = h_el_q2w_sect[ss]->GetYaxis()->GetBinLowEdge(1);
      double q2_max = h_el_q2w_sect[ss]->GetYaxis()->GetBinUpEdge(q2_bins);
      double q2_width = (q2_max-q2_min)/q2_bins;
      //double q2_center = q2b*q2_width + q2_min + q2_width/2.0;
      double q2_center = q2_width*(q2b + 1.0/2.0) + q2_min;
      std::cout << " Q2 Bin Center " << q2_center <<std::endl;

      std::cout << " w min " << w_min << " w max " << w_max << " w width " << w_width << std::endl;
      std::cout << " q2 min " << q2_min << " q2 max " << q2_max << " q2 width " << q2_width << std::endl;
	
      TH1F *h_w_per_q2_per_sect = new TH1F(Form("h_w_per_q2_%d_s%d",q2b,ss), Form("h_w_per_q2_%d_s%d",q2b,ss), w_bins, w_min, w_max);
      TH1F *h_w_per_q2_per_sect_temp = new TH1F(Form("h_w_per_q2_%d_s%d_temp",q2b,ss), Form("h_w_per_q2_%d_s%d_temp",q2b,ss), w_bins, w_min, w_max);

      std::vector<double> v_wb_center;
      std::vector<double> v_wb_err;
      std::vector<double> v_cs_value;      
      std::vector<double> v_rate_lum_norm;      
      std::vector<double> v_rate_nel_norm;
      std::vector<double> v_cs_err;
      std::vector<double> v_nel_norm_err;
      std::vector<double> v_lum_norm_err;

      for( int wb = 0; wb <= w_bins; wb++ ){

	// 1 - get counts , bin center
	double w_count = h_el_q2w_sect[ss]->GetBinContent( wb+1, q2b+1 );
	double w_count_err = h_el_q2w_sect[ss]->GetBinError( wb+1, q2b+1 );
	double w_center = h_w_per_q2_per_sect->GetBinCenter(wb+1);
	if( w_center <= 1.1 ) continue;

	//std::cout << " w bin " << wb << " center " << w_center << " content " << w_count <<  " err " << w_count_err << std::endl;
	h_w_per_q2_per_sect->SetBinContent(wb+1, w_count );       
	//h_w_per_q2_per_sect->Draw();

	double nu = (pow(w_center,2) - pow(proton_mass,2) + q2_center)/(2*proton_mass); 
	double del_energy = beam_energy - nu; 
	double theta  =  2.0*asin( sqrt( q2_center/(4.0 * beam_energy * del_energy )) );
	double epsi = 1.0/(1.0 + 2.0*(1.0 + pow(nu,2)/q2_center)*pow(tan(theta/2.0),2));
	double pho_f = (alpha/pow(pi,2)) * (del_energy/beam_energy ) * (pow(w_center,2) - pow(proton_mass,2))/(2.0*proton_mass)/q2_center * (1.0/(1.0-epsi));
	double sector_scaler = 6.0; //remember this is calculated per sector so in effect we are dividing by 6
	double q2bin_scaler = q2_bins; 
	double delta_phi_bins = 1.0;//
	
	//add later
	double accp = 1.0;
	double rad_corr = 1.0;	
	double bin_cent_corr = 1.0;

	double cs = w_count * (1.0/(q2_width * w_width * lum_factor_ug * accp * rad_corr * bin_cent_corr * pho_f * delta_phi_bins));// * sector_scaler ;
	double csErr = sqrt(w_count) * (1.0/(q2_width * w_width * lum_factor_ug * accp * rad_corr * bin_cent_corr *pho_f ));
	double accp_err=0.0;// change later
	double cs_err = sqrt( (csErr/accp)*(csErr/accp) + (cs*accp_err)*(cs*accp_err));

	// normalized values
	double lum_norm =  w_count * (1.0/(lum_factor_ug)); 
	double nel_norm =  w_count * (1.0/(el_charge*n_el_ug)); 

	double lum_err = 0;//sqrt(lum_factor_ug);
	double nel_err = 0;///sqrt(n_el_ug);
	
	double lum_norm_err = 0.0;//sqrt( ((w_count)/accp)*((w_count)/accp) + (cs*lum_err)*(cs*lum_err)); 
	double nel_norm_err = 0.0;//sqrt( ((w_count)/accp)*((w_count)/accp) + (cs*nel_err)*(cs*nel_err)); 

	//pow( sqrt(w_count) * (1.0/( accp * rad_corr * bin_cent_corr * lum_factor_tot * q2_width*w_width ) ),2);
	std::cout << " >> ------------------------------------------- << " << std::endl;
	std::cout << " W BIN NUM    " <<  wb << std::endl;
	std::cout << " W bin center " << h_w_per_q2_per_sect->GetBinCenter(wb+1) << " Q2 bin center " << q2_center << std::endl;
	std::cout << " RAW COUNT    " << w_count << std::endl;
	std::cout << " PHOTON FLUX  " << pho_f << std::endl;
	std::cout << " DELTA W      " << w_width << std::endl;
	std::cout << " DELTA Q2     " << q2_width << std::endl;
	std::cout << " ACCP         " << accp << std::endl;
	std::cout << " RAD CORR     " << rad_corr << std::endl;
	std::cout << " BIN CORR     " << bin_cent_corr << std::endl;
	std::cout << " CS           " << cs << std::endl;
	std::cout << " CS ERR       " << cs_err << std::endl;
	std::cout << " LUM NORM     " << lum_norm << std::endl;
	std::cout << " LUM NORMERR  " << lum_norm_err << std::endl;
	std::cout << " NEL NORM     " << nel_norm << std::endl;
	std::cout << " NEL NORMERR  " << nel_norm_err << std::endl;
       
	v_wb_center.push_back(w_center);
	v_cs_value.push_back( cs );
	v_wb_err.push_back(0.0);
	v_cs_err.push_back(cs_err);

	// plots for monitoring
	v_rate_lum_norm.push_back(lum_norm);
	v_rate_nel_norm.push_back(nel_norm);

	v_nel_norm_err.push_back(nel_norm_err);
	v_lum_norm_err.push_back(lum_norm_err);
	
	h_w_per_q2_per_sect_temp->SetBinContent(wb, cs );
      }
      v_h_w_per_q2_per_sect_temp.push_back(h_w_per_q2_per_sect_temp);
      v_g_cs.push_back( new TGraphErrors(v_wb_center.size(), &v_wb_center[0],&v_cs_value[0], &v_wb_err[0], &v_cs_err[0]) );

      v_g_rate_lum_norm.push_back( new TGraphErrors(v_wb_center.size(), &v_wb_center[0], &v_rate_lum_norm[0], &v_wb_err[0], &v_lum_norm_err[0]) );
      v_g_rate_nel_norm.push_back( new TGraphErrors(v_wb_center.size(), &v_wb_center[0], &v_rate_nel_norm[0], &v_wb_err[0], &v_nel_norm_err[0]) );

      v_wb_temp.push_back(v_wb_center);
      v_cs_value_temp.push_back(v_cs_value);
      v_wb_err_temp.push_back(v_wb_err);
      v_cs_err_temp.push_back(v_cs_err);

    }

    //store information for later
    h_el_q2w_cs[ss] = v_h_w_per_q2_per_sect_temp;
    m_wb_center[ss] = v_wb_temp;
    m_cs_value[ss] = v_cs_value_temp;
    m_wb_err[ss] = v_wb_err_temp;
    m_cs_err[ss] = v_cs_err_temp;
    m_g_cs[ss] = v_g_cs;
    m_g_ln[ss] = v_g_rate_lum_norm;
    m_g_neln[ss] = v_g_rate_nel_norm;

    
    c1->SaveAs(Form("h_w_per_q2_sector_%d.pdf",ss));   
  }
  c0->SaveAs(Form("h_el_q2w_final_per_sect_%s.pdf",config));

  for( int ss = 0; ss < n_sectors; ss++ ){
    TCanvas *c_cs_sect = new TCanvas(Form("c_cs_sect%d",ss),Form("c_cs_sect%d",ss),900,900);
    c_cs_sect->Divide(4,5);
    std::cout << " Making Plot for Sector 1: Number of Q2 bins - " << m_g_cs[ss].size() << std::endl;
    for( int q2b = 0; q2b < m_g_cs[ss].size(); q2b++ ){
      c_cs_sect->cd(q2b+1);
      std::cout << " generating plot for sector " << ss+1 << "q2b range " << q2b*0.25+1.0 << " to  "<< (q2b+1.0)*0.25 + 1.0 << std::endl;
      double q2min = q2b*0.25 + 1.0;
      double q2max = (q2b+1.0)*0.25 + 1.0;
     
      m_g_cs[ss][q2b]->SetTitle(Form("Sector %d %f < Q2 < %f", ss+1, q2min, q2max));
      m_g_cs[ss][q2b]->GetXaxis()->SetTitle("W (GeV)");
      m_g_cs[ss][q2b]->GetYaxis()->SetTitle("Norm. Yields");
      m_g_cs[ss][q2b]->SetMarkerSize(2);
      //m_g_cs[ss][q2b]->Draw("AP+");      
      m_g_cs[ss][q2b]->GetXaxis()->SetLimits(0.98, 2.4);
      //m_g_cs[ss][q2b]->GetHistogram()->SetMinimum(0.0);
      //m_g_cs[ss][q2b]->GetHistogram()->SetMaximum(0.18);
      //m_g_cs[ss][q2b]->Draw("AP+");      
      c_cs_sect->Update();
    }
    c_cs_sect->SaveAs(Form("h_el_q2s_cs_s%d_%s.pdf",ss,config));
  }



  // PLOTS FOR NORMALIZED YIELDS ( LUM NORM)
  for( int ss = 0; ss < n_sectors; ss++ ){
    TCanvas *c_ln_sect = new TCanvas(Form("c_ln_sect%d",ss),Form("c_ln_sect%d",ss),900,900);
    c_ln_sect->Divide(4,5);
    std::cout << " Making Plot for Sector 1: Number of Q2 bins - " << m_g_ln[ss].size() << std::endl;
    for( int q2b = 0; q2b < m_g_ln[ss].size(); q2b++ ){
      c_ln_sect->cd(q2b+1);
      std::cout << " generating plot for sector " << ss+1 << "q2b range " << q2b*0.25+1.0 << " to  "<< (q2b+1.0)*0.25 + 1.0 << std::endl;
      double q2min = q2b*0.25 + 1.0;
      double q2max = (q2b+1.0)*0.25 + 1.0;
     
      m_g_ln[ss][q2b]->SetTitle(Form("Sector %d %f < Q2 < %f", ss+1, q2min, q2max));
      m_g_ln[ss][q2b]->GetXaxis()->SetTitle("W (GeV)");
      m_g_ln[ss][q2b]->GetYaxis()->SetTitle("Lum. Norm. Yields");
      m_g_ln[ss][q2b]->SetMarkerStyle(20);
      m_g_ln[ss][q2b]->SetMarkerSize(0.75);
      //m_g_ln[ss][q2b]->Draw("AP+");      
      //m_g_ln[ss][q2b]->GetXaxis()->SetLimits(0.98, 2.4);
      //m_g_ln[ss][q2b]->GetHistogram()->SetMinimum(0.0);
      //m_g_ln[ss][q2b]->GetHistogram()->SetMaximum(0.18);
      //m_g_ln[ss][q2b]->Draw("AP+");      
      //c_ln_sect->Update();
    }
    c_ln_sect->SaveAs(Form("h_el_q2s_lumnorm_s%d_%s.pdf",ss,config));
  }

  // PLOTS FOR NORMALIZED YIELDS ( NEL NORM )
  for( int ss = 0; ss < n_sectors; ss++ ){
    TCanvas *c_neln_sect = new TCanvas(Form("c_neln_sect%d",ss),Form("c_neln_sect%d",ss),900,900);
    c_neln_sect->Divide(4,5);
    std::cout << " Making Plot for Sector 1: Number of Q2 bins - " << m_g_neln[ss].size() << std::endl;
    for( int q2b = 0; q2b < m_g_neln[ss].size(); q2b++ ){
      c_neln_sect->cd(q2b+1);
      std::cout << " generating plot for sector " << ss+1 << "q2b range " << q2b*0.25+1.0 << " to  "<< (q2b+1.0)*0.25 + 1.0 << std::endl;
      double q2min = q2b*0.25 + 1.0;
      double q2max = (q2b+1.0)*0.25 + 1.0;
     
      m_g_neln[ss][q2b]->SetTitle(Form("Sector %d %f < Q2 < %f", ss+1, q2min, q2max));
      m_g_neln[ss][q2b]->GetXaxis()->SetTitle("W (GeV)");
      m_g_neln[ss][q2b]->GetYaxis()->SetTitle("N_{el} Norm. Yields");
      m_g_neln[ss][q2b]->SetMarkerStyle(20);
      m_g_neln[ss][q2b]->SetMarkerSize(0.75);
      //m_g_neln[ss][q2b]->Draw("AP+");      
      //m_g_ln[ss][q2b]->GetXaxis()->SetLimits(0.98, 2.4);
      //m_g_ln[ss][q2b]->GetHistogram()->SetMinimum(0.0);
      //m_g_ln[ss][q2b]->GetHistogram()->SetMaximum(0.18);
      //m_g_neln[ss][q2b]->Draw("AP+");      
      //c_neln_sect->Update();
    }
    c_neln_sect->SaveAs(Form("h_el_q2s_nelnorm_s%d_%s.pdf",ss,config));
  }
  
	  
  std::vector< std::vector<double> > v_cs_model =   GetModelValues("model_inclusive_q2_1-6_w_1p125-2p115_delta0p03.txt");
  std::vector< std::vector<double> > v_wcenter_model = GetModelWCenterValues("model_inclusive_q2_1-6_w_1p125-2p115_delta0p03.txt");
  std::vector< TGraphErrors* > v_g_model;
  //first loop over the q2 bins 
  double w_start = 1;
  double w_bins = (double)h_el_q2w_sect[0]->GetXaxis()->GetNbins();
  double w_min = h_el_q2w_sect[0]->GetXaxis()->GetBinLowEdge(1);
  double w_max = h_el_q2w_sect[0]->GetXaxis()->GetBinUpEdge(w_bins);
  double w_width = (w_max-w_min)/w_bins;
  
  TCanvas *c_model_data = new TCanvas("c_model_data","c_model_data",900,900);

  for(int q2 =0; q2<v_cs_model.size(); q2++ ){
    //then the individual W bins
    std::vector<double> temp_w_bin_num;
    std::vector<double> temp_model_val;
    std::vector<double> temp_w_bin_err;
    std::vector<double> temp_model_err;    
    for( int ww = v_cs_model[q2].size(); ww > 0; ww-- ){
      double cs_model =  v_cs_model[q2][ww];
      double w_center_model = v_wcenter_model[q2][ww];
      std::cout << " >>> filling MODEL Q2 bin " << q2 << " W bin " << ww << " W bin center " << w_center_model << " CS Model " << cs_model << std::endl;
      double w_bin_center = 1;      
      // get bin info
      temp_w_bin_num.push_back(w_center_model);
      temp_model_val.push_back(cs_model);
      temp_w_bin_err.push_back(0.0);
      temp_model_err.push_back(0.0);      
    }
    TGraphErrors *g_model = new TGraphErrors(temp_w_bin_num.size(), &temp_w_bin_num[0], &temp_model_val[0], &temp_w_bin_err[0], &temp_model_err[0] );
    g_model->SetTitle(Form("Model Q2B %d",q2));
    v_g_model.push_back( g_model );
  }

  ////////////////////////////////////////////////////////////////////
  // draw model with data now
  for( int ss = 0; ss < n_sectors; ss++ ){
    TCanvas *c_modeldata_sect = new TCanvas(Form("c_modeldata_sect%d",ss),Form("c_modeldata_sect%d",ss),900,900);
    c_modeldata_sect->Divide(4,5);
    for( int q2b = 0; q2b < v_g_model.size(); q2b++ ){
      c_modeldata_sect->cd(q2b+1);

      TMultiGraph *mg_cs_result = new TMultiGraph();
      m_g_cs[ss][q2b]->SetMarkerColor(kBlack);
      m_g_cs[ss][q2b]->SetMarkerStyle(20);
      m_g_cs[ss][q2b]->SetMarkerSize(0.6);

      v_g_model[q2b]->SetMarkerColor(kRed);
      v_g_model[q2b]->SetMarkerStyle(22);
      v_g_model[q2b]->SetMarkerSize(0.4);

      
      mg_cs_result->Add(m_g_cs[ss][q2b]);
      mg_cs_result->Add(v_g_model[q2b]);
      mg_cs_result->SetTitle(Form("Cross Section S%d",ss));
      mg_cs_result->Draw("APE");

      //continue adding plotting style code here


    }
  }

  //draw overlapping cross sections per sector
  

  return 0;
}
