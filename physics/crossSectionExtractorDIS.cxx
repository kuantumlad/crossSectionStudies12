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

  TCanvas *c0 = new TCanvas("c0","c0",900,900);
  c0->Divide(3,2);
  for( int ss = 0; ss < 6; ss++ ){
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

    for( int q2b = 1; q2b <= 20; q2b++){
      c1->cd(q2b);
      double w_bins = (double)h_el_q2w_sect[ss]->GetXaxis()->GetNbins();
      double w_min = h_el_q2w_sect[ss]->GetXaxis()->GetBinLowEdge(1);
      double w_max = h_el_q2w_sect[ss]->GetXaxis()->GetBinUpEdge(w_bins);
      double w_width = (w_max-w_min)/w_bins;

      double q2_bins = (double)h_el_q2w_sect[ss]->GetYaxis()->GetNbins();
      double q2_min = h_el_q2w_sect[ss]->GetYaxis()->GetBinLowEdge(1);
      double q2_max = h_el_q2w_sect[ss]->GetYaxis()->GetBinUpEdge(q2_bins);
      double q2_width = (q2_max-q2_min)/q2_bins;
      double q2_center = q2b*q2_width + q2_min;
      std::cout << " Q2 Bin Center " << q2_center <<std::endl;

      TH1F *h_w_per_q2_per_sect = new TH1F(Form("h_w_per_q2_%d_s%d",q2b,ss), Form("h_w_per_q2_%d_s%d",q2b,ss), w_bins, w_min, w_max);
      TH1F *h_w_per_q2_per_sect_temp = new TH1F(Form("h_w_per_q2_%d_s%d_temp",q2b,ss), Form("h_w_per_q2_%d_s%d_temp",q2b,ss), w_bins, w_min, w_max);

      std::vector<double> v_wb_center;
      std::vector<double> v_wb_err;
      std::vector<double> v_cs_value;      
      std::vector<double> v_cs_err;

      for( int wb = 1; wb <= w_bins; wb++ ){

	// 1 - get counts , bin center
	double w_count = h_el_q2w_sect[ss]->GetBinContent( wb, q2b );
	double w_count_err = h_el_q2w_sect[ss]->GetBinError( wb, q2b );
	double w_center = h_w_per_q2_per_sect->GetBinCenter(wb);
	if( w_center <= 1.1 ) continue;

	//std::cout << " w bin " << wb << " center " << w_center << " content " << w_count <<  " err " << w_count_err << std::endl;
	h_w_per_q2_per_sect->SetBinContent(wb, w_count );       
	h_w_per_q2_per_sect->Draw();

	double nu = (pow(w_center,2) - pow(proton_mass,2) + q2_center)/(2*proton_mass); 
	double del_energy = beam_energy - nu; 
	double theta  =  2.0*asin( sqrt( q2_center/(4.0 * beam_energy * del_energy )) );
	double epsi = 1.0/(1.0 + 2.0*(1.0 + pow(nu,2)/q2_center)*pow(tan(theta/2.0),2));
	double pho_f = (alpha/pow(pi,2)) * (del_energy/beam_energy ) * (pow(w_center,2) - pow(proton_mass,2))/(2.0*proton_mass)/q2_center * (1.0/(1.0-epsi));

	//add later
	double accp = 1.0;
	double rad_corr = 1.0;	
	double bin_cent_corr = 1.0;

	double cs = w_count * (1.0/(q2_width * w_width * lum_factor_ug * accp * rad_corr * bin_cent_corr * pho_f ));
	double csErr = sqrt(w_count) * (1.0/(q2_width * w_width * lum_factor_ug * accp * rad_corr * bin_cent_corr *pho_f ));
	double accp_err=0.0;// change later
	double cs_err = sqrt( (csErr/accp)*(csErr/accp) + (cs*accp_err)*(cs*accp_err));

	//pow( sqrt(w_count) * (1.0/( accp * rad_corr * bin_cent_corr * lum_factor_tot * q2_width*w_width ) ),2);
	std::cout << " >> ------------------------------------------- << " << std::endl;
	std::cout << " W BIN NUM    " <<  wb << std::endl;
	std::cout << " W bin center " << h_w_per_q2_per_sect->GetBinCenter(wb) << "Q2 bin center " << q2_center << std::endl;
	std::cout << " RAW COUNT    " << w_count << std::endl;
	std::cout << " PHOTON FLUX  " << pho_f << std::endl;
	std::cout << " DELTA W      " << w_width << std::endl;
	std::cout << " DELTA Q2     " << q2_width << std::endl;
	std::cout << " ACCP         " << accp << std::endl;
	std::cout << " RAD CORR     " << rad_corr << std::endl;
	std::cout << " BIN CORR     " << bin_cent_corr << std::endl;
	std::cout << " CS           " << cs << std::endl;
	std::cout << " CS ERR       " << cs_err << std::endl;
	
	

	v_wb_center.push_back(w_center);
	v_cs_value.push_back( cs );
	v_wb_err.push_back(0.0);
	v_cs_err.push_back(cs_err);

	h_w_per_q2_per_sect_temp->SetBinContent(wb, cs );
      }
      v_h_w_per_q2_per_sect_temp.push_back(h_w_per_q2_per_sect_temp);
      v_g_cs.push_back( new TGraphErrors(v_wb_center.size(), &v_wb_center[0],&v_cs_value[0], &v_wb_err[0], &v_cs_err[0]) );
      
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
    
    c1->SaveAs(Form("h_w_per_q2_sector_%d.pdf",ss));   
  }
  c0->SaveAs(Form("h_el_q2w_final_per_sect_%s.pdf",config));

  for( int ss = 0; ss < 6; ss++ ){
    TCanvas *c_cs_sect = new TCanvas(Form("c_cs_sect%d",ss),Form("c_cs_sect%d",ss),900,900);
    c_cs_sect->Divide(4,5);
    std::cout << " Making Plot for Sector 1: Number of Q2 bins - " << m_g_cs[ss].size() << std::endl;
    for( int q2b = 0; q2b < m_g_cs[ss].size(); q2b++ ){
      c_cs_sect->cd(q2b+1); 
      m_g_cs[ss][q2b]->SetTitle(Form("Sector %d Q2 %d",ss,q2b));
      m_g_cs[ss][q2b]->GetXaxis()->SetTitle("W (GeV)");
      m_g_cs[ss][q2b]->GetYaxis()->SetTitle("CS (nB)");
      m_g_cs[ss][q2b]->SetMarkerSize(2);
      m_g_cs[ss][q2b]->Draw("AP+");      
    }
    c_cs_sect->SaveAs(Form("h_el_q2s_cs_s%d_%s.pdf",ss,config));
  }
	  
    

  return 0;
}
