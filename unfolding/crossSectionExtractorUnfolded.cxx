#include <iostream>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TMath.h>
#include <TH2F.h>
#include <vector>
#include <map>
#include <TLine.h>

TH1D* getSectorCrossSection(TH1D *h_in, TH1D *h_in_unfolded, std::map<int, std::vector<double> > accp_corr, std::map<int, std::vector<double> > accp_corr_err ){

  int phi_bins = 72;//h_in->GetNbinsX();


  //std::cout << " Bin with max entries: " << temp_bin << " with entries: " << temp_entries << std::endl;
  int delta_bin = 0;
  int bin_min = 1.0;//temp_bin;// - delta_bin;
  int bin_max = 1.0;//temp_bin + delta_bin;
  
  std::cout << " Extending bins by " << delta_bin << " bins on each side; bin min: " << bin_min << " bin max: " << bin_max << std::endl;
  
  TH1D *h_theta_proj = h_in_unfolded; //h2_in->ProjectionY(Form("proj_%s",h2_in->GetTitle()),bin_min, bin_max);
  
  std::vector<double> accp_corr_phibin = accp_corr[37];//bin_min-1];
  std::vector<double> accp_corr_phibin_err = accp_corr_err[37];//bin_min-1];

  //h_theta_proj->Rebin(5);
  
  double smudge = 1.0;
  //double lum_run2391 = 6091588.274652874; //609158.862236702; //6.0915886E6;
  double lum_run2587 = 10469161.819; //2113774.4627232; //1164585.2;//1288395.597;//824013.67435; //713825.6009; //897613.632;// 89637.465;//47362.916709528225;// 171137.8663359962; //211835.29;// 16663.43;//
  double lum_sim = 1.0;
  double bin_phi_size = ( 2.0 * 3.141592658 ) / (double)phi_bins; //h_in->GetNbinsX();
  double bin_theta_size = 0.0175;//(h_theta_proj->GetBinCenter(2) - h_theta_proj->GetBinCenter(1)) * ( 3.1415/180.0);
  double theta_max = h_theta_proj->GetBinCenter(h_theta_proj->GetNbinsX()) + (h_theta_proj->GetBinCenter(2) - h_theta_proj->GetBinCenter(1))/2.0;

  std::cout << " >> Bin info " << std::endl;
  std::cout << " Number of Phi Bins " << phi_bins << std::endl;//h_in_unfolded->GetNbinsX() << std::endl;
  std::cout << " Number of Theta Bins " << h_theta_proj->GetNbinsX() << std::endl;  
  std::cout << " >> theta max " << theta_max << std::endl;
  std::cout << " >> phi bin size " <<bin_phi_size << " theta bin size " << bin_theta_size << std::endl;
  std::cout << " new theta bin count " << h_theta_proj->GetNbinsX() << std::endl;
  TH1D *h_temp = new TH1D(Form("h_cs_%s",h_in->GetTitle()), Form("h_cs_%s",h_in->GetTitle()), h_theta_proj->GetNbinsX(), 0.0, theta_max);
  h_temp->GetXaxis()->SetTitle("#theta [deg]");
  h_temp->GetYaxis()->SetTitle("mb/sr");

  std::cout << " creating CS histogram with " << h_theta_proj->GetNbinsX() << std::endl;
  for( int b = 1; b <= h_theta_proj->GetNbinsX(); b++ ){
    double accp_theta_bin = accp_corr_phibin[b];
    double accp_theta_bin_err =  accp_corr_phibin_err[b];
    if( accp_theta_bin <= 0.001 || accp_theta_bin == 0 ) accp_theta_bin = 1.0; // prevent division by 0   
    std::cout << b <<  " >> bin center " << h_temp->GetBinCenter(b) << " proj bin value " << h_theta_proj->GetBinContent(b) << " accp corr " <<  accp_theta_bin << std::endl;
    double sin_angle = TMath::Sin( h_theta_proj->GetBinCenter(b) * 3.1415/180.0 );
    delta_bin = 1.0;
    double newbincontent= h_theta_proj->GetBinContent(b) / ( accp_theta_bin * lum_run2587 * (bin_phi_size*delta_bin*1.0) * bin_theta_size * sin_angle );
    int bincontent = h_theta_proj->GetBinContent(b);

    //double nev_theta_err = pow(sqrt(h_theta_proj->GetBinContent(b)) * (1.0/( accp_theta_bin * lum_run2587 * (bin_phi_size*delta_bin*1.0) * bin_theta_size * sin_angle ) ),2);
    double nev_theta_err = newbincontent / (accp_theta_bin * lum_run2587 * (bin_phi_size*delta_bin*1.0) * bin_theta_size * sin_angle );
    double acc_err = pow( accp_theta_bin_err/accp_theta_bin, 2) * pow(newbincontent, 2) * (1.0/accp_theta_bin);
    //double acc_err = pow( accp_theta_bin_err * ( (1.0/accp_theta_bin)*newbincontent ), 2);
    double newbincontent_error = newbincontent*sqrt( nev_theta_err + acc_err );

    //double newbincontent_error = newbincontent * (sqrt(bincontent)/bincontent);    
    std::cout << b <<  " >> bin center " << h_temp->GetBinCenter(b) << " proj bin value " << h_theta_proj->GetBinContent(b) << " cross section value " << newbincontent << " error " << newbincontent_error << " acceptance value " << accp_theta_bin << std::endl;
    h_temp->SetBinContent(b,newbincontent);
    h_temp->SetBinError(b,newbincontent_error);
  }

  return h_temp;
}


std::map<int, std::vector<double> > GetAcceptanceFactors(const char* fileName ){

  std::map<int, std::vector<double> > acceptance_values;

  double acceptance;
  int theta_bin = 0;
  int phi_bin=0;

  std::string parentDirectory = "/w/hallb-scifs17exp/clas12/bclary/CLAS12/electron_studies/physics/parameters/";
	
  string line;
  ifstream readFromAcceptance(parentDirectory+fileName);
  if( readFromAcceptance.is_open() ){
    std::cout << " OPENED FILES " << std::endl;	
    for( std::string line; getline( readFromAcceptance, line); ){
      //std::cout << line << std::endl;
      
      std::stringstream stream(line);
      std::string accp_val;
      std::cout << " GETTING VALUES FOR PHI BIN " << phi_bin << std::endl;
      std::vector<double> theta_accp_val;
      while( stream >> accp_val ){       
	theta_accp_val.push_back(std::atof(accp_val.c_str()));
	//std::cout << " accp_val " << std::atof(accp_val.c_str());
      }
      acceptance_values[phi_bin]=theta_accp_val;
      //std::cout << "\n";
      phi_bin++;
      theta_accp_val.clear();    
    }
  }
  else{
    std::cout << " ERROR OPENING FILE " << std::endl;
  }

  return acceptance_values;
}


int crossSectionExtractorUnfolded(const char* infile, int run){


  TFile *fIn = new TFile(infile,"");
  TFile *fOut = new TFile(Form("final_elastic_cs_unfolded_run%d.root",run),"RECREATE");
  
  if( fIn->IsZombie() ){
    std::cout << " Input file "<< fIn << " doesn't exist" << std::endl;
    std::cout << " bye " << std::endl;
    return 0;
  }

  TH1D *h_model = new TH1D("h_elastic_model_bin_theta1deg","h_elastic_model_bin_theta1deg",30,0.0, 30.0);
  double model_value;
  double model_val;
  double model_theta=5.5;
  int model_theta_bin=6;
  std::string model = "model_cross_section.txt";
  string line;
  ifstream readFromEModel(model);//parentDirectory+"w_cut_limits_run"+std::to_string(run)+".txt");
  if( readFromEModel.is_open() ){
    std::cout << " OPENED FILES " << std::endl;	
    while(readFromEModel >> model_value ) {//std::getline (readFromWCut, line) ){      
      ////readFromEModel >> model_val;
      if ( model_theta < 9 ) {
	std::cout << " dont plot model theta " << std::endl;
      }
      else{
	h_model->SetBinContent(model_theta_bin,model_value);
      }
      std::cout << " MODEL BIN " << model_theta_bin << " MODEL THETA " << model_theta <<  " >> MODEL VALUE " << model_value << std::endl;
      model_theta+=1.0;
      model_theta_bin+=1;
    
    }
  }

  //get acceptance correction map first
  std::map<int, std::vector<double> > accp_corr = GetAcceptanceFactors("elastic_theta_phi_acceptance_2587_ftm06sm06.txt");
  std::map<int, std::vector<double> > accp_corr_err = GetAcceptanceFactors("elastic_theta_phi_acceptance_error_2587_ftm06sm06.txt");

  std::vector< TH1D* > h_el_phi_sect_final;
  std::vector< TH1D* > h_el_phi_sect_final2;
  
  for( int s = 0; s < 6; s++ ){
    std::cout << " Sector " << s <<std::endl;
    h_el_phi_sect_final.push_back( (TH1D*)fIn->Get(Form("h_elastic_cs_unfolded_s%d",s+1) ) );
    h_el_phi_sect_final2.push_back( (TH1D*)fIn->Get(Form("h_elastic_cs_unfolded_s%d",s+1) ) );
    //std::cout << " Title of input histgram " << h_el_phi_sect_final[s]->GetTitle() << std::endl;
  }

  std::vector<TH1D*> cs_results;
  for( int s=0; s < 5; s++ ){
    TH1D *h_temp_cs = getSectorCrossSection( h_el_phi_sect_final2[s], h_el_phi_sect_final[s], accp_corr, accp_corr_err );
    cs_results.push_back(h_temp_cs);    
  }


  TCanvas *cs_out = new TCanvas("cs_out","cs_out",800,1200);
  cs_out->Divide(2,3);
  for( int s = 0; s < 5; s++ ){
    cs_out->cd(s+1);
    //gPad->SetLogy();
    h_model->SetLineColor(kRed);
    h_model->GetXaxis()->SetTitle("#theta [deg]");
    h_model->GetXaxis()->CenterTitle();    
    h_model->Draw();
    h_model->SetAxisRange(14,26,"X");
    h_model->Draw();
    cs_results[s]->Draw("same");
    cs_results[s]->SetAxisRange(14,26,"X");
    cs_results[s]->Draw("same");
  }
  cs_out->SaveAs(Form("cs_results_allsectors_r%d.pdf",run));

  // get difference in model and measured
  std::vector<TGraphErrors*> g_result_diff;
  std::vector<TGraphErrors*> g_cs_result;
  std::vector<TGraphErrors*> g_cs_model;
  for( int s = 0; s<cs_results.size(); s++ ){
    std::vector<double> bins;
    std::vector<double> bin_error_x;
    std::vector<double> bin_error_y;
    std::vector<double> model;
    std::vector<double> model_error;
    std::vector<double> data;
    std::vector<double> data_diff;
    std::vector<double> data_diff_err_x;   
    std::vector<double> data_diff_err_y;   
    std::cout << " comparing model and data for sector " << s << std::endl;
    for( int b=1; b<=h_model->GetNbinsX(); b++){
      if ( b <= 14 || b >= 24 ) continue;
      double bin_center_model = h_model->GetBinCenter(b);
      double bin_center_data = cs_results[s]->GetBinCenter(b);
      double model_result = h_model->GetBinContent(b);
      double data_result = cs_results[s]->GetBinContent(b);
      double data_err_y = cs_results[s]->GetBinError(b);
      double data_err_x = 1.0;
      double data_ratio = data_result/model_result;
      double data_ratio_err = data_ratio * data_err_y;
      std::cout << " model result " << model_result << " data result " << data_result << " data error " << data_err_y << std::endl;
      bins.push_back(bin_center_data);
      bin_error_x.push_back(0.0);
      bin_error_y.push_back(data_err_y);
      model_error.push_back(0.0);
      if ( model_result == 0.0 ) continue;
      model.push_back((model_result));
      data.push_back((data_result));
      data_diff.push_back( data_ratio );
      data_diff_err_y.push_back( data_ratio_err );			  
      data_diff_err_x.push_back( 0.0 );			        
      std::cout << " bin " << b << " bin center " << bin_center_model << " bin center data " << bin_center_data << " diff " << data_result/model_result << std::endl;
    }
    g_result_diff.push_back( new TGraphErrors(bins.size(), &(bins[0]), &(data_diff[0]), &(data_diff_err_x[0]), &(data_diff_err_y[0]) ) );
    g_cs_result.push_back( new TGraphErrors(bins.size(), &(bins[0]), &(data[0]), &(bin_error_x[0]), &(bin_error_y[0]) ) );
    g_cs_model.push_back( new TGraphErrors(bins.size(), &(bins[0]), &(model[0]), &(model_error[0]) ) );
  }
  
  TCanvas *c_diff = new TCanvas("c_diff","c_diff",400,600);
  c_diff->Divide(2,3);
  for( int s = 0; s < g_result_diff.size(); s++ ){
    c_diff->cd(s+1);
    g_result_diff[s]->SetTitle(Form("Ratio of Cross Section Model to Data Sector %d",s+1));
    g_result_diff[s]->GetXaxis()->SetTitle("#theta [deg]");
    g_result_diff[s]->GetXaxis()->CenterTitle();
    g_result_diff[s]->SetMarkerStyle(20);
    g_result_diff[s]->SetMarkerSize(0.5);
    g_result_diff[s]->Draw("AP");
    g_result_diff[s]->GetXaxis()->SetLimits(10.0, 26.5);//SetRangeUser(0.0,26.5);
    g_result_diff[s]->GetHistogram()->SetMaximum(1.5);   // along          
    g_result_diff[s]->GetHistogram()->SetMinimum(0.50);  //   Y     
    g_result_diff[s]->Draw("AP");
    c_diff->Update();
  }
  c_diff->SaveAs(Form("g_ratio_model_data_r%d.pdf",run));
  
  std::vector<TMultiGraph*> v_mg_cs;
  TCanvas *c_result = new TCanvas("c_result","c_result",400,600);
  c_result->Divide(2,3);

  for( int s = 0; s < g_cs_result.size(); s++ ){    
    c_result->cd(s+1);
    c_result->SetGrid();
    v_mg_cs.push_back( new TMultiGraph() );
    g_cs_result[s]->SetTitle(Form("Elastic Model Cross Section Sector %d",s+1));
    g_cs_result[s]->GetXaxis()->SetTitle("#theta [deg]");
    g_cs_result[s]->GetXaxis()->CenterTitle();
    g_cs_result[s]->SetMarkerStyle(20);
    g_cs_model[s]->SetMarkerStyle(22);
    //g_result_diff[s]->SetMarkerSize(3);
    //g_cs_result[s]->Draw("AP");

    //g_cs_model[s]->Draw("AP");
    g_cs_result[s]->SetMarkerSize(0.5);
    g_cs_model[s]->SetMarkerSize(0.5);
    g_cs_model[s]->SetMarkerColor(kRed);
    v_mg_cs[s]->Add(g_cs_result[s]);
    //v_mg_cs[s]->Add(g_cs_model[s]);
    v_mg_cs[s]->SetTitle(Form("Cross Section Sector %d",s+1));
    v_mg_cs[s]->GetXaxis()->SetTitle("#theta [deg]");
    v_mg_cs[s]->GetXaxis()->CenterTitle();
    v_mg_cs[s]->GetYaxis()->SetTitle("#sigma [mb/sr]");
    v_mg_cs[s]->GetYaxis()->CenterTitle();
    v_mg_cs[s]->Draw("APE");
    v_mg_cs[s]->GetXaxis()->SetLimits(13.0, 22.0);
    v_mg_cs[s]->GetHistogram()->SetMaximum(1.9); 
    v_mg_cs[s]->GetHistogram()->SetMinimum(0.001); 
    v_mg_cs[s]->Draw("APE");
    h_model->Draw("SAME C");

    TLegend *legend = new TLegend(0.6, 0.7, 0.89, 0.89);
    legend->AddEntry(g_cs_result[s],"Data");
    //legend->AddEntry(g_cs_model[s],"Model");
    legend->AddEntry(h_model,"Model");  
    legend->SetBorderSize(0);
    legend->Draw();

  }
  c_result->SaveAs(Form("g_cs_result_r%d.pdf",run));

  readFromEModel.close();
  fOut->Write();

  return 0;

}
