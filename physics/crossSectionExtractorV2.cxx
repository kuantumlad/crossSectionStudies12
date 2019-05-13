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

std::vector<TH1D*> getPhiBinCrossSection(TH2F *h2_in, std::map<int, std::vector<double> > accp_corr, std::map<int, std::vector<double> > accp_corr_err){

  int n_phi_bins = h2_in->GetNbinsX();
  std::vector<TH1D*> h_cs_per_phi_bin;
    


  for( int phib = 1; phib < n_phi_bins; phib++ ){

    TH1D *h_theta_proj = h2_in->ProjectionY(Form("proj_%s",h2_in->GetTitle()),phib, phib);
    std::vector<double> accp_corr_phibin = accp_corr[phib-1];
    std::vector<double> accp_corr_phibin_err = accp_corr_err[phib-1];
    
    double lum_run2587 = 10469161.819; //2113774.4627232; //1164585.2;//1288395.597;//824013.67435; //713825.6009; //897613.632;// 89637.465;//47362.916709528225;// 171137.8663359962; //211835.29;// 16663.43;//
    double lum_sim = 1.0;
    double bin_phi_size = ( 2.0 * 3.141592658 ) / (double)h2_in->GetNbinsX();
    double bin_theta_size = 0.0175;//(h_theta_proj->GetBinCenter(2) - h_theta_proj->GetBinCenter(1)) * ( 3.1415/180.0);
    double theta_max = h_theta_proj->GetBinCenter(h_theta_proj->GetNbinsX()) + (h_theta_proj->GetBinCenter(2) - h_theta_proj->GetBinCenter(1))/2.0;

    std::cout << " >> Bin info " << std::endl;
    std::cout << " Number of Phi Bins " << h2_in->GetNbinsX() << std::endl;
    std::cout << " Number of Theta Bins " << h_theta_proj->GetNbinsX() << std::endl;  
    std::cout << " >> theta max " << theta_max << std::endl;
    std::cout << " >> phi bin size " <<bin_phi_size << " theta bin size " << bin_theta_size << std::endl;
    std::cout << " new theta bin count " << h_theta_proj->GetNbinsX() << std::endl;
    TH1D *h_temp = new TH1D(Form("h_cs_%d",phib), Form("h_cs_%d",phib), 30, 0.0, 30.0);
    h_temp->GetXaxis()->SetTitle("#theta [deg]");
    h_temp->GetYaxis()->SetTitle("mb/sr");

    
    std::cout << " creating CS histogram with " << h_theta_proj->GetNbinsX() << std::endl;
    for( int b = 1; b <= h_theta_proj->GetNbinsX(); b++ ){
      double accp_theta_bin = accp_corr_phibin[b];
      double accp_theta_bin_err = accp_corr_phibin_err[b];
      if( accp_theta_bin <= 0.001 || accp_theta_bin == 0 ) accp_theta_bin = 1.0; // prevent division by 0   
      std::cout << b <<  " >> bin center " << h_temp->GetBinCenter(b) << " proj bin value " << h_theta_proj->GetBinContent(b) << " accp corr " <<  accp_theta_bin << std::endl;
      double sin_angle = TMath::Sin( h_theta_proj->GetBinCenter(b) * 3.1415/180.0 );
      double delta_bin = 1.0;
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

    h_cs_per_phi_bin.push_back(h_temp);
  }
  
  return h_cs_per_phi_bin;
}



TH1D* getSectorCrossSection(TH1F *h_in, TH2F *h2_in, std::map<int, std::vector<double> > accp_corr, std::map<int, std::vector<double> > accp_corr_err ){

  int phi_bins = h_in->GetNbinsX();
  std::cout << " Number of phi bins: " << phi_bins << std::endl;
  std::cout << " Finding phi bin with most events " << std::endl;
  int temp_entries=-1;
  int temp_bin = -1;
  for( int b = 0; b <= phi_bins; b++ ){
    int b_entries = h_in->GetBinContent(b);
    std::cout << " bin entries " << b_entries << std::endl;

    if( b_entries > temp_entries ){
      temp_entries=b_entries;
      temp_bin = b;
    }   
  }


  std::cout << " Bin with max entries: " << temp_bin << " with entries: " << temp_entries << std::endl;
  int delta_bin = 0;
  int bin_min = temp_bin;// - delta_bin;
  int bin_max = temp_bin + delta_bin;
  
  std::cout << " Extending bins by " << delta_bin << " bins on each side; bin min: " << bin_min << " bin max: " << bin_max << std::endl;
  
  TH1D *h_theta_proj = h2_in->ProjectionY(Form("proj_%s",h2_in->GetTitle()),bin_min, bin_max);
  std::vector<double> accp_corr_phibin = accp_corr[bin_min-1];
  std::vector<double> accp_corr_phibin_err = accp_corr_err[bin_min-1];

  //h_theta_proj->Rebin(5);
 
  //double lum_run2391 = 6091588.274652874; //609158.862236702; //6.0915886E6;
  double lum_run2587 = 10469161.819; //2113774.4627232; //1164585.2;//1288395.597;//824013.67435; //713825.6009; //897613.632;// 89637.465;//47362.916709528225;// 171137.8663359962; //211835.29;// 16663.43;//
  double lum_sim = 1.0;
  double bin_phi_size = ( 2.0 * 3.141592658 ) / (double)h_in->GetNbinsX();
  double bin_theta_size = 0.0175;//(h_theta_proj->GetBinCenter(2) - h_theta_proj->GetBinCenter(1)) * ( 3.1415/180.0);
  double theta_max = h_theta_proj->GetBinCenter(h_theta_proj->GetNbinsX()) + (h_theta_proj->GetBinCenter(2) - h_theta_proj->GetBinCenter(1))/2.0;

  std::cout << " >> Bin info " << std::endl;
  std::cout << " Number of Phi Bins " << h2_in->GetNbinsX() << std::endl;
  std::cout << " Number of Theta Bins " << h_theta_proj->GetNbinsX() << std::endl;  
  std::cout << " >> theta max " << theta_max << std::endl;
  std::cout << " >> phi bin size " <<bin_phi_size << " theta bin size " << bin_theta_size << std::endl;
  std::cout << " new theta bin count " << h_theta_proj->GetNbinsX() << std::endl;
  TH1D *h_temp = new TH1D(Form("h_cs_%s",h_in->GetTitle()), Form("h_cs_%s",h_in->GetTitle()), h_theta_proj->GetNbinsX(), 0.0, 30.0);
  h_temp->GetXaxis()->SetTitle("#theta [deg]");
  h_temp->GetYaxis()->SetTitle("mb/sr");

  std::cout << " creating CS histogram with " << h_theta_proj->GetNbinsX() << std::endl;
  for( int b = 1; b <= h_theta_proj->GetNbinsX(); b++ ){
    double accp_theta_bin = accp_corr_phibin[b];
    double accp_theta_bin_err = accp_corr_phibin_err[b];
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


std::vector<double> GetBinCenteringCorr(const char* fileName){

  std::vector< double > bin_center_corr;
  double model_theta = 9.5;
  int model_theta_bin = 0;
  double bin_corr;
  double bin_center;
  
  std::string parentDir = "/work/clas12/bclary/CLAS12/electron_studies/bin_centering/";
  string line;
  ifstream readFromBinCorr(parentDir+fileName);
  if( readFromBinCorr.is_open() ){
    std::cout << " OPENED FILES " << std::endl;	
    while(readFromBinCorr >> bin_center >> bin_corr ) {
      //if ( model_theta < 9 ) {/
      //std::cout << " dont plot model theta " << std::endl;
      //}
      //else{
      bin_center_corr.push_back(bin_corr);
      //}
      std::cout << " MODEL THETA " << model_theta <<  " >> BIN CORRECTION " << bin_corr << std::endl;
      model_theta+=1.0;
      model_theta_bin+=1;    
    }
  }

  return bin_center_corr;
}


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
    while(readFromRadCorr >> rad_corr ) {
      //if ( model_theta < 9 ) {/
      //std::cout << " dont plot model theta " << std::endl;
      //}
      //else{
      v_rad_corr.push_back(rad_corr);
      //}
      std::cout << " MODEL THETA " << model_theta <<  " >> BIN CORRECTION " << rad_corr << std::endl;
      model_theta+=1.0;
      model_theta_bin+=1;    
    }
  }

  return v_rad_corr;
}


int crossSectionExtractorV2(const char* infile, int run){

  TFile *fIn = new TFile(infile,"");
  TFile *fOut = new TFile(Form("cs_test_run%dV2.root",run),"RECREATE");
  
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

  //get bin centering corrections
  std::vector<double> bin_center_corr = GetBinCenteringCorr("bin_corr_values.txt");

  // get radiative corrections 
  std::vector<double> rad_corr = GetRadiativeCorr("elastrc.dat");
  

  std::vector< TH1F* > h_el_p_sect_final;
  std::vector< TH1F* > h_el_theta_sect_final;
  std::vector< TH1F* > h_el_phi_sect_final;
  std::vector< TH2F* > h_el_ptheta_sect_final;
  std::vector< TH2F* > h_el_phitheta_sect_final;

  TH2F *h_el_phitheta_final = (TH2F*)fIn->Get(Form("/kinematics/h_el_phitheta_final_run%d",run));
  
  for( int s = 1; s <= 6; s++ ){
    h_el_p_sect_final.push_back( (TH1F*)fIn->Get(Form("/kinematics/h_el_p_s%d_final",s) ) );
    h_el_theta_sect_final.push_back( (TH1F*)fIn->Get(Form("/kinematics/h_el_theta_s%d_final",s) ) );
    h_el_phi_sect_final.push_back( (TH1F*)fIn->Get(Form("/kinematics/h_el_phi_s%d_final",s) ) );

    h_el_ptheta_sect_final.push_back( (TH2F*)fIn->Get(Form("/kinematics/h_el_ptheta_s%d_final",s) ) );
    h_el_phitheta_sect_final.push_back( (TH2F*)fIn->Get(Form("/kinematics/h_el_phitheta_s%d_final",s) ) );
  }
  

  int phi_bins = h_el_phi_sect_final[0]->GetNbinsX();
  std::cout << " Number of phi bins: " << phi_bins << std::endl;

  int temp_entries=-1;
  int temp_bin = -1;
  for( int b = 0; b <= phi_bins; b++ ){
    int b_entries = h_el_phi_sect_final[0]->GetBinContent(b);
    std::cout << " bin entries " << b_entries << std::endl;

    if( b_entries > temp_entries ){
      temp_entries=b_entries;
      temp_bin = b;
    }   
  }

  std::cout << " Bin with max entries: " << temp_bin << " with entries: " << temp_entries << std::endl;
  int delta_bin = 0;
  int bin_min = temp_bin - delta_bin;
  int bin_max = temp_bin + delta_bin;
  std::cout << " min x " << h_el_phi_sect_final[0]->GetBinCenter(bin_min)  << std::endl;

  std::cout << " Extending bins by 2 bins on each side; bin min: " << bin_min << " bin max: " << bin_max << std::endl;
  
  TH1D *h_theta_proj = h_el_phitheta_sect_final[0]->ProjectionY("proj_s1_2391",39,39);//bin_min, bin_max);
  std::vector<double> accp_corr_phibin = accp_corr[38]; //38 because starts at 39 -1 
  std::vector<double> accp_corr_phibin_err = accp_corr_err[38]; //38 because starts at 38 = 39 -1 

  TCanvas *c0 = new TCanvas("c0","c0",800,800);
  gStyle->SetOptStat(0);
  c0->cd(1);
  h_el_phitheta_sect_final[0]->SetTitle("#phi vs #theta (S1)");
  h_el_phitheta_sect_final[0]->GetXaxis()->SetTitle("#phi [deg]");
  h_el_phitheta_sect_final[0]->GetXaxis()->CenterTitle();
  h_el_phitheta_sect_final[0]->GetYaxis()->SetTitle("#theta [deg]");
  h_el_phitheta_sect_final[0]->GetYaxis()->CenterTitle();
  h_el_phitheta_sect_final[0]->Draw("colz");
  c0->Update();
  double phi_center = h_el_phi_sect_final[0]->GetBinCenter(bin_min);
  double phi_min_center = phi_center - (h_el_phi_sect_final[0]->GetBinCenter(bin_min-1) - phi_center)/2.0;
  double phi_max_center = phi_center + (h_el_phi_sect_final[0]->GetBinCenter(bin_min-1) - phi_center)/2.0;
  std::cout << " phi min " << phi_min_center << " phi max " << phi_max_center << std::endl;

  TLine *phi_min = new TLine(phi_min_center, 0.0, phi_min_center, 30.0);
  TLine *phi_max = new TLine(phi_max_center, 0.0, phi_max_center, 30.0);
  phi_min->SetLineColor(kRed);
  phi_min->Draw();
  phi_max->SetLineColor(kRed);
  phi_max->Draw();
  c0->SaveAs(Form("h_phitheta_s1_%d.pdf",run));

  TCanvas *c1 = new TCanvas("c1","c1",1600,800);
  gStyle->SetOptStat(0);
  c1->Divide(2,1);
  c1->cd(1);
  h_el_phitheta_sect_final[0]->SetTitle("#phi vs #theta (S1)");
  h_el_phitheta_sect_final[0]->GetXaxis()->SetTitle("#phi [deg]");
  h_el_phitheta_sect_final[0]->GetXaxis()->CenterTitle();
  h_el_phitheta_sect_final[0]->GetYaxis()->SetTitle("#theta [deg]");
  h_el_phitheta_sect_final[0]->GetYaxis()->CenterTitle();
  h_el_phitheta_sect_final[0]->Draw("colz");
  c1->Update();

  phi_min->SetLineColor(kRed);
  phi_min->Draw();
  phi_max->SetLineColor(kRed);
  phi_max->Draw();
  
  c1->cd(2);
  //h_theta_proj->Rebin(10);

  double lum_run2391 = 6091588.62236702; //6.0915886E6;
  double lum_run2587 = 10469161.819;//4762209.0;//10469161.819; //2113774.4627232; //1164585.2;//1288395.597;// 824013.67435; //713825.6009;//897613.632;//47362.916709528225;// 1711378.663359962;// this value is for 2476run //211835.29;// 16663.43;
  double bin_phi_size = (h_el_phi_sect_final[0]->GetBinCenter(bin_max) - h_el_phi_sect_final[0]->GetBinCenter(bin_max-1)) * 3.1415/180.0;
  double bin_theta_size = (h_theta_proj->GetBinCenter(2) - h_theta_proj->GetBinCenter(1)) * (3.1415/180.0);
  double theta_max = h_theta_proj->GetBinCenter(h_theta_proj->GetNbinsX()) +  (h_theta_proj->GetBinCenter(2) - h_theta_proj->GetBinCenter(1))/2.0;

  std::cout << " >> theta max " << theta_max << std::endl;
  std::cout << " >> phi bin size " <<bin_phi_size << " theta bin size " << bin_theta_size << std::endl;
  std::cout << " new theta bin count " << h_theta_proj->GetNbinsX() << std::endl;
  TH1D *h_temp = new TH1D("h_Temp","h_Temp", h_theta_proj->GetNbinsX(), 0.0, 30.0);

  for( int b = 1; b <= h_theta_proj->GetNbinsX(); b++ ){
    double accp_theta_bin = accp_corr_phibin[b];
    double accp_theta_bin_err = accp_corr_phibin_err[b];
    if( accp_theta_bin <= 0.001 || accp_theta_bin == 0 ) accp_theta_bin = 1.0; // prevent division by 0
    std::cout << b <<  " >> bin center " << h_temp->GetBinCenter(b) << " proj bin value " << h_theta_proj->GetBinContent(b) << " acceptance correction " << accp_theta_bin << std::endl;
    double sin_angle = TMath::Sin( h_theta_proj->GetBinCenter(b) * 3.1415/180.0 );
    delta_bin=1.0;
    double newbincontent= h_theta_proj->GetBinContent(b) / (accp_theta_bin * lum_run2587 * (bin_phi_size*(delta_bin*1.0)) * bin_theta_size * sin_angle );

    double nev_theta_err = pow(sqrt(h_theta_proj->GetBinContent(b)) * (1.0/( accp_theta_bin * lum_run2587 * (bin_phi_size*delta_bin*1.0) * bin_theta_size * sin_angle ) ),2);
    double acc_err = pow( sqrt(accp_theta_bin_err) * ( (1.0/accp_theta_bin)*newbincontent ),2);
    double newbincontent_error = newbincontent*sqrt( nev_theta_err + acc_err );

    std::cout << " cross section " << newbincontent << std::endl;
    h_temp->SetBinContent(b,newbincontent);
    //h_temp->SetBinContentError(b, newbincontent_error);
    
  }

  h_model->SetLineColor(kRed);
  h_model->SetTitle("Cross Section (S1)");
  h_model->GetXaxis()->SetTitle("#theta [deg]");
  h_model->GetXaxis()->CenterTitle();
  h_model->Draw();
  h_temp->Draw("same");
  c1->SaveAs(Form("h_phitheta_csresult_s1_%d.pdf",run));


  //////////////////////////////////////////////////////////////////
  // cross section at phi bin with most events / sector  
  std::vector<TH1D*> cs_results;
  for( int s=0; s<6; s++ ){
    TH1D *h_temp_cs = getSectorCrossSection( h_el_phi_sect_final[s], h_el_phitheta_sect_final[s], accp_corr, accp_corr_err );
    cs_results.push_back(h_temp_cs);
  }

  //////////////////////////////////////////////////////////////////
  // cross section for each bin in phi
  std::vector<TH1D*> cs_results_per_phi_bin = getPhiBinCrossSection(h_el_phitheta_final, accp_corr, accp_corr_err);



  TCanvas *c_raw = new TCanvas("c_raw","c_raw",800,800);
  c_raw->cd(1);
  gStyle->SetOptStat(0);
  h_theta_proj->GetXaxis()->SetTitle("#theta [deg]");
  h_theta_proj->GetXaxis()->CenterTitle();
  h_theta_proj->SetTitle("Raw #theta Counts");
  h_theta_proj->Draw();
  c_raw->SaveAs(Form("h_theta_proj_raw_%d.pdf",run));

  TCanvas *cs_out = new TCanvas("cs_out","cs_out",800,1200);
  cs_out->Divide(2,3);
  for( int s = 0; s < 6; s++ ){
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

      //corrections
      double bin_center_correction = bin_center_corr[b-10]; // minus 14 + 6 because start histogram with bin center of 15.5 
      double radiative_correction = rad_corr[b-6];

      std::cout << " model result " << model_result << " data result " << data_result << " data error " << data_err_y << std::endl;
      std::cout << " center bin " << bin_center_data << "  bin centering corection  " << bin_center_correction << std::endl;
      std::cout << " center bin " << bin_center_data << " radiative correction " << radiative_correction << std::endl;

      double final_data_result = data_result/(bin_center_correction*radiative_correction);

      bins.push_back(bin_center_data);
      bin_error_x.push_back(0.0);
      bin_error_y.push_back(data_err_y);
      model_error.push_back(0.0);
      if ( model_result == 0.0 ) continue;
      model.push_back((model_result));
      data.push_back((final_data_result));
      data_diff.push_back( data_ratio );
      data_diff_err_y.push_back( data_ratio_err );			  
      data_diff_err_x.push_back( 0.0 );			        
      std::cout << " bin " << b << " bin center " << bin_center_model << " bin center data " << bin_center_data << " diff " << final_data_result/model_result << std::endl;
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

  std::vector<TMultiGraph*> v_mg_cs_logy;

  TCanvas *c_result_log = new TCanvas("c_result_log","c_result_log",400,600);
  c_result_log->Divide(2,3);

  for( int s = 0; s < g_cs_result.size(); s++ ){    
    c_result_log->cd(s+1);
    gPad->SetLogy();
    c_result_log->SetGrid();
    v_mg_cs_logy.push_back( new TMultiGraph() );
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
    v_mg_cs_logy[s]->Add(g_cs_result[s]);
    //v_mg_cs_logy[s]->Add(g_cs_model[s]);
    v_mg_cs_logy[s]->SetTitle(Form("Cross Section Sector %d",s+1));
    v_mg_cs_logy[s]->GetXaxis()->SetTitle("#theta [deg]");
    v_mg_cs_logy[s]->GetXaxis()->CenterTitle();
    v_mg_cs_logy[s]->GetYaxis()->SetTitle("log( #sigma ) [mb/sr]");
    v_mg_cs_logy[s]->GetYaxis()->CenterTitle();
    v_mg_cs_logy[s]->Draw("APE");
    v_mg_cs_logy[s]->GetXaxis()->SetLimits(13.0, 22.0);
    v_mg_cs_logy[s]->GetHistogram()->SetMaximum(1.9); 
    v_mg_cs_logy[s]->GetHistogram()->SetMinimum(0.1); 
    v_mg_cs_logy[s]->Draw("APE");
    h_model->Draw("SAME C");
  
    TLegend *legend = new TLegend(0.6, 0.7, 0.89, 0.89);
    legend->AddEntry(g_cs_result[s],"Data");
    //legend->AddEntry(g_cs_model[s],"Model");
    legend->AddEntry(h_model,"Model");
    legend->SetBorderSize(0);
    legend->Draw();

  }
  c_result_log->SaveAs(Form("g_cs_result_log_r%d.pdf",run));
  
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // elastic cross section per phi bin 
  TCanvas *ca = new TCanvas("ca","ca",900,900);
  ca->Divide(12,6);
  for( int s = 0; s < cs_results_per_phi_bin.size(); s++ ){
    ca->cd(s+1);    
    h_model->SetLineColor(kRed);
    h_model->GetXaxis()->SetTitle("#theta [deg]");
    h_model->GetXaxis()->CenterTitle();    
    h_model->Draw();
    h_model->SetAxisRange(14,26,"X");
    h_model->Draw();
    cs_results_per_phi_bin[s]->Draw("same");
    cs_results_per_phi_bin[s]->SetAxisRange(14,26,"X");
    cs_results_per_phi_bin[s]->Draw("same");
  }
  ca->SaveAs(Form("cs_results_per_phi_bin_r%d.pdf",run));

  // get difference in model and measured
  std::vector<TGraphErrors*> g_result_diff_per_phi_bin;
  std::vector<TGraphErrors*> g_cs_result_per_phi_bin;
  std::vector<TGraphErrors*> g_cs_model_per_phi_bin;
  for( int pb = 0; pb<cs_results_per_phi_bin.size(); pb++ ){
    std::vector<double> bins;
    std::vector<double> bin_error_x;
    std::vector<double> bin_error_y;
    std::vector<double> model;
    std::vector<double> model_error;
    std::vector<double> data;
    std::vector<double> data_diff;
    std::vector<double> data_diff_err_x;   
    std::vector<double> data_diff_err_y;   
    std::cout << " comparing model and data for phi bin " << pb << std::endl;
    for( int b=1; b<=h_model->GetNbinsX(); b++){
      if ( b <= 14 || b >= 24 ) continue;
      double bin_center_model = h_model->GetBinCenter(b);
      double bin_center_data = cs_results_per_phi_bin[pb]->GetBinCenter(b);
      double model_result = h_model->GetBinContent(b);
      double data_result = cs_results_per_phi_bin[pb]->GetBinContent(b);
      double data_err_y = cs_results_per_phi_bin[pb]->GetBinError(b);
      double data_err_x = 1.0;
      double data_ratio = data_result/model_result;
      double data_ratio_err = data_ratio * data_err_y;

      // corrections
      double bin_center_correction = bin_center_corr[b-10]; // minus 14 + 6 because start histogram with bin center of 15.5 
      double radiative_correction = rad_corr[b-6]; //shift to the get the correct index 

      std::cout << " model result " << model_result << " data result " << data_result << " data error " << data_err_y << std::endl;
      std::cout << " center bin " << bin_center_data << "  bin centering corection  " << bin_center_correction << std::endl;
      std::cout << " center bin " << bin_center_data << " radiative correction " << radiative_correction << std::endl;

      double final_data_result = data_result/(bin_center_correction*radiative_correction);

      bins.push_back(bin_center_data);
      bin_error_x.push_back(0.0);
      bin_error_y.push_back(data_err_y);
      model_error.push_back(0.0);
      if ( model_result == 0.0 ) continue;
      model.push_back((model_result));
      data.push_back((final_data_result));
      data_diff.push_back( data_ratio );
      data_diff_err_y.push_back( data_ratio_err );			  
      data_diff_err_x.push_back( 0.0 );			        
      std::cout << " bin " << b << " bin center " << bin_center_model << " bin center data " << bin_center_data << " diff " << data_result/model_result << std::endl;
    }
    g_result_diff_per_phi_bin.push_back( new TGraphErrors(bins.size(), &(bins[0]), &(data_diff[0]), &(data_diff_err_x[0]), &(data_diff_err_y[0]) ) );
    g_cs_result_per_phi_bin.push_back( new TGraphErrors(bins.size(), &(bins[0]), &(data[0]), &(bin_error_x[0]), &(bin_error_y[0]) ) );
    g_cs_model_per_phi_bin.push_back( new TGraphErrors(bins.size(), &(bins[0]), &(model[0]), &(model_error[0]) ) );
  }

  TCanvas *c_diff_per_phi_bin = new TCanvas("c_diff_per_phi_bin","c_diff_per_phi_bin",400,600);
  c_diff_per_phi_bin->Divide(12,6);
  for( int pb = 0; pb < g_result_diff.size(); pb++ ){
    c_diff->cd(pb+1);
    g_result_diff_per_phi_bin[pb]->SetTitle(Form("Ratio of Cross Section Model to Data #phi Bin %d",pb+1));
    g_result_diff_per_phi_bin[pb]->GetXaxis()->SetTitle("#theta [deg]");
    g_result_diff_per_phi_bin[pb]->GetXaxis()->CenterTitle();
    g_result_diff_per_phi_bin[pb]->SetMarkerStyle(20);
    g_result_diff_per_phi_bin[pb]->SetMarkerSize(0.5);
    g_result_diff_per_phi_bin[pb]->Draw("AP");
    g_result_diff_per_phi_bin[pb]->GetXaxis()->SetLimits(10.0, 26.5);//SetRangeUser(0.0,26.5);
    g_result_diff_per_phi_bin[pb]->GetHistogram()->SetMaximum(1.5);   // along          
    g_result_diff_per_phi_bin[pb]->GetHistogram()->SetMinimum(0.50);  //   Y     
    g_result_diff_per_phi_bin[pb]->Draw("AP");
    c_diff_per_phi_bin->Update();
  }
  c_diff_per_phi_bin->SaveAs(Form("g_ratio_model_data_r%d_per_phi_bin.pdf",run));
  
  std::vector<TMultiGraph*> v_mg_cs_logy_per_phi_bin;

  TCanvas *c_result_log_per_phi_bin = new TCanvas("c_result_log_per_phi_bin","c_result_log_per_phi_bin",400,600);
  c_result_log_per_phi_bin->Divide(12,6);

  for( int pb = 0; pb < g_cs_result_per_phi_bin.size(); pb++ ){    
    std::cout << " plotting phi bi " << pb << std::endl;
    c_result_log_per_phi_bin->cd(pb+1);
    gPad->SetLogy();
    c_result_log_per_phi_bin->SetGrid();
    v_mg_cs_logy_per_phi_bin.push_back( new TMultiGraph() );
    g_cs_result_per_phi_bin[pb]->SetTitle(Form("Elastic Model Cross Section #phi Bin %d",pb+1));
    g_cs_result_per_phi_bin[pb]->GetXaxis()->SetTitle("#theta [deg]");
    g_cs_result_per_phi_bin[pb]->GetXaxis()->CenterTitle();
    g_cs_result_per_phi_bin[pb]->SetMarkerStyle(20);
    g_cs_model_per_phi_bin[pb]->SetMarkerStyle(22);
    //g_result_diff[s]->SetMarkerSize(3);
    //g_cs_result[s]->Draw("AP");

    //g_cs_model[s]->Draw("AP");
    g_cs_result_per_phi_bin[pb]->SetMarkerSize(0.5);
    g_cs_model_per_phi_bin[pb]->SetMarkerSize(0.5);
    g_cs_model_per_phi_bin[pb]->SetMarkerColor(kRed);
    v_mg_cs_logy_per_phi_bin[pb]->Add(g_cs_result_per_phi_bin[pb]);
    v_mg_cs_logy_per_phi_bin[pb]->Add(g_cs_model_per_phi_bin[pb]);
    v_mg_cs_logy_per_phi_bin[pb]->SetTitle(Form("Cross Section Sector %d",pb+1));
    v_mg_cs_logy_per_phi_bin[pb]->GetXaxis()->SetTitle("#theta [deg]");
    v_mg_cs_logy_per_phi_bin[pb]->GetXaxis()->CenterTitle();
    v_mg_cs_logy_per_phi_bin[pb]->GetYaxis()->SetTitle("log( #sigma ) [mb/sr]");
    v_mg_cs_logy_per_phi_bin[pb]->GetYaxis()->CenterTitle();
    v_mg_cs_logy_per_phi_bin[pb]->Draw("APE");
    v_mg_cs_logy_per_phi_bin[pb]->GetXaxis()->SetLimits(13.0, 22.0);
    v_mg_cs_logy_per_phi_bin[pb]->GetHistogram()->SetMaximum(1.9); 
    v_mg_cs_logy_per_phi_bin[pb]->GetHistogram()->SetMinimum(0.1); 
    v_mg_cs_logy_per_phi_bin[pb]->Draw("APE");
    h_model->Draw("SAME C");
  
    TLegend *legend = new TLegend(0.6, 0.7, 0.89, 0.89);
    legend->AddEntry(g_cs_result_per_phi_bin[pb],"Data");
    //legend->AddEntry(g_cs_model[s],"Model");
    legend->AddEntry(h_model,"Model");
    legend->SetBorderSize(0);
    legend->Draw();

  }
  c_result_log_per_phi_bin->SaveAs(Form("g_cs_result_log_r%d_per_phi_bin.pdf",run));
  

  readFromEModel.close();
  fOut->Write();



  return 0;
}
