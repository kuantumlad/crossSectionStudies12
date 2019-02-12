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
double Ebeam = 2.22193;

int ele_sector[BUFFER];
TLorentzVector ele[BUFFER]; 
TLorentzVector prot[BUFFER];
TLorentzVector mc_ele[BUFFER];
TLorentzVector mc_prot[BUFFER];

int select_ele;

double kin_W(TLorentzVector ele);
double kin_Q2(TLorentzVector ele);
double kin_x(TLorentzVector ele);
double kin_y(TLorentzVector ele);


int selector_mc(const char* inFile, const char* outputfile, int run ){

  Char_t *inTree="out_tree";
  Char_t *mcTree="mc_tree";
  
  TFile *f; 
  Char_t tmpstr[80];
  Double_t fraction;

  f = new TFile(inFile,"");   // Input File

  if(f->IsZombie()){   // Check if TFile exists!
    cout<<"Input file " << inFile << "doesn't exist!" << endl;
    cout<<"Exit program" << endl;
    return 0;
  }

  cout << "Reading from File: " << inFile << endl;

  TTree *anaTree=(TTree *) f->Get(inTree);
  TTree *mcanaTree=(TTree *) f->Get(mcTree);

  if(anaTree==0 || mcanaTree==0){  // Check if TTree exists!
    cout << "Tree " << inTree << " doesn't exist!!!" << endl;
    cout <<"Exit program" << endl;
    return 0;
  }

  double Pival = 3.141592658;

  vector<int> *p4_mc_pid = 0;
  vector<double> *p4_mc_px  = 0;
  vector<double> *p4_mc_py = 0;
  vector<double> *p4_mc_pz = 0;
  vector<double> *p4_mc_vx = 0;
  vector<double> *p4_mc_vy = 0;
  vector<double> *p4_mc_vz = 0;
  vector<int> *sectorE = 0;
  vector<double> *p4_ele_px = 0;
  vector<double> *p4_ele_py = 0;
  vector<double> *p4_ele_pz = 0;
  vector<double> *p4_ele_E = 0;
  vector<double> *p4_ele_vx = 0;
  vector<double> *p4_ele_vy = 0;
  vector<double> *p4_ele_vz = 0;
  vector<double> *p4_prot_px = 0;
  vector<double> *p4_prot_py = 0;
  vector<double> *p4_prot_pz = 0;
  vector<double> *p4_prot_E = 0;
  vector<double> *p4_prot_vx = 0;
  vector<double> *p4_prot_vy = 0;
  vector<double> *p4_prot_vz = 0;

  mcanaTree->SetBranchAddress("p4_mc_pid",&p4_mc_pid);
  mcanaTree->SetBranchAddress("p4_mc_px",&p4_mc_px);
  mcanaTree->SetBranchAddress("p4_mc_py",&p4_mc_py);
  mcanaTree->SetBranchAddress("p4_mc_pz",&p4_mc_pz);
  mcanaTree->SetBranchAddress("p4_mc_vx",&p4_mc_vx);
  mcanaTree->SetBranchAddress("p4_mc_vy",&p4_mc_vy);
  mcanaTree->SetBranchAddress("p4_mc_vz",&p4_mc_vz);
  anaTree->SetBranchAddress("sectorE",&sectorE);
  anaTree->SetBranchAddress("p4_ele_px",&p4_ele_px);
  anaTree->SetBranchAddress("p4_ele_py",&p4_ele_py);
  anaTree->SetBranchAddress("p4_ele_pz",&p4_ele_pz);
  anaTree->SetBranchAddress("p4_ele_E",&p4_ele_E);
  anaTree->SetBranchAddress("p4_ele_vx",&p4_ele_vx);
  anaTree->SetBranchAddress("p4_ele_vy",&p4_ele_vy);
  anaTree->SetBranchAddress("p4_ele_vz",&p4_ele_vz);
  anaTree->SetBranchAddress("p4_prot_px",&p4_prot_px);
  anaTree->SetBranchAddress("p4_prot_py",&p4_prot_py);
  anaTree->SetBranchAddress("p4_prot_pz",&p4_prot_pz);
  anaTree->SetBranchAddress("p4_prot_E",&p4_prot_E);
  anaTree->SetBranchAddress("p4_prot_vx",&p4_prot_vx);
  anaTree->SetBranchAddress("p4_prot_vy",&p4_prot_vy);
  anaTree->SetBranchAddress("p4_prot_vz",&p4_prot_vz);

  double el_theta_max = 30.0;
  double pr_theta_max = 50.0;

  TFile *mc_out = new TFile(outputfile,"RECREATE");

  TH1F *hist_all_electron_p; TH1F *hist_all_electron_theta; TH1F *hist_all_electron_phi;
  TH2F *hist_all_electron_p_vs_theta; TH2F *hist_all_electron_p_vs_phi; TH2F *hist_all_electron_theta_vs_phi;

  TH1F *hist_all_proton_p; TH1F *hist_all_proton_theta; TH1F *hist_all_proton_phi;
  TH2F *hist_all_proton_p_vs_theta; TH2F *hist_all_proton_p_vs_phi; TH2F *hist_all_proton_theta_vs_phi;

  hist_all_electron_p = new TH1F("hist_mc_all_electron_p", "electron momentum", 500,0,Ebeam+0.3);   
  hist_all_electron_p->GetXaxis()->SetTitle("p /GeV");
  hist_all_electron_p->GetYaxis()->SetTitle("counts");
  hist_all_electron_theta = new TH1F("hist_mc_all_electron_theta", "electron #Theta", 140,0,el_theta_max);   
  hist_all_electron_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_electron_theta->GetYaxis()->SetTitle("counts");
  hist_all_electron_phi = new TH1F("hist_mc_all_electron_phi", "electron #phi", 73,-180,180);   
  hist_all_electron_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_electron_phi->GetYaxis()->SetTitle("counts");
  hist_all_electron_p_vs_theta = new TH2F("hist_mc_all_electron_p_vs_theta", "electron p vs #Theta", 140,0,el_theta_max,500,0,Ebeam+0.3);   
  hist_all_electron_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_electron_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_all_electron_p_vs_phi = new TH2F("hist_all_electron_p_vs_phi", "electron p vs #phi", 180,-180,180, 500,0,Ebeam+0.3);   
  hist_all_electron_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_electron_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_all_electron_theta_vs_phi = new TH2F("hist_mc_all_electron_theta_vs_phi", "electron #Theta vs #phi", 73, -180,180, 30 , 0, el_theta_max);   
  hist_all_electron_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_electron_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  hist_all_proton_p = new TH1F("hist_mc_all_proton_p", "proton momentum", 500,0,Ebeam+0.3);   
  hist_all_proton_p->GetXaxis()->SetTitle("p /GeV");
  hist_all_proton_p->GetYaxis()->SetTitle("counts");
  hist_all_proton_theta = new TH1F("hist_mc_all_proton_theta", "proton #Theta", 140,0,el_theta_max);   
  hist_all_proton_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_proton_theta->GetYaxis()->SetTitle("counts");
  hist_all_proton_phi = new TH1F("hist_mc_all_proton_phi", "proton #phi", 180,-180,180);   
  hist_all_proton_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_proton_phi->GetYaxis()->SetTitle("counts");
  hist_all_proton_p_vs_theta = new TH2F("hist_mc_all_proton_p_vs_theta", "proton p vs #Theta", 140,0,pr_theta_max,500,0,Ebeam+0.3);   
  hist_all_proton_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_proton_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_all_proton_p_vs_phi = new TH2F("hist_mc_all_proton_p_vs_phi", "proton p vs #phi", 180,-180,180, 500,0,Ebeam+0.3);   
  hist_all_proton_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_proton_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_all_proton_theta_vs_phi = new TH2F("hist_mc_all_proton_theta_vs_phi", "proton #Theta vs #phi", 180,-180,180, 140,0,pr_theta_max);   
  hist_all_proton_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_proton_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  mc_out->mkdir("particle_histograms_selected");				
  mc_out->cd ("particle_histograms_selected");

  TH1F *hist_electron_p; TH1F *hist_electron_theta; TH1F *hist_electron_phi; 
  TH2F *hist_electron_p_vs_theta; TH2F *hist_electron_p_vs_phi; TH2F *hist_electron_theta_vs_phi;
  TH1F *hist_proton_p; TH1F *hist_proton_theta; TH1F *hist_proton_phi; 
  TH2F *hist_proton_p_vs_theta; TH2F *hist_proton_p_vs_phi; TH2F *hist_proton_theta_vs_phi;

  hist_electron_p = new TH1F("hist_mc_electron_p", "electron momentum", 500,0,Ebeam+0.3);   
  hist_electron_p->GetXaxis()->SetTitle("p /GeV");
  hist_electron_p->GetYaxis()->SetTitle("counts");
  hist_electron_theta = new TH1F("hist_mc_electron_theta", "electron #Theta", 140,0,40);   
  hist_electron_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_electron_theta->GetYaxis()->SetTitle("counts");
  hist_electron_phi = new TH1F("hist_mc_electron_phi", "electron #phi", 73,-180,180);   
  hist_electron_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_electron_phi->GetYaxis()->SetTitle("counts");
  hist_electron_p_vs_theta = new TH2F("hist_mc_electron_p_vs_theta", "electron p vs #Theta", 140,0,40,500,0,Ebeam+0.3);   
  hist_electron_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_electron_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_electron_p_vs_phi = new TH2F("hist_mc_electron_p_vs_phi", "electron p vs #phi", 180,-180,180, 500,0,Ebeam+0.3);   
  hist_electron_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_electron_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_electron_theta_vs_phi = new TH2F("hist_mc_electron_theta_vs_phi", "electron #Theta vs phi", 73,-180,180, 30,0,el_theta_max);   
  hist_electron_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_electron_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  hist_proton_p = new TH1F("hist_mc_proton_p", "proton momentum", 500,0,Ebeam+0.3);   
  hist_proton_p->GetXaxis()->SetTitle("p /GeV");
  hist_proton_p->GetYaxis()->SetTitle("counts");
  hist_proton_theta = new TH1F("hist_mc_proton_theta", "proton #Theta", 140,0, el_theta_max);   
  hist_proton_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_proton_theta->GetYaxis()->SetTitle("counts");
  hist_proton_phi = new TH1F("hist_mc_proton_phi", "proton #phi", 180,-180,180);   
  hist_proton_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_proton_phi->GetYaxis()->SetTitle("counts");
  hist_proton_p_vs_theta = new TH2F("hist_mc_proton_p_vs_theta", "proton p vs #Theta", 140,0,pr_theta_max,500,0,Ebeam+0.3);   
  hist_proton_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_proton_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_proton_p_vs_phi = new TH2F("hist_mc_proton_p_vs_phi", "proton p vs #phi", 180,-180,180, 500,0,Ebeam+0.3);   
  hist_proton_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_proton_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_proton_theta_vs_phi = new TH2F("hist_mc_proton_theta_vs_phi", "proton #Theta vs phi", 180,-180,180, 140,0,pr_theta_max);   
  hist_proton_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_proton_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");
  
  mc_out->mkdir("mc_kinematics");				
  mc_out->cd ("mc_kinematics");


  // Look at total phase space here as supplied to the detector
  TH1F *hist_W;
  TH1F *hist_Q2;
  TH1F *hist_x;
  TH1F *hist_y;
  TH1F *hist_nu;
  TH1F *hist_t;
  TH1F *hist_cmphi;
  TH1F *hist_cmcostheta;

  //ADDED 
  TH2F *hist_Q2_x;
  TH2F *hist_x_t;
  TH2F *hist_t_phi;
  TH2F *hist_Q2_phi;
  TH2F *hist_Q2_vs_W;
  TH2F *hist_W_vs_phi;

  std::vector< TH1F* > h_el_w_sect;
  std::vector< TH1F* > h_el_theta_sect;
  std::vector< TH2F* > h_el_q2w_sect;
  std::vector< TH2F* > h_el_ptheta_sect;

  std::vector< TH1F* > h_el_p_sect_final;
  std::vector< TH1F* > h_el_theta_sect_final;
  std::vector< TH1F* > h_el_phi_sect_final;

  std::vector< TH2F* > h_el_ptheta_sect_final;
  std::vector< TH2F* > h_el_phitheta_sect_final;

  TH2F *h_el_phitheta_final = new TH2F(Form("h_mc_el_phitheta_final_run%d",run), Form("h_mc_el_phitheta_final_run%d",run), 300, -180.0, 180.0, 200, 0.0, el_theta_max );

  for( int s = 0; s < 7; s++ ){
    h_el_w_sect.push_back( new TH1F(Form("h_mc_el_w_s%d",s),Form("h_mc_el_w_s%d",s), 100, 0.9, 1.2) );
    h_el_theta_sect.push_back( new TH1F(Form("h_mc_el_theta_s%d",s),Form("h_mc_el_theta_s%d",s), 100, 0.0, 40.0) );
   
    h_el_q2w_sect.push_back( new TH2F(Form("h_mc_el_q2w_s%d",s),Form("h_mc_el_q2w_s%d",s), 200, 0.9, 1.2, 200, 0.0, 2.5) );
    h_el_ptheta_sect.push_back( new TH2F(Form("h_mc_el_ptheta_s%d",s),Form("h_mc_el_ptheta_s%d",s), 200, 0.0, 2.5, 200, 0.0, 40.0) );
    
    h_el_p_sect_final.push_back( new TH1F(Form("h_mc_el_p_s%d_final",s),Form("h_mc_el_p_s%d_final",s),100, 0.0, 2.5) );
    h_el_theta_sect_final.push_back( new TH1F(Form("h_mc_el_theta_s%d_final",s),Form("h_mc_el_theta_s%d_final",s),30, 0.0, el_theta_max) );
    h_el_phi_sect_final.push_back( new TH1F(Form("h_mc_el_phi_s%d_final",s),Form("h_mc_el_phi_s%d_final",s), 73, -180.0, 180.0 ) );

    h_el_ptheta_sect_final.push_back( new TH2F(Form("h_mc_el_ptheta_s%d_final",s),Form("h_mc_el_ptheta_s%d_final",s), 200, 0.0, 2.5, 200, 0.0, el_theta_max) );
    h_el_phitheta_sect_final.push_back( new TH2F(Form("h_mc_el_phitheta_s%d_final",s),Form("h_mc_el_phitheta_s%d_final",s), 73, -180.0, 180.0, 30, 0.0, el_theta_max) );

  }

  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  int ele_count;
  int prot_count;
  int process_Events = -1;            // process all events

  for(Int_t k=0; k < mcanaTree->GetEntriesFast();k++){    
      
    mcanaTree->GetEntry(k);
    
    if(process_Events > 0 && k == process_Events) break;

    /// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// progress:

    if(k % 100000 == 0){

      double events = mcanaTree->GetEntriesFast();
      double percent = k/(events/100);
      
      printf("Analysing event number %i of %.00f (%.01f percent)\n", k, events, percent);
    }

    for( int i = 0; i < p4_mc_pid->size(); i++ ){
      //std::cout << " pid " << p4_mc_pid->at(i) << std::endl;
      int mc_pid = p4_mc_pid->at(i); 
      if( mc_pid == 11 ){ 
	double p4_mc_E = sqrt(p4_mc_px->at(i)*p4_mc_px->at(i) + p4_mc_py->at(i)*p4_mc_py->at(i) + p4_mc_pz->at(i)*p4_mc_pz->at(i) + 0.000511*0.000511);
	mc_ele[i].SetPxPyPzE(p4_mc_px->at(i),p4_mc_py->at(i),p4_mc_pz->at(i),p4_mc_E);
      }
      if( mc_pid == 2212 ){ 
	double p4_mc_E = sqrt(p4_mc_px->at(i)*p4_mc_px->at(i) + p4_mc_py->at(i)*p4_mc_py->at(i) + p4_mc_pz->at(i)*p4_mc_pz->at(i) + 0.938*0.938);      
	mc_prot[i].SetPxPyPzE(p4_mc_px->at(i),p4_mc_py->at(i),p4_mc_pz->at(i),p4_mc_E);
      }

      if( mc_ele[i].Theta()*180/Pival > 5.0 ){ //ignore events in the central detector < 5 deg
	if(mc_ele[i].P() > 0) hist_all_electron_p->Fill(mc_ele[i].P());
	if(mc_ele[i].Theta() > 0) hist_all_electron_theta->Fill(mc_ele[i].Theta()*180/Pival);
	if(mc_ele[i].Phi() != 0) hist_all_electron_phi->Fill(mc_ele[i].Phi()*180/Pival);
	if(mc_ele[i].P() > 0) hist_all_electron_p_vs_theta->Fill(mc_ele[i].Theta()*180/Pival, mc_ele[i].P());
	if(mc_ele[i].P() > 0) hist_all_electron_p_vs_phi->Fill(mc_ele[i].Phi()*180/Pival, mc_ele[i].P());
	if(mc_ele[i].Theta() > 0 && mc_ele[i].Phi() != 0) hist_all_electron_theta_vs_phi->Fill(mc_ele[i].Phi()*180/Pival, mc_ele[i].Theta()*180/Pival);
      }
      
    }
    
    
    
  }




  mc_out->Write();
  mc_out->Close();
  f->Close();
  std::cout << " complete " << std::endl;
  return 0;
}


double kin_W(TLorentzVector ele){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  TLorentzVector fCM = fGamma + target;
  return fCM.M();
}

double kin_Q2(TLorentzVector ele){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  return -fGamma.M2();
}

double kin_x(TLorentzVector ele){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  return -fGamma.M2()/(2*0.938272*fGamma.E());
}


double kin_y(TLorentzVector ele){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  return fGamma.E()/Ebeam;
}

double kin_nu(TLorentzVector ele){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  return fGamma.E();
}

double kin_pT(TLorentzVector ele, TLorentzVector hadron){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  TLorentzVector fCM = fGamma + target;
  TVector3 boost = -fCM.BoostVector();
  TLorentzVector hadronClone(hadron); 
  hadronClone.Boost(boost);
  return hadronClone.Perp(); 
}

double kin_eta(TLorentzVector ele, TLorentzVector hadron){ 
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  TLorentzVector fCM = fGamma + target;
  TVector3 boost = -fCM.BoostVector();
  TLorentzVector hadronClone(hadron); 
  hadronClone.Boost(boost);
  return hadronClone.Rapidity(); 
}

double kin_z(TLorentzVector ele, TLorentzVector hadron){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;  
  return hadron.T()/fGamma.T();
}

double kin_t(TLorentzVector ele, TLorentzVector hadron){

  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  TLorentzVector transfer = fGamma - hadron;
  return transfer.M2();
}

double kin_cmphi(TLorentzVector ele, TLorentzVector hadron){

  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);

  TLorentzVector fGamma = beam - ele;
  TLorentzVector fCM = fGamma + target;

  TVector3 boost = -fCM.BoostVector();  // get boost vector to CM frame

  TLorentzVector gammaClone(fGamma);    // clone particles
  TLorentzVector eleClone(ele); 
  TLorentzVector hadronClone(hadron); 

  gammaClone.Boost(boost);     // boost particles to CM frame
  eleClone.Boost(boost);
  hadronClone.Boost(boost);

  TVector3 v0, v1;
  v0 = gammaClone.Vect().Cross(eleClone.Vect());
  v1 = gammaClone.Vect().Cross(hadronClone.Vect());    // hadronClone.Vect().Cross(gammaClone.Vect());

  Double_t c0, c1, c2, c3;
  c0 = v0.Dot(hadronClone.Vect());
  c1 = v0.Dot(v1);
  c2 = v0.Mag();
  c3 = v1.Mag();

  return (c0/TMath::Abs(c0)) * TMath::ACos(c1 /(c2*c3));
}

double kin_cmcostheta(TLorentzVector ele, TLorentzVector hadron){

  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);

  TLorentzVector fGamma = beam - ele;
  TLorentzVector fCM = fGamma + target;

  TVector3 boost = -fCM.BoostVector();

  TLorentzVector hadronClone(hadron); 
  hadronClone.Boost(boost);

  return TMath::Cos(hadronClone.Theta()); 
}


/////////////////////////////////////////////////////////////////////
///  missing mass and energy

TLorentzVector miss_X(TLorentzVector ele, TLorentzVector hadron){

  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector miss   = beam + target - ele - hadron;
  return miss;
}

double kin_mismass(TLorentzVector ele, TLorentzVector hadron){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector miss   = beam + target - ele - hadron;
  return sqrt(miss.E()*miss.E() - miss.Px()*miss.Px() - miss.Py()*miss.Py() - miss.Pz()*miss.Pz());
}

double kin_mismass2(TLorentzVector ele, TLorentzVector hadron){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector miss   = beam + target - ele - hadron;
  return (miss.E()*miss.E() - miss.Px()*miss.Px() - miss.Py()*miss.Py() - miss.Pz()*miss.Pz());
}
