/// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///  ROOT Macro to select physics reactions from filtered clas 6 and clas 12 data  (simple version for inklusive electron scattering)
///
///  Stefan Diehl  (sdiehl@jlab.org)
///
///  To execute run ROOT and compile macro with .L selector_simple.cxx++
///
///  Then call function:  selector_simple("output/2391_junemap.root", "output/selected_2391_junemap.root")
///   
/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h> 
#include "Riostream.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TMath.h"
#include <vector>
#include <TAxis.h>
#include <TLorentzRotation.h>
#include<vector>
#include <algorithm>
#include <functional>
#include <TCutG.h>
using namespace std;

#include "Math/GenVector/PxPyPzE4D.h"
#include "Math/GenVector/PtEtaPhiE4D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

using namespace ROOT::Math;

/// //////////////////////////////////////////////////////////////////////////////////////////////////////
/// settings:

//double Ebeam = 6.42313;
//double Ebeam = 2.22193;
double Ebeam = 10.594;

int process_Events = -1;            // process all events
//int process_Events = 500000;     // process given number of events


/// /////////////////////////////////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////////////////////////////////
/// common variables
///

Float_t Pival = 3.14159265359;

Float_t m_e = 0.000511;
Float_t m_p = 0.93827;
Float_t m_n = 0.9396;
Float_t m_pip = 0.1396;
Float_t m_pim = 0.1396;
Float_t m_Kp = 0.4937;
Float_t m_Km = 0.4937;
Float_t c = 299792458;

// vectors with components of the Lorentzvector to fill into the tree

int helicity;
double fcup;
vector<double> *p4_ele_px = 0;
vector<double> *p4_ele_py = 0;
vector<double> *p4_ele_pz = 0;
vector<double> *p4_ele_E = 0;
vector<double> *p4_prot_px = 0;
vector<double> *p4_prot_py = 0;
vector<double> *p4_prot_pz = 0;
vector<double> *p4_prot_E = 0;
vector<double> *p4_neutr_px = 0;
vector<double> *p4_neutr_py = 0;
vector<double> *p4_neutr_pz = 0;
vector<double> *p4_neutr_E = 0;
vector<double> *p4_pip_px = 0;
vector<double> *p4_pip_py = 0;
vector<double> *p4_pip_pz = 0;
vector<double> *p4_pip_E = 0;
vector<double> *p4_pim_px = 0;
vector<double> *p4_pim_py = 0;
vector<double> *p4_pim_pz = 0;
vector<double> *p4_pim_E = 0;
vector<double> *p4_Kp_px = 0;
vector<double> *p4_Kp_py = 0;
vector<double> *p4_Kp_pz = 0;
vector<double> *p4_Kp_E = 0;
vector<double> *p4_Km_px = 0;
vector<double> *p4_Km_py = 0;
vector<double> *p4_Km_pz = 0;
vector<double> *p4_Km_E = 0;
vector<double> *p4_phot_px = 0;
vector<double> *p4_phot_py = 0;
vector<double> *p4_phot_pz = 0;
vector<double> *p4_phot_E = 0;
vector<double> *prot_beta_final = 0;
vector<double> *Kp_beta_final = 0;
vector<double> *Km_beta_final = 0;
vector<double> *el_sector = 0;
vector<double> *prot_sector = 0;
vector<double> *Kp_sector = 0;
vector<double> *Km_sector = 0;

/// varibales for particles:

const static int BUFFER = 20;

TLorentzVector ele[BUFFER]; 
TLorentzVector prot[BUFFER]; 
TLorentzVector neutr[BUFFER]; 
TLorentzVector pip[BUFFER]; 
TLorentzVector pim[BUFFER]; 
TLorentzVector Kp[BUFFER]; 
TLorentzVector Km[BUFFER]; 
TLorentzVector phot[BUFFER];


//added for beta vs p
double prot_beta[BUFFER];
double kaonP_beta[BUFFER];
double kaonM_beta[BUFFER];

int sect_el[BUFFER];
int sect_pr[BUFFER];
int sect_kp[BUFFER];
int sect_km[BUFFER];

///  counting variables:

double ele_count;
double prot_count;
double neutr_count;
double pip_count;
double pim_count;
double Kp_count;
double Km_count;
double phot_count;


// kinematic variables
double W, Q2, nu, x, y;
double t1, cmphi1, cmcostheta1, pt1, eta1, z1;
double M_e_p_X_miss, M_e_p_X_miss2;

double W_out, Q2_out, nu_out, x_out, y_out;
double t1_out, cmphi1_out, cmcostheta1_out, pt1_out, eta1_out, z1_out;
double perp_mntm, missing_e, missing_mm2, epX_mm2, epX_mm;
 

double W_2, Q2_2, nu_2, x_2, y_2;
double t1_2, cmphi1_2, cmcostheta1_2, pt1_2, eta1_2, z1_2;


///////////////////////////////////////////////////////////////////////////////////////////////////////
///  selected particles:

Int_t evcat; 
Double_t E_ele; Double_t px_ele; Double_t py_ele; Double_t pz_ele;
Double_t E_prot_1; Double_t px_prot_1; Double_t py_prot_1; Double_t pz_prot_1;

//ADDED 
Double_t E_kaonP_1; Double_t px_kaonP_1; Double_t py_kaonP_1; Double_t pz_kaonP_1;
Double_t E_kaonM_1; Double_t px_kaonM_1; Double_t py_kaonM_1; Double_t pz_kaonM_1;

Double_t prot_beta_out; Double_t kaonP_beta_out; Double_t kaonM_beta_out;
Int_t el_sect; Int_t pr_sect; Int_t kp_sect; Int_t km_sect;

Int_t channel;

/// ////////////////////////////////////////////////////////////////////////////////////////////////////
///  index of selected particle:

int select_ele;
int select_prot_1; 

//NEW
int select_kaonP_1;
int select_kaonM_1;

/// ////////////////////////////////////////////////////////////////////////////////////////////////////
///  input and output files:

TFile *out;

char name[200];
char title[200];

/// /////////////////////////////////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////////////////////////////////
/// define functions

double kin_W(TLorentzVector ele);
double kin_Q2(TLorentzVector ele);
double kin_x(TLorentzVector ele);
double kin_y(TLorentzVector ele);
double kin_nu(TLorentzVector ele);

double kin_pT(TLorentzVector ele, TLorentzVector hadron);
double kin_eta(TLorentzVector ele, TLorentzVector hadron);
double kin_z(TLorentzVector ele, TLorentzVector hadron);
double kin_t(TLorentzVector ele, TLorentzVector hadron);

double kin_cmphi(TLorentzVector ele, TLorentzVector hadron);
double kin_cmcostheta(TLorentzVector ele, TLorentzVector hadron);

double kin_mismass(TLorentzVector ele, TLorentzVector hadron);
double kin_mismass2(TLorentzVector ele, TLorentzVector hadron);

TLorentzVector miss_X(TLorentzVector ele, TLorentzVector hadron);

double alpha_e_X(TLorentzVector ele, TLorentzVector hadron);
double alpha_p1p2(TLorentzVector particle1, TLorentzVector particle2);

double Phi_Mass(TLorentzVector kaonP, TLorentzVector kaonM);

//ADDED 
double kin_epkXMass(TLorentzVector el, TLorentzVector pr, TLorentzVector kp);
double kin_epKpKmXMass(TLorentzVector el, TLorentzVector pr, TLorentzVector kp, TLorentzVector km);



/// /////////////////////////////////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// main
///

Int_t selector_simple( Char_t *inFile, Char_t *outputfile)
{
		
  Char_t *inTree="out_tree";
  cout << "Initalize the input tree ... " << endl;

  /// ///////////////////////////////////////////////////////////////////////////
  ///  get input file:
  /// ///////////////////////////////////////////////////////////////////////////

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

  if(anaTree==0){  // Check if TTree exists!
    cout << "Tree " << inTree << " doesn't exist!!!" << endl;
    cout <<"Exit program" << endl;
    return 0;
  }
    
    
  /// /////////////////////////////////////////////////////////////////////////////    
  ///  get branches from input file:
  /// ///////////////////////////////////////////////////////////////////////////// 
    
  anaTree->SetBranchAddress("helicity", &helicity);
  anaTree->SetBranchAddress("fcup", &fcup);
  anaTree->SetBranchAddress("p4_ele_px", &p4_ele_px);
  anaTree->SetBranchAddress("p4_ele_py", &p4_ele_py);
  anaTree->SetBranchAddress("p4_ele_pz", &p4_ele_pz);
  anaTree->SetBranchAddress("p4_ele_E", &p4_ele_E);
  anaTree->SetBranchAddress("p4_prot_px", &p4_prot_px);
  anaTree->SetBranchAddress("p4_prot_py", &p4_prot_py);
  anaTree->SetBranchAddress("p4_prot_pz", &p4_prot_pz);
  anaTree->SetBranchAddress("p4_prot_E", &p4_prot_E);
  anaTree->SetBranchAddress("p4_neutr_px", &p4_neutr_px);
  anaTree->SetBranchAddress("p4_neutr_py", &p4_neutr_py);
  anaTree->SetBranchAddress("p4_neutr_pz", &p4_neutr_pz);
  anaTree->SetBranchAddress("p4_neutr_E", &p4_neutr_E);
  anaTree->SetBranchAddress("p4_pip_px", &p4_pip_px);
  anaTree->SetBranchAddress("p4_pip_py", &p4_pip_py);
  anaTree->SetBranchAddress("p4_pip_pz", &p4_pip_pz);
  anaTree->SetBranchAddress("p4_pip_E", &p4_pip_E);
  anaTree->SetBranchAddress("p4_pim_px", &p4_pim_px);
  anaTree->SetBranchAddress("p4_pim_py", &p4_pim_py);
  anaTree->SetBranchAddress("p4_pim_pz", &p4_pim_pz);
  anaTree->SetBranchAddress("p4_pim_E", &p4_pim_E);
  anaTree->SetBranchAddress("p4_Kp_px", &p4_Kp_px);
  anaTree->SetBranchAddress("p4_Kp_py", &p4_Kp_py);
  anaTree->SetBranchAddress("p4_Kp_pz", &p4_Kp_pz);
  anaTree->SetBranchAddress("p4_Kp_E", &p4_Kp_E);
  anaTree->SetBranchAddress("p4_Km_px", &p4_Km_px);
  anaTree->SetBranchAddress("p4_Km_py", &p4_Km_py);
  anaTree->SetBranchAddress("p4_Km_pz", &p4_Km_pz);
  anaTree->SetBranchAddress("p4_Km_E", &p4_Km_E);
  anaTree->SetBranchAddress("p4_phot_px", &p4_phot_px);
  anaTree->SetBranchAddress("p4_phot_py", &p4_phot_py);
  anaTree->SetBranchAddress("p4_phot_pz", &p4_phot_pz);
  anaTree->SetBranchAddress("p4_phot_E", &p4_phot_E);  
  /*  anaTree->SetBranchAddress("prot_beta_final",&prot_beta_final);
  anaTree->SetBranchAddress("Kp_beta_final",&Kp_beta_final);
  anaTree->SetBranchAddress("Km_beta_final",&Km_beta_final);
  anaTree->SetBranchAddress("el_sector",&el_sector);
  anaTree->SetBranchAddress("prot_sector",&prot_sector);
  anaTree->SetBranchAddress("Kp_sector",&Kp_sector);
  anaTree->SetBranchAddress("Km_sector",&Km_sector);
  anaTree->SetBranchAddress("el_sector",&el_sector);
  anaTree->SetBranchAddress("prot_sector",&prot_sector);
  anaTree->SetBranchAddress("Kp_sector",&Kp_sector);
  anaTree->SetBranchAddress("Km_sector",&Km_sector);
  */
  /// /////////////////////////////////////////////////////////////////////////////    
  ///  create output tree:
  /// /////////////////////////////////////////////////////////////////////////////

  out = new TFile(outputfile, "RECREATE");

  // categroy 1:  e p -->  e p X

  TTree out_tree1("events_epX","events_epX");
  
  out_tree1.Branch("W", &W);
  out_tree1.Branch("Q2", &Q2);
  out_tree1.Branch("x", &x);
  out_tree1.Branch("y", &y);
  out_tree1.Branch("nu", &nu);
  out_tree1.Branch("minus_t", &t1);
  out_tree1.Branch("cmphi", &cmphi1);
  out_tree1.Branch("cmcostheta", &cmcostheta1);
  out_tree1.Branch("pt", &pt1);
  out_tree1.Branch("eta", &eta1);
  out_tree1.Branch("z", &z1);
  out_tree1.Branch("missmass", &M_e_p_X_miss);
  out_tree1.Branch("missmass2", &M_e_p_X_miss2);
  out_tree1.Branch("helicity", &helicity);
  out_tree1.Branch("fcup", &fcup);
  out_tree1.Branch("E_ele", &E_ele);
  out_tree1.Branch("px_ele", &px_ele);
  out_tree1.Branch("py_ele", &py_ele);
  out_tree1.Branch("pz_ele", &pz_ele);
  out_tree1.Branch("E_prot", &E_prot_1);
  out_tree1.Branch("px_prot", &px_prot_1);
  out_tree1.Branch("py_prot", &py_prot_1);
  out_tree1.Branch("pz_prot", &pz_prot_1);
    
  out_tree1.Branch("channel",&channel);

  /// ///////////////////////////////////////////////////////////////
  ///  create histograms:
  /// ///////////////////////////////////////////////////////////////

  out->mkdir("event_information");				
  out->cd ("event_information");

  TH1F *hist_helicity;
  TH1F *hist_faraday_cup;

  hist_helicity = new TH1F("hist_helicity", "helicity", 5, -2.5, 2.5);   
  hist_helicity->GetXaxis()->SetTitle("helicity");
  hist_helicity->GetYaxis()->SetTitle("counts");
  hist_faraday_cup = new TH1F("hist_faraday_cup", "faraday cup", 200, 0, 200);   
  hist_faraday_cup->GetXaxis()->SetTitle("faraday cup");
  hist_faraday_cup->GetYaxis()->SetTitle("counts");

  out->mkdir("particle_histograms_all");				
  out->cd ("particle_histograms_all");

  TH1F *hist_all_electron_p; TH1F *hist_all_electron_theta; TH1F *hist_all_electron_phi;
  TH2F *hist_all_electron_p_vs_theta; TH2F *hist_all_electron_p_vs_phi; TH2F *hist_all_electron_theta_vs_phi;

  TH1F *hist_all_proton_p; TH1F *hist_all_proton_theta; TH1F *hist_all_proton_phi;
  TH2F *hist_all_proton_p_vs_theta; TH2F *hist_all_proton_p_vs_phi; TH2F *hist_all_proton_theta_vs_phi;

  hist_all_electron_p = new TH1F("hist_all_electron_p", "electron momentum", 500,0,Ebeam+0.3);   
  hist_all_electron_p->GetXaxis()->SetTitle("p /GeV");
  hist_all_electron_p->GetYaxis()->SetTitle("counts");
  hist_all_electron_theta = new TH1F("hist_all_electron_theta", "electron #Theta", 140,0,140);   
  hist_all_electron_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_electron_theta->GetYaxis()->SetTitle("counts");
  hist_all_electron_phi = new TH1F("hist_all_electron_phi", "electron #phi", 180,-180,180);   
  hist_all_electron_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_electron_phi->GetYaxis()->SetTitle("counts");
  hist_all_electron_p_vs_theta = new TH2F("hist_all_electron_p_vs_theta", "electron p vs #Theta", 140,0,140,500,0,Ebeam+0.3);   
  hist_all_electron_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_electron_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_all_electron_p_vs_phi = new TH2F("hist_all_electron_p_vs_phi", "electron p vs #phi", 180,-180,180, 500,0,Ebeam+0.3);   
  hist_all_electron_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_electron_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_all_electron_theta_vs_phi = new TH2F("hist_all_electron_theta_vs_phi", "electron #Theta vs #phi", 180,-180,180, 140,0,140);   
  hist_all_electron_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_electron_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  hist_all_proton_p = new TH1F("hist_all_proton_p", "proton momentum", 500,0,Ebeam+0.3);   
  hist_all_proton_p->GetXaxis()->SetTitle("p /GeV");
  hist_all_proton_p->GetYaxis()->SetTitle("counts");
  hist_all_proton_theta = new TH1F("hist_all_proton_theta", "proton #Theta", 140,0,140);   
  hist_all_proton_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_proton_theta->GetYaxis()->SetTitle("counts");
  hist_all_proton_phi = new TH1F("hist_all_proton_phi", "proton #phi", 180,-180,180);   
  hist_all_proton_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_proton_phi->GetYaxis()->SetTitle("counts");
  hist_all_proton_p_vs_theta = new TH2F("hist_all_proton_p_vs_theta", "proton p vs #Theta", 140,0,140,500,0,Ebeam+0.3);   
  hist_all_proton_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_proton_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_all_proton_p_vs_phi = new TH2F("hist_all_proton_p_vs_phi", "proton p vs #phi", 180,-180,180, 500,0,Ebeam+0.3);   
  hist_all_proton_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_proton_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_all_proton_theta_vs_phi = new TH2F("hist_all_proton_theta_vs_phi", "proton #Theta vs #phi", 180,-180,180, 140,0,140);   
  hist_all_proton_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_proton_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  out->mkdir("particle_histograms_selected");				
  out->cd ("particle_histograms_selected");

  TH1F *hist_electron_p; TH1F *hist_electron_theta; TH1F *hist_electron_phi; 
  TH2F *hist_electron_p_vs_theta; TH2F *hist_electron_p_vs_phi; TH2F *hist_electron_theta_vs_phi;
  TH1F *hist_proton_p; TH1F *hist_proton_theta; TH1F *hist_proton_phi; 
  TH2F *hist_proton_p_vs_theta; TH2F *hist_proton_p_vs_phi; TH2F *hist_proton_theta_vs_phi;

  hist_electron_p = new TH1F("hist_electron_p", "electron momentum", 500,0,Ebeam+0.3);   
  hist_electron_p->GetXaxis()->SetTitle("p /GeV");
  hist_electron_p->GetYaxis()->SetTitle("counts");
  hist_electron_theta = new TH1F("hist_electron_theta", "electron #Theta", 140,0,140);   
  hist_electron_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_electron_theta->GetYaxis()->SetTitle("counts");
  hist_electron_phi = new TH1F("hist_electron_phi", "electron #phi", 180,-180,180);   
  hist_electron_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_electron_phi->GetYaxis()->SetTitle("counts");
  hist_electron_p_vs_theta = new TH2F("hist_electron_p_vs_theta", "electron p vs #Theta", 140,0,140,500,0,Ebeam+0.3);   
  hist_electron_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_electron_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_electron_p_vs_phi = new TH2F("hist_electron_p_vs_phi", "electron p vs #phi", 180,-180,180, 500,0,Ebeam+0.3);   
  hist_electron_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_electron_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_electron_theta_vs_phi = new TH2F("hist_electron_theta_vs_phi", "electron #Theta vs phi", 180,-180,180, 140,0,140);   
  hist_electron_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_electron_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  hist_proton_p = new TH1F("hist_proton_p", "proton momentum", 500,0,Ebeam+0.3);   
  hist_proton_p->GetXaxis()->SetTitle("p /GeV");
  hist_proton_p->GetYaxis()->SetTitle("counts");
  hist_proton_theta = new TH1F("hist_proton_theta", "proton #Theta", 140,0,140);   
  hist_proton_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_proton_theta->GetYaxis()->SetTitle("counts");
  hist_proton_phi = new TH1F("hist_proton_phi", "proton #phi", 180,-180,180);   
  hist_proton_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_proton_phi->GetYaxis()->SetTitle("counts");
  hist_proton_p_vs_theta = new TH2F("hist_proton_p_vs_theta", "proton p vs #Theta", 140,0,140,500,0,Ebeam+0.3);   
  hist_proton_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_proton_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_proton_p_vs_phi = new TH2F("hist_proton_p_vs_phi", "proton p vs #phi", 180,-180,180, 500,0,Ebeam+0.3);   
  hist_proton_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_proton_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_proton_theta_vs_phi = new TH2F("hist_proton_theta_vs_phi", "proton #Theta vs phi", 180,-180,180, 140,0,140);   
  hist_proton_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_proton_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");
  
  out->mkdir("kinematics");				
  out->cd ("kinematics");


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

  std::vector<TH1F*> h_w_cut;
  std::vector<TH1F*> h_q2_cut;
  std::vector<TH1F*> h_x_cut;
  std::vector<TH1F*> h_y_cut;
  std::vector<TH1F*> h_nu_cut;
  std::vector<TH1F*> h_t_cut;
  std::vector<TH1F*> h_cmphi_cut;
  std::vector<TH1F*> h_cmcostheta_cut;

  std::vector<TH2F*> h_q2_x_cut;
  std::vector<TH2F*> h_x_t_cut;
  std::vector<TH2F*> h_t_phi_cut;
  std::vector<TH2F*> h_q2_phi_cut;
  std::vector<TH2F*> h_q2_w_cut;
  std::vector<TH2F*> h_w_phi_cut;
  std::vector<TH2F*> h_q2_t_cut;

  for( int c = 0; c < 12; c++ ){
    h_w_cut.push_back(new TH1F(Form("h_w_cutlvl_%d",c), Form("h_w_cutlvl_%d",c), 100, 0.0, 5.0));
    h_q2_cut.push_back(new TH1F(Form("h_q2_cutlvl_%d",c), Form("h_q2_cutlvl_%d",c), 100, 0.0, 10.0));
    h_x_cut.push_back(new TH1F(Form("h_x_cutlvl_%d",c), Form("h_x_cutlvl_%d",c), 100, 0.0, 1.10));
    h_y_cut.push_back(new TH1F(Form("h_y_cutlvl_%d",c), Form("h_y_cutlvl_%d",c), 100, 0.0, 1.0));
    h_nu_cut.push_back(new TH1F(Form("h_nu_cutlvl_%d",c), Form("h_nu_cutlvl_%d",c), 100, 0.0, 10.0));
    h_t_cut.push_back(new TH1F(Form("h_t_cutlvl_%d",c), Form("h_t_cutlvl_%d",c), 100, 0.0, 5.0));
    h_cmphi_cut.push_back(new TH1F(Form("h_cmphi_cutlvl_%d",c), Form("h_cmphi_cutlvl_%d",c), 200, -180.0, 180.0));
    h_cmcostheta_cut.push_back(new TH1F(Form("h_cmcostheta_cutlvl_%d",c), Form("h_cmcostheta_cutlvl_%d",c), 100, -1.0, 1.0));


    h_q2_x_cut.push_back( new TH2F(Form("h_q2_x_cutlvl_%d",c), Form("h_q2_x_cutlvl_%d",c), 200, 0.0, 1.1, 200, 0.0, 10.5));
    h_x_t_cut.push_back( new TH2F(Form("h_x_t_cutlvl_%d",c), Form("h_x_t_cutlvl_%d",c), 200, 0.0, 1.1, 200, 0.0, 5.0));
    h_t_phi_cut.push_back( new TH2F(Form("h_t_phi_cutlvl_%d",c), Form("h_t_phi_cutlvl_%d",c), 200, -180.0, 180.0, 200, 0.0, 5.0));
    h_q2_phi_cut.push_back( new TH2F(Form("h_q2_phi_cutlvl_%d",c), Form("h_q2_phi_cutlvl_%d",c), 200, -180.0, 180.0, 200, 0.0, 10.5));
    h_q2_w_cut.push_back( new TH2F(Form("h_q2_w_cutlvl_%d",c), Form("h_q2_w_cutlvl_%d",c), 200, 1.0, 5.0, 200, 0.0, 12.0));
    h_w_phi_cut.push_back( new TH2F(Form("h_w_phi_cutlvl_%d",c), Form("h_w_phi_cutlvl_%d",c), 200, -180.0, 180.0, 200, 0.0, 5.0));
    h_q2_t_cut.push_back( new TH2F(Form("h_q2_t_cutlvl_%d",c), Form("h_q2_t_cutlvl_%d",c), 200, 0.0, 12.5, 200, 0.0, 10.0));
    
  }

  
  TH2F *h_el_ptheta;
  TH2F *h_el_pphi;
  TH2F *h_el_phitheta;

  TH2F *h_pr_ptheta;
  TH2F *h_pr_pphi;
  TH2F *h_pr_phitheta;

  h_el_ptheta = new TH2F("h_el_ptheta","Electron #theta vs p",75, 0.0, Ebeam, 75, 0.0, 65 );
  h_el_pphi = new TH2F("h_el_pphi","Electron #phi vs p",75, 0.0, Ebeam, 200, -180.0, 180.0 );
  h_el_phitheta = new TH2F("h_el_phitheta","Electron #phi vs #theta",200, -180.0, 180.0, 75, 0.0, 65 );

  h_pr_ptheta = new TH2F("h_pr_ptheta","Proton #theta vs p",75, 0.0, Ebeam, 95, 0.0, 120 );
  h_pr_pphi = new TH2F("h_pr_pphi","Proton #phi vs p",75, 0.0, Ebeam, 150, -180.0, 180.0 );
  h_pr_phitheta = new TH2F("h_pr_phitheta","Proton #phi vs #theta", 200, -180.0, 180.0, 95, 0.0, 120 );


  TH2F *h_el_ptheta_final;
  TH2F *h_el_pphi_final;
  TH2F *h_el_phitheta_final;

  TH2F *h_pr_ptheta_final;
  TH2F *h_pr_pphi_final;
  TH2F *h_pr_phitheta_final;

  h_el_ptheta_final = new TH2F("h_el_ptheta_final","Electron #theta vs p",75, 0.0, Ebeam, 75, 0.0, 65 );
  h_el_pphi_final = new TH2F("h_el_pphi_final","Electron #phi vs p",75, 0.0, Ebeam, 200, -180.0, 180.0 );
  h_el_phitheta_final = new TH2F("h_el_phitheta_final","Electron #phi vs #theta",200, -180.0, 180.0, 75, 0.0, 65 );

  h_pr_ptheta_final = new TH2F("h_pr_ptheta_final","Proton #theta vs p",75, 0.0, Ebeam, 95, 0.0, 120 );
  h_pr_pphi_final = new TH2F("h_pr_pphi_final","Proton #phi vs p",75, 0.0, Ebeam, 150, -180.0, 180.0 );
  h_pr_phitheta_final = new TH2F("h_pr_phitheta_final","Proton #phi vs #theta", 200, -180.0, 180.0, 95, 0.0, 120 );
 
  hist_W = new TH1F("W", "W all sectors", 700, 0.5, Ebeam-2);   
  hist_W->GetXaxis()->SetTitle("W /GeV");
  hist_W->GetYaxis()->SetTitle("counts");

  hist_Q2 = new TH1F("Q2", "Q^{2} all sectors", 800, 0, Ebeam+2);   
  hist_Q2->GetXaxis()->SetTitle("Q^{2} /GeV^{2}");
  hist_Q2->GetYaxis()->SetTitle("counts");

  hist_Q2_vs_W = new TH2F("Q2_vs_W", "Q^{2} vs W all sectors", 700, 0.5, Ebeam-2, 800, 0, Ebeam+2);   
  hist_Q2_vs_W->GetXaxis()->SetTitle("W /GeV");
  hist_Q2_vs_W->GetYaxis()->SetTitle("Q^{2}");

  hist_W_vs_phi = new TH2F("W_vs_phi", "W vs phi all sectors", 180, -180, 180, 700, 0.5, Ebeam-2);   
  hist_W_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_W_vs_phi->GetYaxis()->SetTitle("W /GeV");

  hist_t = new TH1F("t", "t #pi^{0} all sectors", 50, 0, 2*Ebeam);   
  hist_t->GetXaxis()->SetTitle("t / GeV^{2}");
  hist_t->GetYaxis()->SetTitle("counts");

  hist_x = new TH1F("x", "x all sectors", 500, 0, 1.25);   
  hist_x->GetXaxis()->SetTitle("x");
  hist_x->GetYaxis()->SetTitle("counts");

  hist_y = new TH1F("y", "y all sectors", 440, 0, 1.1);   
  hist_y->GetXaxis()->SetTitle("y");
  hist_y->GetYaxis()->SetTitle("counts");

  hist_nu = new TH1F("nu", "#nu all sectors", 500, 0, Ebeam);   
  hist_nu->GetXaxis()->SetTitle("#nu /GeV");
  hist_nu->GetYaxis()->SetTitle("counts");

  hist_cmphi = new TH1F("cmphi", "#phi_{CM} for p of e p --> e p #pi^{0}", 180, -180, 180);   
  hist_cmphi->GetXaxis()->SetTitle("#phi_{CM} /deg");
  hist_cmphi->GetYaxis()->SetTitle("counts");
  hist_cmcostheta = new TH1F("cmcostheta_p", "cos(#Theta_{CM})", 100, -1, 1);   
  hist_cmcostheta->GetXaxis()->SetTitle("cos(#Theta_{CM})");
  hist_cmcostheta->GetYaxis()->SetTitle("counts");

   //ADDED
  hist_Q2_x = new TH2F("hist_Q2_x","Q^{2} vs Bjorken x",200, 0.0, 1.1, 200, 0.0, Ebeam);  
  h_mass_phi = new TH1F("h_mass_phi","Mass of K^{+}K^{-}",15,0.85,1.8);

  hist_x_t = new TH2F("hist_x_t","Xb vs -t", 100, 0.0, 4.0, 100, 0.0, 1.1);
  hist_x_t->GetXaxis()->SetTitle("-t [GeV]");
  hist_x_t->GetYaxis()->SetTitle("Xb");
  hist_t_phi = new TH2F("hist_t_phi","#phi vs -t", 200, -180.0, 180.0, 100, 0.0, 4.0);
  hist_t_phi->GetXaxis()->SetTitle("#phi [deg]");
  hist_t_phi->GetYaxis()->SetTitle("-t [GeV]");

  hist_Q2_phi = new TH2F("hist_q2_phi","#phi vs Q^{2}",200, -180.0, 180.0, 100, 0.0, Ebeam);
  hist_Q2_phi->GetXaxis()->SetTitle("#phi [deg]");
  hist_Q2_phi->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");  

  out->mkdir("missing_mass");				
  out->cd ("missing_mass");

  TH1F *hist_e_p_X_mismass;
  TH1F *hist_e_p_X_mismass2;

  hist_e_p_X_mismass = new TH1F("hist_e_p_X_mismass", "e p --> e p X missing mass", 200, 0.0, 2.5);   
  hist_e_p_X_mismass->GetXaxis()->SetTitle("missing mass /GeV");
  hist_e_p_X_mismass->GetYaxis()->SetTitle("counts");
  hist_e_p_X_mismass2 = new TH1F("hist_e_p_X_mismass2", "e p --> e p X missing mass squared", 800, -1, 4);   
  hist_e_p_X_mismass2->GetXaxis()->SetTitle("missing mass /GeV");
  hist_e_p_X_mismass2->GetYaxis()->SetTitle("counts");

  //ADDED
 
  TH1F *hist_eX_mass;
  TH1F *hist_epX_mass;

  out->mkdir("statistics");				
  out->cd ("statistics");

  TH1F *hist_electron_count;
  TH1F *hist_proton_count;
  TH1F *hist_neutron_count;
  TH1F *hist_pip_count;
  TH1F *hist_pim_count;
  TH1F *hist_Kp_count;
  TH1F *hist_Km_count;
  TH1F *hist_photon_count;

  hist_electron_count = new TH1F("hist_electron_count", "electron count per event", 6, -0.5, 5.5);   
  hist_electron_count->GetXaxis()->SetTitle("number of electrons per event");
  hist_electron_count->GetYaxis()->SetTitle("counts");
  hist_proton_count = new TH1F("hist_proton_count", "proton count per event", 11, -0.5, 10.5);   
  hist_proton_count->GetXaxis()->SetTitle("number of protons per event");
  hist_proton_count->GetYaxis()->SetTitle("counts");
  hist_neutron_count = new TH1F("hist_neutron_count", "neutron count per event", 11, -0.5, 10.5);   
  hist_neutron_count->GetXaxis()->SetTitle("number of neutrons per event");
  hist_neutron_count->GetYaxis()->SetTitle("counts");
  hist_pip_count = new TH1F("hist_pip_count", "pip count per event", 11, -0.5, 10.5);   
  hist_pip_count->GetXaxis()->SetTitle("number of pips per event");
  hist_pip_count->GetYaxis()->SetTitle("counts");
  hist_pim_count = new TH1F("hist_pim_count", "pim count per event", 11, -0.5, 10.5);   
  hist_pim_count->GetXaxis()->SetTitle("number of pims per event");
  hist_pim_count->GetYaxis()->SetTitle("counts");
  hist_Kp_count = new TH1F("hist_Kp_count", "Kp count per event", 11, -0.5, 10.5);   
  hist_Kp_count->GetXaxis()->SetTitle("number of Kps per event");
  hist_Kp_count->GetYaxis()->SetTitle("counts");
  hist_Km_count = new TH1F("hist_Km_count", "Km count per event", 11, -0.5, 10.5);   
  hist_Km_count->GetXaxis()->SetTitle("number of Kms per event");
  hist_Km_count->GetYaxis()->SetTitle("counts");
  hist_photon_count = new TH1F("hist_photon_count", "photon count per event", 11, -0.5, 10.5);   
  hist_photon_count->GetXaxis()->SetTitle("number of photons per event");
  hist_photon_count->GetYaxis()->SetTitle("counts");


  /// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///  start of the event loop     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  cout << "Analysing Tree: " << inTree << endl;
  cout << "Event Loop starting ... " << endl;
  cout << endl;
  
  int channel_counter = 0;
  for(Int_t k=0; k < anaTree->GetEntriesFast();k++){    

    anaTree->GetEntry(k);

    if(process_Events > 0 && k == process_Events) break;

    /// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// progress:

    if(k % 100000 == 0){

      double events = anaTree->GetEntriesFast();
      double percent = k/(events/100);
      
      printf("Analysing event number %i of %.00f (%.01f percent)\n", k, events, percent);
    }

    //pr_sect=sect_pr[select_prot_1];//
    //kp_sect=sect_kp[select_kaonP_1];//
    //km_sect=sect_kp[select_kaonM_1];

    //std::cout << " pr sector " << pr_sect << std::endl;
	//	std::cout << " kp sector " << kp_sect << std::endl;
    //	std::cout << " km sector " << km_sect << std::endl;


    /// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// assign number of particles of each type PER EVENT:
    // (SIZE OF ARRAYS)
    ele_count = p4_ele_px->size();
    prot_count = p4_prot_px->size();
    neutr_count = p4_neutr_px->size();
    pip_count = p4_pip_px->size();
    pim_count = p4_pim_px->size();
    Kp_count = p4_Kp_px->size();
    Km_count = p4_Km_px->size();
    phot_count = p4_phot_px->size();

    /// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// initalisatize variables:

    TLorentzVector beam(0,0,Ebeam,Ebeam);
    TLorentzVector target(0,0,0,0.93827);

    ///  index of selected particle for the different categories:
    // out tree variables
    select_ele = 0;  
    select_prot_1 = 0;
    evcat = 0; 
    E_ele = 0; px_ele = 0; py_ele = 0; pz_ele = 0;
    E_prot_1 = 0; px_prot_1 = 0; py_prot_1 = 0; pz_prot_1 = 0;
    //ADDED
    E_kaonP_1 = 0; px_kaonP_1 = 0; py_kaonP_1 = 0; pz_kaonP_1 = 0;
    E_kaonM_1 = 0; px_kaonM_1 = 0; py_kaonM_1 = 0; pz_kaonM_1 = 0;

    W = 0; Q2 = 0; nu = 0; x = 0; y = 0;
    t1 = 0; cmphi1 = 0; cmcostheta1 = 0; pt1 = 0; eta1 = 0; z1 = 0;

    M_e_p_X_miss = 0; M_e_p_X_miss2 = 0;

    el_sect=0; pr_sect=0; kp_sect=0; km_sect=0;
    
    /// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Assign tree components to Lorentz vectors:
    /// //////////////////////////////////////////////////////////////////

    for(Int_t i = 0; i < BUFFER; i++){ 
      if(ele_count > i){ele[i].SetPxPyPzE(p4_ele_px->at(i),p4_ele_py->at(i),p4_ele_pz->at(i),p4_ele_E->at(i));}
      else{ele[i].SetPxPyPzE(0, 0, 0, 0);}
      if(prot_count > i){prot[i].SetPxPyPzE(p4_prot_px->at(i),p4_prot_py->at(i),p4_prot_pz->at(i),p4_prot_E->at(i));}
      else{prot[i].SetPxPyPzE(0, 0, 0, 0);}
      if(neutr_count > i){neutr[i].SetPxPyPzE(p4_neutr_px->at(i),p4_neutr_py->at(i),p4_neutr_pz->at(i),p4_neutr_E->at(i));}
      else{neutr[i].SetPxPyPzE(0, 0, 0, 0);}
      if(pip_count > i){pip[i].SetPxPyPzE(p4_pip_px->at(i),p4_pip_py->at(i),p4_pip_pz->at(i),p4_pip_E->at(i));}
      else{pip[i].SetPxPyPzE(0, 0, 0, 0);}
      if(pim_count > i){pim[i].SetPxPyPzE(p4_pim_px->at(i),p4_pim_py->at(i),p4_pim_pz->at(i),p4_pim_E->at(i));}
      else{pim[i].SetPxPyPzE(0, 0, 0, 0);}
      if(Kp_count > i){Kp[i].SetPxPyPzE(p4_Kp_px->at(i),p4_Kp_py->at(i),p4_Kp_pz->at(i),p4_Kp_E->at(i));}
      else{Kp[i].SetPxPyPzE(0, 0, 0, 0);}
      if(Km_count > i){Km[i].SetPxPyPzE(p4_Km_px->at(i),p4_Km_py->at(i),p4_Km_pz->at(i),p4_Km_E->at(i));}
      else{Km[i].SetPxPyPzE(0, 0, 0, 0);}
      if(phot_count > i){phot[i].SetPxPyPzE(p4_phot_px->at(i),p4_phot_py->at(i),p4_phot_pz->at(i),p4_phot_E->at(i));}
      else{phot[i].SetPxPyPzE(0, 0, 0, 0);}      
    }


    // Assign beta from tree to array
    for( Int_t i = 0; i < BUFFER; i++ ){
      //if( prot_count > i ){ prot_beta[i]=prot_beta_final->at(i); }
      //if( Kp_count > i ){ kaonP_beta[i]=Kp_beta_final->at(i); }
      //if( Km_count > i ){ kaonM_beta[i]=Km_beta_final->at(i); }

      // and assign DC sector for hadrons from tree to array
      //if( prot_count > i ){ sect_pr[i]=prot_sector->at(i); }
      //if( Kp_count > i ){ sect_kp[i]=Kp_sector->at(i); }
      //if( Km_count > i ){ sect_km[i]=Km_sector->at(i); }

    }


    /// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Build event:
    /// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    if(ele_count == 1){    // one detected electron is the basic trigger condition for all topologies 

      select_ele = 0;   // if more than one electron is detected (< 0.2 percent of the events), take the one with the hihest momentum
          
      W  = kin_W(ele[select_ele]);
      Q2 = kin_Q2(ele[select_ele]);
      x  = kin_x(ele[select_ele]);
      y  = kin_y(ele[select_ele]);
      nu = kin_nu(ele[select_ele]);
      
      E_ele  = ele[select_ele].E();
      px_ele = ele[select_ele].Px();
      py_ele = ele[select_ele].Py();
      pz_ele = ele[select_ele].Pz();
            
      
      
      // fill kinematics:
      
      if(W > 0) hist_W->Fill(W);
      if(Q2 > 0) hist_Q2->Fill(Q2);
      //if(W > 0 && Q2 > 0) hist_Q2_vs_W->Fill(W, Q2);
      if(ele[select_ele].Phi() != 0 && W > 0) hist_W_vs_phi->Fill(ele[select_ele].Phi()*180/Pival, W);
      if(x > 0) hist_x->Fill(x);
      if(y > 0) hist_y->Fill(y);
      if(nu > 0) hist_nu->Fill(nu);
      

      if( Q2 > 0 && x > 0 )h_q2_x_cut[0]->Fill(x,Q2);
      if( Q2 > 0 && W > 0 )h_q2_w_cut[0]->Fill(W,Q2);

      if( Q2 > 1 && x > 0 )h_q2_x_cut[1]->Fill(x,Q2);
      if( Q2 > 1 && x > 0 && W > 2 )h_q2_x_cut[2]->Fill(x,Q2);
      if( Q2 > 1 && x > 0 && W > 2 )h_q2_w_cut[1]->Fill(W,Q2);

      if(W > 2 && Q2 > 1) hist_Q2_vs_W->Fill(W, Q2);
      if( Q2 > 1 && x > 0 && W > 2 ) hist_Q2_x->Fill(x,Q2);

      if( Q2 > 1 && W > 0 ) h_w_cut[0]->Fill(W);
      if( Q2 > 1 ) h_q2_cut[0]->Fill(Q2);
      if( Q2 > 1 && x > 0 ) h_x_cut[0]->Fill(x);
      if( Q2 > 1 && y > 0 ) h_y_cut[0]->Fill(y);
      if( Q2 > 1 && nu > 0 ) h_nu_cut[0]->Fill(nu);

      if( W > 2 ) h_w_cut[1]->Fill(W);
      if( W > 2 && Q2 > 0 ) h_q2_cut[1]->Fill(Q2);
      if( W > 2 && x > 0 ) h_x_cut[1]->Fill(x);
      if( W > 2 && y > 0 ) h_y_cut[1]->Fill(y);
      if( W > 2 && nu > 0 ) h_nu_cut[1]->Fill(nu);

      TLorentzVector eX_missing_mass = beam + target - ele[select_ele];
      h_eX_mass_cut[0]->Fill(eX_missing_mass.M());
      if ( Q2 > 1 && W >2 ) h_eX_mass_cut[1]->Fill( eX_missing_mass.M() );
    
      ///////////////////////////////////////////////////////////
      // category 1: e p X
      
      if(prot_count == 1){	
	evcat = 1;
	select_prot_1 = 0; // FASTEST PROTON AT INDEX 0, SAME FOR OTHER PARTICLES
      
	if(prot_count > select_prot_1){	  
	  
	  TLorentzVector X1 = miss_X(ele[select_ele], prot[select_prot_1]);
	  
	  t1 = -kin_t(ele[select_ele], X1);
	  cmphi1 = kin_cmphi(ele[select_ele], prot[select_prot_1]);
	  cmcostheta1 = kin_cmcostheta(ele[select_ele], prot[select_prot_1]); 
	  pt1 = kin_pT(ele[select_ele], prot[select_prot_1]);
	  eta1 = kin_eta(ele[select_ele], prot[select_prot_1]);
	  z1 = kin_z(ele[select_ele], prot[select_prot_1]); 
	  
	  M_e_p_X_miss = kin_mismass(ele[select_ele], prot[select_prot_1]);
	  M_e_p_X_miss2 = kin_mismass2(ele[select_ele], prot[select_prot_1]);
	  
	  E_prot_1 = prot[select_prot_1].E();
	  px_prot_1 = prot[select_prot_1].Px();
	  py_prot_1 = prot[select_prot_1].Py();
	  pz_prot_1 = prot[select_prot_1].Pz();

	  if( x > 0 && t1 > 0 ) hist_x_t->Fill(t1,x);
	  //if( cmphi1 > 0 && t1 > 0 ) hist_x_t->Fill(cmphi1,t1);
	  if( cmphi1 > 0 && Q2 > 0) hist_Q2_phi->Fill(cmphi1*180/Pival,Q2);
	  
	  
	  if( Q2 > 1 && W > 0 ) h_w_cut[2]->Fill(W);
	  if( Q2 > 1 ) h_q2_cut[2]->Fill(Q2);
	  if( Q2 > 1 && x > 0 ) h_x_cut[2]->Fill(x);
	  if( Q2 > 1 && y > 0 ) h_y_cut[2]->Fill(y);
	  if( Q2 > 1 && nu > 0 ) h_nu_cut[2]->Fill(nu);
	  if( Q2 > 1 && t1 > 0 ) h_t_cut[0]->Fill(t1);
	  if( Q2 > 1 && cmphi1 != 0 ) h_cmphi_cut[0]->Fill(cmphi1*180/Pival);
	  if( Q2 > 1 && cmcostheta1 != 0 ) h_cmcostheta_cut[0]->Fill(cmcostheta1);
	  if( Q2 > 1 && W > 0 ) h_q2_t_cut[0]->Fill(Q2,t1);

	  
	  if( W > 2 ) h_w_cut[3]->Fill(W);
	  if( W > 2 && Q2 > 0 ) h_q2_cut[3]->Fill(Q2);
	  if( W > 2 && x > 0 ) h_x_cut[3]->Fill(x);
	  if( W > 2 && y > 0 ) h_y_cut[3]->Fill(y);
	  if( W > 2 && nu > 0 ) h_nu_cut[3]->Fill(nu);
	  if( W > 2 && t1 > 0 ) h_t_cut[1]->Fill(t1);
	  if( W > 2 && cmphi1 != 0 ) h_cmphi_cut[1]->Fill(cmphi1*180/Pival);
	  if( W > 2 && cmcostheta1 != 0 ) h_cmcostheta_cut[1]->Fill(cmcostheta1);
	  if( W > 2 && t1 > 0 ) h_q2_t_cut[1]->Fill(Q2,t1);

	  
	  if(t1 > 0) hist_t->Fill(t1);
	  if(cmphi1 != 0) hist_cmphi->Fill(cmphi1*180/Pival);
	  if(cmcostheta1 != 0) hist_cmcostheta->Fill(cmcostheta1);
	  
	  h_epX_mass_cut[0]->Fill(X1.M());
	  if ( Q2 > 1 && W > 2 ) h_epX_mass_cut[1]->Fill(X1.M());
	  if( Q2 > 1 && W > 2 ) h_q2_t_cut[2]->Fill(Q2,t1);


	  if( Q2 > 0 && x > 0 )h_q2_x_cut[3]->Fill(x,Q2);
	  if( Q2 > 0 && W > 0 )h_q2_w_cut[2]->Fill(W,Q2);
	  
	  if( Q2 > 1 && x > 0 )h_q2_x_cut[4]->Fill(x,Q2);
	  if( Q2 > 1 && x > 0 && W > 2 )h_q2_x_cut[5]->Fill(x,Q2);
	  if( Q2 > 1 && W > 2 )h_q2_w_cut[3]->Fill(W,Q2);

	  if( x > 0 && t1 > 0 ) h_x_t_cut[0]->Fill(x,t1);
	  if( x > 0 && t1 > 0 && Q2 > 1 ) h_x_t_cut[1]->Fill(x,t1);
	  if( x > 0 && t1 > 0 && W > 2 ) h_x_t_cut[2]->Fill(x,t1);
	  if( x > 0 && t1 > 0 && Q2 > 1 && W > 2 ) h_x_t_cut[3]->Fill(x,t1);

	  if( cmphi1 != 0 && t1 > 0 ) h_t_phi_cut[0]->Fill(cmphi1*180/Pival,t1);
	  if( cmphi1 != 0 && t1 > 0 && Q2 > 1 ) h_t_phi_cut[1]->Fill(cmphi1*180/Pival,t1);
	  if( cmphi1 != 0 && t1 > 0 && W > 2 ) h_t_phi_cut[2]->Fill(cmphi1*180/Pival,t1);
	  if( cmphi1 != 0 && t1 > 0 && Q2 > 1 && W > 2 ) h_t_phi_cut[3]->Fill(cmphi1*180/Pival,t1);

	  if( Q2 > 0 && cmphi1 != 0 ) h_q2_phi_cut[0]->Fill(cmphi1*180/Pival, Q2);
	  if( cmphi1 != 0 && Q2 > 1 ) h_q2_phi_cut[1]->Fill(cmphi1*180/Pival, Q2);
	  if( cmphi1 != 0 && Q2 > 0 && W > 2 ) h_q2_phi_cut[2]->Fill(cmphi1*180/Pival, Q2);
	  if( cmphi1 != 0 && Q2 > 1 && W > 2 ) h_q2_phi_cut[3]->Fill(cmphi1*180/Pival, Q2);

	  if( W > 0 && cmphi1 != 0 ) h_q2_phi_cut[0]->Fill(cmphi1*180/Pival, Q2);
	  if( cmphi1 != 0 && W > 2 ) h_q2_phi_cut[1]->Fill(cmphi1*180/Pival, Q2);
	  if( cmphi1 != 0 && W > 0 && Q2 > 1 ) h_q2_phi_cut[2]->Fill(cmphi1*180/Pival, Q2);
	  if( cmphi1 != 0 && W > 2 && Q2 > 1 ) h_q2_phi_cut[3]->Fill(cmphi1*180/Pival, Q2);
	 
	}
      }
    
    
  
    
    
      // }   // end of event builder condition
      /// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
      /// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///  fill the particle historgrams:
      /// ///////////////////////////////////////////////////////////////

      // fill event properties


      for(Int_t i = 0; i < BUFFER; i++){ 

	if(ele[i].P() > 0) hist_all_electron_p->Fill(ele[i].P());
	if(ele[i].Theta() > 0) hist_all_electron_theta->Fill(ele[i].Theta()*180/Pival);
	if(ele[i].Phi() != 0) hist_all_electron_phi->Fill(ele[i].Phi()*180/Pival);
	if(ele[i].P() > 0) hist_all_electron_p_vs_theta->Fill(ele[i].Theta()*180/Pival, ele[i].P());
	if(ele[i].P() > 0) hist_all_electron_p_vs_phi->Fill(ele[i].Phi()*180/Pival, ele[i].P());
	if(ele[i].Theta() > 0 && ele[i].Phi() != 0) hist_all_electron_theta_vs_phi->Fill(ele[i].Phi()*180/Pival, ele[i].Theta()*180/Pival);

	if(prot[i].P() > 0) hist_all_proton_p->Fill(prot[i].P());
	if(prot[i].Theta() > 0) hist_all_proton_theta->Fill(prot[i].Theta()*180/Pival);
	if(prot[i].Phi() != 0) hist_all_proton_phi->Fill(prot[i].Phi()*180/Pival);
	if(prot[i].P() > 0) hist_all_proton_p_vs_theta->Fill(prot[i].Theta()*180/Pival, prot[i].P());
	if(prot[i].P() > 0) hist_all_proton_p_vs_phi->Fill(prot[i].Phi()*180/Pival, prot[i].P());
	if(prot[i].Theta() > 0 && prot[i].Phi() != 0) hist_all_proton_theta_vs_phi->Fill(prot[i].Phi()*180/Pival, prot[i].Theta()*180/Pival);

      }


      if(ele[select_ele].P() > 0) hist_electron_p->Fill(ele[select_ele].P());
      if(ele[select_ele].Theta() > 0) hist_electron_theta->Fill(ele[select_ele].Theta()*180/Pival);
      if(ele[select_ele].Phi() != 0) hist_electron_phi->Fill(ele[select_ele].Phi()*180/Pival);
      if(ele[select_ele].P() > 0) hist_electron_p_vs_theta->Fill(ele[select_ele].Theta()*180/Pival, ele[select_ele].P());
      if(ele[select_ele].P() > 0) hist_electron_p_vs_phi->Fill(ele[select_ele].Phi()*180/Pival, ele[select_ele].P());
      if(ele[select_ele].Theta() > 0 && ele[select_ele].Phi() != 0) hist_electron_theta_vs_phi->Fill(ele[select_ele].Phi()*180/Pival, ele[select_ele].Theta()*180/Pival);

      if(prot[select_prot_1].P() > 0) hist_proton_p->Fill(prot[select_prot_1].P());
      if(prot[select_prot_1].Theta() > 0) hist_proton_theta->Fill(prot[select_prot_1].Theta()*180/Pival);
      if(prot[select_prot_1].Phi() != 0) hist_proton_phi->Fill(prot[select_prot_1].Phi()*180/Pival);
      if(prot[select_prot_1].P() > 0) hist_proton_p_vs_theta->Fill(prot[select_prot_1].Theta()*180/Pival, prot[select_prot_1].P());
      if(prot[select_prot_1].P() > 0) hist_proton_p_vs_phi->Fill(prot[select_prot_1].Phi()*180/Pival, prot[select_prot_1].P());
      if(prot[select_prot_1].Theta() > 0 && prot[select_prot_1].Phi() != 0) hist_proton_theta_vs_phi->Fill(prot[select_prot_1].Phi()*180/Pival, prot[select_prot_1].Theta()*180/Pival);


      // fill missing mass and missing energy plots:

      if(M_e_p_X_miss > 0 && prot_count > 0) hist_e_p_X_mismass->Fill(M_e_p_X_miss);
      if(M_e_p_X_miss2 != 0 && prot_count > 0) hist_e_p_X_mismass2->Fill(M_e_p_X_miss2);



      /// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///  fill the statistics histograms:
      /// ///////////////////////////////////////////////////////////////

      hist_electron_count->Fill(ele_count);
      hist_proton_count->Fill(prot_count);
      hist_neutron_count->Fill(neutr_count);
      hist_pip_count->Fill(pip_count);
      hist_pim_count->Fill(pim_count);
      hist_Kp_count->Fill(Kp_count);
      hist_Km_count->Fill(Km_count);
      hist_photon_count->Fill(phot_count);
    }


    /// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
  }  // end of event loop
  /// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////

  cout << endl;
  cout << "Tree successfully analysed!" << endl;
  cout << "Writing the output file ... " << endl;
  out->Write(); // Saving Histograms
  cout << "Histograms saved in File: " << outputfile << endl;
  out->Close(); // Closing Output File
  f->Close();   // Closing Input File
  cout << "... Completed!" << endl;


  
  return 1;
    

  /// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}   /// end of main
/// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////
/// angles between particles

double alpha_e_X(TLorentzVector ele, TLorentzVector hadron){

  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector miss   = beam + target - ele - hadron;

  return ele.Vect().Angle(miss.Vect());
}

double alpha_p1p2(TLorentzVector particle1, TLorentzVector particle2){

  return particle1.Vect().Angle(particle2.Vect());
}


///////////////////////////////////////////////////////////////
///  kinematics

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

double Phi_Mass( TLorentzVector kaonP, TLorentzVector kaonM ){

  TLorentzVector lv_phi(0,0,0,0);
  lv_phi =  kaonP + kaonM;
   
  return lv_phi.M();
}


double kin_epkXMass( TLorentzVector el, TLorentzVector pr, TLorentzVector kp){

  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);

  TLorentzVector temp = beam + target - el - pr - kp;

  return temp.M();

}

double kin_epKpKmXMass( TLorentzVector el, TLorentzVector pr, TLorentzVector kp, TLorentzVector km){

  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);

  //cout << " >> " <<  beam.E() << " " << beam.Pz() << " " << el.E() << " " << pr.E() << " " << kp.E() << " " << km.E() << endl;

  TLorentzVector temp = beam + target - el - pr - kp - km; 

  return temp.M2();

}

