/// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///  ROOT Macro to filter clas 12 data
///
///  Stefan Diehl  (sdiehl@jlab.org)
///
///  To execute run ROOT and compile macro with .L filter_clas12.cxx++
///
///  Then call function:  filter_clas12("trains_v2/skim4_4013.root", "output/skim4_4013.root", 4013)
///                       filter_clas12("trains_v2/skim4_4014.root", "output/skim4_4014.root", 4014)
///
///                       filter_clas12("trains_v2/skim4_5030.root", "output/inbending_test.root", 5030)
///                       filter_clas12("trains_v2/skim4_5532.root", "output/outbending_test.root", 5532)
///
///                       filter_clas12("SIDIS_sim/MC_DIS_F18_0.root", "output/MC_DIS_F18_sim_0.root", 11)
///                       filter_clas12("SIDIS_sim/clas12_pip_sim_inbending.root", "output/recon.root", 11)
//
///  FD eid cuts:  0 = negative charge,  1 = PID default, 2 = PID + charge + EC inner vs outer,  3 = PID + charge + EC sampling fraction,  4 = PID + charge + PCAL hit position fiducial,
///                5 = PID + charge + DC fiducial region 1,  6 = PID + charge + DC fiducial region 3,  7 = PID + charge + DC z vertex, 
///                8 = all cuts passed, except the one based on the shown histogram, 9 = all cuts,  10 = inverse PID cut
///
///  FD hadron ID cuts:  0 = PID default,  1 = charge,  2 = DC fiducial region 1, 3 = DC fiducial region 3
///                      4 = beta vs p,  5 = delta beta,  6 = tofmass, 7 = maximum probability,  8 = delta vz, 9 = all  
///                      --> cuts 3 to 8 contain charge and DC fiducial cut as basis, for neutrons, the DC fiducial cut is not applied
///                      --> For negative pions and Kaons also a a cut, requireing, that the particle is not a correctly identified electron is applied to cut 1 - 7
///                      --> For negative pion also a EC inner vs outer cut is applied to reject electrons to cut 1 - 7
///
///                      Proton = + 0,  Neutron = + 10,  Pip = + 20,  Pim = + 30,  Kp = + 40,  Km = + 50
///
///  CD hadron ID:  0 = PID default,   1 = charge,   2 = PID + charge + beta vs p,   3 = PID + charge + maximum probability,   4 = PID + charge + z vertex,  5 = all cuts
///
///
///  Photon ID:  0 = PID default,  1 = charge,  2 = beta cut,  3 = PCAL hit position fiducial cut,  4 = all cuts
///
/// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
#include "TH1F.h"
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
#include <TChain.h>
#include <TCutG.h>
#include <fstream>
using namespace std;

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

using namespace ROOT::Math;


/// //////////////////////////////////////////////////////////////////////////////////////////////////////
/// settings:

int process_Events = -1;
//int process_Events = 2000000;

//float Ebeam = 6.42313;
//float Ebeam = 2.22193;
//float Ebeam = 10.594;
float Ebeam = 2.211;//10.6;

/// select the field polarity (relevant for fiducial cuts, sampling fraction and vertex cuts)
/// only set 1 as true!

bool inbending  = true;
bool outbending = false;

///

bool simulation  = true;

///

bool use_own_PID_electron = false;   
bool use_own_PID_proton =   false;
bool use_own_PID_neutron =  false;
bool use_own_PID_pip =      false;
bool use_own_PID_pim =      false;
bool use_own_PID_Kp =       false;
bool use_own_PID_Km =       false;
bool use_own_PID_photon =   false;

bool use_own_PID_CD_proton =  false;
bool use_own_PID_CD_neutron = false;
bool use_own_PID_CD_pip =     false;
bool use_own_PID_CD_pim =     false;
bool use_own_PID_CD_Kp =      false;
bool use_own_PID_CD_Km =      false;

bool use_FT = false;
bool use_FD = true;
bool use_CD = true;

// only cut on beta vs p for hadrons, do not apply delta beta and tofmass cut in addition?

bool cut_maximum_probability_prot = true; 
bool cut_maximum_probability_pip  = true;  
bool cut_maximum_probability_pim  = true;  
bool cut_maximum_probability_Kp   = true;  
bool cut_maximum_probability_Km   = true;  
bool population_weighting = false;

bool cut_beta_vs_p_prot = false;
bool cut_beta_vs_p_pip  = false;
bool cut_beta_vs_p_pim  = false;
bool cut_beta_vs_p_Kp   = false; 
bool cut_beta_vs_p_Km   = false;  
bool cut_deltabeta_prot = false;
bool cut_deltabeta_pip  = false;
bool cut_deltabeta_pim  = false;
bool cut_deltabeta_Kp   = false;
bool cut_deltabeta_Km   = false;
bool cut_tofmass_prot = false;
bool cut_tofmass_pip  = false;
bool cut_tofmass_pim  = false;
bool cut_tofmass_Kp   = false;
bool cut_tofmass_Km   = false;

bool cut_beta_vs_p_only_neutron = true;

// CD:

bool CD_cut_maximum_probability_prot = true; 
bool CD_cut_maximum_probability_pip  = true; 
bool CD_cut_maximum_probability_pim  = true; 
bool CD_cut_maximum_probability_Kp   = true; 
bool CD_cut_maximum_probability_Km   = true;
bool population_weighting_CD = false;

bool CD_cut_beta_vs_p_prot = false;
bool CD_cut_beta_vs_p_pip = false;
bool CD_cut_beta_vs_p_pim = false;
bool CD_cut_beta_vs_p_Kp = false;
bool CD_cut_beta_vs_p_Km = false;

// apply time of flight correction

bool correct_FTOF = false;
bool correct_CTOF = false;

// use default, not corrected beta?

bool use_own_beta_charged = false;
bool use_own_beta_neutrals = false;
bool use_own_beta_charged_FT = false;
bool use_own_beta_neutrals_FT = false;

bool apply_correction_electron = false;
bool apply_correction_proton = false;
bool apply_correction_neutron = false;
bool apply_correction_pip = false;
bool apply_correction_pim = false;
bool apply_correction_Kp = false;
bool apply_correction_Km = false;
bool apply_correction_photon = false;

bool show_FD_eid_statistics = true;
bool show_charged_ID_statistics = true;
bool show_photon_ID_statistics = true;

bool fill_electron_pid_histograms = true;
bool fill_hadron_pid_histograms = true;
bool fill_photon_pid_histograms = true;


const static int BUFFER = 60;   // increased from 40


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


/// //////////////////////////////////////////////////////////////////////////////////////////////////////
/// List of CLAS12 hipo branches

   // REC::Event

   vector<int>            *vNRUN     = 0;     // Run Number
   vector<int>            *vNEVENT   = 0;     // Event Number
   vector<float>          *vEVNTime  = 0;     // Enet time
   vector<short>          *vTYPE     = 0;     // Data or MC
   vector<long long int>  *vTRG      = 0;     // Trigger
   vector<float>          *vBCG      = 0;     // FCUP scaler
   vector<float>          *vSTTime   = 0;     // Event Start Time (ns)
   vector<float>          *vRFTime   = 0;     // RF Time (ns)
   vector<short>          *vHelic    = 0;     // helicity

   // REC::Particle

   vector<int>    *vpart_pid     = 0;   
   vector<short>  *vpart_charge  = 0; 
   vector<float>  *vpart_px      = 0;  
   vector<float>  *vpart_py      = 0;  
   vector<float>  *vpart_pz      = 0; 
   vector<float>  *vpart_vx      = 0;   
   vector<float>  *vpart_vy      = 0;  
   vector<float>  *vpart_vz      = 0;
   vector<float>  *vpart_beta    = 0;
   vector<short>  *vpart_status  = 0;     // indicated FT (1000 <= status < 2000), FD (2000 <= status < 3000) or CD (status >= 40000)  
 
   // REC::Calorimeter

   vector<int>    *vCal_pindex   = 0;   
   vector<short>  *vCal_detector = 0;   
   vector<short>  *vCal_sector   = 0; 
   vector<short>  *vCal_layer    = 0;   
   vector<float>  *vCal_energy   = 0;   
   vector<float>  *vCal_time     = 0; 
   vector<float>  *vCal_path     = 0;   
   vector<float>  *vCal_x        = 0; 
   vector<float>  *vCal_y 	 = 0;   
   vector<float>  *vCal_z 	 = 0; 
   vector<float>  *vCal_lu 	 = 0; 
   vector<float>  *vCal_lv 	 = 0;   
   vector<float>  *vCal_lw 	 = 0; 


   // REC::Cherenkov

   vector<int>    *vCC_pindex   = 0;   
   vector<short>  *vCC_detector = 0;   
   vector<short>  *vCC_sector   = 0; 
   vector<float>  *vCC_nphe     = 0; 
   vector<float>  *vCC_time     = 0; 
   vector<float>  *vCC_path     = 0;   
   vector<float>  *vCC_theta    = 0;   
   vector<float>  *vCC_phi      = 0; 

   // REC::ForwardTagger

   vector<int>    *vFT_pindex   = 0;   
   vector<short>  *vFT_detector = 0;   
   vector<float>  *vFT_energy   = 0;   
   vector<float>  *vFT_time     = 0; 
   vector<float>  *vFT_path     = 0;   
   vector<float>  *vFT_x        = 0; 
   vector<float>  *vFT_y        = 0;   
   vector<float>  *vFT_z        = 0; 
   vector<float>  *vFT_dx       = 0; 
   vector<float>  *vFT_dy       = 0;  
   vector<float>  *vFT_radius   = 0;    // cluster radius
   vector<float>  *vFT_size     = 0;    // cluster size 

   // REC::Scintillator

   vector<int>    *vSC_pindex    = 0;   
   vector<short>  *vSC_detector  = 0;  
   vector<short>  *vSC_sector    = 0;   
   vector<short>  *vSC_layer     = 0;   
   vector<int>    *vSC_component = 0;   
   vector<float>  *vSC_energy    = 0;   
   vector<float>  *vSC_time      = 0; 
   vector<float>  *vSC_path      = 0;   

   // REC::Track
 
   vector<int>    *vTRK_pindex     = 0;   
   vector<short>  *vTRK_detector   = 0;  
   vector<short>  *vTRK_sector     = 0;  
   vector<float>  *vTRK_chi2       = 0;  
   vector<int>    *vTRK_NDF        = 0;  
   vector<float>  *vTRK_chi2_nomm  = 0;  
   vector<int>    *vTRK_NDF_nomm   = 0;  

   // REC::Traj

   vector<int>    *vTraj_pindex = 0;   
   vector<int>    *vTraj_detID  = 0;   
   vector<float>  *vTraj_x      = 0;    
   vector<float>  *vTraj_y      = 0;    
   vector<float>  *vTraj_z      = 0;   
   vector<float>  *vTraj_cx     = 0;    
   vector<float>  *vTraj_cy     = 0;    
   vector<float>  *vTraj_cz     = 0;    

   // MC::Particle  and MC::Lund

   vector<int>    *vMC_Header_helicity  = 0; 
   vector<int>    *vMC_Event_npart  = 0;     
   vector<float>  *vMC_Event_ebeam  = 0;    
   vector<float>  *vMC_Event_weight = 0; 

   vector<int>    *vMC_Particle_pid = 0;     
   vector<float>  *vMC_Particle_px  = 0;    
   vector<float>  *vMC_Particle_py  = 0;     
   vector<float>  *vMC_Particle_pz  = 0;     
   vector<float>  *vMC_Particle_vx  = 0;    
   vector<float>  *vMC_Particle_vy  = 0;
   vector<float>  *vMC_Particle_vz  = 0;

   vector<int>    *vMC_Lund_pid  = 0; 
   vector<float>  *vMC_Lund_mass = 0;  
   vector<float>  *vMC_Lund_E    = 0;      
   vector<float>  *vMC_Lund_px   = 0;    
   vector<float>  *vMC_Lund_py   = 0;     
   vector<float>  *vMC_Lund_pz   = 0;     
   vector<float>  *vMC_Lund_vx   = 0;    
   vector<float>  *vMC_Lund_vy   = 0;
   vector<float>  *vMC_Lund_vz   = 0;


/// /////////////////////////////////////////////////////////////////////////////////////
///  output file:

  TFile *out;

  char name[200];
  char title[200];


/// //////////////////////////////////////////////////////////////////////////////////////////////////////
/// general particle properties and additional particle bank properties for non identified particles:

  short TYPE, Helic;
  int NRUN, NEVENT;
  long long int TRG;
  float EVNTime, BCG, STTime, RFTime;

/// //////////////////////////////////////////////////////////////////////////////////////////////////////
/// general particle properties and additional particle bank properties for non identified particles:

  int Npart;

  int eventNumber[BUFFER];
  float part_px[BUFFER], part_py[BUFFER], part_pz[BUFFER], part_beta[BUFFER];
  float part_p[BUFFER], part_theta[BUFFER], part_phi[BUFFER];
  float part_vx[BUFFER], part_vy[BUFFER], part_vz[BUFFER]; 
  short part_charge[BUFFER], part_status[BUFFER];
  int   part_pid[BUFFER];

  int ele_detect[BUFFER], prot_detect[BUFFER], neutr_detect[BUFFER], pip_detect[BUFFER], pim_detect[BUFFER], Kp_detect[BUFFER], Km_detect[BUFFER], phot_detect[BUFFER]; 

  /// ///////////////////////////////
  /// MC:

  float MC_Ebeam, MC_weight; 
  int MC_helicity, MC_Npart;

  int partMC_pid[BUFFER];
  float partMC_px[BUFFER], partMC_py[BUFFER], partMC_pz[BUFFER];
  float partMC_p[BUFFER], partMC_theta[BUFFER], partMC_phi[BUFFER];
  float partMC_vx[BUFFER], partMC_vy[BUFFER], partMC_vz[BUFFER]; 

  int partLUND_pid[BUFFER];
  float partLUND_mass[BUFFER], partLUND_E[BUFFER];
  float partLUND_px[BUFFER], partLUND_py[BUFFER], partLUND_pz[BUFFER];
  float partLUND_p[BUFFER], partLUND_theta[BUFFER], partLUND_phi[BUFFER];
  float partLUND_vx[BUFFER], partLUND_vy[BUFFER], partLUND_vz[BUFFER]; 

  /// //////////////////////////////

  int   part_Cal_PCAL_sector[BUFFER], part_Cal_ECin_sector[BUFFER], part_Cal_ECout_sector[BUFFER];  
  float part_Cal_PCAL_energy[BUFFER], part_Cal_ECin_energy[BUFFER], part_Cal_ECout_energy[BUFFER], part_Cal_energy_total[BUFFER]; 
  float part_Cal_PCAL_time[BUFFER], part_Cal_ECin_time[BUFFER], part_Cal_ECout_time[BUFFER];   
  float part_Cal_PCAL_path[BUFFER], part_Cal_ECin_path[BUFFER], part_Cal_ECout_path[BUFFER];   
  float part_Cal_PCAL_x[BUFFER], part_Cal_PCAL_y[BUFFER], part_Cal_PCAL_z[BUFFER]; 
  float part_Cal_ECin_x[BUFFER], part_Cal_ECin_y[BUFFER], part_Cal_ECin_z[BUFFER]; 
  float part_Cal_ECout_x[BUFFER], part_Cal_ECout_y[BUFFER], part_Cal_ECout_z[BUFFER]; 
  float part_Cal_PCAL_lu[BUFFER], part_Cal_PCAL_lv[BUFFER], part_Cal_PCAL_lw[BUFFER]; 
  float part_Cal_ECin_lu[BUFFER], part_Cal_ECin_lv[BUFFER], part_Cal_ECin_lw[BUFFER]; 
  float part_Cal_ECout_lu[BUFFER], part_Cal_ECout_lv[BUFFER], part_Cal_ECout_lw[BUFFER]; 

  int part_CC_HTCC_sector[BUFFER], part_CC_HTCC_nphe[BUFFER];   
  float part_CC_HTCC_time[BUFFER], part_CC_HTCC_path[BUFFER];   
  float part_CC_HTCC_theta[BUFFER], part_CC_HTCC_phi[BUFFER]; 

  int part_CC_LTCC_sector[BUFFER], part_CC_LTCC_nphe[BUFFER];   
  float part_CC_LTCC_time[BUFFER], part_CC_LTCC_path[BUFFER];   
  float part_CC_LTCC_theta[BUFFER], part_CC_LTCC_phi[BUFFER]; 

  int   part_FTOF_layer[BUFFER];
  int   part_FTOF_sector_layer1[BUFFER], part_FTOF_sector_layer2[BUFFER], part_FTOF_sector_layer3[BUFFER];
  int   part_FTOF_component_layer1[BUFFER], part_FTOF_component_layer2[BUFFER], part_FTOF_component_layer3[BUFFER];   
  float part_FTOF_energy[BUFFER], part_FTOF_time[BUFFER], part_FTOF_path[BUFFER];  
  float part_FTOF_energy_layer1[BUFFER], part_FTOF_time_layer1[BUFFER], part_FTOF_path_layer1[BUFFER];  
  float part_FTOF_energy_layer3[BUFFER], part_FTOF_time_layer3[BUFFER], part_FTOF_path_layer3[BUFFER];  

  int   part_CTOF_component[BUFFER]; 
  float part_CTOF_energy[BUFFER], part_CTOF_time[BUFFER], part_CTOF_path[BUFFER];   

  int   part_CND_component[BUFFER]; 
  float part_CND_energy[BUFFER], part_CND_time[BUFFER], part_CND_path[BUFFER];   
  
  float part_FT_energy[BUFFER], part_FT_radius[BUFFER], part_FTHODO_time[BUFFER], part_FTHODO_path[BUFFER] , part_FTCAL_time[BUFFER], part_FTCAL_path[BUFFER];  
  float part_FTCAL_x[BUFFER], part_FTCAL_y[BUFFER], part_FTCAL_z[BUFFER];
  float part_FTTRK_x[BUFFER], part_FTTRK_y[BUFFER], part_FTTRK_z[BUFFER];
  float part_FTHODO_x[BUFFER], part_FTHODO_y[BUFFER], part_FTHODO_z[BUFFER];

  int part_DC_sector[BUFFER], part_DC_index[BUFFER];
  float part_DC_c1x[BUFFER], part_DC_c1y[BUFFER], part_DC_c1z[BUFFER];
  float part_DC_c2x[BUFFER], part_DC_c2y[BUFFER], part_DC_c2z[BUFFER];
  float part_DC_c3x[BUFFER], part_DC_c3y[BUFFER], part_DC_c3z[BUFFER];


/// //////////////////////////////////////////////////////////////////////////////////////////////////////
/// partcile variables for identified particles:

  int e_count, p_count, n_count, pip_count, pim_count, Kp_count, Km_count, g_count;
  int e_MCcount, p_MCcount, n_MCcount, pip_MCcount, pim_MCcount, Kp_MCcount, Km_MCcount, g_MCcount;  
  int e_ind[BUFFER], p_ind[BUFFER], n_ind[BUFFER], pip_ind[BUFFER], pim_ind[BUFFER], Kp_ind[BUFFER], Km_ind[BUFFER], g_ind[BUFFER];

  TLorentzVector p4_ele[BUFFER];
  TLorentzVector p4_ele_raw[BUFFER];
  TLorentzVector p4_prot[BUFFER];
  TLorentzVector p4_neutr[BUFFER];
  TLorentzVector p4_pip[BUFFER];
  TLorentzVector p4_pim[BUFFER];
  TLorentzVector p4_Kp[BUFFER];
  TLorentzVector p4_Km[BUFFER];
  TLorentzVector p4_phot[BUFFER];

  float e_vx[BUFFER], e_vy[BUFFER], e_vz[BUFFER], e_beta[BUFFER], e_FTOF_sec[BUFFER], e_PCAL_sec[BUFFER];
  float p_vx[BUFFER], p_vy[BUFFER], p_vz[BUFFER], p_beta[BUFFER], p_FTOF_sec[BUFFER], p_PCAL_sec[BUFFER];
  float n_vx[BUFFER], n_vy[BUFFER], n_vz[BUFFER], n_beta[BUFFER];
  float pip_vx[BUFFER], pip_vy[BUFFER], pip_vz[BUFFER], pip_beta[BUFFER], pip_FTOF_sec[BUFFER];
  float pim_vx[BUFFER], pim_vy[BUFFER], pim_vz[BUFFER], pim_beta[BUFFER], pim_FTOF_sec[BUFFER];
  float Kp_vx[BUFFER], Kp_vy[BUFFER], Kp_vz[BUFFER], Kp_beta[BUFFER], Kp_FTOF_sec[BUFFER];
  float Km_vx[BUFFER], Km_vy[BUFFER], Km_vz[BUFFER], Km_beta[BUFFER], Km_FTOF_sec[BUFFER];
  float g_vx[BUFFER], g_vy[BUFFER], g_vz[BUFFER], g_sec[BUFFER];


// variables for reconstructed particles

TLorentzVector neutral_iter[28];  // possible pair combinations for 8 photons
TLorentzVector pi0;
TLorentzVector eta;

int neutral_ind[28][2];
int pi0_ind[2];
int eta_ind[2];

int pi0_g_ind[28];

double mass_neutral_iter[28];
double mass_neutral_iter2[28];
double alpha_gg[28];

int select_pi0;
int select_eta;

// cut variables

  bool electron_sector_cut[BUFFER];


// kinematic variables

  float W, Q2, nu, x, y, t_pi0, t_pip, t_pim, cmphi_p, cmcostheta_p, cmphi_pi0, cmcostheta_pi0, cmphi_pip, cmcostheta_pip, cmphi_pim, cmcostheta_pim;


// vectors with components of the Lorentzvector to fill into the output tree 
// --> vectors are used to optimize the file size

  int helicity;
  double fcup;
  vector<int> sectorE;
  vector<int> electron_event_number;
  vector<double> p4_ele_px;
  vector<double> p4_ele_py;
  vector<double> p4_ele_pz;
  vector<double> p4_ele_E;
  vector<double> p4_prot_px;
  vector<double> p4_prot_py;
  vector<double> p4_prot_pz;
  vector<double> p4_prot_E;
  vector<double> p4_neutr_px; 
  vector<double> p4_neutr_py; 
  vector<double> p4_neutr_pz; 
  vector<double> p4_neutr_E; 
  vector<double> p4_pip_px; 
  vector<double> p4_pip_py; 
  vector<double> p4_pip_pz; 
  vector<double> p4_pip_E; 
  vector<double> p4_pim_px; 
  vector<double> p4_pim_py;
  vector<double> p4_pim_pz;
  vector<double> p4_pim_E;
  vector<double> p4_Kp_px; 
  vector<double> p4_Kp_py;
  vector<double> p4_Kp_pz;
  vector<double> p4_Kp_E;
  vector<double> p4_Km_px;
  vector<double> p4_Km_py;
  vector<double> p4_Km_pz;
  vector<double> p4_Km_E;
  vector<double> p4_phot_px;
  vector<double> p4_phot_py;
  vector<double> p4_phot_pz;
  vector<double> p4_phot_E;
  vector<int> ele_det;
  vector<int> prot_det;
  vector<int> neutr_det;
  vector<int> pip_det;
  vector<int> pim_det;
  vector<int> Kp_det;
  vector<int> Km_det;
  vector<int> phot_det;


// cut statistics:
int event;
long neg_part_count;
long pos_part_count;
long neut_part_count;

long Track_Quality_pass;

long FD_eid_default_PID_pass;
long FD_eid_charge_pass;
long FD_eid_CC_nphe_pass;
long FD_eid_EC_outer_vs_EC_inner_pass;
long FD_eid_EC_sampling_fraction_pass;
long FD_eid_EC_hit_position_fiducial_pass;
long FD_eid_DC_hit_position_region1_fiducial_pass;
long FD_eid_DC_hit_position_region2_fiducial_pass;
long FD_eid_DC_hit_position_region3_fiducial_pass;
long FD_eid_DC_z_vertex_pass;
long FD_eid_all_pass;

long FD_protid_default_PID_pass;
long FD_protid_charge_pass;
long FD_protid_DC_hit_position_region1_fiducial_pass;
long FD_protid_DC_hit_position_region2_fiducial_pass;
long FD_protid_DC_hit_position_region3_fiducial_pass;
long FD_protid_beta_pass;
long FD_protid_delta_beta_pass;
long FD_protid_tofmass_pass;
long FD_protid_maximum_probability_pass;
long FD_protid_delta_vz_pass;
long FD_protid_all_pass;

long FD_neutrid_default_PID_pass;
long FD_neutrid_charge_pass;
long FD_neutrid_beta_pass;
long FD_neutrid_delta_beta_pass;
long FD_neutrid_tofmass_pass;
long FD_neutrid_delta_vz_pass;
long FD_neutrid_all_pass;

long FD_pipid_default_PID_pass;
long FD_pipid_charge_pass;
long FD_pipid_DC_hit_position_region1_fiducial_pass;
long FD_pipid_DC_hit_position_region2_fiducial_pass;
long FD_pipid_DC_hit_position_region3_fiducial_pass;
long FD_pipid_beta_pass;
long FD_pipid_delta_beta_pass;
long FD_pipid_tofmass_pass;
long FD_pipid_maximum_probability_pass;
long FD_pipid_delta_vz_pass;
long FD_pipid_all_pass;

long FD_pimid_default_PID_pass;
long FD_pimid_charge_pass;
long FD_pimid_DC_hit_position_region1_fiducial_pass;
long FD_pimid_DC_hit_position_region2_fiducial_pass;
long FD_pimid_DC_hit_position_region3_fiducial_pass;
long FD_pimid_beta_pass;
long FD_pimid_delta_beta_pass;
long FD_pimid_tofmass_pass;
long FD_pimid_maximum_probability_pass;
long FD_pimid_delta_vz_pass;
long FD_pimid_all_pass;

long FD_Kpid_default_PID_pass;
long FD_Kpid_charge_pass;
long FD_Kpid_DC_hit_position_region1_fiducial_pass;
long FD_Kpid_DC_hit_position_region2_fiducial_pass;
long FD_Kpid_DC_hit_position_region3_fiducial_pass;
long FD_Kpid_beta_pass;
long FD_Kpid_delta_beta_pass;
long FD_Kpid_tofmass_pass;
long FD_Kpid_maximum_probability_pass;
long FD_Kpid_delta_vz_pass;
long FD_Kpid_all_pass;

long FD_Kmid_default_PID_pass;
long FD_Kmid_charge_pass;
long FD_Kmid_DC_hit_position_region1_fiducial_pass;
long FD_Kmid_DC_hit_position_region2_fiducial_pass;
long FD_Kmid_DC_hit_position_region3_fiducial_pass;
long FD_Kmid_beta_pass;
long FD_Kmid_delta_beta_pass;
long FD_Kmid_tofmass_pass;
long FD_Kmid_maximum_probability_pass;
long FD_Kmid_delta_vz_pass;
long FD_Kmid_all_pass;

long FD_photid_default_PID_pass;
long FD_photid_charge_pass;
long FD_photid_beta_pass;
long FD_photid_EC_sampling_fraction_pass;
long FD_photid_EC_hit_position_fiducial_pass;
long FD_photid_all_pass;

/// FT

long FT_eid_charge_pass;
long FT_eid_PID_pass;
long FT_eid_FTCAL_fiducial_pass;
long FT_eid_FTTRK_fiducial_pass;
long FT_eid_FTHODO_fiducial_pass;
long FT_eid_energy_vs_radius_pass;
long FT_eid_all_pass;

long FT_photid_charge_pass;
long FT_photid_PID_pass;
long FT_photid_FTCAL_fiducial_pass;
long FT_photid_beta_pass;
long FT_photid_all_pass;

/// CD

long CD_protid_default_PID_pass;
long CD_protid_charge_pass;
long CD_protid_beta_pass;
long CD_protid_maximum_probability_pass;
long CD_protid_delta_vz_pass;
long CD_protid_all_pass;

long CD_neutrid_default_PID_pass;
long CD_neutrid_charge_pass;
long CD_neutrid_beta_pass;
long CD_neutrid_maximum_probability_pass;
long CD_neutrid_delta_vz_pass;
long CD_neutrid_all_pass;

long CD_pipid_default_PID_pass;
long CD_pipid_charge_pass;
long CD_pipid_beta_pass;
long CD_pipid_maximum_probability_pass;
long CD_pipid_delta_vz_pass;
long CD_pipid_all_pass;

long CD_pimid_default_PID_pass;
long CD_pimid_charge_pass;
long CD_pimid_beta_pass;
long CD_pimid_maximum_probability_pass;
long CD_pimid_delta_vz_pass;
long CD_pimid_all_pass;

long CD_Kpid_default_PID_pass;
long CD_Kpid_charge_pass;
long CD_Kpid_beta_pass;
long CD_Kpid_maximum_probability_pass;
long CD_Kpid_delta_vz_pass;
long CD_Kpid_all_pass;

long CD_Kmid_default_PID_pass;
long CD_Kmid_charge_pass;
long CD_Kmid_beta_pass;
long CD_Kmid_maximum_probability_pass;
long CD_Kmid_delta_vz_pass;
long CD_Kmid_all_pass;


// cut selectors variables

bool Track_Quality_check[BUFFER];

/// /////////////////////////////////////////////////////////////
/// Forward detector

bool FD_eid_default_PID_check[BUFFER];
bool FD_eid_charge_check[BUFFER];
bool FD_eid_CC_nphe_check[BUFFER];
bool FD_eid_EC_outer_vs_EC_inner_check[BUFFER];
bool FD_eid_EC_sampling_fraction_check[BUFFER];
bool FD_eid_EC_hit_position_fiducial_check[BUFFER];
bool FD_eid_DC_hit_position_region1_fiducial_check[BUFFER];
bool FD_eid_DC_hit_position_region2_fiducial_check[BUFFER];
bool FD_eid_DC_hit_position_region3_fiducial_check[BUFFER];
bool FD_eid_DC_z_vertex_check[BUFFER];
bool FD_eid_all_check[BUFFER];

bool FD_protid_default_PID_check[BUFFER];
bool FD_protid_charge_check[BUFFER];
bool FD_protid_DC_hit_position_region1_fiducial_check[BUFFER];
bool FD_protid_DC_hit_position_region2_fiducial_check[BUFFER];
bool FD_protid_DC_hit_position_region3_fiducial_check[BUFFER];
bool FD_protid_beta_check[BUFFER];
bool FD_protid_delta_beta_check[BUFFER];
bool FD_protid_tofmass_check[BUFFER];
bool FD_protid_maximum_probability_check[BUFFER];
bool FD_protid_delta_vz_check[BUFFER];
bool FD_protid_all_check[BUFFER];

bool FD_neutrid_default_PID_check[BUFFER];
bool FD_neutrid_charge_check[BUFFER];
bool FD_neutrid_beta_check[BUFFER];
bool FD_neutrid_delta_beta_check[BUFFER];
bool FD_neutrid_tofmass_check[BUFFER];
bool FD_neutrid_delta_vz_check[BUFFER];
bool FD_neutrid_all_check[BUFFER];

bool FD_pipid_default_PID_check[BUFFER];
bool FD_pipid_charge_check[BUFFER];
bool FD_pipid_DC_hit_position_region1_fiducial_check[BUFFER];
bool FD_pipid_DC_hit_position_region2_fiducial_check[BUFFER];
bool FD_pipid_DC_hit_position_region3_fiducial_check[BUFFER];
bool FD_pipid_beta_check[BUFFER];
bool FD_pipid_delta_beta_check[BUFFER];
bool FD_pipid_tofmass_check[BUFFER];
bool FD_pipid_maximum_probability_check[BUFFER];
bool FD_pipid_delta_vz_check[BUFFER];
bool FD_pipid_all_check[BUFFER];

bool FD_pimid_default_PID_check[BUFFER];
bool FD_pimid_charge_check[BUFFER];
bool FD_pimid_ele_reject_check[BUFFER];
bool FD_pimid_EC_outer_vs_EC_inner_check[BUFFER];
bool FD_pimid_DC_hit_position_region1_fiducial_check[BUFFER];
bool FD_pimid_DC_hit_position_region2_fiducial_check[BUFFER];
bool FD_pimid_DC_hit_position_region3_fiducial_check[BUFFER];
bool FD_pimid_beta_check[BUFFER];
bool FD_pimid_delta_beta_check[BUFFER];
bool FD_pimid_tofmass_check[BUFFER];
bool FD_pimid_maximum_probability_check[BUFFER];
bool FD_pimid_delta_vz_check[BUFFER];
bool FD_pimid_all_check[BUFFER];

bool FD_Kpid_default_PID_check[BUFFER];
bool FD_Kpid_charge_check[BUFFER];
bool FD_Kpid_DC_hit_position_region1_fiducial_check[BUFFER];
bool FD_Kpid_DC_hit_position_region2_fiducial_check[BUFFER];
bool FD_Kpid_DC_hit_position_region3_fiducial_check[BUFFER];
bool FD_Kpid_beta_check[BUFFER];
bool FD_Kpid_delta_beta_check[BUFFER];
bool FD_Kpid_tofmass_check[BUFFER];
bool FD_Kpid_maximum_probability_check[BUFFER];
bool FD_Kpid_delta_vz_check[BUFFER];
bool FD_Kpid_all_check[BUFFER];

bool FD_Kmid_default_PID_check[BUFFER];
bool FD_Kmid_charge_check[BUFFER];
bool FD_Kmid_ele_reject_check[BUFFER];
bool FD_Kmid_EC_outer_vs_EC_inner_check[BUFFER];
bool FD_Kmid_DC_hit_position_region1_fiducial_check[BUFFER];
bool FD_Kmid_DC_hit_position_region2_fiducial_check[BUFFER];
bool FD_Kmid_DC_hit_position_region3_fiducial_check[BUFFER];
bool FD_Kmid_beta_check[BUFFER];
bool FD_Kmid_delta_beta_check[BUFFER];
bool FD_Kmid_tofmass_check[BUFFER];
bool FD_Kmid_maximum_probability_check[BUFFER];
bool FD_Kmid_delta_vz_check[BUFFER];
bool FD_Kmid_all_check[BUFFER];

bool FD_photid_default_PID_check[BUFFER];
bool FD_photid_charge_check[BUFFER];
bool FD_photid_beta_check[BUFFER];
bool FD_photid_EC_sampling_fraction_check[BUFFER];
bool FD_photid_EC_outer_vs_EC_inner_check[BUFFER];
bool FD_photid_EC_hit_position_fiducial_check[BUFFER];
bool FD_photid_all_check[BUFFER];


/// //////////////////////////////////////////////////////////
/// FT

bool FT_eid_charge_check[BUFFER];
bool FT_eid_PID_check[BUFFER];
bool FT_eid_FTCAL_fiducial_check[BUFFER];
bool FT_eid_FTTRK_fiducial_check[BUFFER];
bool FT_eid_FTHODO_fiducial_check[BUFFER];
bool FT_eid_energy_vs_radius_check[BUFFER];
bool FT_eid_all_check[BUFFER];

bool FT_photid_charge_check[BUFFER];
bool FT_photid_PID_check[BUFFER];
bool FT_photid_FTCAL_fiducial_check[BUFFER];
bool FT_photid_beta_check[BUFFER];
bool FT_photid_all_check[BUFFER];

/// //////////////////////////////////////////////////////////////
/// Central detector

bool CD_protid_default_PID_check[BUFFER];
bool CD_protid_charge_check[BUFFER];
bool CD_protid_beta_check[BUFFER];
bool CD_protid_maximum_probability_check[BUFFER];
bool CD_protid_delta_vz_check[BUFFER];
bool CD_protid_all_check[BUFFER];

bool CD_neutrid_default_PID_check[BUFFER];
bool CD_neutrid_charge_check[BUFFER];
bool CD_neutrid_beta_check[BUFFER];
bool CD_neutrid_maximum_probability_check[BUFFER];
bool CD_neutrid_delta_vz_check[BUFFER];
bool CD_neutrid_all_check[BUFFER];

bool CD_pipid_default_PID_check[BUFFER];
bool CD_pipid_charge_check[BUFFER];
bool CD_pipid_beta_check[BUFFER];
bool CD_pipid_maximum_probability_check[BUFFER];
bool CD_pipid_delta_vz_check[BUFFER];
bool CD_pipid_all_check[BUFFER];

bool CD_pimid_default_PID_check[BUFFER];
bool CD_pimid_charge_check[BUFFER];
bool CD_pimid_beta_check[BUFFER];
bool CD_pimid_maximum_probability_check[BUFFER];
bool CD_pimid_delta_vz_check[BUFFER];
bool CD_pimid_all_check[BUFFER];

bool CD_Kpid_default_PID_check[BUFFER];
bool CD_Kpid_charge_check[BUFFER];
bool CD_Kpid_beta_check[BUFFER];
bool CD_Kpid_maximum_probability_check[BUFFER];
bool CD_Kpid_delta_vz_check[BUFFER];
bool CD_Kpid_all_check[BUFFER];

bool CD_Kmid_default_PID_check[BUFFER];
bool CD_Kmid_charge_check[BUFFER];
bool CD_Kmid_beta_check[BUFFER];
bool CD_Kmid_maximum_probability_check[BUFFER];
bool CD_Kmid_delta_vz_check[BUFFER];
bool CD_Kmid_all_check[BUFFER];


/// ////////////////////////////////////////////////////////////////////////////////////////////////////
/// PID histograms:

/// a) electron ID:

  const int FD_eid_cuts = 16;

  // CC cuts

  TH2F *hist_HTCC_theta_vs_phi[FD_eid_cuts];
  TH1F *hist_HTCC_nphe[FD_eid_cuts];
  TH2F *hist_HTCC_nphe_vs_sampling_fraction[FD_eid_cuts];

  // EC cuts

  TH2F *hist_EC_PCAL_vs_EC_ECAL[FD_eid_cuts];
  TH2F *hist_EC_outer_vs_EC_inner[FD_eid_cuts];

  TH2F *hist_EC_total_sampling_fraction_sec1[FD_eid_cuts];
  TH2F *hist_EC_total_sampling_fraction_sec2[FD_eid_cuts];
  TH2F *hist_EC_total_sampling_fraction_sec3[FD_eid_cuts];
  TH2F *hist_EC_total_sampling_fraction_sec4[FD_eid_cuts];
  TH2F *hist_EC_total_sampling_fraction_sec5[FD_eid_cuts];
  TH2F *hist_EC_total_sampling_fraction_sec6[FD_eid_cuts];

  TH2F *hist_EC_PCAL_sampling_fraction_sec1[FD_eid_cuts];
  TH2F *hist_EC_PCAL_sampling_fraction_sec2[FD_eid_cuts];
  TH2F *hist_EC_PCAL_sampling_fraction_sec3[FD_eid_cuts];
  TH2F *hist_EC_PCAL_sampling_fraction_sec4[FD_eid_cuts];
  TH2F *hist_EC_PCAL_sampling_fraction_sec5[FD_eid_cuts];
  TH2F *hist_EC_PCAL_sampling_fraction_sec6[FD_eid_cuts];

  TH2F *hist_EC_ECAL_sampling_fraction_sec1[FD_eid_cuts];
  TH2F *hist_EC_ECAL_sampling_fraction_sec2[FD_eid_cuts];
  TH2F *hist_EC_ECAL_sampling_fraction_sec3[FD_eid_cuts];
  TH2F *hist_EC_ECAL_sampling_fraction_sec4[FD_eid_cuts];
  TH2F *hist_EC_ECAL_sampling_fraction_sec5[FD_eid_cuts];
  TH2F *hist_EC_ECAL_sampling_fraction_sec6[FD_eid_cuts];

  TH2F *hist_EC_PCAL_hit_position[FD_eid_cuts];
  TH2F *hist_EC_inner_hit_position[FD_eid_cuts];
  TH2F *hist_EC_outer_hit_position[FD_eid_cuts];

  // DC cuts

  TH2F *hist_DC_hit_position_region1[FD_eid_cuts];
  TH2F *hist_DC_hit_position_region2[FD_eid_cuts];
  TH2F *hist_DC_hit_position_region3[FD_eid_cuts];

  TH1F *hist_DC_z_vertex_sec1[FD_eid_cuts];
  TH1F *hist_DC_z_vertex_sec2[FD_eid_cuts];
  TH1F *hist_DC_z_vertex_sec3[FD_eid_cuts];
  TH1F *hist_DC_z_vertex_sec4[FD_eid_cuts];
  TH1F *hist_DC_z_vertex_sec5[FD_eid_cuts];
  TH1F *hist_DC_z_vertex_sec6[FD_eid_cuts];


/// b) FTOF + additional cuts for charged hadrons:

  const int FD_hid_count = 60;

  TH2F *hist_DC_hit_position_region1_hadron[FD_hid_count];
  TH2F *hist_DC_hit_position_region2_hadron[FD_hid_count];
  TH2F *hist_DC_hit_position_region3_hadron[FD_hid_count];
  TH2F *hist_EC_outer_vs_EC_inner_hadron[FD_hid_count];

  TH2F *hist_beta_vs_p[FD_hid_count];
  TH2F *hist_beta_vs_p_sec1[FD_hid_count];
  TH2F *hist_beta_vs_p_sec2[FD_hid_count];
  TH2F *hist_beta_vs_p_sec3[FD_hid_count];
  TH2F *hist_beta_vs_p_sec4[FD_hid_count];
  TH2F *hist_beta_vs_p_sec5[FD_hid_count];
  TH2F *hist_beta_vs_p_sec6[FD_hid_count];

  TH2F *hist_delta_beta_vs_p[FD_hid_count];
  TH1F *hist_delta_beta[FD_hid_count];
  TH2F *hist_tofmass_vs_p[FD_hid_count];
  TH1F *hist_tofmass[FD_hid_count];
  TH1F *hist_delta_vz[FD_hid_count];


/// CTOF + additional cuts for charged hadrons:

  const int FD_hid_CD_count = 60;

  TH2F *hist_CD_beta_vs_p[FD_hid_CD_count];
  TH1F *hist_CD_delta_vz[FD_hid_CD_count];


/// c) photon ID:

  const int FD_photid_count = 8;

  TH2F *hist_beta_vs_p_phot[FD_photid_count];
  TH1F *hist_beta_phot[FD_photid_count];
  TH2F *hist_EC_sampling_fraction_phot[FD_photid_count];
  TH2F *hist_EC_PCAL_vs_EC_ECAL_phot[FD_photid_count];
  TH2F *hist_EC_PCAL_hit_position_phot[FD_photid_count];


////////////////////////////////////////////////////////////////////
///  FT plots

  const int FT_pid_count = 20;

  TH2F *hist_FT_FTCAL_energy_vs_radius[FT_pid_count];
  TH2F *hist_FT_FTCAL_hit_position[FT_pid_count];
  TH2F *hist_FT_FTTRK_hit_position[FT_pid_count];
  TH2F *hist_FT_FTHODO_hit_position[FT_pid_count];
  TH1F *hist_FT_beta[FT_pid_count];


/// /////////////////////////////////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////////////////////////////////
/// define functions

/// electron ID:

// CC cuts

TH2F *create_hist_HTCC_theta_vs_phi(int cutnum);
TH1F *create_hist_HTCC_nphe(int cutnum);
TH2F *create_hist_HTCC_nphe_vs_sampling_fraction(int cutnum);

// EC cuts

TH2F *create_hist_EC_PCAL_vs_EC_ECAL(int cutnum);
TH2F *create_hist_EC_outer_vs_EC_inner(int cutnum);

TH2F *create_hist_EC_total_sampling_fraction_sec1(int cutnum);
TH2F *create_hist_EC_total_sampling_fraction_sec2(int cutnum);
TH2F *create_hist_EC_total_sampling_fraction_sec3(int cutnum);
TH2F *create_hist_EC_total_sampling_fraction_sec4(int cutnum);
TH2F *create_hist_EC_total_sampling_fraction_sec5(int cutnum);
TH2F *create_hist_EC_total_sampling_fraction_sec6(int cutnum);

TH2F *create_hist_EC_PCAL_sampling_fraction_sec1(int cutnum);
TH2F *create_hist_EC_PCAL_sampling_fraction_sec2(int cutnum);
TH2F *create_hist_EC_PCAL_sampling_fraction_sec3(int cutnum);
TH2F *create_hist_EC_PCAL_sampling_fraction_sec4(int cutnum);
TH2F *create_hist_EC_PCAL_sampling_fraction_sec5(int cutnum);
TH2F *create_hist_EC_PCAL_sampling_fraction_sec6(int cutnum);

TH2F *create_hist_EC_ECAL_sampling_fraction_sec1(int cutnum);
TH2F *create_hist_EC_ECAL_sampling_fraction_sec2(int cutnum);
TH2F *create_hist_EC_ECAL_sampling_fraction_sec3(int cutnum);
TH2F *create_hist_EC_ECAL_sampling_fraction_sec4(int cutnum);
TH2F *create_hist_EC_ECAL_sampling_fraction_sec5(int cutnum);
TH2F *create_hist_EC_ECAL_sampling_fraction_sec6(int cutnum);

TH2F *create_hist_EC_PCAL_hit_position(int cutnum);
TH2F *create_hist_EC_inner_hit_position(int cutnum);
TH2F *create_hist_EC_outer_hit_position(int cutnum);

// DC cuts

TH2F *create_hist_DC_hit_position_region1(int cutnum);
TH2F *create_hist_DC_hit_position_region2(int cutnum);
TH2F *create_hist_DC_hit_position_region3(int cutnum);

TH1F *create_hist_DC_z_vertex_sec1(int cutnum);
TH1F *create_hist_DC_z_vertex_sec2(int cutnum);
TH1F *create_hist_DC_z_vertex_sec3(int cutnum);
TH1F *create_hist_DC_z_vertex_sec4(int cutnum);
TH1F *create_hist_DC_z_vertex_sec5(int cutnum);
TH1F *create_hist_DC_z_vertex_sec6(int cutnum);

// TOF + others for charged hadrons

TH2F *create_DC_hit_position_region1_hadron(int cutnum);
TH2F *create_DC_hit_position_region2_hadron(int cutnum);
TH2F *create_DC_hit_position_region3_hadron(int cutnum);
TH2F *create_hist_EC_outer_vs_EC_inner_hadron(int cutnum);

TH2F *create_hist_beta_vs_p(int cutnum);
TH2F *create_hist_beta_vs_p_sec1(int cutnum);
TH2F *create_hist_beta_vs_p_sec2(int cutnum);
TH2F *create_hist_beta_vs_p_sec3(int cutnum);
TH2F *create_hist_beta_vs_p_sec4(int cutnum);
TH2F *create_hist_beta_vs_p_sec5(int cutnum);
TH2F *create_hist_beta_vs_p_sec6(int cutnum);

TH2F *create_hist_delta_beta_vs_p(int cutnum);
TH1F *create_hist_delta_beta(int cutnum);
TH2F *create_hist_tofmass_vs_p(int cutnum);
TH1F *create_hist_tofmass(int cutnum);
TH1F *create_hist_delta_vz(int cutnum);

// CD charged hadrons

TH2F *create_hist_CD_beta_vs_p(int cutnum);
TH1F *create_hist_CD_delta_vz(int cutnum);

// TOF for photons

TH2F *create_hist_beta_vs_p_phot(int cutnum);
TH1F *create_hist_beta_phot(int cutnum);
TH2F *create_hist_EC_sampling_fraction_phot(int cutnum);
TH2F *create_hist_EC_PCAL_vs_EC_ECAL_phot(int cutnum);
TH2F *create_hist_EC_PCAL_hit_position_phot(int cutnum);

////////////////////////////////////////////////////////////////
// FT plots

TH2F *create_hist_FT_FTCAL_energy_vs_radius(int cutnum);
TH2F *create_hist_FT_FTCAL_hit_position(int cutnum);
TH2F *create_hist_FT_FTTRK_hit_position(int cutnum);
TH2F *create_hist_FT_FTHODO_hit_position(int cutnum);
TH1F *create_hist_FT_beta(int cutnum);


// get event properties from the EVNT bank

void get_event_properties(void);

// raw particle assignment

void assign_particles(void);

// particle selection:

void select_electron(int run);
void select_proton(int run);
void select_neutron(int run);
void select_pip(int run);
void select_pim(int run);
void select_Kplus(int run);
void select_Kminus(int run);
void select_photon(int run);

// momentum correction:

TLorentzVector correct_lepton_negative(double thetaeld, double phield, double pel, int secte);
TVector3 correct_hadron_positive(double thetahd, double phihd, double ph, int secth);
TVector3 correct_hadron_negative(double thetahd, double phihd, double ph, int secth);

void correct_electron(void);
void correct_proton(void);
void correct_neutron(void);
void correct_pip(void);
void correct_pim(void);
void correct_Kplus(void);
void correct_Kminus(void);
void correct_photon(void);

// reconstruct neutrals from photons

void create_neutrals(void);

// fill output variables:

void fill_output_vector_electron(void);
void fill_output_vector_proton(void);
void fill_output_vector_neutron(void);
void fill_output_vector_pip(void);
void fill_output_vector_pim(void);
void fill_output_vector_Kplus(void);
void fill_output_vector_Kminus(void);
void fill_output_vector_photon(void);
void fill_output_vector_MC(void);

// define cuts:

bool basic_FTOF_cut(int j);

// PID checks:

bool Track_Quality_cut(int j);

bool ele_default_PID_cut(int j);
bool ele_charge_cut(int j);
bool CC_nphe_cut(int j);
bool EC_outer_vs_EC_inner_cut(int j);
bool EC_sampling_fraction_cut(int j);
bool EC_hit_position_fiducial_cut(int j);
bool DC_hit_position_region1_fiducial_cut(int j);
bool DC_hit_position_region2_fiducial_cut(int j);
bool DC_hit_position_region3_fiducial_cut(int j);
bool DC_z_vertex_cut(int j);

bool DC_hit_position_region1_fiducial_cut_hadrons_positive(int j);
bool DC_hit_position_region2_fiducial_cut_hadrons_positive(int j);
bool DC_hit_position_region3_fiducial_cut_hadrons_positive(int j);

bool DC_hit_position_region1_fiducial_cut_hadrons_negative(int j);
bool DC_hit_position_region2_fiducial_cut_hadrons_negative(int j);
bool DC_hit_position_region3_fiducial_cut_hadrons_negative(int j);

bool prot_default_PID_cut(int j);
bool prot_charge_cut(int j);
bool prot_beta_cut(int j, int run);
bool prot_delta_beta_cut(int j, int run);
bool prot_tofmass_cut(int j, int run);
bool prot_delta_vz_cut(int j);

bool neutr_default_PID_cut(int j);
bool neutr_charge_cut(int j);
bool neutr_beta_cut(int j, int run);
bool neutr_delta_beta_cut(int j, int run);
bool neutr_tofmass_cut(int j, int run);
bool neutr_delta_vz_cut(int j);

bool pip_default_PID_cut(int j);
bool pip_charge_cut(int j);
bool pip_beta_cut(int j, int run);
bool pip_delta_beta_cut(int j, int run);
bool pip_tofmass_cut(int j, int run);
bool pip_delta_vz_cut(int j);

bool pim_default_PID_cut(int j);
bool pim_charge_cut(int j);
bool pim_ele_reject_cut(int j);
bool pim_EC_outer_vs_EC_inner_cut(int j);
bool pim_beta_cut(int j, int run);
bool pim_delta_beta_cut(int j, int run);
bool pim_tofmass_cut(int j, int run);
bool pim_delta_vz_cut(int j);

bool Kp_default_PID_cut(int j);
bool Kp_charge_cut(int j);
bool Kp_beta_cut(int j, int run);
bool Kp_delta_beta_cut(int j, int run);
bool Kp_tofmass_cut(int j, int run);
bool Kp_delta_vz_cut(int j);

bool Km_default_PID_cut(int j);
bool Km_charge_cut(int j);
bool Km_ele_reject_cut(int j);
bool Km_EC_outer_vs_EC_inner_cut(int j);
bool Km_beta_cut(int j, int run);
bool Km_delta_beta_cut(int j, int run);
bool Km_tofmass_cut(int j, int run);
bool Km_delta_vz_cut(int j);

bool maximum_probability_cut(int j, int hypothesis, double conflvl, double anticonflvl, int run);

bool phot_default_PID_cut(int j);
bool phot_charge_cut(int j);
bool phot_beta_cut(int j, int run);
bool phot_EC_sampling_fraction_cut(int j);
bool phot_EC_outer_vs_EC_inner_cut(int j);
bool phot_EC_hit_position_fiducial_cut(int j);

// FT

bool FT_eid_charge_cut(int j);
bool FT_eid_PID_cut(int j);
bool FT_eid_FTCAL_fiducial_cut(int j);
bool FT_eid_FTTRK_fiducial_cut(int j);
bool FT_eid_FTHODO_fiducial_cut(int j);
bool FT_eid_energy_vs_radius_cut(int j);

bool FT_photid_charge_cut(int j);
bool FT_photid_PID_cut(int j);
bool FT_photid_FTCAL_fiducial_cut(int j);
bool FT_photid_beta_cut(int j, int run);

// CD

bool CD_prot_default_PID_cut(int j);
bool CD_prot_charge_cut(int j);
bool CD_prot_beta_cut(int j, int run);
bool CD_prot_delta_vz_cut(int j);

bool CD_neutr_default_PID_cut(int j);
bool CD_neutr_charge_cut(int j);
bool CD_neutr_beta_cut(int j, int run);
bool CD_neutr_delta_vz_cut(int j);

bool CD_pip_default_PID_cut(int j);
bool CD_pip_charge_cut(int j);
bool CD_pip_beta_cut(int j, int run);
bool CD_pip_delta_vz_cut(int j);

bool CD_pim_default_PID_cut(int j);
bool CD_pim_charge_cut(int j);
bool CD_pim_beta_cut(int j, int run);
bool CD_pim_delta_vz_cut(int j);

bool CD_Kp_default_PID_cut(int j);
bool CD_Kp_charge_cut(int j);
bool CD_Kp_beta_cut(int j, int run);
bool CD_Kp_delta_vz_cut(int j);

bool CD_Km_default_PID_cut(int j);
bool CD_Km_charge_cut(int j);
bool CD_Km_beta_cut(int j, int run);
bool CD_Km_delta_vz_cut(int j);

bool CD_maximum_probability_cut(int j, int hypothesis, double conflvl, double anticonflvl, int run);

// SC timing correction

double electron_sct_corr(int j, int run);
double hadron_sct_corr(int j, int run);
double hadron_sct_corr_central(int j, int run);
bool good_sc_paddle(int j);

// additional calculation functions:

double GetTheta(int j);
double GetPhi(int j);
TVector3 GetUVWVector(int j);
double GetTOFmass2(int j, int run);
double GetTOFmass2_CD(int j, int run);
double Getdvz(int j);
double Beta_charged(int j, int run);
double Beta_charged_central(int j, int run);
double Beta_charged_FT(int j, int run);
double Beta_neutral(int j, int run);
double Beta_neutral_FT(int j, int run);
double Get_Starttime(int j, int run);

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


/// /////////////////////////////////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// main
///

Int_t filter_clas12( Char_t *inFile, Char_t *outputfile, int run)
{
		
     const Char_t *inTree="clas12";
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
///  create output-txtfiles:
/// /////////////////////////////////////////////////////////////////////////////
    
     //ofstream outputFile_electrons_elastic("output/electrons_elastic.txt");
     //ofstream outputFile_ep_elastic("output/ep_elastic.txt");
     //ofstream outputFile_epip_missing_neutron("output/e_pip_X.txt");

/// /////////////////////////////////////////////////////////////////////////////    
///  create output-file and tree for saving histograms:
/// /////////////////////////////////////////////////////////////////////////////

    out = new TFile(outputfile, "RECREATE");

    TTree out_tree("out_tree","out_tree");
    out_tree.Branch("helicity", &helicity);
    out_tree.Branch("fcup", &fcup);
    out_tree.Branch("sectorE", &sectorE);
    out_tree.Branch("eventNumber",&electron_event_number);
    out_tree.Branch("p4_ele_px", &p4_ele_px);
    out_tree.Branch("p4_ele_py", &p4_ele_py);
    out_tree.Branch("p4_ele_pz", &p4_ele_pz);
    out_tree.Branch("p4_ele_E", &p4_ele_E);
    out_tree.Branch("p4_prot_px", &p4_prot_px);
    out_tree.Branch("p4_prot_py", &p4_prot_py);
    out_tree.Branch("p4_prot_pz", &p4_prot_pz);
    out_tree.Branch("p4_prot_E", &p4_prot_E);
    out_tree.Branch("p4_neutr_px", &p4_neutr_px);
    out_tree.Branch("p4_neutr_py", &p4_neutr_py);
    out_tree.Branch("p4_neutr_pz", &p4_neutr_pz);
    out_tree.Branch("p4_neutr_E", &p4_neutr_E);
    out_tree.Branch("p4_pip_px", &p4_pip_px);
    out_tree.Branch("p4_pip_py", &p4_pip_py);
    out_tree.Branch("p4_pip_pz", &p4_pip_pz);
    out_tree.Branch("p4_pip_E", &p4_pip_E);
    out_tree.Branch("p4_pim_px", &p4_pim_px);
    out_tree.Branch("p4_pim_py", &p4_pim_py);
    out_tree.Branch("p4_pim_pz", &p4_pim_pz);
    out_tree.Branch("p4_pim_E", &p4_pim_E);
    out_tree.Branch("p4_Kp_px", &p4_Kp_px);
    out_tree.Branch("p4_Kp_py", &p4_Kp_py);
    out_tree.Branch("p4_Kp_pz", &p4_Kp_pz);
    out_tree.Branch("p4_Kp_E", &p4_Kp_E);
    out_tree.Branch("p4_Km_px", &p4_Km_px);
    out_tree.Branch("p4_Km_py", &p4_Km_py);
    out_tree.Branch("p4_Km_pz", &p4_Km_pz);
    out_tree.Branch("p4_Km_E", &p4_Km_E);
    out_tree.Branch("p4_phot_px", &p4_phot_px);
    out_tree.Branch("p4_phot_py", &p4_phot_py);
    out_tree.Branch("p4_phot_pz", &p4_phot_pz);
    out_tree.Branch("p4_phot_E", &p4_phot_E);  
    out_tree.Branch("ele_det", &ele_det);
    out_tree.Branch("prot_det", &prot_det);
    out_tree.Branch("neutr_det", &neutr_det);
    out_tree.Branch("pip_det", &pip_det);
    out_tree.Branch("pim_det", &pim_det);
    out_tree.Branch("Kp_det", &Kp_det);
    out_tree.Branch("Km_det", &Km_det);
    out_tree.Branch("phot_det", &phot_det);

    out_tree.Branch("gen_event_helicity", &MC_helicity);
    out_tree.Branch("gen_event_npart", &MC_Npart);
    out_tree.Branch("gen_event_ebeam", &MC_Ebeam);
    out_tree.Branch("gen_event_weight", &MC_weight);
    out_tree.Branch("gen_pid", &vMC_Particle_pid);
    out_tree.Branch("gen_px", &vMC_Particle_px);
    out_tree.Branch("gen_py", &vMC_Particle_py);
    out_tree.Branch("gen_pz", &vMC_Particle_pz);
    out_tree.Branch("gen_vx", &vMC_Particle_vx);
    out_tree.Branch("gen_vy", &vMC_Particle_vy);
    out_tree.Branch("gen_vz", &vMC_Particle_vz);


/// ////////////////////////////////////////////////////////////////////////////
///  assign the tree branches to the variables:
/// ////////////////////////////////////////////////////////////////////////////

   anaTree->SetBranchAddress("REC_Event_NRUN", &vNRUN);
   anaTree->SetBranchAddress("REC_Event_NEVENT", &vNEVENT);
   anaTree->SetBranchAddress("REC_Event_EVNTime", &vEVNTime);
   anaTree->SetBranchAddress("REC_Event_TYPE", &vTYPE);
   anaTree->SetBranchAddress("REC_Event_TRG", &vTRG);
   anaTree->SetBranchAddress("REC_Event_BCG", &vBCG);
   anaTree->SetBranchAddress("REC_Event_STTime", &vSTTime);
   anaTree->SetBranchAddress("REC_Event_RFTime", &vRFTime);
   anaTree->SetBranchAddress("REC_Event_Helic", &vHelic);

   anaTree->SetBranchAddress("REC_Particle_pid", &vpart_pid);
   anaTree->SetBranchAddress("REC_Particle_charge", &vpart_charge);   
   anaTree->SetBranchAddress("REC_Particle_px", &vpart_px);
   anaTree->SetBranchAddress("REC_Particle_py", &vpart_py);
   anaTree->SetBranchAddress("REC_Particle_pz", &vpart_pz);
   anaTree->SetBranchAddress("REC_Particle_vx", &vpart_vx);
   anaTree->SetBranchAddress("REC_Particle_vy", &vpart_vy);
   anaTree->SetBranchAddress("REC_Particle_vz", &vpart_vz);
   anaTree->SetBranchAddress("REC_Particle_beta", &vpart_beta);
   anaTree->SetBranchAddress("REC_Particle_status", &vpart_status);

   anaTree->SetBranchAddress("REC_Calorimeter_pindex", &vCal_pindex);
   anaTree->SetBranchAddress("REC_Calorimeter_detector", &vCal_detector);   
   anaTree->SetBranchAddress("REC_Calorimeter_sector", &vCal_sector);
   anaTree->SetBranchAddress("REC_Calorimeter_layer", &vCal_layer);  
   anaTree->SetBranchAddress("REC_Calorimeter_energy", &vCal_energy);
   anaTree->SetBranchAddress("REC_Calorimeter_time", &vCal_time);   
   anaTree->SetBranchAddress("REC_Calorimeter_path", &vCal_path);
   anaTree->SetBranchAddress("REC_Calorimeter_x", &vCal_x);   
   anaTree->SetBranchAddress("REC_Calorimeter_y", &vCal_y);
   anaTree->SetBranchAddress("REC_Calorimeter_z", &vCal_z);    
   anaTree->SetBranchAddress("REC_Calorimeter_lu", &vCal_lu);   
   anaTree->SetBranchAddress("REC_Calorimeter_lv", &vCal_lv);
   anaTree->SetBranchAddress("REC_Calorimeter_lw", &vCal_lw);   

   anaTree->SetBranchAddress("REC_Cherenkov_pindex", &vCC_pindex);
   anaTree->SetBranchAddress("REC_Cherenkov_detector", &vCC_detector);   
   anaTree->SetBranchAddress("REC_Cherenkov_sector", &vCC_sector);
   anaTree->SetBranchAddress("REC_Cherenkov_nphe", &vCC_nphe);  
   anaTree->SetBranchAddress("REC_Cherenkov_time", &vCC_time);   
   anaTree->SetBranchAddress("REC_Cherenkov_path", &vCC_path);
   anaTree->SetBranchAddress("REC_Cherenkov_theta", &vCC_theta);
   anaTree->SetBranchAddress("REC_Cherenkov_phi", &vCC_phi); 

   anaTree->SetBranchAddress("REC_ForwardTagger_pindex", &vFT_pindex);
   anaTree->SetBranchAddress("REC_ForwardTagger_detector", &vFT_detector);   
   anaTree->SetBranchAddress("REC_ForwardTagger_energy", &vFT_energy);
   anaTree->SetBranchAddress("REC_ForwardTagger_time", &vFT_time);   
   anaTree->SetBranchAddress("REC_ForwardTagger_path", &vFT_path);
   anaTree->SetBranchAddress("REC_ForwardTagger_x", &vFT_x);   
   anaTree->SetBranchAddress("REC_ForwardTagger_y", &vFT_y);
   anaTree->SetBranchAddress("REC_ForwardTagger_z", &vFT_z);    
   anaTree->SetBranchAddress("REC_ForwardTagger_radius", &vFT_radius);
   anaTree->SetBranchAddress("REC_ForwardTagger_size", &vFT_size); 

   anaTree->SetBranchAddress("REC_Scintillator_pindex", &vSC_pindex);
   anaTree->SetBranchAddress("REC_Scintillator_detector", &vSC_detector);   
   anaTree->SetBranchAddress("REC_Scintillator_sector", &vSC_sector);
   anaTree->SetBranchAddress("REC_Scintillator_layer", &vSC_layer);  
   anaTree->SetBranchAddress("REC_Scintillator_component", &vSC_component); 
   anaTree->SetBranchAddress("REC_Scintillator_energy", &vSC_energy);
   anaTree->SetBranchAddress("REC_Scintillator_time", &vSC_time);   
   anaTree->SetBranchAddress("REC_Scintillator_path", &vSC_path);
 
   anaTree->SetBranchAddress("REC_Track_pindex", &vTRK_pindex);
   anaTree->SetBranchAddress("REC_Track_detector", &vTRK_detector);  
   anaTree->SetBranchAddress("REC_Track_sector", &vTRK_sector);   

   anaTree->SetBranchAddress("REC_Traj_pindex", &vTraj_pindex); 
   anaTree->SetBranchAddress("REC_Traj_detId", &vTraj_detID);   
   anaTree->SetBranchAddress("REC_Traj_x", &vTraj_x);   
   anaTree->SetBranchAddress("REC_Traj_y", &vTraj_y);
   anaTree->SetBranchAddress("REC_Traj_z", &vTraj_z); 
   anaTree->SetBranchAddress("REC_Traj_cx", &vTraj_cx);   
   anaTree->SetBranchAddress("REC_Traj_cy", &vTraj_cy);
   anaTree->SetBranchAddress("REC_Traj_cz", &vTraj_cz);   

if(simulation == true){

  //anaTree->SetBranchAddress("MC_Header_helicity", &vMC_Header_helicity);
  // anaTree->SetBranchAddress("MC_Event_npart", &vMC_Event_npart);
  //anaTree->SetBranchAddress("MC_Event_ebeam", &vMC_Event_ebeam);
  // anaTree->SetBranchAddress("MC_Event_weight", &vMC_Event_weight);

   anaTree->SetBranchAddress("MC_Particle_pid", &vMC_Particle_pid);  
   anaTree->SetBranchAddress("MC_Particle_px", &vMC_Particle_px);
   anaTree->SetBranchAddress("MC_Particle_py", &vMC_Particle_py);
   anaTree->SetBranchAddress("MC_Particle_pz", &vMC_Particle_pz);
   anaTree->SetBranchAddress("MC_Particle_vx", &vMC_Particle_vx);
   anaTree->SetBranchAddress("MC_Particle_vy", &vMC_Particle_vy);
   anaTree->SetBranchAddress("MC_Particle_vz", &vMC_Particle_vz);

   /*anaTree->SetBranchAddress("MC_Lund_pid", &vMC_Lund_pid); 
   anaTree->SetBranchAddress("MC_Lund_mass", &vMC_Lund_mass);
   anaTree->SetBranchAddress("MC_Lund_E", &vMC_Lund_E); 
   anaTree->SetBranchAddress("MC_Lund_px", &vMC_Lund_px);
   anaTree->SetBranchAddress("MC_Lund_py", &vMC_Lund_py);
   anaTree->SetBranchAddress("MC_Lund_pz", &vMC_Lund_pz);
   anaTree->SetBranchAddress("MC_Lund_vx", &vMC_Lund_vx);
   anaTree->SetBranchAddress("MC_Lund_vy", &vMC_Lund_vy);
   anaTree->SetBranchAddress("MC_Lund_vz", &vMC_Lund_vz);
   */
}


/// ///////////////////////////////////////////////////////////////
///  reset cut statistics
/// ///////////////////////////////////////////////////////////////

neg_part_count = 0;
pos_part_count = 0;
neut_part_count = 0;

FD_eid_default_PID_pass = 0;
FD_eid_charge_pass = 0;
FD_eid_EC_outer_vs_EC_inner_pass = 0;
FD_eid_EC_sampling_fraction_pass = 0;
FD_eid_EC_hit_position_fiducial_pass = 0;
FD_eid_DC_hit_position_region1_fiducial_pass = 0;
FD_eid_DC_hit_position_region2_fiducial_pass = 0;
FD_eid_DC_hit_position_region3_fiducial_pass = 0;
FD_eid_DC_z_vertex_pass = 0;
FD_eid_all_pass = 0;

FD_protid_default_PID_pass = 0;
FD_protid_charge_pass = 0;
FD_protid_DC_hit_position_region1_fiducial_pass = 0;
FD_protid_DC_hit_position_region2_fiducial_pass = 0;
FD_protid_DC_hit_position_region3_fiducial_pass = 0;
FD_protid_beta_pass = 0;
FD_protid_delta_beta_pass = 0;
FD_protid_tofmass_pass = 0;
FD_protid_maximum_probability_pass = 0;
FD_protid_delta_vz_pass = 0;
FD_protid_all_pass = 0;
FD_neutrid_default_PID_pass = 0;
FD_neutrid_charge_pass = 0;
FD_neutrid_beta_pass = 0;
FD_neutrid_delta_beta_pass = 0;
FD_neutrid_tofmass_pass = 0;
FD_neutrid_delta_vz_pass = 0;
FD_neutrid_all_pass = 0;
FD_pipid_default_PID_pass = 0;
FD_pipid_charge_pass = 0;
FD_pipid_DC_hit_position_region1_fiducial_pass = 0;
FD_pipid_DC_hit_position_region2_fiducial_pass = 0;
FD_pipid_DC_hit_position_region3_fiducial_pass = 0;
FD_pipid_beta_pass = 0;
FD_pipid_delta_beta_pass = 0;
FD_pipid_tofmass_pass = 0;
FD_pipid_maximum_probability_pass = 0;
FD_pipid_delta_vz_pass = 0;
FD_pipid_all_pass = 0;
FD_pimid_default_PID_pass = 0;
FD_pimid_charge_pass = 0;
FD_pimid_DC_hit_position_region1_fiducial_pass = 0;
FD_pimid_DC_hit_position_region2_fiducial_pass = 0;
FD_pimid_DC_hit_position_region3_fiducial_pass = 0;
FD_pimid_beta_pass = 0;
FD_pimid_delta_beta_pass = 0;
FD_pimid_tofmass_pass = 0;
FD_pimid_maximum_probability_pass = 0;
FD_pimid_delta_vz_pass = 0;
FD_pimid_all_pass = 0;
FD_Kpid_default_PID_pass = 0;
FD_Kpid_charge_pass = 0;
FD_Kpid_DC_hit_position_region1_fiducial_pass = 0;
FD_Kpid_DC_hit_position_region2_fiducial_pass = 0;
FD_Kpid_DC_hit_position_region3_fiducial_pass = 0;
FD_Kpid_beta_pass = 0;
FD_Kpid_delta_beta_pass = 0;
FD_Kpid_tofmass_pass = 0;
FD_Kpid_maximum_probability_pass = 0;
FD_Kpid_delta_vz_pass = 0;
FD_Kpid_all_pass = 0;
FD_Kmid_default_PID_pass = 0;
FD_Kmid_charge_pass = 0;
FD_Kmid_DC_hit_position_region1_fiducial_pass = 0;
FD_Kmid_DC_hit_position_region2_fiducial_pass = 0;
FD_Kmid_DC_hit_position_region3_fiducial_pass = 0;
FD_Kmid_beta_pass = 0;
FD_Kmid_delta_beta_pass = 0;
FD_Kmid_tofmass_pass = 0;
FD_Kmid_maximum_probability_pass = 0;
FD_Kmid_delta_vz_pass = 0;
FD_Kmid_all_pass = 0;

CD_protid_default_PID_pass = 0;
CD_protid_charge_pass = 0;
CD_protid_beta_pass = 0;
CD_protid_maximum_probability_pass = 0;
CD_protid_delta_vz_pass = 0;
CD_protid_all_pass = 0;
CD_neutrid_default_PID_pass = 0;
CD_neutrid_charge_pass = 0;
CD_neutrid_beta_pass = 0;
CD_neutrid_delta_vz_pass = 0;
CD_neutrid_all_pass = 0;
CD_pipid_default_PID_pass = 0;
CD_pipid_charge_pass = 0;
CD_pipid_beta_pass = 0;
CD_pipid_maximum_probability_pass = 0;
CD_pipid_delta_vz_pass = 0;
CD_pipid_all_pass = 0;
CD_pimid_default_PID_pass = 0;
CD_pimid_charge_pass = 0;
CD_pimid_beta_pass = 0;
CD_pimid_maximum_probability_pass = 0;
CD_pimid_delta_vz_pass = 0;
CD_pimid_all_pass = 0;
CD_Kpid_default_PID_pass = 0;
CD_Kpid_charge_pass = 0;
CD_Kpid_beta_pass = 0;
CD_Kpid_maximum_probability_pass = 0;
CD_Kpid_delta_vz_pass = 0;
CD_Kpid_all_pass = 0;
CD_Kmid_default_PID_pass = 0;
CD_Kmid_charge_pass = 0;
CD_Kmid_beta_pass = 0;
CD_Kmid_maximum_probability_pass = 0;
CD_Kmid_delta_vz_pass = 0;
CD_Kmid_all_pass = 0;

FD_photid_default_PID_pass = 0;
FD_photid_charge_pass = 0;
FD_photid_beta_pass = 0;
FD_photid_EC_sampling_fraction_pass = 0;
FD_photid_EC_hit_position_fiducial_pass = 0;
FD_photid_all_pass = 0;

FT_eid_charge_pass = 0;
FT_eid_PID_pass = 0;
FT_eid_FTCAL_fiducial_pass = 0;
FT_eid_FTTRK_fiducial_pass = 0;
FT_eid_FTHODO_fiducial_pass = 0;
FT_eid_energy_vs_radius_pass = 0;
FT_eid_all_pass = 0;

FT_photid_charge_pass = 0;
FT_photid_PID_pass = 0;
FT_photid_FTCAL_fiducial_pass = 0;
FT_photid_beta_pass = 0;
FT_photid_all_pass = 0;


/// ///////////////////////////////////////////////////////////////
///  create histograms
/// ///////////////////////////////////////////////////////////////

out->mkdir("event_information");				
out->cd ("event_information");

  TH1F *hist_run_number;
  TH1F *hist_number_of_events;
  TH1F *hist_event_type;
  TH1F *hist_trigger;
  TH1F *hist_helicity;
  TH1F *hist_eventtime;
  TH1F *hist_faraday_cup;
  TH1F *hist_event_starttime;
  TH1F *hist_RF_time;
  TH1F *hist_status;

  hist_run_number = new TH1F("hist_run_number", "run number", 10000, 0, 10000);   
  hist_run_number->GetXaxis()->SetTitle("run number");
  hist_run_number->GetYaxis()->SetTitle("counts");
  hist_number_of_events = new TH1F("hist_number_of_events", "number of events", 1000, 0, 100000000);   
  hist_number_of_events->GetXaxis()->SetTitle("number of events");
  hist_number_of_events->GetYaxis()->SetTitle("counts");
  hist_event_type = new TH1F("hist_event_type", "event type", 11, -0.5, 10.5);   
  hist_event_type->GetXaxis()->SetTitle("event type");
  hist_event_type->GetYaxis()->SetTitle("counts");
  hist_trigger = new TH1F("hist_trigger", "trigger", 1000, 0, 1000);   
  hist_trigger->GetXaxis()->SetTitle("trigger");
  hist_trigger->GetYaxis()->SetTitle("counts");
  hist_helicity = new TH1F("hist_helicity", "helicity", 5, -2.5, 2.5);   
  hist_helicity->GetXaxis()->SetTitle("helicity");
  hist_helicity->GetYaxis()->SetTitle("counts");
  hist_eventtime = new TH1F("hist_eventtime", "event time", 1000, 0, 10000);   
  hist_eventtime->GetXaxis()->SetTitle("event time /ns");
  hist_eventtime->GetYaxis()->SetTitle("counts");
  hist_faraday_cup = new TH1F("hist_faraday_cup", "faraday cup", 500, 0, 5000);   
  hist_faraday_cup->GetXaxis()->SetTitle("faraday cup");
  hist_faraday_cup->GetYaxis()->SetTitle("counts");
  hist_event_starttime = new TH1F("hist_event_starttime", "event starttime", 600, 0, 300);   
  hist_event_starttime->GetXaxis()->SetTitle("event starttime /ns");
  hist_event_starttime->GetYaxis()->SetTitle("counts");
  hist_RF_time = new TH1F("hist_RF_time", "RF time", 800, -50, 150);
  hist_RF_time->GetXaxis()->SetTitle("RF time /ns");
  hist_RF_time->GetYaxis()->SetTitle("counts");
  hist_status = new TH1F("hist_status", "status", 600, 0, 6000);
  hist_status->GetXaxis()->SetTitle("status");
  hist_status->GetYaxis()->SetTitle("counts");



out->mkdir("FD_PID_electron_FTOF_plots");				
out->cd ("FD_PID_electron_FTOF_plots");

  TH2F *hist_FTOF_hit_position;

  hist_FTOF_hit_position = new TH2F("hist_FTOF_hit_position", "FTOF hit position", 500, -500, 500, 500, -500, 500);   
  hist_FTOF_hit_position->GetXaxis()->SetTitle("x /cm");
  hist_FTOF_hit_position->GetYaxis()->SetTitle("y /cm");


out->mkdir("FD_PID_electron_HTCC_plots");				
out->cd ("FD_PID_electron_HTCC_plots");

for(Int_t i = 0; i < 11; i++){ 

  create_hist_HTCC_theta_vs_phi(i);
  create_hist_HTCC_nphe(i);
  create_hist_HTCC_nphe_vs_sampling_fraction(i);
}


out->mkdir("FD_PID_electron_EC_plots");				
out->cd ("FD_PID_electron_EC_plots");

for(Int_t i = 0; i < 11; i++){ 
  create_hist_EC_PCAL_vs_EC_ECAL(i);
  create_hist_EC_outer_vs_EC_inner(i);

  create_hist_EC_total_sampling_fraction_sec1(i);
  create_hist_EC_total_sampling_fraction_sec2(i);
  create_hist_EC_total_sampling_fraction_sec3(i);
  create_hist_EC_total_sampling_fraction_sec4(i);
  create_hist_EC_total_sampling_fraction_sec5(i);
  create_hist_EC_total_sampling_fraction_sec6(i);

  create_hist_EC_PCAL_sampling_fraction_sec1(i);
  create_hist_EC_PCAL_sampling_fraction_sec2(i);
  create_hist_EC_PCAL_sampling_fraction_sec3(i);
  create_hist_EC_PCAL_sampling_fraction_sec4(i);
  create_hist_EC_PCAL_sampling_fraction_sec5(i);
  create_hist_EC_PCAL_sampling_fraction_sec6(i);

  create_hist_EC_ECAL_sampling_fraction_sec1(i);
  create_hist_EC_ECAL_sampling_fraction_sec2(i);
  create_hist_EC_ECAL_sampling_fraction_sec3(i);
  create_hist_EC_ECAL_sampling_fraction_sec4(i);
  create_hist_EC_ECAL_sampling_fraction_sec5(i);
  create_hist_EC_ECAL_sampling_fraction_sec6(i);

  create_hist_EC_PCAL_hit_position(i);
  create_hist_EC_inner_hit_position(i);
  create_hist_EC_outer_hit_position(i);
}


out->mkdir("FD_PID_electron_DC_plots");				
out->cd ("FD_PID_electron_DC_plots");

for(Int_t i = 0; i < 11; i++){ 

  create_hist_DC_hit_position_region1(i);
  create_hist_DC_hit_position_region2(i);
  create_hist_DC_hit_position_region3(i);

  create_hist_DC_z_vertex_sec1(i);
  create_hist_DC_z_vertex_sec2(i);
  create_hist_DC_z_vertex_sec3(i);
  create_hist_DC_z_vertex_sec4(i);
  create_hist_DC_z_vertex_sec5(i);
  create_hist_DC_z_vertex_sec6(i);
}

TH2F *hist_DC_hit_position_region2_cut5a;
hist_DC_hit_position_region2_cut5a = new TH2F("DC_hit_position_region2_cut_05a", "DC_hit_position_region2_cut_05a", 1000,-500,500, 1000,-500,500);   
hist_DC_hit_position_region2_cut5a->GetXaxis()->SetTitle("x /cm");
hist_DC_hit_position_region2_cut5a->GetYaxis()->SetTitle("y /cm");




out->mkdir("FD_PID_hadron_DC_fiducial_plot");				
out->cd ("FD_PID_hadron_DC_fiducial_plot");

for(Int_t i = 0; i < 60; i++){ create_DC_hit_position_region1_hadron(i);}
for(Int_t i = 0; i < 60; i++){ create_DC_hit_position_region2_hadron(i);}
for(Int_t i = 0; i < 60; i++){ create_DC_hit_position_region3_hadron(i);}

TH2F *hist_DC_hit_position_region2_hadron_cut_02a;
TH2F *hist_DC_hit_position_region2_hadron_cut_12a;
TH2F *hist_DC_hit_position_region2_hadron_cut_22a;
TH2F *hist_DC_hit_position_region2_hadron_cut_32a;
TH2F *hist_DC_hit_position_region2_hadron_cut_42a;
TH2F *hist_DC_hit_position_region2_hadron_cut_52a;

hist_DC_hit_position_region2_hadron_cut_02a = new TH2F("DC_hit_position_region2_hadron_cut_02a", "DC_hit_position_region2_hadron_cut_02a", 1000,-450,450, 1000,-450,450);   
hist_DC_hit_position_region2_hadron_cut_02a->GetXaxis()->SetTitle("x /cm");
hist_DC_hit_position_region2_hadron_cut_02a->GetYaxis()->SetTitle("y /cm");
hist_DC_hit_position_region2_hadron_cut_12a = new TH2F("DC_hit_position_region2_hadron_cut_12a", "DC_hit_position_region2_hadron_cut_12a", 1000,-450,450, 1000,-450,450);  
hist_DC_hit_position_region2_hadron_cut_12a->GetXaxis()->SetTitle("x /cm");
hist_DC_hit_position_region2_hadron_cut_12a->GetYaxis()->SetTitle("y /cm");
hist_DC_hit_position_region2_hadron_cut_22a = new TH2F("DC_hit_position_region2_hadron_cut_22a", "DC_hit_position_region2_hadron_cut_22a", 1000,-450,450, 1000,-450,450);
hist_DC_hit_position_region2_hadron_cut_22a->GetXaxis()->SetTitle("x /cm");
hist_DC_hit_position_region2_hadron_cut_22a->GetYaxis()->SetTitle("y /cm");
hist_DC_hit_position_region2_hadron_cut_32a = new TH2F("DC_hit_position_region2_hadron_cut_32a", "DC_hit_position_region2_hadron_cut_32a", 1000,-450,450, 1000,-450,450); 
hist_DC_hit_position_region2_hadron_cut_32a->GetXaxis()->SetTitle("x /cm");
hist_DC_hit_position_region2_hadron_cut_32a->GetYaxis()->SetTitle("y /cm");
hist_DC_hit_position_region2_hadron_cut_42a = new TH2F("DC_hit_position_region2_hadron_cut_42a", "DC_hit_position_region2_hadron_cut_42a", 1000,-450,450, 1000,-450,450);  
hist_DC_hit_position_region2_hadron_cut_42a->GetXaxis()->SetTitle("x /cm");
hist_DC_hit_position_region2_hadron_cut_42a->GetYaxis()->SetTitle("y /cm");
hist_DC_hit_position_region2_hadron_cut_52a = new TH2F("DC_hit_position_region2_hadron_cut_52a", "DC_hit_position_region2_hadron_cut_52a", 1000,-450,450, 1000,-450,450);  
hist_DC_hit_position_region2_hadron_cut_52a->GetXaxis()->SetTitle("x /cm");
hist_DC_hit_position_region2_hadron_cut_52a->GetYaxis()->SetTitle("y /cm");


out->mkdir("FD_PID_hadron_EC_plots");				
out->cd ("FD_PID_hadron_EC_plots");

for(Int_t i = 0; i < 60; i++){ create_hist_EC_outer_vs_EC_inner_hadron(i);}


out->mkdir("FD_PID_hadron_beta_plots");				
out->cd ("FD_PID_hadron_beta_plots");

for(Int_t i = 0; i < 60; i++){
  create_hist_beta_vs_p(i);
  create_hist_beta_vs_p_sec1(i);
  create_hist_beta_vs_p_sec2(i);
  create_hist_beta_vs_p_sec3(i);
  create_hist_beta_vs_p_sec4(i);
  create_hist_beta_vs_p_sec5(i);
  create_hist_beta_vs_p_sec6(i);
}   

for(Int_t i = 0; i < 60; i++){ create_hist_delta_beta_vs_p(i);}

for(Int_t i = 0; i < 60; i++){ create_hist_delta_beta(i);}


out->mkdir("FD_PID_hadron_beta_per_paddle");				
out->cd ("FD_PID_hadron_beta_per_paddle");

TH2F *hist_beta_vs_p_positive_layer1_comp[23]; 

for(Int_t i = 0; i < 23; i++){
  sprintf(name,"hist_beta_vs_p_positive_layer1_comp%01d", i+1);  
  sprintf(title," #beta vs p positive layer 1 component %01d", i+1);
  hist_beta_vs_p_positive_layer1_comp[i] = new TH2F(name, title, 600,0,Ebeam, 350, 0.0, 1.4);   
  hist_beta_vs_p_positive_layer1_comp[i]->GetXaxis()->SetTitle("p /GeV"); 
  hist_beta_vs_p_positive_layer1_comp[i]->GetYaxis()->SetTitle("#beta");
}

TH2F *hist_beta_vs_p_negative_layer1_comp[23];

for(Int_t i = 0; i < 23; i++){
  sprintf(name,"hist_beta_vs_p_negative_layer1_comp%01d", i+1);  
  sprintf(title," #beta vs p negative layer 1 component %01d", i+1);
  hist_beta_vs_p_negative_layer1_comp[i] = new TH2F(name, title, 600,0,Ebeam, 350, 0.0, 1.4);   
  hist_beta_vs_p_negative_layer1_comp[i]->GetXaxis()->SetTitle("p /GeV");
  hist_beta_vs_p_negative_layer1_comp[i]->GetYaxis()->SetTitle("#beta");
}

TH2F *hist_beta_vs_p_positive_layer2_comp[62]; 

for(Int_t i = 0; i < 62; i++){
  sprintf(name,"hist_beta_vs_p_positive_layer2_comp%01d", i+1);  
  sprintf(title," #beta vs p positive layer 2 component %01d", i+1);
  hist_beta_vs_p_positive_layer2_comp[i] = new TH2F(name, title, 600,0,Ebeam, 350, 0.0, 1.4);   
  hist_beta_vs_p_positive_layer2_comp[i]->GetXaxis()->SetTitle("p /GeV"); 
  hist_beta_vs_p_positive_layer2_comp[i]->GetYaxis()->SetTitle("#beta");
}

TH2F *hist_beta_vs_p_negative_layer2_comp[62];

for(Int_t i = 0; i < 62; i++){
  sprintf(name,"hist_beta_vs_p_negative_layer2_comp%01d", i+1);  
  sprintf(title," #beta vs p negative layer 2 component %01d", i+1);
  hist_beta_vs_p_negative_layer2_comp[i] = new TH2F(name, title, 600,0,Ebeam, 350, 0.0, 1.4);   
  hist_beta_vs_p_negative_layer2_comp[i]->GetXaxis()->SetTitle("p /GeV");
  hist_beta_vs_p_negative_layer2_comp[i]->GetYaxis()->SetTitle("#beta");
}


TH2F *hist_beta_vs_p_positive_layer3_comp[5]; 

for(Int_t i = 0; i < 5; i++){
  sprintf(name,"hist_beta_vs_p_positive_layer3_comp%01d", i+1);  
  sprintf(title," #beta vs p positive layer 3 component %01d", i+1);
  hist_beta_vs_p_positive_layer3_comp[i] = new TH2F(name, title, 600,0,Ebeam, 350, 0.0, 1.4);   
  hist_beta_vs_p_positive_layer3_comp[i]->GetXaxis()->SetTitle("p /GeV"); 
  hist_beta_vs_p_positive_layer3_comp[i]->GetYaxis()->SetTitle("#beta");
}

TH2F *hist_beta_vs_p_negative_layer3_comp[5];

for(Int_t i = 0; i < 5; i++){
  sprintf(name,"hist_beta_vs_p_negative_layer3_comp%01d", i+1);  
  sprintf(title," #beta vs p negative layer 3 component %01d", i+1);
  hist_beta_vs_p_negative_layer3_comp[i] = new TH2F(name, title, 600,0,Ebeam, 350, 0.0, 1.4);   
  hist_beta_vs_p_negative_layer3_comp[i]->GetXaxis()->SetTitle("p /GeV");
  hist_beta_vs_p_negative_layer3_comp[i]->GetYaxis()->SetTitle("#beta");
}


out->mkdir("FD_PID_hadron_time_per_paddle");				
out->cd ("FD_PID_hadron_time_per_paddle");

TH2F *hist_time_vs_p_positive_comp[62]; 
TH2F *hist_time_vs_p_negative_comp[62]; 

for(Int_t i = 0; i < 62; i++){
  sprintf(name,"hist_time_vs_p_positive_comp%01d", i+1);  sprintf(title,"flight time vs p positive component %01d", i+1);
  hist_time_vs_p_positive_comp[i] = new TH2F(name, title, 600,0,Ebeam, 500, 0, 50);   
  hist_time_vs_p_positive_comp[i]->GetXaxis()->SetTitle("p /GeV"); hist_time_vs_p_positive_comp[i]->GetYaxis()->SetTitle("(TOF time - start time) /ns");
  sprintf(name,"hist_time_vs_p_negative_comp%01d", i+1);  sprintf(title,"flight time vs p negative component %01d", i+1);
  hist_time_vs_p_negative_comp[i] = new TH2F(name, title, 600,0,Ebeam, 500, 0, 50);   
  hist_time_vs_p_negative_comp[i]->GetXaxis()->SetTitle("p /GeV"); hist_time_vs_p_negative_comp[i]->GetYaxis()->SetTitle("(TOF time - start time) /ns");
}


out->mkdir("FD_PID_hadron_tofmass_plots");				
out->cd ("FD_PID_hadron_tofmass_plots");

for(Int_t i = 0; i < 60; i++){ create_hist_tofmass_vs_p(i);}
for(Int_t i = 0; i < 60; i++){ create_hist_tofmass(i);}


out->mkdir("FD_PID_hadron_vertex_plots");				
out->cd ("FD_PID_hadron_vertex_plots");

for(Int_t i = 0; i < 60; i++){ create_hist_delta_vz(i);}


// CD hadrons

out->mkdir("CD_PID_hadron_beta_plots");				
out->cd ("CD_PID_hadron_beta_plots");

for(Int_t i = 0; i < 6; i++){ create_hist_CD_beta_vs_p(i); }
for(Int_t i = 10; i < 16; i++){ create_hist_CD_beta_vs_p(i); }  
for(Int_t i = 20; i < 26; i++){ create_hist_CD_beta_vs_p(i); }  
for(Int_t i = 30; i < 36; i++){ create_hist_CD_beta_vs_p(i); }  
for(Int_t i = 40; i < 46; i++){ create_hist_CD_beta_vs_p(i); }   
for(Int_t i = 50; i < 56; i++){ create_hist_CD_beta_vs_p(i); }  


out->mkdir("CD_PID_hadron_vertex_plots");				
out->cd ("CD_PID_hadron_vertex_plots");

for(Int_t i = 0; i < 6; i++){ create_hist_CD_delta_vz(i);}
for(Int_t i = 10; i < 16; i++){ create_hist_CD_delta_vz(i);}
for(Int_t i = 20; i < 26; i++){ create_hist_CD_delta_vz(i);}
for(Int_t i = 30; i < 36; i++){ create_hist_CD_delta_vz(i);}
for(Int_t i = 40; i < 46; i++){ create_hist_CD_delta_vz(i);}
for(Int_t i = 50; i < 56; i++){ create_hist_CD_delta_vz(i);}


out->mkdir("FD_PID_neutrals_plots");				
out->cd ("FD_PID_neutrals_plots");

for(Int_t i = 0; i <= 4; i++){
  create_hist_beta_vs_p_phot(i);
  create_hist_beta_phot(i);
  create_hist_EC_sampling_fraction_phot(i);
  create_hist_EC_PCAL_vs_EC_ECAL_phot(i);
  create_hist_EC_PCAL_hit_position_phot(i);
}


out->mkdir("FT_PID_plots");				
out->cd ("FT_PID_plots");

// elctrons
for(Int_t i = 0; i < 8; i++){
  create_hist_FT_FTCAL_energy_vs_radius(i);
  create_hist_FT_FTCAL_hit_position(i);
  create_hist_FT_FTTRK_hit_position(i);
  create_hist_FT_FTHODO_hit_position(i);
  create_hist_FT_beta(i);
}

// photons
for(Int_t i = 10; i < 16; i++){
  create_hist_FT_FTCAL_energy_vs_radius(i);
  create_hist_FT_FTCAL_hit_position(i);
  create_hist_FT_FTTRK_hit_position(i);
  create_hist_FT_FTHODO_hit_position(i);
  create_hist_FT_beta(i);
}


/// //////////////////////////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////////////////////////
/// FTOF monitoring
///

out->mkdir("FTOF_monitoring_layer1");				
out->cd ("FTOF_monitoring_layer1");

TH1F *hist_deltaT_proton_layer1[6][23]; 
TH1F *hist_deltaT_pip_layer1[6][23]; 
TH1F *hist_deltaT_pim_layer1[6][23]; 
TH1F *hist_deltaT_pion_combined_layer1[6][23]; 
TH1F *hist_deltaT_all_combined_layer1[6][23]; 

for(Int_t i = 0; i < 6; i++){
  for(Int_t j = 0; j < 23; j++){
    sprintf(name,"hist_deltat_proton_layer1_sec%01d_comp%01d", i+1, j+1);  
    sprintf(title," #Delta T proton layer 1 sector %01d component %01d", i+1, j+1);
    hist_deltaT_proton_layer1[i][j] = new TH1F(name, title, 400,-4,4);   
    hist_deltaT_proton_layer1[i][j]->GetXaxis()->SetTitle("#Delta T /ns"); 
    hist_deltaT_proton_layer1[i][j]->GetYaxis()->SetTitle("counts");
  }
}

for(Int_t i = 0; i < 6; i++){
  for(Int_t j = 0; j < 23; j++){
    sprintf(name,"hist_deltat_pip_layer1_sec%01d_comp%01d", i+1, j+1);  
    sprintf(title," #Delta T pip layer 1 sector %01d component %01d", i+1, j+1);
    hist_deltaT_pip_layer1[i][j] = new TH1F(name, title, 400,-4,4);   
    hist_deltaT_pip_layer1[i][j]->GetXaxis()->SetTitle("#Delta T /ns"); 
    hist_deltaT_pip_layer1[i][j]->GetYaxis()->SetTitle("counts");
  }
}

for(Int_t i = 0; i < 6; i++){
  for(Int_t j = 0; j < 23; j++){
    sprintf(name,"hist_deltat_pim_layer1_sec%01d_comp%01d", i+1, j+1);  
    sprintf(title," #Delta T pim layer 1 sector %01d component %01d", i+1, j+1);
    hist_deltaT_pim_layer1[i][j] = new TH1F(name, title, 400,-4,4);   
    hist_deltaT_pim_layer1[i][j]->GetXaxis()->SetTitle("#Delta T /ns"); 
    hist_deltaT_pim_layer1[i][j]->GetYaxis()->SetTitle("counts");
  }
}

for(Int_t i = 0; i < 6; i++){
  for(Int_t j = 0; j < 23; j++){
    sprintf(name,"hist_deltat_pion_combined_layer1_sec%01d_comp%01d", i+1, j+1);  
    sprintf(title," #Delta T pion combined layer 1 sector %01d component %01d", i+1, j+1);
    hist_deltaT_pion_combined_layer1[i][j] = new TH1F(name, title, 400,-4,4);   
    hist_deltaT_pion_combined_layer1[i][j]->GetXaxis()->SetTitle("#Delta T /ns"); 
    hist_deltaT_pion_combined_layer1[i][j]->GetYaxis()->SetTitle("counts");
  }
}

for(Int_t i = 0; i < 6; i++){
  for(Int_t j = 0; j < 23; j++){
    sprintf(name,"hist_deltat_all_combined_layer1_sec%01d_comp%01d", i+1, j+1);  
    sprintf(title," #Delta T all combined layer 1 sector %01d component %01d", i+1, j+1);
    hist_deltaT_all_combined_layer1[i][j] = new TH1F(name, title, 400,-4,4);   
    hist_deltaT_all_combined_layer1[i][j]->GetXaxis()->SetTitle("#Delta T /ns"); 
    hist_deltaT_all_combined_layer1[i][j]->GetYaxis()->SetTitle("counts");
  }
}



out->mkdir("FTOF_monitoring_layer2");				
out->cd ("FTOF_monitoring_layer2"); 

TH1F *hist_deltaT_proton_layer2[6][62]; 
TH1F *hist_deltaT_pip_layer2[6][62]; 
TH1F *hist_deltaT_pim_layer2[6][62]; 
TH1F *hist_deltaT_pion_combined_layer2[6][62]; 
TH1F *hist_deltaT_all_combined_layer2[6][62]; 


for(Int_t i = 0; i < 6; i++){
  for(Int_t j = 0; j < 62; j++){
    sprintf(name,"hist_deltat_proton_layer2_sec%01d_comp%01d", i+1, j+1);  
    sprintf(title," #Delta T proton layer 2 sector %01d component %01d", i+1, j+1);
    hist_deltaT_proton_layer2[i][j] = new TH1F(name, title, 400,-4,4);   
    hist_deltaT_proton_layer2[i][j]->GetXaxis()->SetTitle("#Delta T /ns"); 
    hist_deltaT_proton_layer2[i][j]->GetYaxis()->SetTitle("counts");
  }
}

for(Int_t i = 0; i < 6; i++){
  for(Int_t j = 0; j < 62; j++){
    sprintf(name,"hist_deltat_pip_layer2_sec%01d_comp%01d", i+1, j+1);  
    sprintf(title," #Delta T pip layer 2 sector %01d component %01d", i+1, j+1);
    hist_deltaT_pip_layer2[i][j] = new TH1F(name, title, 400,-4,4);   
    hist_deltaT_pip_layer2[i][j]->GetXaxis()->SetTitle("#Delta T /ns"); 
    hist_deltaT_pip_layer2[i][j]->GetYaxis()->SetTitle("counts");
  }
}

for(Int_t i = 0; i < 6; i++){
  for(Int_t j = 0; j < 62; j++){
    sprintf(name,"hist_deltat_pim_layer2_sec%01d_comp%01d", i+1, j+1);  
    sprintf(title," #Delta T pim layer 2 sector %01d component %01d", i+1, j+1);
    hist_deltaT_pim_layer2[i][j] = new TH1F(name, title, 400,-4,4);   
    hist_deltaT_pim_layer2[i][j]->GetXaxis()->SetTitle("#Delta T /ns"); 
    hist_deltaT_pim_layer2[i][j]->GetYaxis()->SetTitle("counts");
  }
}

for(Int_t i = 0; i < 6; i++){
  for(Int_t j = 0; j < 62; j++){
    sprintf(name,"hist_deltat_pion_combined_layer2_sec%01d_comp%01d", i+1, j+1);  
    sprintf(title," #Delta T pion combined layer 2 sector %01d component %01d", i+1, j+1);
    hist_deltaT_pion_combined_layer2[i][j] = new TH1F(name, title, 400,-4,4);   
    hist_deltaT_pion_combined_layer2[i][j]->GetXaxis()->SetTitle("#Delta T /ns"); 
    hist_deltaT_pion_combined_layer2[i][j]->GetYaxis()->SetTitle("counts");
  }
}

for(Int_t i = 0; i < 6; i++){
  for(Int_t j = 0; j < 62; j++){
    sprintf(name,"hist_deltat_all_combined_layer2_sec%01d_comp%01d", i+1, j+1);  
    sprintf(title," #Delta T all combined layer 2 sector %01d component %01d", i+1, j+1);
    hist_deltaT_all_combined_layer2[i][j] = new TH1F(name, title, 400,-4,4);   
    hist_deltaT_all_combined_layer2[i][j]->GetXaxis()->SetTitle("#Delta T /ns"); 
    hist_deltaT_all_combined_layer2[i][j]->GetYaxis()->SetTitle("counts");
  }
}



out->mkdir("FTOF_monitoring_layer3");				
out->cd ("FTOF_monitoring_layer3");


TH1F *hist_deltaT_proton_layer3[6][5]; 
TH1F *hist_deltaT_pip_layer3[6][5]; 
TH1F *hist_deltaT_pim_layer3[6][5]; 
TH1F *hist_deltaT_pion_combined_layer3[6][5]; 
TH1F *hist_deltaT_all_combined_layer3[6][5]; 


for(Int_t i = 0; i < 6; i++){
  for(Int_t j = 0; j < 5; j++){
    sprintf(name,"hist_deltat_proton_layer3_sec%01d_comp%01d", i+1, j+1);  
    sprintf(title," #Delta T proton layer 3 sector %01d component %01d", i+1, j+1);
    hist_deltaT_proton_layer3[i][j] = new TH1F(name, title, 1000,-10,10);   
    hist_deltaT_proton_layer3[i][j]->GetXaxis()->SetTitle("#Delta T /ns"); 
    hist_deltaT_proton_layer3[i][j]->GetYaxis()->SetTitle("counts");
  }
}

for(Int_t i = 0; i < 6; i++){
  for(Int_t j = 0; j < 5; j++){
    sprintf(name,"hist_deltat_pip_layer3_sec%01d_comp%01d", i+1, j+1);  
    sprintf(title," #Delta T pip layer 3 sector %01d component %01d", i+1, j+1);
    hist_deltaT_pip_layer3[i][j] = new TH1F(name, title, 1000,-10,10);   
    hist_deltaT_pip_layer3[i][j]->GetXaxis()->SetTitle("#Delta T /ns"); 
    hist_deltaT_pip_layer3[i][j]->GetYaxis()->SetTitle("counts");
  }
}

for(Int_t i = 0; i < 6; i++){
  for(Int_t j = 0; j < 5; j++){
    sprintf(name,"hist_deltat_pim_layer3_sec%01d_comp%01d", i+1, j+1);  
    sprintf(title," #Delta T pim layer 3 sector %01d component %01d", i+1, j+1);
    hist_deltaT_pim_layer3[i][j] = new TH1F(name, title, 1000,-10,10);   
    hist_deltaT_pim_layer3[i][j]->GetXaxis()->SetTitle("#Delta T /ns"); 
    hist_deltaT_pim_layer3[i][j]->GetYaxis()->SetTitle("counts");
  }
}

for(Int_t i = 0; i < 6; i++){
  for(Int_t j = 0; j < 5; j++){
    sprintf(name,"hist_deltat_pion_combined_layer3_sec%01d_comp%01d", i+1, j+1);  
    sprintf(title," #Delta T pion combined layer 3 sector %01d component %01d", i+1, j+1);
    hist_deltaT_pion_combined_layer3[i][j] = new TH1F(name, title, 1000,-10,10);   
    hist_deltaT_pion_combined_layer3[i][j]->GetXaxis()->SetTitle("#Delta T /ns"); 
    hist_deltaT_pion_combined_layer3[i][j]->GetYaxis()->SetTitle("counts");
  }
}

for(Int_t i = 0; i < 6; i++){
  for(Int_t j = 0; j < 5; j++){
    sprintf(name,"hist_deltat_all_combined_layer3_sec%01d_comp%01d", i+1, j+1);  
    sprintf(title," #Delta T all combined layer 3 sector %01d component %01d", i+1, j+1);
    hist_deltaT_all_combined_layer3[i][j] = new TH1F(name, title, 1000,-10,10);   
    hist_deltaT_all_combined_layer3[i][j]->GetXaxis()->SetTitle("#Delta T /ns"); 
    hist_deltaT_all_combined_layer3[i][j]->GetYaxis()->SetTitle("counts");
  }
}


out->mkdir("CTOF_monitoring");				
out->cd ("CTOF_monitoring");

TH1F *hist_deltaT_proton_CTOF[48]; 
TH1F *hist_deltaT_pip_CTOF[48]; 
TH1F *hist_deltaT_pim_CTOF[48]; 
TH1F *hist_deltaT_pion_combined_CTOF[48]; 
TH1F *hist_deltaT_all_combined_CTOF[48]; 


for(Int_t j = 0; j < 48; j++){
  sprintf(name,"hist_deltat_proton_CTOF_comp%01d", j+1);  
  sprintf(title," #Delta T proton CTOF component %01d", j+1);
  hist_deltaT_proton_CTOF[j] = new TH1F(name, title, 400,-4,4);   
  hist_deltaT_proton_CTOF[j]->GetXaxis()->SetTitle("#Delta T /ns"); 
  hist_deltaT_proton_CTOF[j]->GetYaxis()->SetTitle("counts");
}

for(Int_t j = 0; j < 48; j++){
  sprintf(name,"hist_deltat_pip_CTOF_comp%01d", j+1);  
  sprintf(title," #Delta T pip CTOF component %01d", j+1);
  hist_deltaT_pip_CTOF[j] = new TH1F(name, title, 400,-4,4);   
  hist_deltaT_pip_CTOF[j]->GetXaxis()->SetTitle("#Delta T /ns"); 
  hist_deltaT_pip_CTOF[j]->GetYaxis()->SetTitle("counts");
}

for(Int_t j = 0; j < 48; j++){
  sprintf(name,"hist_deltat_pim_CTOF_comp%01d",j+1);  
  sprintf(title," #Delta T pim CTOF component %01d", j+1);
  hist_deltaT_pim_CTOF[j] = new TH1F(name, title, 400,-4,4);   
  hist_deltaT_pim_CTOF[j]->GetXaxis()->SetTitle("#Delta T /ns"); 
  hist_deltaT_pim_CTOF[j]->GetYaxis()->SetTitle("counts");
}

for(Int_t j = 0; j < 48; j++){
  sprintf(name,"hist_deltat_pion_combined_CTOF_comp%01d", j+1);  
  sprintf(title," #Delta T pion combined CTOF component %01d", j+1);
  hist_deltaT_pion_combined_CTOF[j] = new TH1F(name, title, 400,-4,4);   
  hist_deltaT_pion_combined_CTOF[j]->GetXaxis()->SetTitle("#Delta T /ns"); 
  hist_deltaT_pion_combined_CTOF[j]->GetYaxis()->SetTitle("counts");
}

for(Int_t j = 0; j < 48; j++){
  sprintf(name,"hist_deltat_all_combined_CTOF_comp%01d", j+1);  
  sprintf(title," #Delta T all combined CTOF component %01d", j+1);
  hist_deltaT_all_combined_CTOF[j] = new TH1F(name, title, 400,-4,4);   
  hist_deltaT_all_combined_CTOF[j]->GetXaxis()->SetTitle("#Delta T /ns"); 
  hist_deltaT_all_combined_CTOF[j]->GetYaxis()->SetTitle("counts");
}




out->mkdir("TOFmass");				
out->cd ("TOFmass");

TH1F *hist_TOFmass_FTOF_positive; 
TH1F *hist_TOFmass_FTOF_negative; 
TH1F *hist_TOFmass_CTOF_positive; 
TH1F *hist_TOFmass_CTOF_negative; 

hist_TOFmass_FTOF_positive = new TH1F("hist_TOFmass_FTOF_positive", "TOF mass from the forward detector for hadron with positive charge", 1200,-0.4,2.0);   
hist_TOFmass_FTOF_positive->GetXaxis()->SetTitle("m_{TOF}^{2} /GeV^{2}"); 
hist_TOFmass_FTOF_positive->GetYaxis()->SetTitle("counts");
hist_TOFmass_FTOF_negative = new TH1F("hist_TOFmass_FTOF_negative", "TOF mass from the forward detector for hadron with negative charge", 1200,-0.4,2.0);   
hist_TOFmass_FTOF_negative->GetXaxis()->SetTitle("m_{TOF}^{2} /GeV^{2}"); 
hist_TOFmass_FTOF_negative->GetYaxis()->SetTitle("counts");

hist_TOFmass_CTOF_positive = new TH1F("hist_TOFmass_CTOF_positive", "TOF mass from the central detector for hadron with positive charge", 1200,-0.4,2.0);   
hist_TOFmass_CTOF_positive->GetXaxis()->SetTitle("m_{TOF}^{2} /GeV^{2}"); 
hist_TOFmass_CTOF_positive->GetYaxis()->SetTitle("counts");
hist_TOFmass_CTOF_negative = new TH1F("hist_TOFmass_CTOF_negative", "TOF mass from the central detector for hadron with negative charge", 1200,-0.4,2.0);   
hist_TOFmass_CTOF_negative->GetXaxis()->SetTitle("m_{TOF}^{2} /GeV^{2}"); 
hist_TOFmass_CTOF_negative->GetYaxis()->SetTitle("counts");



/// //////////////////////////////////////////////////////////////////////////////
/// Fiducial cut plots

out->mkdir("EC_fiducial cut plots");
out->cd ("EC_fiducial cut plots");

  TH2F *hist_electron_sampfrac_vs_u_coord_sec1;
  TH2F *hist_electron_sampfrac_vs_u_coord_sec2;
  TH2F *hist_electron_sampfrac_vs_u_coord_sec3;
  TH2F *hist_electron_sampfrac_vs_u_coord_sec4;
  TH2F *hist_electron_sampfrac_vs_u_coord_sec5;
  TH2F *hist_electron_sampfrac_vs_u_coord_sec6;
  TH2F *hist_electron_sampfrac_vs_u_coord;

  TH2F *hist_electron_sampfrac_vs_v_coord_sec1;
  TH2F *hist_electron_sampfrac_vs_v_coord_sec2;
  TH2F *hist_electron_sampfrac_vs_v_coord_sec3;
  TH2F *hist_electron_sampfrac_vs_v_coord_sec4;
  TH2F *hist_electron_sampfrac_vs_v_coord_sec5;
  TH2F *hist_electron_sampfrac_vs_v_coord_sec6;
  TH2F *hist_electron_sampfrac_vs_v_coord;

  TH2F *hist_electron_sampfrac_vs_w_coord_sec1;
  TH2F *hist_electron_sampfrac_vs_w_coord_sec2;
  TH2F *hist_electron_sampfrac_vs_w_coord_sec3;
  TH2F *hist_electron_sampfrac_vs_w_coord_sec4;
  TH2F *hist_electron_sampfrac_vs_w_coord_sec5;
  TH2F *hist_electron_sampfrac_vs_w_coord_sec6;
  TH2F *hist_electron_sampfrac_vs_w_coord;

  TH3F *hist_electron_sampfrac_vs_u_coord_vs_p_sec1;
  TH3F *hist_electron_sampfrac_vs_u_coord_vs_p_sec2;
  TH3F *hist_electron_sampfrac_vs_u_coord_vs_p_sec3;
  TH3F *hist_electron_sampfrac_vs_u_coord_vs_p_sec4;
  TH3F *hist_electron_sampfrac_vs_u_coord_vs_p_sec5;
  TH3F *hist_electron_sampfrac_vs_u_coord_vs_p_sec6;

  TH3F *hist_electron_sampfrac_vs_v_coord_vs_p_sec1;
  TH3F *hist_electron_sampfrac_vs_v_coord_vs_p_sec2;
  TH3F *hist_electron_sampfrac_vs_v_coord_vs_p_sec3;
  TH3F *hist_electron_sampfrac_vs_v_coord_vs_p_sec4;
  TH3F *hist_electron_sampfrac_vs_v_coord_vs_p_sec5;
  TH3F *hist_electron_sampfrac_vs_v_coord_vs_p_sec6;

  TH3F *hist_electron_sampfrac_vs_w_coord_vs_p_sec1;
  TH3F *hist_electron_sampfrac_vs_w_coord_vs_p_sec2;
  TH3F *hist_electron_sampfrac_vs_w_coord_vs_p_sec3;
  TH3F *hist_electron_sampfrac_vs_w_coord_vs_p_sec4;
  TH3F *hist_electron_sampfrac_vs_w_coord_vs_p_sec5;
  TH3F *hist_electron_sampfrac_vs_w_coord_vs_p_sec6;

  hist_electron_sampfrac_vs_u_coord_vs_p_sec1 = new TH3F("hist_electron_sampfrac_vs_u_coord_vs_p_sec1", "electron E/p vs u coordinate vs p sector 1", 450, 0, 450, 50, 0, 0.5, 100, 0.0 , Ebeam+1);   
  hist_electron_sampfrac_vs_u_coord_vs_p_sec1->GetXaxis()->SetTitle("u /cm");
  hist_electron_sampfrac_vs_u_coord_vs_p_sec1->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_u_coord_vs_p_sec1->GetZaxis()->SetTitle("p /GeV");
  hist_electron_sampfrac_vs_u_coord_vs_p_sec2 = new TH3F("hist_electron_sampfrac_vs_u_coord_vs_p_sec2", "electron E/p vs u coordinate vs p sector 2", 450, 0, 450, 50, 0, 0.5, 100, 0.0 , Ebeam+1);   
  hist_electron_sampfrac_vs_u_coord_vs_p_sec2->GetXaxis()->SetTitle("u /cm");
  hist_electron_sampfrac_vs_u_coord_vs_p_sec2->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_u_coord_vs_p_sec2->GetZaxis()->SetTitle("p /GeV");
  hist_electron_sampfrac_vs_u_coord_vs_p_sec3 = new TH3F("hist_electron_sampfrac_vs_u_coord_vs_p_sec3", "electron E/p vs u coordinate vs p sector 3", 450, 0, 450, 50, 0, 0.5, 100, 0.0 , Ebeam+1);   
  hist_electron_sampfrac_vs_u_coord_vs_p_sec3->GetXaxis()->SetTitle("u /cm");
  hist_electron_sampfrac_vs_u_coord_vs_p_sec3->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_u_coord_vs_p_sec3->GetZaxis()->SetTitle("p /GeV");
  hist_electron_sampfrac_vs_u_coord_vs_p_sec4 = new TH3F("hist_electron_sampfrac_vs_u_coord_vs_p_sec4", "electron E/p vs u coordinate vs p sector 4", 450, 0, 450, 50, 0, 0.5, 100, 0.0 , Ebeam+1);   
  hist_electron_sampfrac_vs_u_coord_vs_p_sec4->GetXaxis()->SetTitle("u /cm");
  hist_electron_sampfrac_vs_u_coord_vs_p_sec4->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_u_coord_vs_p_sec4->GetZaxis()->SetTitle("p /GeV");
  hist_electron_sampfrac_vs_u_coord_vs_p_sec5 = new TH3F("hist_electron_sampfrac_vs_u_coord_vs_p_sec5", "electron E/p vs u coordinate vs p sector 5", 450, 0, 450, 50, 0, 0.5, 100, 0.0 , Ebeam+1);   
  hist_electron_sampfrac_vs_u_coord_vs_p_sec5->GetXaxis()->SetTitle("u /cm");
  hist_electron_sampfrac_vs_u_coord_vs_p_sec5->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_u_coord_vs_p_sec5->GetZaxis()->SetTitle("p /GeV");
  hist_electron_sampfrac_vs_u_coord_vs_p_sec6 = new TH3F("hist_electron_sampfrac_vs_u_coord_vs_p_sec6", "electron E/p vs u coordinate vs p sector 6", 450, 0, 450, 50, 0, 0.5, 100, 0.0 , Ebeam+1);   
  hist_electron_sampfrac_vs_u_coord_vs_p_sec6->GetXaxis()->SetTitle("u /cm");
  hist_electron_sampfrac_vs_u_coord_vs_p_sec6->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_u_coord_vs_p_sec6->GetZaxis()->SetTitle("p /GeV");

  hist_electron_sampfrac_vs_v_coord_vs_p_sec1 = new TH3F("hist_electron_sampfrac_vs_v_coord_vs_p_sec1", "electron E/p vs v coordinate vs p sector 1", 450, 0, 450, 50, 0, 0.5, 100, 0.0 , Ebeam+1);   
  hist_electron_sampfrac_vs_v_coord_vs_p_sec1->GetXaxis()->SetTitle("v /cm");
  hist_electron_sampfrac_vs_v_coord_vs_p_sec1->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_v_coord_vs_p_sec1->GetZaxis()->SetTitle("p /GeV");
  hist_electron_sampfrac_vs_v_coord_vs_p_sec2 = new TH3F("hist_electron_sampfrac_vs_v_coord_vs_p_sec2", "electron E/p vs v coordinate vs p sector 2", 450, 0, 450, 50, 0, 0.5, 100, 0.0 , Ebeam+1);   
  hist_electron_sampfrac_vs_v_coord_vs_p_sec2->GetXaxis()->SetTitle("v /cm");
  hist_electron_sampfrac_vs_v_coord_vs_p_sec2->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_v_coord_vs_p_sec2->GetZaxis()->SetTitle("p /GeV");
  hist_electron_sampfrac_vs_v_coord_vs_p_sec3 = new TH3F("hist_electron_sampfrac_vs_v_coord_vs_p_sec3", "electron E/p vs v coordinate vs p sector 3", 450, 0, 450, 50, 0, 0.5, 100, 0.0 , Ebeam+1);   
  hist_electron_sampfrac_vs_v_coord_vs_p_sec3->GetXaxis()->SetTitle("v /cm");
  hist_electron_sampfrac_vs_v_coord_vs_p_sec3->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_v_coord_vs_p_sec3->GetZaxis()->SetTitle("p /GeV");
  hist_electron_sampfrac_vs_v_coord_vs_p_sec4 = new TH3F("hist_electron_sampfrac_vs_v_coord_vs_p_sec4", "electron E/p vs v coordinate vs p sector 4", 450, 0, 450, 50, 0, 0.5, 100, 0.0 , Ebeam+1);   
  hist_electron_sampfrac_vs_v_coord_vs_p_sec4->GetXaxis()->SetTitle("v /cm");
  hist_electron_sampfrac_vs_v_coord_vs_p_sec4->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_v_coord_vs_p_sec4->GetZaxis()->SetTitle("p /GeV");
  hist_electron_sampfrac_vs_v_coord_vs_p_sec5 = new TH3F("hist_electron_sampfrac_vs_v_coord_vs_p_sec5", "electron E/p vs v coordinate vs p sector 5", 450, 0, 450, 50, 0, 0.5, 100, 0.0 , Ebeam+1);   
  hist_electron_sampfrac_vs_v_coord_vs_p_sec5->GetXaxis()->SetTitle("v /cm");
  hist_electron_sampfrac_vs_v_coord_vs_p_sec5->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_v_coord_vs_p_sec5->GetZaxis()->SetTitle("p /GeV");
  hist_electron_sampfrac_vs_v_coord_vs_p_sec6 = new TH3F("hist_electron_sampfrac_vs_v_coord_vs_p_sec6", "electron E/p vs v coordinate vs p sector 6", 450, 0, 450, 50, 0, 0.5, 100, 0.0 , Ebeam+1);   
  hist_electron_sampfrac_vs_v_coord_vs_p_sec6->GetXaxis()->SetTitle("v /cm");
  hist_electron_sampfrac_vs_v_coord_vs_p_sec6->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_v_coord_vs_p_sec6->GetZaxis()->SetTitle("p /GeV");

  hist_electron_sampfrac_vs_w_coord_vs_p_sec1 = new TH3F("hist_electron_sampfrac_vs_w_coord_vs_p_sec1", "electron E/p vs w coordinate vs p sector 1", 450, 0, 450, 50, 0, 0.5, 100, 0.0 , Ebeam+1);   
  hist_electron_sampfrac_vs_w_coord_vs_p_sec1->GetXaxis()->SetTitle("w /cm");
  hist_electron_sampfrac_vs_w_coord_vs_p_sec1->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_w_coord_vs_p_sec1->GetZaxis()->SetTitle("p /GeV");
  hist_electron_sampfrac_vs_w_coord_vs_p_sec2 = new TH3F("hist_electron_sampfrac_vs_w_coord_vs_p_sec2", "electron E/p vs w coordinate vs p sector 2", 450, 0, 450, 50, 0, 0.5, 100, 0.0 , Ebeam+1);   
  hist_electron_sampfrac_vs_w_coord_vs_p_sec2->GetXaxis()->SetTitle("w /cm");
  hist_electron_sampfrac_vs_w_coord_vs_p_sec2->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_w_coord_vs_p_sec2->GetZaxis()->SetTitle("p /GeV");
  hist_electron_sampfrac_vs_w_coord_vs_p_sec3 = new TH3F("hist_electron_sampfrac_vs_w_coord_vs_p_sec3", "electron E/p vs w coordinate vs p sector 3", 450, 0, 450, 50, 0, 0.5, 100, 0.0 , Ebeam+1);   
  hist_electron_sampfrac_vs_w_coord_vs_p_sec3->GetXaxis()->SetTitle("w /cm");
  hist_electron_sampfrac_vs_w_coord_vs_p_sec3->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_w_coord_vs_p_sec3->GetZaxis()->SetTitle("p /GeV");
  hist_electron_sampfrac_vs_w_coord_vs_p_sec4 = new TH3F("hist_electron_sampfrac_vs_w_coord_vs_p_sec4", "electron E/p vs w coordinate vs p sector 4", 450, 0, 450, 50, 0, 0.5, 100, 0.0 , Ebeam+1);   
  hist_electron_sampfrac_vs_w_coord_vs_p_sec4->GetXaxis()->SetTitle("w /cm");
  hist_electron_sampfrac_vs_w_coord_vs_p_sec4->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_w_coord_vs_p_sec4->GetZaxis()->SetTitle("p /GeV");
  hist_electron_sampfrac_vs_w_coord_vs_p_sec5 = new TH3F("hist_electron_sampfrac_vs_w_coord_vs_p_sec5", "electron E/p vs w coordinate vs p sector 5", 450, 0, 450, 50, 0, 0.5, 100, 0.0 , Ebeam+1);   
  hist_electron_sampfrac_vs_w_coord_vs_p_sec5->GetXaxis()->SetTitle("w /cm");
  hist_electron_sampfrac_vs_w_coord_vs_p_sec5->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_w_coord_vs_p_sec5->GetZaxis()->SetTitle("p /GeV");
  hist_electron_sampfrac_vs_w_coord_vs_p_sec6 = new TH3F("hist_electron_sampfrac_vs_w_coord_vs_p_sec6", "electron E/p vs w coordinate vs p sector 6", 450, 0, 450, 50, 0, 0.5, 100, 0.0 , Ebeam+1);   
  hist_electron_sampfrac_vs_w_coord_vs_p_sec6->GetXaxis()->SetTitle("w /cm");
  hist_electron_sampfrac_vs_w_coord_vs_p_sec6->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_w_coord_vs_p_sec6->GetZaxis()->SetTitle("p /GeV");


  hist_electron_sampfrac_vs_u_coord_sec1 = new TH2F("hist_electron_sampfrac_vs_u_coord_sec1", "electron E/p vs u coordinate sector 1", 900, 0, 450, 50, 0, 0.5);   
  hist_electron_sampfrac_vs_u_coord_sec1->GetXaxis()->SetTitle("u /cm");
  hist_electron_sampfrac_vs_u_coord_sec1->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_u_coord_sec2 = new TH2F("hist_electron_sampfrac_vs_u_coord_sec2", "electron E/p vs u coordinate sector 2", 900, 0, 450, 50, 0, 0.5);   
  hist_electron_sampfrac_vs_u_coord_sec2->GetXaxis()->SetTitle("u /cm");
  hist_electron_sampfrac_vs_u_coord_sec2->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_u_coord_sec3 = new TH2F("hist_electron_sampfrac_vs_u_coord_sec3", "electron E/p vs u coordinate sector 3", 900, 0, 450, 50, 0, 0.5);   
  hist_electron_sampfrac_vs_u_coord_sec3->GetXaxis()->SetTitle("u /cm");
  hist_electron_sampfrac_vs_u_coord_sec3->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_u_coord_sec4 = new TH2F("hist_electron_sampfrac_vs_u_coord_sec4", "electron E/p vs u coordinate sector 4", 900, 0, 450, 50, 0, 0.5);   
  hist_electron_sampfrac_vs_u_coord_sec4->GetXaxis()->SetTitle("u /cm");
  hist_electron_sampfrac_vs_u_coord_sec4->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_u_coord_sec5 = new TH2F("hist_electron_sampfrac_vs_u_coord_sec5", "electron E/p vs u coordinate sector 5", 900, 0, 450, 50, 0, 0.5);   
  hist_electron_sampfrac_vs_u_coord_sec5->GetXaxis()->SetTitle("u /cm");
  hist_electron_sampfrac_vs_u_coord_sec5->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_u_coord_sec6 = new TH2F("hist_electron_sampfrac_vs_u_coord_sec6", "electron E/p vs u coordinate sector 6", 900, 0, 450, 50, 0, 0.5);   
  hist_electron_sampfrac_vs_u_coord_sec6->GetXaxis()->SetTitle("u /cm");
  hist_electron_sampfrac_vs_u_coord_sec6->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_u_coord = new TH2F("hist_electron_sampfrac_vs_u_coord", "electron E/p vs u coordinate", 900, 0, 450, 50, 0, 0.5);   
  hist_electron_sampfrac_vs_u_coord->GetXaxis()->SetTitle("u /cm");
  hist_electron_sampfrac_vs_u_coord->GetYaxis()->SetTitle("E/p");

  hist_electron_sampfrac_vs_v_coord_sec1 = new TH2F("hist_electron_sampfrac_vs_v_coord_sec1", "electron E/p vs v coordinate sector 1", 900, 0, 450, 50, 0, 0.5);   
  hist_electron_sampfrac_vs_v_coord_sec1->GetXaxis()->SetTitle("v /cm");
  hist_electron_sampfrac_vs_v_coord_sec1->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_v_coord_sec2 = new TH2F("hist_electron_sampfrac_vs_v_coord_sec2", "electron E/p vs v coordinate sector 2", 900, 0, 450, 50, 0, 0.5);   
  hist_electron_sampfrac_vs_v_coord_sec2->GetXaxis()->SetTitle("v /cm");
  hist_electron_sampfrac_vs_v_coord_sec2->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_v_coord_sec3 = new TH2F("hist_electron_sampfrac_vs_v_coord_sec3", "electron E/p vs v coordinate sector 3", 900, 0, 450, 50, 0, 0.5);   
  hist_electron_sampfrac_vs_v_coord_sec3->GetXaxis()->SetTitle("v /cm");
  hist_electron_sampfrac_vs_v_coord_sec3->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_v_coord_sec4 = new TH2F("hist_electron_sampfrac_vs_v_coord_sec4", "electron E/p vs v coordinate sector 4", 900, 0, 450, 50, 0, 0.5);   
  hist_electron_sampfrac_vs_v_coord_sec4->GetXaxis()->SetTitle("v /cm");
  hist_electron_sampfrac_vs_v_coord_sec4->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_v_coord_sec5 = new TH2F("hist_electron_sampfrac_vs_v_coord_sec5", "electron E/p vs v coordinate sector 5", 900, 0, 450, 50, 0, 0.5);   
  hist_electron_sampfrac_vs_v_coord_sec5->GetXaxis()->SetTitle("v /cm");
  hist_electron_sampfrac_vs_v_coord_sec5->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_v_coord_sec6 = new TH2F("hist_electron_sampfrac_vs_v_coord_sec6", "electron E/p vs v coordinate sector 6", 900, 0, 450, 50, 0, 0.5);   
  hist_electron_sampfrac_vs_v_coord_sec6->GetXaxis()->SetTitle("v /cm");
  hist_electron_sampfrac_vs_v_coord_sec6->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_v_coord = new TH2F("hist_electron_sampfrac_vs_v_coord", "electron E/p vs v coordinate", 900, 0, 450, 50, 0, 0.5);   
  hist_electron_sampfrac_vs_v_coord->GetXaxis()->SetTitle("v /cm");
  hist_electron_sampfrac_vs_v_coord->GetYaxis()->SetTitle("E/p");

  hist_electron_sampfrac_vs_w_coord_sec1 = new TH2F("hist_electron_sampfrac_vs_w_coord_sec1", "electron E/p vs w coordinate sector 1", 900, 0, 450, 50, 0, 0.5);   
  hist_electron_sampfrac_vs_w_coord_sec1->GetXaxis()->SetTitle("w /cm");
  hist_electron_sampfrac_vs_w_coord_sec1->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_w_coord_sec2 = new TH2F("hist_electron_sampfrac_vs_w_coord_sec2", "electron E/p vs w coordinate sector 2", 900, 0, 450, 50, 0, 0.5);   
  hist_electron_sampfrac_vs_w_coord_sec2->GetXaxis()->SetTitle("w /cm");
  hist_electron_sampfrac_vs_w_coord_sec2->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_w_coord_sec3 = new TH2F("hist_electron_sampfrac_vs_w_coord_sec3", "electron E/p vs w coordinate sector 3", 900, 0, 450, 50, 0, 0.5);   
  hist_electron_sampfrac_vs_w_coord_sec3->GetXaxis()->SetTitle("w /cm");
  hist_electron_sampfrac_vs_w_coord_sec3->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_w_coord_sec4 = new TH2F("hist_electron_sampfrac_vs_w_coord_sec4", "electron E/p vs w coordinate sector 4", 900, 0, 450, 50, 0, 0.5);   
  hist_electron_sampfrac_vs_w_coord_sec4->GetXaxis()->SetTitle("w /cm");
  hist_electron_sampfrac_vs_w_coord_sec4->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_w_coord_sec5 = new TH2F("hist_electron_sampfrac_vs_w_coord_sec5", "electron E/p vs w coordinate sector 5", 900, 0, 450, 50, 0, 0.5);   
  hist_electron_sampfrac_vs_w_coord_sec5->GetXaxis()->SetTitle("w /cm");
  hist_electron_sampfrac_vs_w_coord_sec5->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_w_coord_sec6 = new TH2F("hist_electron_sampfrac_vs_w_coord_sec6", "electron E/p vs w coordinate sector 6", 900, 0, 450, 50, 0, 0.5);   
  hist_electron_sampfrac_vs_w_coord_sec6->GetXaxis()->SetTitle("w /cm");
  hist_electron_sampfrac_vs_w_coord_sec6->GetYaxis()->SetTitle("E/p");
  hist_electron_sampfrac_vs_w_coord = new TH2F("hist_electron_sampfrac_vs_w_coord", "electron E/p vs w coordinate", 900, 0, 450, 50, 0, 0.5);   
  hist_electron_sampfrac_vs_w_coord->GetXaxis()->SetTitle("w /cm");
  hist_electron_sampfrac_vs_w_coord->GetYaxis()->SetTitle("E/p");


/// //////////////////////////////////////////////////////////////////////////////
/// Fiducial cut plots

out->mkdir("HTCC_Nphe");				
out->cd ("HTCC_Nphe");

  TH2F *hist_HTCC_Nphe_vs_momentum;
  TH2F *hist_HTCC_Nphe_vs_momentum_sec1;
  TH2F *hist_HTCC_Nphe_vs_momentum_sec2;
  TH2F *hist_HTCC_Nphe_vs_momentum_sec3;
  TH2F *hist_HTCC_Nphe_vs_momentum_sec4;
  TH2F *hist_HTCC_Nphe_vs_momentum_sec5;
  TH2F *hist_HTCC_Nphe_vs_momentum_sec6;

  TH2F *hist_HTCC_Nphe_vs_momentum_prot;
  TH2F *hist_HTCC_Nphe_vs_momentum_prot_sec1;
  TH2F *hist_HTCC_Nphe_vs_momentum_prot_sec2;
  TH2F *hist_HTCC_Nphe_vs_momentum_prot_sec3;
  TH2F *hist_HTCC_Nphe_vs_momentum_prot_sec4;
  TH2F *hist_HTCC_Nphe_vs_momentum_prot_sec5;
  TH2F *hist_HTCC_Nphe_vs_momentum_prot_sec6;

  TH2F *hist_HTCC_Nphe_vs_momentum_pip;
  TH2F *hist_HTCC_Nphe_vs_momentum_pip_sec1;
  TH2F *hist_HTCC_Nphe_vs_momentum_pip_sec2;
  TH2F *hist_HTCC_Nphe_vs_momentum_pip_sec3;
  TH2F *hist_HTCC_Nphe_vs_momentum_pip_sec4;
  TH2F *hist_HTCC_Nphe_vs_momentum_pip_sec5;
  TH2F *hist_HTCC_Nphe_vs_momentum_pip_sec6;

  TH2F *hist_HTCC_Nphe_vs_beta;
  TH2F *hist_HTCC_Nphe_vs_beta_sec1;
  TH2F *hist_HTCC_Nphe_vs_beta_sec2;
  TH2F *hist_HTCC_Nphe_vs_beta_sec3;
  TH2F *hist_HTCC_Nphe_vs_beta_sec4;
  TH2F *hist_HTCC_Nphe_vs_beta_sec5;
  TH2F *hist_HTCC_Nphe_vs_beta_sec6;

  TH2F *hist_HTCC_Nphe_vs_beta_prot;
  TH2F *hist_HTCC_Nphe_vs_beta_prot_sec1;
  TH2F *hist_HTCC_Nphe_vs_beta_prot_sec2;
  TH2F *hist_HTCC_Nphe_vs_beta_prot_sec3;
  TH2F *hist_HTCC_Nphe_vs_beta_prot_sec4;
  TH2F *hist_HTCC_Nphe_vs_beta_prot_sec5;
  TH2F *hist_HTCC_Nphe_vs_beta_prot_sec6;

  TH2F *hist_HTCC_Nphe_vs_beta_pip;
  TH2F *hist_HTCC_Nphe_vs_beta_pip_sec1;
  TH2F *hist_HTCC_Nphe_vs_beta_pip_sec2;
  TH2F *hist_HTCC_Nphe_vs_beta_pip_sec3;
  TH2F *hist_HTCC_Nphe_vs_beta_pip_sec4;
  TH2F *hist_HTCC_Nphe_vs_beta_pip_sec5;
  TH2F *hist_HTCC_Nphe_vs_beta_pip_sec6;


  hist_HTCC_Nphe_vs_momentum = new TH2F("hist_HTCC_Nphe_vs_momentum", "HTCC Nphe vs momentum", 500, 0, Ebeam+0.5, 400, 0, 200);   
  hist_HTCC_Nphe_vs_momentum->GetXaxis()->SetTitle("p /GeV");
  hist_HTCC_Nphe_vs_momentum->GetYaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_momentum_sec1 = new TH2F("hist_HTCC_Nphe_vs_momentum_sec1", "HTCC Nphe vs momentum sector 1", 500, 0, Ebeam+0.5, 400, 0, 200);  
  hist_HTCC_Nphe_vs_momentum_sec1->GetXaxis()->SetTitle("p /GeV");
  hist_HTCC_Nphe_vs_momentum_sec1->GetYaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_momentum_sec2 = new TH2F("hist_HTCC_Nphe_vs_momentum_sec2", "HTCC Nphe vs momentum sector 2", 500, 0, Ebeam+0.5, 400, 0, 200);    
  hist_HTCC_Nphe_vs_momentum_sec2->GetXaxis()->SetTitle("p /GeV");
  hist_HTCC_Nphe_vs_momentum_sec2->GetYaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_momentum_sec3 = new TH2F("hist_HTCC_Nphe_vs_momentum_sec3", "HTCC Nphe vs momentum sector 3", 500, 0, Ebeam+0.5, 400, 0, 200);    
  hist_HTCC_Nphe_vs_momentum_sec3->GetXaxis()->SetTitle("p /GeV");
  hist_HTCC_Nphe_vs_momentum_sec3->GetYaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_momentum_sec4 = new TH2F("hist_HTCC_Nphe_vs_momentum_sec4", "HTCC Nphe vs momentum sector 4", 500, 0, Ebeam+0.5, 400, 0, 200);  
  hist_HTCC_Nphe_vs_momentum_sec4->GetXaxis()->SetTitle("p /GeV");
  hist_HTCC_Nphe_vs_momentum_sec4->GetYaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_momentum_sec5 = new TH2F("hist_HTCC_Nphe_vs_momentum_sec5", "HTCC Nphe vs momentum sector 5", 500, 0, Ebeam+0.5, 400, 0, 200);    
  hist_HTCC_Nphe_vs_momentum_sec5->GetXaxis()->SetTitle("p /GeV");
  hist_HTCC_Nphe_vs_momentum_sec5->GetYaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_momentum_sec6 = new TH2F("hist_HTCC_Nphe_vs_momentum_sec6", "HTCC Nphe vs momentum sector 6", 500, 0, Ebeam+0.5, 400, 0, 200);    
  hist_HTCC_Nphe_vs_momentum_sec6->GetXaxis()->SetTitle("p /GeV");
  hist_HTCC_Nphe_vs_momentum_sec6->GetYaxis()->SetTitle("Nphe");

  hist_HTCC_Nphe_vs_momentum_prot = new TH2F("hist_HTCC_Nphe_vs_momentum_prot", "HTCC Nphe vs momentum for protons", 500, 0, Ebeam-0.5, 120, 0, 60);   
  hist_HTCC_Nphe_vs_momentum_prot->GetXaxis()->SetTitle("p /GeV");
  hist_HTCC_Nphe_vs_momentum_prot->GetYaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_momentum_prot_sec1 = new TH2F("hist_HTCC_Nphe_vs_momentum_prot_sec1", "HTCC Nphe vs momentum for protons sector 1", 500, 0, Ebeam-0.5, 120, 0, 60);    
  hist_HTCC_Nphe_vs_momentum_prot_sec1->GetXaxis()->SetTitle("p /GeV");
  hist_HTCC_Nphe_vs_momentum_prot_sec1->GetYaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_momentum_prot_sec2 = new TH2F("hist_HTCC_Nphe_vs_momentum_prot_sec2", "HTCC Nphe vs momentum for protons sector 2", 500, 0, Ebeam-0.5, 120, 0, 60);    
  hist_HTCC_Nphe_vs_momentum_prot_sec2->GetXaxis()->SetTitle("p /GeV");
  hist_HTCC_Nphe_vs_momentum_prot_sec2->GetYaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_momentum_prot_sec3 = new TH2F("hist_HTCC_Nphe_vs_momentum_prot_sec3", "HTCC Nphe vs momentum for protons sector 3", 500, 0, Ebeam-0.5, 120, 0, 60);    
  hist_HTCC_Nphe_vs_momentum_prot_sec3->GetXaxis()->SetTitle("p /GeV");
  hist_HTCC_Nphe_vs_momentum_prot_sec3->GetYaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_momentum_prot_sec4 = new TH2F("hist_HTCC_Nphe_vs_momentum_prot_sec4", "HTCC Nphe vs momentum for protons sector 4", 500, 0, Ebeam-0.5, 120, 0, 60);   
  hist_HTCC_Nphe_vs_momentum_prot_sec4->GetXaxis()->SetTitle("p /GeV");
  hist_HTCC_Nphe_vs_momentum_prot_sec4->GetYaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_momentum_prot_sec5 = new TH2F("hist_HTCC_Nphe_vs_momentum_prot_sec5", "HTCC Nphe vs momentum for protons sector 5", 500, 0, Ebeam-0.5, 120, 0, 60);    
  hist_HTCC_Nphe_vs_momentum_prot_sec5->GetXaxis()->SetTitle("p /GeV");
  hist_HTCC_Nphe_vs_momentum_prot_sec5->GetYaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_momentum_prot_sec6 = new TH2F("hist_HTCC_Nphe_vs_momentum_prot_sec6", "HTCC Nphe vs momentum for protons sector 6", 500, 0, Ebeam-0.5, 120, 0, 60);   
  hist_HTCC_Nphe_vs_momentum_prot_sec6->GetXaxis()->SetTitle("p /GeV");
  hist_HTCC_Nphe_vs_momentum_prot_sec6->GetYaxis()->SetTitle("Nphe");

  hist_HTCC_Nphe_vs_momentum_pip = new TH2F("hist_HTCC_Nphe_vs_momentum_pip", "HTCC Nphe vs momentum for pip", 500, 0, Ebeam-0.5, 120, 0, 60);    
  hist_HTCC_Nphe_vs_momentum_pip->GetXaxis()->SetTitle("p /GeV");
  hist_HTCC_Nphe_vs_momentum_pip->GetYaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_momentum_pip_sec1 = new TH2F("hist_HTCC_Nphe_vs_momentum_pip_sec1", "HTCC Nphe vs momentum for pip sector 1", 500, 0, Ebeam-0.5, 120, 0, 60);  
  hist_HTCC_Nphe_vs_momentum_pip_sec1->GetXaxis()->SetTitle("p /GeV");
  hist_HTCC_Nphe_vs_momentum_pip_sec1->GetYaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_momentum_pip_sec2 = new TH2F("hist_HTCC_Nphe_vs_momentum_pip_sec2", "HTCC Nphe vs momentum for pip sector 2", 500, 0, Ebeam-0.5, 120, 0, 60);   
  hist_HTCC_Nphe_vs_momentum_pip_sec2->GetXaxis()->SetTitle("p /GeV");
  hist_HTCC_Nphe_vs_momentum_pip_sec2->GetYaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_momentum_pip_sec3 = new TH2F("hist_HTCC_Nphe_vs_momentum_pip_sec3", "HTCC Nphe vs momentum for pip sector 3", 500, 0, Ebeam-0.5, 120, 0, 60);   
  hist_HTCC_Nphe_vs_momentum_pip_sec3->GetXaxis()->SetTitle("p /GeV");
  hist_HTCC_Nphe_vs_momentum_pip_sec3->GetYaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_momentum_pip_sec4 = new TH2F("hist_HTCC_Nphe_vs_momentum_pip_sec4", "HTCC Nphe vs momentum for pip sector 4", 500, 0, Ebeam-0.5, 120, 0, 60);   
  hist_HTCC_Nphe_vs_momentum_pip_sec4->GetXaxis()->SetTitle("p /GeV");
  hist_HTCC_Nphe_vs_momentum_pip_sec4->GetYaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_momentum_pip_sec5 = new TH2F("hist_HTCC_Nphe_vs_momentum_pip_sec5", "HTCC Nphe vs momentum for pip sector 5", 500, 0, Ebeam-0.5, 120, 0, 60);    
  hist_HTCC_Nphe_vs_momentum_pip_sec5->GetXaxis()->SetTitle("p /GeV");
  hist_HTCC_Nphe_vs_momentum_pip_sec5->GetYaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_momentum_pip_sec6 = new TH2F("hist_HTCC_Nphe_vs_momentum_pip_sec6", "HTCC Nphe vs momentum for pip sector 6", 500, 0, Ebeam-0.5, 120, 0, 60);    
  hist_HTCC_Nphe_vs_momentum_pip_sec6->GetXaxis()->SetTitle("p /GeV");
  hist_HTCC_Nphe_vs_momentum_pip_sec6->GetYaxis()->SetTitle("Nphe");


  hist_HTCC_Nphe_vs_beta = new TH2F("hist_HTCC_Nphe_vs_beta", "HTCC Nphe vs beta", 300, 0, 300, 700, 0.5, 1.2); 
  hist_HTCC_Nphe_vs_beta->GetXaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_beta->GetYaxis()->SetTitle("#beta");
  hist_HTCC_Nphe_vs_beta_sec1 = new TH2F("hist_HTCC_Nphe_vs_beta_sec1", "HTCC Nphe vs beta sector 1", 300, 0, 300, 700, 0.5, 1.2);  
  hist_HTCC_Nphe_vs_beta_sec1->GetXaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_beta_sec1->GetYaxis()->SetTitle("#beta");
  hist_HTCC_Nphe_vs_beta_sec2 = new TH2F("hist_HTCC_Nphe_vs_beta_sec2", "HTCC Nphe vs beta sector 2", 300, 0, 300, 700, 0.5, 1.2);     
  hist_HTCC_Nphe_vs_beta_sec2->GetXaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_beta_sec2->GetYaxis()->SetTitle("#beta");
  hist_HTCC_Nphe_vs_beta_sec3 = new TH2F("hist_HTCC_Nphe_vs_beta_sec3", "HTCC Nphe vs beta sector 3", 300, 0, 300, 700, 0.5, 1.2);     
  hist_HTCC_Nphe_vs_beta_sec3->GetXaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_beta_sec3->GetYaxis()->SetTitle("#beta");
  hist_HTCC_Nphe_vs_beta_sec4 = new TH2F("hist_HTCC_Nphe_vs_beta_sec4", "HTCC Nphe vs beta sector 4", 300, 0, 300, 700, 0.5, 1.2);   
  hist_HTCC_Nphe_vs_beta_sec4->GetXaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_beta_sec4->GetYaxis()->SetTitle("#beta");
  hist_HTCC_Nphe_vs_beta_sec5 = new TH2F("hist_HTCC_Nphe_vs_beta_sec5", "HTCC Nphe vs beta sector 5", 300, 0, 300, 700, 0.5, 1.2);   
  hist_HTCC_Nphe_vs_beta_sec5->GetXaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_beta_sec5->GetYaxis()->SetTitle("#beta");
  hist_HTCC_Nphe_vs_beta_sec6 = new TH2F("hist_HTCC_Nphe_vs_beta_sec6", "HTCC Nphe vs beta sector 6", 300, 0, 300, 700, 0.5, 1.2);    
  hist_HTCC_Nphe_vs_beta_sec6->GetXaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_beta_sec6->GetYaxis()->SetTitle("#beta");

  hist_HTCC_Nphe_vs_beta_prot = new TH2F("hist_HTCC_Nphe_vs_beta_prot", "HTCC Nphe vs beta for protons", 50, 0, 50, 700, 0.5, 1.2);   
  hist_HTCC_Nphe_vs_beta_prot->GetXaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_beta_prot->GetYaxis()->SetTitle("#beta");
  hist_HTCC_Nphe_vs_beta_prot_sec1 = new TH2F("hist_HTCC_Nphe_vs_beta_prot_sec1", "HTCC Nphe vs beta for protons sector 1", 50, 0, 50, 700, 0.5, 1.2);   
  hist_HTCC_Nphe_vs_beta_prot_sec1->GetXaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_beta_prot_sec1->GetYaxis()->SetTitle("#beta");
  hist_HTCC_Nphe_vs_beta_prot_sec2 = new TH2F("hist_HTCC_Nphe_vs_beta_prot_sec2", "HTCC Nphe vs beta for protons sector 2", 50, 0, 50, 700, 0.5, 1.2);     
  hist_HTCC_Nphe_vs_beta_prot_sec2->GetXaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_beta_prot_sec2->GetYaxis()->SetTitle("#beta");
  hist_HTCC_Nphe_vs_beta_prot_sec3 = new TH2F("hist_HTCC_Nphe_vs_beta_prot_sec3", "HTCC Nphe vs beta for protons sector 3", 50, 0, 50, 700, 0.5, 1.2);     
  hist_HTCC_Nphe_vs_beta_prot_sec3->GetXaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_beta_prot_sec3->GetYaxis()->SetTitle("#beta");
  hist_HTCC_Nphe_vs_beta_prot_sec4 = new TH2F("hist_HTCC_Nphe_vs_beta_prot_sec4", "HTCC Nphe vs beta for protons sector 4", 50, 0, 50, 700, 0.5, 1.2);   
  hist_HTCC_Nphe_vs_beta_prot_sec4->GetXaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_beta_prot_sec4->GetYaxis()->SetTitle("#beta");
  hist_HTCC_Nphe_vs_beta_prot_sec5 = new TH2F("hist_HTCC_Nphe_vs_beta_prot_sec5", "HTCC Nphe vs beta for protons sector 5", 50, 0, 50, 700, 0.5, 1.2);     
  hist_HTCC_Nphe_vs_beta_prot_sec5->GetXaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_beta_prot_sec5->GetYaxis()->SetTitle("#beta");
  hist_HTCC_Nphe_vs_beta_prot_sec6 = new TH2F("hist_HTCC_Nphe_vs_beta_prot_sec6", "HTCC Nphe vs beta for protons sector 6", 50, 0, 50, 700, 0.5, 1.2);    
  hist_HTCC_Nphe_vs_beta_prot_sec6->GetXaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_beta_prot_sec6->GetYaxis()->SetTitle("#beta");

  hist_HTCC_Nphe_vs_beta_pip = new TH2F("hist_HTCC_Nphe_vs_beta_pip", "HTCC Nphe vs beta for pip", 100, 0, 50, 700, 0.5, 1.2);    
  hist_HTCC_Nphe_vs_beta_pip->GetXaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_beta_pip->GetYaxis()->SetTitle("#beta");
  hist_HTCC_Nphe_vs_beta_pip_sec1 = new TH2F("hist_HTCC_Nphe_vs_beta_pip_sec1", "HTCC Nphe vs beta for pip sector 1", 50, 0, 50, 700, 0.5, 1.2);  
  hist_HTCC_Nphe_vs_beta_pip_sec1->GetXaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_beta_pip_sec1->GetYaxis()->SetTitle("#beta");
  hist_HTCC_Nphe_vs_beta_pip_sec2 = new TH2F("hist_HTCC_Nphe_vs_beta_pip_sec2", "HTCC Nphe vs beta for pip sector 2", 50, 0, 50, 700, 0.5, 1.2);   
  hist_HTCC_Nphe_vs_beta_pip_sec2->GetXaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_beta_pip_sec2->GetYaxis()->SetTitle("#beta");
  hist_HTCC_Nphe_vs_beta_pip_sec3 = new TH2F("hist_HTCC_Nphe_vs_beta_pip_sec3", "HTCC Nphe vs beta for pip sector 3", 50, 0, 50, 700, 0.5, 1.2);   
  hist_HTCC_Nphe_vs_beta_pip_sec3->GetXaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_beta_pip_sec3->GetYaxis()->SetTitle("#beta");
  hist_HTCC_Nphe_vs_beta_pip_sec4 = new TH2F("hist_HTCC_Nphe_vs_beta_pip_sec4", "HTCC Nphe vs beta for pip sector 4", 50, 0, 50, 700, 0.5, 1.2);   
  hist_HTCC_Nphe_vs_beta_pip_sec4->GetXaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_beta_pip_sec4->GetYaxis()->SetTitle("#beta");
  hist_HTCC_Nphe_vs_beta_pip_sec5 = new TH2F("hist_HTCC_Nphe_vs_beta_pip_sec5", "HTCC Nphe vs beta for pip sector 5", 50, 0, 50, 700, 0.5, 1.2);    
  hist_HTCC_Nphe_vs_beta_pip_sec5->GetXaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_beta_pip_sec5->GetYaxis()->SetTitle("#beta");
  hist_HTCC_Nphe_vs_beta_pip_sec6 = new TH2F("hist_HTCC_Nphe_vs_beta_pip_sec6", "HTCC Nphe vs beta for pip sector 6", 50, 0, 50, 700, 0.5, 1.2);    
  hist_HTCC_Nphe_vs_beta_pip_sec6->GetXaxis()->SetTitle("Nphe");
  hist_HTCC_Nphe_vs_beta_pip_sec6->GetYaxis()->SetTitle("#beta");




out->mkdir("electron_kinematics_sector");				
out->cd ("electron_kinematics_sector");

  TH1F *hist_ele_theta_sec1;
  TH1F *hist_ele_theta_sec2;
  TH1F *hist_ele_theta_sec3;
  TH1F *hist_ele_theta_sec4;
  TH1F *hist_ele_theta_sec5;
  TH1F *hist_ele_theta_sec6;

  TH1F *hist_ele_W_sec1;
  TH1F *hist_ele_W_sec2;
  TH1F *hist_ele_W_sec3;
  TH1F *hist_ele_W_sec4;
  TH1F *hist_ele_W_sec5;
  TH1F *hist_ele_W_sec6;

  TH1F *hist_ele_Q2_sec1;
  TH1F *hist_ele_Q2_sec2;
  TH1F *hist_ele_Q2_sec3;
  TH1F *hist_ele_Q2_sec4;
  TH1F *hist_ele_Q2_sec5;
  TH1F *hist_ele_Q2_sec6;

  TH2F *hist_ele_W_vs_Q2_sec1;
  TH2F *hist_ele_W_vs_Q2_sec2;
  TH2F *hist_ele_W_vs_Q2_sec3;
  TH2F *hist_ele_W_vs_Q2_sec4;
  TH2F *hist_ele_W_vs_Q2_sec5;
  TH2F *hist_ele_W_vs_Q2_sec6;

  hist_ele_theta_sec1 = new TH1F("hist_ele_theta_sec1", "electron theta sector 1", 400, 0, 50);  
  hist_ele_theta_sec1->GetXaxis()->SetTitle("#Theta /deg");
  hist_ele_theta_sec1->GetYaxis()->SetTitle("counts");
  hist_ele_theta_sec2 = new TH1F("hist_ele_theta_sec2", "electron theta sector 2", 400, 0, 50);  
  hist_ele_theta_sec2->GetXaxis()->SetTitle("#Theta /deg");
  hist_ele_theta_sec2->GetYaxis()->SetTitle("counts");
  hist_ele_theta_sec3 = new TH1F("hist_ele_theta_sec3", "electron theta sector 3", 400, 0, 50);  
  hist_ele_theta_sec3->GetXaxis()->SetTitle("#Theta /deg");
  hist_ele_theta_sec3->GetYaxis()->SetTitle("counts");
  hist_ele_theta_sec4 = new TH1F("hist_ele_theta_sec4", "electron theta sector 4", 400, 0, 50);  
  hist_ele_theta_sec4->GetXaxis()->SetTitle("#Theta /deg");
  hist_ele_theta_sec4->GetYaxis()->SetTitle("counts");
  hist_ele_theta_sec5 = new TH1F("hist_ele_theta_sec5", "electron theta sector 5", 400, 0, 50);  
  hist_ele_theta_sec5->GetXaxis()->SetTitle("#Theta /deg");
  hist_ele_theta_sec5->GetYaxis()->SetTitle("counts");
  hist_ele_theta_sec6 = new TH1F("hist_ele_theta_sec6", "electron theta sector 6", 400, 0, 50);  
  hist_ele_theta_sec6->GetXaxis()->SetTitle("#Theta /deg");
  hist_ele_theta_sec6->GetYaxis()->SetTitle("counts");

  hist_ele_W_sec1 = new TH1F("hist_ele_W_sec1", "electron W sector 1", 500, 0, 5);  
  hist_ele_W_sec1->GetXaxis()->SetTitle("W /GeV");
  hist_ele_W_sec1->GetYaxis()->SetTitle("counts");
  hist_ele_W_sec2 = new TH1F("hist_ele_W_sec2", "electron W sector 2", 500, 0, 5);  
  hist_ele_W_sec2->GetXaxis()->SetTitle("W /GeV");
  hist_ele_W_sec2->GetYaxis()->SetTitle("counts");
  hist_ele_W_sec3 = new TH1F("hist_ele_W_sec3", "electron W sector 3", 500, 0, 5);  
  hist_ele_W_sec3->GetXaxis()->SetTitle("W /GeV");
  hist_ele_W_sec3->GetYaxis()->SetTitle("counts");
  hist_ele_W_sec4 = new TH1F("hist_ele_W_sec4", "electron W sector 4", 500, 0, 5);  
  hist_ele_W_sec4->GetXaxis()->SetTitle("W /GeV");
  hist_ele_W_sec4->GetYaxis()->SetTitle("counts");
  hist_ele_W_sec5 = new TH1F("hist_ele_W_sec5", "electron W sector 5", 500, 0, 5);  
  hist_ele_W_sec5->GetXaxis()->SetTitle("W /GeV");
  hist_ele_W_sec5->GetYaxis()->SetTitle("counts");
  hist_ele_W_sec6 = new TH1F("hist_ele_W_sec6", "electron W sector 6", 500, 0, 5);  
  hist_ele_W_sec6->GetXaxis()->SetTitle("W /GeV");
  hist_ele_W_sec6->GetYaxis()->SetTitle("counts");

  hist_ele_Q2_sec1 = new TH1F("hist_ele_Q2_sec1", "electron Q^{2} sector 1", 600, 0, 12);  
  hist_ele_Q2_sec1->GetXaxis()->SetTitle("Q^{2} /GeV^{2}");
  hist_ele_Q2_sec1->GetYaxis()->SetTitle("counts");
  hist_ele_Q2_sec2 = new TH1F("hist_ele_Q2_sec2", "electron Q^{2} sector 2", 600, 0, 12);  
  hist_ele_Q2_sec2->GetXaxis()->SetTitle("Q^{2} /GeV^{2}");
  hist_ele_Q2_sec2->GetYaxis()->SetTitle("counts");
  hist_ele_Q2_sec3 = new TH1F("hist_ele_Q2_sec3", "electron Q^{2} sector 3", 600, 0, 12);  
  hist_ele_Q2_sec3->GetXaxis()->SetTitle("Q^{2} /GeV^{2}");
  hist_ele_Q2_sec3->GetYaxis()->SetTitle("counts");
  hist_ele_Q2_sec4 = new TH1F("hist_ele_Q2_sec4", "electron Q^{2} sector 4", 600, 0, 12);  
  hist_ele_Q2_sec4->GetXaxis()->SetTitle("Q^{2} /GeV^{2}");
  hist_ele_Q2_sec4->GetYaxis()->SetTitle("counts");
  hist_ele_Q2_sec5 = new TH1F("hist_ele_Q2_sec5", "electron Q^{2} sector 5", 600, 0, 12);  
  hist_ele_Q2_sec5->GetXaxis()->SetTitle("Q^{2} /GeV^{2}");
  hist_ele_Q2_sec5->GetYaxis()->SetTitle("counts");
  hist_ele_Q2_sec6 = new TH1F("hist_ele_Q2_sec6", "electron Q^{2} sector 6", 600, 0, 12);  
  hist_ele_Q2_sec6->GetXaxis()->SetTitle("Q^{2} /GeV^{2}");
  hist_ele_Q2_sec6->GetYaxis()->SetTitle("counts");

  hist_ele_W_vs_Q2_sec1 = new TH2F("hist_ele_W_vs_Q2_sec1", "electron W vs Q^{2} sector 1", 500, 0, 5, 600, 0, 12);  
  hist_ele_W_vs_Q2_sec1->GetXaxis()->SetTitle("W /GeV");
  hist_ele_W_vs_Q2_sec1->GetYaxis()->SetTitle("Q^{2} /GeV^{2}");
  hist_ele_W_vs_Q2_sec2 = new TH2F("hist_ele_W_vs_Q2_sec2", "electron W vs Q^{2} sector 2", 500, 0, 5, 600, 0, 12);  
  hist_ele_W_vs_Q2_sec2->GetXaxis()->SetTitle("W /GeV");
  hist_ele_W_vs_Q2_sec2->GetYaxis()->SetTitle("Q^{2} /GeV^{2}");
  hist_ele_W_vs_Q2_sec3 = new TH2F("hist_ele_W_vs_Q2_sec3", "electron W vs Q^{2} sector 3", 500, 0, 5, 600, 0, 12);  
  hist_ele_W_vs_Q2_sec3->GetXaxis()->SetTitle("W /GeV");
  hist_ele_W_vs_Q2_sec3->GetYaxis()->SetTitle("Q^{2} /GeV^{2}");
  hist_ele_W_vs_Q2_sec4 = new TH2F("hist_ele_W_vs_Q2_sec4", "electron W vs Q^{2} sector 4", 500, 0, 5, 600, 0, 12);  
  hist_ele_W_vs_Q2_sec4->GetXaxis()->SetTitle("W /GeV");
  hist_ele_W_vs_Q2_sec4->GetYaxis()->SetTitle("Q^{2} /GeV^{2}");
  hist_ele_W_vs_Q2_sec5 = new TH2F("hist_ele_W_vs_Q2_sec5", "electron W vs Q^{2} sector 5", 500, 0, 5, 600, 0, 12);  
  hist_ele_W_vs_Q2_sec5->GetXaxis()->SetTitle("W /GeV");
  hist_ele_W_vs_Q2_sec5->GetYaxis()->SetTitle("Q^{2} /GeV^{2}");
  hist_ele_W_vs_Q2_sec6 = new TH2F("hist_ele_W_vs_Q2_sec6", "electron W vs Q^{2} sector 6", 500, 0, 5, 600, 0, 12);  
  hist_ele_W_vs_Q2_sec6->GetXaxis()->SetTitle("W /GeV");
  hist_ele_W_vs_Q2_sec6->GetYaxis()->SetTitle("Q^{2} /GeV^{2}");



/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

out->mkdir("W_bins");				
out->cd ("W_bins");

  char name[100];
  char title[100];

  TH1F *hist_W_binned_sec_01[17]; 
  TH1F *hist_W_binned_sec_02[17]; 
  TH1F *hist_W_binned_sec_03[17]; 
  TH1F *hist_W_binned_sec_04[17]; 
  TH1F *hist_W_binned_sec_05[17]; 
  TH1F *hist_W_binned_sec_06[17];

  for(Int_t i = 0; i < 17; i++){
    sprintf(name,"W_sec1_thetabin_%02d", i);
    sprintf(title,"W sector 1 thetabin %02d", i);
    hist_W_binned_sec_01[i] = new TH1F(name, title, 1000, 0.0, Ebeam-0.2);   
    hist_W_binned_sec_01[i]->GetXaxis()->SetTitle("W");
    hist_W_binned_sec_01[i]->GetYaxis()->SetTitle("counts");
    sprintf(name,"W_sec2_thetabin_%02d", i);
    sprintf(title,"W sector 2 thetabin %02d", i);
    hist_W_binned_sec_02[i] = new TH1F(name, title, 1000, 0.0, Ebeam-0.2);   
    hist_W_binned_sec_02[i]->GetXaxis()->SetTitle("W");
    hist_W_binned_sec_02[i]->GetYaxis()->SetTitle("counts");
    sprintf(name,"W_sec3_thetabin_%02d", i);
    sprintf(title,"W sector 3 thetabin %02d", i);
    hist_W_binned_sec_03[i] = new TH1F(name, title, 1000, 0.0, Ebeam-0.2);   
    hist_W_binned_sec_03[i]->GetXaxis()->SetTitle("W");
    hist_W_binned_sec_03[i]->GetYaxis()->SetTitle("counts");
    sprintf(name,"W_sec4_thetabin_%02d", i);
    sprintf(title,"W sector 4 thetabin %02d", i);
    hist_W_binned_sec_04[i] = new TH1F(name, title, 1000, 0.0, Ebeam-0.2);   
    hist_W_binned_sec_04[i]->GetXaxis()->SetTitle("W");
    hist_W_binned_sec_04[i]->GetYaxis()->SetTitle("counts");
    sprintf(name,"W_sec5_thetabin_%02d", i);
    sprintf(title,"W sector 5 thetabin %02d", i);
    hist_W_binned_sec_05[i] = new TH1F(name, title, 1000, 0.0, Ebeam-0.2);   
    hist_W_binned_sec_05[i]->GetXaxis()->SetTitle("W");
    hist_W_binned_sec_05[i]->GetYaxis()->SetTitle("counts");
    sprintf(name,"W_sec6_thetabin_%02d", i);
    sprintf(title,"W sector 6 thetabin %02d", i);
    hist_W_binned_sec_06[i] = new TH1F(name, title, 1000, 0.0, Ebeam-0.2);   
    hist_W_binned_sec_06[i]->GetXaxis()->SetTitle("W");
    hist_W_binned_sec_06[i]->GetYaxis()->SetTitle("counts");
  }



/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

out->mkdir("momentum_correction_electron");				
out->cd ("momentum_correction_electron");


  TH2F *hist_all_electron_p_vs_theta_sec1;
  TH2F *hist_all_electron_p_vs_theta_sec2;
  TH2F *hist_all_electron_p_vs_theta_sec3;
  TH2F *hist_all_electron_p_vs_theta_sec4;
  TH2F *hist_all_electron_p_vs_theta_sec5;
  TH2F *hist_all_electron_p_vs_theta_sec6;

  hist_all_electron_p_vs_theta_sec1 = new TH2F("hist_all_electron_p_vs_theta_sec_01", "electron p vs #Theta sector 1", 200,0,37,500,0,Ebeam+1);   
  hist_all_electron_p_vs_theta_sec1->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_electron_p_vs_theta_sec1->GetYaxis()->SetTitle("p /GeV");
  hist_all_electron_p_vs_theta_sec2 = new TH2F("hist_all_electron_p_vs_theta_sec_02", "electron p vs #Theta sector 2", 200,0,37,500,0,Ebeam+1);   
  hist_all_electron_p_vs_theta_sec2->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_electron_p_vs_theta_sec2->GetYaxis()->SetTitle("p /GeV");
  hist_all_electron_p_vs_theta_sec3 = new TH2F("hist_all_electron_p_vs_theta_sec_03", "electron p vs #Theta sector 3", 200,0,37,500,0,Ebeam+1);   
  hist_all_electron_p_vs_theta_sec3->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_electron_p_vs_theta_sec3->GetYaxis()->SetTitle("p /GeV");
  hist_all_electron_p_vs_theta_sec4 = new TH2F("hist_all_electron_p_vs_theta_sec_04", "electron p vs #Theta sector 4", 200,0,37,500,0,Ebeam+1);   
  hist_all_electron_p_vs_theta_sec4->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_electron_p_vs_theta_sec4->GetYaxis()->SetTitle("p /GeV");
  hist_all_electron_p_vs_theta_sec5 = new TH2F("hist_all_electron_p_vs_theta_sec_05", "electron p vs #Theta sector 5", 200,0,37,500,0,Ebeam+1);   
  hist_all_electron_p_vs_theta_sec5->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_electron_p_vs_theta_sec5->GetYaxis()->SetTitle("p /GeV");
  hist_all_electron_p_vs_theta_sec6 = new TH2F("hist_all_electron_p_vs_theta_sec_06", "electron p vs #Theta sector 6", 200,0,37,500,0,Ebeam+1);   
  hist_all_electron_p_vs_theta_sec6->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_electron_p_vs_theta_sec6->GetYaxis()->SetTitle("p /GeV");


  TH2F *hist_theta_vs_phi;

  //Float_t bins_theta[] = {5, 7, 8 , 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 23, 28, 35};   // 2383

  //Float_t bins_theta[] = {10, 12, 13 , 14, 15, 16, 17, 18, 19, 20, 22, 24, 26, 29, 32, 37};  // 2549

  //Float_t bins_theta[] = {5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 18, 21, 25, 30};  // 2395 + 2397 (s = -1, t = +0.6)

  //Float_t bins_theta[] = {10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 24, 26, 29, 32, 37};  // 3050

  //Float_t bins_theta[] = {5, 7, 8 , 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 23, 28, 35};   // 2391

  Float_t bins_theta[] = {10, 12, 13 , 14, 15, 16, 17, 18, 19, 20, 22, 24, 26, 29, 32, 37};  // 2587


  Int_t  bincount_theta = sizeof(bins_theta)/sizeof(Float_t) - 1;

  Float_t bins_phi[] = {-180, -177, -174, -171, -168, -165, -162, -159, -156, -153, -150, -147, -144, -141, -138, -135, -132, -129, -126, -123, 
                         -120, -117, -114, -111, -108, -105, -102,  -99,  -96,  -93,  -90,  -87,  -84,  -81,  -78,  -75,  -72,  -69,  -66,  -63, 
                          -60,  -57,  -54,  -51,  -48,  -45,  -42,  -39,  -36,  -33,  -30,  -27,  -24,  -21,  -18,  -15,  -12,   -9,   -6,   -3, 
                            0,    3,    6,    9,   12,   15,   18,   21,   24,   27,   30,   33,   36,   39,   42,   45,   48,   51,   54,   57,  
                           60,   63,   66,   69,   72,   75,   78,   81,   84,   87,   90,   93,   96,   99,  102,  105,  108,  111,  114,  117, 
                          120,  123,  126,  129,  132,  135,  138,  141,  144,  147,  150,  153,  156,  159,  162,  165,  168,  171,  174,  177, 180};

  Int_t  bincount_phi = sizeof(bins_phi)/sizeof(Float_t) - 1;

  hist_theta_vs_phi = new TH2F("theta_vs_phi", "#Theta vs #Phi", bincount_phi, bins_phi, bincount_theta, bins_theta);   
  hist_theta_vs_phi->GetXaxis()->SetTitle("#Phi /deg");
  hist_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");


  TH1F *hist_W_raw;
  TH1F *hist_W_raw_sec_01; 
  TH1F *hist_W_raw_sec_02; 
  TH1F *hist_W_raw_sec_03; 
  TH1F *hist_W_raw_sec_04; 
  TH1F *hist_W_raw_sec_05; 
  TH1F *hist_W_raw_sec_06;
  TH2F *hist_W_vs_phi_raw;

  hist_W_raw = new TH1F("W_raw", "W raw all sectors", 500, 0.0, Ebeam-0.2);   
  hist_W_raw->GetXaxis()->SetTitle("W /GeV");
  hist_W_raw->GetYaxis()->SetTitle("counts");
  hist_W_raw_sec_01 = new TH1F("W_raw_sec_01", "W raw sector 1", 500, 0.0, Ebeam-0.2);   
  hist_W_raw_sec_01->GetXaxis()->SetTitle("W");
  hist_W_raw_sec_01->GetYaxis()->SetTitle("counts");
  hist_W_raw_sec_02 = new TH1F("W_raw_sec_02", "W raw sector 2", 500, 0.0, Ebeam-0.2);   
  hist_W_raw_sec_02->GetXaxis()->SetTitle("W");
  hist_W_raw_sec_02->GetYaxis()->SetTitle("counts");
  hist_W_raw_sec_03 = new TH1F("W_raw_sec_03", "W raw sector 3", 500, 0.0, Ebeam-0.2);   
  hist_W_raw_sec_03->GetXaxis()->SetTitle("W");
  hist_W_raw_sec_03->GetYaxis()->SetTitle("counts");
  hist_W_raw_sec_04 = new TH1F("W_raw_sec_04", "W raw sector 4", 500, 0.0, Ebeam-0.2);   
  hist_W_raw_sec_04->GetXaxis()->SetTitle("W");
  hist_W_raw_sec_04->GetYaxis()->SetTitle("counts");
  hist_W_raw_sec_05 = new TH1F("W_raw_sec_05", "W raw sector 5", 500, 0.0, Ebeam-0.2);   
  hist_W_raw_sec_05->GetXaxis()->SetTitle("W");
  hist_W_raw_sec_05->GetYaxis()->SetTitle("counts");
  hist_W_raw_sec_06 = new TH1F("W_raw_sec_06", "W raw sector 6", 500, 0.0, Ebeam-0.2);   
  hist_W_raw_sec_06->GetXaxis()->SetTitle("W");
  hist_W_raw_sec_06->GetYaxis()->SetTitle("counts");
  hist_W_vs_phi_raw = new TH2F("W_vs_phi_raw", "W vs phi raw all sectors", 60, -180, 180, 500, 0.0, Ebeam-0.2);   
  hist_W_vs_phi_raw->GetXaxis()->SetTitle("#phi /deg");
  hist_W_vs_phi_raw->GetYaxis()->SetTitle("W /GeV");

  TH1F *hist_W_corr;
  TH1F *hist_W_corr_sec_01; 
  TH1F *hist_W_corr_sec_02; 
  TH1F *hist_W_corr_sec_03; 
  TH1F *hist_W_corr_sec_04; 
  TH1F *hist_W_corr_sec_05; 
  TH1F *hist_W_corr_sec_06;
  TH2F *hist_W_vs_phi_corr;

  hist_W_corr = new TH1F("W_corr", "W corr all sectors", 500, 0.0, Ebeam-0.2);   
  hist_W_corr->GetXaxis()->SetTitle("W /GeV");
  hist_W_corr->GetYaxis()->SetTitle("counts");
  hist_W_corr_sec_01 = new TH1F("W_corr_sec_01", "W corr sector 1", 500, 0.0, Ebeam-0.2);   
  hist_W_corr_sec_01->GetXaxis()->SetTitle("W");
  hist_W_corr_sec_01->GetYaxis()->SetTitle("counts");
  hist_W_corr_sec_02 = new TH1F("W_corr_sec_02", "W corr sector 2", 500, 0.0, Ebeam-0.2);   
  hist_W_corr_sec_02->GetXaxis()->SetTitle("W");
  hist_W_corr_sec_02->GetYaxis()->SetTitle("counts");
  hist_W_corr_sec_03 = new TH1F("W_corr_sec_03", "W corr sector 3", 500, 0.0, Ebeam-0.2);   
  hist_W_corr_sec_03->GetXaxis()->SetTitle("W");
  hist_W_corr_sec_03->GetYaxis()->SetTitle("counts");
  hist_W_corr_sec_04 = new TH1F("W_corr_sec_04", "W corr sector 4", 500, 0.0, Ebeam-0.2);   
  hist_W_corr_sec_04->GetXaxis()->SetTitle("W");
  hist_W_corr_sec_04->GetYaxis()->SetTitle("counts");
  hist_W_corr_sec_05 = new TH1F("W_corr_sec_05", "W corr sector 5", 500, 0.0, Ebeam-0.2);   
  hist_W_corr_sec_05->GetXaxis()->SetTitle("W");
  hist_W_corr_sec_05->GetYaxis()->SetTitle("counts");
  hist_W_corr_sec_06 = new TH1F("W_corr_sec_06", "W corr sector 6", 500, 0.0, Ebeam-0.2);   
  hist_W_corr_sec_06->GetXaxis()->SetTitle("W");
  hist_W_corr_sec_06->GetYaxis()->SetTitle("counts");
  hist_W_vs_phi_corr = new TH2F("W_vs_phi_corr", "W vs phi corr all sectors", 60, -180, 180, 500, 0.0, Ebeam-0.2);   
  hist_W_vs_phi_corr->GetXaxis()->SetTitle("#phi /deg");
  hist_W_vs_phi_corr->GetYaxis()->SetTitle("W /GeV");



  TH1F *hist_delta_P[bincount_theta][bincount_phi];   // phi is equally binned in 3 deg bins from -180 to + 180, theta goes from 5 to 41 in 2 deg bins

  for(Int_t i = 0; i < bincount_theta; i++){
    for(Int_t j = 0; j < bincount_phi; j++){
      char name[100];
      char title[100];
      sprintf(name,"Pcorr_div_pmeas_thetabin_%02d_phibin_%02d", i, j);
      sprintf(title,"P_{corr}/P_{meas} #Theta bin %02d #phi bin %02d", i, j);
      hist_delta_P[i][j] = new TH1F(name, title, 200, 0.8, 1.2);   
      hist_delta_P[i][j]->GetXaxis()->SetTitle("x = P_{corr}/P_{meas}");
      hist_delta_P[i][j]->GetYaxis()->SetTitle("counts");
    }
  }

  TH1F *hist_delta_P_corr[bincount_theta][bincount_phi];   // phi is equally binned in 3 deg bins from -180 to + 180, theta has variable binning

  for(Int_t i = 0; i < bincount_theta; i++){
    for(Int_t j = 0; j < bincount_phi; j++){
      char name[100];
      char title[100];
      sprintf(name,"corr_Pcorr_div_pmeas_thetabin_%02d_phibin_%02d", i, j);
      sprintf(title,"P_{corr}/P_{meas} after correction #Theta bin %02d #phi bin %02d", i, j);
      hist_delta_P_corr[i][j] = new TH1F(name, title, 200, 0.8, 1.2);   
      hist_delta_P_corr[i][j]->GetXaxis()->SetTitle("x = P_{corr}/P_{meas}");
      hist_delta_P_corr[i][j]->GetYaxis()->SetTitle("counts");
    }
  }



out->mkdir("theta_correction_electron");				
out->cd ("theta_correction_electron");


  TH2F *hist_p_vs_phi;

  //Float_t bins_p[] = {0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7, 3.9, 4.1, 4.3, 4.5,
  //                    4.7, 4.9, 5.1, 5.3, 5.5, 5.7, 5.9, 6.1, 6.3, 6.5, 6.7, 6.9, 7.1, 7.3, 7.5, 7.7, 7.9, 8.1, 8.3, 8.5, 8.7,
  //                    8.9, 9.1, 9.3, 9.5, 9.7, 9.9, 10.1, 10.3, 10.5, 10.7};

  Float_t bins_p[] = {0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7, 3.9, 4.1, 4.3, 4.5,
                      4.7, 4.9, 5.1, 5.3, 5.5, 5.7, 5.9, 6.1, 6.3, 6.5, 6.7, 6.9, 7.1, 7.3, 7.5, 7.7, 7.9, 8.1, 8.3, 8.5, 8.7,
                      8.9, 9.1, 9.3, 9.5, 9.7, 9.9, 10.1, 10.3, 10.5, 10.7};

  Int_t  bincount_p = sizeof(bins_p)/sizeof(Float_t) - 1;

  hist_p_vs_phi = new TH2F("p_vs_phi", "P vs #Phi", bincount_phi, bins_phi, bincount_p, bins_p);   
  hist_p_vs_phi->GetXaxis()->SetTitle("#Phi /deg");
  hist_p_vs_phi->GetYaxis()->SetTitle("P /GeV");


  TH1F *hist_delta_Theta[bincount_p][bincount_phi];   // phi is equally binned in 3 deg bins from -180 to + 180, p has 0.2 GeV binning

  for(Int_t i = 0; i < bincount_p; i++){
    for(Int_t j = 0; j < bincount_phi; j++){
      char name[100];
      char title[100];
      sprintf(name,"Pcorr_div_pmeas_pbin_%02d_phibin_%02d", i, j);
      sprintf(title,"#Theta_{corr}/#Theta_{meas} P bin %02d #phi bin %02d", i, j);
      hist_delta_Theta[i][j] = new TH1F(name, title, 200, 0.8, 1.2);   
      hist_delta_Theta[i][j]->GetXaxis()->SetTitle("#Delta #Theta /deg");
      hist_delta_Theta[i][j]->GetYaxis()->SetTitle("counts");
    }
  }


/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

out->mkdir("Calorimeter_energy_depositions");				
out->cd ("Calorimeter_energy_depositions");

  TH1F *hist_CAL_Edep_electron;
  TH1F *hist_CAL_Edep_proton;
  TH1F *hist_CAL_Edep_pip;
  TH1F *hist_CAL_Edep_pim;
  TH1F *hist_CAL_Edep_pip_diffsecele;
  TH1F *hist_CAL_Edep_pim_diffsecele;

  TH1F *hist_PCAL_Edep_electron;
  TH1F *hist_PCAL_Edep_proton;
  TH1F *hist_PCAL_Edep_pip;
  TH1F *hist_PCAL_Edep_pim;
  TH1F *hist_PCAL_Edep_pip_diffsecele;
  TH1F *hist_PCAL_Edep_pim_diffsecele;

  TH1F *hist_ECin_Edep_electron;
  TH1F *hist_ECin_Edep_proton;
  TH1F *hist_ECin_Edep_pip;
  TH1F *hist_ECin_Edep_pim;
  TH1F *hist_ECin_Edep_pip_diffsecele;
  TH1F *hist_ECin_Edep_pim_diffsecele;

  TH1F *hist_ECout_Edep_electron;
  TH1F *hist_ECout_Edep_proton;
  TH1F *hist_ECout_Edep_pip;
  TH1F *hist_ECout_Edep_pim;
  TH1F *hist_ECout_Edep_pip_diffsecele;
  TH1F *hist_ECout_Edep_pim_diffsecele;

  TH1F *hist_ECAL_Edep_electron;
  TH1F *hist_ECAL_Edep_proton;
  TH1F *hist_ECAL_Edep_pip;
  TH1F *hist_ECAL_Edep_pim;
  TH1F *hist_ECAL_Edep_pip_diffsecele;
  TH1F *hist_ECAL_Edep_pim_diffsecele;

  hist_CAL_Edep_electron = new TH1F("hist_CAL_Edep_electron", "Energy deposition of electrons in the complete Calorimeter", 500, 0, 2.5);   
  hist_CAL_Edep_electron->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_CAL_Edep_electron->GetYaxis()->SetTitle("counts");
  hist_CAL_Edep_proton = new TH1F("hist_CAL_Edep_proton", "Energy deposition of protons in the complete Calorimeter", 500, 0, 1);   
  hist_CAL_Edep_proton->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_CAL_Edep_proton->GetYaxis()->SetTitle("counts");
  hist_CAL_Edep_pip = new TH1F("hist_CAL_Edep_pip", "Energy deposition of pip in the complete Calorimeter", 600, 0, 1.5);   
  hist_CAL_Edep_pip->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_CAL_Edep_pip->GetYaxis()->SetTitle("counts");
  hist_CAL_Edep_pim = new TH1F("hist_CAL_Edep_pim", "Energy deposition of pim in the complete Calorimeter", 600, 0, 1.5);   
  hist_CAL_Edep_pim->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_CAL_Edep_pim->GetYaxis()->SetTitle("counts");
  hist_CAL_Edep_pip_diffsecele = new TH1F("hist_CAL_Edep_pip_diffsecele", "Energy deposition of pip (in a different sector than electrons) in the complete Calorimeter", 600, 0, 1.5);   
  hist_CAL_Edep_pip_diffsecele->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_CAL_Edep_pip_diffsecele->GetYaxis()->SetTitle("counts");
  hist_CAL_Edep_pim_diffsecele = new TH1F("hist_CAL_Edep_pim_diffsecele", "Energy deposition of pim (in a different sector than electrons) in the complete Calorimeter", 600, 0, 1.5);   
  hist_CAL_Edep_pim_diffsecele->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_CAL_Edep_pim_diffsecele->GetYaxis()->SetTitle("counts");

  hist_PCAL_Edep_electron = new TH1F("hist_PCAL_Edep_electron", "Energy deposition of electrons in PCAL", 400, 0, 2);   
  hist_PCAL_Edep_electron->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_PCAL_Edep_electron->GetYaxis()->SetTitle("counts");
  hist_PCAL_Edep_proton = new TH1F("hist_PCAL_Edep_proton", "Energy deposition of protons in PCAL", 600, 0, 0.6);   
  hist_PCAL_Edep_proton->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_PCAL_Edep_proton->GetYaxis()->SetTitle("counts");
  hist_PCAL_Edep_pip = new TH1F("hist_PCAL_Edep_pip", "Energy deposition of pip in PCAL", 500, 0, 1);   
  hist_PCAL_Edep_pip->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_PCAL_Edep_pip->GetYaxis()->SetTitle("counts");
  hist_PCAL_Edep_pim = new TH1F("hist_PCAL_Edep_pim", "Energy deposition of pim in PCAL", 500, 0, 1);   
  hist_PCAL_Edep_pim->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_PCAL_Edep_pim->GetYaxis()->SetTitle("counts");
  hist_PCAL_Edep_pip_diffsecele = new TH1F("hist_PCAL_Edep_pip_diffsecele", "Energy deposition of pip (in a different sector than electrons) in PCAL", 500, 0, 1);   
  hist_PCAL_Edep_pip_diffsecele->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_PCAL_Edep_pip_diffsecele->GetYaxis()->SetTitle("counts");
  hist_PCAL_Edep_pim_diffsecele = new TH1F("hist_PCAL_Edep_pim_diffsecele", "Energy deposition of pim (in a different sector than electrons) in PCAL", 500, 0, 1);   
  hist_PCAL_Edep_pim_diffsecele->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_PCAL_Edep_pim_diffsecele->GetYaxis()->SetTitle("counts");

  hist_ECin_Edep_electron = new TH1F("hist_ECin_Edep_electron", "Energy deposition of electrons in ECin", 400, 0, 2);   
  hist_ECin_Edep_electron->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_ECin_Edep_electron->GetYaxis()->SetTitle("counts");
  hist_ECin_Edep_proton = new TH1F("hist_ECin_Edep_proton", "Energy deposition of protons in ECin", 500, 0, 1);   
  hist_ECin_Edep_proton->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_ECin_Edep_proton->GetYaxis()->SetTitle("counts");
  hist_ECin_Edep_pip = new TH1F("hist_ECin_Edep_pip", "Energy deposition of pip in ECin", 500, 0, 1);   
  hist_ECin_Edep_pip->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_ECin_Edep_pip->GetYaxis()->SetTitle("counts");
  hist_ECin_Edep_pim = new TH1F("hist_ECin_Edep_pim", "Energy deposition of pim in ECin", 500, 0, 1);   
  hist_ECin_Edep_pim->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_ECin_Edep_pim->GetYaxis()->SetTitle("counts");
  hist_ECin_Edep_pip_diffsecele = new TH1F("hist_ECin_Edep_pip_diffsecele", "Energy deposition of pip (in a different sector than electrons) in ECin", 500, 0, 1);   
  hist_ECin_Edep_pip_diffsecele->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_ECin_Edep_pip_diffsecele->GetYaxis()->SetTitle("counts");
  hist_ECin_Edep_pim_diffsecele = new TH1F("hist_ECin_Edep_pim_diffsecele", "Energy deposition of pim (in a different sector than electrons) in ECin", 500, 0, 1);   
  hist_ECin_Edep_pim_diffsecele->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_ECin_Edep_pim_diffsecele->GetYaxis()->SetTitle("counts");

  hist_ECout_Edep_electron = new TH1F("hist_ECout_Edep_electron", "Energy deposition of electrons in ECout", 400, 0, 2);   
  hist_ECout_Edep_electron->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_ECout_Edep_electron->GetYaxis()->SetTitle("counts");
  hist_ECout_Edep_proton = new TH1F("hist_ECout_Edep_proton", "Energy deposition of protons in ECout", 500, 0, 1);   
  hist_ECout_Edep_proton->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_ECout_Edep_proton->GetYaxis()->SetTitle("counts");
  hist_ECout_Edep_pip = new TH1F("hist_ECout_Edep_pip", "Energy deposition of pip in ECout", 500, 0, 1);   
  hist_ECout_Edep_pip->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_ECout_Edep_pip->GetYaxis()->SetTitle("counts");
  hist_ECout_Edep_pim = new TH1F("hist_ECout_Edep_pim", "Energy deposition of pim in ECout", 500, 0, 1);   
  hist_ECout_Edep_pim->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_ECout_Edep_pim->GetYaxis()->SetTitle("counts");
  hist_ECout_Edep_pip_diffsecele = new TH1F("hist_ECout_Edep_pip_diffsecele", "Energy deposition of pip (in a different sector than electrons) in ECout", 500, 0, 1);   
  hist_ECout_Edep_pip_diffsecele->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_ECout_Edep_pip_diffsecele->GetYaxis()->SetTitle("counts");
  hist_ECout_Edep_pim_diffsecele = new TH1F("hist_ECout_Edep_pim_diffsecele", "Energy deposition of pim (in a different sector than electrons) in ECout", 500, 0, 1);   
  hist_ECout_Edep_pim_diffsecele->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_ECout_Edep_pim_diffsecele->GetYaxis()->SetTitle("counts");

  hist_ECAL_Edep_electron = new TH1F("hist_ECAL_Edep_electron", "Energy deposition of electrons in ECAL", 400, 0, 2);   
  hist_ECAL_Edep_electron->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_ECAL_Edep_electron->GetYaxis()->SetTitle("counts");
  hist_ECAL_Edep_proton = new TH1F("hist_ECAL_Edep_proton", "Energy deposition of protons in ECAL", 500, 0, 1);   
  hist_ECAL_Edep_proton->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_ECAL_Edep_proton->GetYaxis()->SetTitle("counts");
  hist_ECAL_Edep_pip = new TH1F("hist_ECAL_Edep_pip", "Energy deposition of pip in ECAL", 500, 0, 1);   
  hist_ECAL_Edep_pip->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_ECAL_Edep_pip->GetYaxis()->SetTitle("counts");
  hist_ECAL_Edep_pim = new TH1F("hist_ECAL_Edep_pim", "Energy deposition of pim in ECAL", 500, 0, 1);   
  hist_ECAL_Edep_pim->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_ECAL_Edep_pim->GetYaxis()->SetTitle("counts");
  hist_ECAL_Edep_pip_diffsecele = new TH1F("hist_ECAL_Edep_pip_diffsecele", "Energy deposition of pip (in a different sector than electrons) in ECAL", 500, 0, 1);   
  hist_ECAL_Edep_pip_diffsecele->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_ECAL_Edep_pip_diffsecele->GetYaxis()->SetTitle("counts");
  hist_ECAL_Edep_pim_diffsecele = new TH1F("hist_ECAL_Edep_pim_diffsecele", "Energy deposition of pim (in a different sector than electrons) in ECAL", 500, 0, 1);   
  hist_ECAL_Edep_pim_diffsecele->GetXaxis()->SetTitle("E_{dep} /GeV");
  hist_ECAL_Edep_pim_diffsecele->GetYaxis()->SetTitle("counts");



out->mkdir("FT_gamma_energies");				
out->cd ("FT_gamma_energies");

  TH1F *hist_FT_photon_E; 
  TH1F *hist_FT_photon_theta;
  TH2F *hist_FT_photon_E_vs_theta;

  TH1F *hist_FT_photon_E_pi0; 
  TH1F *hist_FT_photon_theta_pi0;
  TH2F *hist_FT_photon_E_vs_theta_pi0;


  hist_FT_photon_E = new TH1F("hist_FT_photon_E", "photon FT energy", 500,0,10.7);   
  hist_FT_photon_E->GetXaxis()->SetTitle("E /GeV");
  hist_FT_photon_E->GetYaxis()->SetTitle("counts");
  hist_FT_photon_theta = new TH1F("hist_FT_ photon_theta", "photon FT #Theta", 50,1,6);   
  hist_FT_photon_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_FT_photon_theta->GetYaxis()->SetTitle("counts");
  hist_FT_photon_E_vs_theta = new TH2F("hist_FT_photon_E_vs_theta", "photon FT E vs #Theta", 50,0,6, 500,0,10.7);   
  hist_FT_photon_E_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_FT_photon_E_vs_theta->GetYaxis()->SetTitle("E /GeV");

  hist_FT_photon_E_pi0 = new TH1F("hist_FT_photon_E_pi0", "photon FT energy pi0", 500,0,10.7);   
  hist_FT_photon_E_pi0->GetXaxis()->SetTitle("E /GeV");
  hist_FT_photon_E_pi0->GetYaxis()->SetTitle("counts");
  hist_FT_photon_theta_pi0 = new TH1F("hist_FT_ photon_theta_pi0", "photon FT #Theta pi0", 50,1,6);   
  hist_FT_photon_theta_pi0->GetXaxis()->SetTitle("#Theta /deg");
  hist_FT_photon_theta_pi0->GetYaxis()->SetTitle("counts");
  hist_FT_photon_E_vs_theta_pi0 = new TH2F("hist_FT_photon_E_vs_theta_pi0", "photon FT E vs #Theta pi0", 50,1,6, 500,0,10.7);   
  hist_FT_photon_E_vs_theta_pi0->GetXaxis()->SetTitle("#Theta /deg");
  hist_FT_photon_E_vs_theta_pi0->GetYaxis()->SetTitle("E /GeV");



/// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

  hist_electron_count = new TH1F("hist_electron_count", "electron count per event", 10, -0.5, 9.5);   
  hist_electron_count->GetXaxis()->SetTitle("number of electrons per event");
  hist_electron_count->GetYaxis()->SetTitle("counts");
  hist_proton_count = new TH1F("hist_proton_count", "proton count per event", 10, -0.5, 9.5);   
  hist_proton_count->GetXaxis()->SetTitle("number of protons per event");
  hist_proton_count->GetYaxis()->SetTitle("counts");
  hist_neutron_count = new TH1F("hist_neutron_count", "neutron count per event", 10, -0.5, 9.5);   
  hist_neutron_count->GetXaxis()->SetTitle("number of neutrons per event");
  hist_neutron_count->GetYaxis()->SetTitle("counts");
  hist_pip_count = new TH1F("hist_pip_count", "pip count per event", 10, -0.5, 9.5);   
  hist_pip_count->GetXaxis()->SetTitle("number of pips per event");
  hist_pip_count->GetYaxis()->SetTitle("counts");
  hist_pim_count = new TH1F("hist_pim_count", "pim count per event", 10, -0.5, 9.5);   
  hist_pim_count->GetXaxis()->SetTitle("number of pims per event");
  hist_pim_count->GetYaxis()->SetTitle("counts");
  hist_Kp_count = new TH1F("hist_Kp_count", "Kp count per event", 10, -0.5, 9.5);   
  hist_Kp_count->GetXaxis()->SetTitle("number of Kps per event");
  hist_Kp_count->GetYaxis()->SetTitle("counts");
  hist_Km_count = new TH1F("hist_Km_count", "Km count per event", 10, -0.5, 9.5);   
  hist_Km_count->GetXaxis()->SetTitle("number of Kms per event");
  hist_Km_count->GetYaxis()->SetTitle("counts");
  hist_photon_count = new TH1F("hist_photon_count", "photon count per event", 10, -0.5, 9.5);   
  hist_photon_count->GetXaxis()->SetTitle("number of photons per event");
  hist_photon_count->GetYaxis()->SetTitle("counts");


out->mkdir("particles_charge");				
out->cd ("particles_charge");

  TH1F *hist_particles_p; TH1F *hist_particles_theta; TH1F *hist_particles_phi;
  TH2F *hist_particles_p_vs_theta; TH2F *hist_particles_p_vs_phi; TH2F *hist_particles_theta_vs_phi;

  TH1F *hist_positives_p; TH1F *hist_positives_theta; TH1F *hist_positives_phi;
  TH2F *hist_positives_p_vs_theta; TH2F *hist_positives_p_vs_phi; TH2F *hist_positives_theta_vs_phi;

  TH1F *hist_negatives_p; TH1F *hist_negatives_theta; TH1F *hist_negatives_phi;
  TH2F *hist_negatives_p_vs_theta; TH2F *hist_negatives_p_vs_phi; TH2F *hist_negatives_theta_vs_phi;

  hist_particles_p = new TH1F("hist_particles_p", "particles momentum", 500,0,Ebeam+1);   
  hist_particles_p->GetXaxis()->SetTitle("p /GeV");
  hist_particles_p->GetYaxis()->SetTitle("counts");
  hist_particles_theta = new TH1F("hist_particles_theta", "particles #Theta", 560,0,140);   
  hist_particles_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_particles_theta->GetYaxis()->SetTitle("counts");
  hist_particles_phi = new TH1F("hist_particles_phi", "particles #phi", 180,-180,180);   
  hist_particles_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_particles_phi->GetYaxis()->SetTitle("counts");
  hist_particles_p_vs_theta = new TH2F("hist_particles_p_vs_theta", "particles p vs #Theta", 560,0,140,500,0,Ebeam+1);   
  hist_particles_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_particles_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_particles_p_vs_phi = new TH2F("hist_particles_p_vs_phi", "particles p vs #phi", 180,-180,180, 500,0,Ebeam+1);   
  hist_particles_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_particles_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_particles_theta_vs_phi = new TH2F("hist_particles_theta_vs_phi", "particles #Theta vs #phi", 180,-180,180, 560,0,140);   
  hist_particles_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_particles_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  hist_positives_p = new TH1F("hist_positives_p", "positives momentum", 500,0,Ebeam+1);   
  hist_positives_p->GetXaxis()->SetTitle("p /GeV");
  hist_positives_p->GetYaxis()->SetTitle("counts");
  hist_positives_theta = new TH1F("hist_positives_theta", "positives #Theta", 560,0,140);   
  hist_positives_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_positives_theta->GetYaxis()->SetTitle("counts");
  hist_positives_phi = new TH1F("hist_positives_phi", "positives #phi", 180,-180,180);   
  hist_positives_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_positives_phi->GetYaxis()->SetTitle("counts");
  hist_positives_p_vs_theta = new TH2F("hist_positives_p_vs_theta", "positives p vs #Theta", 560,0,140,500,0,Ebeam+1);   
  hist_positives_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_positives_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_positives_p_vs_phi = new TH2F("hist_positives_p_vs_phi", "positives p vs #phi", 180,-180,180, 500,0,Ebeam+1);   
  hist_positives_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_positives_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_positives_theta_vs_phi = new TH2F("hist_positives_theta_vs_phi", "positives #Theta vs #phi", 180,-180,180, 560,0,140);   
  hist_positives_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_positives_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  hist_negatives_p = new TH1F("hist_negatives_p", "negatives momentum", 500,0,Ebeam+1);   
  hist_negatives_p->GetXaxis()->SetTitle("p /GeV");
  hist_negatives_p->GetYaxis()->SetTitle("counts");
  hist_negatives_theta = new TH1F("hist_negatives_theta", "negatives #Theta", 560,0,140);   
  hist_negatives_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_negatives_theta->GetYaxis()->SetTitle("counts");
  hist_negatives_phi = new TH1F("hist_negatives_phi", "negatives #phi", 180,-180,180);   
  hist_negatives_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_negatives_phi->GetYaxis()->SetTitle("counts");
  hist_negatives_p_vs_theta = new TH2F("hist_negatives_p_vs_theta", "negatives p vs #Theta", 560,0,140,500,0,Ebeam+1);   
  hist_negatives_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_negatives_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_negatives_p_vs_phi = new TH2F("hist_negatives_p_vs_phi", "negatives p vs #phi", 180,-180,180, 500,0,Ebeam+1);   
  hist_negatives_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_negatives_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_negatives_theta_vs_phi = new TH2F("hist_negatives_theta_vs_phi", "negatives #Theta vs #phi", 180,-180,180, 560,0,140);   
  hist_negatives_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_negatives_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");


out->mkdir("particles_identified_histograms_all");
out->cd ("particles_identified_histograms_all");

  TH1F *hist_all_electron_p; TH1F *hist_all_electron_theta; TH1F *hist_all_electron_phi;
  TH2F *hist_all_electron_p_vs_theta; TH2F *hist_all_electron_p_vs_phi; TH2F *hist_all_electron_theta_vs_phi;

  TH1F *hist_all_proton_p; TH1F *hist_all_proton_theta; TH1F *hist_all_proton_phi;
  TH2F *hist_all_proton_p_vs_theta; TH2F *hist_all_proton_p_vs_phi; TH2F *hist_all_proton_theta_vs_phi;
  TH1F *hist_all_neutron_p; TH1F *hist_all_neutron_theta; TH1F *hist_all_neutron_phi;
  TH2F *hist_all_neutron_p_vs_theta; TH2F *hist_all_neutron_p_vs_phi; TH2F *hist_all_neutron_theta_vs_phi;

  TH1F *hist_all_pip_p; TH1F *hist_all_pip_theta; TH1F *hist_all_pip_phi;
  TH2F *hist_all_pip_p_vs_theta; TH2F *hist_all_pip_p_vs_phi; TH2F *hist_all_pip_theta_vs_phi;
  TH1F *hist_all_pim_p; TH1F *hist_all_pim_theta; TH1F *hist_all_pim_phi;
  TH2F *hist_all_pim_p_vs_theta; TH2F *hist_all_pim_p_vs_phi; TH2F *hist_all_pim_theta_vs_phi;
 
  TH1F *hist_all_Kp_p; TH1F *hist_all_Kp_theta; TH1F *hist_all_Kp_phi;
  TH2F *hist_all_Kp_p_vs_theta; TH2F *hist_all_Kp_p_vs_phi; TH2F *hist_all_Kp_theta_vs_phi;
  TH1F *hist_all_Km_p; TH1F *hist_all_Km_theta; TH1F *hist_all_Km_phi;
  TH2F *hist_all_Km_p_vs_theta; TH2F *hist_all_Km_p_vs_phi; TH2F *hist_all_Km_theta_vs_phi;

  TH1F *hist_all_photon_p; TH1F *hist_all_photon_theta; TH1F *hist_all_photon_phi;
  TH2F *hist_all_photon_p_vs_theta; TH2F *hist_all_photon_p_vs_phi; TH2F *hist_all_photon_theta_vs_phi;


  hist_all_electron_p = new TH1F("hist_all_electron_p", "electron momentum", 500,0,Ebeam+1);   
  hist_all_electron_p->GetXaxis()->SetTitle("p /GeV");
  hist_all_electron_p->GetYaxis()->SetTitle("counts");
  hist_all_electron_theta = new TH1F("hist_all_electron_theta", "electron #Theta", 200,0,50);   
  hist_all_electron_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_electron_theta->GetYaxis()->SetTitle("counts");
  hist_all_electron_phi = new TH1F("hist_all_electron_phi", "electron #phi", 180,-180,180);   
  hist_all_electron_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_electron_phi->GetYaxis()->SetTitle("counts");
  hist_all_electron_p_vs_theta = new TH2F("hist_all_electron_p_vs_theta", "electron p vs #Theta", 200,0,50,500,0,Ebeam+1);   
  hist_all_electron_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_electron_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_all_electron_p_vs_phi = new TH2F("hist_all_electron_p_vs_phi", "electron p vs #phi", 180,-180,180, 500,0,Ebeam+1);   
  hist_all_electron_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_electron_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_all_electron_theta_vs_phi = new TH2F("hist_all_electron_theta_vs_phi", "electron #Theta vs #phi", 180,-180,180, 200,0,50);   
  hist_all_electron_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_electron_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  hist_all_proton_p = new TH1F("hist_all_proton_p", "proton momentum", 500,0,Ebeam+1);   
  hist_all_proton_p->GetXaxis()->SetTitle("p /GeV");
  hist_all_proton_p->GetYaxis()->SetTitle("counts");
  hist_all_proton_theta = new TH1F("hist_all_proton_theta", "proton #Theta", 560,0,140);   
  hist_all_proton_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_proton_theta->GetYaxis()->SetTitle("counts");
  hist_all_proton_phi = new TH1F("hist_all_proton_phi", "proton #phi", 180,-180,180);   
  hist_all_proton_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_proton_phi->GetYaxis()->SetTitle("counts");
  hist_all_proton_p_vs_theta = new TH2F("hist_all_proton_p_vs_theta", "proton p vs #Theta", 560,0,140,500,0,Ebeam+1);   
  hist_all_proton_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_proton_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_all_proton_p_vs_phi = new TH2F("hist_all_proton_p_vs_phi", "proton p vs #phi", 180,-180,180, 500,0,Ebeam+1);   
  hist_all_proton_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_proton_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_all_proton_theta_vs_phi = new TH2F("hist_all_proton_theta_vs_phi", "proton #Theta vs #phi", 180,-180,180, 560,0,140);   
  hist_all_proton_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_proton_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  hist_all_neutron_p = new TH1F("hist_all_neutron_p", "neutron momentum", 500,0,Ebeam+1);   
  hist_all_neutron_p->GetXaxis()->SetTitle("p /GeV");
  hist_all_neutron_p->GetYaxis()->SetTitle("counts");
  hist_all_neutron_theta = new TH1F("hist_all_neutron_theta", "neutron #Theta", 560,0,140);   
  hist_all_neutron_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_neutron_theta->GetYaxis()->SetTitle("counts");
  hist_all_neutron_phi = new TH1F("hist_all_neutron_phi", "neutron #phi", 180,-180,180);   
  hist_all_neutron_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_neutron_phi->GetYaxis()->SetTitle("counts");
  hist_all_neutron_p_vs_theta = new TH2F("hist_all_neutron_p_vs_theta", "neutron p vs #Theta", 560,0,140,500,0,Ebeam+1);   
  hist_all_neutron_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_neutron_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_all_neutron_p_vs_phi = new TH2F("hist_all_neutron_p_vs_phi", "neutron p vs #phi", 180,-180,180, 500,0,Ebeam+1);   
  hist_all_neutron_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_neutron_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_all_neutron_theta_vs_phi = new TH2F("hist_all_neutron_theta_vs_phi", "neutron #Theta vs #phi", 180,-180,180, 560,0,140);   
  hist_all_neutron_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_neutron_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  hist_all_pip_p = new TH1F("hist_all_pip_p", "#pi^{+} momentum", 500,0,Ebeam+1);   
  hist_all_pip_p->GetXaxis()->SetTitle("p /GeV");
  hist_all_pip_p->GetYaxis()->SetTitle("counts");
  hist_all_pip_theta = new TH1F("hist_all_pip_theta", "#pi^{+} #Theta", 560,0,140);   
  hist_all_pip_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_pip_theta->GetYaxis()->SetTitle("counts");
  hist_all_pip_phi = new TH1F("hist_all_pip_phi", "#pi^{+} #phi", 180,-180,180);   
  hist_all_pip_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_pip_phi->GetYaxis()->SetTitle("counts");
  hist_all_pip_p_vs_theta = new TH2F("hist_all_pip_p_vs_theta", "#pi^{+} p vs #Theta", 560,0,140,500,0,Ebeam+1);   
  hist_all_pip_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_pip_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_all_pip_p_vs_phi = new TH2F("hist_all_pip_p_vs_phi", "#pi^{+} p vs #phi", 180,-180,180, 500,0,Ebeam+1);   
  hist_all_pip_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_pip_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_all_pip_theta_vs_phi = new TH2F("hist_all_pip_theta_vs_phi", "#pi^{+} #Theta vs #phi", 180,-180,180, 560,0,140);   
  hist_all_pip_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_pip_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  hist_all_pim_p = new TH1F("hist_all_pim_p", "#pi^{-} momentum", 500,0,Ebeam+1);   
  hist_all_pim_p->GetXaxis()->SetTitle("p /GeV");
  hist_all_pim_p->GetYaxis()->SetTitle("counts");
  hist_all_pim_theta = new TH1F("hist_all_pim_theta", "#pi^{-} #Theta", 560,0,140);   
  hist_all_pim_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_pim_theta->GetYaxis()->SetTitle("counts");
  hist_all_pim_phi = new TH1F("hist_all_pim_phi", "#pi^{-} #phi", 180,-180,180);   
  hist_all_pim_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_pim_phi->GetYaxis()->SetTitle("counts");
  hist_all_pim_p_vs_theta = new TH2F("hist_all_pim_p_vs_theta", "#pi^{-} p vs #Theta", 560,0,140,500,0,Ebeam+1);   
  hist_all_pim_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_pim_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_all_pim_p_vs_phi = new TH2F("hist_all_pim_p_vs_phi", "#pi^{-} p vs #phi", 180,-180,180, 500,0,Ebeam+1);   
  hist_all_pim_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_pim_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_all_pim_theta_vs_phi = new TH2F("hist_all_pim_theta_vs_phi", "#pi^{-} #Theta vs #phi", 180,-180,180, 560,0,140);   
  hist_all_pim_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_pim_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  hist_all_Kp_p = new TH1F("hist_all_Kp_p", "K^{+} momentum", 500,0,Ebeam+1);   
  hist_all_Kp_p->GetXaxis()->SetTitle("p /GeV");
  hist_all_Kp_p->GetYaxis()->SetTitle("counts");
  hist_all_Kp_theta = new TH1F("hist_all_Kp_theta", "K^{+} #Theta", 560,0,140);   
  hist_all_Kp_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_Kp_theta->GetYaxis()->SetTitle("counts");
  hist_all_Kp_phi = new TH1F("hist_all_Kp_phi", "K^{+} #phi", 360,-180,180);   
  hist_all_Kp_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_Kp_phi->GetYaxis()->SetTitle("counts");
  hist_all_Kp_p_vs_theta = new TH2F("hist_all_Kp_p_vs_theta", "K^{+} p vs #Theta", 560,0,140,500,0,Ebeam+1);   
  hist_all_Kp_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_Kp_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_all_Kp_p_vs_phi = new TH2F("hist_all_Kp_p_vs_phi", "K^{+} p vs #phi", 180,-180,180, 500,0,Ebeam+1);   
  hist_all_Kp_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_Kp_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_all_Kp_theta_vs_phi = new TH2F("hist_all_Kp_theta_vs_phi", "K^{+} #Theta vs #phi", 180,-180,180, 560,0,140);   
  hist_all_Kp_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_Kp_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  hist_all_Km_p = new TH1F("hist_all_Km_p", "K^{-} momentum", 500,0,Ebeam+1);   
  hist_all_Km_p->GetXaxis()->SetTitle("p /GeV");
  hist_all_Km_p->GetYaxis()->SetTitle("counts");
  hist_all_Km_theta = new TH1F("hist_all_Km_theta", "K^{-} #Theta", 560,0,140);   
  hist_all_Km_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_Km_theta->GetYaxis()->SetTitle("counts");
  hist_all_Km_phi = new TH1F("hist_all_Km_phi", "K^{-} #phi", 180,-180,180);   
  hist_all_Km_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_Km_phi->GetYaxis()->SetTitle("counts");
  hist_all_Km_p_vs_theta = new TH2F("hist_all_Km_p_vs_theta", "K^{-} p vs #Theta", 560,0,140,500,0,Ebeam+1);   
  hist_all_Km_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_Km_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_all_Km_p_vs_phi = new TH2F("hist_all_Km_p_vs_phi", "K^{-} p vs #phi", 180,-180,180, 500,0,Ebeam+1);   
  hist_all_Km_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_Km_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_all_Km_theta_vs_phi = new TH2F("hist_all_Km_theta_vs_phi", "K^{-} #Theta vs #phi", 180,-180,180, 560,0,140);   
  hist_all_Km_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_Km_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  hist_all_photon_p = new TH1F("hist_all_photon_p", "photon momentum", 500,0,Ebeam+1);   
  hist_all_photon_p->GetXaxis()->SetTitle("p /GeV");
  hist_all_photon_p->GetYaxis()->SetTitle("counts");
  hist_all_photon_theta = new TH1F("hist_all_photon_theta", "photon #Theta", 200,0,50);   
  hist_all_photon_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_photon_theta->GetYaxis()->SetTitle("counts");
  hist_all_photon_phi = new TH1F("hist_all_photon_phi", "photon #phi", 180,-180,180);   
  hist_all_photon_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_photon_phi->GetYaxis()->SetTitle("counts");
  hist_all_photon_p_vs_theta = new TH2F("hist_all_photon_p_vs_theta", "photon p vs #Theta", 200,0,50,500,0,Ebeam+1);   
  hist_all_photon_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_photon_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_all_photon_p_vs_phi = new TH2F("hist_all_photon_p_vs_phi", "photon p vs #phi", 180,-180,180, 500,0,Ebeam+1);   
  hist_all_photon_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_photon_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_all_photon_theta_vs_phi = new TH2F("hist_all_photon_theta_vs_phi", "photon #Theta vs #phi", 180,-180,180, 200,0,50);   
  hist_all_photon_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_photon_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");


out->mkdir("particles_all_correct");
out->cd ("particles_all_correct");

  TH1F *hist_FD_all_corr_electron_p; TH1F *hist_FD_all_corr_electron_theta;
  TH1F *hist_FD_all_corr_proton_p; TH1F *hist_FD_all_corr_proton_theta;
  TH1F *hist_FD_all_corr_neutron_p; TH1F *hist_FD_all_corr_neutron_theta;
  TH1F *hist_FD_all_corr_pip_p; TH1F *hist_FD_all_corr_pip_theta;
  TH1F *hist_FD_all_corr_pim_p; TH1F *hist_FD_all_corr_pim_theta;
  TH1F *hist_FD_all_corr_Kp_p; TH1F *hist_FD_all_corr_Kp_theta;
  TH1F *hist_FD_all_corr_Km_p; TH1F *hist_FD_all_corr_Km_theta;

  TH1F *hist_FD_all_det_electron_p; TH1F *hist_FD_all_det_electron_theta;
  TH1F *hist_FD_all_det_proton_p; TH1F *hist_FD_all_det_proton_theta;
  TH1F *hist_FD_all_det_neutron_p; TH1F *hist_FD_all_det_neutron_theta;
  TH1F *hist_FD_all_det_pip_p; TH1F *hist_FD_all_det_pip_theta;
  TH1F *hist_FD_all_det_pim_p; TH1F *hist_FD_all_det_pim_theta;
  TH1F *hist_FD_all_det_Kp_p; TH1F *hist_FD_all_det_Kp_theta;
  TH1F *hist_FD_all_det_Km_p; TH1F *hist_FD_all_det_Km_theta; 

  TH1F *hist_FD_all_gen_electron_p; TH1F *hist_FD_all_gen_electron_theta;
  TH1F *hist_FD_all_gen_proton_p; TH1F *hist_FD_all_gen_proton_theta;
  TH1F *hist_FD_all_gen_neutron_p; TH1F *hist_FD_all_gen_neutron_theta;
  TH1F *hist_FD_all_gen_pip_p; TH1F *hist_FD_all_gen_pip_theta;
  TH1F *hist_FD_all_gen_pim_p; TH1F *hist_FD_all_gen_pim_theta;
  TH1F *hist_FD_all_gen_Kp_p; TH1F *hist_FD_all_gen_Kp_theta;
  TH1F *hist_FD_all_gen_Km_p; TH1F *hist_FD_all_gen_Km_theta; 

  TH1F *hist_CD_all_corr_electron_p; TH1F *hist_CD_all_corr_electron_theta;
  TH1F *hist_CD_all_corr_proton_p; TH1F *hist_CD_all_corr_proton_theta;
  TH1F *hist_CD_all_corr_neutron_p; TH1F *hist_CD_all_corr_neutron_theta;
  TH1F *hist_CD_all_corr_pip_p; TH1F *hist_CD_all_corr_pip_theta;
  TH1F *hist_CD_all_corr_pim_p; TH1F *hist_CD_all_corr_pim_theta;
  TH1F *hist_CD_all_corr_Kp_p; TH1F *hist_CD_all_corr_Kp_theta;
  TH1F *hist_CD_all_corr_Km_p; TH1F *hist_CD_all_corr_Km_theta;

  TH1F *hist_CD_all_det_electron_p; TH1F *hist_CD_all_det_electron_theta;
  TH1F *hist_CD_all_det_proton_p; TH1F *hist_CD_all_det_proton_theta;
  TH1F *hist_CD_all_det_neutron_p; TH1F *hist_CD_all_det_neutron_theta;
  TH1F *hist_CD_all_det_pip_p; TH1F *hist_CD_all_det_pip_theta;
  TH1F *hist_CD_all_det_pim_p; TH1F *hist_CD_all_det_pim_theta;
  TH1F *hist_CD_all_det_Kp_p; TH1F *hist_CD_all_det_Kp_theta;
  TH1F *hist_CD_all_det_Km_p; TH1F *hist_CD_all_det_Km_theta; 

  TH1F *hist_CD_all_gen_electron_p; TH1F *hist_CD_all_gen_electron_theta;
  TH1F *hist_CD_all_gen_proton_p; TH1F *hist_CD_all_gen_proton_theta;
  TH1F *hist_CD_all_gen_neutron_p; TH1F *hist_CD_all_gen_neutron_theta;
  TH1F *hist_CD_all_gen_pip_p; TH1F *hist_CD_all_gen_pip_theta;
  TH1F *hist_CD_all_gen_pim_p; TH1F *hist_CD_all_gen_pim_theta;
  TH1F *hist_CD_all_gen_Kp_p; TH1F *hist_CD_all_gen_Kp_theta;
  TH1F *hist_CD_all_gen_Km_p; TH1F *hist_CD_all_gen_Km_theta; 

  hist_FD_all_corr_electron_p = new TH1F("hist_FD_all_corr_electron_p", "electron momentum", 125,0,10);   
  hist_FD_all_corr_electron_p->GetXaxis()->SetTitle("p /GeV");
  hist_FD_all_corr_electron_p->GetYaxis()->SetTitle("counts");
  hist_FD_all_corr_electron_theta = new TH1F("hist_FD_all_corr_electron_theta", "electron #Theta", 300,0,150);   
  hist_FD_all_corr_electron_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_FD_all_corr_electron_theta->GetYaxis()->SetTitle("counts");
  hist_FD_all_corr_proton_p = new TH1F("hist_FD_all_corr_proton_p", "proton momentum", 125,0,10);   
  hist_FD_all_corr_proton_p->GetXaxis()->SetTitle("p /GeV");
  hist_FD_all_corr_proton_p->GetYaxis()->SetTitle("counts");
  hist_FD_all_corr_proton_theta = new TH1F("hist_FD_all_corr_proton_theta", "proton #Theta", 300,0,150);   
  hist_FD_all_corr_proton_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_FD_all_corr_proton_theta->GetYaxis()->SetTitle("counts");
  hist_FD_all_corr_neutron_p = new TH1F("hist_FD_all_corr_neutron_p", "neutron momentum", 125,0,10);   
  hist_FD_all_corr_neutron_p->GetXaxis()->SetTitle("p /GeV");
  hist_FD_all_corr_neutron_p->GetYaxis()->SetTitle("counts");
  hist_FD_all_corr_neutron_theta = new TH1F("hist_FD_all_corr_neutron_theta", "neutron #Theta", 300,0,150);   
  hist_FD_all_corr_neutron_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_FD_all_corr_neutron_theta->GetYaxis()->SetTitle("counts");
  hist_FD_all_corr_pip_p = new TH1F("hist_FD_all_corr_pip_p", "#pi^{+} momentum", 125,0,10);   
  hist_FD_all_corr_pip_p->GetXaxis()->SetTitle("p /GeV");
  hist_FD_all_corr_pip_p->GetYaxis()->SetTitle("counts");
  hist_FD_all_corr_pip_theta = new TH1F("hist_FD_all_corr_pip_theta", "#pi^{+} #Theta", 300,0,150);   
  hist_FD_all_corr_pip_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_FD_all_corr_pip_theta->GetYaxis()->SetTitle("counts");
  hist_FD_all_corr_pim_p = new TH1F("hist_FD_all_corr_pim_p", "#pi^{-} momentum", 125,0,10);   
  hist_FD_all_corr_pim_p->GetXaxis()->SetTitle("p /GeV");
  hist_FD_all_corr_pim_p->GetYaxis()->SetTitle("counts");
  hist_FD_all_corr_pim_theta = new TH1F("hist_FD_all_corr_pim_theta", "#pi^{-} #Theta", 300,0,150);   
  hist_FD_all_corr_pim_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_FD_all_corr_pim_theta->GetYaxis()->SetTitle("counts");
  hist_FD_all_corr_Kp_p = new TH1F("hist_FD_all_corr_Kp_p", "K^{+} momentum", 125,0,10);   
  hist_FD_all_corr_Kp_p->GetXaxis()->SetTitle("p /GeV");
  hist_FD_all_corr_Kp_p->GetYaxis()->SetTitle("counts");
  hist_FD_all_corr_Kp_theta = new TH1F("hist_FD_all_corr_Kp_theta", "K^{+} #Theta", 300,0,150);   
  hist_FD_all_corr_Kp_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_FD_all_corr_Kp_theta->GetYaxis()->SetTitle("counts");
  hist_FD_all_corr_Km_p = new TH1F("hist_FD_all_corr_Km_p", "K^{-} momentum", 125,0,10);   
  hist_FD_all_corr_Km_p->GetXaxis()->SetTitle("p /GeV");
  hist_FD_all_corr_Km_p->GetYaxis()->SetTitle("counts");
  hist_FD_all_corr_Km_theta = new TH1F("hist_FD_all_corr_Km_theta", "K^{-} #Theta", 300,0,150);   
  hist_FD_all_corr_Km_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_FD_all_corr_Km_theta->GetYaxis()->SetTitle("counts");

  hist_FD_all_det_electron_p = new TH1F("hist_FD_all_det_electron_p", "electron momentum", 125,0,10);   
  hist_FD_all_det_electron_p->GetXaxis()->SetTitle("p /GeV");
  hist_FD_all_det_electron_p->GetYaxis()->SetTitle("counts");
  hist_FD_all_det_electron_theta = new TH1F("hist_FD_all_det_electron_theta", "electron #Theta", 300,0,150);   
  hist_FD_all_det_electron_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_FD_all_det_electron_theta->GetYaxis()->SetTitle("counts");
  hist_FD_all_det_proton_p = new TH1F("hist_FD_all_det_proton_p", "proton momentum", 125,0,10);   
  hist_FD_all_det_proton_p->GetXaxis()->SetTitle("p /GeV");
  hist_FD_all_det_proton_p->GetYaxis()->SetTitle("counts");
  hist_FD_all_det_proton_theta = new TH1F("hist_FD_all_det_proton_theta", "proton #Theta", 300,0,150);   
  hist_FD_all_det_proton_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_FD_all_det_proton_theta->GetYaxis()->SetTitle("counts");
  hist_FD_all_det_neutron_p = new TH1F("hist_FD_all_det_neutron_p", "neutron momentum", 125,0,10);   
  hist_FD_all_det_neutron_p->GetXaxis()->SetTitle("p /GeV");
  hist_FD_all_det_neutron_p->GetYaxis()->SetTitle("counts");
  hist_FD_all_det_neutron_theta = new TH1F("hist_FD_all_det_neutron_theta", "neutron #Theta", 300,0,150);   
  hist_FD_all_det_neutron_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_FD_all_det_neutron_theta->GetYaxis()->SetTitle("counts");
  hist_FD_all_det_pip_p = new TH1F("hist_FD_all_det_pip_p", "#pi^{+} momentum", 125,0,10);   
  hist_FD_all_det_pip_p->GetXaxis()->SetTitle("p /GeV");
  hist_FD_all_det_pip_p->GetYaxis()->SetTitle("counts");
  hist_FD_all_det_pip_theta = new TH1F("hist_FD_all_det_pip_theta", "#pi^{+} #Theta", 300,0,150);   
  hist_FD_all_det_pip_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_FD_all_det_pip_theta->GetYaxis()->SetTitle("counts");
  hist_FD_all_det_pim_p = new TH1F("hist_FD_all_det_pim_p", "#pi^{-} momentum", 125,0,10);   
  hist_FD_all_det_pim_p->GetXaxis()->SetTitle("p /GeV");
  hist_FD_all_det_pim_p->GetYaxis()->SetTitle("counts");
  hist_FD_all_det_pim_theta = new TH1F("hist_FD_all_det_pim_theta", "#pi^{-} #Theta", 300,0,150);   
  hist_FD_all_det_pim_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_FD_all_det_pim_theta->GetYaxis()->SetTitle("counts");
  hist_FD_all_det_Kp_p = new TH1F("hist_FD_all_det_Kp_p", "K^{+} momentum", 125,0,10);   
  hist_FD_all_det_Kp_p->GetXaxis()->SetTitle("p /GeV");
  hist_FD_all_det_Kp_p->GetYaxis()->SetTitle("counts");
  hist_FD_all_det_Kp_theta = new TH1F("hist_FD_all_det_Kp_theta", "K^{+} #Theta", 300,0,150);   
  hist_FD_all_det_Kp_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_FD_all_det_Kp_theta->GetYaxis()->SetTitle("counts");
  hist_FD_all_det_Km_p = new TH1F("hist_FD_all_det_Km_p", "K^{-} momentum", 125,0,10);   
  hist_FD_all_det_Km_p->GetXaxis()->SetTitle("p /GeV");
  hist_FD_all_det_Km_p->GetYaxis()->SetTitle("counts");
  hist_FD_all_det_Km_theta = new TH1F("hist_FD_all_det_Km_theta", "K^{-} #Theta", 300,0,150);   
  hist_FD_all_det_Km_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_FD_all_det_Km_theta->GetYaxis()->SetTitle("counts");

  hist_FD_all_gen_electron_p = new TH1F("hist_FD_all_gen_electron_p", "electron momentum", 125,0,10);   
  hist_FD_all_gen_electron_p->GetXaxis()->SetTitle("p /GeV");
  hist_FD_all_gen_electron_p->GetYaxis()->SetTitle("counts");
  hist_FD_all_gen_electron_theta = new TH1F("hist_FD_all_gen_electron_theta", "electron #Theta", 300,0,150);   
  hist_FD_all_gen_electron_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_FD_all_gen_electron_theta->GetYaxis()->SetTitle("counts");
  hist_FD_all_gen_proton_p = new TH1F("hist_FD_all_gen_proton_p", "proton momentum", 125,0,10);   
  hist_FD_all_gen_proton_p->GetXaxis()->SetTitle("p /GeV");
  hist_FD_all_gen_proton_p->GetYaxis()->SetTitle("counts");
  hist_FD_all_gen_proton_theta = new TH1F("hist_FD_all_gen_proton_theta", "proton #Theta", 300,0,150);   
  hist_FD_all_gen_proton_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_FD_all_gen_proton_theta->GetYaxis()->SetTitle("counts");
  hist_FD_all_gen_neutron_p = new TH1F("hist_FD_all_gen_neutron_p", "neutron momentum", 125,0,10);   
  hist_FD_all_gen_neutron_p->GetXaxis()->SetTitle("p /GeV");
  hist_FD_all_gen_neutron_p->GetYaxis()->SetTitle("counts");
  hist_FD_all_gen_neutron_theta = new TH1F("hist_FD_all_gen_neutron_theta", "neutron #Theta", 300,0,150);   
  hist_FD_all_gen_neutron_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_FD_all_gen_neutron_theta->GetYaxis()->SetTitle("counts");
  hist_FD_all_gen_pip_p = new TH1F("hist_FD_all_gen_pip_p", "#pi^{+} momentum", 125,0,10);   
  hist_FD_all_gen_pip_p->GetXaxis()->SetTitle("p /GeV");
  hist_FD_all_gen_pip_p->GetYaxis()->SetTitle("counts");
  hist_FD_all_gen_pip_theta = new TH1F("hist_FD_all_gen_pip_theta", "#pi^{+} #Theta", 300,0,150);   
  hist_FD_all_gen_pip_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_FD_all_gen_pip_theta->GetYaxis()->SetTitle("counts");
  hist_FD_all_gen_pim_p = new TH1F("hist_FD_all_gen_pim_p", "#pi^{-} momentum", 125,0,10);   
  hist_FD_all_gen_pim_p->GetXaxis()->SetTitle("p /GeV");
  hist_FD_all_gen_pim_p->GetYaxis()->SetTitle("counts");
  hist_FD_all_gen_pim_theta = new TH1F("hist_FD_all_gen_pim_theta", "#pi^{-} #Theta", 300,0,150);   
  hist_FD_all_gen_pim_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_FD_all_gen_pim_theta->GetYaxis()->SetTitle("counts");
  hist_FD_all_gen_Kp_p = new TH1F("hist_FD_all_gen_Kp_p", "K^{+} momentum", 125,0,10);   
  hist_FD_all_gen_Kp_p->GetXaxis()->SetTitle("p /GeV");
  hist_FD_all_gen_Kp_p->GetYaxis()->SetTitle("counts");
  hist_FD_all_gen_Kp_theta = new TH1F("hist_FD_all_gen_Kp_theta", "K^{+} #Theta", 300,0,150);   
  hist_FD_all_gen_Kp_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_FD_all_gen_Kp_theta->GetYaxis()->SetTitle("counts");
  hist_FD_all_gen_Km_p = new TH1F("hist_FD_all_gen_Km_p", "K^{-} momentum", 125,0,10);   
  hist_FD_all_gen_Km_p->GetXaxis()->SetTitle("p /GeV");
  hist_FD_all_gen_Km_p->GetYaxis()->SetTitle("counts");
  hist_FD_all_gen_Km_theta = new TH1F("hist_FD_all_gen_Km_theta", "K^{-} #Theta", 300,0,150);   
  hist_FD_all_gen_Km_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_FD_all_gen_Km_theta->GetYaxis()->SetTitle("counts");

  hist_CD_all_corr_electron_p = new TH1F("hist_CD_all_corr_electron_p", "electron momentum", 125,0,10);   
  hist_CD_all_corr_electron_p->GetXaxis()->SetTitle("p /GeV");
  hist_CD_all_corr_electron_p->GetYaxis()->SetTitle("counts");
  hist_CD_all_corr_electron_theta = new TH1F("hist_CD_all_corr_electron_theta", "electron #Theta", 300,0,150);   
  hist_CD_all_corr_electron_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_CD_all_corr_electron_theta->GetYaxis()->SetTitle("counts");
  hist_CD_all_corr_proton_p = new TH1F("hist_CD_all_corr_proton_p", "proton momentum", 125,0,10);   
  hist_CD_all_corr_proton_p->GetXaxis()->SetTitle("p /GeV");
  hist_CD_all_corr_proton_p->GetYaxis()->SetTitle("counts");
  hist_CD_all_corr_proton_theta = new TH1F("hist_CD_all_corr_proton_theta", "proton #Theta", 300,0,150);   
  hist_CD_all_corr_proton_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_CD_all_corr_proton_theta->GetYaxis()->SetTitle("counts");
  hist_CD_all_corr_neutron_p = new TH1F("hist_CD_all_corr_neutron_p", "neutron momentum", 125,0,10);   
  hist_CD_all_corr_neutron_p->GetXaxis()->SetTitle("p /GeV");
  hist_CD_all_corr_neutron_p->GetYaxis()->SetTitle("counts");
  hist_CD_all_corr_neutron_theta = new TH1F("hist_CD_all_corr_neutron_theta", "neutron #Theta", 300,0,150);   
  hist_CD_all_corr_neutron_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_CD_all_corr_neutron_theta->GetYaxis()->SetTitle("counts");
  hist_CD_all_corr_pip_p = new TH1F("hist_CD_all_corr_pip_p", "#pi^{+} momentum", 125,0,10);   
  hist_CD_all_corr_pip_p->GetXaxis()->SetTitle("p /GeV");
  hist_CD_all_corr_pip_p->GetYaxis()->SetTitle("counts");
  hist_CD_all_corr_pip_theta = new TH1F("hist_CD_all_corr_pip_theta", "#pi^{+} #Theta", 300,0,150);   
  hist_CD_all_corr_pip_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_CD_all_corr_pip_theta->GetYaxis()->SetTitle("counts");
  hist_CD_all_corr_pim_p = new TH1F("hist_CD_all_corr_pim_p", "#pi^{-} momentum", 125,0,10);   
  hist_CD_all_corr_pim_p->GetXaxis()->SetTitle("p /GeV");
  hist_CD_all_corr_pim_p->GetYaxis()->SetTitle("counts");
  hist_CD_all_corr_pim_theta = new TH1F("hist_CD_all_corr_pim_theta", "#pi^{-} #Theta", 300,0,150);   
  hist_CD_all_corr_pim_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_CD_all_corr_pim_theta->GetYaxis()->SetTitle("counts");
  hist_CD_all_corr_Kp_p = new TH1F("hist_CD_all_corr_Kp_p", "K^{+} momentum", 125,0,10);   
  hist_CD_all_corr_Kp_p->GetXaxis()->SetTitle("p /GeV");
  hist_CD_all_corr_Kp_p->GetYaxis()->SetTitle("counts");
  hist_CD_all_corr_Kp_theta = new TH1F("hist_CD_all_corr_Kp_theta", "K^{+} #Theta", 300,0,150);   
  hist_CD_all_corr_Kp_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_CD_all_corr_Kp_theta->GetYaxis()->SetTitle("counts");
  hist_CD_all_corr_Km_p = new TH1F("hist_CD_all_corr_Km_p", "K^{-} momentum", 125,0,10);   
  hist_CD_all_corr_Km_p->GetXaxis()->SetTitle("p /GeV");
  hist_CD_all_corr_Km_p->GetYaxis()->SetTitle("counts");
  hist_CD_all_corr_Km_theta = new TH1F("hist_CD_all_corr_Km_theta", "K^{-} #Theta", 300,0,150);   
  hist_CD_all_corr_Km_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_CD_all_corr_Km_theta->GetYaxis()->SetTitle("counts");

  hist_CD_all_det_electron_p = new TH1F("hist_CD_all_det_electron_p", "electron momentum", 125,0,10);   
  hist_CD_all_det_electron_p->GetXaxis()->SetTitle("p /GeV");
  hist_CD_all_det_electron_p->GetYaxis()->SetTitle("counts");
  hist_CD_all_det_electron_theta = new TH1F("hist_CD_all_det_electron_theta", "electron #Theta", 300,0,150);   
  hist_CD_all_det_electron_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_CD_all_det_electron_theta->GetYaxis()->SetTitle("counts");
  hist_CD_all_det_proton_p = new TH1F("hist_CD_all_det_proton_p", "proton momentum", 125,0,10);   
  hist_CD_all_det_proton_p->GetXaxis()->SetTitle("p /GeV");
  hist_CD_all_det_proton_p->GetYaxis()->SetTitle("counts");
  hist_CD_all_det_proton_theta = new TH1F("hist_CD_all_det_proton_theta", "proton #Theta", 300,0,150);   
  hist_CD_all_det_proton_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_CD_all_det_proton_theta->GetYaxis()->SetTitle("counts");
  hist_CD_all_det_neutron_p = new TH1F("hist_CD_all_det_neutron_p", "neutron momentum", 125,0,10);   
  hist_CD_all_det_neutron_p->GetXaxis()->SetTitle("p /GeV");
  hist_CD_all_det_neutron_p->GetYaxis()->SetTitle("counts");
  hist_CD_all_det_neutron_theta = new TH1F("hist_CD_all_det_neutron_theta", "neutron #Theta", 300,0,150);   
  hist_CD_all_det_neutron_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_CD_all_det_neutron_theta->GetYaxis()->SetTitle("counts");
  hist_CD_all_det_pip_p = new TH1F("hist_CD_all_det_pip_p", "#pi^{+} momentum", 125,0,10);   
  hist_CD_all_det_pip_p->GetXaxis()->SetTitle("p /GeV");
  hist_CD_all_det_pip_p->GetYaxis()->SetTitle("counts");
  hist_CD_all_det_pip_theta = new TH1F("hist_CD_all_det_pip_theta", "#pi^{+} #Theta", 300,0,150);   
  hist_CD_all_det_pip_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_CD_all_det_pip_theta->GetYaxis()->SetTitle("counts");
  hist_CD_all_det_pim_p = new TH1F("hist_CD_all_det_pim_p", "#pi^{-} momentum", 125,0,10);   
  hist_CD_all_det_pim_p->GetXaxis()->SetTitle("p /GeV");
  hist_CD_all_det_pim_p->GetYaxis()->SetTitle("counts");
  hist_CD_all_det_pim_theta = new TH1F("hist_CD_all_det_pim_theta", "#pi^{-} #Theta", 300,0,150);   
  hist_CD_all_det_pim_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_CD_all_det_pim_theta->GetYaxis()->SetTitle("counts");
  hist_CD_all_det_Kp_p = new TH1F("hist_CD_all_det_Kp_p", "K^{+} momentum", 125,0,10);   
  hist_CD_all_det_Kp_p->GetXaxis()->SetTitle("p /GeV");
  hist_CD_all_det_Kp_p->GetYaxis()->SetTitle("counts");
  hist_CD_all_det_Kp_theta = new TH1F("hist_CD_all_det_Kp_theta", "K^{+} #Theta", 300,0,150);   
  hist_CD_all_det_Kp_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_CD_all_det_Kp_theta->GetYaxis()->SetTitle("counts");
  hist_CD_all_det_Km_p = new TH1F("hist_CD_all_det_Km_p", "K^{-} momentum", 125,0,10);   
  hist_CD_all_det_Km_p->GetXaxis()->SetTitle("p /GeV");
  hist_CD_all_det_Km_p->GetYaxis()->SetTitle("counts");
  hist_CD_all_det_Km_theta = new TH1F("hist_CD_all_det_Km_theta", "K^{-} #Theta", 300,0,150);   
  hist_CD_all_det_Km_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_CD_all_det_Km_theta->GetYaxis()->SetTitle("counts");

  hist_CD_all_gen_electron_p = new TH1F("hist_CD_all_gen_electron_p", "electron momentum", 125,0,10);   
  hist_CD_all_gen_electron_p->GetXaxis()->SetTitle("p /GeV");
  hist_CD_all_gen_electron_p->GetYaxis()->SetTitle("counts");
  hist_CD_all_gen_electron_theta = new TH1F("hist_CD_all_gen_electron_theta", "electron #Theta", 300,0,150);   
  hist_CD_all_gen_electron_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_CD_all_gen_electron_theta->GetYaxis()->SetTitle("counts");
  hist_CD_all_gen_proton_p = new TH1F("hist_CD_all_gen_proton_p", "proton momentum", 125,0,10);   
  hist_CD_all_gen_proton_p->GetXaxis()->SetTitle("p /GeV");
  hist_CD_all_gen_proton_p->GetYaxis()->SetTitle("counts");
  hist_CD_all_gen_proton_theta = new TH1F("hist_CD_all_gen_proton_theta", "proton #Theta", 300,0,150);   
  hist_CD_all_gen_proton_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_CD_all_gen_proton_theta->GetYaxis()->SetTitle("counts");
  hist_CD_all_gen_neutron_p = new TH1F("hist_CD_all_gen_neutron_p", "neutron momentum", 125,0,10);   
  hist_CD_all_gen_neutron_p->GetXaxis()->SetTitle("p /GeV");
  hist_CD_all_gen_neutron_p->GetYaxis()->SetTitle("counts");
  hist_CD_all_gen_neutron_theta = new TH1F("hist_CD_all_gen_neutron_theta", "neutron #Theta", 300,0,150);   
  hist_CD_all_gen_neutron_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_CD_all_gen_neutron_theta->GetYaxis()->SetTitle("counts");
  hist_CD_all_gen_pip_p = new TH1F("hist_CD_all_gen_pip_p", "#pi^{+} momentum", 125,0,10);   
  hist_CD_all_gen_pip_p->GetXaxis()->SetTitle("p /GeV");
  hist_CD_all_gen_pip_p->GetYaxis()->SetTitle("counts");
  hist_CD_all_gen_pip_theta = new TH1F("hist_CD_all_gen_pip_theta", "#pi^{+} #Theta", 300,0,150);   
  hist_CD_all_gen_pip_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_CD_all_gen_pip_theta->GetYaxis()->SetTitle("counts");
  hist_CD_all_gen_pim_p = new TH1F("hist_CD_all_gen_pim_p", "#pi^{-} momentum", 125,0,10);   
  hist_CD_all_gen_pim_p->GetXaxis()->SetTitle("p /GeV");
  hist_CD_all_gen_pim_p->GetYaxis()->SetTitle("counts");
  hist_CD_all_gen_pim_theta = new TH1F("hist_CD_all_gen_pim_theta", "#pi^{-} #Theta", 300,0,150);   
  hist_CD_all_gen_pim_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_CD_all_gen_pim_theta->GetYaxis()->SetTitle("counts");
  hist_CD_all_gen_Kp_p = new TH1F("hist_CD_all_gen_Kp_p", "K^{+} momentum", 125,0,10);   
  hist_CD_all_gen_Kp_p->GetXaxis()->SetTitle("p /GeV");
  hist_CD_all_gen_Kp_p->GetYaxis()->SetTitle("counts");
  hist_CD_all_gen_Kp_theta = new TH1F("hist_CD_all_gen_Kp_theta", "K^{+} #Theta", 300,0,150);   
  hist_CD_all_gen_Kp_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_CD_all_gen_Kp_theta->GetYaxis()->SetTitle("counts");
  hist_CD_all_gen_Km_p = new TH1F("hist_CD_all_gen_Km_p", "K^{-} momentum", 125,0,10);   
  hist_CD_all_gen_Km_p->GetXaxis()->SetTitle("p /GeV");
  hist_CD_all_gen_Km_p->GetYaxis()->SetTitle("counts");
  hist_CD_all_gen_Km_theta = new TH1F("hist_CD_all_gen_Km_theta", "K^{-} #Theta", 300,0,150);   
  hist_CD_all_gen_Km_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_CD_all_gen_Km_theta->GetYaxis()->SetTitle("counts");



out->mkdir("particles_identified_histograms_selected");				
out->cd ("particles_identified_histograms_selected");

  TH1F *hist_electron_p; TH1F *hist_electron_theta; TH1F *hist_electron_phi; 
  TH2F *hist_electron_p_vs_theta; TH2F *hist_electron_p_vs_phi; TH2F *hist_electron_theta_vs_phi;
  TH1F *hist_proton_p; TH1F *hist_proton_theta; TH1F *hist_proton_phi; 
  TH2F *hist_proton_p_vs_theta; TH2F *hist_proton_p_vs_phi; TH2F *hist_proton_theta_vs_phi;
  TH1F *hist_neutron_p; TH1F *hist_neutron_theta; TH1F *hist_neutron_phi;
  TH2F *hist_neutron_p_vs_theta; TH2F *hist_neutron_p_vs_phi; TH2F *hist_neutron_theta_vs_phi;

  TH1F *hist_pip_p; TH1F *hist_pip_theta; TH1F *hist_pip_phi;
  TH2F *hist_pip_p_vs_theta; TH2F *hist_pip_p_vs_phi; TH2F *hist_pip_theta_vs_phi;
  TH1F *hist_pim_p; TH1F *hist_pim_theta; TH1F *hist_pim_phi;
  TH2F *hist_pim_p_vs_theta; TH2F *hist_pim_p_vs_phi;  TH2F *hist_pim_theta_vs_phi;
 
  TH1F *hist_Kp_p; TH1F *hist_Kp_theta; TH1F *hist_Kp_phi;
  TH2F *hist_Kp_p_vs_theta; TH2F *hist_Kp_p_vs_phi; TH2F *hist_Kp_theta_vs_phi;
  TH1F *hist_Km_p; TH1F *hist_Km_theta; TH1F *hist_Km_phi;
  TH2F *hist_Km_p_vs_theta; TH2F *hist_Km_p_vs_phi; TH2F *hist_Km_theta_vs_phi;

  TH1F *hist_photon_p; TH1F *hist_photon_theta; TH1F *hist_photon_phi;
  TH2F *hist_photon_p_vs_theta; TH2F *hist_photon_p_vs_phi; TH2F *hist_photon_theta_vs_phi;

  hist_electron_p = new TH1F("hist_electron_p", "electron momentum", 500,0,Ebeam+1);   
  hist_electron_p->GetXaxis()->SetTitle("p /GeV");
  hist_electron_p->GetYaxis()->SetTitle("counts");
  hist_electron_theta = new TH1F("hist_electron_theta", "electron #Theta", 200,0,50);   
  hist_electron_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_electron_theta->GetYaxis()->SetTitle("counts");
  hist_electron_phi = new TH1F("hist_electron_phi", "electron #phi", 180,-180,180);   
  hist_electron_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_electron_phi->GetYaxis()->SetTitle("counts");
  hist_electron_p_vs_theta = new TH2F("hist_electron_p_vs_theta", "electron p vs #Theta", 200,0,50,500,0,Ebeam+1);   
  hist_electron_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_electron_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_electron_p_vs_phi = new TH2F("hist_electron_p_vs_phi", "electron p vs #phi", 180,-180,180, 500,0,Ebeam+1);   
  hist_electron_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_electron_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_electron_theta_vs_phi = new TH2F("hist_electron_theta_vs_phi", "electron #Theta vs phi", 180,-180,180, 200,0,50);   
  hist_electron_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_electron_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  hist_proton_p = new TH1F("hist_proton_p", "proton momentum", 500,0,Ebeam+1);   
  hist_proton_p->GetXaxis()->SetTitle("p /GeV");
  hist_proton_p->GetYaxis()->SetTitle("counts");
  hist_proton_theta = new TH1F("hist_proton_theta", "proton #Theta", 560,0,140);   
  hist_proton_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_proton_theta->GetYaxis()->SetTitle("counts");
  hist_proton_phi = new TH1F("hist_proton_phi", "proton #phi", 180,-180,180);   
  hist_proton_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_proton_phi->GetYaxis()->SetTitle("counts");
  hist_proton_p_vs_theta = new TH2F("hist_proton_p_vs_theta", "proton p vs #Theta", 560,0,140,500,0,Ebeam+1);   
  hist_proton_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_proton_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_proton_p_vs_phi = new TH2F("hist_proton_p_vs_phi", "proton p vs #phi", 180,-180,180, 500,0,Ebeam+1);   
  hist_proton_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_proton_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_proton_theta_vs_phi = new TH2F("hist_proton_theta_vs_phi", "proton #Theta vs phi", 180,-180,180, 560,0,140);   
  hist_proton_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_proton_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  hist_neutron_p = new TH1F("hist_neutron_p", "neutron momentum", 500,0,Ebeam+1);   
  hist_neutron_p->GetXaxis()->SetTitle("p /GeV");
  hist_neutron_p->GetYaxis()->SetTitle("counts");
  hist_neutron_theta = new TH1F("hist_neutron_theta", "neutron #Theta", 560,0,140);   
  hist_neutron_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_neutron_theta->GetYaxis()->SetTitle("counts");
  hist_neutron_phi = new TH1F("hist_neutron_phi", "neutron #phi", 180,-180,180);   
  hist_neutron_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_neutron_phi->GetYaxis()->SetTitle("counts");
  hist_neutron_p_vs_theta = new TH2F("hist_neutron_p_vs_theta", "neutron p vs #Theta", 560,0,140,500,0,Ebeam+1);   
  hist_neutron_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_neutron_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_neutron_p_vs_phi = new TH2F("hist_neutron_p_vs_phi", "neutron p vs #phi", 180,-180,180, 500,0,Ebeam+1);   
  hist_neutron_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_neutron_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_neutron_theta_vs_phi = new TH2F("hist_neutron_theta_vs_phi", "neutron #Theta vs phi", 180,-180,180, 560,0,140);   
  hist_neutron_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_neutron_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  hist_pip_p = new TH1F("hist_pip_p", "#pi^{+} momentum", 500,0,Ebeam+1);   
  hist_pip_p->GetXaxis()->SetTitle("p /GeV");
  hist_pip_p->GetYaxis()->SetTitle("counts");
  hist_pip_theta = new TH1F("hist_pip_theta", "#pi^{+} #Theta", 560,0,140);   
  hist_pip_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_pip_theta->GetYaxis()->SetTitle("counts");
  hist_pip_phi = new TH1F("hist_pip_phi", "pi^{+} #phi", 180,-180,180);   
  hist_pip_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_pip_phi->GetYaxis()->SetTitle("counts");
  hist_pip_p_vs_theta = new TH2F("hist_pip_p_vs_theta", "#pi^{+} p vs #Theta", 560,0,140,500,0,Ebeam+1);   
  hist_pip_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_pip_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_pip_p_vs_phi = new TH2F("hist_pip_p_vs_phi", "#pi^{+} p vs #phi", 180,-180,180, 500,0,Ebeam+1);   
  hist_pip_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_pip_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_pip_theta_vs_phi = new TH2F("hist_pip_theta_vs_phi", "#pi^{+} #Theta vs phi", 180,-180,180, 560,0,140);   
  hist_pip_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_pip_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  hist_pim_p = new TH1F("hist_pim_p", "#pi^{-} momentum", 500,0,Ebeam+1);   
  hist_pim_p->GetXaxis()->SetTitle("p /GeV");
  hist_pim_p->GetYaxis()->SetTitle("counts");
  hist_pim_theta = new TH1F("hist_pim_theta", "#pi^{-} #Theta", 560,0,140);   
  hist_pim_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_pim_theta->GetYaxis()->SetTitle("counts");
  hist_pim_phi = new TH1F("hist_pim_phi", "#pi^{-} #phi", 180,-180,180);   
  hist_pim_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_pim_phi->GetYaxis()->SetTitle("counts");
  hist_pim_p_vs_theta = new TH2F("hist_pim_p_vs_theta", "#pi^{-} p vs #Theta", 560,0,140,500,0,Ebeam+1);   
  hist_pim_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_pim_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_pim_p_vs_phi = new TH2F("hist_pim_p_vs_phi", "#pi^{-} p vs #phi", 180,-180,180, 500,0,Ebeam+1);   
  hist_pim_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_pim_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_pim_theta_vs_phi = new TH2F("hist_pim_theta_vs_phi", "#pi^{-} #Theta vs phi", 180,-180,180, 560,0,140);   
  hist_pim_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_pim_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  hist_Kp_p = new TH1F("hist_Kp_p", "K^{+} momentum", 500,0,Ebeam+1);   
  hist_Kp_p->GetXaxis()->SetTitle("p /GeV");
  hist_Kp_p->GetYaxis()->SetTitle("counts");
  hist_Kp_theta = new TH1F("hist_Kp_theta", "K^{+} #Theta", 560,0,140);   
  hist_Kp_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_Kp_theta->GetYaxis()->SetTitle("counts");
  hist_Kp_phi = new TH1F("hist_Kp_phi", "K^{+} #phi", 180,-180,180);   
  hist_Kp_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_Kp_phi->GetYaxis()->SetTitle("counts");
  hist_Kp_p_vs_theta = new TH2F("hist_Kp_p_vs_theta", "K^{+} p vs #Theta", 560,0,140,500,0,Ebeam+1);   
  hist_Kp_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_Kp_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_Kp_p_vs_phi = new TH2F("hist_Kp_p_vs_phi", "K^{+} p vs #phi", 180,-180,180, 500,0,Ebeam+1);   
  hist_Kp_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_Kp_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_Kp_theta_vs_phi = new TH2F("hist_Kp_theta_vs_phi", "K^{+} #Theta vs phi", 180,-180,180, 560,0,140);   
  hist_Kp_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_Kp_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  hist_Km_p = new TH1F("hist_Km_p", "K^{-} momentum", 500,0,Ebeam+1);   
  hist_Km_p->GetXaxis()->SetTitle("p /GeV");
  hist_Km_p->GetYaxis()->SetTitle("counts");
  hist_Km_theta = new TH1F("hist_Km_theta", "K^{-} #Theta", 560,0,140);   
  hist_Km_theta->GetXaxis()->SetTitle("theta /deg");
  hist_Km_theta->GetYaxis()->SetTitle("counts");
  hist_Km_phi = new TH1F("hist_Km_phi", "K^{-} #phi", 180,-180,180);   
  hist_Km_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_Km_phi->GetYaxis()->SetTitle("counts");
  hist_Km_p_vs_theta = new TH2F("hist_Km_p_vs_theta", "K^{-} p vs #Theta", 560,0,140,500,0,Ebeam+1);   
  hist_Km_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_Km_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_Km_p_vs_phi = new TH2F("hist_Km_p_vs_phi", "K^{-} p vs #phi", 180,-180,180, 500,0,Ebeam+1);   
  hist_Km_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_Km_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_Km_theta_vs_phi = new TH2F("hist_Km_theta_vs_phi", "K^{-} #Theta vs phi", 180,-180,180, 560,0,140);   
  hist_Km_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_Km_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  hist_photon_p = new TH1F("hist_photon_p", "photon momentum", 500,0,Ebeam+1);   
  hist_photon_p->GetXaxis()->SetTitle("p /GeV");
  hist_photon_p->GetYaxis()->SetTitle("counts");
  hist_photon_theta = new TH1F("hist_photon_theta", "photon #Theta", 200,0,50);   
  hist_photon_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_photon_theta->GetYaxis()->SetTitle("counts");
  hist_photon_phi = new TH1F("hist_photon_phi", "photon #phi", 180,-180,180);   
  hist_photon_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_photon_phi->GetYaxis()->SetTitle("counts");
  hist_photon_p_vs_theta = new TH2F("hist_photon_p_vs_theta", "photon p vs #Theta", 200,0,50,500,0,Ebeam+1);   
  hist_photon_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_photon_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_photon_p_vs_phi = new TH2F("hist_photon_p_vs_phi", "photon p vs #phi", 180,-180,180, 500,0,Ebeam+1);   
  hist_photon_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_photon_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_photon_theta_vs_phi = new TH2F("hist_photon_theta_vs_phi", "photon #Theta vs phi", 180,-180,180, 200,0,50);   
  hist_photon_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_photon_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");


out->mkdir("vertex_plots");				
out->cd ("vertex_plots");

  TH1F *hist_positive_vertex; TH1F *hist_negative_vertex; TH1F *hist_electron_vertex;  TH1F *hist_proton_vertex;
  TH1F *hist_positive_vertex_sec[6]; TH1F *hist_negative_vertex_sec[6]; TH1F *hist_electron_vertex_sec[6];  TH1F *hist_proton_vertex_sec[6];

  TH2F *hist_positive_vertex_vs_theta; TH2F *hist_negative_vertex_vs_theta; TH2F *hist_electron_vertex_vs_theta;  TH2F *hist_proton_vertex_vs_theta;
  TH2F *hist_positive_vertex_vs_theta_sec[6]; TH2F *hist_negative_vertex_vs_theta_sec[6]; TH2F *hist_electron_vertex_vs_theta_sec[6];  TH2F *hist_proton_vertex_vs_theta_sec[6];

  TH2F *hist_positive_vertex_vs_phi; TH2F *hist_negative_vertex_vs_phi; TH2F *hist_electron_vertex_vs_phi;  TH2F *hist_proton_vertex_vs_phi;
  TH2F *hist_positive_vertex_vs_phi_sec[6]; TH2F *hist_negative_vertex_vs_phi_sec[6]; TH2F *hist_electron_vertex_vs_phi_sec[6];  TH2F *hist_proton_vertex_vs_phi_sec[6];

  TH2F *hist_positive_vertex_vs_p; TH2F *hist_negative_vertex_vs_p; TH2F *hist_neutral_vertex_vs_p;  TH2F *hist_electron_vertex_vs_p;  TH2F *hist_proton_vertex_vs_p;
  TH2F *hist_positive_vertex_vs_p_sec[6]; TH2F *hist_negative_vertex_vs_p_sec[6]; TH2F *hist_neutral_vertex_vs_p_sec[6];  
  TH2F *hist_electron_vertex_vs_p_sec[6];  TH2F *hist_proton_vertex_vs_p_sec[6];


  hist_positive_vertex = new TH1F("hist_positive_vertex", "z vertex for particles with positive charge",  1000, -100, 100);   
  hist_positive_vertex->GetXaxis()->SetTitle("v_{z} /cm"); hist_positive_vertex->GetYaxis()->SetTitle("counts");
  hist_negative_vertex = new TH1F("hist_negative_vertex", "z vertex for particles with negative charge",  1000, -100, 100);   
  hist_negative_vertex->GetXaxis()->SetTitle("v_{z} /cm"); hist_negative_vertex->GetYaxis()->SetTitle("counts");
  hist_electron_vertex = new TH1F("hist_electron_vertex", "z vertex for electrons",  1000, -100, 100);   
  hist_electron_vertex->GetXaxis()->SetTitle("v_{z} /cm"); hist_electron_vertex->GetYaxis()->SetTitle("counts");
  hist_proton_vertex = new TH1F("hist_proton_vertex", "z vertex for protons",  1000, -100, 100);   
  hist_proton_vertex->GetXaxis()->SetTitle("v_{z} /cm"); hist_proton_vertex->GetYaxis()->SetTitle("counts");

  hist_positive_vertex_vs_theta = new TH2F("hist_positive_vertex_vs_theta", "z vertex vs #Theta for particles with positive charge", 140, 0, 140, 1000, -100, 100);   
  hist_positive_vertex_vs_theta->GetXaxis()->SetTitle("#Theta /deg"); hist_positive_vertex_vs_theta->GetYaxis()->SetTitle("v_{z} /cm");
  hist_negative_vertex_vs_theta = new TH2F("hist_negative_vertex_vs_theta", "z vertex vs #Theta for particles with negative charge", 140, 0, 140,  1000, -100, 100);   
  hist_negative_vertex_vs_theta->GetXaxis()->SetTitle("#Theta /deg"); hist_negative_vertex_vs_theta->GetYaxis()->SetTitle("v_{z} /cm");
  hist_electron_vertex_vs_theta = new TH2F("hist_electron_vertex_vs_theta", "z vertex vs #Theta for electrons", 140, 0, 140,  1000, -100, 100);   
  hist_electron_vertex_vs_theta->GetXaxis()->SetTitle("#Theta /deg"); hist_electron_vertex_vs_theta->GetYaxis()->SetTitle("v_{z} /cm");
  hist_proton_vertex_vs_theta = new TH2F("hist_proton_vertex_vs_theta", "z vertex vs #Theta for protons", 140, 0, 140,  1000, -100, 100);   
  hist_proton_vertex_vs_theta->GetXaxis()->SetTitle("#Theta /deg"); hist_proton_vertex_vs_theta->GetYaxis()->SetTitle("v_{z} /cm");

  hist_positive_vertex_vs_phi = new TH2F("hist_positive_vertex_vs_phi", "z vertex vs #phi for particles with positive charge", 140, 0, 140, 1000, -100, 100);   
  hist_positive_vertex_vs_phi->GetXaxis()->SetTitle("#phi /deg"); hist_positive_vertex_vs_phi->GetYaxis()->SetTitle("v_{z} /cm");
  hist_negative_vertex_vs_phi = new TH2F("hist_negative_vertex_vs_phi", "z vertex vs #phi for particles with negative charge", 140, 0, 140,  1000, -100, 100);   
  hist_negative_vertex_vs_phi->GetXaxis()->SetTitle("#phi /deg"); hist_negative_vertex_vs_phi->GetYaxis()->SetTitle("v_{z} /cm");
  hist_electron_vertex_vs_phi = new TH2F("hist_electron_vertex_vs_phi", "z vertex vs #phi for electrons", 140, 0, 140,  1000, -100, 100);   
  hist_electron_vertex_vs_phi->GetXaxis()->SetTitle("#phi /deg"); hist_electron_vertex_vs_phi->GetYaxis()->SetTitle("v_{z} /cm");
  hist_proton_vertex_vs_phi = new TH2F("hist_proton_vertex_vs_phi", "z vertex vs #phi for protons", 140, 0, 140,  1000, -100, 100);   
  hist_proton_vertex_vs_phi->GetXaxis()->SetTitle("#phi /deg"); hist_proton_vertex_vs_phi->GetYaxis()->SetTitle("v_{z} /cm");

  hist_positive_vertex_vs_p = new TH2F("hist_positive_vertex_vs_p", "z vertex vs #p for particles with positive charge", 500, 0, Ebeam+1, 1000, -100, 100);   
  hist_positive_vertex_vs_p->GetXaxis()->SetTitle("#p /GeV"); hist_positive_vertex_vs_p->GetYaxis()->SetTitle("v_{z} /cm");
  hist_negative_vertex_vs_p = new TH2F("hist_negative_vertex_vs_p", "z vertex vs #p for particles with negative charge", 500, 0, Ebeam+1,  1000, -100, 100);   
  hist_negative_vertex_vs_p->GetXaxis()->SetTitle("#p /GeV"); hist_negative_vertex_vs_p->GetYaxis()->SetTitle("v_{z} /cm");
  hist_electron_vertex_vs_p = new TH2F("hist_electron_vertex_vs_p", "z vertex vs #p for electrons", 500, 0, Ebeam+1,  1000, -100, 100);   
  hist_electron_vertex_vs_p->GetXaxis()->SetTitle("#p /GeV"); hist_electron_vertex_vs_p->GetYaxis()->SetTitle("v_{z} /cm");
  hist_proton_vertex_vs_p = new TH2F("hist_proton_vertex_vs_p", "z vertex vs #p for protons", 500, 0, Ebeam+1,  1000, -100, 100);   
  hist_proton_vertex_vs_p->GetXaxis()->SetTitle("#p /GeV"); hist_proton_vertex_vs_p->GetYaxis()->SetTitle("v_{z} /cm");


  for(Int_t i = 0; i < 6; i++){

    sprintf(name,"hist_positive_vertex_sec%01d", i+1);  sprintf(title,"z vertex for particles with positive charge in sector %01d", i+1);
    hist_positive_vertex_sec[i] = new TH1F(name, title,  1000, -100, 100);   
    hist_positive_vertex_sec[i]->GetXaxis()->SetTitle("v_{z} /cm"); hist_positive_vertex_sec[i]->GetYaxis()->SetTitle("counts");
    sprintf(name,"hist_negative_vertex_sec%01d", i+1);  sprintf(title,"z vertex for particles with negative charge in sector %01d", i+1);
    hist_negative_vertex_sec[i] = new TH1F(name, title,  1000, -100, 100);   
    hist_negative_vertex_sec[i]->GetXaxis()->SetTitle("v_{z} /cm"); hist_negative_vertex_sec[i]->GetYaxis()->SetTitle("counts");
    sprintf(name,"hist_electron_vertex_sec%01d", i+1);  sprintf(title,"z vertex for electrons in sector %01d", i+1);
    hist_electron_vertex_sec[i] = new TH1F(name, title,  1000, -100, 100);   
    hist_electron_vertex_sec[i]->GetXaxis()->SetTitle("v_{z} /cm"); hist_electron_vertex_sec[i]->GetYaxis()->SetTitle("counts");
    sprintf(name,"hist_proton_vertex_sec%01d", i+1);  sprintf(title,"z vertex for protons in sector %01d", i+1);
    hist_proton_vertex_sec[i] = new TH1F(name, title,  1000, -100, 100);   
    hist_proton_vertex_sec[i]->GetXaxis()->SetTitle("v_{z} /cm"); hist_proton_vertex_sec[i]->GetYaxis()->SetTitle("counts");

    sprintf(name,"hist_positive_vertex_vs_theta_sec%01d", i+1);  sprintf(title,"z vertex vs #Theta for particles with positive charge in sector %01d", i+1);
    hist_positive_vertex_vs_theta_sec[i] = new TH2F(name, title, 140, 0, 140, 1000, -100, 100);   
    hist_positive_vertex_vs_theta_sec[i]->GetXaxis()->SetTitle("#Theta /deg"); hist_positive_vertex_vs_theta_sec[i]->GetYaxis()->SetTitle("v_{z} /cm");
    sprintf(name,"hist_negative_vertex_vs_theta_sec%01d", i+1);  sprintf(title,"z vertex vs #Theta for particles with negative charge in sector %01d", i+1);
    hist_negative_vertex_vs_theta_sec[i] = new TH2F(name, title, 140, 0, 140,  1000, -100, 100);   
    hist_negative_vertex_vs_theta_sec[i]->GetXaxis()->SetTitle("#Theta /deg"); hist_negative_vertex_vs_theta_sec[i]->GetYaxis()->SetTitle("v_{z} /cm");
    sprintf(name,"hist_electrons_vertex_vs_theta_sec%01d", i+1);  sprintf(title,"z vertex vs #Theta for electrons in sector %01d", i+1);
    hist_electron_vertex_vs_theta_sec[i] = new TH2F(name, title, 140, 0, 140,  1000, -100, 100);   
    hist_electron_vertex_vs_theta_sec[i]->GetXaxis()->SetTitle("#Theta /deg"); hist_electron_vertex_vs_theta_sec[i]->GetYaxis()->SetTitle("v_{z} /cm");
    sprintf(name,"hist_protons_vertex_vs_theta_sec%01d", i+1);  sprintf(title,"z vertex vs #Theta for protons in sector %01d", i+1);
    hist_proton_vertex_vs_theta_sec[i] = new TH2F(name, title, 140, 0, 140,  1000, -100, 100);   
    hist_proton_vertex_vs_theta_sec[i]->GetXaxis()->SetTitle("#Theta /deg"); hist_proton_vertex_vs_theta_sec[i]->GetYaxis()->SetTitle("v_{z} /cm");

    sprintf(name,"hist_positive_vertex_vs_phi_sec%01d", i+1);  sprintf(title,"z vertex vs #phi for particles with positive charge in sector %01d", i+1);
    hist_positive_vertex_vs_phi_sec[i] = new TH2F(name, title, 180, -180, 180, 1000, -100, 100);   
    hist_positive_vertex_vs_phi_sec[i]->GetXaxis()->SetTitle("#phi /deg"); hist_positive_vertex_vs_phi_sec[i]->GetYaxis()->SetTitle("v_{z} /cm");
    sprintf(name,"hist_negative_vertex_vs_phi_sec%01d", i+1);  sprintf(title,"z vertex vs #phi for particles with negative charge in sector %01d", i+1);
    hist_negative_vertex_vs_phi_sec[i] = new TH2F(name, title, 180, -180, 180,  1000, -100, 100);   
    hist_negative_vertex_vs_phi_sec[i]->GetXaxis()->SetTitle("#phi /deg"); hist_negative_vertex_vs_phi_sec[i]->GetYaxis()->SetTitle("v_{z} /cm");
    sprintf(name,"hist_electrons_vertex_vs_phi_sec%01d", i+1);  sprintf(title,"z vertex vs #phi for electrons in sector %01d", i+1);
    hist_electron_vertex_vs_phi_sec[i] = new TH2F(name, title, 180, -180, 180,  1000, -100, 100);   
    hist_electron_vertex_vs_phi_sec[i]->GetXaxis()->SetTitle("#phi /deg"); hist_electron_vertex_vs_phi_sec[i]->GetYaxis()->SetTitle("v_{z} /cm");
    sprintf(name,"hist_protons_vertex_vs_phi_sec%01d", i+1);  sprintf(title,"z vertex vs #phi for protons in sector %01d", i+1);
    hist_proton_vertex_vs_phi_sec[i] = new TH2F(name, title, 180, -180, 180,  1000, -100, 100);   
    hist_proton_vertex_vs_phi_sec[i]->GetXaxis()->SetTitle("#phi /deg"); hist_proton_vertex_vs_phi_sec[i]->GetYaxis()->SetTitle("v_{z} /cm");

    sprintf(name,"hist_positive_vertex_vs_p_sec%01d", i+1);  sprintf(title,"z vertex vs p for particles with positive charge in sector %01d", i+1);
    hist_positive_vertex_vs_p_sec[i] = new TH2F(name, title, 500, 0, Ebeam+1, 1000, -100, 100);   
    hist_positive_vertex_vs_p_sec[i]->GetXaxis()->SetTitle("p /GeV"); hist_positive_vertex_vs_p_sec[i]->GetYaxis()->SetTitle("v_{z} /cm");
    sprintf(name,"hist_negative_vertex_vs_p_sec%01d", i+1);  sprintf(title,"z vertex vs p for particles with negative charge in sector %01d", i+1);
    hist_negative_vertex_vs_p_sec[i] = new TH2F(name, title, 500, 0, Ebeam+1,  1000, -100, 100);   
    hist_negative_vertex_vs_p_sec[i]->GetXaxis()->SetTitle("p /GeV"); hist_negative_vertex_vs_p_sec[i]->GetYaxis()->SetTitle("v_{z} /cm");
    sprintf(name,"hist_electrons_vertex_vs_p_sec%01d", i+1);  sprintf(title,"z vertex vs p for electrons in sector %01d", i+1);
    hist_electron_vertex_vs_p_sec[i] = new TH2F(name, title, 500, 0, Ebeam+1,  1000, -100, 100);   
    hist_electron_vertex_vs_p_sec[i]->GetXaxis()->SetTitle("p /GeV"); hist_electron_vertex_vs_p_sec[i]->GetYaxis()->SetTitle("v_{z} /cm");
    sprintf(name,"hist_protons_vertex_vs_p_sec%01d", i+1);  sprintf(title,"z vertex vs p for protons in sector %01d", i+1);
    hist_proton_vertex_vs_p_sec[i] = new TH2F(name, title, 500, 0, Ebeam+1,  1000, -100, 100);   
    hist_proton_vertex_vs_p_sec[i]->GetXaxis()->SetTitle("p /GeV"); hist_proton_vertex_vs_p_sec[i]->GetYaxis()->SetTitle("v_{z} /cm");
  
  }


out->mkdir("neutrals");				
out->cd ("neutrals");

  TH1F *hist_neutral_mass;
  TH1F *hist_neutral_mass2;

  TH1F *hist_pi0_p;
  TH1F *hist_pi0_theta;
  TH1F *hist_pi0_phi;
  TH1F *hist_pi0_mass;
  TH1F *hist_pi0_mass2;

  TH1F *hist_eta_p;
  TH1F *hist_eta_theta;
  TH1F *hist_eta_phi;
  TH1F *hist_eta_mass;
  TH1F *hist_eta_mass2;

  TH2F *hist_gg_openingangle_vs_Egg;


  hist_neutral_mass = new TH1F("hist_neutral_mass", "neutral mass", 525, -0.05, 1.0);   
  hist_neutral_mass->GetXaxis()->SetTitle("mass /GeV");
  hist_neutral_mass->GetYaxis()->SetTitle("counts");
  hist_neutral_mass2 = new TH1F("hist_neutral_mass2", "neutral mass2", 850, -0.05, 0.8);   
  hist_neutral_mass2->GetXaxis()->SetTitle("mass2 /GeV");
  hist_neutral_mass2->GetYaxis()->SetTitle("counts");

  hist_pi0_p = new TH1F("hist_pi0_p", "pi0 momentum", 500,0,Ebeam+1);   
  hist_pi0_p->GetXaxis()->SetTitle("p /GeV");
  hist_pi0_p->GetYaxis()->SetTitle("counts");
  hist_pi0_theta = new TH1F("hist_pi0_theta", "pi0 theta", 60,0,60);   
  hist_pi0_theta->GetXaxis()->SetTitle("theta /deg");
  hist_pi0_theta->GetYaxis()->SetTitle("counts");
  hist_pi0_phi = new TH1F("hist_pi0_phi", "pi0 phi", 180,-180,180);   
  hist_pi0_phi->GetXaxis()->SetTitle("phi /deg");
  hist_pi0_phi->GetYaxis()->SetTitle("counts");
  hist_pi0_mass = new TH1F("hist_pi0_mass", "pi0 mass", 225, -0.05, 0.5);   
  hist_pi0_mass->GetXaxis()->SetTitle("mass /GeV");
  hist_pi0_mass->GetYaxis()->SetTitle("counts");
  hist_pi0_mass2 = new TH1F("hist_pi0_mass2", "pi0 mass2", 850, -0.05, 0.2);   
  hist_pi0_mass2->GetXaxis()->SetTitle("mass2 /GeV");
  hist_pi0_mass2->GetYaxis()->SetTitle("counts");
 
  hist_eta_p = new TH1F("hist_eta_p", "eta momentum", 500,0,Ebeam+1);   
  hist_eta_p->GetXaxis()->SetTitle("p /GeV");
  hist_eta_p->GetYaxis()->SetTitle("counts");
  hist_eta_theta = new TH1F("hist_eta_theta", "eta theta", 60,0,60);   
  hist_eta_theta->GetXaxis()->SetTitle("theta /deg");
  hist_eta_theta->GetYaxis()->SetTitle("counts");
  hist_eta_phi = new TH1F("hist_eta_phi", "eta phi", 180,-180,180);   
  hist_eta_phi->GetXaxis()->SetTitle("phi /deg");
  hist_eta_phi->GetYaxis()->SetTitle("counts");
  hist_eta_mass = new TH1F("hist_eta_mass", "eta mass", 525, -0.05, 1.0);   
  hist_eta_mass->GetXaxis()->SetTitle("mass /GeV");
  hist_eta_mass->GetYaxis()->SetTitle("counts");
  hist_eta_mass2 = new TH1F("hist_eta_mass2", "eta mass2", 650, -0.05, 0.6);   
  hist_eta_mass2->GetXaxis()->SetTitle("mass2 /GeV");
  hist_eta_mass2->GetYaxis()->SetTitle("counts");

  hist_gg_openingangle_vs_Egg = new TH2F("hist_hist_gg_openingangle_vs_Egg", "#gamma#gamma openingangle vs E_{#gamma#gamma}", 500, 0, Ebeam+0.3, 60, 0, 30);   
  hist_gg_openingangle_vs_Egg->GetXaxis()->SetTitle("E_{#gamma#gamma} /GeV");
  hist_gg_openingangle_vs_Egg->GetYaxis()->SetTitle("#alpha_{#gamma#gamma}");


out->mkdir("2photon_inv_mass");				
out->cd ("2photon_inv_mass");

  TH1F *hist_2photon_mass_11; TH1F *hist_2photon_mass_12; TH1F *hist_2photon_mass_13; TH1F *hist_2photon_mass_14; TH1F *hist_2photon_mass_15; TH1F *hist_2photon_mass_16;
  TH1F *hist_2photon_mass_22; TH1F *hist_2photon_mass_23; TH1F *hist_2photon_mass_24; TH1F *hist_2photon_mass_25; TH1F *hist_2photon_mass_26; 
  TH1F *hist_2photon_mass_33; TH1F *hist_2photon_mass_34; TH1F *hist_2photon_mass_35; TH1F *hist_2photon_mass_36;
  TH1F *hist_2photon_mass_44; TH1F *hist_2photon_mass_45; TH1F *hist_2photon_mass_46;
  TH1F *hist_2photon_mass_55; TH1F *hist_2photon_mass_56;
  TH1F *hist_2photon_mass_66;

  TH1F *hist_2photon_mass2_11; TH1F *hist_2photon_mass2_12; TH1F *hist_2photon_mass2_13; TH1F *hist_2photon_mass2_14; TH1F *hist_2photon_mass2_15; TH1F *hist_2photon_mass2_16;
  TH1F *hist_2photon_mass2_22; TH1F *hist_2photon_mass2_23; TH1F *hist_2photon_mass2_24; TH1F *hist_2photon_mass2_25; TH1F *hist_2photon_mass2_26;
  TH1F *hist_2photon_mass2_33; TH1F *hist_2photon_mass2_34; TH1F *hist_2photon_mass2_35; TH1F *hist_2photon_mass2_36;
  TH1F *hist_2photon_mass2_44; TH1F *hist_2photon_mass2_45; TH1F *hist_2photon_mass2_46;
  TH1F *hist_2photon_mass2_55; TH1F *hist_2photon_mass2_56;
  TH1F *hist_2photon_mass2_66;

  hist_2photon_mass_11 = new TH1F("2photon_mass_11", "2photon_mass_11", 50, 0, 0.4);   
  hist_2photon_mass_11->GetXaxis()->SetTitle("mass /GeV");
  hist_2photon_mass_11->GetYaxis()->SetTitle("counts");
  hist_2photon_mass_12 = new TH1F("2photon_mass_12", "2photon_mass_12", 50, 0, 0.4);   
  hist_2photon_mass_12->GetXaxis()->SetTitle("mass /GeV");
  hist_2photon_mass_12->GetYaxis()->SetTitle("counts");
  hist_2photon_mass_13 = new TH1F("2photon_mass_13", "2photon_mass_13", 50, 0, 0.4);   
  hist_2photon_mass_13->GetXaxis()->SetTitle("mass /GeV");
  hist_2photon_mass_13->GetYaxis()->SetTitle("counts");
  hist_2photon_mass_14 = new TH1F("2photon_mass_14", "2photon_mass_14", 50, 0, 0.4);   
  hist_2photon_mass_14->GetXaxis()->SetTitle("mass /GeV");
  hist_2photon_mass_14->GetYaxis()->SetTitle("counts");
  hist_2photon_mass_15 = new TH1F("2photon_mass_15", "2photon_mass_15", 50, 0, 0.4);   
  hist_2photon_mass_15->GetXaxis()->SetTitle("mass /GeV");
  hist_2photon_mass_15->GetYaxis()->SetTitle("counts");
  hist_2photon_mass_16 = new TH1F("2photon_mass_16", "2photon_mass_16", 50, 0, 0.4);   
  hist_2photon_mass_16->GetXaxis()->SetTitle("mass /GeV");
  hist_2photon_mass_16->GetYaxis()->SetTitle("counts");
  hist_2photon_mass_22 = new TH1F("2photon_mass_22", "2photon_mass_22", 50, 0, 0.4);   
  hist_2photon_mass_22->GetXaxis()->SetTitle("mass /GeV");
  hist_2photon_mass_22->GetYaxis()->SetTitle("counts");
  hist_2photon_mass_23 = new TH1F("2photon_mass_23", "2photon_mass_23", 50, 0, 0.4);   
  hist_2photon_mass_23->GetXaxis()->SetTitle("mass /GeV");
  hist_2photon_mass_23->GetYaxis()->SetTitle("counts");
  hist_2photon_mass_24 = new TH1F("2photon_mass_24", "2photon_mass_24", 50, 0, 0.4);   
  hist_2photon_mass_24->GetXaxis()->SetTitle("mass /GeV");
  hist_2photon_mass_24->GetYaxis()->SetTitle("counts");
  hist_2photon_mass_25 = new TH1F("2photon_mass_25", "2photon_mass_25", 50, 0, 0.4);   
  hist_2photon_mass_25->GetXaxis()->SetTitle("mass /GeV");
  hist_2photon_mass_25->GetYaxis()->SetTitle("counts");
  hist_2photon_mass_26 = new TH1F("2photon_mass_26", "2photon_mass_26", 50, 0, 0.4);   
  hist_2photon_mass_26->GetXaxis()->SetTitle("mass /GeV");
  hist_2photon_mass_26->GetYaxis()->SetTitle("counts");
  hist_2photon_mass_33 = new TH1F("2photon_mass_33", "2photon_mass_33", 50, 0, 0.4);   
  hist_2photon_mass_33->GetXaxis()->SetTitle("mass /GeV");
  hist_2photon_mass_33->GetYaxis()->SetTitle("counts");
  hist_2photon_mass_34 = new TH1F("2photon_mass_34", "2photon_mass_34", 50, 0, 0.4);   
  hist_2photon_mass_34->GetXaxis()->SetTitle("mass /GeV");
  hist_2photon_mass_34->GetYaxis()->SetTitle("counts");
  hist_2photon_mass_35 = new TH1F("2photon_mass_35", "2photon_mass_35", 50, 0, 0.4);   
  hist_2photon_mass_35->GetXaxis()->SetTitle("mass /GeV");
  hist_2photon_mass_35->GetYaxis()->SetTitle("counts");
  hist_2photon_mass_36 = new TH1F("2photon_mass_36", "2photon_mass_36", 50, 0, 0.4);   
  hist_2photon_mass_36->GetXaxis()->SetTitle("mass /GeV");
  hist_2photon_mass_36->GetYaxis()->SetTitle("counts");
  hist_2photon_mass_44 = new TH1F("2photon_mass_44", "2photon_mass_44", 50, 0, 0.4);   
  hist_2photon_mass_44->GetXaxis()->SetTitle("mass /GeV");
  hist_2photon_mass_44->GetYaxis()->SetTitle("counts");
  hist_2photon_mass_45 = new TH1F("2photon_mass_45", "2photon_mass_45", 50, 0, 0.4);   
  hist_2photon_mass_45->GetXaxis()->SetTitle("mass /GeV");
  hist_2photon_mass_45->GetYaxis()->SetTitle("counts");
  hist_2photon_mass_46 = new TH1F("2photon_mass_46", "2photon_mass_46", 50, 0, 0.4);   
  hist_2photon_mass_46->GetXaxis()->SetTitle("mass /GeV");
  hist_2photon_mass_46->GetYaxis()->SetTitle("counts");
  hist_2photon_mass_55 = new TH1F("2photon_mass_55", "2photon_mass_55", 50, 0, 0.4);   
  hist_2photon_mass_55->GetXaxis()->SetTitle("mass /GeV");
  hist_2photon_mass_55->GetYaxis()->SetTitle("counts");
  hist_2photon_mass_56 = new TH1F("2photon_mass_56", "2photon_mass_56", 50, 0, 0.4);   
  hist_2photon_mass_56->GetXaxis()->SetTitle("mass /GeV");
  hist_2photon_mass_56->GetYaxis()->SetTitle("counts");
  hist_2photon_mass_66 = new TH1F("2photon_mass_66", "2photon_mass_66", 50, 0, 0.4);   
  hist_2photon_mass_66->GetXaxis()->SetTitle("mass /GeV");
  hist_2photon_mass_66->GetYaxis()->SetTitle("counts");

  hist_2photon_mass2_11 = new TH1F("2photon_mass2_11", "2photon_mass2_11", 250, -0.05, 0.2);   
  hist_2photon_mass2_11->GetXaxis()->SetTitle("mass2 /GeV");
  hist_2photon_mass2_11->GetYaxis()->SetTitle("counts");
  hist_2photon_mass2_12 = new TH1F("2photon_mass2_12", "2photon_mass2_12", 250, -0.05, 0.2);   
  hist_2photon_mass2_12->GetXaxis()->SetTitle("mass2 /GeV");
  hist_2photon_mass2_12->GetYaxis()->SetTitle("counts");
  hist_2photon_mass2_13 = new TH1F("2photon_mass2_13", "2photon_mass2_13", 250, -0.05, 0.2);   
  hist_2photon_mass2_13->GetXaxis()->SetTitle("mass2 /GeV");
  hist_2photon_mass2_13->GetYaxis()->SetTitle("counts");
  hist_2photon_mass2_14 = new TH1F("2photon_mass2_14", "2photon_mass2_14", 250, -0.05, 0.2);   
  hist_2photon_mass2_14->GetXaxis()->SetTitle("mass2 /GeV");
  hist_2photon_mass2_14->GetYaxis()->SetTitle("counts");
  hist_2photon_mass2_15 = new TH1F("2photon_mass2_15", "2photon_mass2_15", 250, -0.05, 0.2);   
  hist_2photon_mass2_15->GetXaxis()->SetTitle("mass2 /GeV");
  hist_2photon_mass2_15->GetYaxis()->SetTitle("counts");
  hist_2photon_mass2_16 = new TH1F("2photon_mass2_16", "2photon_mass2_16", 250, -0.05, 0.2);   
  hist_2photon_mass2_16->GetXaxis()->SetTitle("mass2 /GeV");
  hist_2photon_mass2_16->GetYaxis()->SetTitle("counts");
  hist_2photon_mass2_22 = new TH1F("2photon_mass2_22", "2photon_mass2_22", 250, -0.05, 0.2);   
  hist_2photon_mass2_22->GetXaxis()->SetTitle("mass2 /GeV");
  hist_2photon_mass2_22->GetYaxis()->SetTitle("counts");
  hist_2photon_mass2_23 = new TH1F("2photon_mass2_23", "2photon_mass2_23", 250, -0.05, 0.2);   
  hist_2photon_mass2_23->GetXaxis()->SetTitle("mass2 /GeV");
  hist_2photon_mass2_23->GetYaxis()->SetTitle("counts");
  hist_2photon_mass2_24 = new TH1F("2photon_mass2_24", "2photon_mass2_24", 250, -0.05, 0.2);   
  hist_2photon_mass2_24->GetXaxis()->SetTitle("mass2 /GeV");
  hist_2photon_mass2_24->GetYaxis()->SetTitle("counts");
  hist_2photon_mass2_25 = new TH1F("2photon_mass2_25", "2photon_mass2_25", 250, -0.05, 0.2);   
  hist_2photon_mass2_25->GetXaxis()->SetTitle("mass2 /GeV");
  hist_2photon_mass2_25->GetYaxis()->SetTitle("counts");
  hist_2photon_mass2_26 = new TH1F("2photon_mass2_26", "2photon_mass2_26", 250, -0.05, 0.2);   
  hist_2photon_mass2_26->GetXaxis()->SetTitle("mass2 /GeV");
  hist_2photon_mass2_26->GetYaxis()->SetTitle("counts");
  hist_2photon_mass2_33 = new TH1F("2photon_mass2_33", "2photon_mass2_33", 250, -0.05, 0.2);   
  hist_2photon_mass2_33->GetXaxis()->SetTitle("mass2 /GeV");
  hist_2photon_mass2_33->GetYaxis()->SetTitle("counts");
  hist_2photon_mass2_34 = new TH1F("2photon_mass2_34", "2photon_mass2_34", 250, -0.05, 0.2);   
  hist_2photon_mass2_34->GetXaxis()->SetTitle("mass2 /GeV");
  hist_2photon_mass2_34->GetYaxis()->SetTitle("counts");
  hist_2photon_mass2_35 = new TH1F("2photon_mass2_35", "2photon_mass2_35", 250, -0.05, 0.2);   
  hist_2photon_mass2_35->GetXaxis()->SetTitle("mass2 /GeV");
  hist_2photon_mass2_35->GetYaxis()->SetTitle("counts");
  hist_2photon_mass2_36 = new TH1F("2photon_mass2_36", "2photon_mass2_36", 250, -0.05, 0.2);   
  hist_2photon_mass2_36->GetXaxis()->SetTitle("mass2 /GeV");
  hist_2photon_mass2_36->GetYaxis()->SetTitle("counts");
  hist_2photon_mass2_44 = new TH1F("2photon_mass2_44", "2photon_mass2_44", 250, -0.05, 0.2);   
  hist_2photon_mass2_44->GetXaxis()->SetTitle("mass2 /GeV");
  hist_2photon_mass2_44->GetYaxis()->SetTitle("counts");
  hist_2photon_mass2_45 = new TH1F("2photon_mass2_45", "2photon_mass2_45", 250, -0.05, 0.2);   
  hist_2photon_mass2_45->GetXaxis()->SetTitle("mass2 /GeV");
  hist_2photon_mass2_45->GetYaxis()->SetTitle("counts");
  hist_2photon_mass2_46 = new TH1F("2photon_mass2_46", "2photon_mass2_46", 250, -0.05, 0.2);   
  hist_2photon_mass2_46->GetXaxis()->SetTitle("mass2 /GeV");
  hist_2photon_mass2_46->GetYaxis()->SetTitle("counts");
  hist_2photon_mass2_55 = new TH1F("2photon_mass2_55", "2photon_mass2_55", 250, -0.05, 0.2);   
  hist_2photon_mass2_55->GetXaxis()->SetTitle("mass2 /GeV");
  hist_2photon_mass2_55->GetYaxis()->SetTitle("counts");
  hist_2photon_mass2_56 = new TH1F("2photon_mass2_56", "2photon_mass2_56", 250, -0.05, 0.2);   
  hist_2photon_mass2_56->GetXaxis()->SetTitle("mass2 /GeV");
  hist_2photon_mass2_56->GetYaxis()->SetTitle("counts");
  hist_2photon_mass2_66 = new TH1F("2photon_mass2_66", "2photon_mass2_66", 250, -0.05, 0.2);   
  hist_2photon_mass2_66->GetXaxis()->SetTitle("mass2 /GeV");
  hist_2photon_mass2_66->GetYaxis()->SetTitle("counts");


out->mkdir("kinematics");				
out->cd ("kinematics");

  TH1F *hist_W;
  TH1F *hist_W_sec_01; 
  TH1F *hist_W_sec_02; 
  TH1F *hist_W_sec_03; 
  TH1F *hist_W_sec_04; 
  TH1F *hist_W_sec_05; 
  TH1F *hist_W_sec_06;
  TH1F *hist_Q2;
  TH2F *hist_Q2_vs_W;
  TH2F *hist_Q2_vs_W_sec_01;
  TH2F *hist_Q2_vs_W_sec_02;
  TH2F *hist_Q2_vs_W_sec_03;
  TH2F *hist_Q2_vs_W_sec_04;
  TH2F *hist_Q2_vs_W_sec_05;
  TH2F *hist_Q2_vs_W_sec_06;
  TH2F *hist_W_vs_phi;
  TH1F *hist_x;
  TH1F *hist_y;
  TH1F *hist_nu;
  TH1F *hist_t_pi0;
  TH1F *hist_t_pip;
  TH1F *hist_t_pim;
  TH1F *hist_cmphi_p;
  TH1F *hist_cmcostheta_p;
  TH1F *hist_cmphi_pi0;
  TH1F *hist_cmcostheta_pi0;
  TH1F *hist_cmphi_pip;
  TH1F *hist_cmcostheta_pip;
  TH1F *hist_cmphi_pim;
  TH1F *hist_cmcostheta_pim;

  hist_W = new TH1F("W", "W all sectors", 500, 0.5, Ebeam+0.3);   
  hist_W->GetXaxis()->SetTitle("W /GeV");
  hist_W->GetYaxis()->SetTitle("counts");

  hist_W_sec_01 = new TH1F("W_sec_01", "W sector 1", 500, 0.5, Ebeam+0.3);   
  hist_W_sec_01->GetXaxis()->SetTitle("W");
  hist_W_sec_01->GetYaxis()->SetTitle("counts");
  hist_W_sec_02 = new TH1F("W_sec_02", "W sector 2", 500, 0.5, Ebeam+0.3);   
  hist_W_sec_02->GetXaxis()->SetTitle("W");
  hist_W_sec_02->GetYaxis()->SetTitle("counts");
  hist_W_sec_03 = new TH1F("W_sec_03", "W sector 3", 500, 0.5, Ebeam+0.3);   
  hist_W_sec_03->GetXaxis()->SetTitle("W");
  hist_W_sec_03->GetYaxis()->SetTitle("counts");
  hist_W_sec_04 = new TH1F("W_sec_04", "W sector 4", 500, 0.5, Ebeam+0.3);   
  hist_W_sec_04->GetXaxis()->SetTitle("W");
  hist_W_sec_04->GetYaxis()->SetTitle("counts");
  hist_W_sec_05 = new TH1F("W_sec_05", "W sector 5", 500, 0.5, Ebeam+0.3);   
  hist_W_sec_05->GetXaxis()->SetTitle("W");
  hist_W_sec_05->GetYaxis()->SetTitle("counts");
  hist_W_sec_06 = new TH1F("W_sec_06", "W sector 6", 500, 0.5, Ebeam+0.3);   
  hist_W_sec_06->GetXaxis()->SetTitle("W");
  hist_W_sec_06->GetYaxis()->SetTitle("counts");

  hist_Q2 = new TH1F("Q2", "Q^{2} all sectors", 500, 0, 2 * Ebeam);   
  hist_Q2->GetXaxis()->SetTitle("Q^{2} /GeV^{2}");
  hist_Q2->GetYaxis()->SetTitle("counts");

  hist_Q2_vs_W = new TH2F("Q2_vs_W", "Q^{2} vs W all sectors", 500, 0.5, Ebeam+0.3, 500, 0, 2 * Ebeam);   
  hist_Q2_vs_W->GetXaxis()->SetTitle("W /GeV");
  hist_Q2_vs_W->GetYaxis()->SetTitle("Q^{2}");

  hist_Q2_vs_W_sec_01 = new TH2F("Q2_vs_W_sec_01", "Q^{2} vs W sector 1", 500, 0.5, Ebeam+0.3, 500, 0, 2 * Ebeam);   
  hist_Q2_vs_W_sec_01->GetXaxis()->SetTitle("W");
  hist_Q2_vs_W_sec_01->GetYaxis()->SetTitle("Q^{2}");
  hist_Q2_vs_W_sec_02 = new TH2F("Q2_vs_W_sec_02", "Q^{2} vs W sector 2", 500, 0.5, Ebeam+0.3, 500, 0, 2 * Ebeam);   
  hist_Q2_vs_W_sec_02->GetXaxis()->SetTitle("W");
  hist_Q2_vs_W_sec_02->GetYaxis()->SetTitle("Q^{2}");
  hist_Q2_vs_W_sec_03 = new TH2F("Q2_vs_W_sec_03", "Q^{2} vs W sector 3", 500, 0.5, Ebeam+0.3, 500, 0, 2 * Ebeam);   
  hist_Q2_vs_W_sec_03->GetXaxis()->SetTitle("W");
  hist_Q2_vs_W_sec_03->GetYaxis()->SetTitle("Q^{2}");
  hist_Q2_vs_W_sec_04 = new TH2F("Q2_vs_W_sec_04", "Q^{2} vs W sector 4", 500, 0.5, Ebeam+0.3, 500, 0, 2 * Ebeam);   
  hist_Q2_vs_W_sec_04->GetXaxis()->SetTitle("W");
  hist_Q2_vs_W_sec_04->GetYaxis()->SetTitle("Q^{2}");
  hist_Q2_vs_W_sec_05 = new TH2F("Q2_vs_W_sec_05", "Q^{2} vs W sector 5", 500, 0.5, Ebeam+0.3, 500, 0, 2 * Ebeam);   
  hist_Q2_vs_W_sec_05->GetXaxis()->SetTitle("W");
  hist_Q2_vs_W_sec_05->GetYaxis()->SetTitle("Q^{2}");
  hist_Q2_vs_W_sec_06 = new TH2F("Q2_vs_W_sec_06", "Q^{2} vs W sector 6", 500, 0.5, Ebeam+0.3, 500, 0, 2 * Ebeam);   
  hist_Q2_vs_W_sec_06->GetXaxis()->SetTitle("W");
  hist_Q2_vs_W_sec_06->GetYaxis()->SetTitle("Q^{2}");

  hist_W_vs_phi = new TH2F("W_vs_phi", "W vs phi all sectors", 180, -180, 180, 500, 0.5, Ebeam+0.3);   
  hist_W_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_W_vs_phi->GetYaxis()->SetTitle("W /GeV");

  hist_t_pi0 = new TH1F("t_pi0", "t #pi^{0} all sectors", 500, 0, 2*Ebeam);   
  hist_t_pi0->GetXaxis()->SetTitle("t / GeV^{2}");
  hist_t_pi0->GetYaxis()->SetTitle("counts");

  hist_t_pip = new TH1F("t_pip", "t #pi^{+} all sectors", 500, 0, 2*Ebeam);   
  hist_t_pip->GetXaxis()->SetTitle("t / GeV^{2}");
  hist_t_pip->GetYaxis()->SetTitle("counts");

  hist_t_pim = new TH1F("t_pim", "t #pi^{-} all sectors", 500, 0, 2*Ebeam);   
  hist_t_pim->GetXaxis()->SetTitle("t / GeV^{2}");
  hist_t_pim->GetYaxis()->SetTitle("counts");

  hist_x = new TH1F("x", "x all sectors", 250, -0.1, 1.25);   
  hist_x->GetXaxis()->SetTitle("x");
  hist_x->GetYaxis()->SetTitle("counts");

  hist_y = new TH1F("y", "y all sectors", 250, - 0.1, 1.1);   
  hist_y->GetXaxis()->SetTitle("y");
  hist_y->GetYaxis()->SetTitle("counts");

  hist_nu = new TH1F("nu", "#nu all sectors", 500, -0.1, Ebeam+0.3);   
  hist_nu->GetXaxis()->SetTitle("#nu /GeV");
  hist_nu->GetYaxis()->SetTitle("counts");

  hist_cmphi_p = new TH1F("cmphi_p", "#phi_{CM} for p of e p --> e p #pi^{0}", 90, -180, 180);   
  hist_cmphi_p->GetXaxis()->SetTitle("#phi_{CM} /deg");
  hist_cmphi_p->GetYaxis()->SetTitle("counts");
  hist_cmcostheta_p = new TH1F("cmcostheta_p", "cos(#Theta_{CM})", 100, -1, 1);   
  hist_cmcostheta_p->GetXaxis()->SetTitle("cos(#Theta_{CM})");
  hist_cmcostheta_p->GetYaxis()->SetTitle("counts");

  hist_cmphi_pi0 = new TH1F("cmphi_pi0", "#phi_{CM} for #pi^{0} of e p --> e p #pi^{0}", 90, -180, 180);   
  hist_cmphi_pi0->GetXaxis()->SetTitle("#phi_{CM} /deg");
  hist_cmphi_pi0->GetYaxis()->SetTitle("counts");
  hist_cmcostheta_pi0 = new TH1F("cmcostheta_pi0", "cos(#Theta_{CM})", 100, -1, 1);   
  hist_cmcostheta_pi0->GetXaxis()->SetTitle("cos(#Theta_{CM})");
  hist_cmcostheta_pi0->GetYaxis()->SetTitle("counts");

  hist_cmphi_pip = new TH1F("cmphi_pip", "#phi_{CM} for #pi^{+} of e p --> e p #pi^{+}", 90, -180, 180);   
  hist_cmphi_pip->GetXaxis()->SetTitle("#phi_{CM} /deg");
  hist_cmphi_pip->GetYaxis()->SetTitle("counts");
  hist_cmcostheta_pip = new TH1F("cmcostheta_pip", "cos(#Theta_{CM})", 100, -1, 1);   
  hist_cmcostheta_pip->GetXaxis()->SetTitle("cos(#Theta_{CM})");
  hist_cmcostheta_pip->GetYaxis()->SetTitle("counts");

  hist_cmphi_pim = new TH1F("cmphi_pim", "#phi_{CM} for #pi^{-} of e p --> e p #pi^{-}", 90, -180, 180);   
  hist_cmphi_pim->GetXaxis()->SetTitle("#phi_{CM} /deg");
  hist_cmphi_pim->GetYaxis()->SetTitle("counts");
  hist_cmcostheta_pim = new TH1F("cmcostheta_pim", "cos(#Theta_{CM})", 100, -1, 1);   
  hist_cmcostheta_pim->GetXaxis()->SetTitle("cos(#Theta_{CM})");
  hist_cmcostheta_pim->GetYaxis()->SetTitle("counts");


out->mkdir("missing_mass");				
out->cd ("missing_mass");

  TH1F *hist_e_p_mismass;
  TH1F *hist_e_p_mismass2;
  TH1F *hist_e_pi0_mismass;
  TH1F *hist_e_pip_mismass;
  TH1F *hist_e_pim_mismass;

  TH1F *hist_e_p_mismass_sec[6];
  TH1F *hist_e_p_mismass2_sec[6];
  TH1F *hist_e_pi0_mismass_sec[6];
  TH1F *hist_e_pip_mismass_sec[6];
  TH1F *hist_e_pim_mismass_sec[6];

  hist_e_p_mismass = new TH1F("hist_e_p_mismass", "e p --> e p X missing mass", 400, 0.0, 4.0);   
  hist_e_p_mismass->GetXaxis()->SetTitle("missing mass /GeV");
  hist_e_p_mismass->GetYaxis()->SetTitle("counts");
  hist_e_p_mismass2 = new TH1F("hist_e_p_mismass2", "e p --> e p X missing mass squared", 1000, -0.5, 1);   
  hist_e_p_mismass2->GetXaxis()->SetTitle("missing mass /GeV");
  hist_e_p_mismass2->GetYaxis()->SetTitle("counts");
  hist_e_pi0_mismass = new TH1F("hist_e_pi0_mismass", "e p --> e #pi^{0} X missing mass", 1000, 0.0, 5.0);   
  hist_e_pi0_mismass->GetXaxis()->SetTitle("missing mass /GeV");
  hist_e_pi0_mismass->GetYaxis()->SetTitle("counts");
  hist_e_pip_mismass = new TH1F("hist_e_pip_mismass", "e p --> e #pi^{+} X missing mass", 1000, 0.0, 5.0);   
  hist_e_pip_mismass->GetXaxis()->SetTitle("missing mass /GeV");
  hist_e_pip_mismass->GetYaxis()->SetTitle("counts");
  hist_e_pim_mismass = new TH1F("hist_e_pim_mismass", "e p --> e #pi^{-} X missing mass", 1000, 0.0, 5.0);   
  hist_e_pim_mismass->GetXaxis()->SetTitle("missing mass /GeV");
  hist_e_pim_mismass->GetYaxis()->SetTitle("counts");

  for(Int_t i = 0; i < 6; i++){
    sprintf(name,"hist_e_p_mismass_sec%01d", i+1);
    sprintf(title,"e p --> e p X missing mass sector %01d", i+1);
    hist_e_p_mismass_sec[i] = new TH1F(name, title, 400, 0.0, 4.0);   
    hist_e_p_mismass_sec[i]->GetXaxis()->SetTitle("missing mass /GeV");
    hist_e_p_mismass_sec[i]->GetYaxis()->SetTitle("counts");
    sprintf(name,"hist_e_p_mismass2_sec%01d", i+1);
    sprintf(title,"e p --> e p X missing mass squared sector %01d", i+1);
    hist_e_p_mismass2_sec[i] = new TH1F(name, title, 1000, -0.5, 1);   
    hist_e_p_mismass2_sec[i]->GetXaxis()->SetTitle("missing mass /GeV");
    hist_e_p_mismass2_sec[i]->GetYaxis()->SetTitle("counts");
    sprintf(name,"hist_e_pi0_mismass_sec%01d", i+1);
    sprintf(title,"e p --> e #pi^{0} X missing mass sector %01d", i+1);
    hist_e_pi0_mismass_sec[i] = new TH1F(name, title, 1000, 0.0, 5.0);   
    hist_e_pi0_mismass_sec[i]->GetXaxis()->SetTitle("missing mass /GeV");
    hist_e_pi0_mismass_sec[i]->GetYaxis()->SetTitle("counts");
    sprintf(name,"hist_e_pip_mismass_sec%01d", i+1);
    sprintf(title,"e p --> e #pi^{+} X missing mass sector %01d", i+1);
    hist_e_pip_mismass_sec[i] = new TH1F(name, title, 1000, 0.0, 5.0);   
    hist_e_pip_mismass_sec[i]->GetXaxis()->SetTitle("missing mass /GeV");
    hist_e_pip_mismass_sec[i]->GetYaxis()->SetTitle("counts");
    sprintf(name,"hist_e_pim_mismass_sec%01d", i+1);
    sprintf(title,"e p --> e #pi^{-} X missing mass sector %01d", i+1);
    hist_e_pim_mismass_sec[i] = new TH1F(name, title, 1000, 0.0, 5.0);   
    hist_e_pim_mismass_sec[i]->GetXaxis()->SetTitle("missing mass /GeV");
    hist_e_pim_mismass_sec[i]->GetYaxis()->SetTitle("counts");
  }


out->mkdir("Dalitz_plots");				
out->cd ("Dalitz_plots");

  TH2F *hist_Dalitz_ppip_pippim;
  TH2F *hist_Dalitz_ppim_pippim;
  TH2F *hist_Dalitz_ppip_ppim;

  TH2F *hist_Dalitz2_ppip_pippim;
  TH2F *hist_Dalitz2_ppim_pippim;
  TH2F *hist_Dalitz2_ppip_ppim;

  TH1F *hist_mass_ppip;
  TH1F *hist_mass_ppim;
  TH1F *hist_mass_pippim;

  hist_Dalitz_ppip_pippim = new TH2F("Dalitz_ppip_pippim", "Dalitz p #pi^{+} vs #pi^{+} #pi^{-}", 150, 1.0, 4.0, 150, 0.0, 3.0);   
  hist_Dalitz_ppip_pippim->GetXaxis()->SetTitle("mass p #pi^{+} /GeV");
  hist_Dalitz_ppip_pippim->GetYaxis()->SetTitle("mass #pi^{+} #pi^{-} /GeV");
  hist_Dalitz_ppim_pippim = new TH2F("Dalitz_ppim_pippim", "Dalitz p #pi^{-} vs #pi^{+} #pi^{-}", 150, 1.0, 4.0, 150, 0.0, 3.0);   
  hist_Dalitz_ppim_pippim->GetXaxis()->SetTitle("mass p #pi^{-} /GeV");
  hist_Dalitz_ppim_pippim->GetYaxis()->SetTitle("mass #pi^{+} #pi^{-} /GeV");
  hist_Dalitz_ppip_ppim = new TH2F("Dalitz_ppip_ppim", "Dalitz p #pi^{+} vs p #pi^{-}", 150, 1.0, 4.0, 150, 1.0, 4.0);   
  hist_Dalitz_ppip_ppim->GetXaxis()->SetTitle("mass p #pi^{+} /GeV");
  hist_Dalitz_ppip_ppim->GetYaxis()->SetTitle("mass p #pi^{-} /GeV");

  hist_Dalitz2_ppip_pippim = new TH2F("Dalitz2_ppip_pippim", "Dalitz p #pi^{+} vs #pi^{+} #pi^{-}", 150, 1.0, 12.0, 150, 0.0, 8.0);   
  hist_Dalitz2_ppip_pippim->GetXaxis()->SetTitle("mass p #pi^{+} /GeV");
  hist_Dalitz2_ppip_pippim->GetYaxis()->SetTitle("mass #pi^{+} #pi^{-} /GeV");
  hist_Dalitz2_ppim_pippim = new TH2F("Dalitz2_ppim_pippim", "Dalitz p #pi^{-} vs #pi^{+} #pi^{-}", 150, 1.0, 12.0, 150, 0.0, 8.0);   
  hist_Dalitz2_ppim_pippim->GetXaxis()->SetTitle("mass p #pi^{-} /GeV");
  hist_Dalitz2_ppim_pippim->GetYaxis()->SetTitle("mass #pi^{+} #pi^{-} /GeV");
  hist_Dalitz2_ppip_ppim = new TH2F("Dalitz2_ppip_ppim", "Dalitz p #pi^{+} vs p #pi^{-}", 150, 1.0, 12.0, 150, 1.0, 12.0);   
  hist_Dalitz2_ppip_ppim->GetXaxis()->SetTitle("mass p #pi^{+} /GeV");
  hist_Dalitz2_ppip_ppim->GetYaxis()->SetTitle("mass p #pi^{-} /GeV");

  hist_mass_ppip = new TH1F("mass_ppip", "mass p #pi^{+}", 310, 0.9, 4.0);   
  hist_mass_ppip->GetXaxis()->SetTitle("mass p #pi^{+} /GeV");
  hist_mass_ppip->GetYaxis()->SetTitle("counts");
  hist_mass_ppim = new TH1F("mass_ppim", "mass p #pi^{-}", 310, 0.9, 4.0);   
  hist_mass_ppim->GetXaxis()->SetTitle("mass p #pi^{-} /GeV");
  hist_mass_ppim->GetYaxis()->SetTitle("counts");
  hist_mass_pippim = new TH1F("mass_pippim", "mass #pi^{+} #pi^{-}", 300, 0.0, 3.0);   
  hist_mass_pippim->GetXaxis()->SetTitle("mass #pi^{+} #pi^{-} /GeV");
  hist_mass_pippim->GetYaxis()->SetTitle("count");


out->mkdir("additional_plots");				
out->cd ("additional_plots");

  TH2F *hist_FTOF_phi_vs_sector_positives;
  TH2F *hist_FTOF_phi_vs_sector_negatives;

  hist_FTOF_phi_vs_sector_positives = new TH2F("FTOF_phi_vs_sector_positives", "#phi vs FTOF sector for particles with positive charge", 6, 0.5, 6.5, 180, -180, 180);   
  hist_FTOF_phi_vs_sector_positives->GetXaxis()->SetTitle("FTOF sector");
  hist_FTOF_phi_vs_sector_positives->GetYaxis()->SetTitle("#phi /deg");

  hist_FTOF_phi_vs_sector_negatives = new TH2F("FTOF_phi_vs_sector_negatives", "#phi vs FTOF sector for particles with negative charge", 6, 0.5, 6.5, 180, -180, 180);   
  hist_FTOF_phi_vs_sector_negatives->GetXaxis()->SetTitle("FTOF sector");
  hist_FTOF_phi_vs_sector_negatives->GetYaxis()->SetTitle("#phi /deg");


out->mkdir("sampling_fraction");				
out->cd ("sampling_fraction");

  TH2F *hist_sampling_fraction_vs_E_sec1;
  TH2F *hist_sampling_fraction_vs_E_sec2;
  TH2F *hist_sampling_fraction_vs_E_sec3;
  TH2F *hist_sampling_fraction_vs_E_sec4;
  TH2F *hist_sampling_fraction_vs_E_sec5;
  TH2F *hist_sampling_fraction_vs_E_sec6;

  TH2F *hist_sampling_fraction_vs_p_sec1;
  TH2F *hist_sampling_fraction_vs_p_sec2;
  TH2F *hist_sampling_fraction_vs_p_sec3;
  TH2F *hist_sampling_fraction_vs_p_sec4;
  TH2F *hist_sampling_fraction_vs_p_sec5;
  TH2F *hist_sampling_fraction_vs_p_sec6;

  hist_sampling_fraction_vs_E_sec1 = new TH2F("sampling_fraction_vs_E_sec1", "EC total sampling fraction versus E sec 1", 350,0,3.5, 300,0,0.6);   
  hist_sampling_fraction_vs_E_sec1->GetXaxis()->SetTitle("E /GeV");
  hist_sampling_fraction_vs_E_sec1->GetYaxis()->SetTitle("E/p");
  hist_sampling_fraction_vs_E_sec2 = new TH2F("sampling_fraction_vs_E_sec2", "EC total sampling fraction versus E sec 2", 350,0,3.5, 300,0,0.6);   
  hist_sampling_fraction_vs_E_sec2->GetXaxis()->SetTitle("E /GeV");
  hist_sampling_fraction_vs_E_sec2->GetYaxis()->SetTitle("E/p");
  hist_sampling_fraction_vs_E_sec3 = new TH2F("sampling_fraction_vs_E_sec3", "EC total sampling fraction versus E sec 3", 350,0,3.5, 300,0,0.6);   
  hist_sampling_fraction_vs_E_sec3->GetXaxis()->SetTitle("E /GeV");
  hist_sampling_fraction_vs_E_sec3->GetYaxis()->SetTitle("E/p");
  hist_sampling_fraction_vs_E_sec4 = new TH2F("sampling_fraction_vs_E_sec4", "EC total sampling fraction versus E sec 4", 350,0,3.5, 300,0,0.6);   
  hist_sampling_fraction_vs_E_sec4->GetXaxis()->SetTitle("E /GeV"); 
  hist_sampling_fraction_vs_E_sec4->GetYaxis()->SetTitle("E/p");
  hist_sampling_fraction_vs_E_sec5 = new TH2F("sampling_fraction_vs_E_sec5", "EC total sampling fraction versus E sec 5", 350,0,3.5, 300,0,0.6);   
  hist_sampling_fraction_vs_E_sec5->GetXaxis()->SetTitle("E /GeV");
  hist_sampling_fraction_vs_E_sec5->GetYaxis()->SetTitle("E/p");
  hist_sampling_fraction_vs_E_sec6 = new TH2F("sampling_fraction_vs_E_sec6", "EC total sampling fraction versus E sec 6", 350,0,3.5, 300,0,0.6);   
  hist_sampling_fraction_vs_E_sec6->GetXaxis()->SetTitle("E /GeV");
  hist_sampling_fraction_vs_E_sec6->GetYaxis()->SetTitle("E/p");

  hist_sampling_fraction_vs_p_sec1 = new TH2F("sampling_fraction_vs_p_sec1", "EC total sampling fraction versus p sec 1", 550,0,11, 300,0,0.6);   
  hist_sampling_fraction_vs_p_sec1->GetXaxis()->SetTitle("p /GeV");
  hist_sampling_fraction_vs_p_sec1->GetYaxis()->SetTitle("E/p");
  hist_sampling_fraction_vs_p_sec2 = new TH2F("sampling_fraction_vs_p_sec2", "EC total sampling fraction versus p sec 2", 550,0,11, 300,0,0.6);   
  hist_sampling_fraction_vs_p_sec2->GetXaxis()->SetTitle("p /GeV");
  hist_sampling_fraction_vs_p_sec2->GetYaxis()->SetTitle("E/p");
  hist_sampling_fraction_vs_p_sec3 = new TH2F("sampling_fraction_vs_p_sec3", "EC total sampling fraction versus p sec 3", 550,0,11, 300,0,0.6);   
  hist_sampling_fraction_vs_p_sec3->GetXaxis()->SetTitle("p /GeV");
  hist_sampling_fraction_vs_p_sec3->GetYaxis()->SetTitle("E/p");
  hist_sampling_fraction_vs_p_sec4 = new TH2F("sampling_fraction_vs_p_sec4", "EC total sampling fraction versus p sec 4", 550,0,11, 300,0,0.6);   
  hist_sampling_fraction_vs_p_sec4->GetXaxis()->SetTitle("p /GeV"); 
  hist_sampling_fraction_vs_p_sec4->GetYaxis()->SetTitle("E/p");
  hist_sampling_fraction_vs_p_sec5 = new TH2F("sampling_fraction_vs_p_sec5", "EC total sampling fraction versus p sec 5", 550,0,11, 300,0,0.6);   
  hist_sampling_fraction_vs_p_sec5->GetXaxis()->SetTitle("p /GeV");
  hist_sampling_fraction_vs_p_sec5->GetYaxis()->SetTitle("E/p");
  hist_sampling_fraction_vs_p_sec6 = new TH2F("sampling_fraction_vs_p_sec6", "EC total sampling fraction versus p sec 6", 550,0,11, 300,0,0.6);   
  hist_sampling_fraction_vs_p_sec6->GetXaxis()->SetTitle("p /GeV");
  hist_sampling_fraction_vs_p_sec6->GetYaxis()->SetTitle("E/p");



out->mkdir("MC_LUND");
out->cd ("MC_LUND");


  TH1F *hist_LUND_pid;

  hist_LUND_pid = new TH1F("hist_LUND_pid", "generated PID", 6000, -3000, 3000);   
  hist_LUND_pid->GetXaxis()->SetTitle("PID");
  hist_LUND_pid->GetYaxis()->SetTitle("counts");


/// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///  start of the event loop     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 event = 0;

neg_part_count = 0;
pos_part_count = 0;
neut_part_count = 0;

cout << "Analysing Tree: " << inTree << endl;
cout << "Event Loop starting ... " << endl;
cout << endl;

 int nentries = 500000; //anaTree->GetEntriesFast();

for(Int_t k=0; k<nentries;k++){    

  anaTree->GetEntry(k);
  event=k;
  if(process_Events > 0 && k == process_Events) break;


/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// progress:

  if(k % 10000 == 0){

      double events = anaTree->GetEntriesFast()/100;
      double percent = (k/100)/(events/100);
      
      printf("Analysing event number %i of %.00f (%.01f percent)\n", k, events*100, percent);
  }


/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// initalisatize the particle fourvectors and their components:
///

TLorentzVector beam(0,0,Ebeam,Ebeam);
TLorentzVector target(0,0,0,0.93827);

NRUN = 0; NEVENT = 0; TYPE = 0; TRG = 0; Helic = 0; EVNTime = 0; BCG = 0; STTime = 0; RFTime = 0;

helicity = 0;
fcup = BCG;

p4_ele_px.clear(); p4_prot_px.clear(); p4_neutr_px.clear(); p4_pip_px.clear(); p4_pim_px.clear(); p4_Kp_px.clear(); p4_Km_px.clear(); p4_phot_px.clear();
p4_ele_py.clear(); p4_prot_py.clear(); p4_neutr_py.clear(); p4_pip_py.clear(); p4_pim_py.clear(); p4_Kp_py.clear(); p4_Km_py.clear(); p4_phot_py.clear();
p4_ele_pz.clear(); p4_prot_pz.clear(); p4_neutr_pz.clear(); p4_pip_pz.clear(); p4_pim_pz.clear(); p4_Kp_pz.clear(); p4_Km_pz.clear(); p4_phot_pz.clear();
//added
 sectorE.clear(); electron_event_number.clear();

p4_ele_E.clear(); p4_prot_E.clear(); p4_neutr_E.clear(); p4_pip_E.clear(); p4_pim_E.clear(); p4_Kp_E.clear(); p4_Km_E.clear(); p4_phot_E.clear();
ele_det.clear(); prot_det.clear(); neutr_det.clear(); pip_det.clear(); pim_det.clear(); Kp_det.clear(); Km_det.clear(); phot_det.clear();

e_count = 0; p_count = 0; n_count = 0; pip_count = 0; pim_count = 0; Kp_count = 0; Km_count = 0; g_count = 0; 
e_MCcount = 0; p_MCcount = 0; n_MCcount = 0; pip_MCcount = 0; pim_MCcount = 0; Kp_MCcount = 0; Km_MCcount = 0; g_MCcount = 0; 

W = 0; Q2 = 0; nu = 0; x = 0; y = 0; t_pi0 = 0;  t_pip = 0; t_pim = 0;
cmphi_p = 0; cmcostheta_p = 0; cmphi_pi0 = 0; cmcostheta_pi0 = 0; cmphi_pip = 0; cmcostheta_pip = 0; cmphi_pim = 0; cmcostheta_pim = 0;

MC_helicity = 0; MC_Npart = 0; MC_Ebeam = 0; MC_weight = 0;

for(Int_t i = 0; i < BUFFER; i++){ 

  p4_ele[i].SetPxPyPzE(0,0,0,0);
  p4_ele_raw[i].SetPxPyPzE(0,0,0,0);
  p4_prot[i].SetPxPyPzE(0,0,0,0);
  p4_neutr[i].SetPxPyPzE(0,0,0,0);
  p4_pip[i].SetPxPyPzE(0,0,0,0);
  p4_pim[i].SetPxPyPzE(0,0,0,0);
  p4_Kp[i].SetPxPyPzE(0,0,0,0);
  p4_Km[i].SetPxPyPzE(0,0,0,0);
  p4_phot[i].SetPxPyPzE(0,0,0,0);

  e_ind[i] = -1;  p_ind[i] = -1;  n_ind[i] = -1;  pip_ind[i] = -1;  pim_ind[i] = -1;  Kp_ind[i] = -1;  Km_ind[i] = -1;  g_ind[i] = -1;

  ele_detect[i] = 0; 
  prot_detect[i] = 0; neutr_detect[i] = 0;
  pip_detect[i] = 0; pim_detect[i] = 0;
  Kp_detect[i] = 0; Km_detect[i] = 0;
  phot_detect[i] = 0;

  e_vx[i] = 0;  e_vy[i] = 0;  e_vz[i] = 0;  e_beta[i] = 0, e_FTOF_sec[i] = -1, e_PCAL_sec[i] = -1;
  p_vx[i] = 0;  p_vy[i] = 0;  p_vz[i] = 0;  p_beta[i] = 0, p_FTOF_sec[i] = -1, p_PCAL_sec[i] = -1; 

 eventNumber[i]=-1;

  n_vx[i] = 0;  n_vy[i] = 0;  n_vz[i] = 0;  n_beta[i] = 0; 
  pip_vx[i] = 0;  pip_vy[i] = 0;  pip_vz[i] = 0;  pip_beta[i] = 0, pip_FTOF_sec[i] = -1; 
  pim_vx[i] = 0;  pim_vy[i] = 0;  pim_vz[i] = 0;  pim_beta[i] = 0, pim_FTOF_sec[i] = -1; 
  Kp_vx[i] = 0;  Kp_vy[i] = 0;  Kp_vz[i] = 0;  Kp_beta[i] = 0, Kp_FTOF_sec[i] = -1; 
  Km_vx[i] = 0;  Km_vy[i] = 0;  Km_vz[i] = 0;  Km_beta[i] = 0, Km_FTOF_sec[i] = -1; 
  g_vx[i] = 0;  g_vy[i] = 0;  g_vz[i] = 0;  g_sec[i] = 0; 

  electron_sector_cut[i] = false;

}

for(Int_t i = 0; i < BUFFER; i++){ 

  part_px[i] = 0; part_py[i] = 0; part_pz[i] = 0; part_p[i] = 0; part_beta[i] = 0;
  part_vx[i] = 0; part_vy[i] = 0; part_vz[i] = 0; part_charge[i] = 0; part_pid[i] = 0;
  part_theta[i] = 0; part_phi[i] = 0; part_status[i] = 0;

  partMC_px[i] = 0; partMC_py[i] = 0; partMC_pz[i] = 0; partMC_p[i] = 0; partMC_theta[i] = 0; partMC_phi[i] = 0;
  partMC_vx[i] = 0; partMC_vy[i] = 0; partMC_vz[i] = 0; partMC_pid[i] = 0;
  partLUND_mass[i] = 0; partLUND_E[i] = 0;
  partLUND_px[i] = 0; partLUND_py[i] = 0; partLUND_pz[i] = 0; partLUND_p[i] = 0; partLUND_theta[i] = 0; partLUND_phi[i] = 0;
  partLUND_vx[i] = 0; partLUND_vy[i] = 0; partLUND_vz[i] = 0; partLUND_pid[i] = 0;

  part_Cal_PCAL_sector[i] = 0; part_Cal_ECin_sector[i] = 0; part_Cal_ECout_sector[i] = 0;  
  part_Cal_PCAL_energy[i] = 0; part_Cal_ECin_energy[i] = 0; part_Cal_ECout_energy[i] = 0; part_Cal_energy_total[i] = 0;
  part_Cal_PCAL_time[i] = 0; part_Cal_ECin_time[i] = 0; part_Cal_ECout_time[i] = 0; 
  part_Cal_PCAL_path[i] = 0; part_Cal_ECin_path[i] = 0; part_Cal_ECout_path[i] = 0;  
  part_Cal_PCAL_x[i] = 0; part_Cal_PCAL_y[i] = 0; part_Cal_PCAL_z[i] = 0;
  part_Cal_ECin_x[i] = 0; part_Cal_ECin_y[i] = 0; part_Cal_ECin_z[i] = 0;
  part_Cal_ECout_x[i] = 0; part_Cal_ECout_y[i] = 0; part_Cal_ECout_z[i] = 0;
  part_Cal_PCAL_lu[i] = 0; part_Cal_PCAL_lv[i] = 0; part_Cal_PCAL_lw[i] = 0;
  part_Cal_ECin_lu[i] = 0; part_Cal_ECin_lv[i] = 0; part_Cal_ECin_lw[i] = 0;
  part_Cal_ECout_lu[i] = 0; part_Cal_ECout_lv[i] = 0; part_Cal_ECout_lw[i] = 0;

  part_CC_HTCC_sector[i] = 0; part_CC_HTCC_nphe[i] = 0; part_CC_HTCC_time[i] = 0; part_CC_HTCC_path[i] = 0;   
  part_CC_HTCC_theta[i] = 0; part_CC_HTCC_phi[i] = 0;

  part_CC_LTCC_sector[i] = 0; part_CC_LTCC_nphe[i] = 0; part_CC_LTCC_time[i] = 0; part_CC_LTCC_path[i] = 0;   
  part_CC_LTCC_theta[i] = 0; part_CC_LTCC_phi[i] = 0;

  part_FTOF_sector_layer1[i] = 0; part_FTOF_sector_layer2[i] = 0; part_FTOF_sector_layer3[i] = 0; 
  part_FTOF_component_layer1[i] = 0; part_FTOF_component_layer2[i] = 0; part_FTOF_component_layer3[i] = 0;  
  part_FTOF_energy[i] = 0; part_FTOF_time[i] = 0; part_FTOF_path[i] = 0;
  part_FTOF_energy_layer1[i] = 0; part_FTOF_time_layer1[i] = 0; part_FTOF_path_layer1[i] = 0;  
  part_FTOF_energy_layer3[i] = 0; part_FTOF_time_layer3[i] = 0; part_FTOF_path_layer3[i] = 0;  
  part_FTOF_layer[i] = 0;
  
  part_CTOF_component[i] = 0; 
  part_CTOF_energy[i] = 0; part_CTOF_time[i] = 0; part_CTOF_path[i] = 0;

  part_CND_component[i] = 0; 
  part_CND_energy[i] = 0; part_CND_time[i] = 0; part_CND_path[i] = 0;

  part_FT_energy[i] = 0; part_FT_radius[i] = 0; part_FTHODO_time[i] = 0; part_FTHODO_path[i] = 0; part_FTCAL_time[i] = 0; part_FTCAL_path[i] = 0;  
  part_FTCAL_x[i] = 0; part_FTCAL_y[i] = 0; part_FTCAL_z[i] = 0; part_FTTRK_x[i] = 0; part_FTTRK_y[i] = 0; part_FTTRK_z[i] = 0; 
  part_FTHODO_x[i] = 0; part_FTHODO_y[i] = 0; part_FTHODO_z[i] = 0; 

  part_DC_index[i] = 0; part_DC_sector[i] = 0; 
  part_DC_c1x[i] = 0; part_DC_c1y[i] = 0; part_DC_c1z[i] = 0;
  part_DC_c2x[i] = 0; part_DC_c2y[i] = 0; part_DC_c2z[i] = 0;
  part_DC_c3x[i] = 0; part_DC_c3y[i] = 0; part_DC_c3z[i] = 0;

  FD_eid_default_PID_check[i] = false;
  FD_eid_charge_check[i] = false;
  FD_eid_EC_outer_vs_EC_inner_check[i] = false;
  FD_eid_EC_sampling_fraction_check[i] = false;
  FD_eid_EC_hit_position_fiducial_check[i] = false;
  FD_eid_DC_hit_position_region1_fiducial_check[i] = false;
  FD_eid_DC_hit_position_region2_fiducial_check[i] = false;
  FD_eid_DC_hit_position_region3_fiducial_check[i] = false;
  FD_eid_DC_z_vertex_check[i] = false;
  FD_eid_all_check[i] = false;

  FD_protid_default_PID_check[i] = false;
  FD_protid_charge_check[i] = false;
  FD_protid_DC_hit_position_region1_fiducial_check[i] = false;
  FD_protid_DC_hit_position_region2_fiducial_check[i] = false;
  FD_protid_DC_hit_position_region3_fiducial_check[i] = false;
  FD_protid_beta_check[i] = false;
  FD_protid_delta_beta_check[i] = false;
  FD_protid_tofmass_check[i] = false;
  FD_protid_maximum_probability_check[i] = false;
  FD_protid_delta_vz_check[i] = false;
  FD_protid_all_check[i] = false;
  FD_neutrid_default_PID_check[i] = false;
  FD_neutrid_charge_check[i] = false;
  FD_neutrid_beta_check[i] = false;
  FD_neutrid_delta_beta_check[i] = false;
  FD_neutrid_tofmass_check[i] = false;
  FD_neutrid_delta_vz_check[i] = false;
  FD_neutrid_all_check[i] = false;
  FD_pipid_default_PID_check[i] = false;
  FD_pipid_charge_check[i] = false;
  FD_pipid_DC_hit_position_region1_fiducial_check[i] = false;
  FD_pipid_DC_hit_position_region2_fiducial_check[i] = false;
  FD_pipid_DC_hit_position_region3_fiducial_check[i] = false;
  FD_pipid_beta_check[i] = false;
  FD_pipid_delta_beta_check[i] = false;
  FD_pipid_tofmass_check[i] = false;
  FD_pipid_maximum_probability_check[i] = false;
  FD_pipid_delta_vz_check[i] = false;
  FD_pipid_all_check[i] = false;
  FD_pimid_default_PID_check[i] = false;
  FD_pimid_charge_check[i] = false;
  FD_pimid_DC_hit_position_region1_fiducial_check[i] = false;
  FD_pimid_DC_hit_position_region2_fiducial_check[i] = false;
  FD_pimid_DC_hit_position_region3_fiducial_check[i] = false;
  FD_pimid_beta_check[i] = false;
  FD_pimid_delta_beta_check[i] = false;
  FD_pimid_tofmass_check[i] = false;
  FD_pimid_maximum_probability_check[i] = false;
  FD_pimid_delta_vz_check[i] = false;
  FD_pimid_all_check[i] = false;
  FD_Kpid_default_PID_check[i] = false;
  FD_Kpid_charge_check[i] = false;
  FD_Kpid_DC_hit_position_region1_fiducial_check[i] = false;
  FD_Kpid_DC_hit_position_region2_fiducial_check[i] = false;
  FD_Kpid_DC_hit_position_region3_fiducial_check[i] = false;
  FD_Kpid_beta_check[i] = false;
  FD_Kpid_delta_beta_check[i] = false;
  FD_Kpid_tofmass_check[i] = false;
  FD_Kpid_maximum_probability_check[i] = false;
  FD_Kpid_delta_vz_check[i] = false;
  FD_Kpid_all_check[i] = false;
  FD_Kmid_default_PID_check[i] = false;
  FD_Kmid_charge_check[i] = false;
  FD_Kmid_DC_hit_position_region1_fiducial_check[i] = false;
  FD_Kmid_DC_hit_position_region2_fiducial_check[i] = false;
  FD_Kmid_DC_hit_position_region3_fiducial_check[i] = false;
  FD_Kmid_beta_check[i] = false;
  FD_Kmid_delta_beta_check[i] = false;
  FD_Kmid_tofmass_check[i] = false;
  FD_Kmid_maximum_probability_check[i] = false;
  FD_Kmid_delta_vz_check[i] = false;
  FD_Kmid_all_check[i] = false;

  FD_photid_default_PID_check[i] = false;
  FD_photid_charge_check[i] = false;
  FD_photid_beta_check[i] = false;
  FD_photid_EC_sampling_fraction_check[i] = false;
  FD_photid_EC_hit_position_fiducial_check[i] = false;
  FD_photid_all_check[i] = false;

  // FT

  FT_eid_charge_check[i] = false;
  FT_eid_PID_check[i] = false;
  FT_eid_FTCAL_fiducial_check[i] = false;
  FT_eid_FTTRK_fiducial_check[i] = false;
  FT_eid_FTHODO_fiducial_check[i] = false;
  FT_eid_energy_vs_radius_check[i] = false;
  FT_eid_all_check[i] = false;

  FT_photid_charge_check[i] = false;
  FT_photid_PID_check[i] = false;
  FT_photid_FTCAL_fiducial_check[i] = false;
  FT_photid_beta_check[i] = false;
  FT_photid_all_check[i] = false;

  // CD

  CD_protid_default_PID_check[i] = false;
  CD_protid_charge_check[i] = false;
  CD_protid_beta_check[i] = false;
  CD_protid_maximum_probability_check[i] = false;
  CD_protid_delta_vz_check[i] = false;
  CD_protid_all_check[i] = false;

  CD_neutrid_default_PID_check[i] = false;
  CD_neutrid_charge_check[i] = false;
  CD_neutrid_beta_check[i] = false;
  CD_neutrid_maximum_probability_check[i] = false;
  CD_neutrid_delta_vz_check[i] = false;
  CD_neutrid_all_check[i] = false;

  CD_pipid_default_PID_check[i] = false;
  CD_pipid_charge_check[i] = false;
  CD_pipid_beta_check[i] = false;
  CD_pipid_maximum_probability_check[i] = false;
  CD_pipid_delta_vz_check[i] = false;
  CD_pipid_all_check[i] = false;

  CD_pimid_default_PID_check[i] = false;
  CD_pimid_charge_check[i] = false;
  CD_pimid_beta_check[i] = false;
  CD_pimid_maximum_probability_check[i] = false;
  CD_pimid_delta_vz_check[i] = false;
  CD_pimid_all_check[i] = false;

  CD_Kpid_default_PID_check[i] = false;
  CD_Kpid_charge_check[i] = false;
  CD_Kpid_beta_check[i] = false;
  CD_Kpid_maximum_probability_check[i] = false;
  CD_Kpid_delta_vz_check[i] = false;
  CD_Kpid_all_check[i] = false;

  CD_Kmid_default_PID_check[i] = false;
  CD_Kmid_charge_check[i] = false;
  CD_Kmid_beta_check[i] = false;
  CD_Kmid_maximum_probability_check[i] = false;
  CD_Kmid_delta_vz_check[i] = false;
  CD_Kmid_all_check[i] = false;

}


for(Int_t i = 0; i <= 27; i++){ 
  neutral_iter[i].SetPxPyPzE(0,0,0,0);
  mass_neutral_iter[i] = 0;
  mass_neutral_iter2[i] = 0;
  alpha_gg[i] = 0;

  pi0_g_ind[i] = -1;

  neutral_ind[i][0] = -1;
  neutral_ind[i][1] = -1;

}

pi0.SetPxPyPzE(0,0,0,0);
eta.SetPxPyPzE(0,0,0,0);

pi0_ind[0] = -1;
pi0_ind[1] = -1;
eta_ind[0] = -1;
eta_ind[1] = -1;


/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  get event properties and assign particles to variables
/// 

get_event_properties();

helicity = Helic;

assign_particles();

if(simulation == true && MC_helicity != 0){ helicity = MC_helicity;}

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  do particle ID an momentum correction
/// 

select_electron(run);   // first select good electrons

if(e_count > 0){    // basic trigger condition: at least one good electron

  // continue to select other particles
  select_proton(run);
  select_neutron(run);
  select_pip(run);
  select_pim(run);
  select_Kplus(run);
  select_Kminus(run);
  select_photon(run);

  // transfer uncorrected electron 4 vector

  for(Int_t i = 0; i < e_count; i++){ 
    p4_ele_raw[i].SetPxPyPzE(p4_ele[i].Px(), p4_ele[i].Py(), p4_ele[i].Pz(), p4_ele[i].E());
  }

  // do momentum correction:

  if(apply_correction_electron) correct_electron();
  if(apply_correction_proton)   correct_proton();
  if(apply_correction_neutron)  correct_neutron();
  if(apply_correction_pip)      correct_pip();
  if(apply_correction_pim)      correct_pim();
  if(apply_correction_Kp)       correct_Kplus();
  if(apply_correction_Km)       correct_Kminus();
  if(apply_correction_photon)   correct_photon();

  // reconstruct neutrals from photons

  create_neutrals();

  // fill output tree variables:

  fill_output_vector_electron();
  fill_output_vector_proton();
  fill_output_vector_neutron();
  fill_output_vector_pip();
  fill_output_vector_pim();
  fill_output_vector_Kplus();
  fill_output_vector_Kminus();
  fill_output_vector_photon();

}  // end of basic trigger condition


// write particles to the tree:

out_tree.Fill();



/// ///////////////////////////////////////////////////////////////////////////////////////
/// calculate kinematics
/// ///////////////////////////////////////////////////////////////////////////////////////

  W = kin_W(p4_ele[0]);
  Q2 = kin_Q2(p4_ele[0]);
  x = kin_x(p4_ele[0]);
  y = kin_y(p4_ele[0]);
  nu = kin_nu(p4_ele[0]);

  t_pi0 = kin_t(p4_ele[0], pi0);
  t_pip = kin_t(p4_ele[0], p4_pip[0]);
  t_pim = kin_t(p4_ele[0], p4_pim[0]);

  Double_t cmcostheta_p = kin_cmcostheta(p4_ele[0], p4_prot[0]);
  Double_t cmphi_p = kin_cmphi(p4_ele[0], p4_prot[0]);

  Double_t cmcostheta_pi0 = kin_cmcostheta(p4_ele[0], pi0);
  Double_t cmphi_pi0 = kin_cmphi(p4_ele[0], pi0);

  Double_t cmcostheta_pip = kin_cmcostheta(p4_ele[0], p4_pip[0]);
  Double_t cmphi_pip = kin_cmphi(p4_ele[0], p4_pip[0]);

  Double_t cmcostheta_pim = kin_cmcostheta(p4_ele[0], p4_pim[0]);
  Double_t cmphi_pim = kin_cmphi(p4_ele[0], p4_pim[0]);

/// ////////////////////////////////////////////////////////////////////////////////////////
/// calculate missing masses and two pion masses
/// ////////////////////////////////////////////////////////////////////////////////////////

  TLorentzVector fElectron(p4_ele[0].Px(), p4_ele[0].Py(), p4_ele[0].Pz(), p4_ele[0].E());
  TLorentzVector fProton(p4_prot[0].Px(), p4_prot[0].Py(), p4_prot[0].Pz(), p4_prot[0].E());
  TLorentzVector fNeutron(p4_neutr[0].Px(), p4_neutr[0].Py(), p4_neutr[0].Pz(), p4_neutr[0].E());
  TLorentzVector fPip(p4_pip[0].Px(), p4_pip[0].Py(), p4_pip[0].Pz(), p4_pip[0].E());
  TLorentzVector fPim(p4_pim[0].Px(), p4_pim[0].Py(), p4_pim[0].Pz(), p4_pim[0].E());
  TLorentzVector fPi0(pi0.Px(), pi0.Py(), pi0.Pz(), pi0.E());
  TLorentzVector fEta(eta.Px(), eta.Py(), eta.Pz(), eta.E());

  TLorentzVector f_e_p_miss   = beam + target - fElectron - fProton;
  TLorentzVector f_e_pi0_miss   = beam + target - fElectron - pi0;
  if(pi0.E() < 0.005) f_e_pi0_miss.SetPxPyPzE(0,0,0,0);
  TLorentzVector f_e_pip_miss = beam + target - fElectron - fPip;
  TLorentzVector f_n_pip_miss = beam + target - fNeutron - fPip;
  TLorentzVector f_e_pim_miss = beam + target - fElectron - fPim;

  TLorentzVector f_ppip   = fProton + fPip;
  TLorentzVector f_ppim   = fProton + fPim;
  TLorentzVector f_pippim = fPip + fPim;




/// ////////////////////////////////////////////////////////////////////////////////////////
/// fill output files with tagged events
/// ////////////////////////////////////////////////////////////////////////////////////////

/*
 if(e_count == 1 && W < 2){

   outputFile_electrons_elastic << NEVENT << endl;

 }

 if(e_count == 1 && p_count == 1 && W < 2){

   outputFile_ep_elastic << NEVENT << endl;

 }

 if(e_count == 1 && pip_count == 1){

   outputFile_epip_missing_neutron << NEVENT << endl;

 }
*/

/// ////////////////////////////////////////////////////////////////////////////////////////
/// fill the histograms
/// ////////////////////////////////////////////////////////////////////////////////////////


if(e_count > 0){

  hist_run_number->Fill(NRUN);
  hist_number_of_events->Fill(NEVENT);
  hist_event_type->Fill(TYPE);
  hist_trigger->Fill(TRG);
  hist_helicity->Fill(helicity);
  hist_eventtime->Fill(EVNTime);
  hist_faraday_cup->Fill(BCG);
  hist_event_starttime->Fill(STTime);
  hist_RF_time->Fill(RFTime);

  // fill kinematics:

  if(W > 0) hist_W->Fill(W);
  if(e_FTOF_sec[0] == 1  && W > 0) hist_W_sec_01->Fill(W);
  if(e_FTOF_sec[0] == 2  && W > 0) hist_W_sec_02->Fill(W);
  if(e_FTOF_sec[0] == 3  && W > 0) hist_W_sec_03->Fill(W);
  if(e_FTOF_sec[0] == 4  && W > 0) hist_W_sec_04->Fill(W);
  if(e_FTOF_sec[0] == 5  && W > 0) hist_W_sec_05->Fill(W);
  if(e_FTOF_sec[0] == 6  && W > 0) hist_W_sec_06->Fill(W);
  if(Q2 > 0) hist_Q2->Fill(Q2);
  if(W > 0 && Q2 > 0) hist_Q2_vs_W->Fill(W, Q2);
  if(e_FTOF_sec[0] == 1 && Q2 > 0 && W > 0) hist_Q2_vs_W_sec_01->Fill(W, Q2);
  if(e_FTOF_sec[0] == 2 && Q2 > 0 && W > 0) hist_Q2_vs_W_sec_02->Fill(W, Q2);
  if(e_FTOF_sec[0] == 3 && Q2 > 0 && W > 0) hist_Q2_vs_W_sec_03->Fill(W, Q2);
  if(e_FTOF_sec[0] == 4 && Q2 > 0 && W > 0) hist_Q2_vs_W_sec_04->Fill(W, Q2);
  if(e_FTOF_sec[0] == 5 && Q2 > 0 && W > 0) hist_Q2_vs_W_sec_05->Fill(W, Q2);
  if(e_FTOF_sec[0] == 6 && Q2 > 0 && W > 0) hist_Q2_vs_W_sec_06->Fill(W, Q2);
  if(p4_ele[0].Phi() != 0 && W > 0) hist_W_vs_phi->Fill(p4_ele[0].Phi()*180/Pival, W);
  if(x > 0) hist_x->Fill(x);
  if(y > 0) hist_y->Fill(y);
  if(nu > 0) hist_nu->Fill(nu);
  if(t_pi0 > 0) hist_t_pi0->Fill(-t_pi0);
  if(t_pip > 0) hist_t_pip->Fill(-t_pip);
  if(t_pim > 0) hist_t_pim->Fill(-t_pim);
  if(cmphi_p != 0 && p_count > 0) hist_cmphi_p->Fill(cmphi_p*180/Pival);
  if(cmcostheta_p != 0 && p_count > 0) hist_cmcostheta_p->Fill(cmcostheta_p);
  if(cmphi_pi0 != 0 ) hist_cmphi_pi0->Fill(cmphi_pi0*180/Pival);
  if(cmcostheta_pi0 != 0 ) hist_cmcostheta_pi0->Fill(cmcostheta_pi0);
  if(cmphi_pip != 0 && pip_count > 0) hist_cmphi_pip->Fill(cmphi_pip*180/Pival);
  if(cmcostheta_pip != 0 && pip_count > 0) hist_cmcostheta_pip->Fill(cmcostheta_pip);
  if(cmphi_pim != 0 && pim_count > 0) hist_cmphi_pim->Fill(cmphi_pim*180/Pival);
  if(cmcostheta_pim != 0 && pim_count > 0) hist_cmcostheta_pim->Fill(cmcostheta_pim);

  // fill particle distributions:

  for(Int_t i = 0; i < BUFFER; i++){

    if(part_status[i] > 0) hist_status->Fill(part_status[i]);

    if(part_p[i] > 0) hist_particles_p->Fill(part_p[i]);
    if(part_theta[i] > 0) hist_particles_theta->Fill(part_theta[i]*180/Pival);
    if(part_phi[i] != 0) hist_particles_phi->Fill(part_phi[i]*180/Pival);
    if(part_theta[i] > 0 && part_p[i] > 0) hist_particles_p_vs_theta->Fill(part_theta[i]*180/Pival, part_p[i]);
    if(part_phi[i] != 0 && part_p[i] > 0) hist_particles_p_vs_phi->Fill(part_phi[i]*180/Pival, part_p[i]);
    if(part_theta[i] > 0 && part_phi[i] !=0) hist_particles_theta_vs_phi->Fill(part_phi[i]*180/Pival, part_theta[i]*180/Pival);

    if(part_charge[i] == +1){
      if(part_p[i] > 0) hist_positives_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_positives_theta->Fill(part_theta[i]*180/Pival);
      if(part_phi[i] != 0) hist_positives_phi->Fill(part_phi[i]*180/Pival);
      if(part_theta[i] > 0 && part_p[i] > 0) hist_positives_p_vs_theta->Fill(part_theta[i]*180/Pival, part_p[i]);
      if(part_phi[i] != 0 && part_p[i] > 0) hist_positives_p_vs_phi->Fill(part_phi[i]*180/Pival, part_p[i]);
      if(part_theta[i] > 0 && part_phi[i] !=0) hist_positives_theta_vs_phi->Fill(part_phi[i]*180/Pival, part_theta[i]*180/Pival);
    }
    if(part_charge[i] == -1){
      if(part_p[i] > 0) hist_negatives_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_negatives_theta->Fill(part_theta[i]*180/Pival);
      if(part_phi[i] != 0) hist_negatives_phi->Fill(part_phi[i]*180/Pival);
      if(part_theta[i] > 0 && part_p[i] > 0) hist_negatives_p_vs_theta->Fill(part_theta[i]*180/Pival, part_p[i]);
      if(part_phi[i] != 0 && part_p[i] > 0) hist_negatives_p_vs_phi->Fill(part_phi[i]*180/Pival, part_p[i]);
      if(part_theta[i] > 0 && part_phi[i] !=0) hist_negatives_theta_vs_phi->Fill(part_phi[i]*180/Pival, part_theta[i]*180/Pival);
    }
  }

  if(p4_ele[0].P() > 0) hist_electron_p->Fill(p4_ele[0].P());
  if(p4_ele[0].Theta() > 0) hist_electron_theta->Fill(p4_ele[0].Theta()*180/Pival);
  if(p4_ele[0].Phi() != 0) hist_electron_phi->Fill(p4_ele[0].Phi()*180/Pival);
  if(p4_ele[0].P() > 0) hist_electron_p_vs_theta->Fill(p4_ele[0].Theta()*180/Pival, p4_ele[0].P());
  if(p4_ele[0].P() > 0) hist_electron_p_vs_phi->Fill(p4_ele[0].Phi()*180/Pival, p4_ele[0].P());
  if(p4_ele[0].Theta() > 0 && p4_ele[0].Phi() != 0) hist_electron_theta_vs_phi->Fill(p4_ele[0].Phi()*180/Pival, p4_ele[0].Theta()*180/Pival);

  if(p4_prot[0].P() > 0) hist_proton_p->Fill(p4_prot[0].P());
  if(p4_prot[0].Theta() > 0) hist_proton_theta->Fill(p4_prot[0].Theta()*180/Pival);
  if(p4_prot[0].Phi() != 0) hist_proton_phi->Fill(p4_prot[0].Phi()*180/Pival);
  if(p4_prot[0].P() > 0) hist_proton_p_vs_theta->Fill(p4_prot[0].Theta()*180/Pival, p4_prot[0].P());
  if(p4_prot[0].P() > 0) hist_proton_p_vs_phi->Fill(p4_prot[0].Phi()*180/Pival, p4_prot[0].P());
  if(p4_prot[0].Theta() > 0 && p4_prot[0].Phi() != 0) hist_proton_theta_vs_phi->Fill(p4_prot[0].Phi()*180/Pival, p4_prot[0].Theta()*180/Pival);

  if(p4_neutr[0].P() > 0) hist_neutron_p->Fill(p4_neutr[0].P());
  if(p4_neutr[0].Theta() > 0) hist_neutron_theta->Fill(p4_neutr[0].Theta()*180/Pival);
  if(p4_neutr[0].Phi() != 0) hist_neutron_phi->Fill(p4_neutr[0].Phi()*180/Pival);
  if(p4_neutr[0].P() > 0) hist_neutron_p_vs_theta->Fill(p4_neutr[0].Theta()*180/Pival, p4_neutr[0].P());
  if(p4_neutr[0].P() > 0) hist_neutron_p_vs_phi->Fill(p4_neutr[0].Phi()*180/Pival, p4_neutr[0].P());
  if(p4_neutr[0].Theta() > 0 && p4_neutr[0].Phi() != 0) hist_neutron_theta_vs_phi->Fill(p4_neutr[0].Phi()*180/Pival, p4_neutr[0].Theta()*180/Pival);

  if(p4_pip[0].P() > 0) hist_pip_p->Fill(p4_pip[0].P());
  if(p4_pip[0].Theta() > 0) hist_pip_theta->Fill(p4_pip[0].Theta()*180/Pival);
  if(p4_pip[0].Phi() != 0) hist_pip_phi->Fill(p4_pip[0].Phi()*180/Pival);
  if(p4_pip[0].P() > 0) hist_pip_p_vs_theta->Fill(p4_pip[0].Theta()*180/Pival, p4_pip[0].P());
  if(p4_pip[0].P() > 0) hist_pip_p_vs_phi->Fill(p4_pip[0].Phi()*180/Pival, p4_pip[0].P());
  if(p4_pip[0].Theta() > 0 && p4_pip[0].Phi() != 0) hist_pip_theta_vs_phi->Fill(p4_pip[0].Phi()*180/Pival, p4_pip[0].Theta()*180/Pival);

  if(p4_pim[0].P() > 0) hist_pim_p->Fill(p4_pim[0].P());
  if(p4_pim[0].Theta() > 0) hist_pim_theta->Fill(p4_pim[0].Theta()*180/Pival);
  if(p4_pim[0].Phi() != 0) hist_pim_phi->Fill(p4_pim[0].Phi()*180/Pival);
  if(p4_pim[0].P() > 0) hist_pim_p_vs_theta->Fill(p4_pim[0].Theta()*180/Pival, p4_pim[0].P());
  if(p4_pim[0].P() > 0) hist_pim_p_vs_phi->Fill(p4_pim[0].Phi()*180/Pival, p4_pim[0].P());
  if(p4_pim[0].Theta() > 0 && p4_pim[0].Phi() != 0) hist_pim_theta_vs_phi->Fill(p4_pim[0].Phi()*180/Pival, p4_pim[0].Theta()*180/Pival);

  if(p4_Kp[0].P() > 0) hist_Kp_p->Fill(p4_Kp[0].P());
  if(p4_Kp[0].Theta() > 0) hist_Kp_theta->Fill(p4_Kp[0].Theta()*180/Pival);
  if(p4_Kp[0].Phi() != 0) hist_Kp_phi->Fill(p4_Kp[0].Phi()*180/Pival);
  if(p4_Kp[0].P() > 0) hist_Kp_p_vs_theta->Fill(p4_Kp[0].Theta()*180/Pival, p4_Kp[0].P());
  if(p4_Kp[0].P() > 0) hist_Kp_p_vs_phi->Fill(p4_Kp[0].Phi()*180/Pival, p4_Kp[0].P());
  if(p4_Kp[0].Theta() > 0 && p4_Kp[0].Phi() != 0) hist_Kp_theta_vs_phi->Fill(p4_Kp[0].Phi()*180/Pival, p4_Kp[0].Theta()*180/Pival);

  if(p4_Km[0].P() > 0) hist_Km_p->Fill(p4_Km[0].P());
  if(p4_Km[0].Theta() > 0) hist_Km_theta->Fill(p4_Km[0].Theta()*180/Pival);
  if(p4_Km[0].Phi() != 0) hist_Km_phi->Fill(p4_Km[0].Phi()*180/Pival);
  if(p4_Km[0].P() > 0) hist_Km_p_vs_theta->Fill(p4_Km[0].Theta()*180/Pival, p4_Km[0].P());
  if(p4_Km[0].P() > 0) hist_Km_p_vs_phi->Fill(p4_Km[0].Phi()*180/Pival, p4_Km[0].P());
  if(p4_Km[0].Theta() > 0 && p4_Km[0].Phi() != 0) hist_Km_theta_vs_phi->Fill(p4_Km[0].Phi()*180/Pival, p4_Km[0].Theta()*180/Pival);

  if(p4_phot[0].P() > 0) hist_photon_p->Fill(p4_phot[0].P());
  if(p4_phot[0].Theta() > 0) hist_photon_theta->Fill(p4_phot[0].Theta()*180/Pival);
  if(p4_phot[0].Phi() != 0) hist_photon_phi->Fill(p4_phot[0].Phi()*180/Pival);
  if(p4_phot[0].P() > 0) hist_photon_p_vs_theta->Fill(p4_phot[0].Theta()*180/Pival, p4_phot[0].P());
  if(p4_phot[0].P() > 0) hist_photon_p_vs_phi->Fill(p4_phot[0].Phi()*180/Pival, p4_phot[0].P());
  if(p4_phot[0].Theta() > 0 && p4_phot[0].Phi() != 0) hist_photon_theta_vs_phi->Fill(p4_phot[0].Phi()*180/Pival, p4_phot[0].Theta()*180/Pival);



  for(Int_t i = 0; i < BUFFER; i++){
  
    if(p4_ele[i].P() > 0) hist_all_electron_p->Fill(p4_ele[i].P());
    if(p4_ele[i].Theta() > 0) hist_all_electron_theta->Fill(p4_ele[i].Theta()*180/Pival);
    if(p4_ele[i].Phi() != 0) hist_all_electron_phi->Fill(p4_ele[i].Phi()*180/Pival);
    if(p4_ele[i].P() > 0) hist_all_electron_p_vs_theta->Fill(p4_ele[i].Theta()*180/Pival, p4_ele[i].P());
    if(p4_ele[i].P() > 0) hist_all_electron_p_vs_phi->Fill(p4_ele[i].Phi()*180/Pival, p4_ele[i].P());
    if(p4_ele[i].Theta() > 0 && p4_ele[i].Phi() != 0) hist_all_electron_theta_vs_phi->Fill(p4_ele[i].Phi()*180/Pival, p4_ele[i].Theta()*180/Pival);

    if(p4_prot[i].P() > 0) hist_all_proton_p->Fill(p4_prot[i].P());
    if(p4_prot[i].Theta() > 0) hist_all_proton_theta->Fill(p4_prot[i].Theta()*180/Pival);
    if(p4_prot[i].Phi() != 0) hist_all_proton_phi->Fill(p4_prot[i].Phi()*180/Pival);
    if(p4_prot[i].P() > 0) hist_all_proton_p_vs_theta->Fill(p4_prot[i].Theta()*180/Pival, p4_prot[i].P());
    if(p4_prot[i].P() > 0) hist_all_proton_p_vs_phi->Fill(p4_prot[i].Phi()*180/Pival, p4_prot[i].P());
    if(p4_prot[i].Theta() > 0 && p4_prot[i].Phi() != 0) hist_all_proton_theta_vs_phi->Fill(p4_prot[i].Phi()*180/Pival, p4_prot[i].Theta()*180/Pival);

    if(p4_neutr[i].P() > 0) hist_all_neutron_p->Fill(p4_neutr[i].P());
    if(p4_neutr[i].Theta() > 0) hist_all_neutron_theta->Fill(p4_neutr[i].Theta()*180/Pival);
    if(p4_neutr[i].Phi() != 0) hist_all_neutron_phi->Fill(p4_neutr[i].Phi()*180/Pival);
    if(p4_neutr[i].P() > 0) hist_all_neutron_p_vs_theta->Fill(p4_neutr[i].Theta()*180/Pival, p4_neutr[i].P());
    if(p4_neutr[i].P() > 0) hist_all_neutron_p_vs_phi->Fill(p4_neutr[i].Phi()*180/Pival, p4_neutr[i].P());
    if(p4_neutr[i].Theta() > 0 && p4_neutr[i].Phi() != 0) hist_all_neutron_theta_vs_phi->Fill(p4_neutr[i].Phi()*180/Pival, p4_neutr[i].Theta()*180/Pival);

    if(p4_pip[i].P() > 0) hist_all_pip_p->Fill(p4_pip[i].P());
    if(p4_pip[i].Theta() > 0) hist_all_pip_theta->Fill(p4_pip[i].Theta()*180/Pival);
    if(p4_pip[i].Phi() != 0) hist_all_pip_phi->Fill(p4_pip[i].Phi()*180/Pival);
    if(p4_pip[i].P() > 0) hist_all_pip_p_vs_theta->Fill(p4_pip[i].Theta()*180/Pival, p4_pip[i].P());
    if(p4_pip[i].P() > 0) hist_all_pip_p_vs_phi->Fill(p4_pip[i].Phi()*180/Pival, p4_pip[i].P());
    if(p4_pip[i].Theta() > 0 && p4_pip[i].Phi() != 0) hist_all_pip_theta_vs_phi->Fill(p4_pip[i].Phi()*180/Pival, p4_pip[i].Theta()*180/Pival);

    if(p4_pim[i].P() > 0) hist_all_pim_p->Fill(p4_pim[i].P());
    if(p4_pim[i].Theta() > 0) hist_all_pim_theta->Fill(p4_pim[i].Theta()*180/Pival);
    if(p4_pim[i].Phi() != 0) hist_all_pim_phi->Fill(p4_pim[i].Phi()*180/Pival);
    if(p4_pim[i].P() > 0) hist_all_pim_p_vs_theta->Fill(p4_pim[i].Theta()*180/Pival, p4_pim[i].P());
    if(p4_pim[i].P() > 0) hist_all_pim_p_vs_phi->Fill(p4_pim[i].Phi()*180/Pival, p4_pim[i].P());
    if(p4_pim[i].Theta() > 0 && p4_pim[i].Phi() != 0) hist_all_pim_theta_vs_phi->Fill(p4_pim[i].Phi()*180/Pival, p4_pim[i].Theta()*180/Pival);

    if(p4_Kp[i].P() > 0) hist_all_Kp_p->Fill(p4_Kp[i].P());
    if(p4_Kp[i].Theta() > 0) hist_all_Kp_theta->Fill(p4_Kp[i].Theta()*180/Pival);
    if(p4_Kp[i].Phi() != 0) hist_all_Kp_phi->Fill(p4_Kp[i].Phi()*180/Pival);
    if(p4_Kp[i].P() > 0) hist_all_Kp_p_vs_theta->Fill(p4_Kp[i].Theta()*180/Pival, p4_Kp[i].P());
    if(p4_Kp[i].P() > 0) hist_all_Kp_p_vs_phi->Fill(p4_Kp[i].Phi()*180/Pival, p4_Kp[i].P());
    if(p4_Kp[i].Theta() > 0 && p4_Kp[i].Phi() != 0) hist_all_Kp_theta_vs_phi->Fill(p4_Kp[i].Phi()*180/Pival, p4_Kp[i].Theta()*180/Pival);

    if(p4_Km[i].P() > 0) hist_all_Km_p->Fill(p4_Km[i].P());
    if(p4_Km[i].Theta() > 0) hist_all_Km_theta->Fill(p4_Km[i].Theta()*180/Pival);
    if(p4_Km[i].Phi() != 0) hist_all_Km_phi->Fill(p4_Km[i].Phi()*180/Pival);
    if(p4_Km[i].P() > 0) hist_all_Km_p_vs_theta->Fill(p4_Km[i].Theta()*180/Pival, p4_Km[i].P());
    if(p4_Km[i].P() > 0) hist_all_Km_p_vs_phi->Fill(p4_Km[i].Phi()*180/Pival, p4_Km[i].P());
    if(p4_Km[i].Theta() > 0 && p4_Km[i].Phi() != 0) hist_all_Km_theta_vs_phi->Fill(p4_Km[i].Phi()*180/Pival, p4_Km[i].Theta()*180/Pival);

    if(p4_phot[i].P() > 0) hist_all_photon_p->Fill(p4_phot[i].P());
    if(p4_phot[i].Theta() > 0) hist_all_photon_theta->Fill(p4_phot[i].Theta()*180/Pival);
    if(p4_phot[i].Phi() != 0) hist_all_photon_phi->Fill(p4_phot[i].Phi()*180/Pival);
    if(p4_phot[i].P() > 0) hist_all_photon_p_vs_theta->Fill(p4_phot[i].Theta()*180/Pival, p4_phot[i].P());
    if(p4_phot[i].P() > 0) hist_all_photon_p_vs_phi->Fill(p4_phot[i].Phi()*180/Pival, p4_phot[i].P());
    if(p4_phot[i].Theta() > 0 && p4_phot[i].Phi() != 0) hist_all_photon_theta_vs_phi->Fill(p4_phot[i].Phi()*180/Pival, p4_phot[i].Theta()*180/Pival);

  }


  for(Int_t i = 0; i < BUFFER; i++){
    if(partLUND_pid[i] != 0) hist_LUND_pid->Fill(partLUND_pid[i]);
  }


  /// ///////////////////////////////////////////////////
  /// count generated particles of each type

  for(Int_t i = 0; i < BUFFER; i++){
    if(partLUND_pid[i] == 11)  { e_MCcount   += 1;}  
    if(partLUND_pid[i] == 2212){ p_MCcount   += 1;}  
    if(partLUND_pid[i] == 2112){ n_MCcount   += 1;}  
    if(partLUND_pid[i] == 211) { pip_MCcount += 1;}  
    if(partLUND_pid[i] == -211){ pim_MCcount += 1;}  
    if(partLUND_pid[i] == 321) { Kp_MCcount  += 1;}  
    if(partLUND_pid[i] == -321){ Km_MCcount  += 1;} 
  }

  int e_RECcount = 0;
  int p_RECcount = 0;
  int n_RECcount = 0;
  int pip_RECcount = 0;
  int pim_RECcount = 0;
  int Kp_RECcount = 0;
  int Km_RECcount = 0;

  for(Int_t i = 0; i < BUFFER; i++){
    if(part_pid[i] == 11)  { e_RECcount   += 1;}  
    if(part_pid[i] == 2212){ p_RECcount   += 1;}  
    if(part_pid[i] == 2112){ n_RECcount   += 1;}  
    if(part_pid[i] == 211) { pip_RECcount += 1;}  
    if(part_pid[i] == -211){ pim_RECcount += 1;}  
    if(part_pid[i] == 321) { Kp_RECcount  += 1;}  
    if(part_pid[i] == -321){ Km_RECcount  += 1;} 
  }
  

  for(Int_t i = 0; i < BUFFER; i++){
  
  if(part_theta[i]*180/Pival > 5 && part_theta[i]*180/Pival < 35){

    if(part_pid[i] == 11){
      if(part_p[i] > 0)     hist_FD_all_det_electron_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_FD_all_det_electron_theta->Fill(part_theta[i]*180/Pival);
    }
    if(part_pid[i] == 2212){
      if(part_p[i] > 0)     hist_FD_all_det_proton_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_FD_all_det_proton_theta->Fill(part_theta[i]*180/Pival);
    }
    if(part_pid[i] == 2112){
      if(part_p[i] > 0)     hist_FD_all_det_neutron_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_FD_all_det_neutron_theta->Fill(part_theta[i]*180/Pival);
    }
    if(part_pid[i] == 211){
      if(part_p[i] > 0)     hist_FD_all_det_pip_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_FD_all_det_pip_theta->Fill(part_theta[i]*180/Pival);
    }
    if(part_pid[i] == -211){
      if(part_p[i] > 0)     hist_FD_all_det_pim_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_FD_all_det_pim_theta->Fill(part_theta[i]*180/Pival);
    }
    if(part_pid[i] == 321){
      if(part_p[i] > 0)     hist_FD_all_det_Kp_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_FD_all_det_Kp_theta->Fill(part_theta[i]*180/Pival);
    }
    if(part_pid[i] == -321){
      if(part_p[i] > 0)    hist_FD_all_det_Km_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_FD_all_det_Km_theta->Fill(part_theta[i]*180/Pival);
    }

    if(part_pid[i] == 11 && e_MCcount > 0){
      if(part_p[i] > 0)     hist_FD_all_corr_electron_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_FD_all_corr_electron_theta->Fill(part_theta[i]*180/Pival);
    }
    if(part_pid[i] == 2212 && p_MCcount > 0){
      if(part_p[i] > 0)     hist_FD_all_corr_proton_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_FD_all_corr_proton_theta->Fill(part_theta[i]*180/Pival);
    }
    if(part_pid[i] == 2112 && n_MCcount > 0){
      if(part_p[i] > 0)     hist_FD_all_corr_neutron_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_FD_all_corr_neutron_theta->Fill(part_theta[i]*180/Pival);
    }
    if(part_pid[i] == 211 && pip_MCcount > 0){
      if(part_p[i] > 0)     hist_FD_all_corr_pip_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_FD_all_corr_pip_theta->Fill(part_theta[i]*180/Pival);
    }
    if(part_pid[i] == -211 && pim_MCcount > 0){
      if(part_p[i] > 0)     hist_FD_all_corr_pim_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_FD_all_corr_pim_theta->Fill(part_theta[i]*180/Pival);
    }
    if(part_pid[i] == 321 && Kp_MCcount > 0){
      if(part_p[i] > 0)     hist_FD_all_corr_Kp_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_FD_all_corr_Kp_theta->Fill(part_theta[i]*180/Pival);
    }
    if(part_pid[i] == -321 && Km_MCcount > 0){
      if(part_p[i] > 0)    hist_FD_all_corr_Km_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_FD_all_corr_Km_theta->Fill(part_theta[i]*180/Pival);
    }
  }

  if(partMC_theta[i]*180/Pival > 5 && partMC_theta[i]*180/Pival < 35){

    if(partMC_pid[i] == 11){
      if(partMC_p[i] > 0)     hist_FD_all_gen_electron_p->Fill(partMC_p[i]);
      if(partMC_theta[i] > 0) hist_FD_all_gen_electron_theta->Fill(partMC_theta[i]*180/Pival);
    }
    if(partMC_pid[i] == 2212){
      if(partMC_p[i] > 0)     hist_FD_all_gen_proton_p->Fill(partMC_p[i]);
      if(partMC_theta[i] > 0) hist_FD_all_gen_proton_theta->Fill(partMC_theta[i]*180/Pival);
    }
    if(partMC_pid[i] == 2112){
      if(partMC_p[i] > 0)     hist_FD_all_gen_neutron_p->Fill(partMC_p[i]);
      if(partMC_theta[i] > 0) hist_FD_all_gen_neutron_theta->Fill(partMC_theta[i]*180/Pival);
    }
    if(partMC_pid[i] == 211){
      if(partMC_p[i] > 0)     hist_FD_all_gen_pip_p->Fill(partMC_p[i]);
      if(partMC_theta[i] > 0) hist_FD_all_gen_pip_theta->Fill(partMC_theta[i]*180/Pival);
    }
    if(partMC_pid[i] == -211){
      if(partMC_p[i] > 0)     hist_FD_all_gen_pim_p->Fill(partMC_p[i]);
      if(partMC_theta[i] > 0) hist_FD_all_gen_pim_theta->Fill(partMC_theta[i]*180/Pival);
    }
    if(partMC_pid[i] == 321){
      if(partMC_p[i] > 0)     hist_FD_all_gen_Kp_p->Fill(partMC_p[i]);
      if(partMC_theta[i] > 0) hist_FD_all_gen_Kp_theta->Fill(partMC_theta[i]*180/Pival);
    }
    if(partMC_pid[i] == -321){
      if(partMC_p[i] > 0)    hist_FD_all_gen_Km_p->Fill(partMC_p[i]);
      if(partMC_theta[i] > 0) hist_FD_all_gen_Km_theta->Fill(partMC_theta[i]*180/Pival);
    }
  }

  if(part_theta[i]*180/Pival > 35){

    if(part_pid[i] == 11){
      if(part_p[i] > 0)     hist_CD_all_det_electron_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_CD_all_det_electron_theta->Fill(part_theta[i]*180/Pival);
    }
    if(part_pid[i] == 2212){
      if(part_p[i] > 0)     hist_CD_all_det_proton_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_CD_all_det_proton_theta->Fill(part_theta[i]*180/Pival);
    }
    if(part_pid[i] == 2112){
      if(part_p[i] > 0)     hist_CD_all_det_neutron_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_CD_all_det_neutron_theta->Fill(part_theta[i]*180/Pival);
    }
    if(part_pid[i] == 211){
      if(part_p[i] > 0)     hist_CD_all_det_pip_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_CD_all_det_pip_theta->Fill(part_theta[i]*180/Pival);
    }
    if(part_pid[i] == -211){
      if(part_p[i] > 0)     hist_CD_all_det_pim_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_CD_all_det_pim_theta->Fill(part_theta[i]*180/Pival);
    }
    if(part_pid[i] == 321){
      if(part_p[i] > 0)     hist_CD_all_det_Kp_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_CD_all_det_Kp_theta->Fill(part_theta[i]*180/Pival);
    }
    if(part_pid[i] == -321){
      if(part_p[i] > 0)    hist_CD_all_det_Km_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_CD_all_det_Km_theta->Fill(part_theta[i]*180/Pival);
    }

    if(part_pid[i] == 11 && e_MCcount > 0){
      if(part_p[i] > 0)     hist_CD_all_corr_electron_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_CD_all_corr_electron_theta->Fill(part_theta[i]*180/Pival);
    }
    if(part_pid[i] == 2212 && p_MCcount > 0){
      if(part_p[i] > 0)     hist_CD_all_corr_proton_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_CD_all_corr_proton_theta->Fill(part_theta[i]*180/Pival);
    }
    if(part_pid[i] == 2112 && n_MCcount > 0){
      if(part_p[i] > 0)     hist_CD_all_corr_neutron_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_CD_all_corr_neutron_theta->Fill(part_theta[i]*180/Pival);
    }
    if(part_pid[i] == 211 && pip_MCcount > 0){
      if(part_p[i] > 0)     hist_CD_all_corr_pip_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_CD_all_corr_pip_theta->Fill(part_theta[i]*180/Pival);
    }
    if(part_pid[i] == -211 && pim_MCcount > 0){
      if(part_p[i] > 0)     hist_CD_all_corr_pim_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_CD_all_corr_pim_theta->Fill(part_theta[i]*180/Pival);
    }
    if(part_pid[i] == 321 && Kp_MCcount > 0){
      if(part_p[i] > 0)     hist_CD_all_corr_Kp_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_CD_all_corr_Kp_theta->Fill(part_theta[i]*180/Pival);
    }
    if(part_pid[i] == -321 && Km_MCcount > 0){
      if(part_p[i] > 0)    hist_CD_all_corr_Km_p->Fill(part_p[i]);
      if(part_theta[i] > 0) hist_CD_all_corr_Km_theta->Fill(part_theta[i]*180/Pival);
    }

  }
  if(partMC_theta[i]*180/Pival > 35){

    if(partMC_pid[i] == 11){
      if(partMC_p[i] > 0)     hist_CD_all_gen_electron_p->Fill(partMC_p[i]);
      if(partMC_theta[i] > 0) hist_CD_all_gen_electron_theta->Fill(partMC_theta[i]*180/Pival);
    }
    if(partMC_pid[i] == 2212){
      if(partMC_p[i] > 0)     hist_CD_all_gen_proton_p->Fill(partMC_p[i]);
      if(partMC_theta[i] > 0) hist_CD_all_gen_proton_theta->Fill(partMC_theta[i]*180/Pival);
    }
    if(partMC_pid[i] == 2112){
      if(partMC_p[i] > 0)     hist_CD_all_gen_neutron_p->Fill(partMC_p[i]);
      if(partMC_theta[i] > 0) hist_CD_all_gen_neutron_theta->Fill(partMC_theta[i]*180/Pival);
    }
    if(partMC_pid[i] == 211){
      if(partMC_p[i] > 0)     hist_CD_all_gen_pip_p->Fill(partMC_p[i]);
      if(partMC_theta[i] > 0) hist_CD_all_gen_pip_theta->Fill(partMC_theta[i]*180/Pival);
    }
    if(partMC_pid[i] == -211){
      if(partMC_p[i] > 0)     hist_CD_all_gen_pim_p->Fill(partMC_p[i]);
      if(partMC_theta[i] > 0) hist_CD_all_gen_pim_theta->Fill(partMC_theta[i]*180/Pival);
    }
    if(partMC_pid[i] == 321){
      if(partMC_p[i] > 0)     hist_CD_all_gen_Kp_p->Fill(partMC_p[i]);
      if(partMC_theta[i] > 0) hist_CD_all_gen_Kp_theta->Fill(partMC_theta[i]*180/Pival);
    }
    if(partMC_pid[i] == -321){
      if(partMC_p[i] > 0)    hist_CD_all_gen_Km_p->Fill(partMC_p[i]);
      if(partMC_theta[i] > 0) hist_CD_all_gen_Km_theta->Fill(partMC_theta[i]*180/Pival);
    }

  }

  }



  /// ///////////////////////////////////////////////////////
  /// Calorimeter energy depositions
  ///
  /// e_FTOF_sec[i] e_PCAL_sec[i] pip_FTOF_sec[i] pim_FTOF_sec[i]

  for(Int_t i = 0; i < BUFFER; i++){

    if(part_Cal_energy_total[i] > 0 && FD_eid_all_check[i]) hist_CAL_Edep_electron->Fill(part_Cal_energy_total[i]);
    if(part_Cal_energy_total[i] > 0 && FD_protid_all_check[i]) hist_CAL_Edep_proton->Fill(part_Cal_energy_total[i]);
    if(part_Cal_energy_total[i] > 0 && FD_pipid_all_check[i]) hist_CAL_Edep_pip->Fill(part_Cal_energy_total[i]);
    if(part_Cal_energy_total[i] > 0 && (FD_pimid_default_PID_check[i] && FD_pimid_beta_check[i])) hist_CAL_Edep_pim->Fill(part_Cal_energy_total[i]);
    if(part_Cal_energy_total[i] > 0 && FD_pipid_all_check[i] && (e_FTOF_sec[i] != pip_FTOF_sec[i])) hist_CAL_Edep_pip_diffsecele->Fill(part_Cal_energy_total[i]);
    if(part_Cal_energy_total[i] > 0 && (FD_pimid_default_PID_check[i] && FD_pimid_beta_check[i]) && (e_FTOF_sec[i] != pim_FTOF_sec[i])) hist_CAL_Edep_pim_diffsecele->Fill(part_Cal_energy_total[i]);

    if(part_Cal_PCAL_energy[i] > 0 && FD_eid_all_check[i]) hist_PCAL_Edep_electron->Fill(part_Cal_PCAL_energy[i]);
    if(part_Cal_PCAL_energy[i] > 0 && FD_protid_all_check[i]) hist_PCAL_Edep_proton->Fill(part_Cal_PCAL_energy[i]);
    if(part_Cal_PCAL_energy[i] > 0 && FD_pipid_all_check[i]) hist_PCAL_Edep_pip->Fill(part_Cal_PCAL_energy[i]);
    if(part_Cal_PCAL_energy[i] > 0 && (FD_pimid_default_PID_check[i] && FD_pimid_beta_check[i])) hist_PCAL_Edep_pim->Fill(part_Cal_PCAL_energy[i]);
    if(part_Cal_PCAL_energy[i] > 0 && FD_pipid_all_check[i] && (e_FTOF_sec[i] != pip_FTOF_sec[i])) hist_PCAL_Edep_pip_diffsecele->Fill(part_Cal_PCAL_energy[i]);
    if(part_Cal_PCAL_energy[i] > 0 && (FD_pimid_default_PID_check[i] && FD_pimid_beta_check[i]) && (e_FTOF_sec[i] != pim_FTOF_sec[i])) hist_PCAL_Edep_pim_diffsecele->Fill(part_Cal_PCAL_energy[i]);

    if(part_Cal_ECin_energy[i] > 0 && FD_eid_all_check[i]) hist_ECin_Edep_electron->Fill(part_Cal_ECin_energy[i]);
    if(part_Cal_ECin_energy[i] > 0 && FD_protid_all_check[i]) hist_ECin_Edep_proton->Fill(part_Cal_ECin_energy[i]);
    if(part_Cal_ECin_energy[i] > 0 && FD_pipid_all_check[i]) hist_ECin_Edep_pip->Fill(part_Cal_ECin_energy[i]);
    if(part_Cal_ECin_energy[i] > 0 && (FD_pimid_default_PID_check[i] && FD_pimid_beta_check[i])) hist_ECin_Edep_pim->Fill(part_Cal_ECin_energy[i]);
    if(part_Cal_ECin_energy[i] > 0 && FD_pipid_all_check[i] && (e_FTOF_sec[i] != pip_FTOF_sec[i])) hist_ECin_Edep_pip_diffsecele->Fill(part_Cal_ECin_energy[i]);
    if(part_Cal_ECin_energy[i] > 0 && (FD_pimid_default_PID_check[i] && FD_pimid_beta_check[i]) && (e_FTOF_sec[i] != pim_FTOF_sec[i])) hist_ECin_Edep_pim_diffsecele->Fill(part_Cal_ECin_energy[i]);

    if(part_Cal_ECout_energy[i] > 0 && FD_eid_all_check[i]) hist_ECout_Edep_electron->Fill(part_Cal_ECout_energy[i]);
    if(part_Cal_ECout_energy[i] > 0 && FD_protid_all_check[i]) hist_ECout_Edep_proton->Fill(part_Cal_ECout_energy[i]);
    if(part_Cal_ECout_energy[i] > 0 && FD_pipid_all_check[i]) hist_ECout_Edep_pip->Fill(part_Cal_ECout_energy[i]);
    if(part_Cal_ECout_energy[i] > 0 && (FD_pimid_default_PID_check[i] && FD_pimid_beta_check[i])) hist_ECout_Edep_pim->Fill(part_Cal_ECout_energy[i]);
    if(part_Cal_ECout_energy[i] > 0 && FD_pipid_all_check[i] && (e_FTOF_sec[i] != pip_FTOF_sec[i])) hist_ECout_Edep_pip_diffsecele->Fill(part_Cal_ECout_energy[i]);
    if(part_Cal_ECout_energy[i] > 0 && (FD_pimid_default_PID_check[i] && FD_pimid_beta_check[i]) && (e_FTOF_sec[i] != pim_FTOF_sec[i])) hist_ECout_Edep_pim_diffsecele->Fill(part_Cal_ECout_energy[i]);

    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && FD_eid_all_check[i]) hist_ECAL_Edep_electron->Fill(part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && FD_protid_all_check[i]) hist_ECAL_Edep_proton->Fill(part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && FD_pipid_all_check[i]) hist_ECAL_Edep_pip->Fill(part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && (FD_pimid_default_PID_check[i] && FD_pimid_beta_check[i])) hist_ECAL_Edep_pim->Fill(part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && FD_pipid_all_check[i] && (e_FTOF_sec[i] != pip_FTOF_sec[i])) hist_ECAL_Edep_pip_diffsecele->Fill(part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && (FD_pimid_default_PID_check[i] && FD_pimid_beta_check[i]) && (e_FTOF_sec[i] != pim_FTOF_sec[i])) hist_ECAL_Edep_pim_diffsecele->Fill(part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i]);

  }

  
  /// //////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// FT gamma clusters

  if(p4_phot[0].Theta()*180/Pival > 2.5 && p4_phot[0].Theta()*180/Pival < 4.5){

    if(p4_phot[0].P() > 0) hist_FT_photon_E->Fill(p4_phot[0].P()); 
    if(p4_phot[0].P() > 0) hist_FT_photon_theta->Fill(p4_phot[0].Theta()*180/Pival);
    if(p4_phot[0].P() > 0) hist_FT_photon_E_vs_theta->Fill(p4_phot[0].Theta()*180/Pival, p4_phot[0].P());
    
    if(pi0.P() > 0 && pi0_g_ind[select_pi0] > -1) hist_FT_photon_E_pi0->Fill(p4_phot[pi0_g_ind[select_pi0]].P()); 
    if(pi0.P() > 0 && pi0_g_ind[select_pi0] > -1) hist_FT_photon_theta_pi0->Fill(p4_phot[pi0_g_ind[select_pi0]].Theta()*180/Pival);
    if(pi0.P() > 0 && pi0_g_ind[select_pi0] > -1) hist_FT_photon_E_vs_theta_pi0->Fill(p4_phot[pi0_g_ind[select_pi0]].Theta()*180/Pival, p4_phot[0].P());

  }



  /// ////////////////////////////////////////////////////////

  // fill plots for neutrals

  for(Int_t i = 0; i <= 27; i++){

    if(alpha_gg[i]*180/Pival > 1.5 && neutral_iter[i].E() > 0.5 && alpha_gg[i]*180/Pival > (7.0 - 1.5 * neutral_iter[i].E())){
      if(mass_neutral_iter[i] > 0.0) hist_neutral_mass->Fill(mass_neutral_iter[i]);
      if(mass_neutral_iter2[i] != 0) hist_neutral_mass2->Fill(mass_neutral_iter2[i]);
      if(neutral_iter[i].E() > 0) hist_gg_openingangle_vs_Egg->Fill(neutral_iter[i].E(), alpha_gg[i]*180/Pival);
    }
  }

  if(pi0.P() > 0)     hist_pi0_p->Fill(pi0.P());
  if(pi0.Theta() > 0) hist_pi0_theta->Fill(pi0.Theta()*180/Pival);
  if(pi0.Phi() != 0)  hist_pi0_phi->Fill(pi0.Phi()*180/Pival);
  if(pi0.M() > 0)     hist_pi0_mass->Fill(pi0.M());
  if(pi0.M2() != 0)   hist_pi0_mass2->Fill(pi0.M2());

  if(eta.P() > 0)     hist_eta_p->Fill(eta.P());
  if(eta.Theta() > 0) hist_eta_theta->Fill(eta.Theta()*180/Pival);
  if(eta.Phi() != 0)  hist_eta_phi->Fill(eta.Phi()*180/Pival);
  if(eta.M() > 0.005) hist_eta_mass->Fill(eta.M());
  if(eta.M2() != 0)   hist_eta_mass2->Fill(eta.M2());



  // fill ssector combinations for neutrals

  if(pi0_ind[0] == 1 && pi0_ind[1] == 1) hist_2photon_mass_11->Fill(pi0.M());
  if(pi0_ind[0] == 1 && pi0_ind[1] == 2) hist_2photon_mass_12->Fill(pi0.M());
  if(pi0_ind[0] == 1 && pi0_ind[1] == 3) hist_2photon_mass_13->Fill(pi0.M());
  if(pi0_ind[0] == 1 && pi0_ind[1] == 4) hist_2photon_mass_14->Fill(pi0.M());
  if(pi0_ind[0] == 1 && pi0_ind[1] == 5) hist_2photon_mass_15->Fill(pi0.M());
  if(pi0_ind[0] == 1 && pi0_ind[1] == 6) hist_2photon_mass_16->Fill(pi0.M());
  if(pi0_ind[0] == 2 && pi0_ind[1] == 2) hist_2photon_mass_22->Fill(pi0.M());
  if(pi0_ind[0] == 2 && pi0_ind[1] == 3) hist_2photon_mass_23->Fill(pi0.M());
  if(pi0_ind[0] == 2 && pi0_ind[1] == 4) hist_2photon_mass_24->Fill(pi0.M());
  if(pi0_ind[0] == 2 && pi0_ind[1] == 5) hist_2photon_mass_25->Fill(pi0.M());
  if(pi0_ind[0] == 2 && pi0_ind[1] == 6) hist_2photon_mass_26->Fill(pi0.M());
  if(pi0_ind[0] == 3 && pi0_ind[1] == 3) hist_2photon_mass_33->Fill(pi0.M());
  if(pi0_ind[0] == 3 && pi0_ind[1] == 4) hist_2photon_mass_34->Fill(pi0.M());
  if(pi0_ind[0] == 3 && pi0_ind[1] == 5) hist_2photon_mass_35->Fill(pi0.M());
  if(pi0_ind[0] == 3 && pi0_ind[1] == 6) hist_2photon_mass_36->Fill(pi0.M());
  if(pi0_ind[0] == 4 && pi0_ind[1] == 4) hist_2photon_mass_44->Fill(pi0.M());
  if(pi0_ind[0] == 4 && pi0_ind[1] == 5) hist_2photon_mass_45->Fill(pi0.M());
  if(pi0_ind[0] == 4 && pi0_ind[1] == 6) hist_2photon_mass_46->Fill(pi0.M());
  if(pi0_ind[0] == 5 && pi0_ind[1] == 5) hist_2photon_mass_55->Fill(pi0.M());
  if(pi0_ind[0] == 5 && pi0_ind[1] == 6) hist_2photon_mass_56->Fill(pi0.M());
  if(pi0_ind[0] == 6 && pi0_ind[1] == 6) hist_2photon_mass_66->Fill(pi0.M());


  if(pi0_ind[0] == 1 && pi0_ind[1] == 1) hist_2photon_mass2_11->Fill(pi0.M2());
  if(pi0_ind[0] == 1 && pi0_ind[1] == 2) hist_2photon_mass2_12->Fill(pi0.M2());
  if(pi0_ind[0] == 1 && pi0_ind[1] == 3) hist_2photon_mass2_13->Fill(pi0.M2());
  if(pi0_ind[0] == 1 && pi0_ind[1] == 4) hist_2photon_mass2_14->Fill(pi0.M2());
  if(pi0_ind[0] == 1 && pi0_ind[1] == 5) hist_2photon_mass2_15->Fill(pi0.M2());
  if(pi0_ind[0] == 1 && pi0_ind[1] == 6) hist_2photon_mass2_16->Fill(pi0.M2());
  if(pi0_ind[0] == 2 && pi0_ind[1] == 2) hist_2photon_mass2_22->Fill(pi0.M2());
  if(pi0_ind[0] == 2 && pi0_ind[1] == 3) hist_2photon_mass2_23->Fill(pi0.M2());
  if(pi0_ind[0] == 2 && pi0_ind[1] == 4) hist_2photon_mass2_24->Fill(pi0.M2());
  if(pi0_ind[0] == 2 && pi0_ind[1] == 5) hist_2photon_mass2_25->Fill(pi0.M2());
  if(pi0_ind[0] == 2 && pi0_ind[1] == 6) hist_2photon_mass2_26->Fill(pi0.M2());
  if(pi0_ind[0] == 3 && pi0_ind[1] == 3) hist_2photon_mass2_33->Fill(pi0.M2());
  if(pi0_ind[0] == 3 && pi0_ind[1] == 4) hist_2photon_mass2_34->Fill(pi0.M2());
  if(pi0_ind[0] == 3 && pi0_ind[1] == 5) hist_2photon_mass2_35->Fill(pi0.M2());
  if(pi0_ind[0] == 3 && pi0_ind[1] == 6) hist_2photon_mass2_36->Fill(pi0.M2());
  if(pi0_ind[0] == 4 && pi0_ind[1] == 4) hist_2photon_mass2_44->Fill(pi0.M2());
  if(pi0_ind[0] == 4 && pi0_ind[1] == 5) hist_2photon_mass2_45->Fill(pi0.M2());
  if(pi0_ind[0] == 4 && pi0_ind[1] == 6) hist_2photon_mass2_46->Fill(pi0.M2());
  if(pi0_ind[0] == 5 && pi0_ind[1] == 5) hist_2photon_mass2_55->Fill(pi0.M2());
  if(pi0_ind[0] == 5 && pi0_ind[1] == 6) hist_2photon_mass2_56->Fill(pi0.M2());
  if(pi0_ind[0] == 6 && pi0_ind[1] == 6) hist_2photon_mass2_66->Fill(pi0.M2());


  // fill vertex plots

  for(Int_t i = 0; i < BUFFER; i++){

    if(part_charge[i] == +1 && part_vz[i] != 0)  hist_positive_vertex->Fill(part_vz[i]);
    if(part_charge[i] == -1 && part_vz[i] != 0) hist_negative_vertex->Fill(part_vz[i]);
    if(part_charge[i] == +1 && part_vz[i] != 0 && part_theta[i] > 0)  hist_positive_vertex_vs_theta->Fill(part_theta[i]*180/Pival, part_vz[i]);
    if(part_charge[i] == -1 && part_vz[i] != 0 && part_theta[i] > 0) hist_negative_vertex_vs_theta->Fill(part_theta[i]*180/Pival, part_vz[i]);
    if(part_charge[i] == +1 && part_vz[i] != 0 && part_phi[i] != 0)  hist_positive_vertex_vs_phi->Fill(part_phi[i]*180/Pival, part_vz[i]);
    if(part_charge[i] == -1 && part_vz[i] != 0 && part_phi[i] != 0) hist_negative_vertex_vs_phi->Fill(part_phi[i]*180/Pival, part_vz[i]);
    if(part_charge[i] == +1 && part_vz[i] != 0 && part_p[i] > 0)  hist_positive_vertex_vs_p->Fill(part_p[i], part_vz[i]);
    if(part_charge[i] == -1 && part_vz[i] != 0 && part_p[i] > 0) hist_negative_vertex_vs_p->Fill(part_p[i], part_vz[i]);

    if(part_Cal_PCAL_sector[i] == i+1  && part_charge[i] == +1 && part_vz[i] != 0) hist_positive_vertex_sec[i]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i] == i+1  && part_charge[i] == -1 && part_vz[i] != 0) hist_negative_vertex_sec[i]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i] == i+1  && part_charge[i] == +1 && part_vz[i] != 0 && part_theta[i] > 0) hist_positive_vertex_vs_theta_sec[i]->Fill(part_theta[i]*180/Pival, part_vz[i]);
    if(part_Cal_PCAL_sector[i] == i+1  && part_charge[i] == -1 && part_vz[i] != 0 && part_theta[i] > 0) hist_negative_vertex_vs_theta_sec[i]->Fill(part_theta[i]*180/Pival, part_vz[i]);
    if(part_Cal_PCAL_sector[i] == i+1  && part_charge[i] == +1 && part_vz[i] != 0 && part_phi[i] != 0) hist_positive_vertex_vs_phi_sec[i]->Fill(part_phi[i]*180/Pival, part_vz[i]);
    if(part_Cal_PCAL_sector[i] == i+1  && part_charge[i] == -1 && part_vz[i] != 0 && part_phi[i] != 0) hist_negative_vertex_vs_phi_sec[i]->Fill(part_phi[i]*180/Pival, part_vz[i]);
    if(part_Cal_PCAL_sector[i] == i+1  && part_charge[i] == +1 && part_vz[i] != 0 && part_p[i] > 0) hist_positive_vertex_vs_p_sec[i]->Fill(part_p[i], part_vz[i]);
    if(part_Cal_PCAL_sector[i] == i+1  && part_charge[i] == -1 && part_vz[i] != 0 && part_p[i] > 0) hist_negative_vertex_vs_p_sec[i]->Fill(part_p[i], part_vz[i]);

  }


  for(Int_t i = 0; i < BUFFER; i++){

    if(e_vz[i] != 0) hist_electron_vertex->Fill(e_vz[i]); 
    if(p_vz[i] != 0) hist_proton_vertex->Fill(p_vz[i]);
    if(e_vz[i] != 0 && p4_ele[i].Theta() > 0) hist_electron_vertex_vs_theta->Fill(p4_ele[i].Theta()*180/Pival, e_vz[i]);
    if(p_vz[i] != 0 && p4_prot[i].Theta() > 0) hist_proton_vertex_vs_theta->Fill(p4_prot[i].Theta()*180/Pival, e_vz[i]);
    if(e_vz[i] != 0 && p4_ele[i].Phi() != 0) hist_electron_vertex_vs_phi->Fill(p4_ele[i].Phi()*180/Pival, e_vz[i]);
    if(p_vz[i] != 0 && p4_prot[i].Phi() != 0) hist_proton_vertex_vs_phi->Fill(p4_prot[i].Phi()*180/Pival, e_vz[i]);
    if(e_vz[i] != 0 && p4_ele[i].P() > 0) hist_electron_vertex_vs_p->Fill(p4_ele[i].P(), e_vz[i]);
    if(p_vz[i] != 0 && p4_prot[i].P() > 0) hist_proton_vertex_vs_p->Fill(p4_prot[i].P(), e_vz[i]);

    if(e_PCAL_sec[i] == i+1  && e_vz[i] != 0) hist_electron_vertex_sec[i]->Fill(e_vz[i]); 
    if(p_PCAL_sec[i] == i+1  && p_vz[i] != 0) hist_proton_vertex_sec[i]->Fill(p_vz[i]);
    if(e_PCAL_sec[i] == i+1  && e_vz[i] != 0 && p4_ele[i].Theta() > 0) hist_electron_vertex_vs_theta_sec[i]->Fill(p4_ele[i].Theta()*180/Pival, e_vz[i]);
    if(p_PCAL_sec[i] == i+1  && p_vz[i] != 0 != 0 && p4_prot[i].Theta() > 0) hist_proton_vertex_vs_theta_sec[i]->Fill(p4_prot[i].Theta()*180/Pival, e_vz[i]);
    if(e_PCAL_sec[i] == i+1  && e_vz[i] != 0 && p4_ele[i].Phi() != 0) hist_electron_vertex_vs_phi_sec[i]->Fill(p4_ele[i].Phi()*180/Pival, e_vz[i]);
    if(p_PCAL_sec[i] == i+1  && p_vz[i] != 0 && p4_prot[i].Phi() != 0) hist_proton_vertex_vs_phi_sec[i]->Fill(p4_prot[i].Phi()*180/Pival, e_vz[i]);
    if(e_PCAL_sec[i] == i+1  && e_vz[i] != 0 && p4_ele[i].P() > 0) hist_electron_vertex_vs_p_sec[i]->Fill(p4_ele[i].P(), e_vz[i]);
    if(p_PCAL_sec[i] == i+1  && p_vz[i] != 0 && p4_prot[i].P() > 0) hist_proton_vertex_vs_p_sec[i]->Fill(p4_prot[i].P(), e_vz[i]);

  }


  // fill missing mass plots:

  for(Int_t i = 0; i < 6; i++){
    if(e_FTOF_sec[0] == i+1 && f_e_p_miss.M() != 0 && p_count > 0){
      hist_e_p_mismass_sec[i]->Fill(f_e_p_miss.M());
      hist_e_p_mismass->Fill(f_e_p_miss.M());  
    }
    if(e_FTOF_sec[0] == i+1 && f_e_p_miss.M() != 0 && p_count > 0){
      hist_e_p_mismass2_sec[i]->Fill(f_e_p_miss.M2());
      hist_e_p_mismass2->Fill(f_e_p_miss.M2());
    }
    if(e_FTOF_sec[0] == i+1 && f_e_pi0_miss.M() != 0){ 
      hist_e_pi0_mismass_sec[i]->Fill(f_e_pi0_miss.M());
      hist_e_pi0_mismass->Fill(f_e_pi0_miss.M());
    }
    if(e_FTOF_sec[0] == i+1 && f_e_pip_miss.M() != 0 && pip_count > 0){
      hist_e_pip_mismass_sec[i]->Fill(f_e_pip_miss.M());
      hist_e_pip_mismass->Fill(f_e_pip_miss.M());
    }
    if(e_FTOF_sec[0] == i+1 && f_e_pim_miss.M() != 0 && pim_count > 0){
      hist_e_pim_mismass_sec[i]->Fill(f_e_pim_miss.M());
      hist_e_pim_mismass->Fill(f_e_pim_miss.M());
    }
  }

  // fill Dalitz plots and two pion mass systems:

  if(p_count > 0 && pip_count > 0 && pim_count > 0){

    if(f_ppip.M() != 0 && f_pippim.M() != 0) hist_Dalitz_ppip_pippim->Fill(f_ppip.M(),f_pippim.M());
    if(f_ppim.M() != 0 && f_pippim.M() != 0) hist_Dalitz_ppim_pippim->Fill(f_ppim.M(),f_pippim.M());
    if(f_ppip.M() != 0 && f_ppim.M() != 0)   hist_Dalitz_ppip_ppim->Fill(f_ppip.M(),f_ppim.M());

    if(f_ppip.M2() != 0 && f_pippim.M2() != 0) hist_Dalitz2_ppip_pippim->Fill(f_ppip.M2(),f_pippim.M2());
    if(f_ppim.M2() != 0 && f_pippim.M2() != 0) hist_Dalitz2_ppim_pippim->Fill(f_ppim.M2(),f_pippim.M2());
    if(f_ppip.M2() != 0 && f_ppim.M2() != 0)   hist_Dalitz2_ppip_ppim->Fill(f_ppip.M2(),f_ppim.M2());

  }

  if(p_count > 0 && pip_count > 0 && f_ppip.M() != 0){hist_mass_ppip->Fill(f_ppip.M());}
  if(p_count > 0 && pim_count > 0 && f_ppim.M() != 0){hist_mass_ppim->Fill(f_ppim.M());}
  if(pip_count > 0 && pim_count > 0 && f_pippim.M() != 0){hist_mass_pippim->Fill(f_pippim.M());}


  // Phi angle versus sector

  for(Int_t i = 0; i < BUFFER; i++){
    if(part_charge[i] == +1 && part_phi[i] != 0) hist_FTOF_phi_vs_sector_positives->Fill(part_FTOF_sector_layer2[i], part_phi[i]*180/Pival);
    if(part_charge[i] == -1 && part_phi[i] != 0) hist_FTOF_phi_vs_sector_negatives->Fill(part_FTOF_sector_layer2[i], part_phi[i]*180/Pival);
  }


}



/// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///  fill the PID historgrams:
/// ///////////////////////////////////////////////////////////////

for(Int_t i = 0; i < BUFFER; i++){ 

if(fill_electron_pid_histograms){

// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// all particles with negative charge

  if(FD_eid_charge_check[i]){   // only particles which passed the default pid cut

    if(part_p[i] > 0) hist_HTCC_theta_vs_phi[0]->Fill(part_CC_HTCC_phi[i]*180/Pival, part_CC_HTCC_theta[i]*180/Pival);
    if(part_CC_HTCC_nphe[i] > 0) hist_HTCC_nphe[0]->Fill(part_CC_HTCC_nphe[i]);
    if(part_CC_HTCC_nphe[i] > 0 && part_p[i] > 0 && part_Cal_energy_total[i]/part_p[i] > 0) hist_HTCC_nphe_vs_sampling_fraction[0]->Fill(part_CC_HTCC_nphe[i], part_Cal_energy_total[i]/part_p[i]);

    // EC
    if((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0) hist_EC_PCAL_vs_EC_ECAL[0]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_Cal_ECin_energy[i] > 0 && part_Cal_ECout_energy[i] > 0) hist_EC_outer_vs_EC_inner[0]->Fill(part_Cal_ECin_energy[i], part_Cal_ECout_energy[i]);

    if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec1[0]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec2[0]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec3[0]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec4[0]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec5[0]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec6[0]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);

    if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec1[0]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec2[0]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec3[0]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec4[0]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec5[0]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec6[0]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);

    if(part_Cal_ECin_sector[i] == 1 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec1[0]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 2 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec2[0]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 3 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec3[0]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 4 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec4[0]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 5 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec5[0]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 6 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec6[0]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);

    if(part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0) hist_EC_PCAL_hit_position[0]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
    if(part_Cal_ECin_x[i] != 0 && part_Cal_ECin_y[i] != 0) hist_EC_inner_hit_position[0]->Fill(part_Cal_ECin_x[i], part_Cal_ECin_y[i]);
    if(part_Cal_ECout_x[i] != 0 && part_Cal_ECout_y[i] != 0) hist_EC_outer_hit_position[0]->Fill(part_Cal_ECout_x[i], part_Cal_ECout_y[i]);

    // DC
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1[0]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2[0]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3[0]->Fill(part_DC_c3x[i], part_DC_c3y[i]);

    if(part_Cal_PCAL_sector[i]== 1 && part_vz[i] != 0) hist_DC_z_vertex_sec1[0]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 2 && part_vz[i] != 0) hist_DC_z_vertex_sec2[0]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 3 && part_vz[i] != 0) hist_DC_z_vertex_sec3[0]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 4 && part_vz[i] != 0) hist_DC_z_vertex_sec4[0]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 5 && part_vz[i] != 0) hist_DC_z_vertex_sec5[0]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 6 && part_vz[i] != 0) hist_DC_z_vertex_sec6[0]->Fill(part_vz[i]);


  }


// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// default PID  (neg. charge is fulfilled for all)

  if(FD_eid_default_PID_check[i]){  

    if(part_p[i] > 0) hist_HTCC_theta_vs_phi[1]->Fill(part_CC_HTCC_phi[i]*180/Pival, part_CC_HTCC_theta[i]*180/Pival);
    if(part_CC_HTCC_nphe[i] > 0) hist_HTCC_nphe[1]->Fill(part_CC_HTCC_nphe[i]);
    if(part_CC_HTCC_nphe[i] > 0 && part_p[i] > 0 && part_Cal_energy_total[i]/part_p[i] > 0) hist_HTCC_nphe_vs_sampling_fraction[1]->Fill(part_CC_HTCC_nphe[i], part_Cal_energy_total[i]/part_p[i]);

    // EC
    if((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0) hist_EC_PCAL_vs_EC_ECAL[1]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_Cal_ECin_energy[i] > 0 && part_Cal_ECout_energy[i] > 0) hist_EC_outer_vs_EC_inner[1]->Fill(part_Cal_ECin_energy[i], part_Cal_ECout_energy[i]);

    if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec1[1]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec2[1]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec3[1]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec4[1]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec5[1]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec6[1]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);

    if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec1[1]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec2[1]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec3[1]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec4[1]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec5[1]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec6[1]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);

    if(part_Cal_ECin_sector[i] == 1 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec1[1]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 2 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec2[1]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 3 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec3[1]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 4 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec4[1]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 5 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec5[1]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 6 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec6[1]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);

    if(part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0) hist_EC_PCAL_hit_position[1]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
    if(part_Cal_ECin_x[i] != 0 && part_Cal_ECin_y[i] != 0) hist_EC_inner_hit_position[1]->Fill(part_Cal_ECin_x[i], part_Cal_ECin_y[i]);
    if(part_Cal_ECout_x[i] != 0 && part_Cal_ECout_y[i] != 0) hist_EC_outer_hit_position[1]->Fill(part_Cal_ECout_x[i], part_Cal_ECout_y[i]);

    // DC
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1[1]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2[1]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3[1]->Fill(part_DC_c3x[i], part_DC_c3y[i]);

    if(part_Cal_PCAL_sector[i]== 1 && part_vz[i] != 0) hist_DC_z_vertex_sec1[1]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 2 && part_vz[i] != 0) hist_DC_z_vertex_sec2[1]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 3 && part_vz[i] != 0) hist_DC_z_vertex_sec3[1]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 4 && part_vz[i] != 0) hist_DC_z_vertex_sec4[1]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 5 && part_vz[i] != 0) hist_DC_z_vertex_sec5[1]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 6 && part_vz[i] != 0) hist_DC_z_vertex_sec6[1]->Fill(part_vz[i]);
  }


// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// EC cuts

  if(FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_EC_outer_vs_EC_inner_check[i]){  

    if(part_p[i] > 0) hist_HTCC_theta_vs_phi[2]->Fill(part_CC_HTCC_phi[i]*180/Pival, part_CC_HTCC_theta[i]*180/Pival);
    if(part_CC_HTCC_nphe[i] > 0) hist_HTCC_nphe[2]->Fill(part_CC_HTCC_nphe[i]);
    if(part_CC_HTCC_nphe[i] > 0 && part_p[i] > 0 && part_Cal_energy_total[i]/part_p[i] > 0) hist_HTCC_nphe_vs_sampling_fraction[2]->Fill(part_CC_HTCC_nphe[i], part_Cal_energy_total[i]/part_p[i]);

    // EC
    if((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0) hist_EC_PCAL_vs_EC_ECAL[2]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_Cal_ECin_energy[i] > 0 && part_Cal_ECout_energy[i] > 0) hist_EC_outer_vs_EC_inner[2]->Fill(part_Cal_ECin_energy[i], part_Cal_ECout_energy[i]);

    if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec1[2]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec2[2]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec3[2]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec4[2]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec5[2]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec6[2]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);

    if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec1[2]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec2[2]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec3[2]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec4[2]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec5[2]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec6[2]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);

    if(part_Cal_ECin_sector[i] == 1 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec1[2]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 2 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec2[2]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 3 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec3[2]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 4 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec4[2]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 5 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec5[2]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 6 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec6[2]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);

    if(part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0) hist_EC_PCAL_hit_position[2]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
    if(part_Cal_ECin_x[i] != 0 && part_Cal_ECin_y[i] != 0) hist_EC_inner_hit_position[2]->Fill(part_Cal_ECin_x[i], part_Cal_ECin_y[i]);
    if(part_Cal_ECout_x[i] != 0 && part_Cal_ECout_y[i] != 0) hist_EC_outer_hit_position[2]->Fill(part_Cal_ECout_x[i], part_Cal_ECout_y[i]);

    // DC
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1[2]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2[2]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3[2]->Fill(part_DC_c3x[i], part_DC_c3y[i]);

    if(part_Cal_PCAL_sector[i]== 1 && part_vz[i] != 0) hist_DC_z_vertex_sec1[2]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 2 && part_vz[i] != 0) hist_DC_z_vertex_sec2[2]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 3 && part_vz[i] != 0) hist_DC_z_vertex_sec3[2]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 4 && part_vz[i] != 0) hist_DC_z_vertex_sec4[2]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 5 && part_vz[i] != 0) hist_DC_z_vertex_sec5[2]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 6 && part_vz[i] != 0) hist_DC_z_vertex_sec6[2]->Fill(part_vz[i]);
  }


  //if(FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_EC_outer_vs_EC_inner_check[i] && FD_eid_DC_z_vertex_check[i] && FD_eid_DC_hit_position_region1_fiducial_check[i] && FD_eid_DC_hit_position_region2_fiducial_check[i] && FD_eid_DC_hit_position_region3_fiducial_check[i] && FD_eid_EC_hit_position_fiducial_check[i] && FD_eid_EC_sampling_fraction_check[i]){  

  if(FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_EC_sampling_fraction_check[i]){  

    if(part_p[i] > 0) hist_HTCC_theta_vs_phi[3]->Fill(part_CC_HTCC_phi[i]*180/Pival, part_CC_HTCC_theta[i]*180/Pival);
    if(part_CC_HTCC_nphe[i] > 0) hist_HTCC_nphe[3]->Fill(part_CC_HTCC_nphe[i]);
    if(part_CC_HTCC_nphe[i] > 0 && part_p[i] > 0 && part_Cal_energy_total[i]/part_p[i] > 0) hist_HTCC_nphe_vs_sampling_fraction[3]->Fill(part_CC_HTCC_nphe[i], part_Cal_energy_total[i]/part_p[i]);

    // EC
    if((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0) hist_EC_PCAL_vs_EC_ECAL[3]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_Cal_ECin_energy[i] > 0 && part_Cal_ECout_energy[i] > 0) hist_EC_outer_vs_EC_inner[3]->Fill(part_Cal_ECin_energy[i], part_Cal_ECout_energy[i]);

    if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec1[3]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec2[3]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec3[3]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec4[3]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec5[3]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec6[3]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);

    if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec1[3]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec2[3]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec3[3]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec4[3]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec5[3]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec6[3]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);

    if(part_Cal_ECin_sector[i] == 1 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec1[3]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 2 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec2[3]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 3 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec3[3]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 4 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec4[3]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 5 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec5[3]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 6 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec6[3]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);

    if(part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0) hist_EC_PCAL_hit_position[3]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
    if(part_Cal_ECin_x[i] != 0 && part_Cal_ECin_y[i] != 0) hist_EC_inner_hit_position[3]->Fill(part_Cal_ECin_x[i], part_Cal_ECin_y[i]);
    if(part_Cal_ECout_x[i] != 0 && part_Cal_ECout_y[i] != 0) hist_EC_outer_hit_position[3]->Fill(part_Cal_ECout_x[i], part_Cal_ECout_y[i]);

    // DC
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1[3]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2[3]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3[3]->Fill(part_DC_c3x[i], part_DC_c3y[i]);

    if(part_Cal_PCAL_sector[i]== 1 && part_vz[i] != 0) hist_DC_z_vertex_sec1[3]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 2 && part_vz[i] != 0) hist_DC_z_vertex_sec2[3]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 3 && part_vz[i] != 0) hist_DC_z_vertex_sec3[3]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 4 && part_vz[i] != 0) hist_DC_z_vertex_sec4[3]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 5 && part_vz[i] != 0) hist_DC_z_vertex_sec5[3]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 6 && part_vz[i] != 0) hist_DC_z_vertex_sec6[3]->Fill(part_vz[i]);
  }

  if(FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_EC_hit_position_fiducial_check[i]){   

    if(part_p[i] > 0) hist_HTCC_theta_vs_phi[4]->Fill(part_CC_HTCC_phi[i]*180/Pival, part_CC_HTCC_theta[i]*180/Pival);
    if(part_CC_HTCC_nphe[i] > 0) hist_HTCC_nphe[4]->Fill(part_CC_HTCC_nphe[i]);
    if(part_CC_HTCC_nphe[i] > 0 && part_p[i] > 0 && part_Cal_energy_total[i]/part_p[i] > 0) hist_HTCC_nphe_vs_sampling_fraction[4]->Fill(part_CC_HTCC_nphe[i], part_Cal_energy_total[i]/part_p[i]);

    // EC
    if((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0) hist_EC_PCAL_vs_EC_ECAL[4]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_Cal_ECin_energy[i] > 0 && part_Cal_ECout_energy[i] > 0) hist_EC_outer_vs_EC_inner[4]->Fill(part_Cal_ECin_energy[i], part_Cal_ECout_energy[i]);

    if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec1[4]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec2[4]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec3[4]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec4[4]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec5[4]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec6[4]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);

    if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec1[4]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec2[4]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec3[4]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec4[4]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec5[4]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec6[4]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);

    if(part_Cal_ECin_sector[i] == 1 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec1[4]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 2 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec2[4]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 3 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec3[4]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 4 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec4[4]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 5 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec5[4]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 6 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec6[4]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);

    if(part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0) hist_EC_PCAL_hit_position[4]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
    if(part_Cal_ECin_x[i] != 0 && part_Cal_ECin_y[i] != 0) hist_EC_inner_hit_position[4]->Fill(part_Cal_ECin_x[i], part_Cal_ECin_y[i]);
    if(part_Cal_ECout_x[i] != 0 && part_Cal_ECout_y[i] != 0) hist_EC_outer_hit_position[4]->Fill(part_Cal_ECout_x[i], part_Cal_ECout_y[i]);

    // DC
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1[4]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2[4]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3[4]->Fill(part_DC_c3x[i], part_DC_c3y[i]);

    if(part_Cal_PCAL_sector[i]== 1 && part_vz[i] != 0) hist_DC_z_vertex_sec1[4]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 2 && part_vz[i] != 0) hist_DC_z_vertex_sec2[4]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 3 && part_vz[i] != 0) hist_DC_z_vertex_sec3[4]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 4 && part_vz[i] != 0) hist_DC_z_vertex_sec4[4]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 5 && part_vz[i] != 0) hist_DC_z_vertex_sec5[4]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 6 && part_vz[i] != 0) hist_DC_z_vertex_sec6[4]->Fill(part_vz[i]);
  }


// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DC cuts

  if(FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_DC_hit_position_region1_fiducial_check[i]){   

    if(part_p[i] > 0) hist_HTCC_theta_vs_phi[5]->Fill(part_CC_HTCC_phi[i]*180/Pival, part_CC_HTCC_theta[i]*180/Pival);
    if(part_CC_HTCC_nphe[i] > 0) hist_HTCC_nphe[5]->Fill(part_CC_HTCC_nphe[i]);
    if(part_CC_HTCC_nphe[i] > 0 && part_p[i] > 0 && part_Cal_energy_total[i]/part_p[i] > 0) hist_HTCC_nphe_vs_sampling_fraction[5]->Fill(part_CC_HTCC_nphe[i], part_Cal_energy_total[i]/part_p[i]);

    // EC
    if((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0) hist_EC_PCAL_vs_EC_ECAL[5]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_Cal_ECin_energy[i] > 0 && part_Cal_ECout_energy[i] > 0) hist_EC_outer_vs_EC_inner[5]->Fill(part_Cal_ECin_energy[i], part_Cal_ECout_energy[i]);

    if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec1[5]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec2[5]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec3[5]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec4[5]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec5[5]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec6[5]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);

    if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec1[5]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec2[5]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec3[5]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec4[5]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec5[5]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec6[5]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);

    if(part_Cal_ECin_sector[i] == 1 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec1[5]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 2 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec2[5]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 3 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec3[5]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 4 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec4[5]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 5 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec5[5]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 6 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec6[5]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);

    if(part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0) hist_EC_PCAL_hit_position[5]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
    if(part_Cal_ECin_x[i] != 0 && part_Cal_ECin_y[i] != 0) hist_EC_inner_hit_position[5]->Fill(part_Cal_ECin_x[i], part_Cal_ECin_y[i]);
    if(part_Cal_ECout_x[i] != 0 && part_Cal_ECout_y[i] != 0) hist_EC_outer_hit_position[5]->Fill(part_Cal_ECout_x[i], part_Cal_ECout_y[i]);

    // DC
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1[5]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2[5]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3[5]->Fill(part_DC_c3x[i], part_DC_c3y[i]);

    if(part_Cal_PCAL_sector[i]== 1 && part_vz[i] != 0) hist_DC_z_vertex_sec1[5]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 2 && part_vz[i] != 0) hist_DC_z_vertex_sec2[5]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 3 && part_vz[i] != 0) hist_DC_z_vertex_sec3[5]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 4 && part_vz[i] != 0) hist_DC_z_vertex_sec4[5]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 5 && part_vz[i] != 0) hist_DC_z_vertex_sec5[5]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 6 && part_vz[i] != 0) hist_DC_z_vertex_sec6[5]->Fill(part_vz[i]);
  }

  if(FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_DC_hit_position_region2_fiducial_check[i]){   
    hist_DC_hit_position_region2_cut5a->Fill(part_DC_c2x[i], part_DC_c2y[i]);
  }

  if(FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_DC_hit_position_region3_fiducial_check[i]){   

    if(part_p[i] > 0) hist_HTCC_theta_vs_phi[6]->Fill(part_CC_HTCC_phi[i]*180/Pival, part_CC_HTCC_theta[i]*180/Pival);
    if(part_CC_HTCC_nphe[i] > 0) hist_HTCC_nphe[6]->Fill(part_CC_HTCC_nphe[i]);
    if(part_CC_HTCC_nphe[i] > 0 && part_p[i] > 0 && part_Cal_energy_total[i]/part_p[i] > 0) hist_HTCC_nphe_vs_sampling_fraction[6]->Fill(part_CC_HTCC_nphe[i], part_Cal_energy_total[i]/part_p[i]);

    // EC
    if((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0) hist_EC_PCAL_vs_EC_ECAL[6]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_Cal_ECin_energy[i] > 0 && part_Cal_ECout_energy[i] > 0) hist_EC_outer_vs_EC_inner[6]->Fill(part_Cal_ECin_energy[i], part_Cal_ECout_energy[i]);

    if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec1[6]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec2[6]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec3[6]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec4[6]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec5[6]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec6[6]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);

    if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec1[6]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec2[6]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec3[6]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec4[6]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec5[6]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec6[6]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);

    if(part_Cal_ECin_sector[i] == 1 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec1[6]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 2 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec2[6]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 3 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec3[6]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 4 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec4[6]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 5 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec5[6]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 6 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec6[6]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);

    if(part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0) hist_EC_PCAL_hit_position[6]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
    if(part_Cal_ECin_x[i] != 0 && part_Cal_ECin_y[i] != 0) hist_EC_inner_hit_position[6]->Fill(part_Cal_ECin_x[i], part_Cal_ECin_y[i]);
    if(part_Cal_ECout_x[i] != 0 && part_Cal_ECout_y[i] != 0) hist_EC_outer_hit_position[6]->Fill(part_Cal_ECout_x[i], part_Cal_ECout_y[i]);

    // DC
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1[6]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2[6]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3[6]->Fill(part_DC_c3x[i], part_DC_c3y[i]);

    if(part_Cal_PCAL_sector[i]== 1 && part_vz[i] != 0) hist_DC_z_vertex_sec1[6]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 2 && part_vz[i] != 0) hist_DC_z_vertex_sec2[6]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 3 && part_vz[i] != 0) hist_DC_z_vertex_sec3[6]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 4 && part_vz[i] != 0) hist_DC_z_vertex_sec4[6]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 5 && part_vz[i] != 0) hist_DC_z_vertex_sec5[6]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 6 && part_vz[i] != 0) hist_DC_z_vertex_sec6[6]->Fill(part_vz[i]);
  }

// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// vertex cut

  if(FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_DC_z_vertex_check[i]){  

    if(part_p[i] > 0) hist_HTCC_theta_vs_phi[7]->Fill(part_CC_HTCC_phi[i]*180/Pival, part_CC_HTCC_theta[i]*180/Pival);
    if(part_CC_HTCC_nphe[i] > 0) hist_HTCC_nphe[7]->Fill(part_CC_HTCC_nphe[i]);
    if(part_CC_HTCC_nphe[i] > 0 && part_p[i] > 0 && part_Cal_energy_total[i]/part_p[i] > 0) hist_HTCC_nphe_vs_sampling_fraction[7]->Fill(part_CC_HTCC_nphe[i], part_Cal_energy_total[i]/part_p[i]);

    // EC
    if((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0) hist_EC_PCAL_vs_EC_ECAL[7]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_Cal_ECin_energy[i] > 0 && part_Cal_ECout_energy[i] > 0) hist_EC_outer_vs_EC_inner[7]->Fill(part_Cal_ECin_energy[i], part_Cal_ECout_energy[i]);

    if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec1[7]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec2[7]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec3[7]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec4[7]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec5[7]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec6[7]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);

    if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec1[7]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec2[7]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec3[7]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec4[7]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec5[7]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec6[7]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);

    if(part_Cal_ECin_sector[i] == 1 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec1[7]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 2 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec2[7]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 3 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec3[7]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 4 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec4[7]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 5 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec5[7]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 6 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec6[7]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);

    if(part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0) hist_EC_PCAL_hit_position[7]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
    if(part_Cal_ECin_x[i] != 0 && part_Cal_ECin_y[i] != 0) hist_EC_inner_hit_position[7]->Fill(part_Cal_ECin_x[i], part_Cal_ECin_y[i]);
    if(part_Cal_ECout_x[i] != 0 && part_Cal_ECout_y[i] != 0) hist_EC_outer_hit_position[7]->Fill(part_Cal_ECout_x[i], part_Cal_ECout_y[i]);

    // DC
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1[7]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2[7]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3[7]->Fill(part_DC_c3x[i], part_DC_c3y[i]);

    if(part_Cal_PCAL_sector[i]== 1 && part_vz[i] != 0) hist_DC_z_vertex_sec1[7]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 2 && part_vz[i] != 0) hist_DC_z_vertex_sec2[7]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 3 && part_vz[i] != 0) hist_DC_z_vertex_sec3[7]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 4 && part_vz[i] != 0) hist_DC_z_vertex_sec4[7]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 5 && part_vz[i] != 0) hist_DC_z_vertex_sec5[7]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 6 && part_vz[i] != 0) hist_DC_z_vertex_sec6[7]->Fill(part_vz[i]);
  }


// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// all other

  // CC

  if(FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_EC_outer_vs_EC_inner_check[i] && FD_eid_EC_sampling_fraction_check[i] && FD_eid_EC_hit_position_fiducial_check[i]
     && FD_eid_DC_hit_position_region1_fiducial_check[i] && FD_eid_DC_hit_position_region3_fiducial_check[i] && FD_eid_DC_z_vertex_check[i]){
    if(part_p[i] > 0) hist_HTCC_theta_vs_phi[8]->Fill(part_CC_HTCC_phi[i]*180/Pival, part_CC_HTCC_theta[i]*180/Pival);
    if(part_CC_HTCC_nphe[i] > 0) hist_HTCC_nphe[8]->Fill(part_CC_HTCC_nphe[i]);
    if(part_CC_HTCC_nphe[i] > 0 && part_p[i] > 0 && part_Cal_energy_total[i]/part_p[i] > 0) hist_HTCC_nphe_vs_sampling_fraction[8]->Fill(part_CC_HTCC_nphe[i], part_Cal_energy_total[i]/part_p[i]);
  }

  // EC

  if(FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_EC_sampling_fraction_check[i] && FD_eid_EC_hit_position_fiducial_check[i]
     && FD_eid_DC_hit_position_region1_fiducial_check[i] && FD_eid_DC_hit_position_region3_fiducial_check[i] && FD_eid_DC_z_vertex_check[i]){

    if((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0) hist_EC_PCAL_vs_EC_ECAL[8]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_Cal_ECin_energy[i] > 0 && part_Cal_ECout_energy[i] > 0) hist_EC_outer_vs_EC_inner[8]->Fill(part_Cal_ECin_energy[i], part_Cal_ECout_energy[i]);
  }

  if(FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_EC_outer_vs_EC_inner_check[i] && FD_eid_EC_hit_position_fiducial_check[i]
     && FD_eid_DC_hit_position_region1_fiducial_check[i] && FD_eid_DC_hit_position_region3_fiducial_check[i] && FD_eid_DC_z_vertex_check[i]){
    if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec1[8]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec2[8]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec3[8]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec4[8]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec5[8]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec6[8]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);

    if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec1[8]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec2[8]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec3[8]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec4[8]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec5[8]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec6[8]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);

    if(part_Cal_ECin_sector[i] == 1 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec1[8]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 2 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec2[8]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 3 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec3[8]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 4 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec4[8]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 5 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec5[8]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 6 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec6[8]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
  }

  if(FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_EC_outer_vs_EC_inner_check[i] && FD_eid_EC_sampling_fraction_check[i]
     && FD_eid_DC_hit_position_region1_fiducial_check[i] && FD_eid_DC_hit_position_region3_fiducial_check[i] && FD_eid_DC_z_vertex_check[i]){
    if(part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0) hist_EC_PCAL_hit_position[8]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
    if(part_Cal_ECin_x[i] != 0 && part_Cal_ECin_y[i] != 0) hist_EC_inner_hit_position[8]->Fill(part_Cal_ECin_x[i], part_Cal_ECin_y[i]);
    if(part_Cal_ECout_x[i] != 0 && part_Cal_ECout_y[i] != 0) hist_EC_outer_hit_position[8]->Fill(part_Cal_ECout_x[i], part_Cal_ECout_y[i]);
  }

    // DC

  if(FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_EC_outer_vs_EC_inner_check[i] && FD_eid_EC_sampling_fraction_check[i] 
     && FD_eid_EC_hit_position_fiducial_check[i] && FD_eid_DC_hit_position_region3_fiducial_check[i] && FD_eid_DC_hit_position_region2_fiducial_check[i] && FD_eid_DC_z_vertex_check[i]){

    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1[8]->Fill(part_DC_c1x[i], part_DC_c1y[i]);

  }

  if(FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_EC_outer_vs_EC_inner_check[i] && FD_eid_EC_sampling_fraction_check[i] 
     && FD_eid_EC_hit_position_fiducial_check[i] && FD_eid_DC_hit_position_region1_fiducial_check[i] && FD_eid_DC_hit_position_region3_fiducial_check[i] && FD_eid_DC_z_vertex_check[i]){

    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2[0]->Fill(part_DC_c2x[i], part_DC_c2y[i]);

  }

  if(FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_EC_outer_vs_EC_inner_check[i] && FD_eid_EC_sampling_fraction_check[i] 
     && FD_eid_EC_hit_position_fiducial_check[i] && FD_eid_DC_hit_position_region1_fiducial_check[i] && FD_eid_DC_hit_position_region2_fiducial_check[i] && FD_eid_DC_z_vertex_check[i]){

    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3[8]->Fill(part_DC_c3x[i], part_DC_c3y[i]);

  }

  if(FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_EC_outer_vs_EC_inner_check[i] && FD_eid_EC_sampling_fraction_check[i] 
     && FD_eid_EC_hit_position_fiducial_check[i]&& FD_eid_DC_hit_position_region1_fiducial_check[i] && FD_eid_DC_hit_position_region2_fiducial_check[i] && FD_eid_DC_hit_position_region3_fiducial_check[i]){
    if(part_Cal_PCAL_sector[i]== 1 && part_vz[i] != 0) hist_DC_z_vertex_sec1[8]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 2 && part_vz[i] != 0) hist_DC_z_vertex_sec2[8]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 3 && part_vz[i] != 0) hist_DC_z_vertex_sec3[8]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 4 && part_vz[i] != 0) hist_DC_z_vertex_sec4[8]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 5 && part_vz[i] != 0) hist_DC_z_vertex_sec5[8]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 6 && part_vz[i] != 0) hist_DC_z_vertex_sec6[8]->Fill(part_vz[i]);
  }


// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// all cuts

  if(FD_eid_all_check[i]){ 

    if(part_p[i] > 0) hist_HTCC_theta_vs_phi[9]->Fill(part_CC_HTCC_phi[i]*180/Pival, part_CC_HTCC_theta[i]*180/Pival);
    if(part_CC_HTCC_nphe[i] > 0) hist_HTCC_nphe[9]->Fill(part_CC_HTCC_nphe[i]);
    if(part_CC_HTCC_nphe[i] > 0 && part_p[i] > 0 && part_Cal_energy_total[i]/part_p[i] > 0) hist_HTCC_nphe_vs_sampling_fraction[9]->Fill(part_CC_HTCC_nphe[i], part_Cal_energy_total[i]/part_p[i]);

    // EC
    if((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0) hist_EC_PCAL_vs_EC_ECAL[9]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_Cal_ECin_energy[i] > 0 && part_Cal_ECout_energy[i] > 0) hist_EC_outer_vs_EC_inner[9]->Fill(part_Cal_ECin_energy[i], part_Cal_ECout_energy[i]);

    if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec1[9]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec2[9]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec3[9]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec4[9]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec5[9]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec6[9]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);

    if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec1[9]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec2[9]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec3[9]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec4[9]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec5[9]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec6[9]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);

    if(part_Cal_ECin_sector[i] == 1 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec1[9]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 2 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec2[9]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 3 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec3[9]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 4 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec4[9]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 5 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec5[9]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 6 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec6[9]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);

    if(part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0) hist_EC_PCAL_hit_position[9]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
    if(part_Cal_ECin_x[i] != 0 && part_Cal_ECin_y[i] != 0) hist_EC_inner_hit_position[9]->Fill(part_Cal_ECin_x[i], part_Cal_ECin_y[i]);
    if(part_Cal_ECout_x[i] != 0 && part_Cal_ECout_y[i] != 0) hist_EC_outer_hit_position[9]->Fill(part_Cal_ECout_x[i], part_Cal_ECout_y[i]);

    // DC
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1[9]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2[9]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3[9]->Fill(part_DC_c3x[i], part_DC_c3y[i]);

    if(part_Cal_PCAL_sector[i]== 1 && part_vz[i] != 0) hist_DC_z_vertex_sec1[9]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 2 && part_vz[i] != 0) hist_DC_z_vertex_sec2[9]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 3 && part_vz[i] != 0) hist_DC_z_vertex_sec3[9]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 4 && part_vz[i] != 0) hist_DC_z_vertex_sec4[9]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 5 && part_vz[i] != 0) hist_DC_z_vertex_sec5[9]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 6 && part_vz[i] != 0) hist_DC_z_vertex_sec6[9]->Fill(part_vz[i]);

    if(part_Cal_PCAL_sector[i]== 1 && p4_ele[i].Theta()*180/Pival > 0) hist_ele_theta_sec1->Fill(p4_ele[i].Theta()*180/Pival);
    if(part_Cal_PCAL_sector[i]== 2 && p4_ele[i].Theta()*180/Pival > 0) hist_ele_theta_sec2->Fill(p4_ele[i].Theta()*180/Pival);
    if(part_Cal_PCAL_sector[i]== 3 && p4_ele[i].Theta()*180/Pival > 0) hist_ele_theta_sec3->Fill(p4_ele[i].Theta()*180/Pival);
    if(part_Cal_PCAL_sector[i]== 4 && p4_ele[i].Theta()*180/Pival > 0) hist_ele_theta_sec4->Fill(p4_ele[i].Theta()*180/Pival);
    if(part_Cal_PCAL_sector[i]== 5 && p4_ele[i].Theta()*180/Pival > 0) hist_ele_theta_sec5->Fill(p4_ele[i].Theta()*180/Pival);
    if(part_Cal_PCAL_sector[i]== 6 && p4_ele[i].Theta()*180/Pival > 0) hist_ele_theta_sec6->Fill(p4_ele[i].Theta()*180/Pival);

    if(part_Cal_PCAL_sector[i]== 1 && W > 0) hist_ele_W_sec1->Fill(W);
    if(part_Cal_PCAL_sector[i]== 2 && W > 0) hist_ele_W_sec2->Fill(W);
    if(part_Cal_PCAL_sector[i]== 3 && W > 0) hist_ele_W_sec3->Fill(W);
    if(part_Cal_PCAL_sector[i]== 4 && W > 0) hist_ele_W_sec4->Fill(W);
    if(part_Cal_PCAL_sector[i]== 5 && W > 0) hist_ele_W_sec5->Fill(W);
    if(part_Cal_PCAL_sector[i]== 6 && W > 0) hist_ele_W_sec6->Fill(W);

    if(part_Cal_PCAL_sector[i]== 1 && Q2 > 0) hist_ele_Q2_sec1->Fill(Q2);
    if(part_Cal_PCAL_sector[i]== 2 && Q2 > 0) hist_ele_Q2_sec2->Fill(Q2);
    if(part_Cal_PCAL_sector[i]== 3 && Q2 > 0) hist_ele_Q2_sec3->Fill(Q2);
    if(part_Cal_PCAL_sector[i]== 4 && Q2 > 0) hist_ele_Q2_sec4->Fill(Q2);
    if(part_Cal_PCAL_sector[i]== 5 && Q2 > 0) hist_ele_Q2_sec5->Fill(Q2);
    if(part_Cal_PCAL_sector[i]== 6 && Q2 > 0) hist_ele_Q2_sec6->Fill(Q2);

    if(part_Cal_PCAL_sector[i]== 1 && W > 0 && Q2 > 0) hist_ele_W_vs_Q2_sec1->Fill(W, Q2);
    if(part_Cal_PCAL_sector[i]== 2 && W > 0 && Q2 > 0) hist_ele_W_vs_Q2_sec2->Fill(W, Q2);
    if(part_Cal_PCAL_sector[i]== 3 && W > 0 && Q2 > 0) hist_ele_W_vs_Q2_sec3->Fill(W, Q2);
    if(part_Cal_PCAL_sector[i]== 4 && W > 0 && Q2 > 0) hist_ele_W_vs_Q2_sec4->Fill(W, Q2);
    if(part_Cal_PCAL_sector[i]== 5 && W > 0 && Q2 > 0) hist_ele_W_vs_Q2_sec5->Fill(W, Q2);
    if(part_Cal_PCAL_sector[i]== 6 && W > 0 && Q2 > 0) hist_ele_W_vs_Q2_sec6->Fill(W, Q2);

  }


// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// inverse cut  (neg. charge but no PID 11)

  if(FD_eid_default_PID_check[i] == false && FD_eid_charge_check[i] == true ){  

    if(part_p[i] > 0) hist_HTCC_theta_vs_phi[10]->Fill(part_CC_HTCC_phi[i]*180/Pival, part_CC_HTCC_theta[i]*180/Pival);
    if(part_CC_HTCC_nphe[i] > 0) hist_HTCC_nphe[10]->Fill(part_CC_HTCC_nphe[i]);
    if(part_CC_HTCC_nphe[i] > 0 && part_p[i] > 0 && part_Cal_energy_total[i]/part_p[i] > 0) hist_HTCC_nphe_vs_sampling_fraction[10]->Fill(part_CC_HTCC_nphe[i], part_Cal_energy_total[i]/part_p[i]);

    // EC
    if((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0) hist_EC_PCAL_vs_EC_ECAL[10]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_Cal_ECin_energy[i] > 0 && part_Cal_ECout_energy[i] > 0) hist_EC_outer_vs_EC_inner[10]->Fill(part_Cal_ECin_energy[i], part_Cal_ECout_energy[i]);

    if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec1[10]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec2[10]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec3[10]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec4[10]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec5[10]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_EC_total_sampling_fraction_sec6[10]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);

    if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec1[10]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec2[10]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec3[10]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec4[10]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec5[10]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_EC_PCAL_sampling_fraction_sec6[10]->Fill(part_p[i], part_Cal_PCAL_energy[i]/part_p[i]);

    if(part_Cal_ECin_sector[i] == 1 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec1[10]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 2 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec2[10]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 3 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec3[10]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 4 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec4[10]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 5 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec5[10]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);
    if(part_Cal_ECin_sector[i] == 6 && part_p[i] > 0) hist_EC_ECAL_sampling_fraction_sec6[10]->Fill(part_p[i], (part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])/part_p[i]);

    if(part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0) hist_EC_PCAL_hit_position[10]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
    if(part_Cal_ECin_x[i] != 0 && part_Cal_ECin_y[i] != 0) hist_EC_inner_hit_position[10]->Fill(part_Cal_ECin_x[i], part_Cal_ECin_y[i]);
    if(part_Cal_ECout_x[i] != 0 && part_Cal_ECout_y[i] != 0) hist_EC_outer_hit_position[10]->Fill(part_Cal_ECout_x[i], part_Cal_ECout_y[i]);

    // DC
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1[10]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2[10]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3[10]->Fill(part_DC_c3x[i], part_DC_c3y[i]);

    if(part_Cal_PCAL_sector[i]== 1 && part_vz[i] != 0) hist_DC_z_vertex_sec1[10]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 2 && part_vz[i] != 0) hist_DC_z_vertex_sec2[10]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 3 && part_vz[i] != 0) hist_DC_z_vertex_sec3[10]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 4 && part_vz[i] != 0) hist_DC_z_vertex_sec4[10]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 5 && part_vz[i] != 0) hist_DC_z_vertex_sec5[10]->Fill(part_vz[i]);
    if(part_Cal_PCAL_sector[i]== 6 && part_vz[i] != 0) hist_DC_z_vertex_sec6[10]->Fill(part_vz[i]);


  }

} // end of fill electron histogram selection

/// //////////////////////////////////////////////////////////////////////////////////////////////////////
/// TOF charged hadron plots
/// /////////////////////////////////////////////////////////////////

double beta_charge = Beta_charged(i, run);
double beta_neutr = Beta_neutral(i, run);

if(fill_hadron_pid_histograms){

  // proton

  if(FD_protid_default_PID_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[0]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[0]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[0]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[0]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[0]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[0]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[0]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[0]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[0]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[0]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[0]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[0]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_p*m_p) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[0]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_p*m_p) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[0]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[0]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[0]->Fill(Getdvz(i));
  }
  if(FD_protid_charge_check[i] && FD_protid_DC_hit_position_region1_fiducial_check[i] && FD_protid_DC_hit_position_region2_fiducial_check[i] && FD_protid_DC_hit_position_region3_fiducial_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[1]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[1]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[1]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[1]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] > 0.4 && beta_charge > 0) hist_beta_vs_p[1]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[1]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[1]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[1]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[1]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[1]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[1]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[1]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_p*m_p) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[1]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_p*m_p) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[1]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[1]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[1]->Fill(Getdvz(i));
  }
  if(FD_protid_default_PID_check[i] && FD_protid_charge_check[i] && FD_protid_DC_hit_position_region1_fiducial_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[2]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[2]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[2]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[2]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[2]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[2]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[2]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[2]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[2]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[2]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[2]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[2]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_p*m_p) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[2]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_p*m_p) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[2]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[2]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[2]->Fill(Getdvz(i));
  }
  if(FD_protid_default_PID_check[i] && FD_protid_charge_check[i] && FD_protid_DC_hit_position_region2_fiducial_check[i]){
    hist_DC_hit_position_region2_hadron_cut_02a->Fill(part_DC_c2x[i], part_DC_c2y[i]);
  }
  if(FD_protid_default_PID_check[i] && FD_protid_charge_check[i] && FD_protid_DC_hit_position_region3_fiducial_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[3]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[3]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[3]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[3]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[3]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[3]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[3]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[3]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[3]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[3]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[3]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[3]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_p*m_p) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[3]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_p*m_p) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[3]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[3]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[3]->Fill(Getdvz(i));
  }
  if(FD_protid_default_PID_check[i] && FD_protid_charge_check[i] && FD_protid_DC_hit_position_region1_fiducial_check[i] && FD_protid_DC_hit_position_region2_fiducial_check[i] 
                                    && FD_protid_DC_hit_position_region3_fiducial_check[i] && FD_protid_beta_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[4]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[4]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[4]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[4]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[4]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[4]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[4]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[4]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[4]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[4]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[4]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[4]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_p*m_p) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[4]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_p*m_p) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[4]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[4]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[4]->Fill(Getdvz(i));
  }
  if(FD_protid_default_PID_check[i] && FD_protid_charge_check[i] && FD_protid_DC_hit_position_region1_fiducial_check[i] && FD_protid_DC_hit_position_region2_fiducial_check[i]
                                    && FD_protid_DC_hit_position_region3_fiducial_check[i] && FD_protid_delta_beta_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[5]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[5]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[5]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[5]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[5]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[5]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[5]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[5]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[5]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[5]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[5]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[5]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_p*m_p) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[5]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_p*m_p) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[5]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[5]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[5]->Fill(Getdvz(i));
  }
  if(FD_protid_default_PID_check[i] && FD_protid_charge_check[i] && FD_protid_DC_hit_position_region1_fiducial_check[i] && FD_protid_DC_hit_position_region2_fiducial_check[i] 
                                    && FD_protid_DC_hit_position_region3_fiducial_check[i] && FD_protid_tofmass_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[6]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[6]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[6]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[6]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[6]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[6]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[6]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[6]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[6]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[6]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[6]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[6]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_p*m_p) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[6]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_p*m_p) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[6]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[6]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[6]->Fill(Getdvz(i));
  }
  if(FD_protid_charge_check[i] && FD_protid_DC_hit_position_region1_fiducial_check[i] && FD_protid_DC_hit_position_region2_fiducial_check[i] 
                               && FD_protid_DC_hit_position_region3_fiducial_check[i] && FD_protid_maximum_probability_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[7]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[7]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[7]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[7]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[7]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[7]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[7]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[7]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[7]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[7]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[7]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[7]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_p*m_p) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[7]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_p*m_p) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[7]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[7]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[7]->Fill(Getdvz(i));
  }
  if(FD_protid_default_PID_check[i] && FD_protid_charge_check[i] && FD_protid_DC_hit_position_region1_fiducial_check[i] && FD_protid_DC_hit_position_region2_fiducial_check[i] 
                                    && FD_protid_DC_hit_position_region3_fiducial_check[i] && FD_protid_delta_vz_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[8]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[8]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[8]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[8]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[8]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[8]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[8]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[8]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[8]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[8]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[8]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[8]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_p*m_p) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[8]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_p*m_p) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[8]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[8]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[8]->Fill(Getdvz(i));
  }
  if(FD_protid_all_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[9]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[9]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[9]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[9]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[9]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[9]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[9]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[9]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[9]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[9]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[9]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[9]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_p*m_p) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[9]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_p*m_p) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[9]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[9]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[9]->Fill(Getdvz(i));
  }

  // Neutron

  if(FD_neutrid_default_PID_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[10]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[10]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[10]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[10]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p[10]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec1[10]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec2[10]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec3[10]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec4[10]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec5[10]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec6[10]->Fill(part_p[i], beta_neutr);
    if(part_p[i] >0 && beta_neutr > 0) hist_delta_beta_vs_p[10]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_n*m_n) - beta_neutr);
    if(part_p[i] >0 && beta_neutr > 0) hist_delta_beta[10]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_n*m_n) - beta_neutr);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[10]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[10]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[10]->Fill(Getdvz(i));
  }
  if(FD_neutrid_charge_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[11]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[11]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[11]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[11]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p[11]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec1[11]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec2[11]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec3[11]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec4[11]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec5[11]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec6[11]->Fill(part_p[i], beta_neutr);
    if(part_p[i] >0 && beta_neutr > 0) hist_delta_beta_vs_p[11]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_n*m_n) - beta_neutr);
    if(part_p[i] >0 && beta_neutr > 0) hist_delta_beta[11]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_n*m_n) - beta_neutr);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[11]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[11]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[11]->Fill(Getdvz(i));
  }
  if(FD_neutrid_default_PID_check[i] && FD_neutrid_charge_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[12]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[12]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[12]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[12]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p[12]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec1[12]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec2[12]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec3[12]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec4[12]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec5[12]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec6[12]->Fill(part_p[i], beta_neutr);
    if(part_p[i] >0 && beta_neutr > 0) hist_delta_beta_vs_p[12]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_n*m_n) - beta_neutr);
    if(part_p[i] >0 && beta_neutr > 0) hist_delta_beta[12]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_n*m_n) - beta_neutr);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[12]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[12]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[12]->Fill(Getdvz(i));
  }
  if(FD_neutrid_default_PID_check[i] && FD_neutrid_charge_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[13]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[13]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[13]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[13]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p[13]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec1[13]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec2[13]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec3[13]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec4[13]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec5[13]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec6[13]->Fill(part_p[i], beta_neutr);
    if(part_p[i] >0 && beta_neutr > 0) hist_delta_beta_vs_p[13]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_n*m_n) - beta_neutr);
    if(part_p[i] >0 && beta_neutr > 0) hist_delta_beta[13]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_n*m_n) - beta_neutr);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[13]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[13]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[13]->Fill(Getdvz(i));
  }
  if(FD_neutrid_default_PID_check[i] && FD_neutrid_charge_check[i] && FD_neutrid_beta_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[14]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[14]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[14]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[14]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p[14]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec1[14]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec2[14]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec3[14]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec4[14]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec5[14]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec6[14]->Fill(part_p[i], beta_neutr);
    if(part_p[i] >0 && beta_neutr > 0) hist_delta_beta_vs_p[14]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_n*m_n) - beta_neutr);
    if(part_p[i] >0 && beta_neutr > 0) hist_delta_beta[14]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_n*m_n) - beta_neutr);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[14]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[14]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[14]->Fill(Getdvz(i));
  }
  if(FD_neutrid_default_PID_check[i] && FD_neutrid_charge_check[i] && FD_neutrid_delta_beta_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[15]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[15]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[15]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[15]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p[15]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec1[15]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec2[15]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec3[15]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec4[15]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec5[15]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec6[15]->Fill(part_p[i], beta_neutr);
    if(part_p[i] >0 && beta_neutr > 0) hist_delta_beta_vs_p[15]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_n*m_n) - beta_neutr);
    if(part_p[i] >0 && beta_neutr > 0) hist_delta_beta[15]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_n*m_n) - beta_neutr);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[15]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[15]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[15]->Fill(Getdvz(i));
  }
  if(FD_neutrid_default_PID_check[i] && FD_neutrid_charge_check[i] && FD_neutrid_tofmass_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[16]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[16]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[16]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[16]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p[16]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec1[16]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec2[16]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec3[16]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec4[16]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec5[16]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec6[16]->Fill(part_p[i], beta_neutr);
    if(part_p[i] >0 && beta_neutr > 0) hist_delta_beta_vs_p[16]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_n*m_n) - beta_neutr);
    if(part_p[i] >0 && beta_neutr > 0) hist_delta_beta[16]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_n*m_n) - beta_neutr);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[16]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[16]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[16]->Fill(Getdvz(i));
  }
  if(FD_neutrid_default_PID_check[i] && FD_neutrid_charge_check[i] && FD_neutrid_delta_vz_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[18]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[18]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[18]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[18]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p[18]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec1[18]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec2[18]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec3[18]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec4[18]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec5[18]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec6[18]->Fill(part_p[i], beta_neutr);
    if(part_p[i] >0 && beta_neutr > 0) hist_delta_beta_vs_p[18]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_n*m_n) - beta_neutr);
    if(part_p[i] >0 && beta_neutr > 0) hist_delta_beta[18]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_n*m_n) - beta_neutr);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[18]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[18]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[18]->Fill(Getdvz(i));
  }
  if(FD_neutrid_all_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[19]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[19]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[19]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[19]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p[19]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec1[19]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec2[19]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec3[19]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec4[19]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec5[19]->Fill(part_p[i], beta_neutr);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_sec6[19]->Fill(part_p[i], beta_neutr);
    if(part_p[i] >0 && beta_neutr > 0) hist_delta_beta_vs_p[19]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_n*m_n) - beta_neutr);
    if(part_p[i] >0 && beta_neutr > 0) hist_delta_beta[19]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_n*m_n) - beta_neutr);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[19]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[19]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[19]->Fill(Getdvz(i));
  }

  // pip

  if(FD_pipid_default_PID_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[20]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[20]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[20]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[20]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[20]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[20]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[20]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[20]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[20]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[20]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[20]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[20]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_pip*m_pip) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[20]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_pip*m_pip) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[20]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[20]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[20]->Fill(Getdvz(i));
  }
  if(FD_pipid_charge_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[21]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[21]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[21]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[21]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[21]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[21]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[21]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[21]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[21]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[21]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[21]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[21]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_pip*m_pip) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[21]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_pip*m_pip) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[21]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[21]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[21]->Fill(Getdvz(i));
  }
  if(FD_pipid_default_PID_check[i] && FD_pipid_charge_check[i] && FD_pipid_DC_hit_position_region1_fiducial_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[22]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[22]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[22]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[22]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[22]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[22]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[22]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[22]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[22]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[22]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[22]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[22]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_pip*m_pip) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[22]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_pip*m_pip) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[22]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[22]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[22]->Fill(Getdvz(i));
  }
  if(FD_pipid_default_PID_check[i] && FD_pipid_charge_check[i] && FD_pipid_DC_hit_position_region2_fiducial_check[i]){
    hist_DC_hit_position_region2_hadron_cut_22a->Fill(part_DC_c2x[i], part_DC_c2y[i]);
  }
  if(FD_pipid_default_PID_check[i] && FD_pipid_charge_check[i] && FD_pipid_DC_hit_position_region3_fiducial_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[23]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[23]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[23]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[23]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[23]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[23]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[23]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[23]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[23]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[23]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[23]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[23]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_pip*m_pip) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[23]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_pip*m_pip) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[23]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[23]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[23]->Fill(Getdvz(i));
  }
  if(FD_pipid_default_PID_check[i] && FD_pipid_charge_check[i] && FD_pipid_DC_hit_position_region1_fiducial_check[i] && FD_pipid_DC_hit_position_region2_fiducial_check[i] 
                                   && FD_pipid_DC_hit_position_region3_fiducial_check[i] && FD_pipid_beta_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[24]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[24]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[24]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[24]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[24]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[24]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[24]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[24]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[24]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[24]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[24]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[24]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_pip*m_pip) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[24]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_pip*m_pip) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[24]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[24]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[24]->Fill(Getdvz(i));
  }
  if(FD_pipid_default_PID_check[i] && FD_pipid_charge_check[i] && FD_pipid_DC_hit_position_region1_fiducial_check[i] && FD_pipid_DC_hit_position_region2_fiducial_check[i] 
                                   && FD_pipid_DC_hit_position_region3_fiducial_check[i] && FD_pipid_delta_beta_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[25]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[25]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[25]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[25]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[25]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[25]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[25]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[25]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[25]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[25]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[25]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[25]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_pip*m_pip) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[25]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_pip*m_pip) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[25]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[25]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[25]->Fill(Getdvz(i));
  }
  if(FD_pipid_default_PID_check[i] && FD_pipid_charge_check[i] && FD_pipid_DC_hit_position_region1_fiducial_check[i] && FD_pipid_DC_hit_position_region2_fiducial_check[i] 
                                   && FD_pipid_DC_hit_position_region3_fiducial_check[i] && FD_pipid_tofmass_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[26]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[26]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[26]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[26]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[26]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[26]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[26]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[26]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[26]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[26]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[26]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[26]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_pip*m_pip) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[26]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_pip*m_pip) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[26]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[26]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[26]->Fill(Getdvz(i));
  }
  if(FD_pipid_charge_check[i] && FD_pipid_DC_hit_position_region1_fiducial_check[i] && FD_pipid_DC_hit_position_region2_fiducial_check[i] 
                              && FD_pipid_DC_hit_position_region3_fiducial_check[i] && FD_pipid_maximum_probability_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[27]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[27]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[27]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[27]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[27]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[27]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[27]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[27]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[27]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[27]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[27]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[27]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_pip*m_pip) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[27]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_pip*m_pip) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[27]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[27]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[27]->Fill(Getdvz(i));
  }
  if(FD_pipid_default_PID_check[i] && FD_pipid_charge_check[i] && FD_pipid_DC_hit_position_region1_fiducial_check[i] && FD_pipid_DC_hit_position_region2_fiducial_check[i] 
                                   && FD_pipid_DC_hit_position_region3_fiducial_check[i] && FD_pipid_delta_vz_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[28]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[28]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[28]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[28]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[28]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[28]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[28]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[28]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[28]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[28]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[28]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[28]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_pip*m_pip) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[28]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_pip*m_pip) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[28]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[28]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[28]->Fill(Getdvz(i));
  }
  if(FD_pipid_all_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[29]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[29]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[29]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[29]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[29]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[29]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[29]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[29]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[29]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[29]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[29]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[29]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_pip*m_pip) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[29]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_pip*m_pip) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[29]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[29]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[29]->Fill(Getdvz(i));
  }

  // pim

  if(FD_pimid_default_PID_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[30]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[30]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[30]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[30]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[30]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[30]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[30]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[30]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[30]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[30]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[30]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[30]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_pim*m_pim) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[30]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_pim*m_pim) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[30]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[30]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[30]->Fill(Getdvz(i));
  }
  if(FD_pimid_charge_check[i] && FD_pimid_ele_reject_check[i] && FD_pimid_EC_outer_vs_EC_inner_check[i] && FD_protid_DC_hit_position_region1_fiducial_check[i] && FD_protid_DC_hit_position_region2_fiducial_check[i] && FD_protid_DC_hit_position_region3_fiducial_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[31]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[31]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[31]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[31]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] > 0.3 && beta_charge > 0) hist_beta_vs_p[31]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[31]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[31]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[31]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[31]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[31]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[31]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[31]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_pim*m_pim) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[31]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_pim*m_pim) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[31]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[31]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[31]->Fill(Getdvz(i));
  }
  if(FD_pimid_default_PID_check[i] && FD_pimid_charge_check[i] && FD_pimid_DC_hit_position_region1_fiducial_check[i] && FD_pimid_ele_reject_check[i] && FD_pimid_EC_outer_vs_EC_inner_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[32]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[32]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[32]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[32]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[32]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[32]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[32]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[32]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[32]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[32]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[32]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[32]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_pim*m_pim) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[32]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_pim*m_pim) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[32]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[32]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[32]->Fill(Getdvz(i));
  }
  if(FD_pimid_default_PID_check[i] && FD_pimid_charge_check[i] && FD_pimid_DC_hit_position_region2_fiducial_check[i] && FD_pimid_ele_reject_check[i] && FD_pimid_EC_outer_vs_EC_inner_check[i]){
    hist_DC_hit_position_region2_hadron_cut_32a->Fill(part_DC_c2x[i], part_DC_c2y[i]);
  }
  if(FD_pimid_default_PID_check[i] && FD_pimid_charge_check[i] && FD_pimid_DC_hit_position_region3_fiducial_check[i] && FD_pimid_ele_reject_check[i] && FD_pimid_EC_outer_vs_EC_inner_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[33]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[33]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[33]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[33]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[33]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[33]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[33]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[33]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[33]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[33]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[33]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[33]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_pim*m_pim) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[33]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_pim*m_pim) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[33]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[33]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[33]->Fill(Getdvz(i));
  }
  if(FD_pimid_default_PID_check[i] && FD_pimid_charge_check[i] && FD_pimid_DC_hit_position_region1_fiducial_check[i] && FD_pimid_DC_hit_position_region2_fiducial_check[i] 
                                   && FD_pimid_DC_hit_position_region3_fiducial_check[i] && FD_pimid_ele_reject_check[i] && FD_pimid_EC_outer_vs_EC_inner_check[i] && FD_pimid_beta_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[34]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[34]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[34]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[34]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[34]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[34]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[34]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[34]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[34]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[34]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[34]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[34]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_pim*m_pim) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[34]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_pim*m_pim) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[34]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[34]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[34]->Fill(Getdvz(i));
  }
  if(FD_pimid_default_PID_check[i] && FD_pimid_charge_check[i] && FD_pimid_DC_hit_position_region1_fiducial_check[i] && FD_pimid_DC_hit_position_region2_fiducial_check[i] 
                                   && FD_pimid_DC_hit_position_region3_fiducial_check[i] && FD_pimid_ele_reject_check[i] && FD_pimid_EC_outer_vs_EC_inner_check[i] && FD_pimid_delta_beta_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[35]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[35]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[35]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[35]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[35]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[35]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[35]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[35]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[35]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[35]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[35]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[35]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_pim*m_pim) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[35]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_pim*m_pim) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[35]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[35]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[35]->Fill(Getdvz(i));
  }
  if(FD_pimid_default_PID_check[i] && FD_pimid_charge_check[i] && FD_pimid_DC_hit_position_region1_fiducial_check[i] && FD_pimid_DC_hit_position_region2_fiducial_check[i] 
                                   && FD_pimid_DC_hit_position_region3_fiducial_check[i] && FD_pimid_ele_reject_check[i] && FD_pimid_EC_outer_vs_EC_inner_check[i] && FD_pimid_tofmass_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[36]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[36]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[36]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[36]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[36]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[36]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[36]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[36]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[36]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[36]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[36]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[36]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_pim*m_pim) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[36]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_pim*m_pim) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[36]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[36]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[36]->Fill(Getdvz(i));
  }
  if(FD_pimid_charge_check[i] && FD_pimid_DC_hit_position_region1_fiducial_check[i] && FD_pimid_DC_hit_position_region2_fiducial_check[i] 
                              && FD_pimid_DC_hit_position_region3_fiducial_check[i] && FD_pimid_ele_reject_check[i] && FD_pimid_EC_outer_vs_EC_inner_check[i] && FD_pipid_maximum_probability_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[37]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[37]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[37]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[37]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[37]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[37]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[37]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[37]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[37]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[37]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[37]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[37]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_pim*m_pim) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[37]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_pim*m_pim) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[37]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[37]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[37]->Fill(Getdvz(i));
  }
  if(FD_pimid_default_PID_check[i] && FD_pimid_charge_check[i] && FD_pimid_DC_hit_position_region1_fiducial_check[i] && FD_pimid_DC_hit_position_region2_fiducial_check[i] 
                                   && FD_pimid_DC_hit_position_region3_fiducial_check[i] && FD_pimid_ele_reject_check[i] && FD_pimid_EC_outer_vs_EC_inner_check[i] && FD_pimid_delta_vz_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[38]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[38]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[38]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[38]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[38]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[38]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[38]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[38]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[38]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[38]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[38]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[38]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_pim*m_pim) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[38]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_pim*m_pim) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[38]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[38]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[38]->Fill(Getdvz(i));
  }
  if(FD_pimid_all_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[39]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[39]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[39]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[39]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[39]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[39]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[39]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[39]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[39]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[39]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[39]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[39]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_pim*m_pim) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[39]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_pim*m_pim) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[39]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[39]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[39]->Fill(Getdvz(i));
  }

  // Kp

  if(FD_Kpid_default_PID_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[40]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[40]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[40]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[40]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[40]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[40]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[40]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[40]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[40]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[40]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[40]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[40]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_Kp*m_Kp) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[40]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_Kp*m_Kp) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[40]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[40]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[40]->Fill(Getdvz(i));
  }
  if(FD_Kpid_charge_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[41]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[41]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[41]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[41]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[41]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[41]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[41]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[41]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[41]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[41]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[41]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[41]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_Kp*m_Kp) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[41]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_Kp*m_Kp) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[41]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[41]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[41]->Fill(Getdvz(i));
  }
  if(FD_Kpid_default_PID_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_DC_hit_position_region1_fiducial_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[42]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[42]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[42]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[42]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[42]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[42]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[42]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[42]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[42]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[42]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[42]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[42]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_Kp*m_Kp) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[42]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_Kp*m_Kp) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[42]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[42]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[42]->Fill(Getdvz(i));
  }
  if(FD_Kpid_default_PID_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_DC_hit_position_region2_fiducial_check[i]){
    hist_DC_hit_position_region2_hadron_cut_42a->Fill(part_DC_c2x[i], part_DC_c2y[i]);
  }
  if(FD_Kpid_default_PID_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_DC_hit_position_region3_fiducial_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[43]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[43]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[43]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[43]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[43]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[43]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[43]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[43]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[43]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[43]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[43]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[43]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_Kp*m_Kp) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[43]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_Kp*m_Kp) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[43]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[43]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[43]->Fill(Getdvz(i));
  }
  if(FD_Kpid_default_PID_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_DC_hit_position_region1_fiducial_check[i] && FD_Kpid_DC_hit_position_region2_fiducial_check[i] 
                                  && FD_Kpid_DC_hit_position_region3_fiducial_check[i] && FD_Kpid_beta_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[44]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[44]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[44]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[44]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[44]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[44]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[44]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[44]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[44]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[44]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[44]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[44]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_Kp*m_Kp) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[44]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_Kp*m_Kp) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[44]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[44]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[44]->Fill(Getdvz(i));
  }
  if(FD_Kpid_default_PID_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_DC_hit_position_region1_fiducial_check[i] && FD_Kpid_DC_hit_position_region2_fiducial_check[i] 
                                  && FD_Kpid_DC_hit_position_region3_fiducial_check[i] && FD_Kpid_DC_hit_position_region1_fiducial_check[i] && FD_Kpid_delta_beta_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[45]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[45]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[45]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[45]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[45]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[45]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[45]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[45]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[45]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[45]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[45]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[45]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_Kp*m_Kp) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[45]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_Kp*m_Kp) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[45]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[45]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[45]->Fill(Getdvz(i));
  }
  if(FD_Kpid_default_PID_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_DC_hit_position_region1_fiducial_check[i] && FD_Kpid_DC_hit_position_region2_fiducial_check[i] 
                                  && FD_Kpid_DC_hit_position_region3_fiducial_check[i] && FD_Kpid_tofmass_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[46]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[46]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[46]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[46]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[46]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[46]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[46]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[46]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[46]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[46]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[46]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[46]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_Kp*m_Kp) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[46]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_Kp*m_Kp) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[46]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[46]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[46]->Fill(Getdvz(i));
  }
  if(FD_Kpid_default_PID_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_DC_hit_position_region1_fiducial_check[i] && FD_Kpid_DC_hit_position_region2_fiducial_check[i] 
                                  && FD_Kpid_DC_hit_position_region3_fiducial_check[i]  && FD_Kpid_maximum_probability_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[47]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[47]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[47]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[47]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[47]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[47]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[47]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[47]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[47]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[47]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[47]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[47]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_Kp*m_Kp) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[47]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_Kp*m_Kp) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[47]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[47]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[47]->Fill(Getdvz(i));
  }
  if(FD_Kpid_default_PID_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_DC_hit_position_region1_fiducial_check[i] && FD_Kpid_DC_hit_position_region2_fiducial_check[i] 
                                  && FD_Kpid_DC_hit_position_region3_fiducial_check[i] && FD_Kpid_delta_vz_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[48]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[48]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[48]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[48]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[48]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[48]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[48]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[48]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[48]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[48]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[48]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[48]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_Kp*m_Kp) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[48]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_Kp*m_Kp) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[48]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[48]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[48]->Fill(Getdvz(i));
  }
  if(FD_Kpid_all_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[49]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[49]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[49]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[49]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[49]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[49]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[49]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[49]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[49]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[49]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[49]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[49]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_Kp*m_Kp) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[49]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_Kp*m_Kp) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[49]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[49]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[49]->Fill(Getdvz(i));
  }

  // Km

  if(FD_Kmid_default_PID_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[50]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[50]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[50]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[50]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[50]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[50]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[50]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[50]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[50]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[50]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[50]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[50]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_Km*m_Km) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[50]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_Km*m_Km) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[50]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[50]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[50]->Fill(Getdvz(i));
  }
  if(FD_Kmid_default_PID_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_ele_reject_check[i] && FD_Kmid_EC_outer_vs_EC_inner_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[51]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[51]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[51]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[51]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[51]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[51]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[51]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[51]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[51]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[51]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[51]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[51]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_Km*m_Km) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[51]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_Km*m_Km) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[51]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[51]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[51]->Fill(Getdvz(i));
  }
  if(FD_Kmid_default_PID_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_DC_hit_position_region1_fiducial_check[i] && FD_Kmid_ele_reject_check[i] && FD_Kmid_EC_outer_vs_EC_inner_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[52]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[52]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[52]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[52]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[52]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[52]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[52]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[52]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[52]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[52]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[52]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[52]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_Km*m_Km) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[52]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_Km*m_Km) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[52]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[52]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[52]->Fill(Getdvz(i));
  }
  if(FD_Kmid_default_PID_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_DC_hit_position_region2_fiducial_check[i] && FD_Kmid_ele_reject_check[i] && FD_Kmid_EC_outer_vs_EC_inner_check[i]){
    hist_DC_hit_position_region2_hadron_cut_52a->Fill(part_DC_c2x[i], part_DC_c2y[i]);
  }
  if(FD_Kmid_default_PID_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_DC_hit_position_region3_fiducial_check[i] && FD_Kmid_ele_reject_check[i] && FD_Kmid_EC_outer_vs_EC_inner_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[53]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[53]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[53]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[53]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[53]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[53]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[53]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[53]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[53]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[53]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[53]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[53]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_Km*m_Km) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[53]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_Km*m_Km) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[53]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[53]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[53]->Fill(Getdvz(i));
  }
  if(FD_Kmid_default_PID_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_DC_hit_position_region1_fiducial_check[i] && FD_Kmid_DC_hit_position_region2_fiducial_check[i] 
                                  && FD_Kmid_DC_hit_position_region3_fiducial_check[i] && FD_Kmid_ele_reject_check[i] && FD_Kmid_EC_outer_vs_EC_inner_check[i] && FD_Kmid_beta_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[54]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[54]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[54]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[54]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[54]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[54]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[54]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[54]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[54]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[54]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[54]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[54]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_Km*m_Km) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[54]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_Km*m_Km) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[54]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[54]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[54]->Fill(Getdvz(i));
  }
  if(FD_Kmid_default_PID_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_DC_hit_position_region1_fiducial_check[i] && FD_Kmid_DC_hit_position_region2_fiducial_check[i] 
                                  && FD_Kmid_DC_hit_position_region3_fiducial_check[i] && FD_Kmid_ele_reject_check[i] && FD_Kmid_EC_outer_vs_EC_inner_check[i] && FD_Kmid_delta_beta_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[55]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[55]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[55]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[55]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[55]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[55]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[55]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[55]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[55]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[55]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[55]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[55]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_Km*m_Km) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[55]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_Km*m_Km) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[55]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[55]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[55]->Fill(Getdvz(i));
  }
  if(FD_Kmid_default_PID_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_DC_hit_position_region1_fiducial_check[i] && FD_Kmid_DC_hit_position_region2_fiducial_check[i] 
                                  && FD_Kmid_DC_hit_position_region3_fiducial_check[i] && FD_Kmid_ele_reject_check[i] && FD_Kmid_EC_outer_vs_EC_inner_check[i] && FD_Kmid_tofmass_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[56]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[56]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[56]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[56]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[56]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[56]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[56]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[56]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[56]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[56]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[56]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[56]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_Km*m_Km) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[56]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_Km*m_Km) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[56]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[56]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[56]->Fill(Getdvz(i));
  }
  if(FD_Kmid_default_PID_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_DC_hit_position_region1_fiducial_check[i] && FD_Kmid_DC_hit_position_region2_fiducial_check[i] 
                                  && FD_Kmid_DC_hit_position_region3_fiducial_check[i]  && FD_Kmid_EC_outer_vs_EC_inner_check[i] && FD_Kmid_maximum_probability_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[57]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[57]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[57]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[57]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[57]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[57]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[57]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[57]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[57]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[57]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[57]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[57]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_Km*m_Km) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[57]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_Km*m_Km) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[57]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[57]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[57]->Fill(Getdvz(i));
  }
  if(FD_Kmid_default_PID_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_DC_hit_position_region1_fiducial_check[i] && FD_Kmid_DC_hit_position_region2_fiducial_check[i] 
                                  && FD_Kmid_DC_hit_position_region3_fiducial_check[i] && FD_Kmid_ele_reject_check[i] && FD_Kmid_EC_outer_vs_EC_inner_check[i] && FD_Kmid_delta_vz_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[58]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[58]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[58]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[58]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[58]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[58]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[58]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[58]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[58]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[58]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[58]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[58]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_Km*m_Km) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[58]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_Km*m_Km) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[58]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[58]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[58]->Fill(Getdvz(i));
  }
  if(FD_Kmid_all_check[i]){
    if(part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0) hist_DC_hit_position_region1_hadron[59]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
    if(part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0) hist_DC_hit_position_region2_hadron[59]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
    if(part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0) hist_DC_hit_position_region3_hadron[59]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
    if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_outer_vs_EC_inner_hadron[59]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
    if(part_p[i] >0 && beta_charge > 0) hist_beta_vs_p[59]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec1[59]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec2[59]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec3[59]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec4[59]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec5[59]->Fill(part_p[i], beta_charge);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && beta_charge > 0) hist_beta_vs_p_sec6[59]->Fill(part_p[i], beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta_vs_p[59]->Fill(part_p[i], part_p[i]/sqrt(part_p[i]*part_p[i]+m_Km*m_Km) - beta_charge);
    if(part_p[i] >0 && beta_charge > 0) hist_delta_beta[59]->Fill(part_p[i]/sqrt(part_p[i]*part_p[i]+m_Km*m_Km) - beta_charge);
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass_vs_p[59]->Fill(part_p[i], GetTOFmass2(i, run));
    if(part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_tofmass[59]->Fill(GetTOFmass2(i, run));
    if(Getdvz(i) != 0) hist_delta_vz[59]->Fill(Getdvz(i));
  }



/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////
/// TOF mass

if(part_theta[i]*180/Pival > 5.0 &&  part_theta[i]*180/Pival < 35.0 && part_charge[i] > 0 && part_p[i] >0 && GetTOFmass2(i, run) != 0) hist_TOFmass_FTOF_positive->Fill(GetTOFmass2(i, run));
if(part_theta[i]*180/Pival > 5.0 &&  part_theta[i]*180/Pival < 35.0 && part_charge[i] < 0 && part_p[i] >0 && GetTOFmass2(i, run) != 0 && FD_pimid_ele_reject_check[i] && FD_pimid_EC_outer_vs_EC_inner_check[i]) hist_TOFmass_FTOF_negative->Fill(GetTOFmass2(i, run));
if(part_theta[i]*180/Pival > 35.0 && part_charge[i] > 0 && part_p[i] >0 && GetTOFmass2_CD(i, run) != 0) hist_TOFmass_CTOF_positive->Fill(GetTOFmass2_CD(i, run));
if(part_theta[i]*180/Pival > 35.0 && part_charge[i] < 0 && part_p[i] >0 && GetTOFmass2_CD(i, run) != 0) hist_TOFmass_CTOF_negative->Fill(GetTOFmass2_CD(i, run));

/// ////////////////////////////////////////////////////////////////////////


/// /////////////////////////////////////////////////////////////////////////////////////////
/// Central detector

  double beta_charge_central = Beta_charged_central(i, run);

  // proton

  if(CD_protid_default_PID_check[i]){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[0]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[0]->Fill(Getdvz(i));
  }
  if(CD_protid_charge_check[i]){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[1]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[1]->Fill(Getdvz(i));
  }
  if(CD_protid_charge_check[i] && CD_protid_beta_check[i]){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[2]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[2]->Fill(Getdvz(i));
  }
  if(CD_protid_charge_check[i] && CD_protid_maximum_probability_check[i]){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[3]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[3]->Fill(Getdvz(i));
  }
  if(CD_protid_charge_check[i] && CD_protid_delta_vz_check[i]){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[4]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[4]->Fill(Getdvz(i));
  }
  if(CD_protid_all_check[i]){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[5]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[5]->Fill(Getdvz(i));
  }

  // Neutron

  if(CD_neutrid_default_PID_check[i]){
    if(part_p[i] >0 && beta_neutr > 0) hist_CD_beta_vs_p[10]->Fill(part_p[i], beta_neutr);
    if(Getdvz(i) != 0) hist_CD_delta_vz[10]->Fill(Getdvz(i));
  }
  if(CD_neutrid_charge_check[i]){
    if(part_p[i] >0 && beta_neutr > 0) hist_CD_beta_vs_p[11]->Fill(part_p[i], beta_neutr);
    if(Getdvz(i) != 0) hist_CD_delta_vz[11]->Fill(Getdvz(i));
  }
  if(CD_neutrid_default_PID_check[i] && CD_neutrid_charge_check[i] && CD_neutrid_beta_check[i]){
    if(part_p[i] >0 && beta_neutr > 0) hist_CD_beta_vs_p[12]->Fill(part_p[i], beta_neutr);
    if(Getdvz(i) != 0) hist_CD_delta_vz[12]->Fill(Getdvz(i));
  }
  if(CD_neutrid_default_PID_check[i] && CD_neutrid_charge_check[i] && CD_neutrid_delta_vz_check[i]){
    if(part_p[i] >0 && beta_neutr > 0) hist_CD_beta_vs_p[14]->Fill(part_p[i], beta_neutr);
    if(Getdvz(i) != 0) hist_CD_delta_vz[14]->Fill(Getdvz(i));
  }
  if(CD_neutrid_all_check[i]){
    if(part_p[i] >0 && beta_neutr > 0) hist_CD_beta_vs_p[15]->Fill(part_p[i], beta_neutr);
    if(Getdvz(i) != 0) hist_CD_delta_vz[15]->Fill(Getdvz(i));
  }

  // pip

  if(CD_pipid_default_PID_check[i]){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[20]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[20]->Fill(Getdvz(i));
  }
  if(CD_pipid_charge_check[i]){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[21]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[21]->Fill(Getdvz(i));
  }
  if(CD_pipid_charge_check[i] && CD_pipid_beta_check[i]){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[22]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[22]->Fill(Getdvz(i));
  }
  if(CD_pipid_charge_check[i] && CD_pipid_maximum_probability_check[i]){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[23]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[23]->Fill(Getdvz(i));
  }
  if(CD_pipid_charge_check[i] && CD_pipid_delta_vz_check[i]){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[24]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[24]->Fill(Getdvz(i));
  }
  if(CD_pipid_all_check[i]){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[25]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[25]->Fill(Getdvz(i));
  }

  // pim

  if(CD_pimid_default_PID_check[i]){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[30]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[30]->Fill(Getdvz(i));
  }
  if(CD_pimid_charge_check[i] && (beta_charge_central < 0.9999 || beta_charge_central > 1.0001)){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[31]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[31]->Fill(Getdvz(i));
  }
  if(CD_pimid_charge_check[i] && CD_pimid_beta_check[i] && (beta_charge_central < 0.9999 || beta_charge_central > 1.0001)){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[32]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[32]->Fill(Getdvz(i));
  }
  if(CD_pimid_charge_check[i] && CD_pipid_maximum_probability_check[i] && (beta_charge_central < 0.9999 || beta_charge_central > 1.0001)){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[33]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[33]->Fill(Getdvz(i));
  }
  if(CD_pimid_charge_check[i] && CD_pimid_delta_vz_check[i] && (beta_charge_central < 0.9999 || beta_charge_central > 1.0001)){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[34]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[34]->Fill(Getdvz(i));
  }
  if(CD_pimid_all_check[i]){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[35]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[35]->Fill(Getdvz(i));
  }

  // Kp

  if(CD_Kpid_default_PID_check[i]){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[40]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[40]->Fill(Getdvz(i));
  }
  if(CD_Kpid_charge_check[i]){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[41]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[41]->Fill(Getdvz(i));
  }
  if(CD_Kpid_charge_check[i] && CD_Kpid_beta_check[i]){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[42]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[42]->Fill(Getdvz(i));
  }
  if(CD_Kpid_charge_check[i] && CD_Kpid_maximum_probability_check[i]){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[43]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[43]->Fill(Getdvz(i));
  }
  if(CD_Kpid_charge_check[i] && CD_Kpid_delta_vz_check[i]){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[44]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[44]->Fill(Getdvz(i));
  }
  if(CD_Kpid_all_check[i]){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[45]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[45]->Fill(Getdvz(i));
  }

  // Km

  if(CD_Kmid_default_PID_check[i]){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[50]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[50]->Fill(Getdvz(i));
  }
  if(CD_Kmid_charge_check[i]  && (beta_charge_central < 0.9999 || beta_charge_central > 1.0001)){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[51]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[51]->Fill(Getdvz(i));
  }
  if(CD_Kmid_charge_check[i] && CD_Kmid_beta_check[i] && (beta_charge_central < 0.9999 || beta_charge_central > 1.0001)){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[52]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[52]->Fill(Getdvz(i));
  }
  if(CD_Kmid_charge_check[i] && CD_Kmid_maximum_probability_check[i] && (beta_charge_central < 0.9999 || beta_charge_central > 1.0001)){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[53]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[53]->Fill(Getdvz(i));
  }
  if(CD_Kmid_charge_check[i] && CD_Kmid_delta_vz_check[i] && (beta_charge_central < 0.9999 || beta_charge_central > 1.0001)){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[54]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[54]->Fill(Getdvz(i));
  }
  if(CD_Kmid_all_check[i]){
    if(part_p[i] >0 && beta_charge_central > 0) hist_CD_beta_vs_p[55]->Fill(part_p[i], beta_charge_central);
    if(Getdvz(i) != 0) hist_CD_delta_vz[55]->Fill(Getdvz(i));
  }



} // end of fill hadron pid histogram selection

/// //////////////////////////////////////////////////////////////////////////////////////////////////////
/// photon ID plots
/// /////////////////////////////////////////////////////////////////

if(fill_photon_pid_histograms){

if(FD_photid_default_PID_check[i]){
  if(part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_phot[0]->Fill(part_p[i], beta_neutr);
  if(beta_neutr > 0) hist_beta_phot[0]->Fill(beta_neutr);
  if(part_p[i] > 0) hist_EC_sampling_fraction_phot[0]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
  if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_PCAL_vs_EC_ECAL_phot[0]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
  if(part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0) hist_EC_PCAL_hit_position_phot[0]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
}

if(FD_photid_charge_check[i]){
  if(part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_phot[1]->Fill(part_p[i], beta_neutr);
  if(beta_neutr > 0) hist_beta_phot[1]->Fill(beta_neutr);
  if(part_p[i] > 0) hist_EC_sampling_fraction_phot[1]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
  if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_PCAL_vs_EC_ECAL_phot[1]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
  if(part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0) hist_EC_PCAL_hit_position_phot[1]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
}

if(FD_photid_default_PID_check[i] && FD_photid_charge_check[i] && FD_photid_beta_check[i]){
  if(part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_phot[2]->Fill(part_p[i], beta_neutr);
  if(beta_neutr > 0) hist_beta_phot[2]->Fill(beta_neutr);
  if(part_p[i] > 0) hist_EC_sampling_fraction_phot[2]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
  if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_PCAL_vs_EC_ECAL_phot[2]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
  if(part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0) hist_EC_PCAL_hit_position_phot[2]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
}

if(FD_photid_default_PID_check[i] && FD_photid_charge_check[i] && FD_photid_EC_hit_position_fiducial_check[i]){
  if(part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_phot[3]->Fill(part_p[i], beta_neutr);
  if(beta_neutr > 0) hist_beta_phot[3]->Fill(beta_neutr);
  if(part_p[i] > 0) hist_EC_sampling_fraction_phot[3]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
  if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_PCAL_vs_EC_ECAL_phot[3]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
  if(part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0) hist_EC_PCAL_hit_position_phot[3]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
}

if(FD_photid_all_check[i]){
  if(part_p[i] >0 && beta_neutr > 0) hist_beta_vs_p_phot[4]->Fill(part_p[i], beta_neutr);
  if(beta_neutr > 0) hist_beta_phot[4]->Fill(beta_neutr);
  if(part_p[i] > 0) hist_EC_sampling_fraction_phot[4]->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
  if((part_Cal_ECin_energy[i]+part_Cal_ECout_energy[i])>0 && part_Cal_PCAL_energy[i]>0) hist_EC_PCAL_vs_EC_ECAL_phot[4]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
  if(part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0) hist_EC_PCAL_hit_position_phot[4]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
}

} // end of fill photon ID histogram selection  


/// ////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Fill the FT histograms:

// electrons

if(FT_eid_charge_check[i]){
  hist_FT_FTCAL_energy_vs_radius[0]->Fill(part_FT_energy[i], part_FT_radius[i]);
  hist_FT_FTCAL_hit_position[0]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
  hist_FT_FTTRK_hit_position[0]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
  hist_FT_FTHODO_hit_position[0]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
  hist_FT_beta[0]->Fill(Beta_charged_FT(i, run));
}

if(FT_eid_PID_check[i]){
  hist_FT_FTCAL_energy_vs_radius[1]->Fill(part_FT_energy[i], part_FT_radius[i]);
  hist_FT_FTCAL_hit_position[1]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
  hist_FT_FTTRK_hit_position[1]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
  hist_FT_FTHODO_hit_position[1]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
  hist_FT_beta[1]->Fill(Beta_charged_FT(i, run));
}

if(FT_eid_FTCAL_fiducial_check[i]){
  hist_FT_FTCAL_energy_vs_radius[2]->Fill(part_FT_energy[i], part_FT_radius[i]);
  hist_FT_FTCAL_hit_position[2]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
  hist_FT_FTTRK_hit_position[2]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
  hist_FT_FTHODO_hit_position[2]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
  hist_FT_beta[2]->Fill(Beta_charged_FT(i, run));
}

if(FT_eid_FTTRK_fiducial_check[i]){
  hist_FT_FTCAL_energy_vs_radius[3]->Fill(part_FT_energy[i], part_FT_radius[i]);
  hist_FT_FTCAL_hit_position[3]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
  hist_FT_FTTRK_hit_position[3]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
  hist_FT_FTHODO_hit_position[3]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
  hist_FT_beta[3]->Fill(Beta_charged_FT(i, run));
}

if(FT_eid_FTHODO_fiducial_check[i]){
  hist_FT_FTCAL_energy_vs_radius[4]->Fill(part_FT_energy[i], part_FT_radius[i]);
  hist_FT_FTCAL_hit_position[4]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
  hist_FT_FTTRK_hit_position[4]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
  hist_FT_FTHODO_hit_position[4]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
  hist_FT_beta[4]->Fill(Beta_charged_FT(i, run));
}

if(FT_eid_energy_vs_radius_check[i]){
  hist_FT_FTCAL_energy_vs_radius[5]->Fill(part_FT_energy[i], part_FT_radius[i]);
  hist_FT_FTCAL_hit_position[5]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
  hist_FT_FTTRK_hit_position[5]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
  hist_FT_FTHODO_hit_position[5]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
  hist_FT_beta[5]->Fill(Beta_charged_FT(i, run));
}

if(FT_eid_all_check[i]){
  hist_FT_FTCAL_energy_vs_radius[6]->Fill(part_FT_energy[i], part_FT_radius[i]);
  hist_FT_FTCAL_hit_position[6]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
  hist_FT_FTTRK_hit_position[6]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
  hist_FT_FTHODO_hit_position[6]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
  hist_FT_beta[6]->Fill(Beta_charged_FT(i, run));
}


  hist_FT_FTCAL_energy_vs_radius[7]->Fill(part_FT_energy[i], part_FT_radius[i]);
  hist_FT_FTCAL_hit_position[7]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
  hist_FT_FTTRK_hit_position[7]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
  hist_FT_FTHODO_hit_position[7]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
  hist_FT_beta[7]->Fill(Beta_charged_FT(i, run));

// photons

if(FT_photid_charge_check[i]){
  hist_FT_FTCAL_energy_vs_radius[10]->Fill(part_FT_energy[i], part_FT_radius[i]);
  hist_FT_FTCAL_hit_position[10]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
  hist_FT_FTTRK_hit_position[10]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
  hist_FT_FTHODO_hit_position[10]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
  hist_FT_beta[10]->Fill(Beta_neutral_FT(i, run));
}

if(FT_photid_PID_check[i]){
  hist_FT_FTCAL_energy_vs_radius[11]->Fill(part_FT_energy[i], part_FT_radius[i]);
  hist_FT_FTCAL_hit_position[11]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
  hist_FT_FTTRK_hit_position[11]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
  hist_FT_FTHODO_hit_position[11]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
  hist_FT_beta[11]->Fill(Beta_neutral_FT(i, run));
}

if(FT_photid_FTCAL_fiducial_check[i]){
  hist_FT_FTCAL_energy_vs_radius[12]->Fill(part_FT_energy[i], part_FT_radius[i]);
  hist_FT_FTCAL_hit_position[12]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
  hist_FT_FTTRK_hit_position[12]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
  hist_FT_FTHODO_hit_position[12]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
  hist_FT_beta[12]->Fill(Beta_neutral_FT(i, run));
}

if(FT_photid_beta_check[i]){
  hist_FT_FTCAL_energy_vs_radius[13]->Fill(part_FT_energy[i], part_FT_radius[i]);
  hist_FT_FTCAL_hit_position[13]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
  hist_FT_FTTRK_hit_position[13]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
  hist_FT_FTHODO_hit_position[13]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
  hist_FT_beta[13]->Fill(Beta_neutral_FT(i, run));
}

if(FT_photid_all_check[i]){
  hist_FT_FTCAL_energy_vs_radius[14]->Fill(part_FT_energy[i], part_FT_radius[i]);
  hist_FT_FTCAL_hit_position[14]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
  hist_FT_FTTRK_hit_position[14]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
  hist_FT_FTHODO_hit_position[14]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
  hist_FT_beta[14]->Fill(Beta_neutral_FT(i, run));
}

  hist_FT_FTCAL_energy_vs_radius[15]->Fill(part_FT_energy[i], part_FT_radius[i]);
  hist_FT_FTCAL_hit_position[15]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
  hist_FT_FTTRK_hit_position[15]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
  hist_FT_FTHODO_hit_position[15]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
  hist_FT_beta[15]->Fill(Beta_neutral_FT(i, run));




/// ////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// define calorimeter fiducial cut borders

  if(FD_eid_default_PID_check[i]){ 

    if(part_p[i] > 3 &&  part_p[i] < 6){

      if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_electron_sampfrac_vs_u_coord_sec1->Fill(part_Cal_PCAL_lu[i], part_Cal_energy_total[i]/part_p[i]);
      if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_electron_sampfrac_vs_u_coord_sec2->Fill(part_Cal_PCAL_lu[i], part_Cal_energy_total[i]/part_p[i]);
      if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_electron_sampfrac_vs_u_coord_sec3->Fill(part_Cal_PCAL_lu[i], part_Cal_energy_total[i]/part_p[i]);
      if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_electron_sampfrac_vs_u_coord_sec4->Fill(part_Cal_PCAL_lu[i], part_Cal_energy_total[i]/part_p[i]);
      if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_electron_sampfrac_vs_u_coord_sec5->Fill(part_Cal_PCAL_lu[i], part_Cal_energy_total[i]/part_p[i]);
      if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_electron_sampfrac_vs_u_coord_sec6->Fill(part_Cal_PCAL_lu[i], part_Cal_energy_total[i]/part_p[i]);
      if(part_p[i] > 0) hist_electron_sampfrac_vs_u_coord->Fill(part_Cal_PCAL_lu[i], part_Cal_energy_total[i]/part_p[i]);

      if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_electron_sampfrac_vs_v_coord_sec1->Fill(part_Cal_PCAL_lv[i], part_Cal_energy_total[i]/part_p[i]);
      if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_electron_sampfrac_vs_v_coord_sec2->Fill(part_Cal_PCAL_lv[i], part_Cal_energy_total[i]/part_p[i]);
      if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_electron_sampfrac_vs_v_coord_sec3->Fill(part_Cal_PCAL_lv[i], part_Cal_energy_total[i]/part_p[i]);
      if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_electron_sampfrac_vs_v_coord_sec4->Fill(part_Cal_PCAL_lv[i], part_Cal_energy_total[i]/part_p[i]);
      if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_electron_sampfrac_vs_v_coord_sec5->Fill(part_Cal_PCAL_lv[i], part_Cal_energy_total[i]/part_p[i]);
      if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_electron_sampfrac_vs_v_coord_sec6->Fill(part_Cal_PCAL_lv[i], part_Cal_energy_total[i]/part_p[i]);
      if(part_p[i] > 0) hist_electron_sampfrac_vs_v_coord->Fill(part_Cal_PCAL_lv[i], part_Cal_energy_total[i]/part_p[i]);

      if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_electron_sampfrac_vs_w_coord_sec1->Fill(part_Cal_PCAL_lw[i], part_Cal_energy_total[i]/part_p[i]);
      if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_electron_sampfrac_vs_w_coord_sec2->Fill(part_Cal_PCAL_lw[i], part_Cal_energy_total[i]/part_p[i]);
      if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_electron_sampfrac_vs_w_coord_sec3->Fill(part_Cal_PCAL_lw[i], part_Cal_energy_total[i]/part_p[i]);
      if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_electron_sampfrac_vs_w_coord_sec4->Fill(part_Cal_PCAL_lw[i], part_Cal_energy_total[i]/part_p[i]);
      if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_electron_sampfrac_vs_w_coord_sec5->Fill(part_Cal_PCAL_lw[i], part_Cal_energy_total[i]/part_p[i]);
      if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_electron_sampfrac_vs_w_coord_sec6->Fill(part_Cal_PCAL_lw[i], part_Cal_energy_total[i]/part_p[i]);
      if(part_p[i] > 0) hist_electron_sampfrac_vs_w_coord->Fill(part_Cal_PCAL_lw[i], part_Cal_energy_total[i]/part_p[i]);
    }

    if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_electron_sampfrac_vs_u_coord_vs_p_sec1->Fill(part_Cal_PCAL_lu[i], part_Cal_energy_total[i]/part_p[i], part_p[i]);
    if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_electron_sampfrac_vs_u_coord_vs_p_sec2->Fill(part_Cal_PCAL_lu[i], part_Cal_energy_total[i]/part_p[i], part_p[i]);
    if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_electron_sampfrac_vs_u_coord_vs_p_sec3->Fill(part_Cal_PCAL_lu[i], part_Cal_energy_total[i]/part_p[i], part_p[i]);
    if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_electron_sampfrac_vs_u_coord_vs_p_sec4->Fill(part_Cal_PCAL_lu[i], part_Cal_energy_total[i]/part_p[i], part_p[i]);
    if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_electron_sampfrac_vs_u_coord_vs_p_sec5->Fill(part_Cal_PCAL_lu[i], part_Cal_energy_total[i]/part_p[i], part_p[i]);
    if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_electron_sampfrac_vs_u_coord_vs_p_sec6->Fill(part_Cal_PCAL_lu[i], part_Cal_energy_total[i]/part_p[i], part_p[i]);

    if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_electron_sampfrac_vs_v_coord_vs_p_sec1->Fill(part_Cal_PCAL_lv[i], part_Cal_energy_total[i]/part_p[i], part_p[i]);
    if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_electron_sampfrac_vs_v_coord_vs_p_sec2->Fill(part_Cal_PCAL_lv[i], part_Cal_energy_total[i]/part_p[i], part_p[i]);
    if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_electron_sampfrac_vs_v_coord_vs_p_sec3->Fill(part_Cal_PCAL_lv[i], part_Cal_energy_total[i]/part_p[i], part_p[i]);
    if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_electron_sampfrac_vs_v_coord_vs_p_sec4->Fill(part_Cal_PCAL_lv[i], part_Cal_energy_total[i]/part_p[i], part_p[i]);
    if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_electron_sampfrac_vs_v_coord_vs_p_sec5->Fill(part_Cal_PCAL_lv[i], part_Cal_energy_total[i]/part_p[i], part_p[i]);
    if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_electron_sampfrac_vs_v_coord_vs_p_sec6->Fill(part_Cal_PCAL_lv[i], part_Cal_energy_total[i]/part_p[i], part_p[i]);

    if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_electron_sampfrac_vs_w_coord_vs_p_sec1->Fill(part_Cal_PCAL_lw[i], part_Cal_energy_total[i]/part_p[i], part_p[i]);
    if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_electron_sampfrac_vs_w_coord_vs_p_sec2->Fill(part_Cal_PCAL_lw[i], part_Cal_energy_total[i]/part_p[i], part_p[i]);
    if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_electron_sampfrac_vs_w_coord_vs_p_sec3->Fill(part_Cal_PCAL_lw[i], part_Cal_energy_total[i]/part_p[i], part_p[i]);
    if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_electron_sampfrac_vs_w_coord_vs_p_sec4->Fill(part_Cal_PCAL_lw[i], part_Cal_energy_total[i]/part_p[i], part_p[i]);
    if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_electron_sampfrac_vs_w_coord_vs_p_sec5->Fill(part_Cal_PCAL_lw[i], part_Cal_energy_total[i]/part_p[i], part_p[i]);
    if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_electron_sampfrac_vs_w_coord_vs_p_sec6->Fill(part_Cal_PCAL_lw[i], part_Cal_energy_total[i]/part_p[i], part_p[i]);
  }



  if(part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_momentum->Fill(part_p[i], part_CC_HTCC_nphe[i]);
  if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_momentum_sec1->Fill(part_p[i], part_CC_HTCC_nphe[i]);
  if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_momentum_sec2->Fill(part_p[i], part_CC_HTCC_nphe[i]);
  if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_momentum_sec3->Fill(part_p[i], part_CC_HTCC_nphe[i]);
  if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_momentum_sec4->Fill(part_p[i], part_CC_HTCC_nphe[i]);
  if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_momentum_sec5->Fill(part_p[i], part_CC_HTCC_nphe[i]);
  if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_momentum_sec6->Fill(part_p[i], part_CC_HTCC_nphe[i]);

  if(FD_protid_all_check[i]){
    if(part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_momentum_prot->Fill(part_p[i], part_CC_HTCC_nphe[i]);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_momentum_prot_sec1->Fill(part_p[i], part_CC_HTCC_nphe[i]);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_momentum_prot_sec2->Fill(part_p[i], part_CC_HTCC_nphe[i]);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_momentum_prot_sec3->Fill(part_p[i], part_CC_HTCC_nphe[i]);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_momentum_prot_sec4->Fill(part_p[i], part_CC_HTCC_nphe[i]);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_momentum_prot_sec5->Fill(part_p[i], part_CC_HTCC_nphe[i]);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_momentum_prot_sec6->Fill(part_p[i], part_CC_HTCC_nphe[i]);
  }

  if(FD_pipid_all_check[i]){
    if(part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_momentum_pip->Fill(part_p[i], part_CC_HTCC_nphe[i]);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_momentum_pip_sec1->Fill(part_p[i], part_CC_HTCC_nphe[i]);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_momentum_pip_sec2->Fill(part_p[i], part_CC_HTCC_nphe[i]);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_momentum_pip_sec3->Fill(part_p[i], part_CC_HTCC_nphe[i]);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_momentum_pip_sec4->Fill(part_p[i], part_CC_HTCC_nphe[i]);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_momentum_pip_sec5->Fill(part_p[i], part_CC_HTCC_nphe[i]);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_momentum_pip_sec6->Fill(part_p[i], part_CC_HTCC_nphe[i]);
  }



  if(part_p[i] > 4){

    if(part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_beta->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
    if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_beta_sec1->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
    if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_beta_sec2->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
    if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_beta_sec3->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
    if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_beta_sec4->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
    if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_beta_sec5->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
    if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_beta_sec6->Fill(part_CC_HTCC_nphe[i], part_beta[i]);

    if(FD_protid_all_check[i]){
      if(part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_beta_prot->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
      if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_beta_prot_sec1->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
      if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_beta_prot_sec2->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
      if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_beta_prot_sec3->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
      if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_beta_prot_sec4->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
      if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_beta_prot_sec5->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
      if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_beta_prot_sec6->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
    }

    if(FD_pipid_all_check[i]){
      if(part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_beta_pip->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
      if(part_FTOF_sector_layer2[i] == 1 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_beta_pip_sec1->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
      if(part_FTOF_sector_layer2[i] == 2 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_beta_pip_sec2->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
      if(part_FTOF_sector_layer2[i] == 3 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_beta_pip_sec3->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
      if(part_FTOF_sector_layer2[i] == 4 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_beta_pip_sec4->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
      if(part_FTOF_sector_layer2[i] == 5 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_beta_pip_sec5->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
      if(part_FTOF_sector_layer2[i] == 6 && part_p[i] >0 && part_CC_HTCC_nphe[i] > 0) hist_HTCC_Nphe_vs_beta_pip_sec6->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
    }
  }


  if(FD_eid_charge_check[i] && FD_eid_EC_outer_vs_EC_inner_check[i] && FD_eid_DC_z_vertex_check[i] && FD_eid_DC_hit_position_region1_fiducial_check[i] && FD_eid_DC_hit_position_region2_fiducial_check[i] 
     && FD_eid_DC_hit_position_region3_fiducial_check[i] && FD_eid_EC_hit_position_fiducial_check[i]){  

    if(part_Cal_PCAL_sector[i] == 1 && part_Cal_energy_total[i] > 0) hist_sampling_fraction_vs_E_sec1->Fill(part_Cal_energy_total[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 2 && part_Cal_energy_total[i] > 0) hist_sampling_fraction_vs_E_sec2->Fill(part_Cal_energy_total[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 3 && part_Cal_energy_total[i] > 0) hist_sampling_fraction_vs_E_sec3->Fill(part_Cal_energy_total[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 4 && part_Cal_energy_total[i] > 0) hist_sampling_fraction_vs_E_sec4->Fill(part_Cal_energy_total[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 5 && part_Cal_energy_total[i] > 0) hist_sampling_fraction_vs_E_sec5->Fill(part_Cal_energy_total[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 6 && part_Cal_energy_total[i] > 0) hist_sampling_fraction_vs_E_sec6->Fill(part_Cal_energy_total[i], part_Cal_energy_total[i]/part_p[i]);

    if(part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0) hist_sampling_fraction_vs_p_sec1->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0) hist_sampling_fraction_vs_p_sec2->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0) hist_sampling_fraction_vs_p_sec3->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0) hist_sampling_fraction_vs_p_sec4->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0) hist_sampling_fraction_vs_p_sec5->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);
    if(part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0) hist_sampling_fraction_vs_p_sec6->Fill(part_p[i], part_Cal_energy_total[i]/part_p[i]);

  }


} // end of BUFFER loop for histogram filling



/// ////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Fill beta and time plots per component for FTOF correction
///


for(Int_t i = 0; i < BUFFER; i++){ 

  double beta_charge = Beta_charged(i, run);
  double Tstart = Get_Starttime(i, run);

  for(Int_t j = 0; j < 23; j++){ 
    if(beta_charge > 0 && part_charge[i] == +1 && part_FTOF_component_layer1[i] == j+1) hist_beta_vs_p_positive_layer1_comp[j]->Fill(part_p[i], beta_charge); 
    if(beta_charge > 0 && part_charge[i] == -1 && part_FTOF_component_layer1[i] == j+1) hist_beta_vs_p_negative_layer1_comp[j]->Fill(part_p[i], beta_charge);
  }

  for(Int_t j = 0; j < 62; j++){ 
    if(beta_charge > 0 && part_charge[i] == +1 && part_FTOF_component_layer2[i] == j+1) hist_beta_vs_p_positive_layer2_comp[j]->Fill(part_p[i], beta_charge); 
    if(beta_charge > 0 && part_charge[i] == -1 && part_FTOF_component_layer2[i] == j+1) hist_beta_vs_p_negative_layer2_comp[j]->Fill(part_p[i], beta_charge);
  }

  for(Int_t j = 0; j < 5; j++){ 
    if(beta_charge > 0 && part_charge[i] == +1 && part_FTOF_component_layer3[i] == j+1) hist_beta_vs_p_positive_layer3_comp[j]->Fill(part_p[i], beta_charge); 
    if(beta_charge > 0 && part_charge[i] == -1 && part_FTOF_component_layer3[i] == j+1) hist_beta_vs_p_negative_layer3_comp[j]->Fill(part_p[i], beta_charge);
  }


  for(Int_t j = 0; j < 62; j++){ 
    if(part_FTOF_time[i] != 0 && Get_Starttime(i, run) != 0 && part_charge[i] == +1 && part_FTOF_component_layer2[i] == j+1){
      hist_time_vs_p_positive_comp[j]->Fill(part_p[i], part_FTOF_time[i]-Tstart); 
    }
    if(part_FTOF_time[i] != 0 && Get_Starttime(i, run) != 0 && part_charge[i] == -1 && part_FTOF_component_layer2[i] == j+1){
      hist_time_vs_p_negative_comp[j]->Fill(part_p[i], part_FTOF_time[i]-Tstart);
    }
  }

}


/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////
/// TOF monitoring

for(Int_t i = 0; i < BUFFER; i++){

  double Tstart = 0;
  double momentum = 0;
  double time_layer1 = 0;
  double time_layer2 = 0;
  double time_layer3 = 0;
  double time_CTOF = 0;
  double path_layer1 = 0;
  double path_layer2 = 0;
  double path_layer3 = 0;
  double path_CTOF = 0;
  double component_layer1 = 0;
  double component_layer2 = 0;
  double component_layer3 = 0;
  double component_CTOF = 0;
  double sector_layer1 = 0;
  double sector_layer2 = 0;
  double sector_layer3 = 0;
  double layer = 0;
  double caltime = 0;

  Tstart = Get_Starttime(i, run);
  momentum = part_p[i];
  layer =  part_FTOF_layer[i];

  time_layer1 = part_FTOF_time_layer1[i];
  time_layer2 = part_FTOF_time[i];
  time_layer3 = part_FTOF_time_layer3[i];
  time_CTOF = part_CTOF_time[i];

  path_layer1 = part_FTOF_path_layer1[i];
  path_layer2 = part_FTOF_path[i];
  path_layer3 = part_FTOF_path_layer3[i];
  path_CTOF = part_CTOF_path[i];

  component_layer1 = part_FTOF_component_layer1[i]; 
  component_layer2 = part_FTOF_component_layer2[i]; 
  component_layer3 = part_FTOF_component_layer3[i]; 
  component_CTOF = part_CTOF_component[i];

  sector_layer1 = part_FTOF_sector_layer1[i];
  sector_layer2 = part_FTOF_sector_layer2[i];
  sector_layer3 = part_FTOF_sector_layer3[i];

  
  for(Int_t k = 0; k < 6; k++){ 
    if(sector_layer1 == k+1){
      for(Int_t j = 0; j < 23; j++){
        if(component_layer1 == j+1){
          if(FD_protid_all_check[i] == true && time_layer1 > 0 && Tstart > 0){
            caltime = 10000000 * path_layer1 * sqrt(pow(momentum,2) + pow(m_p,2))/(momentum*c);
            hist_deltaT_proton_layer1[k][j]->Fill(time_layer1-Tstart-caltime); 
            hist_deltaT_all_combined_layer1[k][j]->Fill(time_layer1-Tstart-caltime); 
          }
          if(FD_pipid_all_check[i] == true && time_layer1 > 0 && Tstart > 0){  
            caltime = 10000000 * path_layer1 * sqrt(pow(momentum,2) + pow(m_pip,2))/(momentum*c);
            hist_deltaT_pip_layer1[k][j]->Fill(time_layer1-Tstart-caltime); 
            hist_deltaT_pion_combined_layer1[k][j]->Fill(time_layer1-Tstart-caltime); 
            hist_deltaT_all_combined_layer1[k][j]->Fill(time_layer1-Tstart-caltime); 
          }
          if(FD_pimid_all_check[i] == true && time_layer1 > 0 && Tstart > 0){  
            caltime = 10000000 * path_layer1 * sqrt(pow(momentum,2) + pow(m_pim,2))/(momentum*c);
            hist_deltaT_pim_layer1[k][j]->Fill(time_layer1-Tstart-caltime); 
            hist_deltaT_pion_combined_layer1[k][j]->Fill(time_layer1-Tstart-caltime); 
            hist_deltaT_all_combined_layer1[k][j]->Fill(time_layer1-Tstart-caltime); 
          }
        }
      }
    }
  }

  for(Int_t k = 0; k < 6; k++){ 
    if(sector_layer2 == k+1){
      for(Int_t j = 0; j < 62; j++){ 
        if(component_layer2 == j+1){
          if(FD_protid_all_check[i] == true && time_layer2 > 0 && Tstart > 0){
            caltime =  10000000 * path_layer2 * sqrt(pow(momentum,2) + pow(m_p,2))/(momentum*c);
            hist_deltaT_proton_layer2[k][j]->Fill(time_layer2-Tstart-caltime); 
            hist_deltaT_all_combined_layer2[k][j]->Fill(time_layer2-Tstart-caltime); 
          }
          if(FD_pipid_all_check[i] == true && time_layer2 > 0 && Tstart > 0){  
            caltime = 10000000 * path_layer2 * sqrt(pow(momentum,2) + pow(m_pip,2))/(momentum*c);;
            hist_deltaT_pip_layer2[k][j]->Fill(time_layer2-Tstart-caltime);
            hist_deltaT_pion_combined_layer2[k][j]->Fill(time_layer2-Tstart-caltime); 
            hist_deltaT_all_combined_layer2[k][j]->Fill(time_layer2-Tstart-caltime); 
          }
          if(FD_pimid_all_check[i] == true && time_layer2 > 0 && Tstart > 0){  
            caltime = 10000000 * path_layer2 * sqrt(pow(momentum,2) + pow(m_pim,2))/(momentum*c);
            hist_deltaT_pim_layer2[k][j]->Fill(time_layer2-Tstart-caltime); 
            hist_deltaT_pion_combined_layer2[k][j]->Fill(time_layer2-Tstart-caltime);
            hist_deltaT_all_combined_layer2[k][j]->Fill(time_layer2-Tstart-caltime); 
          }
        }
      }
    }
  }

  for(Int_t k = 0; k < 6; k++){ 
    if(sector_layer3 == k+1){
      for(Int_t j = 0; j < 5; j++){ 
        if(component_layer3 == j+1){
          if(FD_protid_all_check[i] == true && time_layer3 > 0 && Tstart > 0){ 
            caltime =  10000000 * path_layer3 * sqrt(pow(momentum,2) + pow(m_p,2))/(momentum*c);
            hist_deltaT_proton_layer3[k][j]->Fill(time_layer3-Tstart-caltime);
            hist_deltaT_all_combined_layer3[k][j]->Fill(time_layer3-Tstart-caltime);  
          }
          if(FD_pipid_all_check[i] == true && time_layer3 > 0 && Tstart > 0){
            caltime = 10000000 * path_layer3 * sqrt(pow(momentum,2) + pow(m_pip,2))/(momentum*c);
            hist_deltaT_pip_layer3[k][j]->Fill(time_layer3-Tstart-caltime); 
            hist_deltaT_pion_combined_layer3[k][j]->Fill(time_layer3-Tstart-caltime); 
            hist_deltaT_all_combined_layer3[k][j]->Fill(time_layer3-Tstart-caltime); 
          }
          if(FD_pimid_all_check[i] == true && time_layer3 > 0 && Tstart > 0){ 
            caltime = 10000000 * path_layer3 * sqrt(pow(momentum,2) + pow(m_pim,2))/(momentum*c);
            hist_deltaT_pim_layer3[k][j]->Fill(time_layer3-Tstart-caltime);
            hist_deltaT_pion_combined_layer3[k][j]->Fill(time_layer3-Tstart-caltime); 
            hist_deltaT_all_combined_layer3[k][j]->Fill(time_layer3-Tstart-caltime); 
          }
        }
      }
    }
  }

  for(Int_t j = 0; j < 48; j++){ 
    if(component_CTOF == j+1){
      if(CD_protid_all_check[i] == true && time_CTOF > 0 && Tstart > 0){
        caltime =  10000000 * path_CTOF * sqrt(pow(momentum,2) + pow(m_p,2))/(momentum*c);
        hist_deltaT_proton_CTOF[j]->Fill(time_CTOF-Tstart-caltime); 
        hist_deltaT_all_combined_CTOF[j]->Fill(time_CTOF-Tstart-caltime);
      }
      if(CD_pipid_all_check[i] == true && time_CTOF > 0 && Tstart > 0){
        caltime = 10000000 * path_CTOF * sqrt(pow(momentum,2) + pow(m_pip,2))/(momentum*c);
        hist_deltaT_pip_CTOF[j]->Fill(time_CTOF-Tstart-caltime);  
        hist_deltaT_pion_combined_CTOF[j]->Fill(time_CTOF-Tstart-caltime);
        hist_deltaT_all_combined_CTOF[j]->Fill(time_CTOF-Tstart-caltime);
      }
      if(CD_pimid_all_check[i] == true && time_CTOF > 0 && Tstart > 0){ 
        caltime = 10000000 * path_CTOF * sqrt(pow(momentum,2) + pow(m_pim,2))/(momentum*c);
        hist_deltaT_pim_CTOF[j]->Fill(time_CTOF-Tstart-caltime);  
        hist_deltaT_pion_combined_CTOF[j]->Fill(time_CTOF-Tstart-caltime);
        hist_deltaT_all_combined_CTOF[j]->Fill(time_CTOF-Tstart-caltime);
      }
    }
  }

}

/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// W monitoring
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

double W_theta_bin_min[] = {2.5, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35};
double W_theta_bin_max[] = {4.5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37};

for(Int_t i = 0; i < e_count; i++){

  double kinW = kin_W(p4_ele_raw[i]);
 
  for(Int_t j = 0; j < 17; j++){ 
    if(p4_ele[i].Theta()*180/Pival > W_theta_bin_min[j] && p4_ele[i].Theta()*180/Pival < W_theta_bin_max[j]){
      if(e_FTOF_sec[i] == 1 && kinW > 0) hist_W_binned_sec_01[j]->Fill(kinW); 
      if(e_FTOF_sec[i] == 2 && kinW > 0) hist_W_binned_sec_02[j]->Fill(kinW); 
      if(e_FTOF_sec[i] == 3 && kinW > 0) hist_W_binned_sec_03[j]->Fill(kinW); 
      if(e_FTOF_sec[i] == 4 && kinW > 0) hist_W_binned_sec_04[j]->Fill(kinW); 
      if(e_FTOF_sec[i] == 5 && kinW > 0) hist_W_binned_sec_05[j]->Fill(kinW); 
      if(e_FTOF_sec[i] == 6 && kinW > 0) hist_W_binned_sec_06[j]->Fill(kinW); 
    }
  }
}

/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// momentum correction histograms
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////


for(Int_t i = 0; i < e_count; i++){ 


  if(e_FTOF_sec[i] == 1){hist_all_electron_p_vs_theta_sec1->Fill(p4_ele[i].Theta()*180/Pival, p4_ele[i].P());}
  if(e_FTOF_sec[i] == 2){hist_all_electron_p_vs_theta_sec2->Fill(p4_ele[i].Theta()*180/Pival, p4_ele[i].P());}  
  if(e_FTOF_sec[i] == 3){hist_all_electron_p_vs_theta_sec3->Fill(p4_ele[i].Theta()*180/Pival, p4_ele[i].P());}
  if(e_FTOF_sec[i] == 4){hist_all_electron_p_vs_theta_sec4->Fill(p4_ele[i].Theta()*180/Pival, p4_ele[i].P());}
  if(e_FTOF_sec[i] == 5){hist_all_electron_p_vs_theta_sec5->Fill(p4_ele[i].Theta()*180/Pival, p4_ele[i].P());}
  if(e_FTOF_sec[i] == 6){hist_all_electron_p_vs_theta_sec6->Fill(p4_ele[i].Theta()*180/Pival, p4_ele[i].P());}



  double kinW = kin_W(p4_ele_raw[i]);

  if((e_FTOF_sec[i] == 1 || e_FTOF_sec[i] == 2 || e_FTOF_sec[i] == 3 || e_FTOF_sec[i] == 4 || e_FTOF_sec[i] == 5 || e_FTOF_sec[i] == 6) && kinW > 0) hist_W_raw->Fill(kinW); 
  if(e_FTOF_sec[i] == 1 && kinW > 0) hist_W_raw_sec_01->Fill(kinW); 
  if(e_FTOF_sec[i] == 2 && kinW > 0) hist_W_raw_sec_02->Fill(kinW); 
  if(e_FTOF_sec[i] == 3 && kinW > 0) hist_W_raw_sec_03->Fill(kinW); 
  if(e_FTOF_sec[i] == 4 && kinW > 0) hist_W_raw_sec_04->Fill(kinW); 
  if(e_FTOF_sec[i] == 5 && kinW > 0) hist_W_raw_sec_05->Fill(kinW); 
  if(e_FTOF_sec[i] == 6 && kinW > 0) hist_W_raw_sec_06->Fill(kinW); 
  if(kinW > 0) hist_W_vs_phi_raw->Fill(p4_ele_raw[i].Phi()*180/Pival, kinW); 

  double kinW_corr = kin_W(p4_ele[i]);

  if(p4_ele[i].Theta()*180/Pival > 5 && kinW_corr > 0) hist_W_corr->Fill(kinW_corr); 
  if(e_FTOF_sec[i] == 1 && kinW_corr > 0) hist_W_corr_sec_01->Fill(kinW_corr); 
  if(e_FTOF_sec[i] == 2 && kinW_corr > 0) hist_W_corr_sec_02->Fill(kinW_corr); 
  if(e_FTOF_sec[i] == 3 && kinW_corr > 0) hist_W_corr_sec_03->Fill(kinW_corr); 
  if(e_FTOF_sec[i] == 4 && kinW_corr > 0) hist_W_corr_sec_04->Fill(kinW_corr); 
  if(e_FTOF_sec[i] == 5 && kinW_corr > 0) hist_W_corr_sec_05->Fill(kinW_corr); 
  if(e_FTOF_sec[i] == 6 && kinW_corr > 0) hist_W_corr_sec_06->Fill(kinW_corr); 
  if(kinW_corr > 0) hist_W_vs_phi_corr->Fill(p4_ele[i].Phi()*180/Pival, kinW_corr); 

  double kinW_max = 0;

  if(Ebeam < 3) kinW_max = 1.07;
  if(Ebeam > 5 && Ebeam < 8){
    if(e_FTOF_sec[i] == 1) kinW_max = 1.35;
    if(e_FTOF_sec[i] == 2) kinW_max = 1.30;
    if(e_FTOF_sec[i] == 3) kinW_max = 1.30;
    if(e_FTOF_sec[i] == 4) kinW_max = 1.35;
    if(e_FTOF_sec[i] == 5) kinW_max = 1.46;
    if(e_FTOF_sec[i] == 6) kinW_max = 1.45;
  }
  if(Ebeam > 9) kinW_max = 1.6;

  if(kinW < kinW_max){

    hist_theta_vs_phi->Fill(p4_ele_raw[i].Phi()*180/Pival, p4_ele_raw[i].Theta()*180/Pival);

    double pcorr = Ebeam/(1+(2*Ebeam*pow(sin(p4_ele_raw[i].Theta())/2,2))/m_p);

    TAxis *xaxis = hist_theta_vs_phi->GetXaxis();
    TAxis *yaxis = hist_theta_vs_phi->GetYaxis();
    Int_t binx = xaxis->FindBin(p4_ele_raw[i].Phi()*180/Pival);
    Int_t biny = yaxis->FindBin(p4_ele_raw[i].Theta()*180/Pival);

    for(Int_t k = 0; k < bincount_theta; k++){
      for(Int_t j = 0; j < bincount_phi; j++){
        if(j == binx-1 && k == biny-1) hist_delta_P[k][j]->Fill(pcorr/p4_ele_raw[i].P());
      }
    }
  }

  // for cross check:  Residuals after correction

  if(kinW < kinW_max){

    double pcorr = Ebeam/(1+(2*Ebeam*pow(sin(p4_ele[i].Theta())/2,2))/m_p);

    TAxis *xaxis = hist_theta_vs_phi->GetXaxis();
    TAxis *yaxis = hist_theta_vs_phi->GetYaxis();
    Int_t binx = xaxis->FindBin(p4_ele[i].Phi()*180/Pival);
    Int_t biny = yaxis->FindBin(p4_ele[i].Theta()*180/Pival);

    for(Int_t k = 0; k < bincount_theta; k++){
      for(Int_t j = 0; j < bincount_phi; j++){
        if(j == binx-1 && k == biny-1) hist_delta_P_corr[k][j]->Fill(pcorr/p4_ele[i].P());
      }
    }
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////
  /// Theta angle correction for electrons

  if(kinW < kinW_max){

    hist_p_vs_phi->Fill(p4_ele[i].Phi()*180/Pival, p4_ele[i].P());

    double theta = 2 * asin(sqrt((m_p/(2*p4_ele[i].P()))-(m_p/(2*Ebeam))));

    TAxis *xaxis = hist_p_vs_phi->GetXaxis();
    TAxis *yaxis = hist_p_vs_phi->GetYaxis();
    Int_t binx = xaxis->FindBin(p4_ele[i].Phi()*180/Pival);
    Int_t biny = yaxis->FindBin(p4_ele[i].P());

    for(Int_t k = 0; k < bincount_p; k++){
      for(Int_t j = 0; j < bincount_phi; j++){
        if(j == binx-1 && k == biny-1) hist_delta_Theta[k][j]->Fill(theta/p4_ele[i].Theta());
      }
    }

  }

}


/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///  fill the statistics:
/// ///////////////////////////////////////////////////////////////

  hist_electron_count->Fill(e_count);
  hist_proton_count->Fill(p_count);
  hist_neutron_count->Fill(n_count);
  hist_pip_count->Fill(pip_count);
  hist_pim_count->Fill(pim_count);
  hist_Kp_count->Fill(Kp_count);
  hist_Km_count->Fill(Km_count);
  hist_photon_count->Fill(g_count);


/// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
}   // end of event loop
/// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// cut stiatics

if(show_FD_eid_statistics){
cout << endl;
cout << "------------------------------------------------------------------------------------------" << endl;
cout << "electron PID cut statistics " << endl;
cout << endl;
if(FD_eid_default_PID_pass != 0 && FD_eid_charge_pass != 0){
cout << "number of particles which passed default electron reference PID: "  	<< FD_eid_default_PID_pass << "   ( " << 100*FD_eid_default_PID_pass/FD_eid_charge_pass << " \% of neg. particles and " << 100*FD_eid_default_PID_pass/FD_eid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed CC nphe cut (not in chain): " 	<< FD_eid_CC_nphe_pass << "   ( " << 100*FD_eid_CC_nphe_pass/FD_eid_charge_pass << " \% of neg. particles and " << 100*FD_eid_CC_nphe_pass/FD_eid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed EC PCAL vs ECAL cut: " 		<< FD_eid_EC_outer_vs_EC_inner_pass << "   ( " << 100*FD_eid_EC_outer_vs_EC_inner_pass/FD_eid_charge_pass << " \% of neg. particles and " << 100*FD_eid_EC_outer_vs_EC_inner_pass/FD_eid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed EC sampling fraction cut: " 		<< FD_eid_EC_sampling_fraction_pass << "   ( " << 100*FD_eid_EC_sampling_fraction_pass/FD_eid_charge_pass << " \% of neg. particles and " << 100*FD_eid_EC_sampling_fraction_pass/FD_eid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed EC hit position fid. cut: " 		<< FD_eid_EC_hit_position_fiducial_pass << "   ( " << 100*FD_eid_EC_hit_position_fiducial_pass/FD_eid_charge_pass << " \% of neg. particles and " << 100*FD_eid_EC_hit_position_fiducial_pass/FD_eid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed DC region 1 fid. cut: " 	    	<< FD_eid_DC_hit_position_region1_fiducial_pass << "   ( " << 100*FD_eid_DC_hit_position_region1_fiducial_pass/FD_eid_charge_pass << " \% of neg. particles and " << 100*FD_eid_DC_hit_position_region1_fiducial_pass/FD_eid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed DC region 2 fid. cut: " 	    	<< FD_eid_DC_hit_position_region2_fiducial_pass << "   ( " << 100*FD_eid_DC_hit_position_region2_fiducial_pass/FD_eid_charge_pass << " \% of neg. particles and " << 100*FD_eid_DC_hit_position_region2_fiducial_pass/FD_eid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed DC region 3 fid. cut: " 		<< FD_eid_DC_hit_position_region3_fiducial_pass << "   ( " << 100*FD_eid_DC_hit_position_region3_fiducial_pass/FD_eid_charge_pass << " \% of neg. particles and " << 100*FD_eid_DC_hit_position_region3_fiducial_pass/FD_eid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed z vertex cut: " 	   		<< FD_eid_DC_z_vertex_pass << "   ( " << 100*FD_eid_DC_z_vertex_pass/FD_eid_charge_pass << " \% of neg. particles and " << 100*FD_eid_DC_z_vertex_pass/FD_eid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "__________________________________________________________________________________________________________________________________________________________" << endl;
cout << "number of particles which passed all FD_eid cuts: " 	   << FD_eid_all_pass << "   ( " << 100*FD_eid_all_pass/FD_eid_charge_pass << " \% of neg. particles and " << 100*FD_eid_all_pass/FD_eid_default_PID_pass << " \% of reference PID particles)" << endl;
}
else cout << "No electrons detected!" << endl;
cout << endl;
}

if(show_charged_ID_statistics){
cout << "------------------------------------------------------------------------------------------" << endl;
cout << "charged PID cut statistics " << endl;
cout << endl;
cout << "a) Protons " << endl;
if(FD_protid_default_PID_pass != 0  && FD_protid_charge_pass != 0){
cout << "number of particles which passed default proton reference PID: "  	<< FD_protid_default_PID_pass << "   ( " << 100*FD_protid_default_PID_pass/FD_protid_charge_pass << " \% of pos. particles and " << 100*FD_protid_default_PID_pass/FD_protid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed DC region 1 fid. cut: " 	    	<< FD_protid_DC_hit_position_region1_fiducial_pass << "   ( " << 100*FD_protid_DC_hit_position_region1_fiducial_pass/FD_protid_charge_pass << " \% of pos. particles and " << 100*FD_protid_DC_hit_position_region1_fiducial_pass/FD_protid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed DC region 2 fid. cut: " 	    	<< FD_protid_DC_hit_position_region2_fiducial_pass << "   ( " << 100*FD_protid_DC_hit_position_region2_fiducial_pass/FD_protid_charge_pass << " \% of pos. particles and " << 100*FD_protid_DC_hit_position_region2_fiducial_pass/FD_protid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed DC region 3 fid. cut: " 	    	<< FD_protid_DC_hit_position_region3_fiducial_pass << "   ( " << 100*FD_protid_DC_hit_position_region3_fiducial_pass/FD_protid_charge_pass << " \% of pos. particles and " << 100*FD_protid_DC_hit_position_region3_fiducial_pass/FD_protid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed beta cut: " 		<< FD_protid_beta_pass << "   ( " << 100*FD_protid_beta_pass/FD_protid_charge_pass << " \% of pos. particles and " << 100*FD_protid_beta_pass/FD_protid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed delta beta cut: " 		<< FD_protid_delta_beta_pass << "   ( " << 100*FD_protid_delta_beta_pass/FD_protid_charge_pass << " \% of pos. particles and " << 100*FD_protid_delta_beta_pass/FD_protid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed tofmass cut: " 		<< FD_protid_tofmass_pass << "   ( " << 100*FD_protid_tofmass_pass/FD_protid_charge_pass << " \% of pos. particles and " << 100*FD_protid_tofmass_pass/FD_protid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed maximum probability cut: " 		<< FD_protid_maximum_probability_pass << "   ( " << 100*FD_protid_maximum_probability_pass/FD_protid_charge_pass << " \% of pos. particles and " << 100*FD_protid_maximum_probability_pass/FD_protid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed delta vz cut: " 		<< FD_protid_delta_vz_pass << "   ( " << 100*FD_protid_delta_vz_pass/FD_protid_charge_pass << " \% of pos. particles and " << 100*FD_protid_delta_vz_pass/FD_protid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "__________________________________________________________________________________________________________________________________________________________" << endl;
cout << "number of particles which passed all proton ID cuts: " 	   << FD_protid_all_pass << "   ( " << 100*FD_protid_all_pass/FD_protid_charge_pass << " \% of pos. particles and " << 100*FD_protid_all_pass/FD_protid_default_PID_pass << " \% of reference PID particles)" << endl;
}
else cout << "No protons detected!" << endl;
cout << endl;
cout << "b) Neutron " << endl;
if(FD_neutrid_default_PID_pass != 0  && FD_neutrid_charge_pass != 0){ 
cout << "number of particles which passed default neutron reference PID: "  	<< FD_neutrid_default_PID_pass << "   ( " << 100*FD_neutrid_default_PID_pass/FD_neutrid_charge_pass << " \% of neutral particles and " << 100*FD_neutrid_default_PID_pass/FD_neutrid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed beta cut: " 		<< FD_neutrid_beta_pass << "   ( " << 100*FD_neutrid_beta_pass/FD_neutrid_charge_pass << " \% of neutral particles and " << 100*FD_neutrid_beta_pass/FD_neutrid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed delta beta cut: " 		<< FD_neutrid_delta_beta_pass << "   ( " << 100*FD_neutrid_delta_beta_pass/FD_neutrid_charge_pass << " \% of neutral particles and " << 100*FD_neutrid_delta_beta_pass/FD_neutrid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed tofmass cut: " 		<< FD_neutrid_tofmass_pass << "   ( " << 100*FD_neutrid_tofmass_pass/FD_neutrid_charge_pass << " \% of neutral particles and " << 100*FD_neutrid_tofmass_pass/FD_neutrid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed delta vz cut: " 		<< FD_neutrid_delta_vz_pass << "   ( " << 100*FD_neutrid_delta_vz_pass/FD_neutrid_charge_pass << " \% of neutral particles and " << 100*FD_neutrid_delta_vz_pass/FD_neutrid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "__________________________________________________________________________________________________________________________________________________________" << endl;
cout << "number of particles which passed all neutron ID cuts: " 	   << FD_neutrid_all_pass << "   ( " << 100*FD_neutrid_all_pass/FD_neutrid_charge_pass << " \% of neutral particles and " << 100*FD_neutrid_all_pass/FD_neutrid_default_PID_pass << " \% of reference PID particles)" << endl;
}
else cout << "No neutrons detected!" << endl;
cout << endl;
cout << "b) Pip " << endl;
if(FD_pipid_default_PID_pass != 0  && FD_pipid_charge_pass != 0){
cout << "number of particles which passed default pip reference PID: "  	<< FD_pipid_default_PID_pass << "   ( " << 100*FD_pipid_default_PID_pass/FD_pipid_charge_pass << " \% of pos. particles and " << 100*FD_pipid_default_PID_pass/FD_pipid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed DC region 1 fid. cut: " 	    	<< FD_pipid_DC_hit_position_region1_fiducial_pass << "   ( " << 100*FD_pipid_DC_hit_position_region1_fiducial_pass/FD_pipid_charge_pass << " \% of pos. particles and " << 100*FD_pipid_DC_hit_position_region1_fiducial_pass/FD_pipid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed DC region 2 fid. cut: " 	    	<< FD_pipid_DC_hit_position_region2_fiducial_pass << "   ( " << 100*FD_pipid_DC_hit_position_region2_fiducial_pass/FD_pipid_charge_pass << " \% of pos. particles and " << 100*FD_pipid_DC_hit_position_region2_fiducial_pass/FD_pipid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed DC region 3 fid. cut: " 	    	<< FD_pipid_DC_hit_position_region3_fiducial_pass << "   ( " << 100*FD_pipid_DC_hit_position_region3_fiducial_pass/FD_pipid_charge_pass << " \% of pos. particles and " << 100*FD_pipid_DC_hit_position_region3_fiducial_pass/FD_pipid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed beta cut: " 		<< FD_pipid_beta_pass << "   ( " << 100*FD_pipid_beta_pass/FD_pipid_charge_pass << " \% of pos. particles and " << 100*FD_pipid_beta_pass/FD_pipid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed delta beta cut: " 		<< FD_pipid_delta_beta_pass << "   ( " << 100*FD_pipid_delta_beta_pass/FD_pipid_charge_pass << " \% of pos. particles and " << 100*FD_pipid_delta_beta_pass/FD_pipid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed tofmass cut: " 		<< FD_pipid_tofmass_pass << "   ( " << 100*FD_pipid_tofmass_pass/FD_pipid_charge_pass << " \% of pos. particles and " << 100*FD_pipid_tofmass_pass/FD_pipid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed maximum probability cut: " 		<< FD_pipid_maximum_probability_pass << "   ( " << 100*FD_pipid_maximum_probability_pass/FD_pipid_charge_pass << " \% of pos. particles and " << 100*FD_pipid_maximum_probability_pass/FD_pipid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed delta vz cut: " 		<< FD_pipid_delta_vz_pass << "   ( " << 100*FD_pipid_delta_vz_pass/FD_pipid_charge_pass << " \% of pos. particles and " << 100*FD_pipid_delta_vz_pass/FD_pipid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "__________________________________________________________________________________________________________________________________________________________" << endl;
cout << "number of particles which passed all positive pipon ID cuts: " 	   << FD_pipid_all_pass << "   ( " << 100*FD_pipid_all_pass/FD_pipid_charge_pass << " \% of pos. particles and " << 100*FD_pipid_all_pass/FD_pipid_default_PID_pass << " \% of reference PID particles)" << endl;
}
else cout << "No pi plus detected!" << endl;
cout << endl;
cout << "c) Pim " << endl;
if(FD_pimid_default_PID_pass != 0  && FD_pimid_charge_pass != 0){
cout << "number of particles which passed default pim reference PID: "  	<< FD_pimid_default_PID_pass << "   ( " << 100*FD_pimid_default_PID_pass/FD_pimid_charge_pass << " \% of neg. particles and " << 100*FD_pimid_default_PID_pass/FD_pimid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed DC region 1 fid. cut: " 	    	<< FD_pimid_DC_hit_position_region1_fiducial_pass << "   ( " << 100*FD_pimid_DC_hit_position_region1_fiducial_pass/FD_pimid_charge_pass << " \% of neg. particles and " << 100*FD_pimid_DC_hit_position_region1_fiducial_pass/FD_pimid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed DC region 2 fid. cut: " 	    	<< FD_pimid_DC_hit_position_region2_fiducial_pass << "   ( " << 100*FD_pimid_DC_hit_position_region2_fiducial_pass/FD_pimid_charge_pass << " \% of neg. particles and " << 100*FD_pimid_DC_hit_position_region2_fiducial_pass/FD_pimid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed DC region 3 fid. cut: " 	    	<< FD_pimid_DC_hit_position_region3_fiducial_pass << "   ( " << 100*FD_pimid_DC_hit_position_region3_fiducial_pass/FD_pimid_charge_pass << " \% of neg. particles and " << 100*FD_pimid_DC_hit_position_region3_fiducial_pass/FD_pimid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed beta cut: " 		<< FD_pimid_beta_pass << "   ( " << 100*FD_pimid_beta_pass/FD_pimid_charge_pass << " \% of neg. particles and " << 100*FD_pimid_beta_pass/FD_pimid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed delta beta cut: " 		<< FD_pimid_delta_beta_pass << "   ( " << 100*FD_pimid_delta_beta_pass/FD_pimid_charge_pass << " \% of neg. particles and " << 100*FD_pimid_delta_beta_pass/FD_pimid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed tofmass cut: " 		<< FD_pimid_tofmass_pass << "   ( " << 100*FD_pimid_tofmass_pass/FD_pimid_charge_pass << " \% of neg. particles and " << 100*FD_pimid_tofmass_pass/FD_pimid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed maximum probability cut: " 		<< FD_pimid_maximum_probability_pass << "   ( " << 100*FD_pimid_maximum_probability_pass/FD_pimid_charge_pass << " \% of pos. particles and " << 100*FD_pimid_maximum_probability_pass/FD_pimid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed delta vz cut: " 		<< FD_pimid_delta_vz_pass << "   ( " << 100*FD_pimid_delta_vz_pass/FD_pimid_charge_pass << " \% of neg. particles and " << 100*FD_pimid_delta_vz_pass/FD_pimid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "__________________________________________________________________________________________________________________________________________________________" << endl;
cout << "number of particles which passed all negative pion ID cuts: " 	   << FD_pimid_all_pass << "   ( " << 100*FD_pimid_all_pass/FD_pimid_charge_pass << " \% of neg. particles and " << 100*FD_pimid_all_pass/FD_pimid_default_PID_pass << " \% of reference PID particles)" << endl;
}
else cout << "No pi minus detected!" << endl;
cout << endl;
cout << "d) Kp " << endl;
if(FD_Kpid_default_PID_pass != 0  && FD_Kpid_charge_pass != 0){ 
cout << "number of particles which passed default Kp reference PID: "  	<< FD_Kpid_default_PID_pass << "   ( " << 100*FD_Kpid_default_PID_pass/FD_Kpid_charge_pass << " \% of pos. particles and " << 100*FD_Kpid_default_PID_pass/FD_Kpid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed DC region 1 fid. cut: " 	    	<< FD_Kpid_DC_hit_position_region1_fiducial_pass << "   ( " << 100*FD_Kpid_DC_hit_position_region1_fiducial_pass/FD_Kpid_charge_pass << " \% of pos. particles and " << 100*FD_Kpid_DC_hit_position_region1_fiducial_pass/FD_Kpid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed DC region 2 fid. cut: " 	    	<< FD_Kpid_DC_hit_position_region2_fiducial_pass << "   ( " << 100*FD_Kpid_DC_hit_position_region2_fiducial_pass/FD_Kpid_charge_pass << " \% of pos. particles and " << 100*FD_Kpid_DC_hit_position_region2_fiducial_pass/FD_Kpid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed DC region 3 fid. cut: " 	    	<< FD_Kpid_DC_hit_position_region3_fiducial_pass << "   ( " << 100*FD_Kpid_DC_hit_position_region3_fiducial_pass/FD_Kpid_charge_pass << " \% of pos. particles and " << 100*FD_Kpid_DC_hit_position_region3_fiducial_pass/FD_Kpid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed beta cut: " 		<< FD_Kpid_beta_pass << "   ( " << 100*FD_Kpid_beta_pass/FD_Kpid_charge_pass << " \% of pos. particles and " << 100*FD_Kpid_beta_pass/FD_Kpid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed delta beta cut: " 		<< FD_Kpid_delta_beta_pass << "   ( " << 100*FD_Kpid_delta_beta_pass/FD_Kpid_charge_pass << " \% of pos. particles and " << 100*FD_Kpid_delta_beta_pass/FD_Kpid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed tofmass cut: " 		<< FD_Kpid_tofmass_pass << "   ( " << 100*FD_Kpid_tofmass_pass/FD_Kpid_charge_pass << " \% of pos. particles and " << 100*FD_Kpid_tofmass_pass/FD_Kpid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed maximum probability cut: " 		<< FD_Kpid_maximum_probability_pass << "   ( " << 100*FD_Kpid_maximum_probability_pass/FD_Kpid_charge_pass << " \% of pos. particles and " << 100*FD_Kpid_maximum_probability_pass/FD_Kpid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed delta vz cut: " 		<< FD_Kpid_delta_vz_pass << "   ( " << 100*FD_Kpid_delta_vz_pass/FD_Kpid_charge_pass << " \% of pos. particles and " << 100*FD_Kpid_delta_vz_pass/FD_Kpid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "__________________________________________________________________________________________________________________________________________________________" << endl;
cout << "number of particles which passed all positive Kaon ID cuts: " 	   << FD_Kpid_all_pass << "   ( " << 100*FD_Kpid_all_pass/FD_Kpid_charge_pass << " \% of pos. particles and " << 100*FD_Kpid_all_pass/FD_Kpid_default_PID_pass << " \% of reference PID particles)" << endl;
}
else cout << "No K plus detected!" << endl;
cout << endl;
cout << "e) Km " << endl;
if(FD_Kmid_default_PID_pass != 0  && FD_Kmid_charge_pass != 0){ 
cout << "number of particles which passed default Km reference PID: "  	<< FD_Kmid_default_PID_pass << "   ( " << 100*FD_Kmid_default_PID_pass/FD_Kmid_charge_pass << " \% of neg. particles and " << 100*FD_Kmid_default_PID_pass/FD_Kmid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed DC region 1 fid. cut: " 	    	<< FD_Kmid_DC_hit_position_region1_fiducial_pass << "   ( " << 100*FD_Kmid_DC_hit_position_region1_fiducial_pass/FD_Kmid_charge_pass << " \% of neg. particles and " << 100*FD_Kmid_DC_hit_position_region1_fiducial_pass/FD_Kmid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed DC region 2 fid. cut: " 	    	<< FD_Kmid_DC_hit_position_region2_fiducial_pass << "   ( " << 100*FD_Kmid_DC_hit_position_region2_fiducial_pass/FD_Kmid_charge_pass << " \% of neg. particles and " << 100*FD_Kmid_DC_hit_position_region2_fiducial_pass/FD_Kmid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed DC region 3 fid. cut: " 	    	<< FD_Kmid_DC_hit_position_region3_fiducial_pass << "   ( " << 100*FD_Kmid_DC_hit_position_region3_fiducial_pass/FD_Kmid_charge_pass << " \% of neg. particles and " << 100*FD_Kmid_DC_hit_position_region3_fiducial_pass/FD_Kmid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed beta cut: " 		<< FD_Kmid_beta_pass << "   ( " << 100*FD_Kmid_beta_pass/FD_Kmid_charge_pass << " \% of neg. particles and " << 100*FD_Kmid_beta_pass/FD_Kmid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed delta beta cut: " 		<< FD_Kmid_delta_beta_pass << "   ( " << 100*FD_Kmid_delta_beta_pass/FD_Kmid_charge_pass << " \% of neg. particles and " << 100*FD_Kmid_delta_beta_pass/FD_Kmid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed tofmass cut: " 		<< FD_Kmid_tofmass_pass << "   ( " << 100*FD_Kmid_tofmass_pass/FD_Kmid_charge_pass << " \% of neg. particles and " << 100*FD_Kmid_tofmass_pass/FD_Kmid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed maximum probability cut: " 		<< FD_Kmid_maximum_probability_pass << "   ( " << 100*FD_Kmid_maximum_probability_pass/FD_Kmid_charge_pass << " \% of pos. particles and " << 100*FD_Kmid_maximum_probability_pass/FD_Kmid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed delta vz cut: " 		<< FD_Kmid_delta_vz_pass << "   ( " << 100*FD_Kmid_delta_vz_pass/FD_Kmid_charge_pass << " \% of neg. particles and " << 100*FD_Kmid_delta_vz_pass/FD_Kmid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "__________________________________________________________________________________________________________________________________________________________" << endl;
cout << "number of particles which passed all negative Kaon ID cuts: " 	   << FD_Kmid_all_pass << "   ( " << 100*FD_Kmid_all_pass/FD_Kmid_charge_pass << " \% of neg. particles and " << 100*FD_Kmid_all_pass/FD_Kmid_default_PID_pass << " \% of reference PID particles)" << endl;
}
else cout << "No K minus detected!" << endl;
cout << endl;
}

if(show_photon_ID_statistics){
cout << "------------------------------------------------------------------------------------------" << endl;
cout << "photon PID cut statistics " << endl;
cout << endl;
if(FD_photid_default_PID_pass != 0  && FD_photid_charge_pass != 0){ 
cout << "number of particles which passed default photon reference PID: "  	<< FD_photid_default_PID_pass << "   ( " << 100*FD_photid_default_PID_pass/FD_photid_charge_pass << " \% of neutral particles and " << 100*FD_photid_default_PID_pass/FD_photid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed beta cut: " 		<< FD_photid_beta_pass << "   ( " << 100*FD_photid_beta_pass/FD_photid_charge_pass << " \% of neutral particles and " << 100*FD_photid_beta_pass/FD_photid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "number of particles which passed EC hit position fiducial cut: " 		<< FD_photid_EC_hit_position_fiducial_pass << "   ( " << 100*FD_photid_EC_hit_position_fiducial_pass/FD_photid_charge_pass << " \% of neutral particles and " << 100*FD_photid_EC_hit_position_fiducial_pass/FD_photid_default_PID_pass << " \% of reference PID particles)" << endl;
cout << "__________________________________________________________________________________________________________________________________________________________" << endl;
cout << "number of particles which passed all photon cuts: " 	   << FD_photid_all_pass << "   ( " << 100*FD_photid_all_pass/FD_photid_charge_pass << " \% of neutral particles and " << 100*FD_photid_all_pass/FD_photid_default_PID_pass << " \% of reference PID particles)" << endl;
}
else cout << "No photons detected!" << endl;
cout << endl;
}


/// /////////////////////////////////////////////////////////////

    cout << endl;
    cout << "Tree successfully analysed!" << endl;
    cout << "Writing the output file ... " << endl;
    out->Write(); // Saving Histograms
    cout << "Histograms saved in File: " << outputfile << endl;
    out->Close(); // Closing Output File
    f->Close();   // Closing Input File
    cout << "... Completed!" << endl;

    return 1;


/// ///////////////////////////////////////////////////////////////////////////////
}   /// end of main
/// ///////////////////////////////////////////////////////////////////////////////


/// /////////////////////////////////////////////////////////////////////////////////
/// HTCC

TH2F* create_hist_HTCC_theta_vs_phi(const int cutnum){
  char name[100];
  sprintf(name,"HTCC_theta_vs_phi_cut_%02d", cutnum);
  hist_HTCC_theta_vs_phi[cutnum] = new TH2F(name, name, 360,-180,180, 80,0,40);   
  hist_HTCC_theta_vs_phi[cutnum]->GetXaxis()->SetTitle("#Phi_{CC}");
  hist_HTCC_theta_vs_phi[cutnum]->GetYaxis()->SetTitle("#theta_{CC}");
  return(hist_HTCC_theta_vs_phi[cutnum]);
}

TH1F* create_hist_HTCC_nphe(const int cutnum){
  char name[100];
  sprintf(name,"HTCC_nphe_cut_%02d", cutnum);
  hist_HTCC_nphe[cutnum] = new TH1F(name, name, 70, 0, 70);   
  hist_HTCC_nphe[cutnum]->GetXaxis()->SetTitle("nphe");
  hist_HTCC_nphe[cutnum]->GetYaxis()->SetTitle("counts");
  return(hist_HTCC_nphe[cutnum]);
}

TH2F* create_hist_HTCC_nphe_vs_sampling_fraction(const int cutnum){
  char name[100];
  sprintf(name,"HTCC_nphe_vs_sampfrac_cut_%02d", cutnum);
  hist_HTCC_nphe_vs_sampling_fraction[cutnum] = new TH2F(name, name, 70, 0, 70, 60, 0, 0.6);   
  hist_HTCC_nphe_vs_sampling_fraction[cutnum]->GetXaxis()->SetTitle("nphe");
  hist_HTCC_nphe_vs_sampling_fraction[cutnum]->GetYaxis()->SetTitle("EC sampling fraction");
  return(hist_HTCC_nphe_vs_sampling_fraction[cutnum]);
}


// ////////////////////////////////////////////////////////////////////////////////////
// EC cuts

TH2F* create_hist_EC_PCAL_vs_EC_ECAL(const int cutnum){
  char name[100];
  sprintf(name,"EC_PCAL_vs_EC_ECAL_cut_%02d", cutnum);
  hist_EC_PCAL_vs_EC_ECAL[cutnum] = new TH2F(name, name, 1000,0,1.0, 1000,0,1.0);   
  hist_EC_PCAL_vs_EC_ECAL[cutnum]->GetXaxis()->SetTitle("E_{PCAL} /GeV");
  hist_EC_PCAL_vs_EC_ECAL[cutnum]->GetYaxis()->SetTitle("E_{ECAL} /GeV");
  return(hist_EC_PCAL_vs_EC_ECAL[cutnum]);
}

TH2F* create_hist_EC_outer_vs_EC_inner(const int cutnum){
  char name[100];
  sprintf(name,"EC_outer_vs_EC_inner_cut_%02d", cutnum);
  hist_EC_outer_vs_EC_inner[cutnum] = new TH2F(name, name, 1000,0,1.0, 1000,0,1.0);   
  hist_EC_outer_vs_EC_inner[cutnum]->GetXaxis()->SetTitle("E_{EC inner} /GeV");
  hist_EC_outer_vs_EC_inner[cutnum]->GetYaxis()->SetTitle("E_{EC outer} /GeV");
  return(hist_EC_outer_vs_EC_inner[cutnum]);
}


///

TH2F* create_hist_EC_total_sampling_fraction_sec1(const int cutnum){
  char name[100];
  sprintf(name,"EC_total_sampling_fraction_sec1_cut_%02d", cutnum);
  hist_EC_total_sampling_fraction_sec1[cutnum] = new TH2F(name, name, 1100,0,11, 600,0,0.6);   
  hist_EC_total_sampling_fraction_sec1[cutnum]->GetXaxis()->SetTitle("momentum");
  hist_EC_total_sampling_fraction_sec1[cutnum]->GetYaxis()->SetTitle("sampling fraction");
  return(hist_EC_total_sampling_fraction_sec1[cutnum]);
}
TH2F* create_hist_EC_total_sampling_fraction_sec2(const int cutnum){
  char name[100];
  sprintf(name,"EC_total_sampling_fraction_sec2_cut_%02d", cutnum);
  hist_EC_total_sampling_fraction_sec2[cutnum] = new TH2F(name, name, 1100,0,11, 600,0,0.6);   
  hist_EC_total_sampling_fraction_sec2[cutnum]->GetXaxis()->SetTitle("momentum");
  hist_EC_total_sampling_fraction_sec2[cutnum]->GetYaxis()->SetTitle("sampling fraction");
  return(hist_EC_total_sampling_fraction_sec2[cutnum]);
}
TH2F* create_hist_EC_total_sampling_fraction_sec3(const int cutnum){
  char name[100];
  sprintf(name,"EC_total_sampling_fraction_sec3_cut_%02d", cutnum);
  hist_EC_total_sampling_fraction_sec3[cutnum] = new TH2F(name, name, 1100,0,11, 600,0,0.6);   
  hist_EC_total_sampling_fraction_sec3[cutnum]->GetXaxis()->SetTitle("momentum");
  hist_EC_total_sampling_fraction_sec3[cutnum]->GetYaxis()->SetTitle("sampling fraction");
  return(hist_EC_total_sampling_fraction_sec3[cutnum]);
}
TH2F* create_hist_EC_total_sampling_fraction_sec4(const int cutnum){
  char name[100];
  sprintf(name,"EC_total_sampling_fraction_sec4_cut_%02d", cutnum);
  hist_EC_total_sampling_fraction_sec4[cutnum] = new TH2F(name, name, 1100,0,11, 600,0,0.6);   
  hist_EC_total_sampling_fraction_sec4[cutnum]->GetXaxis()->SetTitle("momentum");
  hist_EC_total_sampling_fraction_sec4[cutnum]->GetYaxis()->SetTitle("sampling fraction");
  return(hist_EC_total_sampling_fraction_sec4[cutnum]);
}
TH2F* create_hist_EC_total_sampling_fraction_sec5(const int cutnum){
  char name[100];
  sprintf(name,"EC_total_sampling_fraction_sec5_cut_%02d", cutnum);
  hist_EC_total_sampling_fraction_sec5[cutnum] = new TH2F(name, name, 1100,0,11, 600,0,0.6);   
  hist_EC_total_sampling_fraction_sec5[cutnum]->GetXaxis()->SetTitle("momentum");
  hist_EC_total_sampling_fraction_sec5[cutnum]->GetYaxis()->SetTitle("sampling fraction");
  return(hist_EC_total_sampling_fraction_sec5[cutnum]);
}
TH2F* create_hist_EC_total_sampling_fraction_sec6(const int cutnum){
  char name[100];
  sprintf(name,"EC_total_sampling_fraction_sec6_cut_%02d", cutnum);
  hist_EC_total_sampling_fraction_sec6[cutnum] = new TH2F(name, name, 1100,0,11, 600,0,0.6);   
  hist_EC_total_sampling_fraction_sec6[cutnum]->GetXaxis()->SetTitle("momentum");
  hist_EC_total_sampling_fraction_sec6[cutnum]->GetYaxis()->SetTitle("sampling fraction");
  return(hist_EC_total_sampling_fraction_sec6[cutnum]);
}

TH2F* create_hist_EC_PCAL_sampling_fraction_sec1(const int cutnum){
  char name[100];
  sprintf(name,"EC_PCAL_sampling_fraction_sec1_cut_%02d", cutnum);
  hist_EC_PCAL_sampling_fraction_sec1[cutnum] = new TH2F(name, name, 1100,0,11, 60,0,0.6);   
  hist_EC_PCAL_sampling_fraction_sec1[cutnum]->GetXaxis()->SetTitle("momentum");
  hist_EC_PCAL_sampling_fraction_sec1[cutnum]->GetYaxis()->SetTitle("sampling fraction");
  return(hist_EC_PCAL_sampling_fraction_sec1[cutnum]);
}
TH2F* create_hist_EC_PCAL_sampling_fraction_sec2(const int cutnum){
  char name[100];
  sprintf(name,"EC_PCAL_sampling_fraction_sec2_cut_%02d", cutnum);
  hist_EC_PCAL_sampling_fraction_sec2[cutnum] = new TH2F(name, name, 1100,0,11, 60,0,0.6);   
  hist_EC_PCAL_sampling_fraction_sec2[cutnum]->GetXaxis()->SetTitle("momentum");
  hist_EC_PCAL_sampling_fraction_sec2[cutnum]->GetYaxis()->SetTitle("sampling fraction");
  return(hist_EC_PCAL_sampling_fraction_sec2[cutnum]);
}
TH2F* create_hist_EC_PCAL_sampling_fraction_sec3(const int cutnum){
  char name[100];
  sprintf(name,"EC_PCAL_sampling_fraction_sec3_cut_%02d", cutnum);
  hist_EC_PCAL_sampling_fraction_sec3[cutnum] = new TH2F(name, name, 1100,0,11, 60,0,0.6);   
  hist_EC_PCAL_sampling_fraction_sec3[cutnum]->GetXaxis()->SetTitle("momentum");
  hist_EC_PCAL_sampling_fraction_sec3[cutnum]->GetYaxis()->SetTitle("sampling fraction");
  return(hist_EC_PCAL_sampling_fraction_sec3[cutnum]);
}
TH2F* create_hist_EC_PCAL_sampling_fraction_sec4(const int cutnum){
  char name[100];
  sprintf(name,"EC_PCAL_sampling_fraction_sec4_cut_%02d", cutnum);
  hist_EC_PCAL_sampling_fraction_sec4[cutnum] = new TH2F(name, name, 1100,0,11, 60,0,0.6);   
  hist_EC_PCAL_sampling_fraction_sec4[cutnum]->GetXaxis()->SetTitle("momentum");
  hist_EC_PCAL_sampling_fraction_sec4[cutnum]->GetYaxis()->SetTitle("sampling fraction");
  return(hist_EC_PCAL_sampling_fraction_sec4[cutnum]);
}
TH2F* create_hist_EC_PCAL_sampling_fraction_sec5(const int cutnum){
  char name[100];
  sprintf(name,"EC_PCAL_sampling_fraction_sec5_cut_%02d", cutnum);
  hist_EC_PCAL_sampling_fraction_sec5[cutnum] = new TH2F(name, name, 1100,0,11, 60,0,0.6);   
  hist_EC_PCAL_sampling_fraction_sec5[cutnum]->GetXaxis()->SetTitle("momentum");
  hist_EC_PCAL_sampling_fraction_sec5[cutnum]->GetYaxis()->SetTitle("sampling fraction");
  return(hist_EC_PCAL_sampling_fraction_sec5[cutnum]);
}
TH2F* create_hist_EC_PCAL_sampling_fraction_sec6(const int cutnum){
  char name[100];
  sprintf(name,"EC_PCAL_sampling_fraction_sec6_cut_%02d", cutnum);
  hist_EC_PCAL_sampling_fraction_sec6[cutnum] = new TH2F(name, name, 1100,0,11, 60,0,0.6);   
  hist_EC_PCAL_sampling_fraction_sec6[cutnum]->GetXaxis()->SetTitle("momentum");
  hist_EC_PCAL_sampling_fraction_sec6[cutnum]->GetYaxis()->SetTitle("sampling fraction");
  return(hist_EC_PCAL_sampling_fraction_sec6[cutnum]);
}

TH2F* create_hist_EC_ECAL_sampling_fraction_sec1(const int cutnum){
  char name[100];
  sprintf(name,"EC_ECAL_sampling_fraction_sec1_cut_%02d", cutnum);
  hist_EC_ECAL_sampling_fraction_sec1[cutnum] = new TH2F(name, name, 1100,0,11, 40,0,0.4);   
  hist_EC_ECAL_sampling_fraction_sec1[cutnum]->GetXaxis()->SetTitle("momentum");
  hist_EC_ECAL_sampling_fraction_sec1[cutnum]->GetYaxis()->SetTitle("sampling fraction");
  return(hist_EC_ECAL_sampling_fraction_sec1[cutnum]);
}
TH2F* create_hist_EC_ECAL_sampling_fraction_sec2(const int cutnum){
  char name[100];
  sprintf(name,"EC_ECAL_sampling_fraction_sec2_cut_%02d", cutnum);
  hist_EC_ECAL_sampling_fraction_sec2[cutnum] = new TH2F(name, name, 1100,0,11, 40,0,0.4);   
  hist_EC_ECAL_sampling_fraction_sec2[cutnum]->GetXaxis()->SetTitle("momentum");
  hist_EC_ECAL_sampling_fraction_sec2[cutnum]->GetYaxis()->SetTitle("sampling fraction");
  return(hist_EC_ECAL_sampling_fraction_sec2[cutnum]);
}
TH2F* create_hist_EC_ECAL_sampling_fraction_sec3(const int cutnum){
  char name[100];
  sprintf(name,"EC_ECAL_sampling_fraction_sec3_cut_%02d", cutnum);
  hist_EC_ECAL_sampling_fraction_sec3[cutnum] = new TH2F(name, name, 1100,0,11, 40,0,0.4);   
  hist_EC_ECAL_sampling_fraction_sec3[cutnum]->GetXaxis()->SetTitle("momentum");
  hist_EC_ECAL_sampling_fraction_sec3[cutnum]->GetYaxis()->SetTitle("sampling fraction");
  return(hist_EC_ECAL_sampling_fraction_sec3[cutnum]);
}
TH2F* create_hist_EC_ECAL_sampling_fraction_sec4(const int cutnum){
  char name[100];
  sprintf(name,"EC_ECAL_sampling_fraction_sec4_cut_%02d", cutnum);
  hist_EC_ECAL_sampling_fraction_sec4[cutnum] = new TH2F(name, name, 1100,0,11, 40,0,0.4);   
  hist_EC_ECAL_sampling_fraction_sec4[cutnum]->GetXaxis()->SetTitle("momentum");
  hist_EC_ECAL_sampling_fraction_sec4[cutnum]->GetYaxis()->SetTitle("sampling fraction");
  return(hist_EC_ECAL_sampling_fraction_sec4[cutnum]);
}
TH2F* create_hist_EC_ECAL_sampling_fraction_sec5(const int cutnum){
  char name[100];
  sprintf(name,"EC_ECAL_sampling_fraction_sec5_cut_%02d", cutnum);
  hist_EC_ECAL_sampling_fraction_sec5[cutnum] = new TH2F(name, name, 1100,0,11, 40,0,0.4);   
  hist_EC_ECAL_sampling_fraction_sec5[cutnum]->GetXaxis()->SetTitle("momentum");
  hist_EC_ECAL_sampling_fraction_sec5[cutnum]->GetYaxis()->SetTitle("sampling fraction");
  return(hist_EC_ECAL_sampling_fraction_sec5[cutnum]);
}
TH2F* create_hist_EC_ECAL_sampling_fraction_sec6(const int cutnum){
  char name[100];
  sprintf(name,"EC_ECAL_sampling_fraction_sec6_cut_%02d", cutnum);
  hist_EC_ECAL_sampling_fraction_sec6[cutnum] = new TH2F(name, name, 1100,0,11, 40,0,0.4);   
  hist_EC_ECAL_sampling_fraction_sec6[cutnum]->GetXaxis()->SetTitle("momentum");
  hist_EC_ECAL_sampling_fraction_sec6[cutnum]->GetYaxis()->SetTitle("sampling fraction");
  return(hist_EC_ECAL_sampling_fraction_sec6[cutnum]);
}

///

TH2F* create_hist_EC_PCAL_hit_position(const int cutnum){
  char name[100];
  sprintf(name,"EC_PCAL_hit_position_cut_%02d", cutnum);
  hist_EC_PCAL_hit_position[cutnum] = new TH2F(name, name, 900,-450,450, 900,-450,450);   
  hist_EC_PCAL_hit_position[cutnum]->GetXaxis()->SetTitle("x /cm");
  hist_EC_PCAL_hit_position[cutnum]->GetYaxis()->SetTitle("y /cm");
  return(hist_EC_PCAL_hit_position[cutnum]);
}

TH2F* create_hist_EC_inner_hit_position(const int cutnum){
  char name[100];
  sprintf(name,"EC_inner_hit_position_cut_%02d", cutnum);
  hist_EC_inner_hit_position[cutnum] = new TH2F(name, name, 900,-450,450, 900,-450,450);   
  hist_EC_inner_hit_position[cutnum]->GetXaxis()->SetTitle("x /cm");
  hist_EC_inner_hit_position[cutnum]->GetYaxis()->SetTitle("y /cm");
  return(hist_EC_inner_hit_position[cutnum]);
}

TH2F* create_hist_EC_outer_hit_position(const int cutnum){
  char name[100];
  sprintf(name,"EC_outer_hit_position_cut_%02d", cutnum);
  hist_EC_outer_hit_position[cutnum] = new TH2F(name, name, 900,-450,450, 900,-450,450);   
  hist_EC_outer_hit_position[cutnum]->GetXaxis()->SetTitle("x /cm");
  hist_EC_outer_hit_position[cutnum]->GetYaxis()->SetTitle("y /cm");
  return(hist_EC_outer_hit_position[cutnum]);
}


// ////////////////////////////////////////////////////////////////////////////////////////
// DC cuts

TH2F *create_hist_DC_hit_position_region1(int cutnum){
  char name[100];
  sprintf(name,"DC_hit_position_region1_cut_%02d", cutnum);
  hist_DC_hit_position_region1[cutnum] = new TH2F(name, name, 900,-450,450, 900,-450,450);   
  hist_DC_hit_position_region1[cutnum]->GetXaxis()->SetTitle("x /cm");
  hist_DC_hit_position_region1[cutnum]->GetYaxis()->SetTitle("y /cm");
  return(hist_DC_hit_position_region1[cutnum]);
}

TH2F *create_hist_DC_hit_position_region2(int cutnum){
  char name[100];
  sprintf(name,"DC_hit_position_region2_cut_%02d", cutnum);
  hist_DC_hit_position_region2[cutnum] = new TH2F(name, name, 1000,-500,500, 1000,-500,500);   
  hist_DC_hit_position_region2[cutnum]->GetXaxis()->SetTitle("x /cm");
  hist_DC_hit_position_region2[cutnum]->GetYaxis()->SetTitle("y /cm");
  return(hist_DC_hit_position_region2[cutnum]);
}

TH2F *create_hist_DC_hit_position_region3(int cutnum){
  char name[100];
  sprintf(name,"DC_hit_position_region3_cut_%02d", cutnum);
  hist_DC_hit_position_region3[cutnum] = new TH2F(name, name, 1000,-500,500, 1000,-500,500);   
  hist_DC_hit_position_region3[cutnum]->GetXaxis()->SetTitle("x /cm");
  hist_DC_hit_position_region3[cutnum]->GetYaxis()->SetTitle("y /cm");
  return(hist_DC_hit_position_region3[cutnum]);
}


TH1F *create_hist_DC_z_vertex_sec1(int cutnum){
  char name[100];
  sprintf(name,"DC_z_vertex_sec1_cut_%02d", cutnum);
  hist_DC_z_vertex_sec1[cutnum] = new TH1F(name, name, 400,-40,40);   
  hist_DC_z_vertex_sec1[cutnum]->GetXaxis()->SetTitle("v_{z} /cm");
  hist_DC_z_vertex_sec1[cutnum]->GetYaxis()->SetTitle("counts");
  return(hist_DC_z_vertex_sec1[cutnum]);
}
TH1F *create_hist_DC_z_vertex_sec2(int cutnum){
  char name[100];
  sprintf(name,"DC_z_vertex_sec2_cut_%02d", cutnum);
  hist_DC_z_vertex_sec2[cutnum] = new TH1F(name, name, 400,-40,40);   
  hist_DC_z_vertex_sec2[cutnum]->GetXaxis()->SetTitle("v_{z} /cm");
  hist_DC_z_vertex_sec2[cutnum]->GetYaxis()->SetTitle("counts");
  return(hist_DC_z_vertex_sec2[cutnum]);
}
TH1F *create_hist_DC_z_vertex_sec3(int cutnum){
  char name[100];
  sprintf(name,"DC_z_vertex_sec3_cut_%02d", cutnum);
  hist_DC_z_vertex_sec3[cutnum] = new TH1F(name, name, 400,-40,40);   
  hist_DC_z_vertex_sec3[cutnum]->GetXaxis()->SetTitle("v_{z} /cm");
  hist_DC_z_vertex_sec3[cutnum]->GetYaxis()->SetTitle("counts");
  return(hist_DC_z_vertex_sec3[cutnum]);
}
TH1F *create_hist_DC_z_vertex_sec4(int cutnum){
  char name[100];
  sprintf(name,"DC_z_vertex_sec4_cut_%02d", cutnum);
  hist_DC_z_vertex_sec4[cutnum] = new TH1F(name, name, 400,-40,40);   
  hist_DC_z_vertex_sec4[cutnum]->GetXaxis()->SetTitle("v_{z} /cm");
  hist_DC_z_vertex_sec4[cutnum]->GetYaxis()->SetTitle("counts");
  return(hist_DC_z_vertex_sec4[cutnum]);
}
TH1F *create_hist_DC_z_vertex_sec5(int cutnum){
  char name[100];
  sprintf(name,"DC_z_vertex_sec5_cut_%02d", cutnum);
  hist_DC_z_vertex_sec5[cutnum] = new TH1F(name, name, 400,-40,40);   
  hist_DC_z_vertex_sec5[cutnum]->GetXaxis()->SetTitle("v_{z} /cm");
  hist_DC_z_vertex_sec5[cutnum]->GetYaxis()->SetTitle("counts");
  return(hist_DC_z_vertex_sec5[cutnum]);
}
TH1F *create_hist_DC_z_vertex_sec6(int cutnum){
  char name[100];
  sprintf(name,"DC_z_vertex_sec6_cut_%02d", cutnum);
  hist_DC_z_vertex_sec6[cutnum] = new TH1F(name, name, 400,-40,40);   
  hist_DC_z_vertex_sec6[cutnum]->GetXaxis()->SetTitle("v_{z} /cm");
  hist_DC_z_vertex_sec6[cutnum]->GetYaxis()->SetTitle("counts");
  return(hist_DC_z_vertex_sec6[cutnum]);
}



// /////////////////////////////////////////////////////////////////////////////
// Basic hadron cuts

TH2F *create_DC_hit_position_region1_hadron(int cutnum){
  char name[100];
  sprintf(name,"DC_hit_position_region1_hadron_cut_%02d", cutnum);
  hist_DC_hit_position_region1_hadron[cutnum] = new TH2F(name, name, 900,-450,450, 900,-450,450);   
  hist_DC_hit_position_region1_hadron[cutnum]->GetXaxis()->SetTitle("x /cm");
  hist_DC_hit_position_region1_hadron[cutnum]->GetYaxis()->SetTitle("y /cm");
  return(hist_DC_hit_position_region1_hadron[cutnum]);
}

TH2F *create_DC_hit_position_region2_hadron(int cutnum){
  char name[100];
  sprintf(name,"DC_hit_position_region2_hadron_cut_%02d", cutnum);
  hist_DC_hit_position_region2_hadron[cutnum] = new TH2F(name, name, 900,-450,450, 900,-450,450);   
  hist_DC_hit_position_region2_hadron[cutnum]->GetXaxis()->SetTitle("x /cm");
  hist_DC_hit_position_region2_hadron[cutnum]->GetYaxis()->SetTitle("y /cm");
  return(hist_DC_hit_position_region2_hadron[cutnum]);
}

TH2F *create_DC_hit_position_region3_hadron(int cutnum){
  char name[100];
  sprintf(name,"DC_hit_position_region3_hadron_cut_%02d", cutnum);
  hist_DC_hit_position_region3_hadron[cutnum] = new TH2F(name, name, 900,-450,450, 900,-450,450);   
  hist_DC_hit_position_region3_hadron[cutnum]->GetXaxis()->SetTitle("x /cm");
  hist_DC_hit_position_region3_hadron[cutnum]->GetYaxis()->SetTitle("y /cm");
  return(hist_DC_hit_position_region3_hadron[cutnum]);
}


TH2F* create_hist_EC_outer_vs_EC_inner_hadron(const int cutnum){
  char name[100];
  sprintf(name,"EC_outer_vs_EC_inner_hadron_cut_%02d", cutnum);
  hist_EC_outer_vs_EC_inner_hadron[cutnum] = new TH2F(name, name, 1000,0,1.0, 1000,0,1.0);   
  hist_EC_outer_vs_EC_inner_hadron[cutnum]->GetXaxis()->SetTitle("E_{EC inner} /GeV");
  hist_EC_outer_vs_EC_inner_hadron[cutnum]->GetYaxis()->SetTitle("E_{EC outer} /GeV");
  return(hist_EC_outer_vs_EC_inner_hadron[cutnum]);
}

// /////////////////////////////////////////////////////////////////////////////
// TOF beta cuts

TH2F *create_hist_beta_vs_p(int cutnum){
  char name[100];
  sprintf(name,"beta_vs_momentum_cut_%02d", cutnum);
  hist_beta_vs_p[cutnum] = new TH2F(name, name, 1100,0,11, 1400,0,1.40);   
  hist_beta_vs_p[cutnum]->GetXaxis()->SetTitle("p /GeV");
  hist_beta_vs_p[cutnum]->GetYaxis()->SetTitle("#beta");
  return(hist_beta_vs_p[cutnum]);
}

TH2F *create_hist_beta_vs_p_sec1(int cutnum){
  char name[100];
  sprintf(name,"beta_vs_momentum_sec1_cut_%02d", cutnum);
  hist_beta_vs_p_sec1[cutnum] = new TH2F(name, name, 1100,0,11, 1400,0,1.40);   
  hist_beta_vs_p_sec1[cutnum]->GetXaxis()->SetTitle("p /GeV");
  hist_beta_vs_p_sec1[cutnum]->GetYaxis()->SetTitle("#beta");
  return(hist_beta_vs_p_sec1[cutnum]);
}
TH2F *create_hist_beta_vs_p_sec2(int cutnum){
  char name[100];
  sprintf(name,"beta_vs_momentum_sec2_cut_%02d", cutnum);
  hist_beta_vs_p_sec2[cutnum] = new TH2F(name, name, 1100,0,11, 1400,0,1.40);   
  hist_beta_vs_p_sec2[cutnum]->GetXaxis()->SetTitle("p /GeV");
  hist_beta_vs_p_sec2[cutnum]->GetYaxis()->SetTitle("#beta");
  return(hist_beta_vs_p_sec2[cutnum]);
}
TH2F *create_hist_beta_vs_p_sec3(int cutnum){
  char name[100];
  sprintf(name,"beta_vs_momentum_sec3_cut_%02d", cutnum);
  hist_beta_vs_p_sec3[cutnum] = new TH2F(name, name, 1100,0,11, 1400,0,1.40);   
  hist_beta_vs_p_sec3[cutnum]->GetXaxis()->SetTitle("p /GeV");
  hist_beta_vs_p_sec3[cutnum]->GetYaxis()->SetTitle("#beta");
  return(hist_beta_vs_p_sec3[cutnum]);
}
TH2F *create_hist_beta_vs_p_sec4(int cutnum){
  char name[100];
  sprintf(name,"beta_vs_momentum_sec4_cut_%02d", cutnum);
  hist_beta_vs_p_sec4[cutnum] = new TH2F(name, name, 1100,0,11, 1400,0,1.40);   
  hist_beta_vs_p_sec4[cutnum]->GetXaxis()->SetTitle("p /GeV");
  hist_beta_vs_p_sec4[cutnum]->GetYaxis()->SetTitle("#beta");
  return(hist_beta_vs_p_sec4[cutnum]);
}
TH2F *create_hist_beta_vs_p_sec5(int cutnum){
  char name[100];
  sprintf(name,"beta_vs_momentum_sec5_cut_%02d", cutnum);
  hist_beta_vs_p_sec5[cutnum] = new TH2F(name, name, 1100,0,11, 1400,0,1.40);   
  hist_beta_vs_p_sec5[cutnum]->GetXaxis()->SetTitle("p /GeV");
  hist_beta_vs_p_sec5[cutnum]->GetYaxis()->SetTitle("#beta");
  return(hist_beta_vs_p_sec5[cutnum]);
}
TH2F *create_hist_beta_vs_p_sec6(int cutnum){
  char name[100];
  sprintf(name,"beta_vs_momentum_sec6_cut_%02d", cutnum);
  hist_beta_vs_p_sec6[cutnum] = new TH2F(name, name, 1100,0,11, 1400,0,1.40);   
  hist_beta_vs_p_sec6[cutnum]->GetXaxis()->SetTitle("p /GeV");
  hist_beta_vs_p_sec6[cutnum]->GetYaxis()->SetTitle("#beta");
  return(hist_beta_vs_p_sec6[cutnum]);
}

TH2F *create_hist_delta_beta_vs_p(int cutnum){
  char name[100];
  sprintf(name,"delta_beta_vs_momentum_cut_%02d", cutnum);
  hist_delta_beta_vs_p[cutnum] = new TH2F(name, name, 1100,0,11, 250,-0.5,0.5);   
  hist_delta_beta_vs_p[cutnum]->GetXaxis()->SetTitle("p /GeV");
  hist_delta_beta_vs_p[cutnum]->GetYaxis()->SetTitle("#Delta #beta");
  return(hist_delta_beta_vs_p[cutnum]);
}

TH1F *create_hist_delta_beta(int cutnum){
  char name[100];
  sprintf(name,"delta_beta_cut_%02d", cutnum);
  hist_delta_beta[cutnum] = new TH1F(name, name, 500,-0.5,0.5);   
  hist_delta_beta[cutnum]->GetXaxis()->SetTitle("#Delta #beta");
  hist_delta_beta[cutnum]->GetYaxis()->SetTitle("counts");
  return(hist_delta_beta[cutnum]);
}

TH2F *create_hist_tofmass_vs_p(int cutnum){
  char name[100];
  sprintf(name,"tofmass2_vs_momentum_cut_%02d", cutnum);
  hist_tofmass_vs_p[cutnum] = new TH2F(name, name, 1100,0,11, 1100, -0.2, 2.0);   
  hist_tofmass_vs_p[cutnum]->GetXaxis()->SetTitle("p /GeV");
  hist_tofmass_vs_p[cutnum]->GetYaxis()->SetTitle("m_{TOF}^{2} /GeV^{2}");
  return(hist_tofmass_vs_p[cutnum]);
}

TH1F *create_hist_tofmass(int cutnum){
  char name[100];
  sprintf(name,"tofmass2_cut_%02d", cutnum);
  hist_tofmass[cutnum] = new TH1F(name, name, 1100, -0.2, 2.0);   
  hist_tofmass[cutnum]->GetXaxis()->SetTitle("m_{TOF}^{2} /GeV^{2}");
  hist_tofmass[cutnum]->GetYaxis()->SetTitle("counts");
  return(hist_tofmass[cutnum]);
}

TH1F *create_hist_delta_vz(int cutnum){
  char name[100];
  sprintf(name,"delta_vz_cut_%02d", cutnum);
  hist_delta_vz[cutnum] = new TH1F(name, name, 200, -100, 100);   
  hist_delta_vz[cutnum]->GetXaxis()->SetTitle("#Delta v_{z} /cm");
  hist_delta_vz[cutnum]->GetYaxis()->SetTitle("counts");
  return(hist_delta_vz[cutnum]);
}

// /////////////////////////////////////////////////////////////////////////////
// CTOF beta plots

TH2F *create_hist_CD_beta_vs_p(int cutnum){
  char name[100];
  sprintf(name,"CD_beta_vs_momentum_cut_%02d", cutnum);
  hist_CD_beta_vs_p[cutnum] = new TH2F(name, name, 1100,0,11, 1400,0,1.40);   
  hist_CD_beta_vs_p[cutnum]->GetXaxis()->SetTitle("p /GeV");
  hist_CD_beta_vs_p[cutnum]->GetYaxis()->SetTitle("#beta");
  return(hist_CD_beta_vs_p[cutnum]);
}


TH1F *create_hist_CD_delta_vz(int cutnum){
  char name[100];
  sprintf(name,"CD_delta_vz_cut_%02d", cutnum);
  hist_CD_delta_vz[cutnum] = new TH1F(name, name, 200, -100, 100);   
  hist_CD_delta_vz[cutnum]->GetXaxis()->SetTitle("#Delta v_{z} /cm");
  hist_CD_delta_vz[cutnum]->GetYaxis()->SetTitle("counts");
  return(hist_CD_delta_vz[cutnum]);
}


// /////////////////////////////////////////////////////////////////////////////
// photon ID plots based on ECAL

TH2F *create_hist_beta_vs_p_phot(int cutnum){
  char name[100];
  sprintf(name,"beta_vs_p_phot_cut_%02d", cutnum);
  hist_beta_vs_p_phot[cutnum] = new TH2F(name, name, 300,0,6, 700,0,1.40);   
  hist_beta_vs_p_phot[cutnum]->GetXaxis()->SetTitle("p / GeV");
  hist_beta_vs_p_phot[cutnum]->GetYaxis()->SetTitle("#beta");
  return(hist_beta_vs_p_phot[cutnum]);
}

TH1F *create_hist_beta_phot(int cutnum){
  char name[100];
  sprintf(name,"beta_phot_cut_%02d", cutnum);
  hist_beta_phot[cutnum] = new TH1F(name, name, 750, 0, 1.5);   
  hist_beta_phot[cutnum]->GetXaxis()->SetTitle("#beta");
  hist_beta_phot[cutnum]->GetYaxis()->SetTitle("counts");
  return(hist_beta_phot[cutnum]);
}

TH2F *create_hist_EC_sampling_fraction_phot(int cutnum){
  char name[100];
  sprintf(name,"EC_sampling_fraction_phot_cut_%02d", cutnum);
  hist_EC_sampling_fraction_phot[cutnum] = new TH2F(name, name, 600,0,6, 100,0,1);   
  hist_EC_sampling_fraction_phot[cutnum]->GetXaxis()->SetTitle("momentum");
  hist_EC_sampling_fraction_phot[cutnum]->GetYaxis()->SetTitle("sampling fraction");
  return(hist_EC_sampling_fraction_phot[cutnum]);
}

TH2F* create_hist_EC_PCAL_vs_EC_ECAL_phot(const int cutnum){
  char name[100];
  sprintf(name,"EC_PCAL_vs_EC_ECAL_phot_cut_%02d", cutnum);
  hist_EC_PCAL_vs_EC_ECAL_phot[cutnum] = new TH2F(name, name, 400,0,0.4, 400,0,0.4);   
  hist_EC_PCAL_vs_EC_ECAL_phot[cutnum]->GetXaxis()->SetTitle("E_{PCAL} /GeV");
  hist_EC_PCAL_vs_EC_ECAL_phot[cutnum]->GetYaxis()->SetTitle("E_{ECAL} /GeV");
  return(hist_EC_PCAL_vs_EC_ECAL_phot[cutnum]);
}

TH2F *create_hist_EC_PCAL_hit_position_phot(int cutnum){
  char name[100];
  sprintf(name,"EC_PCAL_hit_position_phot_cut_%02d", cutnum);
  hist_EC_PCAL_hit_position_phot[cutnum] = new TH2F(name, name, 500, -500, 500, 500, -500, 500);   
  hist_EC_PCAL_hit_position_phot[cutnum]->GetXaxis()->SetTitle("PCAL x /cm");
  hist_EC_PCAL_hit_position_phot[cutnum]->GetYaxis()->SetTitle("PCAL y /cm");
  return(hist_EC_PCAL_hit_position_phot[cutnum]);
}

/// ////////////////////////////////////////////////////////////////////////
/// FT PID plots

TH2F *create_hist_FT_FTCAL_energy_vs_radius(int cutnum){
  char name[100];
  sprintf(name,"FT_FTCAL_energy_vs_radius_cut_%02d", cutnum);
  hist_FT_FTCAL_energy_vs_radius[cutnum] = new TH2F(name, name, 600, 0, 6, 100, 0, 10);   
  hist_FT_FTCAL_energy_vs_radius[cutnum]->GetXaxis()->SetTitle("FTCAL energy /GeV");
  hist_FT_FTCAL_energy_vs_radius[cutnum]->GetYaxis()->SetTitle("FTCAL cluster radius /cm");
  return(hist_FT_FTCAL_energy_vs_radius[cutnum]);
}

TH2F *create_hist_FT_FTCAL_hit_position(int cutnum){
  char name[100];
  sprintf(name,"FT_FTCAL_hit_position_cut_%02d", cutnum);
  hist_FT_FTCAL_hit_position[cutnum] = new TH2F(name, name, 40, -20, 20, 40, -20, 20);   
  hist_FT_FTCAL_hit_position[cutnum]->GetXaxis()->SetTitle("FTCAL x /cm");
  hist_FT_FTCAL_hit_position[cutnum]->GetYaxis()->SetTitle("FTCAL y /cm");
  return(hist_FT_FTCAL_hit_position[cutnum]);
}

TH2F *create_hist_FT_FTTRK_hit_position(int cutnum){
  char name[100];
  sprintf(name,"FT_FTTRK_hit_position_cut_%02d", cutnum);
  hist_FT_FTTRK_hit_position[cutnum] = new TH2F(name, name, 40, -20, 20, 40, -20, 20);   
  hist_FT_FTTRK_hit_position[cutnum]->GetXaxis()->SetTitle("FTTRK x /cm");
  hist_FT_FTTRK_hit_position[cutnum]->GetYaxis()->SetTitle("FTTRK y /cm");
  return(hist_FT_FTTRK_hit_position[cutnum]);
}


TH2F *create_hist_FT_FTHODO_hit_position(int cutnum){
  char name[100];
  sprintf(name,"FT_FTHODO_hit_position_cut_%02d", cutnum);
  hist_FT_FTHODO_hit_position[cutnum] = new TH2F(name, name, 40, -20, 20, 40, -20, 20);   
  hist_FT_FTHODO_hit_position[cutnum]->GetXaxis()->SetTitle("FTHODO x /cm");
  hist_FT_FTHODO_hit_position[cutnum]->GetYaxis()->SetTitle("FTHODO y /cm");
  return(hist_FT_FTHODO_hit_position[cutnum]);
}


TH1F *create_hist_FT_beta(int cutnum){
  char name[100];
  sprintf(name,"FT_beta_cut_%02d", cutnum);
  hist_FT_beta[cutnum] = new TH1F(name, name, 130, 0, 1.3);   
  hist_FT_beta[cutnum]->GetXaxis()->SetTitle("#beta");
  hist_FT_beta[cutnum]->GetYaxis()->SetTitle("counts");
  return(hist_FT_beta[cutnum]);
}


/// ///////////////////////////////////////////////////////////////////////////////////////////////////////
/// assign raw particles

void get_event_properties(void){

  if(vNRUN->size() > 0)    NRUN = vNRUN->at(0);
  if(vNEVENT->size() > 0)  NEVENT = vNEVENT->at(0);
  if(vEVNTime->size() > 0) EVNTime = vNEVENT->at(0);
  if(vTYPE->size() > 0)    TYPE = vTYPE->at(0);
  if(vTRG->size() > 0)     TRG = vNEVENT->at(0);
  if(vBCG->size() > 0)     BCG = vNEVENT->at(0);
  if(vSTTime->size() > 0)  STTime = vSTTime->at(0);
  if(vRFTime->size() > 0)  RFTime = vRFTime->at(0);
  if(vHelic->size() > 0)   Helic = vHelic->at(0);

}


/// ///////////////////////////////////////////////////////////////////////////////////////////////////////
/// assign raw particles and detector properties

void assign_particles(void){

  Npart = vpart_px->size();

  for(Int_t i = 0; i < Npart; i++){
    if( i < BUFFER){
      part_px[i] = vpart_px->at(i);
      part_py[i] = vpart_py->at(i);
      part_pz[i] = vpart_pz->at(i);
      part_vx[i] = vpart_vx->at(i);
      part_vy[i] = vpart_vy->at(i);
      part_vz[i] = vpart_vz->at(i);
      part_charge[i] = vpart_charge->at(i);
      part_beta[i] = vpart_beta->at(i);
      part_pid[i] = vpart_pid->at(i);
      part_status[i] = vpart_status->at(i);
      part_p[i] = sqrt(part_px[i]*part_px[i] + part_py[i]*part_py[i] + part_pz[i]*part_pz[i]);
      part_theta[i] = acos(part_pz[i]/part_p[i]);
      part_phi[i] = atan2(part_py[i],part_px[i]);
    }
  } 

  int Cal_Nentries = 0; int CC_Nentries = 0; int FT_Nentries = 0;
  int SC_Nentries = 0; int Traj_Nentries = 0;
  int TRK_Nentries = 0; int TBT_Nentries = 0;

  Cal_Nentries  = vCal_pindex->size();
  CC_Nentries   = vCC_pindex->size();
  FT_Nentries   = vFT_pindex->size();
  SC_Nentries   = vSC_pindex->size();
  TRK_Nentries  = vTRK_pindex->size();
  Traj_Nentries = vTraj_pindex->size();


  // Calorimeter bank (detector = 7  ---  layer:  PCAL = 1, ECin = 4, ECout = 7)

  if(Cal_Nentries > 0){
    for(int i = 0; i < Cal_Nentries; i++){
      if(vCal_pindex->at(i) >= 0 && vCal_pindex->at(i) < BUFFER && vCal_layer->at(i) == 1){
        part_Cal_PCAL_sector[vCal_pindex->at(i)] = vCal_sector->at(i); 
        part_Cal_PCAL_energy[vCal_pindex->at(i)] = vCal_energy->at(i); 
        part_Cal_PCAL_time[vCal_pindex->at(i)] = vCal_time->at(i);
        part_Cal_PCAL_path[vCal_pindex->at(i)] = vCal_path->at(i);
        part_Cal_PCAL_x[vCal_pindex->at(i)] = vCal_x->at(i);
        part_Cal_PCAL_y[vCal_pindex->at(i)] = vCal_y->at(i);
        part_Cal_PCAL_z[vCal_pindex->at(i)] = vCal_z->at(i);
        part_Cal_PCAL_lu[vCal_pindex->at(i)] = vCal_lu->at(i);
        part_Cal_PCAL_lv[vCal_pindex->at(i)] = vCal_lv->at(i);
        part_Cal_PCAL_lw[vCal_pindex->at(i)] = vCal_lw->at(i);
      }
      if(vCal_pindex->at(i) >= 0 && vCal_pindex->at(i) < BUFFER && vCal_layer->at(i) == 4){
        part_Cal_ECin_sector[vCal_pindex->at(i)] = vCal_sector->at(i); 
        part_Cal_ECin_energy[vCal_pindex->at(i)] = vCal_energy->at(i); 
        part_Cal_ECin_time[vCal_pindex->at(i)] = vCal_time->at(i); 
        part_Cal_ECin_path[vCal_pindex->at(i)] = vCal_path->at(i); 
        part_Cal_ECin_x[vCal_pindex->at(i)] = vCal_x->at(i); 
        part_Cal_ECin_y[vCal_pindex->at(i)] = vCal_y->at(i); 
        part_Cal_ECin_z[vCal_pindex->at(i)] = vCal_z->at(i); 
        part_Cal_ECin_lu[vCal_pindex->at(i)] = vCal_lu->at(i); 
        part_Cal_ECin_lv[vCal_pindex->at(i)] = vCal_lv->at(i); 
        part_Cal_ECin_lw[vCal_pindex->at(i)] = vCal_lw->at(i); 
      }
      if(vCal_pindex->at(i) >= 0 && vCal_pindex->at(i) < BUFFER && vCal_layer->at(i) == 7){
        part_Cal_ECout_sector[vCal_pindex->at(i)] = vCal_sector->at(i);  
        part_Cal_ECout_energy[vCal_pindex->at(i)] = vCal_energy->at(i); 
        part_Cal_ECout_time[vCal_pindex->at(i)] = vCal_time->at(i); 
        part_Cal_ECout_path[vCal_pindex->at(i)] = vCal_path->at(i); 
        part_Cal_ECout_x[vCal_pindex->at(i)] = vCal_x->at(i); 
        part_Cal_ECout_y[vCal_pindex->at(i)] = vCal_y->at(i); 
        part_Cal_ECout_z[vCal_pindex->at(i)] = vCal_z->at(i); 
        part_Cal_ECout_lu[vCal_pindex->at(i)] = vCal_lu->at(i); 
        part_Cal_ECout_lv[vCal_pindex->at(i)] = vCal_lv->at(i); 
        part_Cal_ECout_lw[vCal_pindex->at(i)] = vCal_lw->at(i); 
      }
    }
  }

  for(Int_t i = 0; i < BUFFER; i++){
    part_Cal_energy_total[i] = part_Cal_PCAL_energy[i] + part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i];
  }

  // Cherenkov bank (detectors:  HTCC = 15,  LTCC = 16)

  if(CC_Nentries > 0){
    for(int i = 0; i < CC_Nentries; i++){
      if(vCC_pindex->at(i) >= 0 && vCC_pindex->at(i) < BUFFER && vCC_detector->at(i) == 15){ 
        part_CC_HTCC_sector[vCC_pindex->at(i)] = vCC_sector->at(i); 
        part_CC_HTCC_nphe[vCC_pindex->at(i)] = vCC_nphe->at(i);
        part_CC_HTCC_time[vCC_pindex->at(i)] = vCC_time->at(i); 
        part_CC_HTCC_path[vCC_pindex->at(i)] = vCC_path->at(i); 
        part_CC_HTCC_theta[vCC_pindex->at(i)] = vCC_theta->at(i); 
        part_CC_HTCC_phi[vCC_pindex->at(i)] = vCC_phi->at(i);
      }
      if(vCC_pindex->at(i) >= 0 && vCC_pindex->at(i) < BUFFER && vCC_detector->at(i) == 16){
        part_CC_LTCC_sector[vCC_pindex->at(i)] = vCC_sector->at(i); 
        part_CC_LTCC_nphe[vCC_pindex->at(i)] = vCC_nphe->at(i);
        part_CC_LTCC_time[vCC_pindex->at(i)] = vCC_time->at(i); 
        part_CC_LTCC_path[vCC_pindex->at(i)] = vCC_path->at(i);   
        part_CC_LTCC_theta[vCC_pindex->at(i)] = vCC_theta->at(i); 
        part_CC_LTCC_phi[vCC_pindex->at(i)] = vCC_phi->at(i);
      }
    }
  }

  // Forward Tagger bank (detectors: FTCAL = 10, FTTRK = 13, FTHODO = 11)

  if(FT_Nentries > 0){
    for(int i = 0; i < FT_Nentries; i++){
      if(vFT_pindex->at(i) >= 0 && vFT_pindex->at(i) < BUFFER){
        if(vFT_detector->at(i) == 10) part_FT_energy[vFT_pindex->at(i)] = vFT_energy->at(i); 
        if(vFT_detector->at(i) == 10) part_FT_radius[vFT_pindex->at(i)] = vFT_radius->at(i);
        if(vFT_detector->at(i) == 11) part_FTHODO_time[vFT_pindex->at(i)] = vFT_time->at(i); 
        if(vFT_detector->at(i) == 11) part_FTHODO_path[vFT_pindex->at(i)] = vFT_path->at(i); 
        if(vFT_detector->at(i) == 10) part_FTCAL_time[vFT_pindex->at(i)] = vFT_time->at(i); 
        if(vFT_detector->at(i) == 10) part_FTCAL_path[vFT_pindex->at(i)] = vFT_path->at(i); 
        if(vFT_detector->at(i) == 10) part_FTCAL_x[vFT_pindex->at(i)] = vFT_x->at(i); 
        if(vFT_detector->at(i) == 10) part_FTCAL_y[vFT_pindex->at(i)] = vFT_y->at(i); 
        if(vFT_detector->at(i) == 10) part_FTCAL_z[vFT_pindex->at(i)] = vFT_z->at(i); 
        if(vFT_detector->at(i) == 11) part_FTHODO_x[vFT_pindex->at(i)] = vFT_x->at(i); 
        if(vFT_detector->at(i) == 11) part_FTHODO_y[vFT_pindex->at(i)] = vFT_y->at(i); 
        if(vFT_detector->at(i) == 11) part_FTHODO_z[vFT_pindex->at(i)] = vFT_z->at(i); 
        if(vFT_detector->at(i) == 13) part_FTTRK_x[vFT_pindex->at(i)] = vFT_x->at(i); 
        if(vFT_detector->at(i) == 13) part_FTTRK_y[vFT_pindex->at(i)] = vFT_y->at(i); 
        if(vFT_detector->at(i) == 13) part_FTTRK_z[vFT_pindex->at(i)] = vFT_z->at(i); 
      }
    }
  }

  // Scintillator bank (detectors:  CTOF = 4,  CND = 3,  FTOF = 12  ---  for FTOF:  layer = 1,2)

  if(SC_Nentries > 0){
    for(int i = 0; i < SC_Nentries; i++){
      // FTOF
      if(vSC_pindex->at(i) >= 0 && vSC_pindex->at(i) < BUFFER && vSC_detector->at(i) == 12){                  
        if(vSC_layer->at(i) == 1){
          part_FTOF_layer[vSC_pindex->at(i)] = 1;
          part_FTOF_energy_layer1[vSC_pindex->at(i)] = vSC_energy->at(i); 
          part_FTOF_time_layer1[vSC_pindex->at(i)] = vSC_time->at(i); 
          part_FTOF_path_layer1[vSC_pindex->at(i)] = vSC_path->at(i);  
          part_FTOF_sector_layer1[vSC_pindex->at(i)] = vSC_sector->at(i); 
          part_FTOF_component_layer1[vSC_pindex->at(i)] = vSC_component->at(i);
        } 
        if(vSC_layer->at(i) == 2){
          part_FTOF_layer[vSC_pindex->at(i)] = 2;
          part_FTOF_sector_layer2[vSC_pindex->at(i)] = vSC_sector->at(i); 
          part_FTOF_component_layer2[vSC_pindex->at(i)] = vSC_component->at(i); 
        }
        if(vSC_layer->at(i) == 3){
          part_FTOF_layer[vSC_pindex->at(i)] = 3;
          part_FTOF_energy_layer3[vSC_pindex->at(i)] = vSC_energy->at(i); 
          part_FTOF_time_layer3[vSC_pindex->at(i)] = vSC_time->at(i); 
          part_FTOF_path_layer3[vSC_pindex->at(i)] = vSC_path->at(i);  
          part_FTOF_sector_layer3[vSC_pindex->at(i)] = vSC_sector->at(i); 
          part_FTOF_component_layer3[vSC_pindex->at(i)] = vSC_component->at(i); 
        }
        if(vSC_layer->at(i) == 2){
          part_FTOF_energy[vSC_pindex->at(i)] = vSC_energy->at(i); 
          part_FTOF_time[vSC_pindex->at(i)] = vSC_time->at(i); 
          part_FTOF_path[vSC_pindex->at(i)] = vSC_path->at(i);  
        }
      }
      // CTOF
      if(vSC_pindex->at(i) >= 0 && vSC_pindex->at(i) < BUFFER && vSC_detector->at(i) == 4){
        part_CTOF_component[vSC_pindex->at(i)] = vSC_component->at(i); 
        part_CTOF_energy[vSC_pindex->at(i)] = vSC_energy->at(i); 
        part_CTOF_time[vSC_pindex->at(i)] = vSC_time->at(i); 
        part_CTOF_path[vSC_pindex->at(i)] = vSC_path->at(i);  
      }
      // CND
      if(vSC_pindex->at(i) >= 0 && vSC_pindex->at(i) < BUFFER && vSC_detector->at(i) == 3){
        part_CND_component[vSC_pindex->at(i)] = vSC_component->at(i); 
        part_CND_energy[vSC_pindex->at(i)] = vSC_energy->at(i); 
        part_CND_time[vSC_pindex->at(i)] = vSC_time->at(i); 
        part_CND_path[vSC_pindex->at(i)] = vSC_path->at(i);  
      }
    }
  }


  // tracking banks  (detectors: DC = 6, BST = 2,  BMT = 1, FMT = 8) 

  if(TRK_Nentries > 0){
    for(int i = 0; i < TRK_Nentries; i++){
      if(vTRK_pindex->at(i) < BUFFER && vTRK_detector->at(i) == 6){     // DC
        part_DC_sector[vTRK_pindex->at(i)] = vTRK_sector->at(i);
      }
    }
  }


  // trajectory crosses  (detectors: 12 = DC region 1 start,  24 = DC region 2 start,  36 = DC region 3 start ) 

  if(Traj_Nentries > 0){
    for(int i = 0; i < Traj_Nentries; i++){
      if(vTraj_pindex->at(i) >= 0 && vTraj_pindex->at(i) < BUFFER && vTraj_detID->at(i) == 10){    
          part_DC_c1x[vTraj_pindex->at(i)] = vTraj_x->at(i);
          part_DC_c1y[vTraj_pindex->at(i)] = vTraj_y->at(i);
          part_DC_c1z[vTraj_pindex->at(i)] = vTraj_z->at(i);
      }
      if(vTraj_pindex->at(i) >= 0 && vTraj_pindex->at(i) < BUFFER && vTraj_detID->at(i) == 22){  
          part_DC_c2x[vTraj_pindex->at(i)] = vTraj_x->at(i);
          part_DC_c2y[vTraj_pindex->at(i)] = vTraj_y->at(i);
          part_DC_c2z[vTraj_pindex->at(i)] = vTraj_z->at(i);
      }
      if(vTraj_pindex->at(i) >= 0 && vTraj_pindex->at(i) < BUFFER && vTraj_detID->at(i) == 34){  
          part_DC_c3x[vTraj_pindex->at(i)] = vTraj_x->at(i);
          part_DC_c3y[vTraj_pindex->at(i)] = vTraj_y->at(i);
          part_DC_c3z[vTraj_pindex->at(i)] = vTraj_z->at(i);
      }
    }
  }


  /// //////////////////////////////////////////////
  /// MC:
  
  if(simulation == true){

    //if(vMC_Header_helicity->size() > 0)  MC_helicity  = vMC_Header_helicity->at(0);
    //if(vMC_Event_npart->size() > 0)  MC_Npart  = vMC_Event_npart->at(0);
    //if(vMC_Event_ebeam->size() > 0)  MC_Ebeam  = vMC_Event_ebeam->at(0);
    //if(vMC_Event_weight->size() > 0) MC_weight = vMC_Event_weight->at(0);

  int MC_Particle_Nentries = 0;
  int MC_Lund_Nentries = 0;

  MC_Particle_Nentries   = vMC_Particle_pid->size();
  //std::cout << " n particles for mc "  << vMC_Particle_pid->size() << std::endl;
  //MC_Lund_Nentries   = vMC_Lund_pid->size(); 

  if(MC_Particle_Nentries > 0){
    for(int i = 0; i < MC_Particle_Nentries; i++){
      if(i < BUFFER){ 
        partMC_pid[i] = vMC_Particle_pid->at(i);
        partMC_px[i] = vMC_Particle_px->at(i);
        partMC_py[i] = vMC_Particle_py->at(i);
        partMC_pz[i] = vMC_Particle_pz->at(i);
        partMC_vx[i] = vMC_Particle_vx->at(i);
        partMC_vy[i] = vMC_Particle_vy->at(i);
        partMC_vz[i] = vMC_Particle_vz->at(i);
        partMC_p[i] = sqrt(partMC_px[i]*partMC_px[i] + partMC_py[i]*partMC_py[i] + partMC_pz[i]*partMC_pz[i]);
        partMC_theta[i] = acos(partMC_pz[i]/partMC_p[i]);
        partMC_phi[i] = atan2(partMC_py[i],partMC_px[i]);
      }
    }
  }

  /*  if(MC_Lund_Nentries > 0){
    for(int i = 0; i < MC_Lund_Nentries; i++){
      if(i < BUFFER){ 
        partLUND_pid[i] = vMC_Lund_pid->at(i);
        partLUND_mass[i] = vMC_Lund_mass->at(i);
        partLUND_E[i] = vMC_Lund_E->at(i);
        partLUND_px[i] = sqrt(pow(partLUND_E[i],2) - pow(partLUND_mass[i],2)) * vMC_Lund_px->at(i);
        partLUND_py[i] = sqrt(pow(partLUND_E[i],2) - pow(partLUND_mass[i],2)) * vMC_Lund_py->at(i);
        partLUND_pz[i] = sqrt(pow(partLUND_E[i],2) - pow(partLUND_mass[i],2)) * vMC_Lund_pz->at(i);
        partLUND_vx[i] = vMC_Lund_vx->at(i);
        partLUND_vy[i] = vMC_Lund_vy->at(i);
        partLUND_vz[i] = vMC_Lund_vz->at(i);
        partLUND_p[i] = sqrt(pow(partLUND_E[i],2) - pow(partLUND_mass[i],2));
        partLUND_theta[i] = acos(partLUND_pz[i]/partLUND_p[i]);
        partLUND_phi[i] = atan2(partLUND_py[i],partLUND_px[i]);
      }
    }
  }

  */
  }

  /// ////////////////////////////////////////////

}


/// //////////////////////////////////////////////////////////////////////////////////////////////////////
/// particle selection:

void select_electron(int run){

  float mom = 0;
  int check = 0;

  int Npart = vpart_pid->size();

  for(Int_t i = 0; i < Npart; i++){
    if( i < BUFFER){

      /// ////////////////////////////////////////////////////////////////////////////////////////////
      /// Forward detector:

      if(part_status[i] >= 2000 && part_status[i] < 4000){

        // PID checks:

        Track_Quality_check[i] = Track_Quality_cut(i);
        FD_eid_default_PID_check[i] = ele_default_PID_cut(i);
        FD_eid_charge_check[i] = ele_charge_cut(i);
        FD_eid_CC_nphe_check[i] = CC_nphe_cut(i);
        FD_eid_EC_outer_vs_EC_inner_check[i] = EC_outer_vs_EC_inner_cut(i);
        FD_eid_EC_sampling_fraction_check[i] = EC_sampling_fraction_cut(i);
        FD_eid_EC_hit_position_fiducial_check[i] = EC_hit_position_fiducial_cut(i);
        FD_eid_DC_hit_position_region1_fiducial_check[i] = DC_hit_position_region1_fiducial_cut(i);
        FD_eid_DC_hit_position_region2_fiducial_check[i] = DC_hit_position_region2_fiducial_cut(i);
        FD_eid_DC_hit_position_region3_fiducial_check[i] = DC_hit_position_region3_fiducial_cut(i);
        FD_eid_DC_z_vertex_check[i] = DC_z_vertex_cut(i);   

        // count statistics

        if(Track_Quality_check[i]) 													Track_Quality_pass += 1;
        if(FD_eid_default_PID_check[i] == true) 	    				                                                FD_eid_default_PID_pass += 1;
        if(FD_eid_charge_check[i]  == true)					                                                        FD_eid_charge_pass += 1;
        if(FD_eid_CC_nphe_check[i] && FD_eid_charge_check[i] && FD_eid_default_PID_check[i]) 			        		FD_eid_CC_nphe_pass += 1;
        if(FD_eid_EC_outer_vs_EC_inner_check[i] && FD_eid_charge_check[i] && FD_eid_default_PID_check[i]) 				FD_eid_EC_outer_vs_EC_inner_pass += 1;
        if(FD_eid_default_PID_check[i] && FD_eid_EC_sampling_fraction_check[i] && FD_eid_charge_check[i]&& FD_eid_default_PID_check[i])	FD_eid_EC_sampling_fraction_pass += 1;
        if(FD_eid_EC_hit_position_fiducial_check[i] && FD_eid_charge_check[i] && FD_eid_default_PID_check[i] ) 				FD_eid_EC_hit_position_fiducial_pass += 1;
        if(FD_eid_DC_hit_position_region1_fiducial_check[i] && FD_eid_charge_check[i] && FD_eid_default_PID_check[i]) 			FD_eid_DC_hit_position_region1_fiducial_pass += 1;
        if(FD_eid_DC_hit_position_region2_fiducial_check[i] && FD_eid_charge_check[i] && FD_eid_default_PID_check[i]) 			FD_eid_DC_hit_position_region2_fiducial_pass += 1;
        if(FD_eid_DC_hit_position_region3_fiducial_check[i] && FD_eid_charge_check[i] && FD_eid_default_PID_check[i]) 			FD_eid_DC_hit_position_region3_fiducial_pass += 1;
        if(FD_eid_DC_z_vertex_check[i] && FD_eid_charge_check[i] && FD_eid_default_PID_check[i]) 			        	FD_eid_DC_z_vertex_pass += 1;

        if(FD_eid_default_PID_check[i] ){ // && FD_eid_charge_check[i] && FD_eid_EC_outer_vs_EC_inner_check[i] ){ // && FD_eid_EC_sampling_fraction_check[i] && FD_eid_EC_hit_position_fiducial_check[i]
	  //&& FD_eid_DC_hit_position_region1_fiducial_check[i] && FD_eid_DC_hit_position_region2_fiducial_check[i] && FD_eid_DC_hit_position_region3_fiducial_check[i] 
          //                             && FD_eid_DC_z_vertex_check[i]){
          FD_eid_all_pass += 1;
          FD_eid_all_check[i] = true;
        }
      }

      /// ////////////////////////////////////////////////////////////////////////////////////////////
      /// Forward tagger:

      if(part_status[i] >= 1000 && part_status[i] < 2000){
        FT_eid_charge_check[i] = FT_eid_charge_cut(i);
        FT_eid_PID_check[i] = FT_eid_PID_cut(i);
        FT_eid_FTCAL_fiducial_check[i] = FT_eid_FTCAL_fiducial_cut(i);
        FT_eid_FTTRK_fiducial_check[i] = FT_eid_FTTRK_fiducial_cut(i);
        FT_eid_FTHODO_fiducial_check[i] = FT_eid_FTHODO_fiducial_cut(i);
        FT_eid_energy_vs_radius_check[i] = FT_eid_energy_vs_radius_cut(i);

        if(FT_eid_charge_check[i]) 								FT_eid_charge_pass += 1;
        if(FT_eid_PID_check[i]) 								FT_eid_PID_pass += 1;
        if(FT_eid_FTCAL_fiducial_check[i] && FT_eid_charge_check[i] && FT_eid_PID_check[i]) 	FT_eid_FTCAL_fiducial_pass += 1;
        if(FT_eid_FTTRK_fiducial_check[i] && FT_eid_charge_check[i] && FT_eid_PID_check[i]) 	FT_eid_FTTRK_fiducial_pass += 1;
        if(FT_eid_FTHODO_fiducial_check[i] && FT_eid_charge_check[i] && FT_eid_PID_check[i]) 	FT_eid_FTHODO_fiducial_pass += 1;
        if(FT_eid_energy_vs_radius_check[i] && FT_eid_charge_check[i] && FT_eid_PID_check[i]) 	FT_eid_energy_vs_radius_pass += 1;

        if(FT_eid_PID_check[i] ){// && FT_eid_charge_check[i] && FT_eid_FTCAL_fiducial_check[i] && FT_eid_FTTRK_fiducial_check[i] && FT_eid_FTHODO_fiducial_check[i] && FT_eid_energy_vs_radius_check[i]){
          FT_eid_all_pass += 1;
          FT_eid_all_check[i] = true;
        }
      }

      // create vector with electron lorentz vectors if electron PID is succesfull:

      bool selector;

      if(use_FT == false){ FT_eid_all_check[i] = false;}
      if(use_FD == false){ FD_eid_all_check[i] = false;}

      /// ////////////////////////////////////////////////////////////////
      /// pick particle index and sort by momentum
 
      mom = 0;
      check = 0;

      for(int k = 0; k < Npart; k++){
        if(k < BUFFER){

          if(use_own_PID_electron == true){ 
            selector = (FD_eid_all_check[k] && part_status[i] >= 2000 && part_status[i] < 4000)|| (FT_eid_all_check[k] && part_status[i] >= 1000 && part_status[i] < 2000); 
          }
          else{ 
            selector = (FD_eid_default_PID_check[k] && part_status[i] >= 2000 && part_status[i] < 4000) || (FT_eid_PID_check[k]  && part_status[i] >= 1000 && part_status[i] < 2000); 
          }       

          if(selector){

            for(int j = 0; j < i; j++){
              if(k == e_ind[j]) check = -1;
            }

            if(sqrt(vpart_px->at(k)*vpart_px->at(k) + vpart_py->at(k)*vpart_py->at(k) + vpart_pz->at(k)*vpart_pz->at(k)) > mom  && check != -1){   

	      mom = sqrt(vpart_px->at(k)*vpart_px->at(k) + vpart_py->at(k)*vpart_py->at(k) + vpart_pz->at(k)*vpart_pz->at(k)); 
              e_ind[i] = k;
              e_count += 1;
	    }
            check = 0;
	  }
        }
      }
    }
  }

  /// ///////////////////////////////////////////////////////////////////
  /// assign properties

  for(int i = 0; i < BUFFER; i++){
    if(e_ind[i] != -1){
      e_vx[i] = vpart_vx->at(e_ind[i]);
      e_vy[i] = vpart_vy->at(e_ind[i]);
      e_vz[i] = vpart_vz->at(e_ind[i]);
      e_beta[i] = vpart_beta->at(e_ind[i]);
      double p = sqrt(vpart_px->at(e_ind[i])*vpart_px->at(e_ind[i]) + vpart_py->at(e_ind[i])*vpart_py->at(e_ind[i]) + vpart_pz->at(e_ind[i])*vpart_pz->at(e_ind[i]));
      p4_ele[i].SetPxPyPzE(vpart_px->at(e_ind[i]), vpart_py->at(e_ind[i]), vpart_pz->at(e_ind[i]), sqrt(p*p + m_e*m_e));
      e_FTOF_sec[i] = part_FTOF_sector_layer2[e_ind[i]];
      e_PCAL_sec[i] = part_Cal_PCAL_sector[e_ind[i]];
      eventNumber[i] =  event;

      if(part_status[e_ind[i]] >= 1000 && part_status[e_ind[i]] < 2000) ele_detect[i] = 1; 
      if(part_status[e_ind[i]] >= 2000 && part_status[e_ind[i]] < 4000) ele_detect[i] = 2; 
    }
  }


}

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// proton selector
/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void select_proton(int run){

  float mom = 0;
  float mom_new = 0;
  int check = 0;
  int Npart = vpart_pid->size();


  for(Int_t i = 0; i < Npart; i++){
    if( i < BUFFER){

    /// ////////////////////////////////////////////////////////////////////////////////////////////
    /// Forward detector:

    if(part_status[i] >= 2000 && part_status[i] < 4000){

      // PID checks:

      FD_protid_default_PID_check[i] = prot_default_PID_cut(i);
      FD_protid_charge_check[i] = prot_charge_cut(i);
      FD_protid_DC_hit_position_region1_fiducial_check[i] = DC_hit_position_region1_fiducial_cut_hadrons_positive(i);
      FD_protid_DC_hit_position_region2_fiducial_check[i] = DC_hit_position_region2_fiducial_cut_hadrons_positive(i);
      FD_protid_DC_hit_position_region3_fiducial_check[i] = DC_hit_position_region3_fiducial_cut_hadrons_positive(i);
      FD_protid_beta_check[i] = prot_beta_cut(i, run);
      FD_protid_delta_beta_check[i] = prot_delta_beta_cut(i, run);
      FD_protid_tofmass_check[i] = prot_tofmass_cut(i, run);
      FD_protid_maximum_probability_check[i] = maximum_probability_cut(i, 2212, 0.27, 99.73, run); // check if the hypothesis proton is fulfilled for the particle - proton: 2212,  pip: 211,  Kp: 321 
										         // conflvl = required probability, that it is the hypothesis particle in % - 3sigma = 0.27, 2sigma = 4.55, 1sigma = 31.73
                                                                                         // anticonflvl = maximum acceptable probability for another candidate in %
      FD_protid_delta_vz_check[i] = prot_delta_vz_cut(i);

      if(FD_protid_default_PID_check[i] == true) 										FD_protid_default_PID_pass += 1;
      if(FD_protid_charge_check[i] == true) 											FD_protid_charge_pass += 1;
      if(FD_protid_DC_hit_position_region1_fiducial_check[i] && FD_protid_charge_check[i] && FD_protid_default_PID_check[i]) 	FD_protid_DC_hit_position_region1_fiducial_pass += 1;
      if(FD_protid_DC_hit_position_region2_fiducial_check[i] && FD_protid_charge_check[i] && FD_protid_default_PID_check[i]) 	FD_protid_DC_hit_position_region2_fiducial_pass += 1;
      if(FD_protid_DC_hit_position_region3_fiducial_check[i] && FD_protid_charge_check[i] && FD_protid_default_PID_check[i]) 	FD_protid_DC_hit_position_region3_fiducial_pass += 1;
      if(FD_protid_beta_check[i] && FD_protid_charge_check[i] && FD_protid_default_PID_check[i]) 				FD_protid_beta_pass += 1;
      if(FD_protid_delta_beta_check[i] && FD_protid_charge_check[i] && FD_protid_default_PID_check[i]) 				FD_protid_delta_beta_pass += 1;
      if(FD_protid_tofmass_check[i] && FD_protid_charge_check[i] && FD_protid_default_PID_check[i]) 				FD_protid_tofmass_pass += 1;
      if(FD_protid_maximum_probability_check[i] && FD_protid_charge_check[i] && FD_protid_default_PID_check[i]) 		FD_protid_maximum_probability_pass += 1;
      if(FD_protid_delta_vz_check[i] && FD_protid_charge_check[i] && FD_protid_default_PID_check[i]) 				FD_protid_delta_vz_pass += 1;


      if(cut_maximum_probability_prot == true){
        if(FD_protid_default_PID_check[i] && FD_protid_charge_check[i] && FD_protid_DC_hit_position_region1_fiducial_check[i] && FD_protid_DC_hit_position_region2_fiducial_check[i] 
                                          && FD_protid_DC_hit_position_region3_fiducial_check[i] && FD_protid_maximum_probability_check[i]){
          FD_protid_all_pass += 1;
          FD_protid_all_check[i] = true;
        }
      }
      else if(cut_beta_vs_p_prot == true){
        if(FD_protid_default_PID_check[i] && FD_protid_charge_check[i] && FD_protid_DC_hit_position_region1_fiducial_check[i] && FD_protid_DC_hit_position_region2_fiducial_check[i] 
                                          && FD_protid_DC_hit_position_region3_fiducial_check[i] && FD_protid_beta_check[i]){
          FD_protid_all_pass += 1;
          FD_protid_all_check[i] = true;
        }
      }
      else if(cut_beta_vs_p_prot == true && cut_deltabeta_prot == true && cut_tofmass_prot == true){
        if(FD_protid_default_PID_check[i] && FD_protid_charge_check[i] && FD_protid_DC_hit_position_region1_fiducial_check[i] && FD_protid_DC_hit_position_region2_fiducial_check[i] 
                                          && FD_protid_DC_hit_position_region3_fiducial_check[i] && FD_protid_beta_check[i] && FD_protid_tofmass_check[i] && FD_protid_delta_beta_check[i]){
          FD_protid_all_pass += 1;
          FD_protid_all_check[i] = true;
        }
      }
      else{
        if(FD_protid_default_PID_check[i] && FD_protid_charge_check[i] && FD_protid_DC_hit_position_region1_fiducial_check[i] && FD_protid_DC_hit_position_region2_fiducial_check[i] 
                                          && FD_protid_DC_hit_position_region3_fiducial_check[i]){
          FD_protid_all_pass += 1;
          FD_protid_all_check[i] = true;
        }
      }
    }


    /// ////////////////////////////////////////////////////////////////////////////////////////////
    /// Central detector:

    if(part_status[i] >= 4000){

      // PID checks:

      CD_protid_default_PID_check[i] = prot_default_PID_cut(i);
      CD_protid_charge_check[i] = prot_charge_cut(i);
      CD_protid_beta_check[i] = CD_prot_beta_cut(i, run);
      CD_protid_maximum_probability_check[i] = CD_maximum_probability_cut(i, 2212, 0.27, 99.73, run);
      CD_protid_delta_vz_check[i] = CD_prot_delta_vz_cut(i);

      if(CD_protid_default_PID_check[i]) 											CD_protid_default_PID_pass += 1;
      if(CD_protid_charge_check[i]) 												CD_protid_charge_pass += 1;
      if(CD_protid_beta_check[i] && CD_protid_charge_check[i] && CD_protid_default_PID_check[i]) 				CD_protid_beta_pass += 1;
      if(CD_protid_maximum_probability_check[i] && CD_protid_charge_check[i] && CD_protid_default_PID_check[i]) 		CD_protid_maximum_probability_pass += 1;
      if(CD_protid_delta_vz_check[i] && CD_protid_charge_check[i] && CD_protid_default_PID_check[i]) 				CD_protid_delta_vz_pass += 1;

      if(CD_cut_maximum_probability_prot == true){
        if(CD_protid_charge_check[i] && CD_protid_maximum_probability_check[i]){
          CD_protid_all_pass += 1;
          CD_protid_all_check[i] = true;
        }
      }
      else if(CD_cut_beta_vs_p_prot == true){
        if(CD_protid_charge_check[i] && CD_protid_beta_check[i]){
          CD_protid_all_pass += 1;
          CD_protid_all_check[i] = true;
        }
      }
      else{
        if(CD_protid_default_PID_check[i] && CD_protid_charge_check[i]){
          CD_protid_all_pass += 1;
          CD_protid_all_check[i] = true;
        }
      }


      /// ///////////////////////////////////////////////////////////////////////
      /// 2 GeV condition:

      if(Ebeam < 3){
        if(CD_protid_charge_check[i]){
          CD_protid_all_check[i] = true;
        }
      }
    
      /// //////////////////////////////////////////////////////////////////////


    }




    /// ///////////////////////////////////////////////////////////////
    /// Create proton selector

    bool selector;

    if(use_FD == false) FD_protid_all_check[i] = false;
    if(use_CD == false) CD_protid_all_check[i] = false;

 
    /// ////////////////////////////////////////////////////////////////
    /// pick particle index and sort by momentum

    if( i < BUFFER){
      mom = 0;
      check = 0;
      for(int k = 0; k < Npart; k++){
        if(k < BUFFER){
        if(use_own_PID_proton == true){ 
          selector = (FD_protid_all_check[k] && part_status[i] >= 2000 && part_status[i] < 4000) || (CD_protid_all_check[k] && part_status[i] >= 4000); 
        }
        else{ 
          selector = (FD_protid_default_PID_check[k] && part_status[i] >= 2000 && part_status[i] < 4000) || (CD_protid_default_PID_check[k] && part_status[i] >= 4000); 
        }

        /// ///////////////////////////////////////////////////////////////////////
        /// 2 GeV condition:

        if(Ebeam < 3){
          selector = (FD_protid_all_check[k] && part_status[i] >= 2000 && part_status[i] < 4000) || (CD_protid_all_check[k] && part_status[i] >= 4000); 
        }
    
        /// //////////////////////////////////////////////////////////////////////

        if(selector){
          for(int j = 0; j < i; j++){
            if(k == p_ind[j]) check = -1;
          }
          mom_new = sqrt(vpart_px->at(k)*vpart_px->at(k) + vpart_py->at(k)*vpart_py->at(k) + vpart_pz->at(k)*vpart_pz->at(k));
          if(mom_new > mom && mom_new > 0 && check != -1){        
	    mom = sqrt(vpart_px->at(k)*vpart_px->at(k) + vpart_py->at(k)*vpart_py->at(k) + vpart_pz->at(k)*vpart_pz->at(k)); 
            p_ind[i] = k;
            p_count += 1; 
	  }
          check = 0;
}
	}  
      }
    }
  }
}

  /// //////////////////////////////////////////////////////////////////
  /// Assign properties

  for(int i = 0; i < BUFFER; i++){
    if(p_ind[i] != -1){
      p_vx[i] = vpart_vx->at(p_ind[i]);
      p_vy[i] = vpart_vy->at(p_ind[i]);
      p_vz[i] = vpart_vz->at(p_ind[i]);
      p_beta[i] = vpart_beta->at(p_ind[i]);
      double p = sqrt(vpart_px->at(p_ind[i])*vpart_px->at(p_ind[i]) + vpart_py->at(p_ind[i])*vpart_py->at(p_ind[i]) + vpart_pz->at(p_ind[i])*vpart_pz->at(p_ind[i]));
      p4_prot[i].SetPxPyPzE(vpart_px->at(p_ind[i]), vpart_py->at(p_ind[i]), vpart_pz->at(p_ind[i]), sqrt(p*p + m_p*m_p));
      p_FTOF_sec[i] = part_FTOF_sector_layer2[p_ind[i]];
      p_PCAL_sec[i] = part_Cal_PCAL_sector[p_ind[i]];
      if(part_status[p_ind[i]] >= 2000 && part_status[p_ind[i]] < 4000) prot_detect[i] = 2; 
      if(part_status[p_ind[i]] >= 4000)                                 prot_detect[i] = 3; 
    }
  }  

}

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// neutron selector
/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void select_neutron(int run){

  float mom = 0;
  int check = 0;
  int Npart = vpart_pid->size();


  for(Int_t i = 0; i < Npart; i++){
    if( i < BUFFER){

    /// ////////////////////////////////////////////////////////////////////////////////////////////
    /// Forward detector:

    if(part_status[i] >= 2000 && part_status[i] < 4000){

      // PID checks:

      FD_neutrid_default_PID_check[i] = neutr_default_PID_cut(i);
      FD_neutrid_charge_check[i] = neutr_charge_cut(i);
      FD_neutrid_beta_check[i] = neutr_beta_cut(i, run);
      FD_neutrid_delta_beta_check[i] = neutr_delta_beta_cut(i, run);
      FD_neutrid_tofmass_check[i] = neutr_tofmass_cut(i, run);
      FD_neutrid_delta_vz_check[i] = neutr_delta_vz_cut(i);

      if(FD_neutrid_default_PID_check[i] == true) 									FD_neutrid_default_PID_pass += 1;
      if(FD_neutrid_charge_check[i]  == true) 										FD_neutrid_charge_pass += 1;
      if(FD_neutrid_beta_check[i] && FD_neutrid_charge_check[i] && FD_neutrid_default_PID_check[i]) 			FD_neutrid_beta_pass += 1;
      if(FD_neutrid_delta_beta_check[i] && FD_neutrid_charge_check[i] && FD_neutrid_default_PID_check[i]) 		FD_neutrid_delta_beta_pass += 1;
      if(FD_neutrid_tofmass_check[i] && FD_neutrid_charge_check[i] && FD_neutrid_default_PID_check[i]) 			FD_neutrid_tofmass_pass += 1;
      if(FD_neutrid_delta_vz_check[i] && FD_neutrid_charge_check[i] && FD_neutrid_default_PID_check[i]) 		FD_neutrid_delta_vz_pass += 1;

      if(cut_beta_vs_p_only_neutron == true){
        if(part_p[i] > 0.1 && FD_neutrid_default_PID_check[i] && FD_neutrid_charge_check[i] && FD_neutrid_beta_check[i]){
          FD_neutrid_all_pass += 1;
          FD_neutrid_all_check[i] = true;
        } 
      }
      else{
        if(part_p[i] > 0.1 && FD_neutrid_default_PID_check[i] && FD_neutrid_charge_check[i] && FD_neutrid_beta_check[i] && FD_neutrid_delta_beta_check[i] && FD_neutrid_tofmass_check[i]){
          FD_neutrid_all_pass += 1;
          FD_neutrid_all_check[i] = true;
        }
      }
    }

    /// ////////////////////////////////////////////////////////////////////////////////////////////
    /// Central detector:

    if(part_status[i] >= 4000){

      // PID checks:

      CD_neutrid_default_PID_check[i] = neutr_default_PID_cut(i);
      CD_neutrid_charge_check[i] = neutr_charge_cut(i);
      CD_neutrid_beta_check[i] = CD_neutr_beta_cut(i, run);
      CD_neutrid_delta_vz_check[i] = CD_neutr_delta_vz_cut(i);

      if(CD_neutrid_default_PID_check[i]) 									 	CD_neutrid_default_PID_pass += 1;
      if(CD_neutrid_charge_check[i]) 										 	CD_neutrid_charge_pass += 1;
      if(CD_neutrid_beta_check[i] && CD_neutrid_charge_check[i] && CD_neutrid_default_PID_check[i]) 		 	CD_neutrid_beta_pass += 1;
      if(CD_neutrid_delta_vz_check[i] && CD_neutrid_charge_check[i] && CD_neutrid_default_PID_check[i]) 		CD_neutrid_delta_vz_pass += 1;

      if(cut_beta_vs_p_only_neutron == true){
        if(part_p[i] > 0.1 && CD_neutrid_default_PID_check[i] && CD_neutrid_charge_check[i] && CD_neutrid_beta_check[i]){
          CD_neutrid_all_pass += 1;
          CD_neutrid_all_check[i] = true;
        }
      }
      else{
        if(part_p[i] > 0.1 && CD_neutrid_default_PID_check[i] && CD_neutrid_charge_check[i] && CD_neutrid_beta_check[i]){
          CD_neutrid_all_pass += 1;
          CD_neutrid_all_check[i] = true;
        }
      }
    }

    /// ///////////////////////////////////////////////////////////////
    /// Create neutron selector

    bool selector;

    if(use_FD == false) FD_neutrid_all_check[i] = false;
    if(use_CD == false) CD_neutrid_all_check[i] = false;


    /// ////////////////////////////////////////////////////////////////
    /// pick particle index and sort by momentum

    if( i < BUFFER){
      mom = 0;
      check = 0;
      for(int k = 0; k < Npart; k++){
        if(k < BUFFER){
        if(use_own_PID_neutron == true){ 
          selector = (FD_neutrid_all_check[k] && part_status[i] >= 2000 && part_status[i] < 4000) || (CD_neutrid_all_check[k] && part_status[i] >= 4000); 
        }
        else{ 
          selector = (FD_neutrid_default_PID_check[k] && part_status[i] >= 2000 && part_status[i] < 4000) || (CD_neutrid_default_PID_check[k] && part_status[i] >= 4000); 
        }

        if(selector){
          for(int j = 0; j < i; j++){
            if(k == n_ind[j]) check = -1;
          }
          if(sqrt(vpart_px->at(k)*vpart_px->at(k) + vpart_py->at(k)*vpart_py->at(k) + vpart_pz->at(k)*vpart_pz->at(k)) > mom  && check != -1){        
	    mom = sqrt(vpart_px->at(k)*vpart_px->at(k) + vpart_py->at(k)*vpart_py->at(k) + vpart_pz->at(k)*vpart_pz->at(k)); 
            n_ind[i] = k;
            n_count += 1;  
	  }
          check = 0;
}
	}  
      }
    }
  }
}

  /// //////////////////////////////////////////////////////////////////
  /// Assign properties

  for(int i = 0; i < BUFFER; i++){
    if(n_ind[i] != -1){
      n_vx[i] = vpart_vx->at(n_ind[i]);
      n_vy[i] = vpart_vy->at(n_ind[i]);
      n_vz[i] = vpart_vz->at(n_ind[i]);
      n_beta[i] = vpart_beta->at(n_ind[i]);
      double p = sqrt(vpart_px->at(n_ind[i])*vpart_px->at(n_ind[i]) + vpart_py->at(n_ind[i])*vpart_py->at(n_ind[i]) + vpart_pz->at(n_ind[i])*vpart_pz->at(n_ind[i]));
      p4_neutr[i].SetPxPyPzE(vpart_px->at(n_ind[i]), vpart_py->at(n_ind[i]), vpart_pz->at(n_ind[i]), sqrt(p*p + m_n*m_n));
      if(part_status[n_ind[i]] >= 2000 && part_status[n_ind[i]] < 4000) neutr_detect[i] = 2; 
      if(part_status[n_ind[i]] >= 4000)                                 neutr_detect[i] = 3; 
    }
  }  


}

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// pip selector
/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void select_pip(int run){

  float mom = 0;
  int check = 0;
  int Npart = vpart_pid->size();

  for(Int_t i = 0; i < Npart; i++){
    if( i < BUFFER){

    /// ////////////////////////////////////////////////////////////////////////////////////////////
    /// Forward detector:

    if(part_status[i] >= 2000 && part_status[i] < 4000){

      // PID checks:

      FD_pipid_default_PID_check[i] = pip_default_PID_cut(i);
      FD_pipid_charge_check[i] = pip_charge_cut(i);
      FD_pipid_DC_hit_position_region1_fiducial_check[i] = DC_hit_position_region1_fiducial_cut_hadrons_positive(i);
      FD_pipid_DC_hit_position_region2_fiducial_check[i] = DC_hit_position_region2_fiducial_cut_hadrons_positive(i);
      FD_pipid_DC_hit_position_region3_fiducial_check[i] = DC_hit_position_region3_fiducial_cut_hadrons_positive(i);
      FD_pipid_beta_check[i] = pip_beta_cut(i, run);
      FD_pipid_delta_beta_check[i] = pip_delta_beta_cut(i, run);
      FD_pipid_tofmass_check[i] = pip_tofmass_cut(i, run);
      FD_pipid_maximum_probability_check[i] = maximum_probability_cut(i, 211, 0.27, 99.73, run);  // check if the hypothesis proton is fulfilled for the particle - proton: 2212,  pip: 211,  Kp: 321 
										        // conflvl = required probability, that it is the hypothesis particle in % - 3sigma = 0.27, 2sigma = 4.55, 1sigma = 31.73
                                                                                        // anticonflvl = maximum acceptable probability for another candidate in %
      FD_pipid_delta_vz_check[i] = pip_delta_vz_cut(i);

      if(FD_pipid_default_PID_check[i]  == true) 										FD_pipid_default_PID_pass += 1;
      if(FD_pipid_charge_check[i]  == true) 											FD_pipid_charge_pass += 1;
      if(FD_pipid_DC_hit_position_region1_fiducial_check[i] && FD_pipid_charge_check[i] && FD_pipid_default_PID_check[i]) 	FD_pipid_DC_hit_position_region1_fiducial_pass += 1;
      if(FD_pipid_DC_hit_position_region2_fiducial_check[i] && FD_pipid_charge_check[i] && FD_pipid_default_PID_check[i]) 	FD_pipid_DC_hit_position_region2_fiducial_pass += 1;
      if(FD_pipid_DC_hit_position_region3_fiducial_check[i] && FD_pipid_charge_check[i] && FD_pipid_default_PID_check[i]) 	FD_pipid_DC_hit_position_region3_fiducial_pass += 1;
      if(FD_pipid_beta_check[i] && FD_pipid_charge_check[i] && FD_pipid_default_PID_check[i]) 					FD_pipid_beta_pass += 1;
      if(FD_pipid_delta_beta_check[i] && FD_pipid_charge_check[i] && FD_pipid_default_PID_check[i]) 				FD_pipid_delta_beta_pass += 1;
      if(FD_pipid_tofmass_check[i] && FD_pipid_charge_check[i] && FD_pipid_default_PID_check[i]) 				FD_pipid_tofmass_pass += 1;
      if(FD_pipid_maximum_probability_check[i] && FD_pipid_charge_check[i] && FD_pipid_default_PID_check[i]) 			FD_pipid_maximum_probability_pass += 1;
      if(FD_pipid_delta_vz_check[i] && FD_pipid_charge_check[i] && FD_pipid_default_PID_check[i]) 				FD_pipid_delta_vz_pass += 1;

      if(cut_maximum_probability_pip == true){
        if(FD_pipid_default_PID_check[i] && FD_pipid_charge_check[i] && FD_pipid_DC_hit_position_region1_fiducial_check[i] && FD_pipid_DC_hit_position_region2_fiducial_check[i] 
                                         && FD_pipid_DC_hit_position_region3_fiducial_check[i] && FD_pipid_maximum_probability_check[i]){
          FD_pipid_all_pass += 1;
          FD_pipid_all_check[i] = true;
        }
      }
      else if(cut_beta_vs_p_pip == true){
        if(FD_pipid_default_PID_check[i] && FD_pipid_charge_check[i] && FD_pipid_DC_hit_position_region1_fiducial_check[i] && FD_pipid_DC_hit_position_region2_fiducial_check[i] 
                                         && FD_pipid_DC_hit_position_region3_fiducial_check[i] && FD_pipid_beta_check[i]){
          FD_pipid_all_pass += 1;
          FD_pipid_all_check[i] = true;
        }
      }
      else if(cut_beta_vs_p_pip == true && cut_deltabeta_pip == true && cut_tofmass_pip == true){
        if(FD_pipid_default_PID_check[i] && FD_pipid_charge_check[i] && FD_pipid_DC_hit_position_region1_fiducial_check[i] && FD_pipid_DC_hit_position_region2_fiducial_check[i] 
                                         && FD_pipid_DC_hit_position_region3_fiducial_check[i] && FD_pipid_beta_check[i] && FD_pipid_delta_beta_check[i] && FD_pipid_tofmass_check[i]){
          FD_pipid_all_pass += 1;
          FD_pipid_all_check[i] = true;
        }
      }
      else{
        if(FD_pipid_default_PID_check[i] && FD_pipid_charge_check[i] && FD_pipid_DC_hit_position_region1_fiducial_check[i] && FD_pipid_DC_hit_position_region2_fiducial_check[i] 
                                         && FD_pipid_DC_hit_position_region3_fiducial_check[i]){
          FD_pipid_all_pass += 1;
          FD_pipid_all_check[i] = true;
        }
      }
    }
 
    /// ////////////////////////////////////////////////////////////////////////////////////////////
    /// Central detector:

    if(part_status[i] >= 4000){

      // PID checks:

      CD_pipid_default_PID_check[i] = pip_default_PID_cut(i);
      CD_pipid_charge_check[i] = pip_charge_cut(i);
      CD_pipid_beta_check[i] = CD_pip_beta_cut(i, run);
      CD_pipid_maximum_probability_check[i] = CD_maximum_probability_cut(i, 211, 0.27, 99.73, run);
      CD_pipid_delta_vz_check[i] = CD_pip_delta_vz_cut(i);

      if(CD_pipid_default_PID_check[i]) 											CD_pipid_default_PID_pass += 1;
      if(CD_pipid_charge_check[i]) 												CD_pipid_charge_pass += 1;
      if(CD_pipid_beta_check[i] && CD_pipid_charge_check[i] && CD_pipid_default_PID_check[i]) 					CD_pipid_beta_pass += 1;
      if(CD_pipid_maximum_probability_check[i] && CD_pipid_charge_check[i] && CD_pipid_default_PID_check[i]) 			CD_pipid_maximum_probability_pass += 1;
      if(CD_pipid_delta_vz_check[i] && CD_pipid_charge_check[i] && CD_pipid_default_PID_check[i]) 				CD_pipid_delta_vz_pass += 1;

      if(CD_cut_maximum_probability_pip == true){
        if(CD_pipid_charge_check[i] && CD_pipid_maximum_probability_check[i]){
          CD_pipid_all_pass += 1;
          CD_pipid_all_check[i] = true;
        }
      }
      else if(CD_cut_beta_vs_p_pip == true){
        if(CD_pipid_charge_check[i] && CD_pipid_beta_check[i]){
          CD_pipid_all_pass += 1;
          CD_pipid_all_check[i] = true;
        }
      }
      else{
        if(CD_pipid_default_PID_check[i] && CD_pipid_charge_check[i]){
          CD_pipid_all_pass += 1;
          CD_pipid_all_check[i] = true;
        }
      }
    }

    /// ///////////////////////////////////////////////////////////////
    /// Create pip selector

    bool selector;

    if(use_FD == false) FD_pipid_all_check[i] = false;
    if(use_CD == false) CD_pipid_all_check[i] = false;


    /// ////////////////////////////////////////////////////////////////
    /// pick particle index and sort by momentum
  
    if( i < BUFFER){
      mom = 0;
      check = 0;
      for(int k = 0; k < Npart; k++){
        if(k < BUFFER){

        if(use_own_PID_pip == true){ 
          selector = (FD_pipid_all_check[k] && part_status[i] >= 2000 && part_status[i] < 4000) || (CD_pipid_all_check[k] && part_status[i] >= 4000); 
        }
        else{ 
          selector = (FD_pipid_default_PID_check[k] && part_status[i] >= 2000 && part_status[i] < 4000) || (CD_pipid_default_PID_check[k] && part_status[i] >= 4000); 
        }

        if(selector == true){
          for(int j = 0; j < i; j++){
            if(k == pip_ind[j]) check = -1;
          }
          if(sqrt(vpart_px->at(k)*vpart_px->at(k) + vpart_py->at(k)*vpart_py->at(k) + vpart_pz->at(k)*vpart_pz->at(k)) > mom  && check != -1){        
	    mom = sqrt(vpart_px->at(k)*vpart_px->at(k) + vpart_py->at(k)*vpart_py->at(k) + vpart_pz->at(k)*vpart_pz->at(k)); 
            pip_ind[i] = k;
            pip_count += 1;  
	  }
          check = 0;
}
	}  
      }
    }
  }
}


  /// //////////////////////////////////////////////////////////////////
  /// Assign properties

  for(int i = 0; i < BUFFER; i++){
    if(pip_ind[i] != -1){
      pip_vx[i] = vpart_vx->at(pip_ind[i]);
      pip_vy[i] = vpart_vy->at(pip_ind[i]);
      pip_vz[i] = vpart_vz->at(pip_ind[i]);
      pip_beta[i] = vpart_beta->at(pip_ind[i]);
      double p = sqrt(vpart_px->at(pip_ind[i])*vpart_px->at(pip_ind[i]) + vpart_py->at(pip_ind[i])*vpart_py->at(pip_ind[i]) + vpart_pz->at(pip_ind[i])*vpart_pz->at(pip_ind[i]));
      p4_pip[i].SetPxPyPzE(vpart_px->at(pip_ind[i]), vpart_py->at(pip_ind[i]), vpart_pz->at(pip_ind[i]), sqrt(p*p + m_pip*m_pip));
      pip_FTOF_sec[i] = part_FTOF_sector_layer2[pip_ind[i]];
      if(part_status[pip_ind[i]] >= 2000 && part_status[pip_ind[i]] < 4000) pip_detect[i] = 2; 
      if(part_status[pip_ind[i]] >= 4000)                                   pip_detect[i] = 3; 
    }
  }  


}

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// pim selector
/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void select_pim(int run){

  float mom = 0;
  int check = 0;
  int Npart = vpart_pid->size();
 
  for(Int_t i = 0; i < Npart; i++){
    if( i < BUFFER){

    /// ////////////////////////////////////////////////////////////////////////////////////////////
    /// Forward detector:

    if(part_status[i] >= 2000 && part_status[i] < 4000){

      // PID checks:

      FD_pimid_default_PID_check[i] = pim_default_PID_cut(i);
      FD_pimid_charge_check[i] = pim_charge_cut(i);
      FD_pimid_ele_reject_check[i] = pim_ele_reject_cut(i);
      FD_pimid_EC_outer_vs_EC_inner_check[i] = pim_EC_outer_vs_EC_inner_cut(i);
      FD_pimid_DC_hit_position_region1_fiducial_check[i] = DC_hit_position_region1_fiducial_cut_hadrons_negative(i);
      FD_pimid_DC_hit_position_region2_fiducial_check[i] = DC_hit_position_region2_fiducial_cut_hadrons_negative(i);
      FD_pimid_DC_hit_position_region3_fiducial_check[i] = DC_hit_position_region3_fiducial_cut_hadrons_negative(i);
      FD_pimid_beta_check[i] = pim_beta_cut(i, run);
      FD_pimid_delta_beta_check[i] = pim_delta_beta_cut(i, run);
      FD_pimid_tofmass_check[i] = pim_tofmass_cut(i, run);
      FD_pimid_maximum_probability_check[i] = maximum_probability_cut(i, -211, 0.27, 99.73, run); // check if the hypothesis proton is fulfilled for the particle - proton: 2212,  pip: 211,  Kp: 321 
										        // conflvl = required probability, that it is the hypothesis particle in % - 3sigma = 0.27, 2sigma = 4.55, 1sigma = 31.73
                                                                                        // anticonflvl = maximum acceptable probability for another candidate in %
      FD_pimid_delta_vz_check[i] = pim_delta_vz_cut(i);

      if(FD_pimid_default_PID_check[i]  == true) 										FD_pimid_default_PID_pass += 1;
      if(FD_pimid_charge_check[i]  == true) 											FD_pimid_charge_pass += 1;
      if(FD_pimid_DC_hit_position_region1_fiducial_check[i] && FD_pimid_charge_check[i] && FD_pimid_default_PID_check[i]) 	FD_pimid_DC_hit_position_region1_fiducial_pass += 1;
      if(FD_pimid_DC_hit_position_region2_fiducial_check[i] && FD_pimid_charge_check[i] && FD_pimid_default_PID_check[i]) 	FD_pimid_DC_hit_position_region2_fiducial_pass += 1;
      if(FD_pimid_DC_hit_position_region3_fiducial_check[i] && FD_pimid_charge_check[i] && FD_pimid_default_PID_check[i]) 	FD_pimid_DC_hit_position_region3_fiducial_pass += 1;
      if(FD_pimid_beta_check[i] && FD_pimid_charge_check[i] && FD_pimid_default_PID_check[i]) 					FD_pimid_beta_pass += 1;
      if(FD_pimid_delta_beta_check[i] && FD_pimid_charge_check[i] && FD_pimid_default_PID_check[i]) 				FD_pimid_delta_beta_pass += 1;
      if(FD_pimid_tofmass_check[i] && FD_pimid_charge_check[i] && FD_pimid_default_PID_check[i]) 				FD_pimid_tofmass_pass += 1;
      if(FD_pimid_maximum_probability_check[i] && FD_pimid_charge_check[i] && FD_pimid_default_PID_check[i]) 			FD_pimid_maximum_probability_pass += 1;
      if(FD_pimid_delta_vz_check[i] && FD_pimid_charge_check[i] && FD_pimid_default_PID_check[i]) 				FD_pimid_delta_vz_pass += 1;

      if(cut_maximum_probability_pim == true){
        if(FD_pimid_default_PID_check[i] && FD_pimid_charge_check[i] && FD_pimid_ele_reject_check[i] && FD_pimid_EC_outer_vs_EC_inner_check[i] && FD_pimid_DC_hit_position_region1_fiducial_check[i] 
                                         && FD_pimid_DC_hit_position_region2_fiducial_check[i] && FD_pimid_DC_hit_position_region3_fiducial_check[i] && FD_pimid_maximum_probability_check[i]){
          FD_pimid_all_pass += 1;
          FD_pimid_all_check[i] = true;
        }
      }
      else if(cut_beta_vs_p_pim == true){
        if(FD_pimid_default_PID_check[i] && FD_pimid_charge_check[i] && FD_pimid_ele_reject_check[i] && FD_pimid_EC_outer_vs_EC_inner_check[i] && FD_pimid_DC_hit_position_region1_fiducial_check[i] 
                                         && FD_pimid_DC_hit_position_region2_fiducial_check[i] && FD_pimid_DC_hit_position_region3_fiducial_check[i] && FD_pimid_beta_check[i]){
          FD_pimid_all_pass += 1;
          FD_pimid_all_check[i] = true;
        }
      }
      else if(cut_beta_vs_p_pim == true && cut_deltabeta_pim == true && cut_tofmass_pim == true){
        if(FD_pimid_default_PID_check[i] && FD_pimid_charge_check[i] && FD_pimid_ele_reject_check[i] && FD_pimid_EC_outer_vs_EC_inner_check[i] && FD_pimid_DC_hit_position_region1_fiducial_check[i] 
                                         && FD_pimid_DC_hit_position_region2_fiducial_check[i] && FD_pimid_DC_hit_position_region3_fiducial_check[i] && FD_pimid_beta_check[i] && FD_pimid_delta_beta_check[i] 
                                         && FD_pimid_tofmass_check[i]){
          FD_pimid_all_pass += 1;
          FD_pimid_all_check[i] = true;
        }
      }
      else{
        if(FD_pimid_default_PID_check[i] && FD_pimid_charge_check[i] && FD_pimid_ele_reject_check[i] && FD_pimid_EC_outer_vs_EC_inner_check[i] && FD_pimid_DC_hit_position_region1_fiducial_check[i] 
                                         && FD_pimid_DC_hit_position_region2_fiducial_check[i] && FD_pimid_DC_hit_position_region3_fiducial_check[i]){
          FD_pimid_all_pass += 1;
          FD_pimid_all_check[i] = true;
        }
      }
    }

    /// ////////////////////////////////////////////////////////////////////////////////////////////
    /// Central detector:

    if(part_status[i] >= 4000){

      // PID checks:
 
      double beta_charge_central = Beta_charged_central(i, run);

      CD_pimid_default_PID_check[i] = pim_default_PID_cut(i);
      CD_pimid_charge_check[i] = pim_charge_cut(i);
      CD_pimid_beta_check[i] = CD_pim_beta_cut(i, run);
      CD_pimid_maximum_probability_check[i] = CD_maximum_probability_cut(i, -211, 0.27, 99.73, run);
      CD_pimid_delta_vz_check[i] = CD_pim_delta_vz_cut(i);

      if(CD_pimid_default_PID_check[i]) 											CD_pimid_default_PID_pass += 1;
      if(CD_pimid_charge_check[i]) 												CD_pimid_charge_pass += 1;
      if(CD_pimid_beta_check[i] && CD_pimid_charge_check[i] && CD_pimid_default_PID_check[i]) 					CD_pimid_beta_pass += 1;
      if(CD_pimid_maximum_probability_check[i] && CD_pimid_charge_check[i] && CD_pimid_default_PID_check[i]) 			CD_pimid_maximum_probability_pass += 1;
      if(CD_pimid_delta_vz_check[i] && CD_pimid_charge_check[i] && CD_pimid_default_PID_check[i]) 				CD_pimid_delta_vz_pass += 1;

      if(CD_cut_maximum_probability_pim == true){
        if(CD_pimid_charge_check[i] && CD_pimid_maximum_probability_check[i] && (beta_charge_central < 0.9999 || beta_charge_central > 1.0001)){
          CD_pimid_all_pass += 1;
          CD_pimid_all_check[i] = true;
        }
      }
      else if(CD_cut_beta_vs_p_pim == true){
        if(CD_pimid_charge_check[i] && CD_pimid_beta_check[i] && (beta_charge_central < 0.9999 || beta_charge_central > 1.0001)){
          CD_pimid_all_pass += 1;
          CD_pimid_all_check[i] = true;
        }
      }
      else{
        if(CD_pimid_default_PID_check[i] && CD_pimid_charge_check[i]){
          CD_pimid_all_pass += 1;
          CD_pimid_all_check[i] = true;
        }
      }
    }

    /// ///////////////////////////////////////////////////////////////
    /// Create pim selector

    bool selector;

    if(use_FD == false) FD_pimid_all_check[i] = false;
    if(use_CD == false) CD_pimid_all_check[i] = false;


    /// ////////////////////////////////////////////////////////////////
    /// pick particle index and sort by momentum

    if( i < BUFFER){
      mom = 0;
      check = 0;
      for(int k = 0; k < Npart; k++){
        if(k < BUFFER){

        if(use_own_PID_pim == true){ 
          selector = (FD_pimid_all_check[k] && part_status[i] >= 2000 && part_status[i] < 4000) || (CD_pimid_all_check[k] && part_status[i] >= 4000); 
        }
        else{ 
          selector = (FD_pimid_default_PID_check[k] && part_status[i] >= 2000 && part_status[i] < 4000) || (CD_pimid_default_PID_check[k] && part_status[i] >= 4000); 
        }

        if(selector){
          for(int j = 0; j < i; j++){
            if(k == pim_ind[j]) check = -1;
          }
          if(sqrt(vpart_px->at(k)*vpart_px->at(k) + vpart_py->at(k)*vpart_py->at(k) + vpart_pz->at(k)*vpart_pz->at(k)) > mom  && check != -1){        
	    mom = sqrt(vpart_px->at(k)*vpart_px->at(k) + vpart_py->at(k)*vpart_py->at(k) + vpart_pz->at(k)*vpart_pz->at(k)); 
            pim_ind[i] = k;
            pim_count += 1;  
	  }
          check = 0;
}
	}  
      }
    }
  }
}

  /// //////////////////////////////////////////////////////////////////
  /// Assign properties

  for(int i = 0; i < BUFFER; i++){
    if(pim_ind[i] != -1){
      pim_vx[i] = vpart_vx->at(pim_ind[i]);
      pim_vy[i] = vpart_vy->at(pim_ind[i]);
      pim_vz[i] = vpart_vz->at(pim_ind[i]);
      pim_beta[i] = vpart_beta->at(pim_ind[i]);
      double p = sqrt(vpart_px->at(pim_ind[i])*vpart_px->at(pim_ind[i]) + vpart_py->at(pim_ind[i])*vpart_py->at(pim_ind[i]) + vpart_pz->at(pim_ind[i])*vpart_pz->at(pim_ind[i]));
      p4_pim[i].SetPxPyPzE(vpart_px->at(pim_ind[i]), vpart_py->at(pim_ind[i]), vpart_pz->at(pim_ind[i]), sqrt(p*p + m_pim*m_pim));
      pim_FTOF_sec[i] = part_FTOF_sector_layer2[pim_ind[i]];
      if(part_status[pim_ind[i]] >= 2000 && part_status[pim_ind[i]] < 4000) pim_detect[i] = 2; 
      if(part_status[pim_ind[i]] >= 4000)                                   pim_detect[i] = 3; 
    }
  }  

}

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Kp selector
/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void select_Kplus(int run){

  float mom = 0;
  int check = 0;
  int Npart = vpart_pid->size();

  for(Int_t i = 0; i < Npart; i++){
    if( i < BUFFER){

    /// ////////////////////////////////////////////////////////////////////////////////////////////
    /// Forward detector:

    if(part_status[i] >= 2000 && part_status[i] < 4000){

      // PID checks:

      FD_Kpid_default_PID_check[i] = Kp_default_PID_cut(i);
      FD_Kpid_charge_check[i] = Kp_charge_cut(i);
      FD_Kpid_DC_hit_position_region1_fiducial_check[i] = DC_hit_position_region1_fiducial_cut_hadrons_positive(i);
      FD_Kpid_DC_hit_position_region2_fiducial_check[i] = DC_hit_position_region2_fiducial_cut_hadrons_positive(i);
      FD_Kpid_DC_hit_position_region3_fiducial_check[i] = DC_hit_position_region3_fiducial_cut_hadrons_positive(i);
      FD_Kpid_beta_check[i] = Kp_beta_cut(i, run);
      FD_Kpid_delta_beta_check[i] = Kp_delta_beta_cut(i, run);
      FD_Kpid_tofmass_check[i] = Kp_tofmass_cut(i, run);
      FD_Kpid_maximum_probability_check[i] = maximum_probability_cut(i, 321, 0.27, 99.73, run);   // check if the hypothesis proton is fulfilled for the particle - proton: 2212,  pip: 211,  Kp: 321 
										        // conflvl = required probability, that it is the hypothesis particle in % - 3sigma = 0.27, 2sigma = 4.55, 1sigma = 31.73
                                                                                        // anticonflvl = maximum acceptable probability for another candidate in %
      FD_Kpid_delta_vz_check[i] = Kp_delta_vz_cut(i);

      if(FD_Kpid_default_PID_check[i]  == true) 									FD_Kpid_default_PID_pass += 1;
      if(FD_Kpid_charge_check[i]   == true) 										FD_Kpid_charge_pass += 1;
      if(FD_Kpid_DC_hit_position_region1_fiducial_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_default_PID_check[i]) 	FD_Kpid_DC_hit_position_region1_fiducial_pass += 1;
      if(FD_Kpid_DC_hit_position_region2_fiducial_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_default_PID_check[i]) 	FD_Kpid_DC_hit_position_region2_fiducial_pass += 1;
      if(FD_Kpid_DC_hit_position_region3_fiducial_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_default_PID_check[i]) 	FD_Kpid_DC_hit_position_region3_fiducial_pass += 1;
      if(FD_Kpid_beta_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_default_PID_check[i]) 				FD_Kpid_beta_pass += 1;
      if(FD_Kpid_delta_beta_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_default_PID_check[i]) 			FD_Kpid_delta_beta_pass += 1;
      if(FD_Kpid_tofmass_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_default_PID_check[i]) 				FD_Kpid_tofmass_pass += 1;
      if(FD_Kpid_maximum_probability_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_default_PID_check[i]) 		FD_Kpid_maximum_probability_pass += 1;
      if(FD_Kpid_delta_vz_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_default_PID_check[i]) 				FD_Kpid_delta_vz_pass += 1;

      if(cut_maximum_probability_Kp == true){
        if(FD_Kpid_default_PID_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_DC_hit_position_region1_fiducial_check[i] && FD_Kpid_DC_hit_position_region2_fiducial_check[i] 
                                        && FD_Kpid_DC_hit_position_region3_fiducial_check[i] && FD_Kpid_maximum_probability_check[i]){
          FD_Kpid_all_pass += 1;
          FD_Kpid_all_check[i] = true;
        }
      }
      else if(cut_beta_vs_p_Kp == true){
        if(FD_Kpid_default_PID_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_DC_hit_position_region1_fiducial_check[i] && FD_Kpid_DC_hit_position_region2_fiducial_check[i] 
                                        && FD_Kpid_DC_hit_position_region3_fiducial_check[i] && FD_Kpid_beta_check[i]){
          FD_Kpid_all_pass += 1;
          FD_Kpid_all_check[i] = true;
        }
      }
      else if(cut_beta_vs_p_Kp == true && cut_deltabeta_Kp == true && cut_tofmass_Kp == true){
        if(FD_Kpid_default_PID_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_DC_hit_position_region1_fiducial_check[i] && FD_Kpid_DC_hit_position_region2_fiducial_check[i] 
                                        && FD_Kpid_DC_hit_position_region3_fiducial_check[i] && FD_Kpid_beta_check[i] && FD_Kpid_delta_beta_check[i] && FD_Kpid_tofmass_check[i]){
          FD_Kpid_all_pass += 1;
          FD_Kpid_all_check[i] = true;
        }
      }
      else{
        if(FD_Kpid_default_PID_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_DC_hit_position_region1_fiducial_check[i] && FD_Kpid_DC_hit_position_region2_fiducial_check[i] 
                                        && FD_Kpid_DC_hit_position_region3_fiducial_check[i]){
          FD_Kpid_all_pass += 1;
          FD_Kpid_all_check[i] = true;
        }
      }
    }

    /// ////////////////////////////////////////////////////////////////////////////////////////////
    /// Central detector:

    if(part_status[i] >= 4000){

      // PID checks:

      CD_Kpid_default_PID_check[i] = Kp_default_PID_cut(i);
      CD_Kpid_charge_check[i] = Kp_charge_cut(i);
      CD_Kpid_beta_check[i] = CD_Kp_beta_cut(i, run);
      CD_Kpid_maximum_probability_check[i] = CD_maximum_probability_cut(i, 321, 0.27, 99.73, run);
      CD_Kpid_delta_vz_check[i] = CD_Kp_delta_vz_cut(i);

      if(CD_Kpid_default_PID_check[i]) 											CD_Kpid_default_PID_pass += 1;
      if(CD_Kpid_charge_check[i]) 											CD_Kpid_charge_pass += 1;
      if(CD_Kpid_beta_check[i] && CD_Kpid_charge_check[i] && CD_Kpid_default_PID_check[i]) 				CD_Kpid_beta_pass += 1;
      if(CD_Kpid_maximum_probability_check[i] && CD_Kpid_charge_check[i] && CD_Kpid_default_PID_check[i]) 		CD_Kpid_maximum_probability_pass += 1;
      if(CD_Kpid_delta_vz_check[i] && CD_Kpid_charge_check[i] && CD_Kpid_default_PID_check[i]) 				CD_Kpid_delta_vz_pass += 1;

      if(CD_cut_maximum_probability_Kp == true){
        if(CD_Kpid_charge_check[i] && CD_Kpid_maximum_probability_check[i]){
          CD_Kpid_all_pass += 1;
          CD_Kpid_all_check[i] = true;
        }
      }
      else if(CD_cut_beta_vs_p_Kp == true){
        if(CD_Kpid_charge_check[i] && CD_Kpid_beta_check[i]){
          CD_Kpid_all_pass += 1;
          CD_Kpid_all_check[i] = true;
        }
      }
      else{
        if(CD_Kpid_default_PID_check[i] && CD_Kpid_charge_check[i]){
          CD_Kpid_all_pass += 1;
          CD_Kpid_all_check[i] = true;
        }
      }
    }

    /// ///////////////////////////////////////////////////////////////
    /// Create Kp selector

    bool selector;

    if(use_FD == false) FD_Kpid_all_check[i] = false;
    if(use_CD == false) CD_Kpid_all_check[i] = false;


    /// ////////////////////////////////////////////////////////////////
    /// pick particle index and sort by momentum

    if( i < BUFFER){
      mom = 0;
      check = 0;
      for(int k = 0; k < Npart; k++){
        if(k < BUFFER){
        if(use_own_PID_Kp == true){ 
          selector = (FD_Kpid_all_check[k] && part_status[i] >= 2000 && part_status[i] < 4000) || (CD_Kpid_all_check[k] && part_status[i] >= 4000); 
        }
        else{ 
          selector = (FD_Kpid_default_PID_check[k] && part_status[i] >= 2000 && part_status[i] < 4000) || (CD_Kpid_default_PID_check[k] && part_status[i] >= 4000); 
        }

        if(selector){
          for(int j = 0; j < i; j++){
            if(k == Kp_ind[j]) check = -1;
          }
          if(sqrt(vpart_px->at(k)*vpart_px->at(k) + vpart_py->at(k)*vpart_py->at(k) + vpart_pz->at(k)*vpart_pz->at(k)) > mom  && check != -1){        
	    mom = sqrt(vpart_px->at(k)*vpart_px->at(k) + vpart_py->at(k)*vpart_py->at(k) + vpart_pz->at(k)*vpart_pz->at(k)); 
            Kp_ind[i] = k;
            Kp_count += 1;  
	  }
          check = 0;
}
	}  
      }
    }
  }
}

  /// //////////////////////////////////////////////////////////////////
  /// Assign properties

  for(int i = 0; i < BUFFER; i++){
    if(Kp_ind[i] != -1){
      Kp_vx[i] = vpart_vx->at(Kp_ind[i]);
      Kp_vy[i] = vpart_vy->at(Kp_ind[i]);
      Kp_vz[i] = vpart_vz->at(Kp_ind[i]);
      Kp_beta[i] = vpart_beta->at(Kp_ind[i]);
      double p = sqrt(vpart_px->at(Kp_ind[i])*vpart_px->at(Kp_ind[i]) + vpart_py->at(Kp_ind[i])*vpart_py->at(Kp_ind[i]) + vpart_pz->at(Kp_ind[i])*vpart_pz->at(Kp_ind[i]));
      p4_Kp[i].SetPxPyPzE(vpart_px->at(Kp_ind[i]), vpart_py->at(Kp_ind[i]), vpart_pz->at(Kp_ind[i]), sqrt(p*p + m_Kp*m_Kp));
      Kp_FTOF_sec[i] = part_FTOF_sector_layer2[Kp_ind[i]];
      if(part_status[Kp_ind[i]] >= 2000 && part_status[Kp_ind[i]] < 4000) Kp_detect[i] = 2; 
      if(part_status[Kp_ind[i]] >= 4000)                                  Kp_detect[i] = 3; 
    }
  }  


}

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Km selector
/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void select_Kminus(int run){

  float mom = 0;
  int check = 0;
  int Npart = vpart_pid->size();


  for(Int_t i = 0; i < Npart; i++){
    if( i < BUFFER){

    /// ////////////////////////////////////////////////////////////////////////////////////////////
    /// Forward detector:

    if(part_status[i] >= 2000 && part_status[i] < 4000){

      // PID checks:

      FD_Kmid_default_PID_check[i] = Km_default_PID_cut(i);
      FD_Kmid_charge_check[i] = Km_charge_cut(i);
      FD_Kmid_ele_reject_check[i] = Km_ele_reject_cut(i);
      FD_Kmid_EC_outer_vs_EC_inner_check[i] = Km_EC_outer_vs_EC_inner_cut(i);
      FD_Kmid_DC_hit_position_region1_fiducial_check[i] = DC_hit_position_region1_fiducial_cut_hadrons_negative(i);
      FD_Kmid_DC_hit_position_region2_fiducial_check[i] = DC_hit_position_region2_fiducial_cut_hadrons_negative(i);
      FD_Kmid_DC_hit_position_region3_fiducial_check[i] = DC_hit_position_region3_fiducial_cut_hadrons_negative(i);
      FD_Kmid_beta_check[i] = Km_beta_cut(i, run);
      FD_Kmid_delta_beta_check[i] = Km_delta_beta_cut(i, run);
      FD_Kmid_tofmass_check[i] = Km_tofmass_cut(i, run);
      FD_Kmid_maximum_probability_check[i] = maximum_probability_cut(i, -321, 0.27, 99.73, run);  // check if the hypothesis proton is fulfilled for the particle - proton: 2212,  pip: 211,  Kp: 321 
										        // conflvl = required probability, that it is the hypothesis particle in % - 3sigma = 0.27, 2sigma = 4.55, 1sigma = 31.73
                                                                                        // anticonflvl = maximum acceptable probability for another candidate in %
      FD_Kmid_delta_vz_check[i] = Km_delta_vz_cut(i);

      if(FD_Kmid_default_PID_check[i] == true) 										FD_Kmid_default_PID_pass += 1;
      if(FD_Kmid_charge_check[i] == true) 										FD_Kmid_charge_pass += 1;
      if(FD_Kmid_DC_hit_position_region1_fiducial_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_default_PID_check[i]) 	FD_Kmid_DC_hit_position_region1_fiducial_pass += 1;
      if(FD_Kmid_DC_hit_position_region2_fiducial_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_default_PID_check[i]) 	FD_Kmid_DC_hit_position_region2_fiducial_pass += 1;
      if(FD_Kmid_DC_hit_position_region3_fiducial_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_default_PID_check[i]) 	FD_Kmid_DC_hit_position_region3_fiducial_pass += 1;
      if(FD_Kmid_beta_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_default_PID_check[i]) 				FD_Kmid_beta_pass += 1;
      if(FD_Kmid_delta_beta_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_default_PID_check[i]) 			FD_Kmid_delta_beta_pass += 1;
      if(FD_Kmid_tofmass_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_default_PID_check[i]) 				FD_Kmid_tofmass_pass += 1;
      if(FD_Kmid_maximum_probability_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_default_PID_check[i]) 		FD_Kmid_maximum_probability_pass += 1;
      if(FD_Kmid_delta_vz_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_default_PID_check[i]) 				FD_Kmid_delta_vz_pass += 1;

      if(cut_maximum_probability_Km == true){
        if(FD_Kmid_default_PID_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_ele_reject_check[i] && FD_Kmid_DC_hit_position_region1_fiducial_check[i] && FD_Kmid_DC_hit_position_region2_fiducial_check[i] 
                                        && FD_Kmid_DC_hit_position_region3_fiducial_check[i] && FD_Kmid_maximum_probability_check[i]){
          FD_Kmid_all_pass += 1;
          FD_Kmid_all_check[i] = true;
        }
      }
      else if(cut_beta_vs_p_Km == true){
        if(FD_Kmid_default_PID_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_ele_reject_check[i] && FD_Kmid_DC_hit_position_region1_fiducial_check[i] && FD_Kmid_DC_hit_position_region2_fiducial_check[i] 
                                        && FD_Kmid_DC_hit_position_region3_fiducial_check[i] && FD_Kmid_beta_check[i]){
          FD_Kmid_all_pass += 1;
          FD_Kmid_all_check[i] = true;
        }
      }
      else if(cut_beta_vs_p_Km == true && cut_deltabeta_Km == true && cut_tofmass_Km == true){
        if(FD_Kmid_default_PID_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_ele_reject_check[i] && FD_Kmid_DC_hit_position_region1_fiducial_check[i] && FD_Kmid_DC_hit_position_region2_fiducial_check[i] 
                                        && FD_Kmid_DC_hit_position_region3_fiducial_check[i] && FD_Kmid_beta_check[i] && FD_Kmid_delta_beta_check[i] && FD_Kmid_tofmass_check[i]){
          FD_Kmid_all_pass += 1;
          FD_Kmid_all_check[i] = true;
        }
      }
      else{
        if(FD_Kmid_default_PID_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_ele_reject_check[i] && FD_Kmid_DC_hit_position_region1_fiducial_check[i] && FD_Kmid_DC_hit_position_region2_fiducial_check[i]
                                        && FD_Kmid_DC_hit_position_region3_fiducial_check[i]){
          FD_Kmid_all_pass += 1;
          FD_Kmid_all_check[i] = true;
        }
      }
    }

    /// ////////////////////////////////////////////////////////////////////////////////////////////
    /// Central detector:

    if(part_status[i] >= 4000){

      // PID checks:

      double beta_charge_central = Beta_charged_central(i, run);

      CD_Kmid_default_PID_check[i] = Km_default_PID_cut(i);
      CD_Kmid_charge_check[i] = Km_charge_cut(i);
      CD_Kmid_beta_check[i] = CD_Km_beta_cut(i, run);
      CD_Kmid_maximum_probability_check[i] = CD_maximum_probability_cut(i, -321, 0.27, 99.73, run);
      CD_Kmid_delta_vz_check[i] = CD_Km_delta_vz_cut(i);

      if(CD_Kmid_default_PID_check[i]) 											CD_Kmid_default_PID_pass += 1;
      if(CD_Kmid_charge_check[i]) 											CD_Kmid_charge_pass += 1;
      if(CD_Kmid_beta_check[i] && CD_Kmid_charge_check[i] && CD_Kmid_default_PID_check[i]) 				CD_Kmid_beta_pass += 1;
      if(CD_Kmid_maximum_probability_check[i] && CD_Kmid_charge_check[i] && CD_Kmid_default_PID_check[i]) 		CD_Kmid_maximum_probability_pass += 1;
      if(CD_Kmid_delta_vz_check[i] && CD_Kmid_charge_check[i] && CD_Kmid_default_PID_check[i]) 				CD_Kmid_delta_vz_pass += 1;

      if(CD_cut_maximum_probability_Km == true){
        if(CD_Kmid_charge_check[i] && CD_Kmid_maximum_probability_check[i] && (beta_charge_central < 0.9999 || beta_charge_central > 1.0001)){
          CD_Kmid_all_pass += 1;
          CD_Kmid_all_check[i] = true;
        }
      }
      else if(CD_cut_beta_vs_p_Km == true){
        if(CD_Kmid_charge_check[i] && CD_Kmid_beta_check[i] && (beta_charge_central < 0.9999 || beta_charge_central > 1.0001)){
          CD_Kmid_all_pass += 1;
          CD_Kmid_all_check[i] = true;
        }
      }
      else{
        if(CD_Kmid_default_PID_check[i] && CD_Kmid_charge_check[i]){
          CD_Kmid_all_pass += 1;
          CD_Kmid_all_check[i] = true;
        }
      }
    }

    /// ///////////////////////////////////////////////////////////////
    /// Create Km selector

    bool selector;

    if(use_FD == false) FD_Kmid_all_check[i] = false;
    if(use_CD == false) CD_Kmid_all_check[i] = false;


    /// ////////////////////////////////////////////////////////////////
    /// pick particle index and sort by momentum

    if( i < BUFFER){
      mom = 0;
      check = 0;
      for(int k = 0; k < Npart; k++){
        if(k < BUFFER){
        if(use_own_PID_Km == true){ 
          selector = (FD_Kmid_all_check[k] && part_status[i] >= 2000 && part_status[i] < 4000) || (CD_Kmid_all_check[k] && part_status[i] >= 4000); 
        }
        else{ 
          selector = (FD_Kmid_default_PID_check[k] && part_status[i] >= 2000 && part_status[i] < 4000) || (CD_Kmid_default_PID_check[k] && part_status[i] >= 4000); 
        }

        if(selector){
          for(int j = 0; j < i; j++){
            if(k == Km_ind[j]) check = -1;
          }
          if(sqrt(vpart_px->at(k)*vpart_px->at(k) + vpart_py->at(k)*vpart_py->at(k) + vpart_pz->at(k)*vpart_pz->at(k)) > mom  && check != -1){        
	    mom = sqrt(vpart_px->at(k)*vpart_px->at(k) + vpart_py->at(k)*vpart_py->at(k) + vpart_pz->at(k)*vpart_pz->at(k)); 
            Km_ind[i] = k;
            Km_count += 1;  
	  }
          check = 0;
}
	}  
      }
    }
  }
}

  /// //////////////////////////////////////////////////////////////////
  /// Assign properties

  for(int i = 0; i < BUFFER; i++){
    if(Km_ind[i] != -1){
      Km_vx[i] = vpart_vx->at(Km_ind[i]);
      Km_vy[i] = vpart_vy->at(Km_ind[i]);
      Km_vz[i] = vpart_vz->at(Km_ind[i]);
      Km_beta[i] = vpart_beta->at(Km_ind[i]);
      double p = sqrt(vpart_px->at(Km_ind[i])*vpart_px->at(Km_ind[i]) + vpart_py->at(Km_ind[i])*vpart_py->at(Km_ind[i]) + vpart_pz->at(Km_ind[i])*vpart_pz->at(Km_ind[i]));
      p4_Km[i].SetPxPyPzE(vpart_px->at(Km_ind[i]), vpart_py->at(Km_ind[i]), vpart_pz->at(Km_ind[i]), sqrt(p*p + m_Km*m_Km));
      Km_FTOF_sec[i] = part_FTOF_sector_layer2[Km_ind[i]];
      if(part_status[Km_ind[i]] >= 2000 && part_status[Km_ind[i]] < 4000) Km_detect[i] = 2; 
      if(part_status[Km_ind[i]] >= 4000)                                  Km_detect[i] = 3; 
    }
  }  

}

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// photon selector
/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void select_photon(int run){

  float mom = 0;
  int check = 0;
  int Npart = vpart_pid->size();


  for(Int_t i = 0; i < Npart; i++){
    if( i < BUFFER){

    /// ////////////////////////////////////////////////////////////
    /// Forward detector cuts:

    if(part_status[i] >= 2000 && part_status[i] < 4000){

      FD_photid_default_PID_check[i] = phot_default_PID_cut(i);
      FD_photid_charge_check[i] = phot_charge_cut(i);
      FD_photid_beta_check[i] = phot_beta_cut(i, run);
      FD_photid_EC_sampling_fraction_check[i] = phot_EC_sampling_fraction_cut(i);
      FD_photid_EC_outer_vs_EC_inner_check[i] = phot_EC_outer_vs_EC_inner_cut(i);
      FD_photid_EC_hit_position_fiducial_check[i] = phot_EC_hit_position_fiducial_cut(i);

      if(FD_photid_default_PID_check[i]) 										FD_photid_default_PID_pass += 1;
      if(FD_photid_charge_check[i]) 											FD_photid_charge_pass += 1;
      if(FD_photid_beta_check[i] && FD_photid_charge_check[i] && FD_photid_default_PID_check[i]) 			FD_photid_beta_pass += 1;
      if(FD_photid_EC_hit_position_fiducial_check[i] && FD_photid_charge_check[i] && FD_photid_default_PID_check[i]) 	FD_photid_EC_hit_position_fiducial_pass += 1;

      if(FD_photid_default_PID_check[i] && FD_photid_default_PID_check[i] && FD_photid_charge_check[i] && FD_photid_beta_check[i] && FD_photid_EC_hit_position_fiducial_check[i]){
        FD_photid_all_pass += 1;
        FD_photid_all_check[i] = true;
      }
    }

    /// /////////////////////////////////////////////////////////////
    /// Forward tagger cuts:

    if(part_status[i] >= 1000 && part_status[i] < 2000){

      FT_photid_PID_check[i] = FT_photid_PID_cut(i);
      FT_photid_charge_check[i] = FT_photid_charge_cut(i);
      FT_photid_FTCAL_fiducial_check[i] = FT_photid_FTCAL_fiducial_cut(i);
      FT_photid_beta_check[i] = FT_photid_beta_cut(i, run);

      if(FT_photid_PID_check[i])    									FT_photid_PID_pass += 1;
      if(FT_photid_charge_check[i]) 									FT_photid_charge_pass += 1;
      if(FT_photid_FTCAL_fiducial_check[i] && FT_photid_charge_check[i] && FT_photid_PID_check[i])	FT_photid_FTCAL_fiducial_pass += 1;
      if(FT_photid_beta_check[i] && FT_photid_charge_check[i] && FT_photid_PID_check[i]) 		FT_photid_beta_pass += 1;

      if(FT_photid_PID_check[i] && FT_photid_charge_check[i] && FT_photid_FTCAL_fiducial_check[i] && FT_photid_beta_check[i]){
        FT_photid_all_pass += 1;
        FT_photid_all_check[i] = true;
      }
    }

    /// ///////////////////////////////////////////////////////////////
    /// Create photon selector

    bool selector;
    double theta_min = 2.5;

    if(use_FT == false){ FT_photid_all_check[i] = false;}
    if(use_FD == false){ FD_photid_all_check[i] = false;}


    /// ///////////////////////////////////////////////////////////////
    /// Pick photon particle index and sort by momentum
 
    if( i < BUFFER){
      mom = 0;
      check = 0;
      for(int k = 0; k < Npart; k++){
        if(k < BUFFER){
        if(use_own_PID_photon == true){ 
          selector = (FD_photid_all_check[k] && part_status[i] >= 2000 && part_status[i] < 4000)|| (FT_photid_all_check[k] && part_status[i] >= 1000 && part_status[i] < 2000 && part_p[k] > 0.1); 
        }
        else{ 
          selector = (FD_photid_default_PID_check[k] && part_status[i] >= 2000 && part_status[i] < 4000) || (FT_photid_PID_check[k]  && part_status[i] >= 1000 && part_status[i] < 2000); 
        }       

        if(selector){   // photons only properly detected in FT and FD
          for(int j = 0; j < i; j++){
            if(k == g_ind[j]) check = -1;
          }

          //double ectot = part_Cal_PCAL_energy[k] + part_Cal_ECin_energy[k] + part_Cal_ECout_energy[k];
          //double p_corr = ectot/(0.25*(1.029-(0.015/ectot) + (0.00012/pow(ectot,2))));
          //if(p_corr > mom  && check != -1){  
	    //mom = p_corr;

          if(sqrt(vpart_px->at(k)*vpart_px->at(k) + vpart_py->at(k)*vpart_py->at(k) + vpart_pz->at(k)*vpart_pz->at(k)) > mom  && check != -1){  
	    mom = sqrt(vpart_px->at(k)*vpart_px->at(k) + vpart_py->at(k)*vpart_py->at(k) + vpart_pz->at(k)*vpart_pz->at(k)); 
            g_ind[i] = k;
            g_count += 1;  
	  }
          check = 0;
          }
	}  
      }
    }
  }
}

  /// //////////////////////////////////////////////////////////////////
  /// Assign properties

  for(int i = 0; i < BUFFER; i++){
    if(g_ind[i] != -1){
      g_vx[i] = vpart_vx->at(g_ind[i]);
      g_vy[i] = vpart_vy->at(g_ind[i]);
      g_vz[i] = vpart_vz->at(g_ind[i]);

      //double ectot = part_Cal_PCAL_energy[g_ind[i]] + part_Cal_ECin_energy[g_ind[i]] + part_Cal_ECout_energy[g_ind[i]];
      //double p_corr = ectot/(0.25*(1.029-(0.015/ectot) + (0.00012/pow(ectot,2))));
      double p = sqrt(vpart_px->at(g_ind[i])*vpart_px->at(g_ind[i]) + vpart_py->at(g_ind[i])*vpart_py->at(g_ind[i]) + vpart_pz->at(g_ind[i])*vpart_pz->at(g_ind[i]));
      //p4_phot[i].SetPxPyPzE(p_corr * vpart_px->at(g_ind[i])/p, p_corr * vpart_py->at(g_ind[i])/p, p_corr * vpart_pz->at(g_ind[i])/p, p_corr);

      p4_phot[i].SetPxPyPzE(vpart_px->at(g_ind[i]), vpart_py->at(g_ind[i]), vpart_pz->at(g_ind[i]), p);

      g_sec[i] = part_Cal_PCAL_sector[g_ind[i]];
      if(part_status[g_ind[i]] >= 1000 && part_status[g_ind[i]] < 2000) phot_detect[i] = 1; 
      if(part_status[g_ind[i]] >= 2000 && part_status[g_ind[i]] < 4000) phot_detect[i] = 2; 
    }
  }  

}

/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// particle ID cuts for the FD:
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// hit in the FTOF is required for electrons and charged hadrons (basic cut)  --> For now: Selection of good paddels included
///

bool basic_FTOF_cut(int j){

  bool paddel_select = false;
  //if(part_FTOF_component_layer2[j] > 18 && part_FTOF_component_layer2[j] < 55) paddel_select = true;
  paddel_select = true;   // use all paddels

  if((part_FTOF_sector_layer1[j] != 0 || part_FTOF_sector_layer2[j] != 0 || part_FTOF_sector_layer3[j] != 0) && paddel_select) return true;
  else return false;

}


/// ///////////////////////////////////////////////////////////////////////////
/// FD electrons:

// electron has to have a hit in the FTOF, since this determines the start time of the event.  


bool ele_default_PID_cut(int j){
  if(vpart_pid->at(j) == 11) return true;
  else return false;
}

bool ele_charge_cut(int j){
  if(part_charge[j] == -1) return true;
  else return false;
}


bool CC_nphe_cut(int j){

  double nphe_min = 2;

  if(part_CC_HTCC_nphe[j] > nphe_min) return true;
  else return false;
}


// EC cuts

bool EC_outer_vs_EC_inner_cut(int j){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double edep_min;

  if(loose  == true){ edep_min = 0.06;}
  if(medium == true){ edep_min = 0.06;}
  if(tight  == true){ edep_min = 0.1;}

  if(part_Cal_PCAL_energy[j] > edep_min) return true; 
  else return false;

}


bool EC_sampling_fraction_cut(int j){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double p_min = 0.2381 + 0.11905 * Ebeam;  //

  double sigma_range = 3;

  if(tight == true){
    sigma_range = 2;
  }

  if(medium == true){
    sigma_range = 3;
  }

  if(loose == true){
    sigma_range = 4;
  }


  /// //////////////////////////////////////////////////////////////////////////////////////////////////////
  /// a) cut based on sampling fraction versus drift chamber momentum

  // 10.6 GeV (fall inbending)

  double p0mean[] = {0.105631, 0.11551, 0.112799, 0.109937, 0.116249, 0.119057};
  double p1mean[] = {-0.153951, -0.0253273, -0.125718, 0.165414, 0.0768411, 0.0555026};
  double p2mean[] = {0.00860091, 0.00706291, 0.00908884, 0.00499666, 0.00448701, 0.00558927};
  double p3mean[] = {-0.000821675, -0.000711488, -0.000930922, -0.000298311, -0.000455716, -0.000657084};
  double p0sigma[] = {0.0149613, 0.0115116, 0.00580737, 0.0106817, 0.012667, 0.00553471};
  double p1sigma[] = {0.00700773, 0.0116193, 0.0202375, 0.0126958, 0.00892239, 0.0216206};


  if(outbending == true){    // 10.6 GeV (fall outbending)
  
    double p0mean_out[] = {0.105467, 0.115261, 0.127793, 0.113359, 0.112263, 0.113507};
    double p1mean_out[] = {-0.135178, 0.135808, 0.903412, 0.598274, -0.0466815, 0.0550123};
    double p2mean_out[] = {0.00996842, 0.00672508, -0.00035721, 0.00470925, 0.00588451, 0.00923385};
    double p3mean_out[] = {-0.000754536, -0.000515365, -0.000108273, -0.000447278, -0.000358148, -0.00074643};
    double p0sigma_out[] = {0.00683747, 0.0065199, 0.00297734, 0.00759701, 0.0093309, 0.00591988};
    double p1sigma_out[] = {0.0180228, 0.0183979, 0.0250332, 0.0155001, 0.0137594, 0.0215643};
  
    for(Int_t i = 0; i < 6; i++){  
      p0mean[i] = p0mean_out[i];
      p1mean[i] = p1mean_out[i];
      p2mean[i] = p2mean_out[i];
      p3mean[i] = p3mean_out[i];
      p0sigma[i] = p0sigma_out[i];
      p1sigma[i] = p1sigma_out[i];
    }
  }

  double mean = 0;
  double sigma = 0;
  double upper_lim_total = 0;
  double lower_lim_total = 0;

  for(Int_t k = 0; k < 6; k++){  
    if(part_Cal_PCAL_sector[j]-1 == k){
      mean = p0mean[k] *( 1 + part_p[j]/sqrt(pow(part_p[j],2) + p1mean[k])) + p2mean[k] * part_p[j] + p3mean[k] * pow(part_p[j],2);
      sigma = p0sigma[k] + p1sigma[k] / sqrt(part_p[j]);
      upper_lim_total = mean + sigma_range * sigma;
      lower_lim_total = mean - sigma_range * sigma;
    }
  }

  /// /////////////////////////////////////////////////////////////////////////////////////////////////////
  /// b) cut based on sampling fraction versus enrgy deposited in the calorimeter (outdated)
  /*
  // 10.6 GeV
  double p0mean[] = {0.257819, 0.262394, 0.264319, 0.273215, 0.277912, 0.269706};
  double p1mean[] = {0.965033, 0.983467, 0.989758, 1.01724, 1.03319, 1.00817};
  double p2mean[] = {-0.0632092, -0.0290537, -0.0332319, -0.0940031, -0.103213, -0.0420508};
  double p3mean[] = {0.00797259, 0.000159148, 0.00050095, 0.0151809, 0.0168369, 0.000980376};
  double p0sigma[] = {0.0210694, 0.0209485, 0.0137895, 0.016642, 0.017343, 0.0162234};
  double p1sigma[] = {-0.00137114, -0.000904917, 0.00362368, 0.00243129, 0.00180574, 0.0026879};

  double mean = 0;
  double sigma = 0;
  double upper_lim_total = 0;
  double lower_lim_total = 0;

  for(Int_t k = 0; k < 6; k++){  
    if(part_Cal_PCAL_sector[j]-1 == k){
      mean = p0mean[k] *( p1mean[k] + p2mean[k] / part_Cal_energy_total[j] + p3mean[k] / pow(part_Cal_energy_total[j],2));
      sigma = p0sigma[k] + p1sigma[k] / sqrt(part_p[j]);
      upper_lim_total = mean + sigma_range * sigma;
      lower_lim_total = mean - sigma_range * sigma;
    }
  }
  */
  /// ////////////////////////////////////////////////////////////////////////////////////////////////////

  if(part_Cal_energy_total[j]/part_p[j] <= upper_lim_total && part_Cal_energy_total[j]/part_p[j] >= lower_lim_total && part_p[j] > p_min) return true;
  else return false;

}



bool EC_hit_position_fiducial_cut(int j){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

// Cut using the natural directions of the scintillator bars/ fibers:

  double u = part_Cal_PCAL_lu[j];
  double v = part_Cal_PCAL_lv[j];
  double w = part_Cal_PCAL_lw[j];
   
  /// v + w is going from teh side to the back end of the PCAL, u is going from side to side
  /// 1 scintillator bar is 4.5 cm wide. In the outer regions (back) double bars are used.

  double min_u = 0;
  double max_u = 420;    
  double min_v = 0;
  double max_v = 420;
  double min_w = 0;
  double max_w = 420;

  if(tight == true){
    min_u = 9.0;     // 2nd double bar cut completely + a bit more  20
    max_u = 420;     // no cut on outer side of u (done by v and w cuts)
    min_v = 9;       // 2 bars
    max_v = 411.5;   // 2nd bar cut completely
    min_w = 9;       // 2 bars
    max_w = 411.5;   // 2nd bar cut completely
  }

  if(medium == true){
    min_u = 6.75;    // center of 2nd double bar +  a bit more  6.75
    max_u = 420;     // no cut on outer side of u (done by v and w cuts)
    min_v = 6.75;    // 1.5 bars
    max_v = 413.75;  // center of 2nd bar
    min_w = 6.75;    // 1.5 bars
    max_w = 413.75;  // center of 2nd bar
  }

  if(loose == true){
    min_u = 4.5;    // outer double bar cut + a bit more  10
    max_u = 420;    // no cut on outer side of u (done by v and w cuts)
    min_v = 4.5;    // 1 bar
    max_v = 416;    // first bar cut completely
    min_w = 4.5;    // 1 bar
    max_w = 416;    // first bar cut completely
  }

  if(u > min_u && u < max_u && v > min_v && v < max_v && w > min_w && w < max_w) return true;
  else return false;

}


// DC fiducial cuts for the 3 regions

bool DC_hit_position_region1_fiducial_cut(int j){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double add = 0;  // value in cm added to the height and radius of the cut
  if(tight == true){  add = 1.0; }
  if(medium == true){ add = 0.0; }
  if(loose == true){  add = -1.0; }

  double angle = 60; 
  double height = 29;
  double radius = 29;

  double height_inb[]  = {27, 22, 22, 27, 22, 22};
  double height_outb[] = {15, 15, 15, 15, 15, 15};

  int sec = part_DC_sector[j]-1;   

  for(Int_t k = 0; k < 6; k++){  
    if(sec == k && inbending == true){
      height = height_inb[k] + add;
      radius = 29 + add;
    }
    if(sec == k && outbending == true){
      height = height_outb[k] + add;
      radius = 25 + add;
    }
  }

  double x1_rot = part_DC_c1y[j] * sin(sec*60.0*Pival/180) + part_DC_c1x[j] * cos(sec*60.0*Pival/180);
  double y1_rot = part_DC_c1y[j] * cos(sec*60.0*Pival/180) - part_DC_c1x[j] * sin(sec*60.0*Pival/180);

  double slope = 1/tan(0.5*angle*Pival/180);
  double left  = (height - slope * y1_rot);
  double right = (height + slope * y1_rot);

  double radius2_DCr1 = pow(radius,2)-pow(y1_rot,2);    // cut out the inner circle

  if (x1_rot > left && x1_rot > right && pow(x1_rot,2) > radius2_DCr1) return true;
  else return false;

}


bool DC_hit_position_region2_fiducial_cut(int j){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double add = 0;  // value in cm added to the height and radius of the cut
  if(tight == true){  add = 2.0; }
  if(medium == true){ add = 0.0; }
  if(loose == true){  add = -2.0; }

  double angle = 60; 
  double height = 35;
  double radius = 40;

  double height_inb[]  = {40, 34, 34, 40, 34, 34};
  double height_outb[] = {25, 25, 25, 25, 25, 25};

  int sec = part_DC_sector[j]-1;

  for(Int_t k = 0; k < 6; k++){  
    if(sec == k && inbending == true){
      height = height_inb[k] + add;
      radius = 38 + add;
    }
    if(sec == k && outbending == true){
      height = height_outb[k] + add;
      radius = 39 + add;
    }
  }
    
  double x2_rot = part_DC_c2y[j] * sin(sec*60.0*Pival/180) + part_DC_c2x[j] * cos(sec*60.0*Pival/180);
  double y2_rot = part_DC_c2y[j] * cos(sec*60.0*Pival/180) - part_DC_c2x[j] * sin(sec*60.0*Pival/180);

  double slope = 1/tan(0.5*angle*Pival/180);
  double left  = (height - slope * y2_rot);
  double right = (height + slope * y2_rot);

  double radius2_DCr2 = pow(radius,2)-pow(y2_rot,2);    // cut out the inner circle

  //return true;
  if (x2_rot > left && x2_rot > right && pow(x2_rot,2) > radius2_DCr2) return true;
  else return false;

}


bool DC_hit_position_region3_fiducial_cut(int j){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double add = 0;  // value in cm added to the height and radius of the cut
  if(tight == true){  add = 3.0; }
  if(medium == true){ add = 0.0; }
  if(loose == true){  add = -3.0; }

  double angle = 60; 
  double height = 48;
  double radius = 49;

  double height_inb[]  = {47, 39, 39, 47, 39, 39};
  double height_outb[] = {48, 48, 48, 48, 48, 48};

  int sec = part_DC_sector[j]-1;

  for(Int_t k = 0; k < 6; k++){  
    if(sec == k && inbending == true){
      height = height_inb[k] + add;
      radius = 50 + add;
    }
    if(sec == k && outbending == true){
      height = height_outb[k] + add;
      radius = 65 + add;
    }
  }

  double x3_rot = part_DC_c3y[j] * sin(sec*60.0*Pival/180) + part_DC_c3x[j] * cos(sec*60.0*Pival/180);
  double y3_rot = part_DC_c3y[j] * cos(sec*60.0*Pival/180) - part_DC_c3x[j] * sin(sec*60.0*Pival/180);

  double slope = 1/tan(0.5*angle*Pival/180);
  double left  = (height - slope * y3_rot);
  double right = (height + slope * y3_rot);

  double radius2_DCr3 = pow(radius,2)-pow(y3_rot,2);    // cut out the inner circle

  if (x3_rot > left && x3_rot > right && pow(x3_rot,2) > radius2_DCr3) return true;
  else return false;

}


bool DC_z_vertex_cut(int j){
 
  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double add = 0;  // value in cm added to the height and radius of the cut
  if(tight == true){  add = 2.0; }
  if(medium == true){ add = 0.0; }
  if(loose == true){  add = -2.0; }

  double vz_min_sect_inb[] = {-12, -12, -12, -12, -12, -12};
  double vz_max_sect_inb[] = {9, 9, 9, 9, 9, 9};

  double vz_min_sect_outb[] = {-14, -14, -14, -14, -14, -14};
  double vz_max_sect_outb[] = {5, 5, 5, 5, 5, 5};

  double vz_min_sect[6];
  double vz_max_sect[6];

  for(Int_t i = 0; i < 6; i++){
    if(inbending == true){   
      vz_min_sect[i] = vz_min_sect_inb[i] + add;
      vz_max_sect[i] = vz_max_sect_inb[i] + add;
    }
    if(outbending == true){   
      vz_min_sect[i] = vz_min_sect_outb[i] + add;
      vz_max_sect[i] = vz_max_sect_outb[i] + add;
    }
  }

  double vz_min = 0;
  double vz_max = 0;

  for(Int_t k = 0; k < 6; k++){  
    if(part_Cal_PCAL_sector[j]-1 == k){
      vz_min = vz_min_sect[k];
      vz_max = vz_max_sect[k];
    }
  }

  if(part_vz[j] > vz_min && part_vz[j] < vz_max) return true;
  else return false;
}



bool Track_Quality_cut(int j){
  if(part_CC_HTCC_sector[j] > 0 && part_DC_sector[j] > 0 && part_Cal_PCAL_sector[j] > 0) return true;
  else return false;
}


/// ///////////////////////////////////////////////////////////////////
/// FD hadrons:

// a) default:

bool prot_default_PID_cut(int j){
  if(vpart_pid->at(j) == 2212) return true;
  else return false;
}
bool neutr_default_PID_cut(int j){
  if(vpart_pid->at(j) == 2112) return true;
  else return false;
}
bool pip_default_PID_cut(int j){
  if(vpart_pid->at(j) == 211) return true;
  else return false;
}
bool pim_default_PID_cut(int j){
  if(vpart_pid->at(j) == -211) return true;
  else return false;
}
bool Kp_default_PID_cut(int j){
  if(vpart_pid->at(j) == 321) return true;
  else return false;
}
bool Km_default_PID_cut(int j){
  if(vpart_pid->at(j) == -321) return true;
  else return false;
}


// b) charge cuts

bool prot_charge_cut(int j){
  if(part_charge[j] == +1) return true;
  else return false;
}
bool neutr_charge_cut(int j){
  if(part_charge[j] == 0) return true;
  else return false;
}
bool pip_charge_cut(int j){
  if(part_charge[j] == +1) return true;
  else return false;
}
bool pim_charge_cut(int j){
  if(part_charge[j] == -1) return true;
  else return false;
}
bool Kp_charge_cut(int j){
  if(part_charge[j] == +1) return true;
  else return false;
}
bool Km_charge_cut(int j){
  if(part_charge[j] == -1) return true;
  else return false;
}


// b.1) electron rejection cut for negative hadrons

bool pim_ele_reject_cut(int j){
  if(FD_eid_all_check[j] == false) return true;
  else return false;
}
bool Km_ele_reject_cut(int j){
  if(FD_eid_all_check[j] == false) return true;
  else return false;
}

bool pim_EC_outer_vs_EC_inner_cut(int j){
  double edep_max = 0.06;
  if(part_Cal_PCAL_energy[j]  < edep_max) return true; 
  else return false;
}

bool Km_EC_outer_vs_EC_inner_cut(int j){
  double edep_max = 0.06;
  if(part_Cal_PCAL_energy[j] < edep_max) return true; 
  else return false;
}


// DC cuts


bool DC_hit_position_region1_fiducial_cut_hadrons_positive(int j){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double add = 0;  // value in cm added to the height and radius of the cut
  if(tight == true){  add = 1.0; }
  if(medium == true){ add = 0.0; }
  if(loose == true){  add = -1.0; }

  double angle = 60; 
  double height = 0; 
  double radius = 0;

  double height_inb[]  = {16, 16, 16, 16, 16, 16};
  double height_outb[] = {20, 20, 20, 20, 20, 20};

  int sec = part_DC_sector[j]-1;   

  for(Int_t k = 0; k < 6; k++){  
    if(sec == k && inbending == true){
      height = height_inb[k] + add;
      radius = 25 + add;
    }
    if(sec == k && outbending == true){
      height = height_outb[k] + add;
      radius = 30 + add;
    }
  }
    
  double x1_rot = part_DC_c1y[j] * sin(sec*60.0*Pival/180) + part_DC_c1x[j] * cos(sec*60.0*Pival/180);
  double y1_rot = part_DC_c1y[j] * cos(sec*60.0*Pival/180) - part_DC_c1x[j] * sin(sec*60.0*Pival/180);

  double slope = 1/tan(0.5*angle*Pival/180);
  double left  = (height - slope * y1_rot);
  double right = (height + slope * y1_rot);

  double radius2_DCr1_min = pow(radius,2)-pow(y1_rot,2);    // cut out the inner circle
  double radius2_DCr1_max = pow(155,2)-pow(y1_rot,2);       // cut out the outer circle

  if (x1_rot > left && x1_rot > right && pow(x1_rot,2) > radius2_DCr1_min && pow(x1_rot,2) < radius2_DCr1_max) return true;
  else return false;

}

bool DC_hit_position_region2_fiducial_cut_hadrons_positive(int j){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double add = 0;  // value in cm added to the height and radius of the cut
  if(tight == true){  add = 2.0; }
  if(medium == true){ add = 0.0; }
  if(loose == true){  add = -2.0; }

  double angle = 60;
  double height = 0; 
  double radius = 0;

  double height_inb[]  = {30, 30, 30, 30, 30, 30};
  double height_outb[] = {31, 31, 31, 31, 31, 31};

  int sec = part_DC_sector[j]-1;   

  for(Int_t k = 0; k < 6; k++){  
    if(sec == k && inbending == true){
      height = height_inb[k] + add;
      radius = 41 + add;
    }
    if(sec == k && outbending == true){
      height = height_outb[k] + add;
      radius = 50 + add;
    }
  }
    
  double x2_rot = part_DC_c2y[j] * sin(sec*60.0*Pival/180) + part_DC_c2x[j] * cos(sec*60.0*Pival/180);
  double y2_rot = part_DC_c2y[j] * cos(sec*60.0*Pival/180) - part_DC_c2x[j] * sin(sec*60.0*Pival/180);

  double slope = 1/tan(0.5*angle*Pival/180);
  double left  = (height - slope * y2_rot);
  double right = (height + slope * y2_rot);

  double radius2_DCr2_min = pow(radius,2)-pow(y2_rot,2);    // cut out the inner circle
  double radius2_DCr2_max = pow(245,2)-pow(y2_rot,2);       // cut out the outer circle

  //return true;
  if (x2_rot > left && x2_rot > right && pow(x2_rot,2) > radius2_DCr2_min && pow(x2_rot,2) < radius2_DCr2_max) return true;
  else return false;

}

bool DC_hit_position_region3_fiducial_cut_hadrons_positive(int j){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double add = 0;  // value in cm added to the height and radius of the cut
  if(tight == true){  add = 3.0; }
  if(medium == true){ add = 0.0; }
  if(loose == true){  add = -3.0; }

  double angle = 60; 
  double height = 0; 
  double radius = 0;

  double height_inb[]  = {52, 52, 52, 52, 52, 52};
  double height_outb[] = {43, 43, 43, 43, 43, 43};

  int sec = part_DC_sector[j]-1;   

  for(Int_t k = 0; k < 6; k++){  
    if(sec == k && inbending == true){
      height = height_inb[k] + add;
      radius = 71 + add;
    }
    if(sec == k && outbending == true){
      height = height_outb[k] + add;
      radius = 62 + add;
    }
  }
        
  double x3_rot = part_DC_c3y[j] * sin(sec*60.0*Pival/180) + part_DC_c3x[j] * cos(sec*60.0*Pival/180);
  double y3_rot = part_DC_c3y[j] * cos(sec*60.0*Pival/180) - part_DC_c3x[j] * sin(sec*60.0*Pival/180);

  double slope = 1/tan(0.5*angle*Pival/180);
  double left  = (height - slope * y3_rot);
  double right = (height + slope * y3_rot);

  double radius2_DCr3_min = pow(radius,2)-pow(y3_rot,2);    // cut out the inner circle
  double radius2_DCr3_max = pow(355,2)-pow(y3_rot,2);       // cut out the outer circle

  //return true;
  if (x3_rot > left && x3_rot > right && pow(x3_rot,2) > radius2_DCr3_min && pow(x3_rot,2) < radius2_DCr3_max) return true;
  else return false;

}


bool DC_hit_position_region1_fiducial_cut_hadrons_negative(int j){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double add = 0;  // value in cm added to the height and radius of the cut
  if(tight == true){  add = 1.0; }
  if(medium == true){ add = 0.0; }
  if(loose == true){  add = -1.0; }

  double angle = 60; 
  double height = 0; 
  double radius = 0;

  double height_inb[]  = {21, 21, 21, 21, 21, 21};
  double height_outb[] = {15, 15, 15, 15, 15, 15};

  int sec = part_DC_sector[j]-1;   

  for(Int_t k = 0; k < 6; k++){  
    if(sec == k && inbending == true){
      height = height_inb[k] + add;
      radius = 30 + add;
    }
    if(sec == k && outbending == true){
      height = height_outb[k] + add;
      radius = 25 + add;
    }
  }
    
  double x1_rot = part_DC_c1y[j] * sin(sec*60.0*Pival/180) + part_DC_c1x[j] * cos(sec*60.0*Pival/180);
  double y1_rot = part_DC_c1y[j] * cos(sec*60.0*Pival/180) - part_DC_c1x[j] * sin(sec*60.0*Pival/180);

  double slope = 1/tan(0.5*angle*Pival/180);
  double left  = (height - slope * y1_rot);
  double right = (height + slope * y1_rot);

  double radius2_DCr1_min = pow(radius,2)-pow(y1_rot,2);    // cut out the inner circle
  double radius2_DCr1_max = pow(155,2)-pow(y1_rot,2);       // cut out the outer circle

  if (x1_rot > left && x1_rot > right && pow(x1_rot,2) > radius2_DCr1_min && pow(x1_rot,2) < radius2_DCr1_max) return true;
  else return false;

}

bool DC_hit_position_region2_fiducial_cut_hadrons_negative(int j){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double add = 0;  // value in cm added to the height and radius of the cut
  if(tight == true){  add = 2.0; }
  if(medium == true){ add = 0.0; }
  if(loose == true){  add = -2.0; }

  double angle = 60; 
  double height = 0; 
  double radius = 0;

  double height_inb[]  = {31, 31, 31, 31, 31, 31};
  double height_outb[] = {29, 29, 29, 29, 29, 29};

  int sec = part_DC_sector[j]-1;   

  for(Int_t k = 0; k < 6; k++){  
    if(sec == k && inbending == true){
      height = height_inb[k] + add;
      radius = 45 + add;
    }
    if(sec == k && outbending == true){
      height = height_outb[k] + add;
      radius = 43 + add;
    }
  }
    
  double x2_rot = part_DC_c2y[j] * sin(sec*60.0*Pival/180) + part_DC_c2x[j] * cos(sec*60.0*Pival/180);
  double y2_rot = part_DC_c2y[j] * cos(sec*60.0*Pival/180) - part_DC_c2x[j] * sin(sec*60.0*Pival/180);

  double slope = 1/tan(0.5*angle*Pival/180);
  double left  = (height - slope * y2_rot);
  double right = (height + slope * y2_rot);

  double radius2_DCr2_min = pow(radius,2)-pow(y2_rot,2);    // cut out the inner circle
  double radius2_DCr2_max = pow(245,2)-pow(y2_rot,2);       // cut out the outer circle

  //return true;
  if (x2_rot > left && x2_rot > right && pow(x2_rot,2) > radius2_DCr2_min && pow(x2_rot,2) < radius2_DCr2_max) return true;
  else return false;

}

bool DC_hit_position_region3_fiducial_cut_hadrons_negative(int j){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double add = 0;  // value in cm added to the height and radius of the cut
  if(tight == true){  add = 3.0; }
  if(medium == true){ add = 0.0; }
  if(loose == true){  add = -3.0; }

  double angle = 60; 
  double height = 0; 
  double radius = 0;

  double height_inb[]  = {42, 42, 42, 42, 42, 42};
  double height_outb[] = {53, 53, 53, 53, 53, 53};

  int sec = part_DC_sector[j]-1;   

  for(Int_t k = 0; k < 6; k++){  
    if(sec == k && inbending == true){
      height = height_inb[k] + add;
      radius = 51 + add;
    }
    if(sec == k && outbending == true){
      height = height_outb[k] + add;
      radius = 68 + add;
    }
  }
    
  double x3_rot = part_DC_c3y[j] * sin(sec*60.0*Pival/180) + part_DC_c3x[j] * cos(sec*60.0*Pival/180);
  double y3_rot = part_DC_c3y[j] * cos(sec*60.0*Pival/180) - part_DC_c3x[j] * sin(sec*60.0*Pival/180);

  double slope = 1/tan(0.5*angle*Pival/180);
  double left  = (height - slope * y3_rot);
  double right = (height + slope * y3_rot);

  double radius2_DCr3_min = pow(radius,2)-pow(y3_rot,2);    // cut out the inner circle
  double radius2_DCr3_max = pow(355,2)-pow(y3_rot,2);       // cut out the outer circle

  //return true;
  if (x3_rot > left && x3_rot > right && pow(x3_rot,2) > radius2_DCr3_min && pow(x3_rot,2) < radius2_DCr3_max) return true;
  else return false;

}


// c) beta cuts

bool prot_beta_cut(int j, int run){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double sigma_range = 3;
  if(tight == true){  sigma_range = 2.0; }
  if(medium == true){ sigma_range = 3.0; }
  if(loose == true){  sigma_range = 4.0; }

  double prot_mean_p0[] = {1.00236, 0.999663, 0.999094, 0.999417, 1.00021, 1.00034};
  double prot_mean_p1[] = {0.933781, 0.919325, 0.908676, 0.919357, 0.90765, 0.908293};
  double prot_sigma_p0[] = {0.0100589, 0.010606, 0.0101618, 0.0113726, 0.0098664, 0.00972719};
  double prot_sigma_p1[] = {-0.0119585, -0.0125292, -0.0115819, -0.0140509, -0.0117452, -0.0109518};
  double prot_sigma_p2[] = {0.0208288, 0.020747, 0.0202862, 0.0216464, 0.0207485, 0.0202303};

  double mean = 0;
  double sigma = 0;
  double upper_lim = 0;
  double lower_lim = 0;

  for(Int_t k = 0; k < 6; k++){  
    if(part_FTOF_sector_layer2[j]-1 == k){
      mean = prot_mean_p0[k] * part_p[j] / sqrt(pow(part_p[j],2) + prot_mean_p1[k]);
      sigma = prot_sigma_p0[k] + prot_sigma_p1[k]/sqrt(part_p[j]) + prot_sigma_p2[k]/(part_p[j]*part_p[j]);
      upper_lim = mean + sigma_range * sigma;
      lower_lim = mean - sigma_range * sigma;
    }
  }

  if(Beta_charged(j, run) <= upper_lim && Beta_charged(j, run) >= lower_lim) return true;
  else return false;
}

bool neutr_beta_cut(int j, int run){
  if(Beta_neutral(j, run) < 0.95 && Beta_neutral(j, run) > 0) return true;
  else return false;
}

bool pip_beta_cut(int j, int run){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double sigma_range = 3;
  if(tight == true){  sigma_range = 2.0; }
  if(medium == true){ sigma_range = 3.0; }
  if(loose == true){  sigma_range = 4.0; }

  double pip_mean_p0[] = {0.99839, 0.996913, 0.998191, 0.996847, 0.999992, 0.999761};
  double pip_mean_p1[] = {0.0222201, 0.0211859, 0.0194522, 0.0186658, 0.0206589, 0.022149};
  double pip_sigma_p0[] = {0.00429074, 0.00515927, 0.0050134, 0.00562043, 0.00558013, 0.00524988};
  double pip_sigma_p1[] = {0.000723371, -0.000719322, -0.000348278, -0.00176441, -0.00157653, -0.00079676};
  double pip_sigma_p2[] = {0.00530764, 0.0059365, 0.00564377, 0.00663122, 0.00615767, 0.00606232};

  double mean = 0;
  double sigma = 0;
  double upper_lim = 0;
  double lower_lim = 0;

  for(Int_t k = 0; k < 6; k++){  
    if(part_FTOF_sector_layer2[j]-1 == k){
      mean = pip_mean_p0[k] * part_p[j] / sqrt(pow(part_p[j],2) + pip_mean_p1[k]);
      sigma = pip_sigma_p0[k] + pip_sigma_p1[k]/sqrt(part_p[j]) + pip_sigma_p2[k]/(part_p[j]*part_p[j]);
      upper_lim = mean + sigma_range * sigma;
      lower_lim = mean - sigma_range * sigma;
    }
  }

  if(Beta_charged(j, run) <= upper_lim && Beta_charged(j, run) >= lower_lim) return true;
  else return false;
}

bool pim_beta_cut(int j, int run){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double sigma_range = 3;
  if(tight == true){  sigma_range = 2.0; }
  if(medium == true){ sigma_range = 3.0; }
  if(loose == true){  sigma_range = 4.0; }

  double pim_mean_p0[] = {0.998209, 0.99674, 0.997967, 0.99668, 0.999806, 0.999548};
  double pim_mean_p1[] = {0.0205626, 0.0191593, 0.0174171, 0.0170118, 0.0190239, 0.0201642};
  double pim_sigma_p0[] = {0.0048138, 0.00520169, 0.00570919, 0.0062253, 0.00617771, 0.00578081};
  double pim_sigma_p1[] = {-0.000382083, -0.000771668, -0.00179812, -0.00299382, -0.00282061, -0.00191396};
  double pim_sigma_p2[] = {0.00617906, 0.00592537, 0.00671878, 0.00752452, 0.00707839, 0.00691364};

  double mean = 0;
  double sigma = 0;
  double upper_lim = 0;
  double lower_lim = 0;

  for(Int_t k = 0; k < 6; k++){  
    if(part_FTOF_sector_layer2[j]-1 == k){
      mean = pim_mean_p0[k] * part_p[j] / sqrt(pow(part_p[j],2) + pim_mean_p1[k]);
      sigma = pim_sigma_p0[k] + pim_sigma_p1[k]/sqrt(part_p[j]) + pim_sigma_p2[k]/(part_p[j]*part_p[j]);
      upper_lim = mean + sigma_range * sigma;
      lower_lim = mean - sigma_range * sigma;
    }
  }

  if(Beta_charged(j, run) <= upper_lim && Beta_charged(j, run) >= lower_lim) return true;
  else return false;
}

bool Kp_beta_cut(int j, int run){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double sigma_range = 3;
  if(tight == true){  sigma_range = 2.0; }
  if(medium == true){ sigma_range = 3.0; }
  if(loose == true){  sigma_range = 4.0; }

  double Kp_mean_p0[] = {1.00, 1.00, 1.00, 1.00, 1.00, 1.00};         // literature
  double Kp_mean_p1[] = {0.244, 0.244, 0.244, 0.244, 0.244, 0.244};   // literature
  double Kp_sigma_p0[] = {0.00429074, 0.00515927, 0.0050134, 0.00562043, 0.00558013, 0.00524988};           // copied from pip
  double Kp_sigma_p1[] = {0.000723371, -0.000719322, -0.000348278, -0.00176441, -0.00157653, -0.00079676};  // copied from pip
  double Kp_sigma_p2[] = {0.00530764, 0.0059365, 0.00564377, 0.00663122, 0.00615767, 0.00606232};           // copied from pip

  double mean = 0;
  double sigma = 0;
  double upper_lim = 0;
  double lower_lim = 0;

  for(Int_t k = 0; k < 6; k++){  
    if(part_FTOF_sector_layer2[j]-1 == k){
      mean = Kp_mean_p0[k] * part_p[j] / sqrt(pow(part_p[j],2) + Kp_mean_p1[k]);
      sigma = Kp_sigma_p0[k] + Kp_sigma_p1[k]/sqrt(part_p[j]) + Kp_sigma_p2[k]/(part_p[j]*part_p[j]);
      upper_lim = mean + sigma_range * sigma;
      lower_lim = mean - sigma_range * sigma;
    }
  }

  if(Beta_charged(j, run) <= upper_lim && Beta_charged(j, run) >= lower_lim) return true;
  else return false;
}

bool Km_beta_cut(int j, int run){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double sigma_range = 3;
  if(tight == true){  sigma_range = 2.0; }
  if(medium == true){ sigma_range = 3.0; }
  if(loose == true){  sigma_range = 4.0; }

  double Km_mean_p0[] = {1.00, 1.00, 1.00, 1.00, 1.00, 1.00};        // literature
  double Km_mean_p1[] = {0.244, 0.244, 0.244, 0.244, 0.244, 0.244};  // literature
  double Km_sigma_p0[] = {0.0048138, 0.00520169, 0.00570919, 0.0062253, 0.00617771, 0.00578081};            // copied from pim
  double Km_sigma_p1[] = {-0.000382083, -0.000771668, -0.00179812, -0.00299382, -0.00282061, -0.00191396};  // copied from pim
  double Km_sigma_p2[] = {0.00617906, 0.00592537, 0.00671878, 0.00752452, 0.00707839, 0.00691364};          // copied from pim

  double mean = 0;
  double sigma = 0;
  double upper_lim = 0;
  double lower_lim = 0;

  for(Int_t k = 0; k < 6; k++){  
    if(part_FTOF_sector_layer2[j]-1 == k){
      mean = Km_mean_p0[k] * part_p[j] / sqrt(pow(part_p[j],2) + Km_mean_p1[k]);
      sigma = Km_sigma_p0[k] + Km_sigma_p1[k]/sqrt(part_p[j]) + Km_sigma_p1[k]/(part_p[j]*part_p[j]);
      upper_lim = mean + sigma_range * sigma;
      lower_lim = mean - sigma_range * sigma;
    }
  }

  if(Beta_charged(j, run) <= upper_lim && Beta_charged(j, run) >= lower_lim) return true;
  else return false;
}


// d) delta beta cuts

bool prot_delta_beta_cut(int j, int run){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double nsigma = 3;
  if(tight == true){  nsigma = 2.0; }
  if(medium == true){ nsigma = 3.0; }
  if(loose == true){  nsigma = 4.0; }

  double mean = -0.0009413;
  double sigma = 0.009336;
  double delta_beta_min = mean - nsigma * sigma;
  double delta_beta_max = mean + nsigma * sigma;

  if((part_p[j]/sqrt(part_p[j]*part_p[j]+m_p*m_p) - Beta_charged(j, run)) > delta_beta_min && (part_p[j]/sqrt(part_p[j]*part_p[j]+m_p*m_p) - Beta_charged(j, run)) < delta_beta_max) return true;
  else return false;
}

bool neutr_delta_beta_cut(int j, int run){   // maybe usefull for CND

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double nsigma = 3;
  if(tight == true){  nsigma = 2.0; }
  if(medium == true){ nsigma = 3.0; }
  if(loose == true){  nsigma = 4.0; }

  double mean = 0.0;
  double sigma = 0.001;
  double delta_beta_min = mean - nsigma * sigma;
  double delta_beta_max = mean + nsigma * sigma;

  if((part_p[j]/sqrt(part_p[j]*part_p[j]+m_n*m_n) - Beta_charged(j, run)) > delta_beta_min && (part_p[j]/sqrt(part_p[j]*part_p[j]+m_n*m_n) - Beta_charged(j, run)) < delta_beta_max) return true;
  else return false;
}

bool pip_delta_beta_cut(int j, int run){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double nsigma = 3;
  if(tight == true){  nsigma = 2.0; }
  if(medium == true){ nsigma = 3.0; }
  if(loose == true){  nsigma = 4.0; }

  double mean = -0.00397;
  double sigma = 0.008001;
  double delta_beta_min = mean - nsigma * sigma;
  double delta_beta_max = mean + nsigma * sigma;

  if((part_p[j]/sqrt(part_p[j]*part_p[j]+m_pip*m_pip) - Beta_charged(j, run)) > delta_beta_min && (part_p[j]/sqrt(part_p[j]*part_p[j]+m_pip*m_pip) - Beta_charged(j, run)) < delta_beta_max) return true;
  else return false;
}

bool pim_delta_beta_cut(int j, int run){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double nsigma = 3;
  if(tight == true){  nsigma = 2.0; }
  if(medium == true){ nsigma = 3.0; }
  if(loose == true){  nsigma = 4.0; }

  double mean = -0.003268;
  double sigma = 0.007284;
  double delta_beta_min = mean - nsigma * sigma;
  double delta_beta_max = mean + nsigma * sigma;

  if((part_p[j]/sqrt(part_p[j]*part_p[j]+m_pim*m_pim) - Beta_charged(j, run)) > delta_beta_min && (part_p[j]/sqrt(part_p[j]*part_p[j]+m_pim*m_pim) - Beta_charged(j, run)) < delta_beta_max) return true;
  else return false;
}

bool Kp_delta_beta_cut(int j, int run){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double nsigma = 3;
  if(tight == true){  nsigma = 2.0; }
  if(medium == true){ nsigma = 3.0; }
  if(loose == true){  nsigma = 4.0; }

  double mean = -0.0008986;
  double sigma = 0.005239;
  double delta_beta_min = mean - nsigma * sigma;
  double delta_beta_max = mean + nsigma * sigma;

  if((part_p[j]/sqrt(part_p[j]*part_p[j]+m_Kp*m_Kp) - Beta_charged(j, run)) > delta_beta_min && (part_p[j]/sqrt(part_p[j]*part_p[j]+m_Kp*m_Kp) - Beta_charged(j, run)) < delta_beta_max) return true;
  else return false;
}

bool Km_delta_beta_cut(int j, int run){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double nsigma = 3;
  if(tight == true){  nsigma = 2.0; }
  if(medium == true){ nsigma = 3.0; }
  if(loose == true){  nsigma = 4.0; }

  double mean = -0.001214;
  double sigma = 0.006549;
  double delta_beta_min = mean - nsigma * sigma;
  double delta_beta_max = mean + nsigma * sigma;

  if((part_p[j]/sqrt(part_p[j]*part_p[j]+m_Km*m_Km) - Beta_charged(j, run)) > delta_beta_min && (part_p[j]/sqrt(part_p[j]*part_p[j]+m_Km*m_Km) - Beta_charged(j, run)) < delta_beta_max) return true;
  else return false;
}


// e) tofmass cuts

bool prot_tofmass_cut(int j, int run){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double nsigma = 3;
  if(tight == true){  nsigma = 2.0; }
  if(medium == true){ nsigma = 3.0; }
  if(loose == true){  nsigma = 4.0; }

  double mean = 0.8824;
  double sigma = 0.07534;
  double tofmass2_min = mean - nsigma * sigma;
  double tofmass2_max = mean + nsigma * sigma;

  if(GetTOFmass2(j, run) > tofmass2_min && GetTOFmass2(j, run) < tofmass2_max) return true;
  else return false;
}

bool neutr_tofmass_cut(int j, int run){   // maybe usefull fro CND

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double nsigma = 3;
  if(tight == true){  nsigma = 2.0; }
  if(medium == true){ nsigma = 3.0; }
  if(loose == true){  nsigma = 4.0; }

  double mean = 0.8824;
  double sigma = 0.1;
  double tofmass2_min = mean - nsigma * sigma;
  double tofmass2_max = mean + nsigma * sigma;

  if(GetTOFmass2(j, run) > tofmass2_min && GetTOFmass2(j, run) < tofmass2_max) return true;
  else return false;
}

bool pip_tofmass_cut(int j, int run){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double nsigma = 3;
  if(tight == true){  nsigma = 2.0; }
  if(medium == true){ nsigma = 3.0; }
  if(loose == true){  nsigma = 4.0; }

  double mean = 0.01182;
  double sigma = 0.03106;
  double tofmass2_min = mean - nsigma * sigma;
  double tofmass2_max = mean + nsigma * sigma;
  
  if(GetTOFmass2(j, run) < tofmass2_max) return true;
  else return false;
}

bool pim_tofmass_cut(int j, int run){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double nsigma = 3;
  if(tight == true){  nsigma = 2.0; }
  if(medium == true){ nsigma = 3.0; }
  if(loose == true){  nsigma = 4.0; }

  double mean = 0.00872;
  double sigma = 0.01984;
  double tofmass2_min = mean - nsigma * sigma;
  double tofmass2_max = mean + nsigma * sigma;

  if(GetTOFmass2(j, run) < tofmass2_max) return true;
  else return false;
}

bool Kp_tofmass_cut(int j, int run){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double nsigma = 3;
  if(tight == true){  nsigma = 2.0; }
  if(medium == true){ nsigma = 3.0; }
  if(loose == true){  nsigma = 4.0; }

  double mean = 0.2157;
  double sigma = 0.0754;
  double tofmass2_min = mean - nsigma * sigma;
  double tofmass2_max = mean + nsigma * sigma;

  if(GetTOFmass2(j, run) > tofmass2_min && GetTOFmass2(j, run) < tofmass2_max) return true;
  else return false;
}

bool Km_tofmass_cut(int j, int run){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double nsigma = 3;
  if(tight == true){  nsigma = 2.0; }
  if(medium == true){ nsigma = 3.0; }
  if(loose == true){  nsigma = 4.0; }

  double mean = 0.2241;
  double sigma = 0.05122;
  double tofmass2_min = mean - nsigma * sigma;
  double tofmass2_max = mean + nsigma * sigma;

  if(GetTOFmass2(j, run) > tofmass2_min && GetTOFmass2(j, run) < tofmass2_max) return true;
  else return false;
}


// f) maximum probability cut

bool maximum_probability_cut(int j, int hypothesis, double conflvl, double anticonflvl, int run){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double add1 = 0.0;
  if(tight == true){  add1 =  4.73; }
  if(medium == true){ add1 =  0.00; }
  if(loose == true){  add1 = -0.17; }

  conflvl = conflvl + add1;

  double add2 = 0.0;
  if(tight == true){  add2 = -4.73; }
  if(medium == true){ add2 =  0.00; }
  if(loose == true){  add2 = +0.17; }

  anticonflvl = anticonflvl + add2;

  // possible hypotheses which will be tested
  //
  // proton: 2212    pip:  211     Kp:  321
  //                 pim: -211     Km: -321
  //
  // particle variables:
  
  double sector = part_FTOF_sector_layer2[j];
  double charge = part_charge[j];
  double mom = part_p[j];
  double beta = Beta_charged(j, run);

  // //////////////////////////////////////////////////////////////////////////////////////////////////
  // mean value and resolution for beta as a function of p for the different sectors

  double prot_mean_p0[] = {1.00236, 0.999663, 0.999094, 0.999417, 1.00021, 1.00034};
  double prot_mean_p1[] = {0.933781, 0.919325, 0.908676, 0.919357, 0.90765, 0.908293};
  double prot_sigma_p0[] = {0.0100589, 0.010606, 0.0101618, 0.0113726, 0.0098664, 0.00972719};
  double prot_sigma_p1[] = {-0.0119585, -0.0125292, -0.0115819, -0.0140509, -0.0117452, -0.0109518};
  double prot_sigma_p2[] = {0.0208288, 0.020747, 0.0202862, 0.0216464, 0.0207485, 0.0202303};

  double pip_mean_p0[] = {0.99839, 0.996913, 0.998191, 0.996847, 0.999992, 0.999761};
  double pip_mean_p1[] = {0.0222201, 0.0211859, 0.0194522, 0.0186658, 0.0206589, 0.022149};
  double pip_sigma_p0[] = {0.00429074, 0.00515927, 0.0050134, 0.00562043, 0.00558013, 0.00524988};
  double pip_sigma_p1[] = {0.000723371, -0.000719322, -0.000348278, -0.00176441, -0.00157653, -0.00079676};
  double pip_sigma_p2[] = {0.00530764, 0.0059365, 0.00564377, 0.00663122, 0.00615767, 0.00606232};

  double Kp_mean_p0[] = {1.00, 1.00, 1.00, 1.00, 1.00, 1.00};         // literature
  double Kp_mean_p1[] = {0.244, 0.244, 0.244, 0.244, 0.244, 0.244};   // literature
  double Kp_sigma_p0[] = {0.00429074, 0.00515927, 0.0050134, 0.00562043, 0.00558013, 0.00524988};           // copied from pip
  double Kp_sigma_p1[] = {0.000723371, -0.000719322, -0.000348278, -0.00176441, -0.00157653, -0.00079676};  // copied from pip
  double Kp_sigma_p2[] = {0.00530764, 0.0059365, 0.00564377, 0.00663122, 0.00615767, 0.00606232};           // copied from pip

  double pim_mean_p0[] = {0.998209, 0.99674, 0.997967, 0.99668, 0.999806, 0.999548};
  double pim_mean_p1[] = {0.0205626, 0.0191593, 0.0174171, 0.0170118, 0.0190239, 0.0201642};
  double pim_sigma_p0[] = {0.0048138, 0.00520169, 0.00570919, 0.0062253, 0.00617771, 0.00578081};
  double pim_sigma_p1[] = {-0.000382083, -0.000771668, -0.00179812, -0.00299382, -0.00282061, -0.00191396};
  double pim_sigma_p2[] = {0.00617906, 0.00592537, 0.00671878, 0.00752452, 0.00707839, 0.00691364};

  double Km_mean_p0[] = {1.00, 1.00, 1.00, 1.00, 1.00, 1.00};        // literature
  double Km_mean_p1[] = {0.244, 0.244, 0.244, 0.244, 0.244, 0.244};  // literature
  double Km_sigma_p0[] = {0.0048138, 0.00520169, 0.00570919, 0.0062253, 0.00617771, 0.00578081};            // copied from pim
  double Km_sigma_p1[] = {-0.000382083, -0.000771668, -0.00179812, -0.00299382, -0.00282061, -0.00191396};  // copied from pim
  double Km_sigma_p2[] = {0.00617906, 0.00592537, 0.00671878, 0.00752452, 0.00707839, 0.00691364};          // copied from pim


  // //////////////////////////////////////////////////////////////////////////////////////////////////
  // population factors for the different particles (integrated over p)

  // intially no population weighting:

  double popfrac_proton = 1.0; 
  double popfrac_pip = 1.0;
  double popfrac_Kp = 1.0;
  double popfrac_pim = 1.0;
  double popfrac_Km = 1.0;

  // momentum dependent population factor:

  if(population_weighting == true){
    if(inbending == true){   // inbending (torus -1)  
      popfrac_proton = -0.46270 + 1.497000000*pow(part_p[j],1) - 0.9120000000*pow(part_p[j],2) + 0.242400000*pow(part_p[j],3) - 0.0260900000*pow(part_p[j],4) - 0.0009616000000*pow(part_p[j],5) 
                                + 0.000520700*pow(part_p[j],6) - 0.0000484100*pow(part_p[j],7) + 0.000001559*pow(part_p[j],8) + 0.0000000071*pow(part_p[j],9) - 0.0000000008185*pow(part_p[j],10); 
      popfrac_pip =     1.08600 - 0.739500000*pow(part_p[j],1) + 0.3374000000*pow(part_p[j],2) - 0.065160000*pow(part_p[j],3) + 0.0060100000*pow(part_p[j],4) - 0.0002440000000*pow(part_p[j],5) 
                                + 0.000002776*pow(part_p[j],6); 
      popfrac_Kp =      0.08236 - 0.091010000*pow(part_p[j],1) + 0.0784700000*pow(part_p[j],2) - 0.021870000*pow(part_p[j],3) + 0.0023840000*pow(part_p[j],4) - 0.0000307300000*pow(part_p[j],5) 
                                - 0.000011330*pow(part_p[j],6) + 0.0000005536*pow(part_p[j],7); 

      popfrac_pim =     0.95340 + 0.021960000*pow(part_p[j],1) - 0.0405300000*pow(part_p[j],2) + 0.009071000*pow(part_p[j],3) - 0.0006846000*pow(part_p[j],4) + 0.0000155900000*pow(part_p[j],5); 
      popfrac_Km =     -0.18360 + 0.714600000*pow(part_p[j],1) - 0.7519000000*pow(part_p[j],2) + 0.390900000*pow(part_p[j],3) - 0.1085000000*pow(part_p[j],4) + 0.0170200000000*pow(part_p[j],5) 
                                - 0.001513000*pow(part_p[j],6) + 0.0000708600*pow(part_p[j],7) - 0.000001348*pow(part_p[j],8); 
    }
    if(outbending == true){   // outbending (torus +1)
      popfrac_proton = -0.43870 + 1.5410000000*pow(part_p[j],1)  - 0.21180000000*pow(part_p[j],2)  - 1.0700000*pow(part_p[j],3) + 0.95490000*pow(part_p[j],4) - 0.38710000*pow(part_p[j],5) 
                                + 0.0858300000*pow(part_p[j],6)  - 0.00966400000*pow(part_p[j],7)  + 0.0001328*pow(part_p[j],8) +  0.0001076*pow(part_p[j],9) - 0.00001395*pow(part_p[j],10) 
                                + 0.0000007545*pow(part_p[j],11) - 0.00000001582*pow(part_p[j],12);  
      popfrac_pip =     0.90020 - 0.7225000000*pow(part_p[j],1)  + 0.41450000000*pow(part_p[j],2)  - 0.1149000*pow(part_p[j],3) + 0.01728000*pow(part_p[j],4) - 0.00132500*pow(part_p[j],5) 
                                + 0.0000400800*pow(part_p[j],6); 
      popfrac_Kp =      0.17180 - 0.1994000000*pow(part_p[j],1)  + 0.11450000000*pow(part_p[j],2)  - 0.0190200*pow(part_p[j],3) + 0.00040970*pow(part_p[j],4) + 0.00012380*pow(part_p[j],5) 
                                - 0.0000070660*pow(part_p[j],6); 

      popfrac_pim =     0.95950 + 0.0063400000*pow(part_p[j],1) - 0.030890000000*pow(part_p[j],2)  + 0.0052940*pow(part_p[j],3) - 0.00024860*pow(part_p[j],4);
      popfrac_Km =      0.02474 + 0.0132600000*pow(part_p[j],1) + 0.024240000000*pow(part_p[j],2)  - 0.0045600*pow(part_p[j],3) + 0.00022150*pow(part_p[j],4);
    }
  }


  if(charge > 0){
    if(hypothesis <= 0) return false;    // charge does not match with hypothesis

    for(Int_t k = 0; k < 6; k++){  
      if(sector-1 == k){

        double mean_prot = prot_mean_p0[k] * part_p[j] / sqrt(pow(part_p[j],2) + prot_mean_p1[k]);
        double sigma_prot = prot_sigma_p0[k] + prot_sigma_p1[k]/sqrt(part_p[j]) + prot_sigma_p2[k]/(part_p[j]*part_p[j]);
        double prob_prot = popfrac_proton * (1/(sigma_prot*sqrt(2*3.14159))) * exp(-0.5 * pow((beta - mean_prot)/sigma_prot, 2));
        double conf_prot = 100*(1.0 - TMath::Erf(fabs(beta - mean_prot)/sigma_prot/sqrt(2.0))); 

        double mean_pip = pip_mean_p0[k] * part_p[j] / sqrt(pow(part_p[j],2) + pip_mean_p1[k]);
        double sigma_pip = pip_sigma_p0[k] + pip_sigma_p1[k]/sqrt(part_p[j]) + pip_sigma_p2[k]/(part_p[j]*part_p[j]);
        double prob_pip = popfrac_pip * (1/(sigma_pip*sqrt(2*3.14159))) * exp(-0.5 * pow((beta - mean_pip)/sigma_pip, 2));
        double conf_pip = 100*(1.0 - TMath::Erf(fabs(beta - mean_pip)/sigma_pip/sqrt(2.0))); 

        //double mean_Kp = Kp_mean_p0[k] * part_p[j] / sqrt(pow(part_p[j],2) + Kp_mean_p1[k]);
        double mean_Kp = part_p[j] / sqrt(pow(part_p[j],2) + pow(0.493677,2));
        double sigma_Kp = Kp_sigma_p0[k] + Kp_sigma_p1[k]/sqrt(part_p[j]) + Kp_sigma_p2[k]/(part_p[j]*part_p[j]);
        double prob_Kp = popfrac_Kp * (1/(sigma_Kp*sqrt(2*3.14159))) * exp(-0.5 * pow((beta - mean_Kp)/sigma_Kp, 2));
        double conf_Kp = 100*(1.0 - TMath::Erf(fabs(beta - mean_Kp)/sigma_Kp/sqrt(2.0))); 
        //prob_Kp = 0;   // overwrite Kaons
        //conf_Kp = 0;   // overwrite Kaons

        if(prob_prot > prob_pip  && prob_prot > prob_Kp  && hypothesis == 2212 && conf_prot > conflvl  && conf_pip  < anticonflvl  && conf_Kp  < anticonflvl){ 
        //cout << "proton  -  probability: " << prob_prot << "     confidence: " << conf_prot << "     ( mom: " << mom << " beta:" << beta << " mean: " << mean_prot << " sigma: " << sigma_prot <<  " )" << endl;
        //cout << "pip     -  probability: " << prob_pip  << "     confidence: " << conf_pip  << "     ( mom: " << mom << " beta:" << beta << " mean: " << mean_pip  << " sigma: " << sigma_pip  <<  " )" << endl;
        //cout << "Kp      -  probability: " << prob_Kp   << "     confidence: " << conf_Kp   << "     ( mom: " << mom << " beta:" << beta << " mean: " << mean_Kp   << " sigma: " << sigma_Kp   <<  " )" << endl;
        //cout << endl;
        return true;
        }
        if(prob_pip  > prob_prot && prob_pip  > prob_Kp  && hypothesis == 211  && conf_pip  > conflvl  && conf_prot < anticonflvl  && conf_Kp  < anticonflvl){ 
        //cout << "proton  -  probability: " << prob_prot << "     confidence: " << conf_prot << "     ( mom: " << mom << " beta:" << beta << " mean: " << mean_prot << " sigma: " << sigma_prot <<  " )" << endl;
        //cout << "pip     -  probability: " << prob_pip  << "     confidence: " << conf_pip  << "     ( mom: " << mom << " beta:" << beta << " mean: " << mean_pip  << " sigma: " << sigma_pip  <<  " )" << endl;
        //cout << "Kp      -  probability: " << prob_Kp   << "     confidence: " << conf_Kp   << "     ( mom: " << mom << " beta:" << beta << " mean: " << mean_Kp   << " sigma: " << sigma_Kp   <<  " )" << endl;
        //cout << endl;
        return true;
        }
        if(prob_Kp   > prob_prot && prob_Kp   > prob_pip && hypothesis == 321  && conf_Kp   > conflvl  && conf_prot < anticonflvl  && conf_pip < anticonflvl){ 
        //cout << "proton  -  probability: " << prob_prot << "     confidence: " << conf_prot << "     ( mom: " << mom << " beta:" << beta << " mean: " << mean_prot << " sigma: " << sigma_prot <<  " )" << endl;
        //cout << "pip     -  probability: " << prob_pip  << "     confidence: " << conf_pip  << "     ( mom: " << mom << " beta:" << beta << " mean: " << mean_pip  << " sigma: " << sigma_pip  <<  " )" << endl;
        //cout << "Kp      -  probability: " << prob_Kp   << "     confidence: " << conf_Kp   << "     ( mom: " << mom << " beta:" << beta << " mean: " << mean_Kp   << " sigma: " << sigma_Kp   <<  " )" << endl;
        //cout << endl;
        return true;
        }
      }
    }
  }

  if(charge < 0){
    if(hypothesis >= 0) return false;    // charge does not match with hypothesis

    for(Int_t k = 0; k < 6; k++){
      if(sector-1 == k){

        double mean_pim = pim_mean_p0[k] * part_p[j] / sqrt(pow(part_p[j],2) + pim_mean_p1[k]);
        double sigma_pim = pim_sigma_p0[k] + pim_sigma_p1[k]/sqrt(part_p[j]) + pim_sigma_p2[k]/(part_p[j]*part_p[j]);
        double prob_pim = popfrac_pim * (1/(sigma_pim*sqrt(2*3.14159))) * exp(-0.5 * pow((beta - mean_pim)/sigma_pim, 2));
        double conf_pim = 100*(1.0 - TMath::Erf(fabs(beta - mean_pim)/sigma_pim/sqrt(2.0))); 

        //double mean_Km = Km_mean_p0[k] * part_p[j] / sqrt(pow(part_p[j],2) + Km_mean_p1[k]);
        double mean_Km = part_p[j] / sqrt(pow(part_p[j],2) + pow(0.493677,2));
        double sigma_Km = Km_sigma_p0[k] + Km_sigma_p1[k]/sqrt(part_p[j]) + Km_sigma_p2[k]/(part_p[j]*part_p[j]);
        double prob_Km = popfrac_Km * (1/(sigma_Km*sqrt(2*3.14159))) * exp(-0.5 * pow((beta - mean_Km)/sigma_Km, 2));
        double conf_Km = 100*(1.0 - TMath::Erf(fabs(beta - mean_Km)/sigma_Km/sqrt(2.0))); 
        //prob_Km = 0;  // overwrite Kaons
        //conf_Km = 0;  // overwrite Kaons

        if(prob_pim  > prob_Km  && hypothesis == -211  && conf_pim  > conflvl && conf_Km  < anticonflvl){
        //cout << "pim     -  probability: " << prob_pim  << "     confidence: " << conf_pim << "     ( mom: " << mom << " beta:" << beta << " mean: " << mean_pim << " sigma: " << sigma_pim  <<  " )" << endl;
        //cout << "Km      -  probability: " << prob_Km   << "     confidence: " << conf_Km  << "     ( mom: " << mom << " beta:" << beta << " mean: " << mean_Km  << " sigma: " << sigma_Km   <<  " )" << endl;
        //cout << endl;
        return true;
        }
        if(prob_Km   > prob_pim && hypothesis == -321  && conf_Km   > conflvl && conf_pim < anticonflvl){  
        //cout << "pim     -  probability: " << prob_pim  << "     confidence: " << conf_pim << "     ( mom: " << mom << " beta:" << beta << " mean: " << mean_pim << " sigma: " << sigma_pim  <<  " )" << endl;
        //cout << "Km      -  probability: " << prob_Km   << "     confidence: " << conf_Km  << "     ( mom: " << mom << " beta:" << beta << " mean: " << mean_Km  << " sigma: " << sigma_Km   <<  " )" << endl;
        //cout << endl;
        return true;
        }
      }
    }
  }

  return false;  // return false if no particle is clearly identified or if the charge is 0.
}


// g) delta vz cuts

bool prot_delta_vz_cut(int j){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double nsigma = 3;
  if(tight == true){  nsigma = 2.0; }
  if(medium == true){ nsigma = 3.0; }
  if(loose == true){  nsigma = 4.0; }

  double mean = 0.6142;
  double sigma = 5.306;
  double dvz_min = mean - nsigma * sigma;
  double dvz_max = mean + nsigma * sigma;

  if(Getdvz(j) > dvz_min && Getdvz(j) < dvz_max) return true;
  else return false;
}

bool neutr_delta_vz_cut(int j){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double nsigma = 3;
  if(tight == true){  nsigma = 2.0; }
  if(medium == true){ nsigma = 3.0; }
  if(loose == true){  nsigma = 4.0; }

  double mean = 0.9096;
  double sigma = 5.156;
  double dvz_min = mean - nsigma * sigma;
  double dvz_max = mean + nsigma * sigma;

  if(Getdvz(j) > dvz_min && Getdvz(j) < dvz_max) return true;
  else return false;
}

bool pip_delta_vz_cut(int j){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double nsigma = 3;
  if(tight == true){  nsigma = 2.0; }
  if(medium == true){ nsigma = 3.0; }
  if(loose == true){  nsigma = 4.0; }

  double mean = 0.4677;
  double sigma = 5.453;
  double dvz_min = mean - nsigma * sigma;
  double dvz_max = mean + nsigma * sigma;

  if(Getdvz(j) > dvz_min && Getdvz(j) < dvz_max) return true;
  else return false;
}

bool pim_delta_vz_cut(int j){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double nsigma = 3;
  if(tight == true){  nsigma = 2.0; }
  if(medium == true){ nsigma = 3.0; }
  if(loose == true){  nsigma = 4.0; }

  double mean = 0.06017;
  double sigma = 5.564;
  double dvz_min = mean - nsigma * sigma;
  double dvz_max = mean + nsigma * sigma;

  if(Getdvz(j) > dvz_min && Getdvz(j) < dvz_max) return true;
  else return false;
}

bool Kp_delta_vz_cut(int j){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double nsigma = 3;
  if(tight == true){  nsigma = 2.0; }
  if(medium == true){ nsigma = 3.0; }
  if(loose == true){  nsigma = 4.0; }

  double mean = -0.7903;
  double sigma = 5.017;
  double dvz_min = mean - nsigma * sigma;
  double dvz_max = mean + nsigma * sigma;

  if(Getdvz(j) > dvz_min && Getdvz(j) < dvz_max) return true;
  else return false;
}

bool Km_delta_vz_cut(int j){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double nsigma = 3;
  if(tight == true){  nsigma = 2.0; }
  if(medium == true){ nsigma = 3.0; }
  if(loose == true){  nsigma = 4.0; }

  double mean = -2.103;
  double sigma = 5.723;
  double dvz_min = mean - nsigma * sigma;
  double dvz_max = mean + nsigma * sigma;

  if(Getdvz(j) > dvz_min && Getdvz(j) < dvz_max) return true;
  else return false;
}



// /////////////////////////////////////////////////////////////////
// FD photons: 

bool phot_default_PID_cut(int j){
  if(vpart_pid->at(j) == 22) return true;
  else return false;
}

bool phot_charge_cut(int j){
  if(part_charge[j] == 0) return true;
  else return false;
}

bool phot_beta_cut(int j, int run){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

  double p0_neutr_mean = 0.0088669;
  double p0_neutr_sigma = 0.0170422;
  double p1_neutr_sigma = 0.00764302;
  double p2_neutr_sigma = 0.00111672;

  double min, max;

  if(loose == true){   // no additional cut 
    min = 0.9;
    max = 2.0;
  }

  if(medium == true){  // -3 sigma and +4 sigma
    min = (1.0 + p0_neutr_mean/sqrt(part_p[j])) - 3 * (p0_neutr_sigma + p1_neutr_sigma/sqrt(part_p[j]) + p2_neutr_sigma/(part_p[j]*part_p[j]));
    max = (1.0 + p0_neutr_mean/sqrt(part_p[j])) + 4 * (p0_neutr_sigma + p1_neutr_sigma/sqrt(part_p[j]) + p2_neutr_sigma/(part_p[j]*part_p[j]));
  }

  if(tight == true){   // -2 sigma (neutron rejection) and +3 sigma
    min = (1.0 + p0_neutr_mean/sqrt(part_p[j])) - 2 * (p0_neutr_sigma + p1_neutr_sigma/sqrt(part_p[j]) + p2_neutr_sigma/(part_p[j]*part_p[j]));
    max = (1.0 + p0_neutr_mean/sqrt(part_p[j])) + 3 * (p0_neutr_sigma + p1_neutr_sigma/sqrt(part_p[j]) + p2_neutr_sigma/(part_p[j]*part_p[j]));
  }
  
  if(Beta_neutral(j, run) > min && Beta_neutral(j, run) < max && part_p[j] > 0.15) return true;

  else return false;
}

bool phot_EC_sampling_fraction_cut(int j){

  double ECfrac_min = 0;
  double ECfrac_max = 1;

  if(part_Cal_energy_total[j]/part_p[j] > ECfrac_min && part_Cal_energy_total[j]/part_p[j] < ECfrac_max) return true;
  else return false;
}

bool phot_EC_outer_vs_EC_inner_cut(int j){

  double edep_min = 0.01;

  if(part_Cal_ECin_energy[j] + part_Cal_ECout_energy[j] > edep_min) return true; 
  else return false;

}


bool phot_EC_hit_position_fiducial_cut(int j){

  ///////////////////////////
  bool tight = false;
  bool medium = true;
  bool loose = false;
  //////////////////////////

// Cut using the natural directions of the scintillator bars/ fibers:

  double u = part_Cal_PCAL_lu[j];
  double v = part_Cal_PCAL_lv[j];
  double w = part_Cal_PCAL_lw[j];
   
  /// v + w is going from teh side to the back end of the PCAL, u is going from side to side
  /// 1 scintillator bar is 4.5 cm wide. In the outer regions (back) double bars are used.

  double min_u = 0;
  double max_u = 420;    
  double min_v = 0;
  double max_v = 420;
  double min_w = 0;
  double max_w = 420;

  if(tight == true){
    min_u = 9.0;     // 2nd double bar cut completely + a bit more  20
    max_u = 420;     // no cut on outer side of u (done by v and w cuts)
    min_v = 9;       // 2 bars
    max_v = 411.5;   // 2nd bar cut completely
    min_w = 9;       // 2 bars
    max_w = 411.5;   // 2nd bar cut completely
  }

  if(medium == true){
    min_u = 6.75;    // center of 2nd double bar +  a bit more  6.75
    max_u = 420;     // no cut on outer side of u (done by v and w cuts)
    min_v = 6.75;    // 1.5 bars
    max_v = 413.75;  // center of 2nd bar
    min_w = 6.75;    // 1.5 bars
    max_w = 413.75;  // center of 2nd bar
  }

  if(loose == true){
    min_u = 4.5;    // outer double bar cut + a bit more  10
    max_u = 420;    // no cut on outer side of u (done by v and w cuts)
    min_v = 4.5;    // 1 bar
    max_v = 416;    // first bar cut completely
    min_w = 4.5;    // 1 bar
    max_w = 416;    // first bar cut completely
  }

  if(u > min_u && u < max_u && v > min_v && v < max_v && w > min_w && w < max_w) return true;
  else return false;

}


/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Particle ID cuts for the Forward Tagger:
/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// ////////////////////////////////////////////////////////
/// FT electrons:

bool FT_eid_charge_cut(int j){
  if(part_charge[j] == -1) return true;
  else return false;
}

bool FT_eid_PID_cut(int j){
  if(vpart_pid->at(j) == 11) return true;
  else return false;
}

bool FT_eid_FTCAL_fiducial_cut(int j){

  double p = sqrt(part_px[j]*part_px[j] + part_py[j]*part_py[j] + part_pz[j]*part_pz[j]);
  double theta = acos(vpart_pz->at(j)/p);

  if(theta*180/Pival > 2.5 && theta*180/Pival < 4.5) return true;
  else return false;
}

bool FT_eid_FTTRK_fiducial_cut(int j){

  double p = sqrt(part_px[j]*part_px[j] + part_py[j]*part_py[j] + part_pz[j]*part_pz[j]);
  double theta = acos(vpart_pz->at(j)/p);

  if(theta*180/Pival > 2.5 && theta*180/Pival < 4.5) return true;
  else return false;
}

bool FT_eid_FTHODO_fiducial_cut(int j){

  double p = sqrt(part_px[j]*part_px[j] + part_py[j]*part_py[j] + part_pz[j]*part_pz[j]);
  double theta = acos(vpart_pz->at(j)/p);

  if(theta*180/Pival > 2.5 && theta*180/Pival < 4.5) return true;
  else return false;
}

bool FT_eid_energy_vs_radius_cut(int j){

  return true;
}


/// /////////////////////////////////////////////////
/// FT photons


bool FT_photid_charge_cut(int j){
  if(part_charge[j] == 0) return true;
  else return false;
}

bool FT_photid_PID_cut(int j){
  if(vpart_pid->at(j) == 22) return true;
  else return false;
}

bool FT_photid_FTCAL_fiducial_cut(int j){

  double p = sqrt(part_px[j]*part_px[j] + part_py[j]*part_py[j] + part_pz[j]*part_pz[j]);
  double theta = acos(vpart_pz->at(j)/p);

  if(theta*180/Pival > 2.5 && theta*180/Pival < 4.5) return true;
  else return false;
}

bool FT_photid_beta_cut(int j, int run){

  return true;
}



/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Particle ID cuts for the Central Detector:
/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/// beta cuts

bool CD_prot_beta_cut(int j, int run){

  double prot_mean_p0 = 0.0;
  double prot_mean_p1 = 1.04897;
  double prot_mean_p2 = 0.998148;
  double prot_sigma_p0 = 0.0753141;
  double prot_sigma_p1 = 0.0;

  double sigma_range = 2;

  double mean = prot_mean_p0 + prot_mean_p1 * part_p[j] / sqrt(pow(part_p[j],2) + prot_mean_p2);
  double sigma = prot_sigma_p0 + prot_sigma_p1/sqrt(part_p[j]);
  double upper_lim = mean + sigma_range * sigma;
  double lower_lim = mean - sigma_range * sigma;

  if(Beta_charged_central(j, run) <= upper_lim && Beta_charged_central(j, run) >= lower_lim && part_p[j] > 0.25 && Beta_charged_central(j, run) >= 0.2) return true;
  else return false;
}


bool CD_neutr_beta_cut(int j, int run){
  if(Beta_neutral(j, run) < 0.95 && Beta_neutral(j, run) > 0) return true;
  else return false;
}


bool CD_pip_beta_cut(int j, int run){

  double pip_mean_p0 = 0.0;
  double pip_mean_p1 = 1.02181;
  double pip_mean_p2 = 0.0245488;
  double pip_sigma_p0 = 0.10908;
  double pip_sigma_p1 = 0.0;

  double sigma_range = 2;

  double mean = pip_mean_p0 + pip_mean_p1 * part_p[j] / sqrt(pow(part_p[j],2) + pip_mean_p2);
  double sigma = pip_sigma_p0 + pip_sigma_p1/sqrt(part_p[j]);
  double upper_lim = mean + sigma_range * sigma;
  double lower_lim = mean - sigma_range * sigma;

  if(Beta_charged_central(j, run) <= upper_lim && Beta_charged_central(j, run) >= lower_lim && part_p[j] > 0.2) return true;
  else return false;
}

bool CD_pim_beta_cut(int j, int run){

  double pim_mean_p0 = -2.00799;
  double pim_mean_p1 = 3.02637;
  double pim_mean_p2 = 0.0080753;
  double pim_sigma_p0 = 0.0969147;
  double pim_sigma_p1 = 0.0;

  double sigma_range = 2;

  double mean = pim_mean_p0 + pim_mean_p1 * part_p[j] / sqrt(pow(part_p[j],2) + pim_mean_p2);
  double sigma = pim_sigma_p0 + pim_sigma_p1/sqrt(part_p[j]);
  double upper_lim = mean + sigma_range * sigma;
  double lower_lim = mean - sigma_range * sigma;

  if(Beta_charged_central(j, run) <= upper_lim && Beta_charged_central(j, run) >= lower_lim && part_p[j] > 0.2) return true;
  else return false;
}

bool CD_Kp_beta_cut(int j, int run){

  double Kp_mean_p0 = 0.00000;
  double Kp_mean_p1 = 1.00000;
  double Kp_mean_p2 = 0.24372;         // lit. Kaon mass
  double Kp_sigma_p0 = 0.0969147;      // copied from pip
  double Kp_sigma_p1 = 0.0;            // copied from pip

  double sigma_range = 2;

  double mean = Kp_mean_p0 + Kp_mean_p1 * part_p[j] / sqrt(pow(part_p[j],2) + Kp_mean_p2);
  double sigma = Kp_sigma_p0 + Kp_sigma_p1/sqrt(part_p[j]);
  double upper_lim = mean + sigma_range * sigma;
  double lower_lim = mean - sigma_range * sigma;

  if(Beta_charged_central(j, run) <= upper_lim && Beta_charged_central(j, run) >= lower_lim && part_p[j] > 0.2) return true;
  else return false;
}

bool CD_Km_beta_cut(int j, int run){

  double Km_mean_p0 = 0.00000;
  double Km_mean_p1 = 1.00000;
  double Km_mean_p2 = 0.24372;        // lit. Kaon mass
  double Km_sigma_p0 = 0.0969147;     // copied from pim
  double Km_sigma_p1 = 0.0;           // copied from pim

  double sigma_range = 2;

  double mean = Km_mean_p0 + Km_mean_p1 * part_p[j] / sqrt(pow(part_p[j],2) + Km_mean_p2);
  double sigma = Km_sigma_p0 + Km_sigma_p1/sqrt(part_p[j]);
  double upper_lim = mean + sigma_range * sigma;
  double lower_lim = mean - sigma_range * sigma;

  if(Beta_charged_central(j, run) <= upper_lim && Beta_charged_central(j, run) >= lower_lim && part_p[j] > 0.2) return true;
  else return false;
}



// maximum probability

bool CD_maximum_probability_cut(int j, int hypothesis, double conflvl, double anticonflvl, int run){

  // possible hypotheses which will be tested
  //
  // proton: 2212    pip:  211     Kp:  321
  //                 pim: -211     Km: -321
  //
  // particle variables:
  
  double charge = part_charge[j];
  double mom = part_p[j];
  double beta = Beta_charged_central(j, run);

  // //////////////////////////////////////////////////////////////////////////////////////////////////
  // mean value and resolution for beta as a function of p for the different sectors

  double prot_mean_p0 = 0.0;
  double prot_mean_p1 = 1.04897;
  double prot_mean_p2 = 0.998148;
  double prot_sigma_p0 = 0.0753141;
  double prot_sigma_p1 = 0.0;

  double pip_mean_p0 = 0.0;
  double pip_mean_p1 = 1.02181;
  double pip_mean_p2 = 0.0245488;
  double pip_sigma_p0 = 0.10908;
  double pip_sigma_p1 = 0.0;

  double pim_mean_p0 = -2.00799;
  double pim_mean_p1 = 3.02637;
  double pim_mean_p2 = 0.0080753;
  double pim_sigma_p0 = 0.0969147;
  double pim_sigma_p1 = 0.0;

  double Kp_mean_p0 = 0.00000;
  double Kp_mean_p1 = 1.00000;
  double Kp_mean_p2 = 0.24372;         // lit. Kaon mass
  double Kp_sigma_p0 = 0.0969147;      // copied from pip
  double Kp_sigma_p1 = 0.0;            // copied from pip

  double Km_mean_p0 = 0.00000;
  double Km_mean_p1 = 1.00000;
  double Km_mean_p2 = 0.24372;        // lit. Kaon mass
  double Km_sigma_p0 = 0.0969147;     // copied from pim
  double Km_sigma_p1 = 0.0;           // copied from pim

  // //////////////////////////////////////////////////////////////////////////////////////////////////
  // population factors for the different particles (integrated over p)

  // intially no population weighting:

  double popfrac_proton = 1.0; 
  double popfrac_pip = 1.0;
  double popfrac_Kp = 1.0;
  double popfrac_pim = 1.0;
  double popfrac_Km = 1.0;

  // momentum dependent population factor:

  if(population_weighting_CD == true){
    if(outbending == false){   // inbending (torus -1)  
      popfrac_proton =  1.27200 - 6.6060000000*pow(part_p[j],1)  + 12.8100000000*pow(part_p[j],2) - 10.61000*pow(part_p[j],3) + 4.503000000*pow(part_p[j],4) - 1.034000000*pow(part_p[j],5) 
                                + 0.1225000000*pow(part_p[j],6)  - 0.00588400000*pow(part_p[j],7); 
      if(part_p[j] > 5.0) popfrac_proton = 0.001; 
      popfrac_pip =     0.87680 + 1.1150000000*pow(part_p[j],1)  - 3.94000000000*pow(part_p[j],2) + 3.679000*pow(part_p[j],3) - 1.659000000*pow(part_p[j],4) + 0.425800000*pow(part_p[j],5)
                                - 0.0654700000*pow(part_p[j],6)  + 0.00597200000*pow(part_p[j],7) - 0.000298*pow(part_p[j],8) + 0.000006271*pow(part_p[j],9); 
      if(part_p[j] > 6.0) popfrac_pip = 0.998; 
      popfrac_Kp =      0.04908 + 0.0908700000*pow(part_p[j],1)  - 0.032710000000*pow(part_p[j],2) + 0.002760*pow(part_p[j],3); 
      if(part_p[j] > 5.5) popfrac_Kp = 0.001; 
      popfrac_pim =     1.003000 - 0.324900000*pow(part_p[j],1) + 0.184200000*pow(part_p[j],2) - 0.0335700000*pow(part_p[j],3) + 0.00164000*pow(part_p[j],4) + 0.0001136*pow(part_p[j],5)
                                 - 0.000009368*pow(part_p[j],6);
      if(part_p[j] > 3.5) popfrac_pim = 0.999;  
      popfrac_Km =      0.065570 - 0.155300000*pow(part_p[j],1) + 0.652600000*pow(part_p[j],2) - 0.5292000000*pow(part_p[j],3) + 0.15800000*pow(part_p[j],4) - 0.0162400*pow(part_p[j],5);
      if(part_p[j] > 3.5) popfrac_Km = 0.001;  
    }
    if(outbending == true){   // outbending (torus +1)
      popfrac_proton =  0.091220 - 0.774600000*pow(part_p[j],1) + 3.121000000*pow(part_p[j],2) - 3.1810000000*pow(part_p[j],3) + 1.59300000*pow(part_p[j],4) - 0.4622000*pow(part_p[j],5) 
                                 + 0.081690000*pow(part_p[j],6) - 0.008687000*pow(part_p[j],7) + 0.0005144000*pow(part_p[j],8) - 0.00001281*pow(part_p[j],9); 
      if(part_p[j] > 8.0) popfrac_proton = 0.001; 
      popfrac_pip =     1.223000 - 1.400000000*pow(part_p[j],1) + 0.708100000*pow(part_p[j],2) - 0.1232000000*pow(part_p[j],3) + 0.00514000*pow(part_p[j],4) + 0.0005800*pow(part_p[j],5) 
                                 - 0.000013630*pow(part_p[j],6) - 0.000006621*pow(part_p[j],7) + 0.0000003665*pow(part_p[j],8); 
      if(part_p[j] > 7.0) popfrac_pip = 0.998; 
      popfrac_Kp =      -0.06797 + 0.294800000*pow(part_p[j],1) - 0.137800000*pow(part_p[j],2) + 0.0231100000*pow(part_p[j],3) - 0.00133500*pow(part_p[j],4); 
      if(part_p[j] > 6.0) popfrac_Kp = 0.001; 
      popfrac_pim =     0.93350 + 0.1171000000*pow(part_p[j],1)  - 0.23500000000*pow(part_p[j],2) - 0.254100*pow(part_p[j],3) + 0.464600000*pow(part_p[j],4) - 0.255500000*pow(part_p[j],5)
                                + 0.0702800000*pow(part_p[j],6)  - 0.01024000000*pow(part_p[j],7) + 0.000606*pow(part_p[j],8) + 0.000034870*pow(part_p[j],9) - 0.000008172*pow(part_p[j],10)
                                + 0.0000005102*pow(part_p[j],11) - 0.00000001146*pow(part_p[j],12); 
      if(part_p[j] > 4.0) popfrac_pim = 0.999; 
      popfrac_Km =      0.10210 - 0.2916000000*pow(part_p[j],1)  + 0.72860000000*pow(part_p[j],2) - 0.486800*pow(part_p[j],3) + 0.124600000*pow(part_p[j],4) - 0.011080000*pow(part_p[j],5); 
      if(part_p[j] > 4.0) popfrac_Km = 0.001;  
    }
  }


  if(charge > 0){

    if(hypothesis <= 0) return false;    // charge does not match with hypothesis

    double mean_prot = prot_mean_p0 + prot_mean_p1 * part_p[j] / sqrt(pow(part_p[j],2) + prot_mean_p2);
    double sigma_prot = prot_sigma_p0 + prot_sigma_p1/sqrt(part_p[j]);
    double prob_prot = popfrac_proton * (1/(sigma_prot*sqrt(2*3.14159))) * exp(-0.5 * pow((beta - mean_prot)/sigma_prot, 2));
    double conf_prot = 100*(1.0 - TMath::Erf(fabs(beta - mean_prot)/sigma_prot/sqrt(2.0))); 

    double mean_pip = pip_mean_p0 + pip_mean_p1 * part_p[j] / sqrt(pow(part_p[j],2) + pip_mean_p2);
    double sigma_pip = pip_sigma_p0 + pip_sigma_p1/sqrt(part_p[j]);
    double prob_pip = popfrac_pip * (1/(sigma_pip*sqrt(2*3.14159))) * exp(-0.5 * pow((beta - mean_pip)/sigma_pip, 2));
    double conf_pip = 100*(1.0 - TMath::Erf(fabs(beta - mean_pip)/sigma_pip/sqrt(2.0))); 

    double mean_Kp = Kp_mean_p0 + Kp_mean_p1 * part_p[j] / sqrt(pow(part_p[j],2) + Kp_mean_p2);
    double sigma_Kp = Kp_sigma_p0 + Kp_sigma_p1/sqrt(part_p[j]);
    double prob_Kp = popfrac_Kp * (1/(sigma_Kp*sqrt(2*3.14159))) * exp(-0.5 * pow((beta - mean_Kp)/sigma_Kp, 2));
    double conf_Kp = 100*(1.0 - TMath::Erf(fabs(beta - mean_Kp)/sigma_Kp/sqrt(2.0))); 

    if(prob_prot > prob_pip  && prob_prot > prob_Kp  && hypothesis == 2212 && conf_prot > conflvl  && conf_pip  < anticonflvl  && conf_Kp  < anticonflvl && part_p[j] > 0.2 && Beta_charged_central(j, run) >= 0.2){ 
      //cout << "proton  -  probability: " << prob_prot << "     confidence: " << conf_prot << "     ( mom: " << mom << " beta:" << beta << " mean: " << mean_prot << " sigma: " << sigma_prot <<  " )" << endl;
      //cout << "pip     -  probability: " << prob_pip  << "     confidence: " << conf_pip  << "     ( mom: " << mom << " beta:" << beta << " mean: " << mean_pip  << " sigma: " << sigma_pip  <<  " )" << endl;
      //cout << "Kp      -  probability: " << prob_Kp   << "     confidence: " << conf_Kp   << "     ( mom: " << mom << " beta:" << beta << " mean: " << mean_Kp   << " sigma: " << sigma_Kp   <<  " )" << endl;
      //cout << endl;
      return true;
    }
    if(prob_pip  > prob_prot && prob_pip  > prob_Kp  && hypothesis == 211  && conf_pip  > conflvl  && conf_prot < anticonflvl  && conf_Kp  < anticonflvl && part_p[j] > 0.2){ 
      //cout << "proton  -  probability: " << prob_prot << "     confidence: " << conf_prot << "     ( mom: " << mom << " beta:" << beta << " mean: " << mean_prot << " sigma: " << sigma_prot <<  " )" << endl;
      //cout << "pip     -  probability: " << prob_pip  << "     confidence: " << conf_pip  << "     ( mom: " << mom << " beta:" << beta << " mean: " << mean_pip  << " sigma: " << sigma_pip  <<  " )" << endl;
      //cout << "Kp      -  probability: " << prob_Kp   << "     confidence: " << conf_Kp   << "     ( mom: " << mom << " beta:" << beta << " mean: " << mean_Kp   << " sigma: " << sigma_Kp   <<  " )" << endl;
      //cout << endl;
      return true;
    }
    if(prob_Kp   > prob_prot && prob_Kp   > prob_pip && hypothesis == 321  && conf_Kp   > conflvl  && conf_prot < anticonflvl  && conf_pip < anticonflvl && part_p[j] > 0.2){ 
      //cout << "proton  -  probability: " << prob_prot << "     confidence: " << conf_prot << "     ( mom: " << mom << " beta:" << beta << " mean: " << mean_prot << " sigma: " << sigma_prot <<  " )" << endl;
      //cout << "pip     -  probability: " << prob_pip  << "     confidence: " << conf_pip  << "     ( mom: " << mom << " beta:" << beta << " mean: " << mean_pip  << " sigma: " << sigma_pip  <<  " )" << endl;
      //cout << "Kp      -  probability: " << prob_Kp   << "     confidence: " << conf_Kp   << "     ( mom: " << mom << " beta:" << beta << " mean: " << mean_Kp   << " sigma: " << sigma_Kp   <<  " )" << endl;
      //cout << endl;
      return true;
    }
  }

  if(charge < 0){

    if(hypothesis >= 0) return false;    // charge does not match with hypothesis

    double mean_pim = pim_mean_p0 + pim_mean_p1 * part_p[j] / sqrt(pow(part_p[j],2) + pim_mean_p2);
    double sigma_pim = pim_sigma_p0 + pim_sigma_p1/sqrt(part_p[j]);
    double prob_pim = popfrac_pim * (1/(sigma_pim*sqrt(2*3.14159))) * exp(-0.5 * pow((beta - mean_pim)/sigma_pim, 2));
    double conf_pim = 100*(1.0 - TMath::Erf(fabs(beta - mean_pim)/sigma_pim/sqrt(2.0))); 

    double mean_Km = Km_mean_p0 + Km_mean_p1 * part_p[j] / sqrt(pow(part_p[j],2) + Km_mean_p2);
    double sigma_Km = Km_sigma_p0 + Km_sigma_p1/sqrt(part_p[j]);
    double prob_Km = popfrac_Km * (1/(sigma_Km*sqrt(2*3.14159))) * exp(-0.5 * pow((beta - mean_Km)/sigma_Km, 2));
    double conf_Km = 100*(1.0 - TMath::Erf(fabs(beta - mean_Km)/sigma_Km/sqrt(2.0))); 

    if(prob_pim  > prob_Km  && hypothesis == -211  && conf_pim  > conflvl && conf_Km  < anticonflvl && part_p[j] > 0.2){
      //cout << "pim     -  probability: " << prob_pim  << "     confidence: " << conf_pim << "     ( mom: " << mom << " beta:" << beta << " mean: " << mean_pim << " sigma: " << sigma_pim  <<  " )" << endl;
      //cout << "Km      -  probability: " << prob_Km   << "     confidence: " << conf_Km  << "     ( mom: " << mom << " beta:" << beta << " mean: " << mean_Km  << " sigma: " << sigma_Km   <<  " )" << endl;
      //cout << endl;
      return true;
    }
    if(prob_Km   > prob_pim && hypothesis == -321  && conf_Km   > conflvl && conf_pim < anticonflvl && part_p[j] > 0.2){  
      //cout << "pim     -  probability: " << prob_pim  << "     confidence: " << conf_pim << "     ( mom: " << mom << " beta:" << beta << " mean: " << mean_pim << " sigma: " << sigma_pim  <<  " )" << endl;
      //cout << "Km      -  probability: " << prob_Km   << "     confidence: " << conf_Km  << "     ( mom: " << mom << " beta:" << beta << " mean: " << mean_Km  << " sigma: " << sigma_Km   <<  " )" << endl;
      //cout << endl;
      return true;
    }
  }

  return false;  // return false if no particle is clearly identified or if the charge is 0.
}


// dvz cuts

bool CD_prot_delta_vz_cut(int j){

  double mean = 0.7486;
  double sigma = 3.237;
  double dvz_min = mean - 3 * sigma;
  double dvz_max = mean + 3 * sigma;

  if(Getdvz(j) > dvz_min && Getdvz(j) < dvz_max) return true;
  else return false;
}

bool CD_neutr_delta_vz_cut(int j){

  double mean = 2.254;
  double sigma = 2.693;
  double dvz_min = mean - 3 * sigma;
  double dvz_max = mean + 3 * sigma;

  if(Getdvz(j) > dvz_min && Getdvz(j) < dvz_max) return true;
  else return false;
}

bool CD_pip_delta_vz_cut(int j){

  double mean = -0.6183;
  double sigma = 3.684;
  double dvz_min = mean - 3 * sigma;
  double dvz_max = mean + 3 * sigma;

  if(Getdvz(j) > dvz_min && Getdvz(j) < dvz_max) return true;
  else return false;
}

bool CD_pim_delta_vz_cut(int j){

  double mean = -0.5485;
  double sigma = 3.677;
  double dvz_min = mean - 3 * sigma;
  double dvz_max = mean + 3 * sigma;

  if(Getdvz(j) > dvz_min && Getdvz(j) < dvz_max) return true;
  else return false;
}

bool CD_Kp_delta_vz_cut(int j){

  double mean = 1.658;
  double sigma = 2.52;
  double dvz_min = mean - 3 * sigma;
  double dvz_max = mean + 3 * sigma;

  if(Getdvz(j) > dvz_min && Getdvz(j) < dvz_max) return true;
  else return false;
}

bool CD_Km_delta_vz_cut(int j){

  double mean = -1.161;
  double sigma = 2.691;
  double dvz_min = mean - 3 * sigma;
  double dvz_max = mean + 3 * sigma;

  if(Getdvz(j) > dvz_min && Getdvz(j) < dvz_max) return true;
  else return false;
}




/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Momentum corrections:
/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TLorentzVector correct_lepton_negative(TLorentzVector lep_neg_raw){


  // 3050:

  double AA[] = {-1024.3, 435.908, 48.5039, -0.0773582, 13.1506, -27.7467};
  double AB[] = {284.775, -120.385, -12.7129, 0.306666, -2.40695, 11.996};
  double AC[] = {-28.806, 12.2926, 1.25947, -0.0305063, 0.154954, -1.43404};
  double AD[] = {1.25843, -0.548384, -0.0545483, 0.00133536, -0.0030153, 0.0687841};
  double AE[] = {-0.0200686, 0.00901792, 0.000872667, -2.15769e-05, -1.18338e-05, -0.00119579};
  double BA[] = {-12.0807, 7.74428, 1.80618, 0.178173, -0.301548, 0.54574};
  double BB[] = {3.35651, -2.14446, -0.483517, -0.0544012, 0.0544793, -0.214491};
  double BC[] = {-0.339701, 0.219018, 0.0478322, 0.00598467, -0.00279228, 0.0252665};
  double BD[] = {0.0148493, -0.00977274, -0.00206876, -0.000283642, 6.36583e-06, -0.00120819};
  double BE[] = {-0.000236963, 0.000160741, 3.30498e-05, 4.90956e-06, 1.75822e-06, 2.09916e-05};
  double CA[] = {-0.0355374, 0.034334, 0.0168396, -0.0080824, 0.00173718, -0.00243796};
  double CB[] = {0.00987712, -0.00950999, -0.00450376, 0.00263701, -0.000262686, 0.000925225};
  double CC[] = {-0.0010001, 0.000971529, 0.000444998, -0.000304425, 6.11003e-06, -0.000108096};
  double CD[] = {4.37424e-05, -4.336e-05, -1.92187e-05, 1.49378e-05, 6.06691e-07, 5.16721e-06};
  double CE[] = {-6.9847e-07, 7.13304e-07, 3.06578e-07, -2.64993e-07, -2.38277e-08, -8.98978e-08};

  double phi_reg_min[] = {-180, - 124, -64, -6, 54, 120};
  double phi_reg_max[] = {-162, -100, -46, 20, 78, 140};

  double theta_min = 11;
  double theta_max = 21;


  /*
  // 2395 + 2397
  double AA[] = {3.96447, 3.56539, 1.64587, 0, 1.66861, 5.39892};
  double AB[] = {0.283796, -0.253793, -0.0964826, 0, 0.0246473, -0.645162};
  double AC[] = {-0.101252, -0.00427019, 0.00436751, 0, -0.0179275, 0.0306926};
  double AD[] = {0.00652514, 0.00107924, -3.28438e-05, 0, 0.00131993, -0.000439392};
  double AE[] = {-0.000124654, -2.58448e-05, -1.1317e-06, 0, -2.76683e-05, -9.67504e-07};
  double BA[] = {0.0382254, 0.0542664, 0.034586, 0, -0.0205287, -0.0647981};
  double BB[] = {0.0029391, -0.00639047, -0.00609641, 0, 0.000224562, 0.00947448};
  double BC[] = {-0.00119596, 7.61894e-05, 0.000380832, 0, 0.000391467, -0.000437208};
  double BD[] = {7.85484e-05, 1.38945e-05, -9.29705e-06, 0, -3.16525e-05, 5.29236e-06};
  double BE[] = {-1.51203e-06, -3.8653e-07, 7.55498e-08, 0, 6.82443e-07, 3.97512e-08};
  double CA[] = {0.000119612, 0.000278996, 0.000425972, 0, 0.000143804, 0.000233408};
  double CB[] = {8.30727e-06, -3.6505e-05, -8.11078e-05, 0, -3.5762e-06, -3.31042e-05};
  double CC[] = {-3.60433e-06, 1.00778e-06, 5.7759e-06, 0, -2.36587e-06, 1.39742e-06};
  double CD[] = {2.38629e-07, 3.64067e-08, -1.79923e-07, 0, 1.98863e-07, -9.05074e-09};
  double CE[] = {-4.60753e-09, -1.3364e-09, 2.21418e-09, 0, -4.33074e-09, -3.27964e-10};


  double phi_reg_min[] = {-178, - 120, -60, 20, 60, 120};
  double phi_reg_max[] = {-140, -81, -27, 20, 93, 153}; 

  double theta_min = 6.5;
  double theta_max = 23;
  */

/*
if(run == 2383){    // for run 2383 (100% outbending):

  double AA[] = {6.83045, 4.1226, 1.84998, 1.11049, 2.22096, 9.54057};
  double AB[] = {-0.771015, -0.524748, -0.182965, -0.0209766, -0.197437, -1.93703};
  double AC[] = {0.0221712, 0.0307074, 0.0155828, 0.0019323, 0.0115737, 0.167383};
  double AD[] = {0.000689754, -0.00064143, -0.000598563, -8.61622e-05, -0.00027472, -0.00637268};
  double AE[] = {-2.75416e-05, 2.57831e-06, 8.56679e-06, 1.40151e-06, 2.19829e-06, 8.99762e-05};
  double BA[] = {0.0723159, 0.0624633, 0.0420249, -0.0145283, -0.0306753, -0.12591};
  double BB[] = {-0.00968628, -0.0107396, -0.00941917, 0.00325576, 0.00493625, 0.0286776};
  double BC[] = {0.000289564, 0.000643635, 0.000810481, -0.000274498, -0.000276207, -0.00248044};
  double BD[] = {7.8976e-06, -1.39562e-05, -3.08215e-05, 1.01971e-05, 5.79629e-06, 9.4443e-05};
  double BE[] = {-3.30327e-07, 6.62678e-08, 4.35367e-07, -1.41476e-07, -3.28153e-08, -1.33288e-06};
  double CA[] = {0.000222868, 0.000307922, 0.000498406, 0.000440922, 0.000190687, 0.000461586};
  double CB[] = {-2.99233e-05, -5.34312e-05, -0.000112885, -9.8389e-05, -2.95135e-05, -0.000105235};
  double CC[] = {9.07008e-07, 3.23835e-06, 9.79746e-06, 8.32572e-06, 1.54429e-06, 9.11323e-06};
  double CD[] = {2.32709e-08, -7.19333e-08, -3.75518e-07, -3.10079e-07, -2.73518e-08, -3.47474e-07};
  double CE[] = {-9.93027e-10, 3.77415e-10, 5.33749e-09, 4.30429e-09, 5.9649e-11, 4.90934e-09};

  double phi_reg_min[] = {-178, - 117, -57, 3, 63, 123};   // narrow fiducial fit cut!
  double phi_reg_max[] = {-145, -88, -28, 32, 92, 154};    // narrow fiducial fir cut!

  //double phi_reg_min[] = {-180, - 124, -60, -2, 58, 119};    // wide fiducial cut on phi
  //double phi_reg_max[] = {-143,  -82,  -22, 38, 96, 156};    // wide fiducial cut on phi

  double theta_min = 6;
  double theta_max = 28;
}
*/

/*
if(run == 2549){    // 100% inbending (lower statistics(

  double AA[] = {47.1332, -2.04948, 1.18206, -1.62147, -6.29036, 142.2};
  double AB[] = {-15.8686, 1.24107, -0.0301767, 0.630746, 1.81657, -29.9516};
  double AC[] = {1.80461, -0.15149, 0.0024246, -0.0540056, -0.161698, 2.36236};
  double AD[] = {-0.084515, 0.00718902, -7.96367e-05, 0.00199487, 0.0061333, -0.0822112};
  double AE[] = {0.00140862, -0.000114974, 8.83118e-07, -2.69207e-05, -8.37394e-05, 0.00106684};
  double BA[] = {2.92365, -0.221418, -0.0197177, 0.107532, 0.140379, -1.81829};
  double BB[] = {-0.778836, 0.0650796, 0.00594125, -0.0252402, -0.0344905, 0.385599};
  double BC[] = {0.0753388, -0.00675306, -0.000633764, 0.00213721, 0.0030473, -0.0303964};
  double BD[] = {-0.00314442, 0.000292167, 2.84874e-05, -7.76969e-05, -0.000114678, 0.00105723};
  double BE[] = {4.78596e-05, -4.43637e-06, -4.54441e-07, 1.02747e-06, 1.55185e-06, -1.37111e-05};
  double CA[] = {0.0181768, -0.00212726, -0.00111583, -0.00106242, -0.000658469, -0.0288182};
  double CB[] = {-0.00466499, 0.000583096, 0.000308788, 0.000246336, 0.000160361, 0.00608346};
  double CC[] = {0.00043796, -5.76118e-05, -3.09735e-05, -2.05136e-05, -1.4039e-05, -0.000476833};
  double CD[] = {-1.78259e-05, 2.41192e-06, 1.33212e-06, 7.29507e-07, 5.22958e-07, 1.64323e-05};
  double CE[] = {2.65299e-07, -3.58955e-08, -2.06474e-08, -9.37853e-09, -6.98958e-09, -2.09859e-07};

  double phi_reg_min[] = {-155, - 100, -35, 20, 84, 140};   // narrow fiducial fit cut!
  double phi_reg_max[] = {-120, -56, 0, 60, 120, 180};      // narrow fiducial fir cut!

  double theta_min = 10;
  double theta_max = 25;

}
*/

  double theta = lep_neg_raw.Theta()*180/Pival;
  double phi = lep_neg_raw.Phi()*180/Pival;
  double mom = lep_neg_raw.P();

  double pcorr = 0;
  double A = 0;
  double B = 0;
  double C = 0;

  double torus = 100;
  double torratio = torus/100;

  if(theta > theta_min && theta < theta_max){

    if( phi > phi_reg_min[0] && phi < phi_reg_max[0]){
      A = AA[0]*pow(theta,0) + AB[0]*pow(theta,1) + AC[0]*pow(theta,2) + AD[0]*pow(theta,3) + AE[0]*pow(theta,4);
      B = BA[0]*pow(theta,0) + BB[0]*pow(theta,1) + BC[0]*pow(theta,2) + BD[0]*pow(theta,3) + BE[0]*pow(theta,4);
      C = CA[0]*pow(theta,0) + CB[0]*pow(theta,1) + CC[0]*pow(theta,2) + CD[0]*pow(theta,3) + CE[0]*pow(theta,4);
      pcorr = A * pow(phi,0) + B * pow(phi,1) + C * pow(phi,2);
    }

    if( phi > phi_reg_min[1] && phi < phi_reg_max[1]){
      A = AA[1]*pow(theta,0) + AB[1]*pow(theta,1) + AC[1]*pow(theta,2) + AD[1]*pow(theta,3) + AE[1]*pow(theta,4);
      B = BA[1]*pow(theta,0) + BB[1]*pow(theta,1) + BC[1]*pow(theta,2) + BD[1]*pow(theta,3) + BE[1]*pow(theta,4);
      C = CA[1]*pow(theta,0) + CB[1]*pow(theta,1) + CC[1]*pow(theta,2) + CD[1]*pow(theta,3) + CE[1]*pow(theta,4);
      pcorr = A * pow(phi,0) + B * pow(phi,1) + C * pow(phi,2);
    }

    if( phi > phi_reg_min[2] && phi < phi_reg_max[2]){
      A = AA[2]*pow(theta,0) + AB[2]*pow(theta,1) + AC[2]*pow(theta,2) + AD[2]*pow(theta,3) + AE[2]*pow(theta,4);
      B = BA[2]*pow(theta,0) + BB[2]*pow(theta,1) + BC[2]*pow(theta,2) + BD[2]*pow(theta,3) + BE[2]*pow(theta,4);
      C = CA[2]*pow(theta,0) + CB[2]*pow(theta,1) + CC[2]*pow(theta,2) + CD[2]*pow(theta,3) + CE[2]*pow(theta,4);
      pcorr = A * pow(phi,0) + B * pow(phi,1) + C * pow(phi,2);
    }

    if( phi > phi_reg_min[3] && phi < phi_reg_max[3]){
      A = AA[3]*pow(theta,0) + AB[3]*pow(theta,1) + AC[3]*pow(theta,2) + AD[3]*pow(theta,3) + AE[3]*pow(theta,4);
      B = BA[3]*pow(theta,0) + BB[3]*pow(theta,1) + BC[3]*pow(theta,2) + BD[3]*pow(theta,3) + BE[3]*pow(theta,4);
      C = CA[3]*pow(theta,0) + CB[3]*pow(theta,1) + CC[3]*pow(theta,2) + CD[3]*pow(theta,3) + CE[3]*pow(theta,4);
      pcorr = A * pow(phi,0) + B * pow(phi,1) + C * pow(phi,2);
    }

    if( phi > phi_reg_min[4] && phi < phi_reg_max[4]){
      A = AA[4]*pow(theta,0) + AB[4]*pow(theta,1) + AC[4]*pow(theta,2) + AD[4]*pow(theta,3) + AE[4]*pow(theta,4);
      B = BA[4]*pow(theta,0) + BB[4]*pow(theta,1) + BC[4]*pow(theta,2) + BD[4]*pow(theta,3) + BE[4]*pow(theta,4);
      C = CA[4]*pow(theta,0) + CB[4]*pow(theta,1) + CC[4]*pow(theta,2) + CD[4]*pow(theta,3) + CE[4]*pow(theta,4);
      pcorr = A * pow(phi,0) + B * pow(phi,1) + C * pow(phi,2);
    }

    if( phi > phi_reg_min[5] && phi < phi_reg_max[5]){
      A = AA[5]*pow(theta,0) + AB[5]*pow(theta,1) + AC[5]*pow(theta,2) + AD[5]*pow(theta,3) + AE[5]*pow(theta,4);
      B = BA[5]*pow(theta,0) + BB[5]*pow(theta,1) + BC[5]*pow(theta,2) + BD[5]*pow(theta,3) + BE[5]*pow(theta,4);
      C = CA[5]*pow(theta,0) + CB[5]*pow(theta,1) + CC[5]*pow(theta,2) + CD[5]*pow(theta,3) + CE[5]*pow(theta,4);
      pcorr = A * pow(phi,0) + B * pow(phi,1) + C * pow(phi,2);
    }
  }

  double mom_corrected =  (1 + (pcorr-1)*torratio) * mom;

  double px_corrected = mom_corrected * sin(theta*Pival/180) * cos(phi*Pival/180);
  double py_corrected = mom_corrected * sin(theta*Pival/180) * sin(phi*Pival/180);
  double pz_corrected = mom_corrected * cos(theta*Pival/180);


  TLorentzVector lep_neg_corr;
  lep_neg_corr.SetPxPyPzE(px_corrected, py_corrected, pz_corrected, sqrt(mom_corrected*mom_corrected + m_e*m_e)); 
  return lep_neg_corr;
}




TVector3 correct_hadron_positive(double thetahd, double phihd, double ph, int secth){

  double constant, lin, quad, third, fourth;
  double hconst, hlin, hquad, hthird, hfourth; 
  double hnew;
  double phcorr = ph; 
  double nthetahd = thetahd;

  double torcur = 3375.0;
  double tormax = 3375.0;
  double Bratio;

  double phih_real = phihd;

  double pex = phcorr * sin(nthetahd*Pival/180) * cos(phih_real*Pival/180);
  double pey = phcorr * sin(nthetahd*Pival/180) * sin(phih_real*Pival/180);
  double pez = phcorr * cos(nthetahd*Pival/180);

  TVector3 had_pos_corr;
  had_pos_corr.SetXYZ(pex, pey, pez); 
  return had_pos_corr;
}


// electron

void correct_electron(void){

  for(Int_t i = 0; i < BUFFER; i++){ 
    if(p4_ele[i].E() > 0){
      TLorentzVector ele_out = correct_lepton_negative(p4_ele[i]);
      p4_ele[i].SetPxPyPzE(ele_out.Px(), ele_out.Py(), ele_out.Pz(), ele_out.E()); 
    }
  }
}


// positive hadrons

void correct_proton(void){

  for(Int_t i = 0; i < BUFFER; i++){ 
    if(p4_prot[i].E() > 0){
      TVector3 prot_out = correct_hadron_positive(p4_prot[i].Theta()*180/Pival, p4_prot[i].Phi()*180/Pival, p4_prot[i].P(), part_DC_sector[i]);
      p4_prot[i].SetPxPyPzE(prot_out.Px(), prot_out.Py(), prot_out.Pz(), sqrt(prot_out.Mag2() + m_p*m_p)); 
    }
  }
}

void correct_pip(void){

  for(Int_t i = 0; i < BUFFER; i++){ 
    if(p4_pip[i].E() > 0){
      TVector3 pip_out = correct_hadron_positive(p4_pip[i].Theta()*180/Pival, p4_pip[i].Phi()*180/Pival, p4_pip[i].P(), part_DC_sector[i]);
      p4_pip[i].SetPxPyPzE(pip_out.Px(), pip_out.Py(), pip_out.Pz(), sqrt(pip_out.Mag2() + m_pip*m_pip)); 
    }
  }
}

void correct_Kplus(void){

  for(Int_t i = 0; i < BUFFER; i++){ 
    if(p4_pim[i].E() > 0){
      TVector3 Kp_out = correct_hadron_positive(p4_Kp[i].Theta()*180/Pival, p4_Kp[i].Phi()*180/Pival, p4_Kp[i].P(), part_DC_sector[i]);
      p4_Kp[i].SetPxPyPzE(Kp_out.Px(), Kp_out.Py(), Kp_out.Pz(), sqrt(Kp_out.Mag2() + m_Kp*m_Kp)); 
    }
  }
}


// negative hadrons:

void correct_pim(void){

}

void correct_Kminus(void){

}


// neutral hadrons:

void correct_neutron(void){

}


// photons:

void correct_photon(void){

}




/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// fill output tree variables (particles from all detrectors are combined):
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void fill_output_vector_electron(void){
  for(int i = 0; i < BUFFER; i++){
    if(p4_ele[i].E() > 0){
      p4_ele_px.push_back(p4_ele[i].Px());
      p4_ele_py.push_back(p4_ele[i].Py());
      p4_ele_pz.push_back(p4_ele[i].Pz());
      p4_ele_E.push_back(p4_ele[i].E()); 
      ele_det.push_back(ele_detect[i]);      
      sectorE.push_back(e_PCAL_sec[i]);
      electron_event_number.push_back(eventNumber[i]);

    } 
  }
}

void fill_output_vector_proton(void){
  for(int i = 0; i < BUFFER; i++){
    if(p4_prot[i].E() > 0){
      p4_prot_px.push_back(p4_prot[i].Px());
      p4_prot_py.push_back(p4_prot[i].Py());
      p4_prot_pz.push_back(p4_prot[i].Pz());
      p4_prot_E.push_back(p4_prot[i].E());  
      prot_det.push_back(prot_detect[i]); 
    }
  } 
}

void fill_output_vector_neutron(void){
  for(int i = 0; i < BUFFER; i++){
    if(p4_neutr[i].E() > 0){
      p4_neutr_px.push_back(p4_neutr[i].Px());
      p4_neutr_py.push_back(p4_neutr[i].Py());
      p4_neutr_pz.push_back(p4_neutr[i].Pz());
      p4_neutr_E.push_back(p4_neutr[i].E()); 
      neutr_det.push_back(neutr_detect[i]);  
    }
  }
}

void fill_output_vector_pip(void){
  for(int i = 0; i < BUFFER; i++){
    if(p4_pip[i].E() > 0){
      p4_pip_px.push_back(p4_pip[i].Px());
      p4_pip_py.push_back(p4_pip[i].Py());
      p4_pip_pz.push_back(p4_pip[i].Pz());
      p4_pip_E.push_back(p4_pip[i].E()); 
      pip_det.push_back(pip_detect[i]); 
    } 
  }
}

void fill_output_vector_pim(void){
  for(int i = 0; i < BUFFER; i++){
    if(p4_pim[i].E() > 0){
      p4_pim_px.push_back(p4_pim[i].Px());
      p4_pim_py.push_back(p4_pim[i].Py());
      p4_pim_pz.push_back(p4_pim[i].Pz());
      p4_pim_E.push_back(p4_pim[i].E());  
      pim_det.push_back(pim_detect[i]); 
    }
  }
}

void fill_output_vector_Kplus(void){
  for(int i = 0; i < BUFFER; i++){
    if(p4_Kp[i].E() > 0){
      p4_Kp_px.push_back(p4_Kp[i].Px());
      p4_Kp_py.push_back(p4_Kp[i].Py());
      p4_Kp_pz.push_back(p4_Kp[i].Pz());
      p4_Kp_E.push_back(p4_Kp[i].E()); 
      Kp_det.push_back(Kp_detect[i]); 
    } 
  }
}

void fill_output_vector_Kminus(void){
  for(int i = 0; i < BUFFER; i++){
    if(p4_Km[i].E() > 0){
      p4_Km_px.push_back(p4_Km[i].Px());
      p4_Km_py.push_back(p4_Km[i].Py());
      p4_Km_pz.push_back(p4_Km[i].Pz());
      p4_Km_E.push_back(p4_Km[i].E()); 
      Km_det.push_back(Km_detect[i]);  
    }
  }
}

void fill_output_vector_photon(void){
  for(int i = 0; i < BUFFER; i++){
    if(p4_phot[i].E() > 0){
      p4_phot_px.push_back(p4_phot[i].Px());
      p4_phot_py.push_back(p4_phot[i].Py());
      p4_phot_pz.push_back(p4_phot[i].Pz());
      p4_phot_E.push_back(p4_phot[i].E());  
      phot_det.push_back(phot_detect[i]); 
    }
  }
}



/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Additional functions:
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double GetTheta(int j){
    return acos(part_pz[j]/part_p[j]);
}

double GetPhi(int j){
  return 150*Pival/180 + atan2(part_py[j]/part_p[j], part_px[j]/part_p[j]);  
}

TVector3 GetUVWVector(int j){

  double u, v, w, xi, yi, zi;
  double EC_the = 0.4363323;
  double EC_phi;
  double ylow   = -182.974;
  double yhi    = 189.956;
  double tgrho  = 1.95325;
  double sinrho = 0.8901256;
  double cosrho = 0.455715;
  double rot[3][3];
    
  double x = part_Cal_PCAL_x[j];
  double y = part_Cal_PCAL_y[j];
  double z = part_Cal_PCAL_z[j];
    
  double at = atan2(y, x);
  if (at < 0) at += 2*Pival;
      
  double phi = at*180/Pival;
  phi=phi+30.;
  if (phi>=360.) phi=phi-360.;
    
  EC_phi = (int)(phi/60.) * 1.0471975;
    
  rot[0][0] = cos(EC_the)*cos(EC_phi);
  rot[0][1] = -sin(EC_phi);
  rot[0][2] = sin(EC_the)*cos(EC_phi);
  rot[1][0] = cos(EC_the)*sin(EC_phi);
  rot[1][1] = cos(EC_phi);
  rot[1][2] = sin(EC_the)*sin(EC_phi);
  rot[2][0] = -sin(EC_the);
  rot[2][1] = 0.;
  rot[2][2] = cos(EC_the);
    
  yi=x*rot[0][0]+y*rot[1][0]+z*rot[2][0];
  xi=x*rot[0][1]+y*rot[1][1]+z*rot[2][1];
  zi=x*rot[0][2]+y*rot[1][2]+z*rot[2][2];
  
  zi=zi-510.32;
    
  u = (yi-ylow)/sinrho;
  v = (yhi-ylow)/tgrho-xi+(yhi-yi)/tgrho;
  w = ((yhi-ylow)/tgrho+xi+(yhi-yi)/tgrho)/2./cosrho;
   
  TVector3 uvw(u,v,w);
  return uvw;
}


double GetTOFmass2(int j, int run)
{
  double tofmass2 = 0; 
  if(Beta_charged(j, run) > 0) tofmass2 = pow(part_p[j],2) * (1-pow(Beta_charged(j, run),2))/pow(Beta_charged(j, run),2);
  return tofmass2;
}

double GetTOFmass2_CD(int j, int run)
{
  double tofmass2 = 0; 
  if(Beta_charged_central(j, run) > 0) tofmass2 = pow(part_p[j],2) * (1-pow(Beta_charged_central(j, run),2))/pow(Beta_charged_central(j, run),2);
  return tofmass2;
}

double Getdvz(int j)
{
  // find index of electron with highest momentum as best vertex approach
  float mom = 0;
  int eind = -1;
  for(int k = 0; k < BUFFER; k++){
    if(FD_eid_all_check[k]){
      if(part_p[k] > mom){        
	mom = part_p[k]; 
        eind = k;
      }
    }  
  }

  if(eind == -1) return -1000;  // if no electron exists

  double dvz = 0; 
  dvz = part_vz[eind] - part_vz[j];
  return dvz;
}


double Get_Starttime(int j, int run)
{

  // find index of electron with highest momentum in FD or FT as best start time approach

  float mom_FD = 0;
  float mom_FT = 0;
  int eind_FD = -1;
  int eind_FT = -1;

  for(int k = 0; k < BUFFER; k++){
    if(FD_eid_all_check[k] && part_p[k] > mom_FD && part_FTOF_time[k] > 0){        
      mom_FD = part_p[k]; 
      eind_FD = k;
    } 
  }

  for(int k = 0; k < BUFFER; k++){
    if(FT_eid_all_check[k] && part_p[k] > mom_FT && part_FTHODO_time[k] > 0){        
      mom_FT = part_p[k]; 
      eind_FT = k;
    }
  }

  if(eind_FD == -1 && eind_FT == -1) return 0;    // no electron found 
  if(good_sc_paddle(eind_FD) == false)  return 0;                      // reject bad SC paddles 

  double startTime = 0;

  if(part_p[eind_FD] >= part_p[eind_FT] && electron_sct_corr(eind_FD, run) > 0 && part_FTOF_path[eind_FD] > 0){
    startTime = electron_sct_corr(eind_FD, run) - 10000000*part_FTOF_path[eind_FD]/c;                          // take electron in FD for start time calculation
  }
  else if(electron_sct_corr(eind_FT, run) > 0 && part_FTHODO_path[eind_FT] > 0){
    startTime = electron_sct_corr(eind_FT, run) - 10000000*part_FTHODO_path[eind_FT]/c;                        // take electron in FT for start time calculation
  }

  return startTime;

}


double Beta_charged(int j, int run){

  if(use_own_beta_charged == false){ return part_beta[j];}

  // if own beta should be used, ....

  double startTime = Get_Starttime(j, run);

  if(startTime <= 0 || (hadron_sct_corr(j, run)-startTime) <= 0 || hadron_sct_corr(j, run) <= 0) return 0;    // consitency check
  if(good_sc_paddle(j) == false)  return 0;                                   			  // reject bad SC paddles 

  return 10000000*(part_FTOF_path[j]/(hadron_sct_corr(j, run)-startTime))/c;   // for charged particles beta is calculated based on FTOF time

}


double Beta_charged_central(int j, int run){

  if(use_own_beta_charged == false){ return part_beta[j];}

  // if own beta should be used, ....

  double startTime = Get_Starttime(j, run);

  if(startTime <= 0 || (hadron_sct_corr_central(j, run)-startTime) <= 0 || hadron_sct_corr_central(j, run) <= 0) return 0;    // consitency check

  return 10000000*(part_CTOF_path[j]/(hadron_sct_corr_central(j, run)-startTime))/c;   // for charged particles beta is calculated based on FTOF time
}


double Beta_charged_FT(int j, int run){

  if(use_own_beta_charged_FT == false){ return part_beta[j];}

  // if own beta should be used, ....

  double startTime = Get_Starttime(j, run);

  if(startTime <= 0 || (part_FTOF_time[j]-startTime) <= 0 || part_FTOF_time[j] <= 0) return 0;    // consitency check

  return 10000000*(part_FTHODO_path[j]/(part_FTHODO_time[j]-startTime))/c;   // for charged particles beta is calculated based on FTOF time

}


double Beta_neutral(int j, int run){

  if(use_own_beta_neutrals == false){ return part_beta[j];}

  // if own beta should be used, ....

  double startTime = Get_Starttime(j, run);

  if(startTime <= 0 || (part_FTOF_time[j]-startTime) <= 0 || part_FTOF_time[j] <= 0) return 0; // consitency check

  return 10000000*(part_Cal_PCAL_path[j]/(part_Cal_PCAL_time[j]-startTime))/c;     // for neutrals beta is calculated based on EC time

}


double Beta_neutral_FT(int j, int run){

  if(use_own_beta_neutrals_FT == false){ return part_beta[j];}

  // if own beta should be used, ....

  double startTime = Get_Starttime(j, run);

  if(startTime <= 0 || (part_FTOF_time[j]-startTime) <= 0 || part_FTOF_time[j] <= 0) return 0; // consitency check

  return 10000000*(part_FTCAL_path[j]/(part_FTCAL_time[j]-startTime))/c;     // for neutrals beta is calculated based on EC time

}


double electron_sct_corr(int j, int run){   // SC timing correction for electrons

  int component_layer1 = part_FTOF_component_layer1[j]; 
  int component_layer2 = part_FTOF_component_layer2[j]; 
  int component_layer3 = part_FTOF_component_layer3[j]; 

  int sector_layer1 = part_FTOF_sector_layer1[j];
  int sector_layer2 = part_FTOF_sector_layer2[j];
  int sector_layer3 = part_FTOF_sector_layer3[j];;

  if(simulation == true) return part_FTOF_time[j];   // do not correct simulated data

  //if(sector == 1 && paddle == 1 && run >= 1 && run <= 100)	return sc_t[j] + 0;

  return part_FTOF_time[j];

}

double hadron_sct_corr(int j, int run){    // SC timing correction for hadrons

  int component_layer1 = part_FTOF_component_layer1[j]; 
  int component_layer2 = part_FTOF_component_layer2[j]; 
  int component_layer3 = part_FTOF_component_layer3[j]; 

  int sector_layer1 = part_FTOF_sector_layer1[j];
  int sector_layer2 = part_FTOF_sector_layer2[j];
  int sector_layer3 = part_FTOF_sector_layer3[j];

  if(simulation == true)    return part_FTOF_time[j];   // do not correct simulated data
  if(correct_FTOF == false) return part_FTOF_time[j];  

  /// ////////////////////////////////////////////////////////////////////////////////////////////////

  double mean_pion_combined_layer2_sec1[] = { 0.0538753, 0.0894172, 0.0415536, -0.00697086, -0.0256552, 0.0585226, 0.0387975, 0.0168477, 0.0366831, 0.0361079, 0.0315495, 0.0633489, 0.0284496, 0.0046004, 0.0237375, 0.0372987, -0.00869811, 0.0494778, -0.0175471, 0.00121505, 0.00821949, -0.0110225, -0.00702778, -0.00843307, 0.010757, 0.02797, 0.0150352, 0.0127725, 0.00137534, 0.0266265, 0.0312262, 0.0184115, 0.0207408, 0.0356855, 0.0338395, 0.0168999, 0.0441921, 0.0339078, 0.0563256, 0.0483275, 0.0406767, 0.0529402, 0.0510855, 0.0749354, 0.0601476, 0.067262, 0.0633251, 0.0522541, 0.0290266, 0.0843323, 0.0851452, 0.0887319, 0.113017, 0.127794, 0.154692, 0.179721, 0.173495, 0.14365, 0.831491, 0.207532, 0.180673, 0.214689};
  double mean_pion_combined_layer2_sec2[] = { 0.0201007, 0.0465524, 0.0044113, 0.0725965, 0.00931767, 0.0573054, 0.0335023, 0.0177905, 0.0191517, 0.0143014, -0.0506109, 0.0320215, 0.0380263, -0.00331458, 0.00427068, -0.00365863, 0.0156208, 0.0215388, -0.0162937, -0.0260334, -0.0254874, -0.0147608, -0.00136278, 0.00224445, -0.00325317, 0.0227318, -0.000339176, 0.0350471, 0.00436642, 0.0225317, 0.00832329, 0.00801539, 0.024013, 0.0254705, 0.0280339, 0.0326556, 0.0485061, 0.0319243, 0.0285167, 0.0338927, 0.0156294, 0.0347015, 0.0265838, 0.0463763, 0.0337154, 0.0563934, 0.0207176, 0.0624435, 0.0781742, 0.0581556, 0.0817853, 0.0851641, 0.0965995, 0.122388, 0.107749, 0.156704, 0.156899, 0.171214, 0.224282, 0.130809, 0.201327, 0.229981};
  double mean_pion_combined_layer2_sec3[] = { 0.0288409, 0.022084, 0.0568596, 0.0164514, 0.0803507, 0.0404029, 0.0266618, 0.0259344, 0.0864902, -0.00314741, 0.0535332, 0.0232373, 0.00937413, 0.0100241, 0.0214095, 0.00648485, 0.0223418, 0.0260699, 0.0149715, -0.0157239, 0.0147615, -0.011228, -0.0270731, -0.0101031, -0.00793331, 0.0240069, -0.00444922, 0.00261553, -0.00149125, 0.0129531, -0.0010735, 0.00610513, 0.0219764, 0.0100616, 0.000634229, 0.0210686, 0.0186737, 0.0305352, 0.0306624, 0.0225514, 0.0371804, 0.0462444, 0.0261995, 0.031015, 0.0542188, 0.0755648, 0.044206, 0.0492144, 0.0562242, 0.0484466, 0.07956, 0.0891106, 0.0909857, 0.109619, 0.107927, 0.136553, 0.121198, 0.165602, 0.206659, 0.160572, 0.208224, 0.203239};
  double mean_pion_combined_layer2_sec4[] = { -0.0174365, 0.0403755, 0.0621618, 0.100543, 0.0252135, 0.0639297, 0.0314792, 0.0248234, -0.000119261, -0.00698425, 0.0520655, 0.0353832, 0.0372486, 0.0339803, 0.0492965, 0.0258162, 0.041965, 0.0288217, 0.0112785, 0.00615734, 0.012968, -0.0189003, -0.00886905, -0.00585084, -0.0257223, 0.00263944, -0.012168, 0.00739494, -0.0129347, -0.0110806, -0.00072699, 0.00725338, 0.00161277, 0.00947743, 0.0186749, 0.0146464, 0.0207656, 0.014423, 0.0100822, 0.00148744, 0.0338362, 0.0285981, 0.036646, 0.0279241, 0.023769, 0.0466241, 0.0404707, 0.0484161, 0.0370252, 0.0335709, 0.0580888, 0.0806337, 0.0880605, 0.100155, 0.0715631, 0.096832, 0.161302, 0.102406, 0.174471, 0.167655, 0.166398, 0.156839};
  double mean_pion_combined_layer2_sec5[] = { 0.0317948, 0.0565442, 0.0326173, 0.0385375, -0.0237504, 0.0196178, -0.0136122, -0.00860579, 0.0210756, 0.0677887, 0.0118137, -0.00975603, -0.0156249, 0.022647, 0.00889069, -0.01372, 0.00266827, -0.00832111, -0.036579, -0.0109632, 0.00581234, 0.0156149, -0.00328594, 0.00637419, 0.00501017, 0.0153823, 0.0166495, 0.0298507, 0.0303831, 0.0304541, 0.0205921, 0.0361692, 0.035899, 0.0244119, 0.0397405, 0.00487128, 0.0321951, 0.0392289, 0.0478879, 0.0388485, 0.027897, 0.0281919, 0.0155273, 0.0541479, 0.0535899, 0.0762408, 0.0408552, 0.0731012, 0.0955377, 0.087898, 0.102133, 0.111963, 0.095632, 0.1142, 0.151744, 0.179316, 0.184459, 0.215345, 0.228311, 0.217565, 0.301694, 0.276};
  double mean_pion_combined_layer2_sec6[] = { 0.0250528, 0.0435939, 0.157377, 0.08319, 0.0595912, 0.0522503, 0.0558913, 0.0342482, -0.00132545, 0.0464159, -0.00551834, 0.0319466, 0.0497912, 0.0603297, 0.0537583, 0.0549158, 0.00581518, 0.032842, -0.0204257, 0.0380622, 0.000446151, 0.0238682, -0.00408773, 0.00765562, -0.00188372, 0.0219118, 0.0250298, 0.0287166, 0.0433487, 0.0199701, 0.0500815, 0.0331643, 0.0589243, 0.0305802, 0.0493649, 0.0407388, 0.027975, 0.0417421, 0.047593, 0.0701573, 0.0502872, 0.0571604, 0.0407683, 0.0688008, 0.0854068, 0.06495, 0.0856892, 0.0758915, 0.0984645, 0.0791596, 0.107556, 0.130177, 0.117903, 0.14161, 0.121797, 0.154971, 0.149244, 0.194147, 0.234898, 0.184855, 0.215379, 0.239638};


  /// //////////////////////////////////////////////////////////////////////////////////////////////////

  if(sector_layer2 == 1){
    for(int k = 0; k < 62; k++){
      if(component_layer2 == k+1){ return part_FTOF_time[j] - mean_pion_combined_layer2_sec1[k];}
    }
  }
  if(sector_layer2 == 2){
    for(int k = 0; k < 62; k++){
      if(component_layer2 == k+1){ return part_FTOF_time[j] - mean_pion_combined_layer2_sec2[k];}
    }
  }
  if(sector_layer2 == 3){
    for(int k = 0; k < 62; k++){
      if(component_layer2 == k+1){ return part_FTOF_time[j] - mean_pion_combined_layer2_sec3[k];}
    }
  }
  if(sector_layer2 == 4){
    for(int k = 0; k < 62; k++){
      if(component_layer2 == k+1){ return part_FTOF_time[j] - mean_pion_combined_layer2_sec4[k];}
    }
  }
  if(sector_layer2 == 5){
    for(int k = 0; k < 62; k++){
      if(component_layer2 == k+1){ return part_FTOF_time[j] - mean_pion_combined_layer2_sec5[k];}
    }
  }
  if(sector_layer2 == 6){
    for(int k = 0; k < 62; k++){
      if(component_layer2 == k+1){ return part_FTOF_time[j] - mean_pion_combined_layer2_sec6[k];}
    }
  }

  return part_FTOF_time[j];

}


double hadron_sct_corr_central(int j, int run){    // SC timing correction for hadrons

  int component_CTOF = part_CTOF_component[j];

  if(simulation == true) return part_CTOF_time[j];   // do not correct simulated data
  if(correct_CTOF == false) return part_CTOF_time[j];  
    
  double mean_pion_combined_CTOF[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  for(int k = 0; k < 48; k++){
    if(component_CTOF == k+1){
      return part_CTOF_time[j] - mean_pion_combined_CTOF[k]; 
    }
  }

  return part_CTOF_time[j];
}


bool good_sc_paddle(int j){    // reject bad paddles

  int component_layer1 = part_FTOF_component_layer1[j]; 
  int component_layer2 = part_FTOF_component_layer2[j]; 
  int component_layer3 = part_FTOF_component_layer3[j]; 

  int sector_layer1 = part_FTOF_sector_layer1[j];
  int sector_layer2 = part_FTOF_sector_layer2[j];
  int sector_layer3 = part_FTOF_sector_layer3[j];

  int component_CTOF = part_CTOF_component[j];

  if(sector_layer1 == 2 && component_layer1 == 6)   return false;
  if(sector_layer1 == 2 && component_layer1 == 10)  return false;
  if(sector_layer1 == 2 && component_layer1 == 19)  return false;
  if(sector_layer1 == 6 && component_layer1 == 21)  return false;

  if(sector_layer3 == 1 && component_layer3 == 5)   return false;
  if(sector_layer3 == 2 && component_layer3 == 5)   return false;
  if(sector_layer3 == 6 && component_layer3 == 3)   return false;

  if(component_CTOF == 14)                          return false;

  return true;

}


/// ////////////////////////////////////////////////////
/// Kinematics:

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
  v1 = gammaClone.Vect().Cross(hadronClone.Vect());  // hadronClone.Vect().Cross(gammaClone.Vect());

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


/// ///////////////////////////////////////////////////
/// missing particle and mass:

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



/// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// reconstruct neutrals from gammas:
/// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void create_neutrals(void){

float g_thr = 0.1;

if(p4_phot[0].E() > g_thr && p4_phot[1].E() > g_thr){ neutral_iter[0] = p4_phot[0]+p4_phot[1]; pi0_g_ind[0] = 0; alpha_gg[0] = p4_phot[0].Vect().Angle(p4_phot[1].Vect()); neutral_ind[0][0]=g_sec[0]; neutral_ind[0][1]=g_sec[1];}
if(p4_phot[0].E() > g_thr && p4_phot[2].E() > g_thr){ neutral_iter[1] = p4_phot[0]+p4_phot[2]; pi0_g_ind[1] = 0; alpha_gg[1] = p4_phot[0].Vect().Angle(p4_phot[2].Vect()); neutral_ind[1][0]=g_sec[0]; neutral_ind[1][1]=g_sec[2];}
if(p4_phot[1].E() > g_thr && p4_phot[2].E() > g_thr){ neutral_iter[2] = p4_phot[1]+p4_phot[2]; pi0_g_ind[2] = 1; alpha_gg[2] = p4_phot[1].Vect().Angle(p4_phot[2].Vect()); neutral_ind[2][0]=g_sec[1]; neutral_ind[2][1]=g_sec[2];}
if(p4_phot[0].E() > g_thr && p4_phot[3].E() > g_thr){ neutral_iter[3] = p4_phot[0]+p4_phot[3]; pi0_g_ind[3] = 0; alpha_gg[3] = p4_phot[0].Vect().Angle(p4_phot[3].Vect()); neutral_ind[3][0]=g_sec[0]; neutral_ind[3][1]=g_sec[3];}
if(p4_phot[1].E() > g_thr && p4_phot[3].E() > g_thr){ neutral_iter[4] = p4_phot[1]+p4_phot[3]; pi0_g_ind[4] = 1; alpha_gg[4] = p4_phot[1].Vect().Angle(p4_phot[3].Vect()); neutral_ind[4][0]=g_sec[1]; neutral_ind[4][1]=g_sec[3];}
if(p4_phot[2].E() > g_thr && p4_phot[3].E() > g_thr){ neutral_iter[5] = p4_phot[2]+p4_phot[3]; pi0_g_ind[5] = 2; alpha_gg[5] = p4_phot[2].Vect().Angle(p4_phot[3].Vect()); neutral_ind[5][0]=g_sec[2]; neutral_ind[5][1]=g_sec[3];}
if(p4_phot[0].E() > g_thr && p4_phot[4].E() > g_thr){ neutral_iter[6] = p4_phot[0]+p4_phot[4]; pi0_g_ind[6] = 0; alpha_gg[6] = p4_phot[0].Vect().Angle(p4_phot[4].Vect()); neutral_ind[6][0]=g_sec[0]; neutral_ind[6][1]=g_sec[4];}
if(p4_phot[1].E() > g_thr && p4_phot[4].E() > g_thr){ neutral_iter[7] = p4_phot[1]+p4_phot[4]; pi0_g_ind[7] = 1; alpha_gg[7] = p4_phot[1].Vect().Angle(p4_phot[4].Vect()); neutral_ind[7][0]=g_sec[1]; neutral_ind[7][1]=g_sec[4];}
if(p4_phot[2].E() > g_thr && p4_phot[4].E() > g_thr){ neutral_iter[8] = p4_phot[2]+p4_phot[4]; pi0_g_ind[8] = 2; alpha_gg[8] = p4_phot[2].Vect().Angle(p4_phot[4].Vect()); neutral_ind[8][0]=g_sec[2]; neutral_ind[8][1]=g_sec[4];}
if(p4_phot[3].E() > g_thr && p4_phot[4].E() > g_thr){ neutral_iter[9] = p4_phot[3]+p4_phot[4]; pi0_g_ind[9] = 3; alpha_gg[9] = p4_phot[3].Vect().Angle(p4_phot[4].Vect()); neutral_ind[9][0]=g_sec[3]; neutral_ind[9][1]=g_sec[4];}
if(p4_phot[0].E() > g_thr && p4_phot[5].E() > g_thr){ neutral_iter[10] = p4_phot[0]+p4_phot[5]; pi0_g_ind[10] = 0; alpha_gg[10] = p4_phot[0].Vect().Angle(p4_phot[5].Vect());neutral_ind[10][0]=g_sec[0]; neutral_ind[10][1]=g_sec[5];}
if(p4_phot[1].E() > g_thr && p4_phot[5].E() > g_thr){ neutral_iter[11] = p4_phot[1]+p4_phot[5]; pi0_g_ind[11] = 1; alpha_gg[11] = p4_phot[1].Vect().Angle(p4_phot[5].Vect());neutral_ind[11][0]=g_sec[1]; neutral_ind[11][1]=g_sec[5];}
if(p4_phot[2].E() > g_thr && p4_phot[5].E() > g_thr){ neutral_iter[12] = p4_phot[2]+p4_phot[5]; pi0_g_ind[12] = 2; alpha_gg[12] = p4_phot[2].Vect().Angle(p4_phot[5].Vect());neutral_ind[12][0]=g_sec[2]; neutral_ind[12][1]=g_sec[5];}
if(p4_phot[3].E() > g_thr && p4_phot[5].E() > g_thr){ neutral_iter[13] = p4_phot[3]+p4_phot[5]; pi0_g_ind[13] = 3; alpha_gg[13] = p4_phot[3].Vect().Angle(p4_phot[5].Vect());neutral_ind[13][0]=g_sec[3]; neutral_ind[13][1]=g_sec[5];}
if(p4_phot[4].E() > g_thr && p4_phot[5].E() > g_thr){ neutral_iter[14] = p4_phot[4]+p4_phot[5]; pi0_g_ind[14] = 4; alpha_gg[14] = p4_phot[4].Vect().Angle(p4_phot[5].Vect());neutral_ind[14][0]=g_sec[4]; neutral_ind[14][1]=g_sec[5];}
if(p4_phot[0].E() > g_thr && p4_phot[6].E() > g_thr){ neutral_iter[15] = p4_phot[0]+p4_phot[6]; pi0_g_ind[15] = 0; alpha_gg[15] = p4_phot[0].Vect().Angle(p4_phot[6].Vect());neutral_ind[15][0]=g_sec[0]; neutral_ind[15][1]=g_sec[6];}
if(p4_phot[1].E() > g_thr && p4_phot[6].E() > g_thr){ neutral_iter[16] = p4_phot[1]+p4_phot[6]; pi0_g_ind[16] = 1; alpha_gg[16] = p4_phot[1].Vect().Angle(p4_phot[6].Vect());neutral_ind[16][0]=g_sec[1]; neutral_ind[16][1]=g_sec[6];}
if(p4_phot[2].E() > g_thr && p4_phot[6].E() > g_thr){ neutral_iter[17] = p4_phot[2]+p4_phot[6]; pi0_g_ind[17] = 2; alpha_gg[17] = p4_phot[2].Vect().Angle(p4_phot[6].Vect());neutral_ind[17][0]=g_sec[2]; neutral_ind[17][1]=g_sec[6];}
if(p4_phot[3].E() > g_thr && p4_phot[6].E() > g_thr){ neutral_iter[18] = p4_phot[3]+p4_phot[6]; pi0_g_ind[18] = 3; alpha_gg[18] = p4_phot[3].Vect().Angle(p4_phot[6].Vect());neutral_ind[18][0]=g_sec[3]; neutral_ind[18][1]=g_sec[6];}
if(p4_phot[4].E() > g_thr && p4_phot[6].E() > g_thr){ neutral_iter[19] = p4_phot[4]+p4_phot[6]; pi0_g_ind[19] = 4; alpha_gg[19] = p4_phot[4].Vect().Angle(p4_phot[6].Vect());neutral_ind[19][0]=g_sec[4]; neutral_ind[19][1]=g_sec[6];}
if(p4_phot[5].E() > g_thr && p4_phot[6].E() > g_thr){ neutral_iter[20] = p4_phot[5]+p4_phot[6]; pi0_g_ind[20] = 5; alpha_gg[20] = p4_phot[5].Vect().Angle(p4_phot[6].Vect());neutral_ind[20][0]=g_sec[5]; neutral_ind[20][1]=g_sec[6];}
if(p4_phot[0].E() > g_thr && p4_phot[7].E() > g_thr){ neutral_iter[21] = p4_phot[0]+p4_phot[7]; pi0_g_ind[21] = 0; alpha_gg[21] = p4_phot[0].Vect().Angle(p4_phot[7].Vect());neutral_ind[21][0]=g_sec[0]; neutral_ind[21][1]=g_sec[7];}
if(p4_phot[1].E() > g_thr && p4_phot[7].E() > g_thr){ neutral_iter[22] = p4_phot[1]+p4_phot[7]; pi0_g_ind[22] = 1; alpha_gg[22] = p4_phot[1].Vect().Angle(p4_phot[7].Vect());neutral_ind[22][0]=g_sec[1]; neutral_ind[22][1]=g_sec[7];}
if(p4_phot[2].E() > g_thr && p4_phot[7].E() > g_thr){ neutral_iter[23] = p4_phot[2]+p4_phot[7]; pi0_g_ind[23] = 2; alpha_gg[23] = p4_phot[2].Vect().Angle(p4_phot[7].Vect());neutral_ind[23][0]=g_sec[2]; neutral_ind[23][1]=g_sec[7];}
if(p4_phot[3].E() > g_thr && p4_phot[7].E() > g_thr){ neutral_iter[24] = p4_phot[3]+p4_phot[7]; pi0_g_ind[24] = 3; alpha_gg[24] = p4_phot[3].Vect().Angle(p4_phot[7].Vect());neutral_ind[24][0]=g_sec[3]; neutral_ind[24][1]=g_sec[7];}
if(p4_phot[4].E() > g_thr && p4_phot[7].E() > g_thr){ neutral_iter[25] = p4_phot[4]+p4_phot[7]; pi0_g_ind[25] = 4; alpha_gg[25] = p4_phot[4].Vect().Angle(p4_phot[7].Vect());neutral_ind[25][0]=g_sec[4]; neutral_ind[25][1]=g_sec[7];}
if(p4_phot[5].E() > g_thr && p4_phot[7].E() > g_thr){ neutral_iter[26] = p4_phot[5]+p4_phot[7]; pi0_g_ind[26] = 5; alpha_gg[26] = p4_phot[5].Vect().Angle(p4_phot[7].Vect());neutral_ind[26][0]=g_sec[5]; neutral_ind[26][1]=g_sec[7];}
if(p4_phot[6].E() > g_thr && p4_phot[7].E() > g_thr){ neutral_iter[27] = p4_phot[6]+p4_phot[7]; pi0_g_ind[27] = 6; alpha_gg[27] = p4_phot[6].Vect().Angle(p4_phot[7].Vect());neutral_ind[27][0]=g_sec[6]; neutral_ind[27][1]=g_sec[7];}


  for(Int_t i = 0; i <= 27; i++){
    mass_neutral_iter2[i] = neutral_iter[i].E()*neutral_iter[i].E() - neutral_iter[i].Px()*neutral_iter[i].Px() - neutral_iter[i].Py()*neutral_iter[i].Py() - neutral_iter[i].Pz()*neutral_iter[i].Pz();
    mass_neutral_iter[i] = sqrt(mass_neutral_iter2[i]);
  }

  // select the pi0 and eta from the two gamma with the highest momentum which fulfill the condition

  select_pi0 = 0;
  select_eta = 0;

  if(g_count >= 2){
    for(Int_t i = 27; i >= 0; i--){ 

      //if(mass_neutral_iter[i] > 0.05 &&  mass_neutral_iter[i] < 0.218 && mass_neutral_iter2[i] > 0.008 && mass_neutral_iter2[i] < 0.03 && alpha_gg[i]*180/Pival > 1.5 && neutral_iter[i].E() > 0.5 && alpha_gg[i]*180/Pival > (7.0 - 1.5 * neutral_iter[i].E())){

      if(alpha_gg[i]*180/Pival > 1.5 && neutral_iter[i].E() > 0.1 && alpha_gg[i]*180/Pival > (7.0 - 1.5 * neutral_iter[i].E())){
        select_pi0 = i;
        pi0_ind[0] = neutral_ind[i][0];
        pi0_ind[1] = neutral_ind[i][1];
      }
      else{select_pi0 = -1;}
      if(mass_neutral_iter[i] > 0.35 &&  mass_neutral_iter[i] < 0.75 && alpha_gg[i]*180/Pival > 1.5 && neutral_iter[i].E() > 0.5 && alpha_gg[i]*180/Pival > (7.0 - 1.5 * neutral_iter[i].E())){
        select_eta = i;
        eta_ind[0] = neutral_ind[i][0];
        eta_ind[1] = neutral_ind[i][1];
      }
      else{select_eta = -1;}
    }
  }

  if(select_pi0 != -1) pi0 = neutral_iter[select_pi0];
  if(select_eta != -1) eta = neutral_iter[select_eta];

}















