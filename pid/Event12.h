//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Nov 13 09:05:28 2018 by ROOT version 6.14/04
// from TChain clas12/
//////////////////////////////////////////////////////////

#ifndef Event12_h
#define Event12_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class Event12 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<float>   *RUN_config_torus;
   vector<int>     *RUN_config_run;
   vector<int>     *RUN_config_event;
   vector<float>   *RUN_config_solenoid;
   vector<int>     *RAW_scaler_crate;
   vector<int>     *RAW_scaler_slot;
   vector<int>     *RAW_scaler_channel;
   vector<int>     *RAW_scaler_helicity;
   vector<int>     *RAW_scaler_quartet;
   vector<int>     *RAW_scaler_value;
   vector<int>     *REC_Event_NRUN;
   vector<int>     *REC_Event_NEVENT;
   vector<float>   *REC_Event_EVNTime;
   vector<int>     *REC_Event_TYPE;
   vector<int>     *REC_Event_TRG;
   vector<float>   *REC_Event_BCG;
   vector<float>   *REC_Event_STTime;
   vector<float>   *REC_Event_RFTime;
   vector<int>     *REC_Event_Helic;
   vector<int>     *REC_Particle_pid;
   vector<float>   *REC_Particle_px;
   vector<float>   *REC_Particle_py;
   vector<float>   *REC_Particle_pz;
   vector<float>   *REC_Particle_vx;
   vector<float>   *REC_Particle_vy;
   vector<float>   *REC_Particle_vz;
   vector<int>     *REC_Particle_charge;
   vector<float>   *REC_Particle_beta;
   vector<float>   *REC_Particle_chi2pid;
   vector<int>     *REC_Particle_status;
   vector<int>     *REC_Calorimeter_pindex;
   vector<int>     *REC_Calorimeter_detector;
   vector<int>     *REC_Calorimeter_sector;
   vector<int>     *REC_Calorimeter_layer;
   vector<float>   *REC_Calorimeter_energy;
   vector<float>   *REC_Calorimeter_time;
   vector<float>   *REC_Calorimeter_path;
   vector<float>   *REC_Calorimeter_x;
   vector<float>   *REC_Calorimeter_y;
   vector<float>   *REC_Calorimeter_z;
   vector<float>   *REC_Calorimeter_lu;
   vector<float>   *REC_Calorimeter_lv;
   vector<float>   *REC_Calorimeter_lw;
   vector<int>     *REC_Cherenkov_pindex;
   vector<int>     *REC_Cherenkov_detector;
   vector<int>     *REC_Cherenkov_sector;
   vector<float>   *REC_Cherenkov_nphe;
   vector<float>   *REC_Cherenkov_time;
   vector<float>   *REC_Cherenkov_path;
   vector<float>   *REC_Cherenkov_theta;
   vector<float>   *REC_Cherenkov_phi;
   vector<int>     *REC_ForwardTagger_pindex;
   vector<int>     *REC_ForwardTagger_detector;
   vector<float>   *REC_ForwardTagger_energy;
   vector<float>   *REC_ForwardTagger_time;
   vector<float>   *REC_ForwardTagger_path;
   vector<float>   *REC_ForwardTagger_x;
   vector<float>   *REC_ForwardTagger_y;
   vector<float>   *REC_ForwardTagger_z;
   vector<float>   *REC_ForwardTagger_dx;
   vector<float>   *REC_ForwardTagger_dy;
   vector<float>   *REC_ForwardTagger_radius;
   vector<int>     *REC_ForwardTagger_size;
   vector<int>     *REC_Scintillator_pindex;
   vector<int>     *REC_Scintillator_detector;
   vector<int>     *REC_Scintillator_sector;
   vector<int>     *REC_Scintillator_layer;
   vector<int>     *REC_Scintillator_component;
   vector<float>   *REC_Scintillator_energy;
   vector<float>   *REC_Scintillator_time;
   vector<float>   *REC_Scintillator_path;
   vector<int>     *REC_Track_pindex;
   vector<int>     *REC_Track_detector;
   vector<int>     *REC_Track_sector;
   vector<float>   *REC_Track_chi2;
   vector<int>     *REC_Track_NDF;
   vector<float>   *REC_Track_chi2_nomm;
   vector<int>     *REC_Track_NDF_nomm;
   vector<int>     *REC_Traj_pindex;
   vector<int>     *REC_Traj_detId;
   vector<float>   *REC_Traj_x;
   vector<float>   *REC_Traj_y;
   vector<float>   *REC_Traj_z;
   vector<float>   *REC_Traj_cx;
   vector<float>   *REC_Traj_cy;
   vector<float>   *REC_Traj_cz;

   // List of branches
   TBranch        *b_RUN_config_torus;   //!
   TBranch        *b_RUN_config_run;   //!
   TBranch        *b_RUN_config_event;   //!
   TBranch        *b_RUN_config_solenoid;   //!
   TBranch        *b_RAW_scaler_crate;   //!
   TBranch        *b_RAW_scaler_slot;   //!
   TBranch        *b_RAW_scaler_channel;   //!
   TBranch        *b_RAW_scaler_helicity;   //!
   TBranch        *b_RAW_scaler_quartet;   //!
   TBranch        *b_RAW_scaler_value;   //!
   TBranch        *b_REC_Event_NRUN;   //!
   TBranch        *b_REC_Event_NEVENT;   //!
   TBranch        *b_REC_Event_EVNTime;   //!
   TBranch        *b_REC_Event_TYPE;   //!
   TBranch        *b_REC_Event_TRG;   //!
   TBranch        *b_REC_Event_BCG;   //!
   TBranch        *b_REC_Event_STTime;   //!
   TBranch        *b_REC_Event_RFTime;   //!
   TBranch        *b_REC_Event_Helic;   //!
   TBranch        *b_REC_Particle_pid;   //!
   TBranch        *b_REC_Particle_px;   //!
   TBranch        *b_REC_Particle_py;   //!
   TBranch        *b_REC_Particle_pz;   //!
   TBranch        *b_REC_Particle_vx;   //!
   TBranch        *b_REC_Particle_vy;   //!
   TBranch        *b_REC_Particle_vz;   //!
   TBranch        *b_REC_Particle_charge;   //!
   TBranch        *b_REC_Particle_beta;   //!
   TBranch        *b_REC_Particle_chi2pid;   //!
   TBranch        *b_REC_Particle_status;   //!
   TBranch        *b_REC_Calorimeter_pindex;   //!
   TBranch        *b_REC_Calorimeter_detector;   //!
   TBranch        *b_REC_Calorimeter_sector;   //!
   TBranch        *b_REC_Calorimeter_layer;   //!
   TBranch        *b_REC_Calorimeter_energy;   //!
   TBranch        *b_REC_Calorimeter_time;   //!
   TBranch        *b_REC_Calorimeter_path;   //!
   TBranch        *b_REC_Calorimeter_x;   //!
   TBranch        *b_REC_Calorimeter_y;   //!
   TBranch        *b_REC_Calorimeter_z;   //!
   TBranch        *b_REC_Calorimeter_lu;   //!
   TBranch        *b_REC_Calorimeter_lv;   //!
   TBranch        *b_REC_Calorimeter_lw;   //!
   TBranch        *b_REC_Cherenkov_pindex;   //!
   TBranch        *b_REC_Cherenkov_detector;   //!
   TBranch        *b_REC_Cherenkov_sector;   //!
   TBranch        *b_REC_Cherenkov_nphe;   //!
   TBranch        *b_REC_Cherenkov_time;   //!
   TBranch        *b_REC_Cherenkov_path;   //!
   TBranch        *b_REC_Cherenkov_theta;   //!
   TBranch        *b_REC_Cherenkov_phi;   //!
   TBranch        *b_REC_ForwardTagger_pindex;   //!
   TBranch        *b_REC_ForwardTagger_detector;   //!
   TBranch        *b_REC_ForwardTagger_energy;   //!
   TBranch        *b_REC_ForwardTagger_time;   //!
   TBranch        *b_REC_ForwardTagger_path;   //!
   TBranch        *b_REC_ForwardTagger_x;   //!
   TBranch        *b_REC_ForwardTagger_y;   //!
   TBranch        *b_REC_ForwardTagger_z;   //!
   TBranch        *b_REC_ForwardTagger_dx;   //!
   TBranch        *b_REC_ForwardTagger_dy;   //!
   TBranch        *b_REC_ForwardTagger_radius;   //!
   TBranch        *b_REC_ForwardTagger_size;   //!
   TBranch        *b_REC_Scintillator_pindex;   //!
   TBranch        *b_REC_Scintillator_detector;   //!
   TBranch        *b_REC_Scintillator_sector;   //!
   TBranch        *b_REC_Scintillator_layer;   //!
   TBranch        *b_REC_Scintillator_component;   //!
   TBranch        *b_REC_Scintillator_energy;   //!
   TBranch        *b_REC_Scintillator_time;   //!
   TBranch        *b_REC_Scintillator_path;   //!
   TBranch        *b_REC_Track_pindex;   //!
   TBranch        *b_REC_Track_detector;   //!
   TBranch        *b_REC_Track_sector;   //!
   TBranch        *b_REC_Track_chi2;   //!
   TBranch        *b_REC_Track_NDF;   //!
   TBranch        *b_REC_Track_chi2_nomm;   //!
   TBranch        *b_REC_Track_NDF_nomm;   //!
   TBranch        *b_REC_Traj_pindex;   //!
   TBranch        *b_REC_Traj_detId;   //!
   TBranch        *b_REC_Traj_x;   //!
   TBranch        *b_REC_Traj_y;   //!
   TBranch        *b_REC_Traj_z;   //!
   TBranch        *b_REC_Traj_cx;   //!
   TBranch        *b_REC_Traj_cy;   //!
   TBranch        *b_REC_Traj_cz;   //!

   Event12(TTree *tree=0);
   virtual ~Event12();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Int_t    GetEntries(){ return fChain->GetEntries(); } 
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Event12_cxx
Event12::Event12(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
     //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("out_clas12_2GeV_1.root");
     // if (!f || !f->IsOpen()) {
     //    f = new TFile("out_clas12_2GeV_1.root");
     // }
     // f->GetObject("clas12",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("clas12","");
      //chain->Add("out_clas12_2GeV_1.root/clas12");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

Event12::~Event12()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Event12::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Event12::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Event12::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   RUN_config_torus = 0;
   RUN_config_run = 0;
   RUN_config_event = 0;
   RUN_config_solenoid = 0;
   RAW_scaler_crate = 0;
   RAW_scaler_slot = 0;
   RAW_scaler_channel = 0;
   RAW_scaler_helicity = 0;
   RAW_scaler_quartet = 0;
   RAW_scaler_value = 0;
   REC_Event_NRUN = 0;
   REC_Event_NEVENT = 0;
   REC_Event_EVNTime = 0;
   REC_Event_TYPE = 0;
   REC_Event_TRG = 0;
   REC_Event_BCG = 0;
   REC_Event_STTime = 0;
   REC_Event_RFTime = 0;
   REC_Event_Helic = 0;
   REC_Particle_pid = 0;
   REC_Particle_px = 0;
   REC_Particle_py = 0;
   REC_Particle_pz = 0;
   REC_Particle_vx = 0;
   REC_Particle_vy = 0;
   REC_Particle_vz = 0;
   REC_Particle_charge = 0;
   REC_Particle_beta = 0;
   REC_Particle_chi2pid = 0;
   REC_Particle_status = 0;
   REC_Calorimeter_pindex = 0;
   REC_Calorimeter_detector = 0;
   REC_Calorimeter_sector = 0;
   REC_Calorimeter_layer = 0;
   REC_Calorimeter_energy = 0;
   REC_Calorimeter_time = 0;
   REC_Calorimeter_path = 0;
   REC_Calorimeter_x = 0;
   REC_Calorimeter_y = 0;
   REC_Calorimeter_z = 0;
   REC_Calorimeter_lu = 0;
   REC_Calorimeter_lv = 0;
   REC_Calorimeter_lw = 0;
   REC_Cherenkov_pindex = 0;
   REC_Cherenkov_detector = 0;
   REC_Cherenkov_sector = 0;
   REC_Cherenkov_nphe = 0;
   REC_Cherenkov_time = 0;
   REC_Cherenkov_path = 0;
   REC_Cherenkov_theta = 0;
   REC_Cherenkov_phi = 0;
   REC_ForwardTagger_pindex = 0;
   REC_ForwardTagger_detector = 0;
   REC_ForwardTagger_energy = 0;
   REC_ForwardTagger_time = 0;
   REC_ForwardTagger_path = 0;
   REC_ForwardTagger_x = 0;
   REC_ForwardTagger_y = 0;
   REC_ForwardTagger_z = 0;
   REC_ForwardTagger_dx = 0;
   REC_ForwardTagger_dy = 0;
   REC_ForwardTagger_radius = 0;
   REC_ForwardTagger_size = 0;
   REC_Scintillator_pindex = 0;
   REC_Scintillator_detector = 0;
   REC_Scintillator_sector = 0;
   REC_Scintillator_layer = 0;
   REC_Scintillator_component = 0;
   REC_Scintillator_energy = 0;
   REC_Scintillator_time = 0;
   REC_Scintillator_path = 0;
   REC_Track_pindex = 0;
   REC_Track_detector = 0;
   REC_Track_sector = 0;
   REC_Track_chi2 = 0;
   REC_Track_NDF = 0;
   REC_Track_chi2_nomm = 0;
   REC_Track_NDF_nomm = 0;
   REC_Traj_pindex = 0;
   REC_Traj_detId = 0;
   REC_Traj_x = 0;
   REC_Traj_y = 0;
   REC_Traj_z = 0;
   REC_Traj_cx = 0;
   REC_Traj_cy = 0;
   REC_Traj_cz = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RUN_config_torus", &RUN_config_torus, &b_RUN_config_torus);
   fChain->SetBranchAddress("RUN_config_run", &RUN_config_run, &b_RUN_config_run);
   fChain->SetBranchAddress("RUN_config_event", &RUN_config_event, &b_RUN_config_event);
   fChain->SetBranchAddress("RUN_config_solenoid", &RUN_config_solenoid, &b_RUN_config_solenoid);
   fChain->SetBranchAddress("RAW_scaler_crate", &RAW_scaler_crate, &b_RAW_scaler_crate);
   fChain->SetBranchAddress("RAW_scaler_slot", &RAW_scaler_slot, &b_RAW_scaler_slot);
   fChain->SetBranchAddress("RAW_scaler_channel", &RAW_scaler_channel, &b_RAW_scaler_channel);
   fChain->SetBranchAddress("RAW_scaler_helicity", &RAW_scaler_helicity, &b_RAW_scaler_helicity);
   fChain->SetBranchAddress("RAW_scaler_quartet", &RAW_scaler_quartet, &b_RAW_scaler_quartet);
   fChain->SetBranchAddress("RAW_scaler_value", &RAW_scaler_value, &b_RAW_scaler_value);
   fChain->SetBranchAddress("REC_Event_NRUN", &REC_Event_NRUN, &b_REC_Event_NRUN);
   fChain->SetBranchAddress("REC_Event_NEVENT", &REC_Event_NEVENT, &b_REC_Event_NEVENT);
   fChain->SetBranchAddress("REC_Event_EVNTime", &REC_Event_EVNTime, &b_REC_Event_EVNTime);
   fChain->SetBranchAddress("REC_Event_TYPE", &REC_Event_TYPE, &b_REC_Event_TYPE);
   fChain->SetBranchAddress("REC_Event_TRG", &REC_Event_TRG, &b_REC_Event_TRG);
   fChain->SetBranchAddress("REC_Event_BCG", &REC_Event_BCG, &b_REC_Event_BCG);
   fChain->SetBranchAddress("REC_Event_STTime", &REC_Event_STTime, &b_REC_Event_STTime);
   fChain->SetBranchAddress("REC_Event_RFTime", &REC_Event_RFTime, &b_REC_Event_RFTime);
   fChain->SetBranchAddress("REC_Event_Helic", &REC_Event_Helic, &b_REC_Event_Helic);
   fChain->SetBranchAddress("REC_Particle_pid", &REC_Particle_pid, &b_REC_Particle_pid);
   fChain->SetBranchAddress("REC_Particle_px", &REC_Particle_px, &b_REC_Particle_px);
   fChain->SetBranchAddress("REC_Particle_py", &REC_Particle_py, &b_REC_Particle_py);
   fChain->SetBranchAddress("REC_Particle_pz", &REC_Particle_pz, &b_REC_Particle_pz);
   fChain->SetBranchAddress("REC_Particle_vx", &REC_Particle_vx, &b_REC_Particle_vx);
   fChain->SetBranchAddress("REC_Particle_vy", &REC_Particle_vy, &b_REC_Particle_vy);
   fChain->SetBranchAddress("REC_Particle_vz", &REC_Particle_vz, &b_REC_Particle_vz);
   fChain->SetBranchAddress("REC_Particle_charge", &REC_Particle_charge, &b_REC_Particle_charge);
   fChain->SetBranchAddress("REC_Particle_beta", &REC_Particle_beta, &b_REC_Particle_beta);
   fChain->SetBranchAddress("REC_Particle_chi2pid", &REC_Particle_chi2pid, &b_REC_Particle_chi2pid);
   fChain->SetBranchAddress("REC_Particle_status", &REC_Particle_status, &b_REC_Particle_status);
   fChain->SetBranchAddress("REC_Calorimeter_pindex", &REC_Calorimeter_pindex, &b_REC_Calorimeter_pindex);
   fChain->SetBranchAddress("REC_Calorimeter_detector", &REC_Calorimeter_detector, &b_REC_Calorimeter_detector);
   fChain->SetBranchAddress("REC_Calorimeter_sector", &REC_Calorimeter_sector, &b_REC_Calorimeter_sector);
   fChain->SetBranchAddress("REC_Calorimeter_layer", &REC_Calorimeter_layer, &b_REC_Calorimeter_layer);
   fChain->SetBranchAddress("REC_Calorimeter_energy", &REC_Calorimeter_energy, &b_REC_Calorimeter_energy);
   fChain->SetBranchAddress("REC_Calorimeter_time", &REC_Calorimeter_time, &b_REC_Calorimeter_time);
   fChain->SetBranchAddress("REC_Calorimeter_path", &REC_Calorimeter_path, &b_REC_Calorimeter_path);
   fChain->SetBranchAddress("REC_Calorimeter_x", &REC_Calorimeter_x, &b_REC_Calorimeter_x);
   fChain->SetBranchAddress("REC_Calorimeter_y", &REC_Calorimeter_y, &b_REC_Calorimeter_y);
   fChain->SetBranchAddress("REC_Calorimeter_z", &REC_Calorimeter_z, &b_REC_Calorimeter_z);
   fChain->SetBranchAddress("REC_Calorimeter_lu", &REC_Calorimeter_lu, &b_REC_Calorimeter_lu);
   fChain->SetBranchAddress("REC_Calorimeter_lv", &REC_Calorimeter_lv, &b_REC_Calorimeter_lv);
   fChain->SetBranchAddress("REC_Calorimeter_lw", &REC_Calorimeter_lw, &b_REC_Calorimeter_lw);
   fChain->SetBranchAddress("REC_Cherenkov_pindex", &REC_Cherenkov_pindex, &b_REC_Cherenkov_pindex);
   fChain->SetBranchAddress("REC_Cherenkov_detector", &REC_Cherenkov_detector, &b_REC_Cherenkov_detector);
   fChain->SetBranchAddress("REC_Cherenkov_sector", &REC_Cherenkov_sector, &b_REC_Cherenkov_sector);
   fChain->SetBranchAddress("REC_Cherenkov_nphe", &REC_Cherenkov_nphe, &b_REC_Cherenkov_nphe);
   fChain->SetBranchAddress("REC_Cherenkov_time", &REC_Cherenkov_time, &b_REC_Cherenkov_time);
   fChain->SetBranchAddress("REC_Cherenkov_path", &REC_Cherenkov_path, &b_REC_Cherenkov_path);
   fChain->SetBranchAddress("REC_Cherenkov_theta", &REC_Cherenkov_theta, &b_REC_Cherenkov_theta);
   fChain->SetBranchAddress("REC_Cherenkov_phi", &REC_Cherenkov_phi, &b_REC_Cherenkov_phi);
   fChain->SetBranchAddress("REC_ForwardTagger_pindex", &REC_ForwardTagger_pindex, &b_REC_ForwardTagger_pindex);
   fChain->SetBranchAddress("REC_ForwardTagger_detector", &REC_ForwardTagger_detector, &b_REC_ForwardTagger_detector);
   fChain->SetBranchAddress("REC_ForwardTagger_energy", &REC_ForwardTagger_energy, &b_REC_ForwardTagger_energy);
   fChain->SetBranchAddress("REC_ForwardTagger_time", &REC_ForwardTagger_time, &b_REC_ForwardTagger_time);
   fChain->SetBranchAddress("REC_ForwardTagger_path", &REC_ForwardTagger_path, &b_REC_ForwardTagger_path);
   fChain->SetBranchAddress("REC_ForwardTagger_x", &REC_ForwardTagger_x, &b_REC_ForwardTagger_x);
   fChain->SetBranchAddress("REC_ForwardTagger_y", &REC_ForwardTagger_y, &b_REC_ForwardTagger_y);
   fChain->SetBranchAddress("REC_ForwardTagger_z", &REC_ForwardTagger_z, &b_REC_ForwardTagger_z);
   fChain->SetBranchAddress("REC_ForwardTagger_dx", &REC_ForwardTagger_dx, &b_REC_ForwardTagger_dx);
   fChain->SetBranchAddress("REC_ForwardTagger_dy", &REC_ForwardTagger_dy, &b_REC_ForwardTagger_dy);
   fChain->SetBranchAddress("REC_ForwardTagger_radius", &REC_ForwardTagger_radius, &b_REC_ForwardTagger_radius);
   fChain->SetBranchAddress("REC_ForwardTagger_size", &REC_ForwardTagger_size, &b_REC_ForwardTagger_size);
   fChain->SetBranchAddress("REC_Scintillator_pindex", &REC_Scintillator_pindex, &b_REC_Scintillator_pindex);
   fChain->SetBranchAddress("REC_Scintillator_detector", &REC_Scintillator_detector, &b_REC_Scintillator_detector);
   fChain->SetBranchAddress("REC_Scintillator_sector", &REC_Scintillator_sector, &b_REC_Scintillator_sector);
   fChain->SetBranchAddress("REC_Scintillator_layer", &REC_Scintillator_layer, &b_REC_Scintillator_layer);
   fChain->SetBranchAddress("REC_Scintillator_component", &REC_Scintillator_component, &b_REC_Scintillator_component);
   fChain->SetBranchAddress("REC_Scintillator_energy", &REC_Scintillator_energy, &b_REC_Scintillator_energy);
   fChain->SetBranchAddress("REC_Scintillator_time", &REC_Scintillator_time, &b_REC_Scintillator_time);
   fChain->SetBranchAddress("REC_Scintillator_path", &REC_Scintillator_path, &b_REC_Scintillator_path);
   fChain->SetBranchAddress("REC_Track_pindex", &REC_Track_pindex, &b_REC_Track_pindex);
   fChain->SetBranchAddress("REC_Track_detector", &REC_Track_detector, &b_REC_Track_detector);
   fChain->SetBranchAddress("REC_Track_sector", &REC_Track_sector, &b_REC_Track_sector);
   fChain->SetBranchAddress("REC_Track_chi2", &REC_Track_chi2, &b_REC_Track_chi2);
   fChain->SetBranchAddress("REC_Track_NDF", &REC_Track_NDF, &b_REC_Track_NDF);
   fChain->SetBranchAddress("REC_Track_chi2_nomm", &REC_Track_chi2_nomm, &b_REC_Track_chi2_nomm);
   fChain->SetBranchAddress("REC_Track_NDF_nomm", &REC_Track_NDF_nomm, &b_REC_Track_NDF_nomm);
   fChain->SetBranchAddress("REC_Traj_pindex", &REC_Traj_pindex, &b_REC_Traj_pindex);
   fChain->SetBranchAddress("REC_Traj_detId", &REC_Traj_detId, &b_REC_Traj_detId);
   fChain->SetBranchAddress("REC_Traj_x", &REC_Traj_x, &b_REC_Traj_x);
   fChain->SetBranchAddress("REC_Traj_y", &REC_Traj_y, &b_REC_Traj_y);
   fChain->SetBranchAddress("REC_Traj_z", &REC_Traj_z, &b_REC_Traj_z);
   fChain->SetBranchAddress("REC_Traj_cx", &REC_Traj_cx, &b_REC_Traj_cx);
   fChain->SetBranchAddress("REC_Traj_cy", &REC_Traj_cy, &b_REC_Traj_cy);
   fChain->SetBranchAddress("REC_Traj_cz", &REC_Traj_cz, &b_REC_Traj_cz);
   Notify();
}

Bool_t Event12::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Event12::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Event12::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Event12_cxx
