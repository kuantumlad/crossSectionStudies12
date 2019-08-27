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
#include "TStyle.h"
#include "TROOT.h"
#include "TLatex.h"
#include "fstream"

using namespace std;

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

using namespace ROOT::Math;

float eBeam = 7.546;//2.221;


double kin_W(TLorentzVector ele, float Ebeam);
TF1 *fitHisto(TH1D *htemp);

// model is the run numbers
// model == 11 if it is for simulation
int quickFit(const char* input_file, int model){

  std::cout << " Executing quick fit macro to get fit parameters " << std::endl;

  int nWBins = 250;
  int nQ2Bins = 250;
  float lowBorderW = 0.80;
  float highBorderW = 1.7;  
  float lowBorderQ2 = 0.0;
  float highBorderQ2 = 1.0;
  float W, Q2;
  int sectorElectron;
  int sectorProton;

  int nOverlallBins = 250;
  float eBeamLow = 0;
  float eBeamHigh = 2.5;
  float thetaELow = 5;
  float thetaEHigh = 30;
  float thetaPLow = 50;
  float thetaPHigh = 90;

  float pELow = 1.6;
  float pEHigh = 2.2;
  float pPLow = 0.1;
  float pPHigh = 1.0;

  float deltaEenergyLow = -1.2;
  float deltaEenergyHigh = 1.2;

  float deltaThetaELow = -10;
  float deltaThetaEHigh = 10;

  float deltaThetaPLow = -10;
  float deltaThetaPHigh = 14;
	
  float deltaPELow = -0.2;
  float deltaPEHigh = 0.2;

  float deltaPPLow = -0.4;
  float deltaPPHigh = 0.4;

  std::string parentInDirectory = "/w/hallb-scifs17exp/clas12/bclary/CLAS12/electron_studies/";
  TFile f(Form("%s",input_file));
  
  std::string parentOutDirectory = " /w/hallb-scifs17exp/clas12/bclary/CLAS12/electron_studies/physics/parameters/";
  TFile *quick_fit_results = new TFile(Form("%squick_fit%d.root", parentOutDirectory.c_str(), model), "recreate");

  float mProton = 0.98327;
  TTree *anaTree=(TTree *) f.Get("out_tree");
  TTree *simTree=(TTree *) f.Get("mc_tree");


  float part_px[100], part_py[100], part_pz[100], part_E[100];
  float part_vx[100], part_vy[100], part_vz[100];

  float prot_px[100], prot_py[100], prot_pz[100], prot_E[100];
  float prot_vx[100], prot_vy[100], prot_vz[100];

  float ele_vx[100], ele_vy[100], ele_vz[100];
  int sectorE[100];

  vector<float>  *vpart_px      = 0;
  vector<float>  *vpart_py      = 0;
  vector<float>  *vpart_pz      = 0;

  vector<float>  *vpart_vx      = 0;
  vector<float>  *vpart_vy      = 0;
  vector<float>  *vpart_vz      = 0;
  vector<float>  *vpart_E       = 0;
  vector<int>    *v_sectorE     = 0;

  vector<float>  *vprot_px      = 0;
  vector<float>  *vprot_py      = 0;
  vector<float>  *vprot_pz      = 0;
  vector<float>  *vprot_vx      = 0;
  vector<float>  *vprot_vy      = 0;
  vector<float>  *vprot_vz      = 0;
  vector<float>  *vprot_E       = 0;



  TLorentzVector p4_ele[100], p4_ele_gen[100], initialE, initialP, p4_ele_final;
  TLorentzVector p4_proton[100];
  anaTree->SetBranchAddress("p4_ele_px", &vpart_px);
  anaTree->SetBranchAddress("p4_ele_py", &vpart_py);
  anaTree->SetBranchAddress("p4_ele_pz", &vpart_pz);
  anaTree->SetBranchAddress("p4_ele_E", &vpart_E);

  //anaTree->SetBranchAddress("p4_ele_vx", &vpart_vx);
  //anaTree->SetBranchAddress("p4_ele_vy", &vpart_vy);
  //anaTree->SetBranchAddress("p4_ele_vz", &vpart_vz);

  anaTree->SetBranchAddress("p4_prot_px", &vprot_px);
  anaTree->SetBranchAddress("p4_prot_py", &vprot_py);
  anaTree->SetBranchAddress("p4_prot_pz", &vprot_pz);
  anaTree->SetBranchAddress("p4_prot_E",  &vprot_E);

  //anaTree->SetBranchAddress("p4_prot_vx", &vprot_vx);
  //anaTree->SetBranchAddress("p4_prot_vy", &vprot_vy);
  //anaTree->SetBranchAddress("p4_prot_vz", &vprot_vz);

  anaTree->SetBranchAddress("sectorE", &v_sectorE);
  initialE.SetPxPyPzE(0, 0, eBeam, eBeam);
  initialP.SetPxPyPzE(0, 0, 0, mProton);


  // change histogram range before proceeding 
  if( eBeam > 3 ){
    highBorderW = 1.1 ;
  }


  TH1D *wSectorExclusive[6];
  TH1D *phiEphiP[6];

  for( int s = 0; s < 6; s++ ){
    
    wSectorExclusive[s] = new TH1D(Form("wExclusiveS%d", s + 1), Form("wExclusiveS%d", s + 1), nWBins, lowBorderW, highBorderW);
    wSectorExclusive[s]->GetXaxis()->SetTitle("W [GeV]");
    wSectorExclusive[s]->SetTitle(Form("Sector %d, W, e(FD)", s + 1));

    phiEphiP[s] = new TH1D(Form("phiEphiPS%d", s + 1), Form("phiEphiPS%d", s + 1), 360, 150, 210);
    phiEphiP[s]->GetXaxis()->SetTitle("#phi_e - #phi_p [deg]");

  }  

  int NPart = 0;
  int NPartP = 0;      

  double toRD=180.0/TMath::Pi();
  
  int nentries = anaTree->GetEntriesFast();
  if( nentries > 2000000 ){
    nentries = 2000000;
  }


  for(Int_t k=0; k < nentries; k++){
    anaTree->GetEntry(k);
    NPart = vpart_px->size();
    if (NPart > 0 && NPart < 2){
      for(Int_t i = 0; i < NPart; i++){
	part_px[i] = vpart_px->at(i);
	part_py[i] = vpart_py->at(i);
	part_pz[i] = vpart_pz->at(i);
	//part_vx[i] = vpart_vx->at(i);
	//part_vy[i] = vpart_vy->at(i);
	//part_vz[i] = vpart_vz->at(i);
	part_E[i] = vpart_E->at(i);
	sectorE[i] = v_sectorE->at(i);
	sectorElectron = sectorE[i] - 1;
	if (sectorElectron > -1 && sectorElectron < 7){
	  p4_ele[i].SetPxPyPzE(part_px[i], part_py[i], part_pz[i], part_E[i]);
	  W = kin_W(p4_ele[i],  eBeam);
	  wSectorExclusive[sectorElectron]->Fill(W);
	  
	  NPartP = vprot_px->size();
	  if (NPartP > 0 && NPartP < 2){
	    for(Int_t j = 0; j < NPartP; j++){
	      prot_px[j] = vprot_px->at(j);
	      prot_py[j] = vprot_py->at(j);
	      prot_pz[j] = vprot_pz->at(j);
	      prot_E[j]  = vprot_E->at(j);
	      //prot_vx[i] = vprot_vx->at(i);
	      //prot_vy[i] = vprot_vy->at(i);
	      //prot_vz[i] = vprot_vz->at(i);
	      p4_proton[j].SetPxPyPzE(prot_px[j], prot_py[j], prot_pz[j], prot_E[j]);
	      float deltaPhi = p4_ele[i].Phi()*toRD - p4_proton[j].Phi()*toRD;
	      if (deltaPhi < 0) deltaPhi = -deltaPhi;
	      //std::cout <<  " >>  " << deltaPhi << " " << W << std::endl;
	      phiEphiP[sectorElectron]->Fill(deltaPhi);

	    }
	  }	       
	}
      }
    }
  }
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Fit W and delta phi spectrum.
  // The fit parameters will be used in the electronProtonFinal
  // macro.    

  ofstream outputWCuts;
  std::string parentDirectory = "/w/hallb-scifs17exp/clas12/bclary/CLAS12/electron_studies/physics/parameters/";
  std::string f_out_w_name = parentDirectory+"w_cut_limits_run"+std::to_string(model)+".txt";
  std::cout << " Creating W cuts output file " << f_out_w_name << std::endl;
  outputWCuts.open(f_out_w_name);

  for( int s = 0; s < 6; s++ ){    
    TF1 *fittemp=fitHisto(wSectorExclusive[s]);
    std::cout << " Sector s: " << s << ", " << fittemp->GetParameter(1) << ", " << fittemp->GetParameter(2);    
    outputWCuts << s << " " << fittemp->GetParameter(1) << " " << fittemp->GetParameter(2) << std::endl;
  }

  ofstream outputPhiCuts;
  std::string f_out_phi_name = parentDirectory+ "phi_cut_limits_run"+std::to_string(model)+".txt";
  outputPhiCuts.open(f_out_phi_name);

  for( int s = 0; s < 6; s++ ){    
    TF1 *fittemp=fitHisto(phiEphiP[s]);
    std::cout << " Sector s: " << s << ", " << fittemp->GetParameter(1) << ", " << fittemp->GetParameter(2) << std::endl;    
    outputPhiCuts << s << " " << fittemp->GetParameter(1) << " " << fittemp->GetParameter(2) << std::endl;
  }

  outputWCuts.close();
  outputPhiCuts.close();

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Write to the output file
  
  for( int s = 0; s < 6 ; s++ ){
    phiEphiP[s]->Write();
    wSectorExclusive[s]->Write();
  }


  return 0;

}


TF1* fitHisto(TH1D* htemp){
  //////////////////////////////////////////////////
  // Start Andrew's fitting method:
  double xlow,xhigh,histmax;
  int binlow,binhigh,binmax;

  if( eBeam > 3 ){
    htemp->Rebin(6);
  }
  binmax = htemp->GetMaximumBin();
  histmax = htemp->GetMaximum();
  binlow=binmax;
  binhigh=binmax;

  // The 0.65 parameter can be changed, this basically means start at the peak and work your way left and right
  // until you've gotten to 65% of the max amplitude.
  while(htemp->GetBinContent(binhigh++) >= 0.65*histmax&&binhigh<=htemp->GetNbinsX()){};
  while(htemp->GetBinContent(binlow--) >= 0.65*histmax&&binlow>=1){};
    
  xlow = htemp->GetBinLowEdge(binlow);
  xhigh = htemp->GetBinLowEdge(binhigh+1);
    
  htemp->Fit("gaus","","",xlow,xhigh);

  TF1 *ftemp = (TF1*) htemp->GetListOfFunctions()->FindObject("gaus");
  ftemp->SetLineWidth(3);

  return ftemp;

  // End
  /////////////////////////////////////////////

}


double kin_W(TLorentzVector ele, float Ebeam){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  TLorentzVector fCM = fGamma + target;
  return fCM.M();
}
