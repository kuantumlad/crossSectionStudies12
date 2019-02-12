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
  using namespace std;

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

  using namespace ROOT::Math;
float Ebeam = 10.594;
double kin_W(TLorentzVector ele);
double kin_Q2(TLorentzVector ele);
double charge(vector<int> *Crate, vector<int> *Slot, vector<int> *Channel, vector<float> *Value);

TH2F *wQ2 = new TH2F("wQ2", "wQ2", 250, 0, 5, 200, 0, 10);
TH2F *wQ2Gen = new TH2F("wQ2Gen", "wQ2Gen", 250, 0, 5, 200, 0, 10);
TH1F *currentTime = new TH1F("currentTime", "currentTime", 10000, 0, 10000);
TH1F *chargeTime = new TH1F("chargeTime", "chargeTime", 10000, 0, 10000);

TH1F *wInSectorAndQ2[20][6];
TH1F *wGenSectorAndQ2[20][6];

Int_t inclusive(){
  TFile *resultsF = new TFile("results.root", "recreate");
  TFile *f = new TFile("/lustre/expphy/volatile/clas/clase1/markov/12GeV/inclusive/inclusiveDiffPol/farm/filesClara/10.6GeV/InclusiveElastic/data/newCook/4078.root");
  TTree *anaTree=(TTree *) f->Get("out_tree");
  TTree *scalerTree = (TTree*) f->Get("scaler_tree"); 
  TTree *mcTree=(TTree *) f->Get("mc_tree");

 const static int BUFFER = 100;

 float part_px[BUFFER], part_py[BUFFER], part_pz[BUFFER], part_E[BUFFER];
 float part_gen_px[BUFFER], part_gen_py[BUFFER], part_gen_pz[BUFFER];

  int sectorE[BUFFER];
  float W = 0, Q2 = 0;
  for (int s = 0; s < 6; s++){
    for (int qN = 0; qN < 20; qN++){
      wInSectorAndQ2[qN][s] = new TH1F(Form("wInSectorAndQ2S%dQ%d", s + 1, qN + 1), Form("wInSectorAndQ2S%dQ%d", s + 1, qN + 1), 250, 0, 5);
      wGenSectorAndQ2[qN][s] = new TH1F(Form("wGenSectorAndQ2S%dQ%d", s + 1, qN + 1), Form("wGenSectorAndQ2S%dQ%d", s + 1, qN + 1), 250, 0, 5);
    }
  }

  vector<float>  *vsim_px       = 0;
  vector<float>  *vsim_py       = 0;
  vector<float>  *vsim_pz       = 0;
  vector<float>  *vpart_px      = 0;
  vector<float>  *vpart_py      = 0;
  vector<float>  *vpart_pz      = 0;
  vector<float>  *vpart_E       = 0;
  vector<int>    *v_sectorE     = 0;
  vector<float>  *Value         = 0;
  vector<int>    *Crate         = 0;
  vector<int>    *Slot          = 0;
  vector<int>    *Channel       = 0;

  TLorentzVector p4_ele[BUFFER], p4_ele_gen[BUFFER];

  mcTree->SetBranchAddress("p4_mc_px", &vsim_px);
  mcTree->SetBranchAddress("p4_mc_py", &vsim_py);
  mcTree->SetBranchAddress("p4_mc_pz", &vsim_pz);


  anaTree->SetBranchAddress("p4_ele_px", &vpart_px);
  anaTree->SetBranchAddress("p4_ele_py", &vpart_py);
  anaTree->SetBranchAddress("p4_ele_pz", &vpart_pz);
  anaTree->SetBranchAddress("p4_ele_E", &vpart_E);
  anaTree->SetBranchAddress("sectorE", &v_sectorE);

  scalerTree->SetBranchAddress("RAW_scaler_crate", &Crate);
  scalerTree->SetBranchAddress("RAW_scaler_slot", &Slot);
  scalerTree->SetBranchAddress("RAW_scaler_channel", &Channel);
  scalerTree->SetBranchAddress("RAW_scaler_value", &Value);


  int NPart = 0;
  int nEvents = 100000000;
  float chargeIntegrated = 0;
  float current = 0;
  int q2N;
   for(Int_t k=0; k<anaTree->GetEntriesFast();k++){
    anaTree->GetEntry(k);
    NPart = vpart_px->size();
    for(Int_t i = 0; i < NPart; i++){
      if( i < BUFFER){
	part_px[i] = vpart_px->at(i);
        part_py[i] = vpart_py->at(i);
        part_pz[i] = vpart_pz->at(i);
	part_E[i] = vpart_E->at(i);
	sectorE[i] = v_sectorE->at(i);
	p4_ele[i].SetPxPyPzE(part_px[i], part_py[i], part_pz[i], part_E[i]);
	W = kin_W(p4_ele[0]);
	Q2 = kin_Q2(p4_ele[0]);
	wQ2->Fill(W, Q2);
	q2N = int(Q2*2);
	if (q2N > -1 && q2N < 20 && sectorE[i]>0 && sectorE[i] < 7){
	  wInSectorAndQ2[q2N][sectorE[i] - 1]->Fill(W);
	}
      }
    }
    if(k == nEvents) break;
   }


   for(Int_t k=0; k<mcTree->GetEntriesFast();k++){
     mcTree->GetEntry(k);
     NPart = vsim_px->size();
     for(Int_t i = 0; i < NPart; i++){
       if( i < BUFFER){
	 part_gen_px[i] = vsim_px->at(i);
	 part_gen_py[i] = vsim_py->at(i);
	 part_gen_pz[i] = vsim_pz->at(i);
	 p4_ele_gen[i].SetXYZM(part_gen_px[i], part_gen_py[i], part_gen_pz[i], 0.000511);
	 W = kin_W(p4_ele_gen[0]);
	 Q2 = kin_Q2(p4_ele_gen[0]);
	 wQ2Gen->Fill(W, Q2);
	 q2N = int(Q2*2);
	 if (q2N > -1 && q2N < 20){
	   wGenSectorAndQ2[q2N][0]->Fill(W, 1/6.0);
	 }
       }
     }
   }


   int goodCurrentCounter = 0;
   cout <<"scaler entries : " << scalerTree->GetEntriesFast() << endl;
   for(Int_t u=0; u<scalerTree->GetEntriesFast();u++){
     scalerTree->GetEntry(u);
     chargeIntegrated = chargeIntegrated + charge(Crate, Slot, Channel, Value);
     current = charge(Crate, Slot, Channel, Value)/0.03283;
     if (current > 0){
       goodCurrentCounter++;
       currentTime->SetBinContent(goodCurrentCounter+1, current);
       chargeTime->SetBinContent(goodCurrentCounter+1, chargeIntegrated);
     }
   }
   resultsF->cd();
   resultsF->mkdir("overview");
   resultsF->cd("overview");
   wQ2->Write();
   wQ2Gen->Write();
   resultsF->mkdir("chargeCurrent");
   resultsF->cd("chargeCurrent");
   chargeTime->Write();
   currentTime->Write();
   resultsF->mkdir("Inclusive");
   resultsF->cd ("Inclusive");
   for (int s = 0; s < 6; s++){
     for (int qN = 0; qN < 20; qN++){
       wInSectorAndQ2[qN][s]->Write();
       wGenSectorAndQ2[qN][s]->Write();
     }
   }

   resultsF->Close();
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
double charge(vector<int> *Crate, vector<int> *Slot, vector<int> *Channel, vector<float> *Value){
  double cLocal = 0;
  for (int t = 0; t < Value->size(); t++){
    if (Slot->at(t) == 1 && Channel->at(t) == 0 && Crate->at(t)== 64){
      cLocal = 10*0.03283*(Value->at(t)/0.03283-100)/906.2;
    }
  }
  return cLocal;
}
