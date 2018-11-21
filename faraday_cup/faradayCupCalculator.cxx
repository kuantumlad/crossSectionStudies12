#include <iostream>
#include <TMath.h>
#include <vector>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>


int faradayCupCalculator(Char_t *inFile, Char_t *outputfile, int run){

  const Char_t *inTree="clas12";


  vector<int>            *vNRUN     = 0;     // Run Number
  vector<int>            *vNEVENT   = 0;     // Event Number
  vector<float>          *vEVNTime  = 0;     // Enet time
  vector<short>          *vTYPE     = 0;     // Data or MC
  vector<long long int>  *vTRG      = 0;     // Trigger
  vector<float>          *vBCG      = 0;     // FCUP scaler
  vector<float>          *vSTTime   = 0;     // Event Start Time (ns)
  vector<float>          *vRFTime   = 0;     // RF Time (ns)
  vector<short>          *vHelic    = 0;     // helicity
  vector<short>          *vChannel  = 0;     // Channel
  vector<short>          *vSlot     = 0;     // Slot
  vector<int>            *vValue    = 0;     // Value

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

  anaTree->SetBranchAddress("REC_Event_NRUN", &vNRUN);
  anaTree->SetBranchAddress("REC_Event_NEVENT", &vNEVENT);
  anaTree->SetBranchAddress("REC_Event_EVNTime", &vEVNTime);
  anaTree->SetBranchAddress("REC_Event_TYPE", &vTYPE);
  anaTree->SetBranchAddress("REC_Event_TRG", &vTRG);
  anaTree->SetBranchAddress("REC_Event_BCG", &vBCG);
  anaTree->SetBranchAddress("REC_Event_STTime", &vSTTime);
  anaTree->SetBranchAddress("REC_Event_RFTime", &vRFTime);
  anaTree->SetBranchAddress("REC_Event_Helic", &vHelic);
  anaTree->SetBranchAddress("RAW_scaler_slot", &vSlot );
  anaTree->SetBranchAddress("RAW_scaler_channel", &vChannel );
  anaTree->SetBranchAddress("RAW_scaler_value",&vValue );  

  double events = anaTree->GetEntriesFast()/100;
  short TYPE, Helic;
  int NRUN, NEVENT;
  long long int TRG;
  float EVNTime, BCG, STTime, RFTime;
  short channel;
  short slot;
  int value;

  double tot_beam_charge = 0;

  for( int k = 0; k < anaTree->GetEntries(); k++ ){
    
    anaTree->GetEntry(k);
    
    if(k % 10000 == 0){      
      double percent = (k/100)/(events/100);
      printf("Analysing event number %i of %.00f (%.01f percent)\n", k, events*100, percent);
    }
    
    NRUN = 0; NEVENT = 0; TYPE = 0; TRG = 0; Helic = 0; EVNTime = 0; BCG = 0; STTime = 0; RFTime = 0;
    channel=0;
    slot=0;
    value=0;      
    
    //helicity = 0;
    //fcup = 0;
    
    
    if(vNRUN->size() > 0)    NRUN = vNRUN->at(0);
    if(vNEVENT->size() > 0)  NEVENT = vNEVENT->at(0);
    if(vEVNTime->size() > 0) EVNTime = vNEVENT->at(0);
    if(vTYPE->size() > 0)    TYPE = vTYPE->at(0);
    if(vTRG->size() > 0)     TRG = vNEVENT->at(0);
    if(vBCG->size() > 0)     BCG = vNEVENT->at(0);
    if(vSTTime->size() > 0)  STTime = vSTTime->at(0);
    if(vRFTime->size() > 0)  RFTime = vRFTime->at(0);
    if(vHelic->size() > 0)   Helic = vHelic->at(0);

    int scaler_rows = vChannel->size();
    //std::cout << " scaler rows " << scaler_rows << std::endl;
    double beam_charge = 0;
    if( scaler_rows > 0 ){
      for( int i = 0; i < scaler_rows; i++ ){
	int slot = vSlot->at(i);
	int channel = vChannel->at(i);
	if( slot == 0 && channel == 32 ){
	  int fc_scaler = vValue->at(i);
	  double true_freq = fc_scaler/(0.03333 - 0.0005);
	  double beam_current = (true_freq - 100.0)/906.2;
	  beam_charge = beam_current * (0.03333 - 0.0005 );
	  //std::cout << " fc scaler " << fc_scaler << std::endl;
	  //std::cout << " true freq " << true_freq << std::endl;
	  //std::cout << " beam current " << beam_current << std::endl;

	}
      }
      //std::cout << beam_charge << std::endl;
      tot_beam_charge = tot_beam_charge + beam_charge;  
    }
  }  


  std::cout << " TOTAL BEAM CHARGE: " <<  tot_beam_charge << std::endl;

  return 0;

}
