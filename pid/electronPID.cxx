#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <sstream>
#include "TMath.h"
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "Event12.C"

int electronPID(const char *input, const char *output, int run, const char *analysis_type){

  TFile *fIn = new TFile(input);
  TTree *inTree = (TTree*)fIn->Get("clas12");
  


  Event12 event;
  event.Init(inTree);
  int total_entries=event.GetEntries();
  std::cout << " Events to process " << total_entries << std::endl;  


  for( int i=0; i < total_entries;  i++ ){
    event.GetEntry(i);
    
    

    std::vector<float> ev_px = *event.REC_Particle_px;
    std::cout << " >> " << ev_px.size() << std::endl;

  }




  return 0;
}
