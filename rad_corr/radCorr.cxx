#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <string>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TCanvas.h>

TH1D *lundToHisto(const char* h_name, int n_bins, double min, double max, const char* path_to_lunds, int n_lunds ){
  

  TH1D *h_out = new TH1D(h_name,h_name, n_bins, min, max );



  int c_lunds=0;
  std::string in_path = path_to_lunds;

  while( c_lunds < n_lunds ){
    
    int n_part, charge, n0, pid, n1, n2;
    double px, py, pz, e, n3, vx, vy, vz;  
    
    std::cout << " reading from lund " << std::endl;
    std::string lund = path_to_lunds + std::to_string(c_lunds) + ".lund";
    std::cout << " lund file " << lund << std::endl;
    string line;
    ifstream readInLund(lund);//parentDirectory+"w_cut_limits_run"+std::to_string(run)+".txt");
    if( readInLund.is_open() ){
      std::cout << " OPENED FILE " << std::endl;	
      while(readInLund >> n_part ){ // >> charge >> n0 >> pid >> n1 >> n2 >> px >> py >> pz >> e >> n3 >> vx >> vy >> vz  ) {
	std::cout << n_part << std::endl; //" " <<  charge << " " << n0 << " " << pid << " " << n1 << " " << n2 << " " << px << " " << py << " " << pz << " " << e << " " << n3 << " " << vx  << "  " << vy <<" " << vz << std::endl;
	if ( pid == 11 ){
	  h_out->Fill(px);
	}
    
      

      }
    }
  
    c_lunds++;
  }


  return h_out;
}

int radCorr( const char* inFile, const char*outFile, int beam_energy ){

  TH1D *h_test = lundToHisto("tesT", 100, 0.0, 2.3, "/w/hallb-scifs17exp/clas12/bclary/CLAS12/electron_studies/elastic/lund_norad/clas12_2GeV_", 10 );
  

  TCanvas *c1 = new TCanvas("c1","c1",900,900);
  c1->cd(0);
  h_test->Draw();


  return 0;
}
