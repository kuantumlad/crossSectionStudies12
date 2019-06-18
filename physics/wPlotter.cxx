#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <vector>
#include <iostream>
#include <TCanvas.h>
#include <TStyle.h>

void wPlotter(const char* infile, int run){


  TFile *fin = new TFile(infile);
  
  float lowWValueCut[6];
  float highWValueCut[6];
  
  int sector = 0;
  double mean;
  double sigma;

  std::string parentDirectory = "/w/hallb-scifs17exp/clas12/bclary/CLAS12/electron_studies/physics/parameters/";
	
  string line;
  ifstream readFromWCut(parentDirectory+"w_cut_limits_run"+std::to_string(run)+".txt");
  if( readFromWCut.is_open() ){
    std::cout << " OPENED FILES " << std::endl;	
    while(readFromWCut >> sector ) {//std::getline (readFromWCut, line) ){
      readFromWCut >> mean >> sigma;
      std::cout << " >> W CUT PARAMETERS: " << sector << " " << mean << " " << sigma << std::endl;
      lowWValueCut[sector]= mean - 3*sigma;
      highWValueCut[sector]= mean + 3*sigma;
    }
  }
  else{
    std::cout << " ERROR OPENING FILE " << std::endl;
  }

  TH1D *wSector[6];

  TCanvas *c1 = new TCanvas("c1","c1",600,400);
  c1->Divide(3,2);

  for( int s = 0; s < 6; s++ ){
    
    c1->cd(s+1);
    TH1D *w_temp = (TH1D*)fin->Get(Form("wExclusiveS%d", s + 1));
    w_temp->SetTitle(Form("Sector %d",s+1));
    w_temp->GetXaxis()->CenterTitle();
    gStyle->SetOptStat(0);    
    w_temp->SetLineColor(kBlack);
    w_temp->Draw();

    c1->Update();
    TVirtualPad *p1 = c1->GetPad(1);
    double ymax = gPad->GetUymax();
    double ymin = gPad->GetUymin();
    std::cout << ymax << std::endl;
    TLine *l = new TLine(0.938,ymin, 0.938, ymax);
    l->SetLineColor(kRed);
    l->SetLineWidth(1);
    l->SetLineStyle(2);

    TBox *phi_reg = new TBox(lowWValueCut[s], ymin, highWValueCut[s], ymax );
    phi_reg->SetFillColorAlpha(kBlue-4,0.15);
    phi_reg->Draw("same");

    w_temp->Draw("same");
    //l->Draw("same");
    
  }

  c1->SaveAs(Form("wSpectrumRun%d.pdf",run));

}
