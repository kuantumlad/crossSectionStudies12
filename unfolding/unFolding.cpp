#include <iostream>
#include <vector>

#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixTBase.h"
#include "TDecompSVD.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "tPhi.C"

int bins = 200;

void unFolding( const char* tempdir, int max_files ){

  if( max_files < 1 ){ std::cout << ">> ERROR VALUE NEEDS TO BE LARGER THAN 0" << std::endl; exit(1); }
  TChain *fchain = new TChain("tPhi");
  TChain *realdata = new TChain("tPhi");
  for( int i = 1; i < max_files; i++ ){    
    //std::cout << ">> ADDING " << Form("%s/pid_phi_%d.root",tempdir,i) << std::endl;
    fchain->Add(Form("%s/pid_phi_%d.root",tempdir,i));
  }
  for( int i = max_files; i <= 200; i++ ){
    realdata->Add(Form("%s/pid_phi_%d.root", tempdir, i));  
  }

  //SPECIFIY THE BINNING YOU WANT HERE
  TH1D *h_rec = new TH1D("h_rec","h_rec",bins, 0.0, 10.5);
  TH1D *h_gen = new TH1D("h_gen","h_gen",bins, 0.0, 10.5);

  tPhi *treevar = new tPhi(fchain);  
  Long64_t num_entries = fchain->GetEntries();  

  TMatrixD prob_mig(bins,bins);// TMatrixD(20,20);
  std::map<int,int> m_normalize;
  std::vector<double> v_normalize(2*bins);

  std::cout << ">> POPULATING MATRIX NOW..." << std::endl; 
  for( Long64_t nn = 0; nn < num_entries; nn++ ){
    fchain->GetEntry(nn);

    double rec_q2 = -treevar->Q2;
    double gen_q2 = -treevar->mc_Q2;

    //std::cout << " >> " << rec_q2 << " " << gen_q2 << std::endl;
    
    int gen_bin = h_gen->FindBin(gen_q2);
    if( rec_q2 > 0.0  && gen_bin > 0 ){
      int rec_bin = h_rec->FindBin(rec_q2);
      if ( rec_bin > 0 ){
      //std::cout << " >> MATRIX: i " << rec_bin << " " << gen_bin << std::endl;      
      prob_mig(rec_bin,gen_bin)++;
      }
    }
    if( gen_bin > 0.0 ){
      //std::cout << gen_bin << std::endl;
      v_normalize[gen_bin]++; 
    }
  }
  //  std::cout << " >> SIZE OF v_normalize " << v_normalize.size() << std::endl;
  //for( std::vector<double>::iterator it = v_normalize.begin(); it != v_normalize.end(); ++it ){
  //  std::cout << " >> " << *it << std::endl;
  // }

  std::cout << " >> MATRIX BEFORE " << std::endl;
  //prob_mig.Print();
  for( int col = 0; col < bins; col++ ){
    for( int row = 0; row < bins; row++ ){
      if( v_normalize[col] == 0 )
	{
	  prob_mig(row,col) = 0;
	}
      else{
	std::cout << " >> normalize " << prob_mig(row,col) << " / " << v_normalize[col] << " " << prob_mig(row,col)/v_normalize[col] << std::endl;

	prob_mig(row,col) = prob_mig(row,col)/v_normalize[col];
      }
    }
  }
  std::cout << " >> FINAL MATRIX " << std::endl;
  prob_mig.Print();

  ////////////////////////////////////////////////////////////
  //UNFOLD THE OTHER HALF OF THE DATA
  tPhi *datatreevar = new tPhi(realdata);  
  Long64_t num_entries2 = realdata->GetEntries();  

  TH1D *h_gen_q2 = new TH1D("h_gen_q2","h_gen_q2", bins, 0.0, 10.5);
  TH1D *h_rec_q2 = new TH1D("h_rec_q2","h_rec_q2", bins, 0.0, 10.5);
  TH1D *h_rec_tempq2 = new TH1D("h_rec_tempq2","h_rec_tempq2", bins, 0.0, 10.5);

  TH1D *h_unfold_q2 = new TH1D("h_unfold_q2","h_unfold_q2", bins, 0.0, 10.5 );
  std::vector<double> v_unfolded;

  for( int ii = 0; ii < num_entries2; ii++ ){
    realdata->GetEntry(ii);
    double mc_q2 = -datatreevar->mc_Q2;
    double q2 = -datatreevar->Q2;
    //std::cout << q2 << std::endl;
    h_gen_q2->Fill(mc_q2);
    if( q2 > 0.0 ){
      
      h_rec_q2->Fill(q2);
      h_rec_tempq2->Fill(q2);
    }
        
  }
  std::cout << " >> INVERTING MATRIX " << std::endl;
  //////////////////////////////////
  //INVERT MATRIX TO SOLVE y = A^-1 * x
  std::cout << " DETERMINANT BEFORE " << prob_mig.Determinant() << std::endl;
  TDecompSVD svd(prob_mig);
  svd.Decompose();
  TMatrixD inv_prob_mig = svd.Invert();
  
  inv_prob_mig.Print();

  double true_q2 = 0.0;
  
  TH1D *h_unfold_m1 = new TH1D("h_unfold_m1","h_unfold_m1",bins,0.0, 10.5);
  
  for( int row = 0; row < bins; row++ ){
    double q2_bin_value = h_rec_q2->GetBinContent(row);
    //std::cout << " >> BIN CONTENT FOR DATA HISTOGRAMS IS " << q2_bin_value << std::endl;

    for( int col = 0; col < bins; col++ ){      
      double temp = q2_bin_value * inv_prob_mig(row,col);
      //std::cout << q2_bin_value << " " << inv_prob_mig(row,col) << std::endl;
      
      true_q2 = true_q2 + temp;
      //std::cout <<  ">> > > true_q2 " << true_q2 << std::endl;
    
    }
    //std::cout << " >> true q2 in bin " << row << " is " << true_q2 << std::endl;
    h_unfold_q2->SetBinContent(row,true_q2);

    true_q2 = 0.0;
  }
  
  ///ACCEPTANCE
  TH1D *h_rec_temp = h_rec_tempq2;
  h_rec_temp->Divide(h_gen_q2);

  for ( int i = 0; i < bins; i++ ){
    double q2_rec = h_rec_q2->GetBinContent(i);
    double temp = q2_rec * inv_prob_mig(i, i);
    //std::cout << q2_rec << " " << inv_prob_mig(i,i) << std::endl;
    h_unfold_m1->SetBinContent(i, temp );
    }
  
    
  
  TCanvas *c1 = new TCanvas("c1","c1",1600,800);
  c1->Divide(2,1);
  c1->cd(1);
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  h_gen_q2->SetTitle("Unfolded Q^{2} Distribution vs Gen Q^{2}");
  h_gen_q2->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
  h_unfold_q2->SetLineColor(kRed);
  h_rec_q2->SetLineWidth(kBlack);
  h_rec_q2->SetLineStyle(3);
  h_gen_q2->SetLineColor(kBlue+2);
  h_gen_q2->SetLineWidth(1);
  h_gen_q2->Draw();
  h_rec_q2->Draw("h+same");
  h_unfold_q2->Draw("h+same");
  TLegend *l1 = new TLegend(0.7, 0.7, 0.9, 0.9);
  l1->AddEntry(h_gen_q2,"Generated Q^{2} distribution","l");
  l1->AddEntry(h_rec_q2,"Reconstructed Q^{2} distribution","l");
  l1->AddEntry(h_unfold_q2,"Unfolded Q^{2} distribution","l");
  l1->Draw();

  c1->cd(2);
  gPad->SetLogy();
  h_unfold_m1->SetTitle("Standard Acceptance Corrections");
  h_unfold_m1->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
  h_unfold_m1->SetLineColor(kRed);
  h_gen_q2->SetLineColor(kBlue+2);  
  h_unfold_m1->Draw();
  h_gen_q2->Draw("h+same");
  TLegend *l2 = new TLegend(0.7, 0.7, 0.9, 0.9);
  l2->AddEntry(h_gen_q2,"Generated Q^{2} distribution","l");
  l2->AddEntry(h_unfold_m1,"Standard Acceptance corrected Q^{2} distribution","l");
  l2->Draw();
  

  TCanvas *c2 = new TCanvas("c2","c2",800,800); 
  c2->cd();
  h_rec_temp->SetTitle("Acceptance for Q^{2}: REC/GEN");
  h_rec_temp->GetXaxis()->SetTitle("Q^{2} [GeV^2]");
  h_rec_temp->Draw("p0");
  
   c1->SaveAs("/u/home/bclary/CLAS12/phi_analysis/v1/utilities/hlPhysics/h_unfolded_q2.pdf");
  c2->SaveAs("/u/home/bclary/CLAS12/phi_analysis/v1/utilities/hlPhysics/h_recgen_q2.pdf");

  std::cout << " >> COMPLETE " << std::endl;

}
