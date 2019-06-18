#include <TCanvas.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <string>
#include <vector>
#include <map>
#include <sstream>

int checkLund(){

  double el_mass = 0.000511;
  double pr_mass = 0.937282;
  double kp_mass = 0.493;
  double km_mass = 0.493;
  double phi_mass = 1.01946;
  TLorentzVector lv_ebeam(0, 0, 2.221, 2.221 );
  TLorentzVector lv_target(0,0,0,pr_mass);
  
  int el_fract = 8;
  std::string prefix = "test";

  TFile *fOut = new TFile(Form("outCheckLund_fract%d.root",el_fract),"RECREATE");
  TH2D *h_el_gen_w_theta = new TH2D("h_el_gen_w_theta","h_el_gen_w_theta",200, 0.0, 2.5, 200, 7.0, 60);  
  TH1D *h_el_w = new TH1D(Form("h_el_w_f%d",el_fract),Form("h_el_w_f%d",el_fract),200, 0.0, 2.5);
  TH1D *h_el_q2 = new TH1D(Form("h_el_q2_f%d",el_fract),Form("h_el_q2_f%d",el_fract),200, 0.0, 1.9);           
  TH2D *h_el_q2_w = new TH2D("h_el_q2_w","h_el_q2_w",200,0.0, 2.5, 200, 0.0, 1.50);
  TH2D *h_el_p_theta = new TH2D(Form("h_el_p_theta_f%s%d",prefix.c_str(), el_fract),Form("h_el_p_theta_f%s%d",prefix.c_str(),el_fract),100,0.0, lv_ebeam.E()+2, 100, 0.0, 80.0);
  TH2D *h_el_phi_theta = new TH2D(Form("h_el_phi_theta_f%s%d",prefix.c_str() , el_fract),Form("h_el_phi_theta_f%s%d",prefix.c_str(),el_fract), 200, -180.0, 180.0, 100, 0.0, 80.0);

  int f_max = 10;
  int ff = 0;
  
  while( ff < f_max ){

    ifstream infile;
    std::string line;
    //infile.open(Form("/work/clas12/bclary/CLAS12/electron_studies/elastic/lund_rad/elastic_fraction0p%d/clas12_2GeV_%d.lund",el_fract,ff));
    infile.open(Form("/work/clas12/bclary/CLAS12/electron_studies/elastic/lund_rad/elastic_test/clas12_2GeV_%d.lund",ff));
      
    while( std::getline(infile,line) ){
      double n0;
      double n1;
      double n2;
      double n3;
      double n4;
      double n5;
      double n6;
      double n7;
      double n8;
      double n9;
      double n10;
      double n11;
      double n12;
      double n13;
      std::istringstream ss(line);
    
      ss >> n0 >> n1 >> n2 >> n3 >> n4 >> n5 >> n6 >> n7 >> n8 >> n9 >> n10;
      if( n0 != 2 ){
	//std::cout << "n6 " << n6 << " n7 " << n7 << " n8 " << n8 << " " << n9 << std::endl;
      
	TLorentzVector lv_el;
	lv_el.SetPxPyPzE( n6, n7, n8, sqrt( n6*n6 + n7*n7 + n8*n8 + el_mass*el_mass));

	double w = (lv_ebeam + lv_target - lv_el).M();
	double el_theta = lv_el.Theta() * 180.0/TMath::Pi();
	double q2 = 4.0*lv_el.E()*lv_ebeam.E() * sin( lv_el.Theta()/2.0 ) * sin( lv_el.Theta()/2.0 );

	h_el_gen_w_theta->Fill(w, el_theta);
	h_el_w->Fill(w);
	h_el_q2_w->Fill(w,q2);

	h_el_p_theta->Fill( lv_el.P(), lv_el.Theta() * 180.0/TMath::Pi());
	h_el_phi_theta->Fill( lv_el.Phi()* 180.0/TMath::Pi(), lv_el.Theta() * 180.0/TMath::Pi());


      }

    }

    infile.close();
    ff++;
  }
    

  TCanvas *c1 = new TCanvas("c1","c1",900,900);
  c1->cd(0);
  gPad->SetLogz();
  h_el_gen_w_theta->Draw("colz");
  c1->SaveAs(Form("h_el_gen_w_theta_lund_f%s%d.pdf",prefix.c_str(),el_fract));
  
  TCanvas *c2 = new TCanvas("c2","c2",900,900);
  c2->cd(0);
  //gPad->SetLogy();
  h_el_w->Draw();
  c2->SaveAs(Form("h_el_gen_w_lund_f%s%d.pdf",prefix.c_str(),el_fract));

  TCanvas *c3 = new TCanvas("c3","c3",900,900);
  c3->cd(0);
  h_el_p_theta->Draw("colz");
  c3->SaveAs(Form("h_el_gen_p_theta_lund_f%s%d.pdf",prefix.c_str(),el_fract));

  TCanvas *c4 = new TCanvas("c4","c4",900,900);
  c4->cd(0);
  h_el_phi_theta->Draw("colz");
  c4->SaveAs(Form("h_el_gen_w_phi_theta_f%s%d.pdf",prefix.c_str(),el_fract));

  TCanvas *c5 = new TCanvas("c5","c5",900,900);
  c5->cd(0);
  gPad->SetLogz();
  h_el_q2_w->Draw("colz");
  c5->SaveAs("h_el_q2_w.pdf");
  
  TCanvas *c6 = new TCanvas("c6","c6",900,450);
  c6->Divide(2,10);
  c6->cd(1);
  h_el_w->Draw();
  c6->cd(2);
  h_el_q2->Draw();

  //fOut->Write();
  //  fOut->Close();
  

  return 0;
}
