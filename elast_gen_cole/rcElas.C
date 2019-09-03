{

  int b_ntheta = 60;
  double theta_min = 5.0;
  double theta_max = 20.0;

  std::string noRC_file = "elast_gen_7GeV_NOrad_rc_bosted.root";//elas_gen_7GeV_NOradrc.root";///w/hallb-scifs17exp/clas12/bclary/CLAS12/electron_studies/elastic/radcorr/norad/7GeV/elast_gen_7GeV_noradcorr.root";
  TFile *noRC = new TFile(noRC_file.c_str());//"elast_gen_2GeV_10to49_norad.root");//elast_gen_2GeV_0_50_norad.root");//elast_gen_2GeV_norad.root");
  TTree *t1NORC = (TTree*)noRC->Get("h10");
  Float_t TheteNORC;
  Int_t ev;
  t1NORC->SetBranchAddress("Thete",&TheteNORC);

  TH1F *thetaNORCH   = new TH1F("thetaNORCH","thetaNORC",b_ntheta, theta_min, theta_max ); //25,5,30);

  //read all entries and fill the histograms
  Long64_t nentries = t1NORC->GetEntries();
  for (Long64_t i=0;i<nentries;i++) {
    t1NORC->GetEntry(i);
    thetaNORCH->Fill(TheteNORC);
  }

  std::string RC_file = "elast_gen_7GeV_rad_rc_bosted.root";//elas_gen_7GeV_radrc.root";///w/hallb-scifs17exp/clas12/bclary/CLAS12/electron_studies/elastic/radcorr/rad/7GeV/elast_gen_7GeV_radcorr.root";
  TFile *RC = new TFile(RC_file.c_str()); //"elast_gen_2GeV_10to49_rad.root");//elast_gen_2GeV_0_50_rad.root"); //elast_gen_2GeV_rad.root");
  TTree *t1RC = (TTree*)RC->Get("h10");
  Float_t TheteRC;
  Int_t evRC;
  t1RC->SetBranchAddress("Thete",&TheteRC);

  TH1F *thetaRCH   = new TH1F("thetaRCH","thetaRC",b_ntheta, theta_min, theta_max ); //25,5,30);

  //read all entries and fill the histograms                                                                                                                                                                                    
  Long64_t nentriesRC = t1RC->GetEntries();
  for (Long64_t i=0;i<nentriesRC;i++) {
    t1RC->GetEntry(i);    
    thetaRCH->Fill(TheteRC);
  }

  float timeNORC = 0.833;
  float timeRC = 0.528;
  //thetaNORCH->Scale(1/timeNORC);
  thetaRCH->Scale(0.25714);//1/timeRC);
  
  TCanvas *e = new TCanvas("e", "e", 10, 10, 800, 600);
  e->Divide(2, 1);
  e->cd(1);
  thetaNORCH->Draw();
  thetaRCH->SetLineColor(kRed); 
  thetaRCH->Draw("same");
  
  TH1F *rc = new TH1F("rc", "rc",b_ntheta, theta_min, theta_max ); //, 25, 5, 30);
  for (int t = 0; t < 100; t++){
    if (thetaNORCH->GetBinContent(t + 1) > 0) rc->SetBinContent(t + 1, thetaRCH->GetBinContent(t + 1)/thetaNORCH->GetBinContent(t + 1));
  }
  e->cd(2);
  rc->Draw();

  ofstream rcE("elasrc.dat");
  for (int t = 0; t < b_ntheta; t++){
    rcE << rc->GetBinContent(t + 1) << endl;
    std::cout << rc->GetBinCenter(t+1) << " " << rc->GetBinContent(t + 1) <<  std::endl;
  }
}
