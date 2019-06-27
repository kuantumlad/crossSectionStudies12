{
  gROOT->SetStyle("Plain");

  gStyle->SetStatStyle(1001); 
  gStyle->SetStatFont(62);  
  gStyle->SetStatBorderSize(2);

  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadLeftMargin(0.12);   
  //gStyle->SetPadRightMargin(0.15);   

  gStyle->SetOptStat(0000000);//1111111); //Integral, Overflow, Underflow, RMS, Mean, Nent, Name
  //gStyle->SetOptFit(0001);     //probability, Chi2, errors, name/values of parameters

 
  //gStyle->SetMarkerStyle(20);//Marker is a circle
  //gStyle->SetMarkerSize(.7);  
  //gStyle->SetMarkerColor(kGreen+3);
  //gStyle->SetLineColor(kRed);

  //gStyle->SetOptStat(0);

  //gStyle->SetPalette(kRainBow);
  //gStyle->SetPalette(kRainBow);
  gROOT->ForceStyle();

}
