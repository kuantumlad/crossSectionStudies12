#define elForwardDetectorCuts_hh

void loadECCuts(int run){
  
}


bool ele_default_PID_cut(int j){
  if(vpart_pid->at(j) == 11) return true;
  else return false;
}

bool ele_charge_cut(int j){
  if(part_charge[j] == -1) return true;
  else return false;
}


bool CC_nphe_cut(int j){

  double nphe_min = 2;

  if(part_CC_HTCC_nphe[j] > nphe_min) return true;
  else return false;
}


// EC cuts

bool EC_outer_vs_EC_inner_cut(int j){

  double edep_min = 0.06;

  if(part_Cal_PCAL_energy[j] > edep_min) return true; 
  else return false;

}


bool EC_sampling_fraction_cut(int j){

  double p_min = 1.5;

  if(Ebeam > 10) p_min = 1.5;
  if(Ebeam < 10) p_min = 1.0;
  if(Ebeam < 3)  p_min = 0.5;

  double sigma_range = 3;

  /// //////////////////////////////////////////////////////////////////////////////////////////////////////
  /// a) cut based on sampling fraction versus drift chamber momentum

  // 10.6 GeV (4013)
  double p0mean[] = {0.107005, 0.113632, 0.108363, 0.113069, 0.115156, 0.111107};
  double p1mean[] = {-0.346719, 0.125031, -0.0740992, 0.207526, 0.225827, 0.0156116};
  double p2mean[] = {0.00783146, 0.00368394, 0.00829978, 0.00337467, 0.00393947, 0.00910008};
  double p3mean[] = {-0.00071727, -0.000354281, -0.000875159, -0.000123021, -0.000144347, -0.000892991};
  double p0sigma[] = {0.0154995, 0.0247643, 0.00996904, 0.0182595, 0.0193756, 0.0120628};
  double p1sigma[] = {0.00566889, -0.00391144, 0.0138104, 0.0032902, 0.00060427, 0.011209};


  double mean = 0;
  double sigma = 0;
  double upper_lim_total = 0;
  double lower_lim_total = 0;

  for(Int_t k = 0; k < 6; k++){  
    if(part_Cal_PCAL_sector[j]-1 == k){
      mean = p0mean[k] *( 1 + part_p[j]/sqrt(pow(part_p[j],2) + p1mean[k])) + p2mean[k] * part_p[j] + p3mean[k] * pow(part_p[j],2);
      sigma = p0sigma[k] + p1sigma[k] / sqrt(part_p[j]);
      upper_lim_total = mean + sigma_range * sigma;
      lower_lim_total = mean - sigma_range * sigma;
    }
  }

  /// /////////////////////////////////////////////////////////////////////////////////////////////////////
  /// b) cut based on sampling fraction versus enrgy deposited in the calorimeter
  /*
  // 10.6 GeV
  double p0mean[] = {0.257819, 0.262394, 0.264319, 0.273215, 0.277912, 0.269706};
  double p1mean[] = {0.965033, 0.983467, 0.989758, 1.01724, 1.03319, 1.00817};
  double p2mean[] = {-0.0632092, -0.0290537, -0.0332319, -0.0940031, -0.103213, -0.0420508};
  double p3mean[] = {0.00797259, 0.000159148, 0.00050095, 0.0151809, 0.0168369, 0.000980376};
  double p0sigma[] = {0.0210694, 0.0209485, 0.0137895, 0.016642, 0.017343, 0.0162234};
  double p1sigma[] = {-0.00137114, -0.000904917, 0.00362368, 0.00243129, 0.00180574, 0.0026879};

  double mean = 0;
  double sigma = 0;
  double upper_lim_total = 0;
  double lower_lim_total = 0;

  for(Int_t k = 0; k < 6; k++){  
    if(part_Cal_PCAL_sector[j]-1 == k){
      mean = p0mean[k] *( p1mean[k] + p2mean[k] / part_Cal_energy_total[j] + p3mean[k] / pow(part_Cal_energy_total[j],2));
      sigma = p0sigma[k] + p1sigma[k] / sqrt(part_p[j]);
      upper_lim_total = mean + sigma_range * sigma;
      lower_lim_total = mean - sigma_range * sigma;
    }
  }
  */
  /// ////////////////////////////////////////////////////////////////////////////////////////////////////

  if(part_Cal_energy_total[j]/part_p[j] <= upper_lim_total && part_Cal_energy_total[j]/part_p[j] >= lower_lim_total && part_p[j] > p_min) return true;
  else return false;

}



bool EC_hit_position_fiducial_cut(int j){

  int sec_PCAL = part_Cal_PCAL_sector[j]-1;

  double x_PCAL = part_Cal_PCAL_x[j];
  double y_PCAL = part_Cal_PCAL_y[j];

  double x_PCAL_rot = y_PCAL * sin(sec_PCAL*60.0*Pival/180) + x_PCAL * cos(sec_PCAL*60.0*Pival/180);
  double y_PCAL_rot = y_PCAL * cos(sec_PCAL*60.0*Pival/180) - x_PCAL * sin(sec_PCAL*60.0*Pival/180);

  double angle_PCAL = 60;
  double height_PCAL = 47;   // PCAL starts at a hight of 39 + 8  (1.7 PCAL scintillator bars (each is 4.5 cm))  

  double slope_PCAL = 1/tan(0.5*angle_PCAL*Pival/180);
  double left_PCAL  = (height_PCAL - slope_PCAL * y_PCAL_rot);
  double right_PCAL = (height_PCAL + slope_PCAL * y_PCAL_rot);

  double radius2_PCAL = pow(height_PCAL+6,2)-pow(y_PCAL_rot,2);    // cut another 6 cm circle in the inner triangle tip to reject particles influenced by dead metarial

  if(x_PCAL_rot > left_PCAL && x_PCAL_rot > right_PCAL && pow(x_PCAL_rot,2) > radius2_PCAL && x_PCAL_rot < 371) return true;
  else return false;



// Cut using the natural directions of the fibers:
/*
  double u = part_Cal_PCAL_lu[j]/100;
  double v = part_Cal_PCAL_lv[j]/100;
  double w = part_Cal_PCAL_lw[j]/100;
   
  double min_u = 8;
  double max_u = 400;
  double min_v = 8;
  double max_v = 400;
  double min_w = 8;
  double max_w = 400;
   
  if(u > min_u && u < max_u && v > min_v && v < max_v && w > min_w && w < max_w) return true;
  else return false;
*/
}


// DC fiducial cuts for the 3 regions

bool DC_hit_position_region1_fiducial_cut(int j){

  double angle = 60; 
  double height = 19;

  int sec = part_DC_sector[j]-1;
    
  double x1_rot = part_DC_c1y[j] * sin(sec*60.0*Pival/180) + part_DC_c1x[j] * cos(sec*60.0*Pival/180);
  double y1_rot = part_DC_c1y[j] * cos(sec*60.0*Pival/180) - part_DC_c1x[j] * sin(sec*60.0*Pival/180);

  double slope = 1/tan(0.5*angle*Pival/180);
  double left  = (height - slope * y1_rot);
  double right = (height + slope * y1_rot);

  double radius2_DCr1 = pow(29,2)-pow(y1_rot,2);    // cut out the inner circle

  if (x1_rot > left && x1_rot > right && pow(x1_rot,2) > radius2_DCr1) return true;
  else return false;

}


bool DC_hit_position_region2_fiducial_cut(int j){

  double angle = 60; 
  double height = 38;

  int sec = part_DC_sector[j]-1;
    
  double x2_rot = part_DC_c2y[j] * sin(sec*60.0*Pival/180) + part_DC_c2x[j] * cos(sec*60.0*Pival/180);
  double y2_rot = part_DC_c2y[j] * cos(sec*60.0*Pival/180) - part_DC_c2x[j] * sin(sec*60.0*Pival/180);

  double slope = 1/tan(0.5*angle*Pival/180);
  double left  = (height - slope * y2_rot);
  double right = (height + slope * y2_rot);

  double radius2_DCr2 = pow(46,2)-pow(y2_rot,2);    // cut out the inner circle

  //return true;
  if (x2_rot > left && x2_rot > right && pow(x2_rot,2) > radius2_DCr2) return true;
  else return false;

}


bool DC_hit_position_region3_fiducial_cut(int j){

  double angle = 60; 
  double height = 38;

  int sec = part_DC_sector[j]-1;

  double x3_rot = part_DC_c3y[j] * sin(sec*60.0*Pival/180) + part_DC_c3x[j] * cos(sec*60.0*Pival/180);
  double y3_rot = part_DC_c3y[j] * cos(sec*60.0*Pival/180) - part_DC_c3x[j] * sin(sec*60.0*Pival/180);

  double slope = 1/tan(0.5*angle*Pival/180);
  double left  = (height - slope * y3_rot);
  double right = (height + slope * y3_rot);

  double radius2_DCr3 = pow(48,2)-pow(y3_rot,2);    // cut out the inner circle

  if (x3_rot > left && x3_rot > right && pow(x3_rot,2) > radius2_DCr3) return true;
  else return false;

}


bool DC_z_vertex_cut(int j){
 
  //  double vz_min_sect[] = {-12, -12, -12, -14, -12, -11};
  //double vz_max_sect[] = {10, 10, 10, 8, 10, 11};
  double vz_min_sect[] = {-32, -32, -32, -34, -32, -31};                                                                                                                                 
  double vz_max_sect[] = {30, 30, 30, 30, 30, 30};    
  double vz_min = 0;
  double vz_max = 0;

  for(Int_t k = 0; k < 6; k++){  
    if(part_Cal_PCAL_sector[j]-1 == k){
      vz_min = vz_min_sect[k];
      vz_max = vz_max_sect[k];
    }
  }

  if(part_vz[j] > vz_min && part_vz[j] < vz_max) return true;
  else return false;
}



bool Track_Quality_cut(int j){
  if(part_CC_HTCC_sector[j] > 0 && part_DC_sector[j] > 0 && part_Cal_PCAL_sector[j] > 0) return true;
  else return false;
}

