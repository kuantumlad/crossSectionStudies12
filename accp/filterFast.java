import java.io.*;
import java.util.*;
import org.jlab.io.hipo.*;
import org.jlab.io.base.DataEvent;
import org.jlab.io.base.DataBank;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.groot.data.*;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.math.*;
import org.jlab.groot.fitter.*;

public class filterFast{

    private static final double PI = 3.14592658;
    private static final double LOG_SQRT_PI = Math.log(Math.sqrt(PI));
    private static final double I_SQRT_PI = 1 / Math.sqrt(PI);
    public static final int MAX_X = 20; // max value to represent exp(x)
 

    public static void main(String[] input){

	String dir_in = input[0];
	String s_run = input[1];
	int f_min = Integer.valueOf(input[2]);
	int f_max = Integer.valueOf(input[3]);
	int file_counter=f_min;

	int proton_pid = 2212;
	int electron_pid = 11;
	int kaon_plus_pid = 321;
	int kaon_minus_pid = -321;

	double mass_phi = 1.01946;
	double mass_kaon = 0.49368;
        double mass_electron = 0.000511;
	double mass_proton = 0.93827;
	double mass_pion = 0.1395;
	double mass_deut = 1.877;
			
	double radians_to_deg = 180./3.141592658;

	LorentzVector lv_beam = new LorentzVector(0,0,10.6,10.6);
	LorentzVector lv_target = new LorentzVector(0,0,0,mass_proton);

	double test_chi2 = 3.99;
	int test_dof = 3;
	double test_prob=pochisq(test_chi2, test_dof) ;///getChi2Prob(test_chi2,test_dof);

	System.out.println(" >> test chi2 "  + Double.toString(test_chi2) + " dof " + Integer.toString(test_dof) + " test prob " + test_prob );
	
	

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Histograms for Kinematics
	
	H2F h_el_p_theta = new H2F("h_el_p_theta","h_el_p_theta",200, 0.0, 11.0, 200, 0.0, 60.0 );
	H2F h_el_theta_phi = new H2F("h_el_p_theta_phi","h_el_p_theta_phi",200,-180.0, 180.0, 200, 0.0, 60.0 );

	H2F h_pr_p_theta = new H2F("h_pr_p_theta","h_pr_p_theta",200, 0.0, 6.0, 200, 0.0, 80.0 );
	H2F h_pr_theta_phi = new H2F("h_pr_p_theta_phi","h_pr_p_theta_phi",200, -180.0, 180.0, 200, 0.0, 60.0 );

	H2F h_kp_p_theta = new H2F("h_kp_p_theta","h_kp_p_theta",200, 0.0, 6.0, 200, 0.0, 80.0 );
	H2F h_kp_theta_phi = new H2F("h_kp_p_theta_phi","h_kp_p_theta_phi",200, -180.0, 180.0, 200, 0.0, 60.0 );

	H2F h_km_p_theta = new H2F("h_km_p_theta","h_kp_p_theta",200, 0.0, 6.0, 200, 0.0, 40.0 );
	H2F h_km_theta_phi = new H2F("h_km_p_theta_phi","h_kp_p_theta_phi",200, -180.0, 180.0, 200, 0.0, 60.0 );

	H2F h_el_p_theta_mc = new H2F("h_el_p_theta_mc","h_el_p_theta_mc",200, 0.0, 11.0, 200, 0.0, 60.0 );
	H2F h_el_theta_phi_mc = new H2F("h_el_p_theta_phi_mc","h_el_p_theta_phi_mc",200,-180.0, 180.0, 200, 0.0, 60.0 );

	H2F h_pr_p_theta_mc = new H2F("h_pr_p_theta_mc","h_pr_p_theta_mc",200, 0.0, 6.0, 200, 0.0, 80.0 );
	H2F h_pr_theta_phi_mc = new H2F("h_pr_p_theta_phi_mc","h_pr_p_theta_phi_mc",200, -180.0, 180.0, 200, 0.0, 60.0 );

	H2F h_kp_p_theta_mc = new H2F("h_kp_p_theta_mc","h_kp_p_theta_mc",200, 0.0, 6.0, 200, 0.0, 80.0 );
	H2F h_kp_theta_phi_mc = new H2F("h_kp_p_theta_phi_mc","h_kp_p_theta_phi_mc",200, -180.0, 180.0, 200, 0.0, 60.0 );

	H2F h_km_p_theta_mc = new H2F("h_km_p_theta_mc","h_kp_p_theta_mc",200, 0.0, 6.0, 200, 0.0, 40.0 );
	H2F h_km_theta_phi_mc = new H2F("h_km_p_theta_phi_mc","h_kp_p_theta_phi_mc",200, -180.0, 180.0, 200, 0.0, 60.0 );


	// Histogram for the final states
	H2F h_el_p_theta_final = new H2F("h_el_p_theta_final","h_el_p_theta_final",200, 0.0, 11.0, 200, 0.0, 60.0 );
	H2F h_el_theta_phi_final = new H2F("h_el_p_theta_phi_final","h_el_p_theta_phi_final",200,-180.0, 180.0, 200, 0.0, 60.0 );

	H2F h_pr_p_theta_final = new H2F("h_pr_p_theta_final","h_pr_p_theta_final",200, 0.0, 6.0, 200, 0.0, 40.0 );
	H2F h_pr_theta_phi_final = new H2F("h_pr_p_theta_phi_final","h_pr_p_theta_phi_final",200, -180.0, 180.0, 200, 0.0, 60.0 );

	H2F h_kp_p_theta_final = new H2F("h_kp_p_theta_final","h_kp_p_theta_final",200, 0.0, 6.0, 200, 0.0, 40.0 );
	H2F h_kp_theta_phi_final = new H2F("h_kp_p_theta_phi_final","h_kp_p_theta_phi_final",200, -180.0, 180.0, 200, 0.0, 60.0 );

	H2F h_km_p_theta_final = new H2F("h_km_p_theta_final","h_kp_p_theta_final",200, 0.0, 6.0, 200, 0.0, 40.0 );
	H2F h_km_theta_phi_final = new H2F("h_km_p_theta_phi_final","h_kp_p_theta_phi_final",200, -180.0, 180.0, 200, 0.0, 60.0 );
	
	Vector<H2F> h_pos_betap = new Vector<H2F>();
	Vector<H2F> h_neg_betap = new Vector<H2F>();
	for( int s = 0; s<6; s++ ){
	    h_pos_betap.add( new H2F("h_pos_betap_s"+Integer.toString(s),"h_pos_betap_s"+Integer.toString(s), 200, 0.0, 5.0, 200, 0.01, 1.1 ) );
	    h_neg_betap.add( new H2F("h_neg_betap_s"+Integer.toString(s),"h_neg_betap_s"+Integer.toString(s), 200, 0.0, 5.0, 200, 0.01, 1.1 ) );
	}

	H2F h_pr_final_betap = new H2F("h_pr_final_betap","h_pr_final_betap", 100, 0.0, 4.0, 100, 0.01, 1.1);
	H2F h_kp_final_betap = new H2F("h_kp_final_betap","h_kp_final_betap", 100, 0.0, 4.0, 100, 0.01, 1.1);
	H2F h_km_final_betap = new H2F("h_km_final_betap","h_km_final_betap", 100, 0.0, 4.0, 100, 0.01, 1.1);

	H1F h_pr_chi2_before = new H1F("h_pr_chi2_before","h_pr_chi2_before",500, -1.0, 15.0);
	H1F h_kp_chi2_before = new H1F("h_kp_chi2_before","h_kp_chi2_before",500, -1.0, 15.0);
	H1F h_km_chi2_before = new H1F("h_km_chi2_before","h_km_chi2_before",500, -1.0, 15.0);

	H1F h_pr_chi2_prob_before = new H1F("h_pr_chi2_prob_before","h_pr_chi2_prob_before", 200, 0.0, 1.1);
	H1F h_kp_chi2_prob_before = new H1F("h_kp_chi2_prob_before","h_kp_chi2_prob_before", 200, 0.0, 1.1);
	H1F h_km_chi2_prob_before = new H1F("h_km_chi2_prob_before","h_km_chi2_prob_before", 200, 0.0, 1.1);

 	H2F h_pr_chi2_p_before = new H2F("h_pr_chi2_p_before","h_pr_chi2_p_before",500,0.0, 15, 100, 0.0, 5.0);
 	H2F h_kp_chi2_p_before = new H2F("h_kp_chi2_p_before","h_kp_chi2_p_before",500,0.0,15, 100, 0.0, 5.0);
 	H2F h_km_chi2_p_before = new H2F("h_km_chi2_p_before","h_km_chi2_p_before",500,0.0, 15, 100, 0.0, 5.0);

 	H2F h_pr_chi2_theta_before = new H2F("h_pr_chi2_theta_before","h_pr_chi2_theta_before",500,0.0,15, 100, 0.0, 60.0);
 	H2F h_kp_chi2_theta_before = new H2F("h_kp_chi2_theta_before","h_kp_chi2_theta_before",500,0.0, 15, 100, 0.0, 60.0);
 	H2F h_km_chi2_theta_before = new H2F("h_km_chi2_theta_before","h_km_chi2_theta_before",500,0.0, 15, 100, 0.0, 60.0);

 	H2F h_pr_chi2_phi_before = new H2F("h_pr_chi2_phi_before","h_pr_chi2_phi_before",500, 0.0, 15, 180, -180.0, 180.0);
 	H2F h_kp_chi2_phi_before = new H2F("h_kp_chi2_phi_before","h_kp_chi2_phi_before",500, 0.0, 15, 180, -180.0, 180.0);
 	H2F h_km_chi2_phi_before = new H2F("h_km_chi2_phi_before","h_km_chi2_phi_before",500, 0.0, 15, 180, -180.0, 180.0);

	H1F h_pr_chi2_final = new H1F("h_pr_chi2_final","h_pr_chi2_final",500, 0.0, 15.0);
	H1F h_kp_chi2_final = new H1F("h_kp_chi2_final","h_kp_chi2_final",500, 0.0, 15.0);
	H1F h_km_chi2_final = new H1F("h_km_chi2_final","h_km_chi2_final",500, 0.0, 15.0);

	H1F h_pr_chi2_prob_final = new H1F("h_pr_chi2_prob_final","h_pr_chi2_prob_final", 200, 0.0, 1.1);
	H1F h_kp_chi2_prob_final = new H1F("h_kp_chi2_prob_final","h_kp_chi2_prob_final", 200, 0.0, 1.1);
	H1F h_km_chi2_prob_final = new H1F("h_km_chi2_prob_final","h_km_chi2_prob_final", 200, 0.0, 1.1);
		    
	H2F h_pr_chi2_p_final = new H2F("h_pr_chi2_p_final","h_pr_chi2_p_final",500,0.0, 15, 100, 0.0, 5.0);  
	H2F h_kp_chi2_p_final = new H2F("h_kp_chi2_p_final","h_kp_chi2_p_final",500,0.0, 15, 100, 0.0, 5.0);  
	H2F h_km_chi2_p_final = new H2F("h_km_chi2_p_final","h_km_chi2_p_final",500,0.0, 15, 100, 0.0, 5.0);  

	H2F h_pr_chi2_theta_final = new H2F("h_pr_chi2_theta_final","h_pr_chi2_theta_final",500,0.0, 15, 100, 0.0, 60.0);  
	H2F h_kp_chi2_theta_final = new H2F("h_kp_chi2_theta_final","h_kp_chi2_theta_final",500,0.0, 15, 100, 0.0, 60.0);  
	H2F h_km_chi2_theta_final = new H2F("h_km_chi2_theta_final","h_km_chi2_theta_final",500,0.0, 15, 100, 0.0, 60.0);  

	H2F h_pr_chi2_phi_final = new H2F("h_pr_chi2_phi_final","h_pr_chi2_phi_final",500,0.0, 15, 180, -180.0, 180.0);
	H2F h_kp_chi2_phi_final = new H2F("h_kp_chi2_phi_final","h_kp_chi2_phi_final",500,0.0, 15, 180, -180.0, 180.0);    
	H2F h_km_chi2_phi_final = new H2F("h_km_chi2_phi_final","h_km_chi2_phi_final",500,0.0, 15, 180, -180.0, 180.0);    

	H1F h_phi_final_chi2_final_c1 = new H1F("h_phi_final_chi2_final_c1","h_phi_final_chi2_final_c1",75,0.8,1.5);
	H1F h_phi_final_chi2_final_c2 = new H1F("h_phi_final_chi2_final_c2","h_phi_final_chi2_final_c2",75,0.8,1.5);
	H1F h_phi_final_chi2_final_c3 = new H1F("h_phi_final_chi2_final_c3","h_phi_final_chi2_final_c3",75,0.8,1.5);
	H1F h_phi_final_chi2_final_c4 = new H1F("h_phi_final_chi2_final_c4","h_phi_final_chi2_final_c4",75,0.8,1.5);
	
	Vector<H1F> v_pr_chi2 = new Vector<H1F>();
	Vector<H1F> v_kp_chi2 = new Vector<H1F>();
	Vector<H1F> v_km_chi2 = new Vector<H1F>();
	
	for( int i = 0; i < 15; i++ ){
	    v_pr_chi2.add( new H1F("h_pr_chi2_dof"+Integer.toString(i),"h_pr_chi2_dof"+Integer.toString(i),100,0.0, 15.0));
	    v_kp_chi2.add( new H1F("h_kp_chi2_dof"+Integer.toString(i),"h_kp_chi2_dof"+Integer.toString(i),100,0.0, 15.0));
	    v_km_chi2.add( new H1F("h_km_chi2_dof"+Integer.toString(i),"h_km_chi2_dof"+Integer.toString(i),100,0.0, 15.0));
	}


	H2F h_q2w = new H2F("h_q2w","h_q2w",200,0.0, 9.0, 200, 0.0, 6.0);
	H2F h_q2x = new H2F("h_q2x","h_q2x",200,0.0, 1.0, 200, 0.0, 9.0);
	H2F h_q2t = new H2F("h_q2t","hq2t",200,0.0, 5.0, 200, 0.0, 9.0);
	H1F h_t = new H1F("h_t","h_t",70,0.0, 5.0);
	H2F h_xbt = new H2F("h_xbt","h_xbt",70,0.0, 1.0, 70, 0.0, 5.0);

	H2F h_q2w_mc = new H2F("h_q2w_mc","h_q2w_mc",200,0.0, 9.0, 200, 0.0, 6.0);
	H2F h_q2x_mc = new H2F("h_q2x_mc","h_q2x_mc",200,0.0, 1.0, 200, 0.0, 9.0);
	H2F h_q2t_mc = new H2F("h_q2t_mc","hq2t_mc",200,0.0, 5.0, 200, 0.0, 9.0);
	H1F h_t_mc = new H1F("h_t_mc","h_t_mc",70,0.0, 5.0);
	H2F h_xbt_mc = new H2F("h_xbt_mc","h_xbt_mc",70,0.0, 1.0, 70, 0.0, 5.0);

	H2F h_q2w_final = new H2F("h_q2w_final","h_q2w_final",200,0.0, 9.0, 200, 0.0, 6.0);
	H2F h_q2x_final = new H2F("h_q2x_final","h_q2x_final",200,0.0, 1.0, 200, 0.0, 9.0);
	H2F h_q2t_final = new H2F("h_q2t_final","hq2t_final",200,0.0, 5.0, 200, 0.0, 9.0);
	H1F h_t_final = new H1F("h_t_final","h_t_final",70,0.0, 5.0);
	H2F h_xbt_final = new H2F("h_xbt_final","h_xbt_final",70,0.0, 1.0, 70, 0.0, 5.0);
	H1F h_kpkm_theta = new H1F("h_kpkm_costheta","h_kpkm_costheta",100, -1.1, 1.1 );
	H1F h_kpkm_theta_final = new H1F("h_kpkm_costheta_final","h_kpkm_costheta_final",100, -1.1, 1.1 );
			
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Histograms for Missing Mass	
	H1F h_el_eX = new H1F("h_el_eX","h_el_eX",200,-1.0, 6.0);
	H1F h_el_epX = new H1F("h_el_epX","h_el_epX",200,-0.50, 4.0);
	H1F h_el_epkpX = new H1F("h_el_epkpX","h_el_epkpX",200,-0.50, 2.20);
	H1F h_el_epkmX = new H1F("h_el_epkmX","h_el_epkmX",200,-0.50, 2.20);
	H1F h_el_epkpkmX = new H1F("h_el_epkpkmX","h_el_epkpkmX",75,-3.0, 3.0);
	H1F h_el_epkpkmMM2 = new H1F("h_el_epkpkmMM2","h_el_epkpkmMM2",75,-1.0, 6.0);
	H1F h_el_ekpX = new H1F("h_el_ekpX","h_el_ekpX",100,0.0, 4.0);
	H1F h_el_ekpkmX = new H1F("h_el_ekpkmX","h_el_ekpkmX",75,0.50, 3.80);
	H1F h_prelim_phi = new H1F("h_prelim_phi","h_prelim_phi",75,0.8,1.5);
	H1F h_miss_per = new H1F("h_miss_perp","h_miss_perp",75, -0.5, 1.5);
	H1F h_colin_pr = new H1F("h_colin_pr","h_colin_pr",100,-10, 60.0);
	H1F h_colin_kp = new H1F("h_colin_kp","h_colin_kp",100,-10, 60.0);
	H1F h_colin_km = new H1F("h_colin_km","h_colin_km",100,-10, 60.0);
	H1F h_ang_elpr = new H1F("h_ang_elpr","h_ang_elpr",100,-10, 60.0);
	H1F h_ang_kpkm = new H1F("h_ang_kpkm","h_ang_kpkm",100,-10, 60.0);

	H1F h_el_eX_final = new H1F("h_el_eX_final","h_el_eX_final",200,-1.0, 6.0);
	H1F h_el_epX_final = new H1F("h_el_epX_final","h_el_epX_final",200,-0.50, 4.0);
	H1F h_el_epkpX_final = new H1F("h_el_epkpX_final","h_el_epkpX_final",200,-0.50, 2.20);
	H1F h_el_epkmX_final = new H1F("h_el_epkmX_final","h_el_epkmX_final",200,-0.50, 2.20);
	H1F h_el_epkpkmX_final = new H1F("h_el_epkpkmX_final","h_el_epkpkmX_final",75,-3.0, 3.0);
	H1F h_el_epkpkmMM2_final = new H1F("h_el_epkpkmMM2_final","h_el_epkpkmMM2_final",75,-1.0, 6.0);
	H1F h_el_ekpX_final = new H1F("h_el_ekpX_final","h_el_ekpX_final",100,0.0, 4.0);
	H1F h_el_ekpkmX_final = new H1F("h_el_ekpkmX_final","h_el_ekpkmX_final",75,0.50, 3.80);
	H1F h_miss_per_final = new H1F("h_miss_perp_final","h_miss_perp_final",75, -0.5, 1.5);
	H1F h_colin_pr_final = new H1F("h_colin_pr_final","h_colin_pr_final",100,-10, 60.0);
	H1F h_colin_kp_final = new H1F("h_colin_kp_final","h_colin_kp_final",100,-10, 60.0);
	H1F h_colin_km_final = new H1F("h_colin_km_final","h_colin_km_final",100,-10, 60.0);
	H1F h_ang_elpr_final = new H1F("h_ang_elpr_final","h_ang_elpr_final",100,-10, 60.0);
	H1F h_ang_kpkm_final = new H1F("h_ang_kpkm_final","h_ang_kpkm_final",100,-10, 60.0);

	//Histograms for MC
	H1F h_res_q2 = new H1F("h_res_q2","h_res_q2",100, -0.50, 0.50 );
	H1F h_res_w = new H1F("h_res_w","h_res_w",100, -0.50, 0.50 );
	H1F h_res_t = new H1F("h_res_t","h_res_t",100, -0.50, 0.50 );
	H1F h_res_xb = new H1F("h_res_xb","h_res_xb",100, -0.50, 0.50 );

	H2F h2_res_q2 = new H2F("h2_res_q2","h2_res_q2",100, 0.0, 10.75, 100, -0.50, 0.50 );
	H2F h2_res_w = new H2F("h2_res_w","h2_res_w",100, 0.0, 15.0, 100, -0.50, 0.50 );
	H2F h2_res_t = new H2F("h2_res_t","h2_res_t",100, 0.0, 3.0, 100, -0.50, 0.50 );
	H2F h2_res_xb = new H2F("h2_res_xb","h2_res_xb",100, 0.0, 1.0, 100, -0.50, 0.50 );

	/////////// 
	// Histogram for acceptance
	int n_bins = 100;
	double min_q2 = -0.01;
	double max_q2 = 14.0; // using limits in pac proposal paper

	int n_bins_t = 50;
	double min_t = 0;
	double max_t = 3.0;
 	H1F h_accep_gen_q2 = new H1F("h_accep_gen_q2","h_accep_gen_q2", n_bins, min_q2, max_q2 );
  	H1F h_accep_rec_q2 = new H1F("h_accep_rec_q2","h_accep_rec_q2", n_bins, min_q2, max_q2 );
	H1F h_accp_q2 = new H1F("h_accp_q2","h_accp_q2",n_bins,min_q2, max_q2);

	H1F h_accep_gen_t = new H1F("h_accep_gen_t","h_accep_gen_t", n_bins, min_t, max_t );
  	H1F h_accep_rec_t = new H1F("h_accep_rec_t","h_accep_rec_t", n_bins, min_t, max_t );
	H1F h_accp_t = new H1F("h_accp_t","h_accp_t",n_bins,min_t, max_t);
	
	// Run information at 
	// -- /work/clas12/clas12/data/trains/v2/skim4_inclusive/skim4_
	Vector<Integer> run_list = new Vector<Integer>();
	/*run_list.add(5000);
	run_list.add(5001);
	run_list.add(5030);
	run_list.add(5036);
	run_list.add(5038);
	run_list.add(5046);
	run_list.add(5117);*/
	// outbending data below
	//run_list.add(5532);
	//run_list.add(5533);
	run_list.add(5534);
	//run_list.add(5535);
	//run_list.add(5536);
	//run_list.add(5537);
	//run_list.add(5538);

	String dataType="DATA";


	Vector< Vector<LorentzVector> > final_event = new Vector< Vector<LorentzVector > >();
	Vector<Integer> v_helicity = new Vector<Integer>();

	int file_counter_limit=-1;
	if ( dataType == "MC" ){
	    file_counter_limit = f_max;
	}
	else{
	    file_counter_limit = run_list.size();
	}

	while( file_counter < file_counter_limit ){
	    
	    System.out.println("Opening File ");
	    
	    String file_in=null;
	    if (dataType != "MC" ){
		file_in = dir_in + Integer.toString(run_list.get(file_counter)) + ".hipo";
	    }
	    else{
		file_in = dir_in + Integer.toString(file_counter) + ".txt.hipo";
	    }

	    System.out.println(" >> OPENING FILE " + file_in );
	    File fin_temp = new File(file_in);
	    if ( !fin_temp.exists() ){
		file_counter++;
		continue;
	    }

	    if( !fin_temp.exists() ) continue;
	    //if(f.exists() && !f.isDirectory())

	    System.out.println(" Processing File " + file_in );
	    
	    HipoDataSource hiporeader = new HipoDataSource();
	    hiporeader.open(new File(file_in) );
	    DataEvent event = null;

	    int max_event = 10000;/// hiporeader.getSize();
	    int num_ev = 0;
	    System.out.println(" PROCESSING " + Integer.toString(max_event) + " EVENTS " ) ;

	    while( num_ev < max_event ){

		//System.out.println(" event " + Integer.toString(num_ev) );
		//while( hiporeader.hasEvent()){
		event = (DataEvent)hiporeader.gotoEvent(num_ev);
		//DataEvent event = (DataEvent)hiporeader.getNextEvent();
		num_ev++;

		boolean runConfig_pres = event.hasBank("RUN::config");
		boolean recBank_pres  = event.hasBank("REC::Particle");
		boolean rawScalerBank_pres = event.hasBank("RAW::scaler");
		boolean eventBank_pres = event.hasBank("REC::Event");
		boolean configBank_pres = event.hasBank("RUN::config");
		boolean recScint_pres = event.hasBank("REC::Scintillator");
		boolean recTrack_pres = event.hasBank("REC::Track");
		boolean mcBank_pres = event.hasBank("MC::Particle");

		int percent_complete = (int)Math.floor((double)num_ev/(double)max_event * 100.0);
		if( percent_complete%10 == 0 && percent_complete!=0){
		    //System.out.println(Integer.toString(percent_complete) + " complete " );
		}
	
		LorentzVector lv_el_gen = new LorentzVector(0,0,0,0);
		LorentzVector lv_pr_gen = new LorentzVector(0,0,0,0);
		LorentzVector lv_kp_gen = new LorentzVector(0,0,0,0);
		LorentzVector lv_km_gen = new LorentzVector(0,0,0,0);
		
		LorentzVector lv_el_rec = new LorentzVector(0,0,0,0);
		LorentzVector lv_pr_rec = new LorentzVector(0,0,0,0);
		LorentzVector lv_kp_rec = new LorentzVector(0,0,0,0);
		LorentzVector lv_km_rec = new LorentzVector(0,0,0,0);



		if( eventBank_pres ){
		    DataBank eventBank = event.getBank("REC::Event");
		    int helicity = eventBank.getInt("Helic",0);
		    

		    if( recBank_pres && recScint_pres && recTrack_pres ){
			
			//System.out.println(" Event Bank and REC Particle Bank present " );
			
			Vector<LorentzVector> v_el = new Vector<LorentzVector>();
			Vector<LorentzVector> v_pr = new Vector<LorentzVector>();
			Vector<LorentzVector> v_kp = new Vector<LorentzVector>();
			Vector<LorentzVector> v_km = new Vector<LorentzVector>();

			DataBank recbank = event.getBank("REC::Particle");		
			DataBank scintbank = event.getBank("REC::Scintillator");
			DataBank trackbank = event.getBank("REC::Track");

			LorentzVector lv_el = new LorentzVector(0,0,0,0);
			LorentzVector lv_pr = new LorentzVector(0,0,0,0);
			LorentzVector lv_kp = new LorentzVector(0,0,0,0);
			LorentzVector lv_km = new LorentzVector(0,0,0,0);

			double pr_beta_final = 0.0;
			double kp_beta_final = 0.0;
			double km_beta_final = 0.0;

 			double pr_chi2_final = 0.0;
			double kp_chi2_final = 0.0;
			double km_chi2_final = 0.0;

			double pr_chi2 = 0.0;
			double kp_chi2 = 0.0;
			double km_chi2 = 0.0;

			int pr_ndf = 0;
			int kp_ndf = 0;
			int km_ndf = 0;
			

			int rec_pr_index = -1;
			int rec_kp_index = -1;
			int rec_km_index = -1;

			//System.out.println(" >> Number of rec particles " + Integer.toString(event.getBank("REC::Particle").rows()) );
		 	for( int rec_i = 0; rec_i < event.getBank("REC::Particle").rows(); rec_i++ ){ 			    
			    int charge = recbank.getInt("charge",rec_i);
			    if( charge != 0 ){
				int rec_pid = recbank.getInt("pid",rec_i);			    
 				float rec_px = recbank.getFloat("px",rec_i);
				float rec_py = recbank.getFloat("py",rec_i);
				float rec_pz = recbank.getFloat("pz",rec_i);
				float rec_vz = recbank.getFloat("vz",rec_i);
				
				int rec_status = recbank.getInt("status",rec_i);
				float beta_clas12 = recbank.getFloat("beta",rec_i);
				int scint_sector = -1;
				//if( 2000 >= rec_status && rec_status >= 3000 ) continue;
				for( int i = 0; i < scintbank.rows(); i++){			   
				    int pindex = scintbank.getShort("pindex",i);
				    if( pindex == rec_i ){					
					scint_sector = scintbank.getInt("sector",i) - 1;
					break;
				    }
				}

							
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				// EB PID
				boolean el_FD = false;
			 	if( rec_pid == electron_pid ){
				    LorentzVector lv_temp = new LorentzVector(0,0,0,0);
				    double rec_e = Math.sqrt( rec_px*rec_px + rec_py*rec_py + rec_pz*rec_pz + mass_electron*mass_electron );
				    lv_temp.setPxPyPzE(rec_px,rec_py,rec_pz,rec_e);
				    if( true ){// (lv_temp.theta() * radians_to_deg)  > 6.0 ){
					v_el.add(lv_temp);
					h_el_p_theta.fill( lv_temp.p(), lv_temp.theta() * radians_to_deg );
					h_el_theta_phi.fill( lv_temp.phi() * radians_to_deg, lv_temp.theta() * radians_to_deg );
					if( lv_temp.e() > lv_el.e() ){
					    lv_el = lv_temp;
					    lv_el_rec = lv_el;
					    //System.out.println(" >> electron index " + Integer.toString(rec_i) +  " momentum  " + Double.toString(lv_el.p()) );
					}
				    }
				    				
				}				
				else if (rec_pid == proton_pid ){
				    LorentzVector lv_temp = new LorentzVector(0,0,0,0);
				    double rec_e = Math.sqrt( rec_px*rec_px + rec_py*rec_py + rec_pz*rec_pz + mass_proton*mass_proton );
				    lv_temp.setPxPyPzE(rec_px,rec_py,rec_pz,rec_e);
				    v_pr.add(lv_temp);
				    h_pr_p_theta.fill( lv_temp.p(), lv_temp.theta() * radians_to_deg );
				    h_pr_theta_phi.fill( lv_temp.phi() * radians_to_deg, lv_temp.theta() * radians_to_deg );
				    if( lv_temp.e() > lv_pr.e() ){
					lv_pr = lv_temp; pr_beta_final = beta_clas12;
					lv_pr_rec = lv_pr;
					//System.out.println(" >> proton index " + Integer.toString(rec_i) +  " momentum  " + Double.toString(lv_pr.p()) );
					rec_pr_index=rec_i;
				    }
				}
				else if (rec_pid == kaon_plus_pid ){
				    LorentzVector lv_temp = new LorentzVector(0,0,0,0);
				    double rec_e = Math.sqrt( rec_px*rec_px + rec_py*rec_py + rec_pz*rec_pz + mass_kaon*mass_kaon );
				    lv_temp.setPxPyPzE(rec_px,rec_py,rec_pz,rec_e);
				    v_kp.add(lv_temp);
				    h_kp_p_theta.fill( lv_temp.p(), lv_temp.theta() * radians_to_deg );
				    h_kp_theta_phi.fill( lv_temp.phi() * radians_to_deg, lv_temp.theta() * radians_to_deg );
				    if( lv_temp.e() > lv_kp.e() ){
					lv_kp = lv_temp; kp_beta_final = beta_clas12;
					lv_kp_rec = lv_kp;					//System.out.println(" >> kaon plus index " + Integer.toString(rec_i) +  " momentum  " + Double.toString(lv_kp.p()) );
					rec_kp_index=rec_i;
				    }
				}
				else if (rec_pid == kaon_minus_pid){
				    LorentzVector lv_temp = new LorentzVector(0,0,0,0);
				    double rec_e = Math.sqrt( rec_px*rec_px + rec_py*rec_py + rec_pz*rec_pz + mass_kaon*mass_kaon );
				    lv_temp.setPxPyPzE(rec_px,rec_py,rec_pz,rec_e);
				    v_km.add(lv_temp);
				    h_km_p_theta.fill( lv_temp.p(), lv_temp.theta() * radians_to_deg );
				    h_km_theta_phi.fill( lv_temp.phi() * radians_to_deg, lv_temp.theta() * radians_to_deg );
				    if( lv_temp.e() > lv_km.e() ){ 
					lv_km = lv_temp; km_beta_final = beta_clas12;
					lv_km_rec = lv_km;
					//System.out.println(" >> kaon minus index " + Integer.toString(rec_i) +  " momentum  " + Double.toString(lv_km.p()) );
					rec_km_index=rec_i;
				    }
				}			    	      

				if( scint_sector >= 0  && scint_sector < 6){
				    LorentzVector lv_temp = new LorentzVector(0,0,0,0);
				    lv_temp.setPxPyPzE(rec_px,rec_py,rec_pz,0);				    
				    if( charge > 0 ){
					h_pos_betap.get(scint_sector).fill(lv_temp.p(), beta_clas12);
				    }
				    else if( charge < 0 ){
					h_neg_betap.get(scint_sector).fill(lv_temp.p(), beta_clas12);
				    }
				}
						
			    }
			}
			
			for( int i = 0; i < trackbank.rows(); i++ ){
			    int pindex = trackbank.getShort("pindex",i);
			    if( pindex == rec_pr_index ){
 				pr_ndf = trackbank.getShort("NDF",i);
				pr_chi2 = trackbank.getFloat("chi2",i)/(double)pr_ndf;
			    }					
			    else if( pindex == rec_kp_index ){
 				kp_ndf = trackbank.getShort("NDF",i);
				kp_chi2 = trackbank.getFloat("chi2",i)/(double)pr_ndf;
			    }					
			    else if( pindex == rec_km_index ){
 				km_ndf = trackbank.getShort("NDF",i);
				km_chi2 = trackbank.getFloat("chi2",i)/(double)pr_ndf;
			    }					
			}
			
			System.out.println(" >> reduced chi2 values are " + pr_chi2 + " " + kp_chi2 + " " + km_chi2 );

 			if( pr_chi2 > 0 ){
			    h_pr_chi2_before.fill(pr_chi2);
			    h_pr_chi2_p_before.fill(pr_chi2,lv_pr.p());
			    h_pr_chi2_theta_before.fill(pr_chi2,lv_pr.theta()  * radians_to_deg);			    
			    h_pr_chi2_phi_before.fill(pr_chi2, lv_pr.phi() * radians_to_deg);
			    double chi2_prob = getChi2Prob(pr_chi2, pr_ndf);
			    h_pr_chi2_prob_before.fill( chi2_prob );
			    
 			    if ( pr_ndf < 15 && pr_ndf > 0 ){
				v_pr_chi2.get(pr_ndf).fill(pr_chi2);
			    }


			}
			if( kp_chi2 > 0 ){
			    h_kp_chi2_before.fill(kp_chi2);
			    h_kp_chi2_p_before.fill(kp_chi2,lv_kp.p());
			    h_kp_chi2_theta_before.fill(kp_chi2,lv_kp.theta()  * radians_to_deg);
			    h_kp_chi2_phi_before.fill(kp_chi2, lv_kp.phi() * radians_to_deg);
			    double chi2_prob = getChi2Prob(kp_chi2, kp_ndf);
			    h_kp_chi2_prob_before.fill( chi2_prob );
 			    if ( kp_ndf < 15 && kp_ndf > 0 ){
				v_kp_chi2.get(kp_ndf).fill(kp_chi2);
			    }

			}
			if( km_chi2 > 0 ){
			    h_km_chi2_before.fill(km_chi2);
			    h_km_chi2_p_before.fill(km_chi2,lv_km.p());
			    h_km_chi2_theta_before.fill(km_chi2,lv_km.theta()  * radians_to_deg);
			    h_km_chi2_phi_before.fill(km_chi2, lv_km.phi() * radians_to_deg);
			    double chi2_prob = getChi2Prob(km_chi2, km_ndf);
			    h_km_chi2_prob_before.fill( chi2_prob );
 			    if ( km_ndf < 15 && km_ndf > 0 ){
				v_km_chi2.get(km_ndf).fill(km_chi2);
			    }
			}
			    

			



							
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// EB Event Selection
			

			//////////////////////////////////////////////////
			// get most energetic particles here
			// get missing mass LV here
			// cut on Q2 and W
			LorentzVector lv_eX= new LorentzVector(0,0,0,0);
 			LorentzVector lv_epX= new LorentzVector(0,0,0,0);
			LorentzVector lv_epkpX= new LorentzVector(0,0,0,0);
			LorentzVector lv_epkmX= new LorentzVector(0,0,0,0);
			LorentzVector lv_epkpkmX= new LorentzVector(0,0,0,0);
			LorentzVector lv_ekpX= new LorentzVector(0,0,0,0);
			LorentzVector lv_ekpkmX= new LorentzVector(0,0,0,0);

			lv_eX.add(lv_beam);
			lv_eX.add(lv_target);
			lv_eX.sub(lv_el);
			h_el_eX.fill(lv_eX.mass());

			boolean high_e_el = lv_el.p() > 1.5;
			boolean p_pr = lv_pr.p() > 0.5 && lv_pr.p() < 4.0;
			boolean low_theta_kp = lv_kp.theta()  * radians_to_deg < 35;
			boolean low_theta_km = lv_km.theta()  * radians_to_deg < 40;
			boolean low_p_kp  = lv_kp.p() < 2.5;
			boolean low_p_km  = lv_km.p() < 2.5;

			double q2 = 4*lv_beam.e()*lv_el.e()*Math.sin(lv_el.theta()/2.0)*Math.sin(lv_el.theta()/2.0);
			double xb = q2 / (2.0*mass_proton*(lv_beam.e() - lv_el.e()));
			double w = Math.sqrt(-q2 + mass_proton*mass_proton + 2*mass_proton*(lv_beam.e() - lv_el.e()));
			double t = 2.0*mass_proton*(lv_pr.e() - mass_proton);

			boolean w_cut = w > 2.0;
			boolean q2_cut = q2 > 1.0;
			
			if( w_cut && q2_cut ){
			
			    h_q2w.fill(q2,w);
			    h_q2x.fill(xb,q2);
			    h_q2t.fill(t,q2);
			    h_t.fill(t);
			    h_xbt.fill(xb,t);

			    if( v_pr.size() > 0  && high_e_el){
				lv_epX.add(lv_beam);
				lv_epX.add(lv_target);
				lv_epX.sub(lv_el);
				lv_epX.sub(lv_pr);
				h_el_epX.fill(lv_epX.mass2());
			    }

			    if( v_pr.size() > 0 && v_kp.size() > 0 && high_e_el){
				lv_epkpX.add(lv_beam);
				lv_epkpX.add(lv_target);
				lv_epkpX.sub(lv_el);
				lv_epkpX.sub(lv_pr);
				lv_epkpX.sub(lv_kp);
				h_el_epkpX.fill(lv_epkpX.mass2());
			    }
			    
			    if( v_pr.size() > 0 && v_km.size() > 0 && high_e_el){
				lv_epkmX.add(lv_beam);
				lv_epkmX.add(lv_target);
				lv_epkmX.sub(lv_el);
				lv_epkmX.sub(lv_pr);
				lv_epkmX.sub(lv_km);
				h_el_epkmX.fill(lv_epkmX.mass2());
			    }

			    if( v_pr.size() > 0 && v_kp.size() > 0 && v_km.size() > 0 && high_e_el){
				lv_epkpkmX.add(lv_beam);
				lv_epkpkmX.add(lv_target);
				lv_epkpkmX.sub(lv_el);
				lv_epkpkmX.sub(lv_pr);
				lv_epkpkmX.sub(lv_kp);
				lv_epkpkmX.sub(lv_km);
				h_el_epkpkmX.fill(lv_epkpkmX.mass2());
 				h_el_epkpkmMM2.fill(lv_epkpkmX.mass2());
			    }

			    if( v_kp.size() > 0 && high_e_el){
				lv_ekpX.add(lv_beam);
				lv_ekpX.add(lv_target);
				lv_ekpX.sub(lv_kp);
				h_el_ekpX.fill(lv_ekpX.mass2());
			    }

			    if( v_kp.size() > 0 && v_km.size() > 0 && high_e_el ){
				lv_ekpkmX.add(lv_beam);
				lv_ekpkmX.add(lv_target);
				lv_ekpkmX.sub(lv_el);
				lv_ekpkmX.sub(lv_kp);
				lv_ekpkmX.sub(lv_km);
				h_el_ekpkmX.fill(lv_ekpkmX.mass2());
			    }

			    if(  v_pr.size() > 0 && v_kp.size() > 0 && v_km.size() > 0 && high_e_el ){
				if( low_theta_kp && low_theta_km && low_p_kp && low_p_km && p_pr ){

 				    double ang_elpr = Math.acos( lv_pr.vect().dot(lv_el.vect()) / (lv_pr.vect().mag() * lv_el.vect().mag() ) ) * radians_to_deg;
				    double ang_kpkm = Math.acos( lv_kp.vect().dot(lv_km.vect()) / (lv_kp.vect().mag() * lv_km.vect().mag() ) ) * radians_to_deg;

				    double colin_pr_ang = Math.acos( lv_pr.vect().dot(lv_ekpkmX.vect()) / (lv_pr.vect().mag() * lv_ekpkmX.vect().mag() ) ) * radians_to_deg;
				    double colin_km_ang = Math.acos( lv_km.vect().dot(lv_epkpX.vect()) / (lv_km.vect().mag() * lv_epkpX.vect().mag() ) ) * radians_to_deg;
				    double colin_kp_ang = Math.acos( lv_kp.vect().dot(lv_epkmX.vect()) / (lv_kp.vect().mag() * lv_epkmX.vect().mag() ) ) * radians_to_deg;

				    boolean low_epkpkm_mm2 = Math.pow(lv_epkpkmX.mass(),2) < 0.6;
				    boolean miss_perp = Math.sqrt(lv_epkpkmX.px()*lv_epkpkmX.px() + lv_epkpkmX.py()*lv_epkpkmX.py()) < 0.5;
				    boolean colin_pr = colin_pr_ang < 30;
				    boolean colin_km =  colin_km_ang < 20;
				    boolean colin_kp = colin_kp_ang < 20;
				    
				    h_miss_per.fill(Math.sqrt(lv_epkpkmX.px()*lv_epkpkmX.px() + lv_epkpkmX.py()*lv_epkpkmX.py()));
 				    h_colin_pr.fill(colin_pr_ang);
				    h_colin_km.fill(colin_km_ang);
				    h_colin_kp.fill(colin_kp_ang);
				    //System.out.println(" >> final reduced chi2 values are " + pr_chi2 + " " + kp_chi2 + " " + km_chi2 );

				    if( pr_chi2 > 0 ){
					h_pr_chi2_final.fill(pr_chi2);
					h_pr_chi2_p_final.fill(pr_chi2,lv_pr.p());
					h_pr_chi2_theta_final.fill(pr_chi2,lv_pr.theta()  * radians_to_deg);
					h_pr_chi2_phi_final.fill(pr_chi2, lv_pr.phi() * radians_to_deg);
					double chi2_prob = getChi2Prob(pr_chi2, pr_ndf);
					h_pr_chi2_prob_final.fill( chi2_prob );
				    }

				    if( kp_chi2 > 0 ){
					h_kp_chi2_final.fill(kp_chi2);
					h_kp_chi2_p_final.fill(kp_chi2,lv_kp.p());
					h_kp_chi2_theta_final.fill(kp_chi2,lv_kp.theta()  * radians_to_deg);
					h_kp_chi2_phi_final.fill(kp_chi2, lv_kp.phi() * radians_to_deg);
					double chi2_prob = getChi2Prob(kp_chi2, kp_ndf);
					h_kp_chi2_prob_final.fill( chi2_prob );
				    }
				    if( km_chi2 > 0 ){
					h_km_chi2_final.fill(km_chi2);
					h_km_chi2_p_final.fill(km_chi2,lv_km.p());
					h_km_chi2_theta_final.fill(km_chi2,lv_km.theta()  * radians_to_deg);
					h_km_chi2_phi_final.fill(km_chi2, lv_km.phi() * radians_to_deg);
					double chi2_prob = getChi2Prob(km_chi2, km_ndf);
					h_km_chi2_prob_final.fill( chi2_prob );
				    }
			    

				    LorentzVector lv_phi = new LorentzVector(0,0,0,0);			    
				    if( low_epkpkm_mm2 && miss_perp && colin_pr && colin_km ){					
					lv_phi.add(lv_kp);
					lv_phi.add(lv_km);
					h_prelim_phi.fill(lv_phi.mass());

					// include beta vs p of final hadrons
					h_pr_final_betap.fill(lv_pr.p(),pr_beta_final);
					h_kp_final_betap.fill(lv_kp.p(),kp_beta_final);
					h_km_final_betap.fill(lv_km.p(),km_beta_final);

					v_helicity.add(helicity);
					Vector<LorentzVector> final_temp = new Vector<LorentzVector>();
					
					final_temp.add(lv_el);
					final_temp.add(lv_pr);
					final_temp.add(lv_kp);
					final_temp.add(lv_km);
					final_event.add(final_temp);

					h_el_p_theta_final.fill( lv_el.p(), lv_el.theta() * radians_to_deg );
					h_el_theta_phi_final.fill( lv_el.phi() * radians_to_deg, lv_el.theta() * radians_to_deg );
					h_pr_p_theta_final.fill( lv_pr.p(), lv_pr.theta() * radians_to_deg );
					h_pr_theta_phi_final.fill( lv_pr.phi() * radians_to_deg, lv_pr.theta() * radians_to_deg );
					h_kp_p_theta_final.fill( lv_kp.p(), lv_kp.theta() * radians_to_deg );
					h_kp_theta_phi_final.fill( lv_kp.phi() * radians_to_deg, lv_kp.theta() * radians_to_deg );
					h_km_p_theta_final.fill( lv_km.p(), lv_km.theta() * radians_to_deg );
					h_km_theta_phi_final.fill( lv_km.phi() * radians_to_deg, lv_km.theta() * radians_to_deg );
					
					h_q2w_final.fill(q2,w);
					h_q2x_final.fill(xb,q2);
					h_q2t_final.fill(t,q2);
					h_t_final.fill(t);
					h_xbt_final.fill(xb,t);

					h_ang_elpr_final.fill(ang_elpr);
					h_ang_kpkm_final.fill(ang_kpkm);
					
					h_el_epX_final.fill(lv_epX.mass2());
					h_el_epkpX_final.fill(lv_epkpX.mass2());
					h_el_epkpkmX_final.fill(lv_epkpkmX.mass2());
					h_el_epkpkmMM2_final.fill(lv_epkpkmX.mass2());
					h_el_ekpX_final.fill(lv_ekpX.mass());
					h_el_ekpkmX_final.fill(lv_ekpkmX.mass2());
				      
					h_miss_per_final.fill(Math.sqrt(lv_epkpkmX.px()*lv_epkpkmX.px() + lv_epkpkmX.py()*lv_epkpkmX.py()));
					h_colin_pr_final.fill(colin_pr_ang);
					h_colin_km_final.fill(colin_km_ang);
					h_colin_kp_final.fill(colin_kp_ang);
					
					//double cos_theta_kpkm=Math.cos( lv_kp.
					//h_kpkm_theta_final.fill(

					if( dataType == "MC"){
					    h_accep_rec_q2.fill(q2);
					    h_accep_rec_t.fill(t);
					}

					if( pr_chi2 < 20 && pr_chi2 > 0 ){
					    h_phi_final_chi2_final_c1.fill(lv_phi.mass());					
					    if( kp_chi2 < 20 && kp_chi2 > 0 ){
						h_phi_final_chi2_final_c2.fill(lv_phi.mass());					    
						if( km_chi2 < 20 && km_chi2 > 0 ){
						    h_phi_final_chi2_final_c3.fill(lv_phi.mass());
						}
					    }
					}
					
					//System.out.println(" event " + Integer.toString(num_ev-1) );
					//System.out.println(" >> phi event " +  " el momentum  " + Double.toString(lv_el.p()) );
					//System.out.println(" >>           " +  " pr momentum  " + Double.toString(lv_pr.p()) );
					//System.out.println(" >>           " +  " kp momentum  " + Double.toString(lv_kp.p()) );
					//System.out.println(" >>           " +  " km momentum  " + Double.toString(lv_km.p()) );


				    }
				}
			    }
			}
		    		
		    }
		}
		
		if( mcBank_pres ){
		    DataBank mcbank = event.getBank("MC::Particle");
		    
		    for ( int i = 0; i < mcbank.rows(); i++ ){
			

			int rec_pid = mcbank.getInt("pid",i);			    
			float rec_px = mcbank.getFloat("px",i);
			float rec_py = mcbank.getFloat("py",i);
			float rec_pz = mcbank.getFloat("pz",i);
			float rec_vz = mcbank.getFloat("vz",i);

			LorentzVector lv_el_mc = new LorentzVector(0,0,0,0);
			LorentzVector lv_pr_mc = new LorentzVector(0,0,0,0);
			LorentzVector lv_kp_mc = new LorentzVector(0,0,0,0);
			LorentzVector lv_km_mc = new LorentzVector(0,0,0,0);

			if( rec_pid == electron_pid ){
 			    LorentzVector lv_temp = new LorentzVector(0,0,0,0);
			    double rec_e = Math.sqrt( rec_px*rec_px + rec_py*rec_py + rec_pz*rec_pz + mass_electron*mass_electron );
			    lv_temp.setPxPyPzE(rec_px,rec_py,rec_pz,rec_e);
			    lv_el_mc=lv_temp;
			    lv_el_gen = lv_el_mc;
			    h_el_p_theta_mc.fill( lv_temp.p(), lv_temp.theta() * radians_to_deg );
			    h_el_theta_phi_mc.fill( lv_temp.phi() * radians_to_deg, lv_temp.theta() * radians_to_deg );
			}
			else if (rec_pid == proton_pid ){
			    LorentzVector lv_temp = new LorentzVector(0,0,0,0);
			    double rec_e = Math.sqrt( rec_px*rec_px + rec_py*rec_py + rec_pz*rec_pz + mass_proton*mass_proton );
			    lv_temp.setPxPyPzE(rec_px,rec_py,rec_pz,rec_e);
			    lv_pr_mc=lv_temp;
			    lv_pr_gen = lv_pr_mc;

			    h_pr_p_theta_mc.fill( lv_temp.p(), lv_temp.theta() * radians_to_deg );  
			    h_pr_theta_phi_mc.fill( lv_temp.phi() * radians_to_deg, lv_temp.theta() * radians_to_deg );     
			}
			else if ( rec_pid == kaon_plus_pid ){
			    LorentzVector lv_temp = new LorentzVector(0,0,0,0);
			    double rec_e = Math.sqrt( rec_px*rec_px + rec_py*rec_py + rec_pz*rec_pz + mass_kaon*mass_kaon );
			    lv_temp.setPxPyPzE(rec_px,rec_py,rec_pz,rec_e);
			    lv_kp_mc=lv_temp;
			    lv_kp_gen = lv_kp_mc;

			    h_kp_p_theta_mc.fill( lv_temp.p(), lv_temp.theta() * radians_to_deg );   
			    h_kp_theta_phi_mc.fill(lv_temp.phi() * radians_to_deg, lv_temp.theta() * radians_to_deg ); 
			}
			else if ( rec_pid == kaon_minus_pid ){
			    LorentzVector lv_temp = new LorentzVector(0,0,0,0);
			    double rec_e = Math.sqrt( rec_px*rec_px + rec_py*rec_py + rec_pz*rec_pz + mass_kaon*mass_kaon );
			    lv_temp.setPxPyPzE(rec_px,rec_py,rec_pz,rec_e);				    
			    lv_km_mc=lv_temp;
			    lv_km_gen = lv_km_mc;

			    h_km_p_theta_mc.fill( lv_temp.p(), lv_temp.theta() * radians_to_deg );    
			    h_km_theta_phi_mc.fill(lv_temp.phi() * radians_to_deg, lv_temp.theta() * radians_to_deg ); 
			}

			double mc_t = 2.0*mass_proton*(lv_pr_mc.e() - mass_proton);			
			double mc_q2 = 4*lv_beam.e()*lv_el_mc.e()*Math.sin(lv_el_mc.theta()/2.0)*Math.sin(lv_el_mc.theta()/2.0);
			double mc_xb = mc_q2 / (2.0*mass_proton*(lv_beam.e() - lv_el_mc.e()));
			double mc_w = Math.sqrt(-mc_q2 + mass_proton*mass_proton + 2*mass_proton*(lv_beam.e() - lv_el_mc.e()));
			
 			h_q2w_mc.fill(mc_q2,mc_w);
			h_q2x_mc.fill(mc_xb,mc_q2);
			h_q2t_mc.fill(mc_t,mc_q2);
			h_t_mc.fill(mc_t);
			h_xbt_mc.fill(mc_xb,mc_t);
			
			//fill histogram for acceptance, only do it once
			if ( i == 0 ){
			    h_accep_gen_q2.fill(mc_q2);
			    h_accep_gen_t.fill(mc_t);
			}



		    }
		}

		// check out resolutions here
		if( mcBank_pres && recBank_pres ){
		    
		    //rec
		    double rec_t = 2.0*mass_proton*(lv_pr_rec.e() - mass_proton);		    
		    double rec_q2 = 4*lv_beam.e()*lv_el_rec.e()*Math.sin(lv_el_rec.theta()/2.0)*Math.sin(lv_el_rec.theta()/2.0);
		    double rec_xb = rec_q2 / (2.0*mass_proton*(lv_beam.e() - lv_el_rec.e()));
		    double rec_w = Math.sqrt(-rec_q2 + mass_proton*mass_proton + 2*mass_proton*(lv_beam.e() - lv_el_rec.e()));
		    
		    //gen
		    double gen_t = 2.0*mass_proton*(lv_pr_gen.e() - mass_proton);		    
		    double gen_q2 = 4*lv_beam.e()*lv_el_gen.e()*Math.sin(lv_el_gen.theta()/2.0)*Math.sin(lv_el_gen.theta()/2.0);
		    double gen_xb = gen_q2 / (2.0*mass_proton*(lv_beam.e() - lv_el_gen.e()));
		    double gen_w = Math.sqrt(-gen_q2 + mass_proton*mass_proton + 2*mass_proton*(lv_beam.e() - lv_el_gen.e()));
		   
		    double res_q2 = resolution( rec_q2, gen_q2 );
		    double res_w = resolution( rec_w, gen_w );
		    double res_t = resolution(rec_t, gen_t );
		    double res_xb = resolution(rec_xb, gen_xb);
				    		    

			
		    if ( lv_el_gen.e() > 0 && lv_el_rec.e() > 0 ){
			h_res_q2.fill(res_q2);
			h_res_w.fill(res_w);
			h_res_xb.fill( res_xb);

			h2_res_q2.fill( gen_q2, res_q2 );
			h2_res_w.fill( gen_w, res_w );
			h2_res_xb.fill( gen_xb, res_xb);
		    }		    
		    if( lv_pr_gen.e() > 0 && lv_pr_rec.e() > 0){
			h_res_t.fill(res_t);
			h2_res_t.fill( gen_t, res_t);
		    }



		}


	                  		  	    
	    }
	    file_counter++;
	}


	F1D betap_pr = new F1D("betap_pr"," x/sqrt(x*x + (0.938*0.938))",0.10, 4.0);
	F1D betap_kp = new F1D("betap_kp"," x/sqrt(x*x + (0.493*0.493))",0.10, 4.0);
	F1D betap_km = new F1D("betap_km"," x/sqrt(x*x + (0.493*0.493))",0.10, 4.0);
	F1D betap_pip = new F1D("betap_pip"," x/sqrt(x*x + (0.139*0.139))",0.10, 4.0);
	F1D betap_pim = new F1D("betap_pim"," x/sqrt(x*x + (0.139*0.139))",0.10, 4.0);

	betap_pr.setLineColor(2);
	betap_kp.setLineColor(2);
	betap_km.setLineColor(2);
	betap_pip.setLineColor(2);
	betap_pim.setLineColor(2);

	betap_pr.setLineWidth(1);
	betap_kp.setLineWidth(1);
	betap_km.setLineWidth(1);
	betap_pip.setLineWidth(1);
	betap_pim.setLineWidth(1);

	EmbeddedCanvas c_temp = new EmbeddedCanvas();
	c_temp.setSize(600,900);
	c_temp.divide(2,3);
	int s =0;
	for( H2F h_temp :  h_pos_betap ){
	    c_temp.cd(s);
	    
	    s++;
	    c_temp.draw(h_temp);	    
	    c_temp.draw(betap_pr,"same");
	    c_temp.draw(betap_kp,"same");
	    c_temp.draw(betap_pip,"same");
	}
	c_temp.save("/work/clas12/bclary/CLAS12/phi_analysis/pid/h_pos_betap.png");

	EmbeddedCanvas c_temp2 = new EmbeddedCanvas();
	c_temp2.setSize(600,900);
	c_temp2.divide(2,3);
	s=0;
	for( H2F h_temp :  h_neg_betap ){
	    c_temp2.cd(s);
	    s++;
	    c_temp2.draw(h_temp);
	    c_temp2.draw(betap_km,"same");
	    c_temp2.draw(betap_pim,"same");
	}
	c_temp2.save("/work/clas12/bclary/CLAS12/phi_analysis/pid/h_neg_betap.png");

	EmbeddedCanvas c_part = new EmbeddedCanvas();
	c_part.setSize(1000,7500);
	c_part.divide(2,15);
	c_part.cd(0);
	c_part.draw(h_el_p_theta);
	c_part.cd(1);
	c_part.draw(h_el_theta_phi);
	c_part.cd(2);
	c_part.draw(h_pr_p_theta);
	c_part.cd(3);
	c_part.draw(h_pr_theta_phi);
	c_part.cd(4);
	c_part.draw(h_kp_p_theta);
	c_part.cd(5);
	c_part.draw(h_kp_theta_phi);
	c_part.cd(6);
	c_part.draw(h_km_p_theta);
	c_part.cd(7);
	c_part.draw(h_km_theta_phi);
	c_part.cd(8);
	c_part.draw(h_q2w);
	c_part.cd(9);
	c_part.draw(h_q2x);
	c_part.cd(10);
	c_part.draw(h_q2t);
	c_part.cd(11);
	c_part.draw(h_t);
	c_part.cd(12);
 	c_part.draw(h_el_eX);
	c_part.cd(13);
 	c_part.draw(h_el_epX);
	c_part.cd(14);
 	c_part.draw(h_el_epkpX);
	c_part.cd(15);
 	c_part.draw(h_el_epkpkmX);
	c_part.cd(16);
 	c_part.draw(h_el_ekpX);
	c_part.cd(17);
 	c_part.draw(h_el_ekpkmX);
	c_part.cd(18);
	c_part.draw(h_el_epkpkmMM2);
	c_part.cd(19);
	c_part.draw(h_miss_per);
	c_part.cd(20);
	c_part.draw(h_colin_pr);
	c_part.cd(21);
	c_part.draw(h_colin_kp);
	c_part.cd(22);
	c_part.draw(h_colin_km);
 	c_part.cd(23);
	c_part.draw(h_ang_elpr);
	c_part.cd(24);
	c_part.draw(h_ang_elpr);
	c_part.cd(25);
	c_part.draw(h_prelim_phi);
	c_part.cd(26);
	c_part.draw(h_pr_final_betap);
	c_part.cd(27);
	c_part.draw(h_kp_final_betap);
	c_part.cd(28);
	c_part.draw(h_km_final_betap);
	c_part.save("/work/clas12/bclary/CLAS12/phi_analysis/pid/all_particle_kin.png");


	EmbeddedCanvas c_part_final = new EmbeddedCanvas();
	c_part_final.setSize(1000,7500);
	c_part_final.divide(2,15);
	c_part_final.cd(0);
	c_part_final.draw(h_el_p_theta_final);
	c_part_final.cd(1);
	c_part_final.draw(h_el_theta_phi_final);
	c_part_final.cd(2);
	c_part_final.draw(h_pr_p_theta_final);
	c_part_final.cd(3);
	c_part_final.draw(h_pr_theta_phi_final);
	c_part_final.cd(4);
	c_part_final.draw(h_kp_p_theta_final);
	c_part_final.cd(5);
	c_part_final.draw(h_kp_theta_phi_final);
	c_part_final.cd(6);
	c_part_final.draw(h_km_p_theta_final);
	c_part_final.cd(7);
	c_part_final.draw(h_km_theta_phi_final);
	c_part_final.cd(8);
	c_part_final.draw(h_q2w_final);
	c_part_final.cd(9);
	c_part_final.draw(h_q2x_final);
	c_part_final.cd(10);
	c_part_final.draw(h_q2t_final);
	c_part_final.cd(11);
	c_part_final.draw(h_t_final);
	c_part_final.cd(12);
 	c_part_final.draw(h_el_eX_final);
	c_part_final.cd(13);
 	c_part_final.draw(h_el_epX_final);
	c_part_final.cd(14);
 	c_part_final.draw(h_el_epkpX_final);
	c_part_final.cd(15);
 	c_part_final.draw(h_el_epkpkmX_final);
	c_part_final.cd(16);
 	c_part_final.draw(h_el_ekpX_final);
	c_part_final.cd(17);
 	c_part_final.draw(h_el_ekpkmX_final);
	c_part_final.cd(18);
	c_part_final.draw(h_el_epkpkmMM2_final);
	c_part_final.cd(19);
	c_part_final.draw(h_miss_per_final);
	c_part_final.cd(20);
	c_part_final.draw(h_colin_pr_final);
	c_part_final.cd(21);
	c_part_final.draw(h_colin_kp_final);
	c_part_final.cd(22);
	c_part_final.draw(h_colin_km_final);
 	c_part_final.cd(23);
	c_part_final.draw(h_ang_elpr_final);
	c_part_final.cd(24);
	c_part_final.draw(h_ang_elpr_final);
	c_part_final.cd(25);
	c_part_final.draw(h_prelim_phi);
	c_part_final.cd(26);
	c_part_final.draw(h_pr_final_betap);
	c_part_final.cd(27);
	c_part_final.draw(h_kp_final_betap);
	c_part_final.cd(28);
	c_part_final.draw(h_km_final_betap);
	c_part_final.save("/work/clas12/bclary/CLAS12/phi_analysis/pid/all_particle_kin_final.png");


	EmbeddedCanvas c_part_chi2 = new EmbeddedCanvas();
	c_part_chi2.setSize(1800,4200);
	c_part_chi2.divide(3,7);
	c_part_chi2.cd(0);
	c_part_chi2.draw(h_pr_chi2_before);
	c_part_chi2.cd(1);
	c_part_chi2.draw(h_kp_chi2_before);
	c_part_chi2.cd(2);
	c_part_chi2.draw(h_km_chi2_before);
 	c_part_chi2.cd(3);
	c_part_chi2.draw(h_pr_chi2_p_before);
	c_part_chi2.cd(4);
	c_part_chi2.draw(h_pr_chi2_theta_before);
	c_part_chi2.cd(5);
	c_part_chi2.draw(h_pr_chi2_phi_before);
 	c_part_chi2.cd(6);
	c_part_chi2.draw(h_kp_chi2_p_before);
	c_part_chi2.cd(7);
	c_part_chi2.draw(h_kp_chi2_theta_before);
	c_part_chi2.cd(8);
	c_part_chi2.draw(h_kp_chi2_phi_before);
 	c_part_chi2.cd(9);
	c_part_chi2.draw(h_km_chi2_p_before);
	c_part_chi2.cd(10);
	c_part_chi2.draw(h_km_chi2_theta_before);
	c_part_chi2.cd(11);
	c_part_chi2.draw(h_km_chi2_phi_before);
 	c_part_chi2.cd(12);
	c_part_chi2.draw(h_phi_final_chi2_final_c1);
	c_part_chi2.cd(13);
	c_part_chi2.draw(h_phi_final_chi2_final_c2);
	c_part_chi2.cd(14);
	c_part_chi2.draw(h_phi_final_chi2_final_c3);
	c_part_chi2.cd(15);
	c_part_chi2.getPad(15).getAxisY().setLog(true);
	c_part_chi2.draw(h_pr_chi2_prob_before);
	c_part_chi2.cd(16);
	c_part_chi2.getPad(16).getAxisY().setLog(true);
	c_part_chi2.draw(h_kp_chi2_prob_before);
	c_part_chi2.cd(17);
	c_part_chi2.getPad(17).getAxisY().setLog(true);
	c_part_chi2.draw(h_km_chi2_prob_before);
	c_part_chi2.cd(18);
	c_part_chi2.getPad(18).getAxisY().setLog(true);
	c_part_chi2.draw(h_pr_chi2_prob_final);
	c_part_chi2.cd(19);
	c_part_chi2.getPad(19).getAxisY().setLog(true);
	c_part_chi2.draw(h_kp_chi2_prob_final);
	c_part_chi2.cd(20);
	c_part_chi2.getPad(20).getAxisY().setLog(true);
	c_part_chi2.draw(h_km_chi2_prob_final);              
	c_part_chi2.save("/work/clas12/bclary/CLAS12/phi_analysis/pid/all_hadrons_chi2.png");

	EmbeddedCanvas c_pr_chi2 = new EmbeddedCanvas();
	c_pr_chi2.setSize(1800,4200);
	c_pr_chi2.divide(3,7);
	for( int i = 0; i < v_pr_chi2.size(); i++ ){	    
	    c_pr_chi2.cd(i);	    
	    c_pr_chi2.draw(v_pr_chi2.get(i));
	}
	c_pr_chi2.save("/work/clas12/bclary/CLAS12/phi_analysis/pid/pr_chi2_ndf.png");

	EmbeddedCanvas c_kp_chi2 = new EmbeddedCanvas();
	c_kp_chi2.setSize(1800,4200);
	c_kp_chi2.divide(3,7);
	for( int i = 0; i < v_kp_chi2.size(); i++ ){	    
	    c_kp_chi2.cd(i);	    
	    c_kp_chi2.draw(v_kp_chi2.get(i));
	}
	c_kp_chi2.save("/work/clas12/bclary/CLAS12/phi_analysis/pid/kp_chi2_ndf.png");

	EmbeddedCanvas c_km_chi2 = new EmbeddedCanvas();
	c_km_chi2.setSize(1800,4200);
	c_km_chi2.divide(3,7);
	for( int i = 0; i < v_km_chi2.size(); i++ ){	    
	    c_km_chi2.cd(i);	    
	    c_km_chi2.draw(v_km_chi2.get(i));
	}
	c_km_chi2.save("/work/clas12/bclary/CLAS12/phi_analysis/pid/km_chi2_ndf.png");

	EmbeddedCanvas c_part_mc = new EmbeddedCanvas();
	c_part_mc.setSize(1000,7500);
	c_part_mc.divide(2,15);
	c_part_mc.cd(0);
	c_part_mc.draw(h_el_p_theta_mc);
	c_part_mc.cd(1);
	c_part_mc.draw(h_el_theta_phi_mc);
	c_part_mc.cd(2);
	c_part_mc.draw(h_pr_p_theta_mc);
	c_part_mc.cd(3);
	c_part_mc.draw(h_pr_theta_phi_mc);
	c_part_mc.cd(4);
	c_part_mc.draw(h_kp_p_theta_mc);
	c_part_mc.cd(5);
	c_part_mc.draw(h_kp_theta_phi_mc);
	c_part_mc.cd(6);
	c_part_mc.draw(h_km_p_theta_mc);
	c_part_mc.cd(7);
	c_part_mc.draw(h_km_theta_phi_mc);
	c_part_mc.cd(8);
	c_part_mc.draw(h_q2w_mc);
	c_part_mc.cd(9);
	c_part_mc.draw(h_q2x_mc);
	c_part_mc.cd(10);
	c_part_mc.draw(h_q2t_mc);
	c_part_mc.cd(11);
	c_part_mc.draw(h_xbt_mc);	
	c_part_mc.cd(12);
	c_part_mc.draw(h_t_mc);
	c_part_mc.save("/work/clas12/bclary/CLAS12/phi_analysis/pid/all_particle_kin_mc.png");

	////////////
	//resolution and acceptance histograms 
	if (dataType == "MC" ){
	    for( int bin = 0; bin <  h_accep_rec_q2.getXaxis().getNBins(); bin++ ){		
		double gen_value =  h_accep_gen_q2.getBinContent(bin);
		double rec_value =  h_accep_rec_q2.getBinContent(bin);		
		double accp = rec_value/gen_value;
		System.out.println(" >> GEN " + gen_value + " REC " + rec_value + " bin " + bin + " accp for q2 " + accp );
		if ( gen_value == 0 ) continue;
		h_accp_q2.setBinContent(bin, accp );
	    }

	    for( int bin = 0; bin <  h_accep_rec_t.getXaxis().getNBins(); bin++ ){		
		double gen_value =  h_accep_gen_t.getBinContent(bin);
		double rec_value =  h_accep_rec_t.getBinContent(bin);		
		double accp = rec_value/gen_value;
		System.out.println(" >> GEN " + gen_value + " REC " + rec_value + " bin " + bin + " accp for t " + accp );
		if ( gen_value == 0 ) continue;
		h_accp_t.setBinContent(bin, accp );
	    }

	    EmbeddedCanvas c_res = new EmbeddedCanvas();
	    c_res.setSize(1000,7500);
	    c_res.divide(2,15);
	    c_res.cd(0);
	    c_res.draw(h_res_q2);
	    c_res.cd(1);
	    c_res.draw(h_res_w);
	    c_res.cd(2);
	    c_res.draw(h_res_t);
	    c_res.cd(3);
	    c_res.draw(h_res_xb);
	    c_res.cd(4);
	    c_res.draw(h2_res_q2);
	    c_res.cd(5);
	    c_res.draw(h2_res_w);
	    c_res.cd(6);
	    c_res.draw(h2_res_t);
	    c_res.cd(7);
	    c_res.draw(h2_res_xb);
	    c_res.cd(8);
	    c_res.draw( h_accep_gen_q2 );
	    c_res.cd(9);
	    c_res.draw( h_accep_rec_q2 );
	    c_res.cd(10);
	    c_res.draw(h_accp_q2);
	    c_res.cd(11);
	    c_res.draw( h_accep_gen_t );
	    c_res.cd(12);
	    c_res.draw( h_accep_rec_t );
	    c_res.cd(13);
	    c_res.draw(h_accp_t);
	    c_res.save("/work/clas12/bclary/CLAS12/phi_analysis/pid/particles_resolution.png");	    	  
	}
	
	       
	///////////////////////////////////////////////////////////////////////////////////////////
	// write the lorentz vectors of the ep kp km out to txt file
	//CREATE TXT OUTPUT FILE
	BufferedWriter writertxt = null;
	String output_name = "all_final_state_particles.txt";
	try{
	    writertxt = new BufferedWriter( new FileWriter(output_name) );
	}
	catch ( IOException e ){
	    System.out.println(">> ERROR " );
	}

	int event=0;
	for( Vector<LorentzVector> final_particles : final_event ){
	    try{
	    writertxt.write(Integer.toString(v_helicity.get(event)) + " ");
	    }
	    catch(IOException e){
		System.out.println(" error " );
	    }
	    for( LorentzVector lv : final_particles ){
		try{
		    writertxt.write(Double.toString(lv.px()) + " " + Double.toString(lv.py()) + " " + Double.toString(lv.pz()) + " ");
		}
		catch( IOException e){
		    System.out.println(">> ERROR WITH WRITING KINEMATICS TO FILE");
		}	      		
	    }	    
	    try{
		writertxt.write("\n");    
	    }
	    catch( IOException e){	   
		System.out.println(">> ERROR WITH NEW LINE");
	    }
	    event++;
	}
	try{
	    writertxt.close();
	}
	catch( IOException e ){
	    System.out.println(">> ERROR CLOSING FILE");  
	}
	
    }


    public static double resolution( double rec_val, double gen_val ){
	return rec_val - gen_val;
    }


    public static double getChi2Prob( double chi2, int dof ){

	/*double prob = 0.0;
	if( chi2 > 0 ){
	    double lambda = dof/2.0;
	    double norm = gamma(lambda) * Math.pow(2, lambda);
	    prob = Math.pow(chi2,lambda-1) * Math.exp(-0.5*chi2)/norm;			    
	}
	else{
	    prob = 0;
	}
	*/
	double prob = pochisq(chi2,dof);
	System.out.println(" PROB FOR " + Double.toString(chi2)  + "  and DOF  " + Integer.toString(dof)  + " IS  " + Double.toString(prob) );
	return prob;
    }

    /*************************************************************
    Reads in a command line input x and prints Gamma(x) and
    *  log Gamma(x). The Gamma function is defined by
    *  
    *        Gamma(x) = integral( t^(x-1) e^(-t), t = 0 .. infinity)
    *
    *  Uses Lanczos approximation formula. See Numerical Recipes 6.1.
    * from R. Sedgewick and K. Wayne
    *************************************************************/
    static double logGamma(double x) {
	double tmp = (x - 0.5) * Math.log(x + 4.5) - (x + 4.5);
	double ser = 1.0 + 76.18009173    / (x + 0)   - 86.50532033    / (x + 1)
	    + 24.01409822    / (x + 2)   -  1.231739516   / (x + 3)
	    +  0.00120858003 / (x + 4)   -  0.00000536382 / (x + 5);
	return tmp + Math.log(ser * Math.sqrt(2 * Math.PI));
    }
    static double gamma(double x) { return Math.exp(logGamma(x)); }
 
    public static double pochisq(double x, int df) {
        double a, s;
        double e, c, z;
 
        if (x <= 0.0 || df < 1) {
            return 1.0;
        }
        a = 0.5 * x;
        boolean even = (df & 1) == 0;
        double y = 0;
        if (df > 1) {
            y = ex(-a);
        }
        s = (even ? y : (2.0 * poz(-Math.sqrt(x))));
        if (df > 2) {
            x = 0.5 * (df - 1.0);
            z = (even ? 1.0 : 0.5);
            if (a > MAX_X) {
                e = (even ? 0.0 : LOG_SQRT_PI);
                c = Math.log(a);
                while (z <= x) {
                    e = Math.log(z) + e;
                    s += ex(c * z - a - e);
                    z += 1.0;
                }
                return s;
            } else {
                e = (even ? 1.0 : (I_SQRT_PI / Math.sqrt(a)));
                c = 0.0;
                while (z <= x) {
                    e = e * (a / z);
                    c = c + e;
                    z += 1.0;
                }
                return c * y + s;
            }
        } else {
            return s;
        }
    }
 
 
    public static double poz(double z) {
        double y, x, w;
        double Z_MAX = 6.0; // Maximum meaningful z value 
        if (z == 0.0) {
            x = 0.0;
        } else {
            y = 0.5 * Math.abs(z);
            if (y >= (Z_MAX * 0.5)) {
                x = 1.0;
            } else if (y < 1.0) {
                w = y * y;
                x = ((((((((0.000124818987 * w
			    - 0.001075204047) * w + 0.005198775019) * w
			  - 0.019198292004) * w + 0.059054035642) * w
                        - 0.151968751364) * w + 0.319152932694) * w
		      - 0.531923007300) * w + 0.797884560593) * y * 2.0;
            } else {
                y -= 2.0;
                x = (((((((((((((-0.000045255659 * y
				 + 0.000152529290) * y - 0.000019538132) * y
			       - 0.000676904986) * y + 0.001390604284) * y
			     - 0.000794620820) * y - 0.002034254874) * y
			   + 0.006549791214) * y - 0.010557625006) * y
			 + 0.011630447319) * y - 0.009279453341) * y
		       + 0.005353579108) * y - 0.002141268741) * y
		     + 0.000535310849) * y + 0.999936657524;
            }
        }
        return z > 0.0 ? ((x + 1.0) * 0.5) : ((1.0 - x) * 0.5);
    }
 
 
    public static double ex(double x) {
        return (x < -MAX_X) ? 0.0 : Math.exp(x);
    }

}
