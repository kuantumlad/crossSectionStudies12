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

public class accpTest{

    private static final double PI = 3.14592658;
    private static final double LOG_SQRT_PI = Math.log(Math.sqrt(PI));
    private static final double I_SQRT_PI = 1 / Math.sqrt(PI);
    public static final int MAX_X = 20; // max value to represent exp(x)
 
    public static void main(String[] input){

	String dataType="MC";
	String dir_in = input[0];
	String s_run = input[1];
	int f_min = Integer.valueOf(input[2]);
	int f_max = Integer.valueOf(input[3]);
	int file_counter=f_min;
	int file_counter_limit = f_max;

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

	LorentzVector lv_beam = new LorentzVector(0,0,2.221,2.221);
	LorentzVector lv_target = new LorentzVector(0,0,0,mass_proton);
	
	// System.out.println(" >> test chi2 "  + Double.toString(test_chi2) + " dof " + Integer.toString(test_dof) + " test prob " + test_prob );	       
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Histograms for Kinematics
	
	H2F h_el_p_theta = new H2F("h_el_p_theta","h_el_p_theta",200, 0.0, 11.0, 200, 0.0, 60.0 );
	H2F h_el_theta_phi = new H2F("h_el_p_theta_phi","h_el_p_theta_phi",200,-180.0, 180.0, 200, 0.0, 60.0 );

	H2F h_pr_p_theta = new H2F("h_pr_p_theta","h_pr_p_theta",200, 0.0, 6.0, 200, 0.0, 80.0 );
	H2F h_pr_theta_phi = new H2F("h_pr_p_theta_phi","h_pr_p_theta_phi",200, -180.0, 180.0, 200, 0.0, 60.0 );

	// Histogram for the final states
	H2F h_el_p_theta_final = new H2F("h_el_p_theta_final","h_el_p_theta_final",200, 0.0, 11.0, 200, 0.0, 60.0 );
	H2F h_el_theta_phi_final = new H2F("h_el_p_theta_phi_final","h_el_p_theta_phi_final",200,-180.0, 180.0, 200, 0.0, 60.0 );

	H2F h_pr_p_theta_final = new H2F("h_pr_p_theta_final","h_pr_p_theta_final",200, 0.0, 6.0, 200, 0.0, 40.0 );
	H2F h_pr_theta_phi_final = new H2F("h_pr_p_theta_phi_final","h_pr_p_theta_phi_final",200, -180.0, 180.0, 200, 0.0, 60.0 );

	H2F h_pr_final_betap = new H2F("h_pr_final_betap","h_pr_final_betap", 100, 0.0, 4.0, 100, 0.01, 1.1);

	H1F h_pr_chi2_before = new H1F("h_pr_chi2_before","h_pr_chi2_before",500, -1.0, 15.0);

	H1F h_pr_chi2_prob_final = new H1F("h_pr_chi2_prob_final","h_pr_chi2_prob_final", 200, 0.0, 1.1);

	H2F h_el_p_theta_mc = new H2F("h_el_p_theta_mc","h_el_p_theta_mc",200, 0.0, 11.0, 200, 0.0, 60.0 );
	H2F h_el_theta_phi_mc = new H2F("h_el_p_theta_phi_mc","h_el_p_theta_phi_mc",200,-180.0, 180.0, 200, 0.0, 60.0 );

	H2F h_el_w_theta = new H2F("h_el_w_theta","h_el_w_theta",100, 0.0, 2.2, 100, 0.0, 60.0);
	H2F h_el_gen_w_theta = new H2F("h_el_gen_w_theta","h_el_gen_w_theta",100, 0.0, 2.2, 100, 0.0, 60.0);


	H2F h_q2w = new H2F("h_q2w","h_q2w",200,0.0, 9.0, 200, 0.0, 6.0);
	H2F h_q2x = new H2F("h_q2x","h_q2x",200,0.0, 1.0, 200, 0.0, 9.0);

	H1F h_el_eX = new H1F("h_el_eX","h_el_eX",200,-1.0, 6.0);
	H1F h_res_q2 = new H1F("h_res_q2","h_res_q2",100, -0.50, 0.50 );
	H1F h_res_w = new H1F("h_res_w","h_res_w",100, -0.50, 0.50 );

	/////////// 
	// Histogram for acceptance
	int n_bins = 30;
	double min_theta = 0;
	double max_theta = 30.0;
 	H1F h_accep_gen_theta = new H1F("h_accep_gen_theta","h_accep_gen_theta", n_bins, min_theta, max_theta );
  	H1F h_accep_rec_theta = new H1F("h_accep_rec_theta","h_accep_rec_theta", n_bins, min_theta, max_theta );
	H1F h_accp_q2 = new H1F("h_accp_theta","h_accp_theta",n_bins,min_theta, max_theta);

	//get acceptance per sector
	Vector<H1F> h_gen_sect = new Vector<H1F>();
	Vector<H1F> h_rec_sect = new Vector<H1F>();
	Vector<H1F> h_accp_sect = new Vector<H1F>();
	for( int s = 1; s <= 6; s++ ){
	    h_gen_sect.add( new H1F("h_gen_sect"+Integer.toString(s),"h_gen_sect"+Integer.toString(s),n_bins,min_theta,max_theta) );
	    h_rec_sect.add( new H1F("h_rec_sect"+Integer.toString(s),"h_rec_sect"+Integer.toString(s),n_bins,min_theta,max_theta) );
	    h_accp_sect.add( new H1F("h_accp_sect"+Integer.toString(s),"h_accp_sect"+Integer.toString(s),n_bins,min_theta,max_theta) );
	}


	H2F h_gen_phitheta = new H2F("h_gen_phitheta","h_gen_phitheta",74,-180,180,30, 0, 30);
	H2F h_rec_phitheta = new H2F("h_gen_phitheta","h_gen_phitheta",74,-180,180,30, 0, 30);

 	H1F h_gen_theta = new H1F("h_gen_theta","h_gen_theta", n_bins, min_theta, max_theta );
       
	Vector<H1F> h_accp_theta = new Vector<H1F>();
	Vector<H1F> h_rec_theta_s = new Vector<H1F>();
	Vector<H1F> h_gen_theta_s = new Vector<H1F>();
	Vector<H2F>  h_el_w_theta_s = new Vector<H2F>();

	for( int s = 0; s < 6; s++ ){
	    h_rec_theta_s.add( new H1F("h_rec_theta_sector"+Integer.toString(s),"h_rec_theta_sector"+Integer.toString(s), n_bins, min_theta, max_theta) );
	    h_gen_theta_s.add( new H1F("h_gen_theta_sector"+Integer.toString(s),"h_gen_theta_sector"+Integer.toString(s), n_bins, min_theta, max_theta) );
	    h_el_w_theta_s.add( new H2F("h_el_w_theta_sector"+Integer.toString(s),"h_el_w_theta_sector"+Integer.toString(s), 100, 0.0, 2.2, 100, 0.0, 60.0) );
	}

	Vector<Integer> run_list = new Vector<Integer>();
	run_list.add(0);

	while( file_counter < file_counter_limit ){


	    System.out.println("Opening File ");
	    
	    String file_in=null;
	    if (dataType != "MC" ){
		file_in = dir_in + Integer.toString(run_list.get(file_counter)) + ".hipo";
	    }
	    else{
		file_in = dir_in + Integer.toString(file_counter) + ".lund.hipo";
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

	    int max_event = hiporeader.getSize();
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

		int el_cal_sector = -1;		

		if (event.hasBank("REC::Particle")  == true && event.hasBank("REC::Calorimeter") == true&&event.hasBank("REC::Event")){     
		    for( int rec_i = 0; rec_i < event.getBank("REC::Particle").rows(); rec_i++ ){
			DataBank recbank = event.getBank("REC::Particle");		
			DataBank calbank = event.getBank("REC::Calorimeter");
			
			
			int charge = recbank.getInt("charge",rec_i);
			int cal_sector = -1;
			    
			if( charge != 0 ){
			    int rec_pid = recbank.getInt("pid",rec_i);			    
			    float rec_px = recbank.getFloat("px",rec_i);
			    float rec_py = recbank.getFloat("py",rec_i);
			    float rec_pz = recbank.getFloat("pz",rec_i);
			    float rec_vz = recbank.getFloat("vz",rec_i);
			    
			    int rec_status = recbank.getInt("status",rec_i);
			    float beta_clas12 = recbank.getFloat("beta",rec_i);
			    //if( 2000 >= rec_status && rec_status >= 3000 ) continue;
			    for( int i = 0; i < calbank.rows(); i++){			   
				int pindex = calbank.getShort("pindex",i);
				if( pindex == rec_i ){					
				    cal_sector = calbank.getInt("sector",i) - 1;
				    break;
				}
			    }
			    			    
			    Vector<LorentzVector> v_el = new Vector<LorentzVector>();			    
			    LorentzVector lv_el = new LorentzVector(0,0,0,0);
			    
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
					el_cal_sector=cal_sector;
					//System.out.println(" >> electron index " + Integer.toString(rec_i) +  " momentum  " + Double.toString(lv_el.p()) );
				    }
				}				    				
			    }	
			}
		    }

		    //System.out.println( " GOOD EL energy " + Double.toString(lv_el_rec.e()) + " CAL SECTOR " + Integer.toString(el_cal_sector));
		    if ( el_cal_sector >= 0){
			LorentzVector eX = new LorentzVector(0,0,0,0);
			eX.add(lv_beam);
			eX.add(lv_target);
			eX.sub(lv_el_rec);
			double w = eX.mass();

			h_el_w_theta_s.get(el_cal_sector).fill(w, lv_el_rec.theta() * radians_to_deg);
			h_el_w_theta.fill(w, lv_el_rec.theta() * radians_to_deg);
			if ( w < 1.1 ){
			    h_rec_theta_s.get(el_cal_sector).fill(lv_el_rec.theta() * radians_to_deg);
			    h_rec_phitheta.fill(lv_el_rec.phi() * radians_to_deg, lv_el_rec.theta() * radians_to_deg);

			    
			    if( lv_el_rec.phi() * radians_to_deg > -80.0 && lv_el_rec.phi() * radians_to_deg < -20.0 ){
				h_rec_sect.get(5).fill( lv_el_rec.theta() * radians_to_deg );
			    }
			}
		    }
		}

		////////////////////////////////////////////////////////////////////////////
		// get gen bank information
		
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
			    //System.out.println(" filling generated el energy " + Double.toString(lv_temp.e() ) + " theta " +lv_temp.theta() * radians_to_deg);

			    LorentzVector eX = new LorentzVector(0,0,0,0);
			    eX.add(lv_beam);
			    eX.add(lv_target);
			    eX.sub(lv_el_gen);
			    double w = eX.mass();

			    h_el_gen_w_theta.fill(w, lv_el_gen.theta() * radians_to_deg);
			    
			    if ( w < 1.1 ){ //w > 0.9 && w < 1.1 ){
				h_el_p_theta_mc.fill( lv_temp.p(), lv_temp.theta() * radians_to_deg );
				h_el_theta_phi_mc.fill( lv_temp.phi() * radians_to_deg, lv_temp.theta() * radians_to_deg );
				h_gen_theta.fill(lv_temp.theta() * radians_to_deg );
				h_gen_phitheta.fill(lv_temp.phi() * radians_to_deg, lv_temp.theta() * radians_to_deg);

				if( lv_el_gen.phi() * radians_to_deg > -80.0 && lv_el_gen.phi() * radians_to_deg < -20.0 ){
				    h_gen_sect.get(5).fill( lv_el_gen.theta() * radians_to_deg );
				}


			    }
			}
		    }
		}
	    }
	    file_counter++;
	}

	EmbeddedCanvas c_temp = new EmbeddedCanvas();
	c_temp.setSize(600,900);
	c_temp.divide(2,3);
	for( int s = 0; s < 6; s++ ){
	    c_temp.cd(s);	    
	    c_temp.draw(h_rec_theta_s.get(s));
	}
	c_temp.save("h_rec_theta_sector.png");

	EmbeddedCanvas c_temp1 = new EmbeddedCanvas();
	c_temp1.setSize(600,900);
	c_temp1.draw(h_rec_phitheta,"colz");
	c_temp1.save("h_rec_phitheta.png");

	EmbeddedCanvas c_temp2 = new EmbeddedCanvas();
	c_temp2.setSize(600,900);
	c_temp2.draw(h_gen_phitheta,"colz");
	c_temp2.save("h_gen_phitheta.png");

	EmbeddedCanvas c_temp3 = new EmbeddedCanvas();
	c_temp3.setSize(600,900);
	c_temp3.draw(h_gen_theta);
	c_temp3.save("h_gen_theta.png");

	////////////////////////////////////////////////////
	// acceptance stuff

	// check raw counts


	//getBinContent(int bx, int by) {
	// xAxis.getNBins();
	// get theta bins
	 int n_rec_phi_bins = h_rec_phitheta.getXAxis().getNBins();
	 int n_rec_theta_bins = h_rec_phitheta.getYAxis().getNBins();

	//get gen bins
	int n_gen_phi_bins = h_gen_phitheta.getXAxis().getNBins();
	int n_gen_theta_bins = h_gen_phitheta.getYAxis().getNBins();


	Vector<GraphErrors> g_theta_accp_per_phi = new Vector< GraphErrors>();
	Vector<H1F> h_raw_rec_events_per_phi = new Vector<H1F>();
	Vector<H1F> h_raw_gen_events_per_phi = new Vector<H1F>();

	for( int pb = 0; pb < n_rec_phi_bins; pb++ ){	   
	    
	    //g_theta_accp_per_phi.add( new GraphErrors( ) );

	    double[] theta_bin_center;
	    theta_bin_center = new double[n_rec_theta_bins];   
	    
	    double[] theta_bin_center_err;
	    theta_bin_center_err = new double[n_rec_theta_bins];  
	    
	    double[] accp;
	    accp  =  new double[n_rec_theta_bins]; 
	    double[] accp_err ;
	    accp_err =  new double[n_rec_theta_bins];   
	    
	    h_raw_rec_events_per_phi.add( new H1F("h_raw_rec_events_phi_"+Integer.toString(pb),"h_raw_rec_events_phi_"+Integer.toString(pb),40, 0.0, 40.0) );
	    h_raw_gen_events_per_phi.add( new H1F("h_raw_gen_events_phi_"+Integer.toString(pb),"h_raw_gen_events_phi_"+Integer.toString(pb),40, 0.0, 40.0) );

	    for(int tb = 0; tb < n_rec_theta_bins; tb++ ){

		//System.out.println(" phi bin " + Integer.toString(pb) + " theta bin " + Integer.toString(tb) );
		double theta_center = h_rec_phitheta.getDataY(tb);
		double rec_bin_content = h_rec_phitheta.getBinContent(pb,tb);
		double gen_bin_content = h_gen_phitheta.getBinContent(pb,tb);
		double ratio_bin_content = (double)rec_bin_content/(double)gen_bin_content;
		
		theta_bin_center[tb] = theta_center;
		theta_bin_center_err[tb] = 0.0;

		
		accp[tb]=ratio_bin_content;
		if ( gen_bin_content == 0 || gen_bin_content < 1 ) accp[tb]=0.0;

		accp_err[tb] = 0.0;

		//System.out.println(" rec bin content " + Double.toString(rec_bin_content) + " gen bin content"  + Double.toString(gen_bin_content) + " ratio " + Double.toString(ratio_bin_content) );
		if( tb > 9 ){
		    h_raw_rec_events_per_phi.get(pb).setBinContent( tb, rec_bin_content);
		    h_raw_gen_events_per_phi.get(pb).setBinContent( tb, gen_bin_content);
		}


	    }

	    for(int i = 0; i < 18; i++){
		//System.out.println(Integer.toString(i) + " " + Double.toString(theta_bin_center[i]) + " " + Double.toString(accp[i]));
	    }
	    g_theta_accp_per_phi.add( new GraphErrors("accp_theta_phibin_"+Integer.toString(pb), theta_bin_center,  accp, theta_bin_center_err, accp_err) );	    
	}

	EmbeddedCanvas c_temp4 = new EmbeddedCanvas();                                                                                                                                                   
        c_temp4.setSize(10000,9000); 
	c_temp4.divide(10,10);
	for( int pb = 0 ; pb < n_rec_phi_bins; pb++ ){  
	    c_temp4.cd(pb);
	    c_temp4.draw(g_theta_accp_per_phi.get(pb));
	}
	c_temp4.save("g_theta_accp_per_phi_bin.png");

	EmbeddedCanvas c_temp4a = new EmbeddedCanvas();
	c_temp4a.setSize(10000,9000);
	c_temp4a.divide(10,10);     
	for( int pb = 0 ; pb < n_rec_phi_bins; pb++ ){  
	    c_temp4a.cd(pb);
	    h_raw_rec_events_per_phi.get(pb).setLineColor(2);
	    c_temp4a.draw(h_raw_rec_events_per_phi.get(pb));
	    c_temp4a.draw(h_raw_gen_events_per_phi.get(pb),"same");
	    

	}
	c_temp4a.save("h_theta_comp_per_phi_bin.png");
	
		    
	EmbeddedCanvas c_temp5 = new EmbeddedCanvas();
	c_temp5.setSize(1000,1000);
	c_temp5.divide(2,3);
	for( int s = 0; s < 6; s++ ){
	    c_temp5.cd(s);
	    c_temp5.draw(h_el_w_theta_s.get(s));
	}
	c_temp5.save("h_el_w_theta_sector.png");

	EmbeddedCanvas c_temp6 = new EmbeddedCanvas();
	c_temp6.setSize(1000,1000);
	c_temp6.draw(h_el_w_theta);	
	c_temp6.save("h_el_w_theta.png");
    	     
	EmbeddedCanvas c_temp7 = new EmbeddedCanvas();
	c_temp7.setSize(1000,1000);
	c_temp7.draw(h_el_gen_w_theta);	
	c_temp7.save("h_el_gen_w_theta.png");

	//test acceptance for one sector
	EmbeddedCanvas c_temp8 = new EmbeddedCanvas();
	c_temp8.setSize(900,900);
	for( int b = 1; b < n_bins; b++ ){
	    double rec = h_rec_sect.get(5).getBinContent(b);
	    double gen = h_gen_sect.get(5).getBinContent(b);
	    double accp_val = rec/gen;
	    h_accp_sect.get(5).setBinContent(b,accp_val);
	}
	c_temp8.draw(h_accp_sect.get(5));
	c_temp8.save("h_test_accp_sector.png");
	    
	    

    }
}
	
