import java.io.*;
import java.util.*;
import org.jlab.io.hipo.*;
import org.jlab.io.base.DataEvent;
import org.jlab.io.base.DataBank;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.groot.data.*;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.clas.physics.LorentzVector;


public class genCheck{

    public static void main(String[] input){
	

	double el_mass = 0.000511;
	double pr_mass = 0.938;
	double toDeg = 180.0/3.1415926538;

	String file_in = input[0];
	HipoDataSource hiporeader = new HipoDataSource();
	hiporeader.open(new File(file_in) );
	DataEvent event = null;

	double beam = 7.546;
	LorentzVector lv_target = new LorentzVector(0,0,0,0.938);
	LorentzVector lv_ebeam = new LorentzVector(0,0,beam,beam);

	int max_event=hiporeader.getSize();
	
	int num_ev = 0;

	H1F h_el_p = new H1F("h_el_p","h_el_p",200,0.0,beam+2.0);
	H1F h_el_theta = new H1F("h_el_theta","h_el_theta",200,0.0,30.0);
	H1F h_el_phi = new H1F("h_el_phi","h_el_phi",300,-200.0,200.0);
	H2F h_el_thetap = new H2F("h_el_thetap","h_el_thetap",100,0.0,beam+2, 100, 0.0, 30.0);
	H2F h_el_thetaphi = new H2F("h_el_thetaphi","h_el_thetaphi",200,-180.0,180.0, 200, 0.0, 30.0);
	H1F h_el_w = new H1F("h_el_w","h_el_w", 200, 0.0, 2.0);
	H2F h_el_w_q2 = new H2F("h_el_w_q2","h_el_w_q2",200,0.0, beam-2, 200, 0.0, beam-2);	
	
	H1F h_pr_p = new H1F("h_pr_p","h_pr_p",200,0.0,2.0);
	H1F h_pr_theta = new H1F("h_pr_theta","h_pr_theta",200,0.0,80.0);
	H1F h_pr_phi = new H1F("h_pr_phi","h_pr_phi",300,-180.0,180.0);
	H2F h_pr_thetap = new H2F("h_pr_thetap","h_pr_thetap",200,0.0,3.0,200, 0.0, 80.0);


	while( num_ev < max_event ){
	    event = (DataEvent)hiporeader.gotoEvent(num_ev);
	    boolean runConfig_pres = event.hasBank("RUN::config");
	    boolean mcBank_pres  = event.hasBank("REC::Particle");
	    boolean rawScalerBank_pres = event.hasBank("RAW::scaler");
	    num_ev++;

	    if( mcBank_pres ){
		DataBank mcBank = event.getBank("REC::Particle");
		for( int i = 0; i < mcBank.rows(); i++ ){
		    int pid = mcBank.getInt("pid",i);

		    double px = mcBank.getFloat("px",i); 
		    double py = mcBank.getFloat("py",i); 
		    double pz = mcBank.getFloat("pz",i); 
       	    
		    if( pid == 11 ){
			double e = Math.sqrt( px*px + py*py + pz*pz + el_mass*el_mass );
			LorentzVector lv_el = new LorentzVector(px,py,pz,e);
			h_el_p.fill( lv_el.p() );
			h_el_theta.fill( lv_el.theta() * toDeg );
			h_el_phi.fill( lv_el.phi() * toDeg );			
			h_el_thetap.fill( lv_el.p(), lv_el.theta() * toDeg );
			h_el_thetaphi.fill(lv_el.phi() * toDeg,lv_el.theta() * toDeg);

			LorentzVector lv_eX = new LorentzVector(0,0,0,0);
			lv_eX.add(lv_ebeam);
			lv_eX.add(lv_target);
			lv_eX.sub(lv_el);
						
			double q2 = 4.0*beam*lv_el.e() * Math.sin( lv_el.theta()/2.0) * Math.sin( lv_el.theta()/2.0);		       
			double w = lv_eX.mass();

			if ( q2 > 1.0 ) { //&& w > 1.0 ){
			    h_el_w.fill(w);
			    h_el_w_q2.fill(w,q2);			
			}
		    }
		    if( pid == 2212 ){
			double e = Math.sqrt( px*px + py*py + pz*pz + pr_mass*pr_mass );
			LorentzVector lv_pr = new LorentzVector(px,py,pz,e);
			h_pr_p.fill( lv_pr.p() );
			h_pr_theta.fill( lv_pr.theta() * toDeg );
			h_pr_phi.fill( lv_pr.phi() * toDeg );			
			h_pr_thetap.fill( lv_pr.p(), lv_pr.theta() * toDeg );
		    }		    
		}
	    }       	    
	}
    
	System.out.println(" creating plots " );
	EmbeddedCanvas c_el = new EmbeddedCanvas();
	c_el.divide(2,2);
	c_el.setSize(800,800);
	c_el.cd(0);
	h_el_p.setTitle("Momentum");
	c_el.draw(h_el_p);
	c_el.cd(1);
	h_el_theta.setTitle("Theta");
	c_el.draw(h_el_theta);
	c_el.cd(2);
	h_el_phi.setTitle("Phi");
	c_el.draw(h_el_phi);
	c_el.save("h_gen_el_elastic.png");

	EmbeddedCanvas c_ela = new EmbeddedCanvas();  
	c_ela.divide(2,1);
	c_ela.setSize(800,450);
	c_ela.cd(0);	
	h_el_thetap.setTitle("Momentum vs Theta");
	c_ela.getPad(0).getAxisZ().setLog(true);
	c_ela.draw(h_el_thetap);

	c_ela.cd(1);	
	h_el_thetaphi.setTitle("Phi vs Theta");
	c_ela.getPad(1).getAxisZ().setLog(true);
	c_ela.draw(h_el_thetaphi);
	c_ela.save("h2_gen_el_elastic.png");     

	EmbeddedCanvas c_pr = new EmbeddedCanvas();
	c_pr.divide(2,2);
	c_pr.setSize(800,800);
	c_pr.cd(0);
	h_pr_p.setTitle("Momentum");
	c_pr.draw(h_pr_p);
	c_pr.cd(1);
	h_pr_theta.setTitle("Theta");
	c_pr.draw(h_pr_theta);
	c_pr.cd(2);
	h_pr_phi.setTitle("Phi");
	c_pr.draw(h_pr_phi);
	c_pr.cd(3);
	h_pr_thetap.setTitle("Momentum vs Theta");
	c_pr.getPad(3).getAxisZ().setLog(true);
	c_pr.draw(h_pr_thetap);
	c_pr.save("h_rec_pr_elastic.png");

	EmbeddedCanvas c_kin = new EmbeddedCanvas();
	c_kin.setSize(800,425);
	c_kin.divide(2,1);
	c_kin.cd(0);
	//c_kin.getPad(0).getAxisY().setLog(true);
	c_kin.draw(h_el_w);
	c_kin.cd(1);
	c_kin.getPad(1).getAxisZ().setLog(true);
	c_kin.draw(h_el_w_q2);
	c_kin.save("h_gen_el_kin.png");

	

    }

    
}
