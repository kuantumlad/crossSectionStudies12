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

	int max_event=hiporeader.getSize();
	
	int num_ev = 0;

	H1F h_el_p = new H1F("h_el_p","h_el_p",200,0.0,13.0);
	H1F h_el_theta = new H1F("h_el_theta","h_el_theta",200,0.0,30.0);
	H1F h_el_phi = new H1F("h_el_phi","h_el_phi",300,-200.0,200.0);
	H2F h_el_thetap = new H2F("h_el_thetap","h_el_thetap",100,0.0,13.0, 100, 0.0, 30.0);

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
	c_el.cd(3);
	h_el_thetap.setTitle("Momentum vs Theta");
	c_el.draw(h_el_thetap);
	c_el.save("h_gen_el_elastic_peterbostedgen.png");

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
	c_pr.draw(h_pr_thetap);
	c_pr.save("h_gen_pr_elastic_peterbostedgen.png");


	

    }

    
}
