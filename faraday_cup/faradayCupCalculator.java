import java.io.*;
import java.util.*;
import org.jlab.io.hipo.*;
import org.jlab.io.base.DataEvent;
import org.jlab.io.base.DataBank;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.groot.data.*;
import org.jlab.groot.graphics.EmbeddedCanvas;

public class faradayCupCalculator{


    public static void main(String[] input){

	String dir_in = input[0];
	String s_run = input[1];
	int f_min = Integer.valueOf(input[2]);
	int f_max = Integer.valueOf(input[3]);
	int file_counter=f_min;

	double tot_beam_charge=0.0;
	int global_event = 0;
	
	//tot_lum = tot_lum + delta_beam_charge
	//delta_beam_charge = beam_charge2 - beam_charge1
	double beam_charge1 = 0.0;
	double beam_charge2 = 0.0;
	
	Vector<Double> accum_charge = new Vector<Double>();
	Vector<Integer> charge_ev =  new Vector<Integer>();
	Vector<Double> accum_charge_er = new Vector<Double>();
	Vector<Double> charge_ev_er = new Vector<Double>();

	H1F h_bcg = new H1F("h_bcg","h_bcg",1000,0.0,10000);

	double elapsed_time = 0;
	double elapsed_ts = 0;

	while( file_counter <= f_max ){
	    String file_in = dir_in + "out_clas_00"+s_run+".evio."+Integer.toString(file_counter)+".hipo";


		    
	    System.out.println(" Processing File " + file_in );
	    
	    HipoDataSource hiporeader = new HipoDataSource();
	    hiporeader.open(new File(file_in) );
	    DataEvent event = null;

	    int max_event=hiporeader.getSize();
	
	    int num_ev = 0;
	    double start_unix_time = 0;
	    double end_unix_time = 0;
	    double time_stamp_counter=0;
	    while( num_ev < max_event ){
		event = (DataEvent)hiporeader.gotoEvent(num_ev);
		boolean runConfig_pres = event.hasBank("RUN::config");
		boolean recBank_pres  = event.hasBank("REC::Particle");
		boolean rawScalerBank_pres = event.hasBank("RAW::scaler");
		boolean eventBank_pres = event.hasBank("REC::Event");
		boolean configBank_pres = event.hasBank("RUN::config");

		if ( configBank_pres ){
		    if( num_ev == 7 ){
			start_unix_time = event.getBank("RUN::config").getInt("unixtime",0);
		    }
		    end_unix_time = event.getBank("RUN::config").getInt("unixtime",0);
		    time_stamp_counter=time_stamp_counter+event.getBank("RUN::config").getLong("timestamp",0);
		}

 		num_ev++;
		global_event++;

		if( rawScalerBank_pres ){ 		
		    double beam_charge = 0.0;
		    //System.out.println(" present " );
		    DataBank rawBank = event.getBank("RAW::scaler");
		    for(int i=0; i<rawBank.rows(); i++){
			int slot = rawBank.getByte("slot",i);
			int chan = rawBank.getShort("channel",i);
			if(slot==0 && chan==32){
			    int fc_scaler = rawBank.getInt("value",i);
			    double true_freq = fc_scaler/(0.03333 - 0.0005);
			    double beam_current = (true_freq - 100.0)/906.2;
			    beam_charge = 0.03283*(fc_scaler/0.03283-100)/906.2;// beam_current * (0.03333f - 0.0005);			    
			    //System.out.println(">> current " + beam_current);			    
			    //System.out.println( " >> " + fc_scaler + " " + true_freq + " " + beam_current + " " + beam_charge+ " " +tot_beam_charge);
			}
		    }
		    tot_beam_charge=tot_beam_charge+beam_charge;		    		    

		    accum_charge.add(tot_beam_charge);
		    accum_charge_er.add(0.0);
		    charge_ev.add(global_event);
		    charge_ev_er.add(0.0);
		    //System.out.println(" >> " + tot_beam_charge );
		}
		if( eventBank_pres ){
		    DataBank eventBank = event.getBank("REC::Event");
		    double bcg = eventBank.getFloat("BCG",0);
		    //int run = eventBank.getInt("NRUN",0);
		    //System.out.println(">> " + bcg );
		    h_bcg.fill(bcg);
		}
	    }
	    file_counter++;
	    System.out.println(" >> end " + end_unix_time + " start " + start_unix_time );
	    elapsed_time = end_unix_time - start_unix_time;
	    elapsed_ts = time_stamp_counter;

	}

	GraphErrors fc_graph = new GraphErrors();
	for( int i = 0; i < accum_charge.size(); i++ ){
	    fc_graph.addPoint(charge_ev.get(i),  accum_charge.get(i), charge_ev_er.get(i), accum_charge_er.get(i) );
	}

	double target_density = 0.0703; // g/cm^3
	double atomic_mass_hydrogen = 1.00794; // g/mol
	double avogado_const = 6.0221367e23; // Number/mol	
	double target_length = 5.0; //cm
	double cm_to_microbarn = 1e30;
	double el_charge = 1.602177e-19; // Coulomb
	double nC_to_C = 1e-9;
	double ns_to_s = 1e-9;

	System.out.println(" >> target density          " + target_density );
	System.out.println(" >> atomic mass of hydrogen " + atomic_mass_hydrogen );
	System.out.println(" >> avogado constant        " + avogado_const );
	System.out.println(" >> targe length            " + target_length );
	System.out.println(" >> cm to microbarn         " + cm_to_microbarn );
	System.out.println(" >> electron charge         " + el_charge );
	System.out.println(" >> nC to C                 " + nC_to_C );	

	double n_el = (tot_beam_charge * nC_to_C)/el_charge;
	double n_pr =  ( ( target_length * target_density * avogado_const ) / atomic_mass_hydrogen ) ;       	
	double lum_factor = (n_el*n_pr)/cm_to_microbarn;       	

	System.out.println( " >> Number of electron " + n_el );
	System.out.println( " >> Number of protons  " + n_pr );
	System.out.println( " >> Luminosity factor  " + lum_factor );

	System.out.println(" >> elapsed time " + elapsed_time );
	System.out.println(" >> time stamp time " + elapsed_ts);

	double charge_from_time = elapsed_time * 5.5 * nC_to_C;
	double nel_from_time = charge_from_time/el_charge;
	
	System.out.println(" >> charge from time " + charge_from_time );
	System.out.println(" >> nel from time " + nel_from_time );

	System.out.println(" >> TOTAL ACCUMULATED CHARGE " + tot_beam_charge );
	System.out.println(" >> TARGET DENSITY " + target_density );
	System.out.println(" >> TARGET LENGTH " + target_length );
	double tot_lum = tot_beam_charge*lum_factor;
	System.out.println(" >> TOTAL LUMINOSITY " + tot_lum);

	System.out.println(" >> Plotting FC charge now " );	
	EmbeddedCanvas c_fc = new EmbeddedCanvas();
	c_fc.divide(1,1);
	c_fc.cd(0);
	c_fc.setSize(800,800);
	fc_graph.setMarkerSize(1);
	c_fc.draw(fc_graph);
	c_fc.save("fc_graph.png");

	EmbeddedCanvas c_bcg = new EmbeddedCanvas();
	c_bcg.divide(1,1);
	c_bcg.cd(0);
	c_bcg.setSize(800,800);
	h_bcg.setTitleX("BCG");
	c_bcg.draw(h_bcg);
	c_bcg.save("hist_bcg.png");
	
	double bcg_min = h_bcg.getMin();
	double bcg_max = h_bcg.getMax();
	System.out.println(" >> BCG Min: " + bcg_min + " BCG MAX: " + bcg_max );

    }
    

}
