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

	double tot_beam_charge = 0.0;
	int global_event = 0;
	//double offset = 243.0;
	double offset = 209.0;
	double atten=9.8088;

	//tot_lum = tot_lum + delta_beam_charge
	//delta_beam_charge = beam_charge2 - beam_charge1
	//double beam_charge1 = 0.0;
	//double beam_charge2 = 0.0;
	
	Vector<Double> accum_charge = new Vector<Double>();
	Vector<Integer> charge_ev =  new Vector<Integer>();
	Vector<Double> accum_charge_er = new Vector<Double>();
	Vector<Double> charge_ev_er = new Vector<Double>();
	Vector<Double> current_ev = new Vector<Double>();
	Vector<Integer> current_evnt = new Vector<Integer>();
	Vector<Double> current_ev_err = new Vector<Double>();
	Vector<Double> current_evnt_err = new Vector<Double>();


	H1F h_bcg = new H1F("h_bcg","h_bcg",1000,0.0,10000);

	double elapsed_time = 0;
	double elapsed_ts = 0;
	double test_fc_charge=0.0;

	while( file_counter <= f_max ){
	    String file_in = dir_in + "clas_00"+s_run+".evio."+Integer.toString(file_counter)+".hipo";
	    File fin_temp = new File(file_in);

	    if( !fin_temp.exists() ) continue;
	    //if(f.exists() && !f.isDirectory())

	    System.out.println(" Processing File " + file_in );
	    
	    HipoDataSource hiporeader = new HipoDataSource();
	    hiporeader.open(new File(file_in) );
	    DataEvent event = null;

	    int max_event=hiporeader.getSize();
	
	    int num_ev = 0;
	    double start_unix_time = 0;
	    double end_unix_time = 0;
	    double time_stamp_counter=0;
	    int current_run=0;

 	    double test_time1 = ((DataEvent)hiporeader.gotoEvent(0)).getBank("RUN::config").getInt("unixtime",0); 
 	    double test_time2 = ((DataEvent)hiporeader.gotoEvent(max_event-1)).getBank("RUN::config").getInt("unixtime",0);
	    System.out.println(" test time for event 0 " + test_time1);
	    System.out.println(" test time for last event " + test_time2);
	    double delta_time = test_time2 - test_time1;
	    elapsed_ts+=delta_time;
	    System.out.println(" delta time for file " + delta_time );
	    
	    
	    while( num_ev < max_event ){
		//while( hiporeader.hasEvent()){
		event = (DataEvent)hiporeader.gotoEvent(num_ev);
		//DataEvent event = (DataEvent)hiporeader.getNextEvent();

		boolean runConfig_pres = event.hasBank("RUN::config");
		//boolean recBank_pres  = event.hasBank("REC::Particle");
		boolean rawScalerBank_pres = event.hasBank("RAW::scaler");
		//boolean eventBank_pres = event.hasBank("REC::Event");
		boolean configBank_pres = event.hasBank("RUN::config");

		if ( configBank_pres ){
		    current_run = event.getBank("RUN::config").getInt("run",0);
		    //if( num_ev == 7 ){
		    start_unix_time = event.getBank("RUN::config").getInt("unixtime",0);
		    //}
		    end_unix_time = event.getBank("RUN::config").getInt("unixtime",0);
		    time_stamp_counter=time_stamp_counter+event.getBank("RUN::config").getLong("timestamp",0);
		}

 		num_ev++;
		global_event++;

		test_fc_charge=test_fc_charge+charge(event);
		//System.out.println(" >> charge value from charge method " + test_fc_charge);

		if( rawScalerBank_pres ){ 		
		    double beam_charge = 0.0;
		    double beam_charge0 = 0.0;
		    double beam_charge1 = 0.0;
		    double beam_current_ungated = 0.0;
		    //System.out.println(" present " );
		    DataBank rawBank = event.getBank("RAW::scaler");
		    //rawBank.show();
	
		    for(int i=0; i<rawBank.rows(); i++){
			int slot = rawBank.getByte("slot",i);
			int chan = rawBank.getShort("channel",i);			
			offset=0;
			if(slot==0 && chan==32){
			    int fc_scaler = rawBank.getInt("value",i);
			    double true_freq = fc_scaler/(0.03333 - 0.0005);
			    double beam_current = (true_freq - offset)*atten/906.2;
			    //System.out.println(" >> true freeq " + true_freq);
			    beam_current_ungated = beam_current;
			    beam_charge1 = (0.03283*((fc_scaler/0.03283)-offset))/906.2;// beam_current * (0.03333f - 0.0005);			    
			    //System.out.println(" >> UNGATED INFO " );
			    //System.out.println(" >> Slot " + slot + " Channel " + chan + " fc scaler value " + fc_scaler);
			    //System.out.println(">> current " + beam_current*atten);			    
			    //System.out.println( " >> " + fc_scaler + " " + true_freq + " " + beam_current + " " + beam_charge1+ " " +tot_beam_charge);
			}
			//exclude for runs after feb 5th
			if( slot==0 && chan == 0 ){
			    int fc_scaler = rawBank.getInt("value",i);
			    double true_freq = fc_scaler/(0.03333 - 0.0005);
			    double beam_current = (true_freq - offset)*atten/906.2;
			    beam_charge0 = (0.03283*((fc_scaler/0.03283)-offset)*atten)/906.2;// beam_current * (0.03333f - 0.0005);			    
			    //System.out.println(" >> GATED INFO " );
			    //System.out.println(" >> Slot " + slot + " Channel " + chan + " fc scaler value " + fc_scaler);
			    //System.out.println(">> current " + beam_current);			    
			    //System.out.println( " >> " + fc_scaler + " " + true_freq + " " + beam_current + " " + beam_charge0 + " " + tot_beam_charge);

			}
			
		    }
		
		    double gated_diff = beam_charge1 - beam_charge0;
		    //System.out.println(" >> gated diff " + gated_diff );
		    tot_beam_charge=tot_beam_charge+beam_charge0; // commented out when looking at runs before feb 5th ish		    		    
		    
		    //tot_beam_charge=tot_beam_charge+gated_diff;
		    //System.out.println(" tot beam charge " + tot_beam_charge );

		    //accum_charge.add(tot_beam_charge);
		    //////accum_charge_er.add(0.0);
		    //////charge_ev.add(global_event);
		    //////charge_ev_er.add(0.0);
		    if( beam_current_ungated > 0.0 ){
			///////current_evnt.add(global_event);
			///////current_ev.add(beam_current_ungated);
			///////current_evnt_err.add(0.0);
			////////current_ev_err.add(0.0);
		    }
		    //System.out.println(" >> " + tot_beam_charge );
		}
		//if( eventBank_pres ){
		//  DataBank eventBank = event.getBank("REC::Event");
		//  double bcg = eventBank.getFloat("BCG",0);
		    //int run = eventBank.getInt("NRUN",0);
		    //System.out.println(">> " + bcg );
		    /////////h_bcg.fill(bcg);
		///}
	    }
	    file_counter++;
	    System.out.println(" >> end " + end_unix_time + " start " + start_unix_time );
	    elapsed_time = end_unix_time - start_unix_time;
	    //elapsed_ts = time_stamp_counter;

	}

	GraphErrors fc_graph = new GraphErrors();
	for( int i = 0; i < accum_charge.size(); i++ ){
	    fc_graph.addPoint(charge_ev.get(i),  accum_charge.get(i), charge_ev_er.get(i), accum_charge_er.get(i) );
	}

	GraphErrors g_current = new GraphErrors();
	for( int i = 0; i < current_ev.size(); i++ ){
	    g_current.addPoint(current_evnt.get(i), current_ev.get(i), current_evnt_err.get(i), current_ev_err.get(i) );
	}

	double target_density = 0.0701; // g/cm^3
	double atomic_mass_hydrogen = 1.00794; // g/mol
	double avogado_const = 6.0221367E23; // Number/mol	
	double target_length = 5.0; //cm
	double cm_to_microbarn = 1E30;
	double el_charge = 1.602177E-19; // Coulomb
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
	double n_pr = ( ( target_length * target_density * avogado_const ) / atomic_mass_hydrogen ) ;       	
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
	double tot_lum = lum_factor;
	System.out.println(" >> TOTAL LUMINOSITY " + tot_lum);

	
	System.out.println(" >> Writing FC txt file  now " );	


	Writer writer = null;
	try {
	    writer = new BufferedWriter(new OutputStreamWriter( new FileOutputStream("lum_out_r"+s_run+"_"+input[2]+"_"+input[2]+".txt"), "utf-8"));
 	    writer.write( " >> Number of electron " + n_el ); 
	    writer.write(" >> Number of protons  " + n_pr );
	    writer.write( " >> Luminosity factor  " + lum_factor );    
	    writer.write(" >> elapsed time " + elapsed_time );  
	    writer.write(" >> time stamp time " + elapsed_ts); 
	    writer.write(" >> charge from time " + charge_from_time ); 
	    writer.write(" >> nel from time " + nel_from_time );   
	    writer.write(" >> TOTAL ACCUMULATED CHARGE " + tot_beam_charge );
	    writer.write(" >> TARGET DENSITY " + target_density );
	    writer.write(" >> TOTAL LUMINOSITY " + tot_lum);
	} catch (IOException ex) {
	    // Report
	} finally {
	    try {writer.close();} catch (Exception ex) {/*ignore*/}
	}



	System.out.println(" >> Plotting FC charge now " );	
	EmbeddedCanvas c_fc = new EmbeddedCanvas();
	c_fc.divide(1,1);
	c_fc.cd(0);
	c_fc.setSize(800,800);
	fc_graph.setMarkerSize(1);
	fc_graph.setTitleX("Nth Readout");
	fc_graph.setTitleY("Tot. Acc. Charge");
	c_fc.draw(fc_graph);
	c_fc.save("fc_graph_r"+s_run+".png");

	EmbeddedCanvas c_bcg = new EmbeddedCanvas();
	c_bcg.divide(1,1);
	c_bcg.cd(0);
	c_bcg.setSize(800,800);
	h_bcg.setTitleX("BCG");
	c_bcg.draw(h_bcg);
	c_bcg.save("hist_bcg_r"+s_run+".png");

	EmbeddedCanvas c_current = new EmbeddedCanvas();
	c_current.divide(1,1);
	c_current.cd(0);
	c_current.setSize(800,400);
	g_current.setMarkerSize(1);
	g_current.setTitleX("Nth Readout");
	g_current.setTitleY("Current [nA]");
	c_current.draw(g_current);
	c_current.save("g_current_r"+s_run+".png");	
	
	double bcg_min = h_bcg.getMin();
	double bcg_max = h_bcg.getMax();
	System.out.println(" >> BCG Min: " + bcg_min + " BCG MAX: " + bcg_max );

    }

    //used to validate with Nick Markov.
    //variables described online or in logbook
    public static double charge(DataEvent event){
	double cLocal = 0.0;
	if(event.hasBank("RAW::scaler")){
	    DataBank bb = event.getBank("RAW::scaler");
	    for(int ii=0;ii<bb.rows();ii++){
		int slot = bb.getByte("slot",ii);
		int chan = bb.getShort("channel",ii); 
		if(slot==0 && chan==32){
		    cLocal = bb.getInt("value",ii);
		    //System.out.println("clocal " + cLocal);
		    cLocal = 0.03283*(cLocal/0.03283-100)/906.2;
		}
	    }
	}	    
	return cLocal;		
    }
    
}
