import java.io.*;
import java.util.*;
import org.jlab.io.hipo.*;
import org.jlab.io.base.DataEvent;
import org.jlab.io.base.DataBank;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.groot.data.*;
import org.jlab.groot.graphics.EmbeddedCanvas;
import java.lang.Integer; 

public class faradayCupCalculatorH4{


    public static void main(String[] input){

	String dir_in = input[0];
	String s_run = input[1];
	int f_min = Integer.valueOf(input[2]);
	int f_max = Integer.valueOf(input[3]);
	int file_counter=f_min;

	double tot_beam_charge = 0.0;
	int global_event = 0;
	double offset = 243.0;
	//double offset = 209.0;
	double atten=1.0;//9.8088;

	int key_diff_cut = 3;

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

	Vector<Integer> v_global_event = new Vector<Integer>();

	TreeMap<Integer, Double> tmap_global_event_charge = new TreeMap<Integer, Double>();
	TreeMap<Integer, Double> tmap_global_event_gatedcharge = new TreeMap<Integer, Double>();
	TreeMap<Integer, Double> tmap_global_event_ungatedcharge = new TreeMap<Integer, Double>();
	TreeMap<Integer, Double> tmap_unixtime = new TreeMap<Integer, Double>();                                                                                                         
        Vector<Integer> v_last_event = new Vector<Integer>();
                                                                          
	TreeMap<Integer, Double> tmap_beam_current_gated = new TreeMap<Integer,Double>();           
	

	H1F h_bcg = new H1F("h_bcg","h_bcg",1000,0.0,10000);

	double elapsed_time = 0;
	double elapsed_ts = 0;
	double test_fc_charge=0.0;

	File file = new File(dir_in);
	BufferedReader reader = null;
	String st;

	try {
	    reader = new BufferedReader(new FileReader(file));

	while ((st = reader.readLine()) != null) {

	    System.out.println(st); 
	    

	    String file_in = st;
	    // dir_in ;//+ "clas_00"+s_run+".evio."+Zeros+Integer.toString(file_counter)+"-"+Zeros+Integer.toString(file_counter+4)+".hipo";
	    
	    System.out.println(" >> " + file_in );
	    File fin_temp = new File(file_in);

	    if( !fin_temp.exists() ) continue;
	    
	    System.out.println(" Processing File " + file_in );
	    
	    HipoDataSource hiporeader = new HipoDataSource();
	    hiporeader.open(file_in);
	    DataEvent event = null;

	    int max_event=hiporeader.getSize();
	    System.out.println(" Number of events " + Integer.toString(max_event) );

	    int num_ev = 0;
	    double start_unix_time = 0;
	    double end_unix_time = 0;
	    double time_stamp_counter=0;
	    int current_run=0;

 	    double test_time1 = ((DataEvent)hiporeader.gotoEvent(0)).getBank("RUN::config").getInt("unixtime",0); 
 	    double test_time2 = ((DataEvent)hiporeader.gotoEvent(max_event-1)).getBank("RUN::config").getInt("unixtime",0);
	    
	    //System.out.println(" test time for event 0 " + test_time1);
	    //System.out.println(" test time for last event " + test_time2);
	    double delta_time = test_time2 - test_time1;
	    elapsed_ts+=delta_time;
	    //System.out.println(" delta time for file " + delta_time );
	    
	    int start_event = 1000000000;
	    int end_event = 0;
	    double start_charge = 0;//1E9; 
	    double end_charge = 0;
	    int my_counter = 0;

	    while( num_ev < 200){ //max_event ){
	    	event = (DataEvent)hiporeader.gotoEvent(num_ev);


		if(Math.random() < 0.5) {
		    
		    my_counter = my_counter + 5;
		}
		else{
		    my_counter++;
		}

		System.out.println(" >> " + Integer.toString(my_counter) );
		
		boolean runConfig_pres = event.hasBank("RUN::config");
		boolean rawScalerBank_pres = event.hasBank("RUN::scaler");       
		boolean configBank_pres = event.hasBank("RUN::config");

		int config_event = event.getBank("RUN::config").getInt("event",0);

		if ( config_event >  end_event ){
		    end_event = config_event;
		}

		if( config_event < start_event ){
		    start_event = config_event;
		}

		if ( configBank_pres ){
		    current_run = event.getBank("RUN::config").getInt("run",0);
		    
		    start_unix_time = event.getBank("RUN::config").getInt("unixtime",0);
		    end_unix_time = event.getBank("RUN::config").getInt("unixtime",0);
		    time_stamp_counter=time_stamp_counter+event.getBank("RUN::config").getLong("timestamp",0);
		    tmap_unixtime.put(config_event,start_unix_time);		  

		}

 		num_ev++;
		global_event++;
		
		/*
		if( event.hasBank("RAW::scaler") ){
		    DataBank rawScalerBank  = event.getBank("RAW::scaler");
		    rawScalerBank.show();
		    System.out.println(" has raw scaler bank " );


		    //I [nA] = (S [Hz] - offset ) / slope * attenuation;
		    offset = 374.80;
		    atten = 9.80880;

		    System.out.println(" event " + Integer.toString( config_event) );
		    System.out.println(" local event " + Integer.toString(num_ev) );
		    //run 5038 slope = 906.20 offset = 374.80 and atten = 9.80880
		    for(int i=0; i<rawScalerBank.rows(); i++){
			int slot = rawScalerBank.getByte("slot",i);
			int chan = rawScalerBank.getShort("channel",i);			
			//offset=0;
			if(slot==64 && chan==48){
			    System.out.println(" scaler index " + Integer.toString(i));

			    int fc_scaler = rawScalerBank.getInt("value",i);
			    double true_freq = fc_scaler/(0.03333 - 0.0005);
			    double delta_t = start_unix_time-1539254875;
			    System.out.println(" >>unix time " + start_unix_time );

			    double beam_current =  (fc_scaler/delta_t - offset) / 906.20  * atten; //)*atten)/906.2; // ((true_freq - offset)*atten)/906.2;
			    System.out.println(" >> true freeq " + true_freq);
			    //double beam_current_ungated_temp = beam_current;
			    double beam_charge1 = (fc_scaler - (offset*delta_t) )/906.2;// beam_current * (0.03333f - 0.0005);			    

			    System.out.println(" >>>>>>>>>>> UN-GATED INFO <<<<<<<<<<< " );
			    System.out.println(" >> Slot " + slot + " Channel " + chan + " fc scaler value " + fc_scaler);
			    System.out.println(" >> current " + beam_current);			    
			    System.out.println( " >> fcs" + fc_scaler + " truefreq " + true_freq + " bc " + beam_current + " beamcharge " + beam_charge1+ " " +tot_beam_charge);
			}
			if( slot==64 && chan == 16 ){
			    System.out.println(" scaler index " + Integer.toString(i));

			    int fc_scaler = rawScalerBank.getInt("value",i);
			    double true_freq = fc_scaler/(0.03333 - 0.0005);
			    //double beam_current = (true_freq - offset)*atten/906.2;
			    double beam_current = (1.0/(0.03333 - 0.0005))* (fc_scaler - (0.03333 - 0.0005)*offset) / 906.20  * atten; //*(time_stamp_counter-1539254875))*atten)/906.2; // ((true_freq - offset)*atten)/906.2;
			    tmap_beam_current_gated.put(config_event,beam_current);

			    double beam_charge0 = (0.03283*((fc_scaler/0.03283)-offset)*atten)/906.2;// beam_current * (0.03333f - 0.0005);			    
			    System.out.println(" >>>>>>>>>>> GATED INFO <<<<<<<<<<< " );
			    System.out.println(" >> Slot " + slot + " Channel " + chan + " fc scaler value " + fc_scaler);
			    System.out.println(" >> current " + beam_current);			    
			    System.out.println( " >> fcs" + fc_scaler + " truefreq " + true_freq + " bc " + beam_current + " beamcharge " + beam_charge0+ " " +tot_beam_charge);

			}
			
		    }



		}
		
		*/
		if( true ){// rawScalerBank_pres ){ 		
 		    double beam_charge = 0.0;
		    double beam_current_ungated = 0.0;

		    //DataBank rawBank = event.getBank("RUN::scaler");
		    //float live_time = rawBank.getFloat("livetime",0);
		    //System.out.println(" >> event " + Integer.toString(config_event));

		    if( true ){ // live_time >= 0.0 ){
		    
			//double fcup_ug = rawBank.getFloat("fcup",0);
			//double fcup_g = rawBank.getFloat("fcupgated",0);
			//double fcup_lt = fcup_ug/fcup_g;
			//System.out.println(" >> fcup ug " + Double.toString(fcup_ug) + " g " + Double.toString(fcup_g) + " lt " + Double.toString(fcup_lt) );
			
			//put event number and the event's gated diff. sum later
			//tmap_global_event_charge.put(config_event,gated_diff);//tot_beam_charge);
			config_event = my_counter;
			tmap_global_event_gatedcharge.put(config_event, my_counter*3.0 );//fcup_g );
			tmap_global_event_ungatedcharge.put(config_event, my_counter*3.0 );//fcup_ug );			
		    }		   
		}
	    }

	    //System.out.println(" >> end " + end_unix_time + " start " + start_unix_time );
	    elapsed_time = end_unix_time - start_unix_time;
	    System.out.println(" >> Start event " + Double.toString(start_event) + " End event " + Double.toString(end_event) );
	    v_last_event.add(end_event);

  	    //System.out.println(" >> charge at start event " + charge_at_start  + " charge at end event " + charge_at_end );

	}
	
    }
    catch (FileNotFoundException ex) {
        // handle FileNotFoundException
    }
    catch (IOException ex) {
        // handle IOException
    }
    finally {
        if (reader != null) {
            try {
                reader.close();
            }
            catch (IOException ex) {
                // handle IOException
            }
        }
    }
	// plot the accumulated charge by order of event 
	GraphErrors fc_graph = new GraphErrors();
	GraphErrors fcug_graph = new GraphErrors();
	GraphErrors fcg_graph = new GraphErrors();
	H1F h_fcg_event_diff = new H1F("h_fcg_event_diff","h_fcg_event_diff",2000, 0.0, 2000.0);
	H1F h_fcug_event_diff = new H1F("h_fcug_event_diff","h_fcug_event_diff",2000, 0.0, 2000.0);
	

	double acc_beam_charge_sum = 0;
	double acc_charge_g = 0;
	double acc_charge_ug = 0;
	for(Map.Entry<Integer,Double> entry : tmap_global_event_charge.entrySet()) {
	    Integer key = entry.getKey();
	    Double value = entry.getValue();
	    //acc_beam_charge_sum=acc_beam_charge_sum + (Double)value;
	    acc_beam_charge_sum= (Double)value;

	    double key_temp = key.doubleValue();
	    //if( key <= 2406468 ){

	       
	    //	for( double last_event_in_file : v_last_event ){

	    //	    if( Math.abs(last_event_in_file - key_temp) < 500 ){
	    //		System.out.println(" run config event " + Integer.toString(key) + "fc gated diff " + value  + " acc charge " + acc_beam_charge_sum);
	    //	    }
	    //	}
	    if( key < 100 ){
		fc_graph.addPoint( key_temp, acc_beam_charge_sum, 0.0, 0.0 );
	    }
	    
	}
	acc_beam_charge_sum=0;

	System.out.println("gated Entry set values last: ");
	System.out.println(tmap_global_event_gatedcharge.firstEntry());
	System.out.println("gated Entry set values first: ");
	System.out.println(tmap_global_event_gatedcharge.lastEntry());

	
	double true_acc_charge_sum = 0; 
	int sub_key = tmap_global_event_gatedcharge.firstEntry().getKey();
	double sub_val = tmap_global_event_gatedcharge.firstEntry().getValue();

	ArrayList<DataLine> gap_line = new ArrayList<DataLine>();

	for(Map.Entry<Integer,Double> entry :  tmap_global_event_gatedcharge.entrySet()) {
	    Integer key = entry.getKey();
	    Double value = entry.getValue();
	    //acc_beam_charge_sum=acc_beam_charge_sum + (Double)value;
	    acc_beam_charge_sum = (Double)value;
	    acc_charge_g = (Double)value;
	    double key_temp = key.doubleValue();

	    //if( key > 2406466 ){ System.out.println(" gated " + Double.toString(acc_charge_g) + " event " + Integer.toString(key)); }
	    //System.out.println(" gated " + Double.toString(acc_charge_g) + " event " + Integer.toString(key));
	    //if( key <= 2406468 ){	       
	    if( acc_charge_g <= 0 ){
		System.out.println(" >> problem here at event "  + Double.toString(key_temp) );
	    }

	    if ( tmap_global_event_gatedcharge.higherEntry(key) == null ) continue; 

 	    int next_key = tmap_global_event_gatedcharge.higherEntry(key).getKey();
	    int key_diff_above = next_key - key;


	    // dealing with the first key and the previous key to key 0 as being null
	    int prev_key = -100;

	    if( tmap_global_event_gatedcharge.lowerEntry(key) == null ){
		prev_key = tmap_global_event_gatedcharge.firstEntry().getKey();
		System.out.println(" >> There is no previous key to key  " + Integer.toString(key) + " setting previous key to first entry " + Integer.toString(prev_key) );
	    }
	    else{
		prev_key = tmap_global_event_gatedcharge.lowerEntry(key).getKey();	
	    }

	    int key_diff_below = key - Math.abs(prev_key);
	    
	    boolean no_lone_point = false;
	    //trying to skip the points with no adjacent items nearby
	    if( (key_diff_above > key_diff_cut && key_diff_below > key_diff_cut) ||
	        (tmap_global_event_gatedcharge.lowerEntry(key) == null && key_diff_above > key_diff_cut)  )		
		{
		    System.out.println(" >> should be skipping this point with no adjacent neighbors - key is " + Integer.toString(key) );
		    sub_key = tmap_global_event_gatedcharge.higherEntry(key).getKey();
		    no_lone_point=true;
		    //continue;
		}

	    if( key_diff_above > key_diff_cut && no_lone_point == false){
		System.out.println(" current event " + Integer.toString(key) );
		System.out.print("  event lower -> ");
		System.out.println(tmap_global_event_gatedcharge.lowerEntry(key));
		System.out.print("  event higher -> ");
		System.out.println(tmap_global_event_gatedcharge.higherEntry(key));
		
		
		
		

		double subset_charge_sum = value - sub_val; // current - starting value in the subset of continuous events
		System.out.println( " >> subset charge sum " + Double.toString(subset_charge_sum) + " current val " +Double.toString(value) + " - first val " + Double.toString(sub_val) ); 
		true_acc_charge_sum = true_acc_charge_sum + subset_charge_sum;
		System.out.println(" True acc charge " + Double.toString(true_acc_charge_sum) );
		
		sub_val = tmap_global_event_gatedcharge.higherEntry(key).getValue();
		System.out.println(" setting sub_val to next entry value " + Double.toString(sub_val) );
		System.out.println(" no break between current event " + Integer.toString(key) + " and starting event " + Integer.toString(sub_key) );

		// capture points that are good

		DataLine l_mark_below = new DataLine( (double)key, 0.0, (double)key, tmap_global_event_gatedcharge.lastEntry().getValue() );
		gap_line.add(l_mark_below);
		
		DataLine l_mark_above = new DataLine( (double)sub_key, 0.0, (double)sub_key, tmap_global_event_gatedcharge.lastEntry().getValue() );
		gap_line.add(l_mark_above);
		
		sub_key = tmap_global_event_gatedcharge.higherEntry(key).getKey();

	    }
	    System.out.println(" >> filling histogram with difference in key values " + Double.toString((double)key_diff_above));
	    h_fcg_event_diff.fill((double)key_diff_above);
	    
	    
	    if( key < 105 ){
		fcg_graph.addPoint( key_temp, acc_beam_charge_sum, 0.0, 0.0 );
		
	    }
	}
	acc_beam_charge_sum = 0;

	System.out.println("ungated Entry set values last: ");
	System.out.println(tmap_global_event_ungatedcharge.firstEntry());
	System.out.println("ungated Entry set values first: ");
	System.out.println(tmap_global_event_ungatedcharge.lastEntry());


	for(Map.Entry<Integer,Double> entry :  tmap_global_event_ungatedcharge.entrySet()) {
	    Integer key = entry.getKey();
	    Double value = entry.getValue();
	    //acc_beam_charge_sum=acc_beam_charge_sum + (Double)value;
	    acc_beam_charge_sum= (Double)value;
	    acc_charge_ug= (Double)value;

	    double key_temp = key.doubleValue();

	    if ( tmap_global_event_ungatedcharge.higherEntry(key) == null || tmap_global_event_ungatedcharge.lowerEntry(key) == null ) continue;

 	    int next_key = tmap_global_event_ungatedcharge.higherEntry(key).getKey();
	    int prev_key = tmap_global_event_ungatedcharge.lowerEntry(key).getKey();
	    int key_diff = next_key - key;

	    //if( key > 2406466 ){ System.out.println(" ungated " + Double.toString(acc_charge_ug) + " event " + Integer.toString(key)); }
	    //System.out.println(" ungated " + Double.toString(acc_charge_ug) + " event " + Integer.toString(key)); 
	    //if( key <= 2406468 ){
	
	    if( key < 100 ){

		fcug_graph.addPoint( key_temp, acc_beam_charge_sum, 0.0, 0.0 );
		h_fcug_event_diff.fill((double)key_diff);
	    }
		//}
	}


	for(Map.Entry<Integer,Double> entry : tmap_unixtime.entrySet()) {
	    Integer key = entry.getKey();
	    Double value = entry.getValue();

	    //if( key <= 2406468 ){
		//System.out.println(" run config event " + Integer.toString(key) + " unix time " + Double.toString(value));
	    //}
	}



	GraphErrors g_current = new GraphErrors();
	//for(Map.Entry<Integer,Double> entry :  tmap_beam_current_gated.entrySet()) {
	//  Integer key = entry.getKey();
	//   Double value = entry.getValue();
	    //acc_beam_charge_sum=acc_beam_charge_sum + (Double)value;
	// double beam_current_val= (Double)value;
	    //   
	//    double key_temp = key.doubleValue();
	//   if( key <= 2406468 ){

	//	g_current.addPoint( key_temp, beam_current_val, 0.0, 0.0 );
	//   }
	//}

	//for( int i = 0; i < current_ev.size(); i++ ){
	//    g_current.addPoint(current_evnt.get(i), current_ev.get(i), current_evnt_err.get(i), current_ev_err.get(i) );
	//}

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

	double n_el_ug = acc_charge_ug *nC_to_C / el_charge;  
	double n_el_g = acc_charge_g *nC_to_C / el_charge;  

	double n_pr = ( ( target_length * target_density * avogado_const ) / atomic_mass_hydrogen ) ;       	
	double lum_factor = (n_el*n_pr)/cm_to_microbarn;       	



	System.out.println(" >> Total Gated Charge " +  acc_charge_g );
	System.out.println(" >> Total UNGated Charge " +  acc_charge_ug );
	System.out.println( " >> Number of electron " + n_el );
	System.out.println( " >> Number of electron GATED " + n_el_g );
	System.out.println( " >> Number of electron UNGATED " + n_el_ug );
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
	fc_graph.setTitleX("RUN::config event");
	fc_graph.setTitleY("Tot. Acc. Charge");
	for( DataLine ll : gap_line ){
	    c_fc.draw(ll);
	}
	c_fc.draw(fc_graph,"same");
	c_fc.save("fc_graph_r"+s_run+".png");

	EmbeddedCanvas c_fc2 = new EmbeddedCanvas();
	c_fc2.setSize(1600,800);
	c_fc2.divide(2,1);
	c_fc2.cd(0);
	fcg_graph.setMarkerSize(1);
	fcg_graph.setTitleX("RUN::config event");
	fcg_graph.setTitleY("FC Gated ");
	for( DataLine ll : gap_line ){
	    ll.setLineColor(2);
	    c_fc2.draw(ll);
	}
	c_fc2.draw(fcg_graph,"same");

	c_fc2.cd(1);
	fcug_graph.setMarkerSize(1);
	fcug_graph.setTitleX("RUN::config event");
	fcug_graph.setTitleY("FC Ungated ");
	c_fc2.draw(fcug_graph);
	c_fc2.save("fc2_graph_r"+s_run+".png");

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

	EmbeddedCanvas c_diff = new EmbeddedCanvas();
	c_diff.divide(2,1);
	c_diff.setSize(800,400);	
	c_diff.cd(0);
	h_fcug_event_diff.setTitleY("Freq Diff");
	h_fcg_event_diff.setTitleY("Freq Diff");
	c_diff.getPad(0).getAxisY().setLog(true);
	//c_diff.getPad(0).getAxisX().setLog(true);
	c_diff.draw(h_fcug_event_diff);
	c_diff.cd(1);
	c_diff.getPad(1).getAxisY().setLog(true);
	//c_diff.getPad(1).getAxisX().setLog(true);
	c_diff.draw(h_fcg_event_diff);
	c_diff.save("h_diff_event_numbers.png");

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
