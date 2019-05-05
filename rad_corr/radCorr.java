import java.io.*;
import java.util.*;
import org.jlab.clas.reco.io.*;
import org.jlab.clas.physics.*;
import org.jlab.clas.reactions.*;
import org.jlab.clas.physics.*;
import org.jlab.clas.fastmc.*;
import org.jlab.physics.io.LundReader;


public class radCorr{

    public static void main(String[] input){


	String inputFile = input[0];
	
	
	
	
	
	LundReader      reader = new LundReader();
	reader.addFile(inputFile);
	reader.open();

	while(reader.next()==true){
	    PhysicsEvent genEvent = reader.getEvent();
	    
	    if(genEvent.countByPid(11)>0){
		Particle el = genEvent.getParticle(0);//getParticleByPid(11,11);
	
		//double px = el.px();
		System.out.println(el);
		    //System.out.println("==============>  EVENT PRINTOUT");
		    //System.out.println("----->  GENERATED EVENT");
		System.out.println(genEvent.toLundString());
	    }
	}



    }
}
