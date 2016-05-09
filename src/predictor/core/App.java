package predictor.core;

import java.io.IOException;
import java.util.logging.FileHandler;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

import predictor.core.model.secStructures.SSModel;
import predictor.ui.*;

/***
 * 
 * @author Anson
 *
 * Secondary Structure Predictor
 * 
 */
public class App {
	
	private static boolean enableUI = false;
	private static boolean outputFile = true;
	private static boolean outputObject = true;
	private static String predictor = "STRIDE";
	private static final Logger log = Logger.getLogger("main");
	private static SSModel model;

	public static void main(String[] args) {
	    FileHandler fh;  

	    try {
	        fh = new FileHandler("Log.txt");  
	        log.addHandler(fh);
	        SimpleFormatter formatter = new SimpleFormatter();  
	        fh.setFormatter(formatter);
	        log.setUseParentHandlers(false);
	    } catch (SecurityException e) {  
	        e.printStackTrace();  
	    } catch (IOException e) {  
	        e.printStackTrace();
	    }  
		
		
		log.info("Application starts");
		String dir = "res/test/4hhb.pdb";
   		
   		Controller controller = ControllerFactory.createController(enableUI, dir);
   		
   		if(predictor.equalsIgnoreCase("STRIDE")){
				controller.predictSTRIDE();
			}
   		
   		if(outputObject){
   			model = controller.generateOutputObject();
   		}
   		
   		if(enableUI){
   			//TODO: Add display for GUI
   		} else {
   			controller.displayOutput(model);
   		}
   		
   		if(outputFile){
   			controller.generateOutputFile("4hhb", "");
   		}
   	}
}
