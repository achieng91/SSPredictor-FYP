package predictor.ui;

import predictor.core.PredictorController;
import predictor.core.model.secStructures.SSModel;

public abstract class Controller {
	
	protected PredictorController Predictor;
	
	public Controller(String file) {
		System.out.println("PDB file ID: " + file);
		
		Predictor = new PredictorController(file);
		Predictor.STRIDESetup();
	}
	
	public void displayOutput(SSModel model) {
		System.out.println("Generating Output..");
	}
	
	public void generateOutputFile(String fileName, String filePath){
		System.out.println("Generating Output File..");
		Predictor.genOutputFile(fileName, filePath);
	}
	
	public SSModel generateOutputObject(){
		return Predictor.genOutput();
	}
	
	public void predictSTRIDE() {
		Predictor.runSTRIDE();
	}
}
