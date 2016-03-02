package predictor.ui;

import predictor.core.PredictorController;

public abstract class Controller {
	
	protected PredictorController Predictor;
	
	public Controller(String file) {
		System.out.println("PDB file ID: " + file);
		
		Predictor = new PredictorController(file);
		Predictor.STRIDESetup();
	}
	
	public void displayOutput() {
//		System.out.println("Generate Output..");
	}
	
	
	public void predictSTRIDE() {
		Predictor.runSTRIDE();
	}
}
