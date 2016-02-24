package predictor.core;

import predictor.core.model.Model;

import predictor.core.OutputController;
import predictor.util.SSFile;

/***
 * 
 * @author Anson
 *
 * Does initial setup and determine which predictor to use
 */
public class PredictorController {

	protected String filePath;
	protected Model model;
	
	public PredictorController(String filePath) {
		this.filePath = filePath;
	}
	
	/**
	 * Initial setup 
	 * Must be run first before using other functions
	 */
	public void STRIDESetup() {
		InputController inputController = new InputController();
		inputController.getInput(filePath);
		inputController.createModel();
		inputController.addHAtoms();
		inputController.addPhiPsiAngles();
		
		this.model = inputController.getModel();
	}
	
	/**
	 * Run STRIDE predictor
	 */
	public void runSTRIDE() {
		STRIDEController stride = new STRIDEController(this.model);
		stride.predict();
	}
	
	/**
	 * Generate output file
	 */
	public SSFile genOutputFile(String fileType) {
		return new OutputController(this.model).createSSFile(fileType);
	}
	
	/** 
	 * Generate output object
	 */
	public void genOutput() {
		
	}
}
