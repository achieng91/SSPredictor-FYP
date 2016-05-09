package predictor.core;

import predictor.core.model.*;
import predictor.core.model.secStructures.SSModel;

import org.biojava.bio.structure.Structure;

import predictor.core.OutputController;

/***
 * 
 * @author Anson
 *
 * Does initial setup and determine which predictor to use
 */
public class PredictorController {

	protected String filePath;
	protected Model model = new Model();
	protected Molecule mol;
	protected Structure struc;
	
	public PredictorController(String filePath) {
		this.filePath = filePath;
	}
	
	/**
	 * Initial setup 
	 * Must be run first before using other functions
	 */
	public void STRIDESetup() {
		InputController inputController = new InputController(model);
		inputController.getInput(filePath);
		struc = inputController.getStruc();
		
		mol = model.getMolecules().get(0);
		SetupController.addHAtoms(mol);
		SetupController.addPhiPsiAngles(mol);
		
		model = inputController.getModel();
		
	}
	
	/**
	 * Run STRIDE predictor
	 */
	public void runSTRIDE() {
		new Predictor(model).predict();
	}
	
	/**
	 * Generate output file
	 */
	public void genOutputFile(String fileName, String filePath) {
		new OutputController(model).outputFile(fileName, filePath);
	}
	
	/**
	 * Output to standard output
	 * @param model2 
	 */
	public void output(SSModel model2){
		new OutputController(model).output(struc, model2);
	}
	
	/** 
	 * Generate output object
	 */
	public SSModel genOutput() {
		return new OutputController(model).createSSObject();
	}
}
