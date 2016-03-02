package predictor.core;

import predictor.core.model.*;

import predictor.core.OutputController;

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
		
		model = inputController.getModel();
		
	}
	
	/**
	 * Run STRIDE predictor
	 */
	public void runSTRIDE() {
		new STRIDEController(model).predict();
//		Molecule mol = model.getMolecules().get(0);
//		String name =  mol.getName();
//		for(int i=0; i<mol.getChains().size(); i++){
//			Chain c = mol.getChains().get(i);
//			for(int j=0; j<c.getResidues().size(); j++){
//				Residue r = c.getResidues().get(j);
//				System.out.println("REM  |---Residue---|    |--Structure--|   |-Phi-|   |-Psi-|  |-Area-|      " + name);
//				System.out.println("ASG  " + 
//						r.getName() + " " + c.getName() + "  " + r.getResidueSeqNum() + "  " + r.getResidueSeqNum() +
//						"    " + r.getAsn() + String.format("%14s", r.getSSName()) + "      " + r.getPhi() +
//						name);
//			}
//		}
	}
	
	/**
	 * Generate output file
	 */
	public void genOutputFile(String fileName) {
		
	}
	
	/**
	 * Output to standard output
	 */
	public void output(){
		new OutputController(model).output();
	}
	
	/** 
	 * Generate output object
	 */
	public void genOutput() {
		
	}
}
