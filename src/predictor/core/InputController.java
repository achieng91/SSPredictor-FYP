package predictor.core;

import predictor.model.Model;
import predictor.model.Molecule;
import predictor.model.Chain;
import predictor.model.Residue;
import predictor.model.Atom;
import util.FileReader;

import org.biojava.bio.structure.Structure;

public class InputController {
	
	protected Structure struc;
	protected Model model = new Model();
	
	public void getInput(String filePath) {
		struc = FileReader.readFile(filePath);
	}
	
	public void createModel() {
		model.setMolecule(this.struc.getChains());
	}
	
	public void addHAtoms() {
		//TODO: Add H atoms to residues
		
	}
	
	public void setModel(Model model) {
		this.model = model;
	}
	
	public Model getModel() {
		return this.model;
	}
	
	//test inputs
	public static void main(String [] args){
		InputController i = new InputController();
		i.createModel();
	    
	}
}
