package predictor.core;

import org.biojava.bio.structure.Structure;

import predictor.core.model.Model;

import predictor.util.FileReader;

/***
 * 
 * @author Anson
 *
 * Takes input and run initial setup
 */
public class InputController {

	protected Structure struc;
	protected Model model;

	public InputController(Model model){
		this.model = model;
	}
	
	public void getInput(String filePath) {
		struc = FileReader.readFile(filePath);
		model.setMolecule(struc.getChains());
	}

	public void setModel(Model model) {
		this.model = model;
	}

	public Model getModel() {
		return this.model;
	}
	
	public Structure getStruc(){
		return this.struc;
	}
}
