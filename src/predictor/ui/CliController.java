package predictor.ui;

import predictor.core.model.secStructures.SSModel;

public class CliController extends Controller {

	public CliController(String file) {
		super(file);
	}

	public void displayOutput(SSModel model) {
		System.out.println("Generate Output..");
		Predictor.output(model);
		System.out.println("End");
	}
	
	public void generateOutputFile(String fileName, String filePath){
		Predictor.genOutputFile(fileName, filePath);
	}
}
