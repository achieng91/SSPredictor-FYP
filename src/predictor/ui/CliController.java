package predictor.ui;

public class CliController extends Controller {

	public CliController(String file) {
		super(file);
	}

	public void displayOutput() {
		System.out.println("Generate Output..");
		Predictor.output();
		Predictor.genOutputFile("4hhb", "");
		System.out.println("End");
	}
}
