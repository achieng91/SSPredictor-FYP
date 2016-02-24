package predictor.ui;

public class CliController extends Controller {

	public CliController(String file) {
		super(file);
	}

	public void displayOutput(String fileType) {
		Predictor.genOutputFile(fileType);
	}
}
