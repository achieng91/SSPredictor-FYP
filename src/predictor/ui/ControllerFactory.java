package predictor.ui;

public class ControllerFactory {
	public static Controller createController(boolean enableUI, String file) {
		if (enableUI) {
			//TODO: Add GUI Controller
//			return new GuiController();
			return new CliController(file);
		}
		else {
			return new CliController(file);
		}
	}
}
