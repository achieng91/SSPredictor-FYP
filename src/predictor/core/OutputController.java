package predictor.core;

import predictor.util.FileWriter;
import predictor.util.SSFile;

import predictor.core.model.Model;
import predictor.core.model.secStructures.*;

/***
 * 
 * @author Anson
 *
 * Generates output file and output object
 */
public class OutputController {
	
	protected Model model;
	
	public OutputController (Model model) {
		this.model = model;
	}
	
	/**
	 * Create Secondary Structure object
	 * @return Secondary structure model
	 */
	public SSModel createSSObject() {
		//TODO: create SS Object
		return new SSModel();
	}
	
	/**
	 * Generate output file
	 * @param fileType Type of file to generate
	 * @return Generated file based on chosen file type
	 */
	public SSFile createSSFile(String fileType) {
		return FileWriter.createFile(fileType);
	}
}
