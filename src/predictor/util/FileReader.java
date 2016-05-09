package predictor.util;

import java.util.logging.Level;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.io.*;

public abstract class FileReader {
	
	public static Structure readFile(String fileName){
		if(fileName.endsWith(".pdb")) {
			return FileReader.readPDBFile(fileName);
		} else if(fileName.endsWith(".dssp")){
			return FileReader.readDSSPFile(fileName);
		}
			
		return null;
	}
	
	public static Structure readPDBFile (String fileName) {
		// wrapper class for parsing a PDB file.
		PDBFileReader pdbreader = new PDBFileReader();
		Structure struc = null; 

		try{
			// Access to the data of a PDB file.
			struc = pdbreader.getStructure(fileName);	    
		} catch (Exception e){
			java.util.logging.Logger logger = java.util.logging.Logger.getLogger("FileReader");
			logger.log(Level.SEVERE, "Error converting to structure object from file path");
		}
		return struc;
	}
	
	public static Structure readDSSPFile(String fileName) {
		//TODO: Read dssp file
		
		return null;
	}
}
