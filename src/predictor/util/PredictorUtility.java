package predictor.util;

import predictor.core.model.Residue;
import predictor.core.model.Molecule;

/***
 * 
 * @author Anson
 *
 *	Utility class
 */
public final class PredictorUtility {
	
	/**
	 * Find an atom in a residue by its PDB name & return position of first atom found
	 * @param r the residue that contains the atom
	 * @param c PDB name of atom
	 * @return position of atom in residue
	 */
	public static int findAtom(Residue r, String c) {
		for(int i=0; i<r.getAtomList().size(); i++) {
			if(c.equals(r.getAtomList().get(i).getSymbol())) {
				return i;
			}
		}
		return -1;
	}
	
	/**
	 * Find chain in molecule by its name & return position of first chain found
	 * @param m molecule
	 * @param c name of chain
	 * @return position of chain in molecule
	 */
	public static int findChain(Molecule m, String c) {
		for(int i=0; i<m.getChains().size(); i++){
			if(c.equals(m.getChains().get(i).getName())){
				return i;
			}
		}
		
		return -1;
	}
	
	/**
	 * Convert 3 characters residue name to 1
	 * @param name residue name
	 * @return single character residue name
	 */
	public static char processRName(String name){
		switch(name) {
		case "ALA": return 'A';	
		case "ARG": return 'R';
		case "ASN": return 'N';
		case "ASP": return 'D';
		case "ASX": return 'B';
		case "CYS": return 'C';
		case "GLN": return 'Q';
		case "GLU": return 'E';
		case "GLX": return 'Z';
		case "GLY": return 'G';
		case "HIS": return 'H';
		case "ILE": return 'I';
		case "LEU": return 'L';
		case "LYS": return 'K';
		case "MET": return 'M';
		case "PRO": return 'P';
		case "PHE": return 'F';
		case "SER": return 'S';
		case "THR": return 'T';
		case "TRP": return 'W';
		case "TYR": return 'Y';
		case "VAL": return 'V';
		default: return 'X';
		}
	}
}
