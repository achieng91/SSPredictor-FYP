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
}
