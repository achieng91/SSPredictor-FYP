package predictor.util;

import predictor.core.model.Residue;

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
}
