package predictor.core;

import predictor.core.model.Atom;
import predictor.core.model.Bond;
import predictor.core.model.Chain;
import predictor.core.model.Molecule;
import predictor.core.model.Residue;
import predictor.core.model.math.Geometry;
import predictor.util.PredictorUtility;

public abstract class SetupController {

	/**
	 * Add missing hydrogen atoms 
	 */
	public static void addHAtoms(Molecule mol) {
		for(int i=0; i<mol.getChains().size(); i++){
			Chain c = mol.getChains().get(i);
			for(int j=1; j<c.getResidues().size(); j++){
				Residue r = c.getResidues().get(j);

				if(!r.getName().equals("PRO")) {
					int atomDPos = PredictorUtility.findAtom(r, "D");	// Replaces deiterium atoms with hydrogen
					if(atomDPos!=-1) {
						r.getAtomList().get(atomDPos).setSymbol("H");
					}
					
					int atomHPos = PredictorUtility.findAtom(r, "H");
					int atomNPos = PredictorUtility.findAtom(r, "N");
					int atomCPos = PredictorUtility.findAtom(c.getResidues().get(j-1), "C");
					int atomCAPos = PredictorUtility.findAtom(r, "CA");
					if(atomHPos==-1 && atomNPos!=-1 && atomCPos!=-1 && atomCAPos!=-1){
						
						Atom atomN = r.getAtomList().get(atomNPos);
						Atom atomC = c.getResidues().get(j-1).getAtomList().get(atomCPos);
						Atom atomCA = r.getAtomList().get(atomCAPos);
						Atom newH = new Atom();
						double lengthN_C = Geometry.calDist(atomN, atomC);
						double lengthN_CA = Geometry.calDist(atomN, atomCA);
						float[] coordN = atomN.getCoordinates();
						float[] coordC = atomC.getCoordinates();
						float[] coordCA = atomCA.getCoordinates();
						double[] coordH = {0,0,0};
						
						newH.setSymbol("H");

						for(int k=0; k<3; k++) {
							coordH[k] = coordN[k] - ((coordC[k]-coordN[k])/lengthN_C + (coordCA[k]-coordN[k])/lengthN_CA);
						}
						newH.setCoordinates(coordH);
						double lengthN_H = Geometry.calDist(atomN, newH);
						
						for(int k=0; k<3; k++) {
							coordH[k] = coordN[k] + 1.0 * (coordH[k]-coordN[k]) / lengthN_H;
						}
						newH.setCoordinates(coordH);
						newH.setVectorCoord();
						newH.setBond(new Bond(newH, atomN));
						atomN.setBond(new Bond(atomN, newH));
						r.getAtomList().add(atomNPos+1, newH);
					}
				}
			}
		}
	}
	
	public static void addPhiPsiAngles(Molecule mol) {
		for(int i=0; i<mol.getChains().size(); i++){
			Chain c = mol.getChains().get(i);
			for(int j=0; j<c.getResidues().size(); j++){
				c.getResidues().get(j).setPhi(Geometry.calPhi(c, j));
				c.getResidues().get(j).setPsi(Geometry.calPsi(c, j));
			}
		}
	}
}
