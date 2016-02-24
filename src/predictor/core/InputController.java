package predictor.core;

import org.biojava.bio.structure.Structure;

import predictor.core.model.Atom;
import predictor.core.model.Chain;
import predictor.core.model.Model;
import predictor.core.model.Molecule;
import predictor.core.model.Residue;
import predictor.core.model.Bond;
import predictor.core.model.math.Geometry;

import predictor.util.PredictorUtility;
import predictor.util.FileReader;

/***
 * 
 * @author Anson
 *
 * Takes input and run initial setup
 */
public class InputController {

	protected double distN_H = 1.0;
	protected Structure struc;
	protected Model model = new Model();
	protected Molecule mol = new Molecule();

	public void getInput(String filePath) {
		struc = FileReader.readFile(filePath);
	}

	public void createModel() {
		this.model.setMolecule(struc.getChains());
		this.mol = model.getMolecules().get(0);
	}

	/**
	 * Add missing hydrogen atoms 
	 */
	public void addHAtoms() {
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
							coordH[k] = coordN[k] + this.distN_H* (coordH[k]-coordN[k]) / lengthN_H;
						}
						newH.setCoordinates(coordH);
						newH.setBond(new Bond(newH, atomN));
						atomN.setBond(new Bond(atomN, newH));
						r.getAtomList().add(atomNPos+1, newH);
					}
				}
			}
		}
	}
	
	public void addPhiPsiAngles() {
		for(int i=0; i<mol.getChains().size(); i++){
			Chain c = mol.getChains().get(i);
			for(int j=0; j<c.getResidues().size(); j++){
				c.getResidues().get(j).setPhi(Geometry.calPhi(c, j));
				c.getResidues().get(j).setPsi(Geometry.calPsi(c, j));
			}
		}
	}

	public void setModel(Model model) {
		this.model = model;
	}

	public Model getModel() {
		return this.model;
	}

	//test inputs
//	public static void main(String [] args){
//		InputController i = new InputController();
//		i.getInput("res/test/4hhb.pdb");
//		i.createModel();
//
//		i.addHAtoms();
//		i.addPhiPsiAngles();
//		Model m2 = i.getModel();
//		Chain c2 = m2.getMolecules().get(0).getChains().get(0);
//		Residue r2 = c2.getResidues().get(3);
//		for(int z=0;z<r2.getAtomList().size();z++){
//			System.out.println(r2.getAtomList().get(z).getSymbol());
//		}
//		System.out.println();
//		for(int j=0; j<c2.getResidues().size();j++){
//			System.out.println(c2.getResidues().get(j).getName());
//			for(int z=0;z<c2.getResidues().get(j).getAtomList().size();z++){
////				if(c2.getResidues().get(j).getAtomList().get(z).getSymbol()=="H"){
//					System.out.println(j+"--------------------------------");
//					System.out.println(c2.getResidues().get(j).getPhi());
//					System.out.println(c2.getResidues().get(j).getPsi());
////					System.out.println(c2.getResidues().get(j).getAtomList().get(z).getCoordinates()[2]);
////					System.out.println(c2.getResidues().get(j).getAtomList().get(z+1).getCoordinates()[2]);
////					System.out.println(c2.getResidues().get(j).getAtomList().get(z+2).getCoordinates()[2]);
//					System.out.println();
////				}
//			}
//		}
//	}
}
