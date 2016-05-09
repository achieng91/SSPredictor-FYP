package predictor.core;

import java.util.ArrayList;

import predictor.core.model.AsnResidue;
import predictor.core.model.Chain;
import predictor.core.model.HBond;
import predictor.core.model.Residue;
import predictor.core.model.math.PhiPsiMap;

/***
 * Determine Helix structures
 * @author Anson
 *
 */
public class STRIDEHelixController {

	/** 
	 * Check if structure is a helix
	 * @param c
	 * @param hb
	 */
	protected void isHelix(Chain c, ArrayList<HBond> hb){
		double[] prob = new double[c.getResidues().size()];
		
		for(int i=0;i<c.getResidues().size(); i++){
			prob[i] = 0.0;
		}
		
		for(int i=0;i<c.getResidues().size()-5; i++){
			ArrayList<Residue> r = c.getResidues();
			int bondNum = STRIDEController.findPolInt(hb, r.get(i+4), r.get(i));
			if(bondNum != -1){
				prob[i] = hb.get(bondNum).getHBondEnergy() * (0.5*(PhiPsiMap.calProb(r.get(i).getPhi(), r.get(i).getPsi(), "HELIX")+
						PhiPsiMap.calProb(r.get(i+4).getPhi(), r.get(i+4).getPsi(), "HELIX")));
			}
		}
		
		// Check for Alpha Helix
		for(int i=0; i< c.getResidues().size()-5; i++){
			if(Double.compare(prob[i],STRIDEController.THRESHOLD_H1)<0 && Double.compare(prob[i+1],STRIDEController.THRESHOLD_H1)<0){
				ArrayList<Residue> r = c.getResidues();
				
				if(Double.compare(PhiPsiMap.calProb(r.get(i).getPhi(), r.get(i).getPsi(), "HELIX"),STRIDEController.THRESHOLD_H3)>0){
					AsnResidue h0 = new AsnResidue(r.get(i), "ALPHAHELIX", "H");
					r.add(i, h0); r.remove(i+1);
				}
				if(Double.compare(PhiPsiMap.calProb(r.get(i+5).getPhi(), r.get(i+5).getPsi(), "HELIX"),STRIDEController.THRESHOLD_H4)>0){
					AsnResidue h5 = new AsnResidue(r.get(i+5), "ALPHAHELIX", "H");
					r.add(i+5, h5); r.remove(i+6);
				}
				

				AsnResidue h1 = new AsnResidue(r.get(i+1), "ALPHAHELIX", "H");
				r.add(i+1, h1); r.remove(i+2);
				AsnResidue h2 = new AsnResidue(r.get(i+2), "ALPHAHELIX", "H");
				r.add(i+2, h2); r.remove(i+3);
				AsnResidue h3 = new AsnResidue(r.get(i+3), "ALPHAHELIX", "H");
				r.add(i+3, h3); r.remove(i+4);
				AsnResidue h4 = new AsnResidue(r.get(i+4), "ALPHAHELIX", "H");
				r.add(i+4, h4); r.remove(i+5);
			}
		}
		
		// Check for 310 Helix
		for(int i=0; i<c.getResidues().size()-4; i++){
			ArrayList<Residue> r = c.getResidues();
			if(STRIDEController.findBnd(hb, r.get(i+3), r.get(i))!=-1 && STRIDEController.findBnd(hb, r.get(i+4), r.get(i+1))!=-1 &&
					((!r.get(i+1).getAsn().equals("H") && !r.get(i+2).getAsn().equals("H")) || 
							(!r.get(i+2).getAsn().equals("H") && !r.get(i+3).getAsn().equals("H")))){
				AsnResidue h1 = new AsnResidue (r.get(i+1), "310HELIX", "G");
				r.add(i+1, h1); r.remove(i+2);
				AsnResidue h2 = new AsnResidue (r.get(i+2), "310HELIX", "G");
				r.add(i+2, h2); r.remove(i+3);
				AsnResidue h3 = new AsnResidue (r.get(i+3), "310HELIX", "G");
				r.add(i+3, h3); r.remove(i+4);
			}
			
		}
		
		//Check for Pi Helix
		for(int i=0; i<c.getResidues().size()-6; i++){
			ArrayList<Residue> r = c.getResidues();
			if(STRIDEController.findBnd(hb, r.get(i+5), r.get(i))!=-1 && STRIDEController.findBnd(hb, r.get(i+6), r.get(i+1))!=-1 &&
					r.get(i+1).getAsn()=="C" && r.get(i+2).getAsn()=="C" && r.get(i+3).getAsn()=="C" &&
					r.get(i+4).getAsn()=="N" && r.get(i+5).getAsn()=="N"){
				AsnResidue h1 = new AsnResidue (r.get(i+1), "PIHELIX", "I");
				r.add(i+1, h1); r.remove(i+2);
				AsnResidue h2 = new AsnResidue (r.get(i+2), "PIHELIX", "I");
				r.add(i+2, h2); r.remove(i+3);
				AsnResidue h3 = new AsnResidue (r.get(i+3), "PIHELIX", "I");
				r.add(i+3, h3); r.remove(i+4);
				AsnResidue h4 = new AsnResidue (r.get(i+4), "PIHELIX", "I");
				r.add(i+4, h4); r.remove(i+5);
				AsnResidue h5 = new AsnResidue (r.get(i+5), "PIHELIX", "I");
				r.add(i+5, h5); r.remove(i+6);
			}
		}
	}
}
