package predictor.core;

import java.util.ArrayList;

import predictor.core.model.AsnResidue;
import predictor.core.model.Chain;
import predictor.core.model.HBond;
import predictor.core.model.Residue;
import predictor.core.model.math.Geometry;
import predictor.util.PredictorUtility;

/***
 * Determine Turn structures
 * @author Anson
 *
 */
public class STRIDETurnController {
	
	/**
	 * Check for beta turn
	 * @param c current chain
	 */
	protected void isBetaTurn(Chain c){
		int posCA1, posCA4;
		double phi2, phi3, psi2, psi3, range1 = 30.0, range2 = 45.0;
		String turnType;
		
		for(int i=0; i<c.getResidues().size()-4; i++){
			ArrayList<Residue> r = c.getResidues();
			Residue r0 = c.getResidues().get(i);
			Residue r1 = c.getResidues().get(i+1);
			Residue r2 = c.getResidues().get(i+2);
			Residue r3 = c.getResidues().get(i+3);
			
			posCA1 = PredictorUtility.findAtom(r0, "CA");
			posCA4 = PredictorUtility.findAtom(r3, "CA");
			
			double tmpDist = 0;
			if(posCA1!=-1 && posCA4!=-1){
				tmpDist = Geometry.calDist(r0.getAtomList().get(posCA1), r3.getAtomList().get(posCA4));
			}
			
			if(r1.getAsn().matches("H|G|I") || r2.getAsn().matches("H|G|I") || 
					posCA1==-1 || posCA4==-1 || tmpDist > 7.0){
				continue;
			}
			
			phi2 = r1.getPhi();
			psi2 = r1.getPsi();
			phi3 = r2.getPhi();
			psi3 = r2.getPsi();
			
			// Assign type of turns
			if(turnCondition(phi2, -60.0, psi2, -30, phi3, -90.0, psi3, 0, range1, range2)){
				turnType = "TurnI";
			} else if(turnCondition(phi2, 60.0, psi2, 30, phi3, 90.0, psi3, 0, range1, range2)){
				turnType = "TurnI'";
			} else if(turnCondition(phi2, -60.0, psi2, 120, phi3, 80.0, psi3, 0, range1, range2)){
				turnType = "TurnII";
			} else if(turnCondition(phi2, 60.0, psi2,-120, phi3, -80.0, psi3, 0, range1, range2)){
				turnType = "TurnII'";
			} else if(turnCondition(phi2, -60.0, psi2,120, phi3, -90.0, psi3, 0, range1, range2)){
				turnType = "TurnVIa1";
			} else if(turnCondition(phi2, -120.0, psi2,120, phi3, -60.0, psi3, 0, range1, range2)){
				turnType = "TurnVIa2";
			} else if(turnCondition(phi2, -60.0, psi2,-30, phi3, -120.0, psi3, 120, range1, range2)){
				turnType = "TurnVIII";
			} else 
				turnType = "TurnIV";
			
			if(r0.getAsn().equals("C")){
				AsnResidue t0 = new AsnResidue(r0, "TURN", "T");
				t0.setTurnType(turnType);
				r.add(i, t0); r.remove(i+1);
			}
			if(r1.getAsn().equals("C")){
				AsnResidue t1 = new AsnResidue(r1, "TURN", "T");
				t1.setTurnType(turnType);
				r.add(i+1, t1); r.remove(i+2);
			}
			if(r2.getAsn().equals("C")){
				AsnResidue t2 = new AsnResidue(r2, "TURN", "T");
				t2.setTurnType(turnType);
				r.add(i+2, t2); r.remove(i+3);
			}
			if(r3.getAsn().equals("C")){
				AsnResidue t3 = new AsnResidue(r3, "TURN", "T");
				t3.setTurnType(turnType);
				r.add(i+3, t3); r.remove(i+4);
			}
		}
	}
	
	/**
	 * Condition for different types of beta turn
	 * @param phi2
	 * @param phi2S
	 * @param psi2
	 * @param psi2S
	 * @param phi3
	 * @param phi3S
	 * @param psi3
	 * @param psi3S
	 * @param range1
	 * @param range2
	 * @return boolean if there is turn
	 */
	protected boolean turnCondition(double phi2, double phi2S, double psi2, double psi2S, double phi3, double phi3S, double psi3, double psi3S,
			double range1, double range2){
		if((Math.abs(phi2-phi2S)<=range2 && Math.abs(psi2-psi2S)<=range1 && Math.abs(phi3-phi3S)<=range1 && Math.abs(psi3-psi3S)<=range1) ||
		   (Math.abs(phi2-phi2S)<=range1 && Math.abs(psi2-psi2S)<=range2 && Math.abs(phi3-phi3S)<=range1 && Math.abs(psi3-psi3S)<=range1) ||
		   (Math.abs(phi2-phi2S)<=range1 && Math.abs(psi2-psi2S)<=range1 && Math.abs(phi3-phi3S)<=range2 && Math.abs(psi3-psi3S)<=range1) ||
		   (Math.abs(phi2-phi2S)<=range1 && Math.abs(psi2-psi2S)<=range1 && Math.abs(phi3-phi3S)<=range1 && Math.abs(psi3-psi3S)<=range2)) {
			return true;
		}
		
		return false;
	}
	
	/** 
	 * Check for Gamma Turn
	 * @param c
	 * @param hb
	 */
	protected void isGammaTurn(Chain c, ArrayList<HBond> hb){
		double phi2, psi2;
		String turnType;
		
		for(int i=0; i<c.getResidues().size()-2; i++){
			ArrayList<Residue> r = c.getResidues();
			Residue r0 = null;
			Residue r1 = c.getResidues().get(i);
			Residue r2 = c.getResidues().get(i+1);
			Residue r3 = c.getResidues().get(i+2);
			Residue r4 = null;
			try {
				r0 = c.getResidues().get(i-1);
			} catch(Exception IndexOutOfBoundsException){
				r0 = null;
			}
			
			try {
				r4 = c.getResidues().get(i+3);
			} catch(Exception IndexOutOfBoundsException){
				r4 = null;
			}

			
			if(r2.getAsn().matches("H|T|G|I")){
				continue;
			}
			if(r0 != null){
				if(STRIDEController.findBnd(hb, r3, r0)==-1 || (i>0 && STRIDEController.findBnd(hb, r3, r0)!=-1)){
					continue;
				}

			}
			if(r4 != null){
				if((i<c.getResidues().size()-3 && STRIDEController.findBnd(hb, r4, r1)!=-1))
					continue;
			}
			
			phi2 = r2.getPhi();
			psi2 = r2.getPsi();
			
			if(phi2>0.0 && psi2<0.0){
				turnType = "GammaClassic";
			} else if(phi2<0.0 && psi2>0.0){
				turnType = "GammaInv";
			} else {
				continue;
			}
			
			if(r1.getAsn().equals("C")){
				AsnResidue t1 = new AsnResidue(r1, "TURN", "T");
				t1.setTurnType(turnType);
				r.add(i, t1); r.remove(i+1);
			}
			if(r2.getAsn().equals("C")){
				AsnResidue t2 = new AsnResidue(r2, "TURN", "T");
				t2.setTurnType(turnType);
				r.add(i+1, t2); r.remove(i+2);
			}
			if(r3.getAsn().equals("C")){
				AsnResidue t3 = new AsnResidue(r3, "TURN", "T");
				t3.setTurnType(turnType);
				r.add(i+2, t3); r.remove(i+3);
			}
		}
	}
}
