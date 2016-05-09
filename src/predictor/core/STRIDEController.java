package predictor.core;

import predictor.core.model.Chain;
import predictor.core.model.Residue;
import predictor.core.model.math.Geometry;
import predictor.core.model.math.Vector3D;
import predictor.core.model.Atom;
import predictor.core.model.HBond;

import java.util.ArrayList;

/***
 * 
 * @author Anson
 * 
 * Predictor using STRIDE algorithm
 */
public interface STRIDEController {
	
	public final double DISTCUTOFF = 6.0;
	public final double E_M = -2.8, R_M = 3.0;
	public final double C1_E = -0.2, C2_E = 0.2;
	public final double THRESHOLD_H1 = -230.0;
	public final double THRESHOLD_H3 = 0.12;
	public final double THRESHOLD_H4 = 0.06;
	public final double THRESHOLD_E1 = -240.0;
	public final double THRESHOLD_E2 = -310.0;
	
	public void predict();
	public void findHBonds();
	public void findDonor(Chain c);
	public void findAcceptor(Chain c);
	public int noDoubleHBond();
	public void discrPhiPsi();
	
	public static int findBnd(ArrayList<HBond> hb, Residue r1, Residue r2){
		int h;
		if(r1.getNBondDnr()!=0 && r2.getNBondAcc()!=0){
			for(int i=0; i<r1.getNBondDnr(); i++){
				h = r1.getHBondDnr()[i];
				for(int j=0; j<r2.getNBondAcc(); j++){
					if(h == r2.getHBondAcc()[j] && hb.get(h).getExistHBondRose()){
						return h;
					}
				}
			}
		}
		
		return -1;
	}
	
	public static int findBnd310(ArrayList<HBond> hb, Residue r1, Residue r2){
		int h;
		if(r1.getNBondDnr()!=0 && r2.getNBondAcc()!=0){
			for(int i=0; i<r1.getNBondDnr(); i++){
				h = r1.getHBondDnr()[i];
				for(int j=0; j<r2.getNBondAcc(); j++){
					if(h == r2.getHBondAcc()[j]){
						return h;
					}
				}
			}
		}
		
		return -1;
	}
	
	public static int findPolInt(ArrayList<HBond> hb, Residue r1, Residue r2){
		int h;
		
		if(r1.getNBondDnr()!=0 && r2.getNBondAcc()!=0){
			for(int i=0; i<r1.getNBondDnr(); i++){
				h = r1.getHBondDnr()[i];
				for(int j=0; j<r2.getNBondAcc(); j++){
					if(h == r2.getHBondAcc()[j] && hb.get(h).getExistPolarInter()){
						return h;
					}
				}
			}
		}
		
		return -1;
	}
	/**
	 * Calculate Hydrogen Bond Energy
	 * @param CA
	 * @param C
	 * @param O
	 * @param H
	 * @param N
	 * @param hb
	 * @return H Bond Energy
	 */
	public static double calHBondEnergy(Atom CA, Atom C, Atom O, Atom H, Atom N, HBond hb){
		if(hb.getAccDonDist() < R_M) {
			hb.setAccDonDist(R_M);
		}
		hb.setEr((-3.0*E_M*Math.pow(R_M, 8.0)/Math.pow(hb.getAccDonDist(), 8.0)) -
				(-4.0*E_M*Math.pow(R_M, 6.0)/Math.pow(hb.getAccDonDist(), 6.0)));
		
		Vector3D projH = Geometry.calProj(O, C, CA, H);
		
		hb.setTi(Math.abs(180.0 - Geometry.calAngle(projH, O.getPosition(), C.getPosition())));
		hb.setTo(Geometry.calAngle(H.getPosition(), O.getPosition(), projH));
		hb.setP(Geometry.calAngle(N.getPosition(), H.getPosition(), O.getPosition()));

		if(hb.getTi() >= 0.0 && hb.getTi() < 90.0){
			hb.setEt(Math.cos(Math.toRadians(hb.getTo()) * (0.9+0.1*Math.sin(Math.toRadians(2*hb.getTi())))));
		} else if(hb.getTi() >= 90.0 && hb.getTi() < 110.0){
			hb.setEt((0.9/Math.pow(Math.cos(Math.toRadians(110.0)), 6.0)) * Math.cos(Math.toRadians(hb.getTo())) *
					(Math.pow(Math.pow(Math.cos(Math.toRadians(110.0)), 2.0)-Math.pow(Math.cos(Math.toRadians(hb.getTi())), 2.0), 3.0)) );
		} else {
			hb.setEt(0.0);
		}
		
		if(hb.getP() > 90.0 && hb.getP() < 270.0){
			hb.setEp(Math.pow(Math.cos(Math.toRadians(hb.getP())), 2.0));
		} else {
			hb.setEp(0.0);
		}
		
		hb.setHBondEnergy(1000.0 * hb.getEr() * hb.getEt() * hb.getEp());
		
		return hb.getHBondEnergy();
	}
}
