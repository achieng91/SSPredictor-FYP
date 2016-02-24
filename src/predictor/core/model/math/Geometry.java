package predictor.core.model.math;

import predictor.core.model.Atom;
import predictor.core.model.Chain;
import predictor.core.model.Residue;
import predictor.core.model.math.Vector3D;

import predictor.util.PredictorUtility;

/***
 * 
 * @author Anson
 *
 * Geometry class that include various calculations
 */
public abstract class Geometry {
	
	/** 
	 * Calculate dihedral angle between 4 atoms
	 * @param p1 atom 1
	 * @param p2 atom 2
	 * @param p3 atom 3
	 * @param p4 atom 4
	 * @return calculated dihedral angle
	 */
	public static double calDihedralAngle(Vector3D p1, Vector3D p2, Vector3D p3, Vector3D p4) {
		Vector3D v_12 = p1.subtractAndReturn(p2);
		Vector3D v_43 = p4.subtractAndReturn(p3);
		Vector3D v_23 = p2.subtractAndReturn(p3);
		
		Vector3D p = v_23.getCrossProduct(v_12);
		Vector3D x = v_23.getCrossProduct(v_43);
		Vector3D y = v_23.getCrossProduct(x);
		
		double u = x.getScalarProduct(x);
		double v = y.getScalarProduct(y);
		double result = 360.0;
		
		if(u > 0 && v >0) {
			u = p.getScalarProduct(x) / Math.sqrt(u);
			v = p.getScalarProduct(y) / Math.sqrt(v);
			if(u !=0 || v != 0) {
				result = Math.atan2(v, u) * 180/Math.PI;
			}
		}
		return result;
	}
	
	/**
	 * Calculate angle between 3 vectors
	 * @param vA
	 * @param vB
	 * @param vC
	 * @return angle
	 */
	public static double calAngle(Vector3D vA, Vector3D vB, Vector3D vC){
		Vector3D vBA = vA.subtractAndReturn(vB);
		Vector3D vBC = vC.subtractAndReturn(vB);
		
		double A = vBA.getScalarProduct(vBC);
		double lengthBA = vBA.getMagnitude();
		double lengthBC = vBC.getMagnitude();
		double D = A/(lengthBA*lengthBC);
		
		if(D > 0.0 && Math.abs(D-1.0) < 0.000001){
			D -= 0.000001;
		} else if(D < 0.0 && Math.abs(D+1.0) < 0.000001){
			D += 0.000001;
		}
		return 57.2958 * Math.acos(D);
	}
	
	/**
	 * Calculate Phi angle 
	 * @param chain The chain that contains the Phi angle
	 * @param res The residue number in the chain
	 * @return Calculated Phi angle
	 */
	public static double calPhi(Chain chain, int res) {
		if(res==0)
			return 360.0;

		Residue r = chain.getResidues().get(res);
		Residue pr = chain.getResidues().get(res-1);
		
		Vector3D C_prev = pr.getAtomList().get(PredictorUtility.findAtom(pr, "C")).getPosition();
		Vector3D N = r.getAtomList().get(PredictorUtility.findAtom(r, "N")).getPosition();
		Vector3D CA = r.getAtomList().get(PredictorUtility.findAtom(r, "CA")).getPosition();
		Vector3D C = r.getAtomList().get(PredictorUtility.findAtom(r,"C")).getPosition();
		
		return calDihedralAngle(C_prev, N, CA, C);
	}
	
	/**
	 * Calculate Psi angle
	 * @param chain The chain that contains the Psi angle
	 * @param res The residue number in the chain
	 * @return Calculated Psi angle
	 */
	public static double calPsi(Chain chain, int res) {
		if(res == chain.getResidues().size()-1) {
			return 360.0;
		}
		
		Residue r = chain.getResidues().get(res);
		Residue nr = chain.getResidues().get(res+1);
		
		Vector3D N = r.getAtomList().get(PredictorUtility.findAtom(r, "N")).getPosition();
		Vector3D CA = r.getAtomList().get(PredictorUtility.findAtom(r, "CA")).getPosition();
		Vector3D C = r.getAtomList().get(PredictorUtility.findAtom(r, "C")).getPosition();
		Vector3D N_next = nr.getAtomList().get(PredictorUtility.findAtom(nr, "N")).getPosition();
		
		return calDihedralAngle(N, CA, C, N_next);
	}
	
	/**
	 * Calculate distance between 2 atoms
	 * @param a first atom
	 * @param b second atom
	 * @return length between atom a & b
	 */
	public static double calDist(Atom a, Atom b) {
		Vector3D vA = a.getPosition();
		Vector3D vB = b.getPosition();
		
		return vA.subtractAndReturn(vB).getMagnitude();
	}
	
	/**
	 * Calculate projection of atom on plane
	 * @param a
	 * @param b
	 * @param c
	 * @param d
	 * @return Coordinates of atom d projection
	 */
	public static Vector3D calProj(Atom a, Atom b, Atom c, Atom d){
		Vector3D vA = a.getPosition();
		Vector3D vB = b.getPosition();
		Vector3D vC = c.getPosition();
		Vector3D vD = d.getPosition();
		
		Vector3D vBA = vA.subtractAndReturn(vB);
		Vector3D vBC = vC.subtractAndReturn(vB);
		Vector3D vAD = vD.subtractAndReturn(vA);
		double lengthAD = vAD.getMagnitude();
		
		Vector3D vNorm = vBA.getCrossProduct(vBC);
		vNorm.normalize();
		
		double COSNormAD = vNorm.getScalarProduct(vAD) / lengthAD;
		if(COSNormAD > 0.0) {
			COSNormAD = Math.abs(COSNormAD);
		}
		double projADNorm = lengthAD * COSNormAD;
		vNorm.invert();
		vNorm.scale(projADNorm);
		
		Vector3D result = (vAD.subtractAndReturn(vNorm)).addAndReturn(vA);
		return result;
	}
}
