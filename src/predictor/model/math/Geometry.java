package predictor.model.math;

import predictor.model.Chain;
import predictor.model.Residue;
import predictor.model.math.Vector3D;

public class Geometry {
	
	/** 
	 * Calculate dihedral angle between 4 atoms
	 * @param p1 atom 1
	 * @param p2 atom 2
	 * @param p3 atom 3
	 * @param p4 atom 4
	 * @return calculated dihedral angle
	 */
	public double calDihedralAngle(Vector3D p1, Vector3D p2, Vector3D p3, Vector3D p4) {
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
	 * Calculate Phi angle 
	 * @param chain The chain that contains the Phi angle
	 * @param res The residue number in the chain
	 * @return Calculated Phi angle
	 */
	public double calPhi(Chain chain, int res) {
		if(res==0)
			return 360.0;
		
		Residue r = chain.getResidues().get(res);
		Residue pr = chain.getResidues().get(res-1);
		
		Vector3D C_prev = pr.getAtomList().get(this.findAtom(pr, "C")).getPosition();
		Vector3D N = r.getAtomList().get(this.findAtom(r, "N")).getPosition();
		Vector3D CA = r.getAtomList().get(this.findAtom(r, "CA")).getPosition();
		Vector3D C = r.getAtomList().get(this.findAtom(r,"C")).getPosition();
		
		return this.calDihedralAngle(C_prev, N, CA, C);
	}
	
	/**
	 * Calculate Psi angle
	 * @param chain The chain that contains the Psi angle
	 * @param res The residue number in the chain
	 * @return Calculated Psi angle
	 */
	public double calPsi(Chain chain, int res) {
		Residue r = chain.getResidues().get(res);
		Residue nr = chain.getResidues().get(res+1);
		
		Vector3D N = r.getAtomList().get(this.findAtom(r, "N")).getPosition();
		Vector3D CA = r.getAtomList().get(this.findAtom(r, "CA")).getPosition();
		Vector3D C = r.getAtomList().get(this.findAtom(r, "C")).getPosition();
		Vector3D N_next = nr.getAtomList().get(this.findAtom(nr, "N")).getPosition();
		
		return this.calDihedralAngle(N, CA, C, N_next);
	}
	
	/**
	 * Find an atom in a residue by its PDB name
	 * @param r the residue that contains the atom
	 * @param c PDB name of atom
	 * @return position of atom in residue
	 */
	public int findAtom(Residue r, String c) {
		for(int i=0; i<r.getAtomList().size(); i++) {
			if(c.equals(r.getAtomList().get(i).getSymbol())) {
				return i;
			}
		}
		return -1;
	}

}
