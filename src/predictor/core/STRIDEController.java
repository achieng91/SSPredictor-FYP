package predictor.core;

import predictor.core.model.Model;
import predictor.core.model.Molecule;
import predictor.core.model.Chain;
import predictor.core.model.Residue;
import predictor.core.model.Atom;
import predictor.core.model.Donor;
import predictor.core.model.Acceptor;
import predictor.core.model.HBond;

import predictor.core.model.math.Vector3D;
import predictor.core.model.math.Geometry;

import predictor.util.PredictorUtility;
import java.util.ArrayList;

/***
 * 
 * @author Anson
 * 
 * Predictor using STRIDE algorithm
 */
public class STRIDEController {
	
	protected final double DISTCUTOFF = 6.0;
	protected final double E_M = -2.8, R_M = 3.0;

	protected Model model;
	protected Molecule mol;
	protected ArrayList<Donor> donors;
	protected ArrayList<Acceptor> acceptors;
	protected ArrayList<HBond> HBond;
	
	public STRIDEController(Model model){
		this.model = model;
		this.mol = model.getMolecules().get(0);
	}
	
	/** 
	 * Run STRIDE prediction
	 */
	public void predict() {
		this.findHBonds();
	}
	
	/**
	 * Determine Hydrogen bonds
	 */
	protected void findHBonds() {
		for(int i=0; i<mol.getChains().size(); i++){
			this.findDonor(mol.getChains().get(i));
			this.findAcceptor(mol.getChains().get(i));
		}
		
		boolean[] bondedDonor = new boolean[donors.size()];
		boolean[] bondedAcceptor = new boolean[acceptors.size()];
		
		for(int i=0; i<donors.size(); i++){
			bondedDonor[i] = false;
		}
		for(int i=0; i<acceptors.size(); i++){
			bondedAcceptor[i] = false;
		}
		
		for(int i=0; i<donors.size(); i++){
			for(int j=0; j<acceptors.size(); j++){
				Acceptor acc = acceptors.get(j);
				Donor don = donors.get(i);
				if(Math.abs(acc.getResidueNum()-don.getResidueNum())<2 && 
						acc.getChain().getName().equals(don.getChain().getName())){
					continue;
				}
				
				HBond hb = new HBond();
				hb.setAccDonDist(Geometry.calDist(don.getDonorAtomN(), acc.getAcceptorAtomO()));
				if(hb.getAccDonDist() < DISTCUTOFF) {
					if(don.getDonorAtomH()!=null){
						hb.setHBondEnergy(this.calHBondEnergy(
								acc.getAcceptorAtomCA(),
								acc.getAcceptorAtomC(),
								acc.getAcceptorAtomO(),
								don.getDonorAtomH(),
								don.getDonorAtomN(),
								hb)
								);
					}
					
					if(hb.getHBondEnergy()<-10.0 && Math.abs(hb.getEt())>0.000001 && Math.abs(hb.getEp())>0.000001){
						hb.setExistPolarInter(true);
					}
					
					hb.setOHDist(Geometry.calDist(don.getDonorAtomH(), acc.getAcceptorAtomO()));
					hb.setAngNHO(Geometry.calAngle(don.getDonorAtomN().getPosition(), don.getDonorAtomH().getPosition(), 
							acc.getAcceptorAtomO().getPosition()));
					hb.setAngCOH(Geometry.calAngle(acc.getAcceptorAtomC().getPosition(), acc.getAcceptorAtomO().getPosition(), 
							don.getDonorAtomH().getPosition()));
					
					if(hb.getOHDist()<=2.5 && hb.getAngNHO()>=90.0 && hb.getAngCOH()>=90.0 && hb.getAngCOH()<180.0){
						hb.setExistHBondBaker(true);
					}
					
					if(hb.getAccDonDist() <= don.getHBondRadius()+acc.getHBondRadius()){
						hb.setAccAng(Geometry.calAngle(don.getDonorAtomN().getPosition(), acc.getAcceptorAtomO().getPosition(), 
								acc.getAcceptorAtomC().getPosition()));
						
						if(hb.getAccAng() >= 90.0 && hb.getAccAng() <= 180.0){
							hb.setDonAng(Geometry.calAngle(acc.getAcceptorAtomO().getPosition(), don.getDonorAtomN().getPosition(), 
									don.getDonorAtomC().getPosition()));
							
							if(hb.getDonAng() >= 90.0 && hb.getDonAng() <= 180.0){
								hb.setAccDonAng(Math.abs(Geometry.calDihedralAngle(
										don.getDonorAtomCA().getPosition(), 
										don.getDonorAtomN().getPosition(),
										don.getDonorAtomC().getPosition(),
										acc.getAcceptorAtomO().getPosition())));
								
								if(hb.getAccDonAng() > 90.0 && hb.getAccDonAng() < 270.0){
									hb.setAccDonAng(Math.abs(180.0 - hb.getAccDonAng()));
								}
							}
							
							hb.setDonAccAng(Geometry.calDihedralAngle(
									don.getDonorAtomN().getPosition(), 
									acc.getAcceptorAtomO().getPosition(), 
									acc.getAcceptorAtomC().getPosition(),
									acc.getAcceptorAtomCA().getPosition()));
							if(hb.getDonAccAng() > 90.0 && hb.getDonAccAng() < 270.0){
								hb.setDonAccAng(Math.abs(180.0 - hb.getDonAccAng()));
							}
							
							if(hb.getAccDonAng() <= 60.0 && hb.getDonAccAng() <= 90.0){
								hb.setExistHBondRose(true);
							}
						}
					}
				}
				
				if((hb.getExistPolarInter() && hb.getHBondEnergy()<0.0) || hb.getExistHBondRose() || hb.getExistHBondBaker()) {
					hb.setDonor(don);
					hb.setAcceptor(acc);
					bondedDonor[i] = true;
					bondedAcceptor[j] = true;
				}
			}
		}
	}
	
	/**
	 * Find all donors residues & atoms in chain
	 * @param c 
	 */
	protected void findDonor(Chain c){
		for(int i=0; i<c.getResidues().size(); i++){
			Residue r = c.getResidues().get(i);
			Donor d;
			if(i==0){
				d = new Donor(
						c, 
						i, 
						r.getAtomList().get(PredictorUtility.findAtom(r, "CA")),
						r.getAtomList().get(PredictorUtility.findAtom(r, "H"))
						);
			} else {
				d = new Donor(
						c, 
						i, 
						r.getAtomList().get(PredictorUtility.findAtom(r, "N")),
						r.getAtomList().get(PredictorUtility.findAtom(c.getResidues().get(i-1), "C")),
						r.getAtomList().get(PredictorUtility.findAtom(r, "CA")),
						r.getAtomList().get(PredictorUtility.findAtom(r, "H"))
						);
			}
			donors.add(d);
		}
	}
	
	/**
	 * Find all acceptors residues & atoms in chain
	 * @param c
	 */
	protected void findAcceptor(Chain c){
		for(int i=0; i<c.getResidues().size(); i++){
			Residue r = c.getResidues().get(i);
			Acceptor a;
			if(i==c.getResidues().size()-1){
				a = new Acceptor(
						c, 
						i, 
						r.getAtomList().get(PredictorUtility.findAtom(r, "CA"))
						);
			} else {
				a = new Acceptor(
						c, 
						i, 
						r.getAtomList().get(PredictorUtility.findAtom(r, "O")),
						r.getAtomList().get(PredictorUtility.findAtom(r, "C")),
						r.getAtomList().get(PredictorUtility.findAtom(r, "CA"))
						);
			}
			acceptors.add(a);
		}
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
	protected double calHBondEnergy(Atom CA, Atom C, Atom O, Atom H, Atom N, HBond hb){
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
