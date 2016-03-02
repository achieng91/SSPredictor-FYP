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
import predictor.core.model.math.RamaMap;
import predictor.core.model.math.Geometry;

import predictor.core.model.secStructures.*;

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
	protected final double THRESHOLD_H1 = -230.0;
	protected final double THRESHOLD_H3 = 0.12;
	protected final double THRESHOLD_H4 = 0.06;

	protected Model model;
	protected Molecule mol;
	protected ArrayList<Donor> donors = new ArrayList<Donor>();
	protected ArrayList<Acceptor> acceptors = new ArrayList<Acceptor>();
	protected ArrayList<HBond> hBond = new ArrayList<HBond>();
	
	public STRIDEController(Model model){
		this.model = model;
		this.mol = model.getMolecules().get(0);
	}
	
	/** 
	 * Run STRIDE prediction
	 */
	public void predict() {
		this.findHBonds();
		
		for(int i=0; i<mol.getChains().size(); i++){
			isHelix(mol.getChains().get(i), hBond);
			for(int j=0; j<mol.getChains().size(); j++){
				isBetaTurn(mol.getChains().get(i));
//				isGammaTurn(mol.getChains().get(i), hBond);
			}
		}
	}
	
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
			int bondNum = findPolInt(hb, r.get(i+4), r.get(i));
			if(bondNum != -1){
				prob[i] = hb.get(bondNum).getHBondEnergy() * (0.5*(RamaMap.calHelixProb(r.get(i).getPhi(), r.get(i).getPsi())+
						RamaMap.calHelixProb(r.get(i+4).getPhi(), r.get(i+4).getPsi())));
			}
		}
		
		// Check for Alpha Helix
		for(int i=0; i< c.getResidues().size()-5; i++){
			if(prob[i] < THRESHOLD_H1 && prob[i+1] < THRESHOLD_H1){
				ArrayList<Residue> r = c.getResidues();
				
				HelixAlpha h1 = new HelixAlpha(r.get(i+1));
				r.add(i+1, h1); r.remove(i+2);
				HelixAlpha h2 = new HelixAlpha(r.get(i+2));
				r.add(i+2, h2); r.remove(i+3);
				HelixAlpha h3 = new HelixAlpha(r.get(i+3));
				r.add(i+3, h3); r.remove(i+4);
				HelixAlpha h4 = new HelixAlpha(r.get(i+4));
				r.add(i+4, h4); r.remove(i+5);
				
				if(RamaMap.calHelixProb(r.get(i).getPhi(), r.get(i).getPsi()) > THRESHOLD_H3){
					HelixAlpha h0 = new HelixAlpha(r.get(i));
					r.add(i, h0); r.remove(i+1);
				}
				if(RamaMap.calHelixProb(r.get(i+5).getPhi(), r.get(i+5).getPsi()) > THRESHOLD_H4){
					HelixAlpha h5 = new HelixAlpha(r.get(i+5));
					r.add(i+5, h5); r.remove(i+6);
				}
			}
		}
		
		// Check for 310 Helix
		for(int i=0; i<c.getResidues().size()-4; i++){
			ArrayList<Residue> r = c.getResidues();
			if(findBnd(hb, r.get(i+3), r.get(i))!=-1 && findBnd(hb, r.get(i+4), r.get(i+1))!=-1 &&
					(!r.get(i+1).getAsn().equals("H") && !r.get(i+2).getAsn().equals("H") && 
							!r.get(i+3).getAsn().equals("H") && !r.get(i+4).getAsn().equals("H"))){
				HelixThreeTen h1 = new HelixThreeTen (r.get(i+1));
				r.add(i+1, h1); r.remove(i+2);
				HelixThreeTen h2 = new HelixThreeTen (r.get(i+2));
				r.add(i+2, h2); r.remove(i+3);
				HelixThreeTen h3 = new HelixThreeTen (r.get(i+3));
				r.add(i+3, h3); r.remove(i+4);
			}
			
		}
		
		//Check for Pi Helix
		for(int i=0; i<c.getResidues().size()-6; i++){
			ArrayList<Residue> r = c.getResidues();
			if(findBnd(hb, r.get(i+5), r.get(i))!=-1 && findBnd(hb, r.get(i+6), r.get(i+1))!=-1 &&
					r.get(i+1).getAsn()=="C" && r.get(i+2).getAsn()=="C" && r.get(i+3).getAsn()=="C" &&
					r.get(i+4).getAsn()=="N" && r.get(i+5).getAsn()=="N"){
				HelixPi h1 = new HelixPi (r.get(i+1));
				r.add(i+1, h1); r.remove(i+2);
				HelixPi h2 = new HelixPi (r.get(i+2));
				r.add(i+2, h2); r.remove(i+3);
				HelixPi h3 = new HelixPi (r.get(i+3));
				r.add(i+3, h3); r.remove(i+4);
				HelixPi h4 = new HelixPi (r.get(i+4));
				r.add(i+4, h4); r.remove(i+5);
				HelixPi h5 = new HelixPi (r.get(i+5));
				r.add(i+5, h5); r.remove(i+6);
			}
		}
	}
	
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
				turnType = "TurnII";
			} else if(turnCondition(phi2, -60.0, psi2, 120, phi3, 80.0, psi3, 0, range1, range2)){
				turnType = "TurnIII";
			} else if(turnCondition(phi2, 60.0, psi2,-120, phi3, -80.0, psi3, 0, range1, range2)){
				turnType = "TurnIV";
			} else if(turnCondition(phi2, -60.0, psi2,120, phi3, -90.0, psi3, 0, range1, range2)){
				turnType = "TurnV";
			} else if(turnCondition(phi2, -120.0, psi2,120, phi3, -60.0, psi3, 0, range1, range2)){
				turnType = "TurnVI";
			} else if(turnCondition(phi2, -60.0, psi2,-30, phi3, -120.0, psi3, 120, range1, range2)){
				turnType = "TurnVII";
			} else 
				turnType = "TurnVIII";
			
			if(r0.getAsn().equals("C")){
				Turn t0 = new Turn(r0);
				t0.setTurnType(turnType);
				r.add(i, t0); r.remove(i+1);
			}
			if(r1.getAsn().equals("C")){
				Turn t1 = new Turn(r1);
				t1.setTurnType(turnType);
				r.add(i+1, t1); r.remove(i+2);
			}
			if(r2.getAsn().equals("C")){
				Turn t2 = new Turn(r2);
				t2.setTurnType(turnType);
				r.add(i+2, t2); r.remove(i+3);
			}
			if(r3.getAsn().equals("C")){
				Turn t3 = new Turn(r3);
				t3.setTurnType(turnType);
				r.add(i+3, t3); r.remove(i+4);
			}
//			System.out.println("turn");
		}
	}
	
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
				if(findBnd(hb, r3, r0)==-1 || (i>0 && findBnd(hb, r3, r0)!=-1)){
					continue;
				}

			}
			if(r4 != null){
				if((i<c.getResidues().size()-3 && findBnd(hb, r4, r1)!=-1))
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
				Turn t1 = new Turn(r1);
				t1.setTurnType(turnType);
				r.add(i, t1); r.remove(i+1);
			}
			if(r2.getAsn().equals("C")){
				Turn t2 = new Turn(r2);
				t2.setTurnType(turnType);
				r.add(i+1, t2); r.remove(i+2);
			}
			if(r3.getAsn().equals("C")){
				Turn t3 = new Turn(r3);
				t3.setTurnType(turnType);
				r.add(i+2, t3); r.remove(i+3);
			}
		}
	}
	
	/**
	 * Condition for turn
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
	
	protected int findBnd(ArrayList<HBond> hb, Residue r1, Residue r2){
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
	
	protected int findPolInt(ArrayList<HBond> hb, Residue r1, Residue r2){
		int h;
//		System.out.println(r2.getNBondAcc());
		if(r1.getNBondDnr()!=0 && r2.getNBondAcc()!=0){
			for(int i=0; i<r1.getNBondDnr(); i++){
				h = r1.getHBondDnr()[i];
//				System.out.println(h);
				for(int j=0; j<r2.getNBondAcc(); j++){
					if(h == r2.getHBondAcc()[j] && h!=0){
						return h;
					}
				}
			}
		}
		
		return -1;
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
		
		int hc = 0;
		for(int i=0; i<donors.size(); i++){
			for(int j=0; j<acceptors.size(); j++){
				Acceptor acc = acceptors.get(j);
				Donor don = donors.get(i);
				
				if(acc.getAcceptorAtomO() == null || don.getDonorAtomH() == null) {
					continue;
				}
				
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
					
					int posC = PredictorUtility.findChain(mol, don.getChain().getName());
					if(posC != -1){
						if(mol.getChains().get(posC).getResidues().get(don.getResidueNum()).getNBondDnr() < 6){
							mol.getChains().get(posC).getResidues().get(don.getResidueNum())
							.setNBondDnr(mol.getChains().get(posC).getResidues().get(don.getResidueNum()).getNBondDnr()+1);
							mol.getChains().get(posC).getResidues().get(don.getResidueNum())
							.setHBondDnr(hc, mol.getChains().get(posC).getResidues().get(don.getResidueNum()).getNBondDnr());
						} else {
							System.out.println("Error"); 
						}
					}
					int posD = PredictorUtility.findChain(mol, acc.getChain().getName());
					if(posD != -1){
						if(mol.getChains().get(posD).getResidues().get(acc.getResidueNum()).getNBondAcc() < 6){
							mol.getChains().get(posD).getResidues().get(acc.getResidueNum())
							.setNBondAcc(mol.getChains().get(posD).getResidues().get(acc.getResidueNum()).getNBondAcc()+1);
							mol.getChains().get(posD).getResidues().get(acc.getResidueNum())
							.setHBondAcc(hc, mol.getChains().get(posD).getResidues().get(acc.getResidueNum()).getNBondAcc());
						} else {
							System.out.println("Error");
						}
					}
					
					if(posC!=posD && posC!=-1){
						mol.getChains().get(posC).getResidues().get(don.getResidueNum()).setInterChainHBonds(true);
						mol.getChains().get(posD).getResidues().get(acc.getResidueNum()).setInterChainHBonds(true);
					}	
				}
				hc++;
				hBond.add(hb);
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
			int posH = PredictorUtility.findAtom(r, "H");
			if(i==0){
				if(posH == -1){
					d = new Donor(
							c,
							i,
							r.getAtomList().get(PredictorUtility.findAtom(r, "CA"))
							);
				} else {
					d = new Donor(
							c, 
							i, 
							r.getAtomList().get(PredictorUtility.findAtom(r, "CA")),
							r.getAtomList().get(posH)
							);
				}
				donors.add(d);
			} else {
				if(posH == -1){
					d = new Donor(
							c, 
							i, 
							r.getAtomList().get(PredictorUtility.findAtom(r, "N")),
							r.getAtomList().get(PredictorUtility.findAtom(c.getResidues().get(i-1), "C")),
							r.getAtomList().get(PredictorUtility.findAtom(r, "CA"))
							);
				} else {
					d = new Donor(
							c, 
							i, 
							r.getAtomList().get(PredictorUtility.findAtom(r, "N")),
							r.getAtomList().get(PredictorUtility.findAtom(c.getResidues().get(i-1), "C")),
							r.getAtomList().get(PredictorUtility.findAtom(r, "CA")),
							r.getAtomList().get(posH)
							);
				}
				donors.add(d);
			}
			
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
//		System.out.println(hb.getP());
		
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
//			System.out.println("Ep - 1: " + hb.getEp());
		} else {
			hb.setEp(0.0);
//			System.out.println("Ep - 2: " + hb.getEp());
		}
		
		hb.setHBondEnergy(1000.0 * hb.getEr() * hb.getEt() * hb.getEp());
		
		return hb.getHBondEnergy();
	}
}
