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
		
		for(int i=0; i<mol.getChains().size(); i++){
			isHelix(mol.getChains().get(i), HBond);
		}
		for(int i=0; i<mol.getChains().size(); i++){
			for(int j=0; j<mol.getChains().get(i).getResidues().size(); j++){
				System.out.println(mol.getChains().get(i).getResidues().get(j).getAsn());
//				for(int k=0; k<mol.getChains().get(i).getResidues().get(j).getAtomList().size(); k++){
//					System.out.println(mol.getChains().get(i).getResidues().get(j).getAtomList().get(k).getSymbol())
//				}
			} System.out.println();
		}
	}
	
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
		
		for(int i=0; i< c.getResidues().size()-5; i++){
			if(prob[i] < THRESHOLD_H1 && prob[i+1] < THRESHOLD_H1){
				ArrayList<Residue> r = c.getResidues();
				
				HelixAlpha h1 = (HelixAlpha) r.get(i+1);
				r.add(i+1, h1); r.remove(i+2);
				HelixAlpha h2 = (HelixAlpha) r.get(i+2);
				r.add(i+2, h2); r.remove(i+3);
				HelixAlpha h3 = (HelixAlpha) r.get(i+3);
				r.add(i+3, h3); r.remove(i+4);
				HelixAlpha h4 = (HelixAlpha) r.get(i+4);
				r.add(i+4, h4); r.remove(i+5);
				
				if(RamaMap.calHelixProb(r.get(i).getPhi(), r.get(i).getPsi()) > THRESHOLD_H3){
					HelixAlpha h0 = (HelixAlpha) r.get(i);
					r.add(i, h0); r.remove(i+1);
				}
				if(RamaMap.calHelixProb(r.get(i+5).getPhi(), r.get(i+5).getPsi()) > THRESHOLD_H4){
					HelixAlpha h5 = (HelixAlpha) r.get(i+5);
					r.add(i+5, h5); r.remove(i+6);
				}
			}
		}
		
		for(int i=0; i<c.getResidues().size()-4; i++){
			ArrayList<Residue> r = c.getResidues();
			if(findBnd(hb, r.get(i+3), r.get(i))!=-1 && findBnd(hb, r.get(i+4), r.get(i+1))!=-1 &&
					(r.get(i+1).getAsn()!="H" && r.get(i+2).getAsn()!="H" && r.get(i+3).getAsn()!="H" && r.get(i+4).getAsn()!="H")){
				HelixThreeTen h1 = (HelixThreeTen) r.get(i+1);
				r.add(i+1, h1); r.remove(i+2);
				HelixThreeTen h2 = (HelixThreeTen) r.get(i+2);
				r.add(i+2, h2); r.remove(i+3);
				HelixThreeTen h3 = (HelixThreeTen) r.get(i+3);
				r.add(i+3, h3); r.remove(i+4);
			}
		}
		
		for(int i=0; i<c.getResidues().size()-6; i++){
			ArrayList<Residue> r = c.getResidues();
			if(findBnd(hb, r.get(i+5), r.get(i))!=-1 && findBnd(hb, r.get(i+6), r.get(i+1))!=-1 &&
					r.get(i+1).getAsn()=="C" && r.get(i+2).getAsn()=="C" && r.get(i+3).getAsn()=="C" &&
					r.get(i+4).getAsn()=="N" && r.get(i+5).getAsn()=="N"){
				HelixPi h1 = (HelixPi) r.get(i+1);
				r.add(i+1, h1); r.remove(i+2);
				HelixPi h2 = (HelixPi) r.get(i+2);
				r.add(i+2, h2); r.remove(i+3);
				HelixPi h3 = (HelixPi) r.get(i+3);
				r.add(i+3, h3); r.remove(i+4);
				HelixPi h4 = (HelixPi) r.get(i+4);
				r.add(i+4, h4); r.remove(i+5);
				HelixPi h5 = (HelixPi) r.get(i+5);
				r.add(i+5, h5); r.remove(i+6);
			}
		}
	}
	
	protected int findBnd(ArrayList<HBond> hb, Residue r1, Residue r2){
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
	
	protected int findPolInt(ArrayList<HBond> hb, Residue r1, Residue r2){
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
						if(mol.getChains().get(posD).getResidues().get(acc.getResidueNum()).getNBondDnr() < 6){
							mol.getChains().get(posD).getResidues().get(acc.getResidueNum())
							.setNBondDnr(mol.getChains().get(posD).getResidues().get(acc.getResidueNum()).getNBondDnr()+1);
							mol.getChains().get(posD).getResidues().get(acc.getResidueNum())
							.setHBondDnr(hc, mol.getChains().get(posD).getResidues().get(acc.getResidueNum()).getNBondDnr());
						} else {
							System.out.println("Error");
						}
					}
					
					if(posC!=posD && posC!=-1){
						mol.getChains().get(posC).getResidues().get(don.getResidueNum()).setInterChainHBonds(true);
						mol.getChains().get(posD).getResidues().get(acc.getResidueNum()).setInterChainHBonds(true);
					}
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
				donors.add(d);
			} else {
				d = new Donor(
						c, 
						i, 
						r.getAtomList().get(PredictorUtility.findAtom(r, "N")),
						r.getAtomList().get(PredictorUtility.findAtom(c.getResidues().get(i-1), "C")),
						r.getAtomList().get(PredictorUtility.findAtom(r, "CA")),
						r.getAtomList().get(PredictorUtility.findAtom(r, "H"))
						);
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
