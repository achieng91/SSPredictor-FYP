package predictor.core;

import predictor.core.model.Model;
import predictor.core.model.Molecule;
import predictor.core.model.Chain;
import predictor.core.model.Residue;
import predictor.core.model.Atom;
import predictor.core.model.Donor;
import predictor.core.model.Acceptor;
import predictor.core.model.HBond;
import predictor.core.model.Pattern;

import predictor.core.model.math.Vector3D;
import predictor.core.model.math.PhiPsiMap;
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
	protected final double C1_E = -0.2, C2_E = 0.2;
	protected final double THRESHOLD_H1 = -230.0;
	protected final double THRESHOLD_H3 = 0.12;
	protected final double THRESHOLD_H4 = 0.06;
	protected final double THRESHOLD_E1 = -240.0;
	protected final double THRESHOLD_E2 = -310.0;

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
		findHBonds();
//		noDoubleHBond();
		
		for(int i=0; i<mol.getChains().size(); i++){
			isHelix(mol.getChains().get(i), hBond);
			for(int j=0; j<mol.getChains().size(); j++){
				isSheet(mol.getChains().get(i), mol.getChains().get(j), i, j, hBond);
				
				isBetaTurn(mol.getChains().get(i));
				isGammaTurn(mol.getChains().get(i), hBond);
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
				prob[i] = hb.get(bondNum).getHBondEnergy() * (0.5*(PhiPsiMap.calProb(r.get(i).getPhi(), r.get(i).getPsi(), "HELIX")+
						PhiPsiMap.calProb(r.get(i+4).getPhi(), r.get(i+4).getPsi(), "HELIX")));
			}
		}
		
		// Check for Alpha Helix
		for(int i=0; i< c.getResidues().size()-5; i++){
			if(Double.compare(prob[i],THRESHOLD_H1)<0 && Double.compare(prob[i+1],THRESHOLD_H1)<0){
				ArrayList<Residue> r = c.getResidues();
				
				HelixAlpha h1 = new HelixAlpha(r.get(i+1));
				r.add(i+1, h1); r.remove(i+2);
				HelixAlpha h2 = new HelixAlpha(r.get(i+2));
				r.add(i+2, h2); r.remove(i+3);
				HelixAlpha h3 = new HelixAlpha(r.get(i+3));
				r.add(i+3, h3); r.remove(i+4);
				HelixAlpha h4 = new HelixAlpha(r.get(i+4));
				r.add(i+4, h4); r.remove(i+5);
				
				if(Double.compare(PhiPsiMap.calProb(r.get(i).getPhi(), r.get(i).getPsi(), "HELIX"),THRESHOLD_H3)>0){
					HelixAlpha h0 = new HelixAlpha(r.get(i));
					r.add(i, h0); r.remove(i+1);
				}
				if(Double.compare(PhiPsiMap.calProb(r.get(i+5).getPhi(), r.get(i+5).getPsi(), "HELIX"),THRESHOLD_H4)>0){
					HelixAlpha h5 = new HelixAlpha(r.get(i+5));
					r.add(i+5, h5); r.remove(i+6);
				}
			}
		}
		
		// Check for 310 Helix
		for(int i=0; i<c.getResidues().size()-4; i++){
			ArrayList<Residue> r = c.getResidues();
			if(findBnd(hb, r.get(i+3), r.get(i))!=-1 && findBnd(hb, r.get(i+4), r.get(i+1))!=-1 &&
					(!r.get(i+1).getAsn().equals("H") && !r.get(i+2).getAsn().equals("H") || 
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
	 * Check for beta sheets or bridges
	 * @param c1
	 * @param c2
	 * @param cn1
	 * @param cn2
	 * @param hb
	 */
	protected void isSheet(Chain c1, Chain c2, int cn1, int cn2, ArrayList<HBond> hb){
		Residue r1, r3;
		Residue r2 = new Residue();	Residue r4 = new Residue();
		Residue rA = new Residue();	Residue rB = new Residue();
		Residue r1m1 = new Residue(); Residue r3p1 = new Residue();
		ArrayList<Pattern> patN = new ArrayList<Pattern>();
		ArrayList<Pattern> patP = new ArrayList<Pattern>();
		int beg, patCntN=0, patCntP=0;
		String[] antiPar1 = new String[c1.getResidues().size()];
		String[] antiPar2 = new String[c2.getResidues().size()];
		String[] par1 = new String[c1.getResidues().size()];
		String[] par2 = new String[c2.getResidues().size()];
		
		for(int i=0; i<c1.getResidues().size(); i++){
			antiPar1[i] = "C";
			par1[i] = "C";
		}
		for(int i=0; i<c2.getResidues().size(); i++){
			antiPar2[i] = "C";
			par2[i] = "C";
		}
		
		for(int i=0; i<c1.getResidues().size(); i++){
			r1 = c1.getResidues().get(i);
			if((r1.getNBondDnr()==0 && r1.getNBondAcc()==0) || ((cn1!=cn2) && !r1.getInterChainHBonds())) {
				continue;
			}
			
			if(i!=0){
				r1m1 = c1.getResidues().get(i-1);
			}
			if(i+1 < c1.getResidues().size()) {
				rA = c1.getResidues().get(i+1);
			}
			if(i+2 < c1.getResidues().size()){
				r2 = c1.getResidues().get(i+2);
			}
			
			if(cn1 != cn2){
				beg = 0;
			} else {
				beg = i+1;
			}
			
			for(int j=beg; j<c2.getResidues().size(); j++) {
				r3 = c2.getResidues().get(j);
				if((r3.getNBondAcc()==0 && r3.getNBondDnr()==0) || ((cn1!=cn2) && r3.getInterChainHBonds())) {
					continue;
				}
				
				if(j+1<c2.getResidues().size()) {
					r3p1 = c2.getResidues().get(j+1);
				}
				if(j-1>=0){
					rB = c2.getResidues().get(j-1);
				}
				if(j-2>=0){
					r4 = c2.getResidues().get(j-2);
				}
				
				if(cn1!=cn2 || j-i>=3) {
					link(hb, c1, c2, cn1, cn2, r1, r3, r3, r1, r1, r3, THRESHOLD_E1, patN, patCntN, "1331");
				}
				if((i+2)<c1.getResidues().size() && ((cn1!=cn2 && j-2>=0) || (j-2)-(i+2)>=2)){
					link(hb, c1, c2, cn2, cn1, r3, r1, r2, r4, rB, rA, THRESHOLD_E1, patN, patCntN, "3124");
				}
				if(((cn1!=cn2 && (j-1)>=0) || (j-1)-i>4) && ((j-2)>=c1.getResidues().size() || (cn1==cn2 && j-(i-1)<=4) || 
						link(hb, c1, c2, cn1, cn2, r1, r3, r3, rA, null, r3, THRESHOLD_E1, patN, patCntN, "133A")==-1) && (i+1<0 || 
						link(hb, c1, c2, cn1, cn2, r1m1, rB, rB, r1, r1, null, THRESHOLD_E1, patN, patCntN, "1-BB1")==-1)){
					link(hb, c1, c2, cn1, cn2, r1, r3, rB, r1, r1, null, THRESHOLD_E1, patN, patCntN, "13B1");
				}
				if(((i+1)<c1.getResidues().size() && (cn1!=cn2 || j-(i+1)>4)) && ((cn1==cn2 && (j-1)-i<=4) || (cn1!=cn2 && j-1<0) || 
						link(hb, c1, c2, cn1, cn2, r1, r3, rB, r1, r1, null, THRESHOLD_E1, patN, patCntN, "13B1")==-1) && (j+1>=c2.getResidues().size() || 
						link(hb, c1, c2, cn1, cn2, rA, r3p1, r3, rA, rA, null, THRESHOLD_E1, patN, patCntN, "A3+3A")==-1)){
					link(hb, c1, c2, cn1, cn2, r1, r3, r3, rA, null, r3, THRESHOLD_E1, patN, patCntN, "133A");
				}
				
				if((cn1==cn2 && Math.abs(j-i)<=3) || (j+2)>=c2.getResidues().size()){
					continue;
				}
				rB = c2.getResidues().get(j+1);
				r4 = c2.getResidues().get(j+2);
				
				if((i+1)<c1.getResidues().size() && (cn1!=cn2 || Math.abs((i+1)-j)>3)){
					link(hb, c1, c2, cn2, cn1, r3, r1, r2, r3, r3, rA, THRESHOLD_E2, patP, patCntP, "3123");
				}
				if((j+2)<c2.getResidues().size() && (cn1!=cn2 || Math.abs((j+2)-i)>3)){
					link(hb, c1, c2, cn1, cn2, r1, r3, r4, r1, r1, rB, THRESHOLD_E2, patP, patCntP, "1341");
				}
			}
		}
		
		filterAntiPar(patN, patCntN);
		filterPar(patP, patCntP);
		
		mergePatternsAntiPar(patN, patCntN);
		mergePatternsPar(patP, patCntP);
		
		fillAsnAntiPar(antiPar1, antiPar2, c1, cn1, cn2, patN, patCntN);
		fillAsnPar(par1, par2, c1, cn1, cn2, patP, patCntP);
		
		bridge(antiPar1, antiPar2, c1, cn1, cn2, patN, patCntN);
		bridge(par1, par2, c1, cn1, cn2, patP, patCntP);

		for(int i=0; i<c1.getResidues().size(); i++){
			ArrayList<Residue> r = c1.getResidues();
			if(antiPar1[i].equals("N") || par1[i].equals("P")){
				Sheet s = new Sheet(r.get(i));
				r.add(i,s); r.remove(i+1);
			} else if(antiPar1[i].equals("B") || par1[i].equals("B")){
				Bridge b = new Bridge(r.get(i));
				r.add(i,b); r.remove(i+1);
			} else if(antiPar1[i].equals("b") || par1[i].equals("b")){
				Bridge b = new Bridge(r.get(i));
				b.setAsn("b");
				r.add(i,b); r.remove(i+1);
			}
		}
		
		for(int i=0; i<c2.getResidues().size(); i++){
			ArrayList<Residue> r = c2.getResidues();
			if(r.get(i).getAsn().equals("E")){
				continue;
			} else if(antiPar2[i].equals("N") || par2[i].equals("P")){
				Sheet s = new Sheet(r.get(i));
				r.add(i,s); r.remove(i+1);
			} else if(antiPar2[i].equals("B") || par2[i].equals("B")){
				Bridge b = new Bridge(r.get(i));
				r.add(i,b); r.remove(i+1);
			} else if(antiPar2[i].equals("b") || par2[i].equals("b")){
				Bridge b = new Bridge(r.get(i));
				b.setAsn("b");
				r.add(i,b); r.remove(i+1);
			}
		}
	}
	
	protected int link(ArrayList<HBond> hb, Chain c1, Chain c2, int cn1, int cn2, Residue r1_1, Residue r1_2, Residue r2_2, Residue r2_1, 
			Residue cR1, Residue cR2, double threshold, ArrayList<Pattern> pat, int numPat, String text) {
		int bondNumber1, bondNumber2;
		double conf, coeff, prob1, prob2;
		
		bondNumber1 = findPolInt(hb, r1_1, r1_2);
		bondNumber2 = findPolInt(hb, r2_2, r2_1);
		if(bondNumber1==-1 || bondNumber2==-1) {
			return -1;
		}
		if(cR1 == null){
			conf = PhiPsiMap.calProb(cR2.getPhi(), cR2.getPsi(), "SHEET");
		} else if(cR2 == null){
			conf = PhiPsiMap.calProb(cR1.getPhi(), cR1.getPsi() , "SHEET");
		} else {
			conf = 0.5 * (PhiPsiMap.calProb(cR1.getPhi(), cR1.getPsi(), "SHEET") + PhiPsiMap.calProb(cR2.getPhi(), cR2.getPsi(), "SHEET"));
		}
		coeff = 1 + C1_E + C2_E*conf;
		prob1 = hb.get(bondNumber1).getHBondEnergy() * coeff;
		prob2 = hb.get(bondNumber2).getHBondEnergy() * coeff;
		
		if(prob1<threshold && prob2<threshold) {
			Pattern p = pat.get(numPat);
			p.setExistPattern(true);
			p.setHBond1(hb.get(bondNumber1));
			p.setHBond2(hb.get(bondNumber2));
			p.setType(text);
			numPat++;
			
			return 1;
		} else {
			return -1;
		}
	}
	
	protected void bridge(String[] asn1, String[] asn2, Chain c, int cn1, int cn2, ArrayList<Pattern> pat, int NPat){
		int B_Res=0;
		
		for(int i=0; i<NPat; i++){
			if(pat.get(i).getNei1()!=null || pat.get(i).getNei2()!=null){
				continue;
			}
			HBond patHb1 = pat.get(i).getHBond1();
			HBond patHb2 = pat.get(i).getHBond2();
			if(!pat.get(i).getType().equals("1331") && (cn1!=cn2 || 
					Math.abs(patHb1.getDonor().getResidueNum()-patHb1.getAcceptor().getResidueNum())>=3)){
				
				if(patHb1.getDonor().getChain().getName().equals(c.getName())){
					if(asn1[patHb1.getDonor().getResidueNum()].equals("C")){
						asn1[patHb1.getDonor().getResidueNum()] = "B";
					}
					if(asn2[patHb1.getAcceptor().getResidueNum()].equals("C")){
						asn2[patHb1.getAcceptor().getResidueNum()] = "B";
					}
					
				} else if(pat.get(i).getType().equals("3124") && (cn1!=cn2 || 
						Math.abs(patHb1.getDonor().getResidueNum()-patHb1.getAcceptor().getResidueNum())>=2 && 
						Math.abs(patHb2.getDonor().getResidueNum()-patHb2.getAcceptor().getResidueNum())>= 2)){
					
					if(patHb1.getDonor().getChain().getName().equals(c.getName())){
						if(patHb1.getDonor().getResidueNum() > patHb2.getAcceptor().getResidueNum()){
							B_Res = patHb1.getDonor().getResidueNum() - 1;
						} else {
							B_Res = patHb1.getDonor().getResidueNum() + 1;
						}
						
						if(asn1[B_Res].equals("C")){
							asn1[B_Res] = "B";
						}
						
						if(patHb2.getDonor().getResidueNum() > patHb1.getAcceptor().getResidueNum()){
							B_Res = patHb2.getDonor().getResidueNum() - 1;
						} else {
							B_Res = patHb2.getDonor().getResidueNum() + 1;
						}
						
						if(asn2[B_Res].equals("C")){
							asn2[B_Res] = "B";
						} 
					} else {
						if(patHb1.getDonor().getResidueNum() > patHb2.getAcceptor().getResidueNum()){
							B_Res = patHb1.getDonor().getResidueNum() - 1;
						} else {
							B_Res = patHb1.getDonor().getResidueNum() + 1;
						}
						
						if(asn2[B_Res].equals("C")){
							asn2[B_Res] = "B";
						} 
						
						if(patHb2.getDonor().getResidueNum() > patHb1.getAcceptor().getResidueNum()){
							B_Res = patHb2.getDonor().getResidueNum() - 1;
						} else {
							B_Res = patHb2.getDonor().getResidueNum() + 1;
						}
						
						if(asn1[B_Res].equals("C")){
							asn1[B_Res] = "B";
						} 
					}
				} else if((!pat.get(i).getType().equals("3123") || pat.get(i).getType().equals("1341")) &&
							(cn1!=cn2 || Math.abs(patHb1.getDonor().getResidueNum()-patHb1.getAcceptor().getResidueNum())>3 &&
							Math.abs(patHb2.getDonor().getResidueNum()-patHb2.getAcceptor().getResidueNum())>3)) {
					
					if(patHb1.getDonor().getChain().getName().equals(c.getName())){
						if(patHb1.getDonor().getResidueNum()==patHb2.getAcceptor().getResidueNum()){
							if(asn1[patHb1.getDonor().getResidueNum()].equals("C")){
								asn1[patHb1.getDonor().getResidueNum()] = "B";
							}
							
							if(patHb2.getDonor().getResidueNum() > patHb1.getAcceptor().getResidueNum()){
								B_Res = patHb2.getDonor().getResidueNum() - 1;
							} else {
								B_Res = patHb2.getDonor().getResidueNum() + 1;
							}
							
							if(asn2[B_Res].equals("C")){
								asn2[B_Res] = "B";
							} 
						} else if(patHb2.getDonor().getResidueNum()==patHb1.getAcceptor().getResidueNum()){
							if(asn2[patHb2.getDonor().getResidueNum()].equals("C")){
								asn2[patHb2.getDonor().getResidueNum()] = "B";
							}
							
							if(patHb1.getDonor().getResidueNum() > patHb2.getAcceptor().getResidueNum()){
								B_Res = patHb1.getDonor().getResidueNum() - 1;
							} else {
								B_Res = patHb1.getDonor().getResidueNum() + 1;
							}
							
							if(asn1[B_Res].equals("C")){
								asn1[B_Res] = "B";
							} 
						}
					}
				} else if((!pat.get(i).getType().equals("13B1") || pat.get(i).getType().equals("133A")) &&
						(cn1!=cn2 || Math.abs(patHb1.getDonor().getResidueNum()-patHb1.getAcceptor().getResidueNum())>4 &&
						Math.abs(patHb2.getDonor().getResidueNum()-patHb2.getAcceptor().getResidueNum())>4)){
					
					if(patHb1.getDonor().getChain().getName().equals(c.getName())){
						if(patHb1.getDonor().getResidueNum()==patHb2.getAcceptor().getResidueNum()){
							if(asn1[patHb1.getDonor().getResidueNum()].equals("C")){
								asn1[patHb1.getDonor().getResidueNum()] = "B";
							}
							
							if(patHb2.getDonor().getResidueNum() > patHb1.getAcceptor().getResidueNum()){
								B_Res = patHb2.getDonor().getResidueNum() - 1;
							} else {
								B_Res = patHb2.getDonor().getResidueNum() + 1;
							}
							
							if(asn2[B_Res].equals("C")){
								asn2[B_Res] = "B";
							} 
						} else if(patHb2.getDonor().getResidueNum()==patHb1.getAcceptor().getResidueNum()){
							if(asn2[patHb2.getDonor().getResidueNum()].equals("C")){
								asn2[patHb2.getDonor().getResidueNum()] = "b";
							}
							
							if(patHb1.getDonor().getResidueNum() > patHb2.getAcceptor().getResidueNum()){
								B_Res = patHb1.getDonor().getResidueNum() - 1;
							} else {
								B_Res = patHb1.getDonor().getResidueNum() + 1;
							}
							
							if(asn1[B_Res].equals("C")){
								asn1[B_Res] = "b";
							} 
						}
					}
				}
			}
		}
	}
	
	protected void filterAntiPar(ArrayList<Pattern> pat, int NPat){
		int i1A=0, i1D=0, i2A=0, i2D=0, j1D=0, j1A=0, j2A=0, j2D=0;
		String i1ACn="", i1DCn="", i2ACn="", i2DCn="", j1ACn="", j1DCn="", j2ACn="", j2DCn="";
		
		for(int i=0; i<NPat; i++){
			if(!pat.get(i).getExistPattern()){
				continue;
			}
			
			Alias(i1D, i1A, i2D, i2A, i1DCn, i1ACn, i2DCn, i2ACn, pat.get(i));
			
			for(int j=0; j<NPat; j++){
				if(j==i || !pat.get(j).getExistPattern()){
					continue;
				}
				
				Alias(j1D, j1A, j2D, j2A, j1DCn, j1ACn, j2DCn, j2ACn, pat.get(j));
				
				if(j1D==j2A && j2D==j1A && i1D!=i2A && i2D!=i1A && ((j1D==i1D && j1A==i1A) || (j1D==i1A && j1A==i1D) || 
						(j1D==i2A && j1A==i2D) || (j1D==i2D && j1A==i2A))){
					continue;
				}
				if(((i1D<i2A || i2D<i1A) && ((j1A<=i2A && j1A>=i1D && j2D<=i2A && j2D>=i1D && j2DCn.equals(i1DCn) &&
						j2A<=i1A && j2A>=i2D && j1D<=i1A && j1D>=i2D && j1DCn.equals(i2DCn)) || 
						(j2A<= i2A && j2A>=i1D && j1D<=i2A && j1D>=i1D && j1DCn.equals(i1DCn) && j1A<=i1A && j1A>=i2D && 
						j2D<=i1A && j2D>=i2D && j2DCn.equals(i2DCn)))) || ((i1D>i2A || i2D>i1A) && 
						((j1A>=i2A && j1A<=i1D && j2D>=i2A && j2D<=i1D && j2DCn.equals(i1DCn) && j2A>=i1A && j2A<=i2D && 
						j1D>=i1A && j1D<=i2D && j1DCn.equals(i2DCn)) || (j2A>=i2A && j2A<=i1D && j1D>=i2A && j1D<=i1D && 
						j1DCn.equals(i1DCn) && j1A>=i1A && j1A<=i2D && j2D>=i1A && j2D<=i2D && j2DCn.equals(i2DCn))))){
					pat.get(j).setExistPattern(false);
				}
			}
		}
	}
	
	protected void filterPar(ArrayList<Pattern> pat, int NPat){
		int i1A=0, i1D=0, i2A=0, i2D=0, j1A=0, j1D=0, j2A=0, j2D=0;
		String i1ACn="", i1DCn="", i2ACn="", i2DCn="", j1ACn="", j1DCn="", j2ACn="", j2DCn="";
		
		for(int i=0; i<NPat; i++){
			if(!pat.get(i).getExistPattern()){
				continue;
			}
			
			Alias(i1D, i1A, i2D, i2A, i1DCn, i1ACn, i2DCn, i2ACn, pat.get(i));
			
			for(int j=0; j<NPat; j++){
				if(j==i || !pat.get(j).getExistPattern()){
					continue;
				}
				
				Alias(j1D, j1A, j2D, j2A, j1DCn, j1ACn, j2DCn, j2ACn, pat.get(j));
				
				if(((i1A>=i2D && i1D>=i2A) && ((j1A>=i2A && j1A<=i1D && j2D>=i2A && j2D<=i1D && j2DCn.equals(i1DCn) && 
						j2A<=i1A && j2A>=i2D && j1D<=i1A && j1D>=i2D && j1DCn.equals(i2DCn)) || (j2A>=i2A && j2A<=i1D && 
						j1D>=i2A && j1D<=i1D && j1DCn.equals(i1DCn) && j1A<=i1A && j1A>=i2D && j2D<=i1A && j2D>=i2D &&
						j2DCn.equals(i2DCn)))) || (i2A>=i1D && i2D>=i1A && ((j1A<=i2A && j1A>=i1D && j2D<=i2A && j2D>=i1D &&
						j2DCn.equals(i1DCn) && j2A>=i1A && j2A<=i2D && j1D>=i1A && j1D<=i2D && j1DCn.equals(i2DCn)) ||
						(j2A<=i2A && j2A>=i1D && j1D<=i2A && j1D>=i1D && j1DCn.equals(i2DCn) && j1A>=i1A && j1A<=i2D && 
						j2D>=i1A && j2D<=i2D && j2DCn.equals(i2DCn))))){
					pat.get(j).setExistPattern(false);
				}
			}
		}
	}
	
	protected void mergePatternsAntiPar(ArrayList<Pattern> pat, int NPat){
		int dB=0, dW=0, minDB1=0, minDB2=0, minDW1=0, minDW2=0, min=0, lnk1A=0, lnk1D=0;
		int i1A=0, i1D=0, i2A=0, i2D=0, j1A=0, j1D=0, j2A=0, j2D=0;
		String i1ACn="", i1DCn="", i2ACn="", i2DCn="", j1ACn="", j1DCn="", j2ACn="", j2DCn="";
		
		for(int i=0; i<NPat; i++){
			if(!pat.get(i).getExistPattern()){
				continue;
			}
			minDB1 = minDB2 = minDW1 = minDW2 = 1000;
			
			Alias(i1D, i1A, i2D, i2A, i1DCn, i1ACn, i2DCn, i2ACn, pat.get(i));
			
			for(int j=0; j<NPat; j++){
				if(i==j || !pat.get(j).getExistPattern()){
					continue;
				}
				
				Alias(j1D, j1A, j2D, j2A, j1DCn, j1ACn, j2DCn, j2ACn, pat.get(j));
				
				if(near(i1D, j1D, j1A, i1A, j2A, j2D, i2A, i2D, i1DCn, j1DCn, j1ACn, i1ACn, dB, dW)==1 && 
						((dB<minDB1 && dW<=minDW1) || (dB<=minDB1 && dW<minDW1)) && rightSide(j1A, j1D, i1A, i1D, i2A, i2D)==1){
					joinNeighbours(lnk1A, j2D, lnk1D, j2A, pat.get(i).getNei1(), pat.get(j), minDB1, dB, minDW1, dW, min, j);
				}
				
				if( near(i1D, j1A, j1D, i1A, j2D, j2A, i2A, i2D, i1DCn, j1ACn, j1DCn, i1ACn, dB, dW)==1 &&
						((dB<minDB1 && dW<minDW1) || (dB<=minDB1 && dW<minDW1)) && rightSide(j1D, j1A, i1A, i1D, i2A, i2D)==1){
					joinNeighbours(lnk1A, j2A, lnk1D, j2D, pat.get(i).getNei1(), pat.get(j), minDB1, dB, minDW1, dW, min, j);
				}
				
				if(near(i1D, j2D, j2A, i1A, j1A, j1D, i2A, i2D, i1DCn, j2DCn, j2ACn, i1ACn, dB, dW)==1 &&
						((dB<minDB1 && dW<=minDW1) || (dB<=minDB1 && dW<minDW1)) && rightSide(j2A, j2D, i1A, i1D, i2A, i2D)==1){
					joinNeighbours(lnk1A, j1D, lnk1D, j1A, pat.get(i).getNei1(), pat.get(j), minDB1, dB, minDW1, dW, min, j);
				}
				
				if(near(i1D, j2A, j2D, i1A, j1D, j1A, i2A, i2D, i1DCn, j2ACn, j2DCn, i1ACn, dB, dW)==1 &&
						((dB<minDB1 && dW<=minDW1) || (dB<=minDB1 && dW<minDW1)) && rightSide(j2D, j2A, i1A, i1D, i2A, i2D)==1){
					joinNeighbours(lnk1A, j1A, lnk1D, j1D, pat.get(i).getNei1(), pat.get(j), minDB1, dB, minDW1, dW, min, j);
				}
				
				if(near(i1A, j1D, j1A, i1D, j2A, j2D, i2D, i2A, i1ACn, j1DCn, j1ACn, i1DCn, dB, dW)==1 &&
						((dB<minDB1 && dW<=minDW1) || (dB<=minDB1 && dW<minDW1)) && rightSide(j1A, j1D, i1A, i1D, i2A, i2D)==1){
					joinNeighbours(lnk1A, j2A, lnk1D, j2D, pat.get(i).getNei1(), pat.get(j), minDB1, dB, minDW1, dW, min, j);
				}
				
				if(near(i1A, j1A, j1D, i1D, j2D, j2A, i2D, i2A, i1ACn, j1ACn, j1DCn, i1DCn, dB, dW)==1 &&
						((dB<minDB1 && dW<=minDW1) || (dB<=minDB1 && dW<minDW1)) && rightSide(j1A, j1D, i1A, i1D, i2A, i2D)==1){
					joinNeighbours(lnk1A, j2D, lnk1D, j2A, pat.get(i).getNei1(), pat.get(j), minDB1, dB, minDW1, dW, min, j);
				}
				
				if(near(i1A, j2D, j2A, i1D, j1A, j1D, i2D, i2A, i1ACn, j2DCn, j2ACn, i1DCn, dB, dW)==1 &&
						((dB<minDB1 && dW<=minDW1) || (dB<=minDB1 && dW<minDW1)) && rightSide(j2D, j2A, i1A, i1D, i2A, i2D)==1){
					joinNeighbours(lnk1A, j1A, lnk1D, j1D, pat.get(i).getNei1(), pat.get(j), minDB1, dB, minDW1, dW, min, j);
				}
				
				if(near(i1A, j2A, j2D, i1D, j1D, j1A, i2D, i2A, i1ACn, j2ACn, j2DCn, i1DCn, dB, dW)==1 &&
						((dB<minDB1 && dW<=minDW1) || (dB<=minDB1 && dW<minDW1)) && rightSide(j2A, j2D, i1A, i1D, i2A, i2D)==1){
					joinNeighbours(lnk1A, j1D, lnk1D, j2D, pat.get(i).getNei1(), pat.get(j), minDB1, dB, minDW1, dW, min, j);
				}
			}
			
			for(int j=0; j<NPat; j++){
				if(j==min || j==1 || !pat.get(j).getExistPattern()){
					continue;
				}
				
				Alias(j1D, j1A, j2D, j2A, j1DCn, j1ACn, j2DCn, j2ACn, pat.get(j));
				
				if(near(i2D, j1D, j1A, i2A, j2A, j2D, i1A, i1D, i2DCn, j1DCn, j1ACn, i2ACn, dB, dW)==1 &&
						((dB<minDB2 && dW<=minDW2) || (dB<=minDB2 && dW<minDW2)) && rightSide2(lnk1A, lnk1D, j2A, j2D, i1A, i1D, i2A, i2D)==1){
					joinNeighb(pat.get(i).getNei2(), pat.get(j), minDB2, dB, minDW2, dW);
				}
				
				if(near(i2D, j1A, j1D, i2A, j2D, j2A, i1A, i1D, i2DCn, j1ACn, j1DCn, i2ACn, dB, dW)==1 &&
						((dB<minDB2 && dW<=minDW2) || (dB<=minDB2 && dW<minDW2)) && rightSide2(lnk1A, lnk1D, j2D, j2A, i1A, i1D, i2A, i2D)==1){
					joinNeighb(pat.get(i).getNei2(), pat.get(j), minDB2, dB, minDW2, dW);
				}
				
				if(near(i2D, j2D, j2A, i2A, j1A, j1D, i1A, i1D, i2DCn, j2DCn, j2ACn, i2ACn, dB, dW)==1 &&
						((dB<minDB2 && dW<=minDW2) || (dB<=minDB2 && dW<minDW2)) && rightSide2(lnk1A, lnk1D, j1A, j1D, i1A, i1D, i2A, i2D)==1){
					joinNeighb(pat.get(i).getNei2(), pat.get(j), minDB2, dB, minDW2, dW);
				}
				
				if(near(i2D, j2A, j2D, i2A, j1D, j1A, i1A, i1D, i2DCn, j2ACn, j2DCn, i2ACn, dB, dW)==1 &&
						((dB<minDB2 && dW<=minDW2) || (dB<=minDB2 && dW<minDW2)) && rightSide2(lnk1A, lnk1D, j1D, j1A, i1A, i1D, i2A, i2D)==1){
					joinNeighb(pat.get(i).getNei2(), pat.get(j), minDB2, dB, minDW2, dW);
				}
				
				if(near(i2A, j1D, j1A, i2D, j2A, j2D, i1D, i1A, i2ACn, j1DCn, j1ACn, i2DCn, dB, dW)==1 &&
						((dB<minDB2 && dW<=minDW2) || (dB<=minDB2 && dW<minDW2)) && rightSide2(lnk1A, lnk1D, j2A, j2D, i1A, i1D, i2A, i2D)==1){
					joinNeighb(pat.get(i).getNei2(), pat.get(j), minDB2, dB, minDW2, dW);
				}
				
				if(near(i2A, j1A, j1D, i2D, j2D, j2A, i1D, i1A, i2ACn, j1ACn, j1DCn, i2DCn, dB, dW)==1 &&
						((dB<minDB2 && dW<=minDW2) || (dB<=minDB2 && dW<minDW2)) && rightSide2(lnk1A, lnk1D, j2A, j2D, i1A, i1D, i2A, i2D)==1){
					joinNeighb(pat.get(i).getNei2(), pat.get(j), minDB2, dB, minDW2, dW);
				}
				
				if(near(i2A, j2D, j2A, i2D, j1A, j1D, i1D, i1A, i2ACn, j2DCn, j2ACn, i2DCn, dB, dW)==1 &&
						((dB<minDB2 && dW<=minDW2) || (dB<=minDB2 && dW<minDW2)) && rightSide2(lnk1A, lnk1D, j1D, j1A, i1A, i1D, i2A, i2D)==1){
					joinNeighb(pat.get(i).getNei2(), pat.get(j), minDB2, dB, minDW2, dW);
				}
				
				if(near(i2A, j2A, j2D, i2D, j1D, j1A, i1D, i1A, i2ACn, j2ACn, j2DCn, i2DCn, dB, dW)==1 &&
						((dB<minDB2 && dW<=minDW2) || (dB<=minDB2 && dW<minDW2)) && rightSide2(lnk1A, lnk1D, j1A, j1D, i1A, i1D, i2A, i2D)==1){
					joinNeighb(pat.get(i).getNei2(), pat.get(j), minDB2, dB, minDW2, dW);
				}
			}
		}
	}
	
	protected void mergePatternsPar(ArrayList<Pattern> pat, int NPat){
		int dB=0, dW=0, minDB1=0, minDB2=0, minDW1=0, minDW2=0, min=0, lnk1A=0, lnk1D=0;
		int i1A=0, i1D=0, i2A=0, i2D=0, j1A=0, j1D=0, j2A=0, j2D=0;
		String i1ACn="", i1DCn="", i2ACn="", i2DCn="", j1ACn="", j1DCn="", j2ACn="", j2DCn="";
		
		for(int i=0; i<NPat; i++){
			if(!pat.get(i).getExistPattern()){
				continue;
			}
			minDB1 = minDB2 = minDW1 = minDW2 = 1000;
			
			Alias(i1D, i1A, i2D, i2A, i1DCn, i1ACn, i2DCn, i2ACn, pat.get(i));
			
			for(int j=0; j<NPat; j++){
				if(i==j || !pat.get(j).getExistPattern()){
					continue;
				}
				
				Alias(j1D, j1A, j2D, j2A, j1DCn, j1ACn, j2DCn, j2ACn, pat.get(j));
				
				if(nearPar(i1D, j1D, j1A, i1A, j2A, j2D, i2A, i2D, i1DCn, j1DCn, j1ACn, i1ACn, dB, dW)==1 && 
						((dB<minDB1 && dW<=minDW1) || (dB<=minDB1 && dW<minDW1)) && rightSidePar(j1A, j1D, i1A, i1D, i2A, i2D)==1){
					joinNeighbours(lnk1A, j2D, lnk1D, j2A, pat.get(i).getNei1(), pat.get(j), minDB1, dB, minDW1, dW, min, j);
				}
				
				if( nearPar(i1D, j1A, j1D, i1A, j2D, j2A, i2A, i2D, i1DCn, j1ACn, j1DCn, i1ACn, dB, dW)==1 &&
						((dB<minDB1 && dW<minDW1) || (dB<=minDB1 && dW<minDW1)) && rightSidePar(j1D, j1A, i1A, i1D, i2A, i2D)==1){
					joinNeighbours(lnk1A, j2A, lnk1D, j2D, pat.get(i).getNei1(), pat.get(j), minDB1, dB, minDW1, dW, min, j);
				}
				
				if(nearPar(i1D, j2D, j2A, i1A, j1A, j1D, i2A, i2D, i1DCn, j2DCn, j2ACn, i1ACn, dB, dW)==1 &&
						((dB<minDB1 && dW<=minDW1) || (dB<=minDB1 && dW<minDW1)) && rightSidePar(j2A, j2D, i1A, i1D, i2A, i2D)==1){
					joinNeighbours(lnk1A, j1D, lnk1D, j1A, pat.get(i).getNei1(), pat.get(j), minDB1, dB, minDW1, dW, min, j);
				}
				
				if(nearPar(i1D, j2A, j2D, i1A, j1D, j1A, i2A, i2D, i1DCn, j2ACn, j2DCn, i1ACn, dB, dW)==1 &&
						((dB<minDB1 && dW<=minDW1) || (dB<=minDB1 && dW<minDW1)) && rightSidePar(j2D, j2A, i1A, i1D, i2A, i2D)==1){
					joinNeighbours(lnk1A, j1A, lnk1D, j1D, pat.get(i).getNei1(), pat.get(j), minDB1, dB, minDW1, dW, min, j);
				}
				
				if(nearPar(i1A, j1D, j1A, i1D, j2A, j2D, i2D, i2A, i1ACn, j1DCn, j1ACn, i1DCn, dB, dW)==1 &&
						((dB<minDB1 && dW<=minDW1) || (dB<=minDB1 && dW<minDW1)) && rightSidePar(j1A, j1D, i1A, i1D, i2A, i2D)==1){
					joinNeighbours(lnk1A, j2A, lnk1D, j2D, pat.get(i).getNei1(), pat.get(j), minDB1, dB, minDW1, dW, min, j);
				}
				
				if(nearPar(i1A, j1A, j1D, i1D, j2D, j2A, i2D, i2A, i1ACn, j1ACn, j1DCn, i1DCn, dB, dW)==1 &&
						((dB<minDB1 && dW<=minDW1) || (dB<=minDB1 && dW<minDW1)) && rightSidePar(j1A, j1D, i1A, i1D, i2A, i2D)==1){
					joinNeighbours(lnk1A, j2D, lnk1D, j2A, pat.get(i).getNei1(), pat.get(j), minDB1, dB, minDW1, dW, min, j);
				}
				
				if(nearPar(i1A, j2D, j2A, i1D, j1A, j1D, i2D, i2A, i1ACn, j2DCn, j2ACn, i1DCn, dB, dW)==1 &&
						((dB<minDB1 && dW<=minDW1) || (dB<=minDB1 && dW<minDW1)) && rightSidePar(j2D, j2A, i1A, i1D, i2A, i2D)==1){
					joinNeighbours(lnk1A, j1A, lnk1D, j1D, pat.get(i).getNei1(), pat.get(j), minDB1, dB, minDW1, dW, min, j);
				}
				
				if(nearPar(i1A, j2A, j2D, i1D, j1D, j1A, i2D, i2A, i1ACn, j2ACn, j2DCn, i1DCn, dB, dW)==1 &&
						((dB<minDB1 && dW<=minDW1) || (dB<=minDB1 && dW<minDW1)) && rightSidePar(j2A, j2D, i1A, i1D, i2A, i2D)==1){
					joinNeighbours(lnk1A, j1D, lnk1D, j2D, pat.get(i).getNei1(), pat.get(j), minDB1, dB, minDW1, dW, min, j);
				}
			}
			
			for(int j=0; j<NPat; j++){
				if(j==min || j==1 || !pat.get(j).getExistPattern()){
					continue;
				}
				
				Alias(j1D, j1A, j2D, j2A, j1DCn, j1ACn, j2DCn, j2ACn, pat.get(j));
				
				if(nearPar(i2D, j1D, j1A, i2A, j2A, j2D, i1A, i1D, i2DCn, j1DCn, j1ACn, i2ACn, dB, dW)==1 &&
						((dB<minDB2 && dW<=minDW2) || (dB<=minDB2 && dW<minDW2)) && rightSide2(lnk1A, lnk1D, j2A, j2D, i1A, i1D, i2A, i2D)==1){
					joinNeighb(pat.get(i).getNei2(), pat.get(j), minDB2, dB, minDW2, dW);
				}
				
				if(nearPar(i2D, j1A, j1D, i2A, j2D, j2A, i1A, i1D, i2DCn, j1ACn, j1DCn, i2ACn, dB, dW)==1 &&
						((dB<minDB2 && dW<=minDW2) || (dB<=minDB2 && dW<minDW2)) && rightSide2(lnk1A, lnk1D, j2D, j2A, i1A, i1D, i2A, i2D)==1){
					joinNeighb(pat.get(i).getNei2(), pat.get(j), minDB2, dB, minDW2, dW);
				}
				
				if(nearPar(i2D, j2D, j2A, i2A, j1A, j1D, i1A, i1D, i2DCn, j2DCn, j2ACn, i2ACn, dB, dW)==1 &&
						((dB<minDB2 && dW<=minDW2) || (dB<=minDB2 && dW<minDW2)) && rightSide2(lnk1A, lnk1D, j1A, j1D, i1A, i1D, i2A, i2D)==1){
					joinNeighb(pat.get(i).getNei2(), pat.get(j), minDB2, dB, minDW2, dW);
				}
				
				if(nearPar(i2D, j2A, j2D, i2A, j1D, j1A, i1A, i1D, i2DCn, j2ACn, j2DCn, i2ACn, dB, dW)==1 &&
						((dB<minDB2 && dW<=minDW2) || (dB<=minDB2 && dW<minDW2)) && rightSide2(lnk1A, lnk1D, j1D, j1A, i1A, i1D, i2A, i2D)==1){
					joinNeighb(pat.get(i).getNei2(), pat.get(j), minDB2, dB, minDW2, dW);
				}
				
				if(nearPar(i2A, j1D, j1A, i2D, j2A, j2D, i1D, i1A, i2ACn, j1DCn, j1ACn, i2DCn, dB, dW)==1 &&
						((dB<minDB2 && dW<=minDW2) || (dB<=minDB2 && dW<minDW2)) && rightSide2(lnk1A, lnk1D, j2A, j2D, i1A, i1D, i2A, i2D)==1){
					joinNeighb(pat.get(i).getNei2(), pat.get(j), minDB2, dB, minDW2, dW);
				}
				
				if(nearPar(i2A, j1A, j1D, i2D, j2D, j2A, i1D, i1A, i2ACn, j1ACn, j1DCn, i2DCn, dB, dW)==1 &&
						((dB<minDB2 && dW<=minDW2) || (dB<=minDB2 && dW<minDW2)) && rightSide2(lnk1A, lnk1D, j2A, j2D, i1A, i1D, i2A, i2D)==1){
					joinNeighb(pat.get(i).getNei2(), pat.get(j), minDB2, dB, minDW2, dW);
				}
				
				if(nearPar(i2A, j2D, j2A, i2D, j1A, j1D, i1D, i1A, i2ACn, j2DCn, j2ACn, i2DCn, dB, dW)==1 &&
						((dB<minDB2 && dW<=minDW2) || (dB<=minDB2 && dW<minDW2)) && rightSide2(lnk1A, lnk1D, j1D, j1A, i1A, i1D, i2A, i2D)==1){
					joinNeighb(pat.get(i).getNei2(), pat.get(j), minDB2, dB, minDW2, dW);
				}
				
				if(nearPar(i2A, j2A, j2D, i2D, j1D, j1A, i1D, i1A, i2ACn, j2ACn, j2DCn, i2DCn, dB, dW)==1 &&
						((dB<minDB2 && dW<=minDW2) || (dB<=minDB2 && dW<minDW2)) && rightSide2(lnk1A, lnk1D, j1A, j1D, i1A, i1D, i2A, i2D)==1){
					joinNeighb(pat.get(i).getNei2(), pat.get(j), minDB2, dB, minDW2, dW);
				}
			}
		}
	}
	
	protected int near(int r1, int r2, int r3, int r4, int r5, int r6, int r7, int r8, String cn1, String cn2, String cn3, String cn4,
			int dB, int dW){
		int a,b, c1, d1, c, d, nei1, nei2;
		
		if(!cn1.equals(cn2) || !cn3.equals(cn4)){
			return -1;
		}
		
		if(r1>=r2 && r2>=r5 && r7>=r1 && r4<=r3 && r4<=r6 && r8<=r4){
			if(r5==r2){
				nei1 = r2;
			} else {
				nei1 = r2 - 1;
			}
			
			if(r1==r7){
				nei2 = r1;
			} else {
				nei2 = r1 + 1;
			}
			
			a = nei2 - nei1;
			c1 = nei2 - r5;
			
			if(r3==r6) {
				nei1 = r3;
			} else {
				nei1 = r3 + 1;
			}
			
			if(r4==r8){
				nei2 = r4;
			} else {
				nei2 = r4-1;
			}
			
			b = nei1 - nei2;
			d1 = r6 - nei2;
		} else {
			return -1;
		}
		
		c = Math.max(c1, a);
		d = Math.max(d1, b);
		
		if(a>=0 && b>=0 && c>=0 && d>=0 && ((a<=2 && b<=5) || (a<=5 && b<=2))){
			dB = Math.min(a, b);
			dW = Math.max(c, d);
			if(dB <= dW){
				return 1;
			} else {
				return -1;
			}
		}
		return -1;
	}
	
	protected int nearPar(int r1, int r2, int r3, int r4, int r5, int r6, int r7, int r8, String cn1, String cn2, String cn3, String cn4,
			int dB, int dW){
		int a, b, c1, d1, c, d, nei1, nei2;
		
		if(!cn1.equals(cn2) || !cn3.equals(cn4)){
			return -1;
		}
		
		if(r1>=r2 && r2>=r5 && r7>=r1 && r4>=r3 && r4>=r6 && r8>=r4){
			if(r5==r2){
				nei1 = r2;
			} else {
				nei1 = r2 - 1;
			}
			
			if(r1==r7){
				nei2 = r1;
			} else {
				nei2 = r1 + 1;
			}
			
			a = nei2 - nei1;
			c1 = nei2 - r5;
			
			if(r3==r6){
				nei1 = r3;
			} else {
				nei1 = r3 - 1;
			}
			
			if(r4==r8){
				nei2 = r4;
			} else {
				nei2 = r4 + 1;
			}
			
			b = nei2 - nei1;
			d1 = nei2 - r6;
		} else {
			if(r1<=r2 && r2<=r5 && r7<=r1 && r4<=r3 && r4<=r6 && r8<=r4){
				if(r5==r2){
					nei1 = r2;
				} else {
					nei1 = r2 + 1;
				}
				
				if(r1==r7){
					nei2 = r1;
				} else {
					nei2 = r1 - 1;
				}
				
				a = nei1 - nei2;
				c1 = r1 - r7; 
				
				if(r3==r6){
					nei1 = r3;
				} else {
					nei1 = r3 + 1;
				}
				
				if(r4==r8){
					nei2 = r4;
				} else {
					nei2 = r4 - 1;
				}
				
				b = nei1 - nei2;
				d1 = nei1 - r8;
			} else {
				return -1;
			}
			
			c = Math.max(c1, a);
			d = Math.max(d1, b);
			
			if(a>=0 && b>=0 && c>=0 && d>=0 && ((a<=2 && b<=5) || (a<=5 && b<=2))){
				dB = Math.min(a, b);
				dW = Math.max(c, d);
				if(dB <= dW){
					return 1;
				} else {
					return -1;
				}
			}
		}
		return -1;
	}
	
	protected int rightSide(int lnkA, int lnkD, int i1A, int i1D, int i2A, int i2D) {
		if((i1A==i2D && i1D==i2A) || (i1A<i2D && lnkA<=i2D && lnkA<=i1A) || (i1A>i2D && lnkA>=i2D && lnkA>=i1A) ||
				(i1D<i2A && lnkD<=i2A && lnkD<=i1D) || (i1D>i2A && lnkD>=i2A &&lnkD>=i1D)){
			return 1;
		} 
		return -1;
	}
	
	protected int rightSide2(int l_A1, int l_D1, int lnkD, int lnkA, int i1A, int i1D, int i2A, int i2D) {
		if((i2A<i1D && lnkA<=i1D && lnkA<=i2A) || (i2A>i1D && lnkA>=i1D && lnkA>=i2A) || 
				(i2D<i1A && lnkD<=i1A && lnkA<=i2D) || (i2D>i1A && lnkD>=i1A && lnkD>=i2D)){
			return 1;
		} else {
			if(i2A==i1D && i2D==i1A) {
				if((lnkD<=i2D && l_A1<=i2D && lnkA>=i2A && l_D1>=i2A) || (lnkD>=i2D && l_A1>=i2D && lnkA<=i2A && l_D1<=i2A)){
					return -1;
				} else {
					return 1;
				}
			}
		}
		return -1;
	}
	
	protected int rightSidePar(int lnkA, int lnkD, int i1A, int i1D, int i2A, int i2D) {
		if((i1A==i2D && i1D==i2A) || (i1A<i2D && lnkA<i2D && lnkA<=i1A && i1D<=i2A && lnkD<=i2A && lnkD<=i1D) ||
				(i1A>i2D && lnkA>i2D && lnkA>=i1A && i1D>=i2A && lnkD>=i2A && lnkD>=i1D) ||
				(i1D<i2A && lnkD<i2A && lnkD<=i1D && i1A<=i2D && lnkA<=i2D && lnkA<=i1A) ||
				(i1D>i2A && lnkD>i2A && lnkD>=i1D && i1A>=i2D && lnkA>=i2D && lnkA>=i1A)){
			return 1;
		}
		return -1;
	}
	
	protected void joinNeighbours(int lnk1A, int r1, int lnk1D, int r2, Pattern nei, Pattern pat, int minDB1, int dB, int minDW1, int dW, int min, int j){
		lnk1A = r1;
		lnk1D = r2;
		nei = pat;
		minDB1 = dB;
		minDW1 = dW;
		min = j;
	}
	
	protected void joinNeighb(Pattern nei, Pattern pat, int minDB2, int dB, int minDW2, int dW){
		nei = pat;
		minDB2 = dB;
		minDW2 = dW;
	}
	
	protected void Alias(int d1, int a1, int d2, int a2, String d1Cn, String a1Cn, String d2Cn, String a2Cn, Pattern pat){
		d1 = pat.getHBond1().getDonor().getResidueNum();
		a1 = pat.getHBond1().getAcceptor().getResidueNum();
		d2 = pat.getHBond2().getDonor().getResidueNum();
		a2 = pat.getHBond2().getAcceptor().getResidueNum();
		d1Cn = pat.getHBond1().getDonor().getChain().getName();
		a1Cn = pat.getHBond1().getAcceptor().getChain().getName();
		d2Cn = pat.getHBond2().getDonor().getChain().getName();
		a2Cn = pat.getHBond2().getAcceptor().getChain().getName();
	}
	
	protected void fillAsnAntiPar(String[] asn1, String[] asn2, Chain c, int cn1, int cn2, ArrayList<Pattern> pat, int NPat){
		int beg1=0, beg2=0, end1=0, end2=0;
		int b1D=0, b1A=0, b2D=0, b2A=0, e1D=0, e1A=0, e2D=0, e2A=0;
		String b1DCn="", b1ACn="", b2DCn="", b2ACn="", e1DCn="", e1ACn="", e2DCn="", e2ACn="", beg1Cn="", beg2Cn="";
		Pattern currPat = new Pattern();
		Pattern prevPat = new Pattern();
		
		for(int i=0; i<NPat; i++){
			if(pat.get(i).getNei1()!=null && pat.get(i).getNei2()==null){
				currPat = pat.get(i).getNei1();
			} else if(pat.get(i).getNei2()!=null && pat.get(i).getNei1()==null){
				currPat = pat.get(i).getNei2();
			} else {
				continue;
			}
			prevPat = pat.get(i);
			
			while(currPat.getNei1()!=null && currPat.getNei2()!=null){
				if((currPat.getNei1().getNei1()==currPat || currPat.getNei1().getNei2()==currPat) && currPat.getNei1()!=prevPat){
					prevPat = currPat;
					currPat = currPat.getNei1();
				} else if((currPat.getNei2().getNei1()==currPat || currPat.getNei2().getNei2()==currPat) && currPat.getNei2()!=prevPat){
					prevPat = currPat;
					currPat = currPat.getNei2();
				} else {
					break;
				}
			}
			
			Alias(b1D,b1A,b2D,b2A,b1DCn,b1ACn,b2DCn,b2ACn,pat.get(i));
		    Alias(e1D,e1A,e2D,e2A,e1DCn,e1ACn,e2DCn,e2ACn,currPat);
		    
		    if( (cn1 != cn2 || e1D - b2A <  e2D - b2A ) &&
		            ( makeEnds(beg1,b1D,b2A,beg1Cn,b1DCn,end1,e2A,e1D,e2ACn,beg2,e2D,e1A,beg2Cn,e2DCn,
		    		   end2,b1A,b2D,b1ACn,pat,NPat)==1 ||
		              makeEnds(beg1,b1D,b2A,beg1Cn,b1DCn,end1,e1D,e2A,e1DCn,beg2,e1A,e2D,beg2Cn,e1ACn,
		    		   end2,b1A,b2D,b1ACn,pat,NPat)==1 ) ){
		    }else if( ( cn1 != cn2 || e2D - b2A <  e1D - b2A ) && 
		            ( makeEnds(beg1,b1D,b2A,beg1Cn,b1DCn,end1,e1A,e2D,e1ACn,beg2,e1D,e2A,beg2Cn,e1DCn,
		    		   end2,b1A,b2D,b1ACn,pat,NPat)==1 ||
		              makeEnds(beg1,b1D,b2A,beg1Cn,b1DCn,end1,e2D,e1A,e2DCn,beg2,e2A,e1D,beg2Cn,e2ACn,
		    		   end2,b1A,b2D,b1ACn,pat,NPat)==1 ) ){
		    } else if( ( cn1 != cn2 || b2A - e1D < b2A - e2D ) && 
		            ( makeEnds(beg1,b1A,b2D,beg1Cn,b1ACn,end1,e2D,e1A,e2DCn,beg2,e2A,e1D,beg2Cn,e2ACn,
		    		   end2,b1D,b2A,b1DCn,pat,NPat)==1 ||
		              makeEnds(beg1,b1A,b2D,beg1Cn,b1ACn,end1,e1A,e2D,e1ACn,beg2,e1D,e2A,beg2Cn,e1DCn,
		    		   end2,b1D,b2A,b1DCn,pat,NPat)==1 ) ){
		    } else if( ( cn1 != cn2 || b2A - e2D < b2A - e1D ) && 
		            ( makeEnds(beg1,b1A,b2D,beg1Cn,b1ACn,end1,e1D,e2A,e1DCn,beg2,e1A,e2D,beg2Cn,e1ACn,
		    		   end2,b1D,b2A,b1DCn,pat,NPat)==1 ||
		              makeEnds(beg1,b1A,b2D,beg1Cn,b1ACn,end1,e2A,e1D,e2ACn,beg2,e2D,e1A,beg2Cn,e2DCn,
		    		   end2,b1D,b2A,b1DCn,pat,NPat)==1 ) ){
		    } else if( ( cn1 != cn2 || b1D - e2A <  b2D - e2A ) && 
		            ( makeEnds(beg1,e1D,e2A,beg1Cn,e1DCn,end1,b2A,b1D,b2ACn,beg2,b2D,b1A,beg2Cn,b2DCn,
		    		   end2,e1A,e2D,e1ACn,pat,NPat)==1 ||
		              makeEnds(beg1,e1D,e2A,beg1Cn,e1DCn,end1,b1D,b2A,b1DCn,beg2,b1A,b2D,beg2Cn,b1ACn,
		    		   end2,e1A,e2D,e1ACn,pat,NPat)==1 ) ){
		    } else if( ( cn1 != cn2 || b2D - e2A <  b1D - e2A ) && 
		            ( makeEnds(beg1,e1D,e2A,beg1Cn,e1DCn,end1,b1A,b2D,b1ACn,beg2,b1D,b2A,beg2Cn,b1DCn,
		    		   end2,e1A,e2D,e1ACn,pat,NPat)==1 ||
		              makeEnds(beg1,e1D,e2A,beg1Cn,e1DCn,end1,b2D,b1A,b2DCn,beg2,b2A,b1D,beg2Cn,b2ACn,
		    		   end2,e1A,e2D,e1ACn,pat,NPat)==1 ) ){
		    } else if( ( cn1 != cn2 || e2A - b1D < e2A - b2D ) && 
		            ( makeEnds(beg1,e1A,e2D,beg1Cn,e1ACn,end1,b2D,b1A,b2DCn,beg2,b2A,b1D,beg2Cn,b2ACn,
		    		   end2,e1D,e2A,e1DCn,pat,NPat)==1 ||
		              makeEnds(beg1,e1A,e2D,beg1Cn,e1ACn,end1,b1A,b2D,b1ACn,beg2,b1D,b2A,beg2Cn,b1DCn,
		    		   end2,e1D,e2A,e1DCn,pat,NPat)==1 ) ){
		    } else if( ( cn1 != cn2 || e2A - b2D < e2A - b1D ) && 
		            ( makeEnds(beg1,e1A,e2D,beg1Cn,e1ACn,end1,b1D,b2A,b1DCn,beg2,b1A,b2D,beg2Cn,b1ACn,
		    		   end2,e1D,e2A,e1DCn,pat,NPat)==1 ||
		              makeEnds(beg1,e1A,e2D,beg1Cn,e1ACn,end1,b2A,b1D,b2ACn,beg2,b2D,b1A,beg2Cn,b2DCn,
		    		   end2,e1D,e2A,e1DCn,pat,NPat)==1 ) ){
		    } else {
		          continue;
		    }
		    
		    if(beg1Cn.equals(c.getName())){
		    	for(int j=beg1; j<=end1; j++){
		    		asn1[j] = "N";
		    	}
		    	for(int j=beg2; j<=end2; j++){
		    		asn2[j] = "N";
		    	} 
		    }else {
		    	for(int j=beg1; j<=end1; j++){
		    		asn2[j] = "N";
		    	}
		    	for(int j=beg2; j<=end2; j++){
		    		asn1[j] = "N";
		    	}
		    }
		    pat.get(i).setNei1(null);
		    pat.get(i).setNei2(null);
		    currPat.setNei1(null);
		    currPat.setNei2(null);
		}
	}
	
	protected void fillAsnPar(String[] asn1, String[] asn2, Chain c, int cn1, int cn2, ArrayList<Pattern> pat, int NPat){
		int beg1=0, beg2=0, end1=0, end2=0;
		int b1D=0, b1A=0, b2D=0, b2A=0, e1D=0, e1A=0, e2D=0, e2A=0;
		String b1DCn="", b1ACn="", b2DCn="", b2ACn="", e1DCn="", e1ACn="", e2DCn="", e2ACn="", beg1Cn="", beg2Cn="";
		Pattern currPat = new Pattern();
		Pattern prevPat = new Pattern();
		
		for(int i=0; i<NPat; i++){
			
			if(pat.get(i).getNei1()!=null && pat.get(i).getNei2()==null){
				currPat = pat.get(i).getNei1();
			} else if(pat.get(i).getNei2()!=null && pat.get(i).getNei1()==null){
				currPat = pat.get(i).getNei2();
			} else {
				continue;
			}
			prevPat = pat.get(i);
			
			while(currPat.getNei1()!=null && currPat.getNei2()!=null){
				if((currPat.getNei1().getNei1()==currPat || currPat.getNei1().getNei2()==currPat) && currPat.getNei1()!=prevPat){
					prevPat = currPat;
					currPat = currPat.getNei1();
				} else {
					prevPat = currPat;
					currPat = currPat.getNei2();
				}
			}
			
			Alias(b1D,b1A,b2D,b2A,b1DCn,b1ACn,b2DCn,b2ACn,pat.get(i));
		    Alias(e1D,e1A,e2D,e2A,e1DCn,e1ACn,e2DCn,e2ACn,currPat);
		    
		    if( ( cn1 != cn2 || Math.abs(e1D-b2A) < Math.abs(e2D-b2A) ) && 
		            ( makeEnds(beg1,b1D,b2A,beg1Cn,b1DCn,end1,e2A,e1D,e2ACn,beg2,b1A,b2D,beg2Cn,b1ACn,
		    		   end2,e2D,e1A,e2DCn,pat,NPat)==1 ||
		              makeEnds(beg1,b1D,b2A,beg1Cn,b1DCn,end1,e1D,e2A,e1DCn,beg2,b1A,b2D,beg2Cn,b1ACn,
		    		   end2,e1A,e2D,e1ACn,pat,NPat)==1 ) ){
		    } else if( ( cn1 != cn2 || Math.abs(e2D-b2A) < Math.abs(e1D-b2A) ) && 
		            ( makeEnds(beg1,b1D,b2A,beg1Cn,b1DCn,end1,e1A,e2D,e1ACn,beg2,b1A,b2D,beg2Cn,b1ACn,
		    		   end2,e1D,e2A,e1DCn,pat,NPat)==1 ||
		              makeEnds(beg1,b1D,b2A,beg1Cn,b1DCn,end1,e2D,e1A,e2DCn,beg2,b1A,b2D,beg2Cn,b1ACn,
		    		   end2,e2A,e1D,e2ACn,pat,NPat)==1 ) ){
		    } else if( ( cn1 != cn2 || Math.abs(b2A-e1D) < Math.abs(b2A-e2D) ) && 
		            ( makeEnds(beg1,b1A,b2D,beg1Cn,b1ACn,end1,e2D,e1A,e2DCn,beg2,b1D,b2A,beg2Cn,b1DCn,
		    		   end2,e2A,e1D,e2ACn,pat,NPat)==1 ||
		              makeEnds(beg1,b1A,b2D,beg1Cn,b1ACn,end1,e1A,e2D,e1ACn,beg2,b1D,b2A,beg2Cn,b1DCn,
		    		   end2,e1D,e2A,e1DCn,pat,NPat)==1 ) ){
		    } else if( ( cn1 != cn2 || Math.abs(b2A-e2D) < Math.abs(b2A-e1D) ) && 
		            ( makeEnds(beg1,b1A,b2D,beg1Cn,b1ACn,end1,e1D,e2A,e1DCn,beg2,b1D,b2A,beg2Cn,b1DCn,
		    		   end2,e1A,e2D,e1ACn,pat,NPat)==1 ||
		              makeEnds(beg1,b1A,b2D,beg1Cn,b1ACn,end1,e2A,e1D,e2ACn,beg2,b1D,b2A,beg2Cn,b1DCn,
		    		   end2,e2D,e1A,e2DCn,pat,NPat)==1 ) ){
		    } else if( ( cn1 != cn2 || Math.abs(b1D-e2A) < Math.abs(b2D-e2A) ) && 
		            ( makeEnds(beg1,e1D,e2A,beg1Cn,e1DCn,end1,b2A,b1D,b2ACn,beg2,e1A,e2D,beg2Cn,e1ACn,
		    		   end2,b2D,b1A,b2DCn,pat,NPat)==1 ||
		              makeEnds(beg1,e1D,e2A,beg1Cn,e1DCn,end1,b1D,b2A,b1DCn,beg2,e1A,e2D,beg2Cn,e1ACn,
		    		   end2,b1A,b2D,b1ACn,pat,NPat)==1 ) ){
		    } else if( ( cn1 != cn2 || Math.abs(b2D-e2A) < Math.abs(b1D-e2A) ) && 
		            ( makeEnds(beg1,e1D,e2A,beg1Cn,e1DCn,end1,b1A,b2D,b1ACn,beg2,e1A,e2D,beg2Cn,e1ACn,
		    		   end2,b1D,b2A,b1DCn,pat,NPat)==1 ||
		              makeEnds(beg1,e1D,e2A,beg1Cn,e1DCn,end1,b2D,b1A,b2DCn,beg2,e1A,e2D,beg2Cn,e1ACn,
		    		   end2,b2A,b1D,b2ACn,pat,NPat)==1 ) ){
		    } else if( ( cn1 != cn2 || Math.abs(e2A-b1D) < Math.abs(e2A-b2D) ) && 
		            ( makeEnds(beg1,e1A,e2D,beg1Cn,e1ACn,end1,b2D,b1A,b2DCn,beg2,e1D,e2A,beg2Cn,e1DCn,
		    		   end2,b2A,b1D,b2ACn,pat,NPat)==1 ||
		              makeEnds(beg1,e1A,e2D,beg1Cn,e1ACn,end1,b1A,b2D,b1ACn,beg2,e1D,e2A,beg2Cn,e1DCn,
		    		   end2,b1D,b2A,b1DCn,pat,NPat)==1 ) ){
		    } else if( ( cn1 != cn2 || Math.abs(e2A-b2D) < Math.abs(e2A-b1D) ) && 
		            ( makeEnds(beg1,e1A,e2D,beg1Cn,e1ACn,end1,b1D,b2A,b1DCn,beg2,e1D,e2A,beg2Cn,e1DCn,
		    		   end2,b1A,b2D,b1ACn,pat,NPat)==1 ||
		              makeEnds(beg1,e1A,e2D,beg1Cn,e1ACn,end1,b2A,b1D,b2ACn,beg2,e1D,e2A,beg2Cn,e1DCn,
		    		   end2,b2D,b1A,b2DCn,pat,NPat)==1 ) ){
		    } else {
		    	continue;
		    }
		    
		    if(beg1Cn.equals(c.getName())){
		    	for(int j=beg1; j<=end1; j++) {
		    		asn1[j] = "P";
		    	}
		    	for(int j=beg2; j<=end2; j++) {
		    		asn2[j] = "P";
		    	}
		    } else {
		    	for(int j=beg1; j<=end1; j++) {
		    		asn2[j] = "P";
		    	}
		    	for(int j=beg2; j<=end2; j++) {
		    		asn1[j] = "P";
		    	}
		    }
		    pat.get(i).setNei1(null);
		    pat.get(i).setNei2(null);
		    currPat.setNei1(null);
		    currPat.setNei2(null);
		}
	}
	
	protected int makeEnds(int beg1, int resBeg1, int neiBeg1, String beg1Cn, String resBeg1Cn, int end1, int resEnd1, int neiEnd1, 
			String resEnd1Cn, int beg2, int resBeg2, int neiBeg2, String beg2Cn, String resBeg2Cn, int end2, int resEnd2, int neiEnd2,
			String resEnd2Cn, ArrayList<Pattern> pat, int NPat){
		boolean flag1=false, flag2=false;
		
		if( resBeg1 <= neiBeg1 && neiBeg1 <= neiEnd1 && neiEnd1 <= resEnd1 && resBeg2 <= neiBeg2 && neiBeg2 <= neiEnd2 && 
				neiEnd2 <= resEnd2 && resBeg1Cn.equals(resEnd1Cn) && resBeg2Cn .equals(resEnd2Cn) ) {
			beg1 = resBeg1;
			end1 = resEnd1;
			beg2 = resBeg2;
			end2 = resEnd2;
			beg1Cn = resBeg1Cn;
			beg2Cn = resBeg2Cn;
			
			for(int i=0; i<NPat && (!flag1 || !flag2); i++){
				HBond patHb1 = pat.get(i).getHBond1();
				HBond patHb2 = pat.get(i).getHBond2();
				if(((patHb1.getDonor().getResidueNum()==beg1) && patHb1.getAcceptor().getResidueNum()==end2 &&
						patHb1.getDonor().getChain().getName().equals(beg1Cn) && patHb1.getAcceptor().getChain().getName().equals(beg2Cn)) ||
						(patHb1.getAcceptor().getResidueNum()==beg1 && patHb1.getDonor().getResidueNum()==end2 &&
						patHb1.getAcceptor().getChain().getName().equals(beg1Cn) && patHb1.getDonor().getChain().getName().equals(beg2Cn)) &&
						patHb1.getDonor().getResidueNum()==patHb2.getAcceptor().getResidueNum() &&
						patHb2.getDonor().getResidueNum()==patHb1.getAcceptor().getResidueNum()){
							flag1 = true;
						}
				
				if(((patHb1.getDonor().getResidueNum()==beg2) && patHb1.getAcceptor().getResidueNum()==end1 &&
						patHb1.getDonor().getChain().getName().equals(beg2Cn) && patHb1.getAcceptor().getChain().getName().equals(beg1Cn)) ||
						(patHb1.getAcceptor().getResidueNum()==beg2 && patHb1.getDonor().getResidueNum()==end1 &&
						patHb1.getAcceptor().getChain().getName().equals(beg2Cn) && patHb1.getDonor().getChain().getName().equals(beg1Cn)) &&
						patHb1.getDonor().getResidueNum()==patHb2.getAcceptor().getResidueNum() &&
						patHb2.getDonor().getResidueNum()==patHb1.getAcceptor().getResidueNum()){
							flag2 = true;
						}
				}
			if(!flag1){
				if(beg1!=neiBeg1) beg1++;
				if(end2!=neiEnd2) end2--;
			}

			if(!flag2){
				if(end1!=neiEnd1) end1--;
				if(beg2!=neiBeg2) beg2++;
			}
			return 1;
		}
		return -1;
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
			findDonor(mol.getChains().get(i));
			findAcceptor(mol.getChains().get(i));
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
	
	protected int noDoubleHBond(){
		int NExcl=0;
		HBond hI = new HBond();
		HBond hJ = new HBond();

		for(int i=0; i<hBond.size()-1; i++){
			for(int j=i+1; j<hBond.size(); j++){
				System.out.println(i);
				hI = hBond.get(i);
				hJ = hBond.get(j);
				
				if(hI.getDonor()!=null && hJ.getDonor()!=null){
					if(hI.getDonor().getResidueNum()==hJ.getDonor().getResidueNum() && 
						hI.getDonor().getChain().getName().equals(hJ.getDonor().getChain().getName()) &&
						hI.getExistPolarInter() && hJ.getExistPolarInter()){
						
						if(Double.compare(hI.getHBondEnergy(), 5.0*hJ.getHBondEnergy()) < 0){
							hJ.setExistPolarInter(false);
							NExcl++;
						} else if(Double.compare(hJ.getHBondEnergy(), 5.0*hI.getHBondEnergy()) < 0){
							hI.setExistPolarInter(false);
							NExcl++;
						}
					}
				}
			}
		}
		return NExcl;
	}
	
	protected void discrPhiPsi(){
		for(int cn=0; cn<mol.getChains().size(); cn++){
			Chain c = mol.getChains().get(cn);
			for(int res=0; res<c.getResidues().size(); res++){
				Residue r = c.getResidues().get(res);
				
				
			}
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
