package predictor.core;

import java.util.ArrayList;

import predictor.core.model.AsnResidue;
import predictor.core.model.Chain;
import predictor.core.model.HBond;
import predictor.core.model.Pattern;
import predictor.core.model.Residue;
import predictor.core.model.math.PhiPsiMap;

/***
 * Determine Sheet structure
 * @author Anson
 *
 */
public class STRIDESheetController {
	
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
					link(hb, c1, c2, cn1, cn2, r1, r3, r3, r1, r1, r3, STRIDEController.THRESHOLD_E1, patN, patCntN, "1331");
				}
				if((i+2)<c1.getResidues().size() && ((cn1!=cn2 && j-2>=0) || (j-2)-(i+2)>=2)){
					link(hb, c1, c2, cn2, cn1, r3, r1, r2, r4, rB, rA, STRIDEController.THRESHOLD_E1, patN, patCntN, "3124");
				}
				if(((cn1!=cn2 && (j-1)>=0) || (j-1)-i>4) && ((j-2)>=c1.getResidues().size() || (cn1==cn2 && j-(i-1)<=4) || 
						link(hb, c1, c2, cn1, cn2, r1, r3, r3, rA, null, r3, STRIDEController.THRESHOLD_E1, patN, patCntN, "133A")==-1) && (i+1<0 || 
						link(hb, c1, c2, cn1, cn2, r1m1, rB, rB, r1, r1, null, STRIDEController.THRESHOLD_E1, patN, patCntN, "1-BB1")==-1)){
					link(hb, c1, c2, cn1, cn2, r1, r3, rB, r1, r1, null, STRIDEController.THRESHOLD_E1, patN, patCntN, "13B1");
				}
				if(((i+1)<c1.getResidues().size() && (cn1!=cn2 || j-(i+1)>4)) && ((cn1==cn2 && (j-1)-i<=4) || (cn1!=cn2 && j-1<0) || 
						link(hb, c1, c2, cn1, cn2, r1, r3, rB, r1, r1, null, STRIDEController.THRESHOLD_E1, patN, patCntN, "13B1")==-1) && (j+1>=c2.getResidues().size() || 
						link(hb, c1, c2, cn1, cn2, rA, r3p1, r3, rA, rA, null, STRIDEController.THRESHOLD_E1, patN, patCntN, "A3+3A")==-1)){
					link(hb, c1, c2, cn1, cn2, r1, r3, r3, rA, null, r3, STRIDEController.THRESHOLD_E1, patN, patCntN, "133A");
				}
				
				if((cn1==cn2 && Math.abs(j-i)<=3) || (j+2)>=c2.getResidues().size()){
					continue;
				}
				rB = c2.getResidues().get(j+1);
				r4 = c2.getResidues().get(j+2);
				
				if((i+1)<c1.getResidues().size() && (cn1!=cn2 || Math.abs((i+1)-j)>3)){
					link(hb, c1, c2, cn2, cn1, r3, r1, r2, r3, r3, rA, STRIDEController.THRESHOLD_E2, patP, patCntP, "3123");
				}
				if((j+2)<c2.getResidues().size() && (cn1!=cn2 || Math.abs((j+2)-i)>3)){
					link(hb, c1, c2, cn1, cn2, r1, r3, r4, r1, r1, rB, STRIDEController.THRESHOLD_E2, patP, patCntP, "1341");
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
				AsnResidue s = new AsnResidue(r.get(i), "SHEET", "E");
				r.add(i,s); r.remove(i+1);
			} else if(antiPar1[i].equals("B") || par1[i].equals("B")){
				AsnResidue b = new AsnResidue( r.get(i), "BRIDGE", "B");
				r.add(i,b); r.remove(i+1);
			} else if(antiPar1[i].equals("b") || par1[i].equals("b")){
				AsnResidue b = new AsnResidue( r.get(i), "BRIDGE", "B");
				b.setAsn("b");
				r.add(i,b); r.remove(i+1);
			}
		}
		
		for(int i=0; i<c2.getResidues().size(); i++){
			ArrayList<Residue> r = c2.getResidues();
			if(r.get(i).getAsn().equals("E")){
				continue;
			} else if(antiPar2[i].equals("N") || par2[i].equals("P")){
				AsnResidue s = new AsnResidue(r.get(i), "SHEET", "E");
				r.add(i,s); r.remove(i+1);
			} else if(antiPar2[i].equals("B") || par2[i].equals("B")){
				AsnResidue b = new AsnResidue( r.get(i), "BRIDGE", "B");
				r.add(i,b); r.remove(i+1);
			} else if(antiPar2[i].equals("b") || par2[i].equals("b")){
				AsnResidue b = new AsnResidue( r.get(i), "BRIDGE", "B");
				b.setAsn("b");
				r.add(i,b); r.remove(i+1);
			}
		}
	}
	
	protected int link(ArrayList<HBond> hb, Chain c1, Chain c2, int cn1, int cn2, Residue r1_1, Residue r1_2, Residue r2_2, Residue r2_1, 
			Residue cR1, Residue cR2, double threshold, ArrayList<Pattern> pat, int numPat, String text) {
		int bondNumber1, bondNumber2;
		double conf, coeff, prob1, prob2;
		
		bondNumber1 = STRIDEController.findPolInt(hb, r1_1, r1_2);
		bondNumber2 = STRIDEController.findPolInt(hb, r2_2, r2_1);
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
		coeff = 1 + STRIDEController.C1_E + STRIDEController.C2_E*conf;
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
}
