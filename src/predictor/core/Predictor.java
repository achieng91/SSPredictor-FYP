package predictor.core;

import java.util.ArrayList;

import predictor.core.model.Acceptor;
import predictor.core.model.Chain;
import predictor.core.model.Donor;
import predictor.core.model.HBond;
import predictor.core.model.Model;
import predictor.core.model.Molecule;
import predictor.core.model.Residue;
import predictor.core.model.math.Geometry;
import predictor.util.PredictorUtility;

public class Predictor implements STRIDEController{

	public Model model;
	public Molecule mol;
	public ArrayList<Donor> donors = new ArrayList<Donor>();
	public ArrayList<Acceptor> acceptors = new ArrayList<Acceptor>();
	public ArrayList<HBond> hBond = new ArrayList<HBond>();
	
	public Predictor(Model model){
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
			new STRIDEHelixController().isHelix(mol.getChains().get(i), hBond);
			for(int j=0; j<mol.getChains().size(); j++){
				new STRIDESheetController().isSheet(mol.getChains().get(i), mol.getChains().get(j), i, j, hBond);
				
				new STRIDETurnController().isBetaTurn(mol.getChains().get(i));
				new STRIDETurnController().isGammaTurn(mol.getChains().get(i), hBond);
			}
		}
	}
	
	/**
	 * Determine Hydrogen bonds
	 */
	public void findHBonds() {
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
						hb.setHBondEnergy(STRIDEController.calHBondEnergy(
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
							
							
								hb.setDonAccAng(Math.abs(Geometry.calDihedralAngle(
										don.getDonorAtomN().getPosition(), 
										acc.getAcceptorAtomO().getPosition(), 
										acc.getAcceptorAtomC().getPosition(),
										acc.getAcceptorAtomCA().getPosition())));
								
								if(hb.getDonAccAng() > 90.0 && hb.getDonAccAng() < 270.0){
									hb.setDonAccAng(Math.abs(180.0 - hb.getDonAccAng()));
								}
								
								if(hb.getAccDonAng() <= 60.0 && hb.getDonAccAng() <= 90.0){
									hb.setExistHBondRose(true);
								}
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
	public void findDonor(Chain c){
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
	public void findAcceptor(Chain c){
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
	
	public int noDoubleHBond(){
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
	
	public void discrPhiPsi(){
		for(int cn=0; cn<mol.getChains().size(); cn++){
			Chain c = mol.getChains().get(cn);
			for(int res=0; res<c.getResidues().size(); res++){
				Residue r = c.getResidues().get(res);
				
				
			}
		}
	}
}
