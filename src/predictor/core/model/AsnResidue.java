package predictor.core.model;

import predictor.core.model.secStructures.SecStructureRes;

public class AsnResidue extends Residue implements SecStructureRes {

	public AsnResidue(Residue r, String SSName, String asn){
		this.name = r.getName();
		this.SSName = SSName;
		this.asn = asn;
		this.chainSeqNum = r.getResidueSeqNum();
		this.moleculeSeqNum = r.getMoleculeSeqNum();
		this.chain = r.getParent();
		this.atoms = r.getAtomList();
		this.Phi = r.getPhi();
		this.Psi = r.getPsi();
		this.HBondDnr = r.getHBondDnr();
		this.HBondAcc = r.getHBondAcc();
		this.NBondDnr = r.getNBondDnr();
		this.NBondAcc = r.getNBondAcc();
		this.interchainHBonds = r.getInterChainHBonds();	
	}
	
	public AsnResidue(Residue r){
		this.name = r.getName();
		this.SSName = r.getSSName();
		this.asn = r.getAsn();
		this.chainSeqNum = r.getResidueSeqNum();
		this.moleculeSeqNum = r.getMoleculeSeqNum();
		this.chain = r.getParent();
		this.atoms = r.getAtomList();
		this.Phi = r.getPhi();
		this.Psi = r.getPsi();
		this.HBondDnr = r.getHBondDnr();
		this.HBondAcc = r.getHBondAcc();
		this.NBondDnr = r.getNBondDnr();
		this.NBondAcc = r.getNBondAcc();
		this.interchainHBonds = r.getInterChainHBonds();	
	}
}
