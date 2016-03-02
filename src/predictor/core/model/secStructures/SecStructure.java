package predictor.core.model.secStructures;

import predictor.core.model.Residue;

/***
 * 
 * @author Anson
 *
 * Secondary Structure model
 *
 */
public abstract class SecStructure extends Residue {
	
	protected String name = "DEFAULT";
	
	public SecStructure(Residue r){
		this.name = r.getName();
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
	
	public String getName() {
		return this.name;
	}
}
