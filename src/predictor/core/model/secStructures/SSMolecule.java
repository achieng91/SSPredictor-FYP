package predictor.core.model.secStructures;

import java.util.ArrayList;

public class SSMolecule {
	
	protected ArrayList<SecStructure> ss;
	protected int startResidue, endResidue;
	protected String startResidueName, endResidueName;
	protected String startChain, endChain;
	protected String ssType;
	
	public SSMolecule(ArrayList<SecStructure> ss, String ssType, String startChain, String endChain) {
		this.ss = ss;
		this.startResidue = ss.get(0).getResidueSeqNum();
		this.endResidue = ss.get(ss.size()-1).getResidueSeqNum();
		this.ssType = ssType;
		this.startChain = startChain;
		this.endChain = endChain;
	}
	
	public int getStartRes() {
		return this.startResidue;
	}
	
	public int getEndRes() {
		return this.endResidue;
	}
	
	public String getStartResName(){
		return this.startResidueName;
	}
	
	public String getEndResName(){
		return this.endResidueName;
	}
	
	public String getStartChain(){
		return this.startChain;
	}
	
	public String getEndChain(){
		return this.endChain;
	}
	
	public String getSSType(){
		return this.ssType;
	}
}
