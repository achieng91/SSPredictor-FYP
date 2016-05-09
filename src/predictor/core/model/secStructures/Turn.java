package predictor.core.model.secStructures;

import predictor.core.model.AsnResidue;

public class Turn extends SecStructure {
	
	protected String turnType;
	
	public Turn(AsnResidue res) {
		super(res);
		// TODO Auto-generated constructor stub
		
		this.setAsn("T");
		this.setType("TURN");
	}
	
	public void setTurnType(String turnType){
		this.turnType = turnType;
	}
	
	public String getTurnType(){
		return this.turnType;
	}
}
