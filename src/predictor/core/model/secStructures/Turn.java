package predictor.core.model.secStructures;

import predictor.core.model.Residue;

public class Turn extends SecStructure {
	
	protected String turnType;
	
	public Turn(Residue r) {
		super(r);
		// TODO Auto-generated constructor stub
		this.SSName = "TURN";
		this.asn = "T";
	}
	
	public void setTurnType(String turnType){
		this.turnType = turnType;
	}
	
	public String getTurnType(){
		return this.turnType;
	}
}
