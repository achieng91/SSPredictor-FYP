package predictor.core.model;

public class Pattern {

	protected boolean existPattern;
	protected HBond hb1, hb2;
	protected String type;
	Pattern nei1, nei2;
	
	public Pattern getNei1(){
		return this.nei1;
	}
	
	public void setNei1(Pattern nei1){
		this.nei1 = nei1;
	}
	
	public Pattern getNei2(){
		return this.nei2;
	}
	
	public void setNei2(Pattern nei2){
		this.nei2 = nei2;
	}
	
	public boolean getExistPattern(){
		return this.existPattern;
	}
	
	public void setExistPattern(boolean existPattern){
		this.existPattern = existPattern;
	}
	
	public HBond getHBond1(){
		return this.hb1;
	}
	
	public void setHBond1(HBond hb1){
		this.hb1 = hb1;
	}
	
	public HBond getHBond2(){
		return this.hb2;
	}
	
	public void setHBond2(HBond hb2){
		this.hb2 = hb2;
	}
	
	public String getType(){
		return this.type;
	}
	
	public void setType(String type){
		this.type = type;
	}
}
