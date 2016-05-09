package predictor.core.model.secStructures;

import java.util.ArrayList;

import predictor.core.model.AsnResidue;

public class SecStructure {

	protected int startResNum, endResNum, secStructurePos;
	protected String startResName, endResName, secStructureType="Coil", secStructureAsn="C", secStructureChainName;
	protected ArrayList<AsnResidue> res = new ArrayList<AsnResidue>();

	
	public SecStructure(AsnResidue res){
		this.res.add(res);
		this.startResNum = res.getResidueSeqNum();
		this.startResName = res.getName();
		this.secStructureChainName = res.getParent().getName();
		this.secStructureType = res.getSSName();
		this.secStructureAsn = res.getAsn();
	}

	public void setStartResNum(int start){
		this.startResNum = start;
	}
	public int getStartResNum(){
		return this.startResNum;
	}
	
	public void setEndResNum(int end){
		this.endResNum = end;
	}
	public int getEndResNum(){
		return this.getEndResNum();
	}
	
	public void setSecStructurePos(int pos){
		this.secStructurePos = pos;
	}
	public int getSecStructurePos(){
		return this.secStructurePos;
	}
	
	public void setStartResName(String name){
		this.startResName = name;
	}
	public String getStartResName(){
		return this.startResName;
	}
	
	public void setEndResName(String name){
		this.endResName = name;
	}
	public String getEndResName(){
		return this.endResName;
	}
	
	public void setType(String type){
		this.secStructureType = type;
	}
	public String getType(){
		return this.secStructureType;
	}
	
	public void setAsn(String asn){
		this.secStructureAsn = asn;
	}
	public String getAsn(){
		return this.secStructureAsn;
	}
	
	public void setChainName(String name){
		this.secStructureChainName = name;
	}
	public String getChainName(){
		return this.secStructureChainName;
	}
	
	public void setResidues(ArrayList<AsnResidue> res){
		this.res = res;
	}
	public ArrayList<AsnResidue> getResidues(){
		return res;
	}
}
