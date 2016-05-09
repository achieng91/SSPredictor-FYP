package predictor.core.model.secStructures;

import java.util.ArrayList;

public class SSMolecule {
	
	protected ArrayList<SecStructure> ss = new ArrayList<SecStructure>();

	public void setSS(ArrayList<SecStructure> ss){
		this.ss = ss;
	}
	public ArrayList<SecStructure> getSS(){
		return this.ss;
	}
}
