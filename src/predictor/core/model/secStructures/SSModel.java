package predictor.core.model.secStructures;

import java.util.ArrayList;

import predictor.core.model.Chain;
import predictor.core.model.Molecule;
import predictor.core.model.Residue;

public class SSModel {
	
	protected ArrayList<SSMolecule> ssmol;
	
	public void create(Molecule mol){
		for(int i=0; i<mol.getChains().size(); i++){
			Chain c = mol.getChains().get(i);
			ArrayList<SecStructure> ss = new ArrayList<SecStructure>();
			
			for(int j=0; j<c.getResidues().size(); j++){
				Residue r = c.getResidues().get(j);
				if(r.getAsn().equals("H")){
					if(j==0){
						ss.add((HelixAlpha) r);
					} else if(r.getAsn().equals(c.getResidues().get(j-1).getAsn())){
						ss.add((HelixAlpha) r);
					}
				} else if(r.getAsn().equals("G")){
					if(j==0){
						ss.add((HelixThreeTen) r);
					} else if(r.getAsn().equals(c.getResidues().get(j-1).getAsn())){
						ss.add((HelixThreeTen) r);
					}
				} else if(r.getAsn().equals("I")){
					if(j==0){
						ss.add((HelixPi) r);
					} else if(r.getAsn().equals(c.getResidues().get(j-1).getAsn())){
						ss.add((HelixPi) r);
					}
				} else if(r.getAsn().equals("E")){
					if(j==0){
						ss.add((Sheet) r);
					} else if(r.getAsn().equals(c.getResidues().get(j-1).getAsn())){
						ss.add((Sheet) r);
					}
				} else if(r.getAsn().equals("T")){
					if(j==0){
						ss.add((HelixAlpha) r);
					} else if(r.getAsn().equals(c.getResidues().get(j-1).getAsn())){
						Turn t = (Turn) r;
						if(t.getTurnType().equals(((Turn)c.getResidues().get(j-1)).getTurnType())){
							ss.add(t);
						}
					}
				}
				
				if(j<=c.getResidues().size()-1 && !r.getAsn().equals("C")){
					if(!r.getAsn().equals(c.getResidues().get(j+1))){
						ssmol.add(new SSMolecule(ss, r.getSSName(), ss.get(0).getParent().getName(), ss.get(j).getParent().getName()));
						ss = new ArrayList<SecStructure>();
					}
				}
			}
		}
	}
}
