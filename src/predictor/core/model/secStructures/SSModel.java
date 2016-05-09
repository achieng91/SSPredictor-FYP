package predictor.core.model.secStructures;

import java.util.ArrayList;

import predictor.core.model.AsnResidue;
import predictor.core.model.Chain;
import predictor.core.model.Molecule;
import predictor.core.model.Residue;

public class SSModel {
	
	protected SSMolecule ssmol = new SSMolecule();
	
	public void create(Molecule mol){
		for(int i=0; i<mol.getChains().size(); i++){
			Chain c = mol.getChains().get(i);
			ArrayList<SecStructure> ss = new ArrayList<SecStructure>();
			HelixAlpha h = null;
			HelixThreeTen h310 = null;
			HelixPi hPi = null;
			Sheet s = null;
			Bridge b = null;
			Turn t = null;
			SecStructure coil = null;
			
			for(int j=0; j<c.getResidues().size(); j++){
				Residue r = c.getResidues().get(j);
				if(r.getAsn().equals("H")){
					if(j==0){
						h = new HelixAlpha((AsnResidue)r);
					} else if(r.getAsn().equals(c.getResidues().get(j-1).getAsn())){
						if(h==null){
							h = new HelixAlpha((AsnResidue)r);
						} else {
							h.getResidues().add((AsnResidue)r);
						}
					}
				} else if(r.getAsn().equals("G")){
					if(j==0){
						h310 = new HelixThreeTen((AsnResidue)r);
					} else if(r.getAsn().equals(c.getResidues().get(j-1).getAsn())){
						if(h310==null){
							h310 = new HelixThreeTen((AsnResidue)r);
						} else {
							h310.getResidues().add((AsnResidue)r);
						}
					}
				} else if(r.getAsn().equals("I")){
					if(j==0){
						hPi = new HelixPi((AsnResidue)r);
					} else if(r.getAsn().equals(c.getResidues().get(j-1).getAsn())){
						if(hPi==null){
							hPi = new HelixPi((AsnResidue)r);
						} else {
							hPi.getResidues().add((AsnResidue)r);
						}
					}
				} else if(r.getAsn().equals("E")){
					if(j==0){
						s = new Sheet((AsnResidue)r);
					} else if(r.getAsn().equals(c.getResidues().get(j-1).getAsn())){
						if(s==null){
							s = new Sheet((AsnResidue)r);
						} else {
							s.getResidues().add((AsnResidue)r);
						}
					}
				} else if(r.getAsn().equals("T")){

					if(j==0){
						t = new Turn((AsnResidue)r);
						t.setTurnType(r.getTurnType());
					} else if(r.getAsn().equals(c.getResidues().get(j-1).getAsn())){
						if(t==null){
							t = new Turn((AsnResidue)r);
							t.setTurnType(r.getTurnType());
						} else {
							if(r.getTurnType().equals((c.getResidues().get(j-1)).getTurnType())){
								t.getResidues().add((AsnResidue)r);
							} else {
								ss.add(t);
								t = null;
								t = new Turn((AsnResidue)r);
								t.setTurnType(r.getTurnType());
							}
						}
					}
				} else {
					if(j==0){
						coil = new SecStructure((AsnResidue)r);
					} else if(r.getAsn().equals(c.getResidues().get(j-1).getAsn())) {
						if(coil==null){			
							coil = new SecStructure(new AsnResidue(r));
						} else {
							coil.getResidues().add(new AsnResidue(r));
						}
					}
				}
				
				if(j<=c.getResidues().size()-1 && !r.getAsn().equals("C")){
					if(!r.getAsn().equals(c.getResidues().get(j+1))){
						if(h!=null){
							ssmol.getSS().add(h);
							h = null;
						} else if(h310!=null){
							ssmol.getSS().add(h310);
							h310 = null;
						} else if(hPi!=null){
							ssmol.getSS().add(hPi);
							hPi = null;
						} else if(s!=null){
							ssmol.getSS().add(s);
							s = null;
						} else if(b!=null){
							ssmol.getSS().add(b);
							b = null;
						} else if(t!=null){
							ssmol.getSS().add(t);
							t = null;
						} else {
							ssmol.getSS().add(coil);
							coil=null;
						}
					}
				}
			}
		}
	}
	
	public void setSSMol(SSMolecule ssmol){
		this.ssmol = ssmol;
	}
	public SSMolecule getSSMol(){
		return this.ssmol;
	}
}
