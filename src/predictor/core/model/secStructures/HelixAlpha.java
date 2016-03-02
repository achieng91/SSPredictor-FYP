package predictor.core.model.secStructures;

import predictor.core.model.Residue;

public class HelixAlpha extends SecStructure {
	
	public HelixAlpha(Residue r) {
		super(r);
		this.SSName = "ALPHAHELIX";
		this.asn = "H";
	}
}
