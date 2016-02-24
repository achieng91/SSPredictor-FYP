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
	
	public String getName() {
		return this.name;
	}
}
