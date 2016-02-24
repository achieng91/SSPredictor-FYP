package predictor.core.model;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import predictor.core.model.math.Vector3D;

/**
 * @author Xiu Ting
 *
 */
public class Residue extends predictor.core.model.AbstractParticle {

	protected String name;
	protected int chainSeqNum;
	protected int moleculeSeqNum;
	protected Chain chain;
	protected ArrayList<Atom> atoms;
	protected ArrayList<Integer> representativeCharges;
	protected ArrayList<Vector3D> representativePoints;
	protected double Phi;
	protected double Psi;
	
	public Residue() {
		super();
		name = null;
		atoms = new ArrayList<Atom>();
	}
	
	public double getPhi() {
		return this.Phi;
	}
	
	public void setPhi(double Phi) {
		this.Phi = Phi;
	}
	
	public double getPsi() {
		return this.Psi;
	}
	
	public void setPsi(double Psi) {
		this.Psi = Psi;
	}
	
	public String getName() {
		return this.name;
	}

	public void setName(String name) {
		this.name = name;
	}
	
	public Chain getParent(){
		return chain;
	}
	
	public void setParent(Chain chain) {
		this.chain = chain;
	}
	
	public int getResidueSeqNum() {
		return this.chainSeqNum;
	}

	public void setResidueSeqNum(int seqNum) {
		this.chainSeqNum = seqNum;
	}
	
	public int getChainPosition(){
		return chain.position;
	}
	
	public ArrayList<Atom> getAtomList(){
		return atoms;
	}

	public void setAtomList(List<org.biojava.bio.structure.Atom> atoms) {
		AbstractParticle atom;
		for(int i=0;i<atoms.size();i++){
			atom = new Atom();
			((Atom)atom).setParent(this);
			((Atom)atom).setSymbol(atoms.get(i).getName());
			((Atom)atom).setElementSymbol(atoms.get(i).getElement().toString());
			((Atom)atom).setChainSeqNum(chainSeqNum);
			((Atom)atom).setAtomSeqNum(atoms.get(i).getPDBserial());
			((Atom)atom).setCoordinates(atoms.get(i).getCoords());
			((Atom)atom).setVectorCoord();
			this.atoms.add(((Atom)atom));
		}
	}

	public void setAtomHash(HashMap<String, Atom> atomHash, String modelName) {
		for(int i=0;i<atoms.size();i++){
			atomHash.put(modelName+atoms.get(i).atomSeqNum, atoms.get(i));
		}
	}
	
	public String getModelName(){
		return chain.getModelName();
	}
}