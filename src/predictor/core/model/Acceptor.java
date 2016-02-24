package predictor.core.model;

/***
 * 
 * @author Anson
 *
 * Acceptor of H atom in Hydrogen bonding
 */
public class Acceptor {
	
	protected Chain chain;
	protected int res;
	protected Atom dAtomO;
	protected Atom dAtomC;
	protected Atom dAtomCA;
	protected double HBondRadius = 1.90;
	
	public Acceptor(Chain chain, int res, Atom dAtomCA){
		this.chain = chain;
		this.res = res;
		this.dAtomCA = dAtomCA;
	}
	
	public Acceptor(Chain chain, int res, Atom dAtomO, Atom dAtomC, Atom dAtomCA){
		this.chain = chain;
		this.res = res;
		this.dAtomO = dAtomO;
		this.dAtomC = dAtomC;
		this.dAtomCA = dAtomCA;
	}
	
	public void setHBondRadius(double r){
		this.HBondRadius = r;
	}
	
	public double getHBondRadius(){
		return this.HBondRadius;
	}
	
	public void setChain(Chain c){
		this.chain = c;
	}
	
	public Chain getChain(){
		return this.chain;
	}
	
	public void setAcceptorAtomO(Atom O){
		this.dAtomO = O;
	}
	
	public Atom getAcceptorAtomO(){
		return this.dAtomO;
	}
	
	public void setAcceptorAtomC(Atom C){
		this.dAtomC = C;
	}
	
	public Atom getAcceptorAtomC(){
		return this.dAtomC;
	}

	public void setAcceptorAtomCA(Atom CA){
		this.dAtomCA = CA;
	}
	
	public Atom getAcceptorAtomCA(){
		return this.dAtomCA;
	}
	
	public void setResidueNum(int res){
		this.res = res;
	}
	
	public int getResidueNum(){
		return this.res;
	}
}
