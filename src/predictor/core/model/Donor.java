package predictor.core.model;

/***
 * 
 * @author Anson
 *
 * Donor of H atom in Hydrogen bonding
 */
public class Donor {
	
	protected Chain chain;
	protected int res=-1;
	protected Atom dAtomN;
	protected Atom dAtomC;
	protected Atom dAtomCA;
	protected Atom dAtomH;
	protected double HBondRadius = 1.90;
	
	public Donor(Chain chain, int res, Atom dAtomCA) {
		this.chain = chain;
		this.res = res;
		this.dAtomCA = dAtomCA;
	}
	
	public Donor(Chain chain, int res, Atom dAtomCA, Atom dAtomH){
		this.chain = chain;
		this.res = res;
		this.dAtomCA = dAtomCA;
		this.dAtomH = dAtomH;
	}
	
	public Donor(Chain chain, int res, Atom dAtomN, Atom dAtomC, Atom dAtomCA){
		this.chain = chain;
		this.res = res;
		this.dAtomN = dAtomN;
		this.dAtomC = dAtomC;
		this.dAtomCA = dAtomCA;
	}
	
	public Donor(Chain chain, int res, Atom dAtomN, Atom dAtomC, Atom dAtomCA, Atom dAtomH){
		this.chain = chain;
		this.res = res;
		this.dAtomN = dAtomN;
		this.dAtomC = dAtomC;
		this.dAtomCA = dAtomCA;
		this.dAtomH = dAtomH;
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
	
	public void setDonorAtomN(Atom N){
		this.dAtomN = N;
	}
	
	public Atom getDonorAtomN(){
		return this.dAtomN;
	}
	
	public void setDonorAtomC(Atom C){
		this.dAtomC = C;
	}
	
	public Atom getDonorAtomC(){
		return this.dAtomC;
	}

	public void setDonorAtomCA(Atom CA){
		this.dAtomCA = CA;
	}
	
	public Atom getDonorAtomCA(){
		return this.dAtomCA;
	}
	
	public void setDonorAtomH(Atom H){
		this.dAtomH = H;
	}
	
	public Atom getDonorAtomH(){
		return this.dAtomH;
	}
	
	public void setResidueNum(int res){
		this.res = res;
	}
	
	public int getResidueNum(){
		return this.res;
	}
}
