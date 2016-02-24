package predictor.core.model;

public class HBond {
	
	protected Donor donor;
	protected Acceptor acceptor;
	protected double accDonDist, OHDist;
	protected double angleNHO, angleCOH, accAngle, donAngle, accDonAngle, donAccAngle;
	protected double HbondEnergy, Er, Et, Ep, ti, to, p;
	protected boolean existPolarInter=false, existHBondBaker=false, existHBondRose=false;
	
	public void setDonor(Donor d){
		this.donor = d;
	}
	
	public Donor getDonor(){
		return this.donor;
	}
	
	public void setAcceptor(Acceptor a){
		this.acceptor = a;			
	}
	
	public Acceptor getAcceptor(){
		return this.acceptor;
	}
	
	public void setDonAng(double ang){
		this.donAngle = ang;
	}
	
	public double getDonAng(){
		return this.donAngle;
	}
	
	public void setAccDonAng(double ang){
		this.accDonAngle = ang;
	}
	
	public double getAccDonAng(){
		return this.accDonAngle;
	}
	
	public void setDonAccAng(double ang){
		this.donAccAngle = ang;
	}
	
	public double getDonAccAng(){
		return this.donAccAngle;
	}
	
	public void setAccAng(double ang){
		this.accAngle = ang;
	}
	
	public double getAccAng(){
		return this.accAngle;
	}
	
	public void setAngCOH(double ang){
		this.angleCOH = ang;
	}
	
	public double getAngCOH(){
		return this.angleCOH;
	}
	
	public void setAngNHO(double ang){
		this.angleNHO = ang;
	}
	
	public double getAngNHO(){
		return this.angleNHO;
	}
	
	public void setOHDist(double dist){
		this.OHDist = dist;
	}
	
	public double getOHDist(){
		return this.OHDist;
	}
	
	public void setExistPolarInter(boolean b){
		this.existPolarInter = b;
	}
	
	public boolean getExistPolarInter(){
		return this.existPolarInter;
	}
	
	public void setExistHBondBaker(boolean b){
		this.existHBondBaker = b;
	}
	
	public boolean getExistHBondBaker(){
		return this.existHBondBaker;
	}
	
	public void setExistHBondRose(boolean b){
		this.existHBondRose = b;
	}
	
	public boolean getExistHBondRose(){
		return this.existHBondBaker;
	}
	
	public void setAccDonDist(double dist){
		this.accDonDist = dist;
	}
	
	public double getAccDonDist(){
		return this.accDonDist;
	}
	
	public void setHBondEnergy(double e){
		this.HbondEnergy = e;
	}
	
	public double getHBondEnergy() {
		return this.HbondEnergy;
	}
	
	public void setEr(double Er) {
		this.Er = Er;
	}
	
	public double getEr(){
		return this.Er;
	}
	
	public void setEp(double Ep) {
		this.Ep = Ep;
	}
	
	public double getEp(){
		return this.Ep;
	}
	
	public void setEt(double Et) {
		this.Et = Et;
	}
	
	public double getEt(){
		return this.Et;
	}
	
	public void setTi(double ti) {
		this.ti = ti;
	}
	
	public double getTi(){
		return this.ti;
	}
	
	public void setTo(double to){
		this.to = to;
	}
	
	public double getTo(){
		return this.to;
	}
	
	public void setP(double p){
		this.p = p;
	}
	
	public double getP(){
		return this.p;
	}
}
