package predictor.core;

import predictor.util.FileWriter;
import predictor.util.PredictorUtility;
import predictor.core.model.Model;
import predictor.core.model.Molecule;

import org.biojava.bio.structure.Structure;

import predictor.core.model.Chain;
import predictor.core.model.Residue;
import predictor.core.model.secStructures.*;

/***
 * 
 * @author Anson
 *
 * Generates output file and output object
 */
public class OutputController {
	
	protected Model model;
	protected Molecule mol;
	protected String name;
	
	public OutputController (Model model) {
		this.model = model;
		this.mol = model.getMolecules().get(0);
		this.name = mol.getName();
	}
	
	public void output(Structure struc, SSModel ssmodel) {
		System.out.println("REM  --------------------------------------------------------------------  " + name);
		System.out.println("REM                                                                        " + name);
		System.out.println("REM  STRIDE: Knowledge-based secondary structure assignment                " + name);
		System.out.println("REM  Please cite: D.Frishman & P.Argos, Proteins XX, XXX-XXX, 1995         " + name);
		System.out.println("REM                                                                        " + name);
		System.out.println("REM  Residue accessible surface area calculation                           " + name);
		System.out.println("REM  Please cite: F.Eisenhaber & P.Argos, J.Comp.Chem. 14, 1272-1280, 1993 " + name);
		System.out.println("REM               F.Eisenhaber et al., J.Comp.Chem., 1994, submitted       " + name);
		System.out.println("REM                                                                        " + name);
		
		System.out.println("REM  ------------------------ General information -----------------------  " + name);
		System.out.println("REM                                                                        " + name);
		System.out.println("HDR  " + String.format("%-40s", struc.getHeader().get("classification")) + 
				String.format("%-12s", struc.getHeader().get("depDate")) + String.format("%-18s", name) + name);
		
		for(int i=0; i<struc.getCompounds().size(); i++){
			System.out.println("CMP  MOL_ID: " + String.format("%-62s", (i+1)+";") + name);
			System.out.println("CMP   MOLECULE: " + String.format("%-59s", struc.getCompounds().get(i).getMolName()) + name);
			System.out.println("CMP   CHAIN: " + String.format("%-62s", struc.getCompounds().get(i).getChainId()+";") + name);
			System.out.println("CMP   ENGINEERED: " + String.format("%-57s", struc.getCompounds().get(i).getEngineered()+";") + name);
		}
		
		for(int i=0; i<struc.getCompounds().size(); i++){
			System.out.println("SRC  MOL_ID: " + String.format("%-62s", (i+1)+";") + name);
			System.out.println("SRC   ORGANISM_SCIENTIFIC: " + String.format("%-48s", struc.getCompounds().get(i).getOrganismScientific()+";") + name);
		}
		System.out.println("REM                                                                        " + name);
		
		System.out.println("REM  -------------------- Secondary structure summary -------------------  " + name);
		System.out.println("REM                                                                        " + name);		
		for(int i=0; i<mol.getChains().size(); i++){
			int cnt = 1, start = 1;
			Chain c = mol.getChains().get(i);
			System.out.println("CHN  Chain " + c.getName() + "                                                               " + name);
			System.out.println("REM                                                                        " + name);
			
			do{
				// Print markers
				String begin="", content="", end="";
				for(; cnt%50!=0 && cnt<c.getResidues().size(); cnt++){
					if(cnt==start){
						begin = "REM       ";
					} else if((cnt-1)%10==0){
						content += ".";
					} else {
						content += " ";
					}
				}
				if((cnt)%50==0){
					content += " .";
				}
				end = "               " + name;
				System.out.println(begin + String.format("%-50s", content) + end);
				
				// Print residue name
				int j = start;
				int tmp = start;
				begin=""; content=""; end="";
				for(; j<=cnt; j++){
					Residue r = c.getResidues().get(j-1);
					if(j==start){
						begin = "SEQ  " + String.format("%-5s", start);
					}
					content += PredictorUtility.processRName(r.getName());
					
					if(j==cnt){
						end = "   " + String.format("%-12s", cnt) + name;
					}
				}
				System.out.println(begin + String.format("%-50s", content) + end);
				start = j;
				
				// Print residue asn
				begin=""; content=""; end="";
				for(int k=tmp; k<=cnt; k++){
					Residue r = c.getResidues().get(k-1);
					if(k==tmp){
						begin = "STR       ";
					} 
					if(r.getAsn().equals("C")){
						content += " ";
					} else {
						content += r.getAsn();
					}
					if(k==cnt){
						end = "               " + name;
					}
				}
				System.out.println(begin + String.format("%-50s", content) + end);
				System.out.println("REM                                                                        " + name);
				cnt++;
			} while(cnt<=c.getResidues().size());
		}
		System.out.println("REM                                                                        " + name);
//		
//		for(int i=0;i<ssmodel.getSSMol().getSS().size(); i++){
//			SecStructure s = ssmodel.getSSMol().getSS().get(i);
//			System.out.println("LOC  " + s.getAsn()+s.getStartResName() +ssmodel.getSSMol().getSS().size()+ name);
//		}
		System.out.println("REM                                                                        " + name);
		
		
		System.out.println("REM  --------------- Detailed secondary structure assignment-------------  " + name);
		System.out.println("REM                                                                        " + name);
		System.out.println("REM  |---Residue---|    |--Structure--|   |-Phi-|   |-Psi-|                " + name);
		for(int i=0; i<mol.getChains().size(); i++){
			Chain c = mol.getChains().get(i);
			for(int j=0; j<c.getResidues().size(); j++){
				Residue r = c.getResidues().get(j);

				System.out.println("ASG  " + 
						r.getName() + " " + c.getName() + String.format("%5s", r.getResidueSeqNum()) + String.format("%5s",r.getResidueSeqNum()) +
						"    " + r.getAsn() + String.format("%14s", r.getSSName()) +
						"   " + String.format("%7.2f", r.getPhi()) +
						"   " + String.format("%7.2f", r.getPsi()) +
						"                " +
						name);
			}
		}
	}
	
	/**
	 * Output result as txt file
	 * @param fileName
	 * @param filePath
	 */
	public void outputFile(String fileName, String filePath){
		FileWriter.outputAsFile(mol, filePath, fileName, name);
	}
	
	/**
	 * Create Secondary Structure object
	 * @return Secondary structure model
	 */
	public SSModel createSSObject() {
		SSModel ssModel = new SSModel();
		ssModel.create(mol);
		
		return ssModel;
	}
	
}
