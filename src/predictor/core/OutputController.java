package predictor.core;

import predictor.util.FileWriter;

import predictor.core.model.Model;
import predictor.core.model.Molecule;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.Writer;

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
	
	public void output() {
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
					content += processRName(r.getName());
					
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
	
	protected char processRName(String name){
		switch(name) {
		case "ALA": return 'A';	
		case "ARG": return 'R';
		case "ASN": return 'N';
		case "ASP": return 'D';
		case "ASX": return 'B';
		case "CYS": return 'C';
		case "GLN": return 'Q';
		case "GLU": return 'E';
		case "GLX": return 'Z';
		case "GLY": return 'G';
		case "HIS": return 'H';
		case "ILE": return 'I';
		case "LEU": return 'L';
		case "LYS": return 'K';
		case "MET": return 'M';
		case "PRO": return 'P';
		case "PHE": return 'F';
		case "SER": return 'S';
		case "THR": return 'T';
		case "TRP": return 'W';
		case "TYR": return 'Y';
		case "VAL": return 'V';
		default: return 'X';
		}
	}
	
	public void outputFile(String fileName, String filePath){
		//TODO: Shift to FileWriter
		try {
			PrintWriter writer = new PrintWriter(filePath + fileName + ".txt", "UTF-8");
			
			writer.println("REM  --------------------------------------------------------------------  " + name);
			writer.println("REM                                                                        " + name);
			writer.println("REM  STRIDE: Knowledge-based secondary structure assignment                " + name);
			writer.println("REM  Please cite: D.Frishman & P.Argos, Proteins XX, XXX-XXX, 1995         " + name);
			writer.println("REM                                                                        " + name);
			writer.println("REM  Residue accessible surface area calculation                           " + name);
			writer.println("REM  Please cite: F.Eisenhaber & P.Argos, J.Comp.Chem. 14, 1272-1280, 1993 " + name);
			writer.println("REM               F.Eisenhaber et al., J.Comp.Chem., 1994, submitted       " + name);
			writer.println("REM                                                                        " + name);
			
			writer.println("REM  ------------------------ General information -----------------------  " + name);
			writer.println("REM                                                                        " + name);
			
			writer.println("REM  -------------------- Secondary structure summary -------------------  " + name);
			writer.println("REM                                                                        " + name);		
			
			writer.println("REM  --------------- Detailed secondary structure assignment-------------  " + name);
			writer.println("REM                                                                        " + name);
			writer.println("REM  |---Residue---|    |--Structure--|   |-Phi-|   |-Psi-|                " + name);
			for(int i=0; i<mol.getChains().size(); i++){
				Chain c = mol.getChains().get(i);
				for(int j=0; j<c.getResidues().size(); j++){
					Residue r = c.getResidues().get(j);

					writer.println("ASG  " + 
							r.getName() + " " + c.getName() + String.format("%5s", r.getResidueSeqNum()) + String.format("%5s",r.getResidueSeqNum()) +
							"    " + r.getAsn() + String.format("%14s", r.getSSName()) +
							"   " + String.format("%7.2f", r.getPhi()) +
							"   " + String.format("%7.2f", r.getPsi()) +
							"                " +
							name);
				}
			}
			writer.close();
		} catch (IOException ex) {
			  ex.printStackTrace();
		} 
	}
	
	/**
	 * Create Secondary Structure object
	 * @return Secondary structure model
	 */
	public SSModel createSSObject() {
		//TODO: create SS Object
		return new SSModel();
	}
	
}
