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
