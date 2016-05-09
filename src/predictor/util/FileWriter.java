package predictor.util;

import java.io.IOException;
import java.io.PrintWriter;

import predictor.core.model.Chain;
import predictor.core.model.Molecule;
import predictor.core.model.Residue;

public abstract class FileWriter {

	public static void outputAsFile(Molecule mol, String filePath, String fileName, String name) {
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
			for(int i=0; i<mol.getChains().size(); i++){
				int cnt = 1, start = 1;
				Chain c = mol.getChains().get(i);
				writer.println("CHN  Chain " + c.getName() + "                                                               " + name);
				writer.println("REM                                                                        " + name);
				
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
					writer.println(begin + String.format("%-50s", content) + end);
					
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
					writer.println(begin + String.format("%-50s", content) + end);
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
					writer.println(begin + String.format("%-50s", content) + end);
					writer.println("REM                                                                        " + name);
					cnt++;
				} while(cnt<=c.getResidues().size());
			}
			writer.println("REM                                                                        " + name);
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
	
}
