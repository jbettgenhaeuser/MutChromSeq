package mutantHunter;

import java.io.File;
import java.io.IOException;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.TransformerException;

import support.CLI;
import support.CLIParseException;

public class Pileup2XML {

	
	
	/**
	 * 
	 * @param args
	 */
	public static void main(String[] args){
		
		

		CLI cli = new CLI();
		cli.parseOptions(args);
		//iowcaf
		String helpString =		"-i <inputPileup>\t\tInput file in mpileup format\n"+
								"-o <output.xml>\t\t\tThe output file in xml format\n"+
								"-w\t\t\t\tThis pileup belongs to the wildtype\n"+
								"-c <int>\t\t\tminimum coverage to report a SNP\n"+
								"-a\n"+
								"-f";
								

		
		try{
			
			if( !cli.hasOption("i") ){
				throw new CLIParseException("Missing option. -i is mandatory.");
			}
			
			if( !cli.hasOption("o") ){
				throw new CLIParseException("Missing option. -o is mandatory.");
			}
			
			if( !cli.hasOption("f") ){
				throw new CLIParseException("Missing option. -f is mandatory.");
			}
			
			
			File inputPileupFile = new File(cli.getArg("i"));
						
			File outputXMLFile = new File(cli.getArg("o"));
			
			File fastaFile = new File(cli.getArg("f"));
			
			
			boolean isWildtype = false;
			if(cli.hasOption("w")){
				isWildtype = true;
			}
			
			int minCoverage = 5;
			if(cli.hasOption("c")){
				try{
					minCoverage = Integer.parseInt(cli.getArg("c"));
				}catch (NumberFormatException e){System.err.println("Warning: wrong input value for -c. This needs to be in int. Default is used.");};	
			}
			
			double maxAlleleFreqeuncy = 0.01;
			if(cli.hasOption("a")){
				try{
					maxAlleleFreqeuncy = Double.parseDouble(cli.getArg("a"));
				}catch (NumberFormatException e){System.err.println("Warning: wrong input value for -a. This needs to be a float. Default is used.");};	
			}
				
			
						
			
			MutChromSeq hunter = new MutChromSeq();
			
			
			
			hunter.setContigListFromFasta(fastaFile);
			System.err.println("Fasta File Read");
			
			String mutantLine = inputPileupFile.getName();
			if( isWildtype){
				mutantLine = MutChromSeq.wildtype;
			}
			hunter.readPileupFile(mutantLine, inputPileupFile, minCoverage,  maxAlleleFreqeuncy);
			
			
			hunter.exportToXML(outputXMLFile);
			
		}catch(IOException e){
			e.printStackTrace();
			System.out.println(helpString);
		}
		
		
		catch(ParserConfigurationException e){
			e.printStackTrace();
			System.out.println(helpString);
		}
		
		catch (TransformerException e){
			e.printStackTrace();
			System.out.println(helpString);
		}
		
		catch (CLIParseException e){
			e.printStackTrace();
			System.out.println(helpString);
		}
	}
	
	
	
	
	
	
}
