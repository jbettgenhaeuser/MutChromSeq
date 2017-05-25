package mutantHunter;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collections;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerConfigurationException;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import support.BioSequence;
import support.CLI;
import support.CLIParseException;
import support.FastaReader;
import support.MPileupLine;





/**
 * 
 * @version 3.0
 * @author steuernb
 *
 */
public class MutChromSeq {

	public static final double version = 3.0;
	public static final String wildtype = "wildtype";
	
	Hashtable<String, TargetContig> targetContigs;
	
	HashSet<String> mutantLines;

	
	
	/* ********************************************* *
	 * ********       Constructors      ************ * 
	 * ********************************************* */
	
	
	
	
	
	
	public MutChromSeq(  ){
		
		this.targetContigs     = new Hashtable<String, TargetContig>();
		
		this.mutantLines = new HashSet<String>();
	}
	
	public MutChromSeq( File wtXML )throws IOException, ParserConfigurationException, SAXException{
		this.targetContigs = new Hashtable<String, TargetContig>();
		this.mutantLines = new HashSet<String>();
		this.addXML(wtXML, true);
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	/* ********************************************* *
	 * ********         Getters         ************ * 
	 * ********************************************* */
	public static double getVersion(){
		return version;
	}
	
	public Hashtable<String, TargetContig> getContigs(){
		return this.targetContigs;
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	/* ********************************************* *
	 * ******** High Level Filters      ************ * 
	 * ********************************************* */
	  
	 /* filter the list of contigs or set regions in contigs*/
	  
	 

	
	
	
	public void filterContigList(HashSet<String> contigs){
		Hashtable<String,TargetContig> h = new Hashtable<String,TargetContig>();
		for(Enumeration<String> myenum = targetContigs.keys(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			if(contigs.contains(key)){
				h.put(key, targetContigs.get(key));
			}
		}
		this.targetContigs = h;
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	public void setContigListFromFasta(File fastaFile)throws IOException{
		FastaReader fastaReader = new FastaReader(fastaFile);
		for (BioSequence seq = fastaReader.readEntry(); seq != null; seq = fastaReader.readEntry()) {
			TargetContig contig = new TargetContig(seq.getIdentifier(), seq.getLength());
			targetContigs.put(seq.getIdentifier(), contig);
		}
		fastaReader.close();
	}
	 
	 
	
	
	
	
	
	
	/**
	 * 
	 *  Read a pileup file and add coverage and snps to TargetContigs that are already existing. 
	 * 
	 * 
	 * @param mutantLine
	 * 			The name of the line this mapping is from. 	
	 * 	
	 * @param pileup
	 * 			The location of the pileup file.
	 * 
	 * 
	 * @throws IOException
	 */
	public void readPileupFile(String mutantLine, File pileup, int minCoverageForSNP, double maxRefAlleleFreq)throws IOException{
		
		
		BufferedReader in;
		
		//check if the fastq is gzipped
		FileInputStream fis = new FileInputStream(pileup);
		byte[] bytes = new byte[2];
		fis.read(bytes);
		int head = ((int) bytes[0] & 0xff) | ((bytes[1] << 8) & 0xff00);
		boolean gzip = GZIPInputStream.GZIP_MAGIC == head;
		fis.close();
		if(gzip){
			in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(pileup))));
		}else{
			in = new BufferedReader((new FileReader(pileup)));
		}
		
		
		
		
		TargetContig contig = this.targetContigs.elements().nextElement();
		int[]coverage = new int[contig.getLength()];
		
		for(String inputline = in.readLine(); inputline != null;inputline = in.readLine() ){
			String substr = inputline.substring(0, Math.min(2000, inputline.length()));
			String[] split = substr.split("\t");
			
			if(!targetContigs.containsKey(split[0])){
				continue;
			}
			
			
			
			if( !split[0].equalsIgnoreCase(contig.getContigName())){
				System.out.println(split[0]);
				contig.addCoverage(mutantLine, coverage);
				targetContigs.put(contig.getContigName(), contig);
				contig = targetContigs.get(split[0]);
				//System.out.println(contig.getContigName());
				coverage = new int[contig.getLength()];
			}	
			
			
			
			int cov = Integer.parseInt(split[3]);
			int pos = Integer.parseInt(split[1]);
					
			coverage[pos-1] =cov ;
			
			Pattern p = Pattern.compile("[ATGC]");
			Matcher m = p.matcher(split[2]);
			
			if(cov> 5 && cov < 2000 && m.matches() ){
				MPileupLine line = new MPileupLine(inputline);
				if( line.getReferenceAlleleFrequency() < maxRefAlleleFreq){
					SNP snp = new SNP(mutantLine, line);
					contig.addSNP(snp);
				}
			}
			
			
			
		}
	
		
		
		
		
	}
	
	

	
	
	
	
	
	
	public void findCandidatesSnpOnly(Vector<String> mutantLines, File outputFile, 
									  int minCoverageToConsiderSNP, double maxReferenceAlleleFrequency,
									  int minNumberOfTotalMutants,
									  int maxAllowedOccurences)throws IOException{


		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		
		System.err.println("DEBUG: Number of Mutant Lines: "+this.mutantLines.size());
		
		
		System.err.println("DEBUG: finding candidates in " + this.getNumberOfContigs() + " contigs");
		System.err.println("DEBUG: writing candidates to " + outputFile.getAbsolutePath());
		Vector<String> contigNames = new Vector<String>();
		
		
		for(Enumeration<String> myenum = targetContigs.keys(); myenum.hasMoreElements();){
		contigNames.add(myenum.nextElement());
		}
		Collections.sort(contigNames, new Comparator<String>(){ 
										public int compare(String o1, String o2){
											if( !o1.startsWith("contig") || !o2.startsWith("contig")){
												return o1.compareTo(o2);
											}else{
												int i1 = Integer.parseInt(o1.split("_")[1]);
												int i2 = Integer.parseInt(o2.split("_")[1]);
												if(i1<i2){
													return -1;
												}else if(i1 == i2){
													return 0;
												}else{
													return 1;
												}
											}
											
											
										}
							});
		
		
		for(Enumeration<String> myenum = contigNames.elements(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			TargetContig contig = targetContigs.get(key);
			//System.out.println(contig.getContigName());
		
			contig.filterByWildtype();
			contig.filterByPosition( mutantLines, maxAllowedOccurences);
			
			if(contig.isCandidateSnpOnly(mutantLines, minCoverageToConsiderSNP, maxReferenceAlleleFrequency,   minNumberOfTotalMutants) ){
					out.write(contig.getReport());
					System.out.println(contig.getReport());
					out.newLine();
			}
		
		
		}
		out.close();
		
	}
		
			
	
	
	

	
	public void findCandidates( File outputFile, int minWildtypeCoverage, 
			double maxReferenceAlleleFrequency,int minCoverageToConsiderSNP, 
			int minNumberOfZeroCoveragePositions,int minNumberOfTotalMutants,
			boolean filterSynonymous)throws IOException{
		
		Vector<String> mutantLines =this.getMutantLines();
		
		findCandidates(mutantLines, outputFile, minWildtypeCoverage, maxReferenceAlleleFrequency, minCoverageToConsiderSNP, minNumberOfZeroCoveragePositions, minNumberOfTotalMutants, filterSynonymous);
		
	}
	
	public void findCandidates(Vector<String> mutantLines, File outputFile, int minWildtypeCoverage, 
								double maxReferenceAlleleFrequency,int minCoverageToConsiderSNP, 
								int minNumberOfZeroCoveragePositions,int minNumberOfTotalMutants,
								boolean filterSynonymous)throws IOException{
		
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		System.err.println("DEBUG: Number of Mutant Lines: "+this.mutantLines.size());
		
		
		System.err.println("DEBUG: finding candidates in " + this.getNumberOfContigs() + " contigs");
		System.err.println("DEBUG: writing candidates to " + outputFile.getAbsolutePath());
		Vector<String> contigNames = new Vector<String>();
		
		
		for(Enumeration<String> myenum = targetContigs.keys(); myenum.hasMoreElements();){
			contigNames.add(myenum.nextElement());
		}
		Collections.sort(contigNames, new Comparator<String>(){ 
													public int compare(String o1, String o2){
														if( !o1.startsWith("contig") || !o2.startsWith("contig")){
															return o1.compareTo(o2);
														}else{
															int i1 = Integer.parseInt(o1.split("_")[1]);
															int i2 = Integer.parseInt(o2.split("_")[1]);
															if(i1<i2){
																return -1;
															}else if(i1 == i2){
																return 0;
															}else{
																return 1;
															}
														}
														
														
													}
										});
		
		
		for(Enumeration<String> myenum = contigNames.elements(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			TargetContig contig = targetContigs.get(key);
			//System.out.println(contig.getContigName());
			
			contig.filterByWildtype( 0.1);
			if(contig.isCandidate(mutantLines, minWildtypeCoverage, maxReferenceAlleleFrequency, minCoverageToConsiderSNP, minNumberOfZeroCoveragePositions, minNumberOfTotalMutants, filterSynonymous) ){
				out.write(contig.getReport());
				System.out.println(contig.getReport());
				out.newLine();
			}
			
			
		}
		out.close();
		
	}
	
	
	public int getNumberOfContigs(){
		return this.targetContigs.size();
	}
	

	public Vector<String> getMutantLines(){
		Vector<String> v = new Vector<String>();
		for(Iterator<String> iterator = mutantLines.iterator(); iterator.hasNext();){
			v.add(iterator.next());
		}
		Collections.sort(v);
		return v;
	}
	
	public void writeTargetContigStatistics(int maxAllowedOccurences, int minCoverage, double maxReferenceAlleleFrequency, File assemblyFasta, File outputFile)throws IOException{
		
		Vector<String> v = new Vector<String>();
		for(Iterator<String> iterator = mutantLines.iterator(); iterator.hasNext();){
			String mutantLine = iterator.next();
			if(mutantLine.equalsIgnoreCase(this.wildtype)){
				continue;
			}
			v.add(mutantLine);
		}
		
		String[] mutantLineArray = new String[v.size()];
		for(int i = 0; i< v.size(); i++){
			mutantLineArray[i] = v.get(i);
		}
		
		
		
		Hashtable<String,Integer> h = new Hashtable<String,Integer>();
		FastaReader fastaReader = new FastaReader(assemblyFasta);
		for (BioSequence seq = fastaReader.readEntry(); seq != null; seq = fastaReader.readEntry()) {
			h.put(seq.getIdentifier(), seq.getNumberOfNonStandardBases());
		}
		fastaReader.close();
		
		
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		out.write("contigName\tlength\tnumMutantLines\tnumMaskedBases\tnumBasesUsed");
		for( int i = 0; i< mutantLineArray.length; i++){
			out.write("\t" + mutantLineArray[i]);
		}
		out.newLine();
		
		for(Enumeration<String> myenum = this.targetContigs.keys(); myenum.hasMoreElements();){
			String contigName = myenum.nextElement();
			
			
			TargetContig targetContig = this.targetContigs.get(contigName);
			
			
			
			targetContig.filterByPosition(v, maxAllowedOccurences);
			targetContig.filterByWildtype();
			targetContig.deleteSNPsByCoverage(minCoverage);
			targetContig.deleteSNPsByFrequency(maxReferenceAlleleFrequency);
			
			int numMutants = 0;
			for(Enumeration<String> myenum2 = v.elements(); myenum2.hasMoreElements();){
				String mutantLine = myenum2.nextElement();
				if( targetContig.snps.containsKey(mutantLine)  && targetContig.snps.get(mutantLine).size() > 0 ){
					numMutants++;
				}
			}
			
			int numReadBases =  (targetContig.length-h.get(contigName).intValue());
			
			
			
			out.write(contigName + "\t" + targetContig.getLength() + "\t" + numMutants +"\t" + h.get(contigName).intValue() + "\t" +numReadBases  );
			
			
			for( int i = 0; i< mutantLineArray.length; i++){
				if( targetContig.snps.containsKey(mutantLineArray[i])){
					out.write("\t" +  targetContig.snps.get(mutantLineArray[i]).size() );
				}else{
					out.write("\t0");
				}
			}
			
			out.newLine();
			
		}
		
		out.close();
		
		
	}
	
	public  void getNumberOfTargetContigs( int minContigLength, int maxAllowedOccurences, int minCoverage, double maxReferenceAlleleFrequency){
		
		
		
		Vector<String> v = new Vector<String>();
		for(Iterator<String> iterator = mutantLines.iterator(); iterator.hasNext();){
			String mutantLine = iterator.next();
			if(mutantLine.equalsIgnoreCase(this.wildtype)){
				continue;
			}
			v.add(mutantLine);
		}
		
		
		
		
		
		int[] numLines = new int[v.size()+1]; //index is the number of mutantLines mutated, value is the number of contigs.
		
		
		for(Enumeration<String> myenum = this.targetContigs.keys(); myenum.hasMoreElements();){
			String contigName = myenum.nextElement();
			
			
			TargetContig targetContig = this.targetContigs.get(contigName);
			
			if(targetContig.length<minContigLength){
				continue;
			}
			
			targetContig.filterByPosition(v, maxAllowedOccurences);
			targetContig.filterByWildtype();
			targetContig.deleteSNPsByCoverage(minCoverage);
			targetContig.deleteSNPsByFrequency(maxReferenceAlleleFrequency);
			
			int numMutants = 0;
			for(Enumeration<String> myenum2 = v.elements(); myenum2.hasMoreElements();){
				String mutantLine = myenum2.nextElement();
				if( targetContig.snps.containsKey(mutantLine)  && targetContig.snps.get(mutantLine).size() > 0 ){
					numMutants++;
				}
			}
			
			numLines[numMutants] = numLines[numMutants] + 1;
			
		}
		
		
		
		for( int i = 0; i< numLines.length; i++){
			System.out.println("Contigs with " + i + " mutations: "+ numLines[i]);
		}
		
	}
	
	
	
	
	
	

	/* ******************************************************** *
	 * ********  export and import function        ************ * 
	 * ******************************************************** */
	
	
	
	/**
	 * Export the contig list to an xml file.
	 * 
	 * @param xmlFile
	 * @throws IOException
	 * @throws ParserConfigurationException
	 * @throws TransformerConfigurationException
	 * @throws TransformerException
	 */
	public void exportToXML(File xmlFile)throws IOException, ParserConfigurationException, TransformerConfigurationException, TransformerException{
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		DocumentBuilder db = dbf.newDocumentBuilder();
		Document dom = db.newDocument();
		
		Element rootElement = dom.createElement("MutantHunter_old");
		rootElement.setAttribute("version", version+"");
		dom.appendChild(rootElement);
		
		
		
		Element contigsElement = dom.createElement("TargetContigs");
		rootElement.appendChild(contigsElement);
		for(Enumeration<String> myenum1 = targetContigs.keys(); myenum1.hasMoreElements();){
			String key = myenum1.nextElement();
			TargetContig contig = targetContigs.get(key);
			if(!contig.hasNoData()){
				contigsElement.appendChild(contig.getXMLElement(dom));
			
			}
		}
		
		
		DOMSource source = new DOMSource(rootElement) ;
		StreamResult result = new StreamResult(xmlFile);
		Transformer transformer = TransformerFactory.newInstance().newTransformer();
		if (dom.getDoctype() != null) {
		    String systemValue = (new File (dom.getDoctype().getSystemId())).getName();
		    transformer.setOutputProperty(OutputKeys.DOCTYPE_SYSTEM, systemValue);
		}
		transformer.setOutputProperty(OutputKeys.INDENT, "yes"); //without this there is no newline after elements.
		
		transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "2");  //not sure what happens here but this makes the nice whitespace hierarchy in the output txt.
		transformer.transform(source, result);
		
	}
	
	
	
	/**
	 * 
	 * Read an xml file and add the data to this object. This can be used tobuild the list of NLR_Contigs or just to add data vrom additional mutant lines
	 * 
	 * @param xmlFile
	 * 			The location of the file
	 * @param addNewContigs
	 * 			if true, additional contigs that do not already exist will be added.
	 * 			if false, only new information will be added to existing NLR_Contigs.
	 * 
	 * @throws IOException
	 * @throws SAXException
	 * @throws ParserConfigurationException
	 */
	public void addXML(File xmlFile, boolean addNewContigs)throws IOException, SAXException, ParserConfigurationException{
				
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		DocumentBuilder db = dbf.newDocumentBuilder();
		Document dom = db.parse(xmlFile);
		
		
		Element rootElement = dom.getDocumentElement();
		
		
		Element nlrContigsElement = (Element) rootElement.getElementsByTagName("NBLRR_Contigs").item(0);
		if( nlrContigsElement != null){
			NodeList nodeList = nlrContigsElement.getElementsByTagName("TargetContig");
			for( int i = 0; i< nodeList.getLength(); i++){
				Element contigElement = (Element) nodeList.item(i);
				
				
				TargetContig contig = new TargetContig(contigElement);
				
				
				//this is until I get the data mass under control Do not store the coverage
				if(!contig.coverages.keys().nextElement().equalsIgnoreCase(this.wildtype)){
					contig.removeCoverage(contig.coverages.keys().nextElement());
				}
				
				
				
				
				//System.out.println(contig.getContigName());
				if( targetContigs.containsKey(contig.getContigName())){
					targetContigs.get(contig.getContigName()).mergeContig(contig);
				}else{
					if( addNewContigs){
						targetContigs.put(contig.getContigName(), contig);
					}	
				}
				
				
			}
		}
		
		for(Enumeration<String> myenum = this.targetContigs.keys(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			for(Enumeration<String> myenum2 = this.targetContigs.get(key).getMutantLineNames().elements(); myenum2.hasMoreElements();){
				mutantLines.add(myenum2.nextElement());
			}
		}
		
	}
	
	
	
	
	public void addXMLQuick(File xmlFile, boolean addNewContigs)throws IOException{
		BufferedReader in = new BufferedReader(new FileReader(xmlFile));
		int counta = 0;
		
		TargetContig contig = null;
		
		
		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			
			
			
			if(inputline.contains("/TargetContig>")){
				if(contig != null){
					this.targetContigs.put(contig.getContigName(), contig);
					counta++;
					if(counta%10000==0){
						System.out.print(".");
						if(counta%1000000==0){
							System.out.println(counta);
						}
					}
				}	
				contig = null;
			}	
			
			
			if( inputline.contains("<TargetContig ")){
				
				String name = inputline.split("name=\"")[1].split("\"")[0];
				
				if(this.targetContigs.containsKey(name) ){
					contig = this.targetContigs.get(name);
					
				}else if(addNewContigs){
					
					int length = Integer.parseInt(inputline.split("length=\"")[1].split("\"")[0]);
					contig = new TargetContig(name, length);
				}	
				
			}	
			
			if(inputline.trim().startsWith("<SNP ")&& contig != null ){
				String referenceAllele = inputline.split("referenceAllele=\"")[1].split("\"")[0];
				if(!referenceAllele.equalsIgnoreCase("N")){
					int alleleA = Integer.parseInt(inputline.split("alleleA=\"")[1].split("\"")[0] );
					int alleleT = Integer.parseInt(inputline.split("alleleT=\"")[1].split("\"")[0] );
					int alleleG = Integer.parseInt(inputline.split("alleleG=\"")[1].split("\"")[0] );
					int alleleC = Integer.parseInt(inputline.split("alleleC=\"")[1].split("\"")[0] );
					int position = Integer.parseInt(inputline.split("position=\"")[1].split("\"")[0] );
					int coverage = Integer.parseInt(inputline.split("coverage=\"")[1].split("\"")[0] );
					String mutantLine  = inputline.split("mutantLine=\"")[1].split("\"")[0] ;
					this.mutantLines.add(mutantLine);
					SNP snp = new SNP(mutantLine, contig.getContigName(), referenceAllele.toCharArray()[0], position, coverage, alleleA, alleleC, alleleG, alleleT);
					contig.addSNP(snp);
				}	
				
			}	
			
			
		}

		in.close();
		
		System.out.println();
	}
	
	
	
	
	
	
	
	
	
	
	
	
	public static void main(String[] args){

		
		CLI cli = new CLI();
		cli.parseOptions(args);
		
		
		
		
		try{
			
			if( !cli.hasOption("w") || !cli.hasOption("m")  || !cli.hasOption("o")){
				throw new CLIParseException("Missing required option: -w, -m, -b and -o are required.");
			}
			
			File wtFile = new File(cli.getArg("w"));
			Vector<String> mutants = cli.getArgs("m");
			File outputFile = new File(cli.getArg("o"));
			
			
			MutChromSeq mutChromSeq = new MutChromSeq();
			mutChromSeq.addXMLQuick(wtFile, true);
			for(Enumeration<String> myenum = mutants.elements(); myenum.hasMoreElements();){
				String mutant = myenum.nextElement();
				mutChromSeq.addXMLQuick(new File(mutant),false);
			}
			
			Vector<String> v = new Vector<String>();
			for(Iterator<String> iterator = mutChromSeq.mutantLines.iterator(); iterator.hasNext();){
				String s = iterator.next();
				if(!s.equalsIgnoreCase(MutChromSeq.wildtype)){
					v.add(s);
				}
				
			}
			
			
			int minCoverageToConsiderSNP = 15;
			double maxReferenceAlleleFrequency = 0.01;
			int minNumberOfTotalMutants = 5;
			int maxAllowedOccurences = 2;
			
			
			if( cli.hasOption("n")){
				try{
					minNumberOfTotalMutants = Integer.parseInt(cli.getArg("n"));
				}catch(NumberFormatException e){
					throw new CLIParseException("Argument for -n has to be an int");
				}
			}
			
			
			if( cli.hasOption("c")){
				try{
					minCoverageToConsiderSNP = Integer.parseInt(cli.getArg("c"));
				}catch(NumberFormatException e){
					throw new CLIParseException("Argument for -c has to be an int");
				}
			}
			
					
			
			if( cli.hasOption("a")){
				try{
					maxReferenceAlleleFrequency = Double.parseDouble(cli.getArg("a"));
				}catch(NumberFormatException e){
					throw new CLIParseException("Argument for -a has to be an float");
				}
			}
			
			
			if( cli.hasOption("z")){
				try {
					maxAllowedOccurences = Integer.parseInt(cli.getArg("z"));
				} catch (NumberFormatException e) {
					throw new CLIParseException("Argument for -z has to be an int");
				}
			}
			
			
			
			
			
			
			
			mutChromSeq.findCandidatesSnpOnly(v, outputFile ,minCoverageToConsiderSNP, maxReferenceAlleleFrequency,  minNumberOfTotalMutants,maxAllowedOccurences);
			
			
			
			
			
			
			
			
			
			
			
			
		}catch (CLIParseException e){
			e.printStackTrace();
			
			String s = "-w <wt.xml>\t\t\tThe XML file generated with Pileup2XML.jar made from wildtyp\n"+
					   "-m <mt1.xml [mtn.xml]*\t\tThe XML files generated with Pileup2XML.jar made from mutants\n"+
					   "-o <output.txt>\t\t\tOutput file with andidates\n"+
					   "-n <int>\t\t\tMinimum number of mutants to report a contig. Default is 2\n"+
					   "-c <int>\t\t\tMininum coverage for mappings to be regarded. Default is 10\n"+
					   "-a <float>\t\t\tMaximum reference allele frequency to consider a SNP. Default is 0.01\n"+
					   "-z <int>\t\t\tNumber mutant lines that are allowed to have SNV in same position. Default is 2\n";
			
			System.err.println(s);
			
		}
		
		catch (IOException e){
			e.printStackTrace();
		}
		
		
		
		
	}
	
	
	
	
	
	
	
	
}
