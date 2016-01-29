package mutantHunter;

import java.util.Collections;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Vector;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

public class TargetContig {

	String contigName;
	int length;
	Hashtable<String,Hashtable<Integer,SNP>> snps; //<mutantLine, <position, snp>>
	Hashtable<String,int[]> coverages;
	
	Hashtable<String,Vector<int[]>> lowCovRegions;
	
		
	
	public final String wildtype = "wildtype";
	String report;
	
	
	
	
	
	
	
	
	
	
	
	
	/* ********************************************* *
	 * ********       Constructors      ************ * 
	 * ********************************************* */
	
	
	
	
	
	
	
	public TargetContig(String contigName, int length){
		this.contigName = contigName;
		this.snps=new Hashtable<String,Hashtable<Integer,SNP>>();
		this.coverages = new Hashtable<String,int[]>();
		this.length = length;
	}
	
	public TargetContig(Element xmlElement){
		this.contigName = xmlElement.getAttribute("name");
		
		try{
			this.length = Integer.parseInt(xmlElement.getAttribute("length"));
		}catch(NumberFormatException e){}
		this.snps = new Hashtable<String,Hashtable<Integer, SNP>>();
		this.coverages = new Hashtable<String,int[]>();
		
		
		
		
		//SNPs
		Element snpsElement = (Element) xmlElement.getElementsByTagName("SNPs").item(0);
		if( snpsElement != null){
			NodeList nodeListSNPs = snpsElement.getElementsByTagName("SNP");
			for( int i = 0; i< nodeListSNPs.getLength(); i++){
				Element snpElement = (Element) nodeListSNPs.item(i);
				String mutantLine = snpElement.getAttribute("mutantLine");
				int position = Integer.parseInt(snpElement.getAttribute("position"));
				SNP snp = new SNP(snpElement);
				
				Hashtable<Integer, SNP> h = new Hashtable<Integer, SNP>();
				if(this.snps.containsKey(mutantLine)){
					h = this.snps.get(mutantLine);
				}
				h.put(new Integer(position), snp);
				this.snps.put(mutantLine, h);
				
				
			}
		}
		
		
		
		
		
		//coverages
		Element coveragesElement = (Element) xmlElement.getElementsByTagName("Coverages").item(0);
		if(coveragesElement != null){
			NodeList nodeList = coveragesElement.getElementsByTagName("Coverage");
			for( int i = 0; i< nodeList.getLength(); i++){
				Element coverageElement = (Element) nodeList.item(i);
				String mutantLine = coverageElement.getAttribute("name");
				String[] a = coverageElement.getAttribute("value").split(",");
				int[] cov = new int[ a.length];
				for( int j = 0; j< a.length; j++){
					
					cov[j] = Integer.parseInt(a[j]);
				}
				this.coverages.put(mutantLine, cov);
			}
			
		}
		
		
		
		
	}
	
	
	
	
	

	/* ********************************************* *
	 * ********       Getters           ************ * 
	 * ********************************************* */
	
	
	
	
	
	
	public Vector<String> getMutantLineNames(){
		
		//sort the lines
		HashSet<String> hashset = new HashSet<String>();
		for(Enumeration<String> myenum = snps.keys(); myenum.hasMoreElements();){
			hashset.add(myenum.nextElement());
		}
		for(Enumeration<String> myenum = coverages.keys(); myenum.hasMoreElements();){
			hashset.add(myenum.nextElement());
		}
		hashset.remove(wildtype);
		Vector<String> mutantLines = new Vector<String>();
		for(Iterator<String> iterator = hashset.iterator();iterator.hasNext();){
			mutantLines.add(iterator.next());
		}
		Collections.sort(mutantLines);
		
		return mutantLines;
	}
	
	
	
	
	public boolean hasNoData(){
		if( coverages.size() == 0 && snps.size() == 0){
			return true;
		}
		return false;
	}
	
	public int getMedianCoverage(String mutantLine){
		if( this.coverages.get(mutantLine) ==null){
			return 0;
		}
		Vector<Integer> v = new Vector<Integer>();
		for( int i = 0; i< this.coverages.get(mutantLine).length; i++){
			v.add(this.coverages.get(mutantLine)[i]);
		}
		Collections.sort(v);
		if( v.size()>0){
			return v.get(v.size()/2).intValue();
		}else{
			return 0;
		}
		
	}
	
	public String getReport(){
		return this.report;
	}

	
	
	public int getLength(){
		if(this.length == 0){
			return this.coverages.get(this.wildtype).length;
		}
		
		
		return this.length;
	}
	/* ********************************************* *
	 * ********       Add methods       ************ * 
	 * ********************************************* */
	
	
	public void setLength(int length){
		this.length = length;
	}
	
	
	
	/**
	 * Add coverage information for this contig from a mutant line. 
	 * 	 * 
	 * @param mutantLine
	 * @param coverage
	 */
	public void addCoverage(String mutantLine, int[] coverage){
		this.coverages.put(mutantLine, coverage);
	}
	
	
	
	
	/**
	 * Add a snp for this contig between the mutant line and the wildtype.
	 * 
	 * @param mutantLine
	 * @param snp
	 */
	public void addSNP( SNP snp){
		Hashtable<Integer, SNP> h = new Hashtable<Integer, SNP>();
		String mutantLine = snp.getMutantLine();
		if(snps.containsKey(mutantLine)){
			h = snps.get(mutantLine);
		}
		h.put(snp.getPosition(), snp);
		this.snps.put(mutantLine, h);
	}

	
	public boolean mergeContig(TargetContig contig2){
		if( this.contigName.equalsIgnoreCase(contig2.contigName)){
			
			
			
			for(Enumeration<String> myenum = contig2.coverages.keys(); myenum.hasMoreElements();){
				String mutantLine = myenum.nextElement();
				if( !this.coverages.containsKey(mutantLine)){
					this.coverages.put(mutantLine, contig2.coverages.get(mutantLine));
				}
			}
			
			for(Enumeration<String> myenum = contig2.snps.keys(); myenum.hasMoreElements();){
				String mutantLine = myenum.nextElement();
				if(!this.snps.containsKey(mutantLine)){
					this.snps.put(mutantLine, contig2.snps.get(mutantLine));
				}else{
					Hashtable<Integer, SNP> h = this.snps.get(mutantLine);
					for(Enumeration<Integer> myenum2 = contig2.snps.get(mutantLine).keys(); myenum2.hasMoreElements();){
						Integer i = myenum2.nextElement();
						if(!h.containsKey(i)){
							h.put(i, contig2.snps.get(mutantLine).get(i));
						}
					}
					this.snps.put(mutantLine, h);
				}
				
			}
			
			
			
			return true;
		}
			
		
		return false;
	}
	
	
	
	public void removeCoverage(String mutantLine){
		this.coverages.remove(mutantLine);
	}
	
	
	
	/* ********************************************* *
	 * ********       Filters           ************ * 
	 * ********************************************* */
	
	
	public void filterByWildtype(){
		HashSet<Integer> usedPositions=new HashSet<Integer>();
		HashSet<Integer> removePositions=new HashSet<Integer>();
		
		try{
			for(Enumeration<Integer> myenum2 = this.snps.get(this.wildtype).keys(); myenum2.hasMoreElements();){
				Integer position= myenum2.nextElement();
				if(usedPositions.contains(position)){
					removePositions.add(position);
				}
				usedPositions.add(position);
			}
		}catch(NullPointerException e){}	
			
		
		for(Enumeration<String> myenum = this.snps.keys(); myenum.hasMoreElements();){
			String mutantLine = myenum.nextElement();
			if(mutantLine.equalsIgnoreCase(this.wildtype)){
				continue;
			}	
			
			for(Iterator<Integer> iterator = removePositions.iterator(); iterator.hasNext();){
				Integer position = iterator.next();
				this.snps.get(mutantLine).remove(position);
			}
			
			
			
		}
	
	
		
		
		
	}
	
	
	
	/**
	 * Remove every SNP from mutantLines that is also present in wildtype. 
	 * The allowed difference between frequencies is an experimental threshold to cope with collapsed assemblies.
	 * 
	 * 
	 * @param allowedDifferenceBetweenFrequencies
	 * 				only SNPs will be removed there the difference between the reference-allele frequencies is smaller than this value.  
	 */
	public void filterByWildtype(double allowedDifferenceBetweenFrequencies){
		for(Enumeration<String> myenum = this.snps.keys(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			if(key.equalsIgnoreCase(this.wildtype)){
				continue;
			}
			Hashtable<Integer, SNP> h = this.snps.get(key);
			HashSet<Integer> remove = new HashSet<Integer>();
			for(Enumeration<Integer> myenum2 = h.keys(); myenum2.hasMoreElements();){
				int position = myenum2.nextElement();
				if(snps.get(this.wildtype) != null && snps.get(this.wildtype).containsKey(new Integer(position))){
					if( Math.abs(h.get(position).getReferenceAlleleFrequency() - snps.get(this.wildtype).get(position).getReferenceAlleleFrequency() ) <= allowedDifferenceBetweenFrequencies){
						remove.add(position);
					}
				}
			}
			for(Iterator<Integer> iterator = remove.iterator(); iterator.hasNext();){
				snps.get(key).remove(iterator.next());
			}
		}
	}
	
	
	
	
	
	public void filterByPosition(Vector<String> mutantLines, int maxAllowedOccurences){
		Hashtable<Integer,Integer> usedPositions=new Hashtable<Integer,Integer>();
		HashSet<Integer> removePositions=new HashSet<Integer>();
		
		for(Enumeration<String> myenum = mutantLines.elements(); myenum.hasMoreElements();){
			String mutantLine = myenum.nextElement();
			if( this.snps.containsKey(mutantLine)){
				for(Enumeration<Integer> myenum2 = this.snps.get(mutantLine).keys(); myenum2.hasMoreElements();){
					Integer position= myenum2.nextElement();
					int num = 0;
					if(usedPositions.containsKey(position)){
						num = usedPositions.get(position);
					}
					num++;
					usedPositions.put(position, num);
					if(num > maxAllowedOccurences){
						removePositions.add(position);
						//System.out.println(this.contigName + ": removed: "+position + "\t" + num);
					}
				}
			}
			
		}
		for(Enumeration<String> myenum = mutantLines.elements(); myenum.hasMoreElements();){
			String mutantLine = myenum.nextElement();
			if(this.snps.containsKey(mutantLine)){
				for(Iterator<Integer> iterator = removePositions.iterator(); iterator.hasNext();){
					Integer position = iterator.next();
					this.snps.get(mutantLine).remove(position);
				}
			}
			
			
		}
	}
	
	
	
	
	/**
	 * This will interrogate the mutant lines given as argument and remove positions from all mutant lines
	 * if values are within the given range from the mean of mutant lines.
	 * 
	 * 
	 * @param lines
	 * 		The lines that are interrogated
	 * @param allowedDifferenceToMeanFrequency
	 * 		The difference of a reference allele frequency to the mean of interrogated ref allele frequencies.
	 * 		This is the threshold to keep a SNP.
	 */
	public void deleteDuplicateMutantPositions(String[] lines, double allowedDifferenceToMeanFrequency){
		HashSet<Integer> positions = new HashSet<Integer>();
		HashSet<Integer> dupePositions = new HashSet<Integer> ();
		
		for(int i =0; i< lines.length; i++){
			String line = lines[i];
			if( !this.snps.containsKey(line)){
				continue;
			}
			for(Enumeration<Integer> myenum = this.snps.get(line).keys(); myenum.hasMoreElements();){
				int position = myenum.nextElement();
				if( positions.contains(position)){
					dupePositions.add(position);
				}
				positions.add(position);
			}
		}
		
		HashSet<Integer> removeSNPs = new HashSet<Integer>();
		for(Iterator<Integer> iterator = dupePositions.iterator(); iterator.hasNext();){
			Integer position = iterator.next();
			Vector<Double> frequencies = new Vector<Double>();
			for( int i = 0; i< lines.length; i++){
				String line = lines[i];
				if(!snps.containsKey(line)){
					continue;
				}
				
				if(snps.get(line).containsKey(position)){
					frequencies.add(snps.get(line).get(position).getReferenceAlleleFrequency());
				}
			}
			double mean = 0;
			for(Enumeration<Double> myenum = frequencies.elements(); myenum.hasMoreElements();){
				mean = mean + myenum.nextElement().doubleValue();
			}
			boolean add = true;
			for(Enumeration<Double> myenum = frequencies.elements(); myenum.hasMoreElements();){
				double d = myenum.nextElement();
				if( Math.abs(mean/frequencies.size() - d) > allowedDifferenceToMeanFrequency){
					add = false;
				}
			}
			if(add){
				removeSNPs.add(position);
			}
			
		}
		
		
		for(Enumeration<String> myenum = this.snps.keys(); myenum.hasMoreElements();){
			String mutantLine = myenum.nextElement();
			
			for(Iterator<Integer> iterator = removeSNPs.iterator(); iterator.hasNext();){
				this.snps.get(mutantLine).remove(iterator.next());
			}
		}
		
		
	}
	
	
	
	
	
	
	
	/**
	 * 
	 * delete every SNP from every mutant line (excluding the wildtype) where the coverage is below minCoverage
	 * 
	 * 
	 * @param minCoverage
	 */
	public void deleteSNPsByCoverage(int minCoverage){
		for(Enumeration<String> myenum = snps.keys(); myenum.hasMoreElements();){
			String mutantLine = myenum.nextElement();
			if( mutantLine.equalsIgnoreCase(this.wildtype)){
				continue;
			}
			Hashtable<Integer, SNP> h = new Hashtable<Integer, SNP>();
			for(Enumeration<Integer> myenum2 = snps.get(mutantLine).keys(); myenum2.hasMoreElements();){
				Integer key = myenum2.nextElement();
				SNP snp = snps.get(mutantLine).get(key);
				if( snp.getCoverage() >= minCoverage ){
					h.put(key, snp);
				}
			}
			if(h.size() >0){
				snps.put(mutantLine, h);
			}else{
				snps.remove(mutantLine);
			}
			
		}
	}
	
	
	public void deleteSNPsByFrequency(double maxReferenceAlleleFrequency){
		for(Enumeration<String> myenum = snps.keys(); myenum.hasMoreElements();){
			String mutantLine = myenum.nextElement();
			if( mutantLine.equalsIgnoreCase(this.wildtype)){
				continue;
			}
			Hashtable<Integer, SNP> h = new Hashtable<Integer, SNP>();
			for(Enumeration<Integer> myenum2 = snps.get(mutantLine).keys(); myenum2.hasMoreElements();){
				Integer key = myenum2.nextElement();
				SNP snp = snps.get(mutantLine).get(key);
				if( snp.getReferenceAlleleFrequency() <= maxReferenceAlleleFrequency ){
					h.put(key, snp);
				}
			}
			if(h.size() >0){
				snps.put(mutantLine, h);
			}else{
				snps.remove(mutantLine);
			}
			
		}
	}
	
	
	
	
	public String getContigName(){
		return this.contigName;
	}
	

	
	
	
	
	
	
	

	/* ********************************************* *
	 * ********       Analysis          ************ * 
	 * ********************************************* */
	
	
	
	
	public boolean isCandidateSnpOnly(Vector<String> mymutants, int minSNPCoverage, double maxReferenceAlleleFrequency, int minNumberOfMutants){
		this.report = this.getContigName() +"\tlength:"+this.getLength()+"\n";
		
		
		String linereports = "";
		
		
		
		//check snp mutations
		HashSet<String> mutantLinesWithSNPs = new HashSet<String>();
		
		
		for(Enumeration<String> myenum1 =mymutants.elements(); myenum1.hasMoreElements();){
			String mutantLine = myenum1.nextElement();
			linereports = linereports + mutantLine + ": ";
			if( mutantLine.equalsIgnoreCase(wildtype)){
				continue;
			}
			
		
			if( this.snps.containsKey(mutantLine)){
				String snpreport = new String("");
				for(Enumeration<Integer> myenum2 = this.snps.get(mutantLine).keys(); myenum2.hasMoreElements();){
					Integer key = myenum2.nextElement();
					if(this.snps.get(mutantLine).get(key).getReferenceAlleleFrequency() <= maxReferenceAlleleFrequency && this.snps.get(mutantLine).get(key).getCoverage()>= minSNPCoverage){
			
						mutantLinesWithSNPs.add(mutantLine);
				
						snpreport = snpreport + ";SNP(pos:"+key.intValue()+ ",refallelefreq:" + snps.get(mutantLine).get(key).getReferenceAlleleFrequency() +
								    ",cov:" + snps.get(mutantLine).get(key).getCoverage()+
									"," + snps.get(mutantLine).get(key).getReferenceAllele() + "->" + snps.get(mutantLine).get(key).getLargestAlternativeAllele()+")";
			
			
					}
				}
				if( snpreport.length()>0){
					linereports = linereports + "\t" + snpreport.substring(1) ;
				}	
		
			}
			linereports = linereports +"\n";
		
		}
		
		this.report = this.report + "Number Of SNP mutants: " + mutantLinesWithSNPs.size() + "\n";
		
		
		
		this.report = this.report + linereports+"\n";
		
		//final result
		if( mutantLinesWithSNPs.size() >= minNumberOfMutants){
			return true;
		}
		return false;
	}
	
	
	
		public boolean isCandidate(Vector<String> mymutants, 
			int minWildtypeCoverage, 
			double maxReferenceAlleleFrequency, 
			int minCoverageToConsiderSNP,
			int minNumberOfZeroCoveragePositions,
			int minNumberOfTotalMutants, 
			boolean filterSynonymous){

			//initialize report
			this.report = this.getContigName() + "\tlength:"+this.coverages.get(this.wildtype).length + "\n";
			
			
			//check wildtype median coverage
			int wildTypeCoverage = this.getMedianCoverage(this.wildtype);
			if(wildTypeCoverage< minWildtypeCoverage ){
				return false;
			}
			
			
			String linereports = "";
			
			
			
			//check snp mutations
			HashSet<String> mutantLinesWithSNPs = new HashSet<String>();
			HashSet<String> mutantLinesWithDeletions = new HashSet<String> ();
			
			
			for(Enumeration<String> myenum1 =mymutants.elements(); myenum1.hasMoreElements();){
				String mutantLine = myenum1.nextElement();
				linereports = linereports + mutantLine + ": ";
				if( mutantLine.equalsIgnoreCase(wildtype)){
					continue;
				}
			
				int zeroCovPos = this.countZeroCoveragePositions(mutantLine, minWildtypeCoverage);
				if( zeroCovPos >= minNumberOfZeroCoveragePositions){
					linereports = linereports + "\tdeletion mutant(" + zeroCovPos + ")";
					mutantLinesWithDeletions.add(mutantLine);
				}
			
			
				if( this.snps.containsKey(mutantLine)){
					String snpreport = new String("");
					for(Enumeration<Integer> myenum2 = this.snps.get(mutantLine).keys(); myenum2.hasMoreElements();){
						Integer key = myenum2.nextElement();
						if(this.snps.get(mutantLine).get(key).getReferenceAlleleFrequency() <= maxReferenceAlleleFrequency && this.snps.get(mutantLine).get(key).getCoverage()>= minCoverageToConsiderSNP){
				
							mutantLinesWithSNPs.add(mutantLine);
					
							snpreport = snpreport + ";SNP("+key.intValue()+ "," + snps.get(mutantLine).get(key).getReferenceAlleleFrequency() +
										"," + snps.get(mutantLine).get(key).getReferenceAllele() + "->" + snps.get(mutantLine).get(key).getLargestAlternativeAllele()+")";
				
				
						}
					}
					if( snpreport.length()>0){
						linereports = linereports + "\t" + snpreport.substring(1) ;
					}	
			
				}
				linereports = linereports +"\n";
			
			}
			
			this.report = this.report + "Number Of SNP mutants: " + mutantLinesWithSNPs.size() + "\n";
			this.report = this.report + "Number Of Deletion mutants: " + mutantLinesWithDeletions.size() + "\n";
			this.report = this.report + "Wildtype coverage: "+ wildTypeCoverage + "\n";
			
			this.report = this.report + linereports+"\n";
			
			//final result
			if( mutantLinesWithSNPs.size() >= minNumberOfTotalMutants){
			
			return true;
			}
			
			return false;
			
			}
				
				
				
				
		public int countZeroCoveragePositions(String mutantLine, int minWtCoverage){
			
			
			
			int numZero = 0;
			for( int i = 0; i<this.coverages.get(mutantLine).length; i++){
				try{
					if( this.coverages.get(mutantLine)[i] >0 && this.coverages.get(wildtype)!= null && this.coverages.get(wildtype)[i]>minWtCoverage){
						if(this.coverages.get(mutantLine)==null ||this.coverages.get(mutantLine)[i]==0 ){
							numZero++;
						}
					}
				}catch(ArrayIndexOutOfBoundsException e){}
			
			}
			
			
			
			return numZero;
				
		}	
	
	
	
	
	
	
	
	
	
	
	
	/* ********************************************* *
	 * ********        Export           ************ * 
	 * ********************************************* */
	
	

	
	
	public Element getXMLElement(Document dom){
		Element contigElement = dom.createElement("TargetContig");
		
		
		
		
		contigElement.setAttribute("name",this.getContigName());
		contigElement.setAttribute("length",this.getLength()+"");
		
		
		Element snpsElement = dom.createElement("SNPs");
		contigElement.appendChild(snpsElement);
		for(Enumeration<String> myenum2 = this.snps.keys(); myenum2.hasMoreElements();){
			String mutantLine = myenum2.nextElement();
			
			for(Enumeration<Integer> myenum3 = this.snps.get(mutantLine).keys(); myenum3.hasMoreElements();){
				Integer position = myenum3.nextElement();
				SNP snp = this.snps.get(mutantLine).get(position);
				Element snpElement = snp.getXMLElement(dom);
				snpsElement.appendChild(snpElement);
			}
		}
		
		
		
		
		
		Element coveragesElement = dom.createElement("Coverages");
		contigElement.appendChild(coveragesElement);
		for( Enumeration<String> myenum2 = this.coverages.keys(); myenum2.hasMoreElements();){
			String mutantLine = myenum2.nextElement();
			int[] cov = this.coverages.get(mutantLine);
			String s = "";
			for( int i= 0; i< cov.length; i++){
				s = s +"," +cov[i];
			}
			Element coverageElement = dom.createElement("Coverage");
			coverageElement.setAttribute("name", mutantLine);
			coverageElement.setAttribute("value", s.substring(1));
			coveragesElement.appendChild(coverageElement);
		}
		
		
		
		return contigElement;
	}
	
	
	
	
	
	
	
}
