package mutantHunter;

import org.w3c.dom.Document;
import org.w3c.dom.Element;

import support.MPileupLine;

public class SNP {

	String mutantLine;
	String contig;
	char referenceAllele;
	int position;
	int coverage;
	int alleleA;
	int alleleC;
	int alleleG;
	int alleleT;
	
	
	
	
	
	public SNP(String mutantLine, String contig, char referenceAllele,
			int position, int coverage, int alleleA, int alleleC, int alleleG,
			int alleleT) {
		super();
		this.mutantLine = mutantLine;
		this.contig = contig;
		this.referenceAllele = referenceAllele;
		this.position = position;
		this.coverage = coverage;
		this.alleleA = alleleA;
		this.alleleC = alleleC;
		this.alleleG = alleleG;
		this.alleleT = alleleT;
	}

	public SNP(String mutantLine, String mpileupLine){
		this.mutantLine = mutantLine;
		
		MPileupLine line = new MPileupLine(mpileupLine);
		this.position = line.getPosition();
		this.contig = line.getChromosome();
		this.referenceAllele = line.getRefBase();
		this.coverage = line.getCoverage();
		this.alleleA = 0;
		this.alleleT = 0;
		this.alleleG = 0;
		this.alleleC = 0;
		
		if(line.getAlleles().containsKey("A")){
			this.alleleA = line.getAlleles().get("A");
		}
		if(line.getAlleles().containsKey("T")){
			this.alleleT = line.getAlleles().get("T");
		}
		if(line.getAlleles().containsKey("G")){
			this.alleleG = line.getAlleles().get("G");
		}
		if(line.getAlleles().containsKey("C")){
			this.alleleC = line.getAlleles().get("C");
		}
	}
	
	public SNP(String mutantLine, MPileupLine line){
		this.position = line.getPosition();
		this.mutantLine = mutantLine;
		this.contig = line.getChromosome();
		this.referenceAllele = line.getRefBase();
		this.coverage = line.getCoverage();
		this.alleleA = 0;
		this.alleleT = 0;
		this.alleleG = 0;
		this.alleleC = 0;
		
		if(line.getAlleles().containsKey("A")){
			this.alleleA = line.getAlleles().get("A");
		}
		if(line.getAlleles().containsKey("T")){
			this.alleleT = line.getAlleles().get("T");
		}
		if(line.getAlleles().containsKey("G")){
			this.alleleG = line.getAlleles().get("G");
		}
		if(line.getAlleles().containsKey("C")){
			this.alleleC = line.getAlleles().get("C");
		}
		
	}
	
	
	
	public SNP(Element xmlElement){
		this.mutantLine = xmlElement.getAttribute("mutantLine");
		this.contig   = xmlElement.getAttribute("contig");
		this.referenceAllele = xmlElement.getAttribute("referenceAllele").toCharArray()[0];
		this.position = Integer.parseInt(xmlElement.getAttribute("position"));
		this.coverage = Integer.parseInt(xmlElement.getAttribute("coverage"));
		this.alleleA  = Integer.parseInt(xmlElement.getAttribute("alleleA"));
		this.alleleC  = Integer.parseInt(xmlElement.getAttribute("alleleC"));
		this.alleleG  = Integer.parseInt(xmlElement.getAttribute("alleleG"));
		this.alleleT  = Integer.parseInt(xmlElement.getAttribute("alleleT"));
	}
	
	
	
	
	public int getPosition(){
		return this.position;
	}
	
	
	public double getReferenceAlleleFrequency(){
		int cov = alleleA + alleleG + alleleT + alleleC;
		if(referenceAllele=='A'){return (double)alleleA / cov;}
		if(referenceAllele=='T'){return (double)alleleT / cov;}
		if(referenceAllele=='G'){return (double)alleleG / cov;}
		if(referenceAllele=='C'){return (double)alleleC / cov;}
		return 1;
	}
	
	
	
	public Element getXMLElement(Document dom){
		Element snpElement = dom.createElement("SNP");
		
		snpElement.setAttribute("mutantLine", this.mutantLine);
		snpElement.setAttribute("contig", this.contig);
		snpElement.setAttribute("referenceAllele", this.referenceAllele+"");
		snpElement.setAttribute("position", this.position+"");
		snpElement.setAttribute("coverage", this.coverage+"");
		snpElement.setAttribute("alleleA", this.alleleA+"");
		snpElement.setAttribute("alleleT", this.alleleT+"");
		snpElement.setAttribute("alleleG", this.alleleG+"");
		snpElement.setAttribute("alleleC", this.alleleC+"");
		
		
		return snpElement;
	}
	
	public int getCoverage(){
		return this.coverage;
	}
	
	
	
	public char getReferenceAllele(){
		return referenceAllele;
	}
	
	public String getContig(){
		return this.contig;
	}
	public String getMutantLine(){
		return mutantLine;
	}
	
	public char getLargestAlternativeAllele(){
		char myAllele = 'A';
		int alleleCount=alleleA;
		if(this.referenceAllele=='A'){
			myAllele = 'C';
			alleleCount = alleleC;
		}
		if( alleleCount < alleleC && alleleC != getReferenceAllele()  ) {
			myAllele='C';
			alleleCount = alleleC;
		}
		
		
		if( alleleCount < alleleG && alleleG != getReferenceAllele()  ) {
			myAllele='G';
			alleleCount = alleleG;
		}
		
		if( alleleCount < alleleT && alleleT != getReferenceAllele()  ) {
			myAllele='T';
			alleleCount = alleleT;
		}
		
		
		return myAllele;
	}
	
	
	
	public void setMutantLine(String mutantLine){
		this.mutantLine = mutantLine;
	}
	
	
}
