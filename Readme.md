# MutChromSeq

## About
**MutChromSeq** is a method to clone genes in plants. Prerequisites are that the gene must be related to a specific phenotype and that chromosome flow sorting is established in the species.
In a nutshell, you do an EMS mutagenesis screen of your wildtype and look for mutants that have lost the phenotype related to your gene. You need at least 5 or 6 independent mutants. If you have been working on that gene for a while, you most problably know on which chromosome it is. Otherwise find out. Flow sort chromosomes of your wildtype and mutants. Sequence amplified DNA from your chromosomes on Illumina. You need about a HiSeq2000 lane per chromosome. A denovo assembly of the wildtype data will be used as reference. You map the mutants against that reference and will find one or very few genes that have a mutation in every of your mutants. Here you find the pipeline and scripts for the analysis of the data.

## Prerequisites

### JRE 1.6
Make sure you have the Java Runtime Environments 1.6 or higher. Download from [http://java.com](http://java.com)

### BWA
Download from [http://bio-bwa.sourceforge.net/](http://bio-bwa.sourceforge.net/)

### Samtools
Download from [http://samtools.sourceforge.net/](http://samtools.sourceforge.net/)

### Blast+
Download Blast+ from [NCBI](http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

### A *denovo* assembly software.
Use your favourite assembler to assemble the rawdata of your wildtype. Note that the wrong toold might result in contigs representing a consensus of two or more NLRs. For Illumina data we had nice results with [CLC assembly cell](http://www.clcbio.com/products/clc-assembly-cell/) and default parameters. Please note that this one is a comercial software. Feel free to experiment with free software.

### RepeatMasker
Download from [http://www.repeatmasker.org/](http://www.repeatmasker.org/). If you are working on Triticeae, we recommend to use the external library [TREP](http://wheat.pw.usda.gov/ITMI/Repeats/) for masking

### MutChromSeq
Download the jar file from the release of this repository



## Preprocessing

### Clean raw data
Use [sickle](http://bioinformatics.ucdavis.edu/software/) or whatever tool you like to clean your rawdata.

Example

```
sickle pe -t sanger -f read1.fq -r read2.fq -o read1.clean.fq -p read2.clean.fq -s useless.fq
``` 

### De novo assembly of wildtype

Use your favourite assembler to assemble the rawdata of your wildtype. Note that the wrong toold might result in contigs representing a consensus of two or more NLRs. For Illumina data we had nice results with [CLC assembly cell](http://www.clcbio.com/products/clc-assembly-cell/) and default parameters.

### Check for contamination
Make sure your assembly contains contigs from your species only. This is already  pretty big data and you do not want other contamination in there. This step is probably not really necessary but it makes your life easier. A simple blast you can get rid of most of the contamination. 

### Repeatmasking
Mask your assembly. We are not interested in mutations in repeat elements. This will just create noise. But make sure you do not overmask. Some repeat libraries are "contaminated" with real genes.


### Mappings

Map (cleaned) rawdata of your wild type and the mutants against the denovo assembly. We reccomend using [bwa](http://bio-bwa.sourceforge.net/) and [samtools](http://samtools.sourceforge.net/).

**The filename of the pileup file will be used internally as an identifier. Please make sure this filename is unique for each mutant. E.g. pileup-m1.txt, pileup-m2.txt ...**

```
bwa index assembly.fasta
bwa aln assembly.fasta read1.fastq > read1.aln
bwa aln assembly.fasta read2.fastq > read2.aln
bwa sampe assembly.fasta read1.aln read2.aln read1.fastq read2.fastq > raw.sam
samtools view -f2 -Shub -o raw.bam raw.sam
samtools sort raw.bam sorted
samtools rmdup sorted.bam rmdup.bam
samtools index rmdup.bam
samtools faidx assembly.fasta
samtools mpileup -f assembly.fasta -BQ0 rmdup.bam > pileup.txt 
```

## MutChromSeq pipeline

### Parse pileups
Parse individual pileups and convert to an xml format. This has been separated from the actual MutantHunter because the pielups can be parsed in parallel. Download Pileup2XML.jar from MutantHunter release.

```
java -jar Pileup2XML.jar -i pileup.txt -o pileup.xml -f reference.fasta -a 0.01 -c 10

```

**Please not that you have to add the `-w` option for running Pileup2XML on wild type**

parameter | argument        | description
---       |   ---           | ---
**-i**    | *STR*           | The pileup file
**-o**    | *STR* 			 | The output file
**-f**    | *STR*           | The reference fasta file (the denovo assembly of your wildtype)
**-a**    | *float*         | Reference allele frequency to detect a SNP. In an ideal world, this would be 0.0. You might want to allow for sequencing errors or false mappings. Default is 0.01
**-c**    | *int*           | Minimum coverage for position to be regarded. Default is 5
**-w*     |			        | Flag to mark the pileup originating from wildtype data mapping to your assembly.


### MutChromSeq

This will provide the candidate contigs.  You will have to allocate a bit of memory to the java vm. To add e.g. 16 Gb or RAM you would write `java -Xmx16000M -jar ...`

The `-s` parameter allows to load a list of identifiers that refer to sequences in the reference assembly. Further processing is only done for those sequences. If `-s` is ommitted, the entire set of sequences is used. The argument to `-s` is a either text file containing one identifier per line or a TSV file where the first column contains identifiers.

```
java -jar MutChromSeq.jar -w wildtype.pileup.xml -m mutant1.pileup.xml mutant2.plieup.xml [...] -o output.txt -n 6 -c 10 -a 0.01 -z 2 
```
parameter | argument        | description
---       |   ---           | ---
**-w**    | *STR*           | The xml generated from the wildtype pileup using Pileup2XML.jar
**-m**    | *STR* [*STR*]+  | The xml files of each mutant generated using Pileup2XML.jar
**-o**    | *STR*           | Outputfile
**-a**    | *float*         | Maximum reference allele frequency of a SNP to be reported. Default is 0.01
**-c**    | *int*           | Minimum coverage for position to be regarded. Default is 10
**-n**    | *int*           | Minimum number of mutants to report a contig. Default is 2
**-z**    | *int*           | Number mutant lines that are allowed to have SNV in same position. Default is 2
**-s**    | *Str*           | A file with a list of sequence identifiers



## Tips

### Manual validation
Load your bam file in a genome browser, e.g. [Savant](http://www.genomesavant.com/p/home/index) and check your candidate contigs


## Contact
If there are any issues with the tool or if you would like to collaborate with us, please don't hesitate to contact [us](mailto:burkhard.steuernagel@jic.ac.uk).