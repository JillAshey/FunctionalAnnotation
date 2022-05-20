# Functional annotation workflow

## Table of Contents


This project aims to develop a functional genomic annotation workflow for non-model organisms (ie corals). The goal of functional annotation is to identify and tag genes in a refernce genome with known functions of homologous genes in other organisms. The following document is intended as a tutorial in understanding functional annotation in non-model organisms. The genomic information from the coral *Acropora cervicornis* is used in this workflow.

For this functional annotation workflow tutorial, you will need:

- Access to a high performance computing server with the following programs:
	- BLAST+/2.11.0-gompi-2020
	- DIAMOND v2.0.0-GCC-8.3.0
	- InterProScan/5.44-79.0-foss-2018b  
	- Java/11.0.2
- Laptop with access to the Internet and the following programs installed: 
	- Blast2GO Basic 
	- R v4.0.2
	- RStudio v1.3.959

All analyses done on Putnam Lab Node

### Step 1: Obtain sequences of interest. 

In order to conduct functional annotation steps, protein and transcript sequences are needed. There are two main hubs where coral genomic information is stored: [Reef Genomics](http://reefgenomics.org) and [NCBI](https://www.ncbi.nlm.nih.gov). Other researchers store their genomic infomation on their own personal webpages. Genomic information must be downloaded in order to proceed. 

#### i) Identify species to work with. 

For this project, the coral species of interest are *Acropora cervicornis*, *Montastraea cavernosa*, *Montipora capitata*, *Pocillopora acuta*, *Porites lobata*, and *Orbicella favelota*. These species were selected because they are popular corals to work with and the [Putnam Lab](http://putnamlab.com) is currently working with them.

#### ii) Download genomic files for species of interest. 
 
##### [Acropora cervicornis](http://baumslab.org/research/data/)
###### Unpublished, obtained through personal communication with Dr. Iliana Baums at PSU.

Unfortunately, the *A. cervicornis* genomic information is not publically available yet. The Putnam lab got permission to use the version created by Drs. Iliana Baums and Shelia Kitchen. 

To download genome information for other species, go [here](https://github.com/JillAshey/FunctionalAnnotation/blob/main/wget_genomes.md).

This workflow will only use the protein sequences. 

### Step 2: Identify homologous sequences

Homology refers to the similarity of structure or genes in different taxa due to shared ancestry. {example}

Sequence homology is the homology between DNA, RNA, and protein sequences in terms of shared ancestry. Sequence homology is usually inferred by the similarity of nucleotide or amino acid sequences. Strong sequence similarity (or percent homology) provides evidence that two or more sequences are related through shared ancestry. 

The most common tool to compare sequences to various databases is [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) (Basic Local Alignment Search Tool - what a great name). BLAST compares nucleotide or protein sequences of interest (called a query) to their sequence databases to find regions of similarity. If a nucleotide or protein sequence of interest significantly matches a sequence/sequences in the databases, BLAST will tag the sequence of interest with the information about the known sequence(s). For example, if a new transcript sequence is identified in a mouse, BLAST could be used to see if any other animals carry a similar sequence, and if they do, what the biological functions of the sequence are. 

Several databases are available to assess homology, including [NCBI](https://www.ncbi.nlm.nih.gov/protein), [SwissProt](https://www.uniprot.org), [Trembl](https://www.uniprot.org), and [InterPro](https://www.ebi.ac.uk/interpro/). Note that SwissProt and Trembl are both part of the larger Uniprot database. Go [here](https://www.biostars.org/p/77257/) for an explanation about Uniprot.

Below is the code that details how to BLAST to each of these databases (will not be BLASTing with InterPro, it's its own thing).

#### i) NCBI nr database 

Using Diamond BLAST

##### a) On a HPC server, download nr database from [NCBI](ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz). Convert this database to be Diamond-readable. 

The nr (non-redundant) database is a collection of non-identical protein sequences compiled by NCBI. It is updated on a daily basis.

This script is in the sbatch_executables subdirectory in the Putnam Lab shared folder. The original script, created by Erin Chille on August 6, 2020, downloads the most recent nr database in FASTA format from NCBI and uses it to make a Diamond-formatted nr database. This step was updated by Danielle Becker-Polinski on September 24th, 2021 because the scripts were not including the full CPUs to download and a couple other formatting errors. To use this script, run ```download_nr_database.sh```
and ```make_diamond_nr_db.sh``` in that order (located in the `/data/putnamlab/shared/sbatch_executables` directory path). 

##### b) Align query protein sequences against database 

Now that the reference database has been properly generated, the sequences of interest can be aligned against it.

Before aligning, count number of protein sequences:

```
# Count number of protein sequences
zgrep -c ">" /data/putnamlab/jillashey/genome/Acerv/Acerv_assembly_v1.0.protein.fa 
33322
```

Write and run script:

```
nano acerv_diamond_blastp.sh

#!/bin/bash 
#SBATCH --job-name="diamond-blastp-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="acerv_diamond_blastp_out_error"
#SBATCH --output="acerv_diamond_blastp_out"
#SBATCH --exclusive

echo "START" $(date)
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

echo "Updating Acerv annotation" $(date)
diamond blastp -b 2 -d /data/putnamlab/shared/databases/nr.dmnd -q /data/putnamlab/jillashey/genome/Acerv/Acerv_assembly_v1.0.protein.fa -o Acerv_blastp_annot -f 100 -e 0.00001 -k 1 --threads $SLURM_CPUS_ON_NODE --tmpdir /tmp/

echo "Search complete... converting format to XML and tab"

diamond view -a Acerv_blastp_annot.daa -o Acerv_blastp_annot.xml -f 5
diamond view -a Acerv_blastp_annot.daa -o Acerv_blastp_annot.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

echo "STOP" $(date)

sbatch acerv_diamond_blastp.sh
# Submitted batch job 93763
```

Output: .tab and .xml file of aligned sequence info. The .xml file is the important one - it will be used as input for Blast2GO

After script has run successfully (which may take hours to days), check how many hits the protein sequences got: 

```
wc -l Acerv_blastp_annot.tab
30990 Acerv_blastp_annot.tab # ~31000 hits out of 33000
```

##### c) Secure-copy output files to local computer 

```
# From a local terminal window (ie not a remote server)

scp jillashey@bluewaves.uri.edu:xxxxx /path/to/local/computer/directory
```

DIAMOND BLAST results can now be used in further analyses. To see notes for running Diamond w/ the nr database on all species, go [here](https://github.com/JillAshey/FunctionalAnnotation/blob/main/scripts/Diamond-NCBI_BLAST.md).

#### ii) SwissProt database 

##### a) On a HPC server, download swissprot database.

This script is in the shared/sbatch_excutables folder on Putnam node. This script was written by Danielle Becker.

```
#!/bin/bash
#SBATCH --job-name="swiss-ref"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu # CHANGE EMAIL
#SBATCH -D /data/putnamlab/shared/databases/swiss_db

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b

cd databases/swiss_db

echo "Making swissprot database" $date
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
makeblastdb -in uniprot_sprot.fasta -parse_seqids -dbtype prot -out swissprot_20211025
echo "STOP" $(date)

sbatch download_swissprot_database.sh 
# Submitted batch job 94130
```

##### b) Align query protein sequences against database 

Now that the SwissProt reference database has been properly generated, the sequences of interest can be aligned against it.

```
nano swissprot_blast.sh

#!/bin/bash 
#SBATCH --job-name="swissprot-blastp-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="swissprot_blastp_out_error"
#SBATCH --output="swissprot_blastp_out"
#SBATCH --exclusive

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b #load blast module

echo "Blast against swissprot database" $(date)
blastp -max_target_seqs 5 -num_threads 20 -db /data/putnamlab/shared/databases/swiss_db/swissprot_20211022 -query /data/putnamlab/jillashey/genome/Acerv/Acerv_assembly_v1.0.protein.fa -evalue 1e-5 -outfmt '5 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -out acerv_swissprot_protein.out

echo "STOP" $(date)

sbatch swissprot_blast.sh 
# Submitted batch job 94678
```

Output: .xml file of aligned sequence info

##### c) Secure-copy output files to local computer 

```
# From a local terminal window (ie not a remote server)

scp jillashey@bluewaves.uri.edu:xxxxx /path/to/local/computer/directory
```

SwissProt results can now be used in further analyses. To see notes for running blastp w/ the SwissProt database on all species, go [here](https://github.com/JillAshey/FunctionalAnnotation/blob/main/scripts/SwissProt.md).


#### iii) Trembl database 

##### a) On a HPC server, download Trembl database.

This script is in the shared/sbatch_excutables folder on Putnam node. This script was written by Danielle Becker.

```
#!/bin/bash
#SBATCH --job-name="ref"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu # CHANGE EMAIL
#SBATCH -D /data/putnamlab/shared/databases/trembl_db

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b

cd databases/trembl_db

echo "Making trembl database" $date
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
gunzip uniprot_trembl.fasta.gz
makeblastdb -in uniprot_trembl.fasta -parse_seqids -dbtype prot -out trembl_20211025 # CHANGE DATE
echo "STOP" $(date)

sbatch download_trembl_database.sh
# Submitted batch job 94144
```

##### b) Align query protein sequences against database 

Now that the Trembl reference database has been properly generated, the sequences of interest can be aligned against it.

Trembl takes a long time to run. I found that splitting the protein file into multiple files and running each file on its own makes the whole process much faster. 

To split the protein file, use [PyFasta](https://github.com/brentp/pyfasta) - it will be able to split a fasta file into several new files of relatively even size. This software is only available on bluewaves so do this step on bluewaves and then switch back to andromeda. 

```
module load pyfasta/0.5.2

pyfasta split -n 6 Acerv_assembly_v1.0.protein.fa
creating new files:
Acerv_assembly_v1.0.protein.0.fa
Acerv_assembly_v1.0.protein.1.fa
Acerv_assembly_v1.0.protein.2.fa
Acerv_assembly_v1.0.protein.3.fa
Acerv_assembly_v1.0.protein.4.fa
Acerv_assembly_v1.0.protein.5.fa
```

Copy files into trembl folder

```
cp Acerv_assembly_v1.0.protein.*.fa /data/putnamlab/jillashey/annotation/trembl/acerv
```

Run trembl as an array job on ANDROMEDA 

```
nano acerv_trembl_blastp.sh

#!/bin/bash 
#SBATCH --job-name="trembl-blastp-protein"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --mem=120GB
#SBATCH --error="acerv_trembl_blastp_out_error"
#SBATCH --output="acerv_trembl_blastp_out"
#SBATCH --exclusive

echo "START" $(date)

module load BLAST+/2.11.0-gompi-2020b #load blast module

#F=/data/putnamlab/jillashey/annotation/trembl/acerv

echo "Blast against trembl database" $(date)

array1=($(ls Acerv_assembly_v1.0.protein.*.fa))
for i in ${array1[${SLURM_ARRAY_TASK_ID}]}
do
blastp -max_target_seqs 1 \
-num_threads $SLURM_CPUS_ON_NODE \
-db /data/putnamlab/shared/databases/trembl_db/trembl_20220307 \
-query ${i} \
-evalue 1e-5 \
-outfmt '5 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' \
-out ${i}.xml
done

echo "STOP" $(date)

sbatch --array=0-5 acerv_trembl_blastp.sh
#Submitted batch job 135768
```

Trembl results can now be used in further analyses. To see notes for running blastp w/ the Trembl database on all species, go [here](https://github.com/JillAshey/FunctionalAnnotation/blob/main/scripts/Trembl.md).

#### iv) InterPro

[InterProScan](https://www.ebi.ac.uk/interpro/) is a software that provides functional analysis of proteins by using predictive models that look across protein databases to find homologies with query proteins. If a homology is identified in one of the databases, the information about that homology is used to assign the query protein a GO term. Essentially, InterProScan does both steps - finds homologies and assigns GO terms!

Some protein files have * characters in the sequences and InterProScan aint too happy about that. Have to remove the * with sed
 
```
sed -i s/\*//g Acerv_assembly_v1.0.protein.fa
```

i = in-place (edit file in place)
s = substitute 
/replacement_from_reg_exp/replacement_to_text/ = search and replace statement 
\* = what I want to replace
Add nothing for replacement_to_text
g = global (replace all occurances in file)

##### a) On a HPC server, run IPS.

```
nano IPS_acerv.sh

#!/bin/bash
#SBATCH --job-name="InterProScan"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu

cd /data/putnamlab/jillashey/annotation/InterProScan/acerv

echo "START $(date)"

# Load module
# module load InterProScan/5.46-81.0-foss-2019b - version erin had in her code, not on bluewaves
module load InterProScan/5.44-79.0-foss-2018b  
module load Java/11.0.2
java -version

# Run InterProScan
interproscan.sh -version
interproscan.sh -f XML -i Acerv_assembly_v1.0.protein.fa -b acerv.interpro -iprlookup -goterms -pa 
interproscan.sh -mode convert -f GFF3 -i acerv.interpro.xml -b acerv.interpro

echo "DONE $(date)"

sbatch IPS_acerv.sh
```

Submitted batch job 1761429

##### b) Secure-copy output files to local computer 

```
# From a local terminal window (ie not a remote server)

scp jillashey@bluewaves.uri.edu:xxxxx /path/to/local/computer/directory
```

IPS results can now be used in further analyses. To see notes for running IPS on all species, go [here](https://github.com/JillAshey/FunctionalAnnotation/blob/main/InterProScan/InterProScan.md)

### Step 3: Assign gene ontology terms to sequences

After DIAMOND BLAST is completed, analysis can move to assigning gene ontology (GO) terms to sequences. 

The [Gene Ontology](http://geneontology.org) is an extensive consortium that aims to understand gene function and provide functional annotation for organisms across the tree of life. It also maintains a controlled vocabulary of gene and gene attributes across species.

The Gene Ontology has a system to classify genes into terms. The terms are grouped into 3 categories:

1. Molecular funciton - molecular-level activities performed by gene products (ie RNA or proteins). Describes the activities rather than the entities (molecules or protein complexes)
	- Examples: Toll receptor binding, transcription regulator activity
2. Cellular component - cellular anatomy and/or locations where gene products perform function
	- Examples: mitochondrion, ribosome 
3. Biological process - larger biological processes accomplished by interaction of multiple molecular activities
	- Examples: glucose transmembrane transport, DNA repair 

With this information, genes can be described/annotated with multiple terms. These terms are called GO terms. Here is the basic structure of a GO term:

```
GOID: GO:0007165Term: signal transductionOntology: BPDefinition: The cellular process in which a signal is conveyed to    trigger a change in the activity or state of a cell. Signal    transduction begins with reception of a signal (e.g. a ligand    binding to a receptor or receptor activation by a stimulus such as    light), or for signal transduction in the absence of ligand,    signal-withdrawal or the activity of a constitutively active    receptor. Signal transduction ends with regulation of a downstream    cellular process, e.g. regulation of transcription or regulation of    a metabolic process. Signal transduction covers signaling from    receptors located on the surface of the cell and signaling via    molecules located within the cell. For signaling between cells,    signal transduction is restricted to events at and within the    receiving cell.Synonym: GO:0023033Synonym: signaling pathwaySynonym: signalling pathwaySynonym: signaling cascadeSynonym: signalling cascadeSecondary: GO:0023033
```

**Elements**

- GOID - unique 7-digit identifier that is intended to be machine readable
- Term - human-readable name
- Ontology - molecular function (MF), cellular component (CC), biological process (BP)
- Definition - description of what the term means
- Synonym - alternate words closely related in meaning to term, relevant connections
- Secondary - ID created when 2 or more terms are identical in meaning and so are merged into a single term. Secondary ID preserves the excess GO terms

There is so much more information available with these terms, including relationships to other genes/gene products, graphical representation of related GO terms, and much more, but that is beyond the scope of this analysis. 

#### i) Run BLAST2GO to obtain GO terms

[BLAST2GO](https://www.blast2go.com) (B2G) is a bioinformatics tools for functional annotation and analysis of gene or protein sequences. It was originally developed to provide a user-friendly interface for GO annotation and now hosts many different functional annotation tools. The B2G software makes it possible to generate annotation without requiring writing any code. While B2G can do an extensive array of functions, this analysis primarily utilizes the GO mapping and annotation functions. 

##### a) Download BLAST2GO to personal computer and activate the Basic subscription plan. 

The B2G application can be downloaded [here](https://www.blast2go.com/blast2go-pro/download-b2g). B2G is available for Mac, Windows, and Linux systems. 2GB of RAM is recommended. Additionally, Internet connection is required to use most application features.

Register for B2G Basic [here](https://www.blast2go.com/b2g-register-basic). B2G Basic is free and includes the necessary features for this analysis. Registering will generate an activation key, which be put into the B2G software. You must be a part of some research institution to obtain B2G Basic. 

##### b) Load the XML files generated from NCBI, SwissProt, and Trembl databases . 

To load the file, go to File<Load<Load Blast results<Load Blast XML (Legacy)

Once the file is loaded, a table loads with info about those BLAST results (nr, Tags, SeqName, Description, Length, Hits, e-Value, and sim mean). All of the cells should be orange with Tags that say BLASTED. This indicates that these sequences have only been blasted, not mapped or annotated. 

Only one .xml file can be loaded at here! To analyze all .xml files, go through these steps with each one separately.

##### c) Map GO terms

Mapping is the process of retrieving GO terms associated with the Description obtained by the DIAMOND BLAST search. Several mapping {steps} occur:

- BLAST result accessions are used to retrieve gene names from NCBI and those gene names are used to search the GO database. 
- BLAST result accessions are used to retrieve protein information (with GO terms already annotated) through UniProt, which uses the databases SD, UniProt, Swiss-Prot, TrEMBL, RefSeq, GenPept and PDB.
- BLAST result accessions are directly searched in GO database. 

To map results, select the mapping icon (white circle with green circle inside) at the top of the screen. Then select Run Mapping. A box will open up; don't change anything, click run. Depending on the number of BLAST results, mapping could take hours to days. B2G will continue running if the computer is in sleep mode. Mapping status can be checked under the Progress tab in the lower left box. If mapping retrieved a GO term that may be related to a certain sequence, that sequence row will turn light green. Only move forward when the Progress tab is 100% completed.

##### d) Annotate GO terms

Now that GO terms have been retrieved, the annotation process will select GO terms from the GO pool obtained through mapping and assign them to query sequences. Several annotation steps occur:

- For all found GO terms, an annotation rule (AR) is applied. The rule seeks to find the most specific annotations with a certain level of reliability and can be adjusted in the settings prior to running annotation. 
- If a candidate GO term is found, an annotation score is calculated that weights the highest hit similarity of that candidate GO term by its evidence code (from mapping step). The candidate GO term will only be assigned to a sequence if it passes a certain threshold based on calculations above. For more info about the math behind annotation, go [here](http://docs.blast2go.com/user-manual/gene-ontology-annotation/). 

To annotate, select the annot icon (white circle with blue circle inside) at the top of the screen. Then select Run Annotation. A box will open up; don't change anything unless thresholds need to be adjusted. If no changes are necessary, click next through the boxes until the final one and click run. Depending on the mapping results, annotating could take hours to days. B2G will continue running if the computer is in sleep mode. Annotation status can be checked under the Progress tab in the lower left box. If a GO term has been successfully assigned to a sequence, that sequence row will turn blue. Only move forward when the Progress tab is 100% completed.

##### e) Merge with InterProScan results 

A .xml file was generated during the IPS step. This file can be loaded into B2G and merged with the newly mapped and annotated data from NCBI, SwissProt, or Trembl. Once these files are merged, new columns will appear with IPS GO and enzyme code information. 

##### f) Export annotation 

Once the mapping, annotating, and IPS merging are complete, download the table as a .txt file (this is the only option to save in B2G, can change to .csv later if needed). (still need to add specific steps)

To see notes for running B2G with all species, go [here](https://github.com/JillAshey/FunctionalAnnotation/blob/main/Blast2GO/Blast2GO.md). 

### Step 4: Merge all information for full annotation

This should probably be done in R, but I did it the janky way in excel. 

First, open B2G annotation files (should have one for NCBI, SwissProt, and Trembl annotations) check to make sure they all had same number of sequence IDs and that the sequence IDS were in the same order. Once that is checked, open a blank excel sheet and copy/paste the annotation information from all the files into the new file. Make new column names so that it is clear that the info in a specific column is from NCBI, SwissProt, Trembl, or IPS.

Once all annotation info is together in a new sheet, you got your annotation table! There will likely be blank spaces or NAs - that is okay, it just means that that sequence didn't get a hit in one or more of the databases. To evaluate and compare the GO term information from each database, follow the code [here](https://github.com/JillAshey/FunctionalAnnotation/blob/main/RAnalysis/annotation_20211102.Rmd) (still a draft). 


