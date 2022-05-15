## Diamond BLAST 

### Comparison between blastx and blastp 

| Species   | Site        | TotalGenes | BlastxHits | %HitTotal.Gene | TotalProtein | BlastpHits | %HitTotal.Protein |
| --------- | ----------- | ---------- | ---------- | -------------- | ------------ | ---------- | ----------------- |
| Acerv     | Florida     | 33322      | 29515      | 88.5751155     | 33322        | 30990      | 93.0016206        |
| Mcav      | Florida     | 25142      | 3191       | 12.69191       | 25142        | 23521      | 93.5526211        |
| Ofav      | Florida     | 35971      | 10122      | 28.1393345     | 32587        | 32536      | 99.8434959        |
| Mcap      | Hawaii      | 63227      | 55217      | 87.3313616     | 63227        | 55683      | 88.0683885        |
| Pcomp     | Hawaii      | 74728      | 59049      | 79.018574      | 74728        | 58796      | 78.6800128        |
| Pdam.NCBI | Hawaii      | 27287      | 3406       | 12.4821343     | 25183        | 25170      | 99.9483779        |
| Pdam.RG   | Hawaii      | 26077      | 3732       | 14.3114622     | 26077        | 25947      | 99.5014764        |
| Plob      | Hawaii      | 31126      | 4315       | 13.8630084     | 31126        | 28284      | 90.8693697        |
| Apoc      | RhodeIsland | 48184      | 40031      | 83.0794455     | 45867        | 39699      | 86.5524233        |

The original script, created by Erin Chille on August 6, 2020, downloads the most recent nr database in FASTA format from NCBI and uses it to make a Diamond-formatted nr database. This step was updated by Danielle Becker-Polinski on September 24th, 2021 because the scripts were not including the full CPUs to download and a couple other formatting errors. Go to the sbatch_executables subdirectory in the Putnam Lab shared folder and run the scripts, make_diamond_nr_db.sh and make_diamond_nr_db.sh in this order:

```
$ sbatch download_nr_database.sh
Submitted batch job NNN
$ sbatch -d afterok:NNN make_diamond_nr_db.sh
```

Previously, I ran blastp with protein queries against protein db, but here I'm going to try blastx with transcript queries against protein db. 

*All analyses done on Putnam Lab Node*

### Florida species

#### Acerv 

Ran previously, failed before I got any results due to exceeded memory. Rerun on Putnam lab node

##### BLASTx
```
# Check number of genes 
zgrep -c "^>" Acerv_assembly_v1.0.mRNA.fa
33322 

nano acerv_diamond_blastx.sh

#!/bin/bash
#SBATCH --job-name="diamond-blastx"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH -p putnamlab
#SBATCH --mem=100GB
#SBATCH --error="acerv_diamond_blastx_out_error"
#SBATCH --output="acerv_diamond_blastx_out"

echo "START" $(date)
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

echo "Updating Acerv annotation" $(date)
diamond blastx -d /data/putnamlab/shared/databases/nr.dmnd -q Acerv_assembly_v1.0.mRNA.fa -o Acerv_annot -f 100 -b20 --more-sensitive -e 0.00001 -k1

echo "Search complete... converting format to XML and tab"

diamond view -a Acerv_annot.daa -o Acerv_annot.xml -f 5
diamond view -a Acerv_annot.daa -o Acerv_annot.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

echo "STOP" $(date)

sbatch acerv_diamond_blastx.sh

# After script has run:
wc -l Acerv_annot.tab
29515 Acerv_annot.tab # got ~29500 hits out of ~33000 genes

```
Submitted batch job 13114 - finished

##### BLASTp

**20211021 update** - rerunning diamond blast w/ Acerv proteins. Danielle got more hits on Diamond BLAST when she used blastp and used protein sequences instead of using blastx and the CDS sequences.

```
# Count number of protein sequences
zgrep -c ">" /data/putnamlab/jillashey/genome/Acerv/Acerv_assembly_v1.0.protein.fa 
33322

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

# After script has run: 
wc -l Acerv_blastp_annot.tab
30990 Acerv_blastp_annot.tab # ~31000 hits out of 33000
```
Submitted batch job 93763

Going to try with the SwissProt database 

First, I need to make the swissprot db diamond-readable 

```
cd /data/putnamlab/shared/sbatch_executables

```


#### Mcav 

##### BLASTx

```
# Check number of genes 
zgrep -c "^>" Mcavernosa.maker.transcripts.fasta
25142 

nano mcav_diamond_blastx.sh

#!/bin/bash
#SBATCH --job-name="diamond-blastx"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH -p putnamlab
#SBATCH --mem=100GB
#SBATCH --error="mcav_diamond_blastx_out_error"
#SBATCH --output="mcav_diamond_blastx_out"

echo "START" $(date)
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

echo "Updating Mcav annotation" $(date)
diamond blastx -d /data/putnamlab/shared/databases/nr.dmnd -q Mcavernosa.maker.transcripts.fasta -o Mcav_annot -f 100 -b20 --more-sensitive -e 0.00001 -k1

echo "Search complete... converting format to XML and tab"

diamond view -a Mcav_annot.daa -o Mcav_annot.xml -f 5
diamond view -a Mcav_annot.daa -o Mcav_annot.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

echo "STOP" $(date)

sbatch mcav_diamond_blastx.sh

# After script has run:
wc -l Mcav_annot.tab
3191 Mcav_annot.tab # only got ~3000 hits out of ~25000 genes

```

Submitted batch job 13458 - finished

##### BLASTp
**20211021 update** - rerunning diamond blast w/ Mcav proteins. Danielle got more hits on Diamond BLAST when she used blastp and used protein sequences instead of using blastx and the CDS sequences.

```
# Count number of protein sequences
zgrep -c ">" /data/putnamlab/jillashey/genome/Mcav/Mcavernosa_annotation/Mcavernosa.maker.proteins.fasta
25142

nano mcav_diamond_blastp.sh

#!/bin/bash 
#SBATCH --job-name="diamond-blastp-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="mcav_diamond_blastp_out_error"
#SBATCH --output="mcav_diamond_blastp_out"
#SBATCH --exclusive

echo "START" $(date)
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

echo "Updating Mcav annotation" $(date)
diamond blastp -b 2 -d /data/putnamlab/shared/databases/nr.dmnd -q /data/putnamlab/jillashey/genome/Mcav/Mcavernosa_annotation/Mcavernosa.maker.proteins.fasta -o Mcav_blastp_annot -f 100 -e 0.00001 -k 1 --threads $SLURM_CPUS_ON_NODE --tmpdir /tmp/

echo "Search complete... converting format to XML and tab"

diamond view -a Mcav_blastp_annot.daa -o Mcav_blastp_annot.xml -f 5
diamond view -a Mcav_blastp_annot.daa -o Mcav_blastp_annot.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

echo "STOP" $(date)

# After script has run: 
wc -l Mcav_blastp_annot.tab
23521 Mcav_blastp_annot.tab # ~23500 hits out of ~25000 protein seqs, considerably more than blastx
```
Submitted batch job 93764

#### Ofav

##### BLASTx

```
# Check number of genes 
zgrep -c "^>" GCF_002042975.1_ofav_dov_v1_rna.fna.gz
35971 

nano ofav_diamond_blastx.sh

#!/bin/bash
#SBATCH --job-name="diamond-blastx"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH -p putnamlab
#SBATCH --mem=100GB
#SBATCH --error="ofav_diamond_blastx_out_error"
#SBATCH --output="ofav_diamond_blastx_out"

echo "START" $(date)
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

echo "Updating Ofav annotation" $(date)
diamond blastx -d /data/putnamlab/shared/databases/nr.dmnd -q GCF_002042975.1_ofav_dov_v1_rna.fna -o Ofav_annot -f 100 -b20 --more-sensitive -e 0.00001 -k1

echo "Search complete... converting format to XML and tab"

diamond view -a Ofav_annot.daa -o Ofav_annot.xml -f 5
diamond view -a Ofav_annot.daa -o Ofav_annot.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

echo "STOP" $(date)

sbatch ofav_diamond_blastx.sh

# After script has run:
wc -l Ofav_annot.tab
10122 Ofav_annot.tab # only ~10000 hits out of ~36000 genes 

```

Submitted batch job 13460 - finished 

##### BLASTp

```
# Check number of protein seqs 
zgrep -c "^>" GCF_002042975.1_ofav_dov_v1_protein.faa
32587

nano ofav_diamond_blastp.sh

#!/bin/bash 
#SBATCH --job-name="diamond-blastp-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="ofav_diamond_blastp_out_error"
#SBATCH --output="ofav_diamond_blastp_out"
#SBATCH --exclusive

echo "START" $(date)
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

echo "Updating Ofav annotation" $(date)
diamond blastp -b 2 -d /data/putnamlab/shared/databases/nr.dmnd -q /data/putnamlab/jillashey/genome/Ofav/GCF_002042975.1_ofav_dov_v1_protein.faa -o Ofav_blastp_annot -f 100 -e 0.00001 -k 1 --threads $SLURM_CPUS_ON_NODE --tmpdir /tmp/

echo "Search complete... converting format to XML and tab"

diamond view -a Ofav_blastp_annot.daa -o Ofav_blastp_annot.xml -f 5
diamond view -a Ofav_blastp_annot.daa -o Ofav_blastp_annot.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

echo "STOP" $(date)

sbatch ofav_diamond_blastp.sh

# After script has run:
wc -l Ofav_blastp_annot.tab 
32536 Ofav_blastp_annot.tab
```

Submitted batch job 93778

### Hawaii species 

#### Mcap

##### BLASTx

```
# Check number of genes 
zgrep -c "^>" Mcap.mRNA.fa
63227

nano mcap_diamond_blastx.sh

#!/bin/bash
#SBATCH --job-name="diamond-blastx"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH -p putnamlab
#SBATCH --mem=100GB
#SBATCH --error="mcap_diamond_blastx_out_error"
#SBATCH --output="mcap_diamond_blastx_out"

echo "START" $(date)
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

echo "Updating Mcap annotation" $(date)
diamond blastx -d /data/putnamlab/shared/databases/nr.dmnd -q Mcap.mRNA.fa -o Mcap_annot -f 100 -b20 --more-sensitive -e 0.00001 -k1

echo "Search complete... converting format to XML and tab"

diamond view -a Mcap_annot.daa -o Mcap_annot.xml -f 5
diamond view -a Mcap_annot.daa -o Mcap_annot.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

echo "STOP" $(date)

sbatch mcap_diamond_blastx.sh

# After script has run:
wc -l Mcap_annot.tab
55217 Mcap_annot.tab # ~55000 hits out of ~63000 genes 
```

Submitted batch job 14933

##### BLASTp

```
# Check number of protein seqs 
zgrep -c "^>" Mcap.protein.fa
63227

nano mcap_diamond_blastp.sh

#!/bin/bash 
#SBATCH --job-name="diamond-blastp-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="mcap_diamond_blastp_out_error"
#SBATCH --output="mcap_diamond_blastp_out"
#SBATCH --exclusive

echo "START" $(date)
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

echo "Updating Mcap annotation" $(date)
diamond blastp -b 2 -d /data/putnamlab/shared/databases/nr.dmnd -q /data/putnamlab/jillashey/genome/Mcap/Mcap.protein.fa -o Mcap_blastp_annot -f 100 -e 0.00001 -k 1 --threads $SLURM_CPUS_ON_NODE --tmpdir /tmp/

echo "Search complete... converting format to XML and tab"

diamond view -a Mcap_blastp_annot.daa -o Mcap_blastp_annot.xml -f 5
diamond view -a Mcap_blastp_annot.daa -o Mcap_blastp_annot.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

echo "STOP" $(date)

sbatch mcap_diamond_blastp.sh

# After script has run:
wc -l Mcap_blastp_annot.tab
55683 Mcap_blastp_annot.tab
```

Submitted batch job 93780

#### Pcomp

##### BLASTx

```
# Check number of genes 
zgrep -c "^>" Porites_compressa_CDS.fa
74728

nano pcomp_diamond_blastx.sh

#!/bin/bash
#SBATCH --job-name="diamond-blastx"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH -p putnamlab
#SBATCH --mem=100GB
#SBATCH --error="pcomp_diamond_blastx_out_error"
#SBATCH --output="pcomp_diamond_blastx_out"

echo "START" $(date)
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

echo "Updating Pcomp annotation" $(date)
diamond blastx -d /data/putnamlab/shared/databases/nr.dmnd -q Porites_compressa_CDS.fa -o Pcomp_annot -f 100 -b20 --more-sensitive -e 0.00001 -k1

echo "Search complete... converting format to XML and tab"

diamond view -a Pcomp_annot.daa -o Pcomp_annot.xml -f 5
diamond view -a Pcomp_annot.daa -o Pcomp_annot.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

echo "STOP" $(date)

sbatch pcomp_diamond_blastx.sh

# After script has run:
wc -l Pcomp_annot.tab
59049 Pcomp_annot.tab # ~59000 hits out of ~75000 genes 
```

Submitted batch job 14932

##### BLASTp

```
# Check number of protein seqs 
zgrep -c "^>" Porites_compressa_AA.fa
74728

nano pcomp_diamond_blastp.sh

#!/bin/bash 
#SBATCH --job-name="diamond-blastp-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="pcomp_diamond_blastp_out_error"
#SBATCH --output="pcomp_diamond_blastp_out"
#SBATCH --exclusive

echo "START" $(date)
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

echo "Updating Pcomp annotation" $(date)
diamond blastp -b 2 -d /data/putnamlab/shared/databases/nr.dmnd -q /data/putnamlab/jillashey/genome/Pcomp/Porites_compressa_AA.fa -o Pcomp_blastp_annot -f 100 -e 0.00001 -k 1 --threads $SLURM_CPUS_ON_NODE --tmpdir /tmp/

echo "Search complete... converting format to XML and tab"

diamond view -a Pcomp_blastp_annot.daa -o Pcomp_blastp_annot.xml -f 5
diamond view -a Pcomp_blastp_annot.daa -o Pcomp_blastp_annot.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

echo "STOP" $(date)

sbatch pcomp_diamond_blastp.sh

# After script has run:
wc -l Pcomp_blastp_annot.tab
58796 Pcomp_blastp_annot.tab
```

Submitted batch job 93786

#### Pdam

##### BLASTx

###### NCBI

```
# Check number of genes 
zgrep -c "^>" GCF_003704095.1_ASM370409v1_rna.fna
27287

nano pdam_diamond_blastx.sh

#!/bin/bash
#SBATCH --job-name="diamond-blastx"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH -p putnamlab
#SBATCH --mem=100GB
#SBATCH --error="pdam_diamond_blastx_out_error"
#SBATCH --output="pdam_diamond_blastx_out"

echo "START" $(date)
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

echo "Updating Pdam annotation" $(date)
diamond blastx -d /data/putnamlab/shared/databases/nr.dmnd -q GCF_003704095.1_ASM370409v1_rna.fna -o Pdam_annot -f 100 -b20 --more-sensitive -e 0.00001 -k1

echo "Search complete... converting format to XML and tab"

diamond view -a Pdam_annot.daa -o Pdam_annot.xml -f 5
diamond view -a Pdam_annot.daa -o Pdam_annot.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

echo "STOP" $(date)

sbatch pdam_diamond_blastx.sh

# After script has run:
wc -l Pdam_annot.tab
3406 Pdam_annot.tab # only ~3400 hits out of ~27000 genes 
```

Submitted batch job 13728 - finished 

###### Reef Genomics

Running updated diamond database script

```
cd /data/putnamlab/shared/sbatch_executables/
sbatch make_diamond_nr_db.sh 
Submitted batch job 20794
```

```
cd /data/putnamlab/jillashey/annotation/diamond/pdam/ReefGenomics

ln -s /data/putnamlab/jillashey/genome/Pdam/ReefGenomics/pdam_transcripts.fasta .

# Check number of genes 
zgrep -c "^>" pdam_transcripts.fasta
26077

nano pdam_RG_diamond_blastx.sh

#!/bin/bash
#SBATCH --job-name="diamond-blastx"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH -p putnamlab
#SBATCH --mem=100GB
#SBATCH --error="pdam_RG_diamond_blastx_out_error"
#SBATCH --output="pdam_RG_diamond_blastx_out"

echo "START" $(date)
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

echo "Updating Pdam RG annotation" $(date)
diamond blastx -d /data/putnamlab/shared/databases/nr.dmnd -q pdam_transcripts.fasta -o Pdam_RG_annot -f 100 -b20 --more-sensitive -e 0.00001 -k1

echo "Search complete... converting format to XML and tab"

diamond view -a Pdam_RG_annot.daa -o Pdam_RG_annot.xml -f 5
diamond view -a Pdam_RG_annot.daa -o Pdam_RG_annot.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

echo "STOP" $(date)

sbatch pdam_RG_diamond_blastx.sh

# After script has run:
wc -l Pdam_RG_annot.tab
3732 Pdam_RG_annot.tab # only ~3700 hits out of ~26000 genes 
```
Submitted batch job 20795


##### BLASTp

###### NCBI

```
# Check number of protein seqs 
zgrep -c "^>" GCF_003704095.1_ASM370409v1_protein.faa
25183

nano pdam_NCBI_diamond_blastp.sh

#!/bin/bash 
#SBATCH --job-name="diamond-blastp-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="pdam_ncbi_diamond_blastp_out_error"
#SBATCH --output="pdam_ncbi_diamond_blastp_out"
#SBATCH --exclusive

echo "START" $(date)
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

echo "Updating Pdam NCBI annotation" $(date)
diamond blastp -b 2 -d /data/putnamlab/shared/databases/nr.dmnd -q /data/putnamlab/jillashey/genome/Pdam/NCBI/GCF_003704095.1_ASM370409v1_protein.faa -o Pdam_NCBI_blastp_annot -f 100 -e 0.00001 -k 1 --threads $SLURM_CPUS_ON_NODE --tmpdir /tmp/

echo "Search complete... converting format to XML and tab"

diamond view -a Pdam_NCBI_blastp_annot.daa -o Pdam_NCBI_blastp_annot.xml -f 5
diamond view -a Pdam_NCBI_blastp_annot.daa -o Pdam_NCBI_blastp_annot.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

echo "STOP" $(date)

sbatch pdam_NCBI_diamond_blastp.sh

# After script has run:
wc -l Pdam_NCBI_blastp_annot.tab
25170 Pdam_NCBI_blastp_annot.tab
```

Submitted batch job 93826

###### Reef Genomics 

```
# Check number of protein seqs 
zgrep -c "^>" pdam_proteins.fasta
26077

nano pdam_RG_diamond_blastp.sh

#!/bin/bash 
#SBATCH --job-name="diamond-blastp-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="pdam_rg_diamond_blastp_out_error"
#SBATCH --output="pdam_rg_diamond_blastp_out"
#SBATCH --exclusive

echo "START" $(date)
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

echo "Updating Pdam RG annotation" $(date)
diamond blastp -b 2 -d /data/putnamlab/shared/databases/nr.dmnd -q /data/putnamlab/jillashey/genome/Pdam/ReefGenomics/pdam_proteins.fasta -o Pdam_RG_blastp_annot -f 100 -e 0.00001 -k 1 --threads $SLURM_CPUS_ON_NODE --tmpdir /tmp/

echo "Search complete... converting format to XML and tab"

diamond view -a Pdam_RG_blastp_annot.daa -o Pdam_RG_blastp_annot.xml -f 5
diamond view -a Pdam_RG_blastp_annot.daa -o Pdam_RG_blastp_annot.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

echo "STOP" $(date)

sbatch pdam_RG_diamond_blastp.sh

# After script has run:
wc -l Pdam_RG_blastp_annot.tab
25947 Pdam_RG_blastp_annot.tab
```

Submitted batch job 93827

#### Plob

Using plutea genome, proteins, transcripts, etc 

##### BLASTx

```
# Check number of genes 
zgrep -c "^>" plut2v1.1.transcripts.fasta
31126

nano plob_diamond_blastx.sh

#!/bin/bash
#SBATCH --job-name="diamond-blastx"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH -p putnamlab
#SBATCH --mem=100GB
#SBATCH --error="plob_diamond_blastx_out_error"
#SBATCH --output="plob_diamond_blastx_out"

echo "START" $(date)
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

echo "Updating Plob annotation" $(date)
diamond blastx -d /data/putnamlab/shared/databases/nr.dmnd -q plut2v1.1.transcripts.fasta -o Plut_annot -f 100 -b20 --more-sensitive -e 0.00001 -k1

echo "Search complete... converting format to XML and tab"

diamond view -a Plut_annot.daa -o Plut_annot.xml -f 5
diamond view -a Plut_annot.daa -o Plut_annot.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

echo "STOP" $(date)

sbatch plob_diamond_blastx.sh

# After script has run: 
wc -l Plut_annot.tab
4315 Plut_annot.tab # only ~4300 hits out of ~31000 genes 

```

Submitted batch job 13729 - finished 

##### BLASTp

```
# Check number of protein seqs 
zgrep -c "^>" plut2v1.1.proteins.fasta
31126

nano plob_diamond_blastp.sh

#!/bin/bash 
#SBATCH --job-name="diamond-blastp-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="plob_diamond_blastp_out_error"
#SBATCH --output="plob_diamond_blastp_out"
#SBATCH --exclusive

echo "START" $(date)
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

echo "Updating Plob annotation" $(date)
diamond blastp -b 2 -d /data/putnamlab/shared/databases/nr.dmnd -q /data/putnamlab/jillashey/genome/Plutea/plut2v1.1.proteins.fasta -o Plut_blastp_annot -f 100 -e 0.00001 -k 1 --threads $SLURM_CPUS_ON_NODE --tmpdir /tmp/

echo "Search complete... converting format to XML and tab"

diamond view -a Plut_blastp_annot.daa -o Plut_blastp_annot.xml -f 5
diamond view -a Plut_blastp_annot.daa -o Plut_blastp_annot.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

echo "STOP" $(date)

sbatch plob_diamond_blastp.sh

# After script has run:
wc -l Plut_blastp_annot.tab
28284 Plut_blastp_annot.tab
```

Submitted batch job 93787

#### Pacuta

##### BLASTx

```
# Check number of genes 
zgrep -c "^>" braker_v1.codingseq.fasta
38913

nano pacuta_diamond_blastx.sh

#!/bin/bash
#SBATCH --job-name="diamond-blastx"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="pacuta_diamond_blastx_out_error"
#SBATCH --output="pacuta_diamond_blastx_out"

echo "START" $(date)
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

echo "Updating Pacuta annotation" $(date)
diamond blastx -d /data/putnamlab/shared/databases/nr.dmnd -q /data/putnamlab/jillashey/genome/Pacuta/braker_v1.codingseq.fasta -o Pacuta_blastx_annot -f 100 -b20 --more-sensitive -e 0.00001 -k1

echo "Search complete... converting format to XML and tab"

diamond view -a Pacuta_blastx_annot.daa -o Pacuta_blastx_annot.xml -f 5
diamond view -a Pacuta_blastx_annot.daa -o Pacuta_blastx_annot.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

echo "STOP" $(date)

sbatch pacuta_diamond_blastx.sh

# After script has run: 
wc -l Pacuta_annot.tab
XXXXX Pacuta_annot.tab

```

Submitted batch job 104687

##### BLASTp

```
# Check number of protein seqs 
zgrep -c "^>" Pocillopora_acuta_HIv1.genes.pep.faa
38913

nano pacuta_diamond_blastp.sh

#!/bin/bash 
#SBATCH --job-name="diamond-blastp-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="pacuta_diamond_blastp_out_error"
#SBATCH --output="pacuta_diamond_blastp_out"

echo "START" $(date)
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

echo "Updating Pacuta annotation" $(date)
diamond blastp -b 2 -d /data/putnamlab/shared/databases/nr.dmnd -q /data/putnamlab/jillashey/genome/Pacuta/Pocillopora_acuta_HIv1.genes.pep.faa -o Pacuta_blastp_annot -f 100 -e 0.00001 -k 1 --threads $SLURM_CPUS_ON_NODE --tmpdir /tmp/

echo "Search complete... converting format to XML and tab"

diamond view -a Pacuta_blastp_annot.daa -o Pacuta_blastp_annot.xml -f 5
diamond view -a Pacuta_blastp_annot.daa -o Pacuta_blastp_annot.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

echo "STOP" $(date)

sbatch pacuta_diamond_blastp.sh

# After script has run:
wc -l Pacuta_blastp_annot.tab
36945 Pacuta_blastp_annot.tab
```

Submitted batch job 104686

Copy all .tab files onto local computer to run through uniprot 

Copy all .xml files onto local computer to run through B2G. These files are too large to put in github repo 
