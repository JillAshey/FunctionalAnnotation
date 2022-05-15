## Swiss Prot

using protein seqs 

*All analyses done on Putnam Lab Node*


1) Before blasting, create a swiss prot reference database using the code below. Its in the shared/sbatch_excutables folder on Putnam node. This script was written by Danielle Becker.  

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
```

Submitted batch job 94130

I think it ran...it wouldn't let me change the email or date so i just ran it without changing anything.

2) Write and run script that will blast protein seqs against the Swiss Prot database. This script was written by Danielle Becker.  

### Florida species

#### Acerv 

##### BLASTp
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
blastp -max_target_seqs 5 -num_threads 20 -db /data/putnamlab/shared/databases/swiss_db/swissprot_20211022 -query /data/putnamlab/jillashey/genome/Acerv/Acerv_assembly_v1.0.protein.fa -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -out acerv_swissprot_protein.out

echo "STOP" $(date)

sbatch swissprot_blast.sh 

awk '{print $1}' acerv_swissprot_protein.out | sort | uniq | wc -l 
17519
```

Submitted batch job 94142

Get best hit for each swiss gene model (protein)

```
cat acerv_swissprot_protein.out | sort -k1,1 -k2,2 -k3,3r -k4,4r -k11,11 | awk '!seen[$1]++' > acerv_swissprot_protein_besthit.out

wc -l acerv_swissprot_protein_besthit.out
17519 acerv_swissprot_protein_besthit.out
```

Seeing if I can convert blastp output from out --> xml without rerunning whole thing

```
nano swissprot_blast_formatter.sh

#!/bin/bash 
#SBATCH --job-name="blastformat"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="blastformat_out_error"
#SBATCH --output="blastformat_out"

###### change module once running on andreamadokfosd putnam lab node 

echo "START" $(date)
module load BLAST+/2.9.0-gompi-2019b

blast_formatter -archive acerv_swissprot_protein_besthit.out -outfmt "16 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
echo "STOP" $(date)

SBATCH swissprot_blast_formatter.sh
```

Submitted batch job 1934058 --- JOB SUBMITTED ON BLUEWAVES. okay this failed: BLAST query/options error: Invalid input format for BLAST Archive.

Going to rerun ```swissprot_blast.sh``` for acerv, but change the -outfmt 6 to -outfmt 5 so it will give me xml output. Hopefully, I can input the xml file into B2G. 
Submitted batch job 94678

##### BLASTx

In my originial annotations, I did not compare my transcript seqs against the SwissProt database so im going to do that now. maybe it will give me better results in B2G bc im getting lots of b2g results from protein blast swissprot.

```
nano acerv_swissprot_blastx.sh

#!/bin/bash 
#SBATCH --job-name="swissprot-blastx-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="swissprot_blastx_out_error"
#SBATCH --output="swissprot_blastx_out"
#SBATCH --exclusive

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b #load blast module

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b #load blast module

echo "Blast against swissprot database" $(date)
blastx -max_target_seqs 5 -num_threads 20 -db /data/putnamlab/shared/databases/swiss_db/swissprot_20211022 -query /data/putnamlab/jillashey/genome/Acerv/Acerv_assembly_v1.0.mRNA.fa -evalue 1e-5 -outfmt '5 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -out acerv_swissprot_blastx

echo "STOP" $(date)

sbatch acerv_swissprot_blastx.sh 
```

Submitted batch job 95983




#### Mcav 
```
nano mcav_swissprot_blast.sh

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
blastp -max_target_seqs 5 -num_threads 20 -db /data/putnamlab/shared/databases/swiss_db/swissprot_20211022 -query /data/putnamlab/jillashey/genome/Mcav/Mcavernosa_annotation/Mcavernosa.maker.proteins.fasta -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -out mcav_swissprot_protein.out

echo "STOP" $(date)

sbatch mcav_swissprot_blast.sh 
```

Submitted batch job 94169

Get best hit for each swiss gene model (protein)

```
cat mcav_swissprot_protein.out | sort -k1,1 -k2,2 -k3,3r -k4,4r -k11,11 | awk '!seen[$1]++' > mcav_swissprot_protein_besthit.out

wc -l mcav_swissprot_protein_besthit.out
15953 mcav_swissprot_protein_besthit.out
```

Going to rerun ```mcav_swissprot_blast.sh ``` for mcav, but change the -outfmt 6 to -outfmt 5 so it will give me xml output. Hopefully, xml file can be input to B2G. 
Submitted batch job 94987

#### Ofav 
```
nano ofav_swissprot_blast.sh

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
blastp -max_target_seqs 5 -num_threads 20 -db /data/putnamlab/shared/databases/swiss_db/swissprot_20211022 -query /data/putnamlab/jillashey/genome/Ofav/GCF_002042975.1_ofav_dov_v1_protein.faa -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -out ofav_swissprot_protein.out

echo "STOP" $(date)

sbatch ofav_swissprot_blast.sh 
```

Submitted batch job 94171

Get best hit for each swiss gene model (protein)

```
cat ofav_swissprot_protein.out | sort -k1,1 -k2,2 -k3,3r -k4,4r -k11,11 | awk '!seen[$1]++' > ofav_swissprot_protein_besthit.out

wc -l ofav_swissprot_protein_besthit.out
25296 ofav_swissprot_protein_besthit.out
```

Going to rerun ```ofav_swissprot_blast.sh ``` for ofav, but change the -outfmt 6 to -outfmt 5 so it will give me xml output. Hopefully, xml file can be input to B2G. 
Submitted batch job 94989


### Hawaii species

#### Mcap

```
nano mcap_swissprot_blast.sh

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
blastp -max_target_seqs 5 -num_threads 20 -db /data/putnamlab/shared/databases/swiss_db/swissprot_20211022 -query /data/putnamlab/shared/databases/nr.dmnd -q /data/putnamlab/jillashey/genome/Mcap/Mcap.protein.fa -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -out mcap_swissprot_protein.out

echo "STOP" $(date)

sbatch mcap_swissprot_blast.sh 
```

Submitted batch job 94183

Get best hit for each swiss gene model (protein)

```
cat mcap_swissprot_protein.out | sort -k1,1 -k2,2 -k3,3r -k4,4r -k11,11 | awk '!seen[$1]++' > mcap_swissprot_protein_besthit.out

wc -l mcap_swissprot_protein_besthit.out
27180 mcap_swissprot_protein_besthit.out
```

Going to rerun ```mcap_swissprot_blast.sh ``` for mcap, but change the -outfmt 6 to -outfmt 5 so it will give me xml output. Hopefully, xml file can be input to B2G. Submitted batch job 94988


#### Pcomp

```
nano pcomp_swissprot_blast.sh

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
blastp -max_target_seqs 5 -num_threads 20 -db /data/putnamlab/shared/databases/swiss_db/swissprot_20211022 -query /data/putnamlab/jillashey/genome/Pcomp/Porites_compressa_AA.fa -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -out pcomp_swissprot_protein.out

echo "STOP" $(date)

sbatch pcomp_swissprot_blast.sh 
```

 Submitted batch job 94173
 
 Get best hit for each swiss gene model (protein)
 
 ```
cat pcomp_swissprot_protein.out | sort -k1,1 -k2,2 -k3,3r -k4,4r -k11,11 | awk '!seen[$1]++' > pcomp_swissprot_protein_besthit.out

wc -l pcomp_swissprot_protein_besthit.out
35041 pcomp_swissprot_protein_besthit.out
```

Going to rerun ```pcomp_swissprot_blast.sh ``` for pcomp, but change the -outfmt 6 to -outfmt 5 so it will give me xml output. Hopefully, xml file can be input to B2G. Submitted batch job 94990


#### Pdam

###### NCBI 

```
nano pdam_NCBI_swissprot_blast.sh

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
blastp -max_target_seqs 5 -num_threads 20 -db /data/putnamlab/shared/databases/swiss_db/swissprot_20211022 -query /data/putnamlab/jillashey/genome/Pdam/NCBI/GCF_003704095.1_ASM370409v1_protein.faa -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -out pdam_NCBI_swissprot_protein.out

echo "STOP" $(date)

sbatch pdam_NCBI_swissprot_blast.sh 
```

Submitted batch job 94174

Get best hit for each swiss gene model (protein)

```
cat pdam_NCBI_swissprot_protein.out | sort -k1,1 -k2,2 -k3,3r -k4,4r -k11,11 | awk '!seen[$1]++' > pdam_NCBI_swissprot_protein_besthit.out

wc -l pdam_NCBI_swissprot_protein_besthit.out
20369 pdam_NCBI_swissprot_protein_besthit.out
```

Going to rerun ```pdam_NCBI_swissprot_blast.sh ``` for pdam_NCBI, but change the -outfmt 6 to -outfmt 5 so it will give me xml output. Hopefully, xml file can be input to B2G. Submitted batch job 94991


###### Reef Genomics  

```
nano pdam_RG_swissprot_blast.sh

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
blastp -max_target_seqs 5 -num_threads 20 -db /data/putnamlab/shared/databases/swiss_db/swissprot_20211022 -query /data/putnamlab/jillashey/genome/Pdam/ReefGenomics/pdam_proteins.fasta -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -out pdam_RG_swissprot_protein.out

echo "STOP" $(date)

sbatch pdam_RG_swissprot_blast.sh 
```

Submitted batch job 94175

Get best hit for each swiss gene model (protein)

```
cat pdam_RG_swissprot_protein.out | sort -k1,1 -k2,2 -k3,3r -k4,4r -k11,11 | awk '!seen[$1]++' > pdam_RG_swissprot_protein_besthit.out

wc -l pdam_RG_swissprot_protein_besthit.out
15692 pdam_RG_swissprot_protein_besthit.out
```

Going to rerun ```pdam_RG_swissprot_blast.sh ``` for pdam_RG, but change the -outfmt 6 to -outfmt 5 so it will give me xml output. Hopefully, xml file can be input to B2G. Submitted batch job 94992


#### Plob

```
nano plob_swissprot_blast.sh

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
blastp -max_target_seqs 5 -num_threads 20 -db /data/putnamlab/shared/databases/swiss_db/swissprot_20211022 -query /data/putnamlab/jillashey/genome/Plutea/plut2v1.1.proteins.fasta -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -out plob_swissprot_protein.out

echo "STOP" $(date)

sbatch plob_swissprot_blast.sh 
```

Submitted batch job 94176

Get best hit for each swiss gene model (protein)

```
plob]$ cat plob_swissprot_protein.out | sort -k1,1 -k2,2 -k3,3r -k4,4r -k11,11 | awk '!seen[$1]++' > plob_swissprot_protein_besthit.out

wc -l plob_swissprot_protein_besthit.out
21130 plob_swissprot_protein_besthit.out
```

Going to rerun ```plob_swissprot_blast.sh ``` for plob, but change the -outfmt 6 to -outfmt 5 so it will give me xml output. Hopefully, xml file can be input to B2G. Submitted batch job 94994

#### Pacuta 

```
nano pacuta_swissprot_blast.sh

#!/bin/bash 
#SBATCH --job-name="swissprot-blastp-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="swissprot_blastp_out_error"
#SBATCH --output="swissprot_blastp_out"

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b #load blast module

echo "Blast against swissprot database" $(date)
blastp -max_target_seqs 5 -num_threads 20 -db /data/putnamlab/shared/databases/swiss_db/swissprot_20211022 -query /data/putnamlab/jillashey/genome/Pacuta/Pocillopora_acuta_HIv1.genes.pep.faa -evalue 1e-5 -outfmt '5 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -out pacuta_swissprot_protein.xml

echo "STOP" $(date)

sbatch pacuta_swissprot_blast.sh 
```

Submitted batch job 105213


