### P. meandrina functional annotation 

I will be doing the functional annotation for the Pmeand genome. 

Make folders for this species in my own directory

```
cd /data/putnamlab/jillashey/annotation
cd diamond
mkdir pmea
cd ../InterProScan
mkdir pmea
cd ../swiss_prot
mkdir pmea
```


In my genome folder, make a spot for Pmea.

```
cd /data/putnamlab/jillashey/genome/
mkdir Pmea
cd Pmea
```

Download genomic information from P. meadrina from Hawaii

```
# protein fasta 
wget http://cyanophora.rutgers.edu/Pocillopora_meandrina/Pocillopora_meandrina_HIv1.genes.pep.faa.gz
zgrep -c ">" Pocillopora_meandrina_HIv1.genes.pep.faa
31840

# nucleotide fasta 
wget http://cyanophora.rutgers.edu/Pocillopora_meandrina/Pocillopora_meandrina_HIv1.genes.cds.fna.gz

# GFF 
wget http://cyanophora.rutgers.edu/Pocillopora_meandrina/Pocillopora_meandrina_HIv1.genes.gff3.gz

gunzip *
```

#### Diamond BLAST 

Before running BLAST, update the BLAST database on the Putnam lab shared directory. 

```
cd /data/putnamlab/shared/sbatch_executables
sbatch download_nr_database.sh
Submitted batch job 274384



sbatch make_diamond_nr_db.sh
```

Run BLASTp using protein sequences from NCBI

AS OF 8/16/23, HAVING TROUBLE UPDATING THE NR DATABASE ON ANDROMEDA 

Blast against nr database 

```

```

#### SwissProt  

Before running SwissProt, update the SwissProt database on the Putnam lab shared directory. Before running, change the email and the date that the database was downloaded in the code. 

```
cd /data/putnamlab/shared/sbatch_executables
sbatch download_swissprot_database.sh
```

Submitted batch job 274389
ISSUES DOWNLOADING SP DATBASE TOO 

Blast against swissprot database 

```
nano pmea_swissprot_blast.sh

#!/bin/bash 
#SBATCH --job-name="swissprot-blastp-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --mem=200GB
#SBATCH --error="swissprot_blastp_out_error"
#SBATCH --output="swissprot_blastp_out"
#SBATCH --exclusive

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b #load blast module
### ADD NEWER BLAST MODULE

echo "Blast against swissprot database" $(date)
blastp -max_target_seqs 5 -num_threads 20 -db XXXXXXXXX -query XXXXXXXXX -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -out plob_swissprot_protein.out

echo "STOP" $(date)

sbatch pmea_swissprot_blast.sh 
```

#### Interproscan 

`cd /data/putnamlab/jillashey/annotation/InterProScan/pmea`

```
nano IPS_pmea.sh

#!/bin/bash
#SBATCH --job-name="IPS_pacuta"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --mem=200GB
#SBATCH --error="pmea_interproscan_out_error"
#SBATCH --output="pmea_interproscan_out"
#SBATCH --exclusive

cd /data/putnamlab/jillashey/annotation/InterProScan/pmea

echo "START $(date)"

# Load modules
module load InterProScan_data/5.60-92.0-foss-2021b
module load Java/11.0.2
java -version

# Run InterProScan
interproscan.sh -version
interproscan.sh -f XML -i /data/putnamlab/jillashey/genome/Pmea/Pocillopora_meandrina_HIv1.genes.pep.faa -b pmea.interpro -iprlookup -goterms -pa 
interproscan.sh -mode convert -f GFF3 -i pmea.interpro.xml -b pmea.interpro

echo "DONE $(date)"
```

`sbatch IPS_pmea.sh`; Submitted batch job 274397. This might take a long time to run.
