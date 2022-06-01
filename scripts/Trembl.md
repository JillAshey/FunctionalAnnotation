## Trembl

Protein database to blast protein seqs 

*All analyses done on Putnam Lab Node*


1) Before blasting, create a trembl reference database using the code below. Its in the shared/sbatch_excutables folder on Putnam node. This script was written by Danielle Becker.  

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
```

Submitted batch job 94144

I think it ran...it wouldn't let me change the email or date so i just ran it without changing anything.

2) Write and run script that will blast protein seqs against the Trembl database. This script was written by Danielle Becker.  

I was having SO MUCH trouble running trembl. it was just taking forever and it didn't seem like it was even running at all. But Ariana and I were able to figure out that if we break up the protein seq file and run it as an array job, we can run trembl successfully in a reasonable amount of time. Check out code for Mcap [here](https://github.com/JillAshey/FunctionalAnnotation/blob/main/Trembl/Mcap_Trembl.md) where we figured it out!

### Florida species

#### Acerv

Splitting Acerv protein sequence file into multiple files so that trembl doesn't have to run one big file - it can run several smaller files and hopefully that will make it go faster.

On bluewaves (since pyfasta is only on bluewaves and i am doing this over the weekend so I don't want to bother Kevin Bryan :-) ):

```
cd /data/putnamlab/jillashey/genome/Acerv

zgrep -c ">" Acerv_assembly_v1.0.protein.fa
33322
```

[PyFasta](https://github.com/brentp/pyfasta) will be able to split a fasta file into several new files of relatively even size

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
```

Submitted batch job 135768. Started on May 6, ended on May 13. Xml files are between 17-18 MB in size.


#### Mcav

Splitting Mcav protein sequence file into multiple files so that trembl doesn't have to run one big file - it can run several smaller files and hopefully that will make it go faster.

On bluewaves (since pyfasta is only on bluewaves and i am doing this over the weekend so I don't want to bother Kevin Bryan :-) ):

```
cd /data/putnamlab/jillashey/genome/Mcav/Mcavernosa_annotation

zgrep -c ">" Mcavernosa.maker.proteins.fasta
25142
```

[PyFasta](https://github.com/brentp/pyfasta) will be able to split a fasta file into several new files of relatively even size

```
module load pyfasta/0.5.2

pyfasta split -n 4 Mcavernosa.maker.proteins.fasta
creating new files:
Mcavernosa.maker.proteins.0.fasta
Mcavernosa.maker.proteins.1.fasta
Mcavernosa.maker.proteins.2.fasta
Mcavernosa.maker.proteins.3.fasta
```

Copy files into trembl folder

```
cp Mcavernosa.maker.proteins.*.fasta /data/putnamlab/jillashey/annotation/trembl/mcav
```

Run trembl as an array job on ANDROMEDA 

```
nano mcav_trembl_blastp.sh

#!/bin/bash 
#SBATCH --job-name="trembl-blastp-protein"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --mem=120GB
#SBATCH --error="mcav_trembl_blastp_out_error"
#SBATCH --output="mcav_trembl_blastp_out"
#SBATCH --exclusive

echo "START" $(date)

module load BLAST+/2.11.0-gompi-2020b #load blast module

#F=/data/putnamlab/jillashey/annotation/trembl/mcav

echo "Blast against trembl database" $(date)

array1=($(ls Mcavernosa.maker.proteins.*.fasta))
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

sbatch --array=0-3 mcav_trembl_blastp.sh
```

Submitted batch job 138212. Started on May 15, ended on May 23. Xml files are between 21-23 MB in size.

#### Ofav

Splitting Ofav protein sequence file into multiple files so that trembl doesn't have to run one big file - it can run several smaller files and hopefully that will make it go faster.

On bluewaves (since pyfasta is only on bluewaves and i am doing this over the weekend so I don't want to bother Kevin Bryan :-) ):

```
cd /data/putnamlab/jillashey/genome/Ofav

zgrep -c ">" GCF_002042975.1_ofav_dov_v1_protein.faa
32587
```

[PyFasta](https://github.com/brentp/pyfasta) will be able to split a fasta file into several new files of relatively even size

```
module load pyfasta/0.5.2

pyfasta split -n 6 GCF_002042975.1_ofav_dov_v1_protein.faa
creating new files:
GCF_002042975.1_ofav_dov_v1_protein.0.faa
GCF_002042975.1_ofav_dov_v1_protein.1.faa
GCF_002042975.1_ofav_dov_v1_protein.2.faa
GCF_002042975.1_ofav_dov_v1_protein.3.faa
GCF_002042975.1_ofav_dov_v1_protein.4.faa
GCF_002042975.1_ofav_dov_v1_protein.5.faa
```

Copy files into trembl folder

```
cp GCF_002042975.1_ofav_dov_v1_protein.*.faa /data/putnamlab/jillashey/annotation/trembl/ofav
```

Run trembl as an array job on ANDROMEDA 

```
nano ofav_trembl_blastp.sh

#!/bin/bash 
#SBATCH --job-name="trembl-blastp-protein"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --mem=120GB
#SBATCH --error="ofav_trembl_blastp_out_error"
#SBATCH --output="ofav_trembl_blastp_out"
#SBATCH --exclusive

echo "START" $(date)

module load BLAST+/2.11.0-gompi-2020b #load blast module

#F=/data/putnamlab/jillashey/annotation/trembl/ofav

echo "Blast against trembl database" $(date)

array1=($(ls GCF_002042975.1_ofav_dov_v1_protein.*.faa))
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

sbatch --array=0-5 ofav_trembl_blastp.sh
```

Submitted batch job 140694. Started on May 23, ended on May 30. Xml files are between 21-23 MB in size.

### Hawaii species

#### Mcap

code [here](https://github.com/JillAshey/FunctionalAnnotation/blob/main/Trembl/Mcap_Trembl.md)

#### Pacuta

Splitting Pacuta protein sequence file into multiple files so that trembl doesn't have to run one big file - it can run several smaller files and hopefully that will make it go faster.

On bluewaves (since pyfasta is only on bluewaves and i am doing this over the weekend so I don't want to bother Kevin Bryan :-) ):

```
cd /data/putnamlab/jillashey/genome/Pacuta

zgrep -c ">" Pocillopora_acuta_HIv1.genes.pep.faa
38913
```

[PyFasta](https://github.com/brentp/pyfasta) will be able to split a fasta file into several new files of relatively even size

```
module load pyfasta/0.5.2

pyfasta split -n 6 Pocillopora_acuta_HIv1.genes.pep.faa
creating new files:
Pocillopora_acuta_HIv1.genes.pep.0.faa
Pocillopora_acuta_HIv1.genes.pep.1.faa
Pocillopora_acuta_HIv1.genes.pep.2.faa
Pocillopora_acuta_HIv1.genes.pep.3.faa
Pocillopora_acuta_HIv1.genes.pep.4.faa
Pocillopora_acuta_HIv1.genes.pep.5.faa
```

Copy files into trembl folder

```
cp Pocillopora_acuta_HIv1.genes.pep.*.faa /data/putnamlab/jillashey/annotation/trembl/pacuta
```

Run trembl as an array job on ANDROMEDA 

```
nano pacuta_trembl_blastp.sh

#!/bin/bash 
#SBATCH --job-name="trembl-blastp-protein"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --mem=120GB
#SBATCH --error="pacuta_trembl_blastp_out_error"
#SBATCH --output="pacuta_trembl_blastp_out"
#SBATCH --exclusive

echo "START" $(date)

module load BLAST+/2.11.0-gompi-2020b #load blast module

#F=/data/putnamlab/jillashey/annotation/trembl/pacuta

echo "Blast against trembl database" $(date)

array1=($(ls Pocillopora_acuta_HIv1.genes.pep.*.faa))
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

sbatch --array=0-5 pacuta_trembl_blastp.sh
```

Submitted batch job 143162

#### Plob

Splitting Plob protein sequence file into multiple files so that trembl doesn't have to run one big file - it can run several smaller files and hopefully that will make it go faster.

On bluewaves (since pyfasta is only on bluewaves and i am doing this over the weekend so I don't want to bother Kevin Bryan :-) ):

```
cd /data/putnamlab/jillashey/genome/Plutea

zgrep -c ">" plut2v1.1.proteins.fasta
31126
```

[PyFasta](https://github.com/brentp/pyfasta) will be able to split a fasta file into several new files of relatively even size

```
module load pyfasta/0.5.2

pyfasta split -n 4 plut2v1.1.proteins.fasta
creating new files:
plut2v1.1.proteins.0.fasta
plut2v1.1.proteins.1.fasta
plut2v1.1.proteins.2.fasta
plut2v1.1.proteins.3.fasta
```

Copy files into trembl folder

```
cp plut2v1.1.proteins.*.fasta /data/putnamlab/jillashey/annotation/trembl/plob
```

Run trembl as an array job on ANDROMEDA 

```
nano plob_trembl_blastp.sh

#!/bin/bash 
#SBATCH --job-name="trembl-blastp-protein"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --mem=120GB
#SBATCH --error="plob_trembl_blastp_out_error"
#SBATCH --output="plob_trembl_blastp_out"
#SBATCH --exclusive

echo "START" $(date)

module load BLAST+/2.11.0-gompi-2020b #load blast module

#F=/data/putnamlab/jillashey/annotation/trembl/plob

echo "Blast against trembl database" $(date)

array1=($(ls plut2v1.1.proteins.*.fasta))
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

sbatch --array=0-4 plob_trembl_blastp.sh
```

Submitted batch job 143168
