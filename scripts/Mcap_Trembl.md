Splitting Mcap protein sequence file into multiple files so that trembl doesn't have to run one big file - it can run several smaller files and hopefully that will make it go faster. 

On Bluewaves: 

```
cd /data/putnamlab/jillashey/genome/Mcap
```

[PyFasta](https://github.com/brentp/pyfasta) will be able to split a fasta file into several new files of relatively even size

```
module load pyfasta/0.5.2

pyfasta split -n 6 Montipora_capitata_HIv2.genes.pep.faa
```

Copy files into trembl folder 

```
cd /data/putnamlab/jillashey/genome/Mcap
cp Montipora_capitata_HIv2.genes.pep.*.faa /data/putnamlab/jillashey/annotation/trembl/mcap
```

run trembl with multiple files 

```
nano mcap_v2_trembl_blastp.sh

#!/bin/bash 
#SBATCH --job-name="trembl-blastp-protein"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="trembl_blastp_out_error"
#SBATCH --output="trembl_blastp_out"
#SBATCH --exclusive

echo "START" $(date)

module load BLAST+/2.11.0-gompi-2020b #load blast module

#F=/data/putnamlab/jillashey/annotation/trembl/mcap

echo "Blast against trembl database" $(date)

array1=($(ls Montipora_capitata_HIv2.genes.pep.*.faa))
for i in ${array1[@]}
do
blastp -max_target_seqs 1 \
-num_threads 20 \
-db /data/putnamlab/shared/databases/trembl_db/trembl_20220307 \
-query ${i} \
-evalue 1e-5 \
-outfmt '5 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' \
-out ${i}.
done

echo "STOP" $(date)

sbatch --array=1-6 mcap_v2_trembl_blastp.sh
```

Submitted batch job 130986

scancel 130986 on 4/26/22 at 3:30 pm EST

Got this email from Kevin Bryan

I'm sorry I'm just noticing this now, but I don't think this job array is doing what you want, for a few reasons. First, all of the jobs appear to be processing the whole set of files. Second, because of that, they are all overwriting the output files (both the log file and the "Montipora_capitata_HIv2.genes.pep.0.faa." file. Thirdly, you are specifying  -num_threads 20, instead of -num_threads $SLURM_CPUS_ON_NODE, which would scale to the size of the node the job is placed on (24 or 36).

If you want to run this as an array so that all of the files are processed in their own job, you can do this by replacing:

for i in ${array1[@]}

with

for i in ${array1[${SLURM_ARRAY_TASK_ID}]}

And then submit it with sbatch --array=0-5 so the array indexes match.

I suggest you cancel and re-run this job with these modifications, since there is a good chance the output file is corrupted by having several jobs write to it.

It's formatted a little weird in Markdown.

okay so going to edit our script to include modifications Kevin:

```
nano mcap_v2_trembl_blastp.sh

#!/bin/bash 
#SBATCH --job-name="trembl-blastp-protein"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --mem=120GB
#SBATCH --error="trembl_blastp_out_error"
#SBATCH --output="trembl_blastp_out"
#SBATCH --exclusive

echo "START" $(date)

module load BLAST+/2.11.0-gompi-2020b #load blast module

#F=/data/putnamlab/jillashey/annotation/trembl/mcap

echo "Blast against trembl database" $(date)

array1=($(ls Montipora_capitata_HIv2.genes.pep.*.faa))
for i in ${array1[${SLURM_ARRAY_TASK_ID}]}
do
blastp -max_target_seqs 1 \
-num_threads $SLURM_CPUS_ON_NODE \
-db /data/putnamlab/shared/databases/trembl_db/trembl_20220307 \
-query ${i} \
-evalue 1e-5 \
-outfmt '5 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' \
-out ${i}.
done

echo "STOP" $(date)

sbatch --array=0-5 mcap_v2_trembl_blastp.sh
```

Submitted batch job 132577

Omg success!!!!!!!!! wow so the trick is to break the file into chunks and run each chunk individually/in an array. Each chunk was ~9000 sequences and the job took about 9 days to run. Xml files are between 24-25 MB in size.

```
cd/data/putnamlab/jillashey/annotation/trembl/mcap

ls -l -h
total 169M
-rw-r--r-- 1 jillashey putnamlab  873 Apr 26 16:24 mcap_v2_trembl_blastp.sh
-rwxr--r-- 1 jillashey putnamlab 3.6M Apr 26 16:22 Montipora_capitata_HIv2.genes.pep.0.faa
-rw-r--r-- 1 jillashey putnamlab  25M May  4 14:01 Montipora_capitata_HIv2.genes.pep.0.faa.
-rwxr--r-- 1 jillashey putnamlab 3.6M Apr 26 16:22 Montipora_capitata_HIv2.genes.pep.1.faa
-rw-r--r-- 1 jillashey putnamlab  26M May  4 18:55 Montipora_capitata_HIv2.genes.pep.1.faa.
-rwxr--r-- 1 jillashey putnamlab 3.6M Apr 26 16:22 Montipora_capitata_HIv2.genes.pep.2.faa
-rw-r--r-- 1 jillashey putnamlab  25M May  4 04:42 Montipora_capitata_HIv2.genes.pep.2.faa.
-rwxr--r-- 1 jillashey putnamlab 3.6M Apr 26 16:22 Montipora_capitata_HIv2.genes.pep.3.faa
-rw-r--r-- 1 jillashey putnamlab  25M May  4 04:27 Montipora_capitata_HIv2.genes.pep.3.faa.
-rwxr--r-- 1 jillashey putnamlab 3.6M Apr 26 16:22 Montipora_capitata_HIv2.genes.pep.4.faa
-rw-r--r-- 1 jillashey putnamlab  24M May  4 21:27 Montipora_capitata_HIv2.genes.pep.4.faa.
-rwxr--r-- 1 jillashey putnamlab 3.6M Apr 26 16:22 Montipora_capitata_HIv2.genes.pep.5.faa
-rw-r--r-- 1 jillashey putnamlab  25M May  4 01:57 Montipora_capitata_HIv2.genes.pep.5.faa.
-rw-r--r-- 1 jillashey putnamlab  668 May  4 21:27 trembl_blastp_out
-rw-r--r-- 1 jillashey putnamlab  305 Apr 26 16:53 trembl_blastp_out_error
```