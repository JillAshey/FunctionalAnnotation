### P. meandrina functional annotation 

I will be doing the functional annotation for the *Pocillopora meandrina* genome. 

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
ISSUES DOWNLOADING SP DATBASE TOO. Trying again to download on 8/23/23. Submitted batch job 275700. It appears to be running! The only thing I did differently was remove the path information behind the error and the output file. Weird. Appeared to have worked! It look about a minute which is a little weird. But nothing looks weird in the output or error files. Hooray!

In the annotation swiss prot directory, make a folder for Pmea. 

```
cd /data/putnamlab/jillashey/annotation/swiss_prot
mkdir pmea
cd pmea
```

Blast against swissprot database. 

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

blastp -max_target_seqs 5 -num_threads 20 -db /data/putnamlab/shared/databases/swiss_db/swissprot_20230816 -query /data/putnamlab/jillashey/genome/Pmea/Pocillopora_meandrina_HIv1.genes.pep.faa -evalue 1e-5 -outfmt '5 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -out pmea_swissprot_protein.out

echo "STOP" $(date)
```

Run so that `-outfmt` is 5 (XML output file)

Submitted batch job 275702

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

Finally ran. But did not generate an XML output file... there were a lot of errors in both the error and output files. 

```
tail pmea_interproscan_out
15:20:04.173 [responseMonitorJmsContainer-1] ERROR org.hibernate.engine.jdbc.spi.SqlExceptionHelper - Table "STEP_EXECUTION" not found; SQL statement:
select stepexecut0_.id as id1_109_0_, stepexecut0_.time_completed as time_com2_109_0_, stepexecut0_.time_created as time_cre3_109_0_, stepexecut0_.exception_first_chunk as exceptio4_109_0_, stepexecut0_.proportion_completed as proporti5_109_0_, stepexecut0_.time_started_running as time_sta6_109_0_, stepexecut0_.state as state7_109_0_, stepexecut0_.step_instance_id as step_ins9_109_0_, stepexecut0_.time_submitted as time_sub8_109_0_, exceptionc1_.step_execution_id as step_exe1_15_1_, exceptionc1_.exception_chunks as exceptio2_15_1_, exceptionc1_.chunk_index as chunk_in3_1_, stepinstan2_.id as id1_110_2_, stepinstan2_.bottom_model as bottom_m2_110_2_, stepinstan2_.bottom_protein as bottom_p3_110_2_, stepinstan2_.step_id as step_id4_110_2_, stepinstan2_.time_created as time_cre5_110_2_, stepinstan2_.top_model as top_mode6_110_2_, stepinstan2_.top_protein as top_prot7_110_2_, dependsupo3_.step_instance_id as step_ins1_112_3_, stepinstan4_.id as depends_2_112_3_, stepinstan4_.id as id1_110_4_, stepinstan4_.bottom_model as bottom_m2_110_4_, stepinstan4_.bottom_protein as bottom_p3_110_4_, stepinstan4_.step_id as step_id4_110_4_, stepinstan4_.time_created as time_cre5_110_4_, stepinstan4_.top_model as top_mode6_110_4_, stepinstan4_.top_protein as top_prot7_110_4_, parameters5_.step_instance_id as step_ins1_111_5_, parameters5_.parameters as paramete2_111_5_, parameters5_.parameters_key as paramete3_5_ from public.step_execution stepexecut0_ left outer join public.exception_chunk exceptionc1_ on stepexecut0_.id=exceptionc1_.step_execution_id inner join public.step_instance stepinstan2_ on stepexecut0_.step_instance_id=stepinstan2_.id left outer join public.step_instance_step_instance dependsupo3_ on stepinstan2_.id=dependsupo3_.step_instance_id left outer join public.step_instance stepinstan4_ on dependsupo3_.depends_upon_id=stepinstan4_.id left outer join public.step_instance_parameters parameters5_ on stepinstan4_.id=parameters5_.step_instance_id where stepexecut0_.id=? [42102-199]
15:20:04.175 [responseMonitorJmsContainer-1] ERROR org.hibernate.engine.jdbc.spi.SqlExceptionHelper - Table "STEP_EXECUTION" not found; SQL statement:
select stepexecut0_.id as id1_109_0_, stepexecut0_.time_completed as time_com2_109_0_, stepexecut0_.time_created as time_cre3_109_0_, stepexecut0_.exception_first_chunk as exceptio4_109_0_, stepexecut0_.proportion_completed as proporti5_109_0_, stepexecut0_.time_started_running as time_sta6_109_0_, stepexecut0_.state as state7_109_0_, stepexecut0_.step_instance_id as step_ins9_109_0_, stepexecut0_.time_submitted as time_sub8_109_0_, exceptionc1_.step_execution_id as step_exe1_15_1_, exceptionc1_.exception_chunks as exceptio2_15_1_, exceptionc1_.chunk_index as chunk_in3_1_, stepinstan2_.id as id1_110_2_, stepinstan2_.bottom_model as bottom_m2_110_2_, stepinstan2_.bottom_protein as bottom_p3_110_2_, stepinstan2_.step_id as step_id4_110_2_, stepinstan2_.time_created as time_cre5_110_2_, stepinstan2_.top_model as top_mode6_110_2_, stepinstan2_.top_protein as top_prot7_110_2_, dependsupo3_.step_instance_id as step_ins1_112_3_, stepinstan4_.id as depends_2_112_3_, stepinstan4_.id as id1_110_4_, stepinstan4_.bottom_model as bottom_m2_110_4_, stepinstan4_.bottom_protein as bottom_p3_110_4_, stepinstan4_.step_id as step_id4_110_4_, stepinstan4_.time_created as time_cre5_110_4_, stepinstan4_.top_model as top_mode6_110_4_, stepinstan4_.top_protein as top_prot7_110_4_, parameters5_.step_instance_id as step_ins1_111_5_, parameters5_.parameters as paramete2_111_5_, parameters5_.parameters_key as paramete3_5_ from public.step_execution stepexecut0_ left outer join public.exception_chunk exceptionc1_ on stepexecut0_.id=exceptionc1_.step_execution_id inner join public.step_instance stepinstan2_ on stepexecut0_.step_instance_id=stepinstan2_.id left outer join public.step_instance_step_instance dependsupo3_ on stepinstan2_.id=dependsupo3_.step_instance_id left outer join public.step_instance stepinstan4_ on dependsupo3_.depends_upon_id=stepinstan4_.id left outer join public.step_instance_parameters parameters5_ on stepinstan4_.id=parameters5_.step_instance_id where stepexecut0_.id=? [42102-199]
2023-08-23 15:32:37,156 [main] [uk.ac.ebi.interpro.scan.jms.master.AbstractMaster:259] WARN - Master process unable to delete temporary directory /data/putnamlab/jillashey/annotation/InterProScan/pmea/temp/n096.cluster.com_20230821_095649745_gq4y
23/08/2023 15:32:46:656 Welcome to InterProScan-5.60-92.0
23/08/2023 15:32:46:675 Running InterProScan v5 in CONVERT mode... on Linux
For the (-i) option you specified a location which doesn't exist or is not readable:
/data/putnamlab/jillashey/annotation/InterProScan/pmea/pmea.interpro.xml
DONE Wed Aug 23 15:32:54 EDT 2023

tail pmea_interproscan_out_error 
	at org.springframework.jms.listener.AbstractMessageListenerContainer.doInvokeListener(AbstractMessageListenerContainer.java:761)
	at org.springframework.jms.listener.AbstractMessageListenerContainer.invokeListener(AbstractMessageListenerContainer.java:699)
	at org.springframework.jms.listener.AbstractMessageListenerContainer.doExecuteListener(AbstractMessageListenerContainer.java:674)
	at org.springframework.jms.listener.AbstractPollingMessageListenerContainer.doReceiveAndExecute(AbstractPollingMessageListenerContainer.java:318)
	at org.springframework.jms.listener.AbstractPollingMessageListenerContainer.receiveAndExecute(AbstractPollingMessageListenerContainer.java:257)
	at org.springframework.jms.listener.DefaultMessageListenerContainer$AsyncMessageListenerInvoker.invokeListener(DefaultMessageListenerContainer.java:1186)
	at org.springframework.jms.listener.DefaultMessageListenerContainer$AsyncMessageListenerInvoker.executeOngoingLoop(DefaultMessageListenerContainer.java:1176)
	at org.springframework.jms.listener.DefaultMessageListenerContainer$AsyncMessageListenerInvoker.run(DefaultMessageListenerContainer.java:1073)
	at java.base/java.lang.Thread.run(Thread.java:834)
InterProScan analysis failed. Exception thrown by StandaloneBlackBoxMaster. Check the log file for details
```

This is not super informative...not sure what to do now. I'll look into these errors more later. 
