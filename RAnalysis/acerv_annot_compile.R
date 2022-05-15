# Title: Acerv Annotation Compilation
# Author: Jill Ashey
# Date: 12/07/2020

  

## Species: Acropora Cervicornis

#This script takes the results of functional annotation services and combines the results. Nucleotide CDS sequences were annotated using DIAMONDSEARCH BLASTX, resulting in 29515 hits. These hits were used as input into:
#  - Uniprot
#  - Blast2GO

#Additional annotation was provided by
#  - InterProScan





## Load libraries
library(tidyverse)
library(dplyr)



## Diamond BLAST

# Load DIAMOND BLAST results
acerv_blast <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Diamond/Acerv_annot.tab", col_names = FALSE)
colnames(acerv_blast) <- c("seqName", "top_hit", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",  "qlen", "slen")
head(acerv_blast)
dim(acerv_blast) #29515 x 14



## Uniprot

# Uniprot mapping occurred on 11/10/20.
#EMBL/GenBank/DDBJ CDS protein IDs generated from Diamond BLAST were mapped to UniProtKB database IDs.

#Because there were many Diamond BLAST hits for Acerv, I had to break up the tab file into chunks that Uniprot could handle. So there are 15 acerv_Uniprot files to be read and then rbind() together
u1 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/acerv/acerv_Uniprot_1.tab", col_names = TRUE)
colnames(u1) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u1)
dim(u1) # 218 x 12

u2 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/acerv/acerv_Uniprot_2.tab", col_names = TRUE)
colnames(u2) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u2)
dim(u2) # 230 x 12

u3 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/acerv/acerv_Uniprot_3.tab", col_names = TRUE)
colnames(u3) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u3)
dim(u3) # 216 x 12

u4 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/acerv/acerv_Uniprot_4.tab", col_names = TRUE)
colnames(u4) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u4)
dim(u4) # 238 x 12

u5 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/acerv/acerv_Uniprot_5.tab", col_names = TRUE)
colnames(u5) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u5)
dim(u5) # 218 x 12

u6 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/acerv/acerv_Uniprot_6.tab", col_names = TRUE)
colnames(u6) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u6)
dim(u6) # 238 x 12

u7 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/acerv/acerv_Uniprot_7.tab", col_names = TRUE)
colnames(u7) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u7)
dim(u7) # 232 x 12

u8 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/acerv/acerv_Uniprot_8.tab", col_names = TRUE)
colnames(u8) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u8)
dim(u8) # 219 x 12

u9 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/acerv/acerv_Uniprot_9.tab", col_names = TRUE)
colnames(u9) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u9)
dim(u9) # 106 x 12

u10 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/acerv/acerv_Uniprot_10.tab", col_names = TRUE)
colnames(u10) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u10)
dim(u10) # 86 x 12

u11 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/acerv/acerv_Uniprot_11.tab", col_names = TRUE)
colnames(u11) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u11)
dim(u11) # 107 x 12

u12 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/acerv/acerv_Uniprot_12.tab", col_names = TRUE)
colnames(u12) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u12)
dim(u12) # 89 x 12

u13 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/acerv/acerv_Uniprot_13.tab", col_names = TRUE)
colnames(u13) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u13)
dim(u13) # 117 x 12

u14 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/acerv/acerv_Uniprot_14.tab", col_names = TRUE)
colnames(u14) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u14)
dim(u14) # 95 x 12

u15 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/acerv/acerv_Uniprot_15.tab", col_names = TRUE)
colnames(u15) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u15)
dim(u15) # 70 x 12


#Compile the Uniprot files 
uniprot_results <- bind_rows(u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15)
uniprot_results <- unique(uniprot_results)
head(uniprot_results)
dim(uniprot_results) # 1844 x 12
uniprot_results <- filter(uniprot_results, grepl("GO",go_ids)) #Select only gnes with GO terms
dim(uniprot_results) # 1095 x 12



## Blast2GO

#Nucleotide CDS sequences were annotated using DIAMONDSEARCH BLASTX, resulting in 29515. These hits were used as input into Blast2GO to obtain GO terms using the 12/05/2020 obo database.

B2G_results <- read.csv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Blast2GO/acerv_blast2go_table.csv")
B2G_results <- select(B2G_results, c("SeqName", "Description", "Length", "e.Value", "sim.mean", "GO.IDs", "GO.Names"))
colnames(B2G_results) <- c("seqName", "top_hit", "length", "eValue", "simMean", "GO.ID", "GO_names")
head(B2G_results)
dim(B2G_results) # 29515 x 7
B2G_results <- filter(B2G_results, grepl("GO",GO.ID)) #Genes with GO terms - 136
dim(B2G_results) # 2092 x 7
























# Generate a list of GO terms - UniProt
uniprot_GO <- select(uniprot_results, my_list, go_ids)
splitted <- strsplit(as.character(uniprot_GO$go_ids), ";") #split into multiple GO ids
uniprot_GO <- data.frame(v1 = rep.int(uniprot_GO$my_list, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their GO terms in a single row
uniprot_GO <- unique(uniprot_GO)
colnames(uniprot_GO) <- c("gene_id", "GO.ID")
uniprot_GO$GO.ID <- gsub(" ", "", uniprot_GO$GO.ID)
nrow(uniprot_GO) # 2875 total GO terms from Uniprot
length(unique(uniprot_GO$GO.ID)) # 824 unique GO terms from Uniprot
length(unique(uniprot_GO$gene_id)) # 1095 unique gene ids from Uniprot 
# Not technically gene ids. Uniprot has no gene id info--actually ids from Uniprot. When I combine the Uniport and B2G files, the uniprot ids will then be associated with gene ids 







## Blast2GO

#Nucleotide CDS sequences were annotated using DIAMONDSEARCH BLASTX, resulting in 29515. These hits were used as input into Blast2GO to obtain GO terms using the 12/05/2020 obo database.

B2G_results <- read.csv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Blast2GO/acerv_blast2go_table.csv")
B2G_results <- select(B2G_results, c("SeqName", "Description", "Length", "e.Value", "sim.mean", "GO.IDs", "GO.Names"))
colnames(B2G_results) <- c("seqName", "top_hit", "length", "eValue", "simMean", "GO.ID", "GO_names")
head(B2G_results)
dim(B2G_results) # 29515 x 7
B2G_results <- filter(B2G_results, grepl("GO",GO.ID)) #Genes with GO terms - 136
dim(B2G_results) # 2092 x 7

# Generate a list of GO terms - B2G
B2G_results_GO <- select(B2G_results, top_hit, GO.ID)
splitted <- strsplit(as.character(B2G_results_GO$GO.ID), ";") #split into multiple GO ids
B2G_results_GO <- data.frame(v1 = rep.int(B2G_results_GO$top_hit, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their GO terms in a single row
B2G_results_GO <- unique(B2G_results_GO)
colnames(B2G_results_GO) <- c("gene_id", "GO.ID")
B2G_results_GO$GO.ID <- gsub("F:", "", B2G_results_GO$GO.ID)
B2G_results_GO$GO.ID <- gsub("C:", "", B2G_results_GO$GO.ID)
B2G_results_GO$GO.ID <- gsub("P:", "", B2G_results_GO$GO.ID)
B2G_results_GO$GO.ID <- gsub(" ", "", B2G_results_GO$GO.ID)
head(B2G_results_GO)
nrow(B2G_results_GO) # 4822 total GO terms from B2G
length(unique(B2G_results_GO$GO.ID)) # 1130 unique GO terms from B2G
length(unique(B2G_results_GO$gene_id)) # 1458 unique gene ids from B2G



# Find intersections and unique results for each method (Uniprot and B2G)

## GO
# Intersection between GO terms for B2G and Uniprot
BU_GO <- intersect(B2G_results_GO$GO.ID, uniprot_GO$GO.ID)
length(unique(BU_GO)) # 626 similar GO terms between B2G and Uniprot

# Difference in GO terms for B2G and Uniprot - B2G
BUunique_GO <- setdiff(B2G_results_GO$GO.ID, uniprot_GO$GO.ID)
length(unique(BUunique_GO)) # 504 GO terms unique to B2G

# Difference in GO terms for B2G and Uniprot - Uniprot
UBunique <- setdiff(uniprot_GO$GO.ID, B2G_results_GO$GO.ID)
length(unique(UBunique)) # 198 GO terms unique to Uniprot

## gene id
# Intersection between gene id for B2G and Uniprot
BU_gene <- intersect(B2G_results_GO$gene_id, uniprot_GO$gene_id)
length(unique(BU_gene)) # 1042 similar gene ids between B2G and Uniprot

# Difference in gene ids for B2G and Uniprot
BUunique_gene <- setdiff(B2G_results_GO$gene_id, uniprot_GO$gene_id)
length(unique(BUunique_gene)) # 416 gene ids unique to B2G

# Difference in gene ids for B2G and Uniprot
UBunique_gene <- setdiff(uniprot_GO$gene_id, B2G_results_GO$gene_id)
length(unique(UBunique_gene)) # 53 gene ids unique to uniprot



## Merge Annotations - uniprot + b2g
acerv_annot <- left_join(acerv_blast, B2G_results, by="seqName")
acerv_annot <- select(acerv_annot, seqName, top_hit.x, length.x, evalue, bitscore, simMean, GO.ID, GO_names)
colnames(acerv_annot) <- c("gene_id", "top_hit", "length", "evalue", "bitscore", "simMean", "GO.ID", "GO_names")
uniprot_results_GO <- select(uniprot_results_GO, -top_hit)
uniprot_results_GO <- rename(uniprot_results_GO, "top_hit"="my_list")
acerv_annot <- merge(acerv_annot, uniprot_results_GO, by="top_hit", all.x = T)
acerv_annot$GO.IDs <- paste(acerv_annot$GO.ID, acerv_annot$go_ids, sep=';') #generate new column with concatenated GO IDs
acerv_annot$GO_terms <- paste(acerv_annot$GO_names, acerv_annot$gene_ontology, sep=';') #generate new column with concatenated GO IDs
acerv_annot <- select(acerv_annot, c("top_hit", "gene_id", "length.x", "evalue", "bitscore", "simMean", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "ko", "kegg", "GO.IDs", "GO_terms"))
acerv_annot <- rename(acerv_annot, "GO.ID"="GO.IDs")
names(acerv_annot)
head(acerv_annot)
tail(acerv_annot)
dim(acerv_annot) # 29515 x 16
write.csv(acerv_annot, "~/Desktop/acerv_FuncAnn_UniP_B2G.csv")








# IPS
IPS_GO <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOseq/acerv_GOterms.csv", header=TRUE)
IPS_GO <- select(IPS_GO, -X)
colnames(IPS_GO)[1] <-"gene_id"
#IPS_GO$gene_id <- gsub(".m1", "", IPS_GO$gene_id)
#IPS_GO$gene_id <- gsub("model", "TU", IPS_GO$gene_id)
length(unique(IPS_GO$gene_id)) # 12828
IPS_GO <- select(IPS_GO, c("gene_id", "GO_term")) # taking out source and score tho I want to leave them in--just for comparison though
splitted <- strsplit(as.character(IPS_GO$GO), ",") #split into multiple GO ids
IPS_GO <- data.frame(v1 = rep.int(IPS_GO$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their 
colnames(IPS_GO) <- c("gene_id", "GO.ID")
IPS_GO <- unique(IPS_GO)
head(IPS_GO)
nrow(IPS_GO) # 30085 total GO terms from IPS
length(unique(IPS_GO$GO.ID)) # 1945 unique GO terms from IPS
length(unique(IPS_GO$gene_id)) # 12828 unique gene ids from IPS

# From B2G and Uniprot
FuncAnn_UniP_B2G <- read.csv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/final_Annotations/acerv_FuncAnn_UniP_B2G.csv")
length(unique(FuncAnn_UniP_B2G$gene_id)) # 29515
FuncAnn_UniP_B2G_GO <- filter(FuncAnn_UniP_B2G, grepl("GO",GO.ID)) #Select only gnes with GO terms
length(unique(FuncAnn_UniP_B2G_GO$gene_id)) # 2172
FuncAnn_UniP_B2G_GO <- select(FuncAnn_UniP_B2G_GO, c("gene_id", "GO.ID")) # taking out source and score tho I want to leave them in--just for comparison though
splitted <- strsplit(as.character(FuncAnn_UniP_B2G_GO$GO.ID), ";") #split into multiple GO ids
FuncAnn_UniP_B2G_GO <- data.frame(v1 = rep.int(FuncAnn_UniP_B2G_GO$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their 
colnames(FuncAnn_UniP_B2G_GO) <- c("gene_id", "GO.ID")
FuncAnn_UniP_B2G_GO$GO.ID <- gsub("F:", "", FuncAnn_UniP_B2G_GO$GO.ID)
FuncAnn_UniP_B2G_GO$GO.ID <- gsub("C:", "", FuncAnn_UniP_B2G_GO$GO.ID)
FuncAnn_UniP_B2G_GO$GO.ID <- gsub("P:", "", FuncAnn_UniP_B2G_GO$GO.ID)
FuncAnn_UniP_B2G_GO$GO.ID <- gsub(" ", "", FuncAnn_UniP_B2G_GO$GO.ID)
FuncAnn_UniP_B2G_GO <- unique(FuncAnn_UniP_B2G_GO)
head(FuncAnn_UniP_B2G_GO)
nrow(FuncAnn_UniP_B2G_GO) # 8390 total GO terms from B2G + Uniprot
length(unique(FuncAnn_UniP_B2G_GO$GO.ID)) # 1329 unique GO terms from B2G + Uniprot
length(unique(FuncAnn_UniP_B2G_GO$gene_id)) # 2172 unique gene ids from B2G + Uniprot

# Find intersections and unique results for each method (B2G + Uniprot and IPS)
## GO
# Intersection between GO terms for B2G + Uniprot and IPS
IF_GO <- intersect(IPS_GO$GO.ID, FuncAnn_UniP_B2G_GO$GO.ID)
length(unique(IF_GO)) # 747 similar GO terms between B2G + Uniprot and IPS

# Difference in GO terms for B2G + Uniprot and IPS - FuncAnn first
FIunique_GO <- setdiff(FuncAnn_UniP_B2G_GO$GO.ID, IPS_GO$GO.ID)
length(unique(FIunique_GO)) # 582 GO terms unique to B2G + Uniprot

# Difference in GO terms for B2G + Uniprot and IPS - IPS
IFunique_GO <- setdiff(IPS_GO$GO.ID, FuncAnn_UniP_B2G_GO$GO.ID)
length(unique(IFunique_GO)) # 1198 GO terms unique to Uniprot

## gene id
# Intersection between GO terms for B2G + Uniprot and IPS
IF_gene<- intersect(IPS_GO$gene_id, FuncAnn_UniP_B2G_GO$gene_id)
length(unique(IF_gene)) # 1140 similar gene ids between B2G + Uniprot and IPS

# Difference in gene ids for B2G + Uniprot and IPS - B2G + Uniprot
FIunique_gene <- setdiff(FuncAnn_UniP_B2G_GO$gene_id, IPS_GO$gene_id)
length(unique(FIunique_gene)) # 1032 gene ids unique to B2G + Uniprot

# Difference in gene ids for B2G + Uniprot and IPS - IPS
IFunique_gene <- setdiff(IPS_GO$gene_id, FuncAnn_UniP_B2G_GO$gene_id)
length(unique(IFunique_gene)) # 11688 gene ids unique to IPS


# Aggregate results 
agg <- aggregate(IPS_GO$GO.ID, list(IPS_GO$gene_id), paste, collapse = ",")
colnames(agg) <- c("gene_id", "GO.ID")

full_annot <- merge(agg, FuncAnn_UniP_B2G, by = "gene_id", all.x = T)
full_annot$GO.IDs <- paste(full_annot$GO.ID.x, full_annot$GO.ID.y, sep=';') #generate new column with concatenated GO IDs
full_annot <- select(full_annot, -c(X, GO.ID.x, GO.ID.y))
write.csv(full_annot, file = "~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/final_Annotations/acerv_fullAnnot.csv")




























## InterProScan
# IPS <- read.csv("~/Desktop/PutnamLab/FunctionalAnnotation/InterProScan/acerv.interpro.gff3", header = F, sep = "\t")
# colnames(IPS) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
# IPS <- na.omit(IPS)
# IPS <- IPS %>% filter(type=="protein_match")
# IPS_GO <- filter(IPS, grepl("GO",attributes)) #Genes with GO terms
# IPS_GO$GO <- regmatches(IPS_GO$attributes, gregexpr("(?<=Ontology_term=).*", IPS_GO$attributes, perl = TRUE)) #removing everything in attributes up to Ontology_term=
# IPS_GO$GO <- gsub(";.*", "", IPS_GO$GO)
# 
# # Generate a list of GO terms
# IPS_GO <- select(IPS_GO, c("seqid", "GO")) # taking out source and score tho I want to leave them in--just for comparison though
# splitted <- strsplit(as.character(IPS_GO$GO), ",") #split into multiple GO ids
# IPS_GO <- data.frame(v1 = rep.int(IPS_GO$seqid, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their 
# colnames(IPS_GO) <- c("gene_id", "GO.ID")
# IPS_GO <- unique(IPS_GO)
# head(IPS_GO)
# nrow(IPS_GO) # 30096 total GO terms from IPS
# length(unique(IPS_GO$GO)) # 1946 unique GO terms from IPS
# length(unique(IPS_GO$gene_id)) # 12837 unique gene ids from IPS
# 
# 
# # Find intersections and unique results for each method (IPS and B2G)
# # Cannot compare gene ids because B2G_results_GO has different 'gene ids' because I needed it to pair with Uniprot
# 
# ## GO
# # Intersection between GO terms for B2G and IPS
# BI_GO <- intersect(B2G_results_GO$GO.ID, IPS_GO$GO.ID)
# length(unique(BI_GO)) # 660 similar GO terms between B2G and IPS
# 
# # Difference in GO terms for B2G and IPS - B2G
# BIunique_GO <- setdiff(B2G_results_GO$GO.ID, IPS_GO$GO.ID)
# length(unique(BIunique_GO)) # 470 GO terms unique to B2G
# 
# # Difference in GO terms for B2G and IPS - IPS
# IBunique <- setdiff(IPS_GO$GO.ID, B2G_results_GO$GO.ID)
# length(unique(IBunique)) # 1286 GO terms unique to IPS
# # Cannot compare gene ids because B2G_results_GO has different 'gene ids' because I needed it to pair with Uniprot
# 
# # Find intersections and unique results for each method (IPS and Uniprot)
# # Cannot compare gene ids because Uniprot has different 'gene ids' because I needed it to pair with Uniprot
# ## GO
# # Intersection between GO terms for Uniprot and IPS
# UI_GO <- intersect(uniprot_GO$GO.ID, IPS_GO$GO.ID)
# length(unique(UI_GO)) # 548 similar GO terms between Uniprot and IPS
# 
# # Difference in GO terms for Uniprot and IPS - Uniprot
# UIunique_GO <- setdiff(uniprot_GO$GO.ID, IPS_GO$GO.ID)
# length(unique(UIunique_GO)) # 276 GO terms unique to Uniprot
# 
# # Difference in GO terms for B2G and IPS - IPS
# IUunique_GO <- setdiff(IPS_GO$GO.ID, uniprot_GO$GO.ID)
# length(unique(IUunique_GO)) # 1398 GO terms unique to IPS
# 
# # Find intersections and unique results for each method (IPS and acerv_annot)
# # Cannot compare gene ids because Uniprot has different 'gene ids' because I needed it to pair with Uniprot
# 
# 
# 
# ## GO
# # Before comparison, generate a list of GO terms - acerv_annot
# acerv_annot_GO <- select(acerv_annot, gene_id, GO.ID)
# splitted <- strsplit(as.character(acerv_annot_GO$GO.ID), ";") #split into multiple GO ids
# acerv_annot_GO <- data.frame(v1 = rep.int(acerv_annot_GO$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their GO terms in a single row
# colnames(acerv_annot_GO) <- c("gene_id", "GO.ID")
# acerv_annot_GO <- unique(acerv_annot_GO)
# acerv_annot_GO$GO.ID <- gsub("F:", "", acerv_annot_GO$GO.ID)
# acerv_annot_GO$GO.ID <- gsub("C:", "", acerv_annot_GO$GO.ID)
# acerv_annot_GO$GO.ID <- gsub("P:", "", acerv_annot_GO$GO.ID)
# acerv_annot_GO$GO.ID <- gsub(" ", "", acerv_annot_GO$GO.ID)
# acerv_annot_GO <- filter(acerv_annot_GO, grepl("GO",GO.ID)) #Genes with GO terms - 136
# acerv_annot_GO <- unique(acerv_annot_GO)
# colnames(acerv_annot_GO) <- c("gene_id", "GO.ID")
# 
# # Intersection between GO terms for acerv_annot (b2g + uniprot) and IPS
# AI_GO <- intersect(acerv_annot_GO$GO.ID, IPS_GO$GO.ID)
# length(unique(AI_GO)) # 747 similar GO terms between acerv_annot (b2g + uniprot) and IPS
# 
# # Difference in GO terms for acerv_annot (b2g + uniprot) and IPS - acerv_annot (b2g + uniprot)
# AIunique_GO <- setdiff(acerv_annot_GO$GO.ID, IPS_GO$GO.ID)
# length(unique(AIunique_GO)) # 581 GO terms unique to Uniprot
# 
# # Difference in GO terms for acerv_annot (b2g + uniprot) and IPS - IPS
# IAunique_GO <- setdiff(IPS_GO$GO.ID, acerv_annot_GO$GO.ID)
# length(unique(IAunique_GO)) # 1199 GO terms unique to IPS
# 
# 
# ## gene id
# # Intersection between gene_id for acerv_annot (b2g + uniprot) and IPS
# AI_gene <- intersect(acerv_annot_GO$gene_id, IPS_GO$gene_id)
# length(unique(AI_gene)) # 1140 similar gene ids between acerv_annot (b2g + uniprot) and IPS
# 
# # Difference in gene_id terms for acerv_annot (b2g + uniprot) and IPS - acerv_annot (b2g + uniprot)
# AIunique_gene <- setdiff(acerv_annot_GO$gene_id, IPS_GO$gene_id)
# length(unique(AIunique_gene)) # 1032 gene ids terms unique to Uniprot
# 
# # Difference in gene ids for acerv_annot (b2g + uniprot) and IPS - IPS
# IAunique_gene <- setdiff(IPS_GO$gene_id, acerv_annot_GO$gene_id)
# length(unique(IAunique_gene)) # 11697 GO terms unique to IPS
# 
# 
# 
# 
# 
# ## Merge Annotations (again) - merge IPS and acerv_annot
# 
# # Merge IPS gene names + GO terms (and any other info I can bring from original IPS file) with acerv_annot gene names + GO terms (and any other info I can bring along )
# IPS_results <- read.csv("~/Desktop/PutnamLab/FunctionalAnnotation/InterProScan/acerv.interpro.gff3", header = F, sep = "\t")
# colnames(IPS_results) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
# IPS_results <- na.omit(IPS_results)
# IPS_results <- IPS_results %>% filter(type=="protein_match")
# IPS_results_GO <- filter(IPS_results, grepl("GO",attributes)) #Genes with GO terms
# IPS_results_GO$GO <- regmatches(IPS_results_GO$attributes, gregexpr("(?<=Ontology_term=).*", IPS_results_GO$attributes, perl = TRUE)) #removing everything in attributes up to Ontology_term=
# IPS_results_GO$GO <- gsub(";.*", "", IPS_results_GO$GO)
# IPS_results_GO$GO_name <- regmatches(IPS_results_GO$attributes, gregexpr("(?<=signature_desc=).*", IPS_results_GO$attributes, perl = TRUE)) #removing everything in attributes up to signature_desc=
# IPS_results_GO$GO_name <- gsub(";.*", "", IPS_results_GO$GO_name)
# IPS_results_GO$Dbxref <- regmatches(IPS_results_GO$attributes, gregexpr("(?<=Dbxref=).*", IPS_results_GO$attributes, perl = TRUE)) #removing everything in attributes up to Dbxref=
# IPS_results_GO$Dbxref <- gsub(";.*", "", IPS_results_GO$Dbxref)
# length(unique(IPS_results_GO$seqid))
# 
# # selecting only certain cols 
# IPS_results_GO <- select(IPS_results_GO, c("seqid", "GO", "GO_name", "Dbxref"))
# IPS_results_GO <- unique(IPS_results_GO)
# colnames(IPS_results_GO) <- c("gene_id", "GO", "GO_name", "Dbxref")
# 
# # aggregate IPS cols before merging w/ acerv_annot
# agg1 <- aggregate(IPS_results_GO$GO, list(IPS_results_GO$gene_id), paste, collapse = ",")
# colnames(agg1) <- c("gene", "GO.ID")
# agg2 <- aggregate(IPS_results_GO$GO_name, list(IPS_results_GO$gene_id), paste, collapse = ",")
# colnames(agg2) <- c("gene", "GO_name")
# agg3 <- aggregate(IPS_results_GO$Dbxref, list(IPS_results_GO$gene_id), paste, collapse = ",")
# colnames(agg3) <- c("gene", "Dbxref")
# 
# # merge the aggregates of IPS
# merge_all_IPS <- merge(agg1, agg2, all.x = T, by = "gene")
# merge_all_IPS <- merge(merge_all_IPS, agg3, all.x = T, by = "gene")
# colnames(merge_all_IPS) <- c("gene_id", "GO.ID_IPS", "GO.name_IPS", "Dbxref_IPS")
# nrow(merge_all_IPS) # 12837
# 
# # merge all info 
# acerv_annot_bu <- select(acerv_annot, c("top_hit", "gene_id", "uniprotkb_entry", "protein_names", "gene_names", "organism", "ko", "kegg", "GO.ID", "GO_terms"))
# nrow(acerv_annot_bu) # 29215
# full_annot <- merge(acerv_annot_bu, merge_all_IPS, by = "gene_id", all.x = T)
# write.csv(full_annot, "~/Desktop/acerv_FuncAnn_UniP_B2G_IPS.csv")
# 
# 
# 
# 
# 
# ## Find metrics for new annotation - b2g + uniprot (acerv_annot)
# new_Sig_Alingments=nrow(acerv_annot)
# new_Genes_with_GO <- nrow(filter(acerv_annot, grepl("GO",GO.ID))) #Genes with GO terms...2172
# new_Genes_with_Kegg <- nrow(filter(acerv_annot, grepl("K",ko))) #Genes with Kegg terms...36L
# new.avg.Eval <- mean(acerv_annot$evalue) # 3.410318e-08
# new.median.Eval <- median(acerv_annot$evalue) # 6.5e-91
# new.avg.bit <- mean(acerv_annot$bitscore) # 525.0099
# new.median.bit <- median(acerv_annot$bitscore) # 344.4
# # total GO terms
# new_GO <- select(acerv_annot, gene_id, GO.ID)
# splitted <- strsplit(as.character(new_GO$GO.ID), ";") #split into multiple GO ids
# new_GO <- data.frame(v1 = rep.int(new_GO$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their GO terms in a single row
# colnames(new_GO) <- c("gene_id", "GO.ID")
# new_totGO_narm <- filter(new_GO, GO.ID!="NA")
# new_totGO <- nrow(new_totGO_narm) # 10924
# write.csv(new_totGO_narm, "~/Desktop/acerv_new_totGO_narm_UniP_B2G.csv")
# # total unique GO terms
# new_uniqueGO <- unique(new_totGO_narm$GO.ID)
# new_uniqueGO <- length(new_uniqueGO) # 2226
# 
# # total Kegg terms
# new_Kegg <- select(acerv_annot, gene_id, ko)
# splitted <- strsplit(as.character(new_Kegg$ko), ";") #split into multiple Kegg ids
# new.Keggterms <- data.frame(v1 = rep.int(new_Kegg$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their Kegg terms in a single row
# colnames(new.Keggterms) <- c("gene_id", "Kegg.ID")
# new.Keggterms$Kegg.ID <- replace_na(new.Keggterms$Kegg.ID, "NA")
# new_totKegg_narm <- filter(new.Keggterms, Kegg.ID!="NA")
# new_totKegg <- nrow(new_totKegg_narm) # 36
# write.csv(new_totKegg_narm, "~/Desktop/acerv_new_totKegg_narm_UniP_B2G.csv")
# new_uniqueKegg <- unique(new_totKegg_narm$Kegg.ID)
# new_uniqueKegg <- length(new_uniqueKegg) # 32
# 
# 
# 
# ## Find metrics for new annotation - b2g + uniprot (acerv_annot)
# new_Sig_Alingments=nrow(full_annot)
# new_Genes_with_GO <- nrow(filter(acerv_annot, grepl("GO",GO.ID))) #Genes with GO terms...2172
# new_Genes_with_Kegg <- nrow(filter(acerv_annot, grepl("K",ko))) #Genes with Kegg terms...36L
# new.avg.Eval <- mean(acerv_annot$evalue) # 3.410318e-08
# new.median.Eval <- median(acerv_annot$evalue) # 6.5e-91
# new.avg.bit <- mean(acerv_annot$bitscore) # 525.0099
# new.median.bit <- median(acerv_annot$bitscore) # 344.4
# # total GO terms
# new_GO <- select(acerv_annot, gene_id, GO.ID)
# splitted <- strsplit(as.character(new_GO$GO.ID), ";") #split into multiple GO ids
# new_GO <- data.frame(v1 = rep.int(new_GO$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their GO terms in a single row
# colnames(new_GO) <- c("gene_id", "GO.ID")
# new_totGO_narm <- filter(new_GO, GO.ID!="NA")
# new_totGO <- nrow(new_totGO_narm) # 10924
# write.csv(new_totGO_narm, "~/Desktop/acerv_new_totGO_narm_UniP_B2G.csv")
# # total unique GO terms
# new_uniqueGO <- unique(new_totGO_narm$GO.ID)
# new_uniqueGO <- length(new_uniqueGO) # 2226
# 
# # total Kegg terms
# new_Kegg <- select(acerv_annot, gene_id, ko)
# splitted <- strsplit(as.character(new_Kegg$ko), ";") #split into multiple Kegg ids
# new.Keggterms <- data.frame(v1 = rep.int(new_Kegg$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their Kegg terms in a single row
# colnames(new.Keggterms) <- c("gene_id", "Kegg.ID")
# new.Keggterms$Kegg.ID <- replace_na(new.Keggterms$Kegg.ID, "NA")
# new_totKegg_narm <- filter(new.Keggterms, Kegg.ID!="NA")
# new_totKegg <- nrow(new_totKegg_narm) # 36
# write.csv(new_totKegg_narm, "~/Desktop/acerv_new_totKegg_narm_UniP_B2G.csv")
# new_uniqueKegg <- unique(new_totKegg_narm$Kegg.ID)
# new_uniqueKegg <- length(new_uniqueKegg) # 32








