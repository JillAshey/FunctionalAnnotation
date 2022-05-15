# Title: Mcav Annotation Compilation
# Author: Jill Ashey
# Date: 12/09/2020



## Species: M cavernosa

#This script takes the results of functional annotation services and combines the results. Nucleotide CDS sequences were annotated using DIAMONDSEARCH BLASTX, resulting in 3191 hits. These hits were used as input into:
#  - Uniprot
#  - Blast2GO

#Additional annotation was provided by
#  - InterProScan





## Load libraries
library(tidyverse)
library(dplyr)



## Diamond BLAST

# Load DIAMOND BLAST results
mcav_blast <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Diamond/mcav/Mcav_annot.tab", col_names = FALSE)
colnames(mcav_blast) <- c("seqName", "top_hit", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",  "qlen", "slen")
head(mcav_blast)
dim(mcav_blast) #3191 x 14



## Uniprot

# Uniprot mapping occurred on 11/10/20.
#EMBL/GenBank/DDBJ CDS protein IDs generated from Diamond BLAST were mapped to UniProtKB database IDs.

#Because there were many Diamond BLAST hits for Mcav, I had to break up the tab file into chunks that Uniprot could handle. So there are 2 mcav_Uniprot files to be read and then rbind() together
u1 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/mcav/mcav_Uniprot_1.tab", col_names = TRUE)
colnames(u1) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u1)
dim(u1) # 119 x 12

u2 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/mcav/mcav_Uniprot_2.tab", col_names = TRUE)
colnames(u2) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u2)
dim(u2) # 81 x 12

#Compile the Uniprot files 
uniprot_results <- bind_rows(u1,u2)
uniprot_results <- unique(uniprot_results)
head(uniprot_results)
dim(uniprot_results) # 199 x 12
uniprot_results <- filter(uniprot_results, grepl("GO",go_ids)) #Select only gnes with GO terms
dim(uniprot_results) # 122 x 12

# Generate a list of GO terms - UniProt
uniprot_GO <- select(uniprot_results, my_list, go_ids)
splitted <- strsplit(as.character(uniprot_GO$go_ids), ";") #split into multiple GO ids
uniprot_GO <- data.frame(v1 = rep.int(uniprot_GO$my_list, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their GO terms in a single row
uniprot_GO <- unique(uniprot_GO)
colnames(uniprot_GO) <- c("gene_id", "GO.ID")
uniprot_GO$GO.ID <- gsub(" ", "", uniprot_GO$GO.ID)
nrow(uniprot_GO) # 306 total GO terms from Uniprot
length(unique(uniprot_GO$GO.ID)) # 165 unique GO terms from Uniprot
length(unique(uniprot_GO$gene_id)) # 122 unique gene ids from Uniprot 
# Not technically gene ids. Uniprot has no gene id info--actually ids from Uniprot. When I combine the Uniport and B2G files, the uniprot ids will then be associated with gene ids 




## Blast2GO

#Nucleotide CDS sequences were annotated using DIAMONDSEARCH BLASTX, resulting in 3191 These hits were used as input into Blast2GO to obtain GO terms using the 12/05/2020 obo database.

B2G_results <- read.csv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Blast2GO/mcav_blast2go_table.csv")
B2G_results <- select(B2G_results, c("SeqName", "Description", "Length", "e.Value", "sim.mean", "GO.IDs", "GO.Names"))
colnames(B2G_results) <- c("seqName", "top_hit", "length", "eValue", "simMean", "GO.ID", "GO_names")
head(B2G_results)
dim(B2G_results) # 3191 x 7
B2G_results <- filter(B2G_results, grepl("GO",GO.ID)) #Genes with GO terms - 136
dim(B2G_results) # 136 x 7

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
nrow(B2G_results_GO) # 380 total GO terms from B2G
length(unique(B2G_results_GO$GO.ID)) # 180 unique GO terms from B2G
length(unique(B2G_results_GO$gene_id)) # 117 unique gene ids from B2G

# Find intersections and unique results for each method (Uniprot and B2G)
## GO
# Intersection between GO terms for B2G and Uniprot
BU_GO <- intersect(B2G_results_GO$GO.ID, uniprot_GO$GO.ID)
length(unique(BU_GO)) # 134 similar GO terms between B2G and Uniprot

# Difference in GO terms for B2G and Uniprot - B2G
BUunique_GO <- setdiff(B2G_results_GO$GO.ID, uniprot_GO$GO.ID)
length(unique(BUunique_GO)) # 46 GO terms unique to B2G

# Difference in GO terms for B2G and Uniprot - Uniprot
UBunique <- setdiff(uniprot_GO$GO.ID, B2G_results_GO$GO.ID)
length(unique(UBunique)) # 31 GO terms unique to Uniprot



## Merge Annotations - uniprot + b2g
mcav_annot <- left_join(mcav_blast, B2G_results, by="seqName")
mcav_annot <- select(mcav_annot, seqName, top_hit.x, length.x, evalue, bitscore, simMean, GO.ID, GO_names)
colnames(mcav_annot) <- c("gene_id", "top_hit", "length", "evalue", "bitscore", "simMean", "GO.ID", "GO_names")
uniprot_results <- select(uniprot_results, -top_hit)
uniprot_results <- rename(uniprot_results, "top_hit"="my_list")
mcav_annot <- merge(mcav_annot, uniprot_results, by="top_hit", all.x = T)
mcav_annot$GO.IDs <- paste(mcav_annot$GO.ID, mcav_annot$go_ids, sep=';') #generate new column with concatenated GO IDs
mcav_annot$GO_terms <- paste(mcav_annot$GO_names, mcav_annot$gene_ontology, sep=';') #generate new column with concatenated GO IDs
mcav_annot <- select(mcav_annot, c("top_hit", "gene_id", "length.x", "evalue", "bitscore", "simMean", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "ko", "kegg", "GO.IDs", "GO_terms"))
mcav_annot <- rename(mcav_annot, "GO.ID"="GO.IDs")
names(mcav_annot)
head(mcav_annot)
tail(mcav_annot)
dim(mcav_annot) # 3191 x 16
#write.csv(mcav_annot, "~/Desktop/mcav_FuncAnn_UniP_B2G.csv")








# IPS
mcav_IPS <- read.csv("~/Desktop/PutnamLab/FunctionalAnnotation/InterProScan/mcav.interpro.gff3", header = FALSE, sep="\t", skip=4)
colnames(mcav_IPS) <- c("gene_id", "Predict", "id", "start","stop", "pos1", "pos2","pos3", "attr")
mcav_IPS <- select(mcav_IPS, -Predict)
mcav_IPS <- unique(mcav_IPS)
mcav_IPS_GO <- filter(mcav_IPS, grepl("GO:", attr)) # select only rows with GO terms
mcav_IPS_GO$GO.ID <- regmatches(mcav_IPS_GO$attr, gregexpr("(?<=Ontology_term=).*", mcav_IPS_GO$attr, perl = TRUE))
mcav_IPS_GO$GO.ID <- gsub(";.*", "", mcav_IPS_GO$GO.ID)
mcav_IPS_GO <- mcav_IPS_GO %>% mutate(mcav_IPS_GO, length = stop - start)
# Separate out GO terms
mcav_IPS_GO_sub <- select(mcav_IPS_GO, c(gene_id, GO.ID))
splitted <- strsplit(as.character(mcav_IPS_GO_sub$GO.ID), ",") #split into multiple GO ids
mcav_IPS_GO_sub <- data.frame(v1 = rep.int(mcav_IPS_GO_sub$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their 
mcav_IPS_GO_sub <- unique(mcav_IPS_GO_sub)
colnames(mcav_IPS_GO_sub) <- c("gene_id", "GO.ID")
head(mcav_IPS_GO_sub)
length(unique(mcav_IPS_GO_sub$GO.ID)) # 1720 unique GO terms from IPS
length(unique(mcav_IPS_GO_sub$gene_id)) # 10791 unique gene ids from IPS

# Find intersections and unique results for each method (B@G, Uniprot and IPS)
## GO
# Intersection between GO terms for B2G and IPS
BI_GO <- intersect(B2G_results_GO$GO.ID, mcav_IPS_GO_sub$GO.ID)
length(unique(BI_GO)) # 143 similar GO terms between B2G and IPS

# Difference in GO terms for B2G and IPS - B2G
BIunique_GO <- setdiff(B2G_results_GO$GO.ID, mcav_IPS_GO_sub$GO.ID)
length(unique(BIunique_GO)) # 37 GO terms unique to B2G

# Difference in GO terms for B2G and IPS - IPS
IBunique <- setdiff(mcav_IPS_GO_sub$GO.ID, B2G_results_GO$GO.ID)
length(unique(IBunique)) # 1577 GO terms unique to IPS

# Intersection between GO terms for Uniprot and IPS
UI_GO <- intersect(uniprot_GO$GO.ID, mcav_IPS_GO_sub$GO.ID)
length(unique(UI_GO)) # 138 similar GO terms between Uniprot and IPS

# Difference in GO terms for Uniprot and IPS - Uniprot
UIunique_GO <- setdiff(uniprot_GO$GO.ID, mcav_IPS_GO_sub$GO.ID)
length(unique(UIunique_GO)) # 27 GO terms unique to Uniprot

# Difference in GO terms for Uniprot and IPS - IPS
IUunique <- setdiff(mcav_IPS_GO_sub$GO.ID, uniprot_GO$GO.ID)
length(unique(IUunique)) # 1582 GO terms unique to IPS

## Merge annotations - mcav_annot (B2G + Uniprot) with IPS
full_annot <- merge(mcav_IPS_GO, mcav_annot, by = "gene_id", all.x = T)
full_annot$GO.ID <- paste(full_annot$GO.ID.x, full_annot$GO.ID.y, sep=';') #generate new column with concatenated GO IDs
full_annot <- select(full_annot, -c(GO.ID.x, GO.ID.y))
length(unique(full_annot$gene_id)) # 10791 unique gene ids 


write.csv(full_annot, "~/Desktop/mcav_FullAnnot_20210220.csv")






























IPS_GO <- read.csv("~/Desktop/mcav_GO_20210124.csv", header=TRUE)
IPS_GO <- select(IPS_GO, -X)
colnames(IPS_GO)[1] <-"gene_id"
#IPS_GO$gene_id <- gsub(".m1", "", IPS_GO$gene_id)
#IPS_GO$gene_id <- gsub("model", "TU", IPS_GO$gene_id)
length(unique(IPS_GO$gene_id)) # 10859
IPS_GO <- select(IPS_GO, c("gene_id", "GO.ID")) # taking out source and score tho I want to leave them in--just for comparison though
splitted <- strsplit(as.character(IPS_GO$GO.ID), ",") #split into multiple GO ids
IPS_GO <- data.frame(v1 = rep.int(IPS_GO$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their 
splitted <- strsplit(as.character(IPS_GO$GO.ID), ";") #split into multiple GO ids
IPS_GO <- data.frame(v1 = rep.int(IPS_GO$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their 
colnames(IPS_GO) <- c("gene_id", "GO.ID")
IPS_GO <- unique(IPS_GO)
head(IPS_GO)
nrow(IPS_GO) # 24640 total GO terms from IPS
length(unique(IPS_GO$GO.ID)) # 1766 unique GO terms from IPS
length(unique(IPS_GO$gene_id)) # 10859 unique gene ids from IPS

# From B2G and Uniprot
FuncAnn_UniP_B2G <- read.csv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/final_Annotations/mcav_FuncAnn_UniP_B2G.csv")
length(unique(FuncAnn_UniP_B2G$gene_id)) # 3191
FuncAnn_UniP_B2G_GO <- filter(FuncAnn_UniP_B2G, grepl("GO",GO.ID)) #Select only gnes with GO terms
length(unique(FuncAnn_UniP_B2G_GO$gene_id)) # 147
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
nrow(FuncAnn_UniP_B2G_GO) # 560 total GO terms from B2G + Uniprot
length(unique(FuncAnn_UniP_B2G_GO$GO.ID)) # 212 unique GO terms from B2G + Uniprot
length(unique(FuncAnn_UniP_B2G_GO$gene_id)) # 147 unique gene ids from B2G + Uniprot

# Find intersections and unique results for each method (B2G + Uniprot and IPS)
## GO
# Intersection between GO terms for B2G + Uniprot and IPS
IF_GO <- intersect(IPS_GO$GO.ID, FuncAnn_UniP_B2G_GO$GO.ID)
length(unique(IF_GO)) # 212 similar GO terms between B2G + Uniprot and IPS

# Difference in GO terms for B2G + Uniprot and IPS - FuncAnn first
FIunique_GO <- setdiff(FuncAnn_UniP_B2G_GO$GO.ID, IPS_GO$GO.ID)
length(unique(FIunique_GO)) # 0 GO terms unique to B2G + Uniprot

# Difference in GO terms for B2G + Uniprot and IPS - IPS
IFunique_GO <- setdiff(IPS_GO$GO.ID, FuncAnn_UniP_B2G_GO$GO.ID)
length(unique(IFunique_GO)) # 1554 GO terms unique to IPS

## gene id
# Intersection between GO terms for B2G + Uniprot and IPS
IF_gene<- intersect(IPS_GO$gene_id, FuncAnn_UniP_B2G_GO$gene_id)
length(unique(IF_gene)) # 79 similar gene ids between B2G + Uniprot and IPS

# Difference in gene ids for B2G + Uniprot and IPS - B2G + Uniprot
FIunique_gene <- setdiff(FuncAnn_UniP_B2G_GO$gene_id, IPS_GO$gene_id)
length(unique(FIunique_gene)) # 68 gene ids unique to B2G + Uniprot

# Difference in gene ids for B2G + Uniprot and IPS - IPS
IFunique_gene <- setdiff(IPS_GO$gene_id, FuncAnn_UniP_B2G_GO$gene_id)
length(unique(IFunique_gene)) # 10712 gene ids unique to IPS


# Aggregate results 
agg <- aggregate(IPS_GO$GO.ID, list(IPS_GO$gene_id), paste, collapse = ",")
colnames(agg) <- c("gene_id", "GO.ID")

full_annot <- merge(agg, FuncAnn_UniP_B2G, by = "gene_id", all = T)
full_annot$GO.IDs <- paste(full_annot$GO.ID.x, full_annot$GO.ID.y, sep=';') #generate new column with concatenated GO IDs
full_annot <- select(full_annot, -c(X, GO.ID.x, GO.ID.y))
write.csv(full_annot, file = "~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/final_Annotations/mcav_fullAnnot.csv")

















# 
# 
# 
# ## InterProScan
# IPS <- read.csv("~/Desktop/PutnamLab/FunctionalAnnotation/InterProScan/mcav.interpro.gff3", header = F, sep = "\t")
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
# nrow(IPS_GO) # 24262 total GO terms from IPS
# length(unique(IPS_GO$GO)) # 1721 unique GO terms from IPS
# length(unique(IPS_GO$gene_id)) # 10799 unique gene ids from IPS
# 
# 
# 
# 
# # Find intersections and unique results for each method (IPS and B2G)
# # Cannot compare gene ids because B2G_results_GO has different 'gene ids' because I needed it to pair with Uniprot
# 
# ## GO
# # Intersection between GO terms for B2G and IPS
# BI_GO <- intersect(B2G_results_GO$GO.ID, IPS_GO$GO.ID)
# length(unique(BI_GO)) # 143 similar GO terms between B2G and IPS
# 
# # Difference in GO terms for B2G and IPS - B2G
# BIunique_GO <- setdiff(B2G_results_GO$GO.ID, IPS_GO$GO.ID)
# length(unique(BIunique_GO)) # 37 GO terms unique to B2G
# 
# # Difference in GO terms for B2G and IPS - IPS
# IBunique <- setdiff(IPS_GO$GO.ID, B2G_results_GO$GO.ID)
# length(unique(IBunique)) # 1578 GO terms unique to IPS
# # Cannot compare gene ids because B2G_results_GO has different 'gene ids' because I needed it to pair with Uniprot
# 
# # Find intersections and unique results for each method (IPS and Uniprot)
# # Cannot compare gene ids because Uniprot has different 'gene ids' because I needed it to pair with Uniprot
# ## GO
# # Intersection between GO terms for Uniprot and IPS
# UI_GO <- intersect(uniprot_GO$GO.ID, IPS_GO$GO.ID)
# length(unique(UI_GO)) # 138 similar GO terms between Uniprot and IPS
# 
# # Difference in GO terms for Uniprot and IPS - Uniprot
# UIunique_GO <- setdiff(uniprot_GO$GO.ID, IPS_GO$GO.ID)
# length(unique(UIunique_GO)) # 27 GO terms unique to Uniprot
# 
# # Difference in GO terms for B2G and IPS - IPS
# IUunique_GO <- setdiff(IPS_GO$GO.ID, uniprot_GO$GO.ID)
# length(unique(IUunique_GO)) # 1583 GO terms unique to IPS
# 
# # Find intersections and unique results for each method (IPS and acerv_annot)
# # Cannot compare gene ids because Uniprot has different 'gene ids' because I needed it to pair with Uniprot
# 
# 
# 
# ## GO
# # Before comparison, generate a list of GO terms - acerv_annot
# mcav_annot_GO <- select(mcav_annot, gene_id, GO.ID)
# splitted <- strsplit(as.character(mcav_annot_GO$GO.ID), ";") #split into multiple GO ids
# mcav_annot_GO <- data.frame(v1 = rep.int(mcav_annot_GO$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their GO terms in a single row
# colnames(mcav_annot_GO) <- c("gene_id", "GO.ID")
# mcav_annot_GO <- unique(mcav_annot_GO)
# mcav_annot_GO$GO.ID <- gsub("F:", "", mcav_annot_GO$GO.ID)
# mcav_annot_GO$GO.ID <- gsub("C:", "", mcav_annot_GO$GO.ID)
# mcav_annot_GO$GO.ID <- gsub("P:", "", mcav_annot_GO$GO.ID)
# mcav_annot_GO$GO.ID <- gsub(" ", "", mcav_annot_GO$GO.ID)
# mcav_annot_GO <- filter(mcav_annot_GO, grepl("GO",GO.ID)) #Genes with GO terms - 136
# mcav_annot_GO <- unique(mcav_annot_GO)
# colnames(mcav_annot_GO) <- c("gene_id", "GO.ID")
# 
# # Intersection between GO terms for acerv_annot (b2g + uniprot) and IPS
# AI_GO <- intersect(mcav_annot_GO$GO.ID, IPS_GO$GO.ID)
# length(unique(AI_GO)) # 166 similar GO terms between acerv_annot (b2g + uniprot) and IPS
# 
# # Difference in GO terms for acerv_annot (b2g + uniprot) and IPS - acerv_annot (b2g + uniprot)
# AIunique_GO <- setdiff(mcav_annot_GO$GO.ID, IPS_GO$GO.ID)
# length(unique(AIunique_GO)) # 45 GO terms unique to Uniprot
# 
# # Difference in GO terms for acerv_annot (b2g + uniprot) and IPS - IPS
# IAunique_GO <- setdiff(IPS_GO$GO.ID, mcav_annot_GO$GO.ID)
# length(unique(IAunique_GO)) # 1555 GO terms unique to IPS
# 
# ## gene id
# # Intersection between gene_id for acerv_annot (b2g + uniprot) and IPS
# AI_gene <- intersect(mcav_annot_GO$gene_id, IPS_GO$gene_id)
# length(unique(AI_gene)) # 79 similar gene ids between acerv_annot (b2g + uniprot) and IPS
# 
# # Difference in gene_id terms for acerv_annot (b2g + uniprot) and IPS - acerv_annot (b2g + uniprot)
# AIunique_gene <- setdiff(mcav_annot_GO$gene_id, IPS_GO$gene_id)
# length(unique(AIunique_gene)) # 68 gene ids terms unique to Uniprot
# 
# # Difference in gene ids for acerv_annot (b2g + uniprot) and IPS - IPS
# IAunique_gene <- setdiff(IPS_GO$gene_id, mcav_annot_GO$gene_id)
# length(unique(IAunique_gene)) # 10720 GO terms unique to IPS
# 
# 
# 
# 
# 
# 
# ## Merge Annotations (again) - merge IPS and acerv_annot
# 
# # Merge IPS gene names + GO terms (and any other info I can bring from original IPS file) with acerv_annot gene names + GO terms (and any other info I can bring along )
# IPS_results <- read.csv("~/Desktop/PutnamLab/FunctionalAnnotation/InterProScan/mcav.interpro.gff3", header = F, sep = "\t")
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
# length(unique(IPS_results_GO$seqid)) # 10799
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
# nrow(merge_all_IPS) # 10799
# 
# # merge all info 
# mcav_annot_bu <- select(mcav_annot, c("top_hit", "gene_id", "uniprotkb_entry", "protein_names", "gene_names", "organism", "ko", "kegg", "GO.ID", "GO_terms"))
# nrow(mcav_annot_bu) # 3191
# full_annot <- merge(mcav_annot_bu, merge_all_IPS, by = "gene_id", all.x = T)
# write.csv(full_annot, "~/Desktop/mcav_FuncAnn_UniP_B2G_IPS.csv")
# 
# 
# 
# 
# 
# ## Find metrics for new annotation - b2g + uniprot (mcav_annot)
# new_Sig_Alingments=nrow(ofav_annot)
# new_Genes_with_GO <- nrow(filter(ofav_annot, grepl("GO",GO.ID))) #Genes with GO terms...2172
# new_Genes_with_Kegg <- nrow(filter(ofav_annot, grepl("K",ko))) #Genes with Kegg terms...36L
# new.avg.Eval <- mean(ofav_annot$evalue) # 1.79992e-08
# new.median.Eval <- median(ofav_annot$evalue) # 3e-159
# new.avg.bit <- mean(ofav_annot$bitscore) # 784.2934
# new.median.bit <- median(ofav_annot$bitscore) # 572.8
# # total GO terms
# new_GO <- select(ofav_annot, gene_id, GO.ID)
# splitted <- strsplit(as.character(new_GO$GO.ID), ";") #split into multiple GO ids
# new_GO <- data.frame(v1 = rep.int(new_GO$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their GO terms in a single row
# colnames(new_GO) <- c("gene_id", "GO.ID")
# new_totGO_narm <- filter(new_GO, GO.ID!="NA")
# new_totGO <- nrow(new_totGO_narm) # 1224
# write.csv(new_totGO_narm, "~/Desktop/ofav_new_totGO_narm_UniP_B2G.csv")
# # total unique GO terms
# new_uniqueGO <- unique(new_totGO_narm$GO.ID)
# new_uniqueGO <- length(new_uniqueGO) # 498
# 
# # total Kegg terms
# new_Kegg <- select(ofav_annot, gene_id, ko)
# splitted <- strsplit(as.character(new_Kegg$ko), ";") #split into multiple Kegg ids
# new.Keggterms <- data.frame(v1 = rep.int(new_Kegg$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their Kegg terms in a single row
# colnames(new.Keggterms) <- c("gene_id", "Kegg.ID")
# new.Keggterms$Kegg.ID <- replace_na(new.Keggterms$Kegg.ID, "NA")
# new_totKegg_narm <- filter(new.Keggterms, Kegg.ID!="NA")
# new_totKegg <- nrow(new_totKegg_narm) # 0
# #write.csv(new_totKegg_narm, "~/Desktop/mcav_new_totKegg_narm_UniP_B2G.csv")
# #new_uniqueKegg <- unique(new_totKegg_narm$Kegg.ID)
# #new_uniqueKegg <- length(new_uniqueKegg) # 0
# 
# 
# 
# 
# ## Find metrics for new annotation - b2g + uniprot (acerv_annot)
# new_Sig_Alingments=nrow(full_annot)
# new_Genes_with_GO <- nrow(filter(ofav_annot, grepl("GO",GO.ID))) #Genes with GO terms...2172
# new_Genes_with_Kegg <- nrow(filter(ofav_annot, grepl("K",ko))) #Genes with Kegg terms...36L
# new.avg.Eval <- mean(ofav_annot$evalue) # 1.79992e-08
# new.median.Eval <- median(ofav_annot$evalue) # 2.9e-111
# new.avg.bit <- mean(ofav_annot$bitscore) # 571.9534
# new.median.bit <- median(mcav_annot$bitscore) # 412.9
# # total GO terms
# new_GO <- select(mcav_annot, gene_id, GO.ID)
# splitted <- strsplit(as.character(new_GO$GO.ID), ";") #split into multiple GO ids
# new_GO <- data.frame(v1 = rep.int(new_GO$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their GO terms in a single row
# colnames(new_GO) <- c("gene_id", "GO.ID")
# new_totGO_narm <- filter(new_GO, GO.ID!="NA")
# new_totGO <- nrow(new_totGO_narm) # 801
# write.csv(new_totGO_narm, "~/Desktop/mcav_new_totGO_narm_UniP_B2G_IPS.csv")
# # total unique GO terms
# new_uniqueGO <- unique(new_totGO_narm$GO.ID)
# new_uniqueGO <- length(new_uniqueGO) # 382
# 
# # total Kegg terms
# new_Kegg <- select(mcav_annot, gene_id, ko)
# splitted <- strsplit(as.character(new_Kegg$ko), ";") #split into multiple Kegg ids
# new.Keggterms <- data.frame(v1 = rep.int(new_Kegg$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their Kegg terms in a single row
# colnames(new.Keggterms) <- c("gene_id", "Kegg.ID")
# new.Keggterms$Kegg.ID <- replace_na(new.Keggterms$Kegg.ID, "NA")
# new_totKegg_narm <- filter(new.Keggterms, Kegg.ID!="NA")
# new_totKegg <- nrow(new_totKegg_narm) # 0
# write.csv(new_totKegg_narm, "~/Desktop/mcav_new_totKegg_narm_UniP_B2G.csv")
# new_uniqueKegg <- unique(new_totKegg_narm$Kegg.ID)
# new_uniqueKegg <- length(new_uniqueKegg) # 0
# 
# 
# 
# 
