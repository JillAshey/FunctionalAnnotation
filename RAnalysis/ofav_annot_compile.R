# Title: Ofav Annotation Compilation
# Author: Jill Ashey
# Date: 12/09/2020



## Species: O. favelota 

#This script takes the results of functional annotation services and combines the results. Nucleotide CDS sequences were annotated using DIAMONDSEARCH BLASTX, resulting in 10122 hits. These hits were used as input into:
#  - Uniprot
#  - Blast2GO

#Additional annotation was provided by
#  - InterProScan





## Load libraries
library(tidyverse)
library(dplyr)



## Diamond BLAST

# Load DIAMOND BLAST results
ofav_blast <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Diamond/Ofav_annot.tab", col_names = FALSE)
colnames(ofav_blast) <- c("seqName", "top_hit", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",  "qlen", "slen")
head(ofav_blast)
dim(ofav_blast) #10122 x 14



## Uniprot

# Uniprot mapping occurred on 11/10/20.
#EMBL/GenBank/DDBJ CDS protein IDs generated from Diamond BLAST were mapped to UniProtKB database IDs.

#Because there were many Diamond BLAST hits for Acerv, I had to break up the tab file into chunks that Uniprot could handle. So there are 15 acerv_Uniprot files to be read and then rbind() together
u1 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/ofav_Uniprot_1.tab", col_names = TRUE)
colnames(u1) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u1)
dim(u1) # 30 x 12

u2 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/ofav_Uniprot_2.tab", col_names = TRUE)
colnames(u2) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u2)
dim(u2) # 15 x 12

u3 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/ofav_Uniprot_3.tab", col_names = TRUE)
colnames(u3) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u3)
dim(u3) # 18 x 12

u4 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/ofav_Uniprot_4.tab", col_names = TRUE)
colnames(u4) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u4)
dim(u4) # 15 x 12

u5 <- read_tsv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Uniprot/ofav_Uniprot_5.tab", col_names = TRUE)
colnames(u5) <- c("my_list","top_hit", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "length", "go_ids", "gene_ontology", "ko", "kegg")
head(u5)
dim(u5) # 217 x 12

#Compile the Uniprot files 
uniprot_results <- bind_rows(u1,u2,u3,u4,u5)
uniprot_results <- unique(uniprot_results)
head(uniprot_results)
dim(uniprot_results) # 294 x 12
uniprot_results <- filter(uniprot_results, grepl("GO",go_ids)) #Select only gnes with GO terms
dim(uniprot_results) # 178 x 12

# Generate a list of GO terms - UniProt
uniprot_GO <- select(uniprot_results, my_list, go_ids)
splitted <- strsplit(as.character(uniprot_GO$go_ids), ";") #split into multiple GO ids
uniprot_GO <- data.frame(v1 = rep.int(uniprot_GO$my_list, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their GO terms in a single row
uniprot_GO <- unique(uniprot_GO)
colnames(uniprot_GO) <- c("gene_id", "GO.ID")
uniprot_GO$GO.ID <- gsub(" ", "", uniprot_GO$GO.ID)
nrow(uniprot_GO) # 452 total GO terms from Uniprot
length(unique(uniprot_GO$GO.ID)) # 212 unique GO terms from Uniprot
length(unique(uniprot_GO$gene_id)) # 178 unique gene ids from Uniprot 
# Not technically gene ids. Uniprot has no gene id info--actually ids from Uniprot. When I combine the Uniport and B2G files, the uniprot ids will then be associated with gene ids 





## Blast2GO

#Nucleotide CDS sequences were annotated using DIAMONDSEARCH BLASTX, resulting in 3191 These hits were used as input into Blast2GO to obtain GO terms using the 12/05/2020 obo database.

B2G_results <- read.csv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/Blast2GO/ofav_blast2go_table.csv")
B2G_results <- select(B2G_results, c("SeqName", "Description", "Length", "e.Value", "sim.mean", "GO.IDs", "GO.Names"))
colnames(B2G_results) <- c("seqName", "top_hit", "length", "eValue", "simMean", "GO.ID", "GO_names")
head(B2G_results)
dim(B2G_results) # 10122 x 7
B2G_results <- filter(B2G_results, grepl("GO",GO.ID)) #Genes with GO terms - 136
dim(B2G_results) # 228 x 7

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
nrow(B2G_results_GO) # 542 total GO terms from B2G
length(unique(B2G_results_GO$GO.ID)) # 236 unique GO terms from B2G
length(unique(B2G_results_GO$gene_id)) # 186 unique gene ids from B2G

# Find intersections and unique results for each method (Uniprot and B2G)
## GO
# Intersection between GO terms for B2G and Uniprot
BU_GO <- intersect(B2G_results_GO$GO.ID, uniprot_GO$GO.ID)
length(unique(BU_GO)) # 168 similar GO terms between B2G and Uniprot

# Difference in GO terms for B2G and Uniprot - B2G
BUunique_GO <- setdiff(B2G_results_GO$GO.ID, uniprot_GO$GO.ID)
length(unique(BUunique_GO)) # 68 GO terms unique to B2G

# Difference in GO terms for B2G and Uniprot - Uniprot
UBunique <- setdiff(uniprot_GO$GO.ID, B2G_results_GO$GO.ID)
length(unique(UBunique)) # 44 GO terms unique to Uniprot

## gene id
# Intersection between gene id for B2G and Uniprot
BU_gene <- intersect(B2G_results_GO$gene_id, uniprot_GO$gene_id)
length(unique(BU_gene)) # 177 similar gene ids between B2G and Uniprot

# Difference in gene ids for B2G and Uniprot
BUunique_gene <- setdiff(B2G_results_GO$gene_id, uniprot_GO$gene_id)
length(unique(BUunique_gene)) # 9 gene ids unique to B2G

# Difference in gene ids for B2G and Uniprot
UBunique_gene <- setdiff(uniprot_GO$gene_id, B2G_results_GO$gene_id)
length(unique(UBunique_gene)) # 1 gene ids unique to uniprot





## Merge Annotations - uniprot + b2g
ofav_annot <- left_join(ofav_blast, B2G_results, by="seqName")
ofav_annot <- select(ofav_annot, seqName, top_hit.x, length.x, evalue, bitscore, simMean, GO.ID, GO_names)
colnames(ofav_annot) <- c("gene_id", "top_hit", "length", "evalue", "bitscore", "simMean", "GO.ID", "GO_names")
uniprot_results <- select(uniprot_results, -top_hit)
uniprot_results <- rename(uniprot_results, "top_hit"="my_list")
ofav_annot <- merge(ofav_annot, uniprot_results, by="top_hit", all.x = T)
ofav_annot$GO.IDs <- paste(ofav_annot$GO.ID, ofav_annot$go_ids, sep=';') #generate new column with concatenated GO IDs
ofav_annot$GO_terms <- paste(ofav_annot$GO_names, ofav_annot$gene_ontology, sep=';') #generate new column with concatenated GO IDs
ofav_annot <- select(ofav_annot, c("top_hit", "gene_id", "length.x", "evalue", "bitscore", "simMean", "uniprotkb_entry", "status", "protein_names", "gene_names", "organism", "ko", "kegg", "GO.IDs", "GO_terms"))
ofav_annot <- rename(ofav_annot, "GO.ID"="GO.IDs")
ofav_annot <- rename(ofav_annot, "length"="length.x")
names(ofav_annot)
head(ofav_annot)
tail(ofav_annot)
dim(ofav_annot) # 10122 x 15
write.csv(ofav_annot, "~/Desktop/ofav_FuncAnn_UniP_B2G.csv")









# IPS
IPS_GO <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOseq/ofav_GOterms.csv", header=TRUE)
IPS_GO <- select(IPS_GO, -X)
colnames(IPS_GO)[1] <-"gene_id"
#IPS_GO$gene_id <- gsub(".m1", "", IPS_GO$gene_id)
#IPS_GO$gene_id <- gsub("model", "TU", IPS_GO$gene_id)
length(unique(IPS_GO$gene_id)) # 18424
IPS_GO <- select(IPS_GO, c("gene_id", "GO_term")) # taking out source and score tho I want to leave them in--just for comparison though
splitted <- strsplit(as.character(IPS_GO$GO), ",") #split into multiple GO ids
IPS_GO <- data.frame(v1 = rep.int(IPS_GO$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their 
colnames(IPS_GO) <- c("gene_id", "GO.ID")
IPS_GO <- unique(IPS_GO)
head(IPS_GO)
nrow(IPS_GO) # 45666 total GO terms from IPS
length(unique(IPS_GO$GO.ID)) # 2002 unique GO terms from IPS
length(unique(IPS_GO$gene_id)) # 18424 unique gene ids from IPS

# From B2G and Uniprot
FuncAnn_UniP_B2G <- read.csv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/final_Annotations/ofav_FuncAnn_UniP_B2G.csv")
length(unique(FuncAnn_UniP_B2G$gene_id)) # 10122
FuncAnn_UniP_B2G_GO <- filter(FuncAnn_UniP_B2G, grepl("GO",GO.ID)) #Select only gnes with GO terms
length(unique(FuncAnn_UniP_B2G_GO$gene_id)) # 229
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
nrow(FuncAnn_UniP_B2G_GO) # 838 total GO terms from B2G + Uniprot
length(unique(FuncAnn_UniP_B2G_GO$GO.ID)) # 281 unique GO terms from B2G + Uniprot
length(unique(FuncAnn_UniP_B2G_GO$gene_id)) # 229 unique gene ids from B2G + Uniprot

# Find intersections and unique results for each method (B2G + Uniprot and IPS)
## GO
# Intersection between GO terms for B2G + Uniprot and IPS
IF_GO <- intersect(IPS_GO$GO.ID, FuncAnn_UniP_B2G_GO$GO.ID)
length(unique(IF_GO)) # 208 similar GO terms between B2G + Uniprot and IPS

# Difference in GO terms for B2G + Uniprot and IPS - FuncAnn first
FIunique_GO <- setdiff(FuncAnn_UniP_B2G_GO$GO.ID, IPS_GO$GO.ID)
length(unique(FIunique_GO)) # 73 GO terms unique to B2G + Uniprot

# Difference in GO terms for B2G + Uniprot and IPS - IPS
IFunique_GO <- setdiff(IPS_GO$GO.ID, FuncAnn_UniP_B2G_GO$GO.ID)
length(unique(IFunique_GO)) # 1794 GO terms unique to Uniprot

## gene id
# Intersection between GO terms for B2G + Uniprot and IPS
IF_gene<- intersect(IPS_GO$gene_id, FuncAnn_UniP_B2G_GO$gene_id)
length(unique(IF_gene)) # 0 similar gene ids between B2G + Uniprot and IPS

# Difference in gene ids for B2G + Uniprot and IPS - B2G + Uniprot
FIunique_gene <- setdiff(FuncAnn_UniP_B2G_GO$gene_id, IPS_GO$gene_id)
length(unique(FIunique_gene)) # 229 gene ids unique to B2G + Uniprot

# Difference in gene ids for B2G + Uniprot and IPS - IPS
IFunique_gene <- setdiff(IPS_GO$gene_id, FuncAnn_UniP_B2G_GO$gene_id)
length(unique(IFunique_gene)) # 18424 gene ids unique to IPS


# Aggregate results 
agg <- aggregate(IPS_GO$GO.ID, list(IPS_GO$gene_id), paste, collapse = ",")
colnames(agg) <- c("gene_id", "GO.ID")

full_annot <- merge(agg, FuncAnn_UniP_B2G, by = "gene_id", all.x = T)
full_annot$GO.IDs <- paste(full_annot$GO.ID.x, full_annot$GO.ID.y, sep=';') #generate new column with concatenated GO IDs
full_annot <- select(full_annot, -c(X, GO.ID.x, GO.ID.y))
write.csv(full_annot, file = "~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/final_Annotations/ofav_fullAnnot.csv")













# 
# ## InterProScan
# IPS <- read.csv("~/Desktop/PutnamLab/FunctionalAnnotation/InterProScan/ofav.interpro.gff3", header = F, sep = "\t")
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
# nrow(IPS_GO) # 45680 total GO terms from IPS
# length(unique(IPS_GO$GO)) # 2003 unique GO terms from IPS
# length(unique(IPS_GO$gene_id)) # 18437 unique gene ids from IPS
# 
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
# length(unique(BI_GO)) # 178 similar GO terms between B2G and IPS
# 
# # Difference in GO terms for B2G and IPS - B2G
# BIunique_GO <- setdiff(B2G_results_GO$GO.ID, IPS_GO$GO.ID)
# length(unique(BIunique_GO)) # 58 GO terms unique to B2G
# 
# # Difference in GO terms for B2G and IPS - IPS
# IBunique <- setdiff(IPS_GO$GO.ID, B2G_results_GO$GO.ID)
# length(unique(IBunique)) # 1825 GO terms unique to IPS
# # Cannot compare gene ids because B2G_results_GO has different 'gene ids' because I needed it to pair with Uniprot
# 
# # Find intersections and unique results for each method (IPS and Uniprot)
# # Cannot compare gene ids because Uniprot has different 'gene ids' because I needed it to pair with Uniprot
# 
# ## GO
# # Intersection between GO terms for Uniprot and IPS
# UI_GO <- intersect(uniprot_GO$GO.ID, IPS_GO$GO.ID)
# length(unique(UI_GO)) # 163 similar GO terms between Uniprot and IPS
# 
# # Difference in GO terms for Uniprot and IPS - Uniprot
# UIunique_GO <- setdiff(uniprot_GO$GO.ID, IPS_GO$GO.ID)
# length(unique(UIunique_GO)) # 49 GO terms unique to Uniprot
# 
# # Difference in GO terms for B2G and IPS - IPS
# IUunique_GO <- setdiff(IPS_GO$GO.ID, uniprot_GO$GO.ID)
# length(unique(IUunique_GO)) # 1840 GO terms unique to IPS
# 
# # Find intersections and unique results for each method (IPS and acerv_annot)
# # Cannot compare gene ids because Uniprot has different 'gene ids' because I needed it to pair with Uniprot
# 
# 
# 
# 
# ## GO
# # Before comparison, generate a list of GO terms - acerv_annot
# ofav_annot_GO <- select(ofav_annot, gene_id, GO.ID)
# splitted <- strsplit(as.character(ofav_annot_GO$GO.ID), ";") #split into multiple GO ids
# ofav_annot_GO <- data.frame(v1 = rep.int(ofav_annot_GO$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their GO terms in a single row
# colnames(ofav_annot_GO) <- c("gene_id", "GO.ID")
# ofav_annot_GO <- unique(ofav_annot_GO)
# ofav_annot_GO$GO.ID <- gsub("F:", "", ofav_annot_GO$GO.ID)
# ofav_annot_GO$GO.ID <- gsub("C:", "", ofav_annot_GO$GO.ID)
# ofav_annot_GO$GO.ID <- gsub("P:", "", ofav_annot_GO$GO.ID)
# ofav_annot_GO$GO.ID <- gsub(" ", "", ofav_annot_GO$GO.ID)
# ofav_annot_GO <- filter(ofav_annot_GO, grepl("GO",GO.ID)) #Genes with GO terms - 1224
# ofav_annot_GO <- unique(ofav_annot_GO)
# colnames(ofav_annot_GO) <- c("gene_id", "GO.ID")
# 
# # Intersection between GO terms for acerv_annot (b2g + uniprot) and IPS
# OI_GO <- intersect(ofav_annot_GO$GO.ID, IPS_GO$GO.ID)
# length(unique(OI_GO)) # 208 similar GO terms between acerv_annot (b2g + uniprot) and IPS
# 
# # Difference in GO terms for acerv_annot (b2g + uniprot) and IPS - acerv_annot (b2g + uniprot)
# OIunique_GO <- setdiff(ofav_annot_GO$GO.ID, IPS_GO$GO.ID)
# length(unique(OIunique_GO)) # 72 GO terms unique to Uniprot
# 
# # Difference in GO terms for acerv_annot (b2g + uniprot) and IPS - IPS
# IOunique_GO <- setdiff(IPS_GO$GO.ID, ofav_annot_GO$GO.ID)
# length(unique(IOunique_GO)) # 1795 GO terms unique to IPS
# 
# ## gene id
# # Intersection between gene_id for acerv_annot (b2g + uniprot) and IPS
# OI_gene <- intersect(ofav_annot_GO$gene_id, IPS_GO$gene_id)
# length(unique(OI_gene)) # 0 similar gene ids between acerv_annot (b2g + uniprot) and IPS
# 
# # Difference in gene_id terms for acerv_annot (b2g + uniprot) and IPS - acerv_annot (b2g + uniprot)
# OIunique_gene <- setdiff(ofav_annot_GO$gene_id, IPS_GO$gene_id)
# length(unique(OIunique_gene)) # 229 gene ids terms unique to Uniprot
# 
# # Difference in gene ids for acerv_annot (b2g + uniprot) and IPS - IPS
# IOunique_gene <- setdiff(IPS_GO$gene_id, ofav_annot_GO$gene_id)
# length(unique(IOunique_gene)) # 18437 GO terms unique to IPS
# 
# 
# 
# 
# 
# 
# ## Merge Annotations (again) - merge IPS and acerv_annot
# 
# # Merge IPS gene names + GO terms (and any other info I can bring from original IPS file) with acerv_annot gene names + GO terms (and any other info I can bring along )
# IPS_results <- read.csv("~/Desktop/PutnamLab/FunctionalAnnotation/InterProScan/ofav.interpro.gff3", header = F, sep = "\t")
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
# length(unique(IPS_results_GO$seqid)) # 18437
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
# nrow(merge_all_IPS) # 18437
# 
# # merge all info 
# ofav_annot_bu <- select(ofav_annot, c("top_hit", "gene_id", "uniprotkb_entry", "protein_names", "gene_names", "organism", "ko", "kegg", "GO.ID", "GO_terms"))
# nrow(ofav_annot_bu) # 10122
# full_annot <- merge(ofav_annot_bu, merge_all_IPS, by = "gene_id", all.x = T)
# write.csv(full_annot, "~/Desktop/ofav_FuncAnn_UniP_B2G_IPS.csv")
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
# # no evalues or bit scores in full_annot file 
# ## Find metrics for new annotation - b2g + uniprot (acerv_annot) + IPS
# # new_Sig_Alingments=nrow(full_annot)
# # new_Genes_with_GO <- nrow(filter(ofav_annot, grepl("GO",GO.ID))) #Genes with GO terms...2172
# # new_Genes_with_Kegg <- nrow(filter(ofav_annot, grepl("K",ko))) #Genes with Kegg terms...36L
# # new.avg.Eval <- mean(ofav_annot$evalue) # 1.79992e-08
# # new.median.Eval <- median(ofav_annot$evalue) # 2.9e-111
# # new.avg.bit <- mean(ofav_annot$bitscore) # 571.9534
# # new.median.bit <- median(mcav_annot$bitscore) # 412.9
# # # total GO terms
# # new_GO <- select(mcav_annot, gene_id, GO.ID)
# # splitted <- strsplit(as.character(new_GO$GO.ID), ";") #split into multiple GO ids
# # new_GO <- data.frame(v1 = rep.int(new_GO$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their GO terms in a single row
# # colnames(new_GO) <- c("gene_id", "GO.ID")
# # new_totGO_narm <- filter(new_GO, GO.ID!="NA")
# # new_totGO <- nrow(new_totGO_narm) # 801
# # write.csv(new_totGO_narm, "~/Desktop/mcav_new_totGO_narm_UniP_B2G_IPS.csv")
# # # total unique GO terms
# # new_uniqueGO <- unique(new_totGO_narm$GO.ID)
# # new_uniqueGO <- length(new_uniqueGO) # 382
# # 
# # # total Kegg terms
# # new_Kegg <- select(mcav_annot, gene_id, ko)
# # splitted <- strsplit(as.character(new_Kegg$ko), ";") #split into multiple Kegg ids
# # new.Keggterms <- data.frame(v1 = rep.int(new_Kegg$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their Kegg terms in a single row
# # colnames(new.Keggterms) <- c("gene_id", "Kegg.ID")
# # new.Keggterms$Kegg.ID <- replace_na(new.Keggterms$Kegg.ID, "NA")
# # new_totKegg_narm <- filter(new.Keggterms, Kegg.ID!="NA")
# # new_totKegg <- nrow(new_totKegg_narm) # 0
# # write.csv(new_totKegg_narm, "~/Desktop/mcav_new_totKegg_narm_UniP_B2G.csv")
# # new_uniqueKegg <- unique(new_totKegg_narm$Kegg.ID)
# # new_uniqueKegg <- length(new_uniqueKegg) # 0
# 
# 
