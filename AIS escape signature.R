##Calculate AKT senescence signature across TCGA patients
library(GSVA)
library(reshape2)
library(ggplot2)

###Read in gene expression of all patients

tumour <- "LUSC"

##Loading gene expression data from TCGA (RSEM genes normalized)
##and fomatting (e.g. removing genes without gene names and selecting
##tumour samples.
path <- paste("gdac.broadinstitute.org_", tumour, 
                ".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/", 
                tumour, 
                ".rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",
                sep="")
##Load normalized gene expression data from TCGA
exp <- read.delim(path) 
exp <- exp[-1,]

#Select tumour samples and format input
exp_tumour <- exp[,substr(colnames(exp), 14, 16) %in% c("01A", "01B", "01C", "01D")]

colnames(exp_tumour) <- substr(colnames(exp_tumour),9,12)
  
all_genes <- as.character(exp[,1])
  
all_genes <- unlist(strsplit(all_genes, "\\|"))[seq(1, length(all_genes)*2, 2)]

##Removing genes without a gene name and duplicated genes
genes_remove <- c("?", "SLC35E2")
genes_remove_index <- which(all_genes %in% genes_remove)
all_genes <- all_genes[-genes_remove_index]
  
exp_tumour <- exp_tumour[-genes_remove_index,]
  
if(!is.null(nrow(exp_tumour))){
  rownames(exp_tumour) <- all_genes
}
  
expression <- exp_tumour

  
##Loading list of patients with mutation and CNV information in TCGA to 
##select patients with gene expression, mutation and CNV data.

##These objects were obtained as described in Trigos et al. eLife 2019. DOI: 10.7554/eLife.40947
load("patients_with_mut_info.Rdata")
load("patients_with_CNV_info.Rdata")

###Patients with NF1, FADD and CCAR1 mutations obtained from cBioPortal
TCGA_mutations <- list()
TCGA_mutations[["NF1"]] <- #
TCGA_mutations[["FADD"]] <- #
TCGA_mutations[["CCAR1"]] <- #
TCGA_interest <- rbind(cbind(TCGA_mutations[["NF1"]], "NF1"),
                       cbind(TCGA_mutations[["FADD"]], "FADD"),
                       cbind(TCGA_mutations[["CCAR1"]], "CCAR1"))
colnames(TCGA_interest) <- c("Patient", "Mutation")
TCGA_interest <- as.data.frame(TCGA_interest)
TCGA_interest$Patient <- substr(TCGA_interest$Patient, 9, 12)

pat_mut <- patients_with_mut_info[[tumour]]
pat_CNV <- patients_with_CNV_info[[tumour]]
common <- pat_mut[pat_mut %in% pat_CNV]
expression2 <- expression[,colnames(expression) %in% common]


##Loading signatures

##Genes identified from screen. Supplementary table S9.
AIS_escape_genes <- read.delim("Genes_duplex_screen_hits.txt", header=FALSE)

##Table of Human gene ids (associated gene names and entrez)
human_genes <- read.delim("Human_ensembl_(GRCh37.p13).txt", sep=",")

AIS_escape_genes$GeneID <- as.character(human_genes$Associated.Gene.Name[match(AIS_escape_genes$V1, human_genes$EntrezGene.ID)])
AIS_escape_genes[AIS_escape_genes$V1 == "201175","GeneID"] <- "SH3D20"
AIS_escape_genes[AIS_escape_genes$V1 == "83954","GeneID"] <- "FKSG83"

##Define Akt signature as the genes identified from the screen
akt_sig <- list(akt_sig=AIS_escape_genes$GeneID)

expression2 <- as.data.frame(expression2)
expression2 <- apply(expression2, 2, function(x){as.numeric(as.character(x))})
rownames(expression2) <- rownames(expression)

#Calculating ssGSEA scores using GSVA
ssGSEA <- gsva(expression2, akt_sig, method="ssgsea", kcdf="Gaussian", ssgsea.norm=TRUE)

ssGSEA_df <- melt(ssGSEA)
colnames(ssGSEA_df) <- c("Signature", "Patient", "ssGSEA")

ssGSEA_df$Mutation <- "None"
ssGSEA_df[ssGSEA_df$Patient %in% substr(TCGA_mutations[["NF1"]], 9, 12),"Mutation"] <- "NF1"
ssGSEA_df[ssGSEA_df$Patient %in% substr(TCGA_mutations[["FADD"]], 9, 12),"Mutation"] <- "FADD"
ssGSEA_df[ssGSEA_df$Patient %in% substr(TCGA_mutations[["CCAR1"]], 9, 12),"Mutation"] <- "CCAR1"


wilcox.test(ssGSEA_df[ssGSEA_df$Mutation == "None", "ssGSEA"],
                    ssGSEA_df[ssGSEA_df$Mutation != "None", "ssGSEA"], alternative="less")$p.value

ssGSEA_df$Tumour <- "LUSC"

ggplot(ssGSEA_df, aes(x=Mutation, y=ssGSEA))+
  geom_boxplot(aes(fill=Mutation))+
  ggtitle("Akt signature")+
  facet_grid(.~Tumour)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_bw()

