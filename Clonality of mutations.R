###For Keefe. Patients with mutations in these genes (clonal and then subclonal)

library(e1071)

##Genes from the screen. Supplementary table S9
AIS_escape_genes <- read.delim("Genes_duplex_screen_hits.txt", header=FALSE)

##Table of Human gene ids (associated gene names and entrez)
human_genes <- read.delim("Human_ensembl_(GRCh37.p13).txt", sep=",")

AIS_escape_genes$GeneID <- as.character(human_genes$Associated.Gene.Name[match(AIS_escape_genes$V1, human_genes$EntrezGene.ID)])

AIS_escape_genes[AIS_escape_genes$V1 == "201175","GeneID"] <- "SH3D20"

AIS_escape_genes[AIS_escape_genes$V1 == "83954","GeneID"] <- "FKSG83"

##Object from TRACERx available from REVOLVER's github
load("TRACERx.Rdata")

##Collapsing per patient
TRACERx <- TRACERx[,c("patientID", "variantID", "is.clonal")]
TRACERx <- unique(TRACERx)

genes <- unique(TRACERx$variantID)

##Calculating the ratio of clonal vs. non-clonal
ratio_distribution <- vector()
for(gene in genes){
  temp <- TRACERx[TRACERx$variantID == gene,]
  if(nrow(temp) >= 3){  ##at least 3 patients
    ratio <- sum(temp$is.clonal == TRUE)/sum(temp$is.clonal == FALSE)
    ratio_distribution <- c(ratio_distribution, ratio)
    names(ratio_distribution)[length(ratio_distribution)] <- gene
  }
}

ratio_distribution_df <- data.frame(Genes = names(ratio_distribution),
                                    Ratio = ratio_distribution,
                                    InScreen = ifelse(names(ratio_distribution) %in% c(AIS_escape_genes$GeneID), "YES", "No"))

##Excluding cases where all events are clonal
ratio_distribution_df2 <- ratio_distribution_df[is.finite(ratio_distribution_df$Ratio),]
sum(ratio_distribution_df2$InScreen=="YES")
#36 

#Actual values
actual <- ratio_distribution_df2[ratio_distribution_df2$InScreen=="YES",]

##Calculate how the distribution would like by chance
mean_chance <- vector()
skewness_chance <- vector()
chance_data <- vector()
for(i in 1:100000){
  temp <- sample(ratio_distribution_df2$Ratio, 36, replace=FALSE)
  skewness_chance <- c(skewness_chance, skewness(temp))
  mean_chance <- c(mean_chance, mean(temp))
  chance_data <- rbind(chance_data, cbind(temp, i))
  print(i)
}

##Calculating empirical p-value
sum(skewness_chance >= skewness(actual$Ratio))
#2872/100000

mean(actual$Ratio)

mean(mean_chance)

##Plot distributions
plot(density(chance_data[chance_data[,2] == 1, 1]), ylim=c(0,0.5),
     xlab="Clonality/subclonality ratios", main="")
for(i in 1:100){
  lines(density(chance_data[chance_data[,2] == i, 1]))
}
lines(density(actual$Ratio), col="red", lwd=2)


