##Differential expression using DEseq2

library(biomaRt)
library(DESeq2)


run_DEseq2 <- function(input, coding_genes, levels, output_title){
  genes <- rownames(input)
  input <- input[rownames(input) %in% coding_genes$ensembl_gene_id,]
  genes <- genes[genes %in% coding_genes$ensembl_gene_id]
  input <- apply(input, 2, function(x){
    as.numeric(as.character(x))
  })
  rownames(input) <- genes
  
  ##Keep genes with more than 20 reads in either condition
  input_rel_level2 <- input[, 1:3]
  input_rel_level1 <- input[, 4:6]
  
  input_rel_level2 <- apply(input_rel_level2, 1, mean)
  input_rel_level1 <- apply(input_rel_level1, 1, mean)
  
  genes_more_20_level2 <- names(input_rel_level2[input_rel_level2 >= 20])
  genes_more_20_level1 <- names(input_rel_level1[input_rel_level1 >= 20])
  genes_keep <- unique(c(genes_more_20_level2, genes_more_20_level1))
  input_rel_more_20 <- input[rownames(input) %in% genes_keep,]
  
  if(output_title != "AKT_control_all_columns"){
    condition = factor(substr(colnames(input), 1, nchar(colnames(input))-2))
  }else{
    condition = factor(substr(colnames(input), 1, nchar(colnames(input))-1))
  }
  
  dds <- DESeqDataSetFromMatrix(input_rel_more_20, DataFrame(condition), ~ condition)
  
  if(output_title != "AKT_control_all_columns"){
    dds$condition <- relevel(dds$condition, ref = "pBABE.NTC") 
  }else{
    dds$condition <- relevel(dds$condition, ref = "BJ.T") 
  }
  
  dds <- DESeq(dds)
  res <- as.data.frame(results(dds))
  
  input_DE <- data.frame(Genes = rownames(input),
                         input, 
                         res[match(rownames(input), rownames(res)),])
  
  ##Calculate RPKM
  input_RPKM <- input
  input_RPKM <- input_RPKM*10^9
  
  gene_length <- coding_genes[match(rownames(input_RPKM), coding_genes$ensembl_gene_id),]
  
  input_RPKM <- input_RPKM/gene_length$transcript_length
  
  for(i in 1:ncol(input_RPKM)){
    input_RPKM[,i] <- input_RPKM[,i]/sum(input[,i], na.rm=TRUE)
  }
  
  input_RPKM <- input_RPKM[!is.na(input_RPKM[,1]),]
  colnames(input_RPKM) <- paste("RPKM.", colnames(input_RPKM), sep="")
  
  input_DE <- data.frame(input_DE,
                         input_RPKM[match(input_DE$Genes, rownames(input_RPKM)),])
  
  input_DE <- data.frame(input_DE, 
                         gene_length[match(input_DE$Genes, gene_length$ensembl_gene_id),])
  
  if(output_title != "AKT_control_all_columns"){
    input_DE <- input_DE[,c("Genes", "hgnc_symbol", paste(levels[1], 1:3, sep="."),
                            paste(levels[2],  1:3, sep="."),
                            "transcript_length",
                            paste("RPKM", levels[1], 1:3, sep="."),
                            paste("RPKM", levels[2], 1:3, sep="."),
                            "chromosome_name",  "start_position",  "end_position",  "strand",
                            "ensembl_gene_id", "description",
                            "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue","padj")]
    
  }else{
    input_DE <- input_DE[,c("Genes", "hgnc_symbol", "BJ.T1", "BJ.T2", "BJ.T3",          
                            "BJ.T.AKT1",  "BJ.T.AKT2",  "BJ.T.AKT3",
                            "transcript_length",
                            "RPKM.BJ.T1","RPKM.BJ.T2","RPKM.BJ.T3",
                            "RPKM.BJ.T.AKT1","RPKM.BJ.T.AKT2","RPKM.BJ.T.AKT3",
                            "chromosome_name",  "start_position",  "end_position",  "strand",
                            "ensembl_gene_id", "description",
                            "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue","padj")]
  }
  
  write.table(input_DE, file=paste(output_title, "DEseq2.csv", sep=""),
              quote=FALSE, sep=",", row.names=FALSE)
}



##Obtain coding genes
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

genedesc <- getBM(attributes=c('chromosome_name', 'start_position', 'end_position',
                               'strand', 'ensembl_gene_id', 'external_gene_name',
                               'hgnc_symbol', 'description', 'gene_biotype',
                               'transcript_length'),
                  filters = 'ensembl_gene_id', values = rownames(input),
                  mart =ensembl)

coding_genes <- genedesc[genedesc$gene_biotype == "protein_coding",]

##Select the longest canonical transcript
coding_genes <- aggregate(transcript_length ~ chromosome_name+start_position+end_position+strand+ensembl_gene_id+hgnc_symbol+description+gene_biotype, coding_genes, max)

##Load HRAS and pBABE raw counts. These were obtained using featureCounts in Galaxy
##Formatted with genes as rows, and samples as columns
input_HRAS_pBABE <- ###

##Load HRAS and pBABE raw counts. These were obtained using featureCounts in Galaxy
input_AKT_control <- ###


run_DEseq2(input_HRAS_pBABE, coding_genes, c("pBABE.NTC", "HRAS.NTC"), "HRAS_pBABE_all_columns")

run_DEseq2(input_AKT_control, coding_genes, c("BJ.T", "BJ.T.AKT"), "AKT_control_all_columns")
  