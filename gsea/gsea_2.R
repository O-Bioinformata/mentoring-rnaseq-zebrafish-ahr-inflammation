library(clusterProfiler)
library(org.Dr.eg.db)
library(msigdbr)
library(dplyr)

file <- "DEGs_ctCutUncut0_allSamples.csv"
#file <- "DEGs_ctCut1005_allSamples.csv"
#file <- "DEGs_ctCut1000_allSamples.csv"
#file <- "DEGs_ctCu50_allSamples.csv"
#file <- "DEGs_ahCutUncut0_allSamples.csv"
#file <- "DEGs_ahCut1005_allSamples.csv"
#file <- "DEGs_ahCut1000_allSamples.csv"
#file <- "DEGs_ahCut50_allSamples.csv"

df <- read.csv(paste("DEGs/", file, sep = ""))
df <- df[order(df$log2FoldChange, decreasing = TRUE), ]

# Transformar os genes em uma lista nomeada
gene_list <- df$log2FoldChange
names(gene_list) <- df$X

# Obter conjuntos de genes para Danio rerio (por exemplo, C5 collection)
msigdb_genes <- msigdbr(species = "Danio rerio")

# Realizar a análise GSEA com um número fixo de permutações
gsea_result <- GSEA(
  geneList = gene_list,
  TERM2GENE = msigdb_genes %>% select(gs_name, ensembl_gene),
  pvalueCutoff = 0.05,
  verbose = FALSE
)

# Visualizar resultados
gsea_result_df <- as.data.frame(gsea_result)
#print(head(gsea_result_df))
