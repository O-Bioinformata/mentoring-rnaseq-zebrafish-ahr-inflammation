################################################################################





################################################################################
# Instalar e carregar pacotes necessários
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Dr.eg.db")
#BiocManager::install("msigdbr")

library(clusterProfiler)
library(org.Dr.eg.db)
library(msigdbr)
library(dplyr)

# Carregar os dados de genes diferencialmente expressos
# Suponha que seu dataframe se chama 'df' com colunas 'gene' e 'log2FoldChange'
# Ordene os genes pela estatística que deseja usar (e.g., log2FoldChange)

#file <- "DEGs_ctCutUncut0_allSamples.csv"
file <- "DEGs_ctCut1005_allSamples.csv"
#file <- "DEGs_ctCut1000_allSamples.csv"
#file <- "DEGs_ctCu50_allSamples.csv"
#file <- "DEGs_ahCutUncut0_allSamples.csv"
#file <- "DEGs_ahCut1005_allSamples.csv"
#file <- "DEGs_ahCut1000_allSamples.csv"
#file <- "DEGs_ahCut50_allSamples.csv"


df <- read.csv(paste("DEGs/", file, sep = ""))
df <- df[order(df$log2FoldChange, decreasing = TRUE), ]


#df <- genes[order(genes$log2FoldChange, decreasing = TRUE), ]

# Transformar os genes em uma lista nomeada
gene_list <- df$log2FoldChange
names(gene_list) <- df$X

# Obter conjuntos de genes para Danio rerio (por exemplo, C5 collection)
msigdb_genes <- msigdbr(species = "Danio rerio")

# Realizar a análise GSEA
gsea_result <- GSEA(
  geneList = gene_list,
  TERM2GENE = msigdb_genes %>% select(gs_name, ensembl_gene),
  pvalueCutoff = 0.05,
  verbose = FALSE
)

# Visualizar resultados
gsea_result_df <- as.data.frame(gsea_result)
#print(head(gsea_result_df))

write.csv(gsea_result_df, file = paste("results_gsea/GSEA_", file, sep = ""))
save(gsea_result, file = paste("objects_gsea/", gsub(".csv","",file),".RData", sep = ""))





# Plotar resultados significativos
library(enrichplot)
dotplot(gsea_result) + ggtitle("GSEA Dotplot")

# Exibir as 15 principais categorias com um ajuste de cor para o p-value ajustado
dotplot(gsea_result, showCategory = 15, color = "p.adjust") + ggtitle("Top 15 GSEA Dotplot")

gseaplot2(gsea_result, geneSetID = 1, title = "GSEA Plot for Top Pathway")





cnetplot(gsea_result, categorySize = "pvalue", foldChange = gene_list)

heatplot(gsea_result, foldChange = gene_list)

ridgeplot(gsea_result)

upsetplot(gsea_result)

gseatable(gsea_result)














