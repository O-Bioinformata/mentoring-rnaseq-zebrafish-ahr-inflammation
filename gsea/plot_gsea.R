################################################################################





################################################################################

gsea_result <- load("~/GitHub/mentoring-rnaseq-zebrafish-ahr-inflammation/gsea/objects_gsea/DEGs_ctCutUncut0_allSamples.RData")

gsea_result$ID


terms_gsea <- c("GOBP_LEUKOCYTE_PROLIFERATION", "GOBP_T_CELL_PROLIFERATION",
                "GOBP_NEURAL_CREST_CELL_MIGRATION", "GOBP_NEURAL_CREST_CELL_DIFFERENTIATION",
                "GOBP_AMP_METABOLIC_PROCESS", "TSENG_IRS1_TARGETS_DN",
                "GSE9316_IL6_KO_VS_IFNG_KO_INVIVO_EXPANDED_CD4_TCELL_DN",
                "HU_FETAL_RETINA_AMACRINE",
                "GOBP_PURINE_NUCLEOSIDE_MONOPHOSPHATE_BIOSYNTHETIC_PROCESS",
                "REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_ADDITIONAL_CELL_CYCLE_GENES_WHOSE_EXACT_ROLE_IN_THE_P53_PATHWAY_REMAIN_UNCERTAIN")

# Plotar resultados significativos
library(enrichplot)
#dotplot(gsea_result) + ggtitle("GSEA Dotplot")

# Exibir as 15 principais categorias com um ajuste de cor para o p-value ajustado
dotplot(gsea_result, showCategory = terms_gsea, color = "p.adjust") #+ ggtitle("ctCutUncut0")

################################################################################

load(file = "objects_gsea/DEGs_ctCut1005_allSamples.RData")

gsea_result$ID


terms_gsea <- c("GOBP_LEUKOCYTE_PROLIFERATION", "GOBP_T_CELL_PROLIFERATION",
                "GOBP_NEURAL_CREST_CELL_MIGRATION", "GOBP_NEURAL_CREST_CELL_DIFFERENTIATION",
                "GOBP_AMP_METABOLIC_PROCESS", "TSENG_IRS1_TARGETS_DN",
                "GSE9316_IL6_KO_VS_IFNG_KO_INVIVO_EXPANDED_CD4_TCELL_DN",
                "HU_FETAL_RETINA_AMACRINE",
                "GOBP_PURINE_NUCLEOSIDE_MONOPHOSPHATE_BIOSYNTHETIC_PROCESS",
                "REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_ADDITIONAL_CELL_CYCLE_GENES_WHOSE_EXACT_ROLE_IN_THE_P53_PATHWAY_REMAIN_UNCERTAIN")

# Plotar resultados significativos

dotplot(gsea_result, showCategory = terms_gsea, color = "p.adjust") #+ ggtitle("ctCutUncut0")
