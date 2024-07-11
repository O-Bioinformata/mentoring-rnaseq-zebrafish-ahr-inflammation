################################################################################
#
#
#
#
#
################################################################################

DEGs_allSamples <- function(coldata, i, data){
  require(DESeq2)
  # -- coldata
  if(i == "ahCutUncut0_allSamples" | i == "ctCutUncut0_allSamples"){
    coldata_samples <- coldata[coldata$names %in% colnames(data$counts),
                               c("names","Trat_02")]
    
  } else if(i == "ctCut1005_allSamples" | i == "ahCut50_allSamples" |
            i == "ahCut1000_allSamples" | i == "ahCut1005_allSamples" |
            i == "ctCut50_allSamples" | i == "ctCut1000_allSamples"){
    coldata_samples <- coldata[coldata$names %in% colnames(data$counts),
                               c("names","Trat_03")]
  }else{
    cat("verifique a comparação 1")
  }
  
  # -- preparacao do objeto
  if(i == "ahCutUncut0_allSamples" | i == "ctCutUncut0_allSamples"){
    ddsTxi <- DESeqDataSetFromTximport(txi = data,
                                       colData = coldata_samples,
                                       design = ~ Trat_02)
    
  } else if(i == "ctCut1005_allSamples" | i == "ahCut50_allSamples" |
            i == "ahCut1000_allSamples" | i == "ahCut1005_allSamples" |
            i == "ctCut50_allSamples" | i == "ctCut1000_allSamples"){
    ddsTxi <- DESeqDataSetFromTximport(txi = data,
                                       colData = coldata_samples,
                                       design = ~ Trat_03)
  }else{
    cat("verifique a comparação 2")
  }
  
  # -- pre-filtragem
  keep <- rowSums(counts(ddsTxi)) >= 10
  dds <- ddsTxi[keep,]
  
  # -- setting the reference
  if(i == "ahCutUncut0_allSamples" | i == "ctCutUncut0_allSamples"){
    dds$Trat_02 <- relevel(dds$Trat_02, ref = "uncut")
    
  } else if(i == "ahCut50_allSamples" | i == "ahCut1000_allSamples" |
            i == "ctCut50_allSamples" | i == "ctCut1000_allSamples"){
    dds$Trat_03 <- relevel(dds$Trat_03, ref = "0")
    
  } else if(i == "ctCut1005_allSamples" | i == "ahCut1005_allSamples"){
    dds$Trat_03 <- relevel(dds$Trat_03, ref = "5")
  }else{
    cat("verifique a comparação 3")
  }
  
  # -- differential expression
  if(i == "ahCutUncut0_allSamples" | i == "ctCutUncut0_allSamples"){
    dds <- DESeq(dds)
    res_padj05_lfc0 <- results(dds, contrast = c("Trat_02","cut","uncut"), alpha = 0.05)
    
  } else if(i == "ahCut50_allSamples" | i == "ctCut50_allSamples"){
    dds <- DESeq(dds)
    res_padj05_lfc0 <- results(dds, contrast = c("Trat_03","5","0"), alpha = 0.05)
    
  }else if(i == "ahCut1000_allSamples" | i == "ctCut1000_allSamples"){
    dds <- DESeq(dds)
    res_padj05_lfc0 <- results(dds, contrast = c("Trat_03","100","0"), alpha = 0.05)
    
  } else if(i == "ctCut1005_allSamples"){
    dds <- DESeq(dds)
    res_padj05_lfc0 <- results(dds, contrast = c("Trat_03","100","5"), alpha = 0.05)
    
  }else if(i == "ahCut1005_allSamples"){
    dds <- DESeq(dds)
    res_padj05_lfc0 <- results(dds, contrast = c("Trat_03","100","5"), alpha = 0.05)
    
  }else{
    cat("verifique a comparação 4")
  }
  
  # -- save
  if(i == "ahCutUncut0_allSamples"){
    write.csv(res_padj05_lfc0,file = "DEGs_ahCutUncut0_allSamples.csv")
    
  } else if(i == "ahCut50_allSamples"){
    write.csv(res_padj05_lfc0, file = "DEGs_ahCut50_allSamples.csv")
    
  } else if(i == "ahCut1000_allSamples"){
    write.csv(res_padj05_lfc0, file = "DEGS_ahCut1000_allSamples.csv")
    
  } else if(i == "ctCut1005_allSamples"){
    write.csv(res_padj05_lfc0, file = "DEGs_ctCut1005_allSamples.csv")
    
  }else if(i == "ahCut1005_allSamples"){
    write.csv(res_padj05_lfc0, file = "DEGs_ahCut1005_allSamples.csv")
  
  }else if(i == "ctCutUncut0_allSamples"){
    write.csv(res_padj05_lfc0, file = "DEGs_ctCutUncut0_allSamples.csv")
  
  }else if(i == "ctCut50_allSamples"){
    write.csv(res_padj05_lfc0, file = "DEGs_ctCut50_allSamples.csv")
  
  }else if(i == "ctCut1000_allSamples"){
    write.csv(res_padj05_lfc0, file = "DEGs_ctCut1000_allSamples.csv")
    
  } else{
    cat("verifique a comparacao 5")
  }
}

################################################################################
# -- Anotacao -- 
load("EnsDbAnnotation_20230519_atual.RData")
data_protCod <- EnsDbAnnotation[EnsDbAnnotation$gene_biotype=="protein_coding",]
head(data_protCod)
table(data_protCod$gene_biotype) # 30153

# -- colldata -- 
coldata <- readxl::read_xlsx("coldata.xlsx")
str(coldata)
coldata <- as.data.frame(coldata)

for (i in 2:ncol(coldata)){
  coldata[,i] <- as.factor(coldata[,i])
}

str(coldata)

# -- preparando objeto -- ######################################################
load("matrix_salmon_tximport_20230519.RData")

## -- matrizes com allSamples 
data_allSamples <- mat_gse; rm(mat_gse)

head(rownames(data_allSamples$counts)) # sem num da versao
head(rownames(data_protCod)) # sem num da versao

nrow(data_allSamples$counts)

idx <- match(rownames(data_protCod), rownames(data_allSamples$counts))
table(is.na(idx))
idx <- idx[!is.na(idx)]

colnames(data_allSamples$counts)

ahCutUncut0_allSamples <- grep("^AH.*_0_.", colnames(data_allSamples$counts)) # ahCutUncut0_allSamples

ctCutUncut0_allSamples <- grep("^CT.*_0_.", colnames(data_allSamples$counts)) # ctCutUncut0_allSamples

ahCut50_allSamples <- grep("AH_Cut_._.", colnames(data_allSamples$counts)) # ahCut50_allSamples

ctCut50_allSamples <- grep("CT_Cut_._.", colnames(data_allSamples$counts)) # ctCut50_allSamples

ahCut1000_allSamples <- grep("AH_Cut_.*0_.", colnames(data_allSamples$counts)) # ahCut1000_allSamples

ctCut1000_allSamples <- grep("CT_Cut_.*0_.", colnames(data_allSamples$counts)) # ctCut1000_allSamples

ctCut1005_allSamples <- c(grep("CT_Cut_100_.", colnames(data_allSamples$counts)),
                          grep("CT_Cut_5_.", colnames(data_allSamples$counts)))# ctCut1005_allSamples

ahCut1005_allSamples <- c(grep("AH_Cut_100_.", colnames(data_allSamples$counts)),
                          grep("AH_Cut_5_.", colnames(data_allSamples$counts)))# ahCut1005_allSamples

# -- expressao diferencial: all Samples -- #####################################

## -- ahCutUncut0_allSamples -- 
data <- data_allSamples
i="ahCutUncut0_allSamples"
data$abundance <- data$abundance[idx,ahCutUncut0_allSamples]
data$counts <- data$counts[idx,ahCutUncut0_allSamples]
data$length <- data$length[idx,ahCutUncut0_allSamples]

DEGs_allSamples(coldata=coldata, i= i, data= data)

## -- ctCutUncut0_allSamples -- 
data <- data_allSamples
i="ctCutUncut0_allSamples"
data$abundance <- data$abundance[idx,ahCutUncut0_allSamples]
data$counts <- data$counts[idx,ahCutUncut0_allSamples]
data$length <- data$length[idx,ahCutUncut0_allSamples]

DEGs_allSamples(coldata=coldata, i= i, data= data)

## -- ahCut50_allSamples
data <- data_allSamples
i="ahCut50_allSamples"
data$abundance <- data$abundance[idx,ahCut50_allSamples]
data$counts <- data$counts[idx,ahCut50_allSamples]
data$length <- data$length[idx,ahCut50_allSamples]

DEGs_allSamples(coldata=coldata, i= i, data= data)

## -- ctCut50_allSamples
data <- data_allSamples
i="ctCut50_allSamples"
data$abundance <- data$abundance[idx,ahCut50_allSamples]
data$counts <- data$counts[idx,ahCut50_allSamples]
data$length <- data$length[idx,ahCut50_allSamples]

DEGs_allSamples(coldata=coldata, i= i, data= data)

## -- ahCut1000_allSamples
data <- data_allSamples
i="ahCut1000_allSamples"
data$abundance <- data$abundance[idx,ahCut1000_allSamples]
data$counts <- data$counts[idx,ahCut1000_allSamples]
data$length <- data$length[idx,ahCut1000_allSamples]

DEGs_allSamples(coldata=coldata, i= i, data= data)

## -- ctCut1000_allSamples
data <- data_allSamples
i="ctCut1000_allSamples"
data$abundance <- data$abundance[idx,ahCut1000_allSamples]
data$counts <- data$counts[idx,ahCut1000_allSamples]
data$length <- data$length[idx,ahCut1000_allSamples]

DEGs_allSamples(coldata=coldata, i= i, data= data)

## -- ctCut1005_allSamples
data <- data_allSamples
i="ctCut1005_allSamples"
data$abundance <- data$abundance[idx,ctCut1005_allSamples]
data$counts <- data$counts[idx,ctCut1005_allSamples]
data$length <- data$length[idx,ctCut1005_allSamples]

DEGs_allSamples(coldata=coldata, i= i, data= data)

## -- ahCut1005_allSamples
data <- data_allSamples
i="ahCut1005_allSamples"
data$abundance <- data$abundance[idx,ahCut1005_allSamples]
data$counts <- data$counts[idx,ahCut1005_allSamples]
data$length <- data$length[idx,ahCut1005_allSamples]

DEGs_allSamples(coldata=coldata, i= i, data= data)
