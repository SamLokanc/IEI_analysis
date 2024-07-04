# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("TCGAbiolinks")
BiocManager::install("SummarizedExperiment")
BiocManager::install("TCGAutils")
BiocManager::install("maftools")
BiocManager::install("recount")
BiocManager::install("EnhancedVolcano")
BiocManager::install("topGO")
BiocManager::install("sva")
BiocManager::install("EnsDb.Hsapiens.v79")
BiocManager::install("fgsea")
BiocManager::install("goseq")
BiocManager::install("limma")
BiocManager::install("dplyr")
BiocManager::install("ggplot2")
BiocManager::install("data.table")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(TCGAutils)
library(limma)
library(maftools)
library(recount)
library(dplyr)
library(EnhancedVolcano)
library(topGO)
library(sva)
library(EnsDb.Hsapiens.v79)
library(ggplot2)
library(data.table)
library(fgsea)
library(goseq)

#-----GLOBAL VARIABLES-----
gsea_bool <- FALSE
go_bool <- TRUE

p_thresh <- 0.05
fc_thresh <- 1

cond_sub <- "PID"

#set up total samples count file for downstream use
samples.Count <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("tumor", "normal"))

#read list of genes
if(cond_sub != "none"){
  genes <- read.csv(file = "IEI_genes.csv", header = FALSE)
}

#read project codes
proj_codes <- read.csv("proj_codes.csv", header = TRUE)

#----------LOOP START HERE----------
for (i in 1:length(proj_codes$proj)){
  proj <- proj_codes[i, "proj"]
  sample_type <- proj_codes[i, "sample_type"]
  sample_abbr <- proj_codes[i, "sample_abbr"]
  matched_tissue <- proj_codes[i, "matched_tissue"]
  
  dir.create(sprintf("Results/%s", proj), recursive = TRUE)
  dir.create(sprintf("Results/%s/DEA", proj), recursive = TRUE)
  
  if(go_bool){
    dir.create(sprintf("Results/%s/GO", proj), recursive = TRUE)
  }
  if(gsea_bool){
    dir.create(sprintf("Results/%s/GSEA", proj), recursive = TRUE)
  }
  
  #download TARGET data
  query.proj <- GDCquery(
    project = proj,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  samplesDown.proj <- getResults(query.proj, cols = c("cases"))
  samplesDown.proj <- samplesDown.proj[1:1000]
  
  dataSmTP.proj <- TCGAquery_SampleTypes(barcode = samplesDown.proj,
                                           typesample = sample_abbr)
  
  #query findings from TARGET
  query.target.down <- GDCquery(project = proj,
                                data.category = "Transcriptome Profiling",
                                data.type = "Gene Expression Quantification",
                                workflow.type = "STAR - Counts",
                                barcode = c(dataSmTP.proj))
  
  
  #fix any duplicate downloaded data
  query.target.down.2=query.target.down
  tmp=query.target.down.2$results[[1]]
  tmp=tmp[which(!duplicated(tmp$cases)),]
  query.target.down.2$results[[1]]=tmp
  
  #download query
  GDCdownload(query.target.down.2)
  
  #Prepare to preprocess, normalize and filter TARGET data to correct the inputted data
  data <- GDCprepare(query.target.down.2, summarizedExperiment = TRUE)
  
  dataPrep_target <- TCGAanalyze_Preprocessing(object = data,
                                               cor.cut = 0.6)
  
  dataNorm_target <- TCGAanalyze_Normalization(tabDF = dataPrep_target,
                                               geneInfo = geneInfoHT,
                                               method = "gcContent")
  
  dataFilt_target <- TCGAanalyze_Filtering(tabDF = dataNorm_target,
                                           method = "quantile",
                                           qnt.cut = 0.25)
  
  proj.batchCorr.target <- as.data.frame(dataFilt_target)
  
  #download normal tissue samples from the GTEx database with recount2
  proj.recount.gtex <- TCGAquery_recount2(project = "GTEX", tissue = matched_tissue)
  
  #scale gtex data
  tissue <- paste0("GTEX_", gsub(" ","_",matched_tissue))
  proj.eset.gtex <- assays(scale_counts(proj.recount.gtex[[tissue]], round = TRUE))$counts
  
  #quality control
  rse_scaled <- scale_counts(proj.recount.gtex[[tissue]], round = TRUE)
  summary(colSums(assays(rse_scaled)$counts)) / 1e6
  
  #renaming data for downstream use
  proj.eset.target <- proj.batchCorr.target
  
  #removing spaces
  rownames(proj.eset.gtex) <- gsub("\\..*", "", rownames(proj.eset.gtex))
  rownames(proj.eset.target) <- gsub("\\..*", "", rownames(proj.eset.target))
  
  #create summarized data frome of target and gtex data
  dataPrep.proj <- merge(as.data.frame(proj.eset.gtex), 
                           as.data.frame(proj.eset.target), 
                           by = 0, 
                           all = TRUE)
  
  #batch correction between gtex and target
  AnnotationCounts <- matrix(0, ncol(dataPrep.proj), 2)
  colnames(AnnotationCounts) <- c("Condition", "Batch")
  rownames(AnnotationCounts) <- colnames(dataPrep.proj)
  
  AnnotationCounts <- as.data.frame(AnnotationCounts)
  AnnotationCounts$Condition <- colnames(dataPrep.proj)
  
  AnnotationCounts[colnames(proj.eset.gtex), "Batch"] <- "proj.eset.gtex"
  AnnotationCounts[colnames(proj.eset.target), "Batch"] <- "proj.eset.target"
  
  AnnotationCounts <- AnnotationCounts[-1, ]
  
  #remove NA values
  dataPrep.proj <- na.omit(dataPrep.proj)
  
  #make row indices row.names 
  rownames(dataPrep.proj) <- dataPrep.proj$Row.names
  
  #remove row.names column
  dataPrep.proj <- dataPrep.proj[, -which(names(dataPrep.proj) == "Row.names")]
  
  #batch correction
  countsCorrected <- TCGAbatch_Correction(tabDF = dataPrep.proj,
                                          UnpublishedData = TRUE,
                                          AnnotationDF = AnnotationCounts)
  
  #change negative values into zero values for downstream use
  zero.countsCorrected <- countsCorrected
  zero.countsCorrected[zero.countsCorrected < 0] <- 0
  
  #normalize + filter
  dataNorm.proj <- TCGAanalyze_Normalization(tabDF = zero.countsCorrected,
                                               geneInfo = geneInfoHT,
                                               method = "gcContent")
  
  dataFilt.proj <- TCGAanalyze_Filtering(tabDF = dataNorm.proj,
                                           method = "quantile",
                                           qnt.cut = 0.25)
  
  #DEA
  DEG.proj <- TCGAanalyze_DEA(mat1 = dataFilt.proj[,colnames(proj.eset.gtex)],
                                mat2 = dataFilt.proj[,colnames(proj.eset.target)],
                                metadata = FALSE,
                                pipeline = "limma",
                                voom = TRUE,
                                Cond1type = "Normal",
                                Cond2type = "Tumor",
                                method = "glmLRT")
  #change symbol names
  proj.conversion.table <- ensembldb::select(EnsDb.Hsapiens.v79, keys = rownames(DEG.proj), keytype = "GENEID", columns = c("SYMBOL", "GENEID"))
  proj.conversion.inter.DEG <- intersect(proj.conversion.table$GENEID, rownames(DEG.proj))
  
  proj.conversion.table2 <- proj.conversion.table[which(proj.conversion.table$GENEID %in% proj.conversion.inter.DEG),]
  rownames(proj.conversion.table2) <- proj.conversion.table2$GENEID
  
  #convert all genes to their names
  DEG.proj.hgnc <- DEG.proj[proj.conversion.inter.DEG,]
  DEGs.proj.hgnc <- merge(DEG.proj.hgnc, proj.conversion.table2, by = 0)
  rownames(DEGs.proj.hgnc) <- DEGs.proj.hgnc$Row.names
  DEGs.proj.hgnc$Row.names <- NULL
  DEGs.proj.hgnc <- DEGs.proj.hgnc %>% dplyr::filter(!duplicated(SYMBOL))
  rownames(DEGs.proj.hgnc) <- DEGs.proj.hgnc$SYMBOL
  
  #remove unnecessary columns for downstream analysis
  DEGs.proj.hgnc <- subset(DEGs.proj.hgnc, select = -c(AveExpr, t, P.Value, B, GENEID))
  names(DEGs.proj.hgnc)[3] <- "Gene.symbol"
  DEGs.proj.hgnc <- DEGs.proj.hgnc %>% dplyr::select(Gene.symbol, everything())
  
  #write all DEGs to csv file
  file_name = sprintf("Results/%s/DEA/all_DEGs_%s_%s_%s.csv", proj, proj, sample_abbr, tissue)
  write.csv(DEGs.proj.hgnc, file_name)
  
  #subset to select desired genes
  if(cond_sub != "none"){
    genes.DEGs.proj.hgnc <- subset(DEGs.proj.hgnc, Gene.symbol %in% genes$V1)
  } else {
    genes.DEGs.proj.hgnc <- DEGs.proj.hgnc
  }
  #write subset DEGs to csv
  file_name = sprintf("Results/%s/DEA/subset_DEGs_%s_%s_%s.csv", proj, proj, sample_abbr, tissue)
  write.csv(genes.DEGs.proj.hgnc, file_name)
  
  #set file name
  file_name = sprintf("Results/%s/DEA/Volcano_%s_%s_%s.png", proj, proj, sample_abbr, tissue)
  
  #visualize DEA results
  volcano <- EnhancedVolcano(genes.DEGs.proj.hgnc,
                  lab = rownames(genes.DEGs.proj.hgnc),
                  x = 'logFC',
                  y = 'adj.P.Val',
                  ylab = bquote(~-Log[10]~ "(FDR corrected-P values)"),
                  title = sprintf("%s:%s vs. %s",proj,sample_type,tissue),
                  subtitle = sprintf("Differential expression using limma-voom; showing %s related genes", cond_sub),
                  caption = sprintf("FC cutoff, %f; p-value cutoff, %f", fc_thresh, p_thresh),
                  pCutoff = p_thresh,
                  FCcutoff = fc_thresh,
                  pointSize = 6.0,
                  labSize = 6.0,
                  legendPosition = 'none')
  ggsave(file_name, plot = volcano, height = 8, width = 12)
  
  if(go_bool){
    #GO
    #creating gene vector
    genes.DEGs.proj.hgnc05 <- subset(genes.DEGs.proj.hgnc, adj.P.Val < p_thresh)
    
    #filter degs based on adj pval and logfc
    deg <- genes.DEGs.proj.hgnc05[!is.na(genes.DEGs.proj.hgnc05$adj.P.Val) & 
                                    genes.DEGs.proj.hgnc05$adj.P.Val <= p_thresh &
                                    abs(genes.DEGs.proj.hgnc05$logFC) >= fc_thresh, ]
    
    #get gene lists
    assayed_genes <- rownames(genes.DEGs.proj.hgnc)
    de_genes <- rownames(deg)
    gene_vector <- as.integer(assayed_genes %in% de_genes)
    names(gene_vector) <- assayed_genes
    
    #GO object
    pwf <- nullp(gene_vector, "hg19", "geneSymbol")
    go_wall <- goseq(pwf, "hg19", "geneSymbol")
    
    #save full GO analysis
    file_name <- sprintf("Results/%s/GO/file_%s_%s_%s.csv", proj, proj, sample_abbr, tissue)
    write.csv(go_wall, file_name, row.names = FALSE)
    
    #filter the top 50 GO terms based on p-value
    top_go_terms <- go_wall %>%
      slice_min(order = over_represented_pvalue, n = 10) %>%
      bind_rows(go_wall %>%
                  slice_min(order = under_represented_pvalue, n = 10)) %>%
      mutate(hitsPerc = numDEInCat * 100 / numInCat)
    
    #save top 50 GO terms to csv file for easier access later
    file_name <- sprintf("Results/%s/GO/top_file_%s_%s_%s.csv", proj, proj, sample_abbr, tissue)
    write.csv(top_go_terms, file_name, row.names = FALSE)
    
    #plot
    #file_name <- sprintf("Results/%s/GO/diagram_%s_%s_%s.png", proj, proj, sample_abbr, tissue)
    #
    #go <- go_wall %>%
    #  slice_min(order_by = over_represented_pvalue, n = 10) %>%
    #  mutate(hitsPerc = numDEInCat * 100 / numInCat) %>%
    #  ggplot(aes(x = reorder(term, -over_represented_pvalue),
    #             y = hitsPerc,
    #             fill = over_represented_pvalue)) +
    #  geom_col() +
    #  coord_flip() + # Flip coordinates for horizontal bars
    #  expand_limits(y = 0) +
    #  labs(x = "GO term", y = "Hits (%)", fill = "p value") +
    #  theme_classic()
    #
    #ggsave(file_name, plot = go, height = 8, width = 8)
  }
  
  #final sample count updated
  proj.samples.Count <- data.frame(ncol(proj.eset.target), ncol(proj.eset.gtex))
  names(proj.samples.Count) <- c("tumor", "normal")
  proj_str <- paste0(proj,":",sample_abbr, " vs. ", tissue)
  samples.Count[nrow(samples.Count) + 1, ] <- proj.samples.Count
  rownames(samples.Count)[nrow(samples.Count)] <- proj_str
  write.csv(samples.Count, "Results/samples.Count.csv")
  
  #clear data except for gene panel and sample count
  rm(list = setdiff(ls(), c("samples.Count", 
                            "proj_codes", 
                            "genes", 
                            "gsea_bool", 
                            "go_bool", 
                            "p_thresh", 
                            "fc_thresh", 
                            "cond_sub")))
}
