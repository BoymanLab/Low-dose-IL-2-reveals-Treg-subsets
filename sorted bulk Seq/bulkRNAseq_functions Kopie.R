####Introduction####
#Here, we define various functions for bulk RNA Seq data analysis. The functions expand the edgeR package and add Gene Set Enrichment Analysis to the pipeline.

####edgeR utility####


prepDGE <- function(cd = count_data, md = meta, gmd = gene_meta, original_genes = "ensembl_transcript_id", new_genes = "hgnc_symbol"){
  count_data2 <- merge(cd, gmd[c(all_of(original_genes), all_of(new_genes))], by.x = 0, by.y = all_of(original_genes)) %>% dplyr::select(-Row.names)
  count_data2 <- count_data2 %>% 
    group_by(across(all_of(new_genes))) %>% 
    summarize_all(sum) %>% 
    na.omit() %>%
    arrange(get(new_genes)) %>%
    as.data.frame() %>%
    dplyr::filter(get(new_genes) != "") %>%
    column_to_rownames(all_of(new_genes))
  
  gmd <- gmd %>% 
    dplyr::filter(get(new_genes) %in% rownames(count_data2)) %>%
    .[!duplicated(.[all_of(new_genes)]),] %>%
    arrange(get(new_genes))
  
  y <- DGEList(counts = count_data2, group = md$group, genes = gmd)
  
  return(y)
}

plot_MDS <- function(y = y, md = meta, color_by = "tp", sample_id = "sample_id"){
  #plotMDS: "Distances on an MDS plot of a DGEList object correspond to leading log-fold-change between each pair of samples." -edgeR manual
  mds <- plotMDS(y, col = as.numeric(md[,color_by]), plot = FALSE)
  mdsdf <- data.frame(mds[c("x", "y")], "sample_id" = md[,sample_id], "color_by" = md[,color_by])
  mds.p <- ggplot(mdsdf, aes(x = x, y = y))+
    geom_text(aes(label = sample_id, col = color_by), key_glyph = "rect")+
    labs(x = paste0("Leading logFC dim1 (", round(mds$var.explained[1]*100), "%)"), 
         y = paste0("Leading logFC dim2 (", round(mds$var.explained[2]*100), "%)"))+
    scale_color_manual(values = c("firebrick1", "green3", "blue3", "darkviolet"))+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5))
  
  return(mds.p)
}

plot_volcano <- function(toptags = toptags, ylims = c(0,6), xlims = "free", gene_nomen = "hgnc_symbol", labels = NULL){
  toptags$delabel <- NA
  if(!is.null(labels)){
    toptags[toptags[,all_of(gene_nomen)] %in% all_of(labels), "delabel"] <- toptags[toptags[,all_of(gene_nomen)] %in% all_of(labels), all_of(gene_nomen)]
  }else{
    toptags[toptags$FDR < 0.05 & abs(toptags$logFC) > 1.0, "delabel"] <- toptags[toptags$FDR < 0.05 & abs(toptags$logFC) > 1.0, all_of(gene_nomen)]
    toptags[duplicated(toptags$delabel), "delabel"] <- NA #this should not be necessary if the same genes are used to summarize and to plot    
  }
  toptags$col <- "NO"
  toptags$col[toptags$FDR < 0.05 & toptags$logFC > 1.0] <- "UP"
  toptags$col[toptags$FDR < 0.05 & toptags$logFC < -1.0] <- "DOWN"
  
  p1 <- ggplot(toptags, aes(x = logFC, y = -log10(FDR)))+
    geom_point(aes(col = col, size = logCPM))+
    scale_color_manual(breaks = c("DOWN", "NO", "UP"), values = c("blue", "grey", "red"))+
    geom_text_repel(aes(label = delabel), min.segment.length = unit(0, 'lines'), max.overlaps = 10000, max.time = 5, max.iter = 100000)+
    labs(x = "Log2 FC", y = "-log10 FDR", col = "Regulation", size = "Mean logCPM")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_y_continuous(limits = c(0, 6))  
  
  if(length(xlims) == 2){
    p1 <- p1 + scale_x_continuous(limits = xlims)
  }else if(length(xlims) == 1 & xlims == "symmetrical"){
    p1 <- p1 + scale_x_continuous(limits = c(-max(abs(toptags$logFC)), max(abs(toptags$logFC))))
  }
  # if(xlims == "symmetrical"){
  #   p1 <- p1 + scale_x_continuous(limits = c(-max(abs(toptags$logFC)), max(abs(toptags$logFC))))
  # }else if(xlims != "free"){
  #   p1 <- p1 + scale_x_continuous(limits = xlims)
  # }
  
  return(p1)
}

#### GSEA ####
#Here, we define functions for Gene Set Enrichment Analysis
run_gsea <- function(toptags = toptags, gsl = genesetlib, n_perm = 1000, seed = 42){
  toptags$Fdir <- toptags$`F` * sign(toptags$logFC)
  
  res2 <- toptags %>% 
    dplyr::select(hgnc_symbol, Fdir) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(hgnc_symbol) %>% 
    summarize(Fdir = mean(Fdir))
  
  ranks <- deframe(res2)
  set.seed(seed = seed)
  fgseaRes <- fgseaMultilevel(pathways = gsl, stats=ranks, sampleSize = n_perm)
  fgseaResTidy <- fgseaRes %>% 
    as_tibble() %>% 
    arrange(padj) %>%
    mutate("pathway" = factor(pathway, levels = pathway)) %>% 
    mutate("LeadingEdge" = NA)
  
  for(q in 1:nrow(fgseaResTidy)){
    fgseaResTidy$LeadingEdge[q] <- fgseaResTidy$leadingEdge[q] %>% unlist() %>% .[1:6] %>% na.omit() %>% paste(collapse = ", ") %>% unlist()
  }
  fgseaResTidy <- fgseaResTidy %>% dplyr::select(-leadingEdge)
  
  return(fgseaResTidy)
}

gsea_barplot <- function(fgsea_res, top_n = NULL, p_max = 0.05){
  if(!is.null(top_n) & !is.null(p_max)){
    df <- fgsea_res %>% .[1:top_n,] %>% dplyr::filter(padj < p_max)
  }else if(!is.null(p_max)){
    df <- fgsea_res %>% dplyr::filter(padj < p_max)
  }else if(!is.null(top_n)){
    df <- fgsea_res %>% .[1:top_n,]
  }
  ggplot(df, aes(x = pathway, y = -log10(padj)))+
    geom_bar(stat = "identity", fill = "black", col = "black", width = 0.7)+
    theme_classic()+
    coord_flip()+
    geom_text(aes(label = LeadingEdge), hjust = -0.05, col = "black")+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.7)))+
    labs(y = "-log10 FDR")+
    theme(axis.title.y = element_blank())+
    scale_x_discrete(limits = rev)  
}

plotEnrichment2 <- function(pathway, stats, padj,
                            gseaParam=1,
                            ticksSize=0.2) {
  
  rnk <- rank(-stats)
  ord <- order(rnk)
  
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
  statsAdj <- statsAdj / max(abs(statsAdj))
  
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                          returnAllExtremes = TRUE)
  
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))
  
  diff <- (max(tops) - min(bottoms)) / 8
  
  # Getting rid of NOTEs
  x=y=NULL
  g <- ggplot(toPlot, aes(x=x, y=y)) +
    geom_point(color="green", size=0.1) +
    geom_hline(yintercept=max(tops), colour="red", linetype="dashed") +
    geom_hline(yintercept=min(bottoms), colour="red", linetype="dashed") +
    geom_hline(yintercept=0, colour="black") +
    geom_line(color="green") + theme_classic() +
    geom_segment(data=data.frame(x=pathway),
                 mapping=aes(x=x, y=-diff/2,
                             xend=x, yend=diff/2),
                 size=ticksSize) +
    
    theme(panel.border=element_blank(),
          panel.grid.minor=element_blank()) +
    
    labs(x="rank", y="enrichment score")+
    annotate(geom = "text",
             x = length(stats) * 0.8,
             y = max(tops) * 0.8,
             #vjust = 20,
             label = paste0("p = ", signif(padj, 2)))+
    theme(plot.title = element_text(hjust = 0.5))
  
  g
}

plot_gsea <- function(fgsea_res, genesets, toptags = toptags, top_n = NULL, pathways = NULL, plot_title = NULL){
  if(!is.null(top_n)){
    gs_names <- fgsea_res[1:top_n, "pathway"] %>% unlist() %>% unname()
  }else if(!is.null(pathways)){
    gs_names <- pathways
  }else{
    gs_names <- fgsea_res[fgsea_res$padj <= 0.05, "pathway"] %>% unlist() %>% unname()
  }
  
  ranks <- toptags %>% 
    dplyr::mutate("Fdir" = `F` * sign(logFC)) %>%
    dplyr::select("hgnc_symbol", "Fdir") %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(hgnc_symbol) %>% 
    summarize(Fdir = mean(Fdir)) %>%
    deframe()
  
  plotlist <- list()
  for(geneset in gs_names){
    plotlist[[geneset]] <- plotEnrichment2(pathway = genesets[geneset] %>% unlist() %>% unname(), 
                                           stats = ranks, padj = fgsea_res[fgsea_res$pathway == geneset, "padj"])+ labs(title = paste(plot_title, geneset))
  }
  
  return(plotlist)
}

#### Wrapper ####
wrapper <- function(count_data = count_data, md = meta, gmd = gene_meta, gs = genesets, 
                    design = design, contrasts = contrasts, run_GSEA = FALSE,
                    plot_volc = TRUE, plot_GSEA = TRUE, plot_Enrichment = TRUE,
                    write_excel = TRUE, excel_identifier = celltype,
                    volc.x = c(-5.5, 5.5), volc.y = c(0,6)){
  
  output <- ls[] #output list to which we will add during the analysis
  
  y <- prepDGE(cd = count_data, md = md, gmd = gmd)
  
  #filter out genes that are not expressed by at least all samples in one condition
  keep <- filterByExpr(y)
  table(keep)
  y <- y[keep,,keep.lib.sizes=FALSE]
  
  #normalize for lib size
  y <- calcNormFactors(y) 
  
  #extract count data
  count_data <- cpm(y, normalized.lib.sizes=TRUE) %>% as.data.frame()
  rownames(count_data) <- y$genes$hgnc_symbol
  output["count_data"] <- count_data
  
  #plot MDS; this is bugged since a recent update, check possible fixes or workarounds
  # mds.p <- plot_MDS(y, md = meta_ct) + labs(title = paste0("Low-dose IL-2 in SLE: Bulk RNA Seq MDS plot ", celltype), color = "Time point")
  # output["MDS"] <- mds.p
  
  #estimate dispersion
  y <- estimateDisp(y, design)
  
  #fit GLM model
  fit <- glmQLFit(y,design)
  
  #Loop through contrasts to test one by one.
  DGE <- ls[]
  Volcanoplots <- ls[]
  GSEA <- ls[]
  Enrichmentplots <- ls[]
  for(i in 1:length(contrasts)){
    #test GLM model and extract DGE
    QLFTest <- glmQLFTest(fit,contrast = contrasts[[i]])
    toptags <- topTags(QLFTest, n = nrow(QLFTest))$table
    DGE[names(contrasts)[i]] <- toptags
    if(write_excel){write.xlsx(toptags, paste0(path.out, excel_identifier, "_toptags_", "_", names(contrasts)[i], ".xlsx"))}
    
    #use DGE for volcano plot
    volcano_plot <- plot_volcano(toptags = toptags, xlims = volc.x, ylims = volc.y) + labs(title = paste(celltype, names(contrasts)[i]))
    Volcanoplots[names(contrasts)[i]] <- volcano_plot
    if(plot_volc){print(volcano_plot)}
    
    if(run_GSEA){
      #Run GSEA
      fgsea_res <- run_gsea(toptags = toptags, gsl = genesetlib)
      GSEA[names(contrasts)[i]] <- fgsea_res
      if(write_excel){write.xlsx(fgsea_res, paste0(path.out, excel_identifier, "_GSEA", "_", contrast, ".xlsx"))}

      #Use GSEA results for enrichment plots
      enrichment_plots <- plot_gsea(fgsea_res, gsl = genesetlib, toptags = toptags)
      Enrichmentplots[names(contrasts)[i]] <- enrichment_plots
      if(plot_Enrichment){print(enrichment_plots)}
    }
  }
  output["DGE"] <- DGE
  output["Volcanoplots"] <- Volcanoplots
  output["GSEA"] <- GSEA
  output["Enrichmentplots"] <- Enrichmentplots
  
  return(output)
}

boxplots_wrapper <- function(celltype, gois, toptags.unpaired, toptags.paired, norm.cpm){
  
  df <- data.frame(hgnc_symbol = gois[[celltype]])
  for(m in 1:length(toptags.unpaired)){
    toptags.f <- toptags.unpaired[[m]] %>% 
      dplyr::filter(hgnc_symbol %in% gois[[celltype]]) %>% 
      dplyr::select(hgnc_symbol, FDR, entrezgene_id) %>% 
      rename(!!names(toptags.unpaired)[m] := FDR)
    df <- merge(df, toptags.f)
  }
  
  #df <- df %>% dplyr::select(-names(toptags.paired))
  for(m in 1:length(toptags.paired)){
    toptags.f <- toptags.paired[[m]] %>% 
      dplyr::filter(hgnc_symbol %in% gois[[celltype]]) %>% 
      dplyr::select(hgnc_symbol, FDR, entrezgene_id) %>% 
      rename(!!names(toptags.paired)[m] := FDR)
    df <- merge(df, toptags.f)
  }
  
  count_data.ct.f <- norm.cpm[rownames(norm.cpm) %in% gois[[celltype]],] %>% as.data.frame()
  df.f <- merge(df, count_data.ct.f, by.x = "hgnc_symbol", by.y = 0)
  df.f <- df.f %>% tidyr::gather(key = "comparison", value = "FDR", c(names(toptags.unpaired), names(toptags.paired)))
  
  tests <- df.f %>% 
    mutate("p.adj" = signif(FDR, 2)) %>% #round FDR to two significant digits
    mutate("copy" = comparison) %>%
    separate(copy, c("group1", "group2")) %>%
    mutate("xmin" = as.numeric(factor(group1, levels = c("H", "V2", "V3", "V9")))) %>%
    mutate("xmax" = as.numeric(factor(group2, levels = c("H", "V2", "V3", "V9"))))
  
  tests
  tests$y.position <- apply(tests[colnames(count_data.ct.f)], 1, max) #define y-position of p-value
  tests$y.position <- tests$y.position * (1 + 0.11 * (tests$xmax - tests$xmin))
  
  tests$p.adj.signif <- "ns" #add significance levels
  tests[tests$p.adj < 0.05, "p.adj.signif"] <- "*"
  tests[tests$p.adj < 0.01, "p.adj.signif"] <- "**"
  tests[tests$p.adj < 0.001, "p.adj.signif"] <- "***"  
  tests[tests$p.adj < 0.0001, "p.adj.signif"] <- "****"
  tests$p.adj.comb <- tests$p.adj.signif
  tests[tests$p.adj.comb != "ns", "p.adj.comb"] <- tests[tests$p.adj.comb != "ns", "p.adj"] #add column with pvalues for significant genes
  
  df.f <- df.f %>% gather(key = "sample_id", value = "counts", colnames(count_data.ct.f)) #format data frame
  
  meta_ct <- meta %>% dplyr::filter(ct == celltype)
  df.f <- merge(df.f, meta_ct, by = "sample_id") #add meta information (timepoint etc.)
  
  p1 <- ggplot(df.f, aes(x = tp, y = counts))+
    geom_boxplot(aes(col = tp))+ #boxplot
    geom_point(aes(col = tp))+ #dots
    geom_line(aes(group = id), col = "grey")+ #lines
    facet_wrap("hgnc_symbol", scales = "free_y", ncol = round(sqrt(4*length(gois[[celltype]])/3)))+ #one plot per gene
    theme_bw()+ #remove grey background
    scale_color_brewer(palette = "Set1")+ #define colors
    theme(strip.background = element_blank())+ #remove grey background from titles
    stat_pvalue_manual(tests, label = "p.adj.comb", bracket.shorten = 0.05)+ #add pvalues from glmQLFtest
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.11)))+ #make plot bigger to fit in pvalues
    theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5))+ #remove x-axis title
    labs(y = "Normalized Counts per Million", col = "Timepoint", title = paste0("Bulk RNA Seq ", celltype)) #add titles for plot, y-axis and legend
  return(p1)
}
