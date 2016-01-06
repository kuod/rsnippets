#run DESEQ
run_deseq <- function(samples, dir, formula) {
  # only use genes with a count sum greater than 0
  print(dir)
  sums <- apply(count_df[,samples],1,sum)
  minReadBool <- sums > 1
  mainDir <- getwd()
  dir.create(file.path(mainDir, dir))
  dds <- DESeqDataSetFromMatrix(countData=count_df[minReadBool, samples], colData=design_mat[samples, ], formula )
  colnames(dds) <- colnames(count_df[samples])
  (condition <- factor(colData(dds)$CP))
  dds <- DESeq(dds)
  res <- results(dds)
  order_res <- res[order(res$padj), ]
  annotation <- as.data.frame(getBM(mart=mart, attributes=c("external_gene_name","ensembl_gene_id","description","chromosome_name","start_position",             ↪"end_position"), values=row.names(order_res)))
  row.names(annotation) <- annotation$ensembl_gene_id
  rld <- rlogTransformation(dds)
  #### Plots ###
  #QC
  filePath <- file.path(mainDir, dir,'/')
  (mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])
  sampleDists <- as.matrix(dist(t(assay(rld))))

  #Euclidean distance matrix of similarity between samples
  pdf(paste(filePath, paste(dir,"qc_dist_samples.pdf",sep='_'), sep=''))
  heatmap.2(as.matrix(sampleDists), key=F, trace="none",
            col=colorpanel(100, "black", "white"),
            ColSideColors=mycols[condition], RowSideColors=mycols[condition],
            margin=c(10, 10), main="Sample Distance Matrix")
  dev.off()
  #PCA using top 500 genes
  pdf(paste(filePath, paste(dir,"qc_pca_samples.pdf",sep='_'), sep=''))
  rld_pca(rld, colors=mycols, intgroup="CP", ntop=500)
  dev.off()
  #Histogram of pvalues
  pdf(paste(filePath, paste(dir,'pvalhist.pdf', sep='_'), sep=''))
  hist(res$padj, breaks=40)
  dev.off()
  #standard maplot
  pdf(paste(filePath, paste(dir,'maplot.pdf', sep='_'), sep=''))
  plotMA(res)
  dev.off()
  #dispersion plot
  pdf(paste(filePath, paste(dir, 'disp.pdf', sep='_'), sep=''))
  plotDispEsts(dds)
  dev.off()

  output_df <- merge(as.data.frame(order_res), as.data.frame(counts(dds, normalized=TRUE)), by=0)
  row.names(output_df) <- output_df$Row.names
  output_df <- merge(as.data.frame(output_df), annotation, by=0)
  row.names(output_df) <- output_df$Row.names
  output_df <- subset(output_df, select= -c(Row.names, ensembl_gene_id))
  output_df <- subset(output_df, select= -c(Row.names))
  output_df <- output_df[order(output_df["padj"]),]
  write.table(as.data.frame(output_df), file=paste(filePath, paste(dir,paste('results','csv',sep='.'),sep="_"), sep=""),quote=FALSE, row.names=TRUE,  sep='\t')
  return(order_res)
}
