deseqFunc <- function(data, condition, name)
{
  # Runs DESeq differential gene analyses on a count matrix
  #
  # Args:
  #   data: A data frame with rows as genes, samples as columns and counts 
  #   condition: A vector with the same length as the ncol(data) 
  # Returns:
  #   Plots
  #     _plot1          dispersion plot
  #     _plot2          MA plot (volcano plot)
  #     _plot3          histogram of p-values
  #   Data
  #     _output.txt     genes, pvals, log-fold change calculated from DESeq
  #     _sigOutput.txt  only genes with an adjusted pvalue < 0.1
  ### data looks like this
  #       before  before  after  after
  # gene
  # gene
  
  ### condition looks like this
  # before  before  after after
  
  dispersionName <- paste(name,"_plot1",sep="")
  maName <- paste(name,"_plot2",sep="")
  pvalName <- paste(name,"_plot3",sep="")
  library('DESeq')
  cds <- newCountDataSet(data, condition)
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
    
  png(filename=paste(dispersionName, ".png", sep=""),width=800, height=650, units="px", res=100)
  plotDispEsts(cds)
  dev.off()
  res <- nbinomTest(cds, unique(condition)[1], unique(condition)[2])
  
  png(filename=paste(maName, ".png", sep=""),width=800, height=650, units="px", res=100)
  plotMA(res)
  dev.off()
  
  png(filename=paste(pvalName, ".png", sep=""),width=800, height=650, units="px", res=100)
  hist(res$pval, breaks=200, col="skyblue", border="slateblue", 
       main="P-Value in bins of 100")
  dev.off()
  
  res <- res[order(res$padj), ]
  # Write out all genes
  write.table(res, file=paste(name,"_output.txt",sep=""), quote=FALSE,sep="\t")
  
  # Write out only significant genes
  write.table(res[which(res$padj < 0.1), ], file=paste(name,"_sigOutput.txt",sep=""), quote=FALSE)
}
