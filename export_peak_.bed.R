## export the result and genereate output bed files as an output

export_result<-function(peaks, initial_peak_classification, ...){
  opt<-options(stringsAsFactors = F)
  if(missing(prefix))
    prefix<-paste("initial Peak classification-", Sys.Date(), sep="")

  R_j_S <-paste(prefix, ".Stringent.peak.bed", sep="")
  R_j_W <-paste(prefix, ". Weak.peak.bed", sep="")

  write.table(peak, file=R_j_S, sep = "\t", quote = F)
  write.table(peak, file=R_j_W, sep="\t", quote = F

              #' @param out output directory of the bed files for initial peak classification
              out<-data.frame(peak_chrom=chrom, peak_chrom$start=start, peak_chrom$end=end,
                              peak_chrom$name=name, peak_chrom$pvalue=pvalue)

              write.table(data.frame(out, analyze_result$stringentPeak, row.names = NULL), file=R_j_S,sep = "\t", quote = F, row.names = F)
              write.table(data.frame(out, analyze_result$weakPeak, row.names = NULL), file=R_j_W,sep = "\t", quote = F, row.names = F)
}

