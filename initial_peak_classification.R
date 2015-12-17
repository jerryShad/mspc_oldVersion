################# Initial peak classification #################

@ author: Julaiti Shayiding

@param: source_peak all peak from the source gene replicate
@param: tau_S threshold value for stringent peak
@param: tau_W threshold value for weak peak

initialPeakClassification<-function(source_peak, tau_S, tau_W, outputDir, ...){

  #pvalue =-1xlog10(peak_chrom$score)
  analysisResult<-list()
  for(i in 1: length(souce_peak)){
    peak<-Replicate[i,]
    if(peak$pvalue<= tau_S){
      message(" peak is stringent", peaks)
      analysisResult(peaks)
    }
    else if(peak_chrom$pvalue <= tau_W){
      message(" peak is weak", peaks)
      analysisResults(peaks)
    }
  }
  return(analysisResult)
}
