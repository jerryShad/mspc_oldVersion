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


##################### improve initial peak classification function as below       #####################
@ reference code for initial peak classification


pickPeak <- function(score, threshold, offset=0, sub=FALSE){
  idx <- which(score > threshold)
  if(length(idx) == 0) return(NULL)
  
  idx.dist <- diff(idx)
  change.at <- which(idx.dist > 1)
  idx.start <- c(1, change.at +1)
  idx.end <- c(change.at, length(idx))
  broadPeak <- IRanges::Views(score, idx[idx.start], idx[idx.end])
  peak <- IRanges::viewWhichMaxs(broadPeak)
  ## if the peak is flat we pick the centre
  for(i in 1:length(peak)){
    runs <- as(broadPeak[[i]][(peak[i] - start(broadPeak)[i] + 1):width(broadPeak)[i]], "Rle")
    peakWidth <- runLength(runs)[1]
    if(peakWidth > 1){
      peak[i] <- peak[i] + floor(peakWidth/2) 
    }
  }
  if(sub){
    subPeaks <- mapply(function(s,e, x) {
          if(s == e) p <- s + offset
          else if(isTRUE(all.equal(diff(IRanges::window(x, s, e)), rep(0, e-s)))) 
            p <- s + floor((e - s)/2) + offset
          else if(e-s > 1){
            p <- which(diff(IRanges::window(x, s, e-1)) >= 0 &  
                    diff(IRanges::window(x, s+1, e)) <= 0) + s + offset
            if(x[s] > x[s + 1]) p <- c(s + offset, p)
            if(x[e] > x[e - 1]) p <- c(p, e + offset)
            drop <- which(diff(p) == 1)
            if(length(drop) > 0) p <- p[-(drop + 1)]
          }
          else p <- which.max(x[s:e]) + s -1 + offset
          p
        },idx[idx.start], idx[idx.end], 
        MoreArgs=list(x=score), SIMPLIFY=FALSE)
    ret <- list(peaks=peak+offset, subPeaks=subPeaks)
  }
  else ret <- peak + offset
  
  return(ret)
}
