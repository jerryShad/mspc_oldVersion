##  read peak from gene replicate

readPeak<-function(Replicate, ...){
  if(missing(Replicate)){
    stop("there is no valid gene replicate found", Replicate)
  }
  if(!missing(Replicate)){
    source_peak<-Replicate[1:nrow(Replicate),]
  }
  return(source_peak)
}

## example code
Replicate<-loadReplicate(parentGenomeFolder[2])
peak<-readPeak(Replicate)
print(peak)
