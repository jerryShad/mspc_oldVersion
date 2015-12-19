## data.frame manipulation
## pre-done for p-value conversion

manipulate_df<-function(Replicate, ...){
  if(missing(Replicate$witdth && Replicate$strand)){
    message("gene replicate format is desirable", Replicate)
  }
  if(!missing(Replicate$witdth && Replicate$strand)){
    Replicate$witdth <- NULL
    Replicate$strand <- NULL
  }
  name(Replicate)<-c("chrom","chromStart","chromEnd","name", "score");
  Replicate<-Replicate[,1:5]

  if(missing(Replicate$pvalue)){
    message("replicate gene missing pvalue conversion", Replicate)
    .Replicate <- cbind(Replicate, pvalue = 10^-10/Replicate$score)
  }
  return(.Replicate)
}


## manipulate original data.frame and add p-value column on it
