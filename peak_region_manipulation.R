 ### Intialize
 QC <- data.frame()
  peaks <- list()
  p.samples <- list()

  ##collect start and end co-ordinates (TODO Take centromere into account)
  
  jobs <- rep(0L, length(chr))      
  
  ##################TODO BEGIN:  collect coordinate for each peaks    ##################
  pks<-rep(0l, length(peak))
  pk.region<- cbind(chromStart, chromEnd, pks)    # information about each peak on current chromosome
  
  
  ##################TODO END; ##################
  
  region.full <- cbind(start, end, jobs)
  colnames(region.full) = c("start", "end", "jobs") ##information about each chromosome to be used later

  if(any(is.na(start) | is.na(end)))
  {
    for(i in 1:length(chr))
    {
      ch = chr[i]
      if(any(is.na(start))) ##autofit start location
      {
#        temp <- c(sapply(treatment, function(y){min(y$x[y$chr == ch])}), sapply(control, function(y){min(y$x[y$chr == ch])}))
#        region.full[i, 1] <- min(unlist(temp))
        temp <- sapply(treatment, function(y){min(y$x[y$chr == ch])})
        region.full[i, 1] <- min(unlist(temp))
      }
      if(any(is.na(end))) ##autofit end location
      {
#        temp <- c(sapply(treatment, function(y){max(y$x[y$chr == ch])}), sapply(control, function(y){max(y$x[y$chr == ch])}))
#        region.full[i, 2] <- max(unlist(temp))
        temp <- sapply(treatment, function(y){max(y$x[y$chr == ch])})
        region.full[i, 2] <- max(unlist(temp))
      }
    }
  }
