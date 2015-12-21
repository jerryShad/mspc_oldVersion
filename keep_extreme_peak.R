 ###################  iterate chromosome and keep only extreme peak      ###################
 
  for(j in 1:length(chr))
  {
    ##get peaks on this chr
    ch = chr[j]
    sel <- x$peaks[,1] == ch
    x.sel <- x$peaks[sel,]

    ##threshold
    if(length(threshold) == 1)
    {
      thr <- threshold
    } else {
      thr <- threshold[x.sel$job]
    }
    x.sel <- x.sel[x.sel$PP > thr, ]

    if(is.unsorted(x.sel$start)) {x.sel <- x.sel[order(x.sel),]}##ensure output is sorted

    ##if there are duplicates, then only keep bin with largest PP
    dup <- which(duplicated(x.sel$start))
    dup.all <- sort(c(dup - 1, dup))
    dup.start <- x.sel$start[dup.all]
    for(i in unique(dup.start))
    {
      sel <- dup.all[dup.start == i]
      ##overwrite all PP values with max, then delete later
      temp <- max(x.sel$PP[sel])
      x.sel$PP[sel] <- temp
    }
    if(length(dup)>0) {x.sel <- x.sel[-dup,]}

    ##exit if chromosome has been reduced to no reads
    if(nrow(x.sel) == 0)
    {
      #warning("No enrichment found on ", ch, " at PP > ", threshold)
      unenriched.chr <- c(unenriched.chr, as.character(ch))
      next
    }

    boundaries <- take.union(x.sel[,2:3])

    start <- which(x.sel$start %in% boundaries$start) ##OUTPUT
    end <- which(x.sel$end %in% boundaries$end) ##OUTPUT

    sel <- apply(rbind(cbind(start, end), 0), 1, function(y){y[1]:y[2]}) ##build a list of regions
    sel <- sel[-length(sel)]

    ##combine probabilities - choice of technique
    if(method == "max")
    {
      p <- unlist(lapply(sel, function(y){max(x.sel$PP[y])})) ##build PP vals

    } else if (method == "lowerbound") {

      ##a non overlapping set has probability 1 - product (1-P). Find the set with the highest lower bound via dynamic programming.
      p <- rep(0, length(sel))
      halfwidth <- (x.sel$end[1] - x.sel$start[1])/2

      for(i in 1:length(sel)) ##in each region...
      {
        temp <- x.sel[sel[[i]],c("start","PP")]
        temp$start = floor((temp$start - temp$start[1])/halfwidth + 1)
        Q = rep(1, max(temp$start)) ##can fill gaps with Q = 1-P = 1
        Q[temp$start] = 1 - temp$PP

        best <- c(1,1)
        for(j in 1:length(Q))
        {
          slot <- (j %% 2) + 1
          best[slot] <- min(best[3-slot], Q[j]*best[slot])
        }
        p[i] = 1 - min(best)
      }

    } else {
      stop(paste("method = '", method, "' not supported", sep = ""))
    }

    temp <- data.frame(chr = ch, start = boundaries$start, end = boundaries$end, PP = p)
    output <- rbind(output, temp) ##append
  }

  if(length(unenriched.chr) > 0) {warning("No enrichment found on at PP > ", threshold , " on the following chromosomes: ", paste(unenriched.chr, collapse = ", "))}

  ##output
  if(is.null(output)) {return(NULL)}
  IRanges::RangedData(IRanges(start = output$start, end = output$end), PP = output$PP, space = output$chr)
}
