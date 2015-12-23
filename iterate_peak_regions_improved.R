## iterate peak peak regions

chrom= IRanges(peak$chromStart, peak$chormEnd)
results<-IRanges::as.data.frame(overlap)

 for (i in seq(len=length(chrom))) {
      peak_start <- start(peak$chromStart)[i]
      peak_end <- end(peak$chromEnd)[i]
      peak_IR <- IRanges(peak_start, peak_End)
      
      # how can I scan through the peak region in target gene; how to iterate them
      
      overlap<-findOverlap(peak_IR, target_gene, algorithm=" intervaltree")
    }
    
    
 ## this is how the parameter to be verified:
 ## for instance, we are checking required parameter in findOverlap() function
   type <- match.arg(type)
   select <- match.arg(select)
   algorithm <- match.arg(algorithm)
