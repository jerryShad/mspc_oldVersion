## manipulate gene and split them by chromosome

get_peak_ranges <- function(source_peak, ...){

  chroms = list();
  for (chr in unique(source_peak$chrom)) {
    peaks_chrom = subset(source_peak,chrom == chr);
    chroms[[chr]] = IRanges(start=peaks_chrom$chromStart,end=peaks_chrom$chromEnd);
  }
  return(chroms)
}
