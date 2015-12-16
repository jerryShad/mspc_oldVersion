library(rtracklayer)

loadReplicate <- function(bedFile) {
  if (!file.exists(bedFile)) {
    stop("Can't find BED file: ",bedFile);
  }
  Replicate<-import.bed(bedFile)
  Replicate<-as(Replicate, "data.frame")
  Replicate<-Replicate[,1:7]
  names(Replicate) <- c("chrom","chromStart","chromEnd","witdth", "strand", "name", "score");
  # Check columns.
  for (col in c("chrom","chromStart","chromEnd", "witdth", "strand","name", "score")) {
    if (!col %in% names(Replicate)) {
      stop(sprintf("error in peaks data frame: no column named '%s'",col));
    }
  }
  return(Replicate)
}

################# read parent folder of genome #################

GeneSet<-list.files(path, pattern, full.names = T)
Replicate_List<-lapply(GeneSet, loadReplicate)

#' Need more compactable function for geneset parent folder
#' @example
#' parentGenomeFolder <- list.files("C:/Users/Jvret/Documents/ENCODE_Samples", pattern="bed", full.names=TRUE)
#' temp<-lapply(parentGenomeFolder, loadReplicate)
#' temp
