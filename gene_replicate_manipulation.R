## manipulate gene replicate

# example code: how to ignore the unnecessary column 

##: function for load replicate

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

step 1: load gene replicate fronm parent genome folder

read_gene_replicate <- function (parentGenomeFolder, ...){
  GeneSet<-list.files(path, pattern, full.names = T)
  Replicate_List<-lapply(GeneSet, loadReplicate)
}

step 3: remove unnecessary column in bed files and add p-value column

## example code :
parentGenomeFolder <- list.files("C:/Users/Jvret/Documents/ENCODE_Samples", pattern="bed", full.names=TRUE)
temp<-lapply(parentGenomeFolder, loadReplicate)
temp<- loadReplicate(parentGenomeFolder[1])

# remove witdth, strand column from data.frame tabule objects

temp$witdth <- NULL       # set witdth column as NULL
temp$strand <- NULL       # set witdth column as NULL


# step 4: convert p-value and add p-value column to data.frame tabular object
# step 5: start this gene replicate where p-value column also added; start to do peak classification

