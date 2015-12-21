##############  Fisher' Exact Test    ##############

test_fisher_exact = function(geneset,gpw,alternative="two.sided") {
  # Restrict to only those genes in the genesets.
  gpw = subset(gpw,geneid %in% geneset@all.genes);

  genes = gpw$geneid;
  peaks = gpw$peak;

  results = lapply(ls(geneset@set.gene), function(go_term) {
    go_genes = geneset@set.gene[[go_term]];
    n_go_genes = length(go_genes);

    in_cat = as.numeric(genes %in% go_genes);
    n_go_peak_genes = sum((genes %in% go_genes) & (peaks == 1));
    xt = table(in_cat,peaks);

    pval = 1;
    odds_ratio = 0;
    try({
      fet_result = fisher.test(xt,conf.int=F,alternative=alternative);
      pval = fet_result$p.value;
      odds_ratio = fet_result$estimate;
    },silent=T);

    go_peak_genes = gpw[(gpw$geneid %in% go_genes) & (gpw$peak == 1),]$geneid;
    go_peak_genes_str = paste(go_peak_genes,collapse=", ");

    enr = NA;

    if (odds_ratio > 1) {
      enr = "enriched";
    }

    if (odds_ratio < 1) {
      enr = "depleted";
    }

    if (odds_ratio == 0) {
      enr = "depleted";
    }

    if (is.infinite(odds_ratio)) {
      enr = "enriched";
    }

    data.frame(
      "Geneset ID" = go_term,
      "N Geneset Genes" = n_go_genes,
      "N Geneset Peak Genes" = n_go_peak_genes,
      "Geneset Peak Genes" = go_peak_genes_str,
      "Odds.Ratio" = odds_ratio,
      "Status" = enr,
      "P-value" = pval,
      stringsAsFactors=F
    );
  });

  results = rbind.fill(results);

  results$FDR = p.adjust(results$P.value,method="BH");

  results = results[order(results$P.value),];
  return(results);
}


##############  collect peak analysis result; construct new data.frame for it    ##############

# Collapse results into one table
  results = Reduce(rbind,results_list)

  # Correct for multiple testing
  results$FDR = p.adjust(results$P.value, method="BH");

  # Create enriched/depleted status column
  results$Status = ifelse(results$Effect > 0, 'enriched', 'depleted')

  results = results[order(results$P.value),];

  return(results);
  

##############  restrict the source peak to target genes in genesets    ##############
   # Restrict our genes/weights/peaks to only those genes in the genesets.
  gpw = subset(gpw,geneid %in% geneset@all.genes);

  if (sum(gpw$peak) == 0) {
    stop("Error: no peaks in your data!");
  }
  
  
  
##############  get peak from user' file    ##############


# Get peaks from user's file.
  if (class(peaks) == "data.frame") {
    peakobj = load_peaks(peaks);
  } else if (class(peaks) == "character") {
    if (get_ext(peaks) == "bed") {
      message("Reading BED file: ",peaks);
      peakobj = read_bed(peaks);
    } else {
      message("Reading peaks file: ",peaks);
      peakobj = read_peaks(peaks);
    }
  }

  # Number of peaks in data.
  num_peaks = sum(sapply(peakobj,function(x) length(x)))
  
  
##############  assign the peak to target gene    ############## 
 # purpose: to find out the mathed peak that overlapped with source peak
  # Assign peaks to genes.
  assigned_peaks = assign_peaks(peakobj,ldef,tss);
  peak_genes = unique(assigned_peaks$geneid);
  
  
##############  2nd approach: assign the peak to target gene    ############## 
    # Assign peaks to genes.
  assigned_peaks = assign_peak_segments(peakobj,ldef);
  peak_genes = unique(assigned_peaks$geneid);

  ppg = num_peaks_per_gene(assigned_peaks,ldef,mappa=NULL)

  ppg = calc_peak_gene_overlap(assigned_peaks,ppg)
  
##############  order the peak & calculate FDR    ############## 
    # Order by length.
  gpw = gpw[order(gpw$log10_length),];

  # Calculate prob(all false positives).
  gpw$false_prob = 1 - (1 - (gpw$length/genome_length))^(num_peaks);
  
  
  ##############  create IRanges objects    ############## 
  # Create an IRanges object representing the loci for each gene on that chromosome.
  chroms = list();
  for (chr in unique(d$chrom)) {
    genes_chrom = subset(d,chrom == chr);
    chroms[[chr]] = IRanges(start=genes_chrom$start,end=genes_chrom$end,names=genes_chrom$geneid);
  }

  object@dframe = d;
  object@chrom2iranges = chroms;

  # Store as GRanges object as well, for convenience.
  object@granges = GRanges(
    seqnames=d$chrom,
    ranges=IRanges(d$start,d$end),
    names=d$geneid
  );
  
  ##############  IRanges objects manipulation    ############## 
  gc()
  IR <- IRanges::IRanges(chropmStart = Replicate$chromStart, chromEnd = Replicate$chromEnd)
  gc()
  IRanges::RangedData(IR, chrom = Replicate$chrom, name = Replicate$name, pvalue=Replicate$pvalue)


  ##############  pull out the peak that didn't pass the Fisher' test    ############## 
 # Pull out tests that failed.
  bad_enrich = subset(enrich,is.na(P.value));
  enrich = subset(enrich,!is.na(P.value));


 ##############  write peak analysis result that done by MSPC package and export to the output bed files    ############## 
  # Write results to file.
  if (!is.null(out_name)) {
    filename_analysis = file.path(out_path,sprintf("%s_results.tab",out_name));
    write.table(enrich,file=filename_analysis,row.names=F,quote=F,sep="\t");
    message("Wrote results to: ",filename_analysis);

    filename_peaks = file.path(out_path,sprintf("%s_peaks.tab",out_name));
    write.table(assigned_peaks,file=filename_peaks,row.names=F,quote=F,sep="\t");
    message("Wrote peak-to-gene assignments to: ",filename_peaks);

    filename_opts = file.path(out_path,sprintf("%s_opts.tab",out_name));
    write.table(opts,file=filename_opts,row.names=F,quote=F,sep="\t");
    message("Wrote run options/arguments to: ",filename_opts);

    filename_ppg = file.path(out_path,sprintf("%s_peaks-per-gene.tab",out_name));
    write.table(ppg,file=filename_ppg,row.names=F,quote=F,sep="\t");
    message("Wrote count of peaks per gene to: ",filename_ppg);
