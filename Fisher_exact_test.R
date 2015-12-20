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
  
  
  
  
  
