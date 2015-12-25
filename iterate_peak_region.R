#################   Iterate peak regions    #################
## below code may give me some idea how to iterate peak regions and compare with the region in Intervaltree

 if(length(s$feature)>0){
            step <- step[rownames(s), ]
            step <- as.data.frame(step)
            s <- as.data.frame(s)
            for(i in c("peak", "feature", "start_position", 
                       "end_position", "distancetoFeature", 
                       "shortestDistance")){
                if(any(step[,i]!=s[,i])) stop(paste(i, "is not identical!"))
            }
        }
        
        
  #################   R code block for nested iteration    #################      
  ## learning purpose: apply this iterating data structure to MSPC package
  
  
   if (!is.null(gene.filters))
    {
        disp("Gene filters: ",paste(names(gene.filters),collapse=", "))
        for (gf in names(gene.filters))
        {
            disp("  ",gf,": ")
            for (gfp in names(gene.filters[[gf]]))
            {
                if (length(gene.filters[[gf]][[gfp]])==1 && 
                    is.function(gene.filters[[gf]][[gfp]]))
                    print(gene.filters[[gf]][[gfp]])
                else if (length(gene.filters[[gf]][[gfp]])==1)
                    disp("    ",paste(gfp,gene.filters[[gf]][[gfp]],sep=": "))
                else if (length(gene.filters[[gf]][[gfp]])>1)
                    disp("    ",paste(gfp,paste(gene.filters[[gf]][[gfp]],
                        collapse=", "),sep=": "))
            }
        }
    }
