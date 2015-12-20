
########################## Code block is used for reference; How to iterate over gene from gene sets   #######################
for (i in seq_along(analysisResults)) {
                    val <- analysisResults[[i]]
                    if (comparison(val)=="High vs Normal/Low") {
                        if (analysedSubset(val)=="all"){
                            high_all_p.value<-as.numeric(as.vector(getColumnView(val))[2])      
                        }
                        if (analysedSubset(val)=="males"){
                            high_male_p.value<-as.numeric(as.vector(getColumnView(val))[2])                  
                        }
                        if (analysedSubset(val)=="females"){
                            high_female_p.value<-as.numeric(as.vector(getColumnView(val))[2])                  
                        }
                    }
                    else {
                        if (analysedSubset(val)=="all"){
                            low_all_p.value<-as.numeric(as.vector(getColumnView(val))[2])                  
                        }
                        if (analysedSubset(val)=="males"){
                            low_male_p.value<-as.numeric(as.vector(getColumnView(val))[2])                  
                        }
                        if (analysedSubset(val)=="females"){
                            low_female_p.value<-as.numeric(as.vector(getColumnView(val))[2])                  
                        }
                    }
                }
                
                # High classification p-val is less than threshold and low classification p-val is more that threshold
                if (high_all_p.value < phenotypeThreshold && low_all_p.value >= phenotypeThreshold){
                    direction_all <- "High"
                    all_p.value <- high_all_p.value
                }
                # Low classification p-val is less than threshold and high classification p-val is more that threshold
                else if (high_all_p.value >= phenotypeThreshold && low_all_p.value < phenotypeThreshold){
                    direction_all <- "Low"
                    all_p.value <- low_all_p.value
                }
                
  
  ####################### Code block is used for reference; How to iterate over the peak from source peaks    #######################
  @ iterate peak over source peak data.frame object
                
                 for (i in seq_along(analysisResults)) {
                    val <- analysisResults[[i]]
                    if (analysedSubset(val)=="all"){
                        all_p.value<-as.numeric(as.vector(getColumnView(val))[2]) 
                               
                    }
                    if (analysedSubset(val)=="males"){
                        male_p.value<-as.numeric(as.vector(getColumnView(val))[2])             
                    }
                    if (analysedSubset(val)=="females"){
                        female_p.value<-as.numeric(as.vector(getColumnView(val))[2])                
                    }
                }
