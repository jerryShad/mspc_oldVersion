## MSPC project
@ author: Julaiti Shayiding

@ Description : Project IDEA    <Revised version>

step 1: read the peak from the source replicate genes; and save it in data.frame tabular data objects.
step 2: convert from score to p-value; so make access to each peak with p-value.
step 3: Initial peak classification; make to access to source replicate genes; and grab the peak with its p-value; 
        then use function call "initial peak classification " to classify the peak element to stringent (or weak ) peak;
        
step 4: Before to find out overlapped peak from the reference genes, make sure all peak regions from the source and reference gene
        are placed on the IntervalTree; To do this, I can import "IRanges" library;
        
        IntervalTree function only accept peak regions, so I need to design algorithm how to get detailed overlapped peak information
        for instance, when function report the overlapped regions, I need to make access that peak element with its p-value;
        
step 5: stick on this: function only grab one peak and pass it to peak classification, overlapped peak detection function; 
        so if peak is classified as stringent peak, then pass it to the next functin to find out overlapped regions; 
        if there is no overlapped peak found in reference gene; 
        then peak is discarded and generate output bed files for stringent or weak peak bed file if the peak were classified as stringent or weak peak;
        additinally generate output bed files for discarded peak (discarded peak could be stringent peak, or weak peak)
        
        
        
step 6: based on the result of step 5, function can also report how many overlapped region are found; then further check Is the number of 
        total overlapped regions are enough to proceed Fisher' exact test; If function report sufficient number of overlapped peak for doing 
        Fisher' test, make sure I can access the information of the overlapped peak with its p-value; (look up fisher' test formula)
        
step 7: based on the result of step 6, If fisher' exact test can pass, make sure the function can generate outuput bed files for this peak;
        so output bed files are generated for stringent_confirmed, weak_discarded peak, and so on;

step 8: when process of step 7 is done; readPeak function could grab next peak from source gene replicate,
        and repeat the process in step 3, step 4, step 5, step 6, step 7; Do the same things
        
step 9: for each individual function, I must give test code for how to run that specific function and make sure that result is correct

step 10: Test all function I have implemented and fix the bug if there is any
step 11: write initial documentation for each function and parameters that make sure user can understand the document and easy to use the MSPC package

step 12: attempt to deliever the MSPC package to the Bioconductor reposotory; keep contact with them and refine the code

Step 13: close up the project of MSPC package development;

#### Expected time to refine all MSPC projects: 1 month
###  Expected delievry date : 2016-1-10
          
        
