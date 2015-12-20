## Idea for overlapped peak detection

# key question I am going to ask from R-Studio

1. I gave the peak regions information to let know R that this my input; I also give target reference gene' regions to IntervalTree;

Question 1: ask R to bring me back overlapped peak regions based on my input peak regions from target gene' regions
Question 2: once R report overlapped regions of the peak; I further ask to R which peak overlapped, which chromosome they are come from,
            what's these overlapped peak' pvalue; 
            
Question 3: I am going to ask R-studio to collect all peak regions and place them in IntervalTree; How can I tell to R that iteratively             select the ranges of peak one by one in target_gene, compare the region with source peak resions, and further evaluate where             is the overlapped regions occured, which peak are they ? what's their p-value significant; Do we have sufficient overlapped              regions that asked by user.

Question 4: how R dynamically look up the set of regions and , to find matched overlapped region with source peak region.
