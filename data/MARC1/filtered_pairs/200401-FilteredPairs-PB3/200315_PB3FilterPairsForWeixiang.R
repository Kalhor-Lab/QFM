################################################################################################
################################ LOADING PACKAGES ##############################################
################################################################################################

require(entropy);
require(ape);
require(pvclust);
library(stringdist)
require(RColorBrewer); divcol <- c("RdBu", "PuOr", "PRGn", "PiYG", "BrBG", "RdYIBu"); qualcol <- c("Paired", "Dark2", "Set2"); seqcol <- c("YlOrBr","YlOrRd",  "YlGnBu", "RdPu", "PuRd", "PuBu", "OrRd", "YlGn", "GnBu", "BuPu", "BuGn")
require(beanplot)

################################################################################################
############################# LOADING DATA AND FILTERING #######################################
################################################################################################

setwd('/Users/Reza/OneDrive - Johns Hopkins/Work Documents/Johns Hopkins/Notebooks/SimulationFramework-withWeixiang/200401-FilteredPairs-PB3')

MinRd <- 5000                                   # Minimum total number of reads for a sample to be considered.
PFCOFF_bc <- NULL                                # Read Percentage Cutoff is the the percentage of total reads in a sample that a barcode needs to have to be considered at all.
                                                # This value will be assigned later based on the total number of observed barcodes in each sample.
PFCOFF_pair <- 0.01                              # Read fraction Cutoff is the the fraction of total reads in a barcode that a pair needs to have to be considered at all. 
                # I have found that in the c84m3 OrganLib1 results, 0.003 is the fraction below which one-base shifted variants of the spacer appear.
                # 0.005 was the level at which in the parent most barcodes appeared 100% non-shifted (31 out of 34).
                # I tried defining diversity-based cutoffs, but realized that for the results to be unbiased, the same cutoff should be used for the same barcode in all samples.
PRCOFF_pair <- 2                                # Pairs with equal to or fewer than PRCOFF_pair reads in total will be removed.
timepoints <- c("j1763m3", "E03p5", "E06p5", "E07p5", "E08p5", "E10p5", "E12p5", "E14p5", "E15p5", "E16p5", "P21")
timenames <- c("E0", "E3.5", "E6.5", "E7.5", "E8.5", "E10.5", "E12.5", "E14.5", "E15.5", "E16.5", "Adult")
names(timepoints) <- timenames

spurious_barcodes <- c("ACTCAACGGC", "GAGTCTGGCA", "GCCGAAATCA", "GACCCGACCA", "TAAACTCAAG", "TAATACGACT", "CATACTTAAC", "TCAATACACT", "CCACACGACA", "TATGTCATGC", "CAATCGCACT", "ACACCGGTGT")    # These barcodes have been determined to be spurious based on following reasons and will be removed from each sample:
                                        # ACTCAACGGC belongs to 1763m4 and its presence can only be cross-contamination.
                                        # GAGTCTGGCA belongs to 1763m4 and its presence can only be cross-contamination.
                                        # GCCGAAATCA belongs to 1763m4 and its presence can only be cross-contamination.
                                        # GACCCGACCA belongs to 1763m4 and its presence can only be cross-contamination.
                                        # TAAACTCAAG belongs to 1763m4 and its presence can only be cross-contamination.
                                        # TAATACGACT is observed in very high abundance in empty samples (no-template negative controls).
                                        # CATACTTAAC belongs to 1763m4 and its presence can only be cross-contamination.
                                        # TCAATACACT belongs to 1763m4 and its presence can only be cross-contamination.
                                        # CCACACGACA belongs to 1763m4 and its presence can only be cross-contamination.
                                        # TATGTCATGC belongs to 1763m4 and its presence can only be cross-contamination.
                                        # CAATCGCACT belongs to 1763m4 and its presence can only be cross-contamination.
                                        # ACACCGGTGT is not spurious, but I cannot determine which barcode it is derived from. However, since it is only seen in one sample, it shouldn't alter the results. 
spurious_barcodes <- c(spurious_barcodes, "AAACCCCGGG","AACGCCCTAC","AACTATCGGC","AACTCACCTA","AAGACTTCAT","AAGCCGCGCG","ACATTCGGTT","ACCACTGCTG","ACCCTGGGAC","ACTCCATGTT","ACTCGGTTAC","AGACCCTCGC","AGCACTGTAC","AGCCCAAATC","AGTCTGCCTC","AGTTTCCGAA","ATGCTTAGCT","ATGGCGCCTA","CACAACGCCC","CACTCTCAAG","CATCGTCGTC","CATTGGAGGT","CCCAAAACAC","CCCATCACCC","CCCCTCACTT","CCTCACCCCA","CCTTTACCGC","CGAATCCTTT","CGACAGTTAT","CGCATGATGC","CGCCGTAGTA","CGTAGGGCCC","CGTGTTGTCT","CTACTCGGCC","CTGAGTTTTA","CTGCTATCGA","CTTTTGTCGG","GAAGACCCGC","GACACAGACA","GACCTCCAAT","GATACCCCCA","GCCAAGATGG","GCCAGCCGCT","GCCCCATTCC","GCGAAGTCCC","GCTCTACGCC","GGCACCCTCC","GGCCCCTACA","GGGTGACACG","GTACACAATT","GTCAAATACC","GTGGAGCCTC","TAACTGCTCT","TAACTTATAC","TAGCCATGCA","TCTATCGAGG","TCTCTAGATC","TTAGCTATGT","TTGAGCATAA","TTTGGCACAC")
                                        # These barcodes blown to 1763m7, their presence here can only be cross-contamination.

TEMPLATE_SWITCHING_CORRECTION <- TRUE
max_dist_spacer <- 1;                           # Since spacers for each sample are converted into a consensus separately, it is possible that a specific sample has the reference spacer in the parent only with a PCR or sequencing-based single point mutation which causes it to appear 100% mutated. Such a problem was observed with samples with smaller reads. As a result, the string.dist(method = "hamming") strategy will be re-applied later in this code to count slight variants of the parental sequence as non-mutants.
parent <- "parent-j1763m3-G9b,c"
files1 <- system('ls /Users/Reza/OneDrive\\ -\\ Johns\\ Hopkins/Work\\ Documents/Harvard/GMC_Lab/Notebook_TransgenicMouse/180502_PB3EmbryonicTimeCourse/MiSeq/180502_PB3EmbryonicTimeCourse/3-*adjustment2/*_truepairs.txt', intern = TRUE)                  # These are the processed barcode-spacer pair counts for all samples

PAM_pattern1 <- "GGGT" ; PAM_pattern1_offset <-  0                             # This pattern will be used to extract the length of hgRNA. Offset determines the distance between the pattern and the actual point of interest, PAM here.
PAM_pattern2 <- "GGTT" ; PAM_pattern2_offset <-  -1                            # This pattern will be used to extract the length of hgRNA. Offset determines the distance between the pattern and the actual point of interest, PAM here.Two patterns are being used in case there is a mutation in one patterns region. The two patterns were chosen based on minimal (ideally 0) overlap and non-occurance in the pre-barcode region of the gRNA.
TSS_pattern1 <- "CCGG" ; TSS_pattern1_offset <-  +3                            # This pattern will be used to extract the length of hgRNA, it marks the sequence before TSS. Offset determines the distance between the pattern and the actual point of interest, TSS here.
#TSS_pattern2 <- "AACA" ; TSS_pattern2_offset <-  +7                            # This pattern will be used to extract the length of hgRNA, it marks the sequence before TSS. Offset determines the distance between the pattern and the actual point of interest, TSS here. Two patterns are being used in case there is a mutation in one patterns region.
Post_bc_pattern <- "GAATTC"
barcode_length <- 10
divergence_point <- 12                                                         # Position in the spacer at which the sequence is expected to start to diverge



alldata_bysample <- list()                                                     #All of the data separated by sample in each object within the list
alldata_bybarcode_bysample <- list()                                           #All of the data separated by barcode (i.e., gene) in each object within the list
all_pairs <- vector()                                                          # This vector stores a list of all pairs seen in all samples
allspacers_bybarcode <- list()                                                 # This list stores a list of all pairs seen for each barcode samples
all_barcodes <- vector()

##  Reading the truepairs file for all samples and storing them as a list of dataframes in alldata_bysample
name_base <- paste(Sys.Date(), '_PB3SpacerConsolidationReport', sep = "")
filename_text <- paste(name_base, '.txt', sep = '')
sink(file = filename_text); sink();
for (file1 in files1) {
    sample_name <- tail(unlist(strsplit(file1, split = "/")), n=1); sample_name <- unlist(strsplit(sample_name, split = "_"))[1]
    print(paste("processing", sample_name, Sys.time(), sep = " "))
    sink(file = filename_text, append = TRUE); cat(paste("processing", sample_name, Sys.time(), "\n", sep = " ")); sink()
    
    pairs <- NULL;
    pairs <- read.table(file1, colClasses=c("character", "character", "numeric"))                          # Reading all pairs of the sample
    if (nrow(pairs) < 1 ) { next; }
    pairs <- subset(pairs, !V1 %in% spurious_barcodes)
    if (sum(pairs[,3]) < MinRd) { next; }

    ## Applying length-based filters    (in some rare samples short spacers have been observed which probably means a large deletion extending to the primer)
    if( max(nchar(pairs[,2])) != min(nchar(pairs[,2]))) {
        print(unname(cbind("Removed", pairs[nchar(pairs[,2]) < max(nchar(pairs[,2])),], "because of short spacer")))
        pairs <- pairs[nchar(pairs[,2]) == max(nchar(pairs[,2])),]
    }

    ## Correcting template switching for all samples (the rational is that if the exact same spacer is seen in two barcodes within the same sample, it belongs to only one of those barcodes and that is the barcode where that spacer has a higher count)
    if(TEMPLATE_SWITCHING_CORRECTION) {
        byspacers <- NULL; byspacers <- aggregate(V3~V2, pairs, max)
        for (sp in 1:nrow(byspacers)) {
            pairs <- subset(pairs, V2 != byspacers[sp,1] | (V2 == byspacers[sp,1] & V3 == byspacers[sp,2]) )        # For each spacer sequence, removing all instances, regardless of barcode, that have a submaximum frequency in the sample (since template switching can only shuffle other spacers in the same sample)
        }
    }

    ## [200305] : The "pairs" table has to be sorted properly for below steps to work appropriately (especially the one-base displacement adjustment). Therefore, the following command is added to make sure
    pairs <- pairs[order(pairs$V1, -pairs$V3),]   #sorting pairs on barcode followed by descending frequency because some filtering operations above can disrupt order
    
    pairs$V4 <- NA;                                     # For each entry, this column indicates its percentage among pairs with the same barcode in that sample
    PFCOFF_bc <- NULL; PFCOFF_bc <- 1/length(levels(as.factor(as.vector(pairs[,1]))))^1.75                       # This is reminiscent of how I narrowed the number of barcodes in each sample. But for 2-sequencing_error_adjustment1 the maximum expected number for all samples combined was used. Here it will be tailored to each sample. 
    for(bc in levels(as.factor(as.vector(pairs[,1])))) {
        if ( sum(subset(pairs, V1 == bc)[,3]) / sum(pairs[,3]) > PFCOFF_bc) {
            if (length(which(pairs[,3] <= PRCOFF_pair & pairs[,1] == bc, arr.ind = TRUE)) > 0 ) {pairs <- pairs[-which(pairs[,3] <= PRCOFF_pair & pairs[,1] == bc, arr.ind = TRUE),]}           # Eliminating pairs that don't meat the PRCOFF_pair criteria
            
            # [180515] module for replacing a specific 'N' that appears in some spacers at a specific position
            pairs$V2 <- gsub("ACACCGN", "ACACCGG", pairs$V2)
            pairs$V2 <- gsub("ACACCNG", "ACACCGG", pairs$V2)
            pairs$V2 <- gsub("ACACNGG", "ACACCGG", pairs$V2)
            pairs$V2 <- gsub("ACANCGG", "ACACCGG", pairs$V2)
            pairs$V2 <- gsub("CACCGGN", "CACCGGT", pairs$V2)
            
            # [180315] module for consolidating one-base displacements due to sequencing error, such as "GAAACACCGGTAGCAAACGTTTGGACGTGGGGTTAGAGCTAGAAATAGCA" and "CAAAACACCGGTAGCAAACGTTTGGACGTGGGGTTAGAGCTAGAAATAGC" that create an erroneous mutant spacer
            if (nrow(subset(pairs, V1 == bc)) > 1) {
                for (spacer_obs in subset(pairs, V1 == bc)$V2[length(subset(pairs, V1 == bc)$V2):2]) {
                    spacer_obs_trunc <- substr(spacer_obs, 5, nchar(spacer_obs))
                    matching_spacers <- grep(spacer_obs_trunc, pairs$V2); matching_spacers <- matching_spacers[matching_spacers %in% which(pairs$V1 == bc)]
                    if(length(matching_spacers) > 1 ) {                     # Then the first one would be the parent and last one would be spacer_obs
                        sink(file = filename_text, append = TRUE); cat(paste("Combined ", paste(pairs[matching_spacers[length(matching_spacers)],1:3], collapse = " "), " into ", paste(pairs[matching_spacers[1],1:3], collapse = " "), "\n", sep = "")); sink();
                        pairs[matching_spacers[1],3] <- pairs[matching_spacers[1],3] + pairs[matching_spacers[length(matching_spacers)],3]
                        pairs <- pairs[-matching_spacers[length(matching_spacers)],]
                    }
                }
                
            }
            
            cutoff <- sum(pairs[,3][pairs[,1] == bc]) * PFCOFF_pair
            if (length(which(pairs[,3] < cutoff & pairs[,1] == bc, arr.ind = TRUE)) > 0 ) {pairs <- pairs[-which(pairs[,3] < cutoff & pairs[,1] == bc, arr.ind = TRUE),]}           # Eliminating pairs that don't meat the PFCOFF_pair criteria
            pairs[,4][pairs[,1] == bc] <- round(pairs[,3][pairs[,1] == bc]/(sum(pairs[,3][pairs[,1] == bc]))*100, digits = 2)
            alldata_bybarcode_bysample[[bc]][[sample_name]] <- subset(pairs, V1 == bc)
            allspacers_bybarcode[[bc]] <- append(allspacers_bybarcode[[bc]], as.character(pairs[,2][pairs[,1] == bc])); allspacers_bybarcode[[bc]] <- unique(allspacers_bybarcode[[bc]])
        } else {pairs <- subset(pairs, V1 != bc)}
    }
    
    pairs$V5 <- NA;                                     # For each entry, this column indicates its percentage among all pairs in that sample
    pairs$V5 <- pairs[,3]/sum(pairs[,3])*100
    
    all_pairs <- append(all_pairs, paste(pairs[,1], pairs[,2], sep = "-"))
    all_barcodes <- append(all_barcodes, as.character(pairs[,1]))
    alldata_bysample[[sample_name]] <- pairs;
}
all_pairs <- unique(all_pairs)
all_pairs <- all_pairs[order(all_pairs)]
all_barcodes <- unique(all_barcodes)
all_barcodes <- all_barcodes[order(all_barcodes)]

ext <- NULL;
for(sampl in names(alldata_bysample))  {
    temp <- length(which(!unique(alldata_bysample[[sampl]][,1]) %in% alldata_bysample[[parent]][,1]));
    #print(paste(sampl, temp, sep = ":   "))
    ext <- append(ext, temp)
}
names(ext) <- names(alldata_bysample)


## Determining hgRNA length for each barcode
hgRNA_length <- vector(length= length(all_barcodes), mode = "numeric");                 # This vector will stored the derived length of each hgRNA
names(hgRNA_length) <- all_barcodes;
for (bc in all_barcodes) {
    wt_spacer <- NULL;
    TSS_start <- NULL;
    PAM_start <- NULL;
    if (!is.null(alldata_bybarcode_bysample[[bc]][[parent]])) {                         # If the barcode is observed in the parent, the most frequent spacer of the parent will be used as reference for hgRNA length extraction
        wt_spacer <- alldata_bybarcode_bysample[[bc]][[parent]][1,2]
        TSS_start[1] <- min(gregexpr(TSS_pattern1, wt_spacer)[[1]]) + TSS_pattern1_offset
        PAM_start[1] <- max(gregexpr(PAM_pattern1, wt_spacer)[[1]]) + PAM_pattern1_offset
        #TSS_start[2] <- min(gregexpr(TSS_pattern2, wt_spacer)[[1]]) + TSS_pattern2_offset
        #PAM_start[2] <- max(gregexpr(PAM_pattern2, wt_spacer)[[1]]) + PAM_pattern2_offset
    }
    else {                                                                              # If the barcode is not observed in the parent, the median of length in all observed spacers will be used as hgRNA length reference, without regard to abundance of an individual spacer.
        wt_spacer <- allspacers_bybarcode[[bc]]
        TSS_start[1] <- median(unlist(lapply(gregexpr(TSS_pattern1, wt_spacer), min))) + TSS_pattern1_offset
        PAM_start[1] <- median(unlist(lapply(gregexpr(PAM_pattern1, wt_spacer), max))) + PAM_pattern1_offset
        #TSS_start[2] <- median(unlist(lapply(gregexpr(TSS_pattern2, wt_spacer), min))) + TSS_pattern2_offset
        #PAM_start[2] <- median(unlist(lapply(gregexpr(PAM_pattern2, wt_spacer), max))) + PAM_pattern2_offset
    }
    #if (TSS_start[1] != TSS_start[2] | PAM_start[1] != PAM_start[2]) {print(paste(bc, TSS_start[1], TSS_start[2], PAM_start[1], PAM_start[2], sep = " "))}
    #TSS_start <- min(TSS_start[TSS_start > 0])
    #PAM_start <- max(PAM_start[PAM_start < 38 + TSS_start])                                     # 38 was chosen here because the longer hgRNA is expected to be 35. Should this change, this value should change accordingly.
    hgRNA_length[names(hgRNA_length) == bc] <- PAM_start[1] - TSS_start[1] + 1
    #hgRNA_length[names(hgRNA_length) == bc] <- PAM_start - TSS_start + 1
}

## Fixing the barcodes that have suffered a large deletion which has taken out most of the scaffold
    
    # Custom string.dist function for 5' aligned measurement of distance based on hamming. First mis-match, will be the end of agreement between the strings.
    stringdist.c <- function (a, b) {
                    a <- as.character(a)
                    b <- as.character(b)
                    if (length(a) == 0 || length(b) == 0) {
                        return(numeric(0))
                    }
                    if (max(length(a), length(b))%%min(length(a), length(b)) != 0) {
                        warning(RECYCLEWARNING)
                    }
                    a_split <- unlist(strsplit(a, ''));
                    b_split <- strsplit(b,'');
                    all_distance <- NULL
                    for (b1 in 1:length(b_split)) {
                        if (length(a_split) != length(b_split[[b1]])) {warning("a and b have to be the same length") }
                        d1 <- a_split == b_split[[b1]]
                        d2 <- length(a_split);
                        for (identity in d1) {
                            if(identity) {d2 <- d2 - 1 }
                            else {break;}
                        }
                        all_distance <- append(all_distance, d2);
                    }
                    all_distance
                }

# Filter 1: Fixing the barcodes that have suffered deletion and have a close counterpart in the parental barcodes
all_barcodes -> all_barcodes_prefilter1
all_pairs -> all_pairs_prefilter1
allspacers_bybarcode -> allspacers_bybarcode_prefilter1
alldata_bysample  -> alldata_bysample_prefilter1
alldata_bybarcode_bysample -> alldata_bybarcode_bysample_prefilter1
hgRNA_length -> hgRNA_length_prefilter1

MaxLargeDelLength <- 10
for (bc in all_barcodes[!all_barcodes %in% alldata_bysample[[parent]][,1]]) {
    if (hgRNA_length[bc] < MaxLargeDelLength) {                # If the calculated length suggests that the spacer could be a result of a large deletion
        EcoRI_loc <- NULL; EcoRI_loc <- median(unlist(lapply(gregexpr(Post_bc_pattern, lapply(alldata_bybarcode_bysample[[bc]], "[", 1,2 )), min)))            # Calculates the position of the EcoRI site that comes right after the barcode. If this spacer has a large deletion, the EcoRI site should be seen.
        if (EcoRI_loc > 5) {
            bc_spacer_count <- NULL; bc_spacer_count <- max(unlist(lapply(lapply(alldata_bybarcode_bysample[[bc]], "[", 1,c(2,3) ), "[", 2)))                     # Want to choose the spacer that was observed the most for this barcode as reference, and this command extract its observed count
            bc_spacer <- NULL; bc_spacer <- unlist(lapply(alldata_bybarcode_bysample[[bc]], "[", 1,c(2,3) ))[which(unlist(lapply(alldata_bybarcode_bysample[[bc]], "[", 1,c(2,3) )) == as.character(bc_spacer_count))-1]             # Extracting the reference spacer sequence based on the maximum count measured in the previous step. This is a really dirty way to do it.
            
            barcode_distances <- NULL; barcode_distances <- stringdist.c(bc, unique(alldata_bysample[[parent]][,1]))        # Stores the distance of each parent barcode from the barcode at hand
            spacer_distances <- NULL
            for (bc2 in unique(alldata_bysample[[parent]][,1])) {
                spacer_distances <- append(spacer_distances, stringdist.c(substr(bc_spacer, divergence_point, EcoRI_loc), substr(alldata_bybarcode_bysample[[bc2]][[parent]][1,2], divergence_point, EcoRI_loc) ))
            }
            
            overlap <- barcode_length + (EcoRI_loc - divergence_point + 1) - barcode_distances - spacer_distances               # It stores the total overlap between the current barcode with a parental barcode, including their spacers' divergent sequence regions.  It should be noted, that a base at the junction of truncated spacer barcode may be counted twice as overlap if it matches that particular point of both, but this is ok - perhaps even desirable - given the observatin that deletions often times happen between two bases of the same identity.
            if (length(overlap[overlap == max(overlap)]) == 1) {                                #if there is only one reference barcode from the parent that has the maximum overlapwith the barcode and spacer at hand
                reference_bc <- NULL; reference_bc <- unique(alldata_bysample[[parent]][,1])[which(overlap == max(overlap))]
                if (max(overlap) > 6 ) {                                                    # and if the overlap is more than 7 bases in total.
                    ### correct the barcode in all loaded data
                    print(paste("truncated barcode correction:", bc, reference_bc, max(overlap)))
                    all_barcodes <- all_barcodes[-which(all_barcodes == bc)]
                    
                    for (pos in grep(bc, all_pairs)) { all_pairs[pos] <- paste(reference_bc, strsplit(all_pairs[pos], split = "-")[[1]][2] , sep = "-") }
                    
                    allspacers_bybarcode[[reference_bc]] <- sort(append(allspacers_bybarcode[[reference_bc]], allspacers_bybarcode[[bc]]))
                    allspacers_bybarcode[[bc]] <- NULL
                    
                    hgRNA_length <- hgRNA_length[-which(names(hgRNA_length) == bc)]
                    
                    for(sampl in names(alldata_bybarcode_bysample[[bc]])) {
                        alldata_bysample[[sampl]][which(alldata_bysample[[sampl]][,1] == bc),1] <- reference_bc
                        alldata_bysample[[sampl]] <- alldata_bysample[[sampl]][order(alldata_bysample[[sampl]]$V1),]
                        
                        alldata_bybarcode_bysample[[bc]][[sampl]][,1] <- reference_bc
                        alldata_bybarcode_bysample[[bc]][[sampl]][,ncol(alldata_bybarcode_bysample[[bc]][[sampl]])+1] <- bc
                        if(!is.null(alldata_bybarcode_bysample[[reference_bc]][[sampl]])) {
                            alldata_bybarcode_bysample[[reference_bc]][[sampl]][,ncol(alldata_bybarcode_bysample[[reference_bc]][[sampl]])+1] <- reference_bc
                            alldata_bybarcode_bysample[[reference_bc]][[sampl]] <- rbind(alldata_bybarcode_bysample[[reference_bc]][[sampl]], alldata_bybarcode_bysample[[bc]][[sampl]])        # For these samples the fifth column will be the original barcode
                        } else {
                            alldata_bybarcode_bysample[[reference_bc]][[sampl]] <- alldata_bybarcode_bysample[[bc]][[sampl]]
                        }
                        alldata_bybarcode_bysample[[bc]][[sampl]] <- NULL
                        alldata_bybarcode_bysample[[reference_bc]][[sampl]] <- alldata_bybarcode_bysample[[reference_bc]][[sampl]][order(alldata_bybarcode_bysample[[reference_bc]][[sampl]]$V3, decreasing = TRUE),]
                        alldata_bybarcode_bysample[[reference_bc]][[sampl]][,4] <-alldata_bybarcode_bysample[[reference_bc]][[sampl]][,3]/sum(alldata_bybarcode_bysample[[reference_bc]][[sampl]][,3]) * 100
                    }
                    
                }
            } 
        }
    }
}



# Filter 2: Fixing the barcodes that have suffered deletion and DO NOT have a close counterpart in the parental barcodes
all_barcodes -> all_barcodes_prefilter2
all_pairs -> all_pairs_prefilter2
allspacers_bybarcode -> allspacers_bybarcode_prefilter2
alldata_bysample  -> alldata_bysample_prefilter2
alldata_bybarcode_bysample -> alldata_bybarcode_bysample_prefilter2
hgRNA_length -> hgRNA_length_prefilter2

    # separating non-parental barcodes and sorting them based on number of samples they are observed in, and extracting a consensus (most prevalent) spacer for them.
nonparental_barcodes <- all_barcodes[!all_barcodes %in% alldata_bysample[[parent]][,1]]
nonparental_barcode_dominance <- NULL;
nonparental_barcodes_spacers <- list();
for (bc in nonparental_barcodes) {
    nonparental_barcode_dominance <- append(nonparental_barcode_dominance, length(names(alldata_bybarcode_bysample[[bc]])))
    its_spacers <- NULL; its_spacers <- unlist(lapply(alldata_bybarcode_bysample[[bc]], "[", 1, 2))
    nonparental_barcodes_spacers[[bc]] <- names(sort(table(its_spacers), decreasing = TRUE)[1])
}
nonparental_barcodes <- nonparental_barcodes[order(nonparental_barcode_dominance)]

for (bc in nonparental_barcodes) {
    other_barcodes <- NULL; other_barcodes <- all_barcodes[!all_barcodes %in% alldata_bysample[[parent]][,1]];
    other_barcodes <- other_barcodes[-which(other_barcodes == bc)]              # all non-parental barcodes with the barcode at hand removed
    MaxLargeDelLength <- -1 - (median(unlist(lapply(gregexpr(TSS_pattern1, alldata_bysample[[parent]][,2]), min))))
    if (hgRNA_length[bc] < MaxLargeDelLength) {                # If the calculated length suggests that the spacer could be a result of a large deletion
        EcoRI_loc <- NULL; EcoRI_loc <- median(unlist(lapply(gregexpr(Post_bc_pattern, lapply(alldata_bybarcode_bysample[[bc]], "[", 1,2 )), min)))            # Calculates the position of the EcoRI site that comes right after the barcode. If this spacer has a large deletion, the EcoRI site should be seen.
        if (EcoRI_loc > 5) {
            bc_spacer_count <- NULL; bc_spacer_count <- max(unlist(lapply(lapply(alldata_bybarcode_bysample[[bc]], "[", 1,c(2,3) ), "[", 2)))                     # Want to choose the spacer that was observed the most for this barcode as reference, and this command extract its observed count
            bc_spacer <- NULL; bc_spacer <- unlist(lapply(alldata_bybarcode_bysample[[bc]], "[", 1,c(2,3) ))[which(unlist(lapply(alldata_bybarcode_bysample[[bc]], "[", 1,c(2,3) )) == as.character(bc_spacer_count))-1]             # Extracting the reference spacer sequence based on the maximum count measured in the previous step. This is a really dirty way to do it.
            
            barcode_distances <- NULL; barcode_distances <- stringdist.c(bc, other_barcodes)        # Stores the distance of each parent barcode from the barcode at hand
            spacer_distances <- NULL
            for (bc2 in other_barcodes) {
                spacer_distances <- append(spacer_distances, stringdist.c(substr(bc_spacer, divergence_point, EcoRI_loc), substr(nonparental_barcodes_spacers[[bc2]], divergence_point, EcoRI_loc) ))
            }
            
            overlap <- barcode_length + (EcoRI_loc - divergence_point + 1) - barcode_distances - spacer_distances               # It stores the total overlap between the current barcode with a parental barcode, including their spacers' divergent sequence regions.  It should be noted, that a base at the junction of truncated spacer barcode may be counted twice as overlap if it matches that particular point of both, but this is ok - perhaps even desirable - given the observatin that deletions often times happen between two bases of the same identity.
            if (length(overlap[overlap == max(overlap)]) == 1) {                                #if there is only one reference nonparental barcode that has the maximum overlapwith the barcode and spacer at hand
                reference_bc <- NULL; reference_bc <- other_barcodes[which(overlap == max(overlap))]
                if (max(overlap) > 6 ) {                                                    # and if the overlap is more than 7 bases in total.
                    ### correct the barcode in all loaded data
                    print(paste("truncated barcode correction:", bc, reference_bc, max(overlap)))
                    all_barcodes <- all_barcodes[-which(all_barcodes == bc)]
                    
                    for (pos in grep(bc, all_pairs)) { all_pairs[pos] <- paste(reference_bc, strsplit(all_pairs[pos], split = "-")[[1]][2] , sep = "-") }
                    
                    allspacers_bybarcode[[reference_bc]] <- sort(append(allspacers_bybarcode[[reference_bc]], allspacers_bybarcode[[bc]]))
                    allspacers_bybarcode[[bc]] <- NULL
                    
                    hgRNA_length <- hgRNA_length[-which(names(hgRNA_length) == bc)]
                    
                    for(sampl in names(alldata_bybarcode_bysample[[bc]])) {
                        alldata_bysample[[sampl]][which(alldata_bysample[[sampl]][,1] == bc),1] <- reference_bc
                        alldata_bysample[[sampl]] <- alldata_bysample[[sampl]][order(alldata_bysample[[sampl]]$V1),]
                        
                        alldata_bybarcode_bysample[[bc]][[sampl]][,1] <- reference_bc
                        alldata_bybarcode_bysample[[bc]][[sampl]][,ncol(alldata_bybarcode_bysample[[bc]][[sampl]])+1] <- bc
                        if(!is.null(alldata_bybarcode_bysample[[reference_bc]][[sampl]])) {
                            alldata_bybarcode_bysample[[reference_bc]][[sampl]][,ncol(alldata_bybarcode_bysample[[reference_bc]][[sampl]])+1] <- reference_bc
                            alldata_bybarcode_bysample[[reference_bc]][[sampl]] <- rbind(alldata_bybarcode_bysample[[reference_bc]][[sampl]], alldata_bybarcode_bysample[[bc]][[sampl]])        # For these samples the fifth column will be the original barcode
                        } else {
                            alldata_bybarcode_bysample[[reference_bc]][[sampl]] <- alldata_bybarcode_bysample[[bc]][[sampl]]
                        }
                        alldata_bybarcode_bysample[[bc]][[sampl]] <- NULL
                        alldata_bybarcode_bysample[[reference_bc]][[sampl]] <- alldata_bybarcode_bysample[[reference_bc]][[sampl]][order(alldata_bybarcode_bysample[[reference_bc]][[sampl]]$V3, decreasing = TRUE),]
                        alldata_bybarcode_bysample[[reference_bc]][[sampl]][,4] <-alldata_bybarcode_bysample[[reference_bc]][[sampl]][,3]/sum(alldata_bybarcode_bysample[[reference_bc]][[sampl]][,3]) * 100
                    }
                    
                }
            } 
        }
    }
}

# Filter 3: manual fixing
if (TRUE) {
    all_barcodes -> all_barcodes_prefilter3
    all_pairs -> all_pairs_prefilter3
    allspacers_bybarcode -> allspacers_bybarcode_prefilter3
    alldata_bysample  -> alldata_bysample_prefilter3
    alldata_bybarcode_bysample -> alldata_bybarcode_bysample_prefilter3
    hgRNA_length -> hgRNA_length_prefilter3
    
        # 'GACCCTTCTT' is not spurious, It is most likely from GACCCTTCCT. Eitherway, since it is only seen in one sample, it shouldn't alter the results significantly. 
        # 'TTTCCACAAG' is not spurious, but I cannot determine with certainly which barcode it is derived from. Most likely from TCATCAGGCC. Eitherway, since it is only seen in one sample, it shouldn't alter the results significantly. 
    
    ### correct the manually curated truncated barcodes in all loaded data
    trunc_barcodes <- c('GACCCTTCTT')
    trunc_barcodes_refs <- c('GACCCTTCCT')
    for(i in 1:length(trunc_barcodes)) {
        bc <- reference_bc <- NULL
        bc <- trunc_barcodes[i]
        reference_bc <- trunc_barcodes_refs[i]
        ### correct the barcode in all loaded data
        print(paste("truncated barcode MANUAL correction (Filter 3):", bc, reference_bc))
        all_barcodes <- all_barcodes[-which(all_barcodes == bc)]
        
        for (pos in grep(bc, all_pairs)) { all_pairs[pos] <- paste(reference_bc, strsplit(all_pairs[pos], split = "-")[[1]][2] , sep = "-") }
        
        allspacers_bybarcode[[reference_bc]] <- sort(append(allspacers_bybarcode[[reference_bc]], allspacers_bybarcode[[bc]]))
        allspacers_bybarcode[[bc]] <- NULL
        
        hgRNA_length <- hgRNA_length[-which(names(hgRNA_length) == bc)]
        
        for(sampl in names(alldata_bybarcode_bysample[[bc]])) {
            alldata_bysample[[sampl]][which(alldata_bysample[[sampl]][,1] == bc),1] <- reference_bc
            alldata_bysample[[sampl]] <- alldata_bysample[[sampl]][order(alldata_bysample[[sampl]]$V1),]
            
            alldata_bybarcode_bysample[[bc]][[sampl]][,1] <- reference_bc
            alldata_bybarcode_bysample[[bc]][[sampl]][,ncol(alldata_bybarcode_bysample[[bc]][[sampl]])+1] <- bc
            if(!is.null(alldata_bybarcode_bysample[[reference_bc]][[sampl]])) {
                alldata_bybarcode_bysample[[reference_bc]][[sampl]][,ncol(alldata_bybarcode_bysample[[reference_bc]][[sampl]])+1] <- reference_bc
                if ( ncol(alldata_bybarcode_bysample[[reference_bc]][[sampl]]) > ncol(alldata_bybarcode_bysample[[bc]][[sampl]]) ) {alldata_bybarcode_bysample[[bc]][[sampl]][,ncol(alldata_bybarcode_bysample[[bc]][[sampl]])+1] <- reference_bc}
                alldata_bybarcode_bysample[[reference_bc]][[sampl]] <- rbind(alldata_bybarcode_bysample[[reference_bc]][[sampl]], alldata_bybarcode_bysample[[bc]][[sampl]])        # For these samples the fifth column will be the original barcode
            } else {
                alldata_bybarcode_bysample[[reference_bc]][[sampl]] <- alldata_bybarcode_bysample[[bc]][[sampl]]
            }
            alldata_bybarcode_bysample[[bc]][[sampl]] <- NULL
            alldata_bybarcode_bysample[[reference_bc]][[sampl]] <- alldata_bybarcode_bysample[[reference_bc]][[sampl]][order(alldata_bybarcode_bysample[[reference_bc]][[sampl]]$V3, decreasing = TRUE),]
            alldata_bybarcode_bysample[[reference_bc]][[sampl]][,4] <-alldata_bybarcode_bysample[[reference_bc]][[sampl]][,3]/sum(alldata_bybarcode_bysample[[reference_bc]][[sampl]][,3]) * 100
        }
    }
}

## Consolidating replicates by choosing one or combining them
    # This code does not work on P21 samples
    #If a sample has replicates, as indicated by them having the same name, this code combines those replicates by integrating their reads together and reconstructing the alldata_bysample and alldata_bybarcode_bysample lists
all_barcodes -> all_barcodes_prefilter4
all_pairs -> all_pairs_prefilter4
allspacers_bybarcode -> allspacers_bybarcode_prefilter4
alldata_bysample  -> alldata_bysample_prefilter4
alldata_bybarcode_bysample -> alldata_bybarcode_bysample_prefilter4
hgRNA_length -> hgRNA_length_prefilter4


samplenames_forthis <- names(alldata_bysample)
samplenames_forthis <- samplenames_forthis[-grep("P21", samplenames_forthis)]
samplenames_forthis <- samplenames_forthis[-grep("parent", samplenames_forthis)]
sample_basenames <- unique(unlist(lapply(lapply(strsplit(samplenames_forthis, "-"), "[", 1:3), paste, collapse = "-")))

for(basename in sample_basenames) {
    if (length(grep(basename, names(alldata_bysample))) > 1) {
        replicates <- names(alldata_bysample)[grep(basename, names(alldata_bysample))]
        sample_name <- NULL; sample_name <- paste(basename, "-From", length(replicates), "reps", sep = "")
        pairs <- alldata_bysample[[replicates[1]]];
        for(repl in replicates[2:length(replicates)]) {
            pairs <- rbind(pairs, alldata_bysample[[repl]])
        }
        pairs <- aggregate(formula = V3 ~ V1+V2, data = pairs, FUN = sum)
        pairs <- pairs[order(pairs$V1, -pairs$V3),]
        pairs[,4] <- NA
        for(bc in levels(as.factor(as.vector(pairs[,1])))) {
            pairs[,4][pairs[,1] == bc] <- round(pairs[,3][pairs[,1] == bc]/(sum(pairs[,3][pairs[,1] == bc]))*100, digits = 2)
            alldata_bybarcode_bysample[[bc]][[sample_name]] <- subset(pairs, V1 == bc)
            for(repl in replicates) {alldata_bybarcode_bysample[[bc]][[repl]] <- NULL}
        }
        pairs$V5 <- NA; pairs$V5 <- pairs[,3]/sum(pairs[,3])*100
        alldata_bysample[[sample_name]] <- pairs;
        for(repl in replicates) {alldata_bysample[[repl]] <- NULL}
    }
}


################################################################################################
############################# PROCESSING DATA WHILE MAKING OUTPUTS #############################
################################################################################################

#Writing output pair tables
for(sampl in names(alldata_bysample)) {
    write.table(alldata_bysample[[sampl]][1:3], file = paste(sampl, "_filteredpairs.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}




