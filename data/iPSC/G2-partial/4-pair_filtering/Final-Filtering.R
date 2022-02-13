# This script starts with the sequencing error-corrected files in 3-SP_err_correction and:
#     1) Filters out samples with low coverage (MinRd criterion)
#     2) Filters out IDs with low coverage (PFCOFF_bc_ref and PFCOFF_bc_exp criteria)
#     3) Filters out pairs with low read counts (PFCOFF_pair and PRCOFF_pair criteria)
#     4) Corrects known PCR or sequencing artifacts that are unique to MARC1 samples
#     5) Corrects template-switching during PCR
#     6) Corrects early-cycle PCR mutations that can make a non-mutated spacer appear mutated (max_dist_spacer criterion)
#     7) Corrects one-base displacements due to sequencing error
#     8) Removes spacers with short reads.
#     9) Corrects IDs (barcodes) that have been completely or partially deleted due to a large deletion (***See VARIABLES TO BE SELECTED FOR EACH ANALYSIS***).
#
# The output of the script is
#     1) A table of high-confidence pairs and their counts.
#     2) Mutated fracton of each identifier compared to the sequence in the Founder of that line.
#     3) An allele frequency table that could be readily clustered by samples.
#
# The code creates the following lists within it that can be readily used for further analysis (such as clustering, etc):
#     1) alldata_bybarcode_bysample
#     2) alldata_bysample
library(stringdist)

###################################################################################################################################
################## VARIABLES TO BE SELECTED FOR EACH ANALYSIS #####################################################################
###################################################################################################################################
# If the mouse under analysis is derived from the PB3 line, Founder should be set as "PB3"; if from the PB7 line, the variable should be set as "PB7"
# For EC96Lin12 samples, modify the number in brackets below to 1; otherwise, leave as is.
file.remove(list.files("./", pattern = "*_filteredpairs\\.txt"))
Founder <- c('EC96Lin12', 'EC163Lin11')[1]

# In some experiments, some errors in IDs cannot be resolved by automatic filtering and need to be manually accounted for (orphan barcodes).
# Run this code one time with no values in the following two vectors. Orphan barcodes will be printed in stdout.
# If you can identify the parents of the orphan barcodes, populate these vectors with pairs of orphan barcodes and their real parent barcode.
trunc_barcodes      <- c()          # trunc_barcodes <- c('[orphan_barcode_1]', '[orphan_barcode_2]', ...)
trunc_barcodes_refs <- c()          # trunc_barcodes_refs <- c('[parental_barcode_for_orphan_barcode_1]', '[parental_barcode_for_orphan_barcode_2]', ...)

MinRd <- 20                                   # Minimum total number of reads for a sample to be considered.
PFCOFF_bc_ref <- NULL                               # Read Percentage Cutoff is the the percentage of total reads in a sample that an identifier (bc) needs to have to be included. If not assigned here (is.null) this value will be assigned for each sample based on the total number of observed barcodes in it.
PFCOFF_bc_exp <- 2                            # If not assigned above, PFCOFF_bc_ref will be 1/(number-of-barcodes-in-sample)^PFCOFF_bc_exp . So, if a sample has 10 initially detected barcodes, those with less than 1/(10^PFCOFF_pair) of total read counts will be removed. We typically choose PFCOFF_bc_exp between 1.5 and 2. Higher values are more inclusive. Lower values should be used when there is cross contamination between samples.
PFCOFF_pair <- 0.1                            # Read fraction Cutoff is the the fraction of total reads in a barcode (ID) that a pair needs to have to be considered at all. Spacer sequences that correspond to PCR and sequencing errors tend to be infrequent, constituting less than 0.01 of all spacers for an ID (barcode)
PRCOFF_pair <- 2                                # Pairs with equal to or fewer than PRCOFF_pair reads in total will be removed.
max_dist_spacer <- 1;                                     # Since spacers for each sample are converted into a consensus separately, it is possible that a specific sample has the reference spacer in the parent only with a PCR or sequencing-based single point mutation which causes it to appear 100% mutated. Such a problem was observed with samples with smaller reads. As a result, the string.dist(method = "hamming") strategy will later in this script to count slight variants of the parental sequence as non-mutants.
barcode_length <- 10                                      # Length of identifier sequence
divergence_point <- 12                                    # Position in the spacer at which the sequence is expected to start to diverge
spurious_barcodes <- c("IDSEQUENCE")                      # The barcodes (IDs) that have been determined to be spurious should be added to this vector so that are removed from each sample. For example, if you find an ID is due to cross contamination from another sample, or would like to eliminate it from analysis for any reason, add it to this vector. "IDSEQUENCE" is just a placeholder

######################################################################################################################################################################################
################## IMPORTING REFERENCE hgRNA SEQUENCES AND CLASSIFICATIONS ###################################################################################################################################
######################################################################################################################################################################################
## Loading the mastertable
mastertable <- read.table(file = paste("./INUSE_AllPB-CellBarcodesMasterTable.txt", sep = ""));                   # Based on manuscript supplementary tables
rownames(mastertable) <- mastertable[,1]
mastertable$Number <- 1:nrow(mastertable)
hgRNA_length <- as.numeric(mastertable$Length); names(hgRNA_length) <- mastertable$Identifier                 # This vector will store the derived length of each hgRNA

## Loading barcode class table
barcode_class <- read.table(paste("./INUSE-Cellbarcode_classification.txt", sep = ""), colClasses=c("character", "character"), header = TRUE);        # Based on manuscript supplementary tables
rownames(barcode_class) <- barcode_class[,1]
colnames(barcode_class)[1] <- "category"
barcode_class[barcode_class[,2] == "inactive",1] <- 0
barcode_class[barcode_class[,2] == "slow",1] <- 1
barcode_class[barcode_class[,2] == "mid",1] <- 2
barcode_class[barcode_class[,2] == "fast",1] <- 3
barcode_class <- barcode_class[barcode_class$founder == Founder, c(1, 2)]


######################################################################################################################################################################################
################## LOADING AND FILTERING SAMPLE DATA #################################################################################################################################
######################################################################################################################################################################################

if (Founder == "EC96Lin12") {parent <- "EC96Lin12-PBGeno"} else {parent <- "EC163Lin11-PBGeno"}
PAM_pattern1 <- "GGGT" ; PAM_pattern1_offset <-  0        # This pattern will be used to extract the length of hgRNA. Offset determines the distance between the pattern and the actual point of interest, PAM here.
PAM_pattern2 <- "GGTT" ; PAM_pattern2_offset <-  -1       # This pattern will be used to extract the length of hgRNA. Offset determines the distance between the pattern and the actual point of interest, PAM here.Two patterns are being used in case there is a mutation in one patterns region. The two patterns were chosen based on minimal (ideally 0) overlap and non-occurance in the pre-barcode region of the gRNA.
TSS_pattern1 <- "CCGG" ; TSS_pattern1_offset <-  +3       # This pattern will be used to extract the length of hgRNA, it marks the sequence before TSS. Offset determines the distance between the pattern and the actual point of interest, TSS here.
Post_bc_pattern <- "GAATTC"                               # When this pattern is observed in a spacer, it is indicative of a large deletion. This EcoRI restriction site was used to clone the hgRNA backbone and is thus downstream of the hgRNA and its identifier.

files1 <- system('ls ../3-*correction/*_allpairs.txt', intern = TRUE)                  # These are the processed barcode-spacer pair counts for all samples.
# Removing the file that corresponds to the founder not being used.
if(Founder == "EC163Lin11" &&
   (length(grep("EC96Lin12-PBGeno", files1))) > 0) {
  files1 <- files1[-grep("EC96Lin12-PBGeno", files1)]
} else if (Founder == "EC96Lin12" &&
           (length(grep("EC163Lin11-PBGeno", files1))) > 0) {
  files1 <- files1[-grep("EC163Lin11-PBGeno", files1)]
}


alldata_bysample <- list()                                                     #All of the data separated by sample in each object within the list
alldata_bybarcode_bysample <- list()                                           #All of the data separated by barcode (i.e., gene) in each object within the list
all_pairs <- vector()                                                          # This vector stores a list of all pairs seen in all samples
allspacers_bybarcode <- list()                                                 # This list stores a list of all pairs seen for each barcode samples
all_barcodes <- vector()

##  Reading the truepairs file for all samples and storing them as a list of dataframes in alldata_bysample
parent_sample <- NULL
name_base <- paste(Sys.Date(), '_SpacerFilteringLog', sep = "")
filename_text <- paste(name_base, '.txt', sep = '')
sink(file = filename_text); sink();
for (file1 in files1) {
  message(file1)
  sample_name <- tail(unlist(strsplit(file1, split = "/")), n=1); sample_name <- unlist(strsplit(sample_name, split = "_"))[1];
  if(gregexpr(parent, file1) >= 0) {parent_sample <- sample_name}
  print(paste("processing", sample_name, Sys.time(), sep = " "))
  sink(file = filename_text, append = TRUE); cat(paste("processing", sample_name, Sys.time(), "\n", sep = " ")); sink()
  
  pairs <- NULL;
  pairs <- read.table(file1, colClasses=c("character", "character", "numeric"))                          # Reading all pairs of the sample
  
  ## Applying length-based filters    (in some rare samples short spacers have been observed which probably means a large deletion extending to the primer)
  if( max(nchar(pairs[,2])) != min(nchar(pairs[,2]))) {
    print(unname(cbind("Removed", pairs[nchar(pairs[,2]) < max(nchar(pairs[,2])),], "because of short spacer")))
    pairs <- pairs[nchar(pairs[,2]) == max(nchar(pairs[,2])),]
  }
  
  ## Applying readcount-based filters
  if (nrow(pairs) < 1 ) { print(paste(sample_name, "skipped due to having no data")); next; }
  pairs <- subset(pairs, !V1 %in% spurious_barcodes)
  if (sum(pairs[,3]) < MinRd) { cat(paste("\t\t\t", sample_name, "skipped due to having less than MinRd reads", "\n")); next; }
  
  ## Correcting template switching for all samples (the rational is that if the exact same spacer is seen in two barcodes within the same sample, it belongs to only one of those barcodes and that is the barcode where that spacer has a higher count)
  byspacers <- NULL; byspacers <- aggregate(V3~V2, pairs, max)
  for (sp in 1:nrow(byspacers)) {
    pairs <- subset(pairs, V2 != byspacers[sp,1] | (V2 == byspacers[sp,1] & V3 > byspacers[sp,2]/20) )        # For each spacer sequence, removing all instances, regardless of barcode, that have a submaximum frequency in the sample (since template switching can only shuffle other spacers in the same sample). However, large deletions can lead to the same legitimate spacer in different identifiers. Therefore, to avoid removing these, instead of keeping only the most frequent instance of each spacer, only those with less than 20X that frequency will be removed. Template switching generally leads to spacers at less than 1%, so 20X or 5% is a reasonable cutoff. The effectiveness of this approach (after 190716) was confirmed on bulkPCR-A02-0003pgul_truepairs.txt sample, indicating at least three legitimate spacers that had been deleted by the previous version of this filter.
  }
  
  bc_abundances <- aggregate(V3 ~ V1, pairs, sum);
  for (i in 1:nrow(bc_abundances)) { if (bc_abundances[i,2] <= PRCOFF_pair) {pairs <- pairs[-which(pairs[,1] == bc_abundances[i,1], arr.ind = TRUE),]} }          # Eliminating pairs that don't meat the PRCOFF_pair criteria.
  
  if (nrow(pairs) == 0) {
    next
  }
  
  pairs$V4 <- NA;                                     # For each entry, this column indicates its percentage among pairs with the same barcode in that sample
  if(is.null(PFCOFF_bc_ref)) {PFCOFF_bc <- 1/length(levels(as.factor(as.vector(pairs[,1]))))^PFCOFF_bc_exp} else {PFCOFF_bc <- PFCOFF_bc_ref} #Establishing a cutoff for removing IDs/barcodes with low total read counts. This cutoff is PFCOFF_bc and is calculated based on PFCOFF_bc_ref and PFCOFF_bc_exp. These IDs, which can be a result of cross contamination, are present at lower frequency than IDs that truly belong to the sample. This filtering criteria carries out an adjustment based on total number of IDs so that the cutoff is lower the higher the total number of IDs is.
  
  ## [200305] : The "pairs" table has to be sorted properly for below steps to work appropriately (especially the one-base displacement adjustment). Therefore, the following command is added to make sure
  pairs <- pairs[order(pairs$V1, -pairs$V3),]   #sorting pairs on barcode followed by descending frequency because some filtering operations above can disrupt order
  
  for(bc in levels(as.factor(as.vector(pairs[,1])))) {
    
    if ( sum(subset(pairs, V1 == bc)[,3]) / sum(pairs[,3]) >= PFCOFF_bc) {
      if (length(which(pairs[,3] <= PRCOFF_pair & pairs[,1] == bc, arr.ind = TRUE)) > 0 ) {pairs <- pairs[-which(pairs[,3] <= PRCOFF_pair & pairs[,1] == bc, arr.ind = TRUE),]}           # Eliminating pairs that don't meat the PRCOFF_pair criteria
      
      ## [180608] Applying max_dist_spacer criteria
      if (nrow(subset(pairs, V1 == bc))>1) {
        pairs.bc = subset(pairs, V1 == bc)
        for (sp in nrow(pairs.bc):2) {
          if ( min(stringdist(pairs.bc[sp,2], pairs.bc[(1:(sp-1)),2], method = "hamming")) <= max_dist_spacer ) {
            reference_sp <- which(stringdist(pairs.bc[sp,2], pairs.bc[(1:(sp-1)),2], method = "hamming") == min(stringdist(pairs.bc[sp,2], pairs.bc[(1:(sp-1)),2], method = "hamming")))[1]
            
            # Establishing that the difference is not near cut site
            PAM_pos1 <- min(gregexpr(PAM_pattern1, pairs.bc[sp,2])[[1]])
            PAM_pos2 <- max(gregexpr(PAM_pattern1, pairs.bc[sp,2])[[1]])
            spliced_spacer <- paste( substr(pairs.bc[sp,2], start = 0, stop = PAM_pos1-8) ,
                                     substr(pairs.bc[sp,2], start = PAM_pos2-3, stop = nchar(pairs.bc[sp,2])) , sep = "-")
            spliced_reference <- paste( substr(pairs.bc[reference_sp,2], start = 0, stop = PAM_pos1-8) ,
                                        substr(pairs.bc[reference_sp,2], start = PAM_pos2-3, stop = nchar(pairs.bc[reference_sp,2])) , sep = "-")
            if (stringdist(spliced_spacer, spliced_reference, method = "hamming") > 0) {             # Under this condition, the detected small variation is outside of the cutrange
              sink(file = filename_text, append = TRUE); cat(paste("Combined ", paste(pairs.bc[sp,1:3], collapse = " "), " into ", paste(subset(pairs, V1 == bc)[reference_sp,1:3], collapse = " "), "\tfor max_spacer_dist\n", sep = "")); sink();
              pairs[which(pairs[,2] == pairs.bc[reference_sp,2]),3] <- pairs.bc[reference_sp,3] + pairs.bc[sp,3]
              pairs <- pairs[-which(pairs[,2] == subset(pairs, V1 == bc)[sp,2]),]
            }
          }
        }
      }
      
      # [180315] module for consolidating one-base displacements due to sequencing error, such as "GAAACACCGGTAGCAAACGTTTGGACGTGGGGTTAGAGCTAGAAATAGCA" and "CAAAACACCGGTAGCAAACGTTTGGACGTGGGGTTAGAGCTAGAAATAGC" that create an erroneous mutant spacer
      if (nrow(subset(pairs, V1 == bc)) > 1) {
        for (spacer_obs in pairs.bc$V2[length(pairs.bc$V2):2]) {
          spacer_obs_trunc <- substr(spacer_obs, 5, nchar(spacer_obs))
          matching_spacers <- grep(spacer_obs_trunc, pairs$V2);
          matching_spacers <- matching_spacers[matching_spacers %in% which(pairs$V1 == bc)]
          if(length(matching_spacers) > 1 ) {                     # Then the first one would be the parent and last one would be spacer_obs
            sink(file = filename_text, append = TRUE); cat(paste("Combined ", paste(pairs[matching_spacers[length(matching_spacers)],1:3], collapse = " "), " into ", paste(pairs[matching_spacers[1],1:3], collapse = " "), "\tfor one base displacement\n", sep = "")); sink();
            pairs[matching_spacers[1],3] <- pairs[matching_spacers[1],3] + pairs[matching_spacers[length(matching_spacers)],3]
            pairs <- pairs[-matching_spacers[length(matching_spacers)],]
          }
        }
        
      }
      pairs <- pairs[order(pairs$V1, -pairs$V3),]                 # Re-ordering pairs since the above two consolidations can change the abundance numbers and its order
      
      cutoff <- sum(pairs[,3][pairs[,1] == bc]) * PFCOFF_pair
      if (length(which(pairs[,3] < cutoff & pairs[,1] == bc, arr.ind = TRUE)) > 0 ) {   # Eliminating pairs that don't meet the PFCOFF_pair criteria
        pairs <- pairs[-which(pairs[,3] < cutoff & pairs[,1] == bc, arr.ind = TRUE),]}
      pairs[,4][pairs[,1] == bc] <- round(pairs[,3][pairs[,1] == bc]/(sum(pairs[,3][pairs[,1] == bc]))*100, digits = 2)
      if ( nrow(subset(pairs, V1 == bc)) > 0 ) {alldata_bybarcode_bysample[[bc]][[sample_name]] <- subset(pairs, V1 == bc)}
      allspacers_bybarcode[[bc]] <- append(allspacers_bybarcode[[bc]], as.character(pairs[,2][pairs[,1] == bc])); allspacers_bybarcode[[bc]] <- unique(allspacers_bybarcode[[bc]])
    } else {pairs <- subset(pairs, V1 != bc)}
  }
  
  if (nrow(pairs) == 0) {
    next
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

if (is.null(parent_sample)) {       #Printing an error message that no founder sample has been identified
  print("ERROR: No founder sample was identified in ../3-SP_err_correction. Therefore, mutation rates cannot be calculated.")
  print("A founder sample file named 'EC163Lin11-PBGeno_allpairs.txt' or 'EC96Lin12-PBGeno_allpairs.txt' must be present in ../3-SP_err_correction")
}

### Finding the IDs in each sample that were not observed in the Founder/Parent.
Orph1 <- NULL;
for(sampl in names(alldata_bysample))  {
  temp <- length(which(!unique(alldata_bysample[[sampl]][,1]) %in% alldata_bysample[[parent_sample]][,1]));
  Orph1 <- append(Orph1, temp)
}
names(Orph1) <- names(alldata_bysample)

if (sum(Orph1) > 0) {
  for(sampl in names(Orph1)[Orph1 > 0])  {
    nonparental_barcodes <- unique(alldata_bysample[[sampl]][,1])[!unique(alldata_bysample[[sampl]][,1]) %in% alldata_bysample[[parent_sample]][,1]];
    #print(paste("The following identifiers (barcodes) in ", sampl, " do not exist in ", parent, " (the Founder/Parent): ", paste(nonparental_barcodes, collapse = ","), sep = " "))
  }
} else {
  print("All IDs in all samples match the IDs in the Founder. No additional adjustment is necessary")
}

######################################################################################################################################################################################
################## Filtering and adjusting IDs (barcodes) that are not observed in the founder/parent ################################################################################
######################################################################################################################################################################################
# All IDs are expected to descent from a parental ID that was present in the line's founder.
# However, large deletion events can create hgRNAs that appear to have a new identifier sequence.
# The following section of the code addresses these identifiers, attempting to merge them with their original parental identifier in the sample.

# Custom string.dist function for 5' aligned measurement of distance based on hamming. First mismatch, will be the end of agreement between the strings.
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

### Barcode Corrections: Fixing the barcodes that have suffered a large deletion which has taken out most of the scaffold
if (sum(Orph1) > 0) {
  # Filter 1: Fixing the barcodes that have suffered deletion and have a close counterpart in the parental barcodes
  all_barcodes -> all_barcodes_prefilter1
  all_pairs -> all_pairs_prefilter1
  allspacers_bybarcode -> allspacers_bybarcode_prefilter1
  alldata_bysample  -> alldata_bysample_prefilter1
  alldata_bybarcode_bysample -> alldata_bybarcode_bysample_prefilter1
  hgRNA_length -> hgRNA_length_prefilter1
  
  MaxLargeDelLength <- 10
  for (bc in all_barcodes[!all_barcodes %in% alldata_bysample[[parent]][,1]]) {
    EcoRI_loc <- median(unlist(lapply(gregexpr(Post_bc_pattern, lapply(alldata_bybarcode_bysample[[bc]], "[", 1,2 )), min)))            # Calculates the position of the EcoRI site that comes right after the barcode. If this spacer has a large deletion, the EcoRI site should be seen.
    if (EcoRI_loc > 5) {
      bc_spacer_count <- max(unlist(lapply(lapply(alldata_bybarcode_bysample[[bc]], "[", 1,c(2,3) ), "[", 2)))                     # Want to choose the spacer that was observed the most for this barcode as reference, and this command extract its observed count
      bc_spacer <- unlist(lapply(alldata_bybarcode_bysample[[bc]], "[", 1,c(2,3) ))[which(unlist(lapply(alldata_bybarcode_bysample[[bc]], "[", 1,c(2,3) )) == as.character(bc_spacer_count))-1]             # Extracting the reference spacer sequence based on the maximum count measured in the previous step. This is a really dirty way to do it.
      
      barcode_distances <- stringdist.c(bc, unique(alldata_bysample[[parent]][,1]))        # Stores the distance of each parent barcode from the barcode at hand
      spacer_distances <- NULL
      for (bc2 in unique(alldata_bysample[[parent]][,1])) {
        spacer_distances <- append(spacer_distances,
                                   stringdist.c(substr(bc_spacer, divergence_point, EcoRI_loc),
                                                substr(alldata_bybarcode_bysample[[bc2]][[parent]][1,2],
                                                       divergence_point, EcoRI_loc) ))
      }
      
      overlap <- barcode_length + (EcoRI_loc - divergence_point + 1) - barcode_distances - spacer_distances               # It stores the total overlap between the current barcode with a parental barcode, including their spacers' divergent sequence regions.  It should be noted, that a base at the junction of truncated spacer barcode may be counted twice as overlap if it matches that particular point of both, but this is ok - perhaps even desirable - given the observatin that deletions often times happen between two bases of the same identity.
      if (length(overlap[overlap == max(overlap)]) == 1) {                                #if there is only one reference barcode from the parent that has the maximum overlapwith the barcode and spacer at hand
        reference_bc <- unique(alldata_bysample[[parent]][,1])[which(overlap == max(overlap))]
        if (max(overlap) >= 5 ) {                                                    # and if the overlap is more than 7 bases in total.
          ### correct the barcode in all loaded data
          print(paste("truncated barcode correction (Filter1):", bc, reference_bc, max(overlap)))
          all_barcodes <- all_barcodes[-which(all_barcodes == bc)]
          
          for (pos in grep(bc, all_pairs)) { all_pairs[pos] <- paste(reference_bc, strsplit(all_pairs[pos], split = "-")[[1]][2] , sep = "-") }
          
          allspacers_bybarcode[[reference_bc]] <- sort(append(allspacers_bybarcode[[reference_bc]], allspacers_bybarcode[[bc]]))
          allspacers_bybarcode[[bc]] <- NULL
          
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
            alldata_bybarcode_bysample[[bc]][[sampl]] <- NULL                           # Removing that sample from the barcode
            alldata_bybarcode_bysample[[reference_bc]][[sampl]] <- alldata_bybarcode_bysample[[reference_bc]][[sampl]][order(alldata_bybarcode_bysample[[reference_bc]][[sampl]]$V3, decreasing = TRUE),]
            alldata_bybarcode_bysample[[reference_bc]][[sampl]][,4] <-alldata_bybarcode_bysample[[reference_bc]][[sampl]][,3]/sum(alldata_bybarcode_bysample[[reference_bc]][[sampl]][,3]) * 100
          }
          
          if ( length(alldata_bybarcode_bysample[[bc]]) == 0 ) {alldata_bybarcode_bysample[[bc]] <- NULL}         # If all samples within a barcode have been removed in the previous step, the barcode itself should be removed.
          
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
  if (length(nonparental_barcodes) > 0) {
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
          if (max(overlap) >= 5 ) {                                                    # and if the overlap is more than 7 bases in total.
            ### correct the barcode in all loaded data
            print(paste("truncated barcode correction (Filter 2):", bc, reference_bc, max(overlap)))
            all_barcodes <- all_barcodes[-which(all_barcodes == bc)]
            
            for (pos in grep(bc, all_pairs)) { all_pairs[pos] <- paste(reference_bc, strsplit(all_pairs[pos], split = "-")[[1]][2] , sep = "-") }
            
            allspacers_bybarcode[[reference_bc]] <- sort(append(allspacers_bybarcode[[reference_bc]], allspacers_bybarcode[[bc]]))
            allspacers_bybarcode[[bc]] <- NULL
            
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
}

### Finding the IDs in each sample that are not observed in the Founder/Parent and were not accounted for with Filter 1 and Filter 2 above.
Orph2 <- NULL;
for(sampl in names(alldata_bysample))  {
  temp <- length(which(!unique(alldata_bysample[[sampl]][,1]) %in% alldata_bysample[[parent]][,1]));
  Orph2 <- append(Orph2, temp)
}
names(Orph2) <- names(alldata_bysample)

if (sum(Orph2) > 0) {
  orphan_barcodes <- NULL;
  for(sampl in names(Orph2)[Orph2 > 0])  {
    nonparental_barcodes <- unique(alldata_bysample[[sampl]][,1])[!unique(alldata_bysample[[sampl]][,1]) %in% alldata_bysample[[parent_sample]][,1]];
    # print(paste("The following identifiers (barcodes) in ", sampl, " do not exist in ", parent, " (the Founder/Parent) and have not been corrected using automated filters: ", paste(nonparental_barcodes, collapse = ","), sep = " "))
    orphan_barcodes <- append(orphan_barcodes, nonparental_barcodes)
  }
  # print("*********Above IDs have to be manually corrected*********")
  orphan_barcodes <- unique(orphan_barcodes)
} else {
  if (sum(Orph1) > 0) {print("After automated filtering, all IDs in all samples match the IDs in the Founder. No additional adjustment is necessary")}
}

if (sum(Orph2) > 0) {
  # Filter 3: manual fixing
  all_barcodes -> all_barcodes_prefilter3
  all_pairs -> all_pairs_prefilter3
  allspacers_bybarcode -> allspacers_bybarcode_prefilter3
  alldata_bysample  -> alldata_bysample_prefilter3
  alldata_bybarcode_bysample -> alldata_bybarcode_bysample_prefilter3
  hgRNA_length -> hgRNA_length_prefilter3
  
  ### correct the manually curated truncated barcodes in all loaded data
  if ( length(trunc_barcodes) > 0 ) {
    for(i in 1:length(trunc_barcodes)) {
      bc <- trunc_barcodes[i]
      reference_bc <- trunc_barcodes_refs[i]
      ### correct the barcode in all loaded data
      if (!bc %in% all_barcodes) {next} else {
        all_barcodes <- all_barcodes[-which(all_barcodes == bc)]
        for (pos in grep(bc, all_pairs)) { all_pairs[pos] <- paste(reference_bc, strsplit(all_pairs[pos], split = "-")[[1]][2] , sep = "-") }
        
        allspacers_bybarcode[[reference_bc]] <- sort(append(allspacers_bybarcode[[reference_bc]], allspacers_bybarcode[[bc]]))
        allspacers_bybarcode[[bc]] <- NULL
        
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
  }
}

Orph3 <- NULL;
for(sampl in names(alldata_bysample))  {
  temp <- length(which(!unique(alldata_bysample[[sampl]][,1]) %in% alldata_bysample[[parent]][,1]));
  Orph3 <- append(Orph3, temp)
}
names(Orph3) <- names(alldata_bysample)

name_base_err <- paste(Sys.Date(), '_ERROR', sep = "")        # Reporting orphan barcodes to an error file output.
filename_text_err <- paste(name_base_err, '.txt', sep = '')
if (sum(Orph3) > 0) {
  orphan_barcodes <- NULL;
  for(sampl in names(Orph3)[Orph3 > 0])  {
    nonparental_barcodes <- unique(alldata_bysample[[sampl]][,1])[!unique(alldata_bysample[[sampl]][,1]) %in% alldata_bysample[[parent_sample]][,1]];
    print(paste("The following identifiers (barcodes) in ", sampl, " do not exist in ", parent, " (the Founder/Parent) and have not been corrected using automated or manual filters: ", paste(nonparental_barcodes, collapse = ","), sep = " "))
    sink(file = filename_text_err, append = TRUE);
    cat(paste("The following identifiers (barcodes) in ", sampl, " do not exist in ", parent, " (the Founder/Parent) and have not been corrected using automated or manual filters: ", paste(nonparental_barcodes, collapse = ","), sep = " "), "\n")
    sink();
  }
  print("*********Above IDs are orphan and have to be manually corrected. Enter each orphan ID as an element 'trunc_barcodes' in line 32 of this code.*********")
  print("*********Then enter its corresponding reference ID in the corresponding element in 'trunc_barcodes_ref' in line 33.*********")
  print("*********Finally, return to line 417 of this code and rerun 'Filter 3: manual fixing'*********")
  print("*********Warnings: IDs that do not match an ID in the Founder will be excluded from Barcode Table*********")
  sink(file = filename_text_err, append = TRUE);
  cat("Above IDs are orphan and have to be manually corrected. Enter each orphan ID as an element 'trunc_barcodes' in line 32 of this code.", "\n")
  cat("Then enter its corresponding reference ID in the corresponding element in 'trunc_barcodes_ref' in line 33.", "\n")
  cat("Finally, return to line 417 of this code and rerun 'Filter 3: manual fixing'", "\n")
  cat("Warnings: IDs that do not match an ID in the Founder will be excluded from Barcode Table", "\n", "\n", "\n")
  sink();
} else {
  if ( sum(Orph1) > 0 & sum(Orph2) > 0 ) {print("After automated and manual filtering, all IDs in all samples match the IDs in the Founder. No additional adjustment is necessary")}
}

### Writing the filtered pairs into output files
for(sampl in names(alldata_bysample)) {
  table <- NULL; table <- alldata_bysample[[sampl]]
  table[,4] <- round(table[,4], digits = 2)
  table[,5] <- round(table[,5], digits = 3)
  write.table(table, paste(sampl, "_filteredpairs.txt", sep = ''), quote = FALSE , row.names = FALSE, col.names = FALSE, sep = "\t" )
}

#################################################################################################################################################################
################## REPORTING SOME STATISTICS AND GENERAL FEATURES OF THE SAMPLES ################################################################################
#################################################################################################################################################################

######## BREAKDOWN OF BARCODES IN EACH SAMPLE AND EACH MOUSE ########
## Constructing a matrix that determines which barcodes have been observed in which samples. Rows are samples, columns are barcode, 1 denotes presence and 0 denotes absence of a particular barcode in a sample.
barcodesample_matrix <- as.data.frame(matrix(nrow = length(alldata_bysample), ncol = length(all_barcodes), dimnames = list(names(alldata_bysample), all_barcodes)))
barcodesample_matrix[] <- 0
for(bc in all_barcodes) {
  barcodesample_matrix[rownames(barcodesample_matrix) %in% names(alldata_bybarcode_bysample[[bc]]), colnames(barcodesample_matrix) == bc] <- 1
}
print_barcodesample_matrix <- barcodesample_matrix[grep("parent", rownames(barcodesample_matrix), invert = TRUE),]
print_barcodesample_matrix <- print_barcodesample_matrix[,colSums(print_barcodesample_matrix)>0]
write.table(print_barcodesample_matrix, file = paste(Sys.Date(), "_IdentifierSampleMatrix.txt", sep = ""), quote = FALSE, sep = "\t")


######## BREAKDOWN OF MUTATION RATES IN EACH SAMPLE FOR EACH hgRNA########

### Calculating the mutation level of each sample for each barcode
mutations <- list()               # This list stores mutation rate of each hgRNA, first broken down by sample and then by barcode.
samples <- names(alldata_bysample)

for (sampl in names(alldata_bysample)) {
  mutations[[sampl]] <- list()
  for (bc in unique(alldata_bysample[[sampl]][,1])) {
    partial_rate <- NULL;
    this_table <- NULL; this_table <- alldata_bybarcode_bysample[[bc]][[sampl]]                         # This will enable carrying out the operations while taking out every spacer that has been associated with a parent_sample (wt) spacer already.
    if (!is.null(alldata_bybarcode_bysample[[bc]][[parent]]) ) {                                         # This is to skip  barcodes that have only been seen in the samples but not in the parent.
      for (wt in 1:nrow(alldata_bybarcode_bysample[[bc]][[parent]])) {
        partial_rate <- append(partial_rate, round(sum(this_table[,4][stringdist(this_table[,2], alldata_bybarcode_bysample[[bc]][[parent_sample]][wt,2], method = "hamming") <= max_dist_spacer]), digit = 2))
        this_table <- this_table[-which(stringdist(this_table[,2], alldata_bybarcode_bysample[[bc]][[parent_sample]][wt,2], method = "hamming") <= max_dist_spacer),]              # This line removes the spacers that were analyzed, should the parent_sample  have multiple wt spacers, so that one spacer in sampl won't get counted multiple times, each time against a different parent_sample spacer.
      }
      mutations[[sampl]] <- append(mutations[[sampl]], 100 - sum(partial_rate))
      names(mutations[[sampl]])[length(mutations[[sampl]])] <- bc
    }
  }
}

### Printing all mutation levels
name_base <- paste(Sys.Date(), '_hgRNA-mutation-levels.txt', sep = "")
sink(file = name_base);
cat("sample", "\t", "ID", "\t", "class", "\t", "%mutated", "\n")
for(sampl in names(mutations)) {
  for (bc in names(mutations[[sampl]])) {
    cat(sampl, "\t", bc, "\t", barcode_class[bc,2], "\t", mutations[[sampl]][[bc]], "\n")
  }
}
sink()


#***** DID NOT RUN BELOW YET *******************

################## MAKING BARCODE TABLE ###################################################################################################################################
mouse_samples <- c()                  # If wishing to limit the barcode table to a subset of samples, populate this variable with full sample names. Otherwise, leave it empty to include all samples.
# name_pattern <- "";
# mouse_samples <- names(alldata_bysample)[grep(name_pattern, names(alldata_bysample))]

if (is.null(mouse_samples))
  mouse_samples <- names(alldata_bysample)

file1name <- paste(Sys.Date(), "_", 'BarcodeTableAlleleLookup', '.txt', sep = "")
sink(file = file1name); sink();              # Since append will be used later, this makes sure the file is wiped clean every time, if it exists.
names(mouse_samples) <- mouse_samples
mouse_samples <- mouse_samples[(order(names(mouse_samples)))]

barcodes <- NULL;       # Identifiers to be included in Barcode Table

for (mouse_sampl in mouse_samples){
  if (mouse_sampl != parent_sample)
    barcodes <- union(barcodes, alldata_bysample[[mouse_sampl]][,1])
}

barcodes <- barcodes[barcodes %in% alldata_bysample[[parent_sample]][, 1]] # exluding non-parental (orphan barcodes)
barcodes <- sort(barcodes) # Sorting Identifiers/IDs alphabetically (which matches their numbering).

master_barcode_table <- NULL; master_barcode_table <- as.data.frame(matrix(nrow = length(mouse_samples), ncol = 0, dimnames = list( mouse_samples, NULL ) ))               # Table of all barcodes. Each row will be a sample and each column will be a mutant allele. The value in each cell shows the frequency of the corresponding mutant allele in the corresponding sample.
for (bc in barcodes) {
  alldata_bysample <- alldata_bybarcode_bysample[[bc]]
  #deriving parental and all mutant spacer sequences for this bc in relevant samples
  spacers <- NULL; spacers <- alldata_bysample[[parent]][,2]
  for (mouse_sampl in mouse_samples) {
    spacers <- union(spacers, alldata_bysample[[mouse_sampl]][,2])
  }
  parental_spacer_count <- NULL; parental_spacer_count <- nrow(alldata_bysample[[parent]])      # Number of mutant alleles observed for the ID in the founder (expected to be one, but common sequencing errors may lead to more).
  names(spacers) <- paste("ID#", formatC(mastertable[bc,]$Number, width = 3, flag = "0"),"-mut#", formatC(seq(1-parental_spacer_count, 1-parental_spacer_count+length(spacers)-1), width = nchar(length(spacers)), flag = "0"), sep ="" )
  names(spacers)[1:parental_spacer_count] <- paste("ID#", formatC(mastertable[bc,]$Number, width = 3, flag = "0"), "-par#", formatC(1:parental_spacer_count, width = nchar(length(spacers)), flag = "0"), sep = "")
  write.table(as.data.frame(spacers), file = file1name, append = TRUE, quote = FALSE, sep = "\t",  col.names = FALSE)
  
  barcode_table_this_ID <- NULL; barcode_table_this_ID <- as.data.frame(matrix(nrow = length(mouse_samples), ncol = length(spacers), dimnames = list( mouse_samples, names(spacers) ) )) #Barcode table for the identifier currently under analysis. These barcode tables will be combined to obtain the master_barcode_table
  for(sampl_i in 1:length(mouse_samples)) {
    mouse_sampl <- mouse_samples[sampl_i]
    if ( !is.null(nrow(alldata_bysample[[mouse_sampl]])) ) {
      barcode_table_this_ID[sampl_i, ] <- 0                              # if an identifier/hgRNA is seen in a sample, then the frequency of all its unseen mutant alleles in that sample is 0. If the Identifier/hgRNA is absent from a sample altogether, then the frequency of all its mutant alleles is NA.
      for (spacr_r in 1:nrow(alldata_bysample[[mouse_sampl]]) ) {              # instead of going through all alleles for each sampl, going through only those that exist in that sampl
        barcode_table_this_ID[sampl_i, which(spacers == alldata_bysample[[mouse_sampl]][spacr_r,2])] <- alldata_bysample[[mouse_sampl]][spacr_r, 4]
      }
    }
  }
  master_barcode_table <- cbind(master_barcode_table, barcode_table_this_ID);
}

sink(file = file1name, append = TRUE); cat("\n\n\n"); sink();

file2name <- paste(Sys.Date(), "_", 'BarcodeTable', '.txt', sep = "")
write.table(master_barcode_table, file = file2name, append = FALSE, sep = "\t")

barcode_table <- master_barcode_table
barcode_table <- barcode_table[,-grep("par", colnames(barcode_table))]   # Removing the parental alleles from the table.
#This barcode_table can be used for clustering, for example:
#dendrogram <- as.dendrogram(hclust(dist(barcode_table, method = "manhattan"), method = "ward.D2"))
#plot(dendrogram)

