# Script to detect sequencing errors (single-base mismatches) in spacers and comnbining spacers that otherwise have a the same sequence.
# This script will create 3 output files for each input sample in ../1-pair_counting:
#                       1. A [sample]_ambigpairs.txt which contains all pairs from a sample whose ID or spacers have been identified as having sequencing error using the script by Beltman et al, BMC bioinformatics, 2016 below but whose correct spacer or ID cannot be identified.
#                       2. A [sample]_truepairs.txt which contains all pairs from a sample whose ID corresponds to one of the main IDs in ../2-ID_err_correction and whose spacers have been adjusted for sequencing errors.
#                       3. A [sample]_allpairs.txt is the combination of [sample]_truepairs.txt and [sample]_ambigpairs.txt.
#                       3. A [sample]_genotypes.txt describes the IDs in each sample together with their most prominent spacer.
#               Use [sample]_allpairs.txt if using the filtering strategies in the following steps of the pipelne. Use [sample]_truepairs.txt otherwise. [sample]_ambigpairs.txt has no application but record-keeping.
#               Use [sample]_genotypes.txt for purposes of genotyping and deciding the crosses for creating next generations. [sample]_genotypes.txt should not be used for barcoded samples (Cas9+) because it does not contain all the spacers, only the most frequent one.

library(stringdist)
library(VGAM)

#########################################################################################
# PART 1:   This first part of the code takes various tables and sends them to the second part of the code.
#########################################################################################

max_dist_barcode <- 2;                          # Maximum lv distance allowed between sequenced and consensus barcodes for data point to be kept
max_dist_spacer <- 5;                           # Maximum Hamming distance allowed between sequenced and consensus spacers for data point to be kept
write.individual.tables <- FALSE                # if TRUE, all true and error spacers will be written for each spacer in the ./intermediate_files folder.


files1 <- system('ls ../1-*/*_barcode_gRNA_pairing_counts.txt', intern = TRUE)                  # These are the raw barcode-spacer pair counts

for (file1 in files1[351:1104]) {
    if (file.info(file1)[1] == 0) { print(paste(file1, " is empty*******", sep = "")); next }
	if (file.info(file1)[1] == 0) next
    sample_name <- NULL; sample_name <- tail(unlist(strsplit(file1, split = "/")), n=1); sample_name <- unlist(strsplit(sample_name, split = "_"))[1]
    print(paste("processing", sample_name, Sys.time(), sep = " "))
    true_barcodes_file <- NULL; true_barcodes_file <- paste( "../2-ID_err_correction/", sample_name, "_trueID.txt", sep = "" )      # This is the processed barcodes corresponding to sample_name

    raw_pairs <- NULL; raw_pairs <- read.table(file1, colClasses=c("character", "character", "numeric"))                                # If ColClasses are not defined, the barcodes and spacers will be treated as factors down below which will slow the program to a screeching halt.
    true_barcodes <- NULL; true_barcodes <- read.table(true_barcodes_file, colClasses=c("character", "numeric"))

    trbc_pairs <- list()                    # Contains all pairs whose barcodes been assigned a true barcode.
    for (bc in true_barcodes[,1]) {trbc_pairs[[bc]] <- as.data.frame(matrix(nrow=0, ncol=ncol(raw_pairs)))}

    ambig_pairs <- as.data.frame(matrix(nrow=1, ncol=ncol(raw_pairs)+1))                                 #This table will store the processed barcode-spacer counts for those pairs that either the spacer or the barcode have not been properly matched to their true barcode and spacers

    ## Consolidating the barcodes and separating them into different data frames
    print(paste("          ", sample_name, "IDs", Sys.time(), sep = " "))
    for(pair in 1:nrow(raw_pairs)) {
        stringdist_result <- NULL; stringdist_result <- stringdist(raw_pairs[pair, 1], true_barcodes[,1], method = "lv") # This is done here so the computation intense operation of stringdist won't have to be repeated three times in identical fashion
        distance_barcode <- NULL; distance_barcode <- min(stringdist_result)
        consensus_barcode <- NULL; consensus_barcode <- as.character(true_barcodes[which(stringdist_result == min(stringdist_result))[1],1])

        if (distance_barcode <= max_dist_barcode) {
            trbc_pairs[[consensus_barcode]][nrow(trbc_pairs[[consensus_barcode]])+1,] <- raw_pairs[pair,]               # This strategy works slightly faster than rbind
        } else {ambig_pairs <- rbind(ambig_pairs, c( as.character(raw_pairs[pair,1]), as.character(raw_pairs[pair,2]), as.numeric(raw_pairs[pair,3]), "bc"))}
    }

    true_pairs <- as.data.frame(matrix(nrow=1, ncol=ncol(raw_pairs)))                                 #This table will store the processed barcode-spacer counts for those pairs that either the spacer or the barcode have not been properly matched to their true barcode and spacers

    ## Consolidating the spacers for each barcode
    print(paste("          ", sample_name, "spacers", Sys.time(), sep = " "))
    for (bc in names(trbc_pairs)) {
        true_spacers_this_bc <- as.data.frame(matrix(nrow=1, ncol=ncol(raw_pairs) - 1));
        #########################################################################################
        ## Consolidating the spacers for each barcode PART 1: This code is derived from Beltman et al, BMC bioinformatics, 2016. It eliminates sequencing errors based on certain parameters.
        #########################################################################################
        #
        # File: barcodecleanup.R
        #
        # Script to detect spurious barcodes in barcoding deepseq data.
        #
        # The script expects an input file containing a table with tab-separated read counts,
        # organized such that columns contain samples and rows contain barcodes.
        #
        # An example input file with artificial barcode data is given in the file 'testartificialdata.txt'
        # A formal description of the input data format is as follows (with lines uncommented):
        #
        #<samname_1><tab><samname_2><tab><samname_3>...<enter>
        #<barcode_1><tab><readnumber_sam1><tab><readnumber_sam2><tab><readnumber_sam3>...<enter>
        #<barcode_2><tab><readnumber_sam1><tab><readnumber_sam2><tab><readnumber_sam3>...<enter>
        #<barcode_3><tab><readnumber_sam1><tab><readnumber_sam2><tab><readnumber_sam3>...<enter>
        #.
        #.
        #.
        #
        # The input data should be ordered from most frequent to least frequent
        # in terms of total read count over all samples.
        #
        # Moreover, the barcodes containing an 'N' at at least one position should be removed.
        #
        # Required packages are 'stringdist' and 'VGAM'.
        #
        # Author: Joost Beltman
        # Date: August 28, 2015
        #
        #########################################################################################
        #########################################################################################



        #parameters for clean-up approach
        seqthresh = 5 	                	# threshold of number of nucleotide differences between correct mothers and daughters
        fracthresh = 0.2        			# threshold of maximum read count fraction allowed for correct mothers and daughters
        rdthresh  = 0.0   	                # threshold for minimal read count for individual data points (either mother or child #reads)
        totalrdthresh = 0                   # threshold for minimal total read count to consider barcodes at all. RK: to be re-assigned based on data
        slope 	= -2.0 						# slope for log-likelihood cutoff
        offset	= -2.0			 			# offset for log-likelihood cutoff


        # load unfiltered data and adjust format
        spacers.unfiltered <- trbc_pairs[[bc]]                                          # RK: The spacer file is formated differently and needs to be processed first
        spacers.unfiltered <- spacers.unfiltered[,-1]                                   # RK: removing the barcodes from the datafram
        spacers.unfiltered <- aggregate(. ~ V2, data = spacers.unfiltered, FUN = sum)   # RK: After removing the barcodes, spacers are not unique, this will aggregate the rows with the same spacer
        spacers.unfiltered <- spacers.unfiltered[order(-spacers.unfiltered[,2]),]       # Sorting based on the observed values.
        data.unfiltered <- spacers.unfiltered;
        rownames(data.unfiltered) <- spacers.unfiltered[,1]
        data.unfiltered <- subset(data.unfiltered, select = V3)

        #calculate total read counts (sums of rows)
        data.rs = cbind(data.unfiltered, as.data.frame(rowSums(data.unfiltered)))
        names(data.rs)[ncol(data.rs)] = "rs"

        #restrict analysis to data with at least 'totalrdthresh' reads
        data.manyrds = data.rs[data.rs$rs >= totalrdthresh,]

        #generate empty data frames for mothers and children
        mothers = as.data.frame(matrix(nrow=0, ncol=ncol(data.unfiltered)))
        names(mothers) = names(data.unfiltered)

        children = as.data.frame(matrix(nrow=0, ncol=ncol(data.unfiltered)))
        names(children) = names(data.unfiltered)

        #data frames to keep track of errors (set at default of no error)
        is.error = as.data.frame(matrix(nrow=nrow(data.manyrds), ncol=2))
        names(is.error) = c("error", "mother")
        rownames(is.error) = rownames(data.manyrds)
        is.error$error = rep(0, nrow(is.error))
        is.error$mother = rep("none", nrow(is.error))


        repeat {

                #try next mother
                newmoth = data.manyrds[1,]
                newmoth.bc = as.character(rownames(newmoth))

                #select subset of potential children based on total frequency
                data.freq = data.manyrds[data.manyrds$rs < fracthresh*newmoth$rs,]

                #select subset of potential children based on sequence similarity
                potkids = data.freq[stringdist(newmoth.bc, rownames(data.freq), method="hamming") <= seqthresh,]                    #RK: Changed the method from "lv" to "hamming" because we want to give greater weight to insertions and deletions.
                #get just read count data of current mother and potential children (without rowsums)
                rdcnt.moth = newmoth[!grepl("rs", names(newmoth),)]
                rdcnt.child = potkids[!grepl("rs", names(potkids),)]

                #calculate 'success probability' and 'shape parameter' of beta-binomial distribution for all potential children
                psuc = potkids$rs/newmoth$rs
                shape = psuc/10.0

                #data frame to store new kids and doubts
                newkids = as.data.frame(matrix(nrow=0, ncol=ncol(potkids)))
                names(newkids) = names(potkids)

                #loop over all potential children
                if (nrow(potkids) >= 1) {

                        for (i in 1:nrow(potkids)) {

                                #transpose readcounts
                                rdcnt = as.data.frame(t(rbind(rdcnt.moth, rdcnt.child[i,])))
                                names(rdcnt) = c("moth", "child")

                                #find read counts that are above threshold in either mother or child
                                rdcnt.minthresh = subset(rdcnt, moth>rdthresh | child>rdthresh)

                                #if at least one observation we calculate the log-likelihood score
                                nminthresh = nrow(rdcnt.minthresh)
                                if (nminthresh > 0) {

                                        #calculate log-likelihoods for each observation
                                        logllh = vector(length = nrow(rdcnt.minthresh))
                                        for (j in 1:nminthresh) {
                                            if ( packageVersion("VGAM") >= "1.0.7" ) {logllh[j] = loglink(dbetabinom(rdcnt.minthresh[j,2], rdcnt.minthresh[j,1], psuc[i], shape[i]))}
                                            if ( packageVersion("VGAM") < "1.0.7" ) {logllh[j] = loge(dbetabinom(rdcnt.minthresh[j,2], rdcnt.minthresh[j,1], psuc[i], shape[i]))}
                                        }

                                        #calculate mean log-likelihood
                                        meanlogllh = mean(logllh)

                                        #total number of reads for child
                                        nrdcnt.child = sum(rdcnt$child)

                                        #check if log-likelihood above threshold, in that case set flag that this is considered a child
                                        if (meanlogllh > (slope * log10(nrdcnt.child+1) + offset)) {                                        #RK: added the +1 to nrdcnt.child because log10(1) =0, so the slope element becomes completely eliminated if a read has a frequency of  1.

                                                newkids = rbind(newkids, potkids[i,])
                                                is.error[rownames(is.error) == rownames(potkids[i,]),]$error = 1
                                                is.error[rownames(is.error) == rownames(potkids[i,]),]$mother = newmoth.bc

                                        }

                                } #if at least one observation

                        } #loop over potential children
                } #if at least one potential child


                #bind new mother to mothers
                mothers = rbind(mothers, newmoth)

                #bind mother and kids
                newkids = rbind(newmoth, newkids)

                #if there are new children, bind them to the existing ones
                if (nrow(newkids) > 1) {
                        children = rbind(children, newkids[2:nrow(newkids),])
                }

                #remove children
                data.manyrds = data.manyrds[!(rownames(data.manyrds) %in% rownames(newkids)),]

                #continue until no more barcodes
                if (nrow(data.manyrds) == 0) {
                        break
                }

        }

        #chop 4 characters of input file name to construct output names
        #basisrf = substr(inputrf, 0, nchar(inputrf)-4) #chop off last 4 characters, i.e., file extension
        #basisrf = tail(strsplit(basisrf, split = "/")[[1]], n = 1)    #Remove the address in the file
        #basisrf = strsplit(basisrf, split = "_")[[1]][1]    #adjust the name
        basisrf = paste(sample_name, bc, sep = "-")

        #save cleaned data to files
        mothers.save = mothers[!grepl("rs", names(mothers),)]
        if (write.individual.tables) { write.table(mothers.save, paste("intermediate_files/", basisrf, "_truespacer.txt", sep=""), sep="\t", quote=F, row.names=T, col.names = FALSE) }

        if (nrow(children) > 1) {children.order = children[with(children, order(-rs)),]} else {children.order <- children}
        children.save = children.order[!grepl("rs", names(children.order),)]
        if (write.individual.tables) { write.table(children.save, paste("intermediate_files/", basisrf, "_errsp.txt", sep=""), sep="\t", quote=F, row.names=T, col.names = FALSE) }

        true_spacers <- as.data.frame(matrix(nrow = nrow(mothers.save), ncol = 2))
        true_spacers[,1] <- rownames(mothers.save);
        true_spacers[,2] <- mothers.save[,1]

        #########################################################################################
        ## Consolidating the spacers for each barcode PART 2: This part reads the true spacers that were derived in part one as basis for assigning a sequencing-error adjusted spacer to each barcode.
        #########################################################################################


        for(spacer in 1:nrow(trbc_pairs[[bc]])) {
            stringdist_result <- stringdist(trbc_pairs[[bc]][spacer,2], true_spacers[,1], method = "hamming") # This is done here so the computation intense operation of stringdist won't have to be repeated three times in identical fashion
            distance_spacer <- min(stringdist_result)
            consensus_spacer <- as.character(true_spacers[which(stringdist_result == min(stringdist_result))[1],1])

            if (distance_spacer <= max_dist_spacer) {
                true_spacers_this_bc <- rbind(true_spacers_this_bc, c(consensus_spacer, trbc_pairs[[bc]][spacer,3]) );
            } else {ambig_pairs <- rbind(ambig_pairs, c(as.character(trbc_pairs[[bc]][spacer,1]), as.character(trbc_pairs[[bc]][spacer,2]), as.numeric(trbc_pairs[[bc]][spacer,3]), "sp"))}
        }
        if (is.na(true_spacers_this_bc[1,1])) {true_spacers_this_bc <- true_spacers_this_bc[-1,] }          #Removing the row with NA
        true_spacers_this_bc[,2] <- as.numeric(true_spacers_this_bc[,2])
        true_spacers_this_bc <- aggregate(. ~ V1, data = true_spacers_this_bc, FUN = sum)
        true_spacers_this_bc <- cbind(rep(bc, nrow(true_spacers_this_bc)), true_spacers_this_bc)                # Adding the barcode sequence to all the true spacers
        colnames(true_spacers_this_bc) <- colnames(true_pairs);
        true_pairs <- rbind(true_pairs, true_spacers_this_bc)
    }
    if (is.na(true_pairs[1,1])) {true_pairs <- true_pairs[-1,]}
    true_pairs <- true_pairs[order(true_pairs[,1], -true_pairs[,3]),]                       # Reordering the values based on first barcode and second decreasing count
    true_pairs[,3] <- as.numeric(true_pairs[,3])
    write.table(true_pairs, paste(sample_name, "_truepairs.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)

    if (is.na(ambig_pairs[1,1])) {ambig_pairs <- ambig_pairs[-1,]}
    ambig_pairs[,3] <- as.numeric(ambig_pairs[,3])
    ambig_pairs <- ambig_pairs[order(-ambig_pairs[,3]),]                       # Reordering the values based on first barcode and second decreasing count
    write.table(ambig_pairs, paste(sample_name, "_ambigpairs.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)

    all_pairs <- NULL; all_pairs <- rbind(true_pairs, ambig_pairs[,1:ncol(true_pairs)])
    all_pairs <- all_pairs[order(all_pairs[,1], -all_pairs[,3]),]
    write.table(all_pairs, paste(sample_name, "_allpairs.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)

    ### Making and writing the genotype file
    pairs <- NULL; pairs <- true_pairs
    if (("TACACCTGCA" %in% pairs$V1)) {                                             # This manual adjustment should be applied to ALL PB7 line samples (reason in the line below).
        pairs[pairs$V1 == "TACACCTGCA",1] <- "GGCCCCTACA"                           # It turns out that because PB7 line barcode GGCCCCTACA matches the end of the reverse amplification primer, the primer can bind the GGCCCC sequence and make it appear the barcode is TACACCTGCA (end of the actual barcode truncated to the post-barcode sequence of the construct).
        pairs <- aggregate(V3~V1+V2, pairs, sum)
    }
    pairs <- pairs[order(pairs$V1, -pairs$V3),]         #sorting pairs on barcode followed by descending frequency because some filtering operations above can disrupt order
    pairs$V4 <- NA;                                     # For each ID, this column indicates its percentage among all IDs in that sample
    pairs$V5 <- NA;                                     # For each spacer, this column indicates its percentage in that hgRNA in that sample
    PCOFF_bc <- NULL; PCOFF_bc <- 1/min(60,length(levels(as.factor(as.vector(pairs[,1])))))^1.70 * 100                       # Setting up the maximum number of barcodes/IDs that each sample can have. No sample can have more than 60.
    CCOFF_pair <- 10                                    # Read count cutoff is the minimum number of reads that a pair needs to have to be considered.
    pairs <- pairs[-pairs[,3] < CCOFF_pair,]             # Removing pairs with less than CCOFF_pair reads.

    for(bc in levels(as.factor(as.vector(pairs[,1])))) {
        if ( sum(subset(pairs, V1 == bc)[,3]) / sum(pairs[,3]) * 100 >= PCOFF_bc) {
            cutoff <- max(pairs[,3][pairs[,1] == bc])
            SPinID <- round(max(pairs[,3][pairs[,1] == bc])/sum(pairs[,3][pairs[,1] == bc])*100, digits = 1)
            IDinALL <- round(sum(pairs[,3][pairs[,1] == bc])/sum(pairs[,3])*100, digits = 2)
            if (length(which(pairs[,3] < cutoff & pairs[,1] == bc, arr.ind = TRUE)) > 0 ) {pairs <- pairs[-which(pairs[,3] < cutoff & pairs[,1] == bc, arr.ind = TRUE),]}
            pairs[,4][pairs[,1] == bc] <- IDinALL
            pairs[,5][pairs[,1] == bc] <- SPinID
        } else {pairs <- subset(pairs, V1 != bc)}
    }

    pairs$V6 <- NA
    for(bc in levels(as.factor(as.vector(pairs[,1])))) {
        pairs[,6][pairs[,1] == bc] <- which(levels(as.factor(as.vector(pairs[,1]))) == bc)
    }
    pairs <- pairs[!duplicated(pairs$V6),]

    pairs <- pairs[,c(6, 1:5)]
    colnames(pairs) <- c('#', 'ID', 'SP', 'count', 'ID%inAll', 'SP%inID')

    write.table(pairs, paste(sample_name, "_genotypes.txt", sep = ''), quote = FALSE , row.names = FALSE, sep = "\t" )
}


