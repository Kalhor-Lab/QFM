# Script to determine high-confidence identifier sequences, to help address spurios ones that result from sequencing errors.
# in this code and others "barcode" refers to the Identifier.

#load required packages
library(stringdist)
library(VGAM)

#########################################################################################
# PART 1:   This first part of the code takes various tables and sends them to the second part of the code.
#########################################################################################

files1 <- system('ls ../1-*/*_barcode_counts.txt', intern = TRUE)

max_expected_barcodes <- 60;
COF <- 1/(max_expected_barcodes)^2              # Cutoff frequency for barcodes. Any barcode whose total frequency is less than COF*100% of the more frequent barcodes will not be considered. The formulae is based on max_expected barcodes^2 to allow a distribution around the mean of expected minimum amounts. The higher the total number of barcodes, the more variation is their fraction expected to have, thus square of max_expected_barcodes.

for (file1 in files1) {
    message(file1)
    if(file.info(file1)[1] == 0) {print(paste(file1, " is empty ****", sep = "")); next}
        #########################################################################################
        # PART 2: This code is derived from Beltman et al, BMC bioinformatics, 2016. It eliminates sequencing errors based on certain parameters.
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
        inputrf = file1                                 	# name of input file
        seqthresh = 2 	                		      	# threshold of number of nucleotide differences between correct mothers and daughters
        fracthresh = 0.2        				# threshold of maximum read count fraction allowed for correct mothers and daughters
        rdthresh  = 0.0 	                		# threshold for minimal read count for individual data points (either mother or child #reads)
        totalrdthresh = 0                   #RK: to be re-assigned based on data  # threshold for minimal total read count to consider barcodes at all
        slope 	= -2.0 						# slope for log-likelihood cutoff
        offset	= -2.0			 			# offset for log-likelihood cutoff


        # read unfiltered data
        #data.unfiltered = read.table(inputrf, sep="\t", header=TRUE)
        data.unfiltered = read.table(inputrf, sep="\t", header=FALSE, row.names = 1)    # Reza: changed it so it won't require the header row
        if (nrow(data.unfiltered) > 1) {
          for(i in 2:nrow(data.unfiltered)) {
            if (sum(data.unfiltered[1:i-1,1])*COF > data.unfiltered[i,1]) {
              totalrdthresh = data.unfiltered[min(i, max_expected_barcodes),1];
              break; #  Only  barcodes would be considered whose observed frequency is more than COF*100% of all the more frequent barcodes. Though under no circumstance more than max_expected_barcodes will be considered.
            }
          }
        } else {totalrdthresh  = 0 }

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
                potkids = data.freq[stringdist(newmoth.bc, rownames(data.freq), method="lv") <= seqthresh,]

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
                                        if (meanlogllh > (slope * log10(nrdcnt.child) + offset)) {

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
        basisrf = substr(inputrf, 0, nchar(inputrf)-4)                  #chop off last 4 characters, i.e., file extension
        basisrf = tail(strsplit(basisrf, split = "/")[[1]], n = 1)      #Remove the address in the file
        basisrf = strsplit(basisrf, split = "_")[[1]][1]                #adjust the name

        #save cleaned data to files
        mothers.save = mothers[!grepl("rs", names(mothers),)]
        write.table(mothers.save, paste(basisrf, "_trueID.txt", sep=""), sep="\t", quote=F, row.names=T, col.names = FALSE)

        if (nrow(children) > 1) {children.order = children[with(children, order(-rs)),]} else {children.order <- children}
        children.save = children.order[!grepl("rs", names(children.order),)]
        #write.table(children.save, paste(basisrf, "_errID.txt", sep=""), sep="\t", quote=F, row.names=T, col.names = FALSE)
    }
