
# lines = readLines(paste0("../indelphi/indelphi_result_summary_barcode_3cycle_no_cutoff/indelphi_result_summary_barcode", which(table_s1_active$`Identifier (ID)` == "AGGTGTTCAA"), ".txt"))
# seqs = sapply(1:(length(lines)/2), function(i) {
#         lines[i*2-1]
# })
# library(Biostrings)
# spacers = sapply(seqs,function(x) {
#         m = matchPattern(DNAString(toupper("tgcttaccgtaacttgaaagtatttcgatttcttggctttatatatcttgtggaaaggacg")),
#                      DNAString(x),
#                      max.mismatch = 20,
#                      with.indels = T)
#         if (length(m) > 0) {
#                 assertthat::assert_that(length(m) == 1)
#                 substr(x, end(m)+1, end(m)+50)
#         } else {
#                 return(NA)
#         }
# })
# table(is.na(spacers))
# names(spacers) = NULL
# table(sapply(spacers, str_length))
# seqs_len50 = spacers[sapply(spacers, str_length) == 50]
# seqs_len50 = seqs_len50[!is.na(seqs_len50)]
#
# seqs_len40 = spacers[sapply(spacers, str_length) == 40]
# seqs_len40 = seqs_len40[!is.na(seqs_len40)]
# sum(paste0(seqs_len40, toupper("ggcaccgagt")) %in% id_spacer_emb$AGGTGTTCAA$unique_spacer$spacer)
#
# seqs_len40 %in%
#
# msa::msa(inputSeqs = c(toupper("cgatttcttggctttatatatcttgtggaaaggac"), seqs[!ind][1]), type = "dna")
#
# "GAAACACCGGTACCACTCCCGTCGGGAGTAGTGGGTTAGAGCTAGAAATA" %in% spacers
#
# id_spacer_emb$AGGTGTTCAA$unique_spacer
