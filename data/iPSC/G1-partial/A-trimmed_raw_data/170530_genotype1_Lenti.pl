#!/usr/local/bin/perl
use strict;
use warnings;


my $sample = $ARGV[0];
chomp $sample;
$sample = (split(/[\/,_]/, $sample))[4];
print "\n$sample\tstarted at ".(localtime)."\n";


## Setting run-specific parameters
my $UsePairedEnds = 0;										# 1 if you want the pair to be reported only if Read1 and Read2 having a perfectly matching spacer. 0 otherwise. For 1, Read2 should be at least 165 bases in length.
my $operating_system = ('osx','linux')[1];					# 0 for osx and 1 for linux. To be set depending on what operating system this is run on because each system uses a different version of blat.
my $execute_folder;					#This is where the most I/O intensive operations will be done. /dev/shm is a memory based folder in linux that would make I/O much faster. However, this variable can be left equal '' for operation on systems without /dev/shm.
if ($operating_system eq 'linux') {
	$execute_folder = '/dev/shm';							# /dev/shm is a memory based folder in linux that would make I/O much faster.
} else {
	$execute_folder = '.';									# OSX systems do not have /dev/shm. This value (".") uses the current folder for I/O. If your computer does not have an SSD hard-drive, we recommend mounting a RAM-based disk manually and using the path of that folder here. See https://stackoverflow.com/questions/2033362/does-os-x-have-an-equivalent-to-dev-shm for instructions on how to do so.
}
my $execute_subfolder = $sample."1";				#This is the subfolder that will be created in $execute_folder to carry out operations to prevent interference from similar codes.
my $execute_path = $execute_folder. '/' . $execute_subfolder;

my $inputs_folder = '../0-raw_data_Lenti';

my @read1_trim = (0,180);					#Trim coordinates for read1. First value is the start of the region to be kept, last value is the length of the region to be kept.
my @read2_trim = (28, 10);					#Trim coordinates for read1. First value is the start of the region to be kept, last value is the length of the region to be kept.
my @read1_spacer = (87,50);					#Trim coordinates for extracting the gRNA sequence from read1. First value is the start of the region to be kept, last value is the length of the region to be kept.

my $prebarcode_sequence = "CGAGGTCGAGAATTC";																	# Expected sequence to be observed in Read2 before the barcode starts
$prebarcode_sequence = capital($prebarcode_sequence);

my $prespacer_sequence = "atggactatcatatgcttaccgtaacttgaaagtatttcgatttcttggctttatatatcttgtggaaaggacg";					# Expected sequence to be observed in Read1 before the spacer region starts
$prespacer_sequence = capital($prespacer_sequence);
my $prespacer_sequence_short = substr($prespacer_sequence, -12);

#my $postspacer_sequence = "GTTAACCTAA";																				# Expected sequence to be observed in Read1 after the spacer region ends, it is chosen to land on the downstream mutated sites in the hgRNA backbone.
#$postspacer_sequence = capital($postspacer_sequence);




##Creating execute subfolder in execute folder
unless(-e $execute_path) {mkdir($execute_path) or die "Cannot make $execute_path";}

#Moving the blat file to $execute_path
if ($execute_path) {
		my $filename = "$execute_path/blat_$operating_system\_$sample";
		my $make_executable_command = `chmod +x ./blat_$operating_system`;
 		unless (-e $filename) {my $one_time_command = `cp ./blat_$operating_system $execute_path/blat_$operating_system\_$sample`;}
}


my $PRINT_DEFAULT = 0;
my $PRINT = $PRINT_DEFAULT;
my $PAUSE = 0;


opendir(DIR1, "$inputs_folder") or die "Cannot open $inputs_folder/\n";	
my @R1seqs = sort(grep(/$sample.+R1_001.fastq$/, readdir(DIR1)));	print "\nR1seqs are @R1seqs\n" if $PRINT;
rewinddir DIR1;
my @R2seqs = sort(grep(/$sample.+R2_001.fastq$/, readdir(DIR1)));	print "\nR2seqs are @R2seqs\n" if $PRINT;


my %barcodes;								#This hash will store statistics of each barcode that is extracted from the reverse read.
my %barcodesANDspacers;						#This hash will store statistics of each barcode-spacer combination that is extracted from the reverse and forward reads respectively.


foreach (my $i = 0; $i <= $#R1seqs; $i++) {
	my $filenameR1 = (split("_", $R1seqs[$i]))[0];
	my $filenameR2 = (split("_", $R2seqs[$i]))[0];
	if ($filenameR1 ne $filenameR2) { die "FATAL Error E0010: First and second read files, $R1seqs[$i] and $R2seqs[$i], do not correspond to the same sample!\n"; }
	open(R1, "$inputs_folder/$R1seqs[$i]") or die "Cannot open $inputs_folder/$R1seqs[$i]\n";
	open(R2, "$inputs_folder/$R2seqs[$i]") or die "Cannot open $inputs_folder/$R2seqs[$i]\n";
#	open(OUTPUT1, ">$filenameR1\.fastq") || die "Cannot open output $filenameR1\.fastq\n";
#	open(OUTPUT2, ">$filenameR1\.unaligned") || die "Cannot open output $filenameR1\.unaligned\n";
	open(REPORT3, ">$filenameR1\_barcode_counts.txt") || die "Cannot open output $filenameR1\_barcode_counts.txt\n";
	open(REPORT4, ">$filenameR1\_barcode_gRNA_pairing_counts.txt") || die "Cannot open output $filenameR1\_barcode_gRNA_pairing_counts.txt\n";


	my $line = my $valid_amplicon = my $total_amplicon = 0;
	my $name = my $read1 = my $read2 = my $qual1 = my $qual2 = my $barcode = my $spacer = '';

	LOOP1: while (<R1>) {
			$line++;
			chomp; my $file1 = $_;
			LOOP2: while (<R2>) {
					chomp; my $file2 = $_;
		
					if ($line % 4 == 1) {			### First fastq Line
							unless ( (split(' ', $file1))[0] eq (split(' ', $file2))[0] ) {die "FATAL Error E0020: out of sync in reading the forward and reverse files\n";}
							$name = join(':', ((split(/[-, ]/, $file1))[1], (split(':', $file1))[-1] ) );
					}

					elsif ($line % 4 == 2) {		### Second fastq Line
							$read1 = $file1;
							$read2 = $file2;
					}
		
					elsif ($line % 4 == 3) {		### Third fastq line
					}
		
					elsif ($line % 4 == 0) {		### Fourth fastq line
							$total_amplicon++;
							$qual1 = $file1;
							$qual2 = $file2;

							print "\n$filenameR1\t$name\n" if $PRINT;
							
							## Determining the position of the barcode in R2 based on $prebarcode_sequence
							my $barcode_search = index($read2, $prebarcode_sequence);						# Searching for a perfect match to the pre-barcode sequence
							if ($barcode_search > -1) {														# There is a perfect match to the prebarcode_sequence
								$barcode = substr($read2, $barcode_search + length($prebarcode_sequence), $read2_trim[1]);					# Extracting the Barcode region from Read 2.
							} else {																		# There is no exact match to the prebarcode_sequence and a more time consuming but thorough search is needed.
								my @barcode_alignments = blat($read2, $prebarcode_sequence);					# Since I put only one query in, I expect the blat output to have only one line which is stored in the first element of the array. In the line below, I extract that first element into a second array.
								if (@barcode_alignments) {
										my @barcode_alignment = split("\t", $barcode_alignments[0]);	print "prebarcode alignment:\t".join("\t",@barcode_alignment)."\n" if $PRINT;
										if ($barcode_alignment[17] == 1 && $barcode_alignment[0] > 0.8*length($prebarcode_sequence)) { 												# Ensuring that the $prebarcode_sequence (query) has aligned to only one position in $read2 (target) 
												$barcode = substr($read2, $barcode_alignment[16], $read2_trim[1]);		# Extracting the Barcode region from Read 2.
										}
								}
							}
							print "barcode:\t $barcode\n" if $PRINT;
							

							## Determining the position of the spacer in R1 based on $prespacer_sequence
							my $spacer_search = index($read1, $prespacer_sequence_short);							# Searching for a perfect match to the pre-spacer sequence. The short version is being used to speed up the search.
							if ($spacer_search > -1) {														# There is a perfect match to the prespacer_sequence
								$spacer = substr($read1, $spacer_search + length($prespacer_sequence_short)-1, $read1_spacer[1]);			# Extracting the Spacer region from Read 1.
							} else {																		# There is no exact match to the prespacer_sequence and a more time consuming but thorough search is needed. The full prespacer_sequence will be used.							
								my @spacer_alignments1 = blat($read1, $prespacer_sequence);						# Since I put only one query in, I expect the blat output to have only one line which is stored in the first element of the array. In the line below, I extract that first element into a second array.
								#my @spacer_alignments2 = blat($read1, $postspacer_sequence);					# Since I put only one query in, I expect the blat output to have only one line which is stored in the first element of the array. In the line below, I extract that first element into a second array.
								if (@spacer_alignments1) {
										my @spacer_alignment1 = split("\t", $spacer_alignments1[0]);	print "prespacer alignment:\t".join("\t",@spacer_alignment1)."\n" if $PRINT;
										#my @spacer_alignment2 = split("\t", $spacer_alignments2[0]);	print "pstspacer alignment:\t".join("\t",@spacer_alignment2)."\n" if $PRINT;
										if ($spacer_alignment1[17] == 1 && $spacer_alignment1[0] > 0.8*length($prespacer_sequence) ) { 	# Ensuring that the $prespacer_sequence and $postspacer_sequence (queries) have aligned to only one position in $read1 (target) 
												$spacer = substr($read1, $spacer_alignment1[16]-1, $read1_spacer[1]);		# Extracting the spacer region from Read 1. 
										}
								}
							}
							print "spacer:\t $spacer\n" if $PRINT;

							#$read1 = substr($file1, $read1_trim[0], $read1_trim[1]);
							#$read2 = substr($file2, $read2_trim[0], $read2_trim[1]);
							#spacer = substr($file1, $read1_spacer[0], $read1_spacer[1]);
							#barcode = $read2;
							#$qual1 = substr($file1, $read1_trim[0], $read1_trim[1]);
							#$qual2 = substr($file2, $read2_trim[0], $read2_trim[1]);

							## Determining if the spacer sequence observed in Read1 exactly matches that observed in Read 2.
							my $R1R2match = 0;							# Is 0 if Read1 and Read2 don't match in spacer region. Otherwise, it is 1.
							if ($UsePairedEnds) {
								if ( $spacer && index($read2, reverse_complement( substr($spacer,4) )) != -1 ) {$R1R2match = 1;}			# The first five bases are being ommitted from spacer for two reasons: First, they are outside of the cut region. Second, with 165bp Read2s, the first two bases of the spacer  won't be covered in R2.								
							}

				
							if ( $barcode && $spacer && ($UsePairedEnds == $R1R2match) ) {									
									$name = join(':', $name, $barcode);
									#print OUTPUT1 ">$name\n$read1\n+\n$qual1\n";
									$barcodes{$filenameR1}{"$barcode"}++;
									$barcodesANDspacers{$filenameR1}{"$barcode\t$spacer"}++;
									$barcodes{'ALL'}{"$barcode"}++;
									$barcodesANDspacers{'ALL'}{"$barcode\t$spacer"}++;
									$valid_amplicon++;
							} else {#print OUTPUT2 "$filenameR1\t$name\tbarcode:$barcode\tspacer:$spacer\n";
							}
			
							$name = $read1 = $read2 = $qual1 = $qual2 = $barcode = $spacer = '';
					}
					last LOOP2;
			}
			if ($line % 100000 == 0) {print "$line enteries analyzed for $filenameR1 at ".(localtime)."\n"}
			#if ($line % 1000 == 0) {print "$line enteries analyzed for $filenameR1 at ".(localtime)."\n"; last LOOP1;}
	}
	
	if ($line) {
		my $succes_rate = sprintf("%.2f",($valid_amplicon/$total_amplicon*100));
		print "$line total enteries analyzed for $filenameR1 - ".(localtime)."\n";
		print "$succes_rate percent of all reads were identified as hgRNAs for $filenameR1 - ".(localtime)."\n";
	}

	## Printing barcode and spacer reports for each input file
	foreach my $key (sort { $barcodes{$filenameR1}{$b} <=> $barcodes{$filenameR1}{$a} } keys %{$barcodes{$filenameR1}} ) {
			print REPORT3 "$key\t$barcodes{$filenameR1}{$key}\n";
	}

	foreach my $key (sort { $barcodesANDspacers{$filenameR1}{$b} <=> $barcodesANDspacers{$filenameR1}{$a} } keys %{$barcodesANDspacers{$filenameR1}} ) {
			print REPORT4 "$key\t$barcodesANDspacers{$filenameR1}{$key}\n";
	}

#	close OUTPUT1;
#	close OUTPUT2;
	close R1;
	close R2;
	close REPORT3;
	close REPORT4;
}

close DIR1;
my @delete = system("rm -r $execute_path");

END;















### BLAT SUBROUTINE
#Uses the blat compilation in the home folder to produce the alignment.
sub blat {
	open(TARGET, ">$execute_path/.target.tmp.$sample") or die "Cannot open /dev/shm/.target.tmp.$sample\n";
	my $target = shift(@_);
	print TARGET ">target\n$target\n";
	close TARGET;
	
	open(QUERY, ">$execute_path/.query.tmp.$sample") or die "Cannot open $execute_path/.query.tmp.$sample\n";
	my $i = 0;
	while(@_) {
			$i++;
			my $query = shift(@_);
			if ($query) {print QUERY ">query$i\n$query\n";}
	}
	close QUERY;
	
	my @output = `$execute_path/blat_$operating_system\_$sample -t=dna -q=dna -tileSize=6 -stepSize=3 -oneOff=1 -minMatch=1 -minScore=5 -minIdentity=60 -maxGap=3 -noHead -noTrimA -out=psl $execute_path/.target.tmp.$sample $execute_path/.query.tmp.$sample $execute_path/.alignment.psl.$sample`;
	#system('./blat', '.target.tmp', '.query.tmp', '.alignment.psl');
	
	open(ALIGNMENT, "$execute_path/.alignment.psl.$sample") or die "Cannot open $execute_path/.alignment.psl.$sample\n";
	my @result = <ALIGNMENT>;
	close ALIGNMENT;
	chomp @result;
	return(@result);	
}

### ALIGNMENT SUBROUTINE
sub align {
	# scoring scheme
		#from: http://www.ebi.ac.uk/Tools/psa/emboss_water/help/index-nucleotide.html#matrix
		my $MATCH    =  5; # 	+5 for letters that match
		my $MISMATCH = -4; # 	-4 for letters that mismatch
		my $GAP      = -20; # 	-20 for any gap

	# variable transfer
		my $seq1 = $_[0];
		my $seq2 = $_[1];
		
	# initialization
		my @matrix;
		$matrix[0][0]{score}   = 0;
		$matrix[0][0]{pointer} = "none";
		for(my $j = 1; $j <= length($seq1); $j++) {
			$matrix[0][$j]{score}   = 0;
			$matrix[0][$j]{pointer} = "none";
		}
		for (my $i = 1; $i <= length($seq2); $i++) {
			$matrix[$i][0]{score}   = 0;
			$matrix[$i][0]{pointer} = "none";
		}

	# fill
		my $max_i     = 0;
		my $max_j     = 0;
		my $max_score = 0;

		for(my $i = 1; $i <= length($seq2); $i++) {
			for(my $j = 1; $j <= length($seq1); $j++) {
				my ($diagonal_score, $left_score, $up_score);
		
			# calculate match score
				my $letter1 = substr($seq1, $j-1, 1);
				my $letter2 = substr($seq2, $i-1, 1);       
				if ($letter1 eq $letter2) {
					$diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
				}
				else {
					$diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
				}
		
			# calculate gap scores
				$up_score   = $matrix[$i-1][$j]{score} + $GAP;
				$left_score = $matrix[$i][$j-1]{score} + $GAP;
		
				if ($diagonal_score <= 0 and $up_score <= 0 and $left_score <= 0) {
					$matrix[$i][$j]{score}   = 0;
					$matrix[$i][$j]{pointer} = "none";
					next; # terminate this iteration of the loop
				}
		
			# choose best score
				if ($diagonal_score >= $up_score) {
					if ($diagonal_score >= $left_score) {
						$matrix[$i][$j]{score}   = $diagonal_score;
						$matrix[$i][$j]{pointer} = "diagonal";
					}
					else {
						$matrix[$i][$j]{score}   = $left_score;
						$matrix[$i][$j]{pointer} = "left";
					}
				} else {
					if ($up_score >= $left_score) {
						$matrix[$i][$j]{score}   = $up_score;
						$matrix[$i][$j]{pointer} = "up";
					}
					else {
						$matrix[$i][$j]{score}   = $left_score;
						$matrix[$i][$j]{pointer} = "left";
					}
				}
		
			# set maximum score
				if ($matrix[$i][$j]{score} > $max_score) {
					$max_i     = $i;
					$max_j     = $j;
					$max_score = $matrix[$i][$j]{score};
				}
			}
		}

	# trace-back

		my $align1 = "";
		my $align2 = "";

		my $j = $max_j;
		my $i = $max_i;
		
		my $alignment_score = 0;

		while (1) {
			last if $matrix[$i][$j]{pointer} eq "none";
	
			if ($matrix[$i][$j]{pointer} eq "diagonal") {
				$align1 .= substr($seq1, $j-1, 1);
				$align2 .= substr($seq2, $i-1, 1);
				$i--; $j--;
				$alignment_score = $alignment_score + $matrix[$i][$j]{score};
			}
			elsif ($matrix[$i][$j]{pointer} eq "left") {
				$align1 .= substr($seq1, $j-1, 1);
				$align2 .= "-";
				$j--;
				$alignment_score = $alignment_score + $GAP;
			}
			elsif ($matrix[$i][$j]{pointer} eq "up") {
				$align1 .= "-";
				$align2 .= substr($seq2, $i-1, 1);
				$i--;
				$alignment_score = $alignment_score + $GAP;
			}   
		}

		$align1 = reverse $align1;
		$align2 = reverse $align2;

		my @output = ($j, $max_j-1, $i, $max_i-1, $alignment_score, $align1, $align2);	#$max_j-1 and $max_i-1 produce the indexes of the positions where the alignment ends.
		return @output;
}

### REVERSE COMPLEMENT SUBROUTINE
sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
        return $revcomp;
}

### CAPITALIZE SUBROUTINE
sub capital {
        my $dna = shift;

	# convert all letters to capital
        $dna =~ tr/abcdghmnrstuvwxy/ABCDGHMNRSTUVWXY/;
        return $dna;
}

### SMALLIZE SUBROUTINE
sub small {
        my $dna = shift;

	# convert all letters to capital
        $dna =~ tr/ABCDGHMNRSTUVWXY/abcdghmnrstuvwxy/;
        return $dna;
}

### Find minimum SUBROUTINE
sub min {
    my ($min, @vars) = @_;
    for (@vars) {
        $min = $_ if $_ < $min;
    }
    return $min;
}

### Find maximum SUBROUTINE
sub max {
    my ($max, @vars) = @_;
    for (@vars) {
        $max = $_ if $_ > $max;
    }
    return $max;
}
