#! /usr/bin/perl

# EXAMPLE CALL:
# perl sas12_data_cwn.pl my_OLGs.fasta 1.28 1

# PROGRAM: Implementation of the Wei-Zhang method for estimating dN/dS in overlapping genes (OLGs)
# sense-antisense overlap sas12:
# 123123
#   321321

# Copyright (C) 2017 Chase W. Nelson

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# DATE CREATED: November 3, 2017
# AUTHOR: Chase W. Nelson
# CONTACT1: cnelson@amnh.org
# CONTACT2: cwnelson88@gmail.com
# AFFILIATION: Sackler Institute for Comparative Genomics, American Museum of Natural History, New York, NY 10024, USA
# CITATION1: https://github.com/chasewnelson/overlapgenie
# CITATION2: Wei X, Zhang J. 2015. A simple method for estimating the strength of natural selection on overlapping genes. Genome Biol. Evol. 7:381–390.


use strict;   
use warnings; 
use diagnostics;
use autodie;
use Parallel::ForkManager;

#my $procs_per_node = 4;

#########################################################################################
# GATHER INPUT FROM ARGUMENTS
my $fasta_file_name = $ARGV[0];
my $R = $ARGV[1]; # This is transition/transversion ratio =a/(2b)
my $procs_per_node = $ARGV[2];

unless($fasta_file_name =~ /.fa/) { die "\n\n# FASTA file (argument 1) must contain .fa or .fasta extension. TERMINATED\n\n"; }
unless($R > 0) { die "\n\n# TS/TV ratio R (argument 2) must be provided, >0. TERMINATED\n\n"; }
unless($procs_per_node > 0) { die "\n\n# PROCS PER NOTE (argument 3) must be provided, >0. TERMINATED\n\n"; }

# Generate OUTPUT file name from FASTA file prefix
my $results_file_name = "OLG_results.txt";
if($fasta_file_name =~/\.fa/) { 
	$results_file_name = $` . "_OLG_results.txt";
} elsif($fasta_file_name =~/\.txt/) { 
	$results_file_name = $` . "_OLG_results.txt";
} elsif($fasta_file_name =~/\.csv/) { 
	$results_file_name = $` . "_OLG_results.txt";
}

#########################################################################################
# READ IN ALL SEQUENCES FROM FASTA FILE
my $seq = '';
my @seqs_arr;
my $header = '';
my @headers_arr;
my @seq_names_arr;
my $seq_num = 0;
my $last_seq_length;

open(IN_FASTA, "$fasta_file_name") or die "Could not open file $fasta_file_name\n";

print "\nRecording coding sequence data for $fasta_file_name...\n";

while(<IN_FASTA>) {
	chomp;
	if(/>/) {
		if($seq_num == 0) {
			$header = $_;
			$seq_num ++;
		} else {
			$seq = uc($seq);
			$seq =~ tr/U/T/;
			push(@seqs_arr,$seq);
			push(@headers_arr,$header);
			
			if($header =~ /^>([\w\.\-]+)/) {
				push(@seq_names_arr,$1);
			}
			
			$header = $_;
			$seq_num ++;
			
			my $this_seq_length = length($seq);
			
			#print "\nseq $seq_num is of length $this_seq_length\n";
			#print "\nseq: $seq\n";
			
			if($last_seq_length && ($last_seq_length != $this_seq_length)) {
				die "\n\nDIE: The sequences must be aligned, i.e., must be the same length. TERMINATED.\n\n";
			} else {
				$last_seq_length = $this_seq_length;
#				print "\nseq: $seq\n";
				$seq = '';
			}
		}
	} else {
		$seq .= $_;
	}
}

close IN_FASTA;

#print "\nseq: $seq\n";
$seq = uc($seq);
$seq =~ tr/U/T/;
push(@seqs_arr,$seq);
push(@headers_arr,$header);

if($header =~ /^>([\w\.\-]+)/) {
	push(@seq_names_arr,$1);
}

## print "\n";

# We now have an array of all our FASTA sequences, @seqs_arr, and their headers, @headers_arr
# We also have the ts/tv ratio, $R

my $num_seqs = scalar(@seqs_arr);
my $num_pw_comps = ($num_seqs**2 - $num_seqs) / 2;

srand (time|$$);  # introduce a random seed each time 

## # OLD SAMPLE INPUT:
## my $dna33 = 'ATG---GACAGAATTATT---------AGCTCATTAAAACAC------CAAGCATTAATTAGCCATCATAAGGTATTACGTAATACCTACTTTTTATTAGGGCTAACTCTAGCTTTTTCATCTGCTACTGCAACGTTGAGTACAATATACCATTTACCTTCG---CCCGGTTTTATACTTATGATAGCGGGGTTTTATGGCTTAATGTTTCTTACATATCGTTTAGCAAATCATCCAGCGGGTATCTTAGCCACCTTTGCGTTTACTGGTTTTCTAGGTTATTGTTTAGGACCTATTCTGAATTCATTTTTGTATTCAGGCATGGGAGACGTAATTGCTATGGCTCTAGGAGGTACTGCTCTAATATTCTTCTGTTGTTCTGCTTATGTATTAACAACACGCCGTGATATGTCATTTCTAAGTGGTACAATGCTAGTGGGTTGCGTATTATTATTAATAGGCACGCTTGTTAATATGTTTTTGCAATTAACTGTACTACATCTAACTATTAGTGCACTATTTATCTTATTCTCTACAGGTGCCATTTTATGGGAAACTAGTAATATCATTCATGGTGGTGAAACTAATTACATTCGTGCTACAGTAAGTCTTTATGTATCTTTATATAATATTTTTGTCAGTTTGCTA---GGTTTA---GCTGCTAATAAA------GATTAA';
## my $dna36 = 'ATG---GATCGTATTGTCGTTTCCTCTACCTCTGCGCGC---------AGTTCACTGCTGAGTACGCACAAAGTCCTGCGCAATACCTATTTTCTGCTGGCGCTGACCCTGGCTTTCTCCGCCTTAACCGCCACCGTCAGTACGGTCCTGATGCTGCCGGCG---CCGGGCCTGCTGCTGATGCTGGTCGGCTTTTATGGCCTGATGTTTCTGACCTACCGCCTGGCCGATCGCCCGGCCGGTCTGCTGGCCGTATTCGCCCTGACCGGCTTTATGGGCTATACCCTTGGCCCGATCCTCAGCAGCTTTATCGCCAGCGGCGCAGGCGATCTGATCATGCTGGCGCTGGGCGGCACCGCACTGGTCTTCTTCTGCTGCTCCGCCTATGTGCTGACCACCCGTAAGGACATGTCTTTCCTCTCCGGTATGATGATGGCCGGTTTTGTGGTGCTGCTGCTGGCCGTGGTCGCGAATCTGTTCCTGCAGCTGCCGGCGCTCTCTTTAGCCATCAGCGCGCTGTTTATCCTGTTCTCTACCGGCGCTATCCTGTGGGAGACCAGCAATATCATTCACGGCGGAGAAACCAACTACATCCGCGCCACCGTTGGGCTGTATGTCTCGCTGTATAACCTGTTCATCAGCCTGCTCAGCCTGCTGGGATTC---GCCCGCAGCAAC---TGA';
## 
## # We find out positions that match (1) and don't (0)
## my $match_string = &match_two_seqs($dna33,$dna36);
## print "$match_string\n";
## # RESULT: 111111110010111010000000000101110000000000111111000011010010110000110110110010110111111110111010010101110110110111110110010000110110110010111110010000000010110011111110110010010110111010100110111111111010111111110110110110010110011100110110111010010111000110110010110110111010110111000010110110110110100000111010000000111000110110010110000011110110110110110110110010111111110110110110111110010110110110000110111110110110000111100111010100111100110010010010010101001010100111011110011110010010100110001011010110110110110111111010111111110110110110010111110110110111111111110110110111110110111110110110110110010110111110110010111110010110011110011110000000010000000000000000000000000101
## 
## # Find index and length of longest contiguous string of 0's
## my ($longest_mismatch_index,$longest_mismatch_length) = &longest_stretch_chars($match_string,0);
## print "\nThe longest stretch of 0's occurs at index $longest_mismatch_index with length $longest_mismatch_length\n\n";
## # RESULT: "The longest stretch of 0's occurs at index 656 with length 25". VERIFIED.
## 
## # Find index and length of longest contiguous string of TRIPLETS CONTAINING 0's
## my ($longest_triplet_mismatch_index,$longest_triplet_mismatch_length) = &longest_stretch_chars_in_triplets($match_string,0);
## print "\nThe longest stretch of 0's appearing in contiguous codons occurs at index $longest_triplet_mismatch_index with length $longest_triplet_mismatch_length\n";
## my $codon_position = ($longest_triplet_mismatch_index / 3) + 1;
## my $length_codons = ($longest_triplet_mismatch_length / 3);
## print "This mean codon position $codon_position with $length_codons codons\n\n";
## # RESULT: "The longest stretch of 0's appearing in contiguous codons occurs at index 627 with length 57...
## # This mean codon position 210 with 19 codons". VERIFIED.
## 
## my $num_mismatches_in_stretch = &num_zeros_in_substring($match_string,$longest_triplet_mismatch_index,$longest_triplet_mismatch_length);

#########################################################################################
# PARALLELIZE PAIRWISE SEQUENCE COMPARISONS
mkdir("temp_dNdS_comparisons");
chdir("temp_dNdS_comparisons");

my $pm_poly = Parallel::ForkManager->new($procs_per_node);

for(my $seq_index = 0; $seq_index < $num_seqs; $seq_index++) {
	
	for(my $next_seq_index = ($seq_index + 1); $next_seq_index < $num_seqs; $next_seq_index++) {
		$pm_poly->start and next;
		
		# Get the time
		my $time1 = time;
#		$comp_number++;
		
		print "\n$seq_names_arr[$seq_index]\t$seq_names_arr[$next_seq_index]\t";
		
		my $dna1 = $seqs_arr[$seq_index];
		my $dna2 = $seqs_arr[$next_seq_index];
#		print "$dna1\n$dna2\n";
		
		################################################################################
		# SPLIT THE SEQUENCES INTO SUBSTRINGS CONTAINING NO GAPS
		# Determine indices containing gaps
		my %gap_indices_in_pair = &determine_gap_indices_in_pair($dna1,$dna2);
		my @gap_indices_sorted_arr = sort {$a <=> $b} keys(%gap_indices_in_pair);
		my $num_gap_positions = 0;
		foreach(@gap_indices_sorted_arr) {
			$num_gap_positions += $gap_indices_in_pair{$_}; # 1 if gap, 0 if not
		}
		
		# If there are gaps, extract substrings without gaps in reference
		my %dna1_gapless_index_to_substring;
		if($num_gap_positions > 0) {
			%dna1_gapless_index_to_substring = &extract_gapless_substrings($dna1,%gap_indices_in_pair); # returns (@gapless_substring_start_indices,@gapless_substrings);
		} else {
			$dna1_gapless_index_to_substring{0} = $dna1;
		}
		
		my %dna2_gapless_index_to_substring; 
		if($num_gap_positions > 0) {
			%dna2_gapless_index_to_substring = &extract_gapless_substrings($dna2,%gap_indices_in_pair);
		} else {
			$dna2_gapless_index_to_substring{0} = $dna2;
		}
		
		print "\n";
		
#		foreach my $sorted_key (sort {$a <=> $b} keys(%dna1_gapless_index_to_substring)) {
#			print "$sorted_key: $dna1_gapless_index_to_substring{$sorted_key}\n";
#			print "$sorted_key: $dna2_gapless_index_to_substring{$sorted_key}\n\n";
#		}
		
		################################################################################
		# DETERMINE SUBSTRINGS WITH LONGEST CONTIGUOUS MISMATCHES
## 		# Figure out longest mismatch string
## 		my $this_match_string = &match_two_seqs($dna1,$dna2);
## 		
## 		# Find index and length of longest contiguous string of TRIPLETS CONTAINING 0's
## 		my ($this_longest_triplet_mismatch_index,$this_longest_triplet_mismatch_length) = &longest_stretch_chars_in_triplets($this_match_string,0);
## 		my $this_mismatch_codon_position = ($this_longest_triplet_mismatch_index / 3) + 1;
## 		my $this_length_mismatch_codons = ($this_longest_triplet_mismatch_length / 3);
## 		print "$this_mismatch_codon_position\t$this_length_mismatch_codons\t";
## 		
## 		my $num_zeros_in_stretch = &num_zeros_in_substring($this_match_string,$this_longest_triplet_mismatch_index,$this_longest_triplet_mismatch_length);
## 		
## 		print "$num_zeros_in_stretch\t";
		
		################################################################################
		# CALCULATE AND SUM dN/dS PARTS FOR ALL ANALYZABLE SUBSTRINGS, stored in
		# %dna1_gapless_index_to_substring AND %dna2_gapless_index_to_substring
		# M denotes number of mismatches, L denotes number of sites
		my $Mnn_running_total = 0;
		my $Msn_running_total = 0;
		my $Mns_running_total = 0;
		my $Mss_running_total = 0;
		my $Lnn_running_total = 0;
		my $Lsn_running_total = 0;
		my $Lns_running_total = 0;
		my $Lss_running_total = 0;
		
		print "\nCalculating dN/dS for pair...\n";
		
		foreach my $sorted_key (sort {$a <=> $b} keys(%dna1_gapless_index_to_substring)) {
#			print "$sorted_key: $dna1_gapless_index_to_substring{$sorted_key}\n";
#			print "$sorted_key: $dna2_gapless_index_to_substring{$sorted_key}\n\n";
			my $this_dna1_substring = $dna1_gapless_index_to_substring{$sorted_key};
			my $this_dna2_substring = $dna2_gapless_index_to_substring{$sorted_key};
			
			my @this_comp_dnds_parts = &calculate_dnds_parts($this_dna1_substring,$this_dna2_substring,$R);
			# returns ($n_s[0],$n_s[1],$n_s[2],$n_s[3],$N_S[0],$N_S[1],$N_S[2],$N_S[3]);
#			print "\nRESULTS:\n@this_substring_comp_dnds\n";
			# dnds1_by_ns dnds2_by_ns dnds1_by_ss dnds2_by_ss var(dnn) var(dsn) var(dns) var(dss) dnn dsn dns dss Mnn[12] Msn Mns Mss Lnn Lsn Lns Lss";
			$Mnn_running_total += $this_comp_dnds_parts[0];
			$Msn_running_total += $this_comp_dnds_parts[1];
			$Mns_running_total += $this_comp_dnds_parts[2];
			$Mss_running_total += $this_comp_dnds_parts[3];
			$Lnn_running_total += $this_comp_dnds_parts[4];
			$Lsn_running_total += $this_comp_dnds_parts[5];
			$Lns_running_total += $this_comp_dnds_parts[6];
			$Lss_running_total += $this_comp_dnds_parts[7];
		}
		
#		print "Mnn_total = $Mnn_running_total\n";
#		print "Msn_total = $Msn_running_total\n";
#		print "Mns_total = $Mns_running_total\n";
#		print "Mss_total = $Mss_running_total\n";
#		print "Lnn_total = $Lnn_running_total\n";
#		print "Lsn_total = $Lsn_running_total\n";
#		print "Lns_total = $Lns_running_total\n";
#		print "Lss_total = $Lss_running_total\n";
		
		################################################################################
		# CALCULATE dN/dS FROM TOTAL NUMBER OF N AND S DIFFERENCES AND SITES
		my @corrected_dnds = &correct_dnds_results($Mnn_running_total,$Msn_running_total,
			$Mns_running_total,$Mss_running_total,$Lnn_running_total,$Lsn_running_total,
			$Lns_running_total,$Lss_running_total);
		# Returns ($dNN_by_dSN,$dNN_by_dNS,$dNS_by_dSS,$dSN_by_dSS,$var_dn1n2,$var_ds1n2,$var_dn1s2,$var_ds1s2,$dn1n2,$ds1n2,$dn1s2,$ds1s2,$pn1n2,$ps1n2,$pn1s2,$ps1s2,$var_pn1n2,$var_ps1n2,$var_pn1s2,$var_ps1s2,$pn_1,$ps_1,$pn_2,$ps_2,$dn_1,$ds_1,$dn_2,$ds_2)
		
### 		my @this_comp_dnds;
### 		
### 		# Skip comparisons with long gaps || KEEP IN MIND THIS MAY ACTUALLY BE ABOUT NUM. DIFFS, not GAPS
### 		if($dna1 =~ /----------------------------------------/ || $dna2 =~ /----------------------------------------/) { # 40 nts, i.e., longer than 39nt/3=13 codons
### #			print "One sequence contains a stretch of indels >13 codons; skipping comparison.\n"; # $seq_names_arr[$seq_index] has vs. $seq_names_arr[$next_seq_index]\...\n";
### 			print "X";
### 			for(my $i=0; $i<20; $i++) {	
### 				push(@this_comp_dnds,'*');
### 			}
### 		} else {
### 			my $dna1_header = $headers_arr[$seq_index];
### 			my $dna2_header = $headers_arr[$next_seq_index];
### 			
### 			@this_comp_dnds = &calculate_dnds($dna1,$dna2,$R); # this is NOW of size 20
### 		}

		open(THIS_COMP_TEMP_FILE,">>s$seq_index\-s$next_seq_index");
		
		# Print corrected values and variances to temp file
		foreach(@corrected_dnds) {
			print THIS_COMP_TEMP_FILE "$_\t";
		}
		
		# Print raw counts to temp file
		print THIS_COMP_TEMP_FILE "$Mnn_running_total\t$Msn_running_total\t" .
			"$Mns_running_total\t$Mss_running_total\t$Lnn_running_total\t$Lsn_running_total\t" .
			"$Lns_running_total\t$Lss_running_total\t";
		
		my $time2 = time;
		my $time_elapsed = $time2 - $time1;
		
		print THIS_COMP_TEMP_FILE "$time_elapsed";
		
		close THIS_COMP_TEMP_FILE;
		
		$pm_poly->finish;
	}
}

$pm_poly->wait_all_children; # special name, methods within module

#########################################################################################
# SWEEP UP RESULTS FROM TEMP FILES
my %comp_dNdS_ha;

my @temp_comp_FILES_arr = glob "*";
#@temp_comp_FILES_arr = sort {$a <=> $b} @temp_comp_FILES_arr;

if($num_pw_comps != scalar(@temp_comp_FILES_arr)) {
	warn "### WARNING: we expected $num_pw_comps pairwise comparisons but got " . scalar(@temp_comp_FILES_arr) . "\n\n";
} else {
	print "### SUCCESS: $num_pw_comps pairwise sequence comparisons completed\n\n";
}

#print "sorted polymorphic_codons_FILE_arr: @temp_comp_FILES_arr\n";

foreach(@temp_comp_FILES_arr) { # file names, actually
	my $comp_file = $_;
	
	my $seq_1_index;
	my $seq_2_index;
	
	if($comp_file =~ /s(\d+)\-s(\d+)/) {
		$seq_1_index = $1;
		$seq_2_index = $2;
	} else {
		die "\nRegex problem with temp files\n";
	}
	
	open(CURR_COMP_FILE, $comp_file) or die "\n## Cannot open $comp_file. TERMINATED.\n\n";
	while(<CURR_COMP_FILE>) {
		chomp; # Not necessary but hey
		my @this_comp_dnds = split(/\t/,$_,-1);
		$comp_dNdS_ha{$seq_1_index}->{$seq_2_index} = \@this_comp_dnds;
		#push(@polymorphic_codons_arr,$line_arr[0]);	
		#$num_seqs_defined[$comp_file] = $line_arr[1];
	}
	close CURR_COMP_FILE;
	unlink $comp_file;
}

chdir("..");
rmdir("temp_dNdS_comparisons");
# Finish parallelism for polymorphism detection

#########################################################################################
# PRINT RESULTS TO $results_file_name

open(OUTPUT,">>$results_file_name");
# Header
print OUTPUT "seq1\tseq2\tdNN\/dSN\tdNN\/dNS\tdNS\/dSS\tdSN\/dSS\tvar(dNN)\tvar(dSN)\tvar(dNS)\tvar(dSS)\t" . 
	"dNN\tdSN\tdNS\tdSS\tpNN\tpSN\tpNS\tpSS\tvar(pNN)\tvar(pSN)\tvar(pNS)\tvar(pSS)\t".
	"pN_ref\tpS_ref\tpN_alt\tpS_alt\tdN_ref\tdS_ref\tdN_alt\tdS_alt\t".
	"Mnn\tMsn\tMns\tMss\tLnn\tLsn\tLns\tLss\t".
	"time_used\n";

#for (my $i=1; $i<=$size_of_set; $i++) { # if size of set is 1, this happens only once. "SET" must mean a pairwise comparison
## for (my $i=0; $i<$size_of_set; $i++) { # if size of set is 1, this happens only once. "SET" must mean a pairwise comparison
foreach my $i (sort {$a <=> $b} keys %comp_dNdS_ha) {
	
	foreach my $j (sort {$a <=> $b} keys %{$comp_dNdS_ha{$i}}) {
		
		my @this_comp_arr = @{$comp_dNdS_ha{$i}->{$j}};
		
		my $seq_name1 = $seq_names_arr[$i];
		my $seq_name2 = $seq_names_arr[$j];
		
		print OUTPUT "$seq_name1\t$seq_name2\t";
		
		for(my $k=0; $k<scalar(@this_comp_arr); $k++) {
			
			print OUTPUT "$this_comp_arr[$k]\t"; # they didn't know how to do multidimensional arrays?
##			print "$dnds_2genes[$j+20*$i]\t"; # they didn't know how to do multidimensional arrays?
		}
		
		print OUTPUT "\n";
		
	}
	
}

close OUTPUT;

#########################################################################################
# END PROGRAM
print "\n#### PROGRAM COMPLETE ####\n\n";
exit;


#########################################################################################
#########################################################################################
####################################                 ####################################
####################################   SUBROUTINES   ####################################
####################################                 ####################################
#########################################################################################
#########################################################################################


#########################################################################################
sub match_two_seqs {
	my ($seq1,$seq2)=@_;
	
	my $match_string = '';
	
	if(length($seq1) != length($seq2)) {
		$match_string = 'ERROR';
	} else {
		while(length($seq1)>0 && length($seq2)>0) {
			my $seq1_char = substr($seq1,0,1,'');
			my $seq2_char = substr($seq2,0,1,'');
			
			if($seq1_char eq $seq2_char) {
				$match_string .= '1';
			} else {
				$match_string .= '0';
			}
		}
	}
	
	return $match_string;
}


#########################################################################################
sub longest_stretch_chars {
	my ($string,$char) = @_;
	
	my $longest_match_index = -1;
	my $longest_match_length = -1;
	
	my $curr_match_index = -1;
	my $curr_match_length = 0;
	my $within_match_flag = 0;
	
	if(length($string) > 0) {
	
		for(my $i=0; $i<(length($string)); $i++) {
		
			my $this_char = substr($string,$i,1);
			
			if($this_char eq $char) { # IT'S A MATCH
			
				if($within_match_flag == 0) { # not currently within a match, so re-set position
					
					$within_match_flag = 1;
					$curr_match_index = $i;
					$curr_match_length = 1;
					
					if($curr_match_length > $longest_match_length) {
						$longest_match_length = $curr_match_length;
						$longest_match_index = $curr_match_index;
					}
					
				} else { # currently within a match, and we've matched again
					# Already have a $curr_match_index
					$curr_match_length++;
					
					if($curr_match_length > $longest_match_length) {
						$longest_match_length = $curr_match_length;
						$longest_match_index = $curr_match_index;
					}
					
				}
			} else { # NOT A MATCH
			
				if($within_match_flag == 1) { # match is all done; already stored info
					$within_match_flag = 0;
				} # else not currently within a match, so nothing to be done
			}
		}
	}
	
	return ($longest_match_index,$longest_match_length);
}


#########################################################################################
sub longest_stretch_chars_in_triplets {
	my ($string,$char) = @_;
	
	my $longest_match_index = -1;
	my $longest_match_length = -1;
	
	my $curr_match_index = -1;
	my $curr_match_length = 0;
	my $within_match_flag = 0;
	
	if(length($string) > 0) {
	
		for(my $i=0; $i<(length($string)); $i+=3) {
		
			my $this_triplet = substr($string,$i,3);
			
			if($this_triplet =~ /$char/) { # IT CONTAINS A MATCH
			
				if($within_match_flag == 0) { # not currently within a match, so re-set position
					
					$within_match_flag = 1;
					$curr_match_index = $i;
					$curr_match_length = 3;
					
					if($curr_match_length > $longest_match_length) {
						$longest_match_length = $curr_match_length;
						$longest_match_index = $curr_match_index;
					}
					
				} else { # currently within a match, and we've matched again
					# Already have a $curr_match_index
					$curr_match_length+=3;
					
					if($curr_match_length > $longest_match_length) {
						$longest_match_length = $curr_match_length;
						$longest_match_index = $curr_match_index;
					}
					
				}
			} else { # NOT A MATCH
			
				if($within_match_flag == 1) { # match is all done; already stored info
					$within_match_flag = 0;
				} # else not currently within a match, so nothing to be done
			}
		}
	}
	
	return ($longest_match_index,$longest_match_length);
}


#########################################################################################
sub num_zeros_in_substring {
#sub num_chars_in_substring {
#	my ($string,$index,$length,$char) = @_;
	my ($string,$index,$length) = @_;
	
	#&num_chars_in_substring($match_string,$longest_triplet_mismatch_index,$longest_triplet_mismatch_length,0);
	
	my $this_substring = substr($string,$index,$length);
	
#	print "\nthis_substring is $this_substring\n";
	
	
#	my $count = eval "\$this_substring =~ tr/$char//"; # WORKS BUT UNNECESSARY NOW
	my $count = $this_substring =~ tr/0//;
	
#	print "count of $char is $count\n";
#	print "count of 0s is $count\n";

	return $count;
}


#########################################################################################
sub determine_gap_indices_in_pair {
	my ($dna1,$dna2) = @_;
	
	my %indices_with_gaps;
	
	for (my $i=0; $i<length($dna1); $i++) {
		my $dna1_char = substr($dna1,$i,1);
		my $dna2_char = substr($dna2,$i,1);
		
		if ($dna1_char eq '-' || $dna1_char eq 'N' || $dna2_char eq '-' || $dna2_char eq 'N') {
			$indices_with_gaps{$i} = 1;
		} else {
			$indices_with_gaps{$i} = 0;
		}
		
	}
	
	return %indices_with_gaps;
	
}


#########################################################################################
sub extract_gapless_substrings { # DO NOT CALL unless there is AT LEAST ONE GAP POSITION
	my ($dna,%gap_indices) = @_;
	
#	my @gap_indices_keys = keys (%gap_indices);
#	print "\ngap_indices_keys are @gap_indices_keys\n";
	
	# Recall that for sas12, what gets compared is:
#	[1-sense] (ATG)(CAA)(GAG)(TTG)
#	[2-sense] (ATG)(CAG)(GAG)(ATG)
#
#	[1-rvcom] C(AAC)(TCT)(TGC)AT
#	[2-rvcom] C(ATC)(TCC)(TGC)AT
	
#	my @gapless_substring_start_indices;
#	my @gapless_substrings;
	
	my %gapless_index_to_substring;
	
	my $curr_starting_index = 0;
	my $curr_substring = '';
	
	my $within_gapless_stretch_flag = 0;
	
	for (my $i = 0; $i < length($dna); $i+=3) {
		
		if($gap_indices{$i} == 1 || $gap_indices{$i+1} == 1 || $gap_indices{$i+2} == 1) { # there's a gap in this codon
			
			if($within_gapless_stretch_flag == 1) { # we're in a gapless stretch, now ended
				# We've already added previous codon, so just reset parameters
#				print "\nat index $curr_starting_index\, adding string $curr_substring\n";
#				push(@gapless_substrings,$curr_substring);
#				push(@gapless_substring_start_indices,$curr_starting_index);
				$gapless_index_to_substring{$curr_starting_index} = $curr_substring;
				$curr_substring = '';
				$within_gapless_stretch_flag = 0;
			} # else we're not in a gapless stretch yet
			
		} else { # NO GAP HERE

			if($within_gapless_stretch_flag == 1) { # we're already in a gapless stretch, simply add codon
				my $curr_codon = substr($dna,$i,3);
				$curr_substring .= $curr_codon;
				
				# Store substring data if we're at last codon
				if($i == (length($dna)-3)) {
#					print "\nat index $curr_starting_index\, adding string $curr_substring\n";
#					push(@gapless_substrings,$curr_substring);
#					push(@gapless_substring_start_indices,$curr_starting_index);
					$gapless_index_to_substring{$curr_starting_index} = $curr_substring;
				}
				
			} else { # beginning of a gapless stretch, add codon, re-set index, mark flag
				my $curr_codon = substr($dna,$i,3);
				$curr_substring = $curr_codon;
				$curr_starting_index = $i;
				$within_gapless_stretch_flag = 1;
			}
			
		}
	}
	
#	return(@gapless_substring_start_indices, @gapless_substrings);
	return(%gapless_index_to_substring);
	
}


#########################################################################################
sub correct_dnds_results {
	my ($Mnn,$Msn,$Mns,$Mss,$Lnn,$Lsn,$Lns,$Lss) = @_;
	
	# Determine raw values
	my $pn1n2 = '*';
	if ($Lnn > 0) {
		$pn1n2=$Mnn/$Lnn;
	}
	
	my $ps1n2 = '*';
	if ($Lsn > 0) {
		$ps1n2=$Msn/$Lsn;
	}
	
	my $pn1s2 = '*';
	if ($Lns > 0) {
		$pn1s2 = $Mns/$Lns;
	}
	
	my $ps1s2 = '*';
	if ($Lss > 0) {
		$ps1s2=$Mss/$Lss;
	}
	
	# NEW? alt dN/dS can be estimated by dNN/dNS dSN/dSS
	# e.g., pNN is Mnn/Lnn and pSN is Msn/Lsn
	# Thus, we can pool this measure into (Mnn+Msn)/(Lnn+Lsn)
	
	# Jukes-Cantor correct the values
	my $dn1n2;
	if ($pn1n2 >= 3/4){
		$dn1n2='*';
	} else {
		$dn1n2=-3/4*log(1-4/3*$pn1n2);
	}
	
	my $ds1n2;
	if ($ps1n2 >= 3/4){
		$ds1n2='*';
	} else {
		$ds1n2=-3/4*log(1-4/3*$ps1n2);
	}
	
	my $dn1s2;
	if ($pn1s2 >= 3/4){
		$dn1s2='*';
	} else {
		$dn1s2=-3/4*log(1-4/3*$pn1s2);
	}
	
	my $ds1s2;
	if ($ps1s2 >= 3/4){
		$ds1s2='*';
	} else {
		$ds1s2=-3/4*log(1-4/3*$ps1s2);
	}
	
	# Determine modified Nei-Gojobori values (for non-OLG use)
	my $Mn1 = $Mnn+$Mns;
	my $Mn2 = $Mnn+$Msn;
	my $Ms1 = $Msn+$Mss;
	my $Ms2 = $Mns+$Mss;
	
	my $Ln1 = $Lnn+$Lns;
	my $Ln2 = $Lnn+$Lsn;
	my $Ls1 = $Lsn+$Lss;
	my $Ls2 = $Lns+$Lss;
	
	my $pn_1 = '*'; 
	if ($Ln1 > 0) { 
		$pn_1 = $Mn1/$Ln1;
	}
	
	my $ps_1 = '*'; 
	if ($Ls1 > 0) { 
		$ps_1 = $Ms1/$Ls1;
	}
	
	my $pn_2 = '*'; 
	if ($Ln2 > 0) { 
		$pn_2 = $Mn2/$Ln2;
	}
	
	my $ps_2 = '*'; 
	if ($Ls2 > 0) { 
		$ps_2 = $Ms2/$Ls2;
	}
	
	# Jukes-Cantor correct modified Nei-Gojobori values (for non-OLG use)
	my $dn_1;
	if ($pn_1 >= 3/4){
		$dn_1='*';
	} else {
		$dn_1=-3/4*log(1-4/3*$pn_1);
	}
	
	my $dn_2;
	if ($pn_2 >= 3/4){
		$dn_2='*';
	} else {
		$dn_2=-3/4*log(1-4/3*$pn_2);
	}
	
	my $ds_1;
	if ($ps_1 >= 3/4){
		$ds_1='*';
	} else {
		$ds_1=-3/4*log(1-4/3*$ps_1);
	}
	
	my $ds_2;
	if ($ps_2 >= 3/4){
		$ds_2='*';
	} else {
		$ds_2=-3/4*log(1-4/3*$ps_2);
	}
	
	# Reference gene modified Nei-Gojobori dNdS
	my $origin_dnds1 = '*';
	if($ds_1 > 0) {
#		if($ds_1 eq '*') {
#			print "\n\n****AND YET WE GOT INSIDE!********\n\n";
#		}
		
		$origin_dnds1=$dn_1/$ds_1;
	}
	
	# Alternate gene modified Nei-Gojobori dNdS
	my $origin_dnds2 = '*';
	if($ds_2 > 0) {
		$origin_dnds2=$dn_2/$ds_2;
	}
	
	# Ref/Alt dN ratio
	my $dn_ratio = '*';
	if($dn_2 > 0) {
		$dn_ratio = $dn_1/$dn_2;
	}
	
	# Determine p variances for Wei-Zhang OLG dNdS
	my $var_pn1n2 = '*';
	if($Lnn > 0) {
		$var_pn1n2=$pn1n2*(1-$pn1n2)/$Lnn;
	}
	
	my $var_ps1n2 = '*';
	if($Lsn > 0) {
		$var_ps1n2=$ps1n2*(1-$ps1n2)/$Lsn;
	}
	
	my $var_pn1s2 = '*';
	if($Lns > 0) {
		$var_pn1s2=$pn1s2*(1-$pn1s2)/$Lns;
	}
	
	my $var_ps1s2 = '*';
	if($Lss > 0) {
		$var_ps1s2=$ps1s2*(1-$ps1s2)/$Lss;
	}
	
	# Determine d variances for Wei-Zhang OLG dNdS
	my $var_dn1n2 = '*';
	if($pn1n2 != 3/4) {
		$var_dn1n2 = $var_pn1n2/(1-4/3*$pn1n2)**2;
	}
	
	my $var_ds1n2 = '*';
	if($ps1n2 != 3/4) {
		$var_ds1n2 = $var_ps1n2/(1-4/3*$ps1n2)**2;
	}
	
	my $var_dn1s2 = '*';
	if($pn1s2 != 3/4) {
		$var_dn1s2 = $var_pn1s2/(1-4/3*$pn1s2)**2;
	}
	
	my $var_ds1s2 = '*';
	if($ps1s2 != 3/4) {
		$var_ds1s2 = $var_ps1s2/(1-4/3*$ps1s2)**2;
	}
	
	my $dnds_ref_1=0;
	my $dnds_ref_2=0;
	my $dnds_alt_1=0;
	my $dnds_alt_2=0;
	
	# dNN/dSN, i.e., dN/dS for reference gene
	if ($ds1n2!=0){
		$dnds_ref_1=$dn1n2/$ds1n2;
	}
	
	# dNN/dNS, i.e., dN/dS for alternative gene
	if ($dn1s2>0)	{
		$dnds_alt_1=$dn1n2/$dn1s2;
	}
	
	# dNS/dSS and dSN/dSS (alternative dN/dS measures for ref and alt, respectively)
	if ($ds1s2 > 0){
		$dnds_ref_2=$dn1s2/$ds1s2;
		$dnds_alt_2=$ds1n2/$ds1s2;
	}
	
	# dNN/dSN, i.e., dN/dS for reference gene
	my $dNN_by_dSN = '*';
	if($ds1n2 > 0) {
		$dNN_by_dSN = $dn1n2/$ds1n2;
	}
	
	# dNN/dNS, i.e., dN/dS for alternative gene
	my $dNN_by_dNS = '*';
	if($dn1s2 > 0) {
		$dNN_by_dNS = $dn1n2/$dn1s2;
	}
	
	# dNS/dSS for ref
	my $dNS_by_dSS = '*';
	if($ds1s2 > 0) {
		$dNS_by_dSS = $dn1s2/$ds1s2;
		print "\ndn1s2 is $dn1s2 and ds1s2 is $ds1s2\n\n";
	}
	
	# dSN/dSS for alt
	my $dSN_by_dSS = '*';
	if($ds1s2 > 0) {
		$dSN_by_dSS = $ds1n2/$ds1s2;
	}
	
	my @result = ($dNN_by_dSN,$dNN_by_dNS,$dNS_by_dSS,$dSN_by_dSS,$var_dn1n2,$var_ds1n2,$var_dn1s2,$var_ds1s2,$dn1n2,$ds1n2,$dn1s2,$ds1s2,
		$pn1n2,$ps1n2,$pn1s2,$ps1s2,$var_pn1n2,$var_ps1n2,$var_pn1s2,$var_ps1s2,
		$pn_1,$ps_1,$pn_2,$ps_2,$dn_1,$ds_1,$dn_2,$ds_2);

	return @result;

}


#########################################################################################
sub calculate_dnds_parts { # modified from Wei & Zhang (2015)
	my ($dna1,$dna2,$R)=@_;
#	print "Determing number of sites...\n";
	my @N_S=N_S($dna1,$dna2,$R); # NN [0], SN [1], NS [2], and SS [3] SITES (UPPERcase)
	
#	print "Determing number of differences...\n";
	my @n_s=n_s($dna1,$dna2); # NN [0], SN [1], NS [2], and SS [3] DIFFS (lowercase)
	
	my @result = ($n_s[0],$n_s[1],$n_s[2],$n_s[3],$N_S[0],$N_S[1],$N_S[2],$N_S[3]);

	return @result;

}


#########################################################################################
sub N_S { # modified from Wei & Zhang (2015)
	my ($dna1,$dna2,$R)=@_;
	my @result;
	my $N1N2=0;
	my $N1S2=0;
	my $S1N2=0;
	my $S1S2=0;
	my $N1=0;
	my $S1=0;
	my $N2=0;
	my $S2=0;
	my @N_S_1=N_S_dna($dna1,$R); 
	my @N_S_2=N_S_dna($dna2,$R);
	my $length=length($dna1);
##	print "length of dna is $length\n";
	$N1N2=($N_S_1[0]+$N_S_2[0])/2;
	$S1N2=($N_S_1[1]+$N_S_2[1])/2;
	$N1S2=($N_S_1[2]+$N_S_2[2])/2;
	$S1S2=($N_S_1[3]+$N_S_2[3])/2;
	$N1=($N_S_1[4]+$N_S_2[4])/2;
	$S1=($N_S_1[5]+$N_S_2[5])/2;
	$N2=($N_S_1[6]+$N_S_2[6])/2;
	$S2=($N_S_1[7]+$N_S_2[7])/2;
	@result=($N1N2,$S1N2,$N1S2,$S1S2,$N1,$S1,$N2,$S2);
	return @result;
}


#########################################################################################
sub N_S_dna { # modified from Wei & Zhang (2015)
	my ($dna,$R)=@_;
	my $codon1;
	my $codon2;
	my $position1;
	my $position2;
	my $N1N2=0;
	my $N1S2=0;
	my $S1N2=0;
	my $S1S2=0; 
	my $length=length($dna);
	my @partial_result;
	my @result;
	
	for ( my $i=2; $i < ($length-2); $i++){ 
		if ($i%3==1){
			$codon1=substr($dna,$i-1,3);
			$codon2=substr($dna,$i-2,3);
			$codon2=reverse($codon2);
			$codon2=~tr/ACGTacgt/TGCAtgca/;
			$position1=1;
			$position2=0;
			@partial_result=possible_change($codon1,$codon2,$position1,$position2,$R);
		} elsif($i%3==2) {
			$codon1=substr($dna,$i-2,3);
			$codon2=substr($dna,$i,3);
			$codon2=reverse($codon2);
			$codon2=~tr/ACGTacgt/TGCAtgca/;
			$position1=2;
			$position2=2;
			@partial_result=possible_change($codon1,$codon2,$position1,$position2,$R);
		} elsif($i%3==0) {
			$codon1=substr($dna,$i,3);
			$codon2=substr($dna,$i-1,3);
			$codon2=reverse($codon2);
			$codon2=~tr/ACGTacgt/TGCAtgca/;
			$position1=0;
			$position2=1;
			@partial_result=possible_change($codon1,$codon2,$position1,$position2,$R);
		}
		$N1N2=$N1N2+$partial_result[0];
		$S1N2=$S1N2+$partial_result[1];
		$N1S2=$N1S2+$partial_result[2];
		$S1S2=$S1S2+$partial_result[3];
	}
	my $N1=$N1N2+$N1S2;
	my $S1=$S1N2+$S1S2;
	my $N2=$N1N2+$S1N2;
	my $S2=$N1S2+$S1S2;
	@result=($N1N2,$S1N2,$N1S2,$S1S2,$N1,$S1,$N2,$S2);
	return @result;
}

#########################################################################################
# DETERMINE THE NUMBER OF DIFFERENCES
sub n_s { # modified from Wei & Zhang (2015)
	my ($dna1,$dna2)=@_;
	my $n1n2=0;
	my $s1n2=0;
	my $n1s2=0;
	my $s1s2=0;
	my $n1n2_pathway=0;
	my $s1n2_pathway=0;
	my $n1s2_pathway=0;
	my $s1s2_pathway=0;
	my $region1=substr($dna1,0,3);
	my $region2=substr($dna2,0,3);
	my @place_change;
	my @region_change;
	my @partial_result=(0,0,0,0);
	my @result;
	my $pathway_number=0;
	my $length=length($dna1);
	my $i=0;
	
	while ($i < ($length-2)) { # We're still within the DNA sequence
		@region_change = &switch_region($region1,$region2,$i,$dna1,$dna2);
		$region1=$region_change[0];
		$region2=$region_change[1];
		$i=$region_change[2]; #
		
#		print "Determining partial result for i $i\...\n";
		@partial_result = &pathway($region1,$region2,$i,$pathway_number,$n1n2_pathway,$s1n2_pathway,$n1s2_pathway,$s1s2_pathway,$n1n2_pathway,$s1n2_pathway,$n1s2_pathway,$s1s2_pathway);
		
		if ($i%3==0) {
			if ((length($region1)%3)==0) {
				$i=$i+length($region1)-1; 
			} elsif((length($region1)%3)==2) {
				$i=$i+length($region1)-2;
			}
		} elsif($i%3==2) {
			if ( (length($region1)%3)==0){
				$i=$i+length($region1)-2;
			} elsif((length($region1)%3)==1) {
				$i=$i+length($region1)-1;
			}
		}
		
		if ($i<($length-2)){
			$region1=substr($dna1,$i,3);
			$region2=substr($dna2,$i,3);
		}
		
		if ($partial_result[8]>0){
			$n1n2=$n1n2+$partial_result[0]/$partial_result[8];
			$s1n2=$s1n2+$partial_result[1]/$partial_result[8];
			$n1s2=$n1s2+$partial_result[2]/$partial_result[8];
			$s1s2=$s1s2+$partial_result[3]/$partial_result[8];
		}
	}
	
	my $n1=$n1n2+$n1s2;
	my $s1=$s1n2+$s1s2;
	my $n2=$n1n2+$s1n2;
	my $s2=$n1s2+$s1s2;
	@result=($n1n2,$s1n2,$n1s2,$s1s2,$n1,$s1,$n2,$s2);
	return @result;
}

#########################################################################################
sub pathway { # modified from Wei & Zhang (2015)
	my ($region1,$region2,$i,$pathway_number,$n1n2,$s1n2,$n1s2,$s1s2,$n1n2_pathway,$s1n2_pathway,$n1s2_pathway,$s1s2_pathway)=@_;
	my $change_number=change_number($region1,$region2);
	my @change_place=change_place($region1,$region2);
	my $place_region;
	my $base;
	my $new_region1;
	my $codon1;
	my $new_codon1;
	my @partial_result;
	my @result;
	my @set;
	my $codon1_region1;
	my $codon1_region2;
	my $codon2_region1;
	my $codon2_region2;
	my $aa1_region1;
	my $aa2_region1;
	my $aa1_region2;
	my $aa2_region2;
	my $length=length($region1);
	for (my $j=0;$j<$change_number;$j++){
		$place_region=$change_place[$j];
		$new_region1=$region1;
		$base=substr($region2,$place_region,1);
		substr($new_region1,$place_region,1,$base);
		if ($i%3==0){ 
			if ($place_region%3==1){
				$codon1_region1=substr($region1,$place_region-1,3); 
				$codon2_region1=substr($region1,$place_region-2,3);
				$codon2_region1=reverse($codon2_region1); 
				$codon2_region1=~tr/ACGTacgt/TGCAtgca/;
				$codon1_region2=substr($new_region1,$place_region-1,3);
				$codon2_region2=substr($new_region1,$place_region-2,3);
				$codon2_region2=reverse($codon2_region2);
				$codon2_region2=~tr/ACGTacgt/TGCAtgca/;
				$aa1_region1=codon2aa($codon1_region1);
				$aa2_region1=codon2aa($codon2_region1);
				$aa1_region2=codon2aa($codon1_region2);
				$aa2_region2=codon2aa($codon2_region2);
				@partial_result=compare_dna1_dna2($aa1_region1,$aa2_region1,$aa1_region2,$aa2_region2); 
			} elsif($place_region%3==2) {
				$codon1_region1=substr($region1,$place_region-2,3); 
				$codon2_region1=substr($region1,$place_region,3); 
				$codon2_region1=reverse($codon2_region1);
				$codon2_region1=~tr/ACGTacgt/TGCAtgca/;
				$codon1_region2=substr($new_region1,$place_region-2,3);  
				$codon2_region2=substr($new_region1,$place_region,3);    
				$codon2_region2=reverse($codon2_region2); 
				$codon2_region2=~tr/ACGTacgt/TGCAtgca/;
				$aa1_region1=codon2aa($codon1_region1);
				$aa2_region1=codon2aa($codon2_region1);
				$aa1_region2=codon2aa($codon1_region2);
				$aa2_region2=codon2aa($codon2_region2);
				@partial_result=compare_dna1_dna2($aa1_region1,$aa2_region1,$aa1_region2,$aa2_region2);
			} elsif($place_region%3==0) {
				$codon1_region1=substr($region1,$place_region,3);
				$codon2_region1=substr($region1,$place_region-1,3); 
				$codon2_region1=reverse($codon2_region1);
				$codon2_region1=~tr/ACGTacgt/TGCAtgca/;
				$codon1_region2=substr($new_region1,$place_region,3);
				$codon2_region2=substr($new_region1,$place_region-1,3);
				$codon2_region2=reverse($codon2_region2);
				$codon2_region2=~tr/ACGTacgt/TGCAtgca/;
				$aa1_region1=codon2aa($codon1_region1);
				$aa2_region1=codon2aa($codon2_region1);
				$aa1_region2=codon2aa($codon1_region2);
				$aa2_region2=codon2aa($codon2_region2);
				@partial_result=compare_dna1_dna2($aa1_region1,$aa2_region1,$aa1_region2,$aa2_region2);
			}
		}
		
		if ($i%3==2){
			if ($place_region%3==1){
				$codon1_region1=substr($region1,$place_region,3);  
				$codon2_region1=substr($region1,$place_region-1,3);
				$codon2_region1=reverse($codon2_region1);
				$codon2_region1=~tr/ACGTacgt/TGCAtgca/;
				$codon1_region2=substr($new_region1,$place_region,3);
				$codon2_region2=substr($new_region1,$place_region-1,3);
				$codon2_region2=reverse($codon2_region2);
				$codon2_region2=~tr/ACGTacgt/TGCAtgca/;
				$aa1_region1=codon2aa($codon1_region1);
				$aa2_region1=codon2aa($codon2_region1);
				$aa1_region2=codon2aa($codon1_region2);
				$aa2_region2=codon2aa($codon2_region2);
				@partial_result=compare_dna1_dna2($aa1_region1,$aa2_region1,$aa1_region2,$aa2_region2);
			} elsif($place_region%3==2) {
				$codon1_region1=substr($region1,$place_region-1,3);
				$codon2_region1=substr($region1,$place_region-2,3);
				$codon2_region1=reverse($codon2_region1);
				$codon2_region1=~tr/ACGTacgt/TGCAtgca/;
				$codon1_region2=substr($new_region1,$place_region-1,3);
				$codon2_region2=substr($new_region1,$place_region-2,3);
				$codon2_region2=reverse($codon2_region2);
				$codon2_region2=~tr/ACGTacgt/TGCAtgca/;
				$aa1_region1=codon2aa($codon1_region1);
				$aa2_region1=codon2aa($codon2_region1);
				$aa1_region2=codon2aa($codon1_region2);
				$aa2_region2=codon2aa($codon2_region2);
				@partial_result=compare_dna1_dna2($aa1_region1,$aa2_region1,$aa1_region2,$aa2_region2);
			} elsif($place_region%3==0) {
				$codon1_region1=substr($region1,$place_region-2,3);
				$codon2_region1=substr($region1,$place_region,3);
				$codon2_region1=reverse($codon2_region1);
				$codon2_region1=~tr/ACGTacgt/TGCAtgca/;
				$codon1_region2=substr($new_region1,$place_region-2,3);
				$codon2_region2=substr($new_region1,$place_region,3);
				$codon2_region2=reverse($codon2_region2);
				$codon2_region2=~tr/ACGTacgt/TGCAtgca/;
				$aa1_region1=codon2aa($codon1_region1);
				$aa2_region1=codon2aa($codon2_region1);
				$aa1_region2=codon2aa($codon1_region2);
				$aa2_region2=codon2aa($codon2_region2);
				@partial_result=compare_dna1_dna2($aa1_region1,$aa2_region1,$aa1_region2,$aa2_region2);
			}
		}
		
		if (($aa1_region2 ne '_') & ($aa2_region2 ne '_')){
			$n1n2_pathway=$n1n2_pathway+$partial_result[0];   
			$s1n2_pathway=$s1n2_pathway+$partial_result[1];
			$n1s2_pathway=$n1s2_pathway+$partial_result[2];
			$s1s2_pathway=$s1s2_pathway+$partial_result[3];
			if ($change_number >1) {
				$region1=$new_region1;
				@result=pathway($region1,$region2,$i,$pathway_number,$n1n2,$s1n2,$n1s2,$s1s2,$n1n2_pathway,$s1n2_pathway,$n1s2_pathway,$s1s2_pathway);
				$n1n2=$result[0];
				$s1n2=$result[1];
				$n1s2=$result[2];
				$s1s2=$result[3];
				$pathway_number=$result[8];
			} else {
				$pathway_number++;
				$n1n2=$n1n2+$n1n2_pathway;
				$s1n2=$s1n2+$s1n2_pathway;
				$n1s2=$n1s2+$n1s2_pathway;
				$s1s2=$s1s2+$s1s2_pathway;
				$n1n2_pathway=0;
				$s1n2_pathway=0;
				$n1s2_pathway=0;
				$s1s2_pathway=0;
			} 
		} else{
			$n1n2_pathway=0;
			$s1n2_pathway=0;
			$n1s2_pathway=0;
			$s1s2_pathway=0;
		}
	}
	@set=($n1n2,$s1n2,$n1s2, $s1s2,$n1n2_pathway,$s1n2_pathway,$n1s2_pathway,$s1s2_pathway,$pathway_number);
	return @set;
}

#########################################################################################
# Find the next region containing differences and return the start index $i
sub switch_region { # modified from Wei & Zhang (2015)
	# Here's what got called: @region_change = &switch_region($region1,$region2,$i,$dna1,$dna2)
	# $i is the index along the full sequence comparison (dna1/dna2)
	my ($region1,$region2,$i,$dna1,$dna2)=@_; #$i is initialized at 0; regions as first codons in the dnas
	my @result;
	my $change_number;
	my $length=length($dna1);  
	my $new_region1;
	my $new_region2;
	my $j;
	my $new_change_number;
	my @compare_result;
	
	$change_number = change_number($region1,$region2); # Takes two strings and returns the number of mismatches. First call is first codons of $dna1 and $dna2
	
	while (($change_number==0) & ($i<($length-3))) { # While there are no differences between the regions and we are before index of the last triplet's first position (which is $length-3)
		# If we're at position 1 of a reference codon
		if (($i%3)==0) { # $i is 0-based. Thus, $i%3==0 when we are at position 1 of a new reference codon (index 0, 3, 6, ...)
			$i=$i+2; # change $i to position 3 of the current ref codon beginning at index $i, i.e., position 3 (reverse strand) of the next overlapping alternative codon
			
			if ($i<($length-2)) { # if the position 3 this ref (next alt) codon is at or before the last codon in sequence (i.e., before the index of the terminal codon's position 2)
				$region1=substr($dna1,$i,3); # next overlapping alt codon fron $dna1
				$region2=substr($dna2,$i,3); # next overlapping alt codon fron $dna2
			}
		
		# If we're at position 3 on reverse strand (position 1 on ref strand) of an alternative codon
		} elsif (($i%3)==2) { # $i is 0-based. Thus, $i%3==2 when we are at **position 3** of a new alternate (sas12) codon (index 2, 5, 8, ...), also reference position 3
			$i=$i+1; # position 2 of the current alt codon ending at index $i (beginning at $i+2; reverse strand), position 1 of the reference codon
			
			if ($i<($length-2)) { # if position 2 this alt codon (position 1 of reference) is the last ref codon in the sequence or before (i.e., before the index of last ref codon's position 2)
				$region1=substr($dna1,$i,3); # next ref codon fron $dna1
				$region2=substr($dna2,$i,3); # next ref codon fron $dna2
			}
			
		}
		
		# This will cause a break from the loop if sequence differences are detected, i.e., $change_number>0
		$change_number = change_number($region1,$region2);
		
	}
	
	if ($change_number > 0) { # Must be the case, since we broke from the previous loop
		# Calculate the first index of the next codon
		$j=$i+3; # $i was always the first position of a (ref or alt) codon with respect to reference strand coordinates
		
		if ($i%3==0) { # the codon was a reference one
			if ($j<($length-1)) { # if the next codon is before the last index in seq
				$new_region1=substr($dna1,$j,2);
				$new_region2=substr($dna2,$j,2);
				$new_change_number=change_number($new_region1,$new_region2);
				@compare_result=compare_new_region($j,$new_region1,$new_region2,$dna1,$dna2);
			}
		}
		
		if ($i%3==2) { # the codon was an alt one
			if ($j<$length) { # if the next codon is at or before the last index in seq
				$new_region1=substr($dna1,$j,1);
				$new_region2=substr($dna2,$j,1);
				$new_change_number=change_number($new_region1,$new_region2);
				@compare_result=compare_new_region($j,$new_region1,$new_region2,$dna1,$dna2);
			}
		}	
	}
	
	if (@compare_result ne 0) {
		$region1=$region1.$compare_result[0];
		$region2=$region2.$compare_result[1];
	}
	
	@result=($region1,$region2,$i); 
}

#########################################################################################
sub compare_new_region { # modified from Wei & Zhang (2015)
	my ($j,$add_region1,$add_region2,$dna1,$dna2)=@_;
	my $change; 
	my $try_region1;
	my $try_region2;
	my $length=length($dna1);
	my @result;
	$change=change_number($add_region1,$add_region2);  
	
	while (($change>0) & ($j<$length-1)) {#########
		if (($j%3)==0) {  
			$try_region1=substr($dna1,$j+2,1);   
			$try_region2=substr($dna2,$j+2,1);
			$change=change_number($try_region1,$try_region2);
			
			if ($change > 0) {
				$j=$j+2;
			}
			
			$add_region1=$add_region1.$try_region1;
			$add_region2=$add_region2.$try_region2;
		}
		
		if (($j%3)==2) {   
			$try_region1=substr($dna1,$j+1,2);
			$try_region2=substr($dna2,$j+1,2);
			$change=change_number($try_region1,$try_region2);
			
			if ($change> 0) {
				$j=$j+1;
			}
			
			$add_region1=$add_region1.$try_region1;
			$add_region2=$add_region2.$try_region2;
		}
	}
	
	@result=($add_region1,$add_region2);
	return @result;
	
}

#########################################################################################
# Takes two strings and returns the number of mismatches
sub change_number { # modified from Wei & Zhang (2015)
	my ($region1,$region2)=@_;
	my $number=0;
	my $base1;
	my $base2;
	
	for (my $i=0; $i<length($region1);$i++){
		$base1=substr($region1,$i,1);
		$base2=substr($region2,$i,1);
		
		if ($base1 ne $base2) {
			$number++;
		}
		
	}
	return $number;
}

#########################################################################################
sub change_place { # modified from Wei & Zhang (2015)
	my ($region1,$region2)=@_;
	my @result;
	my $base1;
	my $base2;
	for (my $i=0; $i<length($region1);$i++){
		$base1=substr($region1,$i,1);
		$base2=substr($region2,$i,1);
		if ($base1 ne $base2){  
			push(@result,$i);
		}
	}
	return @result;  ##
}

#########################################################################################
sub compare_dna1_dna2 { # modified from Wei & Zhang (2015)
	my ($aa1_dna1,$aa2_dna1,$aa1_dna2,$aa2_dna2)=@_;
	my @result=(0,0,0,0);
	if ($aa1_dna1 eq $aa1_dna2) {
		if ($aa2_dna1 eq $aa2_dna2) {
			$result[3]=1;  
		} else {
			$result[1]=1;
		}
	} elsif ($aa2_dna1 eq $aa2_dna2){
		$result[2]=1;
	} else{
		$result[0]=1;
	}
	return @result;
}

#########################################################################################
sub possible_change { # modified from Wei & Zhang (2015)
	my ($codon1,$codon2,$position1,$position2,$R)=@_;
	my $original_base=substr($codon1,$position1,1);
	my $new_codon1=$codon1;
	my $new_codon2=$codon2;
	my $aa1=codon2aa($codon1);
	my $aa2=codon2aa($codon2);
	my $new_aa1;
	my $new_aa2;
	my $N1N2=0;
	my $N1S2=0;
	my $S1N2=0;
	my $S1S2=0;
	my $sum=0;
	
	if ($original_base ne 'A') {
		substr($new_codon1,$position1,1,'A');
		substr($new_codon2,$position2,1,'T');
		$new_aa1=codon2aa($new_codon1);
		$new_aa2=codon2aa($new_codon2);
		if ($aa1 eq $new_aa1) {
			if ($aa2 eq $new_aa2) {
				if ($original_base eq 'T'){
					$S1S2=$S1S2+$R/(1+$R);
				} else {
					$S1S2=$S1S2+1/(2+2*$R);
				}
				
			} elsif($new_aa1 ne '_') {
				if ($original_base eq 'T'){
					$S1N2=$S1N2+$R/(1+$R);
				} else {
					$S1N2=$S1N2+1/(2+2*$R);
				}
			}
			
		} elsif($new_aa1 ne '_') {
			if ($aa2 eq $new_aa2) {
				if ($original_base eq 'T'){
					$N1S2=$N1S2+$R/(1+$R);
				} else {
					$N1S2=$N1S2+1/(2+2*$R);
				}
				
			} elsif($new_aa2 ne '_') {
				if ($original_base eq 'T'){
					$N1N2=$N1N2+$R/(1+$R);
				} else {
					$N1N2=$N1N2+1/(2+2*$R);
				}
			}
		}
	}
	
	if ($original_base ne 'T') {
		substr($new_codon1,$position1,1,'T');
		substr($new_codon2,$position2,1,'A');
		$new_aa1=codon2aa($new_codon1);
		$new_aa2=codon2aa($new_codon2);
		if ($aa1 eq $new_aa1) {
			if ($aa2 eq $new_aa2) {
				if ($original_base eq 'A') {
					$S1S2=$S1S2+$R/(1+$R);
				} else{
					$S1S2=$S1S2+1/(2+2*$R);
				}
				
			} elsif($new_aa1 ne '_') {
				if ($original_base eq 'A') {
					$S1N2=$S1N2+$R/(1+$R);
				} else{
					$S1N2=$S1N2+1/(2+2*$R);
				}
			}
			
		} elsif($new_aa1 ne '_') {
			if ($aa2 eq $new_aa2) {
				if ($original_base eq 'A') {
					$N1S2=$N1S2+$R/(1+$R);
				} else {
					$N1S2=$N1S2+1/(2+2*$R);
				}
				
			} elsif($new_aa2 ne '_') {
				if ($original_base eq 'A') {
					$N1N2=$N1N2+$R/(1+$R);
				} else {
					$N1N2=$N1N2+1/(2+2*$R);
				}
			}
		}
	}
	
	if ($original_base ne 'G') {
		substr($new_codon1,$position1,1,'G');
		substr($new_codon2,$position2,1,'C');
		$new_aa1=codon2aa($new_codon1);
		$new_aa2=codon2aa($new_codon2);
		
		if ($aa1 eq $new_aa1) {
			if ($aa2 eq $new_aa2) {
				if ($original_base eq 'C') {
					$S1S2=$S1S2+$R/(1+$R);
				} else {
					$S1S2=$S1S2+1/(2+2*$R);
				}
				
			} elsif($new_aa1 ne '_') {
				if ($original_base eq 'C') {
					$S1N2=$S1N2+$R/(1+$R);
				} else {
					$S1N2=$S1N2+1/(2+2*$R);
				}
			}
			
		} elsif($new_aa1 ne '_') {
			if ($aa2 eq $new_aa2) {
				if ($original_base eq 'C') {
					$N1S2=$N1S2+$R/(1+$R);
				} else {
					$N1S2=$N1S2+1/(2+2*$R);
				}
				
			} elsif($new_aa2 ne '_') {
				if ($original_base eq 'C') {
					$N1N2=$N1N2+$R/(1+$R);
				} else {
					$N1N2=$N1N2+1/(2+2*$R);
				}
			}
		}
	}
	
	if ($original_base ne 'C') {
		substr($new_codon1,$position1,1,'C');
		substr($new_codon2,$position2,1,'G');
		$new_aa1=codon2aa($new_codon1);
		$new_aa2=codon2aa($new_codon2);
		
		if ($aa1 eq $new_aa1) {
			if ($aa2 eq $new_aa2) {
				if ($original_base eq 'G') {
					$S1S2=$S1S2+$R/(1+$R);
				} else{
					$S1S2=$S1S2+1/(2+2*$R);
				}
			} elsif($new_aa1 ne '_') {
				if ($original_base eq 'G') {
					$S1N2=$S1N2+$R/(1+$R);
				} else{
					$S1N2=$S1N2+1/(2+2*$R);
				}
			}
		} elsif($new_aa1 ne '_') {
			if ($aa2 eq $new_aa2) {
				if ($original_base eq 'G') {
					$N1S2=$N1S2+$R/(1+$R);
				} else {
					$N1S2=$N1S2+1/(2+2*$R);
				}
				
			} elsif($new_aa2 ne '_') {
				if ($original_base eq 'G') {
					$N1N2=$N1N2+$R/(1+$R);
				} else {
					$N1N2=$N1N2+1/(2+2*$R);
				}
			}
		}
	}
	
	$sum=$N1N2+$N1S2+$S1N2+$S1S2;
	
	if ($sum ne 0) {
		$N1N2=$N1N2/$sum;
		$N1S2=$N1S2/$sum;
		$S1N2=$S1N2/$sum;
		$S1S2=$S1S2/$sum;
	}
	
	my @result=($N1N2,$S1N2,$N1S2,$S1S2);
	return @result;
}

#########################################################################################
sub make_overlap_DNA_set { # modified from Wei & Zhang (2015)
	my($number_of_codons,$size_of_set) = @_;     
	my $dna;
	my @set;
	for (my $i = 0; $i < $size_of_set ; ++$i) {
		$dna = make_overlap_DNA ( $number_of_codons ); 
		push(@set, $dna );
	}
	return @set;
}

#########################################################################################
sub make_mutation_set { # modified from Wei & Zhang (2015)
	my ($distance,$R,@random_DNA)=@_;
	my @set;
	foreach my $dna (@random_DNA){  
		$dna=mutation_and_selection($dna,$distance,$R);  
		push(@set,$dna);
	}
	return @set;
}

#########################################################################################
sub randomposition { # modified from Wei & Zhang (2015)
	my($string)= @_;
	return (int(rand(length($string)-4))+2);  
}

#########################################################################################
sub make_overlap_DNA { # modified from Wei & Zhang (2015)
	my ($number_of_codons)= @_;
	my $dna;
	my $codon;
	my $stop_codon='_';
	my $gapped_codon='-'; # I ADDED THIS
	for(my $i=0;$i<$number_of_codons;++$i){  
		$dna.=random_codon();
	}
	for(my $i=2;$i < (length($dna)-2);){
		$codon=substr($dna,$i,3);
		$codon=reverse($codon);
		$codon=~tr/ACGTacgt/TGCAtgca/;
		
		if (codon2aa($codon) eq $stop_codon || codon2aa($codon) eq $gapped_codon) { #ALTER
			substr($dna,$i,3)="";
			$codon=substr($dna,$i-2,3);
			
			if (codon2aa($codon) eq $stop_codon || codon2aa($codon) eq $gapped_codon) { #ALTER
				substr($dna,$i-2,3)="";
			}
			
		} else {
			$i=$i+3;
		}
	}
	return $dna;
}

#########################################################################################
sub randomnucleotide { # modified from Wei & Zhang (2015)
	my(@nucleotides)=('A','C','G','T');
	return randomelement(@nucleotides);
}

#########################################################################################
sub random_codon { # modified from Wei & Zhang (2015)
	my(@codons)=('TCA','TCC','TCG','TCT','TTC','TTT','TTA','TTG','TAC','TAT','TGC','TGT','TGG','CTA','CTC','CTG','CTT','CCA','CCC','CCG','CCT','CAC','CAT','CAA','CAG','CGA','CGC','CGG','CGT','ATA','ATC','ATT','ATG','ACA','ACC','ACG','ACT','AAC','AAT','AAA','AAG','AGC','AGT','AGA','AGG','GTA','GTC','GTG','GTT','GCA','GCC','GCG','GCT','GAC','GAT','GAA','GAG','GGA','GGC','GGG','GGT');
	return randomelement(@codons);
}

#########################################################################################
sub randomelement { # modified from Wei & Zhang (2015)
	my(@array)= @_;
	return $array[rand @array];
}

#########################################################################################
sub codon2aa { # modified from Wei & Zhang (2015)
	my($codon) = @_;
	$codon = uc $codon;    # uc makes the entire string uppercase
	
	my(%genetic_code) = (
		'TCA' => 'S', # Serine
		'TCC' => 'S', # Serine
		'TCG' => 'S', # Serine
		'TCT' => 'S', # Serine
		'TTC' => 'F', # Phenylalanine
		'TTT' => 'F', # Phenylalanine
		'TTA' => 'L', # Leucine
		'TTG' => 'L', # Leucine
		'TAC' => 'Y', # Tyrosine
		'TAT' => 'Y', # Tyrosine
		'TAA' => '_', # Stop
		'TAG' => '_', # Stop
		'TGC' => 'C', # Cysteine
		'TGT' => 'C', # Cysteine
		'TGA' => '_', # Stop
		'TGG' => 'W', # Tryptophan
		'CTA' => 'L', # Leucine
		'CTC' => 'L', # Leucine
		'CTG' => 'L', # Leucine
		'CTT' => 'L', # Leucine
		'CCA' => 'P', # Proline
		'CCC' => 'P', # Proline
		'CCG' => 'P', # Proline
		'CCT' => 'P', # Proline
		'CAC' => 'H', # Histidine
		'CAT' => 'H', # Histidine
		'CAA' => 'Q', # Glutamine
		'CAG' => 'Q', # Glutamine
		'CGA' => 'R', # Arginine
		'CGC' => 'R', # Arginine
		'CGG' => 'R', # Arginine
		'CGT' => 'R', # Arginine
		'ATA' => 'I', # Isoleucine
		'ATC' => 'I', # Isoleucine
		'ATT' => 'I', # Isoleucine
		'ATG' => 'M', # Methionine
		'ACA' => 'T', # Threonine
		'ACC' => 'T', # Threonine
		'ACG' => 'T', # Threonine
		'ACT' => 'T', # Threonine
		'AAC' => 'N', # Asparagine
		'AAT' => 'N', # Asparagine
		'AAA' => 'K', # Lysine
		'AAG' => 'K', # Lysine
		'AGC' => 'S', # Serine
		'AGT' => 'S', # Serine
		'AGA' => 'R', # Arginine
		'AGG' => 'R', # Arginine
		'GTA' => 'V', # Valine
		'GTC' => 'V', # Valine
		'GTG' => 'V', # Valine
		'GTT' => 'V', # Valine
		'GCA' => 'A', # Alanine
		'GCC' => 'A', # Alanine
		'GCG' => 'A', # Alanine
		'GCT' => 'A', # Alanine
		'GAC' => 'D', # Aspartic Acid
		'GAT' => 'D', # Aspartic Acid
		'GAA' => 'E', # Glutamic Acid
		'GAG' => 'E', # Glutamic Acid
		'GGA' => 'G', # Glycine
		'GGC' => 'G', # Glycine
		'GGG' => 'G', # Glycine
		'GGT' => 'G', # Glycine
	);

	if(exists $genetic_code{$codon}) {
		return $genetic_code{$codon};
	} else {
		return '-';
##		print STDERR "Bad codon \"$codon\"!!\n";
##		exit;
	}
}
