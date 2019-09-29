#! /usr/bin/perl
# /nas3/cnelson/bin/anaconda2/bin/perl
# /usr/bin/perl
#/usr/local/software/PERL/perl-5.26.0/bin/perl

# PROGRAM: SNPGenie for overlapping genes (overlapgenie; OLGenie) using codon-based
# analysis for pNN/pSN/pNS/pSS 

# THIS VERSION COMPUTES all pairwise comparisons among a set of sequences given by the 
# fasta, among 6-mers (minimal overlapping unit), then weighting the results.

# This is the FINAL RELEASE of the phylogeny-agnostic version.

### Need a WARNING when FASTA names are identical or not correct.

#########################################################################################
# EXAMPLE CALL:
#########################################################################################
# OLGenie.pl --fasta_file=all_seqs.fasta --frame=sas11 > OLGenie_log.txt
#########################################################################################

#########################################################################################
## LICENSE
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
#########################################################################################

# AUTHOR: Chase W. Nelson
# Copyright (C) 2019 Chase W. Nelson
# DATE CREATED: April 2019

# CONTACT1: cnelson@amnh.org
# CONTACT2: cwnelson88@gmail.com

# AFFILIATION: Sackler Institute for Comparative Genomics, American Museum of Natural
#     History, New York, NY 10024, USA

# CITATION1: OLGenie, https://github.com/chasewnelson/OLGenie
# CITATION2: Nelson CW, Moncla LH, Hughes AL (2015) SNPGenie: estimating evolutionary 
#	parameters to detect natural selection using pooled next-generation sequencing data. 
#	Bioinformatics 31(22):3709-11, doi: 10.1093/bioinformatics/btv449.

use strict;
use Data::Dumper;
use List::Util qw(max);
use Getopt::Long;

# Get the time
my $time1 = time;
my $local_time1 = localtime;

STDOUT->autoflush(1);

my @nucleotides = qw/A C G T/;


#########################################################################################
# INITIALIZE (OPTIONAL) INPUT VARIABLES
my $fasta_file;
my $frame;
my $output_file;
my $verbose;
my $verbose_messages; # only meant for developing code

my $die_message = "\n\n" . 
	"################################################################################\n" .
	"### OLGenie for analysis of selection in overlapping genes using member pairs.\n" .
	"################################################################################\n\n" .
	"\n################################################################################\n" .
	"### OPTIONS:\n" .
	"################################################################################\n" .
	"\n" .
	"\t--fasta_file (REQUIRED): a FASTA file containing multiple aligned sequences of one coding sequence.\n" .
		"\t\tThe entire coding sequence must be an overlapping gene (OLG), with no non-overlapping codons.\n" .
		"\t\tThe frame must be the frame of the reference gene (ORF1). (See the --frame option.)\n\n" .
	"\t--frame (REQUIRED): the frame of the overlapping gene (OLG) relationship: ss12, ss13, sas11, sas12, or sas13:\n\n" .
		"\t\tSENSE-SENSE:\n" .
		"\t\t\tss12:\n\t\t\tORF1: 1-2-3-1-2-3-1\n\t\t\tORF2: 2-3-1-2-3-1-2\n" .
		"\t\t\tss13:\n\t\t\tORF1: 1-2-3-1-2-3-1\n\t\t\tORF2: 3-1-2-3-1-2-3\n\n ".
		"\t\tSENSE-ANTISENSE:\n" .
		"\t\t\tsas11:\n\t\t\tORF1: 1-2-3-1-2-3-1\n\t\t\tORF2: 1-3-2-1-3-2-1\n" .
		"\t\t\tsas12:\n\t\t\tORF1: 1-2-3-1-2-3-1\n\t\t\tORF2: 2-1-3-2-1-3-2\n" .
		"\t\t\tsas13:\n\t\t\tORF1: 1-2-3-1-2-3-1\n\t\t\tORF2: 3-2-1-3-2-1-3\n\n\n" .
	"\t--output_file (OPTIONAL): name of the TAB-delimited output file to be placed in the working directory\n" .
		"\t\tunless a full path name is given. If not specified, a file will be printed in the working directory\n" .
		"\t\tby the name OLGenie_codon_results.txt (DEFAULT).\n\n" .
	"\t--verbose (OPTIONAL): tell OLGenie to report all unique nonamers (9-mers) overlapping each reference\n" .
		"\t\tcodon, along with their counts, in the output file. May lead to large output files in cases with\n" .
		"\t\tmany and/or divergent sequences.\n\n" .
	"\n################################################################################\n" .
	"### EXAMPLE:\n" .
	"################################################################################\n\n" .
	"\t\$ OLGenie.pl --fasta_file=my_alignment.fasta --frame=ss13 --output_file=OLGenie_codon_results.txt --verbose\n\n" .
	"################################################################################\n\n";


# Get user input, if given. If a Boolean argument is passed, its value is 1; else undef
GetOptions( "fasta_file=s" => \$fasta_file,
			"frame=s" => \$frame,
			"output_file=s" => \$output_file,
			"verbose" => \$verbose,
			"verbose_messages" => \$verbose_messages )
			
			or die $die_message;
			# If an argument is called as a flag, its value is 0; if not called, it's null

unless(-f "$fasta_file" && 
	($frame eq 'ss12' || $frame eq 'ss13' || $frame eq 'sas11' || $frame eq 'sas12' || $frame eq 'sas13')) {
	
	die $die_message;
}

my $fasta_file_short = $fasta_file;
$fasta_file_short =~ s/(.*\/)?(.+\.\w+)/$2/;

unless("$output_file") {
	$output_file = "OLGenie\_codon\_results\.txt";
}

print "\n################################################################################" .
	"\n##                                                                            ##" .
	"\n##                           Naive OLGenie Initiated!                         ##" .
	"\n##                                                                            ##" .
	"\n################################################################################\n";


print "\nOLGenie initiated at local time $local_time1\n";

# Read in the group of sequences from the fasta file
my %header2sequence;
my $seq = '';
#my @seqs_arr;
my $header = '';
my @headers_arr;
my $seq_num = 0;
my $last_seq_length;

open(IN_FASTA, "$fasta_file") or die "Could not open file $fasta_file\n";

print "\n################################################################################";
print "\nRecording coding sequence data for $fasta_file...\n";

while(<IN_FASTA>) {
	chomp;
	if(/>/) {
		if($seq_num == 0) {
			$header = $_;
			$header =~ s/^>//; # get rid of FASTA header indicator
			$header =~ tr/\|/_/; # convert pipes (|) to underscores (_), following IQTree convention
			$header =~ s/\s.*$//; # trim anything including and after the first whitespace
			
			$seq_num ++;
		} else {
			$seq = uc($seq);
			$seq =~ tr/U/T/;
			
			if($header =~ /[^\w^\-^\.]/) {
				die "\n### TAXA NAMES CONTAIN INAPPROPRIATE CHARACTERS IN FASTA FILE: $header.\n" .
					"### Only alphanumeric characters (a-z, A-Z, 0-9), underscores (_), dashes (-), and periods (.) may be used. SCRIPT TERMINATED.\n\n";
			}
			
			$header2sequence{$header} = $seq;
#			push(@seqs_arr, $seq);
			push(@headers_arr, $header);
			
			$header = $_;
			$header =~ s/^>//; # get rid of FASTA header indicator
			$header =~ tr/\|/_/; # convert pipes (|) to underscores (_), following IQTree convention
			$header =~ s/\s.*$//; # trim anything including and after the first whitespace
			
			$seq_num ++;
			
			my $this_seq_length = length($seq);
			
			unless($this_seq_length % 3 == 0) {
				die "\n\n### DIE: Sequences must be a complete set of codons, i.e., the nucleotide length".
					"### must be evenly divisible by 3. Instead, the length is $this_seq_length\. TERMINATED.\n\n";
			}
			
			if($last_seq_length && ($last_seq_length != $this_seq_length)) {
				die "\n\n### DIE: The sequences must be aligned, i.e., must be the same length. TERMINATED.\n\n";
			} else {
				$last_seq_length = $this_seq_length;
				$seq = '';
			}
		}
	} else {
		$seq .= $_;
	}
}

close IN_FASTA;

$seq = uc($seq);
$seq =~ tr/U/T/;

if($header =~ /^[^\w^\-^\.]/) {
	die "\n### TAXA NAMES CONTAIN INAPPROPRIATE CHARACTERS IN FASTA FILE: $header.\n" .
		"### Only alphanumeric characters (a-z, A-Z, 0-9), underscores (_), dashes (-), and periods (.) may be used. SCRIPT TERMINATED.\n\n";
}
$header2sequence{$header} = $seq;

#push(@seqs_arr, $seq);
push(@headers_arr, $header);


##########################################################################################
# Store the sequence index of each FASTA header
#my %header_to_index;
my %header_to_def_len;
foreach my $curr_header (keys %header2sequence) {
	
	# Get seq
	my $seq = $header2sequence{$curr_header};
	
	# Count undefined length
	my $num_N = $seq =~ s/N//g;
	my $num_gap = $seq =~ s/-//g;
	
	my $defined_length = $last_seq_length - $num_N - $num_gap;
	
	$header_to_def_len{$curr_header} = $defined_length;
	
	#print "seq $header has a defined length of $defined_length\n";
}


##########################################################################################
# Print which frame
print "\n################################################################################";
print "\nAnalyzing overlapping gene in frame $frame\...\n";


##########################################################################################
# DETERMINE ALL UNIQUE 9-mers along the alignment
my %unique_9mers;

for(my $i = 0; $i < @headers_arr; $i++) {
	my $curr_header = $headers_arr[$i];
	my $curr_seq = $header2sequence{$curr_header};
	
	# We've already made sure the sequence is a multiple of 3
	for(my $j = 0; $j <= length($curr_seq) - 9; $j += 3) {
		$unique_9mers{$j}->{substr($curr_seq, $j, 9)}++; # increment current 9-mer count
	}
}

#foreach my $site (sort {$a <=> $b} keys %unique_9mers) {
#	print "SITE: " . ($site + 1) . "-" . ($site + 9) . "\n";
#	foreach my $nonamer (sort keys %{$unique_9mers{$site}}) {
#		print "$nonamer\: " . $unique_9mers{$site}->{$nonamer} . "\n";
#	}
#	print "\n";
#}

#exit;

##########################################################################################
# ANALYZE ALL UNIQUE PAIRS, PRINT CODON FILE, AND STORE TOTALS

my %site_diffs_hh;
my %seq2sites;
my %seq_completed;

print "\n\n################################################################################".
	"\nANALYZING ALL UNIQUE SEQUENCE PAIRS\n";

my %codon_data_hh;

# LOOP ALL NONAMER POSITIONS
foreach my $seq_site_index (sort {$a <=> $b} keys %unique_9mers) {
	if($verbose_messages) { print "\n\n###seq_site_index=$seq_site_index\n" }
	
	### NEW PHYLOGENY-NAIVE APPROACH
	my @nonamers_sorted = sort keys %{$unique_9mers{$seq_site_index}};
	my $codon_in_seq = ($seq_site_index + 6) / 3; # start with codon 2 (middle of first nonamer)
	my $comparisons_sum = 0;
	
	my $NN_sites_numerator = 0;
	my $SN_sites_numerator = 0;
	my $NS_sites_numerator = 0;
	my $SS_sites_numerator = 0;
	
	my $NN_diffs_numerator = 0;
	my $SN_diffs_numerator = 0;
	my $NS_diffs_numerator = 0;
	my $SS_diffs_numerator = 0;
	
	my $this_codon_MNV = 'FALSE';
	
	my %ref_codons_to_counts;
	my %alt1_codons_to_counts;
	my %alt2_codons_to_counts;
	
	for(my $nonamer1_index = 0; $nonamer1_index < @nonamers_sorted; $nonamer1_index++) {
		my $nonamer1 = $nonamers_sorted[$nonamer1_index];
		my $nonamer1_count = $unique_9mers{$seq_site_index}->{$nonamer1};
		
		$codon_data_hh{$codon_in_seq}->{nonamers} .= "$nonamer1\:";
		$codon_data_hh{$codon_in_seq}->{nonamer_counts} .= "$nonamer1_count\:";
		
		$ref_codons_to_counts{substr($nonamer1, 3, 3)} += $nonamer1_count; # middle codon
		
		for(my $nonamer2_index = $nonamer1_index; $nonamer2_index < @nonamers_sorted; $nonamer2_index++) { # start with SELF and proceed
			my $nonamer2 = $nonamers_sorted[$nonamer2_index];
			
			if($verbose_messages) { print "\n#nonamer pair: $nonamer1_index vs. $nonamer2_index\n" }
			
			# EXAMINE THE PROVIDED PRODUCT FOR THE CURRENT MEMBER PAIR
			my $nonamer1_len = length($nonamer1);
			my $nonamer2_len = length($nonamer2);
			
			if(($nonamer1_len % 3) != 0) {
				die "\n\nDIE: sequence of $nonamer1 is not a multiple of 3 (complete codon set). TERMINATED.\n\n";
			}
			
			if(($nonamer2_len % 3) != 0) {
				die "\n\nDIE: sequence of in $nonamer2 is not a multiple of 3 (complete codon set). TERMINATED.\n\n";
			}
			
			if($nonamer1_len != $nonamer2_len) {
				die "\n\nDIE: The length of sequence is different in $nonamer1 and $nonamer2\. TERMINATED.\n\n";
			}
			
			my $num_codons = $nonamer1_len / 3;
			
			if($verbose_messages) { print "num_codons=$num_codons\n" }
			
			# STORE CODONS for nonamer1 and nonamer2 here; we've verified they're the same length
			#NONAMER_CODONS: for(my $codon_num = 1; $codon_num <= $num_codons; $codon_num++) { # each ORF1 codon
			my $codon_num = 2; # with the new approach, it's always the central/focal/ref codon
			
			my $codon_index = $codon_num - 1;
			my $site_index = 0 + (3 * $codon_index);
			
			# Extract codons for ORF1
			my $codon_nonamer1_ORF1 = substr($nonamer1, $site_index, 3);
			my $codon_nonamer2_ORF1 = substr($nonamer2, $site_index, 3);
			
			my $AA_nonamer1_ORF1 = get_amino_acid($codon_nonamer1_ORF1);
			my $AA_nonamer2_ORF1 = get_amino_acid($codon_nonamer2_ORF1);
			
			if($codon_num < $num_codons && ($AA_nonamer1_ORF1 eq '*' || $AA_nonamer2_ORF1 eq '*')) {
				print "### WARNING! ORF1, $nonamer1\-$nonamer2 comparison, codon $codon_in_seq encodes a within-frame STOP codon. Wrong frame selection ($frame)?\n";
			}
			
			# IF UNDEFINED CODONS, NOTHING AT ALL CAN BE DETERMINED
			unless($codon_nonamer1_ORF1 =~ /N/ || $codon_nonamer2_ORF1 =~ /N/ || $codon_nonamer1_ORF1 =~ /-/ ||$codon_nonamer2_ORF1 =~ /-/) {
			
				# Extract codons for (OVERLAPPING) ORF2
				
##################################################################################################
##################### SENSE-SENSE:
#####################  ss12:
#####################    ORF1: 1-2-3-1-2-3-1
#####################    ORF2: 2-3-1-2-3-1-2
##################################################################################################
				
				if($frame eq 'ss12') {
					# Last 2 nt of prev codon, first 1 nt of next codon (same strand)
					my $codon_nonamer1_ORF2_prev = substr($nonamer1, ($site_index - 1), 3);
					my $codon_nonamer1_ORF2_next = substr($nonamer1, ($site_index + 2), 3);
					my $codon_nonamer2_ORF2_prev = substr($nonamer2, ($site_index - 1), 3);
					my $codon_nonamer2_ORF2_next = substr($nonamer2, ($site_index + 2), 3);
					# COMEBACK: speed could be improved by re-using last round's 'codon' and '*_next'
					
					if($verbose_messages) {
						print "site_index=$site_index\n";
						print "nonamer1=$nonamer1\n";
						print "nonamer2=$nonamer2\n";
						print "nonamer1_count=$nonamer1_count\n";
						#print "nonamer2=$nonamer2_count\n";s
						print "codon_nonamer1_ORF2_prev=$codon_nonamer1_ORF2_prev\n";
						print "codon_nonamer1_ORF2_next=$codon_nonamer1_ORF2_next\n";
						print "codon_nonamer2_ORF2_prev=$codon_nonamer2_ORF2_prev\n";
						print "codon_nonamer2_ORF2_next=$codon_nonamer2_ORF2_next\n";
					}
					
					my $AA_nonamer1_ORF2_prev = get_amino_acid($codon_nonamer1_ORF2_prev);
					my $AA_nonamer1_ORF2_next = get_amino_acid($codon_nonamer1_ORF2_next);
					my $AA_nonamer2_ORF2_prev = get_amino_acid($codon_nonamer2_ORF2_prev);
					my $AA_nonamer2_ORF2_next = get_amino_acid($codon_nonamer2_ORF2_next);
					
					if($codon_num > 2 && ($AA_nonamer1_ORF2_prev eq '*' || $AA_nonamer1_ORF2_next eq '*' || $AA_nonamer2_ORF2_prev eq '*' || $AA_nonamer2_ORF2_next eq '*')) {
						print "### WARNING! ORF2, $nonamer1\-$nonamer2 comparison, near ORF1 codon $codon_in_seq encodes a within-frame STOP codon. Wrong frame selection ($frame)?\n";
					}
					
					#######################
					# GET NUMBER OF SITES
					if($codon_num > 1 && $codon_num < $num_codons) {
						# both ORF2's PREV and NEXT codons fully defined, 3 sites to examine
						
						# New phylogeny-naive approach
						$alt1_codons_to_counts{$codon_nonamer1_ORF2_prev} += $nonamer1_count;
						$alt1_codons_to_counts{$codon_nonamer2_ORF2_prev} += $nonamer1_count; # ADDED
						$alt2_codons_to_counts{$codon_nonamer1_ORF2_next} += $nonamer1_count;
						$alt2_codons_to_counts{$codon_nonamer2_ORF2_next} += $nonamer1_count; # ADDED
						
						my $nonamer1_num_changes_poss_site1 = 0;
						my $nonamer1_num_changes_NN_site1 = 0;
						my $nonamer1_num_changes_SN_site1 = 0;
						my $nonamer1_num_changes_NS_site1 = 0;
						my $nonamer1_num_changes_SS_site1 = 0;
						
						my $nonamer1_num_changes_poss_site2 = 0;
						my $nonamer1_num_changes_NN_site2 = 0;
						my $nonamer1_num_changes_SN_site2 = 0;
						my $nonamer1_num_changes_NS_site2 = 0;
						my $nonamer1_num_changes_SS_site2 = 0;
						
						my $nonamer1_num_changes_poss_site3 = 0;
						my $nonamer1_num_changes_NN_site3 = 0;
						my $nonamer1_num_changes_SN_site3 = 0;
						my $nonamer1_num_changes_NS_site3 = 0;
						my $nonamer1_num_changes_SS_site3 = 0;
						
						my $nonamer2_num_changes_poss_site1 = 0;
						my $nonamer2_num_changes_NN_site1 = 0;
						my $nonamer2_num_changes_SN_site1 = 0;
						my $nonamer2_num_changes_NS_site1 = 0;
						my $nonamer2_num_changes_SS_site1 = 0;
						
						my $nonamer2_num_changes_poss_site2 = 0;
						my $nonamer2_num_changes_NN_site2 = 0;
						my $nonamer2_num_changes_SN_site2 = 0;
						my $nonamer2_num_changes_NS_site2 = 0;
						my $nonamer2_num_changes_SS_site2 = 0;
						
						my $nonamer2_num_changes_poss_site3 = 0;
						my $nonamer2_num_changes_NN_site3 = 0;
						my $nonamer2_num_changes_SN_site3 = 0;
						my $nonamer2_num_changes_NS_site3 = 0;
						my $nonamer2_num_changes_SS_site3 = 0;
						
						my $NN_diffs = 0;
						my $SN_diffs = 0;
						my $NS_diffs = 0;
						my $SS_diffs = 0;
						
						# Just ORF1 (reference)
						my $nt1_nonamer1_WT = substr($codon_nonamer1_ORF1, 0, 1);
						my $nt2_nonamer1_WT = substr($codon_nonamer1_ORF1, 1, 1);
						my $nt3_nonamer1_WT = substr($codon_nonamer1_ORF1, 2, 1);
						
						my $nt1_nonamer2_WT = substr($codon_nonamer2_ORF1, 0, 1);
						my $nt2_nonamer2_WT = substr($codon_nonamer2_ORF1, 1, 1);
						my $nt3_nonamer2_WT = substr($codon_nonamer2_ORF1, 2, 1);
						
						foreach my $nt (@nucleotides) {
							##################################################################
							# nonamer1
							
							# SITE 1
							# What is each change to ORF1 CODON SITE 1 / ORF2 CODON SITE 2?
							if($nt ne $nt1_nonamer1_WT) { # only one possibility for this site
								 # only one possibility for this site
								
								my $nonamer1_STOP_caused = 0;
								
								my $codon_nonamer1_ORF1_MUT = $codon_nonamer1_ORF1;
								$codon_nonamer1_ORF1_MUT =~ s/[ACGT]([ACGT])([ACGT])/$nt$1$2/;
								
								my $codon_nonamer1_ORF2_prev_MUT = $codon_nonamer1_ORF2_prev;
								$codon_nonamer1_ORF2_prev_MUT =~ s/([ACGT])[ACGT]([ACGT])/$1$nt$2/;
								
								my $nonamer1_ORF1_effect = 'S';
								
								#nonamer1-ORF1
								my $AA_nonamer1_ORF1_MUT = get_amino_acid($codon_nonamer1_ORF1_MUT);
								
								if($AA_nonamer1_ORF1_MUT ne $AA_nonamer1_ORF1) {
									$nonamer1_ORF1_effect = 'N';
								}
								
								if($AA_nonamer1_ORF1 eq '*' || $AA_nonamer1_ORF1_MUT eq '*') {
									$nonamer1_STOP_caused++;
								}
								
								
								my $nonamer1_ORF2_prev_effect = 'S';
								
								#nonamer1-ORF2
								my $AA_nonamer1_ORF2_prev_MUT = get_amino_acid($codon_nonamer1_ORF2_prev_MUT);
								
								if($AA_nonamer1_ORF2_prev_MUT ne $AA_nonamer1_ORF2_prev) {
									$nonamer1_ORF2_prev_effect = 'N';
								}
								
								if($AA_nonamer1_ORF2_prev eq '*' || $AA_nonamer1_ORF2_prev_MUT eq '*') {
									$nonamer1_STOP_caused++;
								}
								
								# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
								unless($nonamer1_STOP_caused > 0) {
									$nonamer1_num_changes_poss_site1++;
									
									if($nonamer1_ORF1_effect eq 'N') {
										if($nonamer1_ORF2_prev_effect eq 'N') {
											$nonamer1_num_changes_NN_site1++;
											
											# ACTUAL DIFF
											if($nt eq $nt1_nonamer2_WT) {
												$NN_diffs++;
												$site_diffs_hh{$codon_num}->{1}->{NN_diffs}++;
											}
											
										} elsif($nonamer1_ORF2_prev_effect eq 'S') {
											$nonamer1_num_changes_NS_site1++;
											
											# ACTUAL DIFF
											if($nt eq $nt1_nonamer2_WT) {
												$NS_diffs++;
												$site_diffs_hh{$codon_num}->{1}->{NS_diffs}++;
											}
										}
									} elsif($nonamer1_ORF1_effect eq 'S') {
										if($nonamer1_ORF2_prev_effect eq 'N') {
											$nonamer1_num_changes_SN_site1++;
											
											# ACTUAL DIFF
											if($nt eq $nt1_nonamer2_WT) {
												$SN_diffs++;
												$site_diffs_hh{$codon_num}->{1}->{SN_diffs}++;
											}
										} elsif($nonamer1_ORF2_prev_effect eq 'S') {
											$nonamer1_num_changes_SS_site1++;
											
											# ACTUAL DIFF
											if($nt eq $nt1_nonamer2_WT) {
												$SS_diffs++;
												$site_diffs_hh{$codon_num}->{1}->{SS_diffs}++;
											}
										}
									}
								}
							} # end member 1 site 1
							
							
							# SITE 2
							# What is each change to ORF1 CODON SITE 2 / ORF2 CODON SITE 3?
							if($nt ne $nt2_nonamer1_WT) { # only one possibility for this site
								
								my $nonamer1_STOP_caused = 0;
								
								my $codon_nonamer1_ORF1_MUT = $codon_nonamer1_ORF1;
								$codon_nonamer1_ORF1_MUT =~ s/([ACGT])[ACGT]([ACGT])/$1$nt$2/;
								
								my $codon_nonamer1_ORF2_prev_MUT = $codon_nonamer1_ORF2_prev;
								$codon_nonamer1_ORF2_prev_MUT =~ s/([ACGT])([ACGT])[ACGT]/$1$2$nt/;
								
								my $nonamer1_ORF1_effect = 'S';
								
								#nonamer1-ORF1
								my $AA_nonamer1_ORF1_MUT = get_amino_acid($codon_nonamer1_ORF1_MUT);
								
								if($AA_nonamer1_ORF1_MUT ne $AA_nonamer1_ORF1) {
									$nonamer1_ORF1_effect = 'N';
								}
								
								if($AA_nonamer1_ORF1 eq '*' || $AA_nonamer1_ORF1_MUT eq '*') {
									$nonamer1_STOP_caused++;
								}
								
								
								my $nonamer1_ORF2_prev_effect = 'S';
								
								#nonamer1-ORF2
								my $AA_nonamer1_ORF2_prev_MUT = get_amino_acid($codon_nonamer1_ORF2_prev_MUT);
								
								if($AA_nonamer1_ORF2_prev_MUT ne $AA_nonamer1_ORF2_prev) {
									$nonamer1_ORF2_prev_effect = 'N';
								}
								
								if($AA_nonamer1_ORF2_prev eq '*' || $AA_nonamer1_ORF2_prev_MUT eq '*') {
									$nonamer1_STOP_caused++;
								}
								
								# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
								unless($nonamer1_STOP_caused > 0) {
									$nonamer1_num_changes_poss_site2++;
									
									if($nonamer1_ORF1_effect eq 'N') {
										if($nonamer1_ORF2_prev_effect eq 'N') {
											$nonamer1_num_changes_NN_site2++;
											
											# ACTUAL DIFF
											if($nt eq $nt2_nonamer2_WT) {
												$NN_diffs++;
												$site_diffs_hh{$codon_num}->{2}->{NN_diffs}++;
											}
											
										} elsif($nonamer1_ORF2_prev_effect eq 'S') {
											$nonamer1_num_changes_NS_site2++;
											
											# ACTUAL DIFF
											if($nt eq $nt2_nonamer2_WT) {
												$NS_diffs++;
												$site_diffs_hh{$codon_num}->{2}->{NS_diffs}++;
											}
										}
									} elsif($nonamer1_ORF1_effect eq 'S') {
										if($nonamer1_ORF2_prev_effect eq 'N') {
											$nonamer1_num_changes_SN_site2++;
											
											# ACTUAL DIFF
											if($nt eq $nt2_nonamer2_WT) {
												$SN_diffs++;
												$site_diffs_hh{$codon_num}->{2}->{SN_diffs}++;
											}
										} elsif($nonamer1_ORF2_prev_effect eq 'S') {
											$nonamer1_num_changes_SS_site2++;
											
											# ACTUAL DIFF
											if($nt eq $nt2_nonamer2_WT) {
												$SS_diffs++;
												$site_diffs_hh{$codon_num}->{2}->{SS_diffs}++;
											}
										}
									}
								}
							} # end member 1 site 2
							
							
							# SITE 3
							# What is each change to ORF1 CODON SITE 3 / ORF2 (next) CODON SITE 1?
							if($nt ne $nt3_nonamer1_WT) { # only one possibility for this site
								
								my $nonamer1_STOP_caused = 0;
								
								my $codon_nonamer1_ORF1_MUT = $codon_nonamer1_ORF1;
								$codon_nonamer1_ORF1_MUT =~ s/([ACGT])([ACGT])[ACGT]/$1$2$nt/;
								
								my $codon_nonamer1_ORF2_next_MUT = $codon_nonamer1_ORF2_next;
								$codon_nonamer1_ORF2_next_MUT =~ s/[ACGT]([ACGT])([ACGT])/$nt$1$2/;
								
								my $nonamer1_ORF1_effect = 'S';
								
								#nonamer1-ORF1
								my $AA_nonamer1_ORF1_MUT = get_amino_acid($codon_nonamer1_ORF1_MUT);
								
								if($AA_nonamer1_ORF1_MUT ne $AA_nonamer1_ORF1) {
									$nonamer1_ORF1_effect = 'N';
								}
								
								if($AA_nonamer1_ORF1 eq '*' || $AA_nonamer1_ORF1_MUT eq '*') {
									$nonamer1_STOP_caused++;
								}
								
								
								my $nonamer1_ORF2_next_effect = 'S';
								
								#nonamer1-ORF2
								my $AA_nonamer1_ORF2_next_MUT = get_amino_acid($codon_nonamer1_ORF2_next_MUT);
								
								if($AA_nonamer1_ORF2_next_MUT ne $AA_nonamer1_ORF2_next) {
									$nonamer1_ORF2_next_effect = 'N';
								}
								
								if($AA_nonamer1_ORF2_next eq '*' || $AA_nonamer1_ORF2_next_MUT eq '*') {
									$nonamer1_STOP_caused++;
								}
								
								# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
								unless($nonamer1_STOP_caused > 0) {
									$nonamer1_num_changes_poss_site3++;
									
									if($nonamer1_ORF1_effect eq 'N') {
										if($nonamer1_ORF2_next_effect eq 'N') {
											$nonamer1_num_changes_NN_site3++;
											
											# ACTUAL DIFF
											if($nt eq $nt3_nonamer2_WT) {
												$NN_diffs++;
												$site_diffs_hh{$codon_num}->{3}->{NN_diffs}++;
											}
											
										} elsif($nonamer1_ORF2_next_effect eq 'S') {
											$nonamer1_num_changes_NS_site3++;
											
											# ACTUAL DIFF
											if($nt eq $nt3_nonamer2_WT) {
												$NS_diffs++;
												$site_diffs_hh{$codon_num}->{3}->{NS_diffs}++;
											}
										}
									} elsif($nonamer1_ORF1_effect eq 'S') {
										if($nonamer1_ORF2_next_effect eq 'N') {
											$nonamer1_num_changes_SN_site3++;
											
											# ACTUAL DIFF
											if($nt eq $nt3_nonamer2_WT) {
												$SN_diffs++;
												$site_diffs_hh{$codon_num}->{3}->{SN_diffs}++;
											}
										} elsif($nonamer1_ORF2_next_effect eq 'S') {
											$nonamer1_num_changes_SS_site3++;
											
											# ACTUAL DIFF
											if($nt eq $nt3_nonamer2_WT) {
												$SS_diffs++;
												$site_diffs_hh{$codon_num}->{3}->{SS_diffs}++;
											}
										}
									}
								}
							} # end member 1 site 3
							
							
							##################################################################
							# nonamer2
							
							# SITE 1
							# What is each change to ORF1 CODON SITE 1 / ORF2 CODON SITE 2?
							if($nt ne $nt1_nonamer2_WT) { # only one possibility for this site
								 # only one possibility for this site
								
								my $nonamer2_STOP_caused = 0;
								
								my $codon_nonamer2_ORF1_MUT = $codon_nonamer2_ORF1;
								$codon_nonamer2_ORF1_MUT =~ s/[ACGT]([ACGT])([ACGT])/$nt$1$2/;
								
								my $codon_nonamer2_ORF2_prev_MUT = $codon_nonamer2_ORF2_prev;
								$codon_nonamer2_ORF2_prev_MUT =~ s/([ACGT])[ACGT]([ACGT])/$1$nt$2/;
								
								my $nonamer2_ORF1_effect = 'S';
								
								#nonamer1-ORF1
								my $AA_nonamer2_ORF1_MUT = get_amino_acid($codon_nonamer2_ORF1_MUT);
								
								if($AA_nonamer2_ORF1_MUT ne $AA_nonamer2_ORF1) {
									$nonamer2_ORF1_effect = 'N';
								}
								
								if($AA_nonamer2_ORF1 eq '*' || $AA_nonamer2_ORF1_MUT eq '*') {
									$nonamer2_STOP_caused++;
								}
								
								
								my $nonamer2_ORF2_prev_effect = 'S';
								
								#nonamer1-ORF2
								my $AA_nonamer2_ORF2_prev_MUT = get_amino_acid($codon_nonamer2_ORF2_prev_MUT);
								
								if($AA_nonamer2_ORF2_prev_MUT ne $AA_nonamer2_ORF2_prev) {
									$nonamer2_ORF2_prev_effect = 'N';
								}
								
								if($AA_nonamer2_ORF2_prev eq '*' || $AA_nonamer2_ORF2_prev_MUT eq '*') {
									$nonamer2_STOP_caused++;
								}
								
								# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
								unless($nonamer2_STOP_caused > 0) {
									$nonamer2_num_changes_poss_site1++;
									
									if($nonamer2_ORF1_effect eq 'N') {
										if($nonamer2_ORF2_prev_effect eq 'N') {
											$nonamer2_num_changes_NN_site1++;
											
											
										} elsif($nonamer2_ORF2_prev_effect eq 'S') {
											$nonamer2_num_changes_NS_site1++;
											
										}
									} elsif($nonamer2_ORF1_effect eq 'S') {
										if($nonamer2_ORF2_prev_effect eq 'N') {
											$nonamer2_num_changes_SN_site1++;
											
										} elsif($nonamer2_ORF2_prev_effect eq 'S') {
											$nonamer2_num_changes_SS_site1++;
											
										}
									}
								}
							} # end member 2 site 1
							
							
							# SITE 2
							# What is each change to ORF1 CODON SITE 2 / ORF2 CODON SITE 3?
							if($nt ne $nt2_nonamer2_WT) { # only one possibility for this site
								
								my $nonamer2_STOP_caused = 0;
								
								my $codon_nonamer2_ORF1_MUT = $codon_nonamer2_ORF1;
								$codon_nonamer2_ORF1_MUT =~ s/([ACGT])[ACGT]([ACGT])/$1$nt$2/;
								
								my $codon_nonamer2_ORF2_prev_MUT = $codon_nonamer2_ORF2_prev;
								$codon_nonamer2_ORF2_prev_MUT =~ s/([ACGT])([ACGT])[ACGT]/$1$2$nt/;
								
								my $nonamer2_ORF1_effect = 'S';
								
								#nonamer1-ORF1
								my $AA_nonamer2_ORF1_MUT = get_amino_acid($codon_nonamer2_ORF1_MUT);
								
								if($AA_nonamer2_ORF1_MUT ne $AA_nonamer2_ORF1) {
									$nonamer2_ORF1_effect = 'N';
								}
								
								if($AA_nonamer2_ORF1 eq '*' || $AA_nonamer2_ORF1_MUT eq '*') {
									$nonamer2_STOP_caused++;
								}
								
								
								my $nonamer2_ORF2_prev_effect = 'S';
								
								#nonamer1-ORF2
								my $AA_nonamer2_ORF2_prev_MUT = get_amino_acid($codon_nonamer2_ORF2_prev_MUT);
								
								if($AA_nonamer2_ORF2_prev_MUT ne $AA_nonamer2_ORF2_prev) {
									$nonamer2_ORF2_prev_effect = 'N';
								}
								
								if($AA_nonamer2_ORF2_prev eq '*' || $AA_nonamer2_ORF2_prev_MUT eq '*') {
									$nonamer2_STOP_caused++;
								}
								
								# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
								unless($nonamer2_STOP_caused > 0) {
									$nonamer2_num_changes_poss_site2++;
									
									if($nonamer2_ORF1_effect eq 'N') {
										if($nonamer2_ORF2_prev_effect eq 'N') {
											$nonamer2_num_changes_NN_site2++;
											
										} elsif($nonamer2_ORF2_prev_effect eq 'S') {
											$nonamer2_num_changes_NS_site2++;
											
										}
									} elsif($nonamer2_ORF1_effect eq 'S') {
										if($nonamer2_ORF2_prev_effect eq 'N') {
											$nonamer2_num_changes_SN_site2++;
											
										} elsif($nonamer2_ORF2_prev_effect eq 'S') {
											$nonamer2_num_changes_SS_site2++;
											
										}
									}
								}
							} # end member 2 site 2
							
							
							# SITE 3
							# What is each change to ORF1 CODON SITE 3 / ORF2 (next) CODON SITE 1?
							if($nt ne $nt3_nonamer2_WT) { # only one possibility for this site
								
								my $nonamer2_STOP_caused = 0;
								
								my $codon_nonamer2_ORF1_MUT = $codon_nonamer2_ORF1;
								$codon_nonamer2_ORF1_MUT =~ s/([ACGT])([ACGT])[ACGT]/$1$2$nt/;
								
								my $codon_nonamer2_ORF2_next_MUT = $codon_nonamer2_ORF2_next;
								$codon_nonamer2_ORF2_next_MUT =~ s/[ACGT]([ACGT])([ACGT])/$nt$1$2/;
								
								my $nonamer2_ORF1_effect = 'S';
								
								#nonamer2-ORF1
								my $AA_nonamer2_ORF1_MUT = get_amino_acid($codon_nonamer2_ORF1_MUT);
								
								if($AA_nonamer2_ORF1_MUT ne $AA_nonamer2_ORF1) {
									$nonamer2_ORF1_effect = 'N';
								}
								
								if($AA_nonamer2_ORF1 eq '*' || $AA_nonamer2_ORF1_MUT eq '*') {
									$nonamer2_STOP_caused++;
								}
								
								
								my $nonamer2_ORF2_next_effect = 'S';
								
								#nonamer2-ORF2
								my $AA_nonamer2_ORF2_next_MUT = get_amino_acid($codon_nonamer2_ORF2_next_MUT);
								
								if($AA_nonamer2_ORF2_next_MUT ne $AA_nonamer2_ORF2_next) {
									$nonamer2_ORF2_next_effect = 'N';
								}
								
								if($AA_nonamer2_ORF2_next eq '*' || $AA_nonamer2_ORF2_next_MUT eq '*') {
									$nonamer2_STOP_caused++;
								}
								
								# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
								unless($nonamer2_STOP_caused > 0) {
									$nonamer2_num_changes_poss_site3++;
									
									if($nonamer2_ORF1_effect eq 'N') {
										if($nonamer2_ORF2_next_effect eq 'N') {
											$nonamer2_num_changes_NN_site3++;
											
										} elsif($nonamer2_ORF2_next_effect eq 'S') {
											$nonamer2_num_changes_NS_site3++;
											
										}
									} elsif($nonamer2_ORF1_effect eq 'S') {
										if($nonamer2_ORF2_next_effect eq 'N') {
											$nonamer2_num_changes_SN_site3++;
											
										} elsif($nonamer2_ORF2_next_effect eq 'S') {
											$nonamer2_num_changes_SS_site3++;
											
										}
									}
								}
							} # end member 2 site 3
						} # end cycling all 4 nucleotides
						
						# TALLY SITE 1
						my $NN_sites_nonamer1_site1 = 'NA';
						my $SN_sites_nonamer1_site1 = 'NA';
						my $NS_sites_nonamer1_site1 = 'NA';
						my $SS_sites_nonamer1_site1 = 'NA';
						
						if($nonamer1_num_changes_poss_site1 > 0) {
							$NN_sites_nonamer1_site1 = $nonamer1_num_changes_NN_site1 / $nonamer1_num_changes_poss_site1;
							$SN_sites_nonamer1_site1 = $nonamer1_num_changes_SN_site1 / $nonamer1_num_changes_poss_site1;
							$NS_sites_nonamer1_site1 = $nonamer1_num_changes_NS_site1 / $nonamer1_num_changes_poss_site1;
							$SS_sites_nonamer1_site1 = $nonamer1_num_changes_SS_site1 / $nonamer1_num_changes_poss_site1;
							
							unless(exists $seq_completed{$nonamer1}) {
								$seq2sites{$nonamer1}->{NN_sites} += $NN_sites_nonamer1_site1;
								$seq2sites{$nonamer1}->{SN_sites} += $SN_sites_nonamer1_site1;
								$seq2sites{$nonamer1}->{NS_sites} += $NS_sites_nonamer1_site1;
								$seq2sites{$nonamer1}->{SS_sites} += $SS_sites_nonamer1_site1;
							}
						}
						
						my $NN_sites_nonamer2_site1 = 'NA';
						my $SN_sites_nonamer2_site1 = 'NA';
						my $NS_sites_nonamer2_site1 = 'NA';
						my $SS_sites_nonamer2_site1 = 'NA';
						
						if($nonamer2_num_changes_poss_site1 > 0) {
							$NN_sites_nonamer2_site1 = $nonamer2_num_changes_NN_site1 / $nonamer2_num_changes_poss_site1;
							$SN_sites_nonamer2_site1 = $nonamer2_num_changes_SN_site1 / $nonamer2_num_changes_poss_site1;
							$NS_sites_nonamer2_site1 = $nonamer2_num_changes_NS_site1 / $nonamer2_num_changes_poss_site1;
							$SS_sites_nonamer2_site1 = $nonamer2_num_changes_SS_site1 / $nonamer2_num_changes_poss_site1;
							
							unless(exists $seq_completed{$nonamer2}) {
								$seq2sites{$nonamer2}->{NN_sites} += $NN_sites_nonamer2_site1;
								$seq2sites{$nonamer2}->{SN_sites} += $SN_sites_nonamer2_site1;
								$seq2sites{$nonamer2}->{NS_sites} += $NS_sites_nonamer2_site1;
								$seq2sites{$nonamer2}->{SS_sites} += $SS_sites_nonamer2_site1;
							}
						}
						
						my $NN_sites_site1 = 'NA';
						my $SN_sites_site1 = 'NA';
						my $NS_sites_site1 = 'NA';
						my $SS_sites_site1 = 'NA';
						
						if($nonamer1_num_changes_poss_site1 > 0 && $nonamer2_num_changes_poss_site1 > 0) {
							$NN_sites_site1 = ($NN_sites_nonamer1_site1 + $NN_sites_nonamer2_site1) / 2;
							$SN_sites_site1 = ($SN_sites_nonamer1_site1 + $SN_sites_nonamer2_site1) / 2;
							$NS_sites_site1 = ($NS_sites_nonamer1_site1 + $NS_sites_nonamer2_site1) / 2;
							$SS_sites_site1 = ($SS_sites_nonamer1_site1 + $SS_sites_nonamer2_site1) / 2;
						} elsif($nonamer1_num_changes_poss_site1 > 0) {
							$NN_sites_site1 = $NN_sites_nonamer1_site1;
							$SN_sites_site1 = $SN_sites_nonamer1_site1;
							$NS_sites_site1 = $NS_sites_nonamer1_site1;
							$SS_sites_site1 = $SS_sites_nonamer1_site1;
						} elsif($nonamer2_num_changes_poss_site1 > 0) {
							$NN_sites_site1 = $NN_sites_nonamer2_site1;
							$SN_sites_site1 = $SN_sites_nonamer2_site1;
							$NS_sites_site1 = $NS_sites_nonamer2_site1;
							$SS_sites_site1 = $SS_sites_nonamer2_site1;
							
						} # else nothing defined, stay NA
						
						# TALLY SITE 2
						my $NN_sites_nonamer1_site2 = 'NA';
						my $SN_sites_nonamer1_site2 = 'NA';
						my $NS_sites_nonamer1_site2 = 'NA';
						my $SS_sites_nonamer1_site2 = 'NA';
						
						if($nonamer1_num_changes_poss_site2 > 0) {
							$NN_sites_nonamer1_site2 = $nonamer1_num_changes_NN_site2 / $nonamer1_num_changes_poss_site2;
							$SN_sites_nonamer1_site2 = $nonamer1_num_changes_SN_site2 / $nonamer1_num_changes_poss_site2;
							$NS_sites_nonamer1_site2 = $nonamer1_num_changes_NS_site2 / $nonamer1_num_changes_poss_site2;
							$SS_sites_nonamer1_site2 = $nonamer1_num_changes_SS_site2 / $nonamer1_num_changes_poss_site2;
							
							unless(exists $seq_completed{$nonamer1}) {
								$seq2sites{$nonamer1}->{NN_sites} += $NN_sites_nonamer1_site2;
								$seq2sites{$nonamer1}->{SN_sites} += $SN_sites_nonamer1_site2;
								$seq2sites{$nonamer1}->{NS_sites} += $NS_sites_nonamer1_site2;
								$seq2sites{$nonamer1}->{SS_sites} += $SS_sites_nonamer1_site2;
							}
						}
						
						my $NN_sites_nonamer2_site2 = 'NA';
						my $SN_sites_nonamer2_site2 = 'NA';
						my $NS_sites_nonamer2_site2 = 'NA';
						my $SS_sites_nonamer2_site2 = 'NA';
						
						if($nonamer2_num_changes_poss_site2 > 0) {
							$NN_sites_nonamer2_site2 = $nonamer2_num_changes_NN_site2 / $nonamer2_num_changes_poss_site2;
							$SN_sites_nonamer2_site2 = $nonamer2_num_changes_SN_site2 / $nonamer2_num_changes_poss_site2;
							$NS_sites_nonamer2_site2 = $nonamer2_num_changes_NS_site2 / $nonamer2_num_changes_poss_site2;
							$SS_sites_nonamer2_site2 = $nonamer2_num_changes_SS_site2 / $nonamer2_num_changes_poss_site2;
							
							unless(exists $seq_completed{$nonamer2}) {
								$seq2sites{$nonamer2}->{NN_sites} += $NN_sites_nonamer2_site2;
								$seq2sites{$nonamer2}->{SN_sites} += $SN_sites_nonamer2_site2;
								$seq2sites{$nonamer2}->{NS_sites} += $NS_sites_nonamer2_site2;
								$seq2sites{$nonamer2}->{SS_sites} += $SS_sites_nonamer2_site2;
							}
						}
						
						my $NN_sites_site2 = 'NA';
						my $SN_sites_site2 = 'NA';
						my $NS_sites_site2 = 'NA';
						my $SS_sites_site2 = 'NA';
						
						if($nonamer1_num_changes_poss_site2 > 0 && $nonamer2_num_changes_poss_site2 > 0) {
							$NN_sites_site2 = ($NN_sites_nonamer1_site2 + $NN_sites_nonamer2_site2) / 2;
							$SN_sites_site2 = ($SN_sites_nonamer1_site2 + $SN_sites_nonamer2_site2) / 2;
							$NS_sites_site2 = ($NS_sites_nonamer1_site2 + $NS_sites_nonamer2_site2) / 2;
							$SS_sites_site2 = ($SS_sites_nonamer1_site2 + $SS_sites_nonamer2_site2) / 2;
						} elsif($nonamer1_num_changes_poss_site2 > 0) {
							$NN_sites_site2 = $NN_sites_nonamer1_site2;
							$SN_sites_site2 = $SN_sites_nonamer1_site2;
							$NS_sites_site2 = $NS_sites_nonamer1_site2;
							$SS_sites_site2 = $SS_sites_nonamer1_site2;
						} elsif($nonamer2_num_changes_poss_site2 > 0) {
							$NN_sites_site2 = $NN_sites_nonamer2_site2;
							$SN_sites_site2 = $SN_sites_nonamer2_site2;
							$NS_sites_site2 = $NS_sites_nonamer2_site2;
							$SS_sites_site2 = $SS_sites_nonamer2_site2;
							
						} # else nothing defined, stay NA
						
						# TALLY SITE 3
						my $NN_sites_nonamer1_site3 = 'NA';
						my $SN_sites_nonamer1_site3 = 'NA';
						my $NS_sites_nonamer1_site3 = 'NA';
						my $SS_sites_nonamer1_site3 = 'NA';
						
						if($nonamer1_num_changes_poss_site3 > 0) {
							$NN_sites_nonamer1_site3 = $nonamer1_num_changes_NN_site3 / $nonamer1_num_changes_poss_site3;
							$SN_sites_nonamer1_site3 = $nonamer1_num_changes_SN_site3 / $nonamer1_num_changes_poss_site3;
							$NS_sites_nonamer1_site3 = $nonamer1_num_changes_NS_site3 / $nonamer1_num_changes_poss_site3;
							$SS_sites_nonamer1_site3 = $nonamer1_num_changes_SS_site3 / $nonamer1_num_changes_poss_site3;
							
							unless(exists $seq_completed{$nonamer1}) {
								$seq2sites{$nonamer1}->{NN_sites} += $NN_sites_nonamer1_site3;
								$seq2sites{$nonamer1}->{SN_sites} += $SN_sites_nonamer1_site3;
								$seq2sites{$nonamer1}->{NS_sites} += $NS_sites_nonamer1_site3;
								$seq2sites{$nonamer1}->{SS_sites} += $SS_sites_nonamer1_site3;
							}
						}
						
						my $NN_sites_nonamer2_site3 = 'NA';
						my $SN_sites_nonamer2_site3 = 'NA';
						my $NS_sites_nonamer2_site3 = 'NA';
						my $SS_sites_nonamer2_site3 = 'NA';
						
						if($nonamer2_num_changes_poss_site3 > 0) {
							$NN_sites_nonamer2_site3 = $nonamer2_num_changes_NN_site3 / $nonamer2_num_changes_poss_site3;
							$SN_sites_nonamer2_site3 = $nonamer2_num_changes_SN_site3 / $nonamer2_num_changes_poss_site3;
							$NS_sites_nonamer2_site3 = $nonamer2_num_changes_NS_site3 / $nonamer2_num_changes_poss_site3;
							$SS_sites_nonamer2_site3 = $nonamer2_num_changes_SS_site3 / $nonamer2_num_changes_poss_site3;
							
							unless(exists $seq_completed{$nonamer2}) {
								$seq2sites{$nonamer2}->{NN_sites} += $NN_sites_nonamer2_site3;
								$seq2sites{$nonamer2}->{SN_sites} += $SN_sites_nonamer2_site3;
								$seq2sites{$nonamer2}->{NS_sites} += $NS_sites_nonamer2_site3;
								$seq2sites{$nonamer2}->{SS_sites} += $SS_sites_nonamer2_site3;
							}
						}
						
						my $NN_sites_site3 = 'NA';
						my $SN_sites_site3 = 'NA';
						my $NS_sites_site3 = 'NA';
						my $SS_sites_site3 = 'NA';
						
						if($nonamer1_num_changes_poss_site3 > 0 && $nonamer2_num_changes_poss_site3 > 0) {
							$NN_sites_site3 = ($NN_sites_nonamer1_site3 + $NN_sites_nonamer2_site3) / 2;
							$SN_sites_site3 = ($SN_sites_nonamer1_site3 + $SN_sites_nonamer2_site3) / 2;
							$NS_sites_site3 = ($NS_sites_nonamer1_site3 + $NS_sites_nonamer2_site3) / 2;
							$SS_sites_site3 = ($SS_sites_nonamer1_site3 + $SS_sites_nonamer2_site3) / 2;
						} elsif($nonamer1_num_changes_poss_site3 > 0) {
							$NN_sites_site3 = $NN_sites_nonamer1_site3;
							$SN_sites_site3 = $SN_sites_nonamer1_site3;
							$NS_sites_site3 = $NS_sites_nonamer1_site3;
							$SS_sites_site3 = $SS_sites_nonamer1_site3;
						} elsif($nonamer2_num_changes_poss_site3 > 0) {
							$NN_sites_site3 = $NN_sites_nonamer2_site3;
							$SN_sites_site3 = $SN_sites_nonamer2_site3;
							$NS_sites_site3 = $NS_sites_nonamer2_site3;
							$SS_sites_site3 = $SS_sites_nonamer2_site3;
							
						} # else nothing defined, stay NA
						
						# SUM THE THREE SITES
						my $NN_sites = $NN_sites_site1 + $NN_sites_site2 + $NN_sites_site3;
						my $SN_sites = $SN_sites_site1 + $SN_sites_site2 + $SN_sites_site3;
						my $NS_sites = $NS_sites_site1 + $NS_sites_site2 + $NS_sites_site3;
						my $SS_sites = $SS_sites_site1 + $SS_sites_site2 + $SS_sites_site3;
						
						# Check if there are multiple variants in these overlapping codons
						my $MNV = 'FALSE';
						if(($NN_diffs + $SN_diffs + $NS_diffs + $SS_diffs) > 1) {
							$MNV = 'TRUE';
							$this_codon_MNV = 'TRUE'; # only needs to be TRUE in one comparison
						}
						
						# NEW NAIVE APPROACH: ADD TO SUMS
						my $comparison_weight = 0; # number of pairwise comparisons involving these nonamers
						if($nonamer1_index == $nonamer2_index) { # if($nonamer1 eq $nonamer2) { # it's a self-comparison
							# Combination
							$comparison_weight = (($nonamer1_count * $nonamer1_count) - $nonamer1_count) / 2; # nonamer 1 and 2 are the same
						} else {
							# m * n
							my $nonamer2_count = $unique_9mers{$seq_site_index}->{$nonamer2};
							$comparison_weight = $nonamer1_count * $nonamer2_count;
						}
						
						$comparisons_sum += $comparison_weight;
						
						$NN_sites_numerator += $NN_sites * $comparison_weight;
						$SN_sites_numerator += $SN_sites * $comparison_weight;
						$NS_sites_numerator += $NS_sites * $comparison_weight;
						$SS_sites_numerator += $SS_sites * $comparison_weight;
						
						$NN_diffs_numerator += $NN_diffs * $comparison_weight;
						$SN_diffs_numerator += $SN_diffs * $comparison_weight;
						$NS_diffs_numerator += $NS_diffs * $comparison_weight;
						$SS_diffs_numerator += $SS_diffs * $comparison_weight;
						
					} # MIDDLE codon (not first or last; two ORF2 codons overlap
					
					
##################################################################################################
##################### SENSE-SENSE:
#####################  ss13:
#####################    ORF1: 1-2-3-1-2-3-1
#####################    ORF2: 3-1-2-3-1-2-3
##################################################################################################
				} elsif($frame eq 'ss13') {
					
					# Last 1 nt of prev codon, first 2 nt of next codon (same strand)
					my $codon_nonamer1_ORF2_prev = substr($nonamer1, ($site_index - 2), 3);
					my $codon_nonamer1_ORF2_next = substr($nonamer1, ($site_index + 1), 3);
					my $codon_nonamer2_ORF2_prev = substr($nonamer2, ($site_index - 2), 3);
					my $codon_nonamer2_ORF2_next = substr($nonamer2, ($site_index + 1), 3);
					
					if($verbose_messages) {
						print "site_index=$site_index\n";
						print "nonamer1=$nonamer1\n";
						print "nonamer2=$nonamer2\n";
						print "nonamer1_count=$nonamer1_count\n";
						#print "nonamer2=$nonamer2_count\n";s
						print "codon_nonamer1_ORF2_prev=$codon_nonamer1_ORF2_prev\n";
						print "codon_nonamer1_ORF2_next=$codon_nonamer1_ORF2_next\n";
						print "codon_nonamer2_ORF2_prev=$codon_nonamer2_ORF2_prev\n";
						print "codon_nonamer2_ORF2_next=$codon_nonamer2_ORF2_next\n";
					}
					
					my $AA_nonamer1_ORF2_prev = get_amino_acid($codon_nonamer1_ORF2_prev);
					my $AA_nonamer1_ORF2_next = get_amino_acid($codon_nonamer1_ORF2_next);
					my $AA_nonamer2_ORF2_prev = get_amino_acid($codon_nonamer2_ORF2_prev);
					my $AA_nonamer2_ORF2_next = get_amino_acid($codon_nonamer2_ORF2_next);
					
					if($codon_num > 2 && ($AA_nonamer1_ORF2_prev eq '*' || $AA_nonamer1_ORF2_next eq '*' || $AA_nonamer2_ORF2_prev eq '*' || $AA_nonamer2_ORF2_next eq '*')) {
						print "### WARNING! ORF2, $nonamer1\-$nonamer2 comparison, near ORF1 codon $codon_in_seq encodes a within-frame STOP codon. Wrong frame selection ($frame)?\n";
					}
					
					#######################
					# GET NUMBER OF SITES
					if($codon_num > 1 && $codon_num < $num_codons) {
						# both ORF2's PREV and NEXT codons fully defined, 3 sites to examine
						
						# New phylogeny-naive approach
						$alt1_codons_to_counts{$codon_nonamer1_ORF2_prev} += $nonamer1_count;
						$alt1_codons_to_counts{$codon_nonamer2_ORF2_prev} += $nonamer1_count; # ADDED
						$alt2_codons_to_counts{$codon_nonamer1_ORF2_next} += $nonamer1_count;
						$alt2_codons_to_counts{$codon_nonamer2_ORF2_next} += $nonamer1_count; # ADDED
					
						my $nonamer1_num_changes_poss_site1 = 0;
						my $nonamer1_num_changes_NN_site1 = 0;
						my $nonamer1_num_changes_SN_site1 = 0;
						my $nonamer1_num_changes_NS_site1 = 0;
						my $nonamer1_num_changes_SS_site1 = 0;
						
						my $nonamer1_num_changes_poss_site2 = 0;
						my $nonamer1_num_changes_NN_site2 = 0;
						my $nonamer1_num_changes_SN_site2 = 0;
						my $nonamer1_num_changes_NS_site2 = 0;
						my $nonamer1_num_changes_SS_site2 = 0;
						
						my $nonamer1_num_changes_poss_site3 = 0;
						my $nonamer1_num_changes_NN_site3 = 0;
						my $nonamer1_num_changes_SN_site3 = 0;
						my $nonamer1_num_changes_NS_site3 = 0;
						my $nonamer1_num_changes_SS_site3 = 0;
						
						my $nonamer2_num_changes_poss_site1 = 0;
						my $nonamer2_num_changes_NN_site1 = 0;
						my $nonamer2_num_changes_SN_site1 = 0;
						my $nonamer2_num_changes_NS_site1 = 0;
						my $nonamer2_num_changes_SS_site1 = 0;
						
						my $nonamer2_num_changes_poss_site2 = 0;
						my $nonamer2_num_changes_NN_site2 = 0;
						my $nonamer2_num_changes_SN_site2 = 0;
						my $nonamer2_num_changes_NS_site2 = 0;
						my $nonamer2_num_changes_SS_site2 = 0;
						
						my $nonamer2_num_changes_poss_site3 = 0;
						my $nonamer2_num_changes_NN_site3 = 0;
						my $nonamer2_num_changes_SN_site3 = 0;
						my $nonamer2_num_changes_NS_site3 = 0;
						my $nonamer2_num_changes_SS_site3 = 0;
						
						my $NN_diffs = 0;
						my $SN_diffs = 0;
						my $NS_diffs = 0;
						my $SS_diffs = 0;
						
						# Just ORF1 (reference)
						my $nt1_nonamer1_WT = substr($codon_nonamer1_ORF1, 0, 1);
						my $nt2_nonamer1_WT = substr($codon_nonamer1_ORF1, 1, 1);
						my $nt3_nonamer1_WT = substr($codon_nonamer1_ORF1, 2, 1);
						
						my $nt1_nonamer2_WT = substr($codon_nonamer2_ORF1, 0, 1);
						my $nt2_nonamer2_WT = substr($codon_nonamer2_ORF1, 1, 1);
						my $nt3_nonamer2_WT = substr($codon_nonamer2_ORF1, 2, 1);
						
						
						foreach my $nt (@nucleotides) {
							##################################################################
							# nonamer1
							
							# SITE 1
							# What is each change to ORF1 CODON SITE 1 / ORF2 CODON SITE 3?
							if($nt ne $nt1_nonamer1_WT) { # only one possibility for this site
								
								my $nonamer1_STOP_caused = 0;
								
								my $codon_nonamer1_ORF1_MUT = $codon_nonamer1_ORF1;
								$codon_nonamer1_ORF1_MUT =~ s/[ACGT]([ACGT])([ACGT])/$nt$1$2/;
								
								my $codon_nonamer1_ORF2_prev_MUT = $codon_nonamer1_ORF2_prev;
								$codon_nonamer1_ORF2_prev_MUT =~ s/([ACGT])([ACGT])[ACGT]/$1$2$nt/;
								
								my $nonamer1_ORF1_effect = 'S';
								
								#nonamer1-ORF1
								my $AA_nonamer1_ORF1_MUT = get_amino_acid($codon_nonamer1_ORF1_MUT);
								
								if($AA_nonamer1_ORF1_MUT ne $AA_nonamer1_ORF1) {
									$nonamer1_ORF1_effect = 'N';
								}
								
								if($AA_nonamer1_ORF1 eq '*' || $AA_nonamer1_ORF1_MUT eq '*') {
									$nonamer1_STOP_caused++;
								}
								
								
								my $nonamer1_ORF2_prev_effect = 'S';
								
								#nonamer1-ORF2
								my $AA_nonamer1_ORF2_prev_MUT = get_amino_acid($codon_nonamer1_ORF2_prev_MUT);
								
								if($AA_nonamer1_ORF2_prev_MUT ne $AA_nonamer1_ORF2_prev) {
									$nonamer1_ORF2_prev_effect = 'N';
								}
								
								if($AA_nonamer1_ORF2_prev eq '*' || $AA_nonamer1_ORF2_prev_MUT eq '*') {
									$nonamer1_STOP_caused++;
								}
								
								# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
								unless($nonamer1_STOP_caused > 0) {
									$nonamer1_num_changes_poss_site1++;
									
									if($nonamer1_ORF1_effect eq 'N') {
										if($nonamer1_ORF2_prev_effect eq 'N') {
											$nonamer1_num_changes_NN_site1++;
											
											# ACTUAL DIFF
											if($nt eq $nt1_nonamer2_WT) {
												$NN_diffs++;
												$site_diffs_hh{$codon_num}->{1}->{NN_diffs}++;
											}
											
										} elsif($nonamer1_ORF2_prev_effect eq 'S') {
											$nonamer1_num_changes_NS_site1++;
											
											# ACTUAL DIFF
											if($nt eq $nt1_nonamer2_WT) {
												$NS_diffs++;
												$site_diffs_hh{$codon_num}->{1}->{NS_diffs}++;
											}
										}
									} elsif($nonamer1_ORF1_effect eq 'S') {
										if($nonamer1_ORF2_prev_effect eq 'N') {
											$nonamer1_num_changes_SN_site1++;
											
											# ACTUAL DIFF
											if($nt eq $nt1_nonamer2_WT) {
												$SN_diffs++;
												$site_diffs_hh{$codon_num}->{1}->{SN_diffs}++;
											}
										} elsif($nonamer1_ORF2_prev_effect eq 'S') {
											$nonamer1_num_changes_SS_site1++;
											
											# ACTUAL DIFF
											if($nt eq $nt1_nonamer2_WT) {
												$SS_diffs++;
												$site_diffs_hh{$codon_num}->{1}->{SS_diffs}++;
											}
										}
									}
								}
							} # end member 1 site 1
							
							
							# SITE 2
							# What is each change to ORF1 CODON SITE 2 / ORF2 CODON SITE 1?
							if($nt ne $nt2_nonamer1_WT) { # only one possibility for this site
								
								my $nonamer1_STOP_caused = 0;
								
								my $codon_nonamer1_ORF1_MUT = $codon_nonamer1_ORF1;
								$codon_nonamer1_ORF1_MUT =~ s/([ACGT])[ACGT]([ACGT])/$1$nt$2/;
								
								my $codon_nonamer1_ORF2_next_MUT = $codon_nonamer1_ORF2_next;
								$codon_nonamer1_ORF2_next_MUT =~ s/[ACGT]([ACGT])([ACGT])/$nt$1$2/;
								
								my $nonamer1_ORF1_effect = 'S';
								
								#nonamer1-ORF1
								my $AA_nonamer1_ORF1_MUT = get_amino_acid($codon_nonamer1_ORF1_MUT);
								
								if($AA_nonamer1_ORF1_MUT ne $AA_nonamer1_ORF1) {
									$nonamer1_ORF1_effect = 'N';
								}
								
								if($AA_nonamer1_ORF1 eq '*' || $AA_nonamer1_ORF1_MUT eq '*') {
									$nonamer1_STOP_caused++;
								}
								
								my $nonamer1_ORF2_next_effect = 'S';
								
								#nonamer1-ORF2
								my $AA_nonamer1_ORF2_next_MUT = get_amino_acid($codon_nonamer1_ORF2_next_MUT);
								
								if($AA_nonamer1_ORF2_next_MUT ne $AA_nonamer1_ORF2_next) {
									$nonamer1_ORF2_next_effect = 'N';
								}
								
								if($AA_nonamer1_ORF2_next eq '*' || $AA_nonamer1_ORF2_next_MUT eq '*') {
									$nonamer1_STOP_caused++;
								}
								
								# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
								unless($nonamer1_STOP_caused > 0) {
									$nonamer1_num_changes_poss_site2++;
									
									if($nonamer1_ORF1_effect eq 'N') {
										if($nonamer1_ORF2_next_effect eq 'N') {
											$nonamer1_num_changes_NN_site2++;
											
											# ACTUAL DIFF
											if($nt eq $nt2_nonamer2_WT) {
												$NN_diffs++;
												$site_diffs_hh{$codon_num}->{2}->{NN_diffs}++;
											}
											
										} elsif($nonamer1_ORF2_next_effect eq 'S') {
											$nonamer1_num_changes_NS_site2++;
											
											# ACTUAL DIFF
											if($nt eq $nt2_nonamer2_WT) {
												$NS_diffs++;
												$site_diffs_hh{$codon_num}->{2}->{NS_diffs}++;
											}
										}
									} elsif($nonamer1_ORF1_effect eq 'S') {
										if($nonamer1_ORF2_next_effect eq 'N') {
											$nonamer1_num_changes_SN_site2++;
											
											# ACTUAL DIFF
											if($nt eq $nt2_nonamer2_WT) {
												$SN_diffs++;
												$site_diffs_hh{$codon_num}->{2}->{SN_diffs}++;
											}
										} elsif($nonamer1_ORF2_next_effect eq 'S') {
											$nonamer1_num_changes_SS_site2++;
											
											# ACTUAL DIFF
											if($nt eq $nt2_nonamer2_WT) {
												$SS_diffs++;
												$site_diffs_hh{$codon_num}->{2}->{SS_diffs}++;
											}
										}
									}
								}
							} # end member 1 site 2
							
							
							# SITE 3
							# What is each change to ORF1 CODON SITE 3 / ORF2 (next) CODON SITE 2?
							if($nt ne $nt3_nonamer1_WT) { # only one possibility for this site
								
								my $nonamer1_STOP_caused = 0;
								
								my $codon_nonamer1_ORF1_MUT = $codon_nonamer1_ORF1;
								$codon_nonamer1_ORF1_MUT =~ s/([ACGT])([ACGT])[ACGT]/$1$2$nt/;
								
								my $codon_nonamer1_ORF2_next_MUT = $codon_nonamer1_ORF2_next;
								$codon_nonamer1_ORF2_next_MUT =~ s/([ACGT])[ACGT]([ACGT])/$1$nt$2/;
								
								my $nonamer1_ORF1_effect = 'S';
								
								#nonamer1-ORF1
								my $AA_nonamer1_ORF1_MUT = get_amino_acid($codon_nonamer1_ORF1_MUT);
								
								if($AA_nonamer1_ORF1_MUT ne $AA_nonamer1_ORF1) {
									$nonamer1_ORF1_effect = 'N';
								}
								
								if($AA_nonamer1_ORF1 eq '*' || $AA_nonamer1_ORF1_MUT eq '*') {
									$nonamer1_STOP_caused++;
								}
								
								
								my $nonamer1_ORF2_next_effect = 'S';
								
								#nonamer1-ORF2
								my $AA_nonamer1_ORF2_next_MUT = get_amino_acid($codon_nonamer1_ORF2_next_MUT);
								
								if($AA_nonamer1_ORF2_next_MUT ne $AA_nonamer1_ORF2_next) {
									$nonamer1_ORF2_next_effect = 'N';
								}
								
								if($AA_nonamer1_ORF2_next eq '*' || $AA_nonamer1_ORF2_next_MUT eq '*') {
									$nonamer1_STOP_caused++;
								}
								
								# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
								unless($nonamer1_STOP_caused > 0) {
									$nonamer1_num_changes_poss_site3++;
									
									if($nonamer1_ORF1_effect eq 'N') {
										if($nonamer1_ORF2_next_effect eq 'N') {
											$nonamer1_num_changes_NN_site3++;
											
											# ACTUAL DIFF
											if($nt eq $nt3_nonamer2_WT) {
												$NN_diffs++;
												$site_diffs_hh{$codon_num}->{3}->{NN_diffs}++;
											}
											
										} elsif($nonamer1_ORF2_next_effect eq 'S') {
											$nonamer1_num_changes_NS_site3++;
											
											# ACTUAL DIFF
											if($nt eq $nt3_nonamer2_WT) {
												$NS_diffs++;
												$site_diffs_hh{$codon_num}->{3}->{NS_diffs}++;
											}
										}
									} elsif($nonamer1_ORF1_effect eq 'S') {
										if($nonamer1_ORF2_next_effect eq 'N') {
											$nonamer1_num_changes_SN_site3++;
											
											# ACTUAL DIFF
											if($nt eq $nt3_nonamer2_WT) {
												$SN_diffs++;
												$site_diffs_hh{$codon_num}->{3}->{SN_diffs}++;
											}
										} elsif($nonamer1_ORF2_next_effect eq 'S') {
											$nonamer1_num_changes_SS_site3++;
											
											# ACTUAL DIFF
											if($nt eq $nt3_nonamer2_WT) {
												$SS_diffs++;
												$site_diffs_hh{$codon_num}->{3}->{SS_diffs}++;
											}
										}
									}
								}
							} # end member 1 site 3
							
							
							##################################################################
							# nonamer2
							
							# SITE 1
							# What is each change to ORF1 CODON SITE 1 / ORF2 CODON SITE 3?
							if($nt ne $nt1_nonamer2_WT) { # only one possibility for this site
								 # only one possibility for this site
								
								my $nonamer2_STOP_caused = 0;
								
								my $codon_nonamer2_ORF1_MUT = $codon_nonamer2_ORF1;
								$codon_nonamer2_ORF1_MUT =~ s/[ACGT]([ACGT])([ACGT])/$nt$1$2/;
								
								my $codon_nonamer2_ORF2_prev_MUT = $codon_nonamer2_ORF2_prev;
								$codon_nonamer2_ORF2_prev_MUT =~ s/([ACGT])([ACGT])[ACGT]/$1$2$nt/;
								
								my $nonamer2_ORF1_effect = 'S';
								
								#nonamer1-ORF1
								my $AA_nonamer2_ORF1_MUT = get_amino_acid($codon_nonamer2_ORF1_MUT);
								
								if($AA_nonamer2_ORF1_MUT ne $AA_nonamer2_ORF1) {
									$nonamer2_ORF1_effect = 'N';
								}
								
								if($AA_nonamer2_ORF1 eq '*' || $AA_nonamer2_ORF1_MUT eq '*') {
									$nonamer2_STOP_caused++;
								}
								
								
								my $nonamer2_ORF2_prev_effect = 'S';
								
								#nonamer1-ORF2
								my $AA_nonamer2_ORF2_prev_MUT = get_amino_acid($codon_nonamer2_ORF2_prev_MUT);
								
								if($AA_nonamer2_ORF2_prev_MUT ne $AA_nonamer2_ORF2_prev) {
									$nonamer2_ORF2_prev_effect = 'N';
								}
								
								if($AA_nonamer2_ORF2_prev eq '*' || $AA_nonamer2_ORF2_prev_MUT eq '*') {
									$nonamer2_STOP_caused++;
								}
								
								# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
								unless($nonamer2_STOP_caused > 0) {
									$nonamer2_num_changes_poss_site1++;
									
									if($nonamer2_ORF1_effect eq 'N') {
										if($nonamer2_ORF2_prev_effect eq 'N') {
											$nonamer2_num_changes_NN_site1++;
											
											
										} elsif($nonamer2_ORF2_prev_effect eq 'S') {
											$nonamer2_num_changes_NS_site1++;
											
										}
									} elsif($nonamer2_ORF1_effect eq 'S') {
										if($nonamer2_ORF2_prev_effect eq 'N') {
											$nonamer2_num_changes_SN_site1++;
											
										} elsif($nonamer2_ORF2_prev_effect eq 'S') {
											$nonamer2_num_changes_SS_site1++;
											
										}
									}
								}
							} # end member 2 site 1
							
							
							# SITE 2
							# What is each change to ORF1 CODON SITE 2 / ORF2 CODON SITE 1?
							if($nt ne $nt2_nonamer2_WT) { # only one possibility for this site
								
								my $nonamer2_STOP_caused = 0;
								
								my $codon_nonamer2_ORF1_MUT = $codon_nonamer2_ORF1;
								$codon_nonamer2_ORF1_MUT =~ s/([ACGT])[ACGT]([ACGT])/$1$nt$2/;
								
								my $codon_nonamer2_ORF2_next_MUT = $codon_nonamer2_ORF2_next;
								$codon_nonamer2_ORF2_next_MUT =~ s/[ACGT]([ACGT])([ACGT])/$nt$1$2/;
								
								my $nonamer2_ORF1_effect = 'S';
								
								#nonamer1-ORF1
								my $AA_nonamer2_ORF1_MUT = get_amino_acid($codon_nonamer2_ORF1_MUT);
								
								if($AA_nonamer2_ORF1_MUT ne $AA_nonamer2_ORF1) {
									$nonamer2_ORF1_effect = 'N';
								}
								
								if($AA_nonamer2_ORF1 eq '*' || $AA_nonamer2_ORF1_MUT eq '*') {
									$nonamer2_STOP_caused++;
								}
								
								
								my $nonamer2_ORF2_next_effect = 'S';
								
								#nonamer1-ORF2
								my $AA_nonamer2_ORF2_next_MUT = get_amino_acid($codon_nonamer2_ORF2_next_MUT);
								
								if($AA_nonamer2_ORF2_next_MUT ne $AA_nonamer2_ORF2_next) {
									$nonamer2_ORF2_next_effect = 'N';
								}
								
								if($AA_nonamer2_ORF2_next eq '*' || $AA_nonamer2_ORF2_next_MUT eq '*') {
									$nonamer2_STOP_caused++;
								}
								
								# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
								unless($nonamer2_STOP_caused > 0) {
									$nonamer2_num_changes_poss_site2++;
									
									if($nonamer2_ORF1_effect eq 'N') {
										if($nonamer2_ORF2_next_effect eq 'N') {
											$nonamer2_num_changes_NN_site2++;
											
										} elsif($nonamer2_ORF2_next_effect eq 'S') {
											$nonamer2_num_changes_NS_site2++;
											
										}
									} elsif($nonamer2_ORF1_effect eq 'S') {
										if($nonamer2_ORF2_next_effect eq 'N') {
											$nonamer2_num_changes_SN_site2++;
											
										} elsif($nonamer2_ORF2_next_effect eq 'S') {
											$nonamer2_num_changes_SS_site2++;
											
										}
									}
								}
							} # end member 2 site 2
							
							
							# SITE 3
							# What is each change to ORF1 CODON SITE 3 / ORF2 (next) CODON SITE 2?
							if($nt ne $nt3_nonamer2_WT) { # only one possibility for this site
								
								my $nonamer2_STOP_caused = 0;
								
								my $codon_nonamer2_ORF1_MUT = $codon_nonamer2_ORF1;
								$codon_nonamer2_ORF1_MUT =~ s/([ACGT])([ACGT])[ACGT]/$1$2$nt/;
								
								my $codon_nonamer2_ORF2_next_MUT = $codon_nonamer2_ORF2_next;
								$codon_nonamer2_ORF2_next_MUT =~ s/([ACGT])[ACGT]([ACGT])/$1$nt$2/;
								
								my $nonamer2_ORF1_effect = 'S';
								
								#nonamer2-ORF1
								my $AA_nonamer2_ORF1_MUT = get_amino_acid($codon_nonamer2_ORF1_MUT);
								
								if($AA_nonamer2_ORF1_MUT ne $AA_nonamer2_ORF1) {
									$nonamer2_ORF1_effect = 'N';
								}
								
								if($AA_nonamer2_ORF1 eq '*' || $AA_nonamer2_ORF1_MUT eq '*') {
									$nonamer2_STOP_caused++;
								}
								
								
								my $nonamer2_ORF2_next_effect = 'S';
								
								#nonamer2-ORF2
								my $AA_nonamer2_ORF2_next_MUT = get_amino_acid($codon_nonamer2_ORF2_next_MUT);
								
								if($AA_nonamer2_ORF2_next_MUT ne $AA_nonamer2_ORF2_next) {
									$nonamer2_ORF2_next_effect = 'N';
								}
								
								if($AA_nonamer2_ORF2_next eq '*' || $AA_nonamer2_ORF2_next_MUT eq '*') {
									$nonamer2_STOP_caused++;
								}
								
								# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
								unless($nonamer2_STOP_caused > 0) {
									$nonamer2_num_changes_poss_site3++;
									
									if($nonamer2_ORF1_effect eq 'N') {
										if($nonamer2_ORF2_next_effect eq 'N') {
											$nonamer2_num_changes_NN_site3++;
											
										} elsif($nonamer2_ORF2_next_effect eq 'S') {
											$nonamer2_num_changes_NS_site3++;
											
										}
									} elsif($nonamer2_ORF1_effect eq 'S') {
										if($nonamer2_ORF2_next_effect eq 'N') {
											$nonamer2_num_changes_SN_site3++;
											
										} elsif($nonamer2_ORF2_next_effect eq 'S') {
											$nonamer2_num_changes_SS_site3++;
											
										}
									}
								}
							} # end member 2 site 3
						} # end cycling all 4 nucleotides
						
						# TALLY SITE 1
						my $NN_sites_nonamer1_site1 = 'NA';
						my $SN_sites_nonamer1_site1 = 'NA';
						my $NS_sites_nonamer1_site1 = 'NA';
						my $SS_sites_nonamer1_site1 = 'NA';
						
						if($nonamer1_num_changes_poss_site1 > 0) {
						
							$NN_sites_nonamer1_site1 = $nonamer1_num_changes_NN_site1 / $nonamer1_num_changes_poss_site1;
							$SN_sites_nonamer1_site1 = $nonamer1_num_changes_SN_site1 / $nonamer1_num_changes_poss_site1;
							$NS_sites_nonamer1_site1 = $nonamer1_num_changes_NS_site1 / $nonamer1_num_changes_poss_site1;
							$SS_sites_nonamer1_site1 = $nonamer1_num_changes_SS_site1 / $nonamer1_num_changes_poss_site1;
							
							unless(exists $seq_completed{$nonamer1}) {
								$seq2sites{$nonamer1}->{NN_sites} += $NN_sites_nonamer1_site1;
								$seq2sites{$nonamer1}->{SN_sites} += $SN_sites_nonamer1_site1;
								$seq2sites{$nonamer1}->{NS_sites} += $NS_sites_nonamer1_site1;
								$seq2sites{$nonamer1}->{SS_sites} += $SS_sites_nonamer1_site1;
							}
						}
						
						my $NN_sites_nonamer2_site1 = 'NA';
						my $SN_sites_nonamer2_site1 = 'NA';
						my $NS_sites_nonamer2_site1 = 'NA';
						my $SS_sites_nonamer2_site1 = 'NA';
						
						if($nonamer2_num_changes_poss_site1 > 0) {
							$NN_sites_nonamer2_site1 = $nonamer2_num_changes_NN_site1 / $nonamer2_num_changes_poss_site1;
							$SN_sites_nonamer2_site1 = $nonamer2_num_changes_SN_site1 / $nonamer2_num_changes_poss_site1;
							$NS_sites_nonamer2_site1 = $nonamer2_num_changes_NS_site1 / $nonamer2_num_changes_poss_site1;
							$SS_sites_nonamer2_site1 = $nonamer2_num_changes_SS_site1 / $nonamer2_num_changes_poss_site1;
							
							unless(exists $seq_completed{$nonamer2}) {
								$seq2sites{$nonamer2}->{NN_sites} += $NN_sites_nonamer2_site1;
								$seq2sites{$nonamer2}->{SN_sites} += $SN_sites_nonamer2_site1;
								$seq2sites{$nonamer2}->{NS_sites} += $NS_sites_nonamer2_site1;
								$seq2sites{$nonamer2}->{SS_sites} += $SS_sites_nonamer2_site1;
							}
						}
						
						my $NN_sites_site1 = 'NA';
						my $SN_sites_site1 = 'NA';
						my $NS_sites_site1 = 'NA';
						my $SS_sites_site1 = 'NA';
						
						if($nonamer1_num_changes_poss_site1 > 0 && $nonamer2_num_changes_poss_site1 > 0) {
							$NN_sites_site1 = ($NN_sites_nonamer1_site1 + $NN_sites_nonamer2_site1) / 2;
							$SN_sites_site1 = ($SN_sites_nonamer1_site1 + $SN_sites_nonamer2_site1) / 2;
							$NS_sites_site1 = ($NS_sites_nonamer1_site1 + $NS_sites_nonamer2_site1) / 2;
							$SS_sites_site1 = ($SS_sites_nonamer1_site1 + $SS_sites_nonamer2_site1) / 2;
						} elsif($nonamer1_num_changes_poss_site1 > 0) {
							$NN_sites_site1 = $NN_sites_nonamer1_site1;
							$SN_sites_site1 = $SN_sites_nonamer1_site1;
							$NS_sites_site1 = $NS_sites_nonamer1_site1;
							$SS_sites_site1 = $SS_sites_nonamer1_site1;
						} elsif($nonamer2_num_changes_poss_site1 > 0) {
							$NN_sites_site1 = $NN_sites_nonamer2_site1;
							$SN_sites_site1 = $SN_sites_nonamer2_site1;
							$NS_sites_site1 = $NS_sites_nonamer2_site1;
							$SS_sites_site1 = $SS_sites_nonamer2_site1;
							
						} # else nothing defined, stay NA
						
						# TALLY SITE 2
						my $NN_sites_nonamer1_site2 = 'NA';
						my $SN_sites_nonamer1_site2 = 'NA';
						my $NS_sites_nonamer1_site2 = 'NA';
						my $SS_sites_nonamer1_site2 = 'NA';
						
						if($nonamer1_num_changes_poss_site2 > 0) {
							$NN_sites_nonamer1_site2 = $nonamer1_num_changes_NN_site2 / $nonamer1_num_changes_poss_site2;
							$SN_sites_nonamer1_site2 = $nonamer1_num_changes_SN_site2 / $nonamer1_num_changes_poss_site2;
							$NS_sites_nonamer1_site2 = $nonamer1_num_changes_NS_site2 / $nonamer1_num_changes_poss_site2;
							$SS_sites_nonamer1_site2 = $nonamer1_num_changes_SS_site2 / $nonamer1_num_changes_poss_site2;
							
							unless(exists $seq_completed{$nonamer1}) {
								$seq2sites{$nonamer1}->{NN_sites} += $NN_sites_nonamer1_site2;
								$seq2sites{$nonamer1}->{SN_sites} += $SN_sites_nonamer1_site2;
								$seq2sites{$nonamer1}->{NS_sites} += $NS_sites_nonamer1_site2;
								$seq2sites{$nonamer1}->{SS_sites} += $SS_sites_nonamer1_site2;
							}
						}
						
						my $NN_sites_nonamer2_site2 = 'NA';
						my $SN_sites_nonamer2_site2 = 'NA';
						my $NS_sites_nonamer2_site2 = 'NA';
						my $SS_sites_nonamer2_site2 = 'NA';
						
						if($nonamer2_num_changes_poss_site2 > 0) {
							$NN_sites_nonamer2_site2 = $nonamer2_num_changes_NN_site2 / $nonamer2_num_changes_poss_site2;
							$SN_sites_nonamer2_site2 = $nonamer2_num_changes_SN_site2 / $nonamer2_num_changes_poss_site2;
							$NS_sites_nonamer2_site2 = $nonamer2_num_changes_NS_site2 / $nonamer2_num_changes_poss_site2;
							$SS_sites_nonamer2_site2 = $nonamer2_num_changes_SS_site2 / $nonamer2_num_changes_poss_site2;
							
							unless(exists $seq_completed{$nonamer2}) {
								$seq2sites{$nonamer2}->{NN_sites} += $NN_sites_nonamer2_site2;
								$seq2sites{$nonamer2}->{SN_sites} += $SN_sites_nonamer2_site2;
								$seq2sites{$nonamer2}->{NS_sites} += $NS_sites_nonamer2_site2;
								$seq2sites{$nonamer2}->{SS_sites} += $SS_sites_nonamer2_site2;
							}
						}
						
						my $NN_sites_site2 = 'NA';
						my $SN_sites_site2 = 'NA';
						my $NS_sites_site2 = 'NA';
						my $SS_sites_site2 = 'NA';
						
						if($nonamer1_num_changes_poss_site2 > 0 && $nonamer2_num_changes_poss_site2 > 0) {
							$NN_sites_site2 = ($NN_sites_nonamer1_site2 + $NN_sites_nonamer2_site2) / 2;
							$SN_sites_site2 = ($SN_sites_nonamer1_site2 + $SN_sites_nonamer2_site2) / 2;
							$NS_sites_site2 = ($NS_sites_nonamer1_site2 + $NS_sites_nonamer2_site2) / 2;
							$SS_sites_site2 = ($SS_sites_nonamer1_site2 + $SS_sites_nonamer2_site2) / 2;
						} elsif($nonamer1_num_changes_poss_site2 > 0) {
							$NN_sites_site2 = $NN_sites_nonamer1_site2;
							$SN_sites_site2 = $SN_sites_nonamer1_site2;
							$NS_sites_site2 = $NS_sites_nonamer1_site2;
							$SS_sites_site2 = $SS_sites_nonamer1_site2;
						} elsif($nonamer2_num_changes_poss_site2 > 0) {
							$NN_sites_site2 = $NN_sites_nonamer2_site2;
							$SN_sites_site2 = $SN_sites_nonamer2_site2;
							$NS_sites_site2 = $NS_sites_nonamer2_site2;
							$SS_sites_site2 = $SS_sites_nonamer2_site2;
							
						} # else nothing defined, stay NA
						
						# TALLY SITE 3
						my $NN_sites_nonamer1_site3 = 'NA';
						my $SN_sites_nonamer1_site3 = 'NA';
						my $NS_sites_nonamer1_site3 = 'NA';
						my $SS_sites_nonamer1_site3 = 'NA';
						
						if($nonamer1_num_changes_poss_site3 > 0) {
							$NN_sites_nonamer1_site3 = $nonamer1_num_changes_NN_site3 / $nonamer1_num_changes_poss_site3;
							$SN_sites_nonamer1_site3 = $nonamer1_num_changes_SN_site3 / $nonamer1_num_changes_poss_site3;
							$NS_sites_nonamer1_site3 = $nonamer1_num_changes_NS_site3 / $nonamer1_num_changes_poss_site3;
							$SS_sites_nonamer1_site3 = $nonamer1_num_changes_SS_site3 / $nonamer1_num_changes_poss_site3;
							
							unless(exists $seq_completed{$nonamer1}) {
								$seq2sites{$nonamer1}->{NN_sites} += $NN_sites_nonamer1_site3;
								$seq2sites{$nonamer1}->{SN_sites} += $SN_sites_nonamer1_site3;
								$seq2sites{$nonamer1}->{NS_sites} += $NS_sites_nonamer1_site3;
								$seq2sites{$nonamer1}->{SS_sites} += $SS_sites_nonamer1_site3;
							}
						}
						
						my $NN_sites_nonamer2_site3 = 'NA';
						my $SN_sites_nonamer2_site3 = 'NA';
						my $NS_sites_nonamer2_site3 = 'NA';
						my $SS_sites_nonamer2_site3 = 'NA';
						
						if($nonamer2_num_changes_poss_site3 > 0) {
							$NN_sites_nonamer2_site3 = $nonamer2_num_changes_NN_site3 / $nonamer2_num_changes_poss_site3;
							$SN_sites_nonamer2_site3 = $nonamer2_num_changes_SN_site3 / $nonamer2_num_changes_poss_site3;
							$NS_sites_nonamer2_site3 = $nonamer2_num_changes_NS_site3 / $nonamer2_num_changes_poss_site3;
							$SS_sites_nonamer2_site3 = $nonamer2_num_changes_SS_site3 / $nonamer2_num_changes_poss_site3;
							
							unless(exists $seq_completed{$nonamer2}) {
								$seq2sites{$nonamer2}->{NN_sites} += $NN_sites_nonamer2_site3;
								$seq2sites{$nonamer2}->{SN_sites} += $SN_sites_nonamer2_site3;
								$seq2sites{$nonamer2}->{NS_sites} += $NS_sites_nonamer2_site3;
								$seq2sites{$nonamer2}->{SS_sites} += $SS_sites_nonamer2_site3;
							}
						}
						
						my $NN_sites_site3 = 'NA';
						my $SN_sites_site3 = 'NA';
						my $NS_sites_site3 = 'NA';
						my $SS_sites_site3 = 'NA';
						
						if($nonamer1_num_changes_poss_site3 > 0 && $nonamer2_num_changes_poss_site3 > 0) {
							$NN_sites_site3 = ($NN_sites_nonamer1_site3 + $NN_sites_nonamer2_site3) / 2;
							$SN_sites_site3 = ($SN_sites_nonamer1_site3 + $SN_sites_nonamer2_site3) / 2;
							$NS_sites_site3 = ($NS_sites_nonamer1_site3 + $NS_sites_nonamer2_site3) / 2;
							$SS_sites_site3 = ($SS_sites_nonamer1_site3 + $SS_sites_nonamer2_site3) / 2;
						} elsif($nonamer1_num_changes_poss_site3 > 0) {
							$NN_sites_site3 = $NN_sites_nonamer1_site3;
							$SN_sites_site3 = $SN_sites_nonamer1_site3;
							$NS_sites_site3 = $NS_sites_nonamer1_site3;
							$SS_sites_site3 = $SS_sites_nonamer1_site3;
						} elsif($nonamer2_num_changes_poss_site3 > 0) {
							$NN_sites_site3 = $NN_sites_nonamer2_site3;
							$SN_sites_site3 = $SN_sites_nonamer2_site3;
							$NS_sites_site3 = $NS_sites_nonamer2_site3;
							$SS_sites_site3 = $SS_sites_nonamer2_site3;
							
						} # else nothing defined, stay NA
						
						# SUM THE TWO SITES
						my $NN_sites = $NN_sites_site1 + $NN_sites_site2 + $NN_sites_site3;
						my $SN_sites = $SN_sites_site1 + $SN_sites_site2 + $SN_sites_site3;
						my $NS_sites = $NS_sites_site1 + $NS_sites_site2 + $NS_sites_site3;
						my $SS_sites = $SS_sites_site1 + $SS_sites_site2 + $SS_sites_site3;
						
						# Check if there are multiple variants in these overlapping codons
						my $MNV = 'FALSE';
						if(($NN_diffs + $SN_diffs + $NS_diffs + $SS_diffs) > 1) {
							$MNV = 'TRUE';
							$this_codon_MNV = 'TRUE'; # only needs to be TRUE in one comparison
						}
						
						# NEW NAIVE APPROACH: ADD TO SUMS
						my $comparison_weight = 0; # number of pairwise comparisons involving these nonamers
						if($nonamer1_index == $nonamer2_index) { # if($nonamer1 eq $nonamer2) { # it's a self-comparison
							# Combination
							$comparison_weight = (($nonamer1_count * $nonamer1_count) - $nonamer1_count) / 2; # nonamer 1 and 2 are the same
						} else {
							# m * n
							my $nonamer2_count = $unique_9mers{$seq_site_index}->{$nonamer2};
							$comparison_weight = $nonamer1_count * $nonamer2_count;
						}
						
						$comparisons_sum += $comparison_weight;
						
						$NN_sites_numerator += $NN_sites * $comparison_weight;
						$SN_sites_numerator += $SN_sites * $comparison_weight;
						$NS_sites_numerator += $NS_sites * $comparison_weight;
						$SS_sites_numerator += $SS_sites * $comparison_weight;
						
						$NN_diffs_numerator += $NN_diffs * $comparison_weight;
						$SN_diffs_numerator += $SN_diffs * $comparison_weight;
						$NS_diffs_numerator += $NS_diffs * $comparison_weight;
						$SS_diffs_numerator += $SS_diffs * $comparison_weight;
						
					} # MIDDLE codon (not first or last; two ORF2 codons overlap
					
##################################################################################################
##################### SENSE-ANTISENSE:
#####################  sas11:
#####################    ORF1: 1-2-3-1-2-3-1
#####################    ORF2: 1-3-2-1-3-2-1
##################################################################################################
				} elsif($frame eq 'sas11') {
					
					# First 1 nt of prev codon, last 2 nt of next codon (opposite strand)
					my $codon_nonamer1_ORF2_prev = revcom(substr($nonamer1, ($site_index - 2), 3)); # remember, looking 3'->5'
					my $codon_nonamer1_ORF2_next = revcom(substr($nonamer1, ($site_index + 1), 3));
					my $codon_nonamer2_ORF2_prev = revcom(substr($nonamer2, ($site_index - 2), 3));
					my $codon_nonamer2_ORF2_next = revcom(substr($nonamer2, ($site_index + 1), 3));
					
					if($verbose_messages) {
						print "site_index=$site_index\n";
						print "nonamer1=$nonamer1\n";
						print "nonamer2=$nonamer2\n";
						print "nonamer1_count=$nonamer1_count\n";
						#print "nonamer2=$nonamer2_count\n";s
						print "codon_nonamer1_ORF2_prev=$codon_nonamer1_ORF2_prev\n";
						print "codon_nonamer1_ORF2_next=$codon_nonamer1_ORF2_next\n";
						print "codon_nonamer2_ORF2_prev=$codon_nonamer2_ORF2_prev\n";
						print "codon_nonamer2_ORF2_next=$codon_nonamer2_ORF2_next\n";
					}
					
					my $AA_nonamer1_ORF2_prev = get_amino_acid($codon_nonamer1_ORF2_prev);
					my $AA_nonamer1_ORF2_next = get_amino_acid($codon_nonamer1_ORF2_next);
					my $AA_nonamer2_ORF2_prev = get_amino_acid($codon_nonamer2_ORF2_prev);
					my $AA_nonamer2_ORF2_next = get_amino_acid($codon_nonamer2_ORF2_next);
					
					if($codon_num > 2 && ($AA_nonamer1_ORF2_prev eq '*' || $AA_nonamer1_ORF2_next eq '*' || $AA_nonamer2_ORF2_prev eq '*' || $AA_nonamer2_ORF2_next eq '*')) {
						print "### WARNING! ORF2, $nonamer1\-$nonamer2 comparison, near ORF1 codon $codon_in_seq encodes a within-frame STOP codon. Wrong frame selection ($frame)?\n";
					}
					
					#######################
					# GET NUMBER OF SITES
					
					if($codon_num > 1 && $codon_num < $num_codons) {
						# both ORF2's PREV and NEXT codons fully defined, 3 sites to examine
						
						# New phylogeny-naive approach
						$alt1_codons_to_counts{$codon_nonamer1_ORF2_prev} += $nonamer1_count;
						$alt1_codons_to_counts{$codon_nonamer2_ORF2_prev} += $nonamer1_count; # ADDED
						$alt2_codons_to_counts{$codon_nonamer1_ORF2_next} += $nonamer1_count;
						$alt2_codons_to_counts{$codon_nonamer2_ORF2_next} += $nonamer1_count; # ADDED
						
						my $nonamer1_num_changes_poss_site1 = 0;
						my $nonamer1_num_changes_NN_site1 = 0;
						my $nonamer1_num_changes_SN_site1 = 0;
						my $nonamer1_num_changes_NS_site1 = 0;
						my $nonamer1_num_changes_SS_site1 = 0;
						
						my $nonamer1_num_changes_poss_site2 = 0;
						my $nonamer1_num_changes_NN_site2 = 0;
						my $nonamer1_num_changes_SN_site2 = 0;
						my $nonamer1_num_changes_NS_site2 = 0;
						my $nonamer1_num_changes_SS_site2 = 0;
						
						my $nonamer1_num_changes_poss_site3 = 0;
						my $nonamer1_num_changes_NN_site3 = 0;
						my $nonamer1_num_changes_SN_site3 = 0;
						my $nonamer1_num_changes_NS_site3 = 0;
						my $nonamer1_num_changes_SS_site3 = 0;
						
						my $nonamer2_num_changes_poss_site1 = 0;
						my $nonamer2_num_changes_NN_site1 = 0;
						my $nonamer2_num_changes_SN_site1 = 0;
						my $nonamer2_num_changes_NS_site1 = 0;
						my $nonamer2_num_changes_SS_site1 = 0;
						
						my $nonamer2_num_changes_poss_site2 = 0;
						my $nonamer2_num_changes_NN_site2 = 0;
						my $nonamer2_num_changes_SN_site2 = 0;
						my $nonamer2_num_changes_NS_site2 = 0;
						my $nonamer2_num_changes_SS_site2 = 0;
						
						my $nonamer2_num_changes_poss_site3 = 0;
						my $nonamer2_num_changes_NN_site3 = 0;
						my $nonamer2_num_changes_SN_site3 = 0;
						my $nonamer2_num_changes_NS_site3 = 0;
						my $nonamer2_num_changes_SS_site3 = 0;
						
						my $NN_diffs = 0;
						my $SN_diffs = 0;
						my $NS_diffs = 0;
						my $SS_diffs = 0;
						
						# Just ORF1 (reference)
						my $nt1_nonamer1_WT = substr($codon_nonamer1_ORF1, 0, 1);
						my $nt2_nonamer1_WT = substr($codon_nonamer1_ORF1, 1, 1);
						my $nt3_nonamer1_WT = substr($codon_nonamer1_ORF1, 2, 1);
						
						my $nt1_nonamer2_WT = substr($codon_nonamer2_ORF1, 0, 1);
						my $nt2_nonamer2_WT = substr($codon_nonamer2_ORF1, 1, 1);
						my $nt3_nonamer2_WT = substr($codon_nonamer2_ORF1, 2, 1);
						
						foreach my $nt (@nucleotides) {
							my $nt_revcom = revcom($nt);
							
							##################################################################
							# nonamer1
							
							# SITE 1
							# What is each change to ORF1 CODON SITE 1 / ORF2 CODON SITE 1?
							if($nt ne $nt1_nonamer1_WT) { # only one possibility for this site
								 # only one possibility for this site
								
								my $nonamer1_STOP_caused = 0;
								
								my $codon_nonamer1_ORF1_MUT = $codon_nonamer1_ORF1;
								$codon_nonamer1_ORF1_MUT =~ s/[ACGT]([ACGT])([ACGT])/$nt$1$2/;
								
								my $codon_nonamer1_ORF2_prev_MUT = $codon_nonamer1_ORF2_prev;
								$codon_nonamer1_ORF2_prev_MUT =~ s/[ACGT]([ACGT])([ACGT])/$nt_revcom$1$2/;
								
								my $nonamer1_ORF1_effect = 'S';
								
								#nonamer1-ORF1
								my $AA_nonamer1_ORF1_MUT = get_amino_acid($codon_nonamer1_ORF1_MUT);
								
								if($AA_nonamer1_ORF1_MUT ne $AA_nonamer1_ORF1) {
									$nonamer1_ORF1_effect = 'N';
								}
								
								if($AA_nonamer1_ORF1 eq '*' || $AA_nonamer1_ORF1_MUT eq '*') {
									$nonamer1_STOP_caused++;
								}
								
								
								my $nonamer1_ORF2_prev_effect = 'S';
								
								#nonamer1-ORF2
								my $AA_nonamer1_ORF2_prev_MUT = get_amino_acid($codon_nonamer1_ORF2_prev_MUT);
								
								if($AA_nonamer1_ORF2_prev_MUT ne $AA_nonamer1_ORF2_prev) {
									$nonamer1_ORF2_prev_effect = 'N';
								}
								
								if($AA_nonamer1_ORF2_prev eq '*' || $AA_nonamer1_ORF2_prev_MUT eq '*') {
									$nonamer1_STOP_caused++;
								}
								
								# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
								unless($nonamer1_STOP_caused > 0) {
									$nonamer1_num_changes_poss_site1++;
									
									if($nonamer1_ORF1_effect eq 'N') {
										if($nonamer1_ORF2_prev_effect eq 'N') {
											$nonamer1_num_changes_NN_site1++;
											
											# ACTUAL DIFF
											if($nt eq $nt1_nonamer2_WT) {
												$NN_diffs++;
												$site_diffs_hh{$codon_num}->{1}->{NN_diffs}++;
											}
											
										} elsif($nonamer1_ORF2_prev_effect eq 'S') {
											$nonamer1_num_changes_NS_site1++;
											
											# ACTUAL DIFF
											if($nt eq $nt1_nonamer2_WT) {
												$NS_diffs++;
												$site_diffs_hh{$codon_num}->{1}->{NS_diffs}++;
											}
										}
									} elsif($nonamer1_ORF1_effect eq 'S') {
										if($nonamer1_ORF2_prev_effect eq 'N') {
											$nonamer1_num_changes_SN_site1++;
											
											# ACTUAL DIFF
											if($nt eq $nt1_nonamer2_WT) {
												$SN_diffs++;
												$site_diffs_hh{$codon_num}->{1}->{SN_diffs}++;
											}
										} elsif($nonamer1_ORF2_prev_effect eq 'S') {
											$nonamer1_num_changes_SS_site1++;
											
											# ACTUAL DIFF
											if($nt eq $nt1_nonamer2_WT) {
												$SS_diffs++;
												$site_diffs_hh{$codon_num}->{1}->{SS_diffs}++;
											}
										}
									}
								}
							} # end member 1 site 1
							
							
							# SITE 2
							# What is each change to ORF1 CODON SITE 2 / ORF2 CODON SITE 3?
							if($nt ne $nt2_nonamer1_WT) { # only one possibility for this site
								
								my $nonamer1_STOP_caused = 0;
								
								my $codon_nonamer1_ORF1_MUT = $codon_nonamer1_ORF1;
								$codon_nonamer1_ORF1_MUT =~ s/([ACGT])[ACGT]([ACGT])/$1$nt$2/;
								
								my $codon_nonamer1_ORF2_next_MUT = $codon_nonamer1_ORF2_next;
								$codon_nonamer1_ORF2_next_MUT =~ s/([ACGT])([ACGT])[ACGT]/$1$2$nt_revcom/;
								
								my $nonamer1_ORF1_effect = 'S';
								
								#nonamer1-ORF1
								my $AA_nonamer1_ORF1_MUT = get_amino_acid($codon_nonamer1_ORF1_MUT);
								
								if($AA_nonamer1_ORF1_MUT ne $AA_nonamer1_ORF1) {
									$nonamer1_ORF1_effect = 'N';
								}
								
								if($AA_nonamer1_ORF1 eq '*' || $AA_nonamer1_ORF1_MUT eq '*') {
									$nonamer1_STOP_caused++;
								}
								
								
								my $nonamer1_ORF2_next_effect = 'S';
								
								#nonamer1-ORF2
								my $AA_nonamer1_ORF2_next_MUT = get_amino_acid($codon_nonamer1_ORF2_next_MUT);
								
								if($AA_nonamer1_ORF2_next_MUT ne $AA_nonamer1_ORF2_next) {
									$nonamer1_ORF2_next_effect = 'N';
								}
								
								if($AA_nonamer1_ORF2_next eq '*' || $AA_nonamer1_ORF2_next_MUT eq '*') {
									$nonamer1_STOP_caused++;
								}
								
								# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
								unless($nonamer1_STOP_caused > 0) {
									$nonamer1_num_changes_poss_site2++;
									
									if($nonamer1_ORF1_effect eq 'N') {
										if($nonamer1_ORF2_next_effect eq 'N') {
											$nonamer1_num_changes_NN_site2++;
											
											# ACTUAL DIFF
											if($nt eq $nt2_nonamer2_WT) {
												$NN_diffs++;
												$site_diffs_hh{$codon_num}->{2}->{NN_diffs}++;
											}
											
										} elsif($nonamer1_ORF2_next_effect eq 'S') {
											$nonamer1_num_changes_NS_site2++;
											
											# ACTUAL DIFF
											if($nt eq $nt2_nonamer2_WT) {
												$NS_diffs++;
												$site_diffs_hh{$codon_num}->{2}->{NS_diffs}++;
											}
										}
									} elsif($nonamer1_ORF1_effect eq 'S') {
										if($nonamer1_ORF2_next_effect eq 'N') {
											$nonamer1_num_changes_SN_site2++;
											
											# ACTUAL DIFF
											if($nt eq $nt2_nonamer2_WT) {
												$SN_diffs++;
												$site_diffs_hh{$codon_num}->{2}->{SN_diffs}++;
											}
										} elsif($nonamer1_ORF2_next_effect eq 'S') {
											$nonamer1_num_changes_SS_site2++;
											
											# ACTUAL DIFF
											if($nt eq $nt2_nonamer2_WT) {
												$SS_diffs++;
												$site_diffs_hh{$codon_num}->{2}->{SS_diffs}++;
											}
										}
									}
								}
							} # end member 1 site 2
							
							
							# SITE 3
							# What is each change to ORF1 CODON SITE 3 / ORF2 (next) CODON SITE 2?
							if($nt ne $nt3_nonamer1_WT) { # only one possibility for this site
								
								my $nonamer1_STOP_caused = 0;
								
								my $codon_nonamer1_ORF1_MUT = $codon_nonamer1_ORF1;
								$codon_nonamer1_ORF1_MUT =~ s/([ACGT])([ACGT])[ACGT]/$1$2$nt/;
								
								my $codon_nonamer1_ORF2_next_MUT = $codon_nonamer1_ORF2_next;
								$codon_nonamer1_ORF2_next_MUT =~ s/([ACGT])[ACGT]([ACGT])/$1$nt_revcom$2/;
								
								my $nonamer1_ORF1_effect = 'S';
								
								#nonamer1-ORF1
								my $AA_nonamer1_ORF1_MUT = get_amino_acid($codon_nonamer1_ORF1_MUT);
								
								if($AA_nonamer1_ORF1_MUT ne $AA_nonamer1_ORF1) {
									$nonamer1_ORF1_effect = 'N';
								}
								
								if($AA_nonamer1_ORF1 eq '*' || $AA_nonamer1_ORF1_MUT eq '*') {
									$nonamer1_STOP_caused++;
								}
								
								
								my $nonamer1_ORF2_next_effect = 'S';
								
								#nonamer1-ORF2
								my $AA_nonamer1_ORF2_next_MUT = get_amino_acid($codon_nonamer1_ORF2_next_MUT);
								
								if($AA_nonamer1_ORF2_next_MUT ne $AA_nonamer1_ORF2_next) {
									$nonamer1_ORF2_next_effect = 'N';
								}
								
								if($AA_nonamer1_ORF2_next eq '*' || $AA_nonamer1_ORF2_next_MUT eq '*') {
									$nonamer1_STOP_caused++;
								}
								
								# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
								unless($nonamer1_STOP_caused > 0) {
									$nonamer1_num_changes_poss_site3++;
									
									if($nonamer1_ORF1_effect eq 'N') {
										if($nonamer1_ORF2_next_effect eq 'N') {
											$nonamer1_num_changes_NN_site3++;
											
											# ACTUAL DIFF
											if($nt eq $nt3_nonamer2_WT) {
												$NN_diffs++;
												$site_diffs_hh{$codon_num}->{3}->{NN_diffs}++;
											}
											
										} elsif($nonamer1_ORF2_next_effect eq 'S') {
											$nonamer1_num_changes_NS_site3++;
											
											# ACTUAL DIFF
											if($nt eq $nt3_nonamer2_WT) {
												$NS_diffs++;
												$site_diffs_hh{$codon_num}->{3}->{NS_diffs}++;
											}
										}
									} elsif($nonamer1_ORF1_effect eq 'S') {
										if($nonamer1_ORF2_next_effect eq 'N') {
											$nonamer1_num_changes_SN_site3++;
											
											# ACTUAL DIFF
											if($nt eq $nt3_nonamer2_WT) {
												$SN_diffs++;
												$site_diffs_hh{$codon_num}->{3}->{SN_diffs}++;
											}
										} elsif($nonamer1_ORF2_next_effect eq 'S') {
											$nonamer1_num_changes_SS_site3++;
											
											# ACTUAL DIFF
											if($nt eq $nt3_nonamer2_WT) {
												$SS_diffs++;
												$site_diffs_hh{$codon_num}->{3}->{SS_diffs}++;
											}
										}
									}
								}
							} # end member 1 site 3
							
							
							##################################################################
							# nonamer2
							
							# SITE 1
							# What is each change to ORF1 CODON SITE 1 / ORF2 CODON SITE 1?
							if($nt ne $nt1_nonamer2_WT) { # only one possibility for this site
								 # only one possibility for this site
								
								my $nonamer2_STOP_caused = 0;
								
								my $codon_nonamer2_ORF1_MUT = $codon_nonamer2_ORF1;
								$codon_nonamer2_ORF1_MUT =~ s/[ACGT]([ACGT])([ACGT])/$nt$1$2/;
								
								my $codon_nonamer2_ORF2_prev_MUT = $codon_nonamer2_ORF2_prev;
								$codon_nonamer2_ORF2_prev_MUT =~ s/[ACGT]([ACGT])([ACGT])/$nt_revcom$1$2/;
								
								my $nonamer2_ORF1_effect = 'S';
								
								#nonamer1-ORF1
								my $AA_nonamer2_ORF1_MUT = get_amino_acid($codon_nonamer2_ORF1_MUT);
								
								if($AA_nonamer2_ORF1_MUT ne $AA_nonamer2_ORF1) {
									$nonamer2_ORF1_effect = 'N';
								}
								
								if($AA_nonamer2_ORF1 eq '*' || $AA_nonamer2_ORF1_MUT eq '*') {
									$nonamer2_STOP_caused++;
								}
								
								
								my $nonamer2_ORF2_prev_effect = 'S';
								
								#nonamer1-ORF2
								my $AA_nonamer2_ORF2_prev_MUT = get_amino_acid($codon_nonamer2_ORF2_prev_MUT);
								
								if($AA_nonamer2_ORF2_prev_MUT ne $AA_nonamer2_ORF2_prev) {
									$nonamer2_ORF2_prev_effect = 'N';
								}
								
								if($AA_nonamer2_ORF2_prev eq '*' || $AA_nonamer2_ORF2_prev_MUT eq '*') {
									$nonamer2_STOP_caused++;
								}
								
								# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
								unless($nonamer2_STOP_caused > 0) {
									$nonamer2_num_changes_poss_site1++;
									
									if($nonamer2_ORF1_effect eq 'N') {
										if($nonamer2_ORF2_prev_effect eq 'N') {
											$nonamer2_num_changes_NN_site1++;
											
											
										} elsif($nonamer2_ORF2_prev_effect eq 'S') {
											$nonamer2_num_changes_NS_site1++;
											
										}
									} elsif($nonamer2_ORF1_effect eq 'S') {
										if($nonamer2_ORF2_prev_effect eq 'N') {
											$nonamer2_num_changes_SN_site1++;
											
										} elsif($nonamer2_ORF2_prev_effect eq 'S') {
											$nonamer2_num_changes_SS_site1++;
											
										}
									}
								}
							} # end member 2 site 1
							
							
							# SITE 2
							# What is each change to ORF1 CODON SITE 2 / ORF2 [next] CODON SITE 3?
							if($nt ne $nt2_nonamer2_WT) { # only one possibility for this site
								
								my $nonamer2_STOP_caused = 0;
								
								my $codon_nonamer2_ORF1_MUT = $codon_nonamer2_ORF1;
								$codon_nonamer2_ORF1_MUT =~ s/([ACGT])[ACGT]([ACGT])/$1$nt$2/;
								
								my $codon_nonamer2_ORF2_next_MUT = $codon_nonamer2_ORF2_next;
								$codon_nonamer2_ORF2_next_MUT =~ s/([ACGT])([ACGT])[ACGT]/$1$2$nt_revcom/;
								
								my $nonamer2_ORF1_effect = 'S';
								
								#nonamer1-ORF1
								my $AA_nonamer2_ORF1_MUT = get_amino_acid($codon_nonamer2_ORF1_MUT);
								
								if($AA_nonamer2_ORF1_MUT ne $AA_nonamer2_ORF1) {
									$nonamer2_ORF1_effect = 'N';
								}
								
								if($AA_nonamer2_ORF1 eq '*' || $AA_nonamer2_ORF1_MUT eq '*') {
									$nonamer2_STOP_caused++;
								}
								
								
								my $nonamer2_ORF2_next_effect = 'S';
								
								#nonamer1-ORF2
								my $AA_nonamer2_ORF2_next_MUT = get_amino_acid($codon_nonamer2_ORF2_next_MUT);
								
								if($AA_nonamer2_ORF2_next_MUT ne $AA_nonamer2_ORF2_next) {
									$nonamer2_ORF2_next_effect = 'N';
								}
								
								if($AA_nonamer2_ORF2_next eq '*' || $AA_nonamer2_ORF2_next_MUT eq '*') {
									$nonamer2_STOP_caused++;
								}
								
								# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
								unless($nonamer2_STOP_caused > 0) {
									$nonamer2_num_changes_poss_site2++;
									
									if($nonamer2_ORF1_effect eq 'N') {
										if($nonamer2_ORF2_next_effect eq 'N') {
											$nonamer2_num_changes_NN_site2++;
											
										} elsif($nonamer2_ORF2_next_effect eq 'S') {
											$nonamer2_num_changes_NS_site2++;
											
										}
									} elsif($nonamer2_ORF1_effect eq 'S') {
										if($nonamer2_ORF2_next_effect eq 'N') {
											$nonamer2_num_changes_SN_site2++;
											
										} elsif($nonamer2_ORF2_next_effect eq 'S') {
											$nonamer2_num_changes_SS_site2++;
											
										}
									}
								}
							} # end member 2 site 2
							
							
							# SITE 3
							# What is each change to ORF1 CODON SITE 3 / ORF2 (next) CODON SITE 2?
							if($nt ne $nt3_nonamer2_WT) { # only one possibility for this site
								
								my $nonamer2_STOP_caused = 0;
								
								my $codon_nonamer2_ORF1_MUT = $codon_nonamer2_ORF1;
								$codon_nonamer2_ORF1_MUT =~ s/([ACGT])([ACGT])[ACGT]/$1$2$nt/;
								
								my $codon_nonamer2_ORF2_next_MUT = $codon_nonamer2_ORF2_next;
								$codon_nonamer2_ORF2_next_MUT =~ s/([ACGT])[ACGT]([ACGT])/$1$nt_revcom$2/;
								
								my $nonamer2_ORF1_effect = 'S';
								
								#nonamer2-ORF1
								my $AA_nonamer2_ORF1_MUT = get_amino_acid($codon_nonamer2_ORF1_MUT);
								
								if($AA_nonamer2_ORF1_MUT ne $AA_nonamer2_ORF1) {
									$nonamer2_ORF1_effect = 'N';
								}
								
								if($AA_nonamer2_ORF1 eq '*' || $AA_nonamer2_ORF1_MUT eq '*') {
									$nonamer2_STOP_caused++;
								}
								
								
								my $nonamer2_ORF2_next_effect = 'S';
								
								#nonamer2-ORF2
								my $AA_nonamer2_ORF2_next_MUT = get_amino_acid($codon_nonamer2_ORF2_next_MUT);
								
								if($AA_nonamer2_ORF2_next_MUT ne $AA_nonamer2_ORF2_next) {
									$nonamer2_ORF2_next_effect = 'N';
								}
								
								if($AA_nonamer2_ORF2_next eq '*' || $AA_nonamer2_ORF2_next_MUT eq '*') {
									$nonamer2_STOP_caused++;
								}
								
								# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
								unless($nonamer2_STOP_caused > 0) {
									$nonamer2_num_changes_poss_site3++;
									
									if($nonamer2_ORF1_effect eq 'N') {
										if($nonamer2_ORF2_next_effect eq 'N') {
											$nonamer2_num_changes_NN_site3++;
											
										} elsif($nonamer2_ORF2_next_effect eq 'S') {
											$nonamer2_num_changes_NS_site3++;
											
										}
									} elsif($nonamer2_ORF1_effect eq 'S') {
										if($nonamer2_ORF2_next_effect eq 'N') {
											$nonamer2_num_changes_SN_site3++;
											
										} elsif($nonamer2_ORF2_next_effect eq 'S') {
											$nonamer2_num_changes_SS_site3++;
											
										}
									}
								}
							} # end member 2 site 3
						} # end cycling all 4 nucleotides
						
						# TALLY SITE 1
						my $NN_sites_nonamer1_site1 = 'NA';
						my $SN_sites_nonamer1_site1 = 'NA';
						my $NS_sites_nonamer1_site1 = 'NA';
						my $SS_sites_nonamer1_site1 = 'NA';
						
						if($nonamer1_num_changes_poss_site1 > 0) {
							$NN_sites_nonamer1_site1 = $nonamer1_num_changes_NN_site1 / $nonamer1_num_changes_poss_site1;
							$SN_sites_nonamer1_site1 = $nonamer1_num_changes_SN_site1 / $nonamer1_num_changes_poss_site1;
							$NS_sites_nonamer1_site1 = $nonamer1_num_changes_NS_site1 / $nonamer1_num_changes_poss_site1;
							$SS_sites_nonamer1_site1 = $nonamer1_num_changes_SS_site1 / $nonamer1_num_changes_poss_site1;
							
							unless(exists $seq_completed{$nonamer1}) {
								$seq2sites{$nonamer1}->{NN_sites} += $NN_sites_nonamer1_site1;
								$seq2sites{$nonamer1}->{SN_sites} += $SN_sites_nonamer1_site1;
								$seq2sites{$nonamer1}->{NS_sites} += $NS_sites_nonamer1_site1;
								$seq2sites{$nonamer1}->{SS_sites} += $SS_sites_nonamer1_site1;
							}
						}
						
						my $NN_sites_nonamer2_site1 = 'NA';
						my $SN_sites_nonamer2_site1 = 'NA';
						my $NS_sites_nonamer2_site1 = 'NA';
						my $SS_sites_nonamer2_site1 = 'NA';
						
						if($nonamer2_num_changes_poss_site1 > 0) {
							$NN_sites_nonamer2_site1 = $nonamer2_num_changes_NN_site1 / $nonamer2_num_changes_poss_site1;
							$SN_sites_nonamer2_site1 = $nonamer2_num_changes_SN_site1 / $nonamer2_num_changes_poss_site1;
							$NS_sites_nonamer2_site1 = $nonamer2_num_changes_NS_site1 / $nonamer2_num_changes_poss_site1;
							$SS_sites_nonamer2_site1 = $nonamer2_num_changes_SS_site1 / $nonamer2_num_changes_poss_site1;
							
							unless(exists $seq_completed{$nonamer2}) {
								$seq2sites{$nonamer2}->{NN_sites} += $NN_sites_nonamer2_site1;
								$seq2sites{$nonamer2}->{SN_sites} += $SN_sites_nonamer2_site1;
								$seq2sites{$nonamer2}->{NS_sites} += $NS_sites_nonamer2_site1;
								$seq2sites{$nonamer2}->{SS_sites} += $SS_sites_nonamer2_site1;
							}
						}
						
						my $NN_sites_site1 = 'NA';
						my $SN_sites_site1 = 'NA';
						my $NS_sites_site1 = 'NA';
						my $SS_sites_site1 = 'NA';
						
						if($nonamer1_num_changes_poss_site1 > 0 && $nonamer2_num_changes_poss_site1 > 0) {
							$NN_sites_site1 = ($NN_sites_nonamer1_site1 + $NN_sites_nonamer2_site1) / 2;
							$SN_sites_site1 = ($SN_sites_nonamer1_site1 + $SN_sites_nonamer2_site1) / 2;
							$NS_sites_site1 = ($NS_sites_nonamer1_site1 + $NS_sites_nonamer2_site1) / 2;
							$SS_sites_site1 = ($SS_sites_nonamer1_site1 + $SS_sites_nonamer2_site1) / 2;
						} elsif($nonamer1_num_changes_poss_site1 > 0) {
							$NN_sites_site1 = $NN_sites_nonamer1_site1;
							$SN_sites_site1 = $SN_sites_nonamer1_site1;
							$NS_sites_site1 = $NS_sites_nonamer1_site1;
							$SS_sites_site1 = $SS_sites_nonamer1_site1;
						} elsif($nonamer2_num_changes_poss_site1 > 0) {
							$NN_sites_site1 = $NN_sites_nonamer2_site1;
							$SN_sites_site1 = $SN_sites_nonamer2_site1;
							$NS_sites_site1 = $NS_sites_nonamer2_site1;
							$SS_sites_site1 = $SS_sites_nonamer2_site1;
							
						} # else nothing defined, stay NA
						
						# TALLY SITE 2
						my $NN_sites_nonamer1_site2 = 'NA';
						my $SN_sites_nonamer1_site2 = 'NA';
						my $NS_sites_nonamer1_site2 = 'NA';
						my $SS_sites_nonamer1_site2 = 'NA';
						
						if($nonamer1_num_changes_poss_site2 > 0) {
							$NN_sites_nonamer1_site2 = $nonamer1_num_changes_NN_site2 / $nonamer1_num_changes_poss_site2;
							$SN_sites_nonamer1_site2 = $nonamer1_num_changes_SN_site2 / $nonamer1_num_changes_poss_site2;
							$NS_sites_nonamer1_site2 = $nonamer1_num_changes_NS_site2 / $nonamer1_num_changes_poss_site2;
							$SS_sites_nonamer1_site2 = $nonamer1_num_changes_SS_site2 / $nonamer1_num_changes_poss_site2;
							
							unless(exists $seq_completed{$nonamer1}) {
								$seq2sites{$nonamer1}->{NN_sites} += $NN_sites_nonamer1_site2;
								$seq2sites{$nonamer1}->{SN_sites} += $SN_sites_nonamer1_site2;
								$seq2sites{$nonamer1}->{NS_sites} += $NS_sites_nonamer1_site2;
								$seq2sites{$nonamer1}->{SS_sites} += $SS_sites_nonamer1_site2;
							}
						}
						
						my $NN_sites_nonamer2_site2 = 'NA';
						my $SN_sites_nonamer2_site2 = 'NA';
						my $NS_sites_nonamer2_site2 = 'NA';
						my $SS_sites_nonamer2_site2 = 'NA';
						
						if($nonamer2_num_changes_poss_site2 > 0) {
							$NN_sites_nonamer2_site2 = $nonamer2_num_changes_NN_site2 / $nonamer2_num_changes_poss_site2;
							$SN_sites_nonamer2_site2 = $nonamer2_num_changes_SN_site2 / $nonamer2_num_changes_poss_site2;
							$NS_sites_nonamer2_site2 = $nonamer2_num_changes_NS_site2 / $nonamer2_num_changes_poss_site2;
							$SS_sites_nonamer2_site2 = $nonamer2_num_changes_SS_site2 / $nonamer2_num_changes_poss_site2;
							
							unless(exists $seq_completed{$nonamer2}) {
								$seq2sites{$nonamer2}->{NN_sites} += $NN_sites_nonamer2_site2;
								$seq2sites{$nonamer2}->{SN_sites} += $SN_sites_nonamer2_site2;
								$seq2sites{$nonamer2}->{NS_sites} += $NS_sites_nonamer2_site2;
								$seq2sites{$nonamer2}->{SS_sites} += $SS_sites_nonamer2_site2;
							}
						}
						
						my $NN_sites_site2 = 'NA';
						my $SN_sites_site2 = 'NA';
						my $NS_sites_site2 = 'NA';
						my $SS_sites_site2 = 'NA';
						
						if($nonamer1_num_changes_poss_site2 > 0 && $nonamer2_num_changes_poss_site2 > 0) {
							$NN_sites_site2 = ($NN_sites_nonamer1_site2 + $NN_sites_nonamer2_site2) / 2;
							$SN_sites_site2 = ($SN_sites_nonamer1_site2 + $SN_sites_nonamer2_site2) / 2;
							$NS_sites_site2 = ($NS_sites_nonamer1_site2 + $NS_sites_nonamer2_site2) / 2;
							$SS_sites_site2 = ($SS_sites_nonamer1_site2 + $SS_sites_nonamer2_site2) / 2;
						} elsif($nonamer1_num_changes_poss_site2 > 0) {
							$NN_sites_site2 = $NN_sites_nonamer1_site2;
							$SN_sites_site2 = $SN_sites_nonamer1_site2;
							$NS_sites_site2 = $NS_sites_nonamer1_site2;
							$SS_sites_site2 = $SS_sites_nonamer1_site2;
						} elsif($nonamer2_num_changes_poss_site2 > 0) {
							$NN_sites_site2 = $NN_sites_nonamer2_site2;
							$SN_sites_site2 = $SN_sites_nonamer2_site2;
							$NS_sites_site2 = $NS_sites_nonamer2_site2;
							$SS_sites_site2 = $SS_sites_nonamer2_site2;
							
						} # else nothing defined, stay NA
						
						# TALLY SITE 3
						my $NN_sites_nonamer1_site3 = 'NA';
						my $SN_sites_nonamer1_site3 = 'NA';
						my $NS_sites_nonamer1_site3 = 'NA';
						my $SS_sites_nonamer1_site3 = 'NA';
						
						if($nonamer1_num_changes_poss_site3 > 0) {
							$NN_sites_nonamer1_site3 = $nonamer1_num_changes_NN_site3 / $nonamer1_num_changes_poss_site3;
							$SN_sites_nonamer1_site3 = $nonamer1_num_changes_SN_site3 / $nonamer1_num_changes_poss_site3;
							$NS_sites_nonamer1_site3 = $nonamer1_num_changes_NS_site3 / $nonamer1_num_changes_poss_site3;
							$SS_sites_nonamer1_site3 = $nonamer1_num_changes_SS_site3 / $nonamer1_num_changes_poss_site3;
							
							unless(exists $seq_completed{$nonamer1}) {
								$seq2sites{$nonamer1}->{NN_sites} += $NN_sites_nonamer1_site3;
								$seq2sites{$nonamer1}->{SN_sites} += $SN_sites_nonamer1_site3;
								$seq2sites{$nonamer1}->{NS_sites} += $NS_sites_nonamer1_site3;
								$seq2sites{$nonamer1}->{SS_sites} += $SS_sites_nonamer1_site3;
							}
						}
						
						my $NN_sites_nonamer2_site3 = 'NA';
						my $SN_sites_nonamer2_site3 = 'NA';
						my $NS_sites_nonamer2_site3 = 'NA';
						my $SS_sites_nonamer2_site3 = 'NA';
						
						if($nonamer2_num_changes_poss_site3 > 0) {
							$NN_sites_nonamer2_site3 = $nonamer2_num_changes_NN_site3 / $nonamer2_num_changes_poss_site3;
							$SN_sites_nonamer2_site3 = $nonamer2_num_changes_SN_site3 / $nonamer2_num_changes_poss_site3;
							$NS_sites_nonamer2_site3 = $nonamer2_num_changes_NS_site3 / $nonamer2_num_changes_poss_site3;
							$SS_sites_nonamer2_site3 = $nonamer2_num_changes_SS_site3 / $nonamer2_num_changes_poss_site3;
							
							unless(exists $seq_completed{$nonamer2}) {
								$seq2sites{$nonamer2}->{NN_sites} += $NN_sites_nonamer2_site3;
								$seq2sites{$nonamer2}->{SN_sites} += $SN_sites_nonamer2_site3;
								$seq2sites{$nonamer2}->{NS_sites} += $NS_sites_nonamer2_site3;
								$seq2sites{$nonamer2}->{SS_sites} += $SS_sites_nonamer2_site3;
							}
						}
						
						my $NN_sites_site3 = 'NA';
						my $SN_sites_site3 = 'NA';
						my $NS_sites_site3 = 'NA';
						my $SS_sites_site3 = 'NA';
						
						if($nonamer1_num_changes_poss_site3 > 0 && $nonamer2_num_changes_poss_site3 > 0) {
							$NN_sites_site3 = ($NN_sites_nonamer1_site3 + $NN_sites_nonamer2_site3) / 2;
							$SN_sites_site3 = ($SN_sites_nonamer1_site3 + $SN_sites_nonamer2_site3) / 2;
							$NS_sites_site3 = ($NS_sites_nonamer1_site3 + $NS_sites_nonamer2_site3) / 2;
							$SS_sites_site3 = ($SS_sites_nonamer1_site3 + $SS_sites_nonamer2_site3) / 2;
						} elsif($nonamer1_num_changes_poss_site3 > 0) {
							$NN_sites_site3 = $NN_sites_nonamer1_site3;
							$SN_sites_site3 = $SN_sites_nonamer1_site3;
							$NS_sites_site3 = $NS_sites_nonamer1_site3;
							$SS_sites_site3 = $SS_sites_nonamer1_site3;
						} elsif($nonamer2_num_changes_poss_site3 > 0) {
							$NN_sites_site3 = $NN_sites_nonamer2_site3;
							$SN_sites_site3 = $SN_sites_nonamer2_site3;
							$NS_sites_site3 = $NS_sites_nonamer2_site3;
							$SS_sites_site3 = $SS_sites_nonamer2_site3;
							
						} # else nothing defined, stay NA
						
						# SUM THE THREE SITES
						my $NN_sites = $NN_sites_site1 + $NN_sites_site2 + $NN_sites_site3;
						my $SN_sites = $SN_sites_site1 + $SN_sites_site2 + $SN_sites_site3;
						my $NS_sites = $NS_sites_site1 + $NS_sites_site2 + $NS_sites_site3;
						my $SS_sites = $SS_sites_site1 + $SS_sites_site2 + $SS_sites_site3;
						
						
						# Check if there are multiple variants in these overlapping codons
						my $MNV = 'FALSE';
						if(($NN_diffs + $SN_diffs + $NS_diffs + $SS_diffs) > 1) {
							$MNV = 'TRUE';
							$this_codon_MNV = 'TRUE'; # only needs to be TRUE in one comparison
						}
						
						# NEW NAIVE APPROACH: ADD TO SUMS
						my $comparison_weight = 0; # number of pairwise comparisons involving these nonamers
						if($nonamer1_index == $nonamer2_index) { # if($nonamer1 eq $nonamer2) { # it's a self-comparison
							# Combination
							$comparison_weight = (($nonamer1_count * $nonamer1_count) - $nonamer1_count) / 2; # nonamer 1 and 2 are the same
						} else {
							# m * n
							my $nonamer2_count = $unique_9mers{$seq_site_index}->{$nonamer2};
							$comparison_weight = $nonamer1_count * $nonamer2_count;
						}
						
						$comparisons_sum += $comparison_weight;
						
						$NN_sites_numerator += $NN_sites * $comparison_weight;
						$SN_sites_numerator += $SN_sites * $comparison_weight;
						$NS_sites_numerator += $NS_sites * $comparison_weight;
						$SS_sites_numerator += $SS_sites * $comparison_weight;
						
						$NN_diffs_numerator += $NN_diffs * $comparison_weight;
						$SN_diffs_numerator += $SN_diffs * $comparison_weight;
						$NS_diffs_numerator += $NS_diffs * $comparison_weight;
						$SS_diffs_numerator += $SS_diffs * $comparison_weight;
						
						
					} # finish MIDDLE (internal) codon (not first or last; two ORF2 codons overlap)
						
				
##################################################################################################
##################### SENSE-ANTISENSE:
#####################  sas12:
#####################    ORF1: 1-2-3-1-2-3-1
#####################    ORF2: 2-1-3-2-1-3-2
##################################################################################################
				} elsif($frame eq 'sas12') {
					
					#print "site_index=$site_index\n";
					#print "nonamer1=$nonamer1\n";
					#print "nonamer2=$nonamer2\n";
					
					# First 2 nt of prev codon, last 1 nt of next codon (opposite strand)
					my $codon_nonamer1_ORF2_prev = revcom(substr($nonamer1, ($site_index - 1), 3));
					my $codon_nonamer1_ORF2_next = revcom(substr($nonamer1, ($site_index + 2), 3));
					my $codon_nonamer2_ORF2_prev = revcom(substr($nonamer2, ($site_index - 1), 3));
					my $codon_nonamer2_ORF2_next = revcom(substr($nonamer2, ($site_index + 2), 3));
					
					if($verbose_messages) {
						print "site_index=$site_index\n";
						print "nonamer1=$nonamer1\n";
						print "nonamer2=$nonamer2\n";
						print "nonamer1_count=$nonamer1_count\n";
						#print "nonamer2=$nonamer2_count\n";s
						print "codon_nonamer1_ORF2_prev=$codon_nonamer1_ORF2_prev\n";
						print "codon_nonamer1_ORF2_next=$codon_nonamer1_ORF2_next\n";
						print "codon_nonamer2_ORF2_prev=$codon_nonamer2_ORF2_prev\n";
						print "codon_nonamer2_ORF2_next=$codon_nonamer2_ORF2_next\n";
					}
					
					my $AA_nonamer1_ORF2_prev = get_amino_acid($codon_nonamer1_ORF2_prev);
					my $AA_nonamer1_ORF2_next = get_amino_acid($codon_nonamer1_ORF2_next);
					my $AA_nonamer2_ORF2_prev = get_amino_acid($codon_nonamer2_ORF2_prev);
					my $AA_nonamer2_ORF2_next = get_amino_acid($codon_nonamer2_ORF2_next);
					
					if($codon_num > 2 && ($AA_nonamer1_ORF2_prev eq '*' || $AA_nonamer1_ORF2_next eq '*' || $AA_nonamer2_ORF2_prev eq '*' || $AA_nonamer2_ORF2_next eq '*')) {
						print "### WARNING! ORF2, $nonamer1\-$nonamer2 comparison, near ORF1 codon $codon_in_seq encodes a within-frame STOP codon. Wrong frame selection ($frame)?\n";
					}
					
					#######################
					# GET NUMBER OF SITES
					
					if($codon_num > 1 && $codon_num < $num_codons) {
						# both ORF2's PREV and NEXT codons fully defined, 3 sites to examine
						
						# New phylogeny-naive approach
						$alt1_codons_to_counts{$codon_nonamer1_ORF2_prev} += $nonamer1_count;
						$alt1_codons_to_counts{$codon_nonamer2_ORF2_prev} += $nonamer1_count; # ADDED
						$alt2_codons_to_counts{$codon_nonamer1_ORF2_next} += $nonamer1_count;
						$alt2_codons_to_counts{$codon_nonamer2_ORF2_next} += $nonamer1_count; # ADDED
						
						my $nonamer1_num_changes_poss_site1 = 0;
						my $nonamer1_num_changes_NN_site1 = 0;
						my $nonamer1_num_changes_SN_site1 = 0;
						my $nonamer1_num_changes_NS_site1 = 0;
						my $nonamer1_num_changes_SS_site1 = 0;
						
						my $nonamer1_num_changes_poss_site2 = 0;
						my $nonamer1_num_changes_NN_site2 = 0;
						my $nonamer1_num_changes_SN_site2 = 0;
						my $nonamer1_num_changes_NS_site2 = 0;
						my $nonamer1_num_changes_SS_site2 = 0;
						
						my $nonamer1_num_changes_poss_site3 = 0;
						my $nonamer1_num_changes_NN_site3 = 0;
						my $nonamer1_num_changes_SN_site3 = 0;
						my $nonamer1_num_changes_NS_site3 = 0;
						my $nonamer1_num_changes_SS_site3 = 0;
						
						my $nonamer2_num_changes_poss_site1 = 0;
						my $nonamer2_num_changes_NN_site1 = 0;
						my $nonamer2_num_changes_SN_site1 = 0;
						my $nonamer2_num_changes_NS_site1 = 0;
						my $nonamer2_num_changes_SS_site1 = 0;
						
						my $nonamer2_num_changes_poss_site2 = 0;
						my $nonamer2_num_changes_NN_site2 = 0;
						my $nonamer2_num_changes_SN_site2 = 0;
						my $nonamer2_num_changes_NS_site2 = 0;
						my $nonamer2_num_changes_SS_site2 = 0;
						
						my $nonamer2_num_changes_poss_site3 = 0;
						my $nonamer2_num_changes_NN_site3 = 0;
						my $nonamer2_num_changes_SN_site3 = 0;
						my $nonamer2_num_changes_NS_site3 = 0;
						my $nonamer2_num_changes_SS_site3 = 0;
						
						my $NN_diffs = 0;
						my $SN_diffs = 0;
						my $NS_diffs = 0;
						my $SS_diffs = 0;
						
						# Just ORF1 (reference)
						my $nt1_nonamer1_WT = substr($codon_nonamer1_ORF1, 0, 1);
						my $nt2_nonamer1_WT = substr($codon_nonamer1_ORF1, 1, 1);
						my $nt3_nonamer1_WT = substr($codon_nonamer1_ORF1, 2, 1);
						
						my $nt1_nonamer2_WT = substr($codon_nonamer2_ORF1, 0, 1);
						my $nt2_nonamer2_WT = substr($codon_nonamer2_ORF1, 1, 1);
						my $nt3_nonamer2_WT = substr($codon_nonamer2_ORF1, 2, 1);
						
						foreach my $nt (@nucleotides) {
							my $nt_revcom = revcom($nt);
							
							##################################################################
							# nonamer1
							
							# SITE 1
							# What is each change to ORF1 CODON SITE 1 / ORF2 CODON SITE 2?
							if($nt ne $nt1_nonamer1_WT) { # only one possibility for this site
								 # only one possibility for this site
								
								my $nonamer1_STOP_caused = 0;
								
								my $codon_nonamer1_ORF1_MUT = $codon_nonamer1_ORF1;
								$codon_nonamer1_ORF1_MUT =~ s/[ACGT]([ACGT])([ACGT])/$nt$1$2/;
								
								my $codon_nonamer1_ORF2_prev_MUT = $codon_nonamer1_ORF2_prev;
								$codon_nonamer1_ORF2_prev_MUT =~ s/([ACGT])[ACGT]([ACGT])/$1$nt_revcom$2/;
								
								my $nonamer1_ORF1_effect = 'S';
								
								#nonamer1-ORF1
								my $AA_nonamer1_ORF1_MUT = get_amino_acid($codon_nonamer1_ORF1_MUT);
								
								if($AA_nonamer1_ORF1_MUT ne $AA_nonamer1_ORF1) {
									$nonamer1_ORF1_effect = 'N';
								}
								
								if($AA_nonamer1_ORF1 eq '*' || $AA_nonamer1_ORF1_MUT eq '*') {
									$nonamer1_STOP_caused++;
								}
								
								my $nonamer1_ORF2_prev_effect = 'S';
								
								#nonamer1-ORF2
								my $AA_nonamer1_ORF2_prev_MUT = get_amino_acid($codon_nonamer1_ORF2_prev_MUT);
								
								if($AA_nonamer1_ORF2_prev_MUT ne $AA_nonamer1_ORF2_prev) {
									$nonamer1_ORF2_prev_effect = 'N';
								}
								
								if($AA_nonamer1_ORF2_prev eq '*' || $AA_nonamer1_ORF2_prev_MUT eq '*') {
									$nonamer1_STOP_caused++;
								}
								
								# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
								unless($nonamer1_STOP_caused > 0) {
									$nonamer1_num_changes_poss_site1++;
									
									if($nonamer1_ORF1_effect eq 'N') {
										if($nonamer1_ORF2_prev_effect eq 'N') {
											$nonamer1_num_changes_NN_site1++;
											
											# ACTUAL DIFF
											if($nt eq $nt1_nonamer2_WT) {
												$NN_diffs++;
												$site_diffs_hh{$codon_num}->{1}->{NN_diffs}++;
											}
											
										} elsif($nonamer1_ORF2_prev_effect eq 'S') {
											$nonamer1_num_changes_NS_site1++;
											
											# ACTUAL DIFF
											if($nt eq $nt1_nonamer2_WT) {
												$NS_diffs++;
												$site_diffs_hh{$codon_num}->{1}->{NS_diffs}++;
											}
										}
									} elsif($nonamer1_ORF1_effect eq 'S') {
										if($nonamer1_ORF2_prev_effect eq 'N') {
											$nonamer1_num_changes_SN_site1++;
											
											# ACTUAL DIFF
											if($nt eq $nt1_nonamer2_WT) {
												$SN_diffs++;
												$site_diffs_hh{$codon_num}->{1}->{SN_diffs}++;
											}
										} elsif($nonamer1_ORF2_prev_effect eq 'S') {
											$nonamer1_num_changes_SS_site1++;
											
											# ACTUAL DIFF
											if($nt eq $nt1_nonamer2_WT) {
												$SS_diffs++;
												$site_diffs_hh{$codon_num}->{1}->{SS_diffs}++;
											}
										}
									}
								}
							} # end member 1 site 1
							
							
							# SITE 2
							# What is each change to ORF1 CODON SITE 2 / ORF2 CODON SITE 3?
							if($nt ne $nt2_nonamer1_WT) { # only one possibility for this site
								
								my $nonamer1_STOP_caused = 0;
								
								my $codon_nonamer1_ORF1_MUT = $codon_nonamer1_ORF1;
								$codon_nonamer1_ORF1_MUT =~ s/([ACGT])[ACGT]([ACGT])/$1$nt$2/;
								
								my $codon_nonamer1_ORF2_prev_MUT = $codon_nonamer1_ORF2_prev;
								$codon_nonamer1_ORF2_prev_MUT =~ s/[ACGT]([ACGT])([ACGT])/$nt_revcom$1$2/;
								
								my $nonamer1_ORF1_effect = 'S';
								
								#nonamer1-ORF1
								my $AA_nonamer1_ORF1_MUT = get_amino_acid($codon_nonamer1_ORF1_MUT);
								
								if($AA_nonamer1_ORF1_MUT ne $AA_nonamer1_ORF1) {
									$nonamer1_ORF1_effect = 'N';
								}
								
								if($AA_nonamer1_ORF1 eq '*' || $AA_nonamer1_ORF1_MUT eq '*') {
									$nonamer1_STOP_caused++;
								}
								
								
								my $nonamer1_ORF2_prev_effect = 'S';
								
								#nonamer1-ORF2
								my $AA_nonamer1_ORF2_prev_MUT = get_amino_acid($codon_nonamer1_ORF2_prev_MUT);
								
								if($AA_nonamer1_ORF2_prev_MUT ne $AA_nonamer1_ORF2_prev) {
									$nonamer1_ORF2_prev_effect = 'N';
								}
								
								if($AA_nonamer1_ORF2_prev eq '*' || $AA_nonamer1_ORF2_prev_MUT eq '*') {
									$nonamer1_STOP_caused++;
								}
								
								# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
								unless($nonamer1_STOP_caused > 0) {
									$nonamer1_num_changes_poss_site2++;
									
									if($nonamer1_ORF1_effect eq 'N') {
										if($nonamer1_ORF2_prev_effect eq 'N') {
											$nonamer1_num_changes_NN_site2++;
											
											# ACTUAL DIFF
											if($nt eq $nt2_nonamer2_WT) {
												$NN_diffs++;
												$site_diffs_hh{$codon_num}->{2}->{NN_diffs}++;
											}
											
										} elsif($nonamer1_ORF2_prev_effect eq 'S') {
											$nonamer1_num_changes_NS_site2++;
											
											# ACTUAL DIFF
											if($nt eq $nt2_nonamer2_WT) {
												$NS_diffs++;
												$site_diffs_hh{$codon_num}->{2}->{NS_diffs}++;
											}
										}
									} elsif($nonamer1_ORF1_effect eq 'S') {
										if($nonamer1_ORF2_prev_effect eq 'N') {
											$nonamer1_num_changes_SN_site2++;
											
											# ACTUAL DIFF
											if($nt eq $nt2_nonamer2_WT) {
												$SN_diffs++;
												$site_diffs_hh{$codon_num}->{2}->{SN_diffs}++;
											}
										} elsif($nonamer1_ORF2_prev_effect eq 'S') {
											$nonamer1_num_changes_SS_site2++;
											
											# ACTUAL DIFF
											if($nt eq $nt2_nonamer2_WT) {
												$SS_diffs++;
												$site_diffs_hh{$codon_num}->{2}->{SS_diffs}++;
											}
										}
									}
								}
							} # end member 1 site 2
							
							
							# SITE 3
							# What is each change to ORF1 CODON SITE 3 / ORF2 (next) CODON SITE 1?
							if($nt ne $nt3_nonamer1_WT) { # only one possibility for this site
								
								my $nonamer1_STOP_caused = 0;
								
								my $codon_nonamer1_ORF1_MUT = $codon_nonamer1_ORF1;
								$codon_nonamer1_ORF1_MUT =~ s/([ACGT])([ACGT])[ACGT]/$1$2$nt/;
								
								my $codon_nonamer1_ORF2_next_MUT = $codon_nonamer1_ORF2_next;
								$codon_nonamer1_ORF2_next_MUT =~ s/([ACGT])([ACGT])[ACGT]/$1$2$nt_revcom/;
								
								my $nonamer1_ORF1_effect = 'S';
								
								#nonamer1-ORF1
								my $AA_nonamer1_ORF1_MUT = get_amino_acid($codon_nonamer1_ORF1_MUT);
								
								if($AA_nonamer1_ORF1_MUT ne $AA_nonamer1_ORF1) {
									$nonamer1_ORF1_effect = 'N';
								}
								
								if($AA_nonamer1_ORF1 eq '*' || $AA_nonamer1_ORF1_MUT eq '*') {
									$nonamer1_STOP_caused++;
								}
								
								
								my $nonamer1_ORF2_next_effect = 'S';
								
								#nonamer1-ORF2
								my $AA_nonamer1_ORF2_next_MUT = get_amino_acid($codon_nonamer1_ORF2_next_MUT);
								
								if($AA_nonamer1_ORF2_next_MUT ne $AA_nonamer1_ORF2_next) {
									$nonamer1_ORF2_next_effect = 'N';
								}
								
								if($AA_nonamer1_ORF2_next eq '*' || $AA_nonamer1_ORF2_next_MUT eq '*') {
									$nonamer1_STOP_caused++;
								}
								
								# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
								unless($nonamer1_STOP_caused > 0) {
									$nonamer1_num_changes_poss_site3++;
									
									if($nonamer1_ORF1_effect eq 'N') {
										if($nonamer1_ORF2_next_effect eq 'N') {
											$nonamer1_num_changes_NN_site3++;
											
											# ACTUAL DIFF
											if($nt eq $nt3_nonamer2_WT) {
												$NN_diffs++;
												$site_diffs_hh{$codon_num}->{3}->{NN_diffs}++;
											}
											
										} elsif($nonamer1_ORF2_next_effect eq 'S') {
											$nonamer1_num_changes_NS_site3++;
											
											# ACTUAL DIFF
											if($nt eq $nt3_nonamer2_WT) {
												$NS_diffs++;
												$site_diffs_hh{$codon_num}->{3}->{NS_diffs}++;
											}
										}
									} elsif($nonamer1_ORF1_effect eq 'S') {
										if($nonamer1_ORF2_next_effect eq 'N') {
											$nonamer1_num_changes_SN_site3++;
											
											# ACTUAL DIFF
											if($nt eq $nt3_nonamer2_WT) {
												$SN_diffs++;
												$site_diffs_hh{$codon_num}->{3}->{SN_diffs}++;
											}
										} elsif($nonamer1_ORF2_next_effect eq 'S') {
											$nonamer1_num_changes_SS_site3++;
											
											# ACTUAL DIFF
											if($nt eq $nt3_nonamer2_WT) {
												$SS_diffs++;
												$site_diffs_hh{$codon_num}->{3}->{SS_diffs}++;
											}
										}
									}
								}
							} # end member 1 site 3
							
							
							
							##################################################################
							# nonamer2
							
							# SITE 1
							# What is each change to ORF1 CODON SITE 1 / ORF2 CODON SITE 2?
							if($nt ne $nt1_nonamer2_WT) { # only one possibility for this site
								 # only one possibility for this site
								
								my $nonamer2_STOP_caused = 0;
								
								my $codon_nonamer2_ORF1_MUT = $codon_nonamer2_ORF1;
								$codon_nonamer2_ORF1_MUT =~ s/[ACGT]([ACGT])([ACGT])/$nt$1$2/;
								
								my $codon_nonamer2_ORF2_prev_MUT = $codon_nonamer2_ORF2_prev;
								$codon_nonamer2_ORF2_prev_MUT =~ s/([ACGT])[ACGT]([ACGT])/$1$nt_revcom$2/;
								
								my $nonamer2_ORF1_effect = 'S';
								
								#nonamer1-ORF1
								my $AA_nonamer2_ORF1_MUT = get_amino_acid($codon_nonamer2_ORF1_MUT);
								
								if($AA_nonamer2_ORF1_MUT ne $AA_nonamer2_ORF1) {
									$nonamer2_ORF1_effect = 'N';
								}
								
								if($AA_nonamer2_ORF1 eq '*' || $AA_nonamer2_ORF1_MUT eq '*') {
									$nonamer2_STOP_caused++;
								}
								
								
								my $nonamer2_ORF2_prev_effect = 'S';
								
								#nonamer1-ORF2
								my $AA_nonamer2_ORF2_prev_MUT = get_amino_acid($codon_nonamer2_ORF2_prev_MUT);
								
								if($AA_nonamer2_ORF2_prev_MUT ne $AA_nonamer2_ORF2_prev) {
									$nonamer2_ORF2_prev_effect = 'N';
								}
								
								if($AA_nonamer2_ORF2_prev eq '*' || $AA_nonamer2_ORF2_prev_MUT eq '*') {
									$nonamer2_STOP_caused++;
								}
								
								# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
								unless($nonamer2_STOP_caused > 0) {
									$nonamer2_num_changes_poss_site1++;
									
									if($nonamer2_ORF1_effect eq 'N') {
										if($nonamer2_ORF2_prev_effect eq 'N') {
											$nonamer2_num_changes_NN_site1++;
											
											
										} elsif($nonamer2_ORF2_prev_effect eq 'S') {
											$nonamer2_num_changes_NS_site1++;
											
										}
									} elsif($nonamer2_ORF1_effect eq 'S') {
										if($nonamer2_ORF2_prev_effect eq 'N') {
											$nonamer2_num_changes_SN_site1++;
											
										} elsif($nonamer2_ORF2_prev_effect eq 'S') {
											$nonamer2_num_changes_SS_site1++;
											
										}
									}
								}
							} # end member 2 site 1
							
							
							# SITE 2
							# What is each change to ORF1 CODON SITE 2 / ORF2 CODON SITE 3?
							if($nt ne $nt2_nonamer2_WT) { # only one possibility for this site
								
								my $nonamer2_STOP_caused = 0;
								
								my $codon_nonamer2_ORF1_MUT = $codon_nonamer2_ORF1;
								$codon_nonamer2_ORF1_MUT =~ s/([ACGT])[ACGT]([ACGT])/$1$nt$2/;
								
								my $codon_nonamer2_ORF2_prev_MUT = $codon_nonamer2_ORF2_prev;
								$codon_nonamer2_ORF2_prev_MUT =~ s/[ACGT]([ACGT])([ACGT])/$nt_revcom$1$2/;
								
								my $nonamer2_ORF1_effect = 'S';
								
								#nonamer1-ORF1
								my $AA_nonamer2_ORF1_MUT = get_amino_acid($codon_nonamer2_ORF1_MUT);
								
								if($AA_nonamer2_ORF1_MUT ne $AA_nonamer2_ORF1) {
									$nonamer2_ORF1_effect = 'N';
								}
								
								if($AA_nonamer2_ORF1 eq '*' || $AA_nonamer2_ORF1_MUT eq '*') {
									$nonamer2_STOP_caused++;
								}
								
								
								my $nonamer2_ORF2_prev_effect = 'S';
								
								#nonamer1-ORF2
								my $AA_nonamer2_ORF2_prev_MUT = get_amino_acid($codon_nonamer2_ORF2_prev_MUT);
								
								if($AA_nonamer2_ORF2_prev_MUT ne $AA_nonamer2_ORF2_prev) {
									$nonamer2_ORF2_prev_effect = 'N';
								}
								
								if($AA_nonamer2_ORF2_prev eq '*' || $AA_nonamer2_ORF2_prev_MUT eq '*') {
									$nonamer2_STOP_caused++;
								}
								
								# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
								unless($nonamer2_STOP_caused > 0) {
									$nonamer2_num_changes_poss_site2++;
									
									if($nonamer2_ORF1_effect eq 'N') {
										if($nonamer2_ORF2_prev_effect eq 'N') {
											$nonamer2_num_changes_NN_site2++;
											
										} elsif($nonamer2_ORF2_prev_effect eq 'S') {
											$nonamer2_num_changes_NS_site2++;
											
										}
									} elsif($nonamer2_ORF1_effect eq 'S') {
										if($nonamer2_ORF2_prev_effect eq 'N') {
											$nonamer2_num_changes_SN_site2++;
											
										} elsif($nonamer2_ORF2_prev_effect eq 'S') {
											$nonamer2_num_changes_SS_site2++;
											
										}
									}
								}
							} # end member 2 site 2
							
							
							# SITE 3
							# What is each change to ORF1 CODON SITE 3 / ORF2 (next) CODON SITE 1?
							if($nt ne $nt3_nonamer2_WT) { # only one possibility for this site
								
								my $nonamer2_STOP_caused = 0;
								
								my $codon_nonamer2_ORF1_MUT = $codon_nonamer2_ORF1;
								$codon_nonamer2_ORF1_MUT =~ s/([ACGT])([ACGT])[ACGT]/$1$2$nt/;
								
								my $codon_nonamer2_ORF2_next_MUT = $codon_nonamer2_ORF2_next;
								$codon_nonamer2_ORF2_next_MUT =~ s/([ACGT])([ACGT])[ACGT]/$1$2$nt_revcom/;
								
								my $nonamer2_ORF1_effect = 'S';
								
								#nonamer2-ORF1
								my $AA_nonamer2_ORF1_MUT = get_amino_acid($codon_nonamer2_ORF1_MUT);
								
								if($AA_nonamer2_ORF1_MUT ne $AA_nonamer2_ORF1) {
									$nonamer2_ORF1_effect = 'N';
								}
								
								if($AA_nonamer2_ORF1 eq '*' || $AA_nonamer2_ORF1_MUT eq '*') {
									$nonamer2_STOP_caused++;
								}
								
								
								my $nonamer2_ORF2_next_effect = 'S';
								
								#nonamer2-ORF2
								my $AA_nonamer2_ORF2_next_MUT = get_amino_acid($codon_nonamer2_ORF2_next_MUT);
								
								if($AA_nonamer2_ORF2_next_MUT ne $AA_nonamer2_ORF2_next) {
									$nonamer2_ORF2_next_effect = 'N';
								}
								
								if($AA_nonamer2_ORF2_next eq '*' || $AA_nonamer2_ORF2_next_MUT eq '*') {
									$nonamer2_STOP_caused++;
								}
								
								# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
								unless($nonamer2_STOP_caused > 0) {
									$nonamer2_num_changes_poss_site3++;
									
									if($nonamer2_ORF1_effect eq 'N') {
										if($nonamer2_ORF2_next_effect eq 'N') {
											$nonamer2_num_changes_NN_site3++;
											
										} elsif($nonamer2_ORF2_next_effect eq 'S') {
											$nonamer2_num_changes_NS_site3++;
											
										}
									} elsif($nonamer2_ORF1_effect eq 'S') {
										if($nonamer2_ORF2_next_effect eq 'N') {
											$nonamer2_num_changes_SN_site3++;
											
										} elsif($nonamer2_ORF2_next_effect eq 'S') {
											$nonamer2_num_changes_SS_site3++;
											
										}
									}
								}
							} # end member 2 site 3
						} # end cycling all 4 nucleotides
						
						
						# TALLY SITE 1
						my $NN_sites_nonamer1_site1 = 'NA';
						my $SN_sites_nonamer1_site1 = 'NA';
						my $NS_sites_nonamer1_site1 = 'NA';
						my $SS_sites_nonamer1_site1 = 'NA';
						
						if($nonamer1_num_changes_poss_site1 > 0) {
							$NN_sites_nonamer1_site1 = $nonamer1_num_changes_NN_site1 / $nonamer1_num_changes_poss_site1;
							$SN_sites_nonamer1_site1 = $nonamer1_num_changes_SN_site1 / $nonamer1_num_changes_poss_site1;
							$NS_sites_nonamer1_site1 = $nonamer1_num_changes_NS_site1 / $nonamer1_num_changes_poss_site1;
							$SS_sites_nonamer1_site1 = $nonamer1_num_changes_SS_site1 / $nonamer1_num_changes_poss_site1;
							
							unless(exists $seq_completed{$nonamer1}) {
								$seq2sites{$nonamer1}->{NN_sites} += $NN_sites_nonamer1_site1;
								$seq2sites{$nonamer1}->{SN_sites} += $SN_sites_nonamer1_site1;
								$seq2sites{$nonamer1}->{NS_sites} += $NS_sites_nonamer1_site1;
								$seq2sites{$nonamer1}->{SS_sites} += $SS_sites_nonamer1_site1;
							}
						}
						
						my $NN_sites_nonamer2_site1 = 'NA';
						my $SN_sites_nonamer2_site1 = 'NA';
						my $NS_sites_nonamer2_site1 = 'NA';
						my $SS_sites_nonamer2_site1 = 'NA';
						
						if($nonamer2_num_changes_poss_site1 > 0) {
							$NN_sites_nonamer2_site1 = $nonamer2_num_changes_NN_site1 / $nonamer2_num_changes_poss_site1;
							$SN_sites_nonamer2_site1 = $nonamer2_num_changes_SN_site1 / $nonamer2_num_changes_poss_site1;
							$NS_sites_nonamer2_site1 = $nonamer2_num_changes_NS_site1 / $nonamer2_num_changes_poss_site1;
							$SS_sites_nonamer2_site1 = $nonamer2_num_changes_SS_site1 / $nonamer2_num_changes_poss_site1;
							
							unless(exists $seq_completed{$nonamer2}) {
								$seq2sites{$nonamer2}->{NN_sites} += $NN_sites_nonamer2_site1;
								$seq2sites{$nonamer2}->{SN_sites} += $SN_sites_nonamer2_site1;
								$seq2sites{$nonamer2}->{NS_sites} += $NS_sites_nonamer2_site1;
								$seq2sites{$nonamer2}->{SS_sites} += $SS_sites_nonamer2_site1;
							}
						}
						
						my $NN_sites_site1 = 'NA';
						my $SN_sites_site1 = 'NA';
						my $NS_sites_site1 = 'NA';
						my $SS_sites_site1 = 'NA';
						
						if($nonamer1_num_changes_poss_site1 > 0 && $nonamer2_num_changes_poss_site1 > 0) {
							$NN_sites_site1 = ($NN_sites_nonamer1_site1 + $NN_sites_nonamer2_site1) / 2;
							$SN_sites_site1 = ($SN_sites_nonamer1_site1 + $SN_sites_nonamer2_site1) / 2;
							$NS_sites_site1 = ($NS_sites_nonamer1_site1 + $NS_sites_nonamer2_site1) / 2;
							$SS_sites_site1 = ($SS_sites_nonamer1_site1 + $SS_sites_nonamer2_site1) / 2;
						} elsif($nonamer1_num_changes_poss_site1 > 0) {
							$NN_sites_site1 = $NN_sites_nonamer1_site1;
							$SN_sites_site1 = $SN_sites_nonamer1_site1;
							$NS_sites_site1 = $NS_sites_nonamer1_site1;
							$SS_sites_site1 = $SS_sites_nonamer1_site1;
						} elsif($nonamer2_num_changes_poss_site1 > 0) {
							$NN_sites_site1 = $NN_sites_nonamer2_site1;
							$SN_sites_site1 = $SN_sites_nonamer2_site1;
							$NS_sites_site1 = $NS_sites_nonamer2_site1;
							$SS_sites_site1 = $SS_sites_nonamer2_site1;
							
						} # else nothing defined, stay NA
						
						
						# TALLY SITE 2
						my $NN_sites_nonamer1_site2 = 'NA';
						my $SN_sites_nonamer1_site2 = 'NA';
						my $NS_sites_nonamer1_site2 = 'NA';
						my $SS_sites_nonamer1_site2 = 'NA';
						
						if($nonamer1_num_changes_poss_site2 > 0) {
							$NN_sites_nonamer1_site2 = $nonamer1_num_changes_NN_site2 / $nonamer1_num_changes_poss_site2;
							$SN_sites_nonamer1_site2 = $nonamer1_num_changes_SN_site2 / $nonamer1_num_changes_poss_site2;
							$NS_sites_nonamer1_site2 = $nonamer1_num_changes_NS_site2 / $nonamer1_num_changes_poss_site2;
							$SS_sites_nonamer1_site2 = $nonamer1_num_changes_SS_site2 / $nonamer1_num_changes_poss_site2;
							
							unless(exists $seq_completed{$nonamer1}) {
								$seq2sites{$nonamer1}->{NN_sites} += $NN_sites_nonamer1_site2;
								$seq2sites{$nonamer1}->{SN_sites} += $SN_sites_nonamer1_site2;
								$seq2sites{$nonamer1}->{NS_sites} += $NS_sites_nonamer1_site2;
								$seq2sites{$nonamer1}->{SS_sites} += $SS_sites_nonamer1_site2;
							}
						}
						
						my $NN_sites_nonamer2_site2 = 'NA';
						my $SN_sites_nonamer2_site2 = 'NA';
						my $NS_sites_nonamer2_site2 = 'NA';
						my $SS_sites_nonamer2_site2 = 'NA';
						
						if($nonamer2_num_changes_poss_site2 > 0) {
							$NN_sites_nonamer2_site2 = $nonamer2_num_changes_NN_site2 / $nonamer2_num_changes_poss_site2;
							$SN_sites_nonamer2_site2 = $nonamer2_num_changes_SN_site2 / $nonamer2_num_changes_poss_site2;
							$NS_sites_nonamer2_site2 = $nonamer2_num_changes_NS_site2 / $nonamer2_num_changes_poss_site2;
							$SS_sites_nonamer2_site2 = $nonamer2_num_changes_SS_site2 / $nonamer2_num_changes_poss_site2;
							
							unless(exists $seq_completed{$nonamer2}) {
								$seq2sites{$nonamer2}->{NN_sites} += $NN_sites_nonamer2_site2;
								$seq2sites{$nonamer2}->{SN_sites} += $SN_sites_nonamer2_site2;
								$seq2sites{$nonamer2}->{NS_sites} += $NS_sites_nonamer2_site2;
								$seq2sites{$nonamer2}->{SS_sites} += $SS_sites_nonamer2_site2;
							}
						}
						
						my $NN_sites_site2 = 'NA';
						my $SN_sites_site2 = 'NA';
						my $NS_sites_site2 = 'NA';
						my $SS_sites_site2 = 'NA';
						
						if($nonamer1_num_changes_poss_site2 > 0 && $nonamer2_num_changes_poss_site2 > 0) {
							$NN_sites_site2 = ($NN_sites_nonamer1_site2 + $NN_sites_nonamer2_site2) / 2;
							$SN_sites_site2 = ($SN_sites_nonamer1_site2 + $SN_sites_nonamer2_site2) / 2;
							$NS_sites_site2 = ($NS_sites_nonamer1_site2 + $NS_sites_nonamer2_site2) / 2;
							$SS_sites_site2 = ($SS_sites_nonamer1_site2 + $SS_sites_nonamer2_site2) / 2;
						} elsif($nonamer1_num_changes_poss_site2 > 0) {
							$NN_sites_site2 = $NN_sites_nonamer1_site2;
							$SN_sites_site2 = $SN_sites_nonamer1_site2;
							$NS_sites_site2 = $NS_sites_nonamer1_site2;
							$SS_sites_site2 = $SS_sites_nonamer1_site2;
						} elsif($nonamer2_num_changes_poss_site2 > 0) {
							$NN_sites_site2 = $NN_sites_nonamer2_site2;
							$SN_sites_site2 = $SN_sites_nonamer2_site2;
							$NS_sites_site2 = $NS_sites_nonamer2_site2;
							$SS_sites_site2 = $SS_sites_nonamer2_site2;
							
						} # else nothing defined, stay NA
						
						
						# TALLY SITE 3
						my $NN_sites_nonamer1_site3 = 'NA';
						my $SN_sites_nonamer1_site3 = 'NA';
						my $NS_sites_nonamer1_site3 = 'NA';
						my $SS_sites_nonamer1_site3 = 'NA';
						
						if($nonamer1_num_changes_poss_site3 > 0) {
							$NN_sites_nonamer1_site3 = $nonamer1_num_changes_NN_site3 / $nonamer1_num_changes_poss_site3;
							$SN_sites_nonamer1_site3 = $nonamer1_num_changes_SN_site3 / $nonamer1_num_changes_poss_site3;
							$NS_sites_nonamer1_site3 = $nonamer1_num_changes_NS_site3 / $nonamer1_num_changes_poss_site3;
							$SS_sites_nonamer1_site3 = $nonamer1_num_changes_SS_site3 / $nonamer1_num_changes_poss_site3;
							
							unless(exists $seq_completed{$nonamer1}) {
								$seq2sites{$nonamer1}->{NN_sites} += $NN_sites_nonamer1_site3;
								$seq2sites{$nonamer1}->{SN_sites} += $SN_sites_nonamer1_site3;
								$seq2sites{$nonamer1}->{NS_sites} += $NS_sites_nonamer1_site3;
								$seq2sites{$nonamer1}->{SS_sites} += $SS_sites_nonamer1_site3;
							}
						}
						
						my $NN_sites_nonamer2_site3 = 'NA';
						my $SN_sites_nonamer2_site3 = 'NA';
						my $NS_sites_nonamer2_site3 = 'NA';
						my $SS_sites_nonamer2_site3 = 'NA';
						
						if($nonamer2_num_changes_poss_site3 > 0) {
							$NN_sites_nonamer2_site3 = $nonamer2_num_changes_NN_site3 / $nonamer2_num_changes_poss_site3;
							$SN_sites_nonamer2_site3 = $nonamer2_num_changes_SN_site3 / $nonamer2_num_changes_poss_site3;
							$NS_sites_nonamer2_site3 = $nonamer2_num_changes_NS_site3 / $nonamer2_num_changes_poss_site3;
							$SS_sites_nonamer2_site3 = $nonamer2_num_changes_SS_site3 / $nonamer2_num_changes_poss_site3;
							
							unless(exists $seq_completed{$nonamer2}) {
								$seq2sites{$nonamer2}->{NN_sites} += $NN_sites_nonamer2_site3;
								$seq2sites{$nonamer2}->{SN_sites} += $SN_sites_nonamer2_site3;
								$seq2sites{$nonamer2}->{NS_sites} += $NS_sites_nonamer2_site3;
								$seq2sites{$nonamer2}->{SS_sites} += $SS_sites_nonamer2_site3;
							}
						}
						
						my $NN_sites_site3 = 'NA';
						my $SN_sites_site3 = 'NA';
						my $NS_sites_site3 = 'NA';
						my $SS_sites_site3 = 'NA';
						
						if($nonamer1_num_changes_poss_site3 > 0 && $nonamer2_num_changes_poss_site3 > 0) {
							$NN_sites_site3 = ($NN_sites_nonamer1_site3 + $NN_sites_nonamer2_site3) / 2;
							$SN_sites_site3 = ($SN_sites_nonamer1_site3 + $SN_sites_nonamer2_site3) / 2;
							$NS_sites_site3 = ($NS_sites_nonamer1_site3 + $NS_sites_nonamer2_site3) / 2;
							$SS_sites_site3 = ($SS_sites_nonamer1_site3 + $SS_sites_nonamer2_site3) / 2;
						} elsif($nonamer1_num_changes_poss_site3 > 0) {
							$NN_sites_site3 = $NN_sites_nonamer1_site3;
							$SN_sites_site3 = $SN_sites_nonamer1_site3;
							$NS_sites_site3 = $NS_sites_nonamer1_site3;
							$SS_sites_site3 = $SS_sites_nonamer1_site3;
						} elsif($nonamer2_num_changes_poss_site3 > 0) {
							$NN_sites_site3 = $NN_sites_nonamer2_site3;
							$SN_sites_site3 = $SN_sites_nonamer2_site3;
							$NS_sites_site3 = $NS_sites_nonamer2_site3;
							$SS_sites_site3 = $SS_sites_nonamer2_site3;
							
						} # else nothing defined, stay NA
						
						# SUM THE THREE SITES
						my $NN_sites = $NN_sites_site1 + $NN_sites_site2 + $NN_sites_site3;
						my $SN_sites = $SN_sites_site1 + $SN_sites_site2 + $SN_sites_site3;
						my $NS_sites = $NS_sites_site1 + $NS_sites_site2 + $NS_sites_site3;
						my $SS_sites = $SS_sites_site1 + $SS_sites_site2 + $SS_sites_site3;
						
						# Check if there are multiple variants in these overlapping codons
						my $MNV = 'FALSE';
						if(($NN_diffs + $SN_diffs + $NS_diffs + $SS_diffs) > 1) {
							$MNV = 'TRUE';
							$this_codon_MNV = 'TRUE'; # only needs to be TRUE in one comparison
						}
						
						# NEW NAIVE APPROACH: ADD TO SUMS
						my $comparison_weight = 0; # number of pairwise comparisons involving these nonamers
						if($nonamer1_index == $nonamer2_index) { # if($nonamer1 eq $nonamer2) { # it's a self-comparison
							# COMBINATION nC2
							$comparison_weight = (($nonamer1_count * $nonamer1_count) - $nonamer1_count) / 2; # nonamer 1 and 2 are the same
						} else {
							# MULTIPLICATION m * n
							my $nonamer2_count = $unique_9mers{$seq_site_index}->{$nonamer2};
							$comparison_weight = $nonamer1_count * $nonamer2_count;
						}
						
						$comparisons_sum += $comparison_weight;
						
						$NN_sites_numerator += $NN_sites * $comparison_weight;
						$SN_sites_numerator += $SN_sites * $comparison_weight;
						$NS_sites_numerator += $NS_sites * $comparison_weight;
						$SS_sites_numerator += $SS_sites * $comparison_weight;
						
						$NN_diffs_numerator += $NN_diffs * $comparison_weight;
						$SN_diffs_numerator += $SN_diffs * $comparison_weight;
						$NS_diffs_numerator += $NS_diffs * $comparison_weight;
						$SS_diffs_numerator += $SS_diffs * $comparison_weight;
						
					} # MIDDLE (internal) codon (not first or last; two ORF2 codons overlap)
						
					
##################################################################################################
##################### SENSE-ANTISENSE:
#####################  sas13:
#####################    ORF1: 1-2-3-1-2-3-1
#####################    ORF2: 3-2-1-3-2-1-3
##################################################################################################
				} elsif($frame eq 'sas13') {
					
					# EXACT OVERLAP (opposite strand). Always 3 sites to examine (first, middle, last codon)
					# FIRST AND LAST CODONS ARE NOT DIFFERENT
					my $codon_nonamer1_ORF2 = revcom(substr($nonamer1, $site_index, 3));
					my $codon_nonamer2_ORF2 = revcom(substr($nonamer2, $site_index, 3));
					
					if($verbose_messages) {
						print "site_index=$site_index\n";
						print "nonamer1=$nonamer1\n";
						print "nonamer2=$nonamer2\n";
						print "nonamer1_count=$nonamer1_count\n";
						#print "nonamer2=$nonamer2_count\n";
						print "codon_nonamer1_ORF2=$codon_nonamer1_ORF2\n";
						print "codon_nonamer2_ORF2=$codon_nonamer2_ORF2\n";
					}
					
					my $AA_nonamer1_ORF2 = get_amino_acid($codon_nonamer1_ORF2);
					my $AA_nonamer2_ORF2 = get_amino_acid($codon_nonamer2_ORF2);
					
					if($codon_num > 2 && ($AA_nonamer1_ORF2 eq '*' || $AA_nonamer2_ORF2 eq '*')) {
						print "### WARNING! ORF2, $nonamer1\-$nonamer2 comparison, near ORF1 codon $codon_in_seq encodes a within-frame STOP codon. Wrong frame selection ($frame)?\n";
					}
					
					#######################
					# GET NUMBER OF SITES
					
					# New phylogeny-naive approach
					$alt1_codons_to_counts{$codon_nonamer1_ORF2} += $nonamer1_count;
					$alt1_codons_to_counts{$codon_nonamer2_ORF2} += $nonamer1_count; # ADDED
					#$alt2_codons_to_counts{$codon_nonamer1_ORF2} += $nonamer1_count; # COMEBACK -- make sure it's counting correctly for this frame
					#$alt1_codons_to_counts{$codon_nonamer1_ORF2} += $nonamer1_count; # these store redundant info for this frame
					# always nonamer1_count because otherwise redundant?
					
					my $nonamer1_num_changes_poss_site1 = 0;
					my $nonamer1_num_changes_NN_site1 = 0;
					my $nonamer1_num_changes_SN_site1 = 0;
					my $nonamer1_num_changes_NS_site1 = 0;
					my $nonamer1_num_changes_SS_site1 = 0;
					
					my $nonamer1_num_changes_poss_site2 = 0;
					my $nonamer1_num_changes_NN_site2 = 0;
					my $nonamer1_num_changes_SN_site2 = 0;
					my $nonamer1_num_changes_NS_site2 = 0;
					my $nonamer1_num_changes_SS_site2 = 0;
					
					my $nonamer1_num_changes_poss_site3 = 0;
					my $nonamer1_num_changes_NN_site3 = 0;
					my $nonamer1_num_changes_SN_site3 = 0;
					my $nonamer1_num_changes_NS_site3 = 0;
					my $nonamer1_num_changes_SS_site3 = 0;
					
					my $nonamer2_num_changes_poss_site1 = 0;
					my $nonamer2_num_changes_NN_site1 = 0;
					my $nonamer2_num_changes_SN_site1 = 0;
					my $nonamer2_num_changes_NS_site1 = 0;
					my $nonamer2_num_changes_SS_site1 = 0;
					
					my $nonamer2_num_changes_poss_site2 = 0;
					my $nonamer2_num_changes_NN_site2 = 0;
					my $nonamer2_num_changes_SN_site2 = 0;
					my $nonamer2_num_changes_NS_site2 = 0;
					my $nonamer2_num_changes_SS_site2 = 0;
					
					my $nonamer2_num_changes_poss_site3 = 0;
					my $nonamer2_num_changes_NN_site3 = 0;
					my $nonamer2_num_changes_SN_site3 = 0;
					my $nonamer2_num_changes_NS_site3 = 0;
					my $nonamer2_num_changes_SS_site3 = 0;
					
					my $NN_diffs = 0;
					my $SN_diffs = 0;
					my $NS_diffs = 0;
					my $SS_diffs = 0;
					
					# Just ORF1 (reference)
					my $nt1_nonamer1_WT = substr($codon_nonamer1_ORF1, 0, 1);
					my $nt2_nonamer1_WT = substr($codon_nonamer1_ORF1, 1, 1);
					my $nt3_nonamer1_WT = substr($codon_nonamer1_ORF1, 2, 1);
					my $nt1_nonamer2_WT = substr($codon_nonamer2_ORF1, 0, 1);
					my $nt2_nonamer2_WT = substr($codon_nonamer2_ORF1, 1, 1);
					my $nt3_nonamer2_WT = substr($codon_nonamer2_ORF1, 2, 1);
					
					foreach my $nt (@nucleotides) {
						my $nt_revcom = revcom($nt);
						
						##################################################################
						# nonamer1
						
						# SITE 1
						# What is each change to ORF1 CODON SITE 1 / ORF2 CODON SITE 3?
						if($nt ne $nt1_nonamer1_WT) { # only one possibility for this site
							 # only one possibility for this site
							
							my $nonamer1_STOP_caused = 0;
							
							my $codon_nonamer1_ORF1_MUT = $codon_nonamer1_ORF1;
							$codon_nonamer1_ORF1_MUT =~ s/[ACGT]([ACGT])([ACGT])/$nt$1$2/;
							
							my $codon_nonamer1_ORF2_MUT = $codon_nonamer1_ORF2;
							$codon_nonamer1_ORF2_MUT =~ s/([ACGT])([ACGT])[ACGT]/$1$2$nt_revcom/;
							
							my $nonamer1_ORF1_effect = 'S';
							
							#nonamer1-ORF1
							my $AA_nonamer1_ORF1_MUT = get_amino_acid($codon_nonamer1_ORF1_MUT);
							
							if($AA_nonamer1_ORF1_MUT ne $AA_nonamer1_ORF1) {
								$nonamer1_ORF1_effect = 'N';
							}
							
							if($AA_nonamer1_ORF1 eq '*' || $AA_nonamer1_ORF1_MUT eq '*') {
								$nonamer1_STOP_caused++;
							}
							
							my $nonamer1_ORF2_effect = 'S';
							
							#nonamer1-ORF2
							my $AA_nonamer1_ORF2_MUT = get_amino_acid($codon_nonamer1_ORF2_MUT);
							
							if($AA_nonamer1_ORF2_MUT ne $AA_nonamer1_ORF2) {
								$nonamer1_ORF2_effect = 'N';
							}
							
							if($AA_nonamer1_ORF2 eq '*' || $AA_nonamer1_ORF2_MUT eq '*') {
								$nonamer1_STOP_caused++;
							}
							
							# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
							unless($nonamer1_STOP_caused > 0) {
								$nonamer1_num_changes_poss_site1++;
								
								if($nonamer1_ORF1_effect eq 'N') {
									if($nonamer1_ORF2_effect eq 'N') {
										$nonamer1_num_changes_NN_site1++;
										
										# ACTUAL DIFF
										if($nt eq $nt1_nonamer2_WT) {
											$NN_diffs++;
											$site_diffs_hh{$codon_num}->{1}->{NN_diffs}++;
										}
										
									} elsif($nonamer1_ORF2_effect eq 'S') {
										$nonamer1_num_changes_NS_site1++;
										
										# ACTUAL DIFF
										if($nt eq $nt1_nonamer2_WT) {
											$NS_diffs++;
											$site_diffs_hh{$codon_num}->{1}->{NS_diffs}++;
										}
									}
								} elsif($nonamer1_ORF1_effect eq 'S') {
									if($nonamer1_ORF2_effect eq 'N') {
										$nonamer1_num_changes_SN_site1++;
										
										# ACTUAL DIFF
										if($nt eq $nt1_nonamer2_WT) {
											$SN_diffs++;
											$site_diffs_hh{$codon_num}->{1}->{SN_diffs}++;
										}
									} elsif($nonamer1_ORF2_effect eq 'S') {
										$nonamer1_num_changes_SS_site1++;
										
										# ACTUAL DIFF
										if($nt eq $nt1_nonamer2_WT) {
											$SS_diffs++;
											$site_diffs_hh{$codon_num}->{1}->{SS_diffs}++;
										}
									}
								}
							}
						} # end member 1 site 1
						
						
						# SITE 2
						# What is each change to ORF1 CODON SITE 2 / ORF2 CODON SITE 2?
						if($nt ne $nt2_nonamer1_WT) { # only one possibility for this site
							
							my $nonamer1_STOP_caused = 0;
							
							my $codon_nonamer1_ORF1_MUT = $codon_nonamer1_ORF1;
							$codon_nonamer1_ORF1_MUT =~ s/([ACGT])[ACGT]([ACGT])/$1$nt$2/;
							
							my $codon_nonamer1_ORF2_MUT = $codon_nonamer1_ORF2;
							$codon_nonamer1_ORF2_MUT =~ s/([ACGT])[ACGT]([ACGT])/$1$nt_revcom$2/;
							
							my $nonamer1_ORF1_effect = 'S';
							
							#nonamer1-ORF1
							my $AA_nonamer1_ORF1_MUT = get_amino_acid($codon_nonamer1_ORF1_MUT);
							
							if($AA_nonamer1_ORF1_MUT ne $AA_nonamer1_ORF1) {
								$nonamer1_ORF1_effect = 'N';
							}
							
							if($AA_nonamer1_ORF1 eq '*' || $AA_nonamer1_ORF1_MUT eq '*') {
								$nonamer1_STOP_caused++;
							}
							
							
							my $nonamer1_ORF2_effect = 'S';
							
							#nonamer1-ORF2
							my $AA_nonamer1_ORF2_MUT = get_amino_acid($codon_nonamer1_ORF2_MUT);
							
							if($AA_nonamer1_ORF2_MUT ne $AA_nonamer1_ORF2) {
								$nonamer1_ORF2_effect = 'N';
							}
							
							if($AA_nonamer1_ORF2 eq '*' || $AA_nonamer1_ORF2_MUT eq '*') {
								$nonamer1_STOP_caused++;
							}
							
							# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
							unless($nonamer1_STOP_caused > 0) {
								$nonamer1_num_changes_poss_site2++;
								
								if($nonamer1_ORF1_effect eq 'N') {
									if($nonamer1_ORF2_effect eq 'N') {
										$nonamer1_num_changes_NN_site2++;
										
										# ACTUAL DIFF
										if($nt eq $nt2_nonamer2_WT) {
											$NN_diffs++;
											$site_diffs_hh{$codon_num}->{2}->{NN_diffs}++;
										}
										
									} elsif($nonamer1_ORF2_effect eq 'S') {
										$nonamer1_num_changes_NS_site2++;
										
										# ACTUAL DIFF
										if($nt eq $nt2_nonamer2_WT) {
											$NS_diffs++;
											$site_diffs_hh{$codon_num}->{2}->{NS_diffs}++;
										}
									}
								} elsif($nonamer1_ORF1_effect eq 'S') {
									if($nonamer1_ORF2_effect eq 'N') {
										$nonamer1_num_changes_SN_site2++;
										
										# ACTUAL DIFF
										if($nt eq $nt2_nonamer2_WT) {
											$SN_diffs++;
											$site_diffs_hh{$codon_num}->{2}->{SN_diffs}++;
										}
									} elsif($nonamer1_ORF2_effect eq 'S') {
										$nonamer1_num_changes_SS_site2++;
										
										# ACTUAL DIFF
										if($nt eq $nt2_nonamer2_WT) {
											$SS_diffs++;
											$site_diffs_hh{$codon_num}->{2}->{SS_diffs}++;
										}
									}
								}
							}
						} # end member 1 site 2
						
						
						# SITE 3
						# What is each change to ORF1 CODON SITE 3 / ORF2 CODON SITE 1?
						if($nt ne $nt3_nonamer1_WT) { # only one possibility for this site
							
							my $nonamer1_STOP_caused = 0;
							
							my $codon_nonamer1_ORF1_MUT = $codon_nonamer1_ORF1;
							$codon_nonamer1_ORF1_MUT =~ s/([ACGT])([ACGT])[ACGT]/$1$2$nt/;
							
							my $codon_nonamer1_ORF2_MUT = $codon_nonamer1_ORF2;
							$codon_nonamer1_ORF2_MUT =~ s/[ACGT]([ACGT])([ACGT])/$nt_revcom$1$2/;
							
							my $nonamer1_ORF1_effect = 'S';
							
							#nonamer1-ORF1
							my $AA_nonamer1_ORF1_MUT = get_amino_acid($codon_nonamer1_ORF1_MUT);
							
							if($AA_nonamer1_ORF1_MUT ne $AA_nonamer1_ORF1) {
								$nonamer1_ORF1_effect = 'N';
							}
							
							if($AA_nonamer1_ORF1 eq '*' || $AA_nonamer1_ORF1_MUT eq '*') {
								$nonamer1_STOP_caused++;
							}
							
							my $nonamer1_ORF2_effect = 'S';
							
							#nonamer1-ORF2
							my $AA_nonamer1_ORF2_MUT = get_amino_acid($codon_nonamer1_ORF2_MUT);
							
							if($AA_nonamer1_ORF2_MUT ne $AA_nonamer1_ORF2) {
								$nonamer1_ORF2_effect = 'N';
							}
							
							if($AA_nonamer1_ORF2 eq '*' || $AA_nonamer1_ORF2_MUT eq '*') {
								$nonamer1_STOP_caused++;
							}
							
							# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
							unless($nonamer1_STOP_caused > 0) {
								$nonamer1_num_changes_poss_site3++;
								
								if($nonamer1_ORF1_effect eq 'N') {
									if($nonamer1_ORF2_effect eq 'N') {
										$nonamer1_num_changes_NN_site3++;
										
										# ACTUAL DIFF
										if($nt eq $nt3_nonamer2_WT) {
											$NN_diffs++;
											$site_diffs_hh{$codon_num}->{3}->{NN_diffs}++;
										}
										
									} elsif($nonamer1_ORF2_effect eq 'S') {
										$nonamer1_num_changes_NS_site3++;
										
										# ACTUAL DIFF
										if($nt eq $nt3_nonamer2_WT) {
											$NS_diffs++;
											$site_diffs_hh{$codon_num}->{3}->{NS_diffs}++;
										}
									}
								} elsif($nonamer1_ORF1_effect eq 'S') {
									if($nonamer1_ORF2_effect eq 'N') {
										$nonamer1_num_changes_SN_site3++;
										
										# ACTUAL DIFF
										if($nt eq $nt3_nonamer2_WT) {
											$SN_diffs++;
											$site_diffs_hh{$codon_num}->{3}->{SN_diffs}++;
										}
									} elsif($nonamer1_ORF2_effect eq 'S') {
										$nonamer1_num_changes_SS_site3++;
										
										# ACTUAL DIFF
										if($nt eq $nt3_nonamer2_WT) {
											$SS_diffs++;
											$site_diffs_hh{$codon_num}->{3}->{SS_diffs}++;
										}
									}
								}
							}
						} # end member 1 site 3
						
						
						##################################################################
						# nonamer2
						
						# SITE 1
						# What is each change to ORF1 CODON SITE 1 / ORF2 CODON SITE 3?
						if($nt ne $nt1_nonamer2_WT) { # only one possibility for this site
							 # only one possibility for this site
							
							my $nonamer2_STOP_caused = 0;
							
							my $codon_nonamer2_ORF1_MUT = $codon_nonamer2_ORF1;
							$codon_nonamer2_ORF1_MUT =~ s/[ACGT]([ACGT])([ACGT])/$nt$1$2/;
							
							my $codon_nonamer2_ORF2_MUT = $codon_nonamer2_ORF2;
							$codon_nonamer2_ORF2_MUT =~ s/([ACGT])([ACGT])[ACGT]/$1$2$nt_revcom/;
							
							my $nonamer2_ORF1_effect = 'S';
							
							#nonamer1-ORF1
							my $AA_nonamer2_ORF1_MUT = get_amino_acid($codon_nonamer2_ORF1_MUT);
							
							if($AA_nonamer2_ORF1_MUT ne $AA_nonamer2_ORF1) {
								$nonamer2_ORF1_effect = 'N';
							}
							
							if($AA_nonamer2_ORF1 eq '*' || $AA_nonamer2_ORF1_MUT eq '*') {
								$nonamer2_STOP_caused++;
							}
							
							my $nonamer2_ORF2_effect = 'S';
							
							#nonamer1-ORF2
							my $AA_nonamer2_ORF2_MUT = get_amino_acid($codon_nonamer2_ORF2_MUT);
							
							if($AA_nonamer2_ORF2_MUT ne $AA_nonamer2_ORF2) {
								$nonamer2_ORF2_effect = 'N';
							}
							
							if($AA_nonamer2_ORF2 eq '*' || $AA_nonamer2_ORF2_MUT eq '*') {
								$nonamer2_STOP_caused++;
							}
							
							# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
							unless($nonamer2_STOP_caused > 0) {
								$nonamer2_num_changes_poss_site1++;
								
								if($nonamer2_ORF1_effect eq 'N') {
									if($nonamer2_ORF2_effect eq 'N') {
										$nonamer2_num_changes_NN_site1++;
										
										
									} elsif($nonamer2_ORF2_effect eq 'S') {
										$nonamer2_num_changes_NS_site1++;
										
									}
								} elsif($nonamer2_ORF1_effect eq 'S') {
									if($nonamer2_ORF2_effect eq 'N') {
										$nonamer2_num_changes_SN_site1++;
										
									} elsif($nonamer2_ORF2_effect eq 'S') {
										$nonamer2_num_changes_SS_site1++;
										
									}
								}
							}
						} # end member 2 site 1
						
						
						# SITE 2
						# What is each change to ORF1 CODON SITE 2 / ORF2 CODON SITE 2?
						if($nt ne $nt2_nonamer2_WT) { # only one possibility for this site
							
							my $nonamer2_STOP_caused = 0;
							
							my $codon_nonamer2_ORF1_MUT = $codon_nonamer2_ORF1;
							$codon_nonamer2_ORF1_MUT =~ s/([ACGT])[ACGT]([ACGT])/$1$nt$2/;
							
							my $codon_nonamer2_ORF2_MUT = $codon_nonamer2_ORF2;
							$codon_nonamer2_ORF2_MUT =~ s/([ACGT])[ACGT]([ACGT])/$1$nt_revcom$2/;
							
							my $nonamer2_ORF1_effect = 'S';
							
							#nonamer1-ORF1
							my $AA_nonamer2_ORF1_MUT = get_amino_acid($codon_nonamer2_ORF1_MUT);
							
							if($AA_nonamer2_ORF1_MUT ne $AA_nonamer2_ORF1) {
								$nonamer2_ORF1_effect = 'N';
							}
							
							if($AA_nonamer2_ORF1 eq '*' || $AA_nonamer2_ORF1_MUT eq '*') {
								$nonamer2_STOP_caused++;
							}
							
							my $nonamer2_ORF2_effect = 'S';
							
							#nonamer1-ORF2
							my $AA_nonamer2_ORF2_MUT = get_amino_acid($codon_nonamer2_ORF2_MUT);
							
							if($AA_nonamer2_ORF2_MUT ne $AA_nonamer2_ORF2) {
								$nonamer2_ORF2_effect = 'N';
							}
							
							if($AA_nonamer2_ORF2 eq '*' || $AA_nonamer2_ORF2_MUT eq '*') {
								$nonamer2_STOP_caused++;
							}
							
							# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
							unless($nonamer2_STOP_caused > 0) {
								$nonamer2_num_changes_poss_site2++;
								
								if($nonamer2_ORF1_effect eq 'N') {
									if($nonamer2_ORF2_effect eq 'N') {
										$nonamer2_num_changes_NN_site2++;
										
									} elsif($nonamer2_ORF2_effect eq 'S') {
										$nonamer2_num_changes_NS_site2++;
										
									}
								} elsif($nonamer2_ORF1_effect eq 'S') {
									if($nonamer2_ORF2_effect eq 'N') {
										$nonamer2_num_changes_SN_site2++;
										
									} elsif($nonamer2_ORF2_effect eq 'S') {
										$nonamer2_num_changes_SS_site2++;
										
									}
								}
							}
						} # end member 2 site 2
						
						
						# SITE 3
						# What is each change to ORF1 CODON SITE 3 / ORF2 CODON SITE 1?
						if($nt ne $nt3_nonamer2_WT) { # only one possibility for this site
							
							my $nonamer2_STOP_caused = 0;
							
							my $codon_nonamer2_ORF1_MUT = $codon_nonamer2_ORF1;
							$codon_nonamer2_ORF1_MUT =~ s/([ACGT])([ACGT])[ACGT]/$1$2$nt/;
							
							my $codon_nonamer2_ORF2_MUT = $codon_nonamer2_ORF2;
							$codon_nonamer2_ORF2_MUT =~ s/[ACGT]([ACGT])([ACGT])/$nt_revcom$1$2/;
							
							my $nonamer2_ORF1_effect = 'S';
							
							#nonamer2-ORF1
							my $AA_nonamer2_ORF1_MUT = get_amino_acid($codon_nonamer2_ORF1_MUT);
							
							if($AA_nonamer2_ORF1_MUT ne $AA_nonamer2_ORF1) {
								$nonamer2_ORF1_effect = 'N';
							}
							
							if($AA_nonamer2_ORF1 eq '*' || $AA_nonamer2_ORF1_MUT eq '*') {
								$nonamer2_STOP_caused++;
							}
							
							my $nonamer2_ORF2_effect = 'S';
							
							#nonamer2-ORF2
							my $AA_nonamer2_ORF2_MUT = get_amino_acid($codon_nonamer2_ORF2_MUT);
							
							if($AA_nonamer2_ORF2_MUT ne $AA_nonamer2_ORF2) {
								$nonamer2_ORF2_effect = 'N';
							}
							
							if($AA_nonamer2_ORF2 eq '*' || $AA_nonamer2_ORF2_MUT eq '*') {
								$nonamer2_STOP_caused++;
							}
							
							# TALLY VIABLE CHANGES, GET NUMBER OF DIFFS
							unless($nonamer2_STOP_caused > 0) {
								$nonamer2_num_changes_poss_site3++;
								
								if($nonamer2_ORF1_effect eq 'N') {
									if($nonamer2_ORF2_effect eq 'N') {
										$nonamer2_num_changes_NN_site3++;
										
									} elsif($nonamer2_ORF2_effect eq 'S') {
										$nonamer2_num_changes_NS_site3++;
										
									}
								} elsif($nonamer2_ORF1_effect eq 'S') {
									if($nonamer2_ORF2_effect eq 'N') {
										$nonamer2_num_changes_SN_site3++;
										
									} elsif($nonamer2_ORF2_effect eq 'S') {
										$nonamer2_num_changes_SS_site3++;
										
									}
								}
							}
						} # end member 2 site 3
					} # end cycling all 4 nucleotides
					
					
					# TALLY SITE 1
					my $NN_sites_nonamer1_site1 = 'NA';
					my $SN_sites_nonamer1_site1 = 'NA';
					my $NS_sites_nonamer1_site1 = 'NA';
					my $SS_sites_nonamer1_site1 = 'NA';
					
					if($nonamer1_num_changes_poss_site1 > 0) {
						$NN_sites_nonamer1_site1 = $nonamer1_num_changes_NN_site1 / $nonamer1_num_changes_poss_site1;
						$SN_sites_nonamer1_site1 = $nonamer1_num_changes_SN_site1 / $nonamer1_num_changes_poss_site1;
						$NS_sites_nonamer1_site1 = $nonamer1_num_changes_NS_site1 / $nonamer1_num_changes_poss_site1;
						$SS_sites_nonamer1_site1 = $nonamer1_num_changes_SS_site1 / $nonamer1_num_changes_poss_site1;
						
						unless(exists $seq_completed{$nonamer1}) {
							$seq2sites{$nonamer1}->{NN_sites} += $NN_sites_nonamer1_site1;
							$seq2sites{$nonamer1}->{SN_sites} += $SN_sites_nonamer1_site1;
							$seq2sites{$nonamer1}->{NS_sites} += $NS_sites_nonamer1_site1;
							$seq2sites{$nonamer1}->{SS_sites} += $SS_sites_nonamer1_site1;
						}
					}
					
					my $NN_sites_nonamer2_site1 = 'NA';
					my $SN_sites_nonamer2_site1 = 'NA';
					my $NS_sites_nonamer2_site1 = 'NA';
					my $SS_sites_nonamer2_site1 = 'NA';
					
					if($nonamer2_num_changes_poss_site1 > 0) {
						$NN_sites_nonamer2_site1 = $nonamer2_num_changes_NN_site1 / $nonamer2_num_changes_poss_site1;
						$SN_sites_nonamer2_site1 = $nonamer2_num_changes_SN_site1 / $nonamer2_num_changes_poss_site1;
						$NS_sites_nonamer2_site1 = $nonamer2_num_changes_NS_site1 / $nonamer2_num_changes_poss_site1;
						$SS_sites_nonamer2_site1 = $nonamer2_num_changes_SS_site1 / $nonamer2_num_changes_poss_site1;
						
						unless(exists $seq_completed{$nonamer2}) {
							$seq2sites{$nonamer2}->{NN_sites} += $NN_sites_nonamer2_site1;
							$seq2sites{$nonamer2}->{SN_sites} += $SN_sites_nonamer2_site1;
							$seq2sites{$nonamer2}->{NS_sites} += $NS_sites_nonamer2_site1;
							$seq2sites{$nonamer2}->{SS_sites} += $SS_sites_nonamer2_site1;
						}
					}
					
					my $NN_sites_site1 = 'NA';
					my $SN_sites_site1 = 'NA';
					my $NS_sites_site1 = 'NA';
					my $SS_sites_site1 = 'NA';
					
					if($nonamer1_num_changes_poss_site1 > 0 && $nonamer2_num_changes_poss_site1 > 0) {
						$NN_sites_site1 = ($NN_sites_nonamer1_site1 + $NN_sites_nonamer2_site1) / 2;
						$SN_sites_site1 = ($SN_sites_nonamer1_site1 + $SN_sites_nonamer2_site1) / 2;
						$NS_sites_site1 = ($NS_sites_nonamer1_site1 + $NS_sites_nonamer2_site1) / 2;
						$SS_sites_site1 = ($SS_sites_nonamer1_site1 + $SS_sites_nonamer2_site1) / 2;
					} elsif($nonamer1_num_changes_poss_site1 > 0) {
						$NN_sites_site1 = $NN_sites_nonamer1_site1;
						$SN_sites_site1 = $SN_sites_nonamer1_site1;
						$NS_sites_site1 = $NS_sites_nonamer1_site1;
						$SS_sites_site1 = $SS_sites_nonamer1_site1;
					} elsif($nonamer2_num_changes_poss_site1 > 0) {
						$NN_sites_site1 = $NN_sites_nonamer2_site1;
						$SN_sites_site1 = $SN_sites_nonamer2_site1;
						$NS_sites_site1 = $NS_sites_nonamer2_site1;
						$SS_sites_site1 = $SS_sites_nonamer2_site1;
						
					} # else nothing defined, stay NA
					
					
					# TALLY SITE 2
					my $NN_sites_nonamer1_site2 = 'NA';
					my $SN_sites_nonamer1_site2 = 'NA';
					my $NS_sites_nonamer1_site2 = 'NA';
					my $SS_sites_nonamer1_site2 = 'NA';
					
					if($nonamer1_num_changes_poss_site2 > 0) {
						$NN_sites_nonamer1_site2 = $nonamer1_num_changes_NN_site2 / $nonamer1_num_changes_poss_site2;
						$SN_sites_nonamer1_site2 = $nonamer1_num_changes_SN_site2 / $nonamer1_num_changes_poss_site2;
						$NS_sites_nonamer1_site2 = $nonamer1_num_changes_NS_site2 / $nonamer1_num_changes_poss_site2;
						$SS_sites_nonamer1_site2 = $nonamer1_num_changes_SS_site2 / $nonamer1_num_changes_poss_site2;
						
						unless(exists $seq_completed{$nonamer1}) {
							$seq2sites{$nonamer1}->{NN_sites} += $NN_sites_nonamer1_site2;
							$seq2sites{$nonamer1}->{SN_sites} += $SN_sites_nonamer1_site2;
							$seq2sites{$nonamer1}->{NS_sites} += $NS_sites_nonamer1_site2;
							$seq2sites{$nonamer1}->{SS_sites} += $SS_sites_nonamer1_site2;
						}
					}
					
					my $NN_sites_nonamer2_site2 = 'NA';
					my $SN_sites_nonamer2_site2 = 'NA';
					my $NS_sites_nonamer2_site2 = 'NA';
					my $SS_sites_nonamer2_site2 = 'NA';
					
					if($nonamer2_num_changes_poss_site2 > 0) {
						$NN_sites_nonamer2_site2 = $nonamer2_num_changes_NN_site2 / $nonamer2_num_changes_poss_site2;
						$SN_sites_nonamer2_site2 = $nonamer2_num_changes_SN_site2 / $nonamer2_num_changes_poss_site2;
						$NS_sites_nonamer2_site2 = $nonamer2_num_changes_NS_site2 / $nonamer2_num_changes_poss_site2;
						$SS_sites_nonamer2_site2 = $nonamer2_num_changes_SS_site2 / $nonamer2_num_changes_poss_site2;
						
						unless(exists $seq_completed{$nonamer2}) {
							$seq2sites{$nonamer2}->{NN_sites} += $NN_sites_nonamer2_site2;
							$seq2sites{$nonamer2}->{SN_sites} += $SN_sites_nonamer2_site2;
							$seq2sites{$nonamer2}->{NS_sites} += $NS_sites_nonamer2_site2;
							$seq2sites{$nonamer2}->{SS_sites} += $SS_sites_nonamer2_site2;
						}
					}
					
					my $NN_sites_site2 = 'NA';
					my $SN_sites_site2 = 'NA';
					my $NS_sites_site2 = 'NA';
					my $SS_sites_site2 = 'NA';
					
					if($nonamer1_num_changes_poss_site2 > 0 && $nonamer2_num_changes_poss_site2 > 0) {
						$NN_sites_site2 = ($NN_sites_nonamer1_site2 + $NN_sites_nonamer2_site2) / 2;
						$SN_sites_site2 = ($SN_sites_nonamer1_site2 + $SN_sites_nonamer2_site2) / 2;
						$NS_sites_site2 = ($NS_sites_nonamer1_site2 + $NS_sites_nonamer2_site2) / 2;
						$SS_sites_site2 = ($SS_sites_nonamer1_site2 + $SS_sites_nonamer2_site2) / 2;
					} elsif($nonamer1_num_changes_poss_site2 > 0) {
						$NN_sites_site2 = $NN_sites_nonamer1_site2;
						$SN_sites_site2 = $SN_sites_nonamer1_site2;
						$NS_sites_site2 = $NS_sites_nonamer1_site2;
						$SS_sites_site2 = $SS_sites_nonamer1_site2;
					} elsif($nonamer2_num_changes_poss_site2 > 0) {
						$NN_sites_site2 = $NN_sites_nonamer2_site2;
						$SN_sites_site2 = $SN_sites_nonamer2_site2;
						$NS_sites_site2 = $NS_sites_nonamer2_site2;
						$SS_sites_site2 = $SS_sites_nonamer2_site2;
						
					} # else nothing defined, stay NA
					
					
					# TALLY SITE 3
					my $NN_sites_nonamer1_site3 = 'NA';
					my $SN_sites_nonamer1_site3 = 'NA';
					my $NS_sites_nonamer1_site3 = 'NA';
					my $SS_sites_nonamer1_site3 = 'NA';
					
					if($nonamer1_num_changes_poss_site3 > 0) {
						$NN_sites_nonamer1_site3 = $nonamer1_num_changes_NN_site3 / $nonamer1_num_changes_poss_site3;
						$SN_sites_nonamer1_site3 = $nonamer1_num_changes_SN_site3 / $nonamer1_num_changes_poss_site3;
						$NS_sites_nonamer1_site3 = $nonamer1_num_changes_NS_site3 / $nonamer1_num_changes_poss_site3;
						$SS_sites_nonamer1_site3 = $nonamer1_num_changes_SS_site3 / $nonamer1_num_changes_poss_site3;
						
						unless(exists $seq_completed{$nonamer1}) {
							$seq2sites{$nonamer1}->{NN_sites} += $NN_sites_nonamer1_site3;
							$seq2sites{$nonamer1}->{SN_sites} += $SN_sites_nonamer1_site3;
							$seq2sites{$nonamer1}->{NS_sites} += $NS_sites_nonamer1_site3;
							$seq2sites{$nonamer1}->{SS_sites} += $SS_sites_nonamer1_site3;
						}
					}
					
					my $NN_sites_nonamer2_site3 = 'NA';
					my $SN_sites_nonamer2_site3 = 'NA';
					my $NS_sites_nonamer2_site3 = 'NA';
					my $SS_sites_nonamer2_site3 = 'NA';
					
					if($nonamer2_num_changes_poss_site3 > 0) {
						$NN_sites_nonamer2_site3 = $nonamer2_num_changes_NN_site3 / $nonamer2_num_changes_poss_site3;
						$SN_sites_nonamer2_site3 = $nonamer2_num_changes_SN_site3 / $nonamer2_num_changes_poss_site3;
						$NS_sites_nonamer2_site3 = $nonamer2_num_changes_NS_site3 / $nonamer2_num_changes_poss_site3;
						$SS_sites_nonamer2_site3 = $nonamer2_num_changes_SS_site3 / $nonamer2_num_changes_poss_site3;
						
						unless(exists $seq_completed{$nonamer2}) {
							$seq2sites{$nonamer2}->{NN_sites} += $NN_sites_nonamer2_site3;
							$seq2sites{$nonamer2}->{SN_sites} += $SN_sites_nonamer2_site3;
							$seq2sites{$nonamer2}->{NS_sites} += $NS_sites_nonamer2_site3;
							$seq2sites{$nonamer2}->{SS_sites} += $SS_sites_nonamer2_site3;
						}
					}
					
					my $NN_sites_site3 = 'NA';
					my $SN_sites_site3 = 'NA';
					my $NS_sites_site3 = 'NA';
					my $SS_sites_site3 = 'NA';
					
					if($nonamer1_num_changes_poss_site3 > 0 && $nonamer2_num_changes_poss_site3 > 0) {
						$NN_sites_site3 = ($NN_sites_nonamer1_site3 + $NN_sites_nonamer2_site3) / 2;
						$SN_sites_site3 = ($SN_sites_nonamer1_site3 + $SN_sites_nonamer2_site3) / 2;
						$NS_sites_site3 = ($NS_sites_nonamer1_site3 + $NS_sites_nonamer2_site3) / 2;
						$SS_sites_site3 = ($SS_sites_nonamer1_site3 + $SS_sites_nonamer2_site3) / 2;
					} elsif($nonamer1_num_changes_poss_site3 > 0) {
						$NN_sites_site3 = $NN_sites_nonamer1_site3;
						$SN_sites_site3 = $SN_sites_nonamer1_site3;
						$NS_sites_site3 = $NS_sites_nonamer1_site3;
						$SS_sites_site3 = $SS_sites_nonamer1_site3;
					} elsif($nonamer2_num_changes_poss_site3 > 0) {
						$NN_sites_site3 = $NN_sites_nonamer2_site3;
						$SN_sites_site3 = $SN_sites_nonamer2_site3;
						$NS_sites_site3 = $NS_sites_nonamer2_site3;
						$SS_sites_site3 = $SS_sites_nonamer2_site3;
						
					} # else nothing defined, stay NA
					
					# SUM THE THREE SITES
					my $NN_sites = $NN_sites_site1 + $NN_sites_site2 + $NN_sites_site3;
					my $SN_sites = $SN_sites_site1 + $SN_sites_site2 + $SN_sites_site3;
					my $NS_sites = $NS_sites_site1 + $NS_sites_site2 + $NS_sites_site3;
					my $SS_sites = $SS_sites_site1 + $SS_sites_site2 + $SS_sites_site3;
					
					# Check if there are multiple variants in these overlapping codons
					my $MNV = 'FALSE';
					if(($NN_diffs + $SN_diffs + $NS_diffs + $SS_diffs) > 1) {
						$MNV = 'TRUE';
						$this_codon_MNV = 'TRUE'; # only needs to be TRUE in one comparison
					}
					
					# NEW NAIVE APPROACH: ADD TO SUMS
					my $comparison_weight = 0; # number of pairwise comparisons involving these nonamers
					if($nonamer1_index == $nonamer2_index) { # if($nonamer1 eq $nonamer2) { # it's a self-comparison
						# Combination
						$comparison_weight = (($nonamer1_count * $nonamer1_count) - $nonamer1_count) / 2; # nonamer 1 and 2 are the same
					} else {
						# m * n
						my $nonamer2_count = $unique_9mers{$seq_site_index}->{$nonamer2};
						$comparison_weight = $nonamer1_count * $nonamer2_count;
					}
					
					$comparisons_sum += $comparison_weight;
					
					$NN_sites_numerator += $NN_sites * $comparison_weight;
					$SN_sites_numerator += $SN_sites * $comparison_weight;
					$NS_sites_numerator += $NS_sites * $comparison_weight;
					$SS_sites_numerator += $SS_sites * $comparison_weight;
					
					$NN_diffs_numerator += $NN_diffs * $comparison_weight;
					$SN_diffs_numerator += $SN_diffs * $comparison_weight;
					$NS_diffs_numerator += $NS_diffs * $comparison_weight;
					$SS_diffs_numerator += $SS_diffs * $comparison_weight;
					
				} # end sas13 
			} # end case in which each member's ORF1 codon is FULLY DEFINED
				
#			} # last codon?
			
			# Note that the member members have been completed in at least one comparison
			$seq_completed{$nonamer1} = 1;
			$seq_completed{$nonamer2} = 1;
			
		} # last pair (member2)
	} # last pair (member1)
	
	if($comparisons_sum > 0) {
		# Store sites
		$codon_data_hh{$codon_in_seq}->{NN_sites} += $NN_sites_numerator / $comparisons_sum;
		$codon_data_hh{$codon_in_seq}->{SN_sites} += $SN_sites_numerator / $comparisons_sum;
		$codon_data_hh{$codon_in_seq}->{NS_sites} += $NS_sites_numerator / $comparisons_sum;
		$codon_data_hh{$codon_in_seq}->{SS_sites} += $SS_sites_numerator / $comparisons_sum;
		
		# Store diffs
		$codon_data_hh{$codon_in_seq}->{NN_diffs} += $NN_diffs_numerator / $comparisons_sum;
		$codon_data_hh{$codon_in_seq}->{SN_diffs} += $SN_diffs_numerator / $comparisons_sum;
		$codon_data_hh{$codon_in_seq}->{NS_diffs} += $NS_diffs_numerator / $comparisons_sum;
		$codon_data_hh{$codon_in_seq}->{SS_diffs} += $SS_diffs_numerator / $comparisons_sum;
	
	} else {
		# Store sites
		$codon_data_hh{$codon_in_seq}->{NN_sites} = 'NA';
		$codon_data_hh{$codon_in_seq}->{SN_sites} = 'NA';
		$codon_data_hh{$codon_in_seq}->{NS_sites} = 'NA';
		$codon_data_hh{$codon_in_seq}->{SS_sites} = 'NA';
		
		# Store diffs
		$codon_data_hh{$codon_in_seq}->{NN_diffs} = 'NA';
		$codon_data_hh{$codon_in_seq}->{SN_diffs} = 'NA';
		$codon_data_hh{$codon_in_seq}->{NS_diffs} = 'NA';
		$codon_data_hh{$codon_in_seq}->{SS_diffs} = 'NA';
	}
	
	$codon_data_hh{$codon_in_seq}->{MNV} = $this_codon_MNV;
	
	# Store majority ref codon
	my $ref_codon_maj = '';
	my $ref_codon_maj_count = 0;
	
	foreach my $ref_codon (keys %ref_codons_to_counts) {
		my $curr_ref_codon_count = $ref_codons_to_counts{$ref_codon};
		
		if($curr_ref_codon_count > $ref_codon_maj_count) {
			$ref_codon_maj = $ref_codon;
			$ref_codon_maj_count = $curr_ref_codon_count;
		}
	}
	
	$codon_data_hh{$codon_in_seq}->{ref_codon_maj} = $ref_codon_maj;
	$codon_data_hh{$codon_in_seq}->{ref_codon_maj_count} = $ref_codon_maj_count;
	
	# Store majority alt codon 1
	my $alt1_codon_maj = '';
	my $alt1_codon_maj_count = 0;
	
	foreach my $alt1_codon (keys %alt1_codons_to_counts) {
		my $curr_alt1_codon_count = $alt1_codons_to_counts{$alt1_codon};
		
		if($curr_alt1_codon_count > $alt1_codon_maj_count) {
			$alt1_codon_maj = $alt1_codon;
			$alt1_codon_maj_count = $curr_alt1_codon_count;
		}
	}
	
	$codon_data_hh{$codon_in_seq}->{alt1_codon_maj} = $alt1_codon_maj;
	$codon_data_hh{$codon_in_seq}->{alt1_codon_maj_count} = $alt1_codon_maj_count;
	
	# Store majority alt codon 2
	if($frame ne 'sas13') { 
		my $alt2_codon_maj = '';
		my $alt2_codon_maj_count = 0;
		
		foreach my $alt2_codon (keys %alt2_codons_to_counts) {
			my $curr_alt2_codon_count = $alt2_codons_to_counts{$alt2_codon};
			
			if($curr_alt2_codon_count > $alt2_codon_maj_count) {
				$alt2_codon_maj = $alt2_codon;
				$alt2_codon_maj_count = $curr_alt2_codon_count;
			}
		}
		
		$codon_data_hh{$codon_in_seq}->{alt2_codon_maj} = $alt2_codon_maj;
		$codon_data_hh{$codon_in_seq}->{alt2_codon_maj_count} = $alt2_codon_maj_count;
		
	} else { # no codon 2 is sas13
		$codon_data_hh{$codon_in_seq}->{alt2_codon_maj} = 'NA';
		$codon_data_hh{$codon_in_seq}->{alt2_codon_maj_count} = 'NA';
	}
	
} # last nonamer


##########################################################################################
# PRINT CODON RESULTS
open(CODON_FILE,">>$output_file");
print CODON_FILE "codon_num\t" .
	"ref_codon_maj\t" .
	"alt_codon1_maj\t" .
	"alt_codon2_maj\t" .
	"nonamers\t" . 
	"nonamer_counts\t" .
	"multiple_variants\t" .
	"NN_sites\tSN_sites\tNS_sites\tSS_sites\tNN_diffs\tSN_diffs\tNS_diffs\tSS_diffs\n"; # either codon might contain N or -

my $NN_sites_total = 0;
my $SN_sites_total = 0;
my $NS_sites_total = 0;
my $SS_sites_total = 0;

my $NN_diffs_total = 0;
my $SN_diffs_total = 0;
my $NS_diffs_total = 0;
my $SS_diffs_total = 0;

my $codon_output_buffer = '';

foreach my $codon_num (sort {$a <=> $b} keys %codon_data_hh) {
	
	# Remove trailing colons (:)
	chop($codon_data_hh{$codon_num}->{nonamers});
	chop($codon_data_hh{$codon_num}->{nonamer_counts});
	
	$codon_output_buffer .= "$codon_num\t" .
		$codon_data_hh{$codon_num}->{ref_codon_maj} . "\t" .
		$codon_data_hh{$codon_num}->{alt1_codon_maj} . "\t" .
		$codon_data_hh{$codon_num}->{alt2_codon_maj} . "\t";
	
	if($verbose) {
		$codon_output_buffer .= $codon_data_hh{$codon_num}->{nonamers} . "\t" .
			$codon_data_hh{$codon_num}->{nonamer_counts} . "\t";
	}
	
	$codon_output_buffer .= $codon_data_hh{$codon_num}->{MNV} . "\t" .
		$codon_data_hh{$codon_num}->{NN_sites} . "\t" .
		$codon_data_hh{$codon_num}->{SN_sites} . "\t" .
		$codon_data_hh{$codon_num}->{NS_sites} . "\t" .
		$codon_data_hh{$codon_num}->{SS_sites} . "\t" .
		$codon_data_hh{$codon_num}->{NN_diffs} . "\t" .
		$codon_data_hh{$codon_num}->{SN_diffs} . "\t" .
		$codon_data_hh{$codon_num}->{NS_diffs} . "\t" .
		$codon_data_hh{$codon_num}->{SS_diffs} . "\n";
	
	if(length($codon_output_buffer) > 50000) {
		print CODON_FILE "$codon_output_buffer";
		$codon_output_buffer = '';
	}
	
	$NN_sites_total += $codon_data_hh{$codon_num}->{NN_sites};
	$SN_sites_total += $codon_data_hh{$codon_num}->{SN_sites};
	$NS_sites_total += $codon_data_hh{$codon_num}->{NS_sites};
	$SS_sites_total += $codon_data_hh{$codon_num}->{SS_sites};
	$NN_diffs_total += $codon_data_hh{$codon_num}->{NN_diffs};
	$SN_diffs_total += $codon_data_hh{$codon_num}->{SN_diffs};
	$NS_diffs_total += $codon_data_hh{$codon_num}->{NS_diffs};
	$SS_diffs_total += $codon_data_hh{$codon_num}->{SS_diffs};
	
}

# CLEAR BUFFER
print CODON_FILE "$codon_output_buffer";
$codon_output_buffer = '';
close CODON_FILE;

print "\n";


#########################################################################################
# PRINT TOTALS TO SCREEN
my $dNN;
my $dSN;
my $dNS;
my $dSS;

if($NN_sites_total > 0) {
	$dNN = $NN_diffs_total / $NN_sites_total;
} else {
	$dNN = 'NaN';
}

if($SN_sites_total > 0) {
	$dSN = $SN_diffs_total / $SN_sites_total;
} else {
	$dSN = 'NaN';
}

if($NS_sites_total > 0) {
	$dNS = $NS_diffs_total / $NS_sites_total;
} else {
	$dNS = 'NaN';
}

if($SS_sites_total > 0) {
	$dSS = $SS_diffs_total / $SS_sites_total;
} else {
	$dSS = 'NaN';
}

my $dNNdSN;
my $dNSdSS;
my $dNNdNS;
my $dSNdSS;

if($dSN > 0 && $dNN ne 'NaN') {
	$dNNdSN = $dNN / $dSN;
} else {
	$dNNdSN = 'NaN';
}

if($dSS > 0 && $dNS ne 'NaN') {
	$dNSdSS = $dNS / $dSS;
} else {
	$dNSdSS = 'NaN';
}

if($dNS > 0 && $dNN ne 'NaN') {
	$dNNdNS = $dNN / $dNS;
} else {
	$dNNdNS = 'NaN';
}

if($dSS > 0 && $dSN ne 'NaN') {
	$dSNdSS = $dSN / $dSS;
} else {
	$dSNdSS = 'NaN';
}


print "\n#########################################################################################\n";
print "##################################### DATASET TOTALS ####################################\n";
print "#########################################################################################\n\n";

# Used sprintf instead of printf because printf changed 'NaN' to 0.

print "(1) Mean numbers of sites and differences:\n";

unless($NN_sites_total eq 'NaN') {
	$NN_sites_total = sprintf("%.3f", $NN_sites_total);
}

unless($SN_sites_total eq 'NaN') {
	$SN_sites_total = sprintf("%.3f", $SN_sites_total);
}

unless($NS_sites_total eq 'NaN') {
	$NS_sites_total = sprintf("%.3f", $NS_sites_total);
}

unless($SS_sites_total eq 'NaN') {
	$SS_sites_total = sprintf("%.3f", $SS_sites_total);
}

unless($NN_diffs_total eq 'NaN') {
	$NN_diffs_total = sprintf("%.3f", $NN_diffs_total);
}

unless($SN_diffs_total eq 'NaN') {
	$SN_diffs_total = sprintf("%.3f", $SN_diffs_total);
}

unless($NS_diffs_total eq 'NaN') {
	$NS_diffs_total = sprintf("%.3f", $NS_diffs_total);
}

unless($SS_diffs_total eq 'NaN') {
	$SS_diffs_total = sprintf("%.3f", $SS_diffs_total);
}

print "NN_sites=$NN_sites_total\n";
print "SN_sites=$SN_sites_total\n";
print "NS_sites=$NS_sites_total\n";
print "SS_sites=$SS_sites_total\n";
print "NN_diffs=$NN_diffs_total\n";
print "SN_diffs=$SN_diffs_total\n";
print "NS_diffs=$NS_diffs_total\n";
print "SS_diffs=$SS_diffs_total\n";
print "\n";
#printf "NN_sites=%.3f\n", $NN_sites_total;
#printf "SN_sites=%.3f\n", $SN_sites_total;
#printf "NS_sites=%.3f\n", $NS_sites_total;
#printf "SS_sites=%.3f\n", $SS_sites_total;
#printf "NN_diffs=%.3f\n", $NN_diffs_total;
#printf "SN_diffs=%.3f\n", $SN_diffs_total;
#printf "NS_diffs=%.3f\n", $NS_diffs_total;
#printf "SS_diffs=%.3f\n", $SS_diffs_total;
#print "\n";

print "(2) Mean substitution rates (between-species) or nucleotide diversities (within-species):\n";

unless($dNN eq 'NaN') {
	$dNN = sprintf("%.6f", $dNN);
}

unless($dSN eq 'NaN') {
	$dSN = sprintf("%.6f", $dSN);
}

unless($dNS eq 'NaN') {
	$dNS = sprintf("%.6f", $dNS);
}

unless($dSS eq 'NaN') {
	$dSS = sprintf("%.6f", $dSS);
}

print "dNN=$dNN\n";
print "dSN=$dSN\n";
print "dNS=$dNS\n";
print "dSS=$dSS\n";
#printf "dNN=%.3f\n", $dNN;
#printf "dSN=%.3f\n", $dSN;
#printf "dNS=%.3f\n", $dNS;
#printf "dSS=%.3f\n", $dSS;
print "\n";

print "(3) dN/dS estimates:\n";
print "For the REFERENCE gene (ORF1):\n";

unless($dNNdSN eq 'NaN') {
	$dNNdSN = sprintf("%.3f", $dNNdSN);
}

unless($dNSdSS eq 'NaN') {
	$dNSdSS = sprintf("%.3f", $dNSdSS);
}

print "dNN/dSN=$dNNdSN\n";
print "dNS/dSS=$dNSdSS\n";
#printf "dNN/dSN=%.3f\n", $dNNdSN;
#printf "dNS/dSS=%.3f\n", $dNSdSS;
print "\n";

print "For the ALTERNATE gene (ORF2):\n";

unless($dNNdNS eq 'NaN') {
	$dNNdNS = sprintf("%.3f", $dNNdNS);
}

unless($dSNdSS eq 'NaN') {
	$dSNdSS = sprintf("%.3f", $dSNdSS);
}

print "dNN/dNS=$dNNdNS\n";
print "dSN/dSS=$dSNdSS\n";
#printf "dNN/dNS=%.3f\n", $dNNdNS;
#printf "dSN/dSS=%.3f\n", $dSNdSS;
print "\n";

print "(4) NOTES:\n";
print "-Values have been rounded to 3 or 6 decimal places.\n";
print "-NaN indicates undefined (Not a Number) values (division by 0).\n";
print "-To test for significance, please bootstrap OLGenie_codon_results.txt using OLGenie_bootstrap.R.\n";
print "\tFor example, at the command line:\n\n";
print "\$ OLGenie_bootstrap.R OLGenie_codon_results.txt > OLGenie_test_results.txt\n\n";
print "\n#########################################################################################\n\n";


##########################################################################################
# AVERAGE NUMBER OF SITES BETWEEN ALL SEQUENCES
# SITES ARE HERE:
#$seq2sites{$member2}->{NN_sites} += $NN_sites_member2_site3;
#$seq2sites{$member2}->{SN_sites} += $SN_sites_member2_site3;
#$seq2sites{$member2}->{NS_sites} += $NS_sites_member2_site3;
#$seq2sites{$member2}->{SS_sites} += $SS_sites_member2_site3;

my @NN_sites;
my @SN_sites;
my @NS_sites;
my @SS_sites;


##########################################################################################
# DIFFERENCES: RAW AND NORMALIZED (EACH SITE IS A MAXIMUM OF 1 DIFFERENCE)
my $NN_diffs_sum = 0;
my $SN_diffs_sum = 0;
my $NS_diffs_sum = 0;
my $SS_diffs_sum = 0;

# Only needed if $verbose_messages is ON
my $NN_diffs_normalized = 0;
my $SN_diffs_normalized = 0;
my $NS_diffs_normalized = 0;
my $SS_diffs_normalized = 0;

foreach my $codon_num (sort {$a <=> $b} keys %site_diffs_hh) {
	
	my $NN_diffs_total = 0;
	my $SN_diffs_total = 0;
	my $NS_diffs_total = 0;
	my $SS_diffs_total = 0;
	
	foreach my $codon_site (sort {$a <=> $b} keys %{$site_diffs_hh{$codon_num}}) {
		my $NN_diffs = $site_diffs_hh{$codon_num}->{$codon_site}->{NN_diffs};
		my $SN_diffs = $site_diffs_hh{$codon_num}->{$codon_site}->{SN_diffs};
		my $NS_diffs = $site_diffs_hh{$codon_num}->{$codon_site}->{NS_diffs};
		my $SS_diffs = $site_diffs_hh{$codon_num}->{$codon_site}->{SS_diffs};
		
		# Add up raw totals
		$NN_diffs_sum += $NN_diffs;
		$SN_diffs_sum += $SN_diffs;
		$NS_diffs_sum += $NS_diffs;
		$SS_diffs_sum += $SS_diffs;
		
		# Now, for normalization
		if($verbose_messages) {
			my $total_diffs = $NN_diffs + $SN_diffs + $NS_diffs + $SS_diffs;
			
			if($total_diffs > 1) {
				my $NN_diffs_corrected = $NN_diffs / $total_diffs;
				my $SN_diffs_corrected = $SN_diffs / $total_diffs;
				my $NS_diffs_corrected = $NS_diffs / $total_diffs;
				my $SS_diffs_corrected = $SS_diffs / $total_diffs;
				
				$NN_diffs_total += $NN_diffs_corrected;
				$SN_diffs_total += $SN_diffs_corrected;
				$NS_diffs_total += $NS_diffs_corrected;
				$SS_diffs_total += $SS_diffs_corrected;
			} else {
				$NN_diffs_total += $NN_diffs;
				$SN_diffs_total += $SN_diffs;
				$NS_diffs_total += $NS_diffs;
				$SS_diffs_total += $SS_diffs;
			}
		}
	}
	
	if($verbose_messages) {
		$NN_diffs_normalized += $NN_diffs_total;
		$SN_diffs_normalized += $SN_diffs_total;
		$NS_diffs_normalized += $NS_diffs_total;
		$SS_diffs_normalized += $SS_diffs_total;
	}
}


#########################################################################################
# PRINT OTHER TOTALS TO SCREEN IF DEVELOPING
if($verbose_messages) {	
	
	my $NN_sites_mean = &mean(@NN_sites);
	my $SN_sites_mean = &mean(@SN_sites);
	my $NS_sites_mean = &mean(@NS_sites);
	my $SS_sites_mean = &mean(@SS_sites);
	
	# LATER, make this output LOOK PRETTY
	print "\n#########################################################################################\n";
	print "### DATA FOR HYPOTHESIS TESTING:\n\n";
	
	# Note these differences may differ from simple sums of pairs IF THERE ARE MULTIPLE VARIANTS IN A CODON
	print "#########################################################################################\n";
	print "# RAW NUMBER OF DIFFERENCES:\n\n"; # \n# (Recommended for maximum power if alignment is certain.)
	print "\tNN\tSN\tNS\tSS\n";
	print "Diffs\t$NN_diffs_sum\t$SN_diffs_sum\t$NS_diffs_sum\t$SS_diffs_sum\n";
	print "Sites-Diffs\t" . ($NN_sites_mean - $NN_diffs_sum) . "\t" .
		($SN_sites_mean - $SN_diffs_sum) . "\t" . 
		($NS_sites_mean - $NS_diffs_sum) . "\t" .
		($SS_sites_mean - $SS_diffs_sum) . "\n";
	print "SITES\t$NN_sites_mean\t$SN_sites_mean\t$NS_sites_mean\t$SS_sites_mean\n";
	
	
	print "\n#########################################################################################\n";
	print "# NORMALIZED NUMBER OF DIFFERENCES (MAX 1 PER SITE):\n# (Recommended if alignment may contain errors or large numbers of differences are observed.)\n\n";
	print "\tNN\tSN\tNS\tSS\n";
	print "Diffs\t$NN_diffs_normalized\t$SN_diffs_normalized\t$NS_diffs_normalized\t$SS_diffs_normalized\n";
	print "Sites-Diffs\t" . ($NN_sites_mean - $NN_diffs_normalized) . "\t" .
		($SN_sites_mean - $SN_diffs_normalized) . "\t" . 
		($NS_sites_mean - $NS_diffs_normalized) . "\t" .
		($SS_sites_mean - $SS_diffs_normalized) . "\n";
	print "SITES\t$NN_sites_mean\t$SN_sites_mean\t$NS_sites_mean\t$SS_sites_mean\n";

	print "\n#########################################################################################\n";

	#####################################################################################
	# Print the above to file for easy reading into R
	open(OLGENIE_FILE, ">>OLGenie\_results\.txt");
	
	print OLGENIE_FILE "measure\t" .
		"NN\tSN\tNS\tSS\n";
	
	print OLGENIE_FILE "diffs_raw\t$NN_diffs_sum\t$SN_diffs_sum\t$NS_diffs_sum\t$SS_diffs_sum\n" .
		"diffs_normalized\t$NN_diffs_normalized\t$SN_diffs_normalized\t$NS_diffs_normalized\t$SS_diffs_normalized\n" .
		"sites\t$NN_sites_mean\t$SN_sites_mean\t$NS_sites_mean\t$SS_sites_mean\n";
	
	close OLGENIE_FILE;
}


#########################################################################################
# Print a completion message to screen
&end_the_program;
exit;


#########################################################################################
#########################################################################################
####################################                 ####################################
####################################   SUBROUTINES   ####################################
####################################                 ####################################
#########################################################################################
#########################################################################################

#########################################################################################
sub median {
    my @values = sort {$a <=> $b} @_;
    my $length = @values;
    if($length%2) { # odd number of elements: return middle
        return $values[int($length/2)];
    } else { # even number of elements: return mean of middle two
        return ($values[int($length/2)-1] + $values[int($length/2)])/2;
    }
}

#########################################################################################
sub mean {
    my @values = @_;
    my $length = @values;
    my $sum;
    
    foreach (@values) {
    	$sum += $_;
    }
    
    if($length > 0) {
 	   return($sum / $length);
    } else {
    	return('NA');
    }
}

#########################################################################################
sub standard_deviation {
    my @values = @_;
    my $length = @values;
	my $mean_of_values = &mean(@values);
	my $sum_squared_deviations;
	
	foreach (@values) {
		$sum_squared_deviations += ($_ - $mean_of_values)**2;
	}
	
	my $variance = ($sum_squared_deviations) / ($length - 1);
	
	return(sqrt($variance));
	
}

#########################################################################################
sub revcom {
    my ($seq) = @_;
	chomp($seq);
	$seq = uc($seq); # uc returns uppercase
	$seq =~ tr/U/T/; # no RNA
	
	my $rev_seq = reverse($seq);
	my $rev_com_seq = $rev_seq;
	$rev_com_seq =~ tr/ACGT/TGCA/;
	
	return($rev_com_seq);
}


#########################################################################################
# Get the amino acid (single-letter code) encoded by a given DNA or RNA codon
# Returns an array with:
#	returned[0] = number of nonsynonymous sites
#	returned[1] = number of synonymous sites
sub get_amino_acid {
	my ($codon) = @_;
	my $amino_acid;
	
	# Establish genetic code for use with synonymous sites; DNA or RNA
	my %code = (
		"AAA" => "K", "AAC" => "N", "AAG" => "K", "AAT" => "N", "AAU" => "N", 
		"ACA" => "T", "ACC" => "T", "ACG" => "T", "ACT" => "T", "ACU" => "T", 
		"AGA" => "R", "AGC" => "S", "AGG" => "R", "AGT" => "S", "AGU" => "S", 
		"ATA" => "I", "ATC" => "I", "ATG" => "M", "ATT" => "I", "AUA" => "I", 
		"AUC" => "I", "AUG" => "M", "AUU" => "I", "CAA" => "Q", "CAC" => "H", 
		"CAG" => "Q", "CAT" => "H", "CAU" => "H", "CCA" => "P", "CCC" => "P", 
		"CCG" => "P", "CCT" => "P", "CCU" => "P", "CGA" => "R", "CGC" => "R", 
		"CGG" => "R", "CGT" => "R", "CGU" => "R", "CTA" => "L", "CTC" => "L", 
		"CTG" => "L", "CTT" => "L", "CUA" => "L", "CUC" => "L", "CUG" => "L", 
		"CUU" => "L", "GAA" => "E", "GAC" => "D", "GAG" => "E", "GAT" => "D", 
		"GAU" => "D", "GCA" => "A", "GCC" => "A", "GCG" => "A", "GCT" => "A", 
		"GCU" => "A", "GGA" => "G", "GGC" => "G", "GGG" => "G", "GGT" => "G", 
		"GGU" => "G", "GTA" => "V", "GTC" => "V", "GTG" => "V", "GTT" => "V", 
		"GUA" => "V", "GUC" => "V", "GUG" => "V", "GUU" => "V", "TAA" => "*", 
		"TAC" => "Y", "TAG" => "*", "TAT" => "Y", "UAA" => "*", "UAC" => "Y", 
		"UAG" => "*", "UAU" => "Y", "TCA" => "S", "TCC" => "S", "TCG" => "S", 
		"TCT" => "S", "UCA" => "S", "UCC" => "S", "UCG" => "S", "UCU" => "S", 
		"TGA" => "*", "TGC" => "C", "TGG" => "W", "TGT" => "C", "UGA" => "*", 
		"UGC" => "C", "UGG" => "W", "UGU" => "C", "TTA" => "L", "TTC" => "F", 
		"TTG" => "L", "TTT" => "F", "UUA" => "L", "UUC" => "F", "UUG" => "L", 
		"UUU" => "F"
	);
	
	$amino_acid = $code{$codon};
	
	if($amino_acid eq '') {
		$amino_acid = '?';
	}
	
	return $amino_acid;
}


#########################################################################################
# End the program by notifying the screen at command line
sub end_the_program {
	my $time2 = time;
	my $local_time2 = localtime;
	
	my $time_diff = ($time2 - $time1);
	my $time_diff_rounded = sprintf("%.2f",$time_diff);
	my $mins_elapsed = ($time_diff / 60);
	my $whole_mins_elapsed = int($mins_elapsed);
	my $whole_mins_in_secs = ($whole_mins_elapsed * 60);
	my $secs_remaining = ($time_diff - $whole_mins_in_secs);
	my $secs_remaining_rounded = sprintf("%.2f",$secs_remaining);
	
	print "OLGenie completed at local time $local_time2. The process took $time_diff_rounded secs, i.e., ".
			"$whole_mins_elapsed mins and $secs_remaining_rounded secs\n";

	print "\n################################################################################".
		"\n##                      OLGenie completed successfully.                       ##".
		"\n##                Please find results in the working directory.               ##\n".
		"################################################################################".
		"\n\n\n"; 
}
