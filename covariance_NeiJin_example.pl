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
#use warnings; 
#use diagnostics;
#use autodie;
use Parallel::ForkManager;
use List::Util qw( min max );

# SET UP EXAMPLE
my $N = 3;
my $m = 120;

my %dNN_hh;

$dNN_hh{1}->{2}->{d} = 0.0009;
$dNN_hh{1}->{2}->{var_d} = 0.00009;

$dNN_hh{1}->{3}->{d} = 0.001;
$dNN_hh{1}->{3}->{var_d} = 0.0001;

$dNN_hh{2}->{3}->{d} = 0.0011;
$dNN_hh{2}->{3}->{var_d} = 0.00011;

my $print_flag = 0;

#my ($new_pi, $new_var_pi) = &var_pairwise_distance_NeiJin($N, $m, \%dNN_hh);
#
#print "\nThe results are: pi=$new_pi\, var_pi=$new_var_pi\n\n";

# Example fron Nei & Jin 1989, Table 1, common chimps only
my %common_chimp;
$common_chimp{1}->{2}->{d} = .0168;
$common_chimp{1}->{3}->{d} = .0220;
$common_chimp{1}->{4}->{d} = .0141;
$common_chimp{1}->{5}->{d} = .0220;

$common_chimp{2}->{3}->{d} = .0048;
$common_chimp{2}->{4}->{d} = .0173;
$common_chimp{2}->{5}->{d} = .0203;

$common_chimp{3}->{4}->{d} = .0122;
$common_chimp{3}->{5}->{d} = .015;

$common_chimp{4}->{5}->{d} = .0072;

$common_chimp{1}->{2}->{mij} = 33;
$common_chimp{1}->{3}->{mij} = 32;
$common_chimp{1}->{4}->{mij} = 34;
$common_chimp{1}->{5}->{mij} = 32;

$common_chimp{2}->{3}->{mij} = 34;
$common_chimp{2}->{4}->{mij} = 32;
$common_chimp{2}->{5}->{mij} = 31;

$common_chimp{3}->{4}->{mij} = 33;
$common_chimp{3}->{5}->{mij} = 32;

$common_chimp{4}->{5}->{mij} = 34;


for (my $i = 1; $i <= 5; $i++) {
	for (my $j = $i+1; $j <= 5; $j++) {
		my $this_d = $common_chimp{$i}->{$j}->{d};
		my $this_m = $common_chimp{$i}->{$j}->{mij};
		
		
#		my $p = (3/4) * (1 - exp(-4/3 * $b));
		my $p = (3/4) * (1 - exp(-4/3 * $this_d)); # SHOULD BE b # COMEBACK
#		my $this_p = (3/4) * (1 - exp(-4 / 3 * $this_d));
		
#		$p = $this_p;
		
		my $this_var_d = 9 * $p * (1 - $p) / ((3 - 4 * $p)**2 * $this_m);
		
		$common_chimp{$i}->{$j}->{var_d} = $this_var_d;
	}
}


#&UPGMA_from_distance_hash(\%common_chimp);

my $common_chimp_distmatrix_ref = &convert_dhh_to_distmatrix(\%common_chimp);
my %common_chimp_distmatrix = %{$common_chimp_distmatrix_ref};

my $common_chimp_tree = &UPGMA_from_distance_hash(\%common_chimp_distmatrix);
# ME MATCH EXAMPLE TREE FROM NEI & JIN 1989 FIG. 3A, TABLE 1 TAXA 1-5

print "common chimp tree: $common_chimp_tree\n";

my ($chimp_pi, $chimp_var_pi) = &var_pairwise_distance_NeiJin(5, 36, \%common_chimp, $common_chimp_tree); # my ($N, $m, $d_hh_ref) # $N = num_taxa; $m = median num. nt in alignment
my $chimp_SE = sqrt($chimp_var_pi);
print "\nThe chimp RESTRICTION ENZYME ANSWERS are: pi=0.0152, var_pi=0.00001866, SE=0.0043\n";
print "The chimp NUCLEOTIDE RESULTS are: pi=$chimp_pi\, var_pi=$chimp_var_pi\, SE_pi=$chimp_SE\n\n"; # VARIANCE not right yet

#my $patristic_distance = &get_patristic_distance($common_chimp_tree, '1', '2');
#print "patristic distance 1-to-2 is: $patristic_distance\n";
#my $patristic_distance = &get_patristic_distance($common_chimp_tree, '1', '3');
#print "patristic distance 1-to-3 is: $patristic_distance\n";
#my $patristic_distance = &get_patristic_distance($common_chimp_tree, '1', '4');
#print "patristic distance 1-to-4 is: $patristic_distance\n";
#my $patristic_distance = &get_patristic_distance($common_chimp_tree, '1', '5');
#print "patristic distance 1-to-5 is: $patristic_distance\n";
#my $patristic_distance = &get_patristic_distance($common_chimp_tree, '2', '3');
#print "patristic distance 2-to-3 is: $patristic_distance\n";
#my $patristic_distance = &get_patristic_distance($common_chimp_tree, '2', '4');
#print "patristic distance 2-to-4 is: $patristic_distance\n";
#my $patristic_distance = &get_patristic_distance($common_chimp_tree, '2', '5');
#print "patristic distance 2-to-5 is: $patristic_distance\n";
#my $patristic_distance = &get_patristic_distance($common_chimp_tree, '3', '4');
#print "patristic distance 3-to-4 is: $patristic_distance\n";
#my $patristic_distance = &get_patristic_distance($common_chimp_tree, '3', '5');
#print "patristic distance 3-to-5 is: $patristic_distance\n";
#my $patristic_distance = &get_patristic_distance($common_chimp_tree, '4', '5');
#print "patristic distance 4-to-5 is: $patristic_distance\n";


# ANOTHER UPGMA EXAMPLE FROM https://www.mygoblet.org/training-portal/materials/upgma-worked-example
my %d_example;
#$d_example{'Apple_78'}->{'B'} = 19;
#$d_example{'Apple_78'}->{'C'} = 27;
#$d_example{'Apple_78'}->{'D'} = 8;
#$d_example{'Apple_78'}->{'E'} = 33;
#$d_example{'Apple_78'}->{'F'} = 18;
#$d_example{'Apple_78'}->{'G'} = 13;

$d_example{'A'}->{'B'} = 19;
$d_example{'A'}->{'C'} = 27;
$d_example{'A'}->{'D'} = 8;
$d_example{'A'}->{'E'} = 33;
$d_example{'A'}->{'F'} = 18;
$d_example{'A'}->{'G'} = 13;

$d_example{'B'}->{'C'} = 31;
$d_example{'B'}->{'D'} = 18;
$d_example{'B'}->{'E'} = 36;
$d_example{'B'}->{'F'} = 1;
$d_example{'B'}->{'G'} = 13;

$d_example{'C'}->{'D'} = 26;
$d_example{'C'}->{'E'} = 41;
$d_example{'C'}->{'F'} = 32;
$d_example{'C'}->{'G'} = 29;

$d_example{'D'}->{'E'} = 31;
$d_example{'D'}->{'F'} = 17;
$d_example{'D'}->{'G'} = 14;

$d_example{'E'}->{'F'} = 35;
$d_example{'E'}->{'G'} = 28;

$d_example{'F'}->{'G'} = 12;

#my $example_tree = &UPGMA_from_distance_hash(\%d_example);
# WE MATCH EXAMPLE FROM https://www.mygoblet.org/training-portal/materials/upgma-worked-example
#print "\nexample tree: $example_tree\n";
#
#$patristic_distance = &get_patristic_distance($example_tree, 'A', 'B');
#print "patristic distance A-to-B is: $patristic_distance\n";
#
#$patristic_distance = &get_patristic_distance($example_tree, 'A', 'G');
#print "patristic distance A-to-G is: $patristic_distance\n";
#
#$patristic_distance = &get_patristic_distance($example_tree, 'E', 'B');
#print "patristic distance B-to-E is: $patristic_distance\n";
#
#$patristic_distance = &get_patristic_distance($example_tree, 'F', 'B');
#print "patristic distance B-to-F is: $patristic_distance\n";
#
#$patristic_distance = &get_patristic_distance($example_tree, 'A', 'D');
#print "patristic distance A-to-D is: $patristic_distance\n";

exit;

#########################################################################################
# Calculate variance of dNN, dSN, dNS, and dSS according to Nei & Jin. Mol. Biol. Evol. 6(3):290-300. 1989.
sub var_pairwise_distance_NeiJin {
	my ($N, $m, $d_hh_ref, $tree) = @_; # $m = median nt in alignment, hh for dNN, dSN, dNS, or dSS
	my %d_hh = %{$d_hh_ref};
#	my %d_hh = %dNN_hh;
#	my ($N, $m, @array) = @_; # $m is the number of nucleotides per sequence
	
#	print "input hash: " . %d_hh . "\n";
#	foreach my $key1 (sort keys %d_hh) {
#		foreach my $key2 (sort keys %{$d_hh{$key1}}) {
#			print "key 1 = $key1 | key 2 = $key2\n";
#		}
#	}
	
	#---------Format of the input-----------
	# $N is the total number of sequences.
	# @array contains 4*N*N numbers: first one could be one of @all_dNN, @all_dNS, @all_dSN, @all_dSS
	# second one could be one of @all_var_dNN, @all_var_dNS, @all_var_dSN, @all_var_dSS and 
	# third and fourth are @all_indices_i and @all_indeces_j.
	# all pairs of dNN_ij are in @all_dNN, and the order of ij should be the same for @all_dNN, @all_dNS, @all_dSN, and all @dSS.  
	# @all_dNN contains results of sequence pair such that seq_index i and j are not equal. 
	
	#---------Requirement of input data-----------
	# The calculation of cov(dNN_ij,dNN_kl) requires the existence of dNN_ij, dNN_kl, 
	# dNN_il, dNN_jk, dNN_ik, dNN_jl.
	# When not all four numbers exist, covariance cannot be calculated.
	# var(pi_NN) can only be calculated when all pairs of sequence satisfy this criteria,
	# so when dNN_ij does not exist, dNN_ix, and dNN_jx should all be removed from the calculation. 
	# This subroutine only calculates var(pi_XX) and pi_XX when dXX (X could be N or S) 
	# from all pairs of sequences exist.
	
	#(1) Regarding @array, I don't understand what this contains. The comments say it contains 
	#	4 * N * N numbers. Is this a multidimensional array? If we have N sequences, then we 
	#	have N-choose-2 = N(N-1)/2 comparisons between sequences. Thus there should be 
	#	N(N-1)/2 numbers within each of @all_dNN, @all_dNS, @all_dSN, @all_dSS, @all_var_dNN, 
	#	@all_var_dNS, @all_var_dSN, and @all_var_dSS. A particular sequence pair will have values 
	#	located at the same index in each of these arrays. But I don't understand how this 
	#	relates to @array.
	#ANSWER:
	#(1) I think you are right. The @array should contain @all_dNN, @all_dNS, @all_dSN, 
	#	@all_dSS, @all_var_dNN, @all_var_dNS, @all_var_dSN, and @all_var_dSS 
	#	as columns for,  and N(N-1)/2 rows are necessary and sufficient. 
	#	I was not sure about needing the indices last time when I tried to do this. 
	#	I think last time I might think something wrong about adding the indices. 
	#	It should be enough to know N, i, j to check the thing. 
	
	#(2) It looks like $m is the number of nucleotides (not codons?) per sequence. Does this 
	#	need to be the same number for all comparisons? If so, this may be a problem, because 
	#	I have changed the actual sites present in each sequence pair by eliminating STOP 
	#	codons and gaps. Perhaps we choose the median number of sites? Not sure.
	#ANSWER:
	#(2) $m is nucleotides in the original covariance calculation.  
	#	I think using should be median is fair. The paper did not say anything about this, 
	#	but this is not a specific problem to overlapping genes but to all indels problem. 
	
	# I don't want it to be an array. I want it to be a hash. So let's do a hash.
#	my %example_hh; # keys are sequence numbers and/or other IDs
#	$example_hh{1}->{2}->{dNN};
#	$example_hh{1}->{2}->{var_dNN};
#	$example_hh{1}->{2}->{dSN};
#	$example_hh{1}->{2}->{var_dSN};
#	$example_hh{1}->{2}->{dNS};
#	$example_hh{1}->{2}->{var_dNS};
#	$example_hh{1}->{2}->{dSS};
#	$example_hh{1}->{2}->{var_dSS};
	
#	my $arrSize = @array;
	
#	my @all_dNN;
	
#	my @all_d = @array[0..$arraySize/4-1];
#	my @all_var_d = @array[$arraySize/4..$arraySize*2/4-1];
#	my @all_indices_1 = @array[$arraySize*2/4..$arraySize*3/4-1];
#	my @all_indeces_2 = @array[$arraySize*3/4..$arraySize-1];
	
	#----------Define parametors to be estimated as output----------
	my $pi_numerator = 0;
	my $pi_counter = 0;
	
#	my $var_pi = 0;
#	my $var_pi_counter = 0;
	
	my $var_pi_numerator = 0;
	my $var_pi_counter = 0;
	
	#-----------Define other parameters----------
	my $d_ij;
	my $d_kl;
	my $d_ik;
	my $d_il;
	my $d_jk;
	my $d_jl;
	
	my $cov_d_ij_kl;
#	my $idx_ij;
#	my $idx_kl;
#	my $idx_ik;
#	my $idx_il;
#	my $idx_jk;
#	my $idx_kl;
	
	#---------------Calculation---------------
#	for (my $seq_i = 1; $seq_i <= $N; $seq_i++) { # each sequence
#		for (my $seq_j = 1; $seq_j <= $N; $seq_j++) { # each sequence
	
	my %all_unique_seq_names;
	foreach my $seq_i (sort keys %d_hh) { # each sequence
		foreach my $seq_j (sort keys %{$d_hh{$seq_i}}) { # each sequence
			$all_unique_seq_names{$seq_i} = 1;
			$all_unique_seq_names{$seq_j} = 1;
		}
	}
	
	my @seq_names_sorted = sort keys %all_unique_seq_names;
	
	if($N != scalar(@seq_names_sorted)) {
		die "\n### WARNING: mismatch in number of sequences. There's a bug.\n";
	}
	
	my %store_covariances;
	
	foreach my $seq_i (@seq_names_sorted) { # each sequence
		foreach my $seq_j (@seq_names_sorted) { # each sequence
			
			if ($seq_i ne $seq_j) { # but not the same sequence
				
				print "\ntesting seq_i = $seq_i vs. seq_j = $seq_j\n";
				
				my @seq_i_j_ordered = sort ($seq_i, $seq_j);
				my $min_ij = $seq_i_j_ordered[0];
				my $max_ij = $seq_i_j_ordered[1];
				
				my $d_ij = $d_hh{$min_ij}->{$max_ij}->{d}; # checked
				my $var_d_ij = $d_hh{$min_ij}->{$max_ij}->{var_d}; # checked
				
				print "\nThe sequences are $seq_i and $seq_j and we have d_ij=$d_ij\, var_d_ij=$var_d_ij\n\n";
				
				#------------------------------Calculate pi and var_pi------------------------
				$pi_numerator += $d_ij;
				$pi_counter++;
				
#				$var_pi += $var_d_ij;
#				$var_pi_counter++;
				
				my $covariance_sum = 0;
				my $covariance_counter = 0;
				
				foreach my $seq_k (@seq_names_sorted) { # each sequence
					foreach my $seq_l (@seq_names_sorted) { # each sequence
						#if ($seq_k ne $seq_l && !($seq_i eq $seq_k && $seq_j eq $seq_l) && !($seq_i eq $seq_l && $seq_j eq $seq_k)) { # but not the same sequence for dkl
						if ($seq_k ne $seq_l) { # but not the same sequence for dkl

							my @sorted_ij = sort($seq_i, $seq_j);
							my @sorted_kl = sort($seq_k, $seq_l);
							my @sorted_ik = sort($seq_i, $seq_k);
							my @sorted_jl = sort($seq_j, $seq_l);
							my @sorted_il = sort($seq_i, $seq_l);
							my @sorted_jk = sort($seq_j, $seq_k);
							
							
							my $min_ij = $sorted_ij[0];
							my $min_kl = $sorted_kl[0];
							my $min_ik = $sorted_ik[0];
							my $min_jl = $sorted_jl[0];
							my $min_il = $sorted_il[0];
							my $min_jk = $sorted_jk[0];
							
							
							my $max_ij = $sorted_ij[1];
							my $max_kl = $sorted_kl[1];
							my $max_ik = $sorted_ik[1];
							my $max_jl = $sorted_jl[1];
							my $max_il = $sorted_il[1];
							my $max_jk = $sorted_jk[1];
							
							
#								print "\nmin of ($seq_i\, $seq_k\) is $min_ik\n" . 
#									"max of ($seq_i\, $seq_k\) is $max_ik\n\n";
							
							# GET PATRISTIC DISTANCES using UPGMA $tree, compute bij values
							my $b_ij = &get_patristic_distance($tree, $min_ij, $max_ij);
							my $b_kl = &get_patristic_distance($tree, $min_kl, $max_kl);
							my $b_ik = &get_patristic_distance($tree, $min_ik, $max_ik);
							my $b_jl = &get_patristic_distance($tree, $min_jl, $max_jl);
							my $b_il = &get_patristic_distance($tree, $min_il, $max_il);
							my $b_jk = &get_patristic_distance($tree, $min_jk, $max_jk);
							
#								print "patristic distance $seq_i to $seq_j is $b_ij\n";
							
#								print "calculating covariances for:\nseq_i=$seq_i\nseq_j=$seq_j\n" .
#									"seq_k=$seq_k\nseq_l=$seq_l\n";
							
							if($seq_i eq '5' && $seq_j eq '4' && $seq_k eq '1' && $seq_l eq '4') {
								#print "bij = b54 = $b_ij\n";
								#print "bkl = b14 = $b_kl\n";
								#print "bik = b51 = $b_ik\n";
								#print "bjl = b44 = $b_jl\n";
								#print "bil = b54 = $b_il\n";
								#print "bjk = b41 = $b_jk\n";
								$print_flag = 1;
							}
							
							# MUST PASS b's not d's
							$cov_d_ij_kl = &cov_calculation_NeiJin($b_ij, $b_kl, $b_ik, $b_jl, $b_il, $b_jk, $m);
							
							#print "\ncovariance of $seq_i / $seq_j and $seq_k / $seq_l is $cov_d_ij_kl\n";
#								print "b_jk was $b_jk\n";
#								if($seq_i eq '1' && $seq_j eq '2' && $seq_k eq '3' && $seq_l eq '4') {
#									print "b_ij=$b_ij, b_kl=$b_kl, b_ik=$b_ik, b_jl=$b_jl, b_il=$b_il, b_jk=$b_jk\n";
#								}
							
							#if($seq_i eq $seq_k && $seq_j eq $seq_l) {
							#	print "\nseq_i=$seq_i\, seq_j=$seq_j\, seq_k=$seq_k\, seq_l=$seq_l\n";
							#	print "covariance=$cov_d_ij_kl\n";
							#}
							
							if($print_flag == 1) {
								print "\nseq_i=$seq_i\, seq_j=$seq_j\, seq_k=$seq_k\, seq_l=$seq_l\n";
								print "covariance=$cov_d_ij_kl\n";
							}
							
							$print_flag = 0;
							
							my $ij_name = "d$seq_i$seq_j";
							my $kl_name = "d$seq_k$seq_l";
							
							$store_covariances{$ij_name}->{$kl_name} = $cov_d_ij_kl;
							
							$covariance_sum += $cov_d_ij_kl;
							$covariance_counter++;
							
						} # k != l
					} # finished all l seqs
				} # finished all k seqs
				
				print "\nI determined $covariance_counter covariances\n";
#				print "\ncovariance of $seq_i and $seq_j is $cov_d_ij_kl\n";
				
				$var_pi_numerator += ($var_d_ij + $covariance_sum);
#				$var_pi_numerator += $covariance_sum;
				$var_pi_counter++;
			} # i != j
		} # finished all j seqs
	} # finished all i seqs
	
	
	print "\nI determined $var_pi_counter variances\n";
	
	print "\nI determined these covariances: \n"; # . %store_covariances . "\n";#$store_covariances{$ij_name}->{$kl_name} = $cov_d_ij_kl;;
	
	foreach my $ij_name (sort keys %store_covariances) {
		foreach my $kl_name (sort keys %{$store_covariances{$ij_name}}) {
			print "covariance $ij_name and $kl_name = $store_covariances{$ij_name}->{$kl_name}\n";
			
		}
	}
	
#	print "42 vs 23 is: " . $store_covariances{d42}->{d23} * 10 . "\n";
	
	# CALCULATE VARIANCE HERE
	# $pi_counter == ($N**2 - $N) / 2 ?
	my $pi = $pi_numerator / (($N**2 - $N));
#	my $pi = $pi_numerator / $pi_counter;
	my $var_pi = $var_pi_numerator / (($N * ($N - 1)) ** 2);
	
#	print "\ndenominator is: " . (($N * ($N - 1)) ** 2) . "\n";
	
	my @result = ($pi, $var_pi);
	return @result;
}
	
# 	# TAKAHATA & TAJIMA 1991 VERSION
# 	my $num_seq_names_Takahata = scalar(@seq_names_sorted);
# 	
# 	for(my $i = 0; $i < $num_seq_names_Takahata; $i++) { # each sequence
#  		for(my $j = $i + 1; $j < $num_seq_names_Takahata; $j++) { # i<j
#  			
#  			if ($i != $j) { # but not the same sequence
#  				
#  				my $seq_i = $seq_names_sorted[$i];
#  				my $seq_j = $seq_names_sorted[$j];
#  				
#  				print "\ntesting seq_i = $seq_i vs. seq_j = $seq_j\n";
#  				
#  #				#----------Find index that correspond to ij pair, and assign number to ij parameters------------
#  #				$idx_ij=find_index($i,$j,\@all_indices_1,\@all_indeces_2);
#  				
#  				#------------------------------------------------------------------------------------------
#  				
#  				# Obtain the values for this pair
#  #				my $dNN_ij = $example_hh{$seq_i}->{$seq_j}->{dNN};
#  #				my $var_dNN_ij = $example_hh{$seq_i}->{$seq_j}->{var_dNN};
#  #				my $dSN_ij = $example_hh{$seq_i}->{$seq_j}->{dSN};
#  #				my $var_dSN_ij = $example_hh{$seq_i}->{$seq_j}->{var_dSN};
#  #				my $dNS_ij = $example_hh{$seq_i}->{$seq_j}->{dNS};
#  #				my $var_dNS_ij = $example_hh{$seq_i}->{$seq_j}->{var_dNS};
#  #				my $dSS_ij = $example_hh{$seq_i}->{$seq_j}->{dSS};
#  #				my $var_dSS_ij = $example_hh{$seq_i}->{$seq_j}->{var_dSS};
#  				
#  #				$d_ij=$all_d[$idx_ij];
#  #				$var_d_ij=$all_var_d[$idx_ij];
#  				
#  				# Let's do for dNN temporarily. Do we loop all four?
#  				my @seq_i_j_ordered = sort ($seq_i, $seq_j);
#  				my $min_ij = $seq_i_j_ordered[0];
#  				my $max_ij = $seq_i_j_ordered[1];
#  				
#  				my $d_ij = $d_hh{$min_ij}->{$max_ij}->{d}; # checked
#  				my $var_d_ij = $d_hh{$min_ij}->{$max_ij}->{var_d}; # checked
#  				
#  				print "\nThe sequences are $seq_i and $seq_j and we have d_ij=$d_ij\, var_d_ij=$var_d_ij\n\n";
#  				
#  				#------------------------------Calculate pi and var_pi------------------------
#  				$pi_numerator += $d_ij;
#  				$pi_counter++;
#  				
#  #				$var_pi += $var_d_ij;
#  #				$var_pi_counter++;
#  				
#  				my $covariance_sum = 0;
#  				my $covariance_counter = 0;
#  				
#  				# dNN again, or dSN? Also EVERY sequence, or start at 2? Complicated.
#  				# k and l must both be larger than i? not sure.
#  #				for (my $seq_k = 1; $seq_k <= $N; $seq_k++) {
#  #					for (my $seq_l = 1; $seq_l <= $N; $seq_l++) {
#  				for(my $k = 0; $k < $num_seq_names_Takahata; $k++) {
#  					for(my $l = $k + 1; $l < $num_seq_names_Takahata; $l++) { # k<l
#  						if ($k != $l) { # but not the same sequence for dkl
#  							
#  							my $seq_k = $seq_names_sorted[$k];
#  							my $seq_l = $seq_names_sorted[$l];
#  							
#  ###							if($seq_i eq $seq_k && $seq_j eq $seq_l) {
#  ###								# IF IT *IS* THE SAME SEQ PAIRS, then Cov(dij, dkl) = Var(dij)
#  ###								$cov_d_ij_kl = $d_ij; # VARIANCE OF THIS, WHICH IS???
#  ###								$covariance_sum += $cov_d_ij_kl;
#  ###								$covariance_counter++;
#  ###							} else {
#  							#----------Find index that correspond to kl pair, and assign number to kl parameters------------
#  	#							$idx_kl = &find_index($k, $l, \@all_indices_1, \@all_indeces_2);
#  	#							$idx_ik = &find_index($i, $k, \@all_indices_1, \@all_indeces_2);
#  	#							$idx_il = &find_index($i, $l, \@all_indices_1, \@all_indeces_2);
#  	#							$idx_jk = &find_index($j, $k, \@all_indices_1, \@all_indeces_2);
#  	#							$idx_jl = &find_index($j, $l, \@all_indices_1, \@all_indeces_2);
#  								
#  							#---------------------------------------------------------------------------------------
#  
#  								my @sorted_ij = sort($seq_i, $seq_j);
#  								my @sorted_kl = sort($seq_k, $seq_l);
#  								my @sorted_ik = sort($seq_i, $seq_k);
#  								my @sorted_jl = sort($seq_j, $seq_l);
#  								my @sorted_il = sort($seq_i, $seq_l);
#  								my @sorted_jk = sort($seq_j, $seq_k);
#  								
#  								
#  								my $min_ij = $sorted_ij[0];
#  								my $min_kl = $sorted_kl[0];
#  								my $min_ik = $sorted_ik[0];
#  								my $min_jl = $sorted_jl[0];
#  								my $min_il = $sorted_il[0];
#  								my $min_jk = $sorted_jk[0];
#  								
#  								
#  								my $max_ij = $sorted_ij[1];
#  								my $max_kl = $sorted_kl[1];
#  								my $max_ik = $sorted_ik[1];
#  								my $max_jl = $sorted_jl[1];
#  								my $max_il = $sorted_il[1];
#  								my $max_jk = $sorted_jk[1];
#  								
#  								
#  #								print "\nmin of ($seq_i\, $seq_k\) is $min_ik\n" . 
#  #									"max of ($seq_i\, $seq_k\) is $max_ik\n\n";
#  								
#  								# GET PATRISTIC DISTANCES using UPGMA $tree, compute bij values
#  								my $b_ij = &get_patristic_distance($tree, $min_ij, $max_ij);
#  								my $b_kl = &get_patristic_distance($tree, $min_kl, $max_kl);
#  								my $b_ik = &get_patristic_distance($tree, $min_ik, $max_ik);
#  								my $b_jl = &get_patristic_distance($tree, $min_jl, $max_jl);
#  								my $b_il = &get_patristic_distance($tree, $min_il, $max_il);
#  								my $b_jk = &get_patristic_distance($tree, $min_jk, $max_jk);
#  								
#  #								print "patristic distance $seq_i to $seq_j is $b_ij\n";
#  								
#  ##								$d_kl = $d_hh{$min_kl}->{$max_kl}->{d};
#  ##								$d_ik = $d_hh{$min_ik}->{$max_ik}->{d};
#  ##								$d_il = $d_hh{$min_il}->{$max_il}->{d};
#  ##								$d_jk = $d_hh{$min_jk}->{$max_jk}->{d};
#  ##								$d_jl = $d_hh{$min_jl}->{$max_jl}->{d};
#  								
#  #								print "calculating covariances for:\nseq_i=$seq_i\nseq_j=$seq_j\n" .
#  #									"seq_k=$seq_k\nseq_l=$seq_l\n";
#  								
#  								# MUST PASS b's not d's
#  								#my ($b_ij, $b_kl, $b_ik, $b_jl, $b_il, $b_jk, $m) = @_; 
#  ##								$cov_d_ij_kl = &cov_calculation_NeiJin($d_ij, $d_kl, $d_ik, $d_jl, $d_il, $d_jk, $m);
#  								
#  								$cov_d_ij_kl = &cov_calculation_NeiJin($b_ij, $b_kl, $b_ik, $b_jl, $b_il, $b_jk, $m);
#  								
#  								#print "\ncovariance of $seq_i / $seq_j and $seq_k / $seq_l is $cov_d_ij_kl\n";
#  #								print "b_jk was $b_jk\n";
#  #								if($seq_i eq '1' && $seq_j eq '2' && $seq_k eq '3' && $seq_l eq '4') {
#  #									print "b_ij=$b_ij, b_kl=$b_kl, b_ik=$b_ik, b_jl=$b_jl, b_il=$b_il, b_jk=$b_jk\n";
#  #								}
#  								
#  								if($seq_i eq $seq_k && $seq_j eq $seq_l) {
#  									print "\nseq_i=$seq_i\, seq_j=$seq_j\, seq_k=$seq_k\, seq_l=$seq_l\n";
#  									print "covariance=$cov_d_ij_kl\n";
#  								}
#  								
#  								my $ij_name = "d$seq_i$seq_j";
#  								my $kl_name = "d$seq_k$seq_l";
#  								
#  								$store_covariances{$ij_name}->{$kl_name} = $cov_d_ij_kl;
#  								
#  								$covariance_sum += $cov_d_ij_kl;
#  								$covariance_counter++;
#  ###							} # not same seq comparisons
#  						} # k != l
#  					} # finished all l seqs
#  				} # finished all k seqs
#  				
#  				print "\nI determined $covariance_counter covariances\n";
#  #				print "\ncovariance of $seq_i and $seq_j is $cov_d_ij_kl\n";
#  				
#  				$var_pi_numerator += ($var_d_ij + $covariance_sum);
#  #				$var_pi_numerator += $covariance_sum;
#  				$var_pi_counter++;
#  			} # i != j
#  		} # finished all j seqs
#  	} # finished all i seqs
# 	
# 	print "\nI determined $var_pi_counter variances\n";
# 	
# 	print "\nI determined these covariances: \n"; # . %store_covariances . "\n";#$store_covariances{$ij_name}->{$kl_name} = $cov_d_ij_kl;;
# 	
# 	foreach my $ij_name (sort keys %store_covariances) {
# 		foreach my $kl_name (sort keys %{$store_covariances{$ij_name}}) {
# 			print "covariance $ij_name and $kl_name = $store_covariances{$ij_name}->{$kl_name}\n";
# 		}
# 	}
# 	
# #	print "42 vs 23 is: " . $store_covariances{d42}->{d23} * 10 . "\n";
# 	
# 	# CALCULATE VARIANCE HERE
# 	# $pi_counter == ($N**2 - $N) / 2 ?
# #	my $pi = $pi_numerator / (($N**2 - $N));
# 	my $pi = 2 * $pi_numerator / ($N * ($N - 1));
# #	my $pi = $pi_numerator / $pi_counter;
# #	my $var_pi = $var_pi_numerator / (($N * ($N - 1)) ** 2);
# 	my $var_pi = 4 / (($N * ($N - 1))**2) * $var_pi_numerator;
# 	
# #	print "\ndenominator is: " . (($N * ($N - 1)) ** 2) . "\n";
# 	
# 	my @result = ($pi, $var_pi);
# 	return @result;
# }

#sub find_index { # Not finished
#	my ($idx_i, $idx_j, $indices1, $indices2)=@_;
#	my $ind_ij
#	#iterate over @$indices1 instead of @indices1
#	return $ind_ij
#}


#########################################################################################
sub cov_calculation_NeiJin {
	#Eqs1. b_ij = d_ij, for all i and j
	
	#Eqs2: b1 = b_ij + b_kl; 
	#	b2 = b_ik + b_jl; 
	#	b3 = b_il + b_jk
	
	#Eq3: b = (b1 - min[b1,b2,b3]) / 2
	
	#Eq4: p = 3/4 * (1 - exp(-4/3b))
	
	#Eq5: Cov(d_ij, d_kl) = 9p(1-p) / ((3-4p)^2m), 
	#	m is the number of nucleotides examined per sequence.
	
	#Eq6: V(pi) = (sum(V(d_ij)) + sum(cov(d_ij,d_kl))) / (n(n-1))^2
	#	all i,j pair i not equal to j
	#	all ij,kl pair, i not equal to j and k not equal to l
	
	my ($b_ij, $b_kl, $b_ik, $b_jl, $b_il, $b_jk, $m) = @_; 
	# $m is the number of nucleotide sites examined
	
###	my ($d_ij, $d_kl, $d_ik, $d_jl, $d_il, $d_jk, $m, $tree) = @_; 
	# if i==k or j==l, some of these will be undefined (0)
	
#	print "\nHere's my array: @_\n\n";
	
#	my $b1 = $d_ij + $d_kl; # FALSE: should be b_ij, b_kl
#	my $b2 = $d_ik + $d_jl; # FALSE: should be b_ik, b_jl
#	my $b3 = $d_il + $d_jk; # FALSE: should be b_il, b_jk
	
	my $b1 = $b_ij + $b_kl;
	my $b2 = $b_ik + $b_jl;
	my $b3 = $b_il + $b_jk;
	
	my $bm = min($b1, $b2, $b3);
	my $b = ($b1 - $bm) / 2;
	my $p = 3/4 * (1 - exp((-4/3) * $b));
#	my $cov_ij_kl = 9 * $p * (1 - $p) / (2 * $m * (3 - 4 * $p));
	my $cov_ij_kl = 9 * $p * (1 - $p) / (((3 - 4 * $p) ** 2) * $m); # SAME as VARIANCE OF SHARED BRANCH
	
#	print "\n3 squared is " . 3 ** 2 . "\n";
	
#	print "b1=$b1 b2=$b2 b3=$b3 bm=$bm b=$b p=$b m=$m cov=$cov_ij_kl\n";
	
	if($print_flag == 1) {
		print "bij = b54 = $b_ij\n";
		print "bkl = b14 = $b_kl\n";
		print "bik = b51 = $b_ik\n";
		print "bjl = b44 = $b_jl\n";
		print "bil = b54 = $b_il\n";
		print "bjk = b41 = $b_jk\n";
		print "b1 = $b1\n";
		print "b2 = $b2\n";
		print "b3 = $b3\n";
		print "bm = $bm\n";
		print "b = $b\n";
		print "p = $p\n";
		print "cov_ij_kl = $cov_ij_kl\n";
	}
	
	return $cov_ij_kl;
}

#########################################################################################
sub UPGMA_from_distance_hash {
	# The provided distance hash must have distance values saved in sorted order, such that
	# the first key is always less than the second key. We shall sort alphabetically 
	# (only using sort), not numerically, since the ID's will likely be FASTA headers.
	
	# EXAMPLE: $common_chimp{1}->{2} = .0168
	
	my ($d_hh_ref) = @_; 
	my %d_hh = %{$d_hh_ref};
	
	# Record all seq IDs by going through hash and storing unique names in LOOP
	my %all_seq_IDs_hash;
	foreach my $seq_ID_1 (sort keys %d_hh) {
		$all_seq_IDs_hash{$seq_ID_1} = 1;
		
		foreach my $seq_ID_2 (sort keys %{$d_hh{$seq_ID_1}}) {
			$all_seq_IDs_hash{$seq_ID_2} = 1;
		}
	}
	
	my @all_seq_IDs = sort keys %all_seq_IDs_hash;
	
	
	# RE-NAME WITH SIMPLE (LETTERS ONLY) NAMES
	my %d_hh_simple_to_original_name;
	my %d_hh_original_to_simple_name;
	
	my $simple_name = 'A';
	my @all_simple_names;
	
	foreach(@all_seq_IDs) {
		my $original_name = $_;
		
		$d_hh_simple_to_original_name{$simple_name} = $original_name;
		$d_hh_original_to_simple_name{$original_name} = $simple_name;
		
		push(@all_simple_names, $simple_name);
		
		$simple_name++;
	}
	
	# RE-NAME VALUES WITHIN DISTANCE MATRIX
	my %d_hh_renamed;
	foreach my $seq_ID_1 (sort keys %d_hh) {
		
		foreach my $seq_ID_2 (sort keys %{$d_hh{$seq_ID_1}}) {
			my $simple_ID_1 = $d_hh_original_to_simple_name{$seq_ID_1};
			my $simple_ID_2 = $d_hh_original_to_simple_name{$seq_ID_2};
			
			$d_hh_renamed{$simple_ID_1}->{$simple_ID_2} = $d_hh{$seq_ID_1}->{$seq_ID_2};
			
			delete $d_hh{$seq_ID_1}->{$seq_ID_2};
			
		}
		
		delete $d_hh{$seq_ID_1};
		
	}
	
	# Build the UPGMA tree

#	my @clusters = @all_seq_IDs;

	my @clusters = @all_simple_names;
	
	# Find first sequence pair
	while(scalar(@clusters) > 1) { # while there are sequences we haven't composited
		my $curr_min_dist = 10000000; # something crazy big
		my @curr_min_cluster_indices;
		
		my @cluster_1_members; # 1 or more
		my @cluster_2_members; # 1 or more
		
#		print "\n";
		
		################
		# FIND CLUSTER 1
		for(my $i = 0; $i < scalar(@clusters); $i++) {
			my $cluster_1 = $clusters[$i];
			
			# First get rid of all ()
			$cluster_1 =~ s/\(//g;
			$cluster_1 =~ s/\)//g;
			
			#print "\ncluster 1 before: $cluster_1\n";
			
			# Get rid of branch lengths
			$cluster_1 =~ s/\:[\.\d]+//g;
				
			#print "\ncluster 1 after: $cluster_1\n";
				
			# Taxa are now separated by commas (,)
			@cluster_1_members = split(",", $cluster_1);
			
#			print "cluster_1_members: @cluster_1_members\n";
			
			################
			# FIND CLUSTER 2
			for(my $j = $i + 1; $j < scalar(@clusters); $j++) {
				my $cluster_2 = $clusters[$j];
				
				# First get rid of all ()
				$cluster_2 =~ s/\(//g;
				$cluster_2 =~ s/\)//g;
				
				#print "\ncluster 2 before: $cluster_2\n";
				# Get rid of branch lengths
				$cluster_2 =~ s/\:[\.\d]+//g;
				#print "\ncluster 2 after: $cluster_2\n";
				
				# Taxa are now separated by commas (,)
				@cluster_2_members = split(",", $cluster_2);
				
#				print "cluster_2_members: @cluster_2_members\n";
				
				####################
				# FOUND TWO CLUSTERS
				# Check to see if it's the cluster pair with the smallest distance thus far
				my $this_num_comparisons;
				my $this_sum_pw_distances;
				
				foreach my $cluster_1_member (@cluster_1_members) {
					foreach my $cluster_2_member (@cluster_2_members) {
						my @member_pair = sort ($cluster_1_member, $cluster_2_member);
						
						my $sorted_cluster_1_member = $member_pair[0];
						my $sorted_cluster_2_member = $member_pair[1];
						
#						print "distance $cluster_1_member to $cluster_2_member is $d_hh{$cluster_1_member}->{$cluster_2_member}\n";
						
#						$this_sum_pw_distances += $d_hh{$sorted_cluster_1_member}->{$sorted_cluster_2_member};
						$this_sum_pw_distances += $d_hh_renamed{$sorted_cluster_1_member}->{$sorted_cluster_2_member};
						$this_num_comparisons++;
					}
				}
				
				# Calculate mean distance between these two clusters
				my $this_mean_dist = $this_sum_pw_distances / $this_num_comparisons;
				
#				print "this_mean_dist is $this_mean_dist\n";
				
				# If it's the smallest seen so far, save it
				if($this_mean_dist < $curr_min_dist) {
					$curr_min_dist = $this_mean_dist;
					@curr_min_cluster_indices = ($i, $j);
				}
			}
		}
		
		# FOUND CLUSTER PAIR WITH MIN DISTANCE, so join them
		my @min_cluster_pair = sort($clusters[$curr_min_cluster_indices[0]], $clusters[$curr_min_cluster_indices[1]]);
		
#		print "cluster pair: @min_cluster_pair\n";
		
		my $new_cluster_1;
		my $new_cluster_2;
		
		# Here, right here, add the branch lengths: we have $curr_min_dist
		# The new branch length will be $curr_min_dist/2, SUBTRACTED FROM half the mean
		# distance between all sequences in the composite OTU being joined. 
		my $half_curr_min_dist = $curr_min_dist / 2;
		
		# Join clusters and add branch lengths
		
		# SOLUTION: ((((B:0.5, F:0.5):5.75, G:6.25):2.0, (A:4.0, D:4.0):4.25):6.25, C:14.5):2.5, E:17.0
		# B:0.5,F:0.5
		# B:0.5,F:0.5	A:4.0,D:4.0
		# (B:0.5,F:0.5):5.75,G:6.25		A:4.0,D:4.0
		# (B:0.5,F:0.5):5.75,G:6.25		A:4.0,D:4.0
		# ((B:0.5,F:0.5):5.75,G:6.25):2.0,(A:4.0,D:4.0):4.25
		
		#((((B:0.5,F:0.5):5.75,G:6.25):2,(A:4,D:4):4.25):6.25,C:14.5):2.5,E:17
		#((((B:0.5,F:0.5):5.75,G:6.25):2,(A:4,D:4):4.25):6.25,C:14.5):2.5,E:17
		
		if($min_cluster_pair[0] =~ /\,/) {
			###################################################
			### cluster 1 is COMPOSITE / cluster 2 is COMPOSITE
			if($min_cluster_pair[1] =~ /\,/) {
				
				my $curr_cluster_1 = $min_cluster_pair[0];
				my $curr_cluster_2 = $min_cluster_pair[1];
				
				# CLUSTER 1
				my $cluster_1_length;
				
				if($curr_cluster_1 =~ /\(/) { # it's got parentheses already
					
					# Reduce all parentheses
					while($curr_cluster_1 =~ /\(([\w\:\.\,\d]+,[\w\:\.\,\d]+)\)/) { # while we've still got branch lengths inside parentheses
						
						my $inside_parentheses = $1;
						
						# Find sum of numbers inside parentheses
						my @numbers_inside;
						while($inside_parentheses =~ /([\.\d]+)/g) {
							push(@numbers_inside, $1);
						}
						
						my $inside_sum = 0;
						foreach(@numbers_inside) {
							$inside_sum += $_;
						}
						
						my $half_inside_sum = $inside_sum / 2;
						
						$curr_cluster_1 =~ s/\($inside_parentheses\)/$half_inside_sum/;
					}
					
					my $remaining_numbers_sum = 0;
					
					# Now deal with all remaining numbers
					while($curr_cluster_1 =~ /[\d\.]+/g) { # iterates over all matches
						$remaining_numbers_sum += $&;
					}
					
					my $cluster_1_to_tip = $remaining_numbers_sum / 2;
					
					# Subtract this from half the mean distance calculated earlier
					$cluster_1_length = $half_curr_min_dist - $cluster_1_to_tip;
					
				} else { # composite of only two simple OTUs, e.g., B:0.5,F:0.5, or unresolved
					
					my $curr_length_sum;
					my $curr_length_count;
					
					while($curr_cluster_1 =~ /\:([\d\.]+)/g) { # iterates over all matches
						$curr_length_sum += $1;
						$curr_length_count++;
					}
					
					my $cluster_1_to_tip = $curr_length_sum / $curr_length_count;
					$cluster_1_length = $half_curr_min_dist - $cluster_1_to_tip;
				}
				
				# CLUSTER 2
				my $cluster_2_length;
				
				if($curr_cluster_2 =~ /\(/) { # it's got parentheses already
					
					# Reduce all parentheses
					while($curr_cluster_2 =~ /\(([\w\:\.\,\d]+,[\w\:\.\,\d]+)\)/) { # while we've still got branch lengths inside parentheses
						
						my $inside_parentheses = $1;
						
						# Find sum of numbers inside parentheses
						my @numbers_inside;
						while($inside_parentheses =~ /([\.\d]+)/g) {
							push(@numbers_inside, $1);
						}
						
						my $inside_sum = 0;
						foreach(@numbers_inside) {
							$inside_sum += $_;
						}
						
						my $half_inside_sum = $inside_sum / 2;
						
						$curr_cluster_2 =~ s/\($inside_parentheses\)/$half_inside_sum/;
					}
					
					my $remaining_numbers_sum = 0;
					# Now deal with all remaining numbers
					while($curr_cluster_2 =~ /[\d\.]+/g) { # iterates over all matches
						$remaining_numbers_sum += $&;
					}
					
					my $cluster_2_to_tip = $remaining_numbers_sum / 2;
					
					# Subtract this from half the mean distance calculated earlier
					$cluster_2_length = $half_curr_min_dist - $cluster_2_to_tip;
					
				} else { # composite of only two simple OTUs, e.g., B:0.5,F:0.5, or unresolved
					
					my $curr_length_sum;
					my $curr_length_count;
					
					while($curr_cluster_2 =~ /\:([\d\.]+)/g) { # iterates over all matches
						$curr_length_sum += $1;
						$curr_length_count++;
					}
					
					my $cluster_2_to_tip = $curr_length_sum / $curr_length_count;
					$cluster_2_length = $half_curr_min_dist - $cluster_2_to_tip;
				}
				
				$new_cluster_1 = "\($min_cluster_pair[0]\):$cluster_1_length";
				$new_cluster_2 = "\($min_cluster_pair[1]\):$cluster_2_length";
		
		
			################################################
			### cluster 1 is COMPOSITE / cluster 2 is SIMPLE
			} else { # surround first in parentheses
				
				# Find length of branch basal to composite cluster 1
				my $curr_cluster_1 = $min_cluster_pair[0];
				
				# Start innermost and work out
				# EXAMPLE:
				# ((B:0.5,F:0.5):5.75,G:6.25):2.0,(A:4.0,D:4.0):4.25 ==> find (B:0.5,F:0.5), replace with 0.5
				# (0.5:5.75,G:6.25):2.0,(A:4.0,D:4.0):4.25 ==> find (0.5:5.75,G:6.25), replace 6.25
				# 6.25:2.0,(A:4.0,D:4.0):4.25 ==> find (A:4.0,D:4.0), replace 4.0
				# 6.25:2.0,4.0:4.25 # THERE ARE NO : inside (), so ADD ALL NUMBERS and divide by 2 ==>
				# 6.25	2.0	4.0	4.25 ==> extract numbers from each element, sum, divide by num elements ==>
				# (6.25 + 2.0 + 4.0 + 4.25) / 2 = 8.25
				# so, for example, 14.5 - 8.25 = 6.25, which is correct for this example
				
				my $cluster_1_length;
				
				# /\(([\w\d\.\-]+\:[\.\d]+)\)/
				if($curr_cluster_1 =~ /\(/) { # it's got parentheses already
					
					# Reduce all parentheses
					while($curr_cluster_1 =~ /\(([\w\:\.\,\d]+,[\w\:\.\,\d]+)\)/) { # while we've still got branch lengths inside parentheses
						
						my $inside_parentheses = $1;
						
						# Find sum of numbers inside parentheses
						my @numbers_inside;
						while($inside_parentheses =~ /([\.\d]+)/g) {
							push(@numbers_inside, $1);
						}
						
						my $inside_sum = 0;
						foreach(@numbers_inside) {
							$inside_sum += $_;
						}
						
						my $half_inside_sum = $inside_sum / 2;
						
						$curr_cluster_1 =~ s/\($inside_parentheses\)/$half_inside_sum/;
					}
					
					my $remaining_numbers_sum = 0;
					
					# Now deal with all remaining numbers
					while($curr_cluster_1 =~ /[\d\.]+/g) { # iterates over all matches
						$remaining_numbers_sum += $&;
					}
					
					my $cluster_1_to_tip = $remaining_numbers_sum / 2;
					
					# Subtract this from half the mean distance calculated earlier
					$cluster_1_length = $half_curr_min_dist - $cluster_1_to_tip;
					
				} else { # composite of only two simple OTUs, e.g., B:0.5,F:0.5, or unresolved
					
					my $curr_length_sum;
					my $curr_length_count;
					
					while($curr_cluster_1 =~ /\:([\d\.]+)/g) { # iterates over all matches
						$curr_length_sum += $1;
						$curr_length_count++;
					}
					
					my $cluster_1_to_tip = $curr_length_sum / $curr_length_count;
					$cluster_1_length = $half_curr_min_dist - $cluster_1_to_tip;
				}
				
				$new_cluster_1 = "\($min_cluster_pair[0]\):$cluster_1_length";
				
				$new_cluster_2 = "$min_cluster_pair[1]\:$half_curr_min_dist";
				
			} # cluster 2 was simple
			
		} else { # no commas inside cluster 1, so it's a simple OTU
			
			###########################
			### cluster 1 is SIMPLE / cluster 2 is COMPOSITE
			if($min_cluster_pair[1] =~ /\,/) {
			
				my $curr_cluster_2 = $min_cluster_pair[1];
				
				# Cluster 2 is difficult
				my $cluster_2_length;
				
				if($curr_cluster_2 =~ /\(/) { # it's got parentheses already
					
					# Reduce all parentheses
					while($curr_cluster_2 =~ /\(([\w\:\.\,\d]+,[\w\:\.\,\d]+)\)/) { # while we've still got branch lengths inside parentheses
						
						my $inside_parentheses = $1;
						
						# Find sum of numbers inside parentheses
						my @numbers_inside;
						while($inside_parentheses =~ /([\.\d]+)/g) {
							push(@numbers_inside, $1);
						}
						
						my $inside_sum = 0;
						foreach(@numbers_inside) {
							$inside_sum += $_;
						}
						
						my $half_inside_sum = $inside_sum / 2;
						
						$curr_cluster_2 =~ s/\($inside_parentheses\)/$half_inside_sum/;
						
					}
					
					my $remaining_numbers_sum = 0;
					
					# Now deal with all remaining numbers
					while($curr_cluster_2 =~ /[\d\.]+/g) { # iterates over all matches
						$remaining_numbers_sum += $&;
					}
					
					my $cluster_2_to_tip = $remaining_numbers_sum / 2;
					
					# Subtract this from half the mean distance calculated earlier
					$cluster_2_length = $half_curr_min_dist - $cluster_2_to_tip;
					
				} else { # composite of only two simple OTUs, e.g., B:0.5,F:0.5, or unresolved
					
					my $curr_length_sum;
					my $curr_length_count;
					
					while($curr_cluster_2 =~ /\:([\d\.]+)/g) { # iterates over all matches
						$curr_length_sum += $1;
						$curr_length_count++;
					}
					
					my $cluster_2_to_tip = $curr_length_sum / $curr_length_count;
					$cluster_2_length = $half_curr_min_dist - $cluster_2_to_tip;
				}
				
				$new_cluster_1 = "$min_cluster_pair[0]\:$half_curr_min_dist";
				$new_cluster_2 = "\($min_cluster_pair[1]\):$cluster_2_length";
				
				
			########################
			### cluster 1 is SIMPLE / cluster 2 is SIMPLE
			} else { # don't enclose: want "B:0.5,F:0.5"
				$new_cluster_1 = "$min_cluster_pair[0]\:$half_curr_min_dist";
				$new_cluster_2 = "$min_cluster_pair[1]\:$half_curr_min_dist";
			}
		}
		
		my $new_cluster = "$new_cluster_1\,$new_cluster_2";
		
		# Add to end of clusters
		push(@clusters, $new_cluster);
		
		# Delete the two clusters we've subsumed, in reverse order based on index
		# second one is always larger in @curr_min_cluster_indices
		splice(@clusters, $curr_min_cluster_indices[1], 1);
		splice(@clusters, $curr_min_cluster_indices[0], 1);
		
	} # while there is more than 1 cluster
	
	# Now replace all names with the originals, in reverse order (e.g., simple name AA before A)
	my $tree = "\($clusters[0]\)\;";
	
	foreach my $simple_name (reverse (@all_simple_names)) {
		my $this_original_name = $d_hh_simple_to_original_name{$simple_name};
		$tree =~ s/$simple_name/$this_original_name/; # should only be one of each, so global not needed
	}
	
#	print "\n\nTREE FOUND IS: $tree\n\n";
	
	return $tree;
	
}


#########################################################################################
sub get_patristic_distance {
	my ($tree, $taxon1, $taxon2) = @_;
	
	if($taxon1 eq $taxon2) {
		return 0;
	}
	
#	print "tree = $tree\n";
	
	# First, record all taxon names to replace with letter codes for ease and consistency
	my %all_taxa_names;
	while($tree =~ /([\w\d\.\-]+)\:/g) {
		$all_taxa_names{$1} = 1;
	}
	
	my @all_taxa_names = sort keys %all_taxa_names;
	
	#print "all the taxa: @all_taxa_names\n";
	
	# RE-NAME WITH SIMPLE (LETTERS ONLY) NAMES
	my %simple_to_original_name;
	my %original_to_simple_name;
	
	my $simple_name = 'A';
	my @all_simple_names;
	
	foreach(@all_taxa_names) {
		my $original_name = $_;
		
		$simple_to_original_name{$simple_name} = $original_name;
		$original_to_simple_name{$original_name} = $simple_name;
		
		push(@all_simple_names, $simple_name);
		
		$simple_name++;
	}
	
	
	# REPLACE TREE WITH SIMPLE NAMES
	
	# First, must be replace in order of reverse length, or we risk replacing substrings
#	print "\nall_taxa_names = @all_taxa_names\n";
	my $original_names_by_length_ref = &sort_array_by_elt_length(\@all_taxa_names);
	my @original_names_by_length = @{$original_names_by_length_ref};
	
#	print "\nall_taxa_names by length = @original_names_by_length\n";
	
	my $tree_simple = $tree;
	
	foreach my $original_name (reverse @original_names_by_length) {
		my $this_simple_name = $original_to_simple_name{$original_name};
		$tree_simple =~ s/$original_name\:/$this_simple_name\:/;
	}
	
#	print "the simple tree version is: $tree_simple\n";
	
	###################
	# GET DISTANCES
	# Let's take the approach of stripping the tree from the outside until we get the 
	# minimum subtree still containing both the taxa of interest. What then?
	
	# EXAMPLE: PRUNE before trying to find A-to-G distance
	# (((((B:0.5,F:0.5):5.75,G:6.25):2,(A:4,D:4):4.25):6.25,C:14.5):2.5,E:17); ==> remove outside parentheses, anything directly inside
	# ((((B:0.5,F:0.5):5.75,G:6.25):2,(A:4,D:4):4.25):6.25,C:14.5) ==> same again
	# (((B:0.5,F:0.5):5.75,G:6.25):2,(A:4,D:4):4.25) ==> this one's rough to regex
	# ((B:0.5,F:0.5):5.75,G:6.25) ==> don't have both A and G anymore, so use the previous tree
	
	# EXAMPLE
	# (((B:0.5,F:0.5):5.75,G:6.25):2,(A:4,D:4):4.25) ==> find A-to-G
	# same elimination and replacement approach looking for clusters of our taxa or x inside ()
	# (((B:0.5,F:0.5):5.75,G:6.25):2,(A:4,D:4):4.25) ==> first cluster doesn't contain; replace with o
	# ((o:5.75,G:6.25):2,(A:4,D:4):4.25) ==> first cluster does contain; add the length of the containing and replace with x (curr sum is 6.25)
	# (x:2,(A:4,D:4):4.25) ==> first cluster does contain; add length of containing and replace with x (curr sum is 10.25)
	# (x:2,x:4.25) ==> first cluster does contain; add length of containing to reach 10.25+2+4.25=16.5, the correct answer
	# Q.E.D.
	
	# HARDER: find B-to-E
	# pruning removes nothing:
	# (((((B:0.5,F:0.5):5.75,G:6.25):2,(A:4,D:4):4.25):6.25,C:14.5):2.5,E:17) ==> first contains B; sum = 0.5
	# ((((x:5.75,G:6.25):2,(A:4,D:4):4.25):6.25,C:14.5):2.5,E:17) ==> first contains x; sum = 0.5+5.75
	# (((x:2,(A:4,D:4):4.25):6.25,C:14.5):2.5,E:17) ==> first contains NOTHING; sum stays 0.5+5.75
	# (((x:2,o:4.25):6.25,C:14.5):2.5,E:17) ==> first contains x; sum = 0.5+5.75+2
	# ((x:6.25,C:14.5):2.5,E:17) ==> first contains x; sum = 0.5+5.75+2+6.25
	# (x:2.5,E:17) ==> first contains x and E; sum = 0.5+5.75+2+6.25+2.5+17 = 34
	# CORRECT ANSWER.
	
	# So above algorithm should work.
	
	my $search1 = $original_to_simple_name{$taxon1};
	my $search2 = $original_to_simple_name{$taxon2};
	
#	print "searching for: $search1 \($taxon1\) and $search2 \($taxon2\)\n\n";
	
#	print "simple tree = $tree_simple\n";
	
#		#################################################
#		# STEP 1: PRUNE TO MINIMUM SUBTREE FROM OUTSIDE
#		my $curr_tree = $tree_simple;
#		$curr_tree =~ s/;//; # get rid of closing semicolon
#		my $previous_subtree = $curr_tree;
#		
#		while($previous_subtree =~ /$search1\:/ && $previous_subtree =~ /$search2\:/) {
#	#		while($previous_subtree =~ /$search2\:/) {
#				# Previous tree is a possibility; save
#				$curr_tree = $previous_subtree;
#				print "curr_tree is $curr_tree\n";
#				
#				### START OF TREE TRIMMING
#				# Trim it and see if it remains viable
#				# Nix first parentheses and anything directly inside
#				$previous_subtree =~ s/^\([\d\:\.\,\w\-]*//;
#				# If there's a lonely leader, end it
#				while($previous_subtree =~ s/^\([\d\:\.\,\w\-]*\)[\d\:\.\,\w\-]+\(//) { # will always be succeeded by a comma if it's not the whole tree
#					# Couldn't get this to work another way at 3AM.
#					print "previous subtree after leader trim = $previous_subtree\n";
#					#exit;
#				}
#				
#				### END OF TREE TRIMMING
#				# Nix last parentheses and anything directly inside
#				$previous_subtree =~ s/[\d\:\.\,\w\-]*\)$//;
#				print "previous subtree = $previous_subtree\n";
#				# If there's a straggling cluster, end it
#				while($previous_subtree =~ s/[\d\:\.\,\w\-]+\([\d\:\.\,\w\-]*\)$//) { # will always be preceded by a comma if it's not the whole tree
#	#				my $match = $&;
#	#				print "match = $match\n";
#	#				my $match = $&;
#	#				$previous_subtree =~ s/$match$//;
#					print "previous subtree after end trim = $previous_subtree\n";
#					#exit;
#				}
#				
#				#print "pruning outermost tree layer resulted in: $previous_subtree\n\n";
#	#		}
#		}
	
	#################################################
	# COLLAPSE CLUSTERS FROM THE INSIDE-OUT
	
	# DIFFICULT EXAMPLE: find B-to-E distance
	# pruning removes nothing:
	# (((((B:0.5,F:0.5):5.75,G:6.25):2,(A:4,D:4):4.25):6.25,C:14.5):2.5,E:17) ==> first contains B; sum = 0.5
	# ((((x:5.75,G:6.25):2,(A:4,D:4):4.25):6.25,C:14.5):2.5,E:17) ==> first contains x; sum = 0.5+5.75
	# (((x:2,(A:4,D:4):4.25):6.25,C:14.5):2.5,E:17) ==> first contains NOTHING; sum stays 0.5+5.75
	# (((x:2,o:4.25):6.25,C:14.5):2.5,E:17) ==> first contains x; sum = 0.5+5.75+2
	# ((x:6.25,C:14.5):2.5,E:17) ==> first contains x; sum = 0.5+5.75+2+6.25
	# (x:2.5,E:17) ==> first contains x and E; sum = 0.5+5.75+2+6.25+2.5+17 = 34
	# CORRECT ANSWER.
	
	# EXAMPLE: A-to-D
	# (((B:0.5,F:0.5):5.75,G:6.25):2,(A:4,D:4):4.25) ==> first cluster doesn't contain; replace with o
	# ((o:5.75,G:6.25):2,(A:4,D:4):4.25) ==> first cluster does't contain; replace with o
	# (o:2,(A:4,D:4):4.25) ==> first cluster does contain BOTH; add length of containing and replace with x (curr sum is 4)
	# TWO OPTIONS: STOP when we've seen both; or else don't add remaining x values if there are no actual values.
	# I think these two are equivalent. Let's try the first.
	# Q.E.D.
	
	my $patristic_distance_sum = 0;
	my $seen_search1 = 0;
	my $seen_search2 = 0;
	### Tree still has our taxa or an x
#	while($tree_simple =~ /x:/ || $tree_simple =~ /$search1\:/ || $tree_simple =~ /$search2\:/) {
	while($tree_simple =~ /$search1\:/ || $tree_simple =~ /$search2\:/) {
		#print "current remaining tree is $tree_simple\n";
		
		if($tree_simple =~ /\(([\w\.\d]+)\:([\w\.\d]+),([\w\.\d]+)\:([\w\.\d]+)\)/) { # while we've still got clusters inside parentheses
			my $first_taxon = $1;
			my $first_length = $2;
			my $second_taxon = $3;
			my $second_length = $4;
			my $entire_match = "$1\:$2,$3\:$4";
			
#			if($first_taxon eq $search1) {
#				$seen_search1 = 1;
#			} elsif($first_taxon eq $search2) {
#				$seen_search2 = 1;
#			}
#			
#			if($second_taxon eq $search1) {
#				$seen_search1 = 1;
#			} elsif($second_taxon eq $search2) {
#				$seen_search2 = 1;
#			}
			
#			print "we matched $entire_match\n";
#			print "we found a cluster with $first_taxon and $second_taxon\n";
			
			my $contained_relevant_taxa = 0;
			
			if($first_taxon eq 'x' || $first_taxon eq $search1 || $first_taxon eq $search2) {
				$patristic_distance_sum += $first_length;
				$contained_relevant_taxa = 1;
			}
			
			if($second_taxon eq 'x' || $second_taxon eq $search1 || $second_taxon eq $search2) {
				$patristic_distance_sum += $second_length;
				$contained_relevant_taxa = 1;
			}
			
			my $replacement_holder;
			
			if($contained_relevant_taxa == 1) {
				$replacement_holder = 'x';
			} else {
				$replacement_holder = 'o';
			}
			
			$tree_simple =~ s/\($entire_match\)/$replacement_holder/;
			
		} else { # we've finished reducing, and must simply sum what remains
			while($tree_simple =~ /([\w\.\d]+)\:([\w\.\d]+)/g) {
				my $this_taxon = $1;
				my $this_length = $2;
				
				print "we found a cluster with $this_taxon\n";
				
				if($this_taxon eq 'x' || $this_taxon eq $search1 || $this_taxon eq $search2) {
					$patristic_distance_sum += $this_length;
				}
				
			}
			last;
		}
		
#		print "new tree is $tree_simple\n";
		
	}
	
	# Check to see if we have two x's in a single cluster
	while($tree_simple =~ /\(x\:([\w\.\d]+),x\:([\w\.\d]+)\)/g) {
	#print "current remaining tree is $tree_simple\n";
	
		my $first_length = $1;
		my $second_length = $2;
		
#		print "we found two clustered x's\n";
		
		$patristic_distance_sum += ($first_length + $second_length);
	}
	
#	print "new tree is $tree_simple\n";
		
		
		
##		if($curr_cluster_2 =~ /\(/) { # it's got parentheses already
##			
##			# Reduce all parentheses
##			while($curr_cluster_2 =~ /\(([\w\:\.\,\d]+,[\w\:\.\,\d]+)\)/) { # while we've still got branch lengths inside parentheses
##				
##				my $inside_parentheses = $1;
##				
##				# Find sum of numbers inside parentheses
##				my @numbers_inside;
##				while($inside_parentheses =~ /([\.\d]+)/g) {
##					push(@numbers_inside, $1);
##				}
##				
##				my $inside_sum = 0;
##				foreach(@numbers_inside) {
##					$inside_sum += $_;
##				}
##				
##				my $half_inside_sum = $inside_sum / 2;
##				
##				$curr_cluster_2 =~ s/\($inside_parentheses\)/$half_inside_sum/;
##				
##			}
##			
##			my $remaining_numbers_sum = 0;
##			
##			# Now deal with all remaining numbers
##			while($curr_cluster_2 =~ /[\d\.]+/g) { # iterates over all matches
##				$remaining_numbers_sum += $&;
##			}
##			
##			my $cluster_2_to_tip = $remaining_numbers_sum / 2;
##			
##			# Subtract this from half the mean distance calculated earlier
##			$cluster_2_length = $half_curr_min_dist - $cluster_2_to_tip;
##			
##		} else { # composite of only two simple OTUs, e.g., B:0.5,F:0.5, or unresolved
##			
##			my $curr_length_sum;
##			my $curr_length_count;
##			
##			while($curr_cluster_2 =~ /\:([\d\.]+)/g) { # iterates over all matches
##				$curr_length_sum += $1;
##				$curr_length_count++;
##			}
##			
##			my $cluster_2_to_tip = $curr_length_sum / $curr_length_count;
##			$cluster_2_length = $half_curr_min_dist - $cluster_2_to_tip;
##		}
##		
##		$new_cluster_1 = "$min_cluster_pair[0]\:$half_curr_min_dist";
##		$new_cluster_2 = "\($min_cluster_pair[1]\):$cluster_2_length";
##		
##		
##	########################
##	### cluster 1 is SIMPLE / cluster 2 is SIMPLE
##	}

#	#################################################
#	# RE-INSERT THE ORIGINAL TAXA NAMES
#	foreach my $simple_name (reverse (@all_simple_names)) {
#		my $this_original_name = $simple_to_original_name{$simple_name};
#		$curr_tree =~ s/$simple_name/$this_original_name/; # should only be one of each, so global not needed
#	}
#	
#	print "comparing $taxon1 to $taxon2 resulted in subtree: $curr_tree\n\n";
	
	return $patristic_distance_sum;
}


#########################################################################################
sub convert_dhh_to_distmatrix {
	my ($d_hh_ref) = @_; 
	my %d_hh = %{$d_hh_ref};
	
	my %distance_only_matrix;
	
	foreach my $key1 (sort keys %d_hh) {
		foreach my $key2 (sort keys %{$d_hh{$key1}}) {
			$distance_only_matrix{$key1}->{$key2} = $d_hh{$key1}->{$key2}->{d};
		}
	}
	
	return \%distance_only_matrix;
}


#########################################################################################
sub sort_array_by_elt_length {
	my ($array_ref) = @_; 
	
#	print "\nmy array_ref input is $array_ref\n";
	
	my @array = @{$array_ref};
	#@array = qw/lalala heh 79jackobsin gah Mitch Chase Nelson/;
	
#	print "\nmy array input is @array\n";
	
	my @ordered_array;
	
	while(scalar(@array) > 0) {
	
		my $curr_min_elt_length = 10000000000; # something crazy big
		my $curr_min_elt;
		my $curr_min_elt_index;
	
		for(my $i = 0; $i < scalar(@array); $i++) {
			my $curr_elt = $array[$i];
			my $length_curr_name = length($curr_elt);
			
			if($length_curr_name < $curr_min_elt_length) {
				$curr_min_elt_length = $length_curr_name;
				$curr_min_elt = $curr_elt;
				$curr_min_elt_index = $i;
			}
			
		}
		
		# Found the name with min length; add to ordered array
		push(@ordered_array, $curr_min_elt);
		
		# Delete from source array
		splice(@array, $curr_min_elt_index, 1);
		
	}
	
#	print "\nmy array output is @ordered_array\n";
	
	return \@ordered_array;
	
}


