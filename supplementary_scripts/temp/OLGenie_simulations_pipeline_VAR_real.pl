#! /nas3/cnelson/bin/anaconda2/bin/perl
# /usr/bin/perl
# /usr/local/software/PERL/perl-5.26.0/bin/perl
# /usr/bin/perl

# Call from any empty directory, like:
# /Users/cwnelson88/Desktop/OLG/simulations/SIMULATIONS

#########################################################################################
# EXAMPLE CALL:
#########################################################################################
# ./OLGenie_controls_pipeline.pl
#########################################################################################

# Copyright (C) 2019 Chase W. Nelson
# DATE CREATED: October 2019
# AUTHOR: Chase W. Nelson
# CONTACT1: cnelson@amnh.org
# AFFILIATION: Sackler Institute for Comparative Genomics, American Museum of Natural 
#	History, New York, NY 10024, USA

# ACKNOWLEDGMENTS: written by C.W.N. with support from a Gerstner Scholars Fellowship from
# 	the Gerstner Family Foundation at the American Museum of Natural History, New York.

use strict;
use Parallel::ForkManager;
#use warnings;
#use List::Util qw(max sum);
#use Getopt::Long;

# Get the time
my $time1 = time;
my $local_time1 = localtime;

#print "\nAnalysis initiated at local time $local_time1\n";
my $original_directory = `pwd`;
chomp($original_directory);
#print "The original working directory is $original_directory\n";

##########################################################################################
# SET GLOBAL PARAMETERS 
my $SIMULATION_SCRIPT_DIR = '/nas3/cnelson/bin'; #'/Users/cwnelson88/Desktop/OLG/simulations/ss13';
my $OLGENIE = '/nas3/cnelson/bin/OLGenie.pl'; # '/Users/cwnelson88/scripts_NGS/overlapgenie_github/OLGenie.pl';
my $OLGENIE_BOOTSTRAP = '/nas3/cnelson/bin/OLGenie_bootstrap.R'; #'/Users/cwnelson88/scripts_NGS/overlapgenie_github_bucket/OLGenie_process_codons_batch.R'; # '/nas3/cnelson/OLG/simulations/OLGenie_process_codons_batch.R'; # or '/Users/cwnelson88/scripts_NGS/overlapgenie_github/OLGenie_process_codons_batch.R';
my $NBOOTSTRAPS = 100;
my $NCPUS = 4; # 8 on Mac, or 60 on HPC

my @frame_types = qw/ss13 sas11 sas12 sas13/; #qw/ss13/; #
my $distance = 0.05;
my $R = 0.5;

if(@ARGV == 2) {
	$NBOOTSTRAPS = $ARGV[0];
	$NCPUS = $ARGV[1];
} elsif(@ARGV == 3) {
	$NBOOTSTRAPS = $ARGV[0];
	$NCPUS = $ARGV[1];
	@frame_types = ($ARGV[2]);
} elsif(@ARGV == 4) {
	$NBOOTSTRAPS = $ARGV[0];
	$NCPUS = $ARGV[1];
	@frame_types = ($ARGV[2]);
	$distance = $ARGV[3];
} elsif(@ARGV == 5) {
	$NBOOTSTRAPS = $ARGV[0];
	$NCPUS = $ARGV[1];
	@frame_types = ($ARGV[2]);
	$distance = $ARGV[3];
	$R = $ARGV[4];
}


##########################################################################################
# Create a directory for each frame type; create 100 replicate datasets

### EXPLORE WHAT REPRODUCES CONTROL SET
my @distances = ($distance); # (0.1, 0.25, 0.5); # 0.005
#(0.0009892, 0.0293707, 0.0467106, 0.0920861, 0.2611777) # min0 quartiles

# FIXED
my @dnds_1_values = (0.1, 0.5, 1, 1.5, 2);
my @dnds_2_values = (0.1, 0.5, 1, 1.5, 2);
my @number_of_codons = (8, 16, 17, 18, 22, 23, 24, 28, 29, 33, 33, 33, 33, 36, 37, 37, 37, 41, 42, 42, 43, 44, 44, 45, 47, 50, 51, 51, 52, 53, 53, 53, 56, 57, 58, 62, 62, 62, 62, 64, 64, 66, 66, 66, 68, 68, 70, 70, 71, 72, 73, 73, 75, 77, 78, 80, 82, 84, 84, 84, 84, 85, 86, 87, 90, 91, 92, 96, 97, 98, 98, 98, 102, 104, 112, 114, 115, 115, 117, 119, 119, 120, 120, 122, 122, 126, 127, 127, 127, 130, 131, 132, 134, 135, 137, 138, 139, 140, 145, 148, 148, 148, 149, 150, 152, 154, 154, 159, 164, 166, 172, 174, 174, 176, 177, 179, 186, 187, 188, 189, 192, 194, 194, 197, 199, 202, 202, 202, 206, 206, 207, 213, 220, 222, 223, 228, 230, 231, 235, 235, 256, 282, 311, 313, 318, 318, 320, 322, 327, 334, 334, 347, 356, 366, 369, 370, 376, 376, 377, 387, 392, 395, 412, 417, 420, 420, 421, 421, 426, 457, 468, 473, 487, 493, 495, 503, 510, 523, 524, 524, 536, 548, 549, 549, 549, 550, 550, 564, 573, 574, 587, 616, 616, 620, 621, 622, 625, 630, 649, 666, 666, 696, 722, 764, 799, 817, 823, 839, 860, 886, 886, 892, 915, 980, 1046, 1118, 1216, 1232, 1246, 1472, 1472, 1581, 1607, 1677, 1692, 2107, 2108, 2146, 2182, 2182, 2227, 2688, 2791, 4369); 
# (50, 90, 100, 200); # INSERT EMPIRICAL VALUES FROM CONTROLS HERE
# (10, 88, 193, 502, 4384) # min0 quartiles
my @number_of_seqs = (6, 6, 6, 6, 6, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 10, 10, 10, 11, 11, 11, 11, 11, 12, 13, 13, 13, 13, 14, 14, 14, 14, 14, 15, 15, 15, 15, 17, 18, 18, 19, 19, 20, 22, 22, 22, 23, 23, 24, 24, 24, 24, 24, 25, 25, 27, 27, 27, 28, 28, 29, 29, 29, 30, 31, 31, 31, 32, 32, 32, 32, 33, 33, 34, 34, 35, 35, 35, 35, 36, 37, 38, 38, 39, 39, 40, 40, 41, 42, 43, 43, 43, 44, 44, 45, 47, 47, 47, 48, 48, 48, 49, 49, 50, 51, 51, 54, 55, 57, 63, 64, 65, 66, 66, 67, 70, 70, 71, 71, 71, 72, 74, 75, 75, 76, 76, 77, 80, 81, 81, 82, 84, 85, 86, 86, 86, 87, 89, 90, 95, 98, 100, 107, 119, 120, 120, 126, 128, 133, 142, 144, 144, 144, 146, 147, 149, 149, 152, 155, 155, 159, 162, 179, 185, 194, 207, 207, 222, 223, 225, 230, 237, 240, 240, 267, 273, 275, 280, 290, 295, 301, 302, 302, 323, 332, 369, 381, 509, 509, 522, 661, 749, 749, 827, 833, 1064, 1107, 1177, 1313, 1360, 1619, 1872, 2747, 3476, 3706, 3748, 3748, 3832, 4778, 5208, 5663, 5671, 5671, 5833, 5833, 6012, 6161, 6812, 7247, 7247, 7306, 7744, 8149, 8820, 8944, 9167, 9954, 11026, 12206);
#(2, 13, 47, 176, 12206); # min0 quartiles
my $num_replicates = 1;
my $num_genes = 234; # 100; # scalar(@number_of_codons)	

foreach my $frame_type (@frame_types) {
	# Create the directory for this frame type
	mkdir("$original_directory\/$frame_type");
	
	foreach my $distance (@distances) {
	
		mkdir("$original_directory\/$frame_type\/distance_$distance");
		
		foreach my $dnds1 (@dnds_1_values) {
			
			foreach my $dnds2 (@dnds_2_values) {
				# Create the directory for this dnds combination
				mkdir("$original_directory\/$frame_type\/distance_$distance\/dnds1_$dnds1\_dnds2_$dnds2");
				
				# Create and analyze 100 replicate simulated datasets for each combination
				for (my $i = 1; $i <= $num_replicates; $i++) {
					mkdir("$original_directory\/$frame_type\/distance_$distance\/dnds1_$dnds1\_dnds2_$dnds2\/replicate_$i");
					chdir("$original_directory\/$frame_type\/distance_$distance\/dnds1_$dnds1\_dnds2_$dnds2\/replicate_$i");
					
					my $pm_poly = Parallel::ForkManager->new($NCPUS);
					
					for (my $gene_num = 1; $gene_num <= $num_genes; $gene_num++) {
						
						$pm_poly->start and next; # this is IT
						
						my $random_seed = srand(time ^ $$ ^ unpack "%32L*", `ps wwaxl | gzip`); # (Programming Perl, p. 955)
						#print "\nRANDOM_SEED: $random_seed\n";
						
						# Select the number of codons (sequence length) and number of sequences randomly.
						my $num_codons = $number_of_codons[int(rand(scalar(@number_of_codons)))];
						my $num_seqs = $number_of_seqs[int(rand(scalar(@number_of_seqs)))];
						
						mkdir("$original_directory\/$frame_type\/distance_$distance\/dnds1_$dnds1\_dnds2_$dnds2\/replicate_$i\/gene_$gene_num\_len_$num_codons\_seqs_$num_seqs");
						chdir("$original_directory\/$frame_type\/distance_$distance\/dnds1_$dnds1\_dnds2_$dnds2\/replicate_$i\/gene_$gene_num\_len_$num_codons\_seqs_$num_seqs");
						
						# LOGFILE FOR TIMES
						open(PIPELINE_LOG, ">$original_directory\/$frame_type\/distance_$distance\/dnds1_$dnds1\_dnds2_$dnds2\/replicate_$i\/gene_$gene_num\_len_$num_codons\_seqs_$num_seqs\/PIPELINE.log");
						
						# Simulate data
						my $outfile_name = "$frame_type\_$distance\_dnds1_$dnds1\_dnds2_$dnds2\_gene_$gene_num\_len_$num_codons\_seqs_$num_seqs\.txt";
						my $outfile = "$original_directory\/$frame_type\/distance_$distance\/dnds1_$dnds1\_dnds2_$dnds2\/replicate_$i\/gene_$gene_num\_len_$num_codons\_seqs_$num_seqs\/$outfile_name";
						my $simulation_logfile = "$original_directory\/$frame_type\/distance_$distance\/dnds1_$dnds1\_dnds2_$dnds2\/replicate_$i\/gene_$gene_num\_len_$num_codons\_seqs_$num_seqs\/SIMULATION.log";
						
						my $simulation_command = "$SIMULATION_SCRIPT_DIR\/simulation_$frame_type\.pl --dnds1=$dnds1 --dnds2=$dnds2 --num_codons=$num_codons --num_seqs=$num_seqs --distance_actual=$distance --R=$R --outfile=$outfile > $simulation_logfile";
						#print "### simulation_command:\n$simulation_command\n";
						my $time1 = time;
						`$simulation_command`;
						my $time2 = time;
						my $time_diff = ($time2 - $time1);
						print PIPELINE_LOG "simulation_command=$simulation_command\, time=$time_diff\n";
						
						# Run OLGenie
						my $OLGenie_command = "$OLGENIE --fasta_file=$outfile --frame=$frame_type --verbose > $original_directory\/$frame_type\/distance_$distance\/dnds1_$dnds1\_dnds2_$dnds2\/replicate_$i\/gene_$gene_num\_len_$num_codons\_seqs_$num_seqs\/OLGenie.out";
						#print "### OLGenie_command:\n$OLGenie_command\n";
						$time1 = time;
						`$OLGenie_command`;
						$time2 = time;
						$time_diff = ($time2 - $time1);
						print PIPELINE_LOG "OLGenie_command=$OLGenie_command\, time=$time_diff\n";
						
						# Perform bootstrapping
						#my $outline = "$frame_type\t$dnds1\t$dnds2\t$i\t$gene_num\t$num_codons\t$num_seqs\t$random_seed";
						my $METADATA = "$frame_type\_$distance\_$dnds1\_$dnds2\_$i\_$gene_num\_$num_codons\_$num_seqs\_$random_seed";
						
						# Now call R; bootstrap with 1,000 replicates
						my $bootstrap_command = "Rscript --quiet --no-restore --no-save $OLGENIE_BOOTSTRAP $original_directory\/$frame_type\/distance_$distance\/dnds1_$dnds1\_dnds2_$dnds2\/replicate_$i\/gene_$gene_num\_len_$num_codons\_seqs_$num_seqs\/OLGenie_codon_results.txt 0 $NBOOTSTRAPS 1 $METADATA > OLGenie_bootstrap_results.txt 2> OLGENIE_BOOTSTRAP.out && touch OLGENIE_BOOTSTRAP.done"; #$NCPUS";
						$time1 = time;
						`$bootstrap_command`;
						$time2 = time;
						$time_diff = ($time2 - $time1);
						print PIPELINE_LOG "bootstrap_command=$bootstrap_command\, time=$time_diff\n";
						
						close PIPELINE_LOG;
						
#						my $this_parameters_results = `$bootstrap_command`;
#						
#						#print "this_parameters_results=$this_parameters_results\n";
#						my @this_parameters_results = split("\n", $this_parameters_results);
#						
#						open(OLGENIE_GENE_RESULTS, ">$original_directory\/$frame_type\/distance_$distance\/dnds1_$dnds1\_dnds2_$dnds2\/replicate_$i\/gene_$gene_num\_len_$num_codons\_seqs_$num_seqs\/OLGenie_gene_results.txt");
#						foreach my $result (@this_parameters_results) {
#							#push(@RESULTS, "$outline\t$_");
#							print OLGENIE_GENE_RESULTS "$result\n"; # "$outline\t$result\n";
#							print "$result\n";
#						}
#						close OLGENIE_GENE_RESULTS;
						
						$pm_poly->finish; # special name
						
					}
					
					$pm_poly->wait_all_children; # special name, methods within module
					
				}
			}
		}
	}
}

print "\n";

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
sub roundup { # to next integer
	my ($value) = @_;
    
    if(int($value) == $value) { # already an integer
    	return $value;
    } else {
    	my $rounded_up_value = int($value+1);
    	return $rounded_up_value
    }
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
	
	print "Analysis completed at local time $local_time2. The process took $time_diff_rounded secs, i.e., ".
			"$whole_mins_elapsed mins and $secs_remaining_rounded secs\n";

	print "\n################################################################################".
		"\n##                      Analysis completed successfully.                      ##".
		"\n##                Please find results in the working directory.               ##\n".
		"################################################################################".
		"\n\n\n"; 
}

