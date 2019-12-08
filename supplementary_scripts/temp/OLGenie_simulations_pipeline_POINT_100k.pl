#! /nas3/cnelson/bin/anaconda2/bin/perl
# /usr/bin/perl
# /usr/bin/perl

# Call from any empty directory, like:
# /Users/cwnelson88/Desktop/OLG/simulations/SIMULATIONS

# HERE WE HAVE CHANGED PARALLELISM TO THE BOOTSTRAP; ONLY ONE REPLICATE EACH!

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
my $SIMULATION_SCRIPT_DIR = '/nas3/cnelson/bin/'; #'/Users/cwnelson88/Desktop/OLG/simulations/ss13';
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
my $R = 0.5; # 0.5 2

# FIXED
my @dnds_1_values = (0.1, 0.5, 1, 1.5, 2);
my @dnds_2_values = (0.1, 0.5, 1, 1.5, 2);
my @number_of_codons = (100000); 
# (50, 90, 100, 200); # INSERT EMPIRICAL VALUES FROM CONTROLS HERE
# (10, 88, 193, 502, 4384) # min0 quartiles
my @number_of_seqs = (1024);
#(2, 13, 47, 176, 12206); # min0 quartiles
my $num_replicates = 1;
my $num_genes = 1; # 100; # scalar(@number_of_codons)	

# Form all possible dnds combos, to be parallelized
my @dnds_combos;
foreach my $dnds1 (@dnds_1_values) {
	foreach my $dnds2 (@dnds_2_values) {
		push (@dnds_combos, "$dnds1\,$dnds2");
	}
}


foreach my $frame_type (@frame_types) {
	# Create the directory for this frame type
	mkdir("$original_directory\/$frame_type");
	
	foreach my $distance (@distances) {
	
		mkdir("$original_directory\/$frame_type\/distance_$distance");
		
		my $pm_poly = Parallel::ForkManager->new($NCPUS); # PARALLEL
		
		foreach my $dnds_combo (@dnds_combos) {
		
			$pm_poly->start and next; # PARALLEL
		
			my @dnds_values = split /,/, $dnds_combo;
			my $dnds1 = $dnds_values[0];
			my $dnds2 = $dnds_values[1];
		
			# Create the directory for this dnds combination
			mkdir("$original_directory\/$frame_type\/distance_$distance\/dnds1_$dnds1\_dnds2_$dnds2");
			
			# Create and analyze ONE replicate simulated datasets for this combination
			for (my $i = 1; $i <= $num_replicates; $i++) {
				mkdir("$original_directory\/$frame_type\/distance_$distance\/dnds1_$dnds1\_dnds2_$dnds2\/replicate_$i");
				chdir("$original_directory\/$frame_type\/distance_$distance\/dnds1_$dnds1\_dnds2_$dnds2\/replicate_$i");
				
				for (my $gene_num = 1; $gene_num <= $num_genes; $gene_num++) {
					
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
					
				}
			}
		
			$pm_poly->finish; # PARALLEL
			
		}
		
		$pm_poly->wait_all_children; # PARALLEL
		
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

