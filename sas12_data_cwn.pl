#! /usr/bin/perl
#
# HUXLEY: /usr/local/software/PERL/perl-5.26.0/bin/perl
# PC: /usr/bin/perl
# overlapping gene selection simulation
# sense-antisense overlap  123123
##########                   321321
use strict;   
use warnings; 
use diagnostics;
use autodie;
use Parallel::ForkManager;

#my $procs_per_node = 4;

####################################   INPUT HERE   ####################################
## my $dna1 = 'CTCTACTCCAAGCGCTTCGCCGTCTTCCTGTCCGAGGTCAGCGAAAGCCGTCTAAAGCAGCTCAATCTCAACCACGAGTGGACGCCCGAGAAGCTTCGACAGAAGCTGCAGCGCAATGCCGCGGGCCGGCTGGAGCTGGCCCTCTGCATGCTGCCGGGTCTGCCCGACACCGTCTTTGAGCTCAGTGAGGTGGAGTCACTCAGGCTGGAGGCCATCTGCGATATCACCTTCCCCCCGGGGCTGTCACAGCTGGTGCACTTGCAGGAGCTCAGCTTGCTCCACTCGCCCGCCAGGCTACCCTTCTCCTTGCAGGTCTTCCTGCGGGACCACCTGAAGGTGATGCGCGTCAAATGCGAGGAGCTCCGCGAGGTGCCGCTTTGGGTGTTTGGGCTGCGGGGCTTGGAGGAGCTGCACCTGGAGGGGCTTTTCCCCCAGGAGCTAGCTCGGGCAGCCACCCTGGAGAGCCTCCGGGAGCTGAAGCAGCTCAAGGTGTTGTCCCTCCGGAGCAACGCCGGGAAGGTGCCAGCCAGTGTGACCGACGTTGCTGGCCACCTGCAGAGGCTCAGCCTGCACAACGATGGGGCCCGTCTGGTTGCCCTGAACAGCCTCAAGAAGCTGGCGGCATTGCGGGAGCTGGAGCTGGTGGCCTGCGGGCTGGAGCGCATCCCCCATGCAGTGTTCAGCCTGGGTGCGCTGCAGGAACTTGACCTCAAGGACAACCACCTGCGCTCCAT';
## if($ARGV[0]) {
## 	$dna1 = $ARGV[0];
## }
## 
## my $dna2 = 'CTTTACTCCAAACGCTTTGCCGTCTTCCTCTCCGAGGTCAGCGAAAGCCGACTCAAGCAACTCAATCTCAACCATGAATGGACGCCAGAGAAACTGCGCCAGAAGCTCCAGCGCAATGGGAGAGGGCGGCTAGAGCTGGCACTGTGCATGCTGCCAGGGCTACCCGACACAGTGTTTGAGCTGATCGAGACAGAGTCTCTCAAGCTGGAGGCCATCTGTGATATCACTTTTCCTCCCACCCTCTCCCAGCTGGTTCACCTGGAAGAGCTCAGCCTGCTCCATTCCCCTGCCAAGCTGCCCTTCTCTTCCTTTGTTTTCTTGAGGGACAGACTGAAGGTGGTGAGGGTCAAGTGTGAGGAGCTTCGAGAAGTACCCCTGTGGGTGTTTGGGCTCCGGAGCCTGGAGGAGCTCCATCTGGAGGGCCTCTTCCCTCCAGAACTTGCTCGAGCCGCTAATCTGGAGAGCCTTAAAGAGCTAAAGCTCCTCAAGTGCCTCTCCTTGAGAAGCAATGCGGGCAAGATGCCCCCCAGTGTGACGGACGTGGCGGGCCACCTGCAGCAGCTCAGTCTCCACAACGATGGGACGAGGCTGCTGACCCTCAATGGACTCAAGAAGCTCACGGCGCTGAGAGAGTTAGAGCTGGTAGGCTGTGGGCTGGAGAGGATACCCCATGCTGTCTTCAGTTTGACTGCGCTTCAGGAGTTAGACCTGAAGGACAACCACCTCCGCTCCAT';
## if($ARGV[1]) {
## 	$dna2 = $ARGV[1];
## }
## 
## my $R=3.61; ######     This is transition transversion ratio =a/(2b)
## if($ARGV[2]) {
## 	$R = $ARGV[2];
## }
####################################   INTPUT ENDS   ####################################


####################################   ALT IN HERE   ####################################
my $fasta_file_name = $ARGV[0];
my $R = $ARGV[1]; ######     This is transition transversion ratio =a/(2b)
my $procs_per_node = $ARGV[2];

unless($fasta_file_name =~ /.fa/) { die "\n\n# FASTA file (argument 1) must contain .fa or .fasta extension. TERMINATED\n\n"; }
unless($R > 0) { die "\n\n# TS/TV ratio R (argument 2) must be provided, >0. TERMINATED\n\n"; }
unless($procs_per_node > 0) { die "\n\n# PROCS PER NOTE (argument 3) must be provided, >0. TERMINATED\n\n"; }

# Generate OUTPUT file name
my $results_file_name = "OLG_results.txt";
if($fasta_file_name =~/\.fa/) { 
	$results_file_name = $` . "_OLG_results.txt";
} elsif($fasta_file_name =~/\.txt/) { 
	$results_file_name = $` . "_OLG_results.txt";
} elsif($fasta_file_name =~/\.csv/) { 
	$results_file_name = $` . "_OLG_results.txt";
}


# Read in the group of sequences from the fasta file
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
####################################   ALT IN ENDS   ####################################

# We now have an array of all our FASTA sequences, @seqs_arr, and their headers, @headers_arr
# We also have the ts/tv ratio, $R

####################################   CALCULATION   ####################################
## my $size_of_set=1;
my $num_seqs = scalar(@seqs_arr);
my $size_of_set = ($num_seqs**2 - $num_seqs) / 2;

## my @dnds_2genes;

srand (time|$$);  # introduce a random seed each time 

## @dnds_2genes=calculate_dnds($dna1,$dna2,$R); # this is of size 20


my @comp_pair;
#my $comp_number = 0;


# PARALLELIZE COMPARISONS
mkdir("temp_dNdS_comparisons");
chdir("temp_dNdS_comparisons");

my $pm_poly = Parallel::ForkManager->new($procs_per_node);

for(my $seq_index = 0; $seq_index < $num_seqs; $seq_index++) {
	
	for(my $next_seq_index = ($seq_index + 1); $next_seq_index < $num_seqs; $next_seq_index++) {
		$pm_poly->start and next;
		
		# Get the time
		my $time1 = time;
#		$comp_number++;
		
		print "\nProcessing pair $seq_names_arr[$seq_index] vs. $seq_names_arr[$next_seq_index]\...\n";
		
		my $dna1 = $seqs_arr[$seq_index];
		my $dna2 = $seqs_arr[$next_seq_index];
#		print "$dna1\n$dna2\n";
		
		my @this_comp_dnds;
		
		# Skip comparisons with long gaps || KEEP IN MIND THIS MAY ACTUALLY BE ABOUT NUM. DIFFS, not GAPS
		if($dna1 =~ /----------------------------------------/ || $dna2 =~ /----------------------------------------/) { # 40 nts, i.e., longer than 39nt/3=13 codons
			print "One sequence contains a stretch of indels >13 codons; skipping comparison.\n"; # $seq_names_arr[$seq_index] has vs. $seq_names_arr[$next_seq_index]\...\n";
			for(my $i=0; $i<20; $i++) {	
				push(@this_comp_dnds,'*');
			}
		} else {
			my $dna1_header = $headers_arr[$seq_index];
			my $dna2_header = $headers_arr[$next_seq_index];
			
			@this_comp_dnds = &calculate_dnds($dna1,$dna2,$R); # this is NOW of size 20
		}
		
##		unless($dna1 =~ /----------------------/) { # 22 nts, i.e., longer than 21nt/3=7 codons
##			print "$seq_names_arr[$seq_index] has a stretch of indels >7 codons; skipping.\n"; #  vs. $seq_names_arr[$next_seq_index]\...\n";
##			if(! $dna2 =~ /----------------------/) { # 22 nts, i.e., longer than 21nt/3=7 codons
##				
##				my $dna1_header = $headers_arr[$seq_index];
##				my $dna2_header = $headers_arr[$next_seq_index];
##				
##				@this_comp_dnds = &calculate_dnds($dna1,$dna2,$R); # this is NOW of size 20
##				
##			} else {
##				
##			}
##		}

		open(THIS_COMP_TEMP_FILE,">>s$seq_index\-s$next_seq_index");
		
		foreach(@this_comp_dnds) {
			print THIS_COMP_TEMP_FILE "$_\t";
		}
		
		my $time2 = time;
		my $time_elapsed = $time2 - $time1;
		
		print THIS_COMP_TEMP_FILE "$time_elapsed";
		
		close THIS_COMP_TEMP_FILE;
		
		$pm_poly->finish; # special name
#		$comp_dNdS_ha{$comp_number} = \@this_comp_dnds;
	}
	
}

$pm_poly->wait_all_children; # special name, methods within module
		
# SWEEP UP RESULTS
my %comp_dNdS_ha;

my @temp_comp_FILES_arr = glob "*";
#@temp_comp_FILES_arr = sort {$a <=> $b} @temp_comp_FILES_arr;

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





open(OUTPUT,">>$results_file_name");
# Header
print OUTPUT "dna1\tdna2\tdnds1_by_ns\tdnds2_by_ns\tdnds1_by_ss\tdnds2_by_ss\tvar(dnn)\tvar(dsn)\tvar(dns)\tvar(dss)\tdnn\tdsn\tdns\tdss\tMnn\tMsn\tMns\tMss\tLnn\tLsn\tLns\tLss\ttime_used\n";

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

####################################   END CALCN.   ####################################

print "\n#### PROGRAM COMPLETE ####\n\n";

exit;


#########################################################################################
#########################################################################################
####################################                 ####################################
####################################   SUBROUTINES   ####################################
####################################                 ####################################
#########################################################################################
#########################################################################################

sub calculate_dnds {
	my ($dna1,$dna2,$R)=@_;
#	print "Determing number of sites...\n";
	my @N_S=N_S($dna1,$dna2,$R); # NN [0], SN [1], NS [2], and SS [3] SITES (UPPERcase)
	
#	print "Determing number of differences...\n";
	my @n_s=n_s($dna1,$dna2); # NN [0], SN [1], NS [2], and SS [3] DIFFS (lowercase)
	
#	print "Calculating and correcting dN/dS measures...\n";
	my @result;
	my $pn1n2=$n_s[0]/$N_S[0];
	my $ps1n2=$n_s[1]/$N_S[1];
	my $pn1s2=$n_s[2]/$N_S[2];
	
	my $ps1s2 = '*';
	if($N_S[3] > 0) {
		$ps1s2=$n_s[3]/$N_S[3];
	}
	
##	my $ps1s2=$n_s[3]/$N_S[3];
##	print "nn\tsn\tns\tss\n";
##	print $n_s[0],"\t",$n_s[1],"\t",$n_s[2],"\t",$n_s[3],"\n";
##	print "NN\tSN\tNS\tSS\n";
##	print $N_S[0],"\t",$N_S[1],"\t",$N_S[2],"\t",$N_S[3],"\n";
	
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
	#print $pn1n2,"\t",$ps1n2,"\t",$pn1s2,"\t",$ps1s2,"\n";
##	print $dn1n2,"\t",$ds1n2,"\t",$dn1s2,"\t",$ds1s2,"\n";
	my $pn_1=$n_s[4]/$N_S[4];
	my $ps_1=$n_s[5]/$N_S[5];
	my $pn_2=$n_s[6]/$N_S[6];
	my $ps_2=$n_s[7]/$N_S[7]; 
	
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
	
	my $origin_dnds1 = '*';
	if($ds_1 > 0) {
		$origin_dnds1=$dn_1/$ds_1;
	}
	
	my $origin_dnds2 = '*';
	if($ds_2 > 0) {
		$origin_dnds2=$dn_2/$ds_2;
	}
	
	my $dn_ratio=$dn_1/$dn_2;
	my $var_pn1n2=$pn1n2*(1-$pn1n2)/$N_S[0];
	my $var_ps1n2=$ps1n2*(1-$ps1n2)/$N_S[1];
	my $var_pn1s2=$pn1s2*(1-$pn1s2)/$N_S[2];
	my $var_ps1s2=$ps1s2*(1-$ps1s2)/$N_S[3];
	my $var_dn1n2=$var_pn1n2/(1-4/3*$pn1n2)**2;
	my $var_ds1n2=$var_ps1n2/(1-4/3*$ps1n2)**2;
	my $var_dn1s2=$var_pn1s2/(1-4/3*$pn1s2)**2;
	my $var_ds1s2=$var_ps1s2/(1-4/3*$ps1s2)**2;
	my $dnds1_1=0;
	my $dnds1_2=0;
	my $dnds2_1=0;
	my $dnds2_2=0;
	if ($ds1n2!=0){
		$dnds1_1=$dn1n2/$ds1n2;
	}
	if ($dn1s2>0)	{
		$dnds2_1=$dn1n2/$dn1s2;
	}
	if ($ds1s2>0){
		$dnds1_2=$dn1s2/$ds1s2;
		$dnds2_2=$ds1n2/$ds1s2;
	}
	if (($dnds1_1!=0)and($dnds1_2!=0)){
		#print "dnds1 is: ",($dnds1_1+$dnds1_2)/2,"\n";
	} elsif ($dnds1_1!=0){
##		print "dnds1 is: ",$dnds1_1,"\n";
	} elsif ($dnds1_2!=0){
##		print "dnds1 is: ",$dnds1_2,"\n";
	} else{
##		print 'Sorry, I could not get the ratio of dnds1';
	}
	if (($dnds2_1!=0)and($dnds2_2!=0)){
		#print "dnds2 is: ",($dnds2_1+$dnds2_2)/2,"\n";
	} elsif ($dnds2_1!=0){
##		print "dnds2 is: ",$dnds2_1,"\n";
	} elsif ($dnds2_2!=0){
##		print "dnds2 is: ",$dnds2_2,"\n";
	} else{
##		print 'Sorry, I could not get the ratio of dnds2';
	}
	
	my $dNN_by_dSN = '*';
	if($ds1n2 > 0) {
		$dNN_by_dSN = $dn1n2/$ds1n2;
	}
	
	my $dNN_by_dNS = '*';
	if($dn1s2 > 0) {
		$dNN_by_dSN = $dn1n2/$dn1s2;
	}
	
	my $dNS_by_dSS = '*';
	if($ds1s2 > 0) {
		$dNS_by_dSS = $dn1s2/$ds1s2;
	}
	
	my $dSN_by_dSS = '*';
	if($ds1s2 > 0) {
		$dSN_by_dSS = $ds1n2/$ds1s2;
	}
	
	@result=($dNN_by_dSN,$dNN_by_dNS,$dNS_by_dSS,$dSN_by_dSS,$var_dn1n2,$var_ds1n2,$var_dn1s2,$var_ds1s2,$dn1n2,$ds1n2,$dn1s2,$ds1s2,$n_s[0],$n_s[1],$n_s[2],$n_s[3],$N_S[0],$N_S[1],$N_S[2],$N_S[3]);
	return @result;
}#$dn1s2/$ds1s2


sub N_S {
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


sub N_S_dna{  
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

# DETERMINE THE NUMBER OF DIFFERENCES
sub n_s {
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
	
	while ($i < ($length-2)) {  
		@region_change=switch_region($region1,$region2,$i,$dna1,$dna2);
		$region1=$region_change[0];
		$region2=$region_change[1];
		$i=$region_change[2]; #
		
#		print "Determining partial result for i $i\...\n";
		@partial_result=pathway($region1,$region2,$i,$pathway_number,$n1n2_pathway,$s1n2_pathway,$n1s2_pathway,$s1s2_pathway,$n1n2_pathway,$s1n2_pathway,$n1s2_pathway,$s1s2_pathway);
		
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

sub pathway{  
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

sub switch_region {  ###
	my ($region1,$region2,$i,$dna1,$dna2)=@_; #$i is where region with mutation starts in $dna1
	my @result;
	my $change_number;
	my $length=length($dna1);  
	my $new_region1;
	my $new_region2;
	my $j;
	my $new_change_number;
	my @compare_result;
	$change_number=change_number($region1,$region2);
	while (($change_number==0) & ($i<($length-3))) {##
		if (($i%3)==0) {
			$i=$i+2;
			if ($i<($length-2)){
				$region1=substr($dna1,$i,3);
				$region2=substr($dna2,$i,3);
			}
		} elsif (($i%3)==2) {
			$i=$i+1;
			if ($i<($length-2)){
				$region1=substr($dna1,$i,3);
				$region2=substr($dna2,$i,3);
			}
		}
		$change_number=change_number($region1,$region2);
	}
	if ($change_number>0){
		$j=$i+3;
		if ($i%3==0){ 
			if ($j<($length-1)){ ####
				$new_region1=substr($dna1,$j,2);
				$new_region2=substr($dna2,$j,2);
				$new_change_number=change_number($new_region1,$new_region2);
				@compare_result=compare_new_region($j,$new_region1,$new_region2,$dna1,$dna2);
			}
		}
		if ($i%3==2) { #
			if ($j<$length){########
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

sub compare_new_region { 
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
			if ($change> 0) {
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

sub change_number {
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

sub change_place {  ###
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

sub compare_dna1_dna2 {
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

sub possible_change {
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
				} else{
					$S1S2=$S1S2+1/(2+2*$R);
				}
			} elsif($new_aa1 ne '_') {
				if ($original_base eq 'T'){
					$S1N2=$S1N2+$R/(1+$R);
				} else{
					$S1N2=$S1N2+1/(2+2*$R);
				}
			}
		} elsif($new_aa1 ne '_') {
			if ($aa2 eq $new_aa2) {
				if ($original_base eq 'T'){
					$N1S2=$N1S2+$R/(1+$R);
				} else{
					$N1S2=$N1S2+1/(2+2*$R);
				}
			} elsif($new_aa2 ne '_') {
				if ($original_base eq 'T'){
					$N1N2=$N1N2+$R/(1+$R);
				} else{
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
				if ($original_base eq 'A'){
					$S1S2=$S1S2+$R/(1+$R);
				} else{
					$S1S2=$S1S2+1/(2+2*$R);
				}
			} elsif($new_aa1 ne '_') {
				if ($original_base eq 'A'){
					$S1N2=$S1N2+$R/(1+$R);
				} else{
					$S1N2=$S1N2+1/(2+2*$R);
				}
			}
		} elsif($new_aa1 ne '_') {
			if ($aa2 eq $new_aa2) {
				if ($original_base eq 'A'){
					$N1S2=$N1S2+$R/(1+$R);
				} else{
					$N1S2=$N1S2+1/(2+2*$R);
				}
			} elsif($new_aa2 ne '_') {
				if ($original_base eq 'A'){
					$N1N2=$N1N2+$R/(1+$R);
				} else{
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
				if ($original_base eq 'C'){
					$S1S2=$S1S2+$R/(1+$R);
				} else{
					$S1S2=$S1S2+1/(2+2*$R);
				}
			} elsif($new_aa1 ne '_') {
				if ($original_base eq 'C'){
					$S1N2=$S1N2+$R/(1+$R);
				} else{
					$S1N2=$S1N2+1/(2+2*$R);
				}
			}
		} elsif($new_aa1 ne '_') {
			if ($aa2 eq $new_aa2) {
				if ($original_base eq 'C'){
					$N1S2=$N1S2+$R/(1+$R);
				} else{
					$N1S2=$N1S2+1/(2+2*$R);
				}
			} elsif($new_aa2 ne '_') {
				if ($original_base eq 'C'){
					$N1N2=$N1N2+$R/(1+$R);
				} else{
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
				if ($original_base eq 'G'){
					$S1S2=$S1S2+$R/(1+$R);
				} else{
					$S1S2=$S1S2+1/(2+2*$R);
				}
			} elsif($new_aa1 ne '_') {
				if ($original_base eq 'G'){
					$S1N2=$S1N2+$R/(1+$R);
				} else{
					$S1N2=$S1N2+1/(2+2*$R);
				}
			}
		} elsif($new_aa1 ne '_') {
			if ($aa2 eq $new_aa2) {
				if ($original_base eq 'G'){
					$N1S2=$N1S2+$R/(1+$R);
				} else{
					$N1S2=$N1S2+1/(2+2*$R);
				}
			} elsif($new_aa2 ne '_') {
				if ($original_base eq 'G'){
					$N1N2=$N1N2+$R/(1+$R);
				} else{
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


sub make_overlap_DNA_set {
	my($number_of_codons,$size_of_set) = @_;     
	my $dna;
	my @set;
	for (my $i = 0; $i < $size_of_set ; ++$i) {
		$dna = make_overlap_DNA ( $number_of_codons ); 
		push(@set, $dna );
	}
	return @set;
}


sub make_mutation_set{
	my ($distance,$R,@random_DNA)=@_;
	my @set;
	foreach my $dna (@random_DNA){  
		$dna=mutation_and_selection($dna,$distance,$R);  
		push(@set,$dna);
	}
	return @set;
}


sub randomposition{
	my($string)= @_;
	return (int(rand(length($string)-4))+2);  
}


sub make_overlap_DNA{ #
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


sub randomnucleotide{    
	my(@nucleotides)=('A','C','G','T');
	return randomelement(@nucleotides);
}


sub random_codon{
	my(@codons)=('TCA','TCC','TCG','TCT','TTC','TTT','TTA','TTG','TAC','TAT','TGC','TGT','TGG','CTA','CTC','CTG','CTT','CCA','CCC','CCG','CCT','CAC','CAT','CAA','CAG','CGA','CGC','CGG','CGT','ATA','ATC','ATT','ATG','ACA','ACC','ACG','ACT','AAC','AAT','AAA','AAG','AGC','AGT','AGA','AGG','GTA','GTC','GTG','GTT','GCA','GCC','GCG','GCT','GAC','GAT','GAA','GAG','GGA','GGC','GGG','GGT');
	return randomelement(@codons);
}


sub randomelement{
	my(@array)= @_;
	return $array[rand @array];
}


sub codon2aa {
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



