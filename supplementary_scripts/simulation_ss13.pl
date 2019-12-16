#! /usr/bin/perl

############################################################################################################
### Supplementary Material for Nelson CW, Ardern Z, Wei X, "OLGenie: Estimating Natural Selection to Predict Functional Overlapping Genes"
### Modified from Wei X, Zhang J, "A Simple Method for Estimating the Strength of Natural Selection on Overlapping Genes", doi: 10.1093/gbe/evu294
### overlapping gene selection simulation, frame ss13

############################################################################################################
### EXAMPLE:
# ./simulation_ss13.pl --dnds1=0.1 --dnds2=0.1 --num_codons=100 --num_seqs=100 --distance=0.05 --R=1 --outfile=sim_output.fasta
############################################################################################################

# second base in ORF1 is first base in ORF2
#select codon, and use Jukes Canter Model to correct, add transition transversion ratio
use strict;   
use warnings; 
use diagnostics;
use autodie;
use List::Util 'shuffle';
use Getopt::Long;

#srand (time|$$);  # introduce a random seed each time
my $random_seed = srand(time ^ $$ ^ unpack "%32L*", `ps wwaxl | gzip`); # (Programming Perl, p. 955)

# Initialize input variables
my $dnds1;
my $dnds2;
my $num_codons; # alignment LENGTH
my $num_seqs;
my $distance = 'NA';
my $distance_actual = 'NA';
my $R;
my $outfile;

my $die_message = "\n### OLGenie **ss13** simulation\n".
		"### PLEASE PROVIDE THE FOLLOWING 7 ARGUMENTS:\n".
		"# --dnds1=\<reference gene dN/dS\>\n" . 
		"# --dnds2=\<alternate gene dN/dS\>\n" . 
		"# --num_codons=\<number of codons\>\n" . 
		"# --num_seqs=\<number of sequences\> (max 1,024)\n" . 
		"# --distance=\<mean pairwise differences per site\> (original script)\n" . 
		"# --distance_actual=\<mean pairwise differences per site\> (a value we will attempt to achieve)\n" . 
		"# --R=\<transition/transversion ratio\>\n" . 
		"# --outfile=\<output file name or path\>\n\n";

# Get user input, if given. If a Boolean argument is passed, its value is 1; else undef
GetOptions( "dnds1=f" => \$dnds1,
			"dnds2=f" => \$dnds2,
			"num_codons=i" => \$num_codons,
			"num_seqs=i" => \$num_seqs,
			"distance=f" => \$distance,
			"distance_actual=f" => \$distance_actual,
			"R=f" => \$R,
			"outfile=s" => \$outfile
			 )
			
			or die $die_message;
			# If an argument is called as a flag, its value is 0; if not called, it's null

# User-defined input
if(! $dnds1 || ! $dnds2 || ! $num_codons || ! $num_seqs || 
	(! $distance && ! $distance_actual ) || ! $R || ! $outfile) { 
	die $die_message;
}

if($distance ne 'NA' && $distance_actual ne 'NA') {
	print "\n### WARNING: both --distance and --distance_actual were used. The second takes priority and overrides the first.\n";
	$distance = 'NA';
}
			
print "### PARAMETERS:\n";
print "random_seed=$random_seed\n";
print "dnds1=$dnds1\n";
print "dnds2=$dnds2\n";
print "num_codons=$num_codons\n";
print "num_seqs=$num_seqs\n";
print "distance=$distance\n";
print "distance_actual=$distance_actual\n";
print "R=$R\n";
print "outfile=$outfile\n";
print "random_seed=$random_seed\n";

# Initialize %distances_expected if we need to achieve a certain distance
if($distance_actual =~ /[0-9]/ && $distance_actual ne 'NA') {
	my %distances_expected = (%{&initialize_distances_expected_hash()});
	#$distances_expected{'ss13'}->{0.5}->{2}->{2} = 0.12354994092148716;
	#print "distances_expected=" . $distances_expected{'ss13'}->{0.5}->{2}->{2} . "\n";
	
	# Correct distance to obtain the actual distance outcome desired
	
	my $R_for_correction;
	my $dnds1_for_correction;
	my $dnds2_for_correction;
	
	# Get closest R
	if($R < 0.5) {
		print "\n### WARNING: R=$R was specified with distance_actual=$distance_actual\.\n" .
			"### However, distance correction is only available for R=\(0.5\|1\|2\)\.\n" .
			"### The correction for R=0.5 will be used as an approximation\.\n";
		$R_for_correction = 0.5;
	} elsif($R == 0.5) {
		$R_for_correction = 0.5;
	} elsif($R < 0.75) {
		print "\n### WARNING: R=$R was specified with distance_actual=$distance_actual\.\n" .
			"### However, distance correction is only available for R=\(0.5\|1\|2\)\.\n" .
			"### The correction for R=0.5 will be used as an approximation\.\n";
		$R_for_correction = 0.5;
	} elsif($R < 1) {
		print "\n### WARNING: R=$R was specified with distance_actual=$distance_actual\.\n" .
			"### However, distance correction is only available for R=\(0.5\|1\|2\)\.\n" .
			"### The correction for R=1 will be used as an approximation\.\n";
		$R_for_correction = 1;
	} elsif($R == 1) {
		$R_for_correction = 1;
	} elsif($R < 1.5) {
		print "\n### WARNING: R=$R was specified with distance_actual=$distance_actual\.\n" .
			"### However, distance correction is only available for R=\(0.5\|1\|2\)\.\n" .
			"### The correction for R=1 will be used as an approximation\.\n";
		$R_for_correction = 1;
	} elsif($R < 2) {
		print "\n### WARNING: R=$R was specified with distance_actual=$distance_actual\.\n" .
			"### However, distance correction is only available for R=\(0.5\|1\|2\)\.\n" .
			"### The correction for R=2 will be used as an approximation\.\n";
		$R_for_correction = 2;
	} elsif($R == 2) {
		$R_for_correction = 2;
	} else {
		print "\n### WARNING: R=$R was specified with distance_actual=$distance_actual\.\n" .
			"### However, distance correction is only available for R=\(0.5\|1\|2\)\.\n" .
			"### The correction for R=2 will be used as an approximation\.\n";
		$R_for_correction = 2;
	} 
	
	# Get closest dnds1
	if($dnds1 < 0.1) {
		print "\n### WARNING: dnds1=$dnds1 was specified with distance_actual=$distance_actual\.\n" .
			"### However, distance correction is only available for dnds1=\(0.1\|0.5\|1\|1.5\|2\)\.\n" .
			"### The correction for dnds1=0.1 will be used as an approximation\.\n";
		$dnds1_for_correction = 0.1;
	} elsif($dnds1 == 0.1) {
		$dnds1_for_correction = 0.1;
	} elsif($dnds1 < 0.3) {
		print "\n### WARNING: dnds1=$dnds1 was specified with distance_actual=$distance_actual\.\n" .
			"### However, distance correction is only available for dnds1=\(0.1\|0.5\|1\|1.5\|2\)\.\n" .
			"### The correction for dnds1=0.1 will be used as an approximation\.\n";
		$dnds1_for_correction = 0.1;
	} elsif($dnds1 < 0.5) {
		print "\n### WARNING: dnds1=$dnds1 was specified with distance_actual=$distance_actual\.\n" .
			"### However, distance correction is only available for dnds1=\(0.1\|0.5\|1\|1.5\|2\)\.\n" .
			"### The correction for dnds1=0.5 will be used as an approximation\.\n";
		$dnds1_for_correction = 0.5;
	} elsif($dnds1 == 0.5) {
		$dnds1_for_correction = 0.5;
	} elsif($dnds1 < 0.75) {
		print "\n### WARNING: dnds1=$dnds1 was specified with distance_actual=$distance_actual\.\n" .
			"### However, distance correction is only available for dnds1=\(0.1\|0.5\|1\|1.5\|2\)\.\n" .
			"### The correction for dnds1=0.5 will be used as an approximation\.\n";
		$dnds1_for_correction = 0.5;
	} elsif($dnds1 < 1) {
		print "\n### WARNING: dnds1=$dnds1 was specified with distance_actual=$distance_actual\.\n" .
			"### However, distance correction is only available for dnds1=\(0.1\|0.5\|1\|1.5\|2\)\.\n" .
			"### The correction for dnds1=1 will be used as an approximation\.\n";
		$dnds1_for_correction = 1;
	} elsif($dnds1 == 1) {
		$dnds1_for_correction = 1;
	} elsif($dnds1 < 1.25) {
		print "\n### WARNING: dnds1=$dnds1 was specified with distance_actual=$distance_actual\.\n" .
			"### However, distance correction is only available for dnds1=\(0.1\|0.5\|1\|1.5\|2\)\.\n" .
			"### The correction for dnds1=1 will be used as an approximation\.\n";
		$dnds1_for_correction = 1;
	} elsif($dnds1 < 1.5) {
		print "\n### WARNING: dnds1=$dnds1 was specified with distance_actual=$distance_actual\.\n" .
			"### However, distance correction is only available for dnds1=\(0.1\|0.5\|1\|1.5\|2\)\.\n" .
			"### The correction for dnds1=1.5 will be used as an approximation\.\n";
		$dnds1_for_correction = 1.5;
	} elsif($dnds1 == 1.5) {
		$dnds1_for_correction = 1.5;
	} elsif($dnds1 < 1.75) {
		print "\n### WARNING: dnds1=$dnds1 was specified with distance_actual=$distance_actual\.\n" .
			"### However, distance correction is only available for dnds1=\(0.1\|0.5\|1\|1.5\|2\)\.\n" .
			"### The correction for dnds1=1.5 will be used as an approximation\.\n";
		$dnds1_for_correction = 1.5;
	} elsif($dnds1 < 2) {
		print "\n### WARNING: dnds1=$dnds1 was specified with distance_actual=$distance_actual\.\n" .
			"### However, distance correction is only available for dnds1=\(0.1\|0.5\|1\|1.5\|2\)\.\n" .
			"### The correction for dnds1=2 will be used as an approximation\.\n";
		$dnds1_for_correction = 2;
	} elsif($dnds1 == 2) {
		$dnds1_for_correction = 2;
	} else {
		print "\n### WARNING: dnds1=$dnds1 was specified with distance_actual=$distance_actual\.\n" .
			"### However, distance correction is only available for dnds1=\(0.1\|0.5\|1\|1.5\|2\)\.\n" .
			"### The correction for dnds1=2 will be used as an approximation\.\n";
		$dnds1_for_correction = 2;
	}
	
	# Get closest dnds2
	if($dnds2 < 0.1) {
		print "\n### WARNING: dnds2=$dnds2 was specified with distance_actual=$distance_actual\.\n" .
			"### However, distance correction is only available for dnds2=\(0.1\|0.5\|1\|1.5\|2\)\.\n" .
			"### The correction for dnds2=0.1 will be used as an approximation\.\n";
		$dnds2_for_correction = 0.1;
	} elsif($dnds2 == 0.1) {
		$dnds2_for_correction = 0.1;
	} elsif($dnds2 < 0.3) {
		print "\n### WARNING: dnds2=$dnds2 was specified with distance_actual=$distance_actual\.\n" .
			"### However, distance correction is only available for dnds2=\(0.1\|0.5\|1\|1.5\|2\)\.\n" .
			"### The correction for dnds2=0.1 will be used as an approximation\.\n";
		$dnds2_for_correction = 0.1;
	} elsif($dnds2 < 0.5) {
		print "\n### WARNING: dnds2=$dnds2 was specified with distance_actual=$distance_actual\.\n" .
			"### However, distance correction is only available for dnds2=\(0.1\|0.5\|1\|1.5\|2\)\.\n" .
			"### The correction for dnds2=0.5 will be used as an approximation\.\n";
		$dnds2_for_correction = 0.5;
	} elsif($dnds2 == 0.5) {
		$dnds2_for_correction = 0.5;
	} elsif($dnds2 < 0.75) {
		print "\n### WARNING: dnds2=$dnds2 was specified with distance_actual=$distance_actual\.\n" .
			"### However, distance correction is only available for dnds2=\(0.1\|0.5\|1\|1.5\|2\)\.\n" .
			"### The correction for dnds2=0.5 will be used as an approximation\.\n";
		$dnds2_for_correction = 0.5;
	} elsif($dnds2 < 1) {
		print "\n### WARNING: dnds2=$dnds2 was specified with distance_actual=$distance_actual\.\n" .
			"### However, distance correction is only available for dnds2=\(0.1\|0.5\|1\|1.5\|2\)\.\n" .
			"### The correction for dnds2=1 will be used as an approximation\.\n";
		$dnds2_for_correction = 1;
	} elsif($dnds2 == 1) {
		$dnds2_for_correction = 1;
	} elsif($dnds2 < 1.25) {
		print "\n### WARNING: dnds2=$dnds2 was specified with distance_actual=$distance_actual\.\n" .
			"### However, distance correction is only available for dnds2=\(0.1\|0.5\|1\|1.5\|2\)\.\n" .
			"### The correction for dnds2=1 will be used as an approximation\.\n";
		$dnds2_for_correction = 1;
	} elsif($dnds2 < 1.5) {
		print "\n### WARNING: dnds2=$dnds2 was specified with distance_actual=$distance_actual\.\n" .
			"### However, distance correction is only available for dnds2=\(0.1\|0.5\|1\|1.5\|2\)\.\n" .
			"### The correction for dnds2=1.5 will be used as an approximation\.\n";
		$dnds2_for_correction = 1.5;
	} elsif($dnds2 == 1.5) {
		$dnds2_for_correction = 1.5;
	} elsif($dnds2 < 1.75) {
		print "\n### WARNING: dnds2=$dnds2 was specified with distance_actual=$distance_actual\.\n" .
			"### However, distance correction is only available for dnds2=\(0.1\|0.5\|1\|1.5\|2\)\.\n" .
			"### The correction for dnds2=1.5 will be used as an approximation\.\n";
		$dnds2_for_correction = 1.5;
	} elsif($dnds2 < 2) {
		print "\n### WARNING: dnds2=$dnds2 was specified with distance_actual=$distance_actual\.\n" .
			"### However, distance correction is only available for dnds2=\(0.1\|0.5\|1\|1.5\|2\)\.\n" .
			"### The correction for dnds2=2 will be used as an approximation\.\n";
		$dnds2_for_correction = 2;
	} elsif($dnds2 == 2) {
		$dnds2_for_correction = 2;
	} else {
		print "\n### WARNING: dnds2=$dnds2 was specified with distance_actual=$distance_actual\.\n" .
			"### However, distance correction is only available for dnds2=\(0.1\|0.5\|1\|1.5\|2\)\.\n" .
			"### The correction for dnds2=2 will be used as an approximation\.\n";
		$dnds2_for_correction = 2;
	}
	
	# Find best correction factor
	my $correction_factor = $distance_actual / $distances_expected{'ss13'}->{$R_for_correction}->{$dnds1_for_correction}->{$dnds2_for_correction};
	
	$distance = (0.05 * $correction_factor); # because the experiments for the observed values were performed with distance == 0.05
	
	print "\n################################################################################\n" .
		"### DISTANCE SUCESSFULLY CORRECTED!\n" .
		"### correction_factor=$correction_factor\n" .
		"### distance=$distance should result in ~d=$distance_actual\n" . 
		"################################################################################\n";
}


# Convert dnds values to array for future flexibility
my @dnds_1_values = ($dnds1);
my @dnds_2_values = ($dnds2);

if($num_seqs > 1024) {
	$num_seqs = 1024;
}

foreach my $dnds_1 (@dnds_1_values){
    foreach my $dnds_2 (@dnds_2_values){
        my $size_of_set = 1;
        my $number_branching_event = 10;
        $distance = $distance / $number_branching_event; # This is the distance between the sister pair.
        my $R=2; ###### This is transition transversion ratio =a/(2b)
        my $factor=0.9;  #### This will correct the simulation bias for neutral site

###--------------------------main part--------------------------------------
        #open(my $fh, '>', "SS13\_Phylogeny_distance_$distance\_dnds1_$dnds_1\_dnds2_$dnds_2\_R_$R.txt");
        open(my $fh, '>', "$outfile");

        if ($dnds_1>1) {
            if($dnds_2>1){
                $factor=$factor/$dnds_2/$dnds_1;
            } else{
                $factor=$factor/$dnds_1;
            }
        } elsif($dnds_2>1) {
            $factor=$factor/$dnds_2;
        }
        $distance=$distance/$factor;
        my @random_DNA=( );

		# Make a seed DNA set that is 1.5 times the needed length, since it will have STOP codons pruned,
		# then prune down to (3 * $num_codons) nucleotides
        @random_DNA = make_overlap_DNA_set((1.5 * $num_codons), $size_of_set); # start bigger because will be pruned
		PRUNE_DNA: for (my $i = 0; $i < @random_DNA; $i++) {
			my $this_random_DNA = $random_DNA[$i];
			
			if(length($this_random_DNA) > (3 * $num_codons)) {
				my $this_random_DNA_trimmed = substr($this_random_DNA, 0, (3 * $num_codons));
				$random_DNA[$i] = $this_random_DNA_trimmed;
			}
		}

        for (my $i=0; $i < $number_branching_event; $i++) {
            my @temp_set =( );
            foreach my $dna (@random_DNA){
                my @return_set = create_branching_set($dna,$distance,$R,$factor,$dnds_1,$dnds_2);
                push(@temp_set, @return_set);
            }
            @random_DNA = @temp_set;
            @temp_set = ();
        }

        my $loopVar = 1;
        
        # Shuffled list of sequence indexes; https://stackoverflow.com/questions/8963228/how-can-i-take-n-elements-at-random-from-a-perl-array
        my @shuffled_indexes = shuffle(0..$#random_DNA);
        my @pick_indexes = @shuffled_indexes[ 0 .. $num_seqs - 1 ];

        foreach my $dna (@random_DNA [ @pick_indexes ] ) {
            print $fh ">sim_seq_$loopVar\n";
            #remove the first bp and add a '-' to the last position to make the frame ss12.
            my $dna_i = mutation_and_selection($dna,$distance,$R,$factor,$dnds_1,$dnds_2);
            #substr($dna_i,1,1) = '';
            print $fh $dna_i;
            print $fh "\n";#add '-'
            $loopVar++;
        }

        close $fh;
    }
}
print "\nDone!\n\n";

exit;




sub make_overlap_DNA_set {
	my($num_codons,$size_of_set) = @_;     
	my $dna;
	my @set;
	for (my $i = 0; $i < $size_of_set ; ++$i) {
		$dna = make_overlap_DNA ( $num_codons ); 
		push(@set, $dna );
	}
	return @set;
}

sub create_branching_set {
    my($seq,$distance,$R,$factor,$dnds_1,$dnds_2) = @_;
    my $dna1;
    my $dna2;
    my @set;
    $dna1 = mutation_and_selection($seq,$distance,$R,$factor,$dnds_1,$dnds_2);
    push(@set, $dna1);
    $dna2 = mutation_and_selection($seq,$distance,$R,$factor,$dnds_1,$dnds_2);
    push(@set, $dna2);
    return @set;
}

#sub make_mutation_set{
#	my ($distance,$R,@random_DNA)=@_;
#	my @set;
#	foreach my $dna (@random_DNA){  
#		$dna=mutation_and_selection($dna,$distance,$R);  
#		push(@set,$dna);
#	}
#	return @set;
#}


sub mutation_and_selection{
	my ($dna,$distance,$R,$factor,$dnds_1,$dnds_2)=@_;
	my $select_base;
	my $newbase;
	my $oldbase;
	my $random_number;
	my $position;  
	my $mutate_time=int(length($dna)*$distance/2);
	for (my $i = 0; $i < $mutate_time; ++$i) {
		$position=randomposition($dna);
		$oldbase=substr($dna,$position,1);
		$random_number=rand(1);
		if ($oldbase eq 'A') {   ######
			if ($random_number < $R/(1+$R)) {
				$newbase='T';
			} elsif ($random_number < (2*$R+1)/(2+2*$R)){
				$newbase='C';
			} else{
				$newbase='G';
			}
		} elsif ($oldbase eq 'T') {
			if ($random_number < $R/(1+$R)) {
				$newbase='A';
			} elsif ($random_number < (2*$R+1)/(2+2*$R)){
				$newbase='C';
			} else{
				$newbase='G';
			}
		} elsif ($oldbase eq 'C') {
			if ($random_number < $R/(1+$R)) {
				$newbase='G';
			} elsif ($random_number < (2*$R+1)/(2+2*$R)){
				$newbase='A';
			} else{
				$newbase='T';
			}
		} elsif ($oldbase eq 'G') {
			if ($random_number < $R/(1+$R)) {
				$newbase='C';
			} elsif ($random_number < (2*$R+1)/(2+2*$R)){
				$newbase='A';
			} else{
				$newbase='T';
			}
		}
		$select_base=selection($dna,$newbase,$position,$factor,$dnds_1,$dnds_2);
		substr($dna,$position,1,$select_base); 
	}
	return $dna;
}


sub selection{
	my ($dna,$newbase,$position,$factor,$dnds_1,$dnds_2)=@_;
	my $codon_original_1;
	my $codon_mutation_1;
	my $codon_original_2;
	my $codon_mutation_2;
	my $stop_codon='_';
	my $original_aa1;
	my $original_aa2;
	my $mutate_aa1;
	my $mutate_aa2;
	my $select_base=substr($dna,$position,1);
	if ($position%3==2){
		$codon_original_1=substr($dna,$position-2,3);
		$codon_original_2=substr($dna,$position-1,3);
		$codon_mutation_1=$codon_original_1;  
		$codon_mutation_2=$codon_original_2;
		substr($codon_mutation_1,2,1,$newbase); 
		substr($codon_mutation_2,1,1,$newbase); 
	} elsif($position%3==0){
		$codon_original_1=substr($dna,$position,3);
		$codon_original_2=substr($dna,$position-2,3);
		$codon_mutation_1=$codon_original_1;  
		$codon_mutation_2=$codon_original_2;
		substr($codon_mutation_1,0,1,$newbase); 
		substr($codon_mutation_2,2,1,$newbase); 
	} elsif($position%3==1){
		$codon_original_1=substr($dna,$position-1,3);
		$codon_original_2=substr($dna,$position,3);
        $codon_mutation_1=$codon_original_1;  
		$codon_mutation_2=$codon_original_2;
		substr($codon_mutation_1,1,1,$newbase); 
		substr($codon_mutation_2,0,1,$newbase); 
	}
	$original_aa1=codon2aa($codon_original_1);
	$original_aa2=codon2aa($codon_original_2);
	$mutate_aa1=codon2aa($codon_mutation_1);
	$mutate_aa2=codon2aa($codon_mutation_2);	
	if (($original_aa1 eq $mutate_aa1)){
		if(($original_aa2 eq $mutate_aa2)){
			if (rand(1)<=$factor){
				$select_base=$newbase;
			}
		} elsif(( rand(1) <= $factor*$dnds_2) & ($mutate_aa2 ne $stop_codon)){
			$select_base=$newbase;
		} 
	} elsif(($mutate_aa1 ne $stop_codon)) { 
		if (($original_aa2 eq $mutate_aa2)){
			if( rand(1) <= $factor*$dnds_1){
				$select_base=$newbase;
			}
		} elsif((rand(1) <= $factor*$dnds_2*$dnds_1) & ($mutate_aa2 ne $stop_codon)){
			$select_base=$newbase;
		}
	}
	return $select_base;
}


sub randomposition{
	my($string)= @_;
	return (int(rand(length($string)-3))+1);  
}


sub make_overlap_DNA{ 
	my ($num_codons)= @_;
	my $dna;
	my $codon;
	my $stop_codon='_';
	for(my $i=0;$i<$num_codons;++$i){  
		$dna.=random_codon();
	}
    my $stop_counter = 1;
    while ($stop_counter > 0) {
        $stop_counter = 0;
        for(my $i=1;$i < (length($dna)-2);){
            $codon=substr($dna,$i,3);
            if (codon2aa($codon) eq $stop_codon) {
                substr($dna,$i,3)="";
                $stop_counter = $stop_counter + 1;
                $codon=substr($dna,$i-1,3);
                if (codon2aa($codon) eq $stop_codon) {
                    substr($dna,$i-1,3)="";
                    $stop_counter = $stop_counter + 1;
                }
            } else {
                $i=$i+3;
            }
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
	}else{
	print STDERR "Bad codon \"$codon\"!!\n";
	exit;
}
}


##########################################################################################
sub initialize_distances_expected_hash {
	my %distances_expected;
	$distances_expected{'ss13'}->{0.5}->{0.1}->{0.1} = 0.002553499079225211;
	$distances_expected{'ss13'}->{0.5}->{0.1}->{0.5} = 0.007405631575385161;
	$distances_expected{'ss13'}->{0.5}->{0.1}->{1} = 0.01346969435506098;
	$distances_expected{'ss13'}->{0.5}->{0.1}->{1.5} = 0.019160371136666018;
	$distances_expected{'ss13'}->{0.5}->{0.1}->{2} = 0.024716619739126547;
	$distances_expected{'ss13'}->{0.5}->{0.5}->{0.1} = 0.0075101989390322565;
	$distances_expected{'ss13'}->{0.5}->{0.5}->{0.5} = 0.01647325860411644;
	$distances_expected{'ss13'}->{0.5}->{0.5}->{1} = 0.02715800990797991;
	$distances_expected{'ss13'}->{0.5}->{0.5}->{1.5} = 0.037712558749170566;
	$distances_expected{'ss13'}->{0.5}->{0.5}->{2} = 0.047773626510093946;
	$distances_expected{'ss13'}->{0.5}->{1}->{0.1} = 0.01376102873310296;
	$distances_expected{'ss13'}->{0.5}->{1}->{0.5} = 0.02733002845879082;
	$distances_expected{'ss13'}->{0.5}->{1}->{1} = 0.04387581351453926;
	$distances_expected{'ss13'}->{0.5}->{1}->{1.5} = 0.059822803391936345;
	$distances_expected{'ss13'}->{0.5}->{1}->{2} = 0.07487910968750217;
	$distances_expected{'ss13'}->{0.5}->{1.5}->{0.1} = 0.01989693470978198;
	$distances_expected{'ss13'}->{0.5}->{1.5}->{0.5} = 0.03804479374004733;
	$distances_expected{'ss13'}->{0.5}->{1.5}->{1} = 0.059889566767553554;
	$distances_expected{'ss13'}->{0.5}->{1.5}->{1.5} = 0.08053053945707048;
	$distances_expected{'ss13'}->{0.5}->{1.5}->{2} = 0.10020835656564256;
	$distances_expected{'ss13'}->{0.5}->{2}->{0.1} = 0.025641521415335718;
	$distances_expected{'ss13'}->{0.5}->{2}->{0.5} = 0.04859309443965315;
	$distances_expected{'ss13'}->{0.5}->{2}->{1} = 0.07531465291437058;
	$distances_expected{'ss13'}->{0.5}->{2}->{1.5} = 0.10011761248588137;
	$distances_expected{'ss13'}->{0.5}->{2}->{2} = 0.12354994092148716;
	$distances_expected{'ss13'}->{1}->{0.1}->{0.1} = 0.0025455908983481835;
	$distances_expected{'ss13'}->{1}->{0.1}->{0.5} = 0.007416585704308192;
	$distances_expected{'ss13'}->{1}->{0.1}->{1} = 0.013364017449634742;
	$distances_expected{'ss13'}->{1}->{0.1}->{1.5} = 0.019261043974890635;
	$distances_expected{'ss13'}->{1}->{0.1}->{2} = 0.02475684992903486;
	$distances_expected{'ss13'}->{1}->{0.5}->{0.1} = 0.007595863807616301;
	$distances_expected{'ss13'}->{1}->{0.5}->{0.5} = 0.01642001754899392;
	$distances_expected{'ss13'}->{1}->{0.5}->{1} = 0.027172198054007755;
	$distances_expected{'ss13'}->{1}->{0.5}->{1.5} = 0.037601698790682;
	$distances_expected{'ss13'}->{1}->{0.5}->{2} = 0.04803170373109444;
	$distances_expected{'ss13'}->{1}->{1}->{0.1} = 0.013787893778724834;
	$distances_expected{'ss13'}->{1}->{1}->{0.5} = 0.027369154044244216;
	$distances_expected{'ss13'}->{1}->{1}->{1} = 0.043887614845288137;
	$distances_expected{'ss13'}->{1}->{1}->{1.5} = 0.05966953694214054;
	$distances_expected{'ss13'}->{1}->{1}->{2} = 0.07491151379725215;
	$distances_expected{'ss13'}->{1}->{1.5}->{0.1} = 0.019752878277867297;
	$distances_expected{'ss13'}->{1}->{1.5}->{0.5} = 0.03801075851861902;
	$distances_expected{'ss13'}->{1}->{1.5}->{1} = 0.05982325871662956;
	$distances_expected{'ss13'}->{1}->{1.5}->{1.5} = 0.08041472422072281;
	$distances_expected{'ss13'}->{1}->{1.5}->{2} = 0.09991302618514517;
	$distances_expected{'ss13'}->{1}->{2}->{0.1} = 0.025522130668651805;
	$distances_expected{'ss13'}->{1}->{2}->{0.5} = 0.04852224386746001;
	$distances_expected{'ss13'}->{1}->{2}->{1} = 0.07517949458326328;
	$distances_expected{'ss13'}->{1}->{2}->{1.5} = 0.10023852358560209;
	$distances_expected{'ss13'}->{1}->{2}->{2} = 0.1235407079122379;
	$distances_expected{'ss13'}->{2}->{0.1}->{0.1} = 0.0025303896544551585;
	$distances_expected{'ss13'}->{2}->{0.1}->{0.5} = 0.007373458575224243;
	$distances_expected{'ss13'}->{2}->{0.1}->{1} = 0.013376065510943258;
	$distances_expected{'ss13'}->{2}->{0.1}->{1.5} = 0.019166535493575134;
	$distances_expected{'ss13'}->{2}->{0.1}->{2} = 0.024729294206261323;
	$distances_expected{'ss13'}->{2}->{0.5}->{0.1} = 0.007528067654960151;
	$distances_expected{'ss13'}->{2}->{0.5}->{0.5} = 0.016402909227826733;
	$distances_expected{'ss13'}->{2}->{0.5}->{1} = 0.027292625751095076;
	$distances_expected{'ss13'}->{2}->{0.5}->{1.5} = 0.03772251048824549;
	$distances_expected{'ss13'}->{2}->{0.5}->{2} = 0.047779411037381564;
	$distances_expected{'ss13'}->{2}->{1}->{0.1} = 0.013682572340261032;
	$distances_expected{'ss13'}->{2}->{1}->{0.5} = 0.02747193351009367;
	$distances_expected{'ss13'}->{2}->{1}->{1} = 0.04393225482627202;
	$distances_expected{'ss13'}->{2}->{1}->{1.5} = 0.059704258725456094;
	$distances_expected{'ss13'}->{2}->{1}->{2} = 0.07495338085279715;
	$distances_expected{'ss13'}->{2}->{1.5}->{0.1} = 0.01972606452986857;
	$distances_expected{'ss13'}->{2}->{1.5}->{0.5} = 0.0381609573802136;
	$distances_expected{'ss13'}->{2}->{1.5}->{1} = 0.05980270138440549;
	$distances_expected{'ss13'}->{2}->{1.5}->{1.5} = 0.08046675627209691;
	$distances_expected{'ss13'}->{2}->{1.5}->{2} = 0.09998959314075238;
	$distances_expected{'ss13'}->{2}->{2}->{0.1} = 0.02552977535284197;
	$distances_expected{'ss13'}->{2}->{2}->{0.5} = 0.04833830002868566;
	$distances_expected{'ss13'}->{2}->{2}->{1} = 0.07532444422306045;
	$distances_expected{'ss13'}->{2}->{2}->{1.5} = 0.10018096747049277;
	$distances_expected{'ss13'}->{2}->{2}->{2} = 0.1233748814173843;
	$distances_expected{'sas11'}->{0.5}->{0.1}->{0.1} = 0.0022768333110837866;
	$distances_expected{'sas11'}->{0.5}->{0.1}->{0.5} = 0.007200018775088351;
	$distances_expected{'sas11'}->{0.5}->{0.1}->{1} = 0.013266948825810613;
	$distances_expected{'sas11'}->{0.5}->{0.1}->{1.5} = 0.019249585756927784;
	$distances_expected{'sas11'}->{0.5}->{0.1}->{2} = 0.024990125698571355;
	$distances_expected{'sas11'}->{0.5}->{0.5}->{0.1} = 0.007166750669700983;
	$distances_expected{'sas11'}->{0.5}->{0.5}->{0.5} = 0.01616542386673126;
	$distances_expected{'sas11'}->{0.5}->{0.5}->{1} = 0.026855260511144546;
	$distances_expected{'sas11'}->{0.5}->{0.5}->{1.5} = 0.037532679544312594;
	$distances_expected{'sas11'}->{0.5}->{0.5}->{2} = 0.04791537250995916;
	$distances_expected{'sas11'}->{0.5}->{1}->{0.1} = 0.013394105534084065;
	$distances_expected{'sas11'}->{0.5}->{1}->{0.5} = 0.02699854301032741;
	$distances_expected{'sas11'}->{0.5}->{1}->{1} = 0.04330411400378343;
	$distances_expected{'sas11'}->{0.5}->{1}->{1.5} = 0.05913812373851546;
	$distances_expected{'sas11'}->{0.5}->{1}->{2} = 0.07426106913384414;
	$distances_expected{'sas11'}->{0.5}->{1.5}->{0.1} = 0.019154169453898647;
	$distances_expected{'sas11'}->{0.5}->{1.5}->{0.5} = 0.03759196240763361;
	$distances_expected{'sas11'}->{0.5}->{1.5}->{1} = 0.05915742645991949;
	$distances_expected{'sas11'}->{0.5}->{1.5}->{1.5} = 0.08000744853963132;
	$distances_expected{'sas11'}->{0.5}->{1.5}->{2} = 0.09932685595967136;
	$distances_expected{'sas11'}->{0.5}->{2}->{0.1} = 0.024955661657272296;
	$distances_expected{'sas11'}->{0.5}->{2}->{0.5} = 0.04777460975152221;
	$distances_expected{'sas11'}->{0.5}->{2}->{1} = 0.07443563279969564;
	$distances_expected{'sas11'}->{0.5}->{2}->{1.5} = 0.09915547787322208;
	$distances_expected{'sas11'}->{0.5}->{2}->{2} = 0.12237235429021477;
	$distances_expected{'sas11'}->{1}->{0.1}->{0.1} = 0.002278217448258623;
	$distances_expected{'sas11'}->{1}->{0.1}->{0.5} = 0.007223309126902647;
	$distances_expected{'sas11'}->{1}->{0.1}->{1} = 0.013371809908552799;
	$distances_expected{'sas11'}->{1}->{0.1}->{1.5} = 0.019103910023989;
	$distances_expected{'sas11'}->{1}->{0.1}->{2} = 0.02489940429536227;
	$distances_expected{'sas11'}->{1}->{0.5}->{0.1} = 0.007157492500208443;
	$distances_expected{'sas11'}->{1}->{0.5}->{0.5} = 0.016129042421645562;
	$distances_expected{'sas11'}->{1}->{0.5}->{1} = 0.02693708768634071;
	$distances_expected{'sas11'}->{1}->{0.5}->{1.5} = 0.037435129773893114;
	$distances_expected{'sas11'}->{1}->{0.5}->{2} = 0.04780261604310584;
	$distances_expected{'sas11'}->{1}->{1}->{0.1} = 0.013200216337639747;
	$distances_expected{'sas11'}->{1}->{1}->{0.5} = 0.026895530235568055;
	$distances_expected{'sas11'}->{1}->{1}->{1} = 0.04330142587343844;
	$distances_expected{'sas11'}->{1}->{1}->{1.5} = 0.05921972447650914;
	$distances_expected{'sas11'}->{1}->{1}->{2} = 0.074395496977429;
	$distances_expected{'sas11'}->{1}->{1.5}->{0.1} = 0.01906871450397401;
	$distances_expected{'sas11'}->{1}->{1.5}->{0.5} = 0.037471175591695646;
	$distances_expected{'sas11'}->{1}->{1.5}->{1} = 0.05912109762107073;
	$distances_expected{'sas11'}->{1}->{1.5}->{1.5} = 0.079802387162011;
	$distances_expected{'sas11'}->{1}->{1.5}->{2} = 0.09925521275627938;
	$distances_expected{'sas11'}->{1}->{2}->{0.1} = 0.02496382998158287;
	$distances_expected{'sas11'}->{1}->{2}->{0.5} = 0.04761467448898306;
	$distances_expected{'sas11'}->{1}->{2}->{1} = 0.0742571345678908;
	$distances_expected{'sas11'}->{1}->{2}->{1.5} = 0.09923450919811587;
	$distances_expected{'sas11'}->{1}->{2}->{2} = 0.12247406202171603;
	$distances_expected{'sas11'}->{2}->{0.1}->{0.1} = 0.0022609684511831307;
	$distances_expected{'sas11'}->{2}->{0.1}->{0.5} = 0.007263477699193108;
	$distances_expected{'sas11'}->{2}->{0.1}->{1} = 0.0132539556580375;
	$distances_expected{'sas11'}->{2}->{0.1}->{1.5} = 0.019192900130862275;
	$distances_expected{'sas11'}->{2}->{0.1}->{2} = 0.024793195351219686;
	$distances_expected{'sas11'}->{2}->{0.5}->{0.1} = 0.007283223214154738;
	$distances_expected{'sas11'}->{2}->{0.5}->{0.5} = 0.016063338173481878;
	$distances_expected{'sas11'}->{2}->{0.5}->{1} = 0.026988129896548366;
	$distances_expected{'sas11'}->{2}->{0.5}->{1.5} = 0.037561285506267764;
	$distances_expected{'sas11'}->{2}->{0.5}->{2} = 0.04764525694605209;
	$distances_expected{'sas11'}->{2}->{1}->{0.1} = 0.013310599471499845;
	$distances_expected{'sas11'}->{2}->{1}->{0.5} = 0.027021579164425796;
	$distances_expected{'sas11'}->{2}->{1}->{1} = 0.04336415767828425;
	$distances_expected{'sas11'}->{2}->{1}->{1.5} = 0.05912716879831455;
	$distances_expected{'sas11'}->{2}->{1}->{2} = 0.07439331376767837;
	$distances_expected{'sas11'}->{2}->{1.5}->{0.1} = 0.019227411858894174;
	$distances_expected{'sas11'}->{2}->{1.5}->{0.5} = 0.03740546684406309;
	$distances_expected{'sas11'}->{2}->{1.5}->{1} = 0.05915837062922837;
	$distances_expected{'sas11'}->{2}->{1.5}->{1.5} = 0.07974174994583488;
	$distances_expected{'sas11'}->{2}->{1.5}->{2} = 0.09933226123461666;
	$distances_expected{'sas11'}->{2}->{2}->{0.1} = 0.02496546444705775;
	$distances_expected{'sas11'}->{2}->{2}->{0.5} = 0.04758834866946604;
	$distances_expected{'sas11'}->{2}->{2}->{1} = 0.07439733778205404;
	$distances_expected{'sas11'}->{2}->{2}->{1.5} = 0.09917220740757293;
	$distances_expected{'sas11'}->{2}->{2}->{2} = 0.12235608626008447;
	$distances_expected{'sas12'}->{0.5}->{0.1}->{0.1} = 0.00794345587198668;
	$distances_expected{'sas12'}->{0.5}->{0.1}->{0.5} = 0.010446468246119259;
	$distances_expected{'sas12'}->{0.5}->{0.1}->{1} = 0.013799135276824921;
	$distances_expected{'sas12'}->{0.5}->{0.1}->{1.5} = 0.017145999044004292;
	$distances_expected{'sas12'}->{0.5}->{0.1}->{2} = 0.020402341913418217;
	$distances_expected{'sas12'}->{0.5}->{0.5}->{0.1} = 0.010568130698290042;
	$distances_expected{'sas12'}->{0.5}->{0.5}->{0.5} = 0.01840971495864523;
	$distances_expected{'sas12'}->{0.5}->{0.5}->{1} = 0.027774084795191195;
	$distances_expected{'sas12'}->{0.5}->{0.5}->{1.5} = 0.03702471854840683;
	$distances_expected{'sas12'}->{0.5}->{0.5}->{2} = 0.046022523295885076;
	$distances_expected{'sas12'}->{0.5}->{1}->{0.1} = 0.013901978350465707;
	$distances_expected{'sas12'}->{0.5}->{1}->{0.5} = 0.02773668403381717;
	$distances_expected{'sas12'}->{0.5}->{1}->{1} = 0.04451406522401335;
	$distances_expected{'sas12'}->{0.5}->{1}->{1.5} = 0.06052137697492862;
	$distances_expected{'sas12'}->{0.5}->{1}->{2} = 0.07610146499273047;
	$distances_expected{'sas12'}->{0.5}->{1.5}->{0.1} = 0.017175104404595323;
	$distances_expected{'sas12'}->{0.5}->{1.5}->{0.5} = 0.03704823784393142;
	$distances_expected{'sas12'}->{0.5}->{1.5}->{1} = 0.06058917907759461;
	$distances_expected{'sas12'}->{0.5}->{1.5}->{1.5} = 0.08286974154514895;
	$distances_expected{'sas12'}->{0.5}->{1.5}->{2} = 0.10401048419332326;
	$distances_expected{'sas12'}->{0.5}->{2}->{0.1} = 0.020330939024488123;
	$distances_expected{'sas12'}->{0.5}->{2}->{0.5} = 0.04609708805539978;
	$distances_expected{'sas12'}->{0.5}->{2}->{1} = 0.07626108639901424;
	$distances_expected{'sas12'}->{0.5}->{2}->{1.5} = 0.10396079703103692;
	$distances_expected{'sas12'}->{0.5}->{2}->{2} = 0.12974397314368125;
	$distances_expected{'sas12'}->{1}->{0.1}->{0.1} = 0.007882942561485537;
	$distances_expected{'sas12'}->{1}->{0.1}->{0.5} = 0.010573233685274372;
	$distances_expected{'sas12'}->{1}->{0.1}->{1} = 0.013891782050940502;
	$distances_expected{'sas12'}->{1}->{0.1}->{1.5} = 0.017190357190931052;
	$distances_expected{'sas12'}->{1}->{0.1}->{2} = 0.020257512486257352;
	$distances_expected{'sas12'}->{1}->{0.5}->{0.1} = 0.010628809296458467;
	$distances_expected{'sas12'}->{1}->{0.5}->{0.5} = 0.018398994797876024;
	$distances_expected{'sas12'}->{1}->{0.5}->{1} = 0.027762805755464463;
	$distances_expected{'sas12'}->{1}->{0.5}->{1.5} = 0.037136680027448105;
	$distances_expected{'sas12'}->{1}->{0.5}->{2} = 0.04602015018697642;
	$distances_expected{'sas12'}->{1}->{1}->{0.1} = 0.013940928870910058;
	$distances_expected{'sas12'}->{1}->{1}->{0.5} = 0.027770860066939084;
	$distances_expected{'sas12'}->{1}->{1}->{1} = 0.04446652768398662;
	$distances_expected{'sas12'}->{1}->{1}->{1.5} = 0.06067944859341203;
	$distances_expected{'sas12'}->{1}->{1}->{2} = 0.07615700323566706;
	$distances_expected{'sas12'}->{1}->{1.5}->{0.1} = 0.017165455737956567;
	$distances_expected{'sas12'}->{1}->{1.5}->{0.5} = 0.037071357925934065;
	$distances_expected{'sas12'}->{1}->{1.5}->{1} = 0.0606984632531126;
	$distances_expected{'sas12'}->{1}->{1.5}->{1.5} = 0.08293050626571176;
	$distances_expected{'sas12'}->{1}->{1.5}->{2} = 0.10397020059579742;
	$distances_expected{'sas12'}->{1}->{2}->{0.1} = 0.020376204373598855;
	$distances_expected{'sas12'}->{1}->{2}->{0.5} = 0.045983781197256934;
	$distances_expected{'sas12'}->{1}->{2}->{1} = 0.07607477398274012;
	$distances_expected{'sas12'}->{1}->{2}->{1.5} = 0.1041295838711888;
	$distances_expected{'sas12'}->{1}->{2}->{2} = 0.1296809372511477;
	$distances_expected{'sas12'}->{2}->{0.1}->{0.1} = 0.007942573983241248;
	$distances_expected{'sas12'}->{2}->{0.1}->{0.5} = 0.010577777855219756;
	$distances_expected{'sas12'}->{2}->{0.1}->{1} = 0.013799947382760481;
	$distances_expected{'sas12'}->{2}->{0.1}->{1.5} = 0.017133226472735928;
	$distances_expected{'sas12'}->{2}->{0.1}->{2} = 0.020240910427509032;
	$distances_expected{'sas12'}->{2}->{0.5}->{0.1} = 0.010612319970049418;
	$distances_expected{'sas12'}->{2}->{0.5}->{0.5} = 0.01831969337853124;
	$distances_expected{'sas12'}->{2}->{0.5}->{1} = 0.027807271572871287;
	$distances_expected{'sas12'}->{2}->{0.5}->{1.5} = 0.03716576299005542;
	$distances_expected{'sas12'}->{2}->{0.5}->{2} = 0.04609893341282696;
	$distances_expected{'sas12'}->{2}->{1}->{0.1} = 0.013924712998245279;
	$distances_expected{'sas12'}->{2}->{1}->{0.5} = 0.02774410618906607;
	$distances_expected{'sas12'}->{2}->{1}->{1} = 0.044438423277719435;
	$distances_expected{'sas12'}->{2}->{1}->{1.5} = 0.06065147265826003;
	$distances_expected{'sas12'}->{2}->{1}->{2} = 0.07601458847073482;
	$distances_expected{'sas12'}->{2}->{1.5}->{0.1} = 0.017180825844686477;
	$distances_expected{'sas12'}->{2}->{1.5}->{0.5} = 0.036963712564722764;
	$distances_expected{'sas12'}->{2}->{1.5}->{1} = 0.060636072338746135;
	$distances_expected{'sas12'}->{2}->{1.5}->{1.5} = 0.08291831996105777;
	$distances_expected{'sas12'}->{2}->{1.5}->{2} = 0.10381061101679581;
	$distances_expected{'sas12'}->{2}->{2}->{0.1} = 0.0204190620183191;
	$distances_expected{'sas12'}->{2}->{2}->{0.5} = 0.04625370871461132;
	$distances_expected{'sas12'}->{2}->{2}->{1} = 0.07602363488178461;
	$distances_expected{'sas12'}->{2}->{2}->{1.5} = 0.1041065713562011;
	$distances_expected{'sas12'}->{2}->{2}->{2} = 0.12953219743565506;
	$distances_expected{'sas13'}->{0.5}->{0.1}->{0.1} = 0.0027436785473042465;
	$distances_expected{'sas13'}->{0.5}->{0.1}->{0.5} = 0.007485086131290487;
	$distances_expected{'sas13'}->{0.5}->{0.1}->{1} = 0.013235739345163255;
	$distances_expected{'sas13'}->{0.5}->{0.1}->{1.5} = 0.018898904110118287;
	$distances_expected{'sas13'}->{0.5}->{0.1}->{2} = 0.024374556361534898;
	$distances_expected{'sas13'}->{0.5}->{0.5}->{0.1} = 0.007411773542217746;
	$distances_expected{'sas13'}->{0.5}->{0.5}->{0.5} = 0.016503135473164815;
	$distances_expected{'sas13'}->{0.5}->{0.5}->{1} = 0.027366743361412975;
	$distances_expected{'sas13'}->{0.5}->{0.5}->{1.5} = 0.038007943725895534;
	$distances_expected{'sas13'}->{0.5}->{0.5}->{2} = 0.04837603283483687;
	$distances_expected{'sas13'}->{0.5}->{1}->{0.1} = 0.013256397227671638;
	$distances_expected{'sas13'}->{0.5}->{1}->{0.5} = 0.027424177872201902;
	$distances_expected{'sas13'}->{0.5}->{1}->{1} = 0.04439315447451625;
	$distances_expected{'sas13'}->{0.5}->{1}->{1.5} = 0.06088038188778819;
	$distances_expected{'sas13'}->{0.5}->{1}->{2} = 0.07670346157895577;
	$distances_expected{'sas13'}->{0.5}->{1.5}->{0.1} = 0.01888196936307576;
	$distances_expected{'sas13'}->{0.5}->{1.5}->{0.5} = 0.0380200145025526;
	$distances_expected{'sas13'}->{0.5}->{1.5}->{1} = 0.060787094928166675;
	$distances_expected{'sas13'}->{0.5}->{1.5}->{1.5} = 0.08234356934078853;
	$distances_expected{'sas13'}->{0.5}->{1.5}->{2} = 0.10272055483239634;
	$distances_expected{'sas13'}->{0.5}->{2}->{0.1} = 0.0244193116681961;
	$distances_expected{'sas13'}->{0.5}->{2}->{0.5} = 0.048477872304676504;
	$distances_expected{'sas13'}->{0.5}->{2}->{1} = 0.07659428985230689;
	$distances_expected{'sas13'}->{0.5}->{2}->{1.5} = 0.10305535287682163;
	$distances_expected{'sas13'}->{0.5}->{2}->{2} = 0.127291174989403;
	$distances_expected{'sas13'}->{1}->{0.1}->{0.1} = 0.002696472143844964;
	$distances_expected{'sas13'}->{1}->{0.1}->{0.5} = 0.007471493414065086;
	$distances_expected{'sas13'}->{1}->{0.1}->{1} = 0.01333416462821556;
	$distances_expected{'sas13'}->{1}->{0.1}->{1.5} = 0.018853108923381233;
	$distances_expected{'sas13'}->{1}->{0.1}->{2} = 0.024358544565655497;
	$distances_expected{'sas13'}->{1}->{0.5}->{0.1} = 0.007475024440666795;
	$distances_expected{'sas13'}->{1}->{0.5}->{0.5} = 0.016403348117518257;
	$distances_expected{'sas13'}->{1}->{0.5}->{1} = 0.02737481744508018;
	$distances_expected{'sas13'}->{1}->{0.5}->{1.5} = 0.0379662070855258;
	$distances_expected{'sas13'}->{1}->{0.5}->{2} = 0.048365634345922455;
	$distances_expected{'sas13'}->{1}->{1}->{0.1} = 0.013204356364487578;
	$distances_expected{'sas13'}->{1}->{1}->{0.5} = 0.02744178185000446;
	$distances_expected{'sas13'}->{1}->{1}->{1} = 0.044358944497285537;
	$distances_expected{'sas13'}->{1}->{1}->{1.5} = 0.0608208988012848;
	$distances_expected{'sas13'}->{1}->{1}->{2} = 0.0766198612120867;
	$distances_expected{'sas13'}->{1}->{1.5}->{0.1} = 0.018934356354425082;
	$distances_expected{'sas13'}->{1}->{1.5}->{0.5} = 0.03808426829888013;
	$distances_expected{'sas13'}->{1}->{1.5}->{1} = 0.06079424510604478;
	$distances_expected{'sas13'}->{1}->{1.5}->{1.5} = 0.08246680605706018;
	$distances_expected{'sas13'}->{1}->{1.5}->{2} = 0.102811461324789;
	$distances_expected{'sas13'}->{1}->{2}->{0.1} = 0.02430187042691854;
	$distances_expected{'sas13'}->{1}->{2}->{0.5} = 0.04850267476441502;
	$distances_expected{'sas13'}->{1}->{2}->{1} = 0.07656661658506057;
	$distances_expected{'sas13'}->{1}->{2}->{1.5} = 0.10285642614632211;
	$distances_expected{'sas13'}->{1}->{2}->{2} = 0.12731204119031883;
	$distances_expected{'sas13'}->{2}->{0.1}->{0.1} = 0.0026864575608146254;
	$distances_expected{'sas13'}->{2}->{0.1}->{0.5} = 0.007479133773118214;
	$distances_expected{'sas13'}->{2}->{0.1}->{1} = 0.013333217811411443;
	$distances_expected{'sas13'}->{2}->{0.1}->{1.5} = 0.018964236117923347;
	$distances_expected{'sas13'}->{2}->{0.1}->{2} = 0.024198530457402497;
	$distances_expected{'sas13'}->{2}->{0.5}->{0.1} = 0.007369439452572063;
	$distances_expected{'sas13'}->{2}->{0.5}->{0.5} = 0.0165382632338518;
	$distances_expected{'sas13'}->{2}->{0.5}->{1} = 0.027361110179141942;
	$distances_expected{'sas13'}->{2}->{0.5}->{1.5} = 0.03804925095655916;
	$distances_expected{'sas13'}->{2}->{0.5}->{2} = 0.04835789217286259;
	$distances_expected{'sas13'}->{2}->{1}->{0.1} = 0.013169002637547039;
	$distances_expected{'sas13'}->{2}->{1}->{0.5} = 0.027365305205155104;
	$distances_expected{'sas13'}->{2}->{1}->{1} = 0.044251874143513024;
	$distances_expected{'sas13'}->{2}->{1}->{1.5} = 0.06072970471200741;
	$distances_expected{'sas13'}->{2}->{1}->{2} = 0.07658162526031544;
	$distances_expected{'sas13'}->{2}->{1.5}->{0.1} = 0.018846627963534368;
	$distances_expected{'sas13'}->{2}->{1.5}->{0.5} = 0.03811529178663164;
	$distances_expected{'sas13'}->{2}->{1.5}->{1} = 0.06079980865792167;
	$distances_expected{'sas13'}->{2}->{1.5}->{1.5} = 0.08239519809208282;
	$distances_expected{'sas13'}->{2}->{1.5}->{2} = 0.10277230002892647;
	$distances_expected{'sas13'}->{2}->{2}->{0.1} = 0.02429807822110538;
	$distances_expected{'sas13'}->{2}->{2}->{0.5} = 0.04836026154734151;
	$distances_expected{'sas13'}->{2}->{2}->{1} = 0.07662744570644713;
	$distances_expected{'sas13'}->{2}->{2}->{1.5} = 0.10291248684898577;
	$distances_expected{'sas13'}->{2}->{2}->{2} = 0.12730618929694906;
	
	return \%distances_expected;
}

