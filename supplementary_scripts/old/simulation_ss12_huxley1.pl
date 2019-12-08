#! /nas3/cnelson/bin/anaconda2/bin/perl
# /usr/bin/perl

# overlapping gene selection simulation
# second base in ORF1 is first base in ORF2
#select codon, and use Jukes Canter Model to correct, add transition transversion ratio
use strict;
use warnings;
use diagnostics;
use autodie;
#srand (time|$$);  # introduce a random seed each time
my $random_seed = srand(time ^ $$ ^ unpack "%32L*", `ps wwaxl | gzip`); # (Programming Perl, p. 955)

# User-defined input
my @dnds_1_values = ($ARGV[0]);
my @dnds_2_values = ($ARGV[1]);
my $number_of_codons = $ARGV[2]; # alignment LENGTH
my $sample_size = $ARGV[3];
my $distance = $ARGV[4];
my $outfile = $ARGV[5];

print("PARAMETERS\n");
print("random_seed=$random_seed\n");
print("dnds_1_values=@dnds_1_values\n");
print("dnds_2_values=@dnds_2_values\n");
print("number_of_codons=$number_of_codons\n");
print "random_seed=$random_seed\n";

if($number_of_sequences > 1024) {
	$number_of_sequences = 1024;
}

######----------------------input--------------------------
#my $dnds_1=1;#shift;
#my @dnds_1_values = (0.1,0.5,1,1.5,2);
#my @dnds_2_values = (0.1,0.5,1,1.5,2);

foreach  my $dnds_1 (@dnds_1_values){
    foreach my $dnds_2 (@dnds_2_values){
        my $size_of_set=1;
        my $number_branching_event = 10;
        #my $sample_size = 2^10;#This is the sample size before, not replaced by branching events
        #my $number_of_codons=100000;
        #my $distance=0.05; #This now is the longest distance between two sequences
        $distance = $distance/$number_branching_event; # This is the distance between the sister pair.
        my $R=0.5; ######     This is transition transversion ratio =a/(2b)
        my $factor=0.9;  ####This will correct the simulation bias for neutral site

###--------------------------main part--------------------------------------
        #open(my $fh, '>', "SS12\_Phylogeny_distance_$distance\_dnds1_$dnds_1\_dnds2_$dnds_2\_R_$R.txt");
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
        
        
        @random_DNA=make_overlap_DNA_set($number_of_codons,$size_of_set);
        
        
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
        foreach my $dna (@random_DNA){
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
print "Done!";

exit;




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

sub create_branching_set {
    my($seq,$distance,$R) = @_;
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
	my ($number_of_codons)= @_;
	my $dna;
	my $codon;
	my $stop_codon='_';
	for(my $i=0;$i<$number_of_codons;++$i){  
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



