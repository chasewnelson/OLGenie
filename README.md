# OLGenie
Perl software for analysis of selection in overlapping genes (OLGs) using sister pairs to estimate dN/dS with the Wei-Zhang method.

## Background
Given the codon triplet nature of the genetic code and the two antiparallel strands of the DNA molecule, a single segment of DNA has the potential to encode 6 reading frames: three in the forward (sense) direction and three in the reverse (antisense) direction. Although all possible frames are seldom if ever used, a substantial fraction of genes in taxa ranging from viruses to humans encode overlapping gene (OLG) pairs, running in either the same (ss; sense-sense) or opposite directions (sas; sense-antisense) (Sabath 2009). We refer to these overlapping phases as ss12, ss13, sas11, sas12, or sas13, where the first number refers to the codon position in the reference gene (ORF1), and the second number refers to the codon position in the alternative (overlapping) gene (ORF2):

<img src="https://github.com/chasewnelson/OLGenie/blob/master/OLGs_all_phases.png" alt="sas12 logo" align="middle" width="400px">

For example, in sas12, genes overlap in sense-antisense relationship such that position 1 of codons in the sense reference gene correspond to position 2 of codons in the reverse strand overlapping gene. In other words, the sense gene's first codon position overlaps the antisense gene's second codon position:
	
<img src="https://github.com/chasewnelson/OLGenie/blob/master/sas12_figure.png" alt="sas12 logo" align="middle" width="2500px">

Unfortunately, standard methods for detecting natural selection using dN/dS do not apply to OLGs, because a mutation that is synonymous in one frame may be nonsynonymous in another, and *vice versa*. Although some methods for detecting natural selection in OLGs have been developed, they are generally computationally intensive and limited in utility (*e.g.*, Wei and Zhang 2015; Sabath et al. 2008). Thus, it is necessary to develop improved approaches for detecting selection in OLGs. 

OLGenie represents a combination of the methods of Hughes et al. (2006), Nelson et al. (2015), and Wei and Zhang (2015), tailored for detecting selection in OLGs. The Wei-Zhang method is a counting method analogous to the modified Nei-Gojobori (1986) method, but for use with OLGs. It considers the effects of mutations in the overlapping frame to determine the numerator and denominator of dN and dS.  For example, dN is usually calculated as the mean number of nonsynonymous (amino acid changing) nucleotide differences per nonsynonymous nucleotide site, and dS is similarly calculated for synonymous (silent) differences and sites. In order to control for the possibility that synonymous sites in the frame of interest may be under selection in the alternative overlapping reading frame, Wei-Zhang instead considers the expanded measures dNN, dNS, dSN, and dSS, where the first subscript refers to the reference (sense) gene (ORF1), and the second to the alternative (sense or antisense) OLG (ORF2). For example, dNN refers to the normalized number of differences which are nonsynonymous in both the reference and the alternative frames (i.e., NN). Using these measures, it is possible to calculate dN/dS for the reference gene as dNN/dSN or dNS/dSS, and to calculate the same for the alternative overlapping gene as dNN/dNS or dSN/dSS, i.e., the subscript in the alternative OLG is held constant to control for OLG effects. 

## Method
Written in Perl with no dependencies, OLGenie performs the following tasks:
			
* Examines a user-provided FASTA alignment of the overlapping region of ONE CODING REGION from the reference gene (i.e., mother gene; larger gene; ORF1) point of view.
* Identifies all well-supported sister-pairs from a user-provided phylogenetic tree. Making the reasonable assumption that the members of each pair are more closely related to one another than to any other sequences in the phylogeny (i.e., the assumption that they are true sisters), it follows that every difference between the two members of a pair has arisen since their last common ancestor, and are thus (statistically) independent of differences between all other pair members.
* Calculates the number of NN, SN, NS, and SS sites and differences for every sister pair.

These tasks are sufficiently fast that no parallelism is necessary at the level of the single gene alignment. For datasets with many genes, the user can implement their own parallelization by running numerous alignments simultaneously.

Several tests can subsequently be used to detect natural selection, e.g., purifying selection for ORF1 (pNN < pSN; pNS < pSS) and ORF2 (pNN < pNS; pSN < pSS). These are implemented in the accompanying R script OLGenie_selection_tests.R:

* For datasets with very few differences between pair members: 
	* Fisher's Exact Test (e.g., [NN\_sitesWithDiffs/SN\_sitesWithDiffs] X [NN\_remainingSites/SN\_remainingSites]. See Zhang et al. 1997.
	* Binomial Test (e.g., for ORF1, expected proportion of ORF1 diffs NN under neutrality is [NN\_sites]/[NN\_sites+SN\_sites])
* For datasets with sufficient differences between sister pairs, a Paired Wilcoxon Test (e.g., that pNN<pSN). See Hughes et al. 2006.

## How to Use
**OLGenie** accepts as input three necessary parameters: 

* **--fasta\_file**: a FASTA file containing multiple aligned sequences of one coding sequence. The entire coding sequence must be an overlapping gene (OLG), with no non-overlapping codons. The frame must be the frame of the reference gene (ORF1) (see the --phase option). It is recommended that the user translate the gene sequences, align at the amino acid level, and then impose the amino acid alignment on the DNA alignment to preserve complete codons. (See align\_codon2aa.pl at <a target="_blank" href="https://github.com/chasewnelson/CHASeq">CHASeq</a>.)
* **--tree\_file**: a text file containing one Newick tree using the exact sequence names (headers) as the FASTA. If multiple trees are present, only the first will be used.
* **--phase**: the phase of the overlapping gene (OLG) relationship: ss12, ss13, sas11, sas12, or sas13 (see above).
* **--prune\_polytomies**: flag that indicates that polytomies should be reduced to the two sequences with most data (fewest gaps).
* **--min\_support**: minimum bootstrap support (0-100) required for a sister pair to be included in the analysis.

Call **OLGenie** as follows:

    OLGenie.pl --fasta_file=<alignment>.fasta --tree_file=<phylogeny>.treefile --min_support=50 --phase=sas12

## References
* Hughes AL, Friedman R, Glenn NL. 2006. The Future of Data Analysis in Evolutionary Genomics. Current Genomics 7:227–234.
* Nei M, Gojobori T. 1986. Simple methods for estimating the numbers of synonymous and nonsynonymous nucleotide substitutions. Mol. Biol. Evol. 3:418–426.
* Nelson CW, Moncla LH, Hughes AL. 2015. SNPGenie: estimating evolutionary parameters to detect natural selection using pooled next-generation sequencing data. Bioinformatics 31:3709–3711.
* Sabath N. 2009. Molecular Evolution of Overlapping Genes. ProQuest Dissertation #3405062.* Sabath N, Landan G, Graur D. 2008. A method for the simultaneous estimation of selection intensities in overlapping genes. PLoS One 3(12): e3996.* Wei X, Zhang J. 2015. A simple method for estimating the strength of natural selection on overlapping genes. Genome Biol. Evol. 7:381–390.
* Zhang J, Kumar S, Nei M. 1997. Small-sample tests of episodic adaptive evolution: a case study of primate lysozymes. Molecular Biology and Evolution 14:1335–1338.