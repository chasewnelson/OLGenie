<img src="https://github.com/chasewnelson/OLGenie/blob/master/OLGenie_logo.png?raw=true" title="OLGenie logo by Mitch Lin" alt="OLGenie logo by Mitch Lin" align="middle">

**OLGenie** is a Perl program for detecting selection with *d*<sub>N</sub>/*d*<sub>S</sub> in overlapping genes (OLGs). It relies on no external dependencies, facilitating maximum portability. Just download and run.

To test the software with the [example data](#examples), execute the program at the Unix command line or Mac Terminal as follows:

### FORMAT:

	OLGenie.pl --fasta_file=<alignment>.fasta --frame=<frame> \
	--output_file=<OLGenie_codon_results>.tsv --verbose > OLGenie_log.txt

Find some real [examples](#examples) below. Check out our <a target="_blank" rel="noopener noreferrer" href="https://github.com/chasewnelson/OLGenie">preprint on bioRxiv</a>.

## <a name="contents"></a>Contents

* [Description](#description)
* [How it Works](#how-it-works)
* [Options](#options)
* [Examples](#examples)
* [Output](#output)
* [Troubleshooting](#troubleshooting)
* [Acknowledgments](#acknowledgments)
* [Citation](#citation)
* [Contact](#contact)
* [References](#references)

## <a name="description"></a>Description
Given the codon triplet basis of the genetic code, and the two antiparallel strands of the DNA molecule, a single segment of DNA has the potential to encode six reading frames: three in the forward (sense) direction and three in the reverse (antisense) direction. This allows for the possibility that two or more genes may overlap the same nucleotide positions in a genome. Indeed, a substantial fraction of genes in taxa ranging from viruses to humans may encode overlapping gene (**OLG**) pairs, running in either the same (**ss**; sense-sense) or opposite (**sas**; sense-antisense) directions (*e.g.*, see Pavesi *et al.* 2018 and Sabath 2009). We use the nomenclature of Wei and Zhang (2015), referring to these overlapping frames as ss12, ss13, sas11, sas12, or sas13, where the first number refers to the codon position in the reference gene (ORF1), and the second number refers to the codon position in the alternative (overlapping) gene (ORF2):

<img src="https://github.com/chasewnelson/OLGenie/blob/master/OLGs_all_frames.png" alt="sas12 logo" align="middle">

For example, in sas12, genes overlap in a sense-antisense relationship such that position 1 of codons in the sense (reference) gene correspond to position 2 of codons in the reverse strand (overlapping; alternative) gene. In other words, the sense gene's first codon position overlaps the antisense gene's second codon position:
	
<img src="https://github.com/chasewnelson/OLGenie/blob/master/sas12_figure.png" alt="sas12 logo" align="middle" width="2500px">

It is common to detect natural selection in a DNA sequence alignment using *d*<sub>N</sub>/*d*<sub>S</sub>, *i.e.*, the ratio of nonsynonymous (changes the amino acid) to synonymous (does not change the amino acid) differences per site. Unfortunately, standard methods for estimating *d*<sub>N</sub>/*d*<sub>S</sub> do not apply to OLGs, because a mutation that is synonymous in one frame may be nonsynonymous in another, and *vice versa*. Although some methods for detecting natural selection in OLGs have been developed, they are generally computationally intensive and limited in utility (*e.g.*, Wei and Zhang 2015; Sabath et al. 2008). Thus, it is necessary to develop improved approaches for detecting selection in OLGs that can be implemented using genomic-scale data.

OLGenie represents a simplification and extension of the method of Wei and Zhang (2015), utilizing the approach of <a target="_blank" href="https://github.com/chasewnelson/SNPGenie">SNPGenie</a> (Nelson et al. 2015), and tailored for detecting selection in OLGs. The method considers the effects of mutations in the overlapping frame to determine the numerator (number of differences) and denominator (number of sites) of *d*<sub>N</sub> and *d*<sub>S</sub>.  For example, *d*<sub>N</sub> is usually calculated as the mean number of nonsynonymous (amino acid changing) nucleotide differences per nonsynonymous nucleotide site, and *d*<sub>S</sub> is similarly calculated for synonymous (silent) differences and sites. In order to control for the possibility that synonymous sites in the frame of interest may be under selection in the alternative overlapping reading frame, Wei-Zhang further considers the expanded measures *d*<sub>NN</sub>, *d*<sub>SN</sub>, *d*<sub>NS</sub>, and *d*<sub>SS</sub>, where the first subscript refers to the reference (sense) gene (ORF1), and the second to the alternative (sense or antisense) OLG (ORF2). For example, *d*<sub>SN</sub> refers to the mean number of differences per site that are synonymous in the reference frame but nonsynonymous in the alternative frame (*i.e.*, SN). Using these measures, it is possible to estimate *d*<sub>N</sub>/*d*<sub>S</sub> for the reference gene using *d*<sub>NN</sub>/*d*<sub>SN</sub> or *d*<sub>NS</sub>/*d*<sub>SS</sub>, and to estimate *d*<sub>N</sub>/*d*<sub>S</sub> for the alternative gene as *d*<sub>NN</sub>/*d*<sub>NS</sub> or *d*<sub>SN</sub>/*d*<sub>SS</sub>, *i.e.*, the subscript in the alternative OLG is held constant to control for OLG effects.

For more details, please refer to our manuscript.

## <a name="how-it-works"></a>How it Works
OLGenie is written in Perl with no dependencies for maximum portability (just download and run). The program examines a user-provided FASTA alignment of one protein-coding gene region from the reference gene point of view. This means that the alignment begins at the first site of a reference gene codon, and ends at the last (third) site of a reference gene codon. 

The choice of which gene to consider the reference gene is arbitrary. Typically, the **reference gene** is the gene whose functional status is known, while the functionality of the **alternate gene** may be in question. Thus, in practice, the reference gene (mother/ORF1 gene) is usually larger than the alternate gene (daughter/ORF2 gene), and the alternate gene is either partially or fully embedded within the reference.

After reading in the user-provided alignment, OLGenie calculates the number of NN, SN, NS, and SS sites and differences, reporting the mean number of each for all pairwise comparisons. This is done separately for each focal reference codon by considering all unique nonamer (9nt) alleles of which the reference codon is the center, of which 6nt constitute a minimum overlapping unit: one reference gene codon and its two overlapping alternate gene codons. (Note that sas13 is unique in that one reference codon overlaps exactly one alternate codon.) These tasks are sufficiently fast that no parallelism is necessary at the level of the single gene alignment. For datasets with many genes, the user can implement their own parallelization by running numerous alignments (genes) simultaneously.

After results are obtained for each focal codon in the alignment, significant deviations from the null expectation (*d*<sub>N</sub> - *d*<sub>S</sub> = 0) may be tested using a *Z*-test, where the standard error is estimated using bootstrapping (focal codon unit). Don't worry — we provide scripts to do it all!

## <a name="options"></a>Options
Call **OLGenie** using the following options: 

* `--fasta_file` (**REQUIRED**): a FASTA file containing multiple aligned sequences of one coding sequence. The entire coding sequence will be analyzed as an overlapping gene (OLG), with no non-overlapping codons, even if only part (or none) of the alignment constitues a true OLG. The frame of the alignment must be the frame of the reference gene (ORF1) (see the --frame option). If the user wishes to align their own sequences, it is recommended to translate the gene sequences, align at the amino acid level, and then impose the amino acid alignment on the DNA alignment to preserve complete codons. (If you need a tool to help with this, see align\_codon2aa.pl at <a target="_blank" href="https://github.com/chasewnelson/EBT">Evolutionary Bioinformatics Toolkit</a>.)
* `--frame` (**REQUIRED**): the frame relationship of the overlapping gene (OLG): ss12, ss13, sas11, sas12, or sas13 (see figures and description above).
* `--output_file` (*OPTIONAL*): name of the TAB-delimited output file to be placed in the working directory unless a full path name is given. If not specified, a file will be printed in the working directory by the name OLGenie_codon_results.txt (DEFAULT).
* `--verbose` (*OPTIONAL*): tell OLGenie to report all unique nonamers (9-mers) overlapping each reference codon, along with their counts, in the output file. May lead to large output files in cases with many and/or divergent sequences. If not specified, verbose output will not be reported (DEFAULT).

## <a name="examples"></a>EXAMPLES

Example input and output files are available in the `EXAMPLE_INPUT` and `EXAMPLE_OUTPUT` directories at this GitHub page, where reproducible examples are numbered (*e.g.*, **output_example1.txt**). Note that, if your input file(s) (*e.g.*, **alignment.fasta**) are not in the working directory (*i.e.*, where your Terminal is currently operating), you will need to specify the full path of the file name (*e.g.*, **/Users/ohta/Desktop/OLGenie\_practice/alignment.fasta**). Also note that, in the examples below, a `\` is used simply to continue the previous command on the line.

### EXAMPLE 1: A SIMPLE RUN
Note that this is a 'real' example and may take up to 60 seconds!

	OLGenie.pl --fasta_file=HIV1_env_BLAST.fa --frame=sas12 > example1.out

### EXAMPLE 2: VERBOSE OUTPUT TO A USER-SPECIFIED FILE
Remember to replace the `--output_file` path with a location that exists on your machine.

	OLGenie.pl --fasta_file=HIV1_env_BLAST.fa --frame=sas12 \
	--output_file=/Users/ohta/Desktop/OLGenie_codon_results_ex2.txt --verbose > example2.out

### EXAMPLE 3: TESTING FOR SIGNIFICANCE WITH BOOTSTRAPPING
Use our script `OLGenie_bootstrap.R`. We provide this script separately so that users can take advantage of the accessible statistical resources offerred by R without having to install Perl modules. Just make sure the R packages `readr` and `boot` have been installed (e.g., by calling `install.packages("readr")` and `install.packages("boot")` at the R console).

Call the script with the following (unnamed) arguments (in this order):

1. The name/path of the file containing the codon results file from the OLGenie analysis.
2. The number of bootstrap replicates to perform (typically 1,000 or 10,000).
3. The number of parallel processes (CPUs) to use when bootstrapping. A typical personal laptop computer can utilize 4-8 CPUs, while a high performance computing cluster might provide access to 10s or 100s.

Thus, the format is:

	OLGenie_bootstrap.R <codon_results>.txt <num_bootstraps> <num_cpus> > <output>.out

For example, try the following using the results from **Example 2**:

	OLGenie_bootstrap.R OLGenie_codons_results_ex2.txt 1000 4 > example3.out

This produces TAB-delimited output, as described in the [bootstrap output](#bootstrap-output) section.

## <a name="output"></a>Output

**OLGenie** will output the following data:

### <a name="standard-output"></a>Standard Output

At the command line (Terminal), OLGenie will first report the data and time, the file and frame relationship used in the analysis, and any warning messages. Following completion of the analysis, OLGenie reports the following summary statistics:

* **Mean numbers of sites and differences**: the total numbers of NN, SN, NS, and SS sites and differences for the entire alignment, obtained by summing the results for all codons.
* **Mean substitution rates (between-species) or nucleotide diversities (within-species):**: OLGenie's estimates of *d*<sub>NN</sub>, *d*<sub>SN</sub>, *d*<sub>NS</sub>, and *d*<sub>SS</sub> for the entire alignment, calculated as (\*\_diffs / \*\_sites) for each site type.
* **dN/dS estimates**: OLGenie's estimates of *d*<sub>N</sub>/*d*<sub>S</sub> for the reference gene (*d*<sub>NN</sub>/*d*<sub>SN</sub>, *d*<sub>NS</sub>/*d*<sub>SS</sub>) and alternate gene (*d*<sub>NN</sub>/*d*<sub>NS</sub> and *d*<sub>SN</sub>/*d*<sub>SS</sub>) for the entire alignment.

### <a name="codon-output-file"></a>Codon Results Output File

OLGenie will report codon-by-codon results in the file **OLGenie\_codon\_results.txt** (or any file specified with the `--output_file` option). The columns contain the following information:

* `codon_num`: the codon position in the alignment, starting at codon 2 and ending at the penultimate codon. The first and last codons are excluded because their values cannot be estimated, as one of their overlapping (alternate gene) codons is unknown, occurring before or after the alignment begins or ends, respectively.
* `ref_codon_maj`: the major (most common) allele for the reference gene codon at this position.
* `alt_codon1_maj`: the major (most common) allele for the alternate gene codon overlapping the beginning (5' side) of the reference codon at this position.
* `alt_codon2_maj`: the major (most common) allele for the alternate gene codon overlapping the end (3' side) of the reference codon at this position. Note that only `alt_codon1_maj` will be reported for the sas13 frame, since OLG codons form one-to-one overlaps in this frame.
* `nonamers`: only included when using the `--verbose` option. This column contains all unique nonamer (9nt) alleles occuring at this position, with the reference focal codon at the center. Different alleles are separated using the colon (`:`) delimiter.
* `nonamer_counts`: only included when using the `--verbose` option. This column contains the counts (number of sequences) having each unique nonamer (9nt) alleles at this position, in the same order given in the `nonamers` column. Values for different alleles are separated using the colon (`:`) delimiter.
* `multiple_variants`: whether the nonamer at this position contains more than one nucleotide variant. If so, the OLGenie method may underestimate *d*<sub>S</sub> at this position. The *d*<sub>N</sub>/*d*<sub>S</sub> ratio will constitute a conservative test of purifying (negative) selection, but positive (Darwinian) selection should be inferred with caution.
* `NN_sites`: the number of sites (*i.e.*, possible nucleotide changes) that are nonsynonymous in both the reference and alternate genes at this reference codon.
* `SN_sites`: the number of sites (*i.e.*, possible nucleotide changes) that are synonymous in the reference gene but nonsynonymous in the alternate gene at this reference codon.
* `NS_sites`: the number of sites (*i.e.*, possible nucleotide changes) that are nonsynonymous in the reference gene but synonymous in the alternate gene at this reference codon.
* `SS_sites`: the number of sites (*i.e.*, possible nucleotide changes) that are synonymous in both the reference and alternate genes at this reference codon.
* `NN_diffs`: the number of differences (*i.e.*, observed nucleotide changes) that are nonsynonymous in both the reference and alternate genes at this reference codon.
* `SN_diffs`: the number of differences (*i.e.*, observed nucleotide changes) that are synonymous in the reference gene but nonsynonymous in the alternate gene at this reference codon.
* `NS_diffs`: the number of differences (*i.e.*, observed nucleotide changes) that are nonsynonymous in the reference gene but synonymous in the alternate gene at this reference codon.
* `SS_diffs`: the number of differences (*i.e.*, observed nucleotide changes) that are synonymous in both the reference and alternate genes at this reference codon.

Note that any desired estimate of *d*<sub>N</sub>, *d*<sub>S</sub>, or their ratio can be obtained for any subregion of the alignment by summing the appropriate numbers of sites and differences and performing the appropriate calculations. For example, to calculate the alternate gene *d*<sub>N</sub>/*d*<sub>S</sub> = *d*<sub>SN</sub>/*d*<sub>SS</sub> ratio for a 25-codon window within an alignment:

1. Calculate *d*<sub>SN</sub> as sum(`SN_diffs`)/sum(`SN_sites`) for those 25 codons;
2. Calculate *d*<sub>SS</sub> as sum(`SS_diffs`)/sum(`SS_sites`) for those 25 codons; and
3. Calculate the *d*<sub>SN</sub>/*d*<sub>SS</sub> value.

### <a name="bootstrap-output"></a>Bootstrap Output

Significant deviations from neutrality (*d*<sub>N</sub> - *d*<sub>S</sub> = 0) can be detected using a *Z*-test, where the standard error *d*<sub>N</sub> - *d*<sub>S</sub> is estimated using bootstrapping (reference codon unit) (Nei and Kumar 2000). Consider using our R script, `OLGenie_bootstrap.R` (see [examples](#examples)). The produces four lines of output, one for each of the four ratios: *d*<sub>NN</sub>/*d*<sub>SN</sub>, *d*<sub>NN</sub>/*d*<sub>NS</sub>, *d*<sub>NS</sub>/*d*<sub>SS</sub>, and *d*<sub>SN</sub>/*d*<sub>SS</sub>. Columns of values are given in the following order (numbered here for clarity, as these headers do not appear in the output):

1. `num_codons`: the total number of codons examined.
2. `NN_sites`: see the description of the [codon output file](#codon-output-file).
3. `SN_sites`: see the description of the [codon output file](#codon-output-file).
4. `NS_sites`: see the description of the [codon output file](#codon-output-file).
5. `SS_sites`: see the description of the [codon output file](#codon-output-file).
6. `NN_diffs`: see the description of the [codon output file](#codon-output-file).
7. `SN_diffs`: see the description of the [codon output file](#codon-output-file).
8. `NS_diffs`: see the description of the [codon output file](#codon-output-file).
9. `SS_diffs`: see the description of the [codon output file](#codon-output-file).
10. `ratio`: the ratio being estimated on this line: dNNdSN denotes *d*<sub>NN</sub>/*d*<sub>SN</sub>; dNNdNS denotes *d*<sub>NN</sub>/*d*<sub>NS</sub>; dNSdSS denotes *d*<sub>NS</sub>/*d*<sub>SS</sub>; and dSNdSS denotes *d*<sub>SN</sub>/*d*<sub>SS</sub>.
11. `site_rich_ratio`: whether this is the most site-rich ratio (**TRUE** or **FALSE**). Note that, for sas12, the more accurate ratios (*d*<sub>NS</sub>/*d*<sub>SS</sub> and *d*<sub>SN</sub>/*d*<sub>SS</sub>) are not the most site-rich.
12. `gene`: whether this line is an estimate of *d*<sub>N</sub>/*d*<sub>S</sub> for the reference gene (**ORF1**) or the alternate gene (**ORF2**).
13. `num_replicates`: number of bootstrap replicates performed.
14. `dN`: the estimate of *d*<sub>N</sub> (numerator of `ratio`).
15. `dS`: the estimate of *d*<sub>S</sub> (denominator of `ratio`).
16. `dNdS`: the estimate of *d*<sub>N</sub>/*d*<sub>S</sub> (value of `ratio`).
17. `dN_m_dS`: the estimate of *d*<sub>N</sub> - *d*<sub>S</sub>.
18. `boot_dN_SE`: the standard error of *d*<sub>N</sub>, estimated by bootstrapping.
19. `boot_dS_SE`: the standard error of *d*<sub>S</sub>, estimated by bootstrapping.
20. `boot_dN_over_dS_SE`: the standard error of *d*<sub>N</sub>/*d*<sub>S</sub>, estimated by bootstrapping.
21. `boot_dN_over_dS_P`: the P value of a deviation from *d*<sub>N</sub>/*d*<sub>S</sub> = 1 (two-sided; *Z*-test).
22. `boot_dN_m_dS_SE`: the standard error of *d*<sub>N</sub> - *d*<sub>S</sub>, estimated by bootstrapping.
23. `boot_dN_m_dS_P`: the P value of a deviation from *d*<sub>N</sub> - *d*<sub>S</sub> = 0 (two-sided; *Z*-test).

## <a name="troubleshooting"></a>Troubleshooting

If you have questions about **OLGenie**, please click on the <a target="_blank" href="https://github.com/chasewnelson/OLGenie/issues">Issues</a> tab at the top of this page and begin a new thread, so that others might benefit from the discussion. Common questions will be addressed in this section.

## <a name="acknowledgments"></a>Acknowledgments
**OLGenie** was written with support from a Gerstner Scholars Fellowship from the Gerstner Family Foundation at the American Museum of Natural History to C.W.N. (2016-2019), and is maintained with support from a 中央研究院 Academia Sinica Postdoctoral Research Fellowship (2019-2021). The logo image was designed by Mitch Lin (2019); copyright-free DNA helix obtained from Pixabay. Thanks to Reed Cartwright, Jim Hussey, Michael Lynch, Sergios Orestis-Kolokotronis, Wen-Hsiung Li, Apurva Narechania, Sally Warring, Jeff Witmer, Meredith Yeager, Jianzhi (George) Zhang, Martine Zilversmit, and the Sackler Institute for Comparative Genomics workgroup for discussion along the way.

## <a name="citation"></a>Citation

When using this software, please refer to and cite:

>Nelson CW, Ardern Z, Wei X. OLGenie: detecting natural selection to identify functional overlapping genes. bioRxiv doi: <a target="_blank" rel="noopener noreferrer" href="https://github.com/chasewnelson/OLGenie">TBA</a>

and this page:

>https://github.com/chasewnelson/OLGenie

## <a name="contact"></a>Contact
If you have questions about **OLGenie**, please click on the <a target="_blank" href="https://github.com/chasewnelson/OLGenie/issues">Issues</a> tab at the top of this page and begin a new thread, so that others might benefit from the discussion.

Other correspondence should be addressed to Chase W. Nelson: 

* cnelson <**AT**> amnh <**DOT**> org

## <a name="references"></a>References
* Nei M, Gojobori T. 1986. <a target="_blank" href="https://academic.oup.com/mbe/article/3/5/418/988012">Simple methods for estimating the numbers of synonymous and nonsynonymous nucleotide substitutions</a>. *Molecular Biology and Evolution* **3**(5):418-426.
* Nei M, Kumar S. 2000. *Molecular Evolution and Phylogenetics*. New York, NY: Oxford University Press.
* Nelson CW, Moncla LH, Hughes AL. 2015. <a target="_blank" href="https://academic.oup.com/bioinformatics/article/31/22/3709/241742">SNPGenie: estimating evolutionary parameters to detect natural selection using pooled next-generation sequencing data</a>. *Bioinformatics* **31**(22):3709-3711.
* Pavesi A, Vianelli A, Chirico N, Bao Y, Blinkova O, Belshaw R, Firth A, Karlin D. 2018. <a target="_blank" href="https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0202513">Overlapping genes and the proteins they encode differ significantly in their sequence composition from non-overlapping genes</a>. *PLoS ONE* **13**:e0202513.
* Sabath N. 2009. <a target="_blank" href="http://nsmn1.uh.edu/dgraur/niv/Sabath_PhD_Thesis.pdf">*Molecular Evolution of Overlapping Genes*</a>. ProQuest Dissertation #3405062.* Sabath N, Landan G, Graur D. 2008. <a target="_blank" href="https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0003996">A method for the simultaneous estimation of selection intensities in overlapping genes</a>. *PLoS One* **3**(12):e3996.* Wei X, Zhang J. 2015. <a target="_blank" href="https://academic.oup.com/gbe/article/7/1/381/604609">A simple method for estimating the strength of natural selection on overlapping genes</a>. *Genome Biology and Evolution* **7**(1):381-390.