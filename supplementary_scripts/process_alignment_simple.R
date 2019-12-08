#####################################################################################################
### Supplementary Material for Nelson CW, Ardern Z, Wei X, "OLGenie: Estimating Natural Selection to Predict Functional Overlapping Genes"
### Process BLAST alignment for OLGenie

library(tidyverse)
library(GenomicFeatures)
library(Biostrings)

input_file_path <- "/ENTER/YOUR/PATH/HERE" # "controls_simple/ss12/hepatitis_B_virus/NOL_X_1628-1813/J7CUJW3F014-Alignment.fasta"
(seqs_DSS <- readDNAStringSet(filepath = input_file_path))

# Extract the sequences as a vector
seqs_nt_vector <- as.data.frame(seqs_DSS)$x

# Remove gaps
seqs_nt_vector_noGap <- str_replace_all(string = seqs_nt_vector, pattern = "-", replacement = "")

# Determine names
seq_nt_names_df <- tibble(seq = seqs_nt_vector_noGap, name = names(seqs_DSS))

# Keep unique (alleles)
seq_nt_distinct_df <- distinct(.data = seq_nt_names_df, seq, .keep_all = TRUE)

# Construct new set of unique sequences as DSS
seqs_DSS_noGap <- DNAStringSet(x = seq_nt_distinct_df$seq)
names(seqs_DSS_noGap) <- seq_nt_distinct_df$name

# Translate
(seqs_aa_AASS <- Biostrings::translate(x = seqs_DSS_noGap, if.fuzzy.codon = 'solve'))
# Warnings are OK; they result from non-complete codons (to be removed)

# Find mid-sequence STOP codons
seqs_aa_vector <- as.data.frame(seqs_aa_AASS)$x
mid_seq_STOP <- stringr::str_detect(string = seqs_aa_vector, pattern = "\\*\\w")
#mid_seq_STOP_count <- stringr::str_count(string = seqs_aa_vector, pattern = "\\*\\w") # optional for counting

# Find sequences that are multiples of 3; a formality, as the query should be a full codon set
complete_codons <- width(seqs_DSS_noGap) %% 3 == 0

# Keep sufficient sequence coverage, which we will define as an exact match to the first hit
sufficient_coverage <- width(seqs_DSS_noGap) == width(seqs_DSS_noGap[1])

# Only keep sequences without mid-sequence STOP,  multiple of 3, exact length
seqs_to_keep <- (! mid_seq_STOP & complete_codons & sufficient_coverage)
#seqs_to_keep <- (mid_seq_STOP < 2 & complete_codons & sufficient_coverage) # optional for counting criterion

# Translate keepers
(seqs_DSS_noGap_noStop <- seqs_DSS_noGap[seqs_to_keep])

# Modify duplicate names with a unique numeric suffix
name_count_table <- table(names(seqs_DSS_noGap_noStop))
name_count_dup <- name_count_table > 1
dup_names <- name_count_table[name_count_table > 1]
name_is_dup <- names(seqs_DSS_noGap_noStop) %in% names(dup_names)

# Generate suffices
name_suffices <- seq(1:sum(name_count_table[name_count_table > 1]))

# Add suffices
new_names <- names(seqs_DSS_noGap_noStop)
new_names[name_is_dup] <- paste0(new_names[name_is_dup], '_', name_suffices)

# ASSIGN MODIFIED NAMES
names(seqs_DSS_noGap_noStop) <- new_names

# Write the processed nucleotide alignment file
output_file_path <- input_file_path
output_file_path <- str_replace(string = output_file_path, pattern = "\\/[\\w\\-\\d]+.fasta", replacement = "\\/seqs_DSS_noGap_noStop_fullCoverage_UNIQUE.fa")
#writeXStringSet(x = seqs_DSS_noGap_noStop, filepath = output_file_path, format = 'fasta')


### NOW RUN OLGENIE AT THE COMMAND LINE ##
#OLGenie.pl --fasta_file=[alignment].fa --frame=[frame_type] --verbose > OLGenie.out

