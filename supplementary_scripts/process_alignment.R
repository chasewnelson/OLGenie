#####################################################################################################
### Supplementary Material for Nelson CW, Ardern Z, Wei X, "OLGenie: Estimating Natural Selection to Predict Functional Overlapping Genes"
### Process BLAST alignment for OLGenie

library(tidyverse)
library(GenomicFeatures)
library(Biostrings)

input_file_path <- "/ENTER/YOUR/PATH/HERE" # e.g, "controls_CAL_min6/ss12/hepatitis_B_virus/NOL_X_1628-1813/J7CUJW3F014-Alignment.fasta"
(seqs_DSS <- readDNAStringSet(filepath = input_file_path))

# Extract the sequences as a vector
seqs_nt_vector <- as.data.frame(seqs_DSS)$x

# Remove gaps
seqs_nt_vector_noGap <- str_replace_all(string = seqs_nt_vector, pattern = "-", replacement = "")

# Save sequences in data frame with names
seq_nt_df <- tibble(seq = seqs_nt_vector_noGap, name = names(seqs_DSS))

# Keep unique (alleles)
seq_nt_df <- distinct(.data = seq_nt_df, seq, .keep_all = TRUE)

# Construct new set unique sequences as DSS
seqs_DSS_noGap <- DNAStringSet(x = seq_nt_df$seq)
names(seqs_DSS_noGap) <- seq_nt_df$name

# Translate
(seqs_aa_AASS <- Biostrings::translate(x = seqs_DSS_noGap, if.fuzzy.codon = 'solve'))
# Warnings are OK; they result from non-complete codons (to be removed)

# Find mid-sequence STOP codons
seqs_aa_vector <- as.data.frame(seqs_aa_AASS)$x
mid_seq_STOP_count <- stringr::str_count(string = seqs_aa_vector, pattern = "\\*\\w")
mid_seq_STOP_count_firstHit <- stringr::str_count(string = seqs_aa_vector[[1]], pattern = "\\*\\w")

# Find sequences that are multiples of 3
complete_codons <- width(seqs_DSS_noGap) %% 3 == 0

# Keep only seqs of sufficient coverage, defined as 70% the length of the first hit, following Hughes et al. (2005)
sufficient_coverage <- width(seqs_DSS_noGap) >= 0.7 * width(seqs_DSS_noGap[1])

# If HIV-1 env example, exclude if certain terms appear in the name
#name_to_exclude <- (stringr::str_detect(string = names(seqs_DSS_noGap), pattern = c("[Ss]imian"))) |
#  (stringr::str_detect(string = names(seqs_DSS_noGap), pattern = c("[Ss]ynthetic construct"))) |
#  (stringr::str_detect(string = names(seqs_DSS_noGap), pattern = c("[Cc]loning"))) |
#  (stringr::str_detect(string = names(seqs_DSS_noGap), pattern = c("T-cell")))

# Only keep sequences with ≤ the number of STOPs of first hit, exact multiple of 3, ≥70% length
seqs_to_keep <- ((mid_seq_STOP_count <= mid_seq_STOP_count_firstHit) & complete_codons & sufficient_coverage) # & (! name_to_exclude)

# Build DSS of keepers ("filtered")
(seqs_DSS_filtered <- seqs_DSS_noGap[seqs_to_keep])

# Translate keepers
(seqs_AASS_filtered <- Biostrings::translate(x = seqs_DSS_filtered, if.fuzzy.codon = 'solve'))

# Modify duplicate names with a unique numeric suffix
name_count_table <- table(names(seqs_AASS_filtered))
name_count_dup <- name_count_table > 1
dup_names <- name_count_table[name_count_dup]
name_is_dup <- names(seqs_AASS_filtered) %in% names(dup_names)

# If there is at least one duplicate name, append suffix
if(sum(name_is_dup) > 0) {
  # Generate suffices
  name_suffices <- seq(1:sum(name_count_table[name_count_table > 1]))
  
  # Add suffices
  new_names <- names(seqs_AASS_filtered)
  new_names[name_is_dup] <- paste0(new_names[name_is_dup], '_', name_suffices)
  
  # ASSIGN MODIFIED NAMES
  names(seqs_DSS_filtered) <- new_names
  names(seqs_AASS_filtered) <- new_names
}

### SAVE NUCLEOTIDE FASTA
output_file_path_nt <- input_file_path
output_file_path_nt <- str_replace(string = output_file_path_nt, pattern = "\\/[\\w\\-\\d]+.fasta", replacement = "\\/seqs_DSS_uniq_goodName_filtered.fa")
writeXStringSet(x = seqs_DSS_filtered, filepath = output_file_path_nt, format = 'fasta')

### SAVE AMINO ACID FASTA
output_file_path_pep <- input_file_path
output_file_path_pep <- str_replace(string = output_file_path_pep, pattern = "\\/[\\w\\-\\d]+.fasta", replacement = "\\/seqs_AASS_uniq_goodName_filtered.fa")
writeXStringSet(x = seqs_AASS_filtered, filepath = output_file_path_pep, format = 'fasta')

##########################################################################################
### NOW RUN MAFFT TO ALIGN AMINO ACIDS AT COMMAND LINE
# mafft seqs_AASS_uniq_goodName_filtered.fa > seqs_AASS_uniq_goodName_filtered_PALmafft.fa 2> seqs_AASS_uniq_goodName_filtered_PALmafft.out

##########################################################################################
### NOW RUN PAL2NAL TO IMPOSE AMINO ACID ALIGNMENT ON NUCLEOTIDES
# pal2nal.pl seqs_AASS_uniq_goodName_filtered_PAL.fa seqs_DSS_uniq_goodName_filtered.fa -output fasta > seqs_DSS_uniq_goodName_filtered_CAL.fa

##########################################################################################
### NOW RUN OLGENIE TO ESTIMATE dN/dS FOR OLG
# OLGenie.pl --fasta_file=seqs_DSS_uniq_goodName_filtered_PALmafft_CAL.fa --frame=sas12 --verbose > OLGenie_uniq_goodName_filtered_PALmafft.out

