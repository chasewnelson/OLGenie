######################################################################################################
### Supplementary Material for Nelson CW, Ardern Z, Wei X, "OLGenie: Estimating Natural Selection to Predict Functional Overlapping Genes"
### Analyze OLGenie control results

library(tidyverse)
library(patchwork)
library(RColorBrewer)


######################################################################################################
## LOAD THE RESULTS
results <- read_tsv("/Users/cwnelson88/Desktop/OLG/tableS2_practice.txt") # Supplementary Table S1


#####################################################################################################
# ONLY ANALYZE ORF2 (ALTERNATE GENE) AND THE MOST ACCURATE RATIO (SITE RICH; dNN/dNS)
results_to_analyze <- filter(results, gene == 'ORF2', site_rich_ratio == TRUE, is.finite(dNdS), is.finite(boot_dN_m_dS_SE))


#####################################################################################################
# Determine empirical P-value cutoffs

# Want to form the ROC by recording the TP and FP rate for every P-value cutoff
alpha_vector <- sort(unique(results_to_analyze$P_value)) # sort all P-values
empirical_P0001 <- 0
empirical_P001 <- 0
empirical_P01 <- 0
empirical_P025 <- 0
empirical_P05 <- 0
empirical_P1 <- 0

for (max_P in sort(alpha_vector)) { # MUST BE SORTED; already is, but re-do anyway
  # Use for P = 0.0001
  if(empirical_P0001 == 0 && max_P > 0.0001) {
    empirical_P0001 <- max_P
  }
  
  # Use for P = 0.001
  if(empirical_P001 == 0 && max_P > 0.001) {
    empirical_P001 <- max_P
  }
  
  # Use for P = 0.01
  if(empirical_P01 == 0 && max_P > 0.01) {
    empirical_P01 <- max_P
  }
  
  # Use for P = 0.025
  if(empirical_P025 == 0 && max_P > 0.025) {
    empirical_P025 <- max_P
  }
  
  # Use for P = 0.05
  if(empirical_P05 == 0 && max_P > 0.05) {
    empirical_P05 <- max_P
  }
  
  # Use for P = 0.1
  if(empirical_P1 == 0 && max_P > 0.1) {
    empirical_P1 <- max_P
  }
}


#####################################################################################################
### ROC ###
#####################################################################################################

#####################################################################################################
# ANALYZE ALL DATA SUBSETS AND P-VALUE CUT-OFFS

# Set up results tables for ROC
FP_vector <- vector()
TN_vector <- vector()
TP_vector <- vector()
FN_vector <- vector()
FP_rate_vector <- vector()
TN_rate_vector <- vector()
TP_rate_vector <- vector()
FN_rate_vector <- vector()
prevalence_vector <- vector()
precision_vector <- vector()
accuracy_vector <- vector()
max_P_vector <- vector()
min_width_vector <- vector()
max_dNdS_vector <- vector()
sample_size_vector <- vector()

min_width_cutoffs <- c(0, 50, 100, 150, 200, 300) # NUCLEOTIDES; similar cutoffs as Schlub et al. 2019
max_dNdS_cutoffs <- c(seq(0.1, 0.9, 0.1), Inf) # max dN/dS values from 0.1 to 0.9

# LOOP ALL P-VALUES, WIDTHS, AND MAX dN/dS VALUES
for (max_P in alpha_vector) {
  
  for (min_width in min_width_cutoffs) { # NUCLEOTIDE, not codon, as Schlub et al. 2019
    
    for (max_dNdS in max_dNdS_cutoffs) {
      
      sample_size <- nrow(results_to_analyze[results_to_analyze$dNdS <= max_dNdS & results_to_analyze$width >= min_width, ])
      
      NEG_count <- nrow(results_to_analyze[results_to_analyze$control_type == 'NEG' & results_to_analyze$dNdS <= max_dNdS & results_to_analyze$width >= min_width, ])
      NEG_sig <- nrow(results_to_analyze[results_to_analyze$control_type == 'NEG' & results_to_analyze$dNdS <= max_dNdS & results_to_analyze$width >= min_width &
                                           results_to_analyze$P_value < max_P, ])
      POS_count <- nrow(results_to_analyze[results_to_analyze$control_type == 'POS' & results_to_analyze$dNdS <= max_dNdS & results_to_analyze$width >= min_width, ])
      POS_sig <- nrow(results_to_analyze[results_to_analyze$control_type == 'POS' & results_to_analyze$dNdS <= max_dNdS & results_to_analyze$width >= min_width &
                                           results_to_analyze$P_value < max_P, ])
      
      prevalence <- POS_count / sample_size
      
      
      # Add to table columns
      FP_vector <- c(FP_vector, NEG_sig)
      TN_vector <- c(TN_vector, (NEG_count - NEG_sig))
      TP_vector <- c(TP_vector, POS_sig)
      FN_vector <- c(FN_vector, (POS_count - POS_sig))
      FP_rate_vector <- c(FP_rate_vector, NEG_sig / NEG_count)
      TN_rate_vector <- c(TN_rate_vector, (NEG_count - NEG_sig) / NEG_count) # specificity
      TP_rate_vector <- c(TP_rate_vector, POS_sig / POS_count) # sensitivity; recall
      FN_rate_vector <- c(FN_rate_vector, (POS_count - POS_sig) / POS_count)
      
      max_P_vector <- c(max_P_vector, max_P)
      min_width_vector <- c(min_width_vector, min_width)
      max_dNdS_vector <- c(max_dNdS_vector, max_dNdS)
      
      sample_size_vector <- c(sample_size_vector, sample_size)
      
      ###############################################################################################
      # ACCURACY & CONFUSION MATRIX
      precision <- POS_sig / (POS_sig + NEG_sig)
      accuracy <- ((NEG_count - NEG_sig) + POS_sig) / sample_size
      
      prevalence_vector <- c(prevalence_vector, prevalence)
      precision_vector <- c(precision_vector, precision)
      accuracy_vector <- c(accuracy_vector, accuracy)
    }
  }
}

# SAVE
ROC_table <- data.frame(max_P = max_P_vector, min_width = min_width_vector, max_dNdS = max_dNdS_vector,
                        sample_size = sample_size_vector,
                        FP = FP_vector, TN = TN_vector, TP = TP_vector, FN = FN_vector,
                        FP_rate = FP_rate_vector, TN_rate = TN_rate_vector, TP_rate = TP_rate_vector, FN_rate = FN_rate_vector,
                        prevalence = prevalence_vector, precision = precision_vector, accuracy = accuracy_vector)
#write_tsv(x = ROC_table, path = "<ENTER/YOUR/PATH/HERE>") # Table SX


#####################################################################################################
# CALCULATE AUC (area under curve) for all conditions
ROC_min_width_vector <- vector()
ROC_max_dNdS_vector <- vector()
ROC_sample_size_vector <- vector()
ROC_num_species_vector <- vector()
ROC_AUC_vector <- vector()

for (this_min_width in min_width_cutoffs) { # NUCLEOTIDE, not codon

  for (this_max_dNdS in max_dNdS_cutoffs) {
    
    this_ROC <- filter(ROC_table, min_width == this_min_width, max_dNdS == this_max_dNdS)
    this_AUC <- sum(this_ROC$TP_rate[-1] * diff(this_ROC$FP_rate)) # the [-1] removed the first element
    
    ROC_min_width_vector <- c(ROC_min_width_vector, this_min_width)
    ROC_max_dNdS_vector <- c(ROC_max_dNdS_vector, this_max_dNdS)
    ROC_sample_size_vector <- c(ROC_sample_size_vector, this_ROC$sample_size[1]) # all same for this combination; take first
    ROC_num_species_vector <- c(ROC_num_species_vector, this_ROC$num_species[1]) # all same for this combination; take first
    ROC_AUC_vector <- c(ROC_AUC_vector, this_AUC)
    
  }
}

AUC_results <- data.frame(min_width = ROC_min_width_vector, max_dNdS = ROC_max_dNdS_vector,
                          sample_size = ROC_sample_size_vector, AUC = ROC_AUC_vector)

# SAVE
#write_tsv(AUC_results, path = "<ENTER/YOUR/PATH/HERE>")

