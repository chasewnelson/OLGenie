############################################################################################################
### Supplementary Material for Nelson CW, Ardern Z, Wei X, "OLGenie: Estimating Natural Selection to Predict Functional Overlapping Genes"
### HIV-1 env, Synplot2 and OLGenie sliding window analyses

library(caTools) # for runmean, for synplot
library(tidyverse)
library(boot) # for bootstrapping, for OLGenie


############################################################################################################
############################################################################################################
### GLOBAL VARIABLES
win_size <- 25
step_size <- 1 # only explicitly needed for OLGenie
NBOOTSTRAPS <- 1000
NCPUS <- 4


############################################################################################################
############################################################################################################
### BOOTSTRAP FUNCTION (dN - dS) for RESAMPLING CODON UNIT
dNdS_diff_boot_fun <- function(codon_results, numerator, denominator, num_replicates, num_cpus) {
  
  # Function for dN
  dN_function <- function(D, indices) {
    dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
    return(dN)
  }
  
  # Function for dN
  dS_function <- function(D, indices) {
    dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
    return(dS)
  }
  
  # Function for dN - dS
  dN_m_dS_function <- function(D, indices) {
    dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
    dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
    dN_m_dS <- dN - dS
    return(dN_m_dS)
  }
  
  # Function for dN/dS
  dN_over_dS_function <- function(D, indices) {
    dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
    dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
    dN_over_dS <- dN / dS
    return(dN_over_dS)
  }
  
  # CREATE FUNCTION FOR dN/dS TO CALCULATE ITS SE
  
  (dN <- sum(as.vector(codon_results[ , paste0(numerator, "_diffs")]), na.rm = TRUE) / sum(as.vector(codon_results[ , paste0(numerator, "_sites")]), na.rm = TRUE))
  (dS <- sum(as.vector(codon_results[ , paste0(denominator, "_diffs")]), na.rm = TRUE) / sum(as.vector(codon_results[ , paste0(denominator, "_sites")]), na.rm = TRUE))
  (dNdS <- dN / dS)
  
  # Run the BOOTSTRAPS
  # boot dN
  (boot_dN <- boot(data = codon_results, R = num_replicates, statistic = dN_function, ncpus = num_cpus)) # parallel = 'multicore', 
  (dN <- boot_dN$t0)
  (boot_dN_SE <- sd(boot_dN$t))
  
  # boot dS
  (boot_dS <- boot(data = codon_results, R = num_replicates, statistic = dS_function, ncpus = num_cpus)) # parallel = 'multicore', 
  (dS <- boot_dS$t0)
  (boot_dS_SE <- sd(boot_dS$t))
  
  # boot dN - dS
  (boot_dN_m_dS <- boot(data = codon_results, R = num_replicates, statistic = dN_m_dS_function, ncpus = num_cpus)) # parallel = 'multicore', 
  (dN_m_dS <- boot_dN_m_dS$t0)
  (boot_dN_m_dS_SE <- sd(boot_dN_m_dS$t))
  (boot_dN_m_dS_Z <- dN_m_dS / boot_dN_m_dS_SE)
  (boot_dN_m_dS_P <- 2 * pnorm(-abs(boot_dN_m_dS_Z)))
  
  # boot dN/dS
  (boot_dN_over_dS <- boot(data = codon_results, R = num_replicates, statistic = dN_over_dS_function, ncpus = num_cpus)) # parallel = 'multicore', 
  (dN_over_dS <- boot_dN_over_dS$t0)
  (boot_dN_over_dS_SE <- sd(boot_dN_over_dS$t))
  (boot_dN_over_dS_Z <- dN_over_dS / boot_dN_over_dS_SE)
  (boot_dN_over_dS_P <- 2 * pnorm(-abs(boot_dN_over_dS_Z)))
  
  return(paste(num_replicates, dN, dS, dNdS, dN_m_dS, boot_dN_SE, boot_dS_SE, boot_dN_over_dS_SE, boot_dN_over_dS_P, boot_dN_m_dS_SE, boot_dN_m_dS_P, sep = "\t"))
}


############################################################################################################
############################################################################################################
### SYNPLOT

### INPUT
#The columns of the plots files are:
#1)  codon number (alignment coordinates)
#2)  obs:    observed number of substitutions (summed over all pairs)
#3)  exp:    expected number of substitutions (summed over all pairs)
#4)  stddev: standard deviation (estimated from the null model)

### CHOOSE ONE (1) SYNPLOT DATASET TO ANALYZE
# Cassan M group
synplot_data <- read_tsv(file = "HIV1_env_CassanM_Synplot2.txt")
# Cassan non-M groups
synplot_data <- read_tsv(file = "HIV1_env_CassanNonM_Synplot2.txt")
# BLAST whole-env
synplot_data <- read_tsv(file = "HIV1_env_BLAST_Synplot2.txt")

# PROCESS Synplot2
synplot_data$P_value <- as.numeric(pnorm(runmean(synplot_data$subs_observed - synplot_data$subs_expected, win_size) / sqrt(runmean(synplot_data$subs_stdev * synplot_data$subs_stdev, win_size)) * sqrt(win_size)))
synplot_data$P_value_inverse <- as.numeric(1 / synplot_data$P_value)
synplot_data$selection_measure <- log(synplot_data$P_value_inverse)

# Add selection type
synplot_data$selection_type <- NA
synplot_data[synplot_data$subs_observed <= synplot_data$subs_expected & is.finite(synplot_data$subs_observed) & is.finite(synplot_data$subs_expected), ]$selection_type <- 'Negative'
synplot_data[synplot_data$subs_observed > synplot_data$subs_expected & is.finite(synplot_data$subs_observed) & is.finite(synplot_data$subs_expected), ]$selection_type <- 'Positive'

synplot_data$analysis <- 'Synplot2'
(synplot_staged <- dplyr::select(synplot_data, codon_num, analysis, selection_type, selection_measure))


############################################################################################################
############################################################################################################
### OLGENIE

### CHOOSE ONE (1) OLGENIE DATASET TO ANALYZE
# Cassan M group dataset
OLGenie_data <- read_tsv(file = "HIV1_env_CassanM_OLGenie_sas12.txt") # ENTER YOUR PATH HERE
# Cassan non-M groups dataset
OLGenie_data <- read_tsv(file = "HIV1_env_CassanNonM_OLGenie_sas12.txt") # ENTER YOUR PATH HERE
# BLAST whole env dataset
OLGenie_data <- read_tsv(file = "HIV1_env_BLAST_OLGenie_sas12.txt") # ENTER YOUR PATH HERE


############################################################################################################
### BOOTSTRAP SLIDING WINDOWS

# Prepare new columns
OLGenie_data$dNNdSN_sw <- NA
OLGenie_data$dNNdSN_sw_P <- NA
OLGenie_data$dNNdNS_sw <- NA
OLGenie_data$dNNdNS_sw_P <- NA
OLGenie_data$dNSdSS_sw <- NA
OLGenie_data$dNSdSS_sw_P <- NA
OLGenie_data$dSNdSS_sw <- NA
OLGenie_data$dSNdSS_sw_P <- NA

for(i in 1:(nrow(OLGenie_data) - win_size + 1)) { # each window starting at row 1
  
  # Extract window; analyze
  # Can add criterion here to keep only sites with a sufficient number of defined codons
  window_codon_data <- OLGenie_data[i:(i + win_size - 1), ]
  
  if(nrow(window_codon_data) > 1) {
    ### BOOTSTRAP EACH RATIO FOR THIS WINDOW. Columns:
    #[1] num_replicates
    #[2] dN
    #[3] dS
    #[4] dNdS
    #[5] dN_m_dS
    #[6] boot_dN_SE
    #[7] boot_dS_SE
    #[8] boot_dN_over_dS_SE
    #[9] boot_dN_over_dS_P
    #[10] boot_dN_m_dS_SE
    #[11] boot_dN_m_dS_P
    
    # dNN/dSN, ORF1
    boot_dNNdSN <- dNdS_diff_boot_fun(window_codon_data, 'NN', 'SN', NBOOTSTRAPS, NCPUS)
    boot_dNNdSN <- str_split(string = boot_dNNdSN, pattern = '\t')[[1]]
    dNNdSN <- as.numeric(boot_dNNdSN[4])
    dNNdSN_P <- as.numeric(boot_dNNdSN[11])
    
    # dNN/dNS, ORF2
    boot_dNNdNS <- dNdS_diff_boot_fun(window_codon_data, 'NN', 'NS', NBOOTSTRAPS, NCPUS)
    boot_dNNdNS <- str_split(string = boot_dNNdNS, pattern = '\t')[[1]]
    dNNdNS <- as.numeric(boot_dNNdNS[4])
    dNNdNS_P <- as.numeric(boot_dNNdNS[11])
    
    # dNS/dSS, ORF1
    boot_dNSdSS <- dNdS_diff_boot_fun(window_codon_data, 'NS', 'SS', NBOOTSTRAPS, NCPUS)
    boot_dNSdSS <- str_split(string = boot_dNSdSS, pattern = '\t')[[1]]
    dNSdSS <- as.numeric(boot_dNSdSS[4])
    dNSdSS_P <- as.numeric(boot_dNSdSS[11])
    
    # dSN/dSS, ORF2
    boot_dSNdSS <- dNdS_diff_boot_fun(window_codon_data, 'SN', 'SS', NBOOTSTRAPS, NCPUS)
    boot_dSNdSS <- str_split(string = boot_dSNdSS, pattern = '\t')[[1]]
    dSNdSS <- as.numeric(boot_dSNdSS[4])
    dSNdSS_P <- as.numeric(boot_dSNdSS[11])
    ######
    
    # Add to table
    OLGenie_data[i, ]$dNNdSN_sw <- dNNdSN
    OLGenie_data[i, ]$dNNdSN_sw_P <- dNNdSN_P
    OLGenie_data[i, ]$dNNdNS_sw <- dNNdNS
    OLGenie_data[i, ]$dNNdNS_sw_P <- dNNdNS_P
    OLGenie_data[i, ]$dNSdSS_sw <- dNSdSS
    OLGenie_data[i, ]$dNSdSS_sw_P <- dNSdSS_P
    OLGenie_data[i, ]$dSNdSS_sw <- dSNdSS
    OLGenie_data[i, ]$dSNdSS_sw_P <- dSNdSS_P
  } # else leave it as NA
} # end last window

# PROCESS OLGenie (choose which ratio you want to analyze by changing 'dSNdSS_sw_P')
OLGenie_data$dSNdSS_sw_P_inverse <- as.numeric(1 / OLGenie_data$dSNdSS_sw_P)
OLGenie_data$selection_measure <- log(OLGenie_data$dSNdSS_sw_P_inverse)
  
# Add selection type
OLGenie_data$selection_type <- NA
OLGenie_data[OLGenie_data$dSNdSS_sw <= 1 & is.finite(OLGenie_data$dSNdSS_sw), ]$selection_type <- 'Purifying'
OLGenie_data[OLGenie_data$dSNdSS_sw > 1 & OLGenie_data$dSNdSS_sw < 1.5 & is.finite(OLGenie_data$dSNdSS_sw), ]$selection_type <- 'Positive'

OLGenie_data$analysis <- 'OLGenie' # or Synplot, or name of ratio
OLGenie_staged <- dplyr::select(OLGenie_data, codon_num, analysis, selection_type, selection_measure)
  
### COMBINE SYNPLOT & OLGENIE
combined_data <- rbind(synplot_staged, OLGenie_staged)
combined_data$central_codon_num <- (2 * combined_data$codon_num + win_size) / 2

# Calculate approximate multiple comparison-corrected P-value using synplot method, i.e., 0.05 / (CDS length/window size)
P_cutoff <- log(1 / ((0.05) / (((length(unique(combined_data$codon_num)) + win_size) / win_size))))

