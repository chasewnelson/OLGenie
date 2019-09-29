#! /usr/local/bin/Rscript --slave --no-restore --no-save

############################################################################################################
# OLGenie process codon results to bootstrap

# Mac
suppressMessages(library(package = tidyverse))
suppressMessages(library(package = boot))

# At command line, call something like:
# Rscript --quiet --no-restore --no-save script_file.R arg1 arg2 arg3 > results.out
# Rscript --slave --quiet --no-restore --no-save OLGenie_process_codons_batch.R OLGenie_codon_results.txt 10000 > OLGenie_results.out

# Gather command-line arguments
argv <- commandArgs(trailingOnly = T)
codon_results_file <- as.character(argv[1])
NBOOTSTRAPS <- as.integer(argv[2]) # <- 10000
NCPUS <- as.integer(argv[3]) # <- 4, 8, or 60


############################################################################################################
# *BASIC* BOOTSTRAP FUNCTION (dN - dS) for CODON UNIT
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
  (boot_dN <- boot(data = codon_results, R = num_replicates, statistic = dN_function, parallel = 'multicore', ncpus = num_cpus))
  (dN <- boot_dN$t0)
  (boot_dN_SE <- sd(boot_dN$t))
  
  # boot dS
  (boot_dS <- boot(data = codon_results, R = num_replicates, statistic = dS_function, parallel = 'multicore', ncpus = num_cpus))
  (dS <- boot_dS$t0)
  (boot_dS_SE <- sd(boot_dS$t))
  
  # boot dN - dS
  (boot_dN_m_dS <- boot(data = codon_results, R = num_replicates, statistic = dN_m_dS_function, parallel = 'multicore', ncpus = num_cpus))
  (dN_m_dS <- boot_dN_m_dS$t0)
  (boot_dN_m_dS_SE <- sd(boot_dN_m_dS$t))
  (boot_dN_m_dS_Z <- dN_m_dS / boot_dN_m_dS_SE)
  (boot_dN_m_dS_P <- 2 * pnorm(-abs(boot_dN_m_dS_Z)))
  
  # boot dN/dS
  (boot_dN_over_dS <- boot(data = codon_results, R = num_replicates, statistic = dN_over_dS_function, parallel = 'multicore', ncpus = num_cpus))
  (dN_over_dS <- boot_dN_over_dS$t0)
  (boot_dN_over_dS_SE <- sd(boot_dN_over_dS$t))
  (boot_dN_over_dS_Z <- dN_over_dS / boot_dN_over_dS_SE)
  (boot_dN_over_dS_P <- 2 * pnorm(-abs(boot_dN_over_dS_Z)))
  
  return(paste(num_replicates, dN, dS, dNdS, dN_m_dS, boot_dN_SE, boot_dS_SE, boot_dN_over_dS_SE, boot_dN_over_dS_P, boot_dN_m_dS_SE, boot_dN_m_dS_P, sep = "\t"))
}

# INPUT
codon_results <- suppressMessages(read_tsv(codon_results_file))

# LEADING SUMMARY COLUMNS:
NN_sites <- sum(codon_results$NN_sites, na.rm = T)
SN_sites <- sum(codon_results$SN_sites, na.rm = T)
NS_sites <- sum(codon_results$NS_sites, na.rm = T)
SS_sites <- sum(codon_results$SS_sites, na.rm = T)
summary_data <- paste(nrow(codon_results),
                      NN_sites, SN_sites, NS_sites, SS_sites, 
                      sum(codon_results$NN_diffs, na.rm = T), sum(codon_results$SN_diffs, na.rm = T), 
                      sum(codon_results$NS_diffs, na.rm = T), sum(codon_results$SS_diffs, na.rm = T), 
                      sep = "\t")

# Determine which ORF1 ratio to use: dNN/dSN or dNS/dSS
dNNdSN_ratio_used <- (min(NN_sites, SN_sites) > min(NS_sites, SS_sites))

# Determine which ORF2 ratio to use: dNN/dNS or dSN/dSS
dNNdNS_ratio_used <- (min(NN_sites, NS_sites) > min(SN_sites, SS_sites))

# BOOTSTRAP EACH RATIO
boot_dNNdSN <- dNdS_diff_boot_fun(codon_results, 'NN', 'SN', NBOOTSTRAPS, NCPUS)
boot_dNNdNS <- dNdS_diff_boot_fun(codon_results, 'NN', 'NS', NBOOTSTRAPS, NCPUS)
boot_dNSdSS <- dNdS_diff_boot_fun(codon_results, 'NS', 'SS', NBOOTSTRAPS, NCPUS)
boot_dSNdSS <- dNdS_diff_boot_fun(codon_results, 'SN', 'SS', NBOOTSTRAPS, NCPUS)

# OUTPUT ALL RESULTS
cat(summary_data, 'dNNdSN', dNNdSN_ratio_used, 'ORF1', boot_dNNdSN, "\n", sep = "\t")
cat(summary_data, 'dNNdNS', dNNdNS_ratio_used, 'ORF2', boot_dNNdNS, "\n", sep = "\t")
cat(summary_data, 'dNSdSS', ! dNNdSN_ratio_used, 'ORF1', boot_dNSdSS, "\n", sep = "\t")
cat(summary_data, 'dSNdSS', ! dNNdNS_ratio_used, 'ORF2', boot_dSNdSS, "\n", sep = "\t")

