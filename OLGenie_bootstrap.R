#! /usr/local/bin/Rscript --slave --no-restore --no-save

############################################################################################################
# OLGenie process codon results to bootstrap

# Mac
suppressMessages(library(package = readr))
suppressMessages(library(package = stringr))
suppressMessages(library(package = boot))

# At command line, call something like:
# Rscript --quiet --no-restore --no-save script_file.R arg1 arg2 arg3 > results.out
# Rscript --slave --quiet --no-restore --no-save OLGenie_process_codons_batch.R OLGenie_codon_results.txt 10000 > OLGenie_results.out

############################################################################################################
############################################################################################################
### GATHER GLOBAL VARIABLES FROM COMMAND LINE
ARGV <- commandArgs(trailingOnly = T)

kill_script <- FALSE

if(! (length(ARGV) >= 3)) {
  kill_script <- TRUE
} else if(! str_detect(string = ARGV[1], pattern = "\\w")) {
  kill_script <- TRUE
} else if(! str_detect(string = ARGV[2], pattern = "\\d")) {
  kill_script <- TRUE
} else if(! str_detect(string = ARGV[3], pattern = "\\d")) {
  kill_script <- TRUE
} # could add more conditions later

if(kill_script) {
  cat("\n\n### WARNING: there must be 3-6 command line arguments, in this order:\n")
  cat("    (1) CODON RESULTS FILE\n")
  cat("    (2) MINIMUM NUMBER OF DEFINED CODONS\n")
  cat("    (3) NUMBER OF BOOTSTRAP REPLICATES\n")
  cat("    (4) NUMBER OF CPUS (OPTIONAL)\n")
  cat("    (5) STRING TO PREPEND TO OUTPUT LINES (OPTIONAL)\n\n")
  quit(save = 'no', status = 1, runLast = TRUE)
} #else {
#  cat("\n\n### RUNNING! ###\n")
#}

CODON_RESULTS_FILE <- as.character(ARGV[1])
MIN_DEFINED_CODONS <- as.integer(ARGV[2]) # 0, 6, or 8 # only used by OLGenie
NBOOTSTRAPS <- as.integer(ARGV[3]) # <- 10000

NCPUS <- 1
if(! is.na(ARGV[4]) && str_detect(string = ARGV[4], pattern = "\\d") && ! str_detect(string = ARGV[4], pattern = "[a-zA-Z]")) {
  NCPUS <- as.integer(ARGV[4]) # <- 4, 8, or 60
}

PREPEND_TO_OUTPUT <- ''
if(! is.na(ARGV[5])) {
  PREPEND_TO_OUTPUT <- as.character(ARGV[5])
}

# EXAMPLES
#CODON_RESULTS_FILE <- "/Users/cwnelson88/Desktop/OLG/simulations/SIMULATIONS/OLGenie_codon_results.txt"
#MIN_DEFINED_CODONS <- 6
#NBOOTSTRAPS <- 100
#NCPUS <- 4
#PREPEND_TO_OUTPUT <- 'METADATA'

# Read in the file
suppressMessages(OLGenie_data <- read_tsv(file = CODON_RESULTS_FILE))


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
codon_results <- suppressMessages(read_tsv(CODON_RESULTS_FILE)) # 1,353 x 16

# FILTER BY MIN NUMBER OF DEFINED CODONS
codon_results <- codon_results[codon_results$num_defined_seqs >= MIN_DEFINED_CODONS, ] # 982 x 16

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

if(PREPEND_TO_OUTPUT != '') {
  summary_data <- paste(PREPEND_TO_OUTPUT, summary_data, sep = "\t")
}

# Determine which ORF1 ratio to use: dNN/dSN or dNS/dSS
dNNdSN_ratio_used <- (min(NN_sites, SN_sites) > min(NS_sites, SS_sites))

# Determine which ORF2 ratio to use: dNN/dNS or dSN/dSS
dNNdNS_ratio_used <- (min(NN_sites, NS_sites) > min(SN_sites, SS_sites))

# BOOTSTRAP EACH RATIO
boot_dNNdSN <- dNdS_diff_boot_fun(codon_results, 'NN', 'SN', NBOOTSTRAPS, NCPUS)
boot_dNNdNS <- dNdS_diff_boot_fun(codon_results, 'NN', 'NS', NBOOTSTRAPS, NCPUS)
boot_dNSdSS <- dNdS_diff_boot_fun(codon_results, 'NS', 'SS', NBOOTSTRAPS, NCPUS)
boot_dSNdSS <- dNdS_diff_boot_fun(codon_results, 'SN', 'SS', NBOOTSTRAPS, NCPUS)

# OUTPUT HEADER
header <- c('num_codons', 'NN_sites', 'SN_sites', 'NS_sites', 'SS_sites', 'NN_diffs', 'SN_diffs', 'NS_diffs', 'SS_diffs', 
            'ratio', 'site_rich_ratio', 'gene', 'num_replicates', 
            'dN', 'dS', 'dNdS', 'dN_m_dS', 'boot_dN_SE', 'boot_dS_SE', 'boot_dN_over_dS_SE', 'boot_dN_over_dS_P', 'boot_dN_m_dS_SE', 'P_value')

if(PREPEND_TO_OUTPUT != '') {
  header <- c('METADATA', header)
}

cat(header, sep = "\t")
cat("\n", sep = '')

# OUTPUT ALL RESULTS
cat(summary_data, 'dNNdSN', dNNdSN_ratio_used, 'ORF1', boot_dNNdSN, sep = "\t")
cat("\n", sep = '')
cat(summary_data, 'dNNdNS', dNNdNS_ratio_used, 'ORF2', boot_dNNdNS, sep = "\t")
cat("\n", sep = '')
cat(summary_data, 'dNSdSS', ! dNNdSN_ratio_used, 'ORF1', boot_dNSdSS, sep = "\t")
cat("\n", sep = '')
cat(summary_data, 'dSNdSS', ! dNNdNS_ratio_used, 'ORF2', boot_dSNdSS, sep = "\t")
cat("\n", sep = '')

