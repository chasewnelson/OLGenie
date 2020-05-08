#! /usr/bin/env Rscript --slave --no-restore --no-save

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
  cat("    (2) MINIMUM NUMBER OF DEFINED CODONS PER CODON POSITION\n")
  cat("    (3) NUMBER OF BOOTSTRAP REPLICATES\n")
  cat("    (4) NUMBER OF CPUS (OPTIONAL; ≥1; DEFAULT=1)\n")
  cat("    (5) MULTIPLE HITS CORRECTION (OPTIONAL; \"NONE\" or \"JC\"; DEFAULT=NONE)\n")
  cat("    (6) STRING TO PREPEND TO OUTPUT LINES (OPTIONAL; DEFAULT=\"\")\n\n")
  quit(save = 'no', status = 1, runLast = TRUE)
} #else {
#  cat("\n\n### RUNNING! ###\n")
#}

CODON_RESULTS_FILE <- as.character(ARGV[1])
MIN_DEFINED_CODONS <- as.integer(ARGV[2]) # 0, 6, or 8 # only used by OLGenie
NBOOTSTRAPS <- as.integer(ARGV[3]) # <- 10000
NCPUS <- 1
CORRECTION <- 'NONE' # "JC"
PREPEND_TO_OUTPUT <- ''

# Produce some helpful warning messages
if(! (NBOOTSTRAPS >= 1 && NBOOTSTRAPS <= 1000000)) {
  cat("### WARNING: NBOOTSTRAPS must be in the range [1,1000000]. Using: 1000.\n")
  NBOOTSTRAPS <- 10000
}

if(! (MIN_DEFINED_CODONS >= 2)) {
  cat("### WARNING: MIN_DEFINED_CODONS must be ≥2. Using: 6.\n")
  MIN_DEFINED_CODONS <- 6
}

if(! is.na(ARGV[4]) && str_detect(string = ARGV[4], pattern = "\\d") && ! str_detect(string = ARGV[4], pattern = "[a-zA-Z]")) {
  NCPUS <- as.integer(ARGV[4]) # <- 4, 8, or 60
}

if(! is.na(ARGV[5]) && ARGV[5] == "JC") {
  CORRECTION <- as.character(ARGV[5])
} #else {
#  cat("### WARNING: unrecognized CORRECTION supplied. Using: \"NONE\".\n")
#}

if(! is.na(ARGV[6])) {
  PREPEND_TO_OUTPUT <- as.character(ARGV[6])
}

# EXAMPLES
#CODON_RESULTS_FILE <- "XXX"
#MIN_DEFINED_CODONS <- 6
#NBOOTSTRAPS <- 100
#NCPUS <- 4
#PREPEND_TO_OUTPUT <- 'METADATA'

# Read in the file
suppressMessages(OLGenie_data <- read_tsv(file = CODON_RESULTS_FILE))


############################################################################################################
############################################################################################################
### BOOTSTRAP FUNCTIONS
############################################################################################################
############################################################################################################

############################################################################################################
# BOOTSTRAP FUNCTION (dN - dS) for CODON UNIT, NO CORRECTION
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
  
  ### NEW: ASL (acheived significance level)
  boot_dN_gt_dS_count <- sum(boot_dN_m_dS$t > 0)
  boot_dN_eq_dS_count <- sum(boot_dN_m_dS$t == 0)
  boot_dN_lt_dS_count <- sum(boot_dN_m_dS$t < 0)
  ASL_dN_gt_dS_P <- boot_dN_lt_dS_count / (boot_dN_gt_dS_count + boot_dN_eq_dS_count + boot_dN_lt_dS_count)
  ASL_dN_lt_dS_P <- boot_dN_gt_dS_count / (boot_dN_gt_dS_count + boot_dN_eq_dS_count + boot_dN_lt_dS_count)
  
  return(paste(num_replicates, dN, dS, dNdS, dN_m_dS, boot_dN_SE, boot_dS_SE, boot_dN_over_dS_SE, boot_dN_over_dS_P, 
               boot_dN_m_dS_SE, boot_dN_m_dS_P, 
               boot_dN_gt_dS_count, boot_dN_eq_dS_count, boot_dN_lt_dS_count, ASL_dN_gt_dS_P, ASL_dN_lt_dS_P,
               sep = "\t"))
}


############################################################################################################
# BOOTSTRAP FUNCTION (dN - dS) for CODON UNIT, ***JUKES-CANTOR CORRECTION***
dNdS_diff_boot_fun_JC <- function(codon_results, numerator, denominator, num_replicates, num_cpus) {
  
  # Function for dN
  dN_function <- function(D, indices) {
    dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
    dN <- -3/4 * log(1 - (4/3) * dN)
    return(dN)
  }
  
  # Function for dN
  dS_function <- function(D, indices) {
    dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
    dS <- -3/4 * log(1 - (4/3) * dS)
    return(dS)
  }
  
  # Function for dN - dS
  dN_m_dS_function <- function(D, indices) {
    dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
    dN <- -3/4 * log(1 - (4/3) * dN)
    dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
    dS <- -3/4 * log(1 - (4/3) * dS)
    dN_m_dS <- dN - dS
    return(dN_m_dS)
  }
  
  # Function for dN/dS
  dN_over_dS_function <- function(D, indices) {
    dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
    dN <- -3/4 * log(1 - (4/3) * dN)
    dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
    dS <- -3/4 * log(1 - (4/3) * dS)
    dN_over_dS <- dN / dS
    return(dN_over_dS)
  }
  
  # CREATE FUNCTION FOR dN/dS TO CALCULATE ITS SE
  dN <- sum(as.vector(codon_results[ , paste0(numerator, "_diffs")]), na.rm = TRUE) / sum(as.vector(codon_results[ , paste0(numerator, "_sites")]), na.rm = TRUE)
  dN <- -3/4 * log(1 - (4/3) * dN)
  dS <- sum(as.vector(codon_results[ , paste0(denominator, "_diffs")]), na.rm = TRUE) / sum(as.vector(codon_results[ , paste0(denominator, "_sites")]), na.rm = TRUE)
  dS <- -3/4 * log(1 - (4/3) * dS)
  dNdS <- dN / dS
  
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
  
  ### NEW: ASL (acheived significance level)
  boot_dN_gt_dS_count <- sum(boot_dN_m_dS$t > 0)
  boot_dN_eq_dS_count <- sum(boot_dN_m_dS$t == 0)
  boot_dN_lt_dS_count <- sum(boot_dN_m_dS$t < 0)
  ASL_dN_gt_dS_P <- boot_dN_lt_dS_count / (boot_dN_gt_dS_count + boot_dN_eq_dS_count + boot_dN_lt_dS_count)
  ASL_dN_lt_dS_P <- boot_dN_gt_dS_count / (boot_dN_gt_dS_count + boot_dN_eq_dS_count + boot_dN_lt_dS_count)
  
  return(paste(num_replicates, dN, dS, dNdS, dN_m_dS, boot_dN_SE, boot_dS_SE, boot_dN_over_dS_SE, boot_dN_over_dS_P, 
               boot_dN_m_dS_SE, boot_dN_m_dS_P, 
               boot_dN_gt_dS_count, boot_dN_eq_dS_count, boot_dN_lt_dS_count, ASL_dN_gt_dS_P, ASL_dN_lt_dS_P,
               sep = "\t"))
}

############################################################################################################
############################################################################################################
# INPUT
codon_results <- suppressMessages(read_tsv(CODON_RESULTS_FILE))

# FILTER BY MIN NUMBER OF DEFINED CODONS
codon_results <- codon_results[codon_results$num_defined_seqs >= MIN_DEFINED_CODONS, ]

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
boot_dNNdSN <- NA
boot_dNNdNS <- NA
boot_dNSdSS <- NA
boot_dSNdSS <- NA

# BOOTSTRAP: run the right kind of calculation for desired correction
if(CORRECTION == "JC") {
  boot_dNNdSN <- dNdS_diff_boot_fun_JC(codon_results, 'NN', 'SN', NBOOTSTRAPS, NCPUS)
  boot_dNNdNS <- dNdS_diff_boot_fun_JC(codon_results, 'NN', 'NS', NBOOTSTRAPS, NCPUS)
  boot_dNSdSS <- dNdS_diff_boot_fun_JC(codon_results, 'NS', 'SS', NBOOTSTRAPS, NCPUS)
  boot_dSNdSS <- dNdS_diff_boot_fun_JC(codon_results, 'SN', 'SS', NBOOTSTRAPS, NCPUS)
} else {
  boot_dNNdSN <- dNdS_diff_boot_fun(codon_results, 'NN', 'SN', NBOOTSTRAPS, NCPUS)
  boot_dNNdNS <- dNdS_diff_boot_fun(codon_results, 'NN', 'NS', NBOOTSTRAPS, NCPUS)
  boot_dNSdSS <- dNdS_diff_boot_fun(codon_results, 'NS', 'SS', NBOOTSTRAPS, NCPUS)
  boot_dSNdSS <- dNdS_diff_boot_fun(codon_results, 'SN', 'SS', NBOOTSTRAPS, NCPUS)
}

# OUTPUT HEADER
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
#[12] boot_dN_gt_dS_count
#[13] boot_dN_eq_dS_count
#[14] boot_dN_lt_dS_count
#[15] ASL_dN_gt_dS_P
#[16] ASL_dN_lt_dS_P
header <- c('num_codons', 'NN_sites', 'SN_sites', 'NS_sites', 'SS_sites', 'NN_diffs', 'SN_diffs', 'NS_diffs', 'SS_diffs', 
            'ratio', 'site_rich_ratio', 'gene', 'num_replicates', 
            'dN', 'dS', 'dNdS', 'dN_m_dS', 'boot_dN_SE', 'boot_dS_SE', 'boot_dN_over_dS_SE', 'boot_dN_over_dS_P', 'boot_dN_m_dS_SE', 'boot_dN_m_dS_P',
            'boot_dN_gt_dS_count', 'boot_dN_eq_dS_count', 'boot_dN_lt_dS_count', 'ASL_dN_gt_dS_P', 'ASL_dN_lt_dS_P')

if(PREPEND_TO_OUTPUT != '') {
  header <- c('METADATA', header, 'ASL_dNdS_P')
} else {
  header <- c(header, 'ASL_dNdS_P')
}

cat(header, sep = "\t")
cat("\n", sep = '')

### Calculate 2-sided ASL P-values
# dNNdSN
boot_dNNdSN_array <- str_split(string = boot_dNNdSN, pattern = '\t')[[1]]
suppressWarnings(boot_dNNdSN_ASL_dN_gt_dS_P <- as.numeric(boot_dNNdSN_array[15]))
suppressWarnings(boot_dNNdSN_ASL_dN_lt_dS_P <- as.numeric(boot_dNNdSN_array[16]))
boot_dNNdSN_ASL_dNdS_P <- 1
if(! is.na(boot_dNNdSN_ASL_dN_gt_dS_P) && ! is.na(boot_dNNdSN_ASL_dN_lt_dS_P) && boot_dNNdSN_ASL_dN_gt_dS_P < boot_dNNdSN_ASL_dN_lt_dS_P) {
  boot_dNNdSN_ASL_dNdS_P <- 2 * boot_dNNdSN_ASL_dN_gt_dS_P
} else if(! is.na(boot_dNNdSN_ASL_dN_gt_dS_P) && ! is.na(boot_dNNdSN_ASL_dN_lt_dS_P)) {
  boot_dNNdSN_ASL_dNdS_P <- 2 * boot_dNNdSN_ASL_dN_lt_dS_P
}
if(boot_dNNdSN_ASL_dNdS_P == 0) {
  boot_dNNdSN_ASL_dNdS_P <- 1 / NBOOTSTRAPS
}

# dNNdNS
boot_dNNdNS_array <- str_split(string = boot_dNNdNS, pattern = '\t')[[1]]
suppressWarnings(boot_dNNdNS_ASL_dN_gt_dS_P <- as.numeric(boot_dNNdNS_array[15]))
suppressWarnings(boot_dNNdNS_ASL_dN_lt_dS_P <- as.numeric(boot_dNNdNS_array[16]))
boot_dNNdNS_ASL_dNdS_P <- 1
if(! is.na(boot_dNNdNS_ASL_dN_gt_dS_P) && ! is.na(boot_dNNdNS_ASL_dN_lt_dS_P) && boot_dNNdNS_ASL_dN_gt_dS_P < boot_dNNdNS_ASL_dN_lt_dS_P) {
  boot_dNNdNS_ASL_dNdS_P <- 2 * boot_dNNdNS_ASL_dN_gt_dS_P
} else if(! is.na(boot_dNNdNS_ASL_dN_gt_dS_P) && ! is.na(boot_dNNdNS_ASL_dN_lt_dS_P)) {
  boot_dNNdNS_ASL_dNdS_P <- 2 * boot_dNNdNS_ASL_dN_lt_dS_P
}
if(boot_dNNdNS_ASL_dNdS_P == 0) {
  boot_dNNdNS_ASL_dNdS_P <- 1 / NBOOTSTRAPS
}

# dNSdSS
boot_dNSdSS_array <- str_split(string = boot_dNSdSS, pattern = '\t')[[1]]
suppressWarnings(boot_dNSdSS_ASL_dN_gt_dS_P <- as.numeric(boot_dNSdSS_array[15]))
suppressWarnings(boot_dNSdSS_ASL_dN_lt_dS_P <- as.numeric(boot_dNSdSS_array[16]))
boot_dNSdSS_ASL_dNdS_P <- 1
if(! is.na(boot_dNSdSS_ASL_dN_gt_dS_P) && ! is.na(boot_dNSdSS_ASL_dN_lt_dS_P) && boot_dNSdSS_ASL_dN_gt_dS_P < boot_dNSdSS_ASL_dN_lt_dS_P) {
  boot_dNSdSS_ASL_dNdS_P <- 2 * boot_dNSdSS_ASL_dN_gt_dS_P
} else if(! is.na(boot_dNSdSS_ASL_dN_gt_dS_P) && ! is.na(boot_dNSdSS_ASL_dN_lt_dS_P)) {
  boot_dNSdSS_ASL_dNdS_P <- 2 * boot_dNSdSS_ASL_dN_lt_dS_P
}
if(boot_dNSdSS_ASL_dNdS_P == 0) {
  boot_dNSdSS_ASL_dNdS_P <- 1 / NBOOTSTRAPS
}

# dSNdSS
boot_dSNdSS_array <- str_split(string = boot_dSNdSS, pattern = '\t')[[1]]
suppressWarnings(boot_dSNdSS_ASL_dN_gt_dS_P <- as.numeric(boot_dSNdSS_array[15]))
suppressWarnings(boot_dSNdSS_ASL_dN_lt_dS_P <- as.numeric(boot_dSNdSS_array[16]))
boot_dSNdSS_ASL_dNdS_P <- 1
if(! is.na(boot_dSNdSS_ASL_dN_gt_dS_P) && ! is.na(boot_dSNdSS_ASL_dN_lt_dS_P) && boot_dSNdSS_ASL_dN_gt_dS_P < boot_dSNdSS_ASL_dN_lt_dS_P) {
  boot_dSNdSS_ASL_dNdS_P <- 2 * boot_dSNdSS_ASL_dN_gt_dS_P
} else if(! is.na(boot_dSNdSS_ASL_dN_gt_dS_P) && ! is.na(boot_dSNdSS_ASL_dN_lt_dS_P)) {
  boot_dSNdSS_ASL_dNdS_P <- 2 * boot_dSNdSS_ASL_dN_lt_dS_P
}
if(boot_dSNdSS_ASL_dNdS_P == 0) {
  boot_dSNdSS_ASL_dNdS_P <- 1 / NBOOTSTRAPS
}


# OUTPUT ALL RESULTS
cat(summary_data, 'dNNdSN', dNNdSN_ratio_used, 'ORF1', boot_dNNdSN, boot_dNNdSN_ASL_dNdS_P, sep = "\t")
cat("\n", sep = '')
cat(summary_data, 'dNNdNS', dNNdNS_ratio_used, 'ORF2', boot_dNNdNS, boot_dNNdNS_ASL_dNdS_P, sep = "\t")
cat("\n", sep = '')
cat(summary_data, 'dNSdSS', ! dNNdSN_ratio_used, 'ORF1', boot_dNSdSS, boot_dNSdSS_ASL_dNdS_P, sep = "\t")
cat("\n", sep = '')
cat(summary_data, 'dSNdSS', ! dNNdNS_ratio_used, 'ORF2', boot_dSNdSS, boot_dSNdSS_ASL_dNdS_P, sep = "\t")
cat("\n", sep = '')

