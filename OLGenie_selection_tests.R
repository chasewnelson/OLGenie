################################################################################################
# OLGenie TESTS FOR SELECTION
################################################################################################
# FISHER's exact test for selection, based on Zhang et al.:
#   Zhang J, Kumar S, Nei M. 1997. Small-sample tests of episodic adaptive evolution: a case 
#      study of primate lysozymes. Molecular Biology and Evolution 14:1335–1338.
#
# WILCOXON paired test for selection using (independent) sister pairs, based on Hughes et al.:
#   Hughes AL, Friedman R, Glenn NL. 2006. The Future of Data Analysis in Evolutionary Genomics. 
#      Current Genomics 7:227–234.
#
# BINOMIAL test for selection, based on plain good sense, but Type I error rate too high.
#


################################################################################################
# MANUALLY ENTER THE PATHS TO YOUR RESULTS FILES BELOW, AS IN:
#OLGenie_results <- read.delim("<PATH_TO>/OLGenie_results.txt") # <-- ENTER HERE, replace <PATH_TO>
#OLGenie_pair_results <- read.delim("<PATH_TO>/OLGenie_pair_results.txt") # <-- ENTER HERE, replace <PATH_TO>

######################################## <-----
######################################## <-----
### *ENTER YOUR OWN FILE PATHS HERE* ### <-----
OLGenie_results <- read.delim("/Users/cwnelson88/Desktop/OLG/ex_20190110/OLGenie_results.txt", header = TRUE)#OLGenie_results <- read.delim("/Users/cwnelson88/Desktop/OLG/simulations/sas12/dist_0.1_dnds1_0.5_dnds2_1_R_1/bootstrap_90/OLGenie_results.txt", header = TRUE)
OLGenie_pair_results <- read.delim("/Users/cwnelson88/Desktop/OLG/ex_20190110/OLGenie_pair_results.txt", header = TRUE)
######################################## <-----
######################################## <-----


################################################################################################
# SIMPLY RUN THE REMAINDER OF THE SCRIPT, COPY & PASTE OUTPUT

CORRECTED <- T # to correct is always best

# SITES
(NN_sites <- OLGenie_results[OLGenie_results$measure == 'sites', 'NN'])
(SN_sites <- OLGenie_results[OLGenie_results$measure == 'sites', 'SN'])
(NS_sites <- OLGenie_results[OLGenie_results$measure == 'sites', 'NS'])
(SS_sites <- OLGenie_results[OLGenie_results$measure == 'sites', 'SS'])

if(CORRECTED == TRUE) {
  # CORRECTED: Un-comment and run the following if you'd prefer to examine CORRECTED differences
  (NN_diffs <- OLGenie_results[OLGenie_results$measure == 'diffs_normalized', 'NN'])
  (SN_diffs <- OLGenie_results[OLGenie_results$measure == 'diffs_normalized', 'SN'])
  (NS_diffs <- OLGenie_results[OLGenie_results$measure == 'diffs_normalized', 'NS'])
  (SS_diffs <- OLGenie_results[OLGenie_results$measure == 'diffs_normalized', 'SS'])
} else {
  # RAW: Un-comment and run the following if you'd prefer to examine RAW differences
  (NN_diffs <- OLGenie_results[OLGenie_results$measure == 'diffs_raw', 'NN'])
  (SN_diffs <- OLGenie_results[OLGenie_results$measure == 'diffs_raw', 'SN'])
  (NS_diffs <- OLGenie_results[OLGenie_results$measure == 'diffs_raw', 'NS'])
  (SS_diffs <- OLGenie_results[OLGenie_results$measure == 'diffs_raw', 'SS'])
}

#results.table <- data.frame(dNNdSN = numeric(), dNSdSS = numeric(), dNNdNS = numeric(), dSNdSS = numeric(),
#                            dNNdSN_Fisher_P = numeric(), dNNdSN_Binomial_P = numeric(), dNNdSN_Wilcoxon_P = numeric(),
#                            dNSdSS_Fisher_P = numeric(), dNSdSS_Binomial_P = numeric(), dNSdSS_Wilcoxon_P = numeric(),
#                            dNNdNS_Fisher_P = numeric(), dNNdNS_Binomial_P = numeric(), dNNdNS_Wilcoxon_P = numeric(),
#                            dSNdSS_Fisher_P = numeric(), dSNdSS_Binomial_P = numeric(), dSNdSS_Wilcoxon_P = numeric())


################################################################################################
# dN/dS POINT ESTIMATES
(dNN <- (NN_diffs / NN_sites))
(dSN <- (SN_diffs / SN_sites))
(dNS <- (NS_diffs / NS_sites))
(dSS <- (SS_diffs / SS_sites))
(dNNdSN <- (NN_diffs / NN_sites) / (SN_diffs / SN_sites)) # dN/dS ORF1, estimate 1
(dNSdSS <- (NS_diffs / NS_sites) / (SS_diffs / SS_sites)) # dN/dS ORF1, estimate 2
(dNNdNS <- (NN_diffs / NN_sites) / (NS_diffs / NS_sites)) # dN/dS ORF2, estimate 1
(dSNdSS <- (SN_diffs / SN_sites) / (SS_diffs / SS_sites)) # dN/dS ORF2, estimate 2
#results.table <- rbind(results.table, data.frame(dNNdSN = dNNdSN, dNSdSS = dNSdSS, dNNdNS = dNNdNS, dSNdSS = dSNdSS))
#temp_row <- data.frame(dNNdSN = dNNdSN, dNSdSS = dNSdSS, dNNdNS = dNNdNS, dSNdSS = dSNdSS)


################################################################################################
# SIMPLY RUN THE FOLLOWING
NN_no_diffs <- NN_sites - NN_diffs
SN_no_diffs <- SN_sites - SN_diffs
NS_no_diffs <- NS_sites - NS_diffs
SS_no_diffs <- SS_sites - SS_diffs


################################################################################################
# ORF1 pNN/pSN
input_table <- c(NN_diffs, SN_diffs, NN_no_diffs, SN_no_diffs)
(input_table <- matrix(input_table, nrow = 2, ncol = 2, byrow = T, dimnames = list(c('DIFFS', 'NO_DIFFS'), c('NN', 'SN'))))

# Three possible tests
## (dNNdSN_Fisher_P <- fisher.test(input_table, alternative = "two.sided"))
## (dNNdSN_Binomial_P <- binom.test(x = round(NN_diffs), n = round(NN_diffs + SN_diffs), p = (NN_sites / (NN_sites + SN_sites))))
(dNNdSN_Wilcoxon_P <- wilcox.test(x = OLGenie_pair_results$pNN, y = OLGenie_pair_results$pSN, paired = T))
#t.test(x = OLGenie_pair_results$pNN, y = OLGenie_pair_results$pSN, paired = T)

################################################################################################
# ORF1 pNS/pSS
input_table <- c(NS_diffs, SS_diffs, NS_no_diffs, SS_no_diffs)
(input_table <- matrix(input_table, nrow = 2, ncol = 2, byrow = T, dimnames = list(c('DIFFS', 'NO_DIFFS'), c('NS', 'SS'))))

# Three possible tests
## (dNSdSS_Fisher_P <- fisher.test(input_table, alternative = "two.sided"))
## (dNSdSS_Binomial_P <- binom.test(x = round(NS_diffs), n = round(NS_diffs + SS_diffs), p = (NS_sites / (NS_sites + SS_sites))))
(dNSdSS_Wilcoxon_P <- wilcox.test(x = OLGenie_pair_results$pNS, y = OLGenie_pair_results$pSS, paired = T))
#t.test(x = OLGenie_pair_results$pNS, y = OLGenie_pair_results$pSS, paired = T)


################################################################################################
#ORF2 pNN/pNS
input_table <- c(NN_diffs, NS_diffs, NN_no_diffs, NS_no_diffs)
(input_table <- matrix(input_table, nrow = 2, ncol = 2, byrow = T, dimnames = list(c('DIFFS', 'NO_DIFFS'), c('NN', 'NS'))))

# Three possible tests
## (dNNdNS_Fisher_P <- fisher.test(input_table, alternative = "two.sided")) # P=0.2949, methinks it is better
## (dNNdNS_Binomial_P <- binom.test(x = round(NN_diffs), n = round(NN_diffs + NS_diffs), p = (NN_sites / (NN_sites + NS_sites))))
(dNNdNS_Wilcoxon_P <- wilcox.test(x = OLGenie_pair_results$pNN, y = OLGenie_pair_results$pNS, paired = T))
#t.test(x = OLGenie_pair_results$pNN, y = OLGenie_pair_results$pNS, paired = T)


################################################################################################
#ORF2 pSN/pSS
input_table <- c(SN_diffs, SS_diffs, SN_no_diffs, SS_no_diffs)
(input_table <- matrix(input_table, nrow = 2, ncol = 2, byrow = T, dimnames = list(c('DIFFS', 'NO_DIFFS'), c('SN', 'SS'))))

# Three possible tests
## (dSNdSS_Fisher_P <- fisher.test(input_table, alternative = "two.sided")) # P=0.2949, methinks it is better
## (dSNdSS_Binomial_P <- binom.test(x = round(SN_diffs), n = round(SN_diffs + SS_diffs), p = (SN_sites / (SN_sites + SS_sites))))
(dSNdSS_Wilcoxon_P <- wilcox.test(x = OLGenie_pair_results$pSN, y = OLGenie_pair_results$pSS, paired = T))
#t.test(x = OLGenie_pair_results$pSN, y = OLGenie_pair_results$pSS, paired = T)

# Number of sister pairs
(num_pairs <- nrow(OLGenie_pair_results))


################################################################################################
# PRINT HEADER AND DATA

# HEADER (un-comment if needed)
#cat('num_pairs',
#    'NN_sites', 'SN_sites', 'NS_sites', 'SS_sites',
#    'CORRECTED',
#    'NN_diffs', 'SN_diffs', 'NS_diffs', 'SS_diffs',
#    'NN_no_diffs', 'SN_no_diffs', 'NS_no_diffs', 'SS_no_diffs',
#    'dNN', 'dSN', 'dNS', 'dSS',
#    'dNNdSN', 'dNSdSS', 'dNNdNS', 'dSNdSS', 
#    'dNNdSN_Fisher_P$p.value', 'dNNdSN_Binomial_P$p.value', 'dNNdSN_Wilcoxon_P$p.value',
#    'dNSdSS_Fisher_P$p.value', 'dNSdSS_Binomial_P$p.value', 'dNSdSS_Wilcoxon_P$p.value',
#    'dNNdNS_Fisher_P$p.value', 'dNNdNS_Binomial_P$p.value', 'dNNdNS_Wilcoxon_P$p.value',
#    'dSNdSS_Fisher_P$p.value', 'dSNdSS_Binomial_P$p.value', 'dSNdSS_Wilcoxon_P$p.value',
#    sep = "\t")

## # DATA
## cat(num_pairs,
##     NN_sites, SN_sites, NS_sites, SS_sites,
##     CORRECTED,
##     NN_diffs, SN_diffs, NS_diffs, SS_diffs,
##     NN_no_diffs, SN_no_diffs, NS_no_diffs, SS_no_diffs,
##     dNN, dSN, dNS, dSS,
##     dNNdSN, dNSdSS, dNNdNS, dSNdSS, 
##     dNNdSN_Fisher_P$p.value, dNNdSN_Binomial_P$p.value, dNNdSN_Wilcoxon_P$p.value,
##     dNSdSS_Fisher_P$p.value, dNSdSS_Binomial_P$p.value, dNSdSS_Wilcoxon_P$p.value,
##     dNNdNS_Fisher_P$p.value, dNNdNS_Binomial_P$p.value, dNNdNS_Wilcoxon_P$p.value,
##     dSNdSS_Fisher_P$p.value, dSNdSS_Binomial_P$p.value, dSNdSS_Wilcoxon_P$p.value,
##     sep = "\t")

# DATA
cat(num_pairs,
    NN_sites, SN_sites, NS_sites, SS_sites,
    CORRECTED,
    NN_diffs, SN_diffs, NS_diffs, SS_diffs,
    NN_no_diffs, SN_no_diffs, NS_no_diffs, SS_no_diffs,
    dNN, dSN, dNS, dSS,
    dNNdSN, dNSdSS, dNNdNS, dSNdSS, 
    dNNdSN_Wilcoxon_P$p.value,
    dNSdSS_Wilcoxon_P$p.value,
    dNNdNS_Wilcoxon_P$p.value,
    dSNdSS_Wilcoxon_P$p.value,
    sep = "\t")


