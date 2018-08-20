overlapgenie.analyze.comparisons <- function(overlapgenie.results) {
  # Execute script to load function overlapgenie.test.means
  # Call function with two arguments: the <OLG>_results.txt file from overlapgenie and a number of bootstraps
  # EXAMPLE CALL: overlapgenie.test.means(my_OLG_results,5000)
  
  # EXAMPLE INPUT:
  #  overlapgenie.results <- read.delim("/Users/cwnelson88/Desktop/my_OLG_results.txt") # nrow(overlapgenie.results)
  
#  dNN.dSN.subset <- overlapgenie.results[overlapgenie.results$dNN!='*' & overlapgenie.results$dSN!='*',] # this works
  
  # Initialized vectors for storage
  dNN.dSN.Zs <- vector()
  dNN.dSN.Ps <- vector()
  
  dNS.dSS.Zs <- vector()
  dNS.dSS.Ps <- vector()
  
  dNN.dNS.Zs <- vector()
  dNN.dNS.Ps <- vector()
  
  dSN.dSS.Zs <- vector()
  dSN.dSS.Ps <- vector()
  
  for(i in 1:nrow(overlapgenie.results)) { # for each bootstrap replicate
    # REF measure 1: dNN/dSN
    if(overlapgenie.results$dNN[i]!='*' & overlapgenie.results$dSN[i]!='*') {
      this.dNN <- as.double(as.character(overlapgenie.results$dNN[i]))
      this.dSN <- as.double(as.character(overlapgenie.results$dSN[i]))
      
      this.dNN.var <- as.double(as.character(overlapgenie.results$var.dNN.[i]))
      this.dSN.var <- as.double(as.character(overlapgenie.results$var.dSN.[i]))
      
      dNN.dSN.Z <- (this.dNN - this.dSN) / sqrt(this.dNN.var + this.dSN.var)
      dNN.dSN.Zs[i] <- dNN.dSN.Z
      
      dNN.dSN.P <- 2*pnorm(-abs(dNN.dSN.Z))
      dNN.dSN.Ps[i] <- dNN.dSN.P
      
    } # else they'll be NA or NaN
    
    # REF measure 2: dNS/dSS
    if(overlapgenie.results$dNS[i]!='*' & overlapgenie.results$dSS[i]!='*') {
      this.dNS <- as.double(as.character(overlapgenie.results$dNS[i]))
      this.dSS <- as.double(as.character(overlapgenie.results$dSS[i]))
      
      this.dNS.var <- as.double(as.character(overlapgenie.results$var.dNS.[i]))
      this.dSS.var <- as.double(as.character(overlapgenie.results$var.dSS.[i]))
      
      dNS.dSS.Z <- (this.dNS - this.dSS) / sqrt(this.dNS.var + this.dSS.var)
      dNS.dSS.Zs[i] <- dNS.dSS.Z
      
      dNS.dSS.P <- 2*pnorm(-abs(dNS.dSS.Z))
      dNS.dSS.Ps[i] <- dNS.dSS.P
      
    } # else they'll be NA or NaN
    
    # ALT measure 1: dNN/dNS
    if(overlapgenie.results$dNN[i]!='*' & overlapgenie.results$dNS[i]!='*') {
      this.dNN <- as.double(as.character(overlapgenie.results$dNN[i]))
      this.dNS <- as.double(as.character(overlapgenie.results$dNS[i]))
      
      this.dNN.var <- as.double(as.character(overlapgenie.results$var.dNN.[i]))
      this.dNS.var <- as.double(as.character(overlapgenie.results$var.dNS.[i]))
      
      dNN.dNS.Z <- (this.dNN - this.dNS) / sqrt(this.dNN.var + this.dNS.var)
      dNN.dNS.Zs[i] <- dNN.dNS.Z
      
      dNN.dNS.P <- 2*pnorm(-abs(dNN.dNS.Z))
      dNN.dNS.Ps[i] <- dNN.dNS.P
      
    } # else they'll be NA or NaN
    
    # ALT measure 1: dSN/dSS
    if(overlapgenie.results$dSN[i]!='*' & overlapgenie.results$dSS[i]!='*') {
      this.dSN <- as.double(as.character(overlapgenie.results$dSN[i]))
      this.dSS <- as.double(as.character(overlapgenie.results$dSS[i]))
      
      this.dSN.var <- as.double(as.character(overlapgenie.results$var.dSN.[i]))
      this.dSS.var <- as.double(as.character(overlapgenie.results$var.dSS.[i]))
      
      dSN.dSS.Z <- (this.dSN - this.dSS) / sqrt(this.dSN.var + this.dSS.var)
      dSN.dSS.Zs[i] <- dSN.dSS.Z
      
      dSN.dSS.P <- 2*pnorm(-abs(dSN.dSS.Z))
      dSN.dSS.Ps[i] <- dSN.dSS.P
      
    } # else they'll be NA or NaN
    
  }
  
  
  dNN.dSN.Ps.corrected <- p.adjust(dNN.dSN.Ps, method = "BY")
  dNS.dSS.Ps.corrected <- p.adjust(dNS.dSS.Ps, method = "BY")
  dNN.dNS.Ps.corrected <- p.adjust(dNN.dNS.Ps, method = "BY")
  dSN.dSS.Ps.corrected <- p.adjust(dSN.dSS.Ps, method = "BY")
  
  ##############################
  # Calculate SEMs, Z statistics, and P values
  
  # Ref gene dN/dS estimate 1: dNN.dSN.statistic, dNN.dSN.slope
#  dNN.dSN.slope.SEM <- sd(dNN.dSN.slopes, na.rm=T)
#  dNN.dSN.Z <- dNN.dSN.slope / dNN.dSN.slope.SEM
#  dNN.dSN.P <- 2*pnorm(-abs(dNN.dSN.Z))
#  cat("\n","REFERENCE dN/dS estimate 1:","\n","dNN/dSN=",dNN.dSN.statistic,"\n",
#      "SLOPE dNN/dSN=",dNN.dSN.slope,"\n",
#      "SE=",dNN.dSN.slope.SEM,"\n",
#      "Z=",dNN.dSN.Z,"\n",
#      "P=",dNN.dSN.P,"\n",
#      sep="")
  
}
