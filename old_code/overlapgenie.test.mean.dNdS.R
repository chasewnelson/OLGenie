overlapgenie.test.means <- function(overlapgenie.results, num.bootstraps) {
  # Execute script to load function overlapgenie.test.means
  # Call function with two arguments: the <OLG>_results.txt file from overlapgenie and a number of bootstraps
  # EXAMPLE CALL: overlapgenie.test.means(my_OLG_results,5000)
  
  # EXAMPLE INPUT:
  #  overlapgenie.results <- read.delim("/Users/cwnelson88/Desktop/my_OLG_results.txt") # nrow(overlapgenie.results)
  
  # Find UNDEFINED values
#  dNN.undefined <- overlapgenie.results$dNN == '*'
#  dSN.undefined <- overlapgenie.results$dSN == '*'
#  dNS.undefined <- overlapgenie.results$dNS == '*'
#  dSS.undefined <- overlapgenie.results$dSS == '*'
  
#  dNN.results <- as.numeric(overlapgenie.results$dNN)
#  dNS.results <- as.numeric(overlapgenie.results$dNS)
  
  
#  dNN.dSN.subset <- subset(overlapgenie.results, !(dNN.undefined || dSN.undefined)) # nrow(dNN.dSN.subset)
#  dNN.dSN.subset <- overlapgenie.results[!(dNN.undefined || dSN.undefined),] # nrow(dNN.dSN.subset)
#  dNN.dSN.subset <- overlapgenie.results[overlapgenie.results$dNN!='*' & overlapgenie.results$dSN!='*',] # this works
  
  ##############################
  # Calculate mean dNN, dSN, dNS, and dSS for all pairwise comparisons
  dNN.mean <- mean(as.double(as.character(overlapgenie.results[overlapgenie.results$dNN != '*',]$dNN)), na.rm=T)
  dSN.mean <- mean(as.double(as.character(overlapgenie.results[overlapgenie.results$dSN != '*',]$dSN)), na.rm=T)
  dNS.mean <- mean(as.double(as.character(overlapgenie.results[overlapgenie.results$dNS != '*',]$dNS)), na.rm=T)
  dSS.mean <- mean(as.double(as.character(overlapgenie.results[overlapgenie.results$dSS != '*',]$dSS)), na.rm=T)
  
#  dNS.dSS.subset <- subset(overlapgenie.results, !(dNS.undefined || dSS.undefined)) # nrow(dNN.dSN.subset)
  
  ##############################
  # Ref gene dN/dS estimate 1
  dNN.dSN.statistic <- dNN.mean/dSN.mean
  
  # Estimate using slope
  dNN.dSN.subset <- overlapgenie.results[overlapgenie.results$dNN!='*' & overlapgenie.results$dSN!='*',] # this works
  dNN.dSN.lm <- lm(as.double(as.character(dNN.dSN.subset$dNN)) ~ as.double(as.character(dNN.dSN.subset$dSN)) - 1) # omitting intercept
  dNN.dSN.slope <- summary(dNN.dSN.lm)$coefficients[1]
  plot(as.double(as.character(dNN.dSN.subset$dNN)) ~ as.double(as.character(dNN.dSN.subset$dSN)), xlab="dSN", ylab="dNN")
  abline(0, 1)
  abline(dNN.dSN.lm, col="red") # regression line (y~x) 
  
  # Ref gene dN/dS estimate 2
  dNS.dSS.statistic <- dNS.mean/dSS.mean
  
  # Estimate using slope
  dNS.dSS.subset <- overlapgenie.results[overlapgenie.results$dNS!='*' & overlapgenie.results$dSS!='*',] # this works
  dNS.dSS.lm <- lm(as.double(as.character(dNS.dSS.subset$dNS)) ~ as.double(as.character(dNS.dSS.subset$dSS)) - 1) # omitting intercept
  dNS.dSS.slope <- summary(dNS.dSS.lm)$coefficients[1]
  
  ##############################
  # Alt gene dN/dS estimate 1
  dNN.dNS.statistic <- dNN.mean/dNS.mean
  
  # Estimate using slope
  dNN.dNS.subset <- overlapgenie.results[overlapgenie.results$dNN!='*' & overlapgenie.results$dNS!='*',] # this works
  dNN.dNS.lm <- lm(as.double(as.character(dNN.dNS.subset$dNN)) ~ as.double(as.character(dNN.dNS.subset$dNS)) - 1) # omitting intercept
  dNN.dNS.slope <- summary(dNN.dNS.lm)$coefficients[1]
  plot(as.double(as.character(dNN.dNS.subset$dNN)) ~ as.double(as.character(dNN.dNS.subset$dNS)), xlab="dNS", ylab="dNN")
  abline(0, 1)
  abline(dNN.dNS.lm, col="red") # regression line (y~x) 
  
  # Alt gene dN/dS estimate 2
  dSN.dSS.statistic <- dSN.mean/dSS.mean
  
  # Estimate using slope
  dSN.dSS.subset <- overlapgenie.results[overlapgenie.results$dSN!='*' & overlapgenie.results$dSS!='*',] # this works
  dSN.dSS.lm <- lm(as.double(as.character(dSN.dSS.subset$dSN)) ~ as.double(as.character(dSN.dSS.subset$dSS)) - 1) # omitting intercept
  dSN.dSS.slope <- summary(dSN.dSS.lm)$coefficients[1]
  
  ##############################
  # PERFORM BOOTSTRAP ANALYSIS
  # Ref gene dN/dS bootstraps
  dNN.dSN.slopes <- vector()
  dNS.dSS.slopes <- vector()
  
  # Alt gene dN/dS bootstraps
  dNN.dNS.slopes <- vector()
  dSN.dSS.slopes <- vector()
  
  for(i in 1:num.bootstraps) { # for each bootstrap replicate
    sample.size <- nrow(overlapgenie.results) # number of rows is sample size
    
    # Store THIS sample of points for which to calculate the slope
    dNN.rand.sample <- vector()
    dSN.rand.sample <- vector()
    dNS.rand.sample <- vector()
    dSS.rand.sample <- vector()
    
    for(j in 1:sample.size) { # create a randomization sample
      # Choose a random row (here, a single dNN/dSN/dNS/dSS observation)
      random.row <- sample(1:sample.size, 1)
      
      # Get values of row
      sampled.dNN <- overlapgenie.results$dNN[random.row]
      sampled.dSN <- overlapgenie.results$dSN[random.row]
      sampled.dNS <- overlapgenie.results$dNS[random.row]
      sampled.dSS <- overlapgenie.results$dSS[random.row]
      
      # Add to bootstrap sample
      if(j == 1) { # first sample
        dNN.rand.sample[j] <- sampled.dNN
        dSN.rand.sample[j] <- sampled.dSN
        dNS.rand.sample[j] <- sampled.dNS
        dSS.rand.sample[j] <- sampled.dSS
      } else {
        dNN.rand.sample <- append(dNN.rand.sample, sampled.dNN)
        dSN.rand.sample <- append(dSN.rand.sample, sampled.dSN)
        dNS.rand.sample <- append(dNS.rand.sample, sampled.dNS)
        dSS.rand.sample <- append(dSS.rand.sample, sampled.dSS)
      }
    } # finished creating sample of points of size N
    
    # Calculate bootstrap statistic: slope through the origin
    # Ref gene dN/dS 1: dNN/dSN
    lm.dNN.dSN <- lm(dNN.rand.sample ~ dSN.rand.sample - 1) # omitting intercept
    sample.slope.dNN.dSN <- summary(lm.dNN.dSN)$coefficients[1]
    
    # Ref gene dN/dS 2: dNS/dSS
    lm.dNS.dSS <- lm(dNS.rand.sample ~ dSS.rand.sample - 1) # omitting intercept
    sample.slope.dNS.dSS <- summary(lm.dNS.dSS)$coefficients[1]
    
    # Alt gene dN/dS 1: dNN/dNS
    lm.dNN.dNS <- lm(dNN.rand.sample ~ dNS.rand.sample - 1) # omitting intercept
    sample.slope.dNN.dNS <- summary(lm.dNN.dNS)$coefficients[1]
    
    # Alt gene dN/dS 2: dSN/dSS
    lm.dSN.dSS <- lm(dSN.rand.sample ~ dSS.rand.sample - 1) # omitting intercept
    sample.slope.dSN.dSS <- summary(lm.dSN.dSS)$coefficients[1]

    dNN.dSN.slopes[i] <- sample.slope.dNN.dSN
    dNS.dSS.slopes[i] <- sample.slope.dNS.dSS
    dNN.dNS.slopes[i] <- sample.slope.dNN.dNS
    dSN.dSS.slopes[i] <- sample.slope.dSN.dSS
    
  } # end this bootstrap replicate, i
  
  ##############################
  # Calculate SEMs, Z statistics, and P values
  
  # Ref gene dN/dS estimate 1: dNN.dSN.statistic, dNN.dSN.slope
  dNN.dSN.slope.SEM <- sd(dNN.dSN.slopes, na.rm=T)
  dNN.dSN.Z <- dNN.dSN.slope / dNN.dSN.slope.SEM
  dNN.dSN.P <- 2*pnorm(-abs(dNN.dSN.Z))
  cat("\n","REFERENCE dN/dS estimate 1:","\n","dNN/dSN=",dNN.dSN.statistic,"\n",
      "SLOPE dNN/dSN=",dNN.dSN.slope,"\n",
      "SE=",dNN.dSN.slope.SEM,"\n",
      "Z=",dNN.dSN.Z,"\n",
      "P=",dNN.dSN.P,"\n",
      sep="")
  
  # Ref gene dN/dS estimate 2: dNS.dSS.statistic, dNS.dSS.slope
  dNS.dSS.slope.SEM <- sd(dNS.dSS.slopes, na.rm=T)
  dNS.dSS.Z <- dNS.dSS.slope / dNS.dSS.slope.SEM
  dNS.dSS.P <- 2*pnorm(-abs(dNS.dSS.Z))
  cat("\n","REFERENCE dN/dS estimate 2:","\n","dNS/dSS=",dNS.dSS.statistic,"\n",
      "SLOPE dNS/dSS=",dNS.dSS.slope,"\n",
      "SE=",dNS.dSS.slope.SEM,"\n",
      "Z=",dNS.dSS.Z,"\n",
      "P=",dNS.dSS.P,"\n",
      sep="")
  
  # Alt gene dN/dS estimate 1: dNN.dNS.statistic, dNS.dNS.slope <- summary(dNN.dNS.lm)$coefficients[1]
  dNN.dNS.slope.SEM <- sd(dNN.dNS.slopes, na.rm=T)
  dNN.dNS.Z <- dNN.dNS.slope / dNN.dNS.slope.SEM
  dNN.dNS.P <- 2*pnorm(-abs(dNN.dNS.Z))
  cat("\n","ALTERNATIVE dN/dS estimate 1:","\n","dNN/dNS=",dNN.dNS.statistic,"\n",
      "SLOPE dNN/dNS=",dNN.dNS.slope,"\n",
      "SE=",dNN.dNS.slope.SEM,"\n",
      "Z=",dNN.dNS.Z,"\n",
      "P=",dNN.dNS.P,"\n",
      sep="")
  
  # Alt gene dN/dS estimate 2: dSN.dSS.statistic, dSN.dSS.slope 
  dSN.dSS.slope.SEM <- sd(dSN.dSS.slopes, na.rm=T)
  dSN.dSS.Z <- dSN.dSS.slope / dSN.dSS.slope.SEM
  dSN.dSS.P <- 2*pnorm(-abs(dSN.dSS.Z))
  cat("\n","ALTERNATIVE dN/dS estimate 2:","\n","dSN/dSS=",dSN.dSS.statistic,"\n",
      "SLOPE dSN/dSS=",dSN.dSS.slope,"\n",
      "SE=",dSN.dSS.slope.SEM,"\n",
      "Z=",dSN.dSS.Z,"\n",
      "P=",dSN.dSS.P,"\n",
      sep="")
  cat("\n")
}
