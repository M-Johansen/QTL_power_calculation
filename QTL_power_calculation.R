#### QTL Power Calculations ####

## A function to calculate the power
## to detect QTLs for a given dataset
## 
#
# # Variables
# Type 1 error (alpha) = 0.01
# Sample size = n
# interval size = cM
# Phenotypic variance: h^2 (integer 0-1)

## Developed based on equations described in Hu & Xu (2008)
# https://www.nature.com/articles/hdy200825

## Parameters: 
#
# sampleSize: size of sample population used
# intervalcM: Interval to calculate power for in centi Morgans
# phenoVariance: ratio of phenotypic trait explained by the QTL
# ## e.g. 0.1 for explaining 10% of the variance, and 1 for explaining 100%
# cross = cross type / mating design
# choose from:
##### "simple" (default)
##### "BC" (backcross)
##### "DH" (double haploid)
##### "RIL" (recombinant inbred line)


QTL_Power_calc <- function(sampleSize, 
                           intervalcM,
                           phenoVariance,
                           cross = "simple")
{
  # Do stop(checks for input data)
  if (sampleSize <= 0) stop("Sample size must be a positive number.")
  if (intervalcM < 0) stop("Interval size must be 0 or greater (in cM).")
  if (phenoVariance < 0 || phenoVariance > 1) stop("Phenotypic variance must be between 0 and 1.")
  
  ## 1 # find recombination rate
  recomb = 0.5*(1-exp(-(intervalcM/sampleSize)))
 
  ## 2 # Calculate the Qx2
  # based on mating strategy:
  ## if cross = "simple"
  if(cross == "simple") {
    Qx2 = ((1-2*recomb)^2)/(2-4*recomb*(1-recomb))
  }
  
  # if cross = backcross ("BC")
  if(cross == "BC") {
    Qx2 = ((1-2*recomb)^2)/(4-8*(1-recomb)*recomb)
  }
  
  # if cross is double haploid (DH)
  if(cross == "DH") {
    Qx2 = 2-(1/(1-2*(1-recomb))*recomb)
  }
  
  # if cross us recombinant inbred line ("RIL")
  if(cross == "RIL") {
    Qx2 = 1 - ((2*recomb)/(1+4*recomb^2))
  }
  
  ## 3 # Calculate the squared QTL effect a^2
  # phenoVariance = 0.1
  effect <- (phenoVariance/(0.5*(1 - phenoVariance)))
  
  ## 4 # Calculate the non-centrality parameters
  nonCent <- sampleSize * Qx2 * (effect/(1 + effect * (0.5 - Qx2)))
  
  ## 5 # Calculate type 1 error (=0.01):
  type_1 <- qf(0.99, 1, sampleSize)
  
  ## 6 # Calculate the Type 2 error, F distribution
  type_2 <- pf(type_1, 1, sampleSize, nonCent)
  
  # 7 # Calculate statistical power to detect QTL
  power_r <- 1-type_2
  print(power_r)
}

#### Examples ##### 
##### To run the function on e.g., 100 samples, with 20 cM interval,
##### explaining 10% of variance, and cross = "simple" 

QTL_Power_calc(sampleSize = 100, 
               intervalcM = 20, 
               phenoVariance = 0.1,
               cross = "simple")
# or:
QTL_Power_calc(100, 20, 0.1)

# or with sample size = 116, cM interval = 5cM and 10% of variance
# and a back cross mating strategy
QTL_Power_calc(116, 5, 0.1, "BC")
