AbundToMatrix <- function (all, sample, cogs, group = "yes", log = "yes") { # default to ln(OR)

# process all bacteria data and return numerical matrix
  
  all_bact <- sapply(all, AllBactMat)
  all_bact <- as.matrix(Reduce(cbind,all_bact)) # convert list to abundance matrix
  rownames(all_bact) <- cogs
  
  # process all bacteria data and return numerical matrix
  
  sample[, c("Func_id", "Func_name") := NULL]
  sample <- as.matrix(sample)
  rownames(sample) <- cogs

# get sum information needed for odds ratio calculations
  
  AllBactCOGSums <- rowSums(all_bact) # vector of totals for each COG in all bacteria
  AllCOGS <- sum(AllBactCOGSums) # total of all COGS in all bacteria
#  SampleOrgSum <- colSums(samp) # vector of all COGS in each sample
#  SampleCOGSum <- rowSums(sample) # vector of totals for each COG in all samples  

# group sample data and add to matrix if necessary

  if(group == "yes") {samp <- cbind(sample, rowSums(sample))} # combine each sample into a single sample and add to sample matrix

# calculate odds ratios

  mat <- apply( # calculate matrix of ln(OR)
    sample, 
    2, 
    CalcOddsRatios,
    AllBactCOGSums = AllBactCOGSums # vector of COG sums for all bacteria
  ) 
}

AllBactMat <- function (file) { # read abund file, remove COG/name columns, return data frame
  DF <- fread(file) # load file into data frame
  drops <- c("Func_id", "Func_name") # identify columns to drop
  DF[,(drops) := NULL] # drop columns
  DF
}

CalcOddsRatios <- function (samp, AllBactCOGSums) { # calculate odds ratios and return matrix on ln(OR)
  SampleOrgSum = sum(samp) # number of COGS in org
  AllCOGSum = sum(AllBactCOGSums) # total number COGs in all bacteria
  data <- sapply(
    1:length(samp), 
    function(x) {
      OR_obj <- oddsratio(c(samp[x],SampleOrgSum), n=c(AllBactCOGSums[x],AllCOGSum))
      OR_obj$OR
    }
  )
}

# all = list of COG abundance files for all bacteria (>1000 bacteria per file) 
# that will be concatenated into a single numerical abundance matrix
# sample = sample COG abundance matrix to be compared to all bacteria
# cogs = vector of COG category names (same length as all and sample)
# group = in addition to individual analysis, group all samples together as single sample (default = "yes")
# log = return matrix of ln(odds ratio) (default = "yes")

# Dependencies
#   data.table
#   bstats