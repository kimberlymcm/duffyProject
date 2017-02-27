# Util functions & data for featuresel_results_github.r

# Real data 
real.stats <- c(1.168, 14, -1.385, 6.6897, -5.5216, 3217.87, 4985, 5, 2, 8, 0.569, 0.018, 0.7417,
               0.0316, 12, 1, 0.67345, 0.001003009)
names(real.stats) <- c('pi', 'seg', 'D', 'theta', 'H', 'iHH', 'fixed', 'singletons', 'doubletons', 
                      'X2', 'H1', 'H2', 'H12', 'H1H2', 'uniq', 'finalFreq', 'ehh.avg', 'sf')

# Takes in a filename and returns cleaned data
readInFile <- function(fileName){
  # Reads in simulation result file and performs some cleaning
  #
  #   Args:
  #     Filename: Expects a file with the summary statistics listed in the names(data) list
  #               Each line of the input file represents 1 simulation & the summary stats
  # 
  #   Returns:
  #     data.frame containing simulations and summary statistics
  # 
  data = read.table(fileName, skipNul = TRUE)
  # s: selection coefficient, pi: # pairwise differences, seg: # segregating sites
  # theta: theta, H: Fay & Wu's H, ehh.L: ehh on left side of simulated region,
  # ehh.R: ehh on right side of simulated region, fixed: # fixed sites, 
  # singletons: # singletons, doubletons: # doubletons, X2: # non-singletons/doubletons, 
  # H1: H1 statistic, H2: H2 statistic, H12: H12 statistic, H1H2: H1H2 statistic, 
  # numHaps: # haplotypes, uniq: # unique haplotypes, finalFreq: final frequency 
  # of selected allele
  names(data) = c('s', 'pi', 'seg', 'D', 'theta', 'H', 'ehh.L', 'iHH', 'ehh.R', 
                  'fixed', 'singletons', 'doubletons', 'X2', 'H1' , 'H2', 'H12', 
                  'H1H2', 'numHaps', 'uniq', 'finalFreq')
  # Data cleaning steps, a few simulations over-wrote each other & we remove those
  data = data[which(data$s != 0),]
  data = data[which(data$finalFreq < 10),]
  # Get average ehh on the sequence ends
  data$ehh.avg = (data$ehh.L + data$ehh.R) / 2
  # Remove unwanted statistics
  data = data[, (!names(data) %in% c('ehh.L', 'ehh.R', 'numHaps'))]
  # sf: # singletons / # fixed sites
  data$sf = data$singletons / data$fixed
  if ((data$finalFreq == 0) || (data$finalFreq == 1)){
    data$seg = data$seg - 1
  }
  return(data)
}
