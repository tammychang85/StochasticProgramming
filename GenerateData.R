simulateDemands = function(){
  
  # static covariates
  x1 = rbinom(1, 1, 0.25)
  x2 = rbinom(1, 1, 0.5)
  x3 = rnorm(1, 900, 200)
  x4 = sample(seq(7.95, 22.95), 1)
  
  # demands of period 1 to 4
  d1 = 0
  d2 = 0
  d3 = 0
  d4 = 0
  while(TRUE){
    d1 = (738 + 1072 * x1 - 403 * x2 + 2.8 * x3 - 55 * x4 + rnorm(1, 0, 1044))
    if(d1 >= 0){break}
  }
  while(TRUE){
    d2 = (-399 - 638 * x1 + 53 * x4 + 0.854 * d1 + rnorm(1, 0, 781))
    if(d2 >= 0){break}
  }
  while(TRUE){
    d3 = (-5 + 0.955 * d2 + rnorm(1, 0, 601))
    if(d3 >= 0){break}
  }
  while(TRUE){
    d4 = (874 - 46 * x4 + 0.516 * d2 + 0.318 * d3 + rnorm(1, 0, 820))
    if(d4 >= 0){break}
  }
  
  return(data.frame(d1, d2, d3, d4, x1, x2, x3, x4))
}
getRealizations = function(realization.size){
  
  realization.df = data.frame()
  for (i in 1:realization.size) {
    sample.i = simulateDemands()
    realization.df = rbind(realization.df, sample.i)
  }
  realization.matrix = rbind(rep(0, realization.size))
  realization.matrix = rbind(realization.matrix, t(data.matrix(realization.df[, 1:4])))
  
  return(list(frame=realization.df, matrix=realization.matrix))
}

realization.size = 50
testing.size = 200
realizations = getRealizations(realization.size)
testingDataSet = getRealizations(testing.size)$matrix
saveRDS(realizations, 'data/realizations.rds')
saveRDS(testingDataSet, 'data/testingDataSet.rds')
