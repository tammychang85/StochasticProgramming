# ---- packages ----
library(DescTools)

# ---- functions ----
getFeatures = function(){
  # static covariates
  x1 = rbinom(1, 1, 0.25)
  x2 = rbinom(1, 1, 0.5)
  x3 = rnorm(1, 900, 200)
  x4 = sample(seq(7.95, 22.95), 1)
  
  return(data.frame(x1, x2, x3, x4))
}
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
getRealizationR = function(realization.size){
  
  realization.df = data.frame()
  for (i in 1:realization.size) {
    sample.i = simulateDemands()
    realization.df = rbind(realization.df, sample.i)
  }
  realization.matrix = rbind(rep(0, realization.size))
  realization.matrix = rbind(realization.matrix, t(data.matrix(realization.df[, 1:4])))
  # realization.matrix = t(realization.matrix)
  
  return(list(frame=realization.df, matrix=realization.matrix))
}
# using all the features to estimate for each period
getEstimatedDemands = function(X0, similar.product.datas, time.period=4){
  
  new.product.demands = list()
  
  ## perform least-squares regression on available data on the n historical demands of similar products
  # period 1
  lm.1 = lm(d1~x1+x2+x3+x4, similar.product.datas)
  estimated.d01 = predict(lm.1, X0) + lm.1$residual
  estimated.d01[which(estimated.d01 < 0)] = 0
  new.product.demands = append(new.product.demands, list(estimated.d01))
  
  # period 2
  lm.2 = lm(d2~x1+x2+x3+x4+d1, similar.product.datas)
  estimated.d02 = c()
  for (d in estimated.d01){
    X0.t = cbind(X0, d1=d)
    estimated.d02 = c(estimated.d02, predict(lm.2, X0.t))
  }
  estimated.d02 = estimated.d02 + lm.2$residual
  estimated.d02[which(estimated.d02 < 0)] = 0
  new.product.demands = append(new.product.demands, list(estimated.d02))
  
  # period 3
  lm.3 = lm(d3~x1+x2+x3+x4+d1+d2, similar.product.datas)
  estimated.d03 = c()
  for (i in seq_along(estimated.d02)){
    X0.t = cbind(X0, d1=estimated.d01[i], d2=estimated.d02[i])
    estimated.d03 = c(estimated.d03, predict(lm.3, X0.t))
  }
  estimated.d03 = estimated.d03 + lm.3$residual
  estimated.d03[which(estimated.d03 < 0)] = 0
  new.product.demands = append(new.product.demands, list(estimated.d03))
  
  # period 4
  lm.4 = lm(d4~x1+x2+x3+x4+d1+d2+d3, similar.product.datas)
  estimated.d04 = c()
  for (i in seq_along(estimated.d03)){
    X0.t = cbind(X0, d1=estimated.d01[i], d2=estimated.d02[i], d3=estimated.d03[i])
    estimated.d04 = c(estimated.d04, predict(lm.4, X0.t))
  }
  estimated.d04 = estimated.d04 + lm.4$residual
  estimated.d04[which(estimated.d04 < 0)] = 0
  new.product.demands = append(new.product.demands, list(estimated.d04))
  
  return(new.product.demands)
}
binDemands = function(demands, bin.num, realization.size, mode='default', probs=c(0.05, 0.95), test=0){
  
  # if to winsorize the demands
  if (mode == 'w'){
    demands = lapply(demands, function(x){Winsorize(x, probs=probs)})
  }
  
  # visulalize the bins if need
  if (test){
    x11(width=70,height=30)
    par(mfrow=c(2,2))
  }
  
  # bin demands into bin.num bins for each period
  binned.demands = list()
  binned.demands.probs = list()
  
  period = 1
  for (demand in demands){
    bins = cut(demand, bin.num)
    bins.median = tapply(demand, bins, median)
    binned.demands = append(binned.demands, list(bins.median))
    binned.demands.probs = append(binned.demands.probs, list(as.vector(table(bins)) / realization.size))
    
    # visulalize the bins if need
    if (test){
      barplot(table(bins), names.arg=lapply(bins.median, round), main=paste0('period', period)) 
      period = period + 1
    }
    
  }
  
  return(list(values=binned.demands, probs=binned.demands.probs))
}
getResidualTree = function(realizations, bin.num, realization.size, mode='default', probs=c(0.05, 0.95), test=0) {
  
  features.x0 = getFeatures()
  demands.x0 = getEstimatedDemands(features.x0, realizations$frame)
  bin.Demands = binDemands(demands.x0, bin.num, realization.size, mode=mode, probs=probs, test=test)
  
  demands = bin.Demands$values
  realization.size = length(demands[[1]])
  probabilities = bin.Demands$probs
  paths.df = data.frame()
  probs.df = data.frame()
  for (i in (1:(length(demands) - 1))) {
    if (i == 1){
      for (j in seq_along(demands[[i]])) {
        curret.demand = demands[[i]][j]
        next.demand = demands[[i + 1]]
        path.i = data.frame(curret.demand, next.demand, row.names=NULL)
        paths.df = rbind(paths.df, path.i)
        
        curret.prob = probabilities[[i]][j]
        next.prob = probabilities[[i + 1]]
        prob.i = data.frame(curret.prob, next.prob, row.names=NULL)
        probs.df = rbind(probs.df, prob.i)
      }
      paths.df = apply(paths.df, 1, toString)
      probs.df = apply(probs.df, 1, toString)
    }
    else {
      current.paths = data.frame()
      current.probs = data.frame()
      for (j in seq_along(paths.df)) {
        curret.period = rep(paths.df[j], realization.size)
        next.period = demands[[i + 1]]
        path.i = data.frame(curret.period, next.period, row.names=NULL)
        current.paths = rbind(current.paths, path.i)
        
        curret.prob = rep(probs.df[j], realization.size)
        next.prob = probabilities[[i + 1]]
        prob.i = data.frame(curret.prob, next.prob, row.names=NULL)
        current.probs = rbind(current.probs, prob.i) 
      }
      paths.df = apply(current.paths, 1, toString)
      probs.df = apply(current.probs, 1, toString)
    }
  }
  
  # convert the paths from string to float
  tree_values = list()
  branch_probabilities = c()
  for (i in seq_along(paths.df)) {
    tree_values = append(tree_values, list(c(0, as.double(strsplit(paths.df[i], ',')[[1]]))))
    branch_probabilities = c(branch_probabilities, prod(as.double(strsplit(probs.df[i], ',')[[1]])))
  }
  tree_values = do.call(cbind, tree_values)
  
  return(list(tree_values=tree_values, branch_probabilities=branch_probabilities))
}


getEstimatedDemands1 = function(X0, similar.product.datas, time.period=4){
  
  new.product.demands = list()
  
  ## perform least-squares regression on available data on the n historical demands of similar products
  # period 1
  lm.1 = lm(d1~x1, similar.product.datas)
  estimated.d01 = predict(lm.1, X0) + lm.1$residual
  estimated.d01[which(estimated.d01 < 0)] = 0
  new.product.demands = append(new.product.demands, list(estimated.d01))
  
  # period 2
  lm.2 = lm(d2~x1, similar.product.datas)
  estimated.d02 = predict(lm.2, X0) + lm.2$residual
  estimated.d02[which(estimated.d02 < 0)] = 0
  new.product.demands = append(new.product.demands, list(estimated.d02))
  
  # period 3
  lm.3 = lm(d3~x1, similar.product.datas)
  estimated.d03 = predict(lm.3, X0) + lm.3$residual
  estimated.d03[which(estimated.d03 < 0)] = 0
  new.product.demands = append(new.product.demands, list(estimated.d03))
  
  # period 4
  lm.4 = lm(d4~x1, similar.product.datas)
  estimated.d04 = predict(lm.4, X0) + lm.4$residual
  estimated.d04[which(estimated.d04 < 0)] = 0
  new.product.demands = append(new.product.demands, list(estimated.d04))
  print(lm.4)
  return(new.product.demands)
}
getEstimatedDemands2 = function(X0, similar.product.datas, time.period=4){
  
  new.product.demands = list()
  
  ## perform least-squares regression on available data on the n historical demands of similar products
  # period 1
  lm.1 = lm(d1~x1+x2, similar.product.datas)
  estimated.d01 = predict(lm.1, X0) + lm.1$residual
  estimated.d01[which(estimated.d01 < 0)] = 0
  new.product.demands = append(new.product.demands, list(estimated.d01))
  
  # period 2
  lm.2 = lm(d2~x1+x2, similar.product.datas)
  estimated.d02 = predict(lm.2, X0) + lm.2$residual
  estimated.d02[which(estimated.d02 < 0)] = 0
  new.product.demands = append(new.product.demands, list(estimated.d02))
  
  # period 3
  lm.3 = lm(d3~x1+x2, similar.product.datas)
  estimated.d03 = predict(lm.3, X0) + lm.3$residual
  estimated.d03[which(estimated.d03 < 0)] = 0
  new.product.demands = append(new.product.demands, list(estimated.d03))
  
  # period 4
  lm.4 = lm(d4~x1+x2, similar.product.datas)
  estimated.d04 = predict(lm.4, X0) + lm.4$residual
  estimated.d04[which(estimated.d04 < 0)] = 0
  new.product.demands = append(new.product.demands, list(estimated.d04))
  print(lm.4)
  return(new.product.demands)
}
getEstimatedDemands3 = function(X0, similar.product.datas, time.period=4){
  
  new.product.demands = list()
  
  ## perform least-squares regression on available data on the n historical demands of similar products
  # period 1
  lm.1 = lm(d1~x1+x2+x3, similar.product.datas)
  estimated.d01 = predict(lm.1, X0) + lm.1$residual
  estimated.d01[which(estimated.d01 < 0)] = 0
  new.product.demands = append(new.product.demands, list(estimated.d01))
  
  # period 2
  lm.2 = lm(d2~x1+x2+x3, similar.product.datas)
  estimated.d02 = predict(lm.2, X0) + lm.2$residual
  estimated.d02[which(estimated.d02 < 0)] = 0
  new.product.demands = append(new.product.demands, list(estimated.d02))
  
  # period 3
  lm.3 = lm(d3~x1+x2+x3, similar.product.datas)
  estimated.d03 = predict(lm.3, X0) + lm.3$residual
  estimated.d03[which(estimated.d03 < 0)] = 0
  new.product.demands = append(new.product.demands, list(estimated.d03))
  
  # period 4
  lm.4 = lm(d4~x1+x2+x3, similar.product.datas)
  estimated.d04 = predict(lm.4, X0) + lm.4$residual
  estimated.d04[which(estimated.d04 < 0)] = 0
  new.product.demands = append(new.product.demands, list(estimated.d04))
  print(lm.4)
  return(new.product.demands)
}
getEstimatedDemands4 = function(X0, similar.product.datas, time.period=4){
  
  new.product.demands = list()
  
  ## perform least-squares regression on available data on the n historical demands of similar products
  # period 1
  lm.1 = lm(d1~x1+x2+x3+x4, similar.product.datas)
  estimated.d01 = predict(lm.1, X0) + lm.1$residual
  estimated.d01[which(estimated.d01 < 0)] = 0
  new.product.demands = append(new.product.demands, list(estimated.d01))
  
  # period 2
  lm.2 = lm(d2~x1+x2+x3+x4, similar.product.datas)
  estimated.d02 = predict(lm.2, X0) + lm.2$residual
  estimated.d02[which(estimated.d02 < 0)] = 0
  new.product.demands = append(new.product.demands, list(estimated.d02))
  
  # period 3
  lm.3 = lm(d3~x1+x2+x3+x4, similar.product.datas)
  estimated.d03 = predict(lm.3, X0) + lm.3$residual
  estimated.d03[which(estimated.d03 < 0)] = 0
  new.product.demands = append(new.product.demands, list(estimated.d03))
  
  # period 4
  lm.4 = lm(d4~x1+x2+x3+x4, similar.product.datas)
  estimated.d04 = predict(lm.4, X0) + lm.4$residual
  estimated.d04[which(estimated.d04 < 0)] = 0
  new.product.demands = append(new.product.demands, list(estimated.d04))
  print(lm.4)
  return(new.product.demands)
}
# used to test the effect of choosing different features
getResidualTree2 = function(realizations, bin.num, realization.size, mode=1) {
  
  features.x0 = getFeatures()
  demands.x0 = NULL
  if (mode == 1){
    demands.x0 = getEstimatedDemands1(features.x0, realizations$frame)
  }
  else if (mode == 2){
    demands.x0 = getEstimatedDemands2(features.x0, realizations$frame)
  }
  else if (mode == 3){
    demands.x0 = getEstimatedDemands3(features.x0, realizations$frame)
  }
  else{
    demands.x0 = getEstimatedDemands4(features.x0, realizations$frame)
  }

  bin.Demands = binDemands(demands.x0, bin.num, realization.size)
  
  demands = bin.Demands$values
  probabilities = bin.Demands$probs
  
  realization.size = length(demands[[1]])
  paths.df = data.frame()
  probs.df = data.frame()
  for (i in (1:(length(demands) - 1))) {
    if (i == 1){
      for (j in seq_along(demands[[i]])) {
        curret.demand = demands[[i]][j]
        next.demand = demands[[i + 1]]
        path.i = data.frame(curret.demand, next.demand, row.names=NULL)
        paths.df = rbind(paths.df, path.i)
        
        curret.prob = probabilities[[i]][j]
        next.prob = probabilities[[i + 1]]
        prob.i = data.frame(curret.prob, next.prob, row.names=NULL)
        probs.df = rbind(probs.df, prob.i)
      }
      paths.df = apply(paths.df, 1, toString)
      probs.df = apply(probs.df, 1, toString)
    }
    else {
      current.paths = data.frame()
      current.probs = data.frame()
      for (j in seq_along(paths.df)) {
        curret.period = rep(paths.df[j], realization.size)
        next.period = demands[[i + 1]]
        path.i = data.frame(curret.period, next.period, row.names=NULL)
        current.paths = rbind(current.paths, path.i)
        
        curret.prob = rep(probs.df[j], realization.size)
        next.prob = probabilities[[i + 1]]
        prob.i = data.frame(curret.prob, next.prob, row.names=NULL)
        current.probs = rbind(current.probs, prob.i) 
      }
      paths.df = apply(current.paths, 1, toString)
      probs.df = apply(current.probs, 1, toString)
    }
  }
  
  # convert the paths from string to float
  tree_values = list()
  branch_probabilities = c()
  for (i in seq_along(paths.df)) {
    tree_values = append(tree_values, list(c(0, as.double(strsplit(paths.df[i], ',')[[1]]))))
    branch_probabilities = c(branch_probabilities, prod(as.double(strsplit(probs.df[i], ',')[[1]])))
  }
  tree_values = do.call(cbind, tree_values)
  
  return(list(tree_values=tree_values, branch_probabilities=branch_probabilities))
}

# ---- start ----
realization.size = 50
realizations = getRealizationR(realization.size)
bin.num = 3
set.seed(1)
residualtree = getResidualTree(realizations, bin.num, realization.size)
set.seed(1)
residualtree2 = getResidualTree(realizations, bin.num, realization.size, mode='w', probs=c(0.01, 0.99))
getResidualTree2(realizations, 2, realization.size, 1)

features.x0 = getFeatures()
demands.x0 = getEstimatedDemands(features.x0, realizations$frame)

bin.Demands = binDemands(demands.x0, 11, realization.size, test=1)
bin.Demands2 = binDemands(demands.x0, 6, realization.size, mode='w', test=1)
bin.Demands3 = binDemands(demands.x0, 6, realization.size, mode='w',probs=c(0.1, 0.90), test=1)

x11(width=70,height=30)
par(mfrow=c(2,2))
for (i in 1:4){
  hist(demands.x0[[i]], main=paste0('Org_period', i), xlab='demand')
}

lapply(demands.x0, summary)
lapply(demands.x0, function(x){table(cut(x, 11))})
cut(demands.x0[[1]], 11)
sum(table(cut(demands.x0[[1]], 11)))
