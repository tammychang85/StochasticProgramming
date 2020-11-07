# ---- packages ----
# to bin the data
library(rbin)
# solver
library(Rglpk)


# ---- functions ----
## generating demands
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

## residual tree
getFeatures = function(){
  # static covariates
  x1 = rbinom(1, 1, 0.25)
  x2 = rbinom(1, 1, 0.5)
  x3 = rnorm(1, 900, 200)
  x4 = sample(seq(7.95, 22.95), 1)
  
  return(data.frame(x1, x2, x3, x4))
}
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
binDemands = function(demands, bin.num, realization.size, test=0){
  
  # bin demands into bin.num bins for each period
  binned.demands = list()
  binned.demands.probs = list()
  
  # if (test){
  #   x11(width=70,height=30)
  #   par(mfrow=c(2,2))
  # }
  
  period = 1
  for (d in demands){
    demand.dt = data.frame(demand=d)
    bin.groups = rbin_equal_length(demand.dt, demand, demand, bin.num)
    breaks = bin.groups$lower_cut[1]
    breaks = c(breaks, bin.groups$upper_cut)
    bins = cut(d, breaks, right = F)
    bins.median = tapply(d, bins, median)
    binned.demands = append(binned.demands, list(bins.median))
    
    bins.count = table(bins)
    binned.demands.probs = append(binned.demands.probs, list(as.vector(table(bins)) / realization.size))
    if (test){
      barplot(bins.count, names.arg=lapply(bins.median, round), main=paste0('period', period)) 
      period = period + 1
    }
    
  }
  
  return(list(values=binned.demands, probs=binned.demands.probs))
}
getResidualTree = function(realizations, bin.num, realization.size) {
  
  features.x0 = getFeatures()
  demands.x0 = getEstimatedDemands(features.x0, realizations$frame)
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


## scenario tree
# receive a vector of nodal structure and return a matrix of tree structure
getTreeStructure = function(node.per.period){
  # how many period does a scenario tree have
  tree.length = length(node.per.period)
  scenario.num = node.per.period[tree.length]
  # build the tree structure
  tree.structure = c()
  node.value = 1
  for(i in node.per.period){
    for(j in 1:i){
      tree.structure = c(tree.structure, rep(node.value, scenario.num / i))
      node.value = node.value + 1
    }
  }
  
  return(matrix(tree.structure, nrow=tree.length, byrow=TRUE))
}
# to build the scenario tree
import::here(buildtree, .from = "BuildScenarioTree.R")
# use modified function of scenario package, able to build tree with only one scenario
getScenarioTree = function(node.per.period, realizations, maxIteration=40000){
  
  treeStructure = getTreeStructure(node.per.period)
  
  return(buildtree(realizations$matrix, treeStructure, jMax=maxIteration))
}


## inear programming
getCostStructure = function(fixedCost=1, holdingCost=0.25, flexible.k=1.5, penalty.k=2){
  
  flexibleCost = fixedCost * flexible.k
  penaltyCost = flexibleCost * penalty.k
  salvageValue = holdingCost
  
  return(c(fixedCost, flexibleCost, holdingCost, penaltyCost, salvageValue))
}
getObjFunction = function(probability, node.per.period, cost.structure, flexible=1) {
  # how many paths does the scenario tree have
  scenario.num = node.per.period[length(node.per.period)]
  period.num = length(node.per.period) - 1
  
  # the objective function
  obj = c()
  
  # coefficient of fixed supplier
  c1qt1 = c()
  for (i in 1:period.num) {
    c1qt1 = c(c1qt1, cost.structure[1] * probability)  }
  # coefficient of inventory
  htIt = c()
  for (i in 1:(period.num-1)) {
    htIt = c(htIt, cost.structure[3] * probability)
  }
  # coefficient of lost
  ptLt = c()
  for (i in 1:(period.num)) {
    ptLt = c(ptLt, cost.structure[4] * probability)
  }
  # coefficient of salvaged value (or inventory, depends on the sign of cost sructure)
  sIT = c()
  sIT = c(sIT, cost.structure[5] * probability)
  
  # results if flexible supplier not allowed
  if (flexible != 1){
    obj = c(c1qt1, htIt, ptLt, sIT)
  }
  # results if flexible supplier allowed
  else {
    # coefficient of flexible supplier
    c2qt2 = c()
    for (i in 1:period.num) {
      c2qt2 = c(c2qt2, cost.structure[2] * probability)
    }
    obj = c(c1qt1, c2qt2, htIt, ptLt, sIT)
  }
  
  return(obj)
}
getConstraintsMatrix = function(obj, node.per.period, flexible=1, mode=1) {
  # how many scenarios does a scenario tree have
  scenario.num = node.per.period[length(node.per.period)]
  period.num = length(node.per.period) - 1
  # total coefficient number of the objective function
  coefficient.num = length(obj)
  # number of available suppliers
  supplier.num = 0
  # constraints vector to return
  constraints = c()
  
  # if flexible supplier allowed
  if (flexible == 1){
    supplier.num = 2
    if(scenario.num != 1){
      ## constraits for fixed supplier of all period + flexible supplier of period 1
      # should have (period.num + 1) * (scenario.num - 1) constraints
      for(t in 1:(period.num + 1)){
        for(k in 1:(scenario.num - 1)){
          constraint = rep(0, coefficient.num)
          # coefficient of qt1 & q12 , k=1
          constraint[(t - 1) * scenario.num + 1] = 1
          # coefficient of  qt1 & q12 , k=2~scenario.num
          constraint[(t - 1) * scenario.num + 1 + k] = -1
          constraints = c(constraints, constraint)
        }
      }
      
      ## constraits for flexible supplier after period 1
      # for scenario tree
      if (mode == 1){
        # should have node.per.period[t] * ((scenario.num/node.per.period[t]) - 1) constraints
        for (t in 1:(period.num - 1)){
          for(i in 1:node.per.period[t + 1]){
            # how many scenarios does a node have at time t
            node.scenario.num = scenario.num / node.per.period[t + 1]
            # if a node has more than one scenario
            if ((node.scenario.num - 1) != 0){
              for(k in 1:(node.scenario.num-1)){
                constraint = rep(0, coefficient.num)
                constraint[scenario.num * (period.num + t) + node.scenario.num * (i - 1) + 1] = 1
                constraint[scenario.num * (period.num + t) + node.scenario.num * (i - 1) + 1 + k] = -1
                constraints = c(constraints, constraint)
              } 
            }
          }
        }
      }
      
      # for residual tree
      else{
        for (i in 1:(period.num - 1)) {
          group.num = 2^i
          path.num.group = scenario.num / group.num
          for (group in 1:group.num) {
            for (path in 1:(path.num.group - 1)) {
              constraint = rep(0, coefficient.num)
              constraint[((period.num + 1 + (i - 1)) * scenario.num) + (1 + path.num.group * (group - 1))] = 1
              constraint[((period.num + 1 + (i - 1)) * scenario.num) + ((path + 1) + path.num.group * (group - 1))] = -1
              constraints = c(constraints, constraint)
            }
          }
        }
      }
      
    }
    
    # constraits for invetory and lost, should have scenario.num *  period.num constraints
    for (i in 1:(period.num)){
      # number of coefficients of suppliers
      supplier.coef.num = scenario.num * period.num * supplier.num
      for (j in 1:(scenario.num)){
        # constraints for fixed & flexible supplier and inventory & lost of first period
        if ((i == 1)){
          I = rep(0, coefficient.num)
          # fiexd supplier
          I[j] = 1
          # flexible supplier
          I[j + period.num * scenario.num] = 1
          # inventory of t=1
          I[supplier.coef.num + j] = -1
          # lost of t=1
          I[supplier.coef.num + j + ((period.num - 1) * scenario.num)] = 1
          constraints = c(constraints, I)
        }
        
        # constraints for inventory & lost between the first and last period (exclude the two periods)
        else if (i != period.num){
          I = rep(0, coefficient.num)
          # fiexd supplier
          I[((i - 1) * scenario.num) + j] = 1
          # flexible supplier
          I[((i - 1) * scenario.num) + j + period.num * scenario.num] = 1
          # invetory of t-1
          I[supplier.coef.num + ((i - 2) * scenario.num) + j] = 1
          # inventory of t
          I[supplier.coef.num + ((i - 1) * scenario.num) + j] = -1
          # lost of t
          I[supplier.coef.num + ((i - 1) * scenario.num) + j + ((period.num - 1) * scenario.num)] = 1
          constraints = c(constraints, I)
        }
        
        # constraints for inventory & lost of the last period
        else{
          I = rep(0, coefficient.num)
          # fixed supplier
          I[((i - 1) * scenario.num) + j] = 1
          # flexible supplier
          I[((i - 1) * scenario.num) + j + period.num * scenario.num] = 1
          # invetory of T - 1
          I[supplier.coef.num + ((i - 2) * scenario.num) + j] = 1
          # invetory of T
          I[supplier.coef.num + ((2 * period.num - 1) * scenario.num) + j] = -1
          # lsot of T
          I[supplier.coef.num + (2 * (period.num - 1) * scenario.num) + j] = 1
          constraints = c(constraints, I)
        } 
      }
    }
  }
  
  # if flexible supplier not allowed
  else {
    supplier.num = 1
    if(scenario.num != 1){
      # constraits for fixed supplier of all period
      # should have (period.num) * (scenario.num - 1) constraints
      for(t in 1:(period.num)){
        for(k in 1:(scenario.num - 1)){
          constraint = rep(0, coefficient.num)
          # coefficient of constraintt1-1
          constraint[(t - 1) * scenario.num + 1] = 1
          # coefficient of constraintt1-i+1
          constraint[(t - 1) * scenario.num + 1 + k] = -1
          constraints = c(constraints, constraint)
        }
      }
    }
    
    # constraits for invetory and lost, should have scenario.num *  period.num constraints
    for (i in 1:(period.num)){
      # coefficient number of suppliers
      supplier.coef.num = scenario.num * period.num * supplier.num
      for (j in 1:(scenario.num)){
        # constraints for fixed & flexible supplier and inventory & lost of first period
        if ((i == 1)){
          I = rep(0, coefficient.num)
          # fiexd supplier
          I[j] = 1
          # invetory
          I[supplier.coef.num + j] = -1
          # lost
          I[supplier.coef.num + j + ((period.num - 1) * scenario.num)] = 1
          constraints = c(constraints, I)
        }
        # constraints for inventory & lost between the first and last period (exclude the two periods)
        else if (i != period.num){
          I = rep(0, coefficient.num)
          # fiexd supplier
          I[((i - 1) * scenario.num) + j] = 1
          # invetory of t - 1
          I[supplier.coef.num + ((i - 2) * scenario.num) + j] = 1
          # invetory of t
          I[supplier.coef.num + ((i - 1) * scenario.num) + j] = -1
          # lost of t
          I[supplier.coef.num + ((i - 1) * scenario.num) + j + ((period.num - 1) * scenario.num)] = 1
          constraints = c(constraints, I)
        }
        # constraints for inventory & lost of the last period
        else{
          I = rep(0, coefficient.num)
          # fiexd supplier
          I[((i - 1) * scenario.num) + j] = 1
          # invetory of T - 1
          I[supplier.coef.num + ((i - 2) * scenario.num) + j] = 1
          # invetory of T
          I[supplier.coef.num + ((2 * period.num - 1) * scenario.num) + j] = -1
          # lost of T
          I[supplier.coef.num + (2 * (period.num - 1) * scenario.num) + j] = 1
          constraints = c(constraints, I)
        } 
      }
    }
  }
  
  return(matrix(constraints, ncol=coefficient.num, byrow=TRUE))
}
# return the right hand side of the constraints
getRHS = function(constraints.matrix, node.per.period, scenario.tree) {
  
  # how many scenarios does a scenario tree have
  scenario.num = node.per.period[length(node.per.period)]
  period.num = length(node.per.period) - 1
  
  # the right hand side of constraints for suppliers
  rhs = rep(0, (dim(constraints.matrix)[1]) - (scenario.num * period.num))
  for(i in 1:period.num){
    rhs = c(rhs, scenario.tree$tree_values[(i + 1), ])
  }
  
  return(rhs)
}
# optimize and return the results for scenario tree
optimize = function(scenario.tree, node.per.period, cost.structure, flexible=1, mode=1){
  # obejct function for optimization
  obj = getObjFunction(scenario.tree$branch_probabilities, node.per.period, cost.structure, flexible)
  # constraints matrix for optimization
  constraints.matrix = constraints.matrix = getConstraintsMatrix(obj, node.per.period, flexible, mode)
  # directions of the constraints
  dir = rep('==' ,dim(constraints.matrix)[1])
  # the right hand side of the constraints
  rhs = getRHS(constraints.matrix, node.per.period, scenario.tree)
  
  # optimize and return
  return(Rglpk_solve_LP(obj, constraints.matrix, dir, rhs, max=FALSE))
}


## seeing and checking the optimization model 
getColumnName = function(period.num, scenario.num, flexible=1) {
  
  col.name = c()
  # qt1
  for(i in 1:period.num){
    for(j in 1:scenario.num){
      time = as.character(i)
      scenario = as.character(j)
      col.name = c(col.name, paste0('fi', time, scenario))
    }
  }
  if (flexible == 1){
    # qt2
    for(i in 1:period.num){
      for(j in 1:scenario.num){
        time = as.character(i)
        scenario = as.character(j)
        col.name = c(col.name, paste0('fl', time, scenario))
      }
    } 
  }
  # invetory
  for(i in 1:(period.num-1)){
    for(j in 1:scenario.num){
      time = as.character(i)
      scenario = as.character(j)
      col.name = c(col.name, paste0('I', time, scenario))
    }
  }
  # lost
  for(i in 1:(period.num)){
    for(j in 1:scenario.num){
      time = as.character(i)
      scenario = as.character(j)
      col.name = c(col.name, paste0('L', time, scenario))
    }
  }
  # salvaged value
  for(j in 1:scenario.num){
    time = period.num
    scenario = as.character(j)
    col.name = c(col.name, paste0('S', time, scenario))
  }
  return(col.name)
}
checkObjectFunction = function(obj, node.per.period, flexible=1){
  # how many scenarios does the scenario tree have
  scenario.num = node.per.period[length(node.per.period)]
  # how many period of demands does the tree have
  period.num = length(node.per.period) - 1
  #
  m = matrix(obj, nrow=1)
  col.name = getColumnName(period.num, scenario.num, flexible)
  colnames(m) = col.name
  
  return(m)
}
checkConstraintMetrix = function(constraints.matrix, node.per.period, flexible=1){
  
  # how many paths does the scenario tree have
  scenario.num = node.per.period[length(node.per.period)]
  # how many period of demands does the tree have
  period.num = length(node.per.period) - 1
  
  col.name = getColumnName(period.num, scenario.num, flexible)
  row.name = c()
  
  # constraints of suppliers if total sceanrio num is more than one
  if(scenario.num != 1){
    # constraints of fixed supplier, should have period.num * (scenario.num - 1)
    for(i in 1:period.num){
      for (j in 1:(scenario.num-1)) {
        time.count = as.character(i)
        scenario.count = as.character(j)
        row.name = c(row.name, paste0('fixed', time.count, '-', scenario.count)) 
      }
    }
    
    # constraints of flexible supplier if allowed
    if (flexible == 1){
      # should have node.per.period[t] * ((scenario.num/node.per.period[t])-1)
      for(i in 1:(period.num)){
        # how many scenarios does a node have at time i
        node.scenario.num = node.per.period[i] * ((scenario.num / node.per.period[i]) - 1)
        if (node.scenario.num != 0){
          for(s in 1:node.scenario.num){
            time.count = as.character(i)
            scenario.count = as.character(s)
            row.name = c(row.name, paste0('flexible', time.count, '-', scenario.count))
          }
        }
      }
    }
  }
  
  # constraints of invetory and lost, should have period.num * scenario.num
  for (i in 1:period.num){
    for (j in 1:scenario.num){
      timeCount = as.character(i)
      count = as.character(j)
      row.name = c(row.name, paste0('I&L', timeCount, '-', count))
    } 
  }
  
  # add the column and row names
  colnames(constraints.matrix) = col.name
  rownames(constraints.matrix) = row.name
  View(constraints.matrix)
  
  return(constraints.matrix)
}
checkRHS = function(RHS, constraints.matrix){
  return(matrix(RHS, ncol=1, byrow=TRUE))
}


## parsing optimizaton results
# parse the solution form optimization model and return the order policy matrix
getOrderPolicy = function(node.per.period, solutions, flexible=1) {
  
  scenario.num = node.per.period[length(node.per.period)]
  period.num = length(node.per.period)-1
  s.num = 0 # number of available supplier
  
  # set supplier number depends on if flexible suppler allowed
  if (flexible == 1){
    s.num = 2
  }
  else{
    s.num = 1
  }
  
  # extract only solutions of order policy
  policy = solutions$solution[1:(s.num * scenario.num * period.num)]
  
  
  # make a matrix of fixed policy
  fixed.policy = policy[1:(period.num * scenario.num)]
  fixed.policy.matrix = rep(NA, period.num)
  index = 1
  for (i in 1:period.num){
    # record policy and remove duplicated ones
    fixed.policy.matrix[i] = unique(fixed.policy[index:(i * scenario.num)])
    index = i * scenario.num + 1
  }
  
  # make a matrix of flexible policy (values equal to 0 if not flexibe supplier allowed)
  flexible.policy.matrix = matrix(0, nrow=period.num, ncol=scenario.num)
  if (flexible == 1){
    flexible.costs = list()
    flexible.policy = policy[(period.num * scenario.num + 1):length(policy)]
    index = 1
    for (i in 1:period.num){
      flexible.costs[[i]] = flexible.policy[index:(i * scenario.num)]
      index = i * scenario.num + 1
    }
    for(i in seq_along(flexible.costs)){
      for(j in 1:scenario.num){
        flexible.policy.matrix[i, j] = flexible.costs[[i]][j]
      }
    } 
  }
  
  # combine fixed and flexible policies
  final.policy = cbind(flexible.policy.matrix, fixed.policy.matrix)
  
  return(final.policy)
}
# decide which flexible options to apply in each perid for a given scenario and return the policy matrix
decideFlexiblePolicy = function(scenario.tree, testingData, policy) {
  # record the decision to the final strategy for each period
  scenario.num = dim(scenario.tree$tree_values)[2]
  possible.scenarios = scenario.tree$tree_values[-length(testingData), , drop=F]
  possible.flexible.policy = policy[,-dim(policy)[2], drop=F]
  
  # decide decison for each period of the scenario
  for (i in 1:(length(testingData)-1)) {
    possible.outcomes.i = possible.scenarios[i, ]
    # get demand distance between testing data and scenario tree for period j 
    distance = as.matrix(dist(matrix(c(testingData[i], possible.outcomes.i))))[,1]
    # exclude the distance with itself
    distance = distance[2:length(distance)]
    # find the closest one to choose the decision
    min.distance = min(distance)
    min.distance.index = which(distance==min.distance)[1]
    chosen.outcome = possible.outcomes.i[min.distance.index]
    
    # choose policy
    possible.scenarios.index = possible.scenarios[i,] %in% chosen.outcome
    possible.scenarios =  possible.scenarios[, possible.scenarios.index, drop=F]
    possible.flexible.policy = possible.flexible.policy[, possible.scenarios.index, drop=F]
    if(dim(possible.scenarios)[2] == 1){break}
  }
  
  # only return one if have duplicated ones
  if(dim(possible.flexible.policy)[2] > 1){
    possible.flexible.policy = possible.flexible.policy[,1]
  }
  if(!is.null(dim(possible.flexible.policy))){
    possible.flexible.policy = possible.flexible.policy[,1]
  }
  
  return(possible.flexible.policy)
}


## for costs
# calculate cost for a given testing scneario
getCost = function(scenario.tree, testingData, cost.structure, policy) {
  #
  fixed.policy = policy[, dim(policy)[2]]
  flexible.policy = decideFlexiblePolicy(scenario.tree, testingData, policy)
  # record cost of each period
  costs = c()
  inventory = 0
  
  # get cost for each period
  for (i in 1:(length(testingData)-1)) {
    # available items of period i
    quantity = inventory + fixed.policy[i] + flexible.policy[i]
    # calculate inventory or lost 
    difference = quantity - testingData[i+1]
    # cost of period i
    cost.i = 0
    
    # if not the last period
    if (i != (length(testingData)-1)){
      # if have unmet demands
      if (difference < 0){
        # cost of period i
        cost.i = fixed.policy[i] * cost.structure[1] + flexible.policy[i] * cost.structure[2] +  abs(difference) * cost.structure[4]
        inventory = 0
      }
      # if have invetory
      else{
        cost.i = fixed.policy[i] * cost.structure[1] + flexible.policy[i] * cost.structure[2] + difference * cost.structure[3]
        inventory = difference
      }
      
    }
    else{
      # if have unmet demands
      if (difference < 0){
        cost.i = fixed.policy[i] * cost.structure[1] + flexible.policy[i] * cost.structure[2] + abs(difference) * cost.structure[4]
      }
      # if have invetory
      else{
        cost.i = fixed.policy[i] * cost.structure[1] + flexible.policy[i] * cost.structure[2] + difference * cost.structure[5]
      }
    }
    costs = c(costs, cost.i)
  }
  
  
  return(list(costs, flexible.policy, fixed.policy))
}
# calculate the average cost of a given set of testing scenarios
getAverageCost = function(scenario.tree, testingDataSet, cost.structure, policy) {
  # record costs for all testing data
  costs = rep(0, 4)
  flexible.order = rep(0, 4)
  fixed.order = rep(0, 4)
  
  # get the total cost for all testing data
  for ( i in 1: (dim(testingDataSet)[2])){
    # get the cost of testing data i
    result.i = getCost(scenario.tree, testingDataSet[, i], cost.structure, policy)
    costs = costs + result.i[[1]]
    flexible.order = flexible.order + result.i[[2]]
    fixed.order = fixed.order + result.i[[3]]
  }
  costs = costs / dim(testingDataSet)[2]
  flexible.order = flexible.order / dim(testingDataSet)[2]
  fixed.order = fixed.order / dim(testingDataSet)[2]
  
  return(list(cost=costs, flexible.order=flexible.order, fixed.order=fixed.order))
}


## coordinate functions above to get the results of simulation for a given tree
# mode： 1 for scenario tree; 0 for residual tree
getTestResults = function(tree, testingDataSet, node.per.period, cost.structure, flexible=1, mode=1) {
  
  # get optimal solutions
  solutions = optimize(tree, node.per.period, cost.structure, flexible, mode)
  # parse order policy
  policy = getOrderPolicy(node.per.period, solutions, flexible)
  # results on testing datas
  results = getAverageCost(tree, testingDataSet, cost.structure, policy)

  return(results)
}

# ---- main ----

## simply test the functions
# read the data
realizations = readRDS('data/realizations.rds')
testingDataSet = readRDS('data/testingDataSet.rds')

st.node.test = c(1, 2, 4, 8, 8)
st = getScenarioTree(st.node.test, realizations)

rt.node.test = c(1, 2, 4, 8, 16)
bin.num = 2
rt = getResidualTree(realizations, bin.num, realization.size)

cost.structure = getCostStructure()
st.results = getTestResults(st, testingDataSet, st.node.test, cost.structure, mode=1)
rt.results = getTestResults(rt, testingDataSet, rt.node.test, cost.structure, mode=0)
sum(st.results[[1]])
sum(rt.results[[1]])


## compare secnario and residual trees with diferents numbers of bin
# read the data
realizations = readRDS('data/realizations.rds')
testingDataSet = readRDS('data/testingDataSet.rds')

# tree sturctures
st.node = c(1, 2, 4, 8, 16)
rt2.node = c(1, 2, 4, 8, 16)
rt4.node = c(1, 4, 16, 64, 256)

# start
rounds = 2
flexible.costs = c(1.5, 3, 6)
high.cost.structure = rep(list(0), length(flexible.costs)) # cost structure with high peanlty
low.cost.structure = rep(list(0), length(flexible.costs)) # cost structure with low peanlty
for (i in seq_along(flexible.costs)) {
  high.cost.structure[[i]] = getCostStructure(flexible.k=flexible.costs[i], penalty.k=4.5)
  low.cost.structure[[i]] = getCostStructure(flexible.k=flexible.costs[i], penalty.k=1.5)
}

# scenario tree
st.results.high = rep(list(rep(0, 3)), 3)
st.results.low = rep(list(rep(0, 3)), 3)
for (round in 1:rounds) {
  
  print(paste0('---round ', round, '---'))
  
  print('get scenario tree')
  st = getScenarioTree(st.node, realizations)

  print('get results')
  for (i in seq_along(flexible.costs)) {
    st.results.high.r = getTestResults(st, testingDataSet, st.node, high.cost.structure[[i]], mode=1)
    st.results.high.r = lapply(st.results.high.r, sum)
    for (j in seq_along(st.results.high.r)) {
      st.results.high[[i]][j] = st.results.high[[i]][j] + st.results.high.r[[j]]
    }
    
    st.results.low.r = getTestResults(st, testingDataSet, st.node, low.cost.structure[[i]], mode=1)
    st.results.low.r = lapply(st.results.low.r, sum)
    for (j in seq_along(st.results.low.r)) {
      st.results.low[[i]][j] = st.results.low[[i]][j] + st.results.low.r[[j]]
    }
  }
  
  if (round == rounds){
    print('avearge results')
    st.results.high = lapply(st.results.high, function(x){x / rounds})
    st.results.low = lapply(st.results.low, function(x){x / rounds})
    
    saveRDS(st.results.high, 'results/stResultsHigh.rds')
    saveRDS(st.results.low, 'results/stResultsLow.rds')
  }  
}

# residual tree, bin = 2
rt2.results.high = rep(list(rep(0, 3)), 3)
rt2.results.low = rep(list(rep(0, 3)), 3)
for (round in 1:rounds) {
  
  print(paste0('---round ', round, '---'))
  
  print('get residual tree, bin = 2')
  rt = getResidualTree(realizations, 2, 50)
  
  print('get results')
  for (i in seq_along(flexible.costs)) {
    rt2.results.high.r = getTestResults(rt, testingDataSet, rt2.node, high.cost.structure[[i]], mode=0)
    rt2.results.high.r = lapply(rt2.results.high.r, sum)
    for (j in seq_along(rt2.results.high.r)) {
      rt2.results.high[[i]][j] = rt2.results.high[[i]][j] + rt2.results.high.r[[j]]
    }
    
    rt2.results.low.r = getTestResults(rt, testingDataSet, rt2.node, low.cost.structure[[i]], mode=0)
    rt2.results.low.r = lapply(rt2.results.low.r, sum)
    for (j in seq_along(rt2.results.low.r)) {
      rt2.results.low[[i]][j] = rt2.results.low[[i]][j] + rt2.results.low.r[[j]]
    }
  }
  
  if (round == rounds){
    print('avearge results')
    rt2.results.high = lapply(rt2.results.high, function(x){x / rounds})
    rt2.results.low = lapply(rt2.results.low, function(x){x / rounds})
    
    saveRDS(rt2.results.high, 'results/rt2ResultsHigh.rds')
    saveRDS(rt2.results.low, 'results/rt2ResultsLow.rds')
  }   
}

# residual tree, bin = 4
rt4.results.high = rep(list(rep(0, 3)), 3)
rt4.results.low = rep(list(rep(0, 3)), 3)
for (round in 1:rounds) {
  
  print(paste0('---round ', round, '---'))
  
  print('get residual tree, bin = 4')
  rt = getResidualTree(realizations, 4, 50)
  
  print('get results')
  for (i in seq_along(flexible.costs)) {
    rt4.results.high.r = getTestResults(rt, testingDataSet, rt4.node, high.cost.structure[[i]], mode=0)
    rt4.results.high.r = lapply(rt4.results.high.r, sum)
    for (j in seq_along(rt4.results.high.r)) {
      rt4.results.high[[i]][j] = rt4.results.high[[i]][j] + rt4.results.high.r[[j]]
    }
    
    rt4.results.low.r = getTestResults(rt, testingDataSet, rt4.node, low.cost.structure[[i]], mode=0)
    rt4.results.low.r = lapply(rt4.results.low.r, sum)
    for (j in seq_along(rt4.results.low.r)) {
      rt4.results.low[[i]][j] = rt4.results.low[[i]][j] + rt4.results.low.r[[j]]
    }
  }
  
  if (round == rounds){
    print('avearge results')
    rt4.results.high = lapply(rt4.results.high, function(x){x / rounds})
    rt4.results.low = lapply(rt4.results.low, function(x){x / rounds})
    
    saveRDS(rt4.results.high, 'results/rt4ResultsHigh.rds')
    saveRDS(rt4.results.low, 'results/rt4ResultsLow.rds')
  }   
}

# visulize the results
ratio.st.rt2.high = c()
ratio.st.rt2.low = c()
ratio.st.rt4.high = c()
ratio.st.rt4.low = c()
for (i in seq_along(flexible.costs)) {
  ratio.st.rt2.high =c(ratio.st.rt2.high, (st.results.high[[i]][1] / rt2.results.high[[i]][1]))
  ratio.st.rt2.low =c(ratio.st.rt2.low, (st.results.low[[i]][1] / rt2.results.low[[i]][1]))
  
  ratio.st.rt4.high =c(ratio.st.rt4.high, (st.results.high[[i]][1] / rt4.results.high[[i]][1]))
  ratio.st.rt4.low =c(ratio.st.rt4.low, (st.results.low[[i]][1] / rt4.results.low[[i]][1]))
}

png('graphs/AllFeatures.png')
x11(width=70,height=30)
par(mfrow=c(1,2))
plot(1:length(flexible.costs), ratio.st.rt2.high, type='b',lty=2, lwd=2, col='blue',
     xlab='flexible cosst', ylab='cost ratio', xaxt='n', yaxt='n', ylim=c(0.5, 1.2), main='high penalty')
axis(1, at=1:length(flexible.costs), labels=flexible.costs, legend())
axis(2, at=seq(0.6, 1.2, 0.05))
legend('bottomright', legend=c('st / rt (bin=2)', 'st / rt (bin=4)'), col=c('blue', 'red'), text.col=c('blue', 'red'), lty=2, lwd=2, cex = 0.85)
lines(1:length(flexible.costs), ratio.st.rt4.high, type='b', lty=2, lwd=2, col='red')

plot(1:length(flexible.costs), ratio.st.rt2.low, type='b',lty=2, lwd=2, col='blue',
     xlab='flexible cosst', ylab='cost ratio', xaxt='n', yaxt='n', ylim=c(0.5, 1.2), main='low penalty')
axis(1, at=1:length(flexible.costs), labels=flexible.costs)
axis(2, at=seq(0.6, 1.2, 0.05))
legend('bottomright', legend=c('st / rt (bin=2)', 'st / rt (bin=4)'), col=c('blue', 'red'), text.col=c('blue', 'red'), lty=2, lwd=2, cex = 0.85)
lines(1:length(flexible.costs), ratio.st.rt4.low, type='b', lty=2, lwd=2, col='red')
dev.off()
