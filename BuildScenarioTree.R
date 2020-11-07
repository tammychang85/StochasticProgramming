buildtree <- function(x, treeStruct, lambda_0 = 10, lambda_f = 0.01,
                      e_0 = 0.5, e_f = 0.05, jMax = 40000, plot = FALSE){
  
  
  ## -------- SET UP FUNCTIONS --------##
  getEucDist <- function(scenario, member){
    EucDist <- sqrt(sum((scenario - member)^2))
    return(EucDist)
  }
  step_size <- function(e_0, e_f, jMax, j){
    e_0 * ( (e_f / e_0) ^ (j / jMax))
  }
  getLambda_j <- function(lambda_f, lambda_0, jMax, j){
    lambda_j <- lambda_0 * ( (lambda_f/lambda_0) ^ (j / jMax))
    return(lambda_j)
  }
  getH <- function(O, lambda_j){
    h <- exp(-O/lambda_j)
    return(h)
  }
  
  ## -------- INITIALIZE TREE NODES -------- ##
  numScenarios <- ncol(treeStruct)
  tree <- x[,sample(ncol(x), numScenarios), drop=F]    # select x at random
  for (i in 1:max(treeStruct)){ #max(treeStruct) gives the number of nodes in the tree
    tree[which(treeStruct == i)] <- mean(tree[which(treeStruct == i)])
  }
  
  ## -------- ITERATE NODES -------- ##
  for (j in 1:jMax){
    randMember <- x[,sample(ncol(x), 1)] # Sample a single member from the x
    Distances <- apply(tree, 2, getEucDist, member = randMember)
    O <- rank(Distances)
    for (node in 1:max(treeStruct)){
      adapt <- getH(O[which(treeStruct==node, arr.ind=TRUE)[,2]],getLambda_j(lambda_f, lambda_0, jMax, j))
      diff <- randMember[which(treeStruct==node, arr.ind = TRUE)[,1]] - tree[which(treeStruct==node)][1]
      delta_node <- step_size(e_0, e_f, jMax, j) * mean(adapt * diff)
      tree[which(treeStruct==node)] <- tree[which(treeStruct==node)] + delta_node
    }
  }
  
  
  ## -------- CALCULATE BRANCH PROBABILITIES -------- ##
  branchProbs <- vector("integer", ncol(tree))
  for (m in 1:ncol(x)){
    dist <- apply(tree, 2, getEucDist, member = x[,m])
    branchProbs[which.min(dist)] <- branchProbs[which.min(dist)] + 1 / ncol(x)
  }
  
  ## -------- OUTPUT -------- ##
  if (plot) {
    matplot(x, type = "l", col = "grey", lty = 2, ylab = "Disturbance", xlab = "Time step")
    matlines(tree, pch = 3, lty = 1)
  }
  output <- list(treeStruct,tree, branchProbs)
  names(output) <- c("tree_structure", "tree_values", "branch_probabilities")
  return(output)
}