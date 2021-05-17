library("igraph")
library ("Matrix")

#' @title Maps initial expression scores to an interaction network
#'
#' @description This combines gene expression of various conditions with interaction network
#'
#' @param network species specific protein protein interaction list (typically first two columns:protein1 & protein2) from STRING database or species specific gene interaction list from any desired database.
#' @param expr_mat gene expression measures (absolute log2 fold changes or negative log p-values) or mutation profiles or any input values in a dataframe or matrix format. First column of expr_mat should contain gene id and other columns corresponding to other conditions (samples) that has to be propagated on the network.
#' @param type Type of graph that is created from the protein interaction list. By default: laplacian matrix is generated where each edge is normalized by the degree of interacting proteins.
#' @param merge.by string specifying a character vector in expr_mat (usually the first column name) that is matching with gene or protein names in interaction list
#' @param global TRUE indicates that all the input values of a condition should be mapped if the corresponding nodes are present in network. The nodes which doesnot have input values are assigned with 0 as their starting value.
#'
#' @export
#' @return list of laplacian matrix from the protein interaction list, intensity (initial scores) matrix matching to the nodes in network and node ids as gene_names
#'
#' @examples
#' data (expr_mat)
#' data (network)
#'
#' # To map the p-values from differential expression analysis to the STRING interaction network
#' netmap <- network_mapping(network = network, expr_mat = expr_mat, type = "laplacian", merge.by = "gene_id", global = TRUE)
#'
network_mapping <- function(network, expr_mat, type = "laplacian", merge.by = 'gene_id', global = TRUE)
{
  # Getting the unique nodes in the network from the input edgelist
  interaction_list <- cbind(as.vector(network[,1]), as.vector(network[,2]))
  ppi <- sort(unique(c(as.vector(network[,1]), as.vector(network[,2]))))
  ppi <- as.data.frame(ppi)
  colnames(ppi) <- merge.by

  # Making the expression matrix to the number of nodes in network (retaining the values only for the nodes present in the network and )
  expr_start <- merge(ppi, expr_mat, by=merge.by, all.x = global)
  expr_start <- expr_start[!duplicated(expr_start[,merge.by]),]
  mat_intensities <- apply (expr_start[,-1],2,as.numeric)
  mat_intensities[is.na(mat_intensities)] <- 0
  select_prot <- as.matrix(sort(expr_start[,c(merge.by)]))
  all_prot <- as.matrix(ppi)

  g <- construct_network(interaction_list, type)
  S <- g[all_prot %in% select_prot, all_prot %in% select_prot]
  return(list('G' = S, 'mat_intensities' = mat_intensities, "gene_names"=as.character(expr_start[,1])))
}

#' @title Builds a graph using a list of interacting genes or proteins
#'
#' @description construct_network is called within network_mapping for generating a laplacian matrix of graph from the protein interaction list
#'
#' @param edges a matrix or a data.frame with two columns listing pairs of interacting proteins or genes (column1 corresponds to protein1 and column2 corresponds to protein2)
#' @param type Type of graph that is created from the protein interaction list. By default: laplacian matrix is generated where each edge is normalized by the degree of interacting proteins. If the proteinA interacts with proteinB it gives 1 & upon normalization 1 will be divided by the square root of product of degrees of proteinA and proteinB.
#'
#' @export
#' @return Graph created from the edge list (interacting proteins)
#' @example
#'
construct_network <- function(edges, type="laplacian")
{
  nodes <- unique(sort(c(edges[,1], edges[,2])))
  size <- length(nodes)
  cat(paste("Network with:", size, "nodes\n"))
  #Construct interaction matrix

  indices_s <- cbind(match(edges[,1], nodes), match(edges[,2], nodes))
  indices_g <- cbind(match(edges[,2], nodes), match(edges[,1], nodes))

  indices <- rbind(indices_s, indices_g)

  mat <- Matrix(0, nrow = size, ncol= size, sparse = TRUE)
  mat[indices] <- 1
  diag(mat) <- 0


  if (type == "laplacian")
  {
    g <- graph.adjacency(as.matrix(mat), mode="undirected")
    L <- graph.laplacian(g, normalized=TRUE)
    G <- abs(L)
    return(G)
  }
  if (type == "adjacency")
  {
    G <- graph.adjacency(as.matrix(mat), mode="undirected")
    return(G)
  }
}

#' @title function for identifing active genes and their neighbours through propagating expression scores on an interaction network
#'
#' @description Using the propagate algebra function initial gene expression scores are spread on the network for every given condition. This propagation function will
#'
#' @param net igraph object respresenting the laplacian matrix (degree normalized adjacency matrix) of graph generated from protein interaction list. It is a symmetric matrix with equal number of rows and columns which is equal to the number of nodes in graph. "G" object from network_mapping output list.
#' @param mat_intensities matrix of mapped initial scores (gene expression values of each condition) corresponding to nodes in graph "G".
#' @param conditions character vector representing the sample names. This will used as column names of smooth matrix (output from this function).
#' @param iter an integer indicating number of runs of propagation of initial scores on network. At iter "t" every node propagates the scores received in the t-1 iteration to its neighbors. The node scores will get converged over the progression of iterations. The converged smooth matrix is found by taking the norm of smooth mat at t subtracted with smooth mat at t-1, the norm is exprected to be a very small value (in the ranges of 10^-6).
#' @param alpha ranging between 0 and 1, the fraction of intial scores that has to be diffused to the adjacent nodes. This is a tunable parameter. If 0.5 then 50 percent of intial score of the nodes will be spread to its neighbours, if it is 1 then 100 percent of initial scores will be spread to the adjacent nodes (complete loss of initial information).
#' @param network_type string specifying the type of graph that is created from the protein interaction list. By default: "laplacian".
#' @param stop_val if specified the scores are propagated in the network till the norm of smooth mat at t subtracted with smooth mat at t-1 reaches this value. By default:0, therefore it returns a matrix of smooth scores at the specified iteration.
#'
#' @return a matrix of propagated node scores at the specified alpha and iterations.
#' @export
#'
#' @examples
#' # To propagate the scores on the network with 20 runs and alpha of 0.5
#' # netmap is a list generated from the network_mapping function
#' netmap$smoothmat <- network_smoothing(net = netmap$G, mat_intensities = netmap$mat_intensities, conditions = colnames (netmap$mat_intensities), iter = 20, alpha = 0.5, network_type = "laplacian")
#'
#' # To generate propagated scores at 1-30 continuous iterations with alpha of 0.5
#' netmap$SmoothMats<- lapply (1:30, function (x) network_smoothing(net = netmap$G, mat_intensities = netmap$mat_intensities, conditions = colnames (netmap$mat_intensities), iter = x, alpha = 0.5, network_type="laplacian"))
#'
#' # To see the difference between the propagated matrices
#' norm (netmap$SmoothMats[[29]]- netmap$SmoothMats[[30]])
#'
#' # With this example the smoothed scores on network got converged at 30 iterations, therefore these scores could be used for furthur exploration
#' netmap$ConvergedSmoothMat <- netmap$SmoothMats[[30]]
#'
#' @references Vanunu et al. Associating Genes and Protein Complexes with Disease via Network Propagation. PLoS Comput Biol 2010, 6(1):e1000641.
#'
#' Hofree et al. Network-based stratification of tumor mutations. Nature Methods 2013, 10,1108–1115.
#'
#' Ruffalo et al. Network-Based Integration of Disparate Omic Data To Identify "Silent Players" in Cancer. PLoS Comput Biol 2015, 11(12): e1004595.
#'

network_smoothing <- function(net, mat_intensities, conditions, iter=1, alpha, network_type="laplacian", stop_val = 0)
{
  if (dim(mat_intensities)[1] == dim(net)[1])
  {
    if (stop_val != 0)
    {
      delta = 0
      i = 0
      while(delta < stop_val)
      {
        t_mat_intensities <- mat_intensities
        mat_intensities <- apply(t_mat_intensities, 2, function(x) as.vector(propagate_algebra(net, x, alpha, i)))
        delta <- mean((mat_intensities - t_mat_intensities)^2, na.rm = TRUE)
        cat(paste(delta, 'delta \n'))
        cat(paste(i, 'i \n'))
        i = i + 1
      }
      smooth_net <- apply(mat_intensities, 2 , I)
      smooth_net <- as.matrix(smooth_net)
      colnames(smooth_net) <- conditions
      smooth_net
    }
    else
    {
      smooth_net <- apply(mat_intensities, 2, function(x) as.vector(propagate_algebra(net, x, alpha, iter)))
      smooth_net <- apply(smooth_net, 2 , I)
      smooth_net <- as.matrix(smooth_net)
      colnames(smooth_net) <- conditions
      smooth_net
    }
  }
  else
  {
    cat("Please ensure correspondance between network and intensity matrix.")
  }

}

#' @title propagates node scores on a graph
#'
#' @description This is the actual function for diffusing the scores on the network. The sharing of initial scores to adjacent nodes in network is constrained by two parameters alpha and iter.
#'
#' @param W igraph object respresenting the laplacian matrix (degree normalized adjacency matrix) of graph generated from protein interaction list. It is a symmetric matrix with equal number of rows and columns which is equal to the number of nodes in graph.
#' @param z vector of mapped initial scores (gene expression values of each condition) corresponding to nodes in graph "G".
#' @param alpha ranging between 0 and 1, the fraction of intial scores that has to be diffused to the adjacent nodes. This is a tunable parameter. If 0.5 then 50 percent of intial score of the nodes will be spread to its neighbours, if it is 1 then 100 percent of initial scores will be spread to the adjacent nodes (complete loss of initial information).
#' @param iter an integer indicating number of runs of propagation of initial scores on network. At iter "t" every node propagates the scores received in the t-1 iteration to its neighbors. The node scores will get converged over the progression of iterations. The converged smooth matrix is found by taking the norm of smooth mat at t subtracted with smooth mat at t-1, the norm is exprected to be a very small value (in the ranges of 10^-6).
#'
#' @export
#' @return vector of propagated node scores at the specified alpha and iterations
#' @example
#'
propagate_algebra <- function(W, z, alpha, iter)
{
  F <- z
  W_p <- W * alpha
  diag(W_p) <- 0
  for (i in c(1:iter))
  {
    F <- (W_p %*% F) + (1 - alpha) * z
  }
  F
}

#' @title Builds an edge weighted graph using a list of interacting genes or proteins
#'
#' @description construct_weighted_network is called within network_mapping_weighted for generating a edge weighted laplacian matrix of graph from the protein interaction list
#'
#' @param edges a matrix or a data.frame with two columns listing pairs of interacting proteins or genes (column1 corresponds to protein1 and column2 corresponds to protein2)
#' @param type Type of graph that is created from the protein interaction list. By default: laplacian matrix is generated where each edge is normalized by the degree of interacting proteins. If the proteinA interacts with proteinB it gives 1 & upon normalization 1 will be divided by the square root of product of degrees of proteinA and proteinB.
#' @param weights numeric vector representing the weights of the edges. If the source of protein interactions is STRING, then combined score can be used as weights by dividing them with 1000. Typically weights gives the confidence of the interaction between nodes and ranges from 0 to 1, 1 being the 100 percent confidence of that the two proteins are interacting in the species.
#'
#' @return Edge weighted graph created from the protein interaction list.
#' @export
#'
#' @examples
#'
construct_weighted_network <- function(edges, type="laplacian", weights)
{
  nodes <- unique(sort(c(edges[,1], edges[,2])))
  size <- length(nodes)
  cat(paste("Network with:", size, "nodes\n"))

  # Creating indices for interactions
  indices_s <- cbind(match(edges[,1], nodes), match(edges[,2], nodes))
  indices_g <- cbind(match(edges[,2], nodes), match(edges[,1], nodes))
  indices <- rbind(indices_s, indices_g)

  # creating sparse square matrix to the size of network and whether there is interaction assigning 1 (essentially making adjacency mat of network)
  mat <- Matrix(0, nrow = size, ncol= size, sparse = TRUE)
  mat[indices] <- weights
  diag(mat) <- 0

  # Laplacian gives the propability of observing an edge between the same end-points in a random network with the same node degrees
  # if there is interaction between i&j, then (i,j) is 1, -w/sqrt(d[i] d[j])
  if (type == "laplacian")
  {
    g <- graph.adjacency(as.matrix(mat), mode="undirected", weighted = TRUE)
    L <- graph.laplacian(g, normalized=TRUE)
    G <- abs(L)
    return(G)
  }

  if (type == "adjacency")
  {
    G <- graph.adjacency(as.matrix(mat), mode="undirected", weighted = TRUE)
    return(G)
  }
}

#' @title Maps initial expression scores to an edge weighted interaction network
#'
#' @description This combines gene expression of various conditions with an edge weighted interaction network
#'
#' @param network species specific protein protein interaction list (typically first two columns:protein1 & protein2) from STRING database or species specific gene interaction list from any desired database.
#' @param expr_mat gene expression measures (absolute log2 fold changes or negative log p-values) or mutation profiles or any input values in a dataframe or matrix format. First column of expr_mat should contain gene id and other columns corresponding to other conditions (samples) that has to be propagated on the network.
#' @param weights numeric vector representing the weights of the edges. If the source of protein interactions is STRING, then combined score can be used as weights by dividing them with 1000. Typically weights gives the confidence of the interaction between nodes and ranges from 0 to 1, 1 being the 100 percent confidence of that the two proteins are interacting in the species.
#' @param type Type of graph that is created from the protein interaction list. By default: laplacian matrix is generated where each edge is normalized by the degree of interacting proteins. If the proteinA interacts with proteinB it gives 1 & upon normalization 1 will be divided by the square root of product of degrees of proteinA and proteinB.
#' @param merge.by string specifying a character vector in expr_mat (usually the first column name) that is matching with gene or protein names in interaction list
#' @param global TRUE indicates that all the input values of a condition should be mapped if the corresponding nodes are present in network. The nodes which doesnot have input values are assigned with 0 as their starting value.
#'
#' @return list of laplacian matrix from the protein interaction list, intensity (initial scores) matrix matching to the nodes in network and node ids as gene_names
#' @export
#'
#' @examples
#' # To construct edge weighted network using combined scores from STRING and mapping the expression values on to network
#' netmapWeighted <- network_mapping_weighted(network = network, expr_mat = expr_mat, weights = as.numeric(network$combined_score_weight), type = "laplacian", merge.by = "gene_id", global = TRUE)
#'
network_mapping_weighted <- function(network, expr_mat, weights, type = "laplacian", merge.by = 'gene_id', global = TRUE)
{
  # Getting the unique nodes in the network from the input edgelist
  interaction_list <- cbind(as.vector(network[,1]), as.vector(network[,2]))
  ppi <- sort(unique(c(as.vector(network[,1]), as.vector(network[,2]))))
  ppi <- as.data.frame(ppi)
  colnames(ppi) <- merge.by

  # Making the expression matrix to the size of nodes in network (retaining the values only for the nodes present in the network and )
  expr_start <- merge(ppi, expr_mat, by=merge.by, all.x = global)
  expr_start <- expr_start[!duplicated(expr_start[,merge.by]),]
  mat_intensities <- apply (expr_start[,-1],2,as.numeric)
  mat_intensities[is.na(mat_intensities)] <- 0
  select_prot <- as.matrix(sort(expr_start[,c(merge.by)]))
  all_prot <- as.matrix(ppi)

  # Construct network with the given edgelist matrix
  g <- construct_weighted_network(interaction_list, type, weights)
  S <- g[all_prot %in% select_prot, all_prot %in% select_prot]
  return(list('G' = S, 'mat_intensities' = mat_intensities, "gene_names"=as.character(expr_start[,1])))
}

### This is to account for scoring bias introduced on the hub nodes by smoothing
#' @title Network propagated scores are corrected by subtracted with mean inital propagated scores.
#'
#' @description Network propagation suffers from gene-specific biases created by the connections in the network ('topology bias'). For example, genes with many neighbors will generally tend to accumulate higher scores independent of the gene expression data. Therefore, network topology bias correction on the propagated scores is essential.
#' @param net igraph object respresenting the laplacian matrix (degree normalized adjacency matrix) of graph generated from protein interaction list. It is a symmetric matrix with equal number of rows and columns which is equal to the number of nodes in graph. "G" object from network_mapping output list.
#' @param mat_intensities matrix of mapped initial scores (gene expression values of each condition) corresponding to nodes in graph "G".
#' @param smooth_mat matrix of propagated scores for each gene for various conditions (converged smoothed matrix).
#' @param conditions character vector representing the sample names. This will used as column names of smooth matrix (output from this function).
#' @param alpha ranging between 0 and 1, the fraction of intial scores that has to be diffused to the adjacent nodes. This is a tunable parameter. If 0.5 then 50 percent of intial score of the nodes will be spread to its neighbours, if it is 1 then 100 percent of initial scores will be spread to the adjacent nodes (complete loss of initial information).
#' @param iter an integer value where the propagated scores has converged.
#'
#' @return a matrix with network bias corrected propagated scores for each condition.
#' @export
#'
#' @examples
#' # To correct for gene-sepcific biases on the propagated scores due to network structure
#' netmap$CorrectedSmoothMat <- network_bias_correction(net = netmap$G, mat_intensities = netmap$mat_intensities, smooth_mat = netmap$ConvergedSmoothMat, conditions = colnames (netmap$mat_intensities),alpha = 0.5, iter = 30)
#'
network_bias_correction <- function (net, mat_intensities, smooth_mat, conditions, alpha, iter)
{
  # This works only if the avg.value for conditions are above zero
  # To remove the gene id column in mapped matrix & also to remove na in the mapped mat (the genes that are in network but not in the input mat then assigned with NA)
  #mat_intensities <- as.matrix (mat_intensities[,-1])
  #mat_intensities [is.na(mat_intensities)]<- 0
  MatColMeans <- colMeans(as.matrix(mat_intensities))
  # To assign the mean value to all the nodes in network
  # replicate colmeans to the number of nodes in network
  t <- list()
  for (i in 1:length (MatColMeans)){
    t[[i]]<- rep(MatColMeans[i],nrow(mat_intensities))
  }
  MatColMeans <- do.call(cbind,t)
  # Smooth the mean values similarly as done for the real data
  SmoothCorrectionFactors <- network_smoothing(net, MatColMeans, conditions ,iter, alpha)

  # Subtracting the real scores with the avg.intial smooth scores on all nodes (Avg. value for each condition)
  network_bias_corrected_Mat <- smooth_mat - SmoothCorrectionFactors
  return (network_bias_corrected_Mat)
}
