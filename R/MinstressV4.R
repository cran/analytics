#' @title Better Starting Configuration For Non-Metric MDS
#' @description \code{Minstress} is a heuristic to find better non-metric MDS solutions,
#' by finding better starting configurations, instead of just using a random one.
#' @details This function performs several iterations, each using a different starting seed,
#' and in turn each one of those iterations performs non-metric MDS many times (typically, thousands or more)
#' in an attempt to find the best seed (which induces a particular initial configuration) of them all.
#' @param x a data frame containing numeric values only
#' @param p the size of the population of seeds (any positive integer)
#' @param s the number of seeds we sample (any positive integer)
#' @param k the number of dimensions wanted (any positive integer)
#' @param iter a positive integer specifying the number of iterations.
#' @param pb a Boolean variable declaring if one wants to display a pogress bar (default: False)
#' @param m a string specifying the distance method (default: 'euclidean')
#' @return A list informing about dimensionality, minimum STRESS level found, and best seed found.
#' One can then use the best seed found to perform non-metric MDS with a better initial configuration (generally).
#' @author Albert Dorador
#' @export
#' @import cluster
#' @import tcltk
#' @import MASS
#' @import stats
#' @examples
#'
#' require(MASS)
#'
#' swiss.x <- as.data.frame(swiss[, -1])
#' Minstress(swiss.x, 1e5, 50, 2, iter = 3)
#'
#' # Comparing without using Minstress (for such a low value of s, difference is minimal)
#' swiss.x <- as.matrix(swiss[, -1])
#' swiss.dist <- dist(swiss.x)
#' swiss.mds <- isoMDS(swiss.dist)
#'

Minstress <- function(x, p, s, k, iter = 5, pb = F, m = 'euclidean'){
  stopifnot(is.numeric(c(p, s, k)))
  if (m == "gower"){
    Dist <- as.matrix(daisy(x, metric = "gower"))
  } else {
    Dist <- as.matrix(dist(x, method = m))
  }
  n <- nrow(x)
  best.seeds <- integer(iter)
  Stress.vector <- numeric(iter)
  if (pb == T){
    progbar <- tkProgressBar(title = "Progress bar", min = 0, max = iter, width = 300)
    setTkProgressBar(progbar, 0, label = paste(0,"% done"))
  }
  for (it in 1:iter){
    ind <- sample.int(p, s)
    STRESS <- numeric(length(ind))
    for (i in ind){
      set.seed(i)
      init <- scale(matrix(runif(n*k), ncol = k), scale = FALSE)
      nmmds.out <- isoMDS(Dist, y = init, k = k, maxit = 100, trace = F)
      STRESS[which(ind == i)] <- nmmds.out$stress
    }
    if (pb == T) setTkProgressBar(progbar, it, label = paste(round(it/iter*100, 0),"% done"))
    best.seeds[it] <- ind[which.min(STRESS)]
    Stress.vector[it] <- min(STRESS)
  }
  if (pb == T){
    Sys.sleep(0.5)
    close(progbar)
  }
  solution <- list('Dimensionality'= k, 'Minimum found'= min(Stress.vector),
                   'Best seed'= best.seeds[which.min(Stress.vector)])
  return(solution)
}
