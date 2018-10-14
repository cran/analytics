#' @title Give Column Means of a Matrix-like Object, Based on a Grouping Variable
#' @description Compute column (weighted) means across rows of a numeric matrix-like object for each level of a grouping variable.
#' @details This function is a wrapper for base function \code{rowsum} which allows one to compute the (weighted) mean instead of the sum,
#' while handling integer overflow.
#'
#' Note: although data frames ara allowed, keep in mind that data frames do not allow duplicate row names. Hence if you have a dataframe with more than 1 group, you may want to use the function as.matrix() to convert it to an object of class matrix
#'
#' To compute the mean over all the rows of a matrix (i.e. a single group) use colMeans, which should be even faster.
#' @param M a matrix, data frame or vector of numeric data. Missing values are allowed. A numeric vector will be treated as a column vector.
#' @param group a vector or factor giving the grouping, with one element per row of M. Default: rownames of M.
#' @param w a vector giving the weights that must be applied to each of the stacked blocks of an original object
#' @param reord if TRUE, then the result will be in order of sort(unique(group)), if FALSE (the default), it will be in the order that groups were encountered.
#' @param na_rm logical (TRUE or FALSE). Should NA (including NaN) values be discarded?
#' @param big is your object big and integer overflow is likely? If TRUE, then M is multiplied by 1.0 to ensure values are of type double (perhaps taking more RAM).
#' @param ... other arguments to be passed to or from methods.
#' @return A matrix-like object containing the means by group. There will be one row per unique value of group.
#' If object supplied in fact (explicitly) had just one group, base function
#' \code{colMeans} is called for maximum efficiency and a numeric vector containing
#' the mean of each column is returned.
#'
#' @author Albert Dorador
#' @export
#' @seealso
#' \code{\link[base]{rowsum}}
#' @examples
#'
#' A <- matrix(1:8, ncol = 2)
#' rownames(A) <- c("A", "B", "A", "B")
#' rowmean(A)
#'
#' B <- matrix(1:40, ncol = 2)
#' gr <- rep(1:5, 4)
#' B.mean <- rowmean(B, group = gr)
#' sum(B.mean[, 1])*4 == sum(B[, 1]) #basic sanity check
#' sum(B.mean[, 2])*4 == sum(B[, 2]) #basic sanity check
#'
#' dfB <- as.data.frame(B)
#' gr <- rep(1:5, 4)
#' dfB.mean <- rowmean(dfB, group = gr)
#'
#' numbers <- rnorm(1e7, mean = 3)
#' C <- matrix(numbers, ncol = 5)
#' gr <- rep(1:20, 1e5)
#' rowmean(C, group = gr) # Handles Big Data fast
#'
#' vec <- 1:10
#' gr <- rep(1:2, 5)
#' rowmean(vec, gr)
#'
#' onegroup = matrix(1:40, ncol = 2)
#' gr = rep(1,20)
#' rowmean(onegroup, gr)[1] == mean(onegroup[,1])
#' rowmean(onegroup, gr)[2] == mean(onegroup[,2])
#'
#' numbers <- rnorm(30, mean = 3)
#' D <- matrix(numbers, ncol = 3)
#' num_blocks <- 2
#' gr <- rep(1:5, num_blocks)
#' rownames(D) <- gr
#' rowmean(D, w = c(0.1,0.9))
#' rowmean(D, w = c(0,1))
#' rowmean(D, w = c(0.5,0.5))
#' rowmean(D)
#'

rowmean <- function(M, group=rownames(M), w = FALSE, reord=FALSE, na_rm=FALSE, big=TRUE, ...){
  if (is.null(group))
    stop("group cannot be NULL. Possibly: M has no rownames.")
  group <- as.factor(group)
  gr.count <- unlist(lapply(split(group, f = group), length))

  if (big) {
    M <- M*1.0 #to prevent integer overflow
  }

  if (length(w) == 1) {
    if (length(gr.count) == 1) { #in case there's really just 1 group
      return(colMeans(x = M, na.rm = na_rm))
    }
    M.sum <- rowsum(M, group, reorder = reord, na.rm = na_rm, ...)
    return(M.sum / gr.count)

  } else { #assuming we deal with stacked blocks of an original object
    if (sum(w) != 1)
      stop("weights must sum to 1.")
    aux <- NA
    flag <- FALSE
    i <- 1
    while (!flag & (i <= nrow(M))) {
      aux[i] <- rownames(M)[i]
      flag <- any(duplicated(aux)) #find first repeated rowname
      i <- i + 1
    }
    len_block <- i-2

    full_w <- rep(w, each = len_block)
    M.w <- M * full_w
    M.sum <- rowsum(M.w, group, reorder = reord, na.rm = na_rm, ...)
    return(M.sum)
  }
}

