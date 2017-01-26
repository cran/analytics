#' @title Takes Outliers Off
#' @description \code{offliers} Finds the existing outliers after fitting a linear model, according to various criteria.
#' @details Criteria available: Cook's Distance, DFBetas, and COVRATIO.
#' DFFits has not been included because it is conceptually equivalent to Cook's Distance; in fact
#' there's a closed-form formula to convert one value to the other. See references.
#' The user can select any combination of those, and take off the outliers in the intersection (default) or union.
#' @references Cook & Weisberg 1982, "Residuals and Influence in Regression".
#' \url{http://conservancy.umn.edu/handle/11299/37076}
#' @param dataset an object containing the data used in the linear regression model, typically a data frame.
#' @param mod an object of class "lm" (the model fitted by the user).
#' @param CD a Boolean variable, indicating whether the Cook's Distance criterion is to be used.
#' @param DFB a Boolean variable, indicating whether the DFBetas criterion is to be used.
#' @param COVR a Boolean variable, indicating whether the COVRATIO criterion is to be used.
#' @param pctg a real number between 0 and 100, indicating the maximum percentage of original observations to be removed.
#' @param intersection a Boolean variable, indicating whether the intersection or the union of the outliers detected
#' must be considered (in case more than one criterion has been selected).
#' @return a list containing, according to each criterion selected, which observations have been identified as outliers,
#' how many they are, what percentage of the total number of observations they represent, and, if more than one criterion
#' has been selected, the final outliers, quantity and percentage.
#'
#' @author Albert Dorador
#' @export
#' @import stats
#' @examples
#'
#' require(graphics)
#' ## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
#' ## Page 9: Plant Weight Data.
#' ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
#' trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
#' group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
#' weight <- c(ctl, trt)
#' lm.D9 <- lm(weight ~ group)
#' db <- data.frame(weight, group)
#'
#' offliers(db,lm.D9)
#' offliers(db,lm.D9, CD = FALSE, DFB = TRUE, COVR = TRUE)
#' offliers(db,lm.D9, CD = TRUE, DFB = TRUE, COVR = TRUE, intersection = FALSE)
#' offliers(db,lm.D9, CD = TRUE, DFB = TRUE, COVR = TRUE, pctg = 10, intersection = FALSE)
#'

offliers = function(dataset, mod, CD = TRUE, DFB = FALSE, COVR = FALSE, pctg = 100, intersection = TRUE){
  stopifnot(class(pctg)=="numeric" || class(pctg) == "integer" && pctg <=100 && pctg >=0)
  n <- nrow(dataset)
  outmax <- floor(n*pctg/100)
  if (outmax == 0) return("You have requested to remove 0 outliers")
  sqrt.n <- sqrt(n)
  model <- mod
  k <- length(model$coefficients)-1
  if (CD == TRUE){
    cd <- na.omit(cooks.distance(model))
    cutoff.cd <- 4/(n-k-1)
    outliers.cd <- cd > cutoff.cd
    outliers.cd <- which(outliers.cd == TRUE)
    num.outliers.cd <- length(outliers.cd)
    if (num.outliers.cd > 0){
      top.outliers.cd <- order(cd[outliers.cd], decreasing = TRUE)[1:(min(num.outliers.cd,outmax))]
      outliers.cd <- outliers.cd[top.outliers.cd]
      num.outliers.cd <- length(outliers.cd)
    }
    prop.outliers.cd <- num.outliers.cd/n
  }
  if (DFB == TRUE){
    outliers.dfb <- integer(n) #store memory for potentially n outliers
    dfb <- dfbetas(model)
    cutoff.dfb <- 2/sqrt.n
    for (i in 1:n){
      for (j in 1:(k+1)){
        if (abs(dfb[i,j]) > cutoff.dfb){
          outliers.dfb[i] <- i
          break #as soon as an observation fails the test for 1 of the variables, signal it as outlier
        }
      }
    }
    outliers.dfb <- outliers.dfb[outliers.dfb > 0]
    num.outliers.dfb <- length(outliers.dfb)
    if (num.outliers.dfb > 0){
      outliers.dfb.values <- rowSums(abs(dfb[outliers.dfb, ]), na.rm = T) #sum by rows in subset of dfb 2D array
      top.outliers.dfb <- order(outliers.dfb.values, decreasing = TRUE)[1:(min(num.outliers.dfb, outmax))]
      outliers.dfb <- outliers.dfb[top.outliers.dfb]
      num.outliers.dfb <- length(outliers.dfb)
    }
    prop.outliers.dfb <- num.outliers.dfb/n
  }
  if (COVR == TRUE){
    covr <- covratio(model)
    cutoff.covr <- 3*(k+1)/n
    outliers.covr <- abs(covr-1) > cutoff.covr #|covr-1| > 3(k+1)/n <-> covr > 1+3(k+1)/n or covr < 1-3(k+1)/n
    outliers.covr <- which(outliers.covr == TRUE)
    num.outliers.covr <- length(outliers.covr)
    if(num.outliers.covr > 0){
      top.outliers.covr <- order(covr[outliers.covr], decreasing = TRUE)[1:(min(num.outliers.covr,outmax))]
      outliers.covr <- outliers.covr[top.outliers.covr]
      num.outliers.covr <- length(outliers.covr)
    }
    prop.outliers.covr <- num.outliers.covr/n
  }
  if (intersection == FALSE){ #avoids code duplication
    variant <- union
  }else{
    variant <- intersect
  }

  # There are 3C1 + 3C2 + 3C3 = 7 mutually exclusive cases

  # Case 1
  if (CD == TRUE && DFB == TRUE && COVR == TRUE){
    final.outliers <- sort(Reduce(variant, list(outliers.cd, outliers.dfb, outliers.covr)))
    num.final.outliers <- length(final.outliers)
    prop.final.outliers <- num.final.outliers/n
    if (num.final.outliers > 0){
      data.offliers <- dataset[-final.outliers, ]
    }else{
      data.offliers <- dataset
    }
    solution <- list(new.dataset = data.offliers,
                     "Outliers according to Cook's Distance" = sort(as.vector(outliers.cd)),
                     "# outliers according to Cook's Distance" = num.outliers.cd,
                     "% of outliers according to Cook's Distance" = prop.outliers.cd*100,
                     "Outliers according to DFBetas" = sort(outliers.dfb),
                     "# outliers according to DFBetas" = num.outliers.dfb,
                     "% of outliers according to DFBetas" = prop.outliers.dfb*100,
                     "Outliers according to COVRATIO" = sort(as.vector(outliers.covr)),
                     "# outliers according to COVRATIO" = num.outliers.covr,
                     "% of outliers according to COVRATIO" = prop.outliers.covr*100,
                     "Final outliers identified" = sort(final.outliers),
                     "# final outliers" = num.final.outliers,
                     "% of final outliers" = prop.final.outliers*100)
    return(solution)
  }
  # Case 2
  if (CD == TRUE && DFB == TRUE && COVR == FALSE){
    final.outliers <- sort(Reduce(variant, list(outliers.cd, outliers.dfb)))
    num.final.outliers <- length(final.outliers)
    prop.final.outliers <- num.final.outliers/n
    if (num.final.outliers > 0){
      data.offliers <- dataset[-final.outliers, ]
    }else{
      data.offliers <- dataset
    }
    solution <- list(new.dataset = data.offliers,
                     "Outliers according to Cook's Distance" = sort(as.vector(outliers.cd)),
                     "# outliers according to Cook's Distance" = num.outliers.cd,
                     "% of outliers according to Cook's Distance" = prop.outliers.cd*100,
                     "Outliers according to DFBetas" = sort(outliers.dfb),
                     "# outliers according to DFBetas" = num.outliers.dfb,
                     "% of outliers according to DFBetas" = prop.outliers.dfb*100,
                     "Final outliers identified" = sort(final.outliers),
                     "# final outliers" = num.final.outliers,
                     "% of final outliers" = prop.final.outliers*100)
    return(solution)
  }
  # Case 3
  if (CD == TRUE && DFB == FALSE && COVR == TRUE){
    final.outliers <- sort(Reduce(variant, list(outliers.cd, outliers.covr)))
    num.final.outliers <- length(final.outliers)
    prop.final.outliers <- num.final.outliers/n
    if (num.final.outliers > 0){
      data.offliers <- dataset[-final.outliers, ]
    }else{
      data.offliers <- dataset
    }
    solution <- list(new.dataset = data.offliers,
                     "Outliers according to Cook's Distance" = sort(as.vector(outliers.cd)),
                     "# outliers according to Cook's Distance" = num.outliers.cd,
                     "% of outliers according to Cook's Distance" = prop.outliers.cd*100,
                     "Outliers according to COVRATIO" = sort(as.vector(outliers.covr)),
                     "# outliers according to COVRATIO" = num.outliers.covr,
                     "% of outliers according to COVRATIO" = prop.outliers.covr*100,
                     "Final outliers identified" = sort(final.outliers),
                     "# final outliers" = num.final.outliers,
                     "% of final outliers" = prop.final.outliers*100)
    return(solution)
  }
  # Case 4
  if (CD == FALSE && DFB == TRUE && COVR == TRUE){
    final.outliers <- sort(Reduce(variant, list(outliers.dfb, outliers.covr)))
    num.final.outliers <- length(final.outliers)
    prop.final.outliers <- num.final.outliers/n
    if (num.final.outliers > 0){
      data.offliers <- dataset[-final.outliers, ]
    }else{
      data.offliers <- dataset
    }
    solution <- list(new.dataset = data.offliers,
                     "Outliers according to DFBetas" = sort(outliers.dfb),
                     "# outliers according to DFBetas" = num.outliers.dfb,
                     "% of outliers according to DFBetas" = prop.outliers.dfb*100,
                     "Outliers according to COVRATIO" = sort(as.vector(outliers.covr)),
                     "# outliers according to COVRATIO" = num.outliers.covr,
                     "% of outliers according to COVRATIO" = prop.outliers.covr*100,
                     "Final outliers identified" = sort(final.outliers),
                     "# final outliers" = num.final.outliers,
                     "% of final outliers" = prop.final.outliers*100)
    return(solution)
  }
  # Case 5
  if (CD == TRUE && DFB == FALSE && COVR == FALSE){
    if (num.outliers.cd > 0){
      data.offliers <- dataset[-outliers.cd, ]
    }else{
      data.offliers <- dataset
    }
    solution <- list(new.dataset = data.offliers,
                     "Outliers according to Cook's Distance" = sort(as.vector(outliers.cd)),
                     "# outliers according to Cook's Distance" = num.outliers.cd,
                     "% of outliers according to Cook's Distance" = prop.outliers.cd*100)
    return(solution)
  }
  # Case 6
  if (CD == FALSE && DFB == TRUE && COVR == FALSE){
    if (num.outliers.dfb > 0){
      data.offliers <- dataset[-outliers.dfb, ]
    }else{
      data.offliers <- dataset
    }
    solution <- list(new.dataset = data.offliers,
                     "Outliers according to DFBetas" = sort(outliers.dfb),
                     "# outliers according to DFBetas" = num.outliers.dfb,
                     "% of outliers according to DFBetas" = prop.outliers.dfb*100)
    return(solution)
  }
  # Case 7
  if (CD == FALSE && DFB == FALSE && COVR == TRUE){
    if (num.outliers.covr > 0){
      data.offliers <- dataset[-outliers.covr, ]
    }else{
      data.offliers <- dataset
    }
    solution <- list(new.dataset = data.offliers,
                     "Outliers according to COVRATIO" = sort(as.vector(outliers.covr)),
                     "# outliers according to COVRATIO" = num.outliers.covr,
                     "% of outliers according to COVRATIO" = prop.outliers.covr*100)
    return(solution)
  }
}
