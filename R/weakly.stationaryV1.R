#' @title Testing for Weak Stationarity in a Time Series
#' @description Performs a series of statistical tests aimed at detecting non-stationarity.
#' @details This function offers a great deal of customization: diverse significance levels,
#' multiple tests specialized in certain aspects of (weak) stationarity, as well as
#' handy predefined sets of parameters providing a more or less strict diagnostic:
#' "neutral", "strict" and "loose" modes. By including this possibility, the technical burden
#' on the user is made lighter. Mode "strict" includes two tests for constant mean
#' (basic & Mann-Kendall), two tests for constant variance (McLeod-Li & Breusch-Pagan tests),
#' the Priestley-Subba Rao (PSR) test for nonstationarity across time, and three tests for
#' weak dependence (Phillips-Perron, Augmented Dickey-Fuller, and KPSS tests), which test
#' weak stationarity if and only if the underlying data generating process is assumed to be an AR(p).
#' Mode "loose" just performs the basic test for constant mean (a linear model that includes a trend
#' whose statistical significance is determined using robust regression if the Durbin Watson
#' test detects serial correlation in the residuals),
#' and the Breusch-Pagan test (on the previous auxiliary linear model's residuals) for constant
#' variance. Mode "neutral" (the default) provides all the default parameter options.
#' Significance levels also differ across modes.
#' This function differentiates two significance levels: general (\code{signific_gen})
#' and specific to the Phillips-Perron and Augmented Dickey-Fuller tests (\code{signific_pp.df}).
#' In mode "strict", \code{signific_gen} is 0.1, and \code{signific_pp.df} is 0.01.
#' In mode "loose", \code{signific_gen} is 0.01, and \code{signific_pp.df} is irrelevant.
#' In mode "neutral", both significance levels are set to 0.05.
#'
#' @param tseries a 1-D or 2-D array. In the latter case, the time series to be evaluated
#' must be placed in the 2nd dimension (columns). If that's not your case, transpose it.
#' @param signific_gen significance level for all tests except Phillips-Perron and Augmented Dickey-Fuller.
#' @param signific_pp.df significance level for the Phillips-Perron and Augmented Dickey-Fuller tests.
#' @param MK if TRUE, the Mann-Kendall test for constant mean is executed, instead of a faster basic test. Default is FALSE.
#' @param BP if TRUE, the Breusch-Pagan test for constant (residual) variance is executed
#' on the residuals of an auxiliary linear model that includes a time variable, instead of
#' the McLeod-Li test. The default is TRUE.
#' @param PSR if TRUE, the Priestley-Subba Rao test for nonstationarity across time is
#' executed. The default is TRUE.
#' @param weak.dep if TRUE, then the Phillips-Perron, Augmented Dickey-Fuller, and
#' KPSS tests for weak stationarity (assuming an AR(p)) are performed.
#' @param mode one of "neutral", "strict", "loose". Case insenstive. The default is "neutral".
#' @return if a 1-D array is supplied, then a Boolean is returned indicating whether
#' the time series supplied is weakly stationary (TRUE) or not (FALSE).
#' If a 2-D array is supplied, then a vector of Booleans is returned indicating whether
#' each individual time series supplied is weakly stationary (TRUE) or not (FALSE).
#' @author Albert Dorador
#' @export
#' @import urca
#' @importFrom car durbinWatsonTest
#' @importFrom robust lmRob
#' @importFrom trend mk.test
#' @importFrom TSA McLeod.Li.test
#' @importFrom lmtest bptest
#' @importFrom fractal stationarity
#' @examples
#' x1 <- rnorm(1e2)
#' weakly.stationary(tseries = x1)
#' weakly.stationary(tseries = x1, signific_gen = 0.025)
#' weakly.stationary(tseries = x1, signific_pp.df = 0.1)
#' weakly.stationary(tseries = x1, PSR = FALSE)
#' weakly.stationary(tseries = x1, weak.dep = TRUE)
#' weakly.stationary(tseries = x1, mode = "strict")
#' weakly.stationary(tseries = x1, mode = "loose")
#'
#' require(stats)
#' set.seed(123)
#' x2 <- arima.sim(n = 1e2, list(ar = 0.4))
#' weakly.stationary(tseries = x2)
#' weakly.stationary(tseries = x2, signific_gen = 0.01)
#' weakly.stationary(tseries = x2, PSR = FALSE)
#' weakly.stationary(tseries = x2, weak.dep = TRUE)
#' weakly.stationary(tseries = x2, mode = "strict")
#' weakly.stationary(tseries = x2, mode = "loose")
#'

weakly.stationary <- function(tseries, signific_gen = 0.05, signific_pp.df = 0.05,
                              MK = FALSE, BP = TRUE, PSR = TRUE,
                              weak.dep = FALSE, mode = "neutral"){
  # Note:2) MK has bad time complexity 3) mode: "strict", "loose"
  if (!signific_gen %in% c(0.01, 0.025, 0.05, 0.1))
    stop("Signific_gen currently supports: 0.01, 0.025, 0.05, 0.1")

  if (!signific_pp.df %in% c(0.01, 0.05, 0.1))
    stop("Signific_pp.df currently supports: 0.01, 0.05, 0.1")

  mode <- tolower(mode) #make this parameter case-insensitive

  if (!mode %in% c("loose", "neutral", "strict"))
    stop("mode must be one of: 'loose', 'neutral', 'strict'.")


  if (mode == "strict") { #Trend(lm&MK)+Var(McLi&BP)+PSR+weak.dep
    signific_gen <- 0.1
    signific_pp.df <- 0.01
    weak.dep <- TRUE
  }
  if (mode == "loose") { #Trend(lm)+Var(BP)
    signific_gen <- 0.01
    PSR <- FALSE
  }

  Dim <- dim(tseries)

  if (length(Dim) == 2) {
    n.series <- Dim[2]
    time <- 1:Dim[1]

    #preallocate memory
    results.vec <- logical(n.series)

    for (j in 1:n.series) {
      series.j <- tseries[,j]

      # (1) Testing for trend

      mod <- lm(series.j ~ time)

      if (!MK && mode != "strict") {
        mod.res <- as.vector(mod$residuals)
        DW <- durbinWatsonTest(mod.res)
        if (DW > 1.85) { #Don't reject H0 of serial uncorrelation
          s.mod <- summary(mod)
        } else { #Residuals are serially correlated
          mod.R <- lmRob(series.j ~ time)
          s.mod <- summary(mod.R)
        }
        stat_in_trend <- !(s.mod$coefficients[8] < signific_gen) #0.01 is looser than 0.05

      } else if (MK && mode != "strict") { #if Mann-Kendall test is requested
        series.j <- as.ts(series.j)
        mk <- mk.test(series.j)
        stat_in_trend <- !(mk$pvalg < signific_gen) #0.01 is looser than 0.05

      } else { # mode is "strict"
        mod.res <- as.vector(mod$residuals)
        DW <- durbinWatsonTest(mod.res)
        if (DW > 1.85) { #Don't reject H0 of serial uncorrelation
          s.mod <- summary(mod)
        } else { #Residuals are serially correlated
          mod.R <- lmRob(series.j ~ time)
          s.mod <- summary(mod.R)
        }
        stat_in_trend.t <- !(s.mod$coefficients[8] < signific_gen)
        series.j <- as.ts(series.j)
        mk <- mk.test(series.j)
        stat_in_trend.MK <- !(mk$pvalg < signific_gen)

        stat_in_trend <- stat_in_trend.t && stat_in_trend.MK
      }

      if (!stat_in_trend) {
        results.vec[j] <- FALSE
        next #skip to next iteration
      }

      # (2) Testing for non-constant variance

      if (!BP && mode != "strict") {
        # a) McLeod-Li test (Ljung-Box test on squared data)
        McLi <- McLeod.Li.test(y = series.j, plot = FALSE) #0.01 is looser than 0.05
        mcli <- !(any(McLi$p.values < signific_gen)) #!(Non-constant variance)
        stat_in_var <- mcli

      } else if (BP && mode != "strict") {
        # b) Breusch-Pagan test
        bp <- bptest(mod) #0.01 is looser than 0.05
        bp <- !(bp$p.value[[1]] < signific_gen) #!(Non-constant variance (of the error))
        stat_in_var <- bp

      } else { # mode is "strict"
        # a) McLeod-Li test (Ljung-Box test on squared data)
        McLi <- McLeod.Li.test(y = series.j, plot = FALSE)
        mcli <- !(any(McLi$p.values < signific_gen)) #!(Non-constant variance)

        # b) Breusch-Pagan test
        bp <- bptest(mod)
        bp <- !(bp$p.value[[1]] < signific_gen) #!(Non-constant variance (of the error))

        stat_in_var <- mcli && bp
      }

      if (!stat_in_var) {
        results.vec[j] <- FALSE
        next
      }

      if (PSR && !weak.dep) {
        sink("NUL") #starts suppressing output until sink() is called
        psr.info <- stationarity(series.j, significance = signific_gen) #Ho: stat.
        sink()
        psr <- attributes(psr.info)$stationarity[[1]] #0.01 is looser than 0.05
        results.vec[j] <- psr
      }

      if (PSR && weak.dep) { #battery of weak dependence tests: PP, ADF, KPSS
        # PSR
        sink("NUL") #starts suppressing output until sink() is called
        psr.info <- stationarity(series.j, significance = signific_gen) #Ho: stat.
        sink()
        psr <- attributes(psr.info)$stationarity[[1]]

        # Weak dep
        # (a) Phillips-Perron
        a <- summary(ur.pp(series.j, type = "Z-tau", model = "constant"))
        z.tau <- attributes(a)$teststat[1] #0.05 is looser than 0.01
        pp.crit.val <- ifelse(signific_pp.df == 0.01, attributes(a)$cval[1],
                              ifelse(signific_pp.df == 0.05, attributes(a)$cval[2],
                                     attributes(a)$cval[3]))
        PP <- z.tau < pp.crit.val #reject H0 of unit root

        # (b) Augmented Dickey-Fuller
        found <- FALSE
        lag.order <- 0
        while (lag.order <= 20 && !found) {
          ADF.info <- ur.df(series.j, lags = lag.order)
          ADF.res <- as.vector(ADF.info@res)
          DW <- durbinWatsonTest(ADF.res)
          if (DW > 1.85) { #Verbeek 2004, p.185; don't reject H0
            found <- TRUE #stop
          } else {
            lag.order <- lag.order + 1
          }
        }
        b <- summary(ur.df(series.j, lags = lag.order))
        adf.stat <- attributes(b)$teststat[1] #0.05 is looser than 0.01
        adf.crit.val <- ifelse(signific_pp.df == 0.01, attributes(b)$cval[1],
                               ifelse(signific_pp.df == 0.05, attributes(b)$cval[2],
                                      attributes(b)$cval[3]))
        ADF <- adf.stat < adf.crit.val #reject H0 of unit root

        # (c) KPSS
        c <- summary(ur.kpss(series.j, type = "mu", use.lag = lag.order))
        kpss.stat <- attributes(c)$teststat[1] #0.01 is looser than 0.05
        kpss.crit.val <- ifelse(signific_gen == 0.1, attributes(c)$cval[1],
                                ifelse(signific_gen == 0.05, attributes(c)$cval[2],
                                       ifelse(signific_gen == 0.025, attributes(c)$cval[3],
                                              attributes(c)$cval[4])))
        KPSS <- kpss.stat < kpss.crit.val # Do NOT reject H0 of weakly dependent

        weakly.dependent <- PP && ADF && KPSS

        results.vec[j] <- psr && weakly.dependent
      }

      if (!PSR && weak.dep) {
        # (a) Phillips-Perron
        a <- summary(ur.pp(series.j, type = "Z-tau", model = "constant"))
        z.tau <- attributes(a)$teststat[1] #0.05 is looser than 0.01
        pp.crit.val <- ifelse(signific_pp.df == 0.01, attributes(a)$cval[1],
                              ifelse(signific_pp.df == 0.05, attributes(a)$cval[2],
                                     attributes(a)$cval[3]))
        PP <- z.tau < pp.crit.val #reject H0 of unit root

        # (b) Augmented Dickey-Fuller
        found <- FALSE
        lag.order <- 0
        while (lag.order <= 20 && !found) {
          ADF.info <- ur.df(series.j, lags = lag.order)
          ADF.res <- as.vector(ADF.info@res)
          DW <- durbinWatsonTest(ADF.res)
          if (DW > 1.85) { #Verbeek 2004, p.185; don't reject H0
            found <- TRUE #stop
          } else {
            lag.order <- lag.order + 1
          }
        }
        b <- summary(ur.df(series.j, lags = lag.order))
        adf.stat <- attributes(b)$teststat[1] #0.05 is looser than 0.01
        adf.crit.val <- ifelse(signific_pp.df == 0.01, attributes(b)$cval[1],
                               ifelse(signific_pp.df == 0.05, attributes(b)$cval[2],
                                      attributes(b)$cval[3]))
        ADF <- adf.stat < adf.crit.val #reject H0 of unit root

        # (c) KPSS
        c <- summary(ur.kpss(series.j, type = "mu", use.lag = lag.order))
        kpss.stat <- attributes(c)$teststat[1] #0.01 is looser than 0.05
        kpss.crit.val <- ifelse(signific_gen == 0.1, attributes(c)$cval[1],
                                ifelse(signific_gen == 0.05, attributes(c)$cval[2],
                                       ifelse(signific_gen == 0.025, attributes(c)$cval[3],
                                              attributes(c)$cval[4])))
        KPSS <- kpss.stat < kpss.crit.val # Do NOT reject H0 of weakly dependent

        weakly.dependent <- PP && ADF && KPSS

        results.vec[j] <- weakly.dependent
      }

      if (!PSR && !weak.dep) {
        results.vec[j] <- TRUE #if reached this far, must have passed trend & var
      }
    }
    return(results.vec)

  } else { #only 1 time series

    # (1) Testing for trend

    time <- 1:length(tseries)
    mod <- lm(tseries ~ time)

    if (!MK && mode != "strict") {
      mod.res <- as.vector(mod$residuals)
      DW <- durbinWatsonTest(mod.res)
      if (DW > 1.85) { #Don't reject H0 of serial uncorrelation
        s.mod <- summary(mod)
      } else { #Residuals are serially correlated
        mod.R <- lmRob(tseries ~ time)
        s.mod <- summary(mod.R)
      }
      stat_in_trend <- !(s.mod$coefficients[8] < signific_gen) #0.01 is looser than 0.05

    } else if (MK && mode != "strict") { #if Mann-Kendall test is requested
      tseries <- as.ts(tseries)
      mk <- mk.test(tseries)
      stat_in_trend <- !(mk$pvalg < signific_gen) #0.01 is looser than 0.05

    } else { # mode is "strict"
      mod.res <- as.vector(mod$residuals)
      DW <- durbinWatsonTest(mod.res)
      if (DW > 1.85) { #Don't reject H0 of serial uncorrelation
        s.mod <- summary(mod)
      } else { #Residuals are serially correlated
        mod.R <- lmRob(tseries ~ time)
        s.mod <- summary(mod.R)
      }
      stat_in_trend.t <- !(s.mod$coefficients[8] < signific_gen)
      tseries <- as.ts(tseries)
      mk <- mk.test(tseries)
      stat_in_trend.MK <- !(mk$pvalg < signific_gen)

      stat_in_trend <- stat_in_trend.t && stat_in_trend.MK
    }

    if (!stat_in_trend) {
      return(FALSE)
    }

    # (2) Testing for non-constant variance

    if (!BP && mode != "strict") {
      # a) McLeod-Li test (Ljung-Box test on squared data)
      McLi <- McLeod.Li.test(y = tseries, plot = FALSE) #0.01 is looser than 0.05
      mcli <- !(any(McLi$p.values < signific_gen)) #!(Non-constant variance)
      stat_in_var <- mcli

    } else if (BP && mode != "strict") {
      # b) Breusch-Pagan test
      bp <- bptest(mod) #0.01 is looser than 0.05
      bp <- !(bp$p.value[[1]] < signific_gen) #!(Non-constant variance (of the error))
      stat_in_var <- bp

    } else { # mode is "strict"
      # a) McLeod-Li test (Ljung-Box test on squared data)
      McLi <- McLeod.Li.test(y = tseries, plot = FALSE)
      mcli <- !(any(McLi$p.values < signific_gen)) #!(Non-constant variance)

      # b) Breusch-Pagan test
      bp <- bptest(mod)
      bp <- !(bp$p.value[[1]] < signific_gen) #!(Non-constant variance (of the error))

      stat_in_var <- mcli && bp
    }

    if (!stat_in_var) {
      return(FALSE)
    }

    if (PSR && !weak.dep) {
      sink("NUL") #starts suppressing output until sink() is called
      psr.info <- stationarity(tseries, significance = signific_gen) #Ho: stat.
      sink()
      psr <- attributes(psr.info)$stationarity[[1]] #0.01 is looser than 0.05
      return(psr)
    }

    if (PSR && weak.dep) { #battery of weak dependence tests: PP, ADF, KPSS
      # PSR
      sink("NUL") #starts suppressing output until sink() is called
      psr.info <- stationarity(tseries, significance = signific_gen) #Ho: stat.
      sink()
      psr <- attributes(psr.info)$stationarity[[1]]

      # Weak dep
      # (a) Phillips-Perron
      a <- summary(ur.pp(tseries, type = "Z-tau", model = "constant"))
      z.tau <- attributes(a)$teststat[1] #0.05 is looser than 0.01
      pp.crit.val <- ifelse(signific_pp.df == 0.01, attributes(a)$cval[1],
                            ifelse(signific_pp.df == 0.05, attributes(a)$cval[2],
                                   attributes(a)$cval[3]))
      PP <- z.tau < pp.crit.val #reject H0 of unit root

      # (b) Augmented Dickey-Fuller
      found <- FALSE
      lag.order <- 0
      while (lag.order <= 20 && !found) {
        ADF.info <- ur.df(tseries, lags = lag.order)
        ADF.res <- as.vector(ADF.info@res)
        DW <- durbinWatsonTest(ADF.res)
        if (DW > 1.85) { #Verbeek 2004, p.185; don't reject H0
          found <- TRUE #stop
        } else {
          lag.order <- lag.order + 1
        }
      }
      b <- summary(ur.df(tseries, lags = lag.order))
      adf.stat <- attributes(b)$teststat[1] #0.05 is looser than 0.01
      adf.crit.val <- ifelse(signific_pp.df == 0.01, attributes(b)$cval[1],
                             ifelse(signific_pp.df == 0.05, attributes(b)$cval[2],
                                    attributes(b)$cval[3]))
      ADF <- adf.stat < adf.crit.val #reject H0 of unit root

      # (c) KPSS
      c <- summary(ur.kpss(tseries, type = "mu", use.lag = lag.order))
      kpss.stat <- attributes(c)$teststat[1] #0.01 is looser than 0.05
      kpss.crit.val <- ifelse(signific_gen == 0.1, attributes(c)$cval[1],
                              ifelse(signific_gen == 0.05, attributes(c)$cval[2],
                                     ifelse(signific_gen == 0.025, attributes(c)$cval[3],
                                            attributes(c)$cval[4])))
      KPSS <- kpss.stat < kpss.crit.val # Do NOT reject H0 of weakly dependent

      weakly.dependent <- PP && ADF && KPSS

      return(psr && weakly.dependent)
    }

    if (!PSR && weak.dep) {
      # (a) Phillips-Perron
      a <- summary(ur.pp(tseries, type = "Z-tau", model = "constant"))
      z.tau <- attributes(a)$teststat[1] #0.05 is looser than 0.01
      pp.crit.val <- ifelse(signific_pp.df == 0.01, attributes(a)$cval[1],
                            ifelse(signific_pp.df == 0.05, attributes(a)$cval[2],
                                   attributes(a)$cval[3]))
      PP <- z.tau < pp.crit.val #reject H0 of unit root

      # (b) Augmented Dickey-Fuller
      found <- FALSE
      lag.order <- 0
      while (lag.order <= 20 && !found) {
        ADF.info <- ur.df(tseries, lags = lag.order)
        ADF.res <- as.vector(ADF.info@res)
        DW <- durbinWatsonTest(ADF.res)
        if (DW > 1.85) { #Verbeek 2004, p.185; don't reject H0
          found <- TRUE #stop
        } else {
          lag.order <- lag.order + 1
        }
      }
      b <- summary(ur.df(tseries, lags = lag.order))
      adf.stat <- attributes(b)$teststat[1] #0.05 is looser than 0.01
      adf.crit.val <- ifelse(signific_pp.df == 0.01, attributes(b)$cval[1],
                             ifelse(signific_pp.df == 0.05, attributes(b)$cval[2],
                                    attributes(b)$cval[3]))
      ADF <- adf.stat < adf.crit.val #reject H0 of unit root

      # (c) KPSS
      c <- summary(ur.kpss(tseries, type = "mu", use.lag = lag.order))
      kpss.stat <- attributes(c)$teststat[1] #0.01 is looser than 0.05
      kpss.crit.val <- ifelse(signific_gen == 0.1, attributes(c)$cval[1],
                              ifelse(signific_gen == 0.05, attributes(c)$cval[2],
                                     ifelse(signific_gen == 0.025, attributes(c)$cval[3],
                                            attributes(c)$cval[4])))
      KPSS <- kpss.stat < kpss.crit.val # Do NOT reject H0 of weakly dependent

      weakly.dependent <- PP && ADF && KPSS

      return(weakly.dependent)
    }

    if (!PSR && !weak.dep) {
      return(TRUE) #if reached this far, must have passed trend & var
    }
  }
}
