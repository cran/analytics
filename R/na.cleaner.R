#' @title Missing Value Imputation
#' @description Missing value imputation based on different methods. Can handle continuous and
#' categorical variables.
#' @details Each of the available methods in this function may be the best choice for a particular
#' dataset, but since it is impossible to know which one it is in each particular case, Mode "all"
#' might be a good, robust choice.
#' For categorical variables, the only mode implemented is knn, so argument Mode really refers
#' only to the continuous variables.
#' @param dataset a matrix or data frame. May have continuous and/or categorical variables.
#' @param t1 the threshold value in interval 0-1 beyond which a record is deemed as having a high
#' \% of NAs. Default: 0.5.
#' @param t2 the threshold value in interval 0-1 beyond which a variable is deemed as having a
#' high \% of NAs. Default: 0.5.
#' @param auto If TRUE (the default), it will eliminate those records and/or variables deemed as
#' having a high \% of NAs. If FALSE, one handpicks which records/variables will be deleted.
#' @param maxDel1 the proportion in interval 0-1 of records that can at most be deleted.
#' Default: 0.2.
#' @param maxDel2 the proportion in interval 0-1 of variables that can at most be deleted.
#' Default: 0.3.
#' @param Mode a string specifying the imputation method to be used, among "mean" (default),
#' "median", "mean&lm", "median&lm", "knn".
#' @param neigh the neighbours to be used in knn, both for continuous and categorical variables.
#' Default: interval 3-7.
#' For each value in neigh, knn is run, and then in the case of continuous variables, the outcome
#' of those runs
#' are averaged out. In the case of categorical variables, the imputed value is the most common
#' imputed value across runs.
#' @return the original dataset with imputed missing values.
#'
#' @author Albert Dorador
#' @export
#' @importFrom VIM kNN
#' @seealso
#' \code{\link[VIM]{kNN}}
#' \code{\link[analytics]{rowmean}}
#' @examples
#'
#' mtcars_mod <- mtcars
#' set.seed(1)
#' mtcars_mod <- as.data.frame(lapply(mtcars_mod, function(cc) cc[ sample(c(TRUE, NA),
#' prob = c(0.6, 0.4), size = length(cc), replace = TRUE) ]))
#' rownames(mtcars_mod) <- rownames(mtcars)
#'
#' # Compare methods
#' kNN_dt <- na.cleaner(dataset = mtcars_mod, Mode = "kNN")
#' mean_lm_dt <- na.cleaner(dataset = mtcars_mod, Mode = "mean&lm")
#' median_dt <- na.cleaner(dataset = mtcars_mod, Mode = "median")
#' all_dt <- na.cleaner(dataset = mtcars_mod, Mode = "all")
#' dev_kNN <- norm(as.matrix(mtcars[-c(4,6,8,13,18,20), -6])-as.matrix(kNN_dt))
#' dev_m_ml <- norm(as.matrix(mtcars[-c(4,6,8,13,18,20), -6])-as.matrix(mean_lm_dt))
#' dev_md <- norm(as.matrix(mtcars[-c(4,6,8,13,18,20), -6])-as.matrix(median_dt))
#' dev_all <- norm(as.matrix(mtcars[-c(4,6,8,13,18,20), -6])-as.matrix(all_dt))
#'
#' iris_mod <- iris
#' set.seed(5)
#' iris_mod <- as.data.frame(lapply(iris_mod, function(cc) cc[ sample(c(TRUE, NA),
#' prob = c(0.6, 0.4), size = length(cc), replace = TRUE) ]))
#' rownames(iris_mod) <- rownames(iris)
#' na.cleaner(dataset = iris_mod, neigh = 1, Mode = "all")
#'

na.cleaner <- function(dataset, t1 = 0.5, t2 = 0.5, auto = TRUE, maxDel1 = 0.2, maxDel2 = 0.3, Mode = "mean", neigh = 3:7){
  nrows0 <- nrow(dataset)
  ncols0 <- ncol(dataset)
  rownames(dataset) <- 1:nrows0
  colnames(dataset) <- 1:ncols0
  bNAs <- is.na(dataset)
  if (!any(bNAs))
    stop("There are no missing values in the original dataset provided.")

  if (is.null(rownames(dataset))){
    brow_names <- FALSE
  } else {
    brow_names <- TRUE
    row_names <- rownames(dataset)
  }

  if (is.null(colnames(dataset))){
    bcol_names <- FALSE
  } else {
    bcol_names <- TRUE
    col_names <- colnames(dataset)
  }

  NAs_by_obs <- rowSums(bNAs) # 1 by m
  NAs_by_obs_prop <- rowMeans(bNAs)
  NAs_by_var <- colSums(bNAs) # 1 by n
  NAs_by_var_prop <- colMeans(bNAs)

  info_list <- list("NAs by record" = NAs_by_obs,
                    "NAs by record, proportion" = NAs_by_obs_prop,
                    "NAs by variable" = NAs_by_var,
                    "NAs by variable, proportion" = NAs_by_var_prop)
  print(info_list)

  quick_NAs_by_obs_prop <- which(NAs_by_obs_prop > t1)
  unnamed_quick_NAs_by_obs_prop <- unname(quick_NAs_by_obs_prop)
  quick_NAs_by_variable_prop <- which(NAs_by_var_prop > t2)
  unnamed_quick_NAs_by_variable_prop <- unname(quick_NAs_by_variable_prop)

  quick_info_list <- list("Index of records with high % of NAs" = unnamed_quick_NAs_by_obs_prop,
                          "Index of Variables with high % of NAs" = unnamed_quick_NAs_by_variable_prop)

  print(quick_info_list)

  new_dataset <- dataset

  # If auto enabled, get rid of any rows and/or columns up to maxDel1 and/or maxDel2
  if (auto) {

    sorted_NAs_by_obs <- sort(x = quick_NAs_by_obs_prop, decreasing = TRUE)
    sorted_NAs_by_variable <- sort(x = quick_NAs_by_variable_prop, decreasing = TRUE)
    num_maxDel1 <- floor(maxDel1*nrows0)
    num_maxDel2 <- floor(maxDel2*ncols0)
    num_high_NAs_obs <- length(quick_NAs_by_obs_prop)
    del_rows <- NULL
    if (num_high_NAs_obs > 0){
      num_to_delete_1 <- min(num_high_NAs_obs, num_maxDel1)
      del_rows <- as.integer(names(sorted_NAs_by_obs[1:num_to_delete_1]))
      new_dataset <- new_dataset[-del_rows,]
    }
    num_high_NAs_var <- length(quick_NAs_by_variable_prop)
    del_cols <- NULL
    if (num_high_NAs_var > 0){
      num_to_delete_2 <- min(num_high_NAs_var, num_maxDel2)
      del_cols <- as.integer(names(sorted_NAs_by_variable[1:num_to_delete_2]))
      new_dataset <- new_dataset[,-del_cols,drop = FALSE]
    }
  } else {
    # If manual, first we ask the user if she wants to get rid of any rows or columns

    prompt <- readline("Please provide comma-separatedly the records you want to delete - insert 0 if none \n")

    del_rows <- NULL
    if (prompt != "0"){
      del_rows <- as.integer(strsplit(prompt, ",")[[1]])
      new_dataset <- new_dataset[-del_rows,]
    }

    prompt <- readline("Please provide comma-separatedly the variables you want to delete - insert 0 if none \n")

    del_cols <- NULL
    if (prompt != "0"){
      del_cols <- as.integer(strsplit(prompt, ",")[[1]])
      new_dataset <- new_dataset[,-del_cols,drop = FALSE]
    }
  }

  # After possibly sanitizing the dataset, let's input missing values according to Mode
  # and only if there are any NAs (left) in the dataset
  if (anyNA(new_dataset)){

    cont_var_idx <- which(!sapply(new_dataset, is.factor))
    there_are_cont <- length(cont_var_idx > 0)
    categ_var_idx <- which(sapply(new_dataset, is.factor))
    there_are_categ <- length(categ_var_idx > 0)

    # A) Continuous variables first #
    # ----------------------------- #

    if (there_are_cont){
      cont.new_dataset <- new_dataset[, cont_var_idx, drop = FALSE]

      bNAs <- is.na(cont.new_dataset) #compute again because numrows and/or cols may have changed

      Mode <- tolower(Mode)
      ncols <- ncol(cont.new_dataset)

      if (Mode == "all"){
        ALL <- TRUE
        Mode <- "mean"
      } else {ALL <- FALSE}

      if (Mode == "mean"){
        var_means <- colMeans(cont.new_dataset, na.rm = TRUE)
        cont.new_dataset_final <- cont.new_dataset
        for (j in 1:ncols){
          jcolNAs <- bNAs[,j]
          cont.new_dataset_final[,j][jcolNAs] <- var_means[j]
        }
        if (ALL){Mode <- "median"; cont.new_dataset_1 <- cont.new_dataset_final}
      }

      if (Mode == "median"){
        var_medians <- apply(X = cont.new_dataset, MARGIN = 2, FUN = median, na.rm = TRUE)
        cont.new_dataset_final <- cont.new_dataset
        for (j in 1:ncols){
          jcolNAs <- bNAs[,j]
          cont.new_dataset_final[,j][jcolNAs] <- var_medians[j]
        }
        if (ALL){Mode <- "mean&lm"; cont.new_dataset_2 <- cont.new_dataset_final}
      }
      if (Mode == "mean&lm"){
        var_means <- colMeans(cont.new_dataset, na.rm = TRUE)
        cont.new_dataset_final <- cont.new_dataset
        for (j in 1:ncols){
          rng <- j:ncols
          for (jj in rng){
            if (jj < ncols){
              jjcolNAs <- bNAs[,jj+1]
              cont.new_dataset_final[,jj+1][jjcolNAs] <- var_means[jj+1]
            }
          }
          jcolNAs <- bNAs[,j]
          y <- cont.new_dataset_final[!jcolNAs,j] # we want to predict column j (excluding NAs)...
          X <- as.matrix(cont.new_dataset_final[!jcolNAs,-j]) # ... with the rest of columns
          model <- lm(y ~ X)
          Xpred <- cont.new_dataset_final[jcolNAs,-j] #the records that have an NA in column j
          Xpred <- cbind(rep(1,sum(jcolNAs)), Xpred) #attach 1s for Beta0
          ypred <- as.matrix(Xpred) %*% as.matrix(model$coefficients)
          cont.new_dataset_final[,j][jcolNAs] <- ypred
          print(paste0("Adjusted R^2 for variable ", j, ": ",summary(model)$adj.r.squared))
        }
        if (ALL){Mode <- "median&lm"; cont.new_dataset_3 <- cont.new_dataset_final}
      }
      if (Mode == "median&lm"){
        var_medians <- apply(X = cont.new_dataset, MARGIN = 2, FUN = median, na.rm = TRUE)
        cont.new_dataset_final <- cont.new_dataset
        for (j in 1:ncols){
          rng <- j:ncols
          for (jj in rng){
            if (jj < ncols){
              jjcolNAs <- bNAs[,jj+1]
              cont.new_dataset_final[,jj+1][jjcolNAs] <- var_medians[jj+1]
            }
          }
          jcolNAs <- bNAs[,j]
          y <- cont.new_dataset_final[!jcolNAs,j] # we want to predict column j (excluding NAs)...
          X <- as.matrix(cont.new_dataset_final[!jcolNAs,-j]) # ... with the rest of columns
          model <- lm(y ~ X)
          Xpred <- cont.new_dataset_final[jcolNAs,-j] #the records that have an NA in column j
          Xpred <- cbind(rep(1,sum(jcolNAs)), Xpred) #attach 1s for Beta0
          ypred <- as.matrix(Xpred) %*% as.matrix(model$coefficients)
          cont.new_dataset_final[,j][jcolNAs] <- ypred
          print(paste0("Adjusted R^2 for variable ", j, ": ",summary(model)$adj.r.squared))
        }
        if (ALL){Mode <- "knn"; cont.new_dataset_4 <- cont.new_dataset_final}
      }
      if (Mode == "knn"){
        cont.new_dataset_final <- cont.new_dataset

        # if more than 1 run of kNN
        if (length(neigh) > 1){
          datalist <- list()
          counter <- 0
          for (i in neigh) {
            counter <- counter + 1
            datalist[[counter]] <- kNN(data = cont.new_dataset_final, k = i, imp_var = FALSE)
          }
          cont.Tot_kNN_dataset <- do.call(cbind, datalist)

          group <- rep(1:ncols, length(neigh))

          # Average across kNN repetitions using a different #neigh
          cont.new_dataset_final <- t(rowmean(t(cont.Tot_kNN_dataset), group = group))

        } else {
          cont.new_dataset_final <- kNN(data = cont.new_dataset_final, k = neigh, imp_var = FALSE)
        }

        if (length(del_rows) > 0){
          rownames(cont.new_dataset_final) <- (1:nrows0)[-del_rows]
        } else {
          rownames(cont.new_dataset_final) <- 1:nrows0
        }
        if (length(del_cols) > 0){
          colnames(cont.new_dataset_final) <- (1:ncols0)[-del_cols]
        } else {
          colnames(cont.new_dataset_final) <- 1:ncols
        }
        if (ALL){cont.new_dataset_5 <- cont.new_dataset_final}
      }

      if (ALL){
        cont.Total_new_dataset <- rbind(as.matrix(cont.new_dataset_1),
                                        as.matrix(cont.new_dataset_2),
                                        as.matrix(cont.new_dataset_3),
                                        as.matrix(cont.new_dataset_4),
                                        as.matrix(cont.new_dataset_5))
        cont.new_dataset_final <- rowmean(cont.Total_new_dataset)
      }
    }

    # B) Categorical variables next #
    # ----------------------------- #

    if (there_are_categ){
      if (there_are_cont) {
        cat.new_dataset <- new_dataset
        cat.new_dataset[,cont_var_idx] <- cont.new_dataset_final #We take all vars and AFTER imputing NAs
        colnames(cat.new_dataset) <- colnames(new_dataset)
      } else {
        cat.new_dataset <- new_dataset
      }

      # if more than 1 run of knn
      if (length(neigh) > 1){
        a.mode <- function(x) {
          ux <- unique(x)
          ux[which.max(tabulate(match(x, ux)))]
        }

        datalist <- list()

        counter <- 0
        for (i in neigh) {
          counter <- counter + 1
          datalist[[counter]] <- kNN(data = cat.new_dataset, variable = colnames(cat.new_dataset)[categ_var_idx], k = i, imp_var = FALSE)
        }
        cat.Tot_kNN_dataset <- do.call(cbind, datalist)

        categ_var_idx_2 <- which(sapply(cat.Tot_kNN_dataset, is.factor))
        only_cat_vars <- cat.Tot_kNN_dataset[,categ_var_idx_2,drop = FALSE]
        cat.ncols <- length(categ_var_idx)
        group <- rep(1:cat.ncols, length(neigh))

        # Find mode by variable, across kNN repetitions using a different #neigh
        cat.new_dataset_final <- t(aggregate(t(only_cat_vars), list(group), a.mode))
        cat.new_dataset_final <- as.data.frame(cat.new_dataset_final[-1,]) #remove 1st row 'cause is just title

      } else {
        cat.new_dataset_final <- kNN(data = cat.new_dataset, variable = colnames(cat.new_dataset)[categ_var_idx], k = neigh, imp_var = FALSE)
        categ_var_idx_2 <- which(sapply(cat.new_dataset_final, is.factor))
        cat.new_dataset_final <- cat.new_dataset_final[,categ_var_idx_2,drop = FALSE]
      }

    }

    if (there_are_cont & there_are_categ){

      # Join the continuous and categorical parts of the dataset in original order
      new_dataset_final <- cbind(cont.new_dataset_final, cat.new_dataset_final)
      new_dataset_final <- new_dataset_final[, c(cont_var_idx, categ_var_idx)]
    } else if (there_are_cont){
      new_dataset_final <- cont.new_dataset_final
    } else if (there_are_categ){
      new_dataset_final <- cat.new_dataset_final
    }

    new_dataset <- new_dataset_final # to make it compatible with a df without NAs

  }

  if (brow_names){
    if (length(del_rows) > 0){
      rownames(new_dataset) <- row_names[-del_rows]
    } else {
      rownames(new_dataset) <- row_names
    }

  } else {
    rownames(new_dataset) <- NULL
  }

  if (bcol_names){
    if (length(del_cols) > 0){
      colnames(new_dataset) <- col_names[-del_cols]
    } else {
      colnames(new_dataset) <- col_names
    }

  } else {
    colnames(new_dataset) <- NULL
  }

  return(new_dataset)
}
