
# (k) = exp(0.095*Temp − 0.00014 × Temp^2) × (1 − exp[−1.21 × Precip]) from Yasso07 model
# (see also Brian’s Nature paper for more info)
# NEW: fixed so that map can be in units of mm, rather than m.
yasso_k <- function(mat, map) {
  exp(0.095 * mat - 0.00014 * mat^2) * (1 - exp(-1.21 * map/1000))
}

includeQuadraticTerms <- function(x, exclude=NULL) {
  for(p in colnames(x)) {
    if(!p %in% exclude) x[,paste0(p, "_2")] <- x[,p]^2
  }
  x
}

includeInteractions <- function(x, exclude=NULL) {
  f <- as.formula(~. * .)
  y <- x[,which(!colnames(x) %in% exclude)]
  cbind(x, as.data.frame(model.matrix.lm(f, y, na.action = "na.pass")[,-(1:ncol(x))]))
}

# "exclude": character vector of terms NOT to expand
expandTerms <- function(x, exclude=NULL, quadratic=TRUE, interactions=TRUE) {
  if(!quadratic && !interactions) {
    return(x)
  }
  if(!quadratic && interactions) {
    return(includeInteractions(x, exclude))
  }
  if(quadratic && !interactions) {
    return(includeQuadraticTerms(x, exclude))
  }
  if(quadratic && interactions) {
    return(cbind(
      includeInteractions(x, exclude=exclude),
      includeQuadraticTerms(x, exclude=exclude)[,-(1:ncol(x))]
    ))
  }
}

scaleToReference <- function(x, ref) {
  if(ncol(x) != ncol(ref)) {
    stop("Input data frame and 'ref' must have the same number of columns.")
  }
  means <- apply(ref, 2, mean, na.rm=TRUE)
  sds <- apply(ref, 2, sd, na.rm=TRUE)
  sweep(sweep(x, 2, means, FUN = "-"), 2, sds, FUN = "/")
}

scaleAndExpandPredictors <- function(surface, physeq, models=NULL, exclude="k") {
  pred_vars <- c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality", "k")
  # Expand predictors to match those of models, if provided
  if(!is.null(models)) {
    covars <- rownames(coef(models[[1]]))
    if(any(grepl(":", covars))) inclInteractions <- TRUE else inclInteractions <- FALSE
    if(any(grepl("_2$", covars))) inclQuadraticTerms <- TRUE else inclQuadraticTerms <- FALSE
    message("Based on covariate names of the first model, assuming that inclQuadraticTerms = ", inclQuadraticTerms,
            ", and inclInteractions = ", inclInteractions)
    if(inclQuadraticTerms && inclInteractions) {
      surface <- expandTerms(surface, exclude=exclude)
      predictors <- expandTerms(as(sample_data(physeq)[,pred_vars], "data.frame"), exclude=exclude)
    } else if(inclQuadraticTerms && !inclInteractions) {
      surface <- includeQuadraticTerms(surface, exclude=exclude)
      predictors <- includeQuadraticTerms(as(sample_data(physeq)[,pred_vars], "data.frame"), exclude=exclude)
    } else if(!inclQuadraticTerms && inclInteractions) {
      surface <- includeInteractions(surface, exclude=exclude)
      predictors <- includeInteractions(as(sample_data(physeq)[,pred_vars], "data.frame"), exclude=exclude)
    } else {
      surface <- surface
      predictors <- as(sample_data(physeq)[,pred_vars], "data.frame")
    }
    # If models are not provided...
  } else {
    warning("No models provided. Assuming that all quadratic terms and interactions should be included.")
    surface <- expandTerms(surface, exclude=exclude)
    predictors <- expandTerms(as(sample_data(physeq)[,pred_vars], "data.frame"), exclude=exclude)
  }
  surface_scaled <- scaleToReference(surface, predictors)
  return(surface_scaled)
}


# getNicheMetrics2D0 <- function(physeq, models, thresholds, x1, x2, r_climate,
#                               n_bins_per_axis=100, constant_values=NULL,
#                               hull=c("2D", "4D", "none"), inflate=0, ncores=1,
#                               progressbar=TRUE, dev_version=FALSE,
#                               exclude="k",
#                               pred_vars = c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality", "k")) {
#   # pred_vars <- c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality", "k")
#   if(!x1 %in% pred_vars | !x2 %in% pred_vars) {
#     stop("`x1` and `x2` must be among `pred_vars`")
#   }
#   hull <- match.arg(hull, c("2D", "4D", "none"))
#   n_gradient_steps <- n_bins_per_axis + 1
#   n_taxa <- length(taxa_names(physeq))
#   obs_env <- as(sample_data(physeq), "data.frame")[, pred_vars]
#   obs_env <- obs_env[complete.cases(obs_env),]
#   env_medians <- apply(obs_env, 2, function(x) median(x, na.rm=TRUE))
#   if(is.null(constant_values)) {
#     env_constant <- env_medians
#   } else {
#     if(length(constant_values) != length(pred_vars)) {
#       stop(paste0("`constant_values` must be length of ", length(pred_vars)))
#     }
#     env_constant <- constant_values
#     env_constant[which(is.na(constant_values))] <- env_medians[which(is.na(constant_values))]
#     # env_constant <- c(env_constant, yasso_k(env_constant[1], env_constant[3]))
#     names(env_constant) <- pred_vars
#   }
#   message("Using as constant unless otherwise specified: ", paste0(pred_vars, " ", env_constant, ", "))
#
#   quantize <- function(x, n) {
#     seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length.out=n)
#   }
#   surface_var <- expand.grid(
#     quantize(obs_env[,x1], n_gradient_steps),
#     quantize(obs_env[,x2], n_gradient_steps)
#   )
#   names(surface_var) <- c(x1, x2)
#   for(p in pred_vars) {
#     if(!p %in% colnames(surface_var)) {
#       surface_var[,p] <- env_constant[p]
#     }
#   }
#   surface_var[,"k"] <- yasso_k(surface_var[,"mat_celsius"], surface_var[,"map_mm"])
#   surface_var <- dplyr::select(surface_var, all_of(pred_vars))
#
#   # covars <- rownames(coef(models[[1]]))
#   # if(any(grepl(":", covars))) inclInteractions <- TRUE else inclInteractions <- FALSE
#   # if(any(grepl("_2$", covars))) inclQuadraticTerms <- TRUE else inclQuadraticTerms <- FALSE
#   # message("Based on covariate names of the first model, assuming that inclQuadraticTerms = ", inclQuadraticTerms,
#   #         ", and inclInteractions = ", inclInteractions)
#   #
#   # if(inclQuadraticTerms && inclInteractions) {
#   #   surface_var <- expandTerms(surface_var, exclude=exclude)
#   #   predictors <- expandTerms(as(sample_data(physeq)[,pred_vars], "data.frame"), exclude=exclude)
#   # } else if(inclQuadraticTerms && !inclInteractions) {
#   #   surface_var <- includeQuadraticTerms(surface_var, exclude=exclude)
#   #   predictors <- includeQuadraticTerms(as(sample_data(physeq)[,pred_vars], "data.frame"), exclude=exclude)
#   # } else if(!inclQuadraticTerms && inclInteractions) {
#   #   surface_var <- includeInteractions(surface_var, exclude=exclude)
#   #   predictors <- includeInteractions(as(sample_data(physeq)[,pred_vars], "data.frame"), exclude=exclude)
#   # } else {
#   #   surface_var <- surface_var
#   #   predictors <- as(sample_data(physeq)[,pred_vars], "data.frame")
#   # }
#   # surface_var_scaled <- scaleToReference(surface_var, predictors)
#
#   # covars <- rownames(coef(models[[1]]))[-1]
#   # surface_var <- expandTerms(surface_var, "data.frame")
#   # predictors <- as(sample_data(physeq)[,covars], "data.frame")
#   # keep_cols <- which(colnames(surface_var) %in% covars)
#   # surface_var <- surface_var[,keep_cols]
#   # surface_var_scaled <- scaleToReference(surface_var, predictors)
#
#   covars <- covarNamesFromLognets(models)
#   surface_var <- includeQuadraticTerms(surface_var)
#   surface_var <- surface_var[,match(covars,colnames(surface_var))]
#   predictors <- as(sample_data(physeq)[,covars], "data.frame")
#   predictors <- predictors[,match(covars, colnames(predictors))]
#   surface_var_scaled <- scaleToReference(surface_var, predictors)
#
#   registerDoParallel(ncores)
#   surface_suitability <- foreach(i = seq_along(taxa_names(physeq)), .combine=cbind) %dopar% {
#     tryCatch(
#       as.vector(predict(models[[i]], newx = as.matrix(surface_var_scaled), type = "response", s = "lambda.min")),
#       error = function(e) return(NA))
#   }
#   stopImplicitCluster()
#   surface_presence <- do.call(cbind, lapply(seq_along(thresholds),
#                                             function(i) surface_suitability[,i] >= thresholds[[i]]))
#
#   colnames(surface_suitability) <- colnames(surface_presence) <- taxa_names(physeq)
#
#   inhull_list <- list()
#
#   # Options for constraining by a convex hull
#
#   # If constraining by the convex hull around observations for each pair of climate axes
#   if(identical(hull, "4D")) {
#     for(i in 1:3) {
#       for(j in (i+1):4) {
#         # cvhull <- obs_env[chull(obs_env[,i], obs_env[,j]),]
#         # inhull_list <- c(inhull_list, list(inpolygon(surface_var[,i], surface_var[,j], cvhull[,i], cvhull[,j], boundary=TRUE)))
#         x <- obs_env[,i]
#         y <- obs_env[,j]
#         ch_ind <- chull(x, y)
#         infl_x <- x[ch_ind] + inflate * (x[ch_ind] - mean(x))
#         infl_y <- y[ch_ind] + inflate * (y[ch_ind] - mean(y))
#         inhull_list <- c(inhull_list, list(inpolygon(surface_var[,i], surface_var[,j], infl_x, infl_y, boundary=TRUE)))
#       }
#     }
#     inhull <- rowSums(do.call(cbind, inhull_list)) == 6
#   }
#
#   # If constraining by the convex hull around observations for the two climate axes of interest
#   if(identical(hull, "2D")) {
#     x <- obs_env[,x1]
#     y <- obs_env[,x2]
#     ch_ind <- chull(x, y)
#     infl_x <- x[ch_ind] + inflate * (x[ch_ind] - mean(x))
#     infl_y <- y[ch_ind] + inflate * (y[ch_ind] - mean(y))
#     inhull <- inpolygon(surface_var[,x1], surface_var[,x2], infl_x, infl_y, boundary=TRUE)
#   }
#
#   # Otherwise, inhull is NA
#   if(identical(hull, "none")) {
#     inhull <- NA
#   }
#
#   # Break here for dev version
#   if(dev_version) {
#     message("Returning dev version. Only suitability, n_niches, and inhull returned.")
#     return(cbind(surface_var,
#                  mean_suitability = rowMeans(surface_suitability, na.rm=TRUE),
#                  n_niches = rowSums(surface_presence, na.rm=TRUE),
#                  inhull = inhull))
#   }
#
#   surface_ordinal <- expand.grid(
#     seq_len(n_gradient_steps),
#     seq_len(n_gradient_steps)
#   )
#   surface_ordinal_vec <- apply(surface_ordinal, 1, paste, collapse='.')
#
#   surface_n_limits <- rep(NA, nrow(surface_var))
#   if(progressbar) {
#     message("Progress for step 2 of 2:")
#     pb <- txtProgressBar(max=n_gradient_steps^2, style=3)
#   }
#   for(i1 in seq_len(n_gradient_steps)) {
#     for(i2 in seq_len(n_gradient_steps)) {
#       ticker <- (i1-1)*n_gradient_steps + i2
#       this_ord <- c(i1, i2)
#       this_ind <- which(surface_ordinal_vec == paste0(this_ord, collapse="."))
#       this_val <- surface_presence[this_ind,]
#       if(!identical(hull, "none")) { # If a hull is used,
#         if(!inhull[this_ind]) {      # and this cell isn't within the hull,
#           surface_n_limits[this_ind] <- NA # then skip it.
#           if(progressbar) setTxtProgressBar(pb, ticker)
#           next
#         }
#       }
#       offsets <- rbind(
#         c(-1, 0),
#         c(0, -1),
#         c(1, 0),
#         c(0, 1)
#       )
#       neighbors_ord <- sweep(offsets, 2, c(i1, i2), '+')
#       neighbors_ord <- neighbors_ord[which(apply(neighbors_ord, 1, function(x) all(x %in% seq_len(n_gradient_steps)))),]
#       neighbors_vec <- apply(neighbors_ord, 1, paste, collapse='.')
#       neighbors_ind <- which(surface_ordinal_vec %in% neighbors_vec)
#       neighbors_val <- foreach(k = neighbors_ind) %do% {
#         return(surface_presence[k,])
#       }
#       crossings <- lapply(neighbors_val, function(x) x==0 & this_val==1) # "inner" boundary
#       isboundary <- rowSums(do.call(cbind, crossings)) > 0
#       surface_n_limits[this_ind] <- sum(isboundary, na.rm=TRUE)
#       if(progressbar) setTxtProgressBar(pb, ticker)
#     }
#   }
#   if(progressbar) close(pb)
#
#   # Get the mean suitability and total number of niches within each climate voxel
#   surface_mean_suitability <- rowMeans(surface_suitability, na.rm=TRUE)
#   surface_n_niches <- rowSums(surface_presence, na.rm=TRUE)
#
#   # # Bring in climate raster
#   # df_climate <- as.matrix(r_climate)[,pred_vars]
#   # int <- cbind(
#   #   findInterval(df_climate[,1], quantize(obs_env[,1], n_gradient_steps)),
#   #   findInterval(df_climate[,2], quantize(obs_env[,2], n_gradient_steps)),
#   #   findInterval(df_climate[,3], quantize(obs_env[,3], n_gradient_steps)),
#   #   findInterval(df_climate[,4], quantize(obs_env[,4], n_gradient_steps))
#   # )
#   # int_vec <- apply(int, 1, paste0, collapse=".")
#   # int_vec[grep("NA", int_vec)] <- NA
#   #
#   # # Get the total area within each climate voxel
#   # area <- as.vector(area(r_climate))
#   # int_area <- cbind.data.frame(int = int_vec, area) %>%
#   #   group_by(int) %>%
#   #   summarise(area = sum(area)) %>%
#   #   mutate(int = as.character(int))
#   # surface_area <- int_area$area[match(surface_ordinal_vec, int_area$int)]
#   # surface_area[which(is.na(surface_area))] <- 0
#
#   # Calculate Sorensen-like climate sensitivity
#   surface_sorensen <- surface_n_limits / surface_n_niches
#   # # As raster map:
#   # r_sorensen <- raster(r_climate)
#   # r_sorensen[] <- surface_sorensen[match(int_vec, surface_ordinal_vec)]
#
#   # # Since metrics are only defined within intervals, discard the ends.
#   # # In this case, only the upper ends remain to be discarded.
#   # discard_surface <- which(rowSums(apply(surface_ordinal, 2, function(x) x == n_gradient_steps)) > 0)
#   # surface_area[discard_surface] <- 0
#
#   # Collect variables
#   # out_df <- cbind(surface_var, mean_suitability = surface_mean_suitability, n_niches = surface_n_niches,
#   #                 n_limits = surface_n_limits, sorensen = surface_sorensen, area = surface_area, inhull = inhull)
#   out_df <- cbind(surface_var, mean_suitability = surface_mean_suitability, n_niches = surface_n_niches,
#                   n_limits = surface_n_limits, sorensen = surface_sorensen, inhull = inhull)
#   return(out_df)
# }




# n_taxa <- 100
# spp_subset_ind <- sort(sample(1:length(taxa_names(neon_dob_prevalent)), size=n_taxa))
# neon_dob_prevalent_sppsub <- prune_taxa(taxa_names(neon_dob_prevalent)[spp_subset_ind],
#                                         neon_dob_prevalent)
# physeq <- permuteClimate(neon_dob_prevalent_sppsub, randomseed=1010101)
# models <-  foreach(i = seq_along(spp_subset_ind)) %dopar% {
#   y01 <- as.vector(otu_table(physeq)[,i]) > 0
#   set.seed(1010101)
#   tryCatch(cv.glmnet(as.matrix(predictors_train_std2)[complete_records2,],
#                      y01[complete_records2], family = "binomial",
#                      alpha = 0, type.measure = "default", standardize=FALSE),
#            error = function(e) return(NA))
#   }
# thresholds <- findThresholds(
#   physeq, models,
#   pred_vars = terms2)
# x1 <- "mat_celsius"
# x2 <- "map_mm"
# r_climate <- r_present_northam
# n_bins_per_axis <- 10
# hull <- "2D"
# inflate <- 0
# ncores <- 8
# progressbar <- TRUE
# exclude <- "k"
# pred_vars <- c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality",
#                "mat_celsius_2", "map_mm_2", "soilInCaClpH", "nitrogenPercent", "organicCPercent")

# Not sure if this works yet. It's supposed to reduce runtime by only
# getting lognet predictions for the climate cells that are within the
# hull, but in a test, it didn't reduce time nor did it return the
# same output.

getNicheMetrics2D2 <- function(physeq, models, thresholds, x1, x2,
                              x1_range=NULL, x2_range=NULL, constant_values=NULL,
                              r_climate, n_bins_per_axis=50,
                              hull=c("2D", "4D", "none"), inflate=0, ncores=1,
                              progressbar=TRUE, dev_version=FALSE,
                              pred_vars = c("mat_celsius", "temp_seasonality", "map_mm", "mat_celsius_2", "map_mm_2", "soilInCaClpH", "organicCPercent")) {

  if(!x1 %in% pred_vars | !x2 %in% pred_vars) {
    stop("`x1` and `x2` must be among `pred_vars`")
  }
  hull <- match.arg(hull, c("2D", "4D", "none"))
  n_gradient_steps <- n_bins_per_axis + 1
  n_taxa <- length(taxa_names(physeq))
  obs_env <- as(sample_data(physeq), "data.frame")[, pred_vars]
  obs_env <- obs_env[complete.cases(obs_env),]

  env_medians <- apply(obs_env, 2, function(x) median(x, na.rm=TRUE))
  if(is.null(constant_values)) {
    env_constant <- env_medians
  } else {
    if(length(constant_values) != length(pred_vars)) {
      stop("`constant_values` must have same length as `pred_vars`.
           To use the median of observed values for a variable, enter NA for
           that variable.")
    }
    env_constant <- constant_values
    env_constant[which(is.na(constant_values))] <- env_medians[which(is.na(constant_values))]
    names(env_constant) <- pred_vars
  }
  message("Using as constant unless otherwise specified: ", paste0(pred_vars, " ", env_constant, ", "))

  if(is.null(x1_range)) {
    x1_range <- range(obs_env[,x1], na.rm=TRUE)
    message("Using as `x1_range`: ", paste0(x1_range, collapse=", "))
  }
  if(is.null(x2_range)) {
    x2_range <- range(obs_env[,x2], na.rm=TRUE)
    message("Using as `x2_range`: ", paste0(x2_range, collapse=", "))
  }

  quantize <- function(x, n) {
    seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length.out=n)
  }
  surface_var <- tidyr::expand_grid(
    # quantize(obs_env[,x1], n_gradient_steps),
    # quantize(obs_env[,x2], n_gradient_steps)
    quantize(x1_range, n_gradient_steps),
    quantize(x2_range, n_gradient_steps)
  )
  names(surface_var) <- c(x1, x2)
  for(p in pred_vars) {
    if(!p %in% colnames(surface_var)) {
      surface_var[,p] <- env_constant[p]
    }
  }

  covars <- covarNamesFromLognets(models)
  surface_var <- includeQuadraticTerms(surface_var)
  surface_var <- surface_var[,match(covars,colnames(surface_var))]
  predictors <- as(sample_data(physeq)[,covars], "data.frame")
  predictors <- predictors[,match(covars, colnames(predictors))]
  surface_var_scaled <- scaleToReference(surface_var, predictors)

  # Options for constraining by a convex hull

  # # If constraining by the convex hull around observations for each pair of climate axes
  # inhull_list <- list()
  # if(identical(hull, "4D")) {
  #   for(i in 1:3) {
  #     for(j in (i+1):4) {
  #       # cvhull <- obs_env[chull(obs_env[,i], obs_env[,j]),]
  #       # inhull_list <- c(inhull_list, list(inpolygon(surface_var[,i], surface_var[,j], cvhull[,i], cvhull[,j], boundary=TRUE)))
  #       x <- obs_env[,i]
  #       y <- obs_env[,j]
  #       ch_ind <- chull(x, y)
  #       infl_x <- x[ch_ind] + inflate * (x[ch_ind] - mean(x))
  #       infl_y <- y[ch_ind] + inflate * (y[ch_ind] - mean(y))
  #       inhull_list <- c(inhull_list, list(inpolygon(surface_var[,i], surface_var[,j], infl_x, infl_y, boundary=TRUE)))
  #     }
  #   }
  #   inhull <- rowSums(do.call(cbind, inhull_list)) == 6
  # }

  # If constraining by the convex hull around observations for the two climate axes of interest
  if(identical(hull, "2D")) {
    x <- obs_env[,x1]
    y <- obs_env[,x2]
    ch_ind <- chull(x, y)
    infl_x <- x[ch_ind] + inflate * (x[ch_ind] - mean(x))
    infl_y <- y[ch_ind] + inflate * (y[ch_ind] - mean(y))
    inhull <- pracma::inpolygon(surface_var[,x1,drop=TRUE], surface_var[,x2,drop=TRUE], infl_x, infl_y, boundary=TRUE)
  }

  # Otherwise, inhull is NA
  if(identical(hull, "none")) {
    inhull <- NA
  }

  surface_var_scaled <- surface_var_scaled[inhull,]

  surface_suitability <- matrix(NA, nrow = nrow(surface_var), ncol = length(taxa_names(physeq)))

  registerDoParallel(ncores)
  temp <- foreach(i = seq_along(taxa_names(physeq)), .combine=cbind) %dopar% {
    tryCatch(
      as.vector(predict(models[[i]], newx = as.matrix(surface_var_scaled), type = "response", s = "lambda.min")),
      error = function(e) return(NA))
  }
  stopImplicitCluster()
  surface_suitability[inhull,] <- temp

  surface_presence <- do.call(cbind, lapply(seq_along(thresholds),
                                            function(i) surface_suitability[,i] >= thresholds[[i]]))

  colnames(surface_suitability) <- colnames(surface_presence) <- taxa_names(physeq)

  # Break here for dev version
  if(dev_version) {
    message("Returning dev version. Only suitability, n_niches, and inhull returned.")
    return(cbind(surface_var,
                 mean_suitability = rowMeans(surface_suitability, na.rm=TRUE),
                 n_niches = rowSums(surface_presence, na.rm=TRUE),
                 inhull = inhull))
  }

  surface_n_limits <- rep(NA, nrow(surface_var))
  if(progressbar) {
    message("Progress for step 2 of 2:")
    pb <- txtProgressBar(max=n_gradient_steps^2, style=3)
  }
  for(i1 in seq_len(n_gradient_steps)) {
    for(i2 in seq_len(n_gradient_steps)) {
      this_ind <- (i1-1)*n_gradient_steps + i2
      this_ord <- c(i1, i2)
      this_val <- surface_presence[this_ind,]
      if(!identical(hull, "none")) { # If a hull is used,
        if(!inhull[this_ind]) {      # and this cell isn't within the hull,
          surface_n_limits[this_ind] <- NA # then skip it.
          if(progressbar) setTxtProgressBar(pb, this_ind)
          next
        }
      }
      offsets <- rbind(
        c(-1, 0),
        c(0, -1),
        c(1, 0),
        c(0, 1)
      )
      neighbors_ord <- sweep(offsets, 2, this_ord, '+')
      neighbors_ord <- neighbors_ord[which(apply(neighbors_ord, 1, function(x) all(x %in% seq_len(n_gradient_steps)))),]
      neighbors_ind <- apply(neighbors_ord, 1, function(x) (x[1]-1)*n_gradient_steps + x[2])
      neighbors_val <- foreach(k = neighbors_ind) %do% {
        return(surface_presence[k,])
      }
      crossings <- lapply(neighbors_val, function(x) x==0 & this_val==1) # "inner" boundary
      isboundary <- rowSums(do.call(cbind, crossings)) > 0
      surface_n_limits[this_ind] <- sum(isboundary, na.rm=TRUE)
      if(progressbar) setTxtProgressBar(pb, this_ind)
    }
  }
  if(progressbar) close(pb)

  # Get the mean suitability and total number of niches within each climate voxel
  surface_mean_suitability <- rowMeans(surface_suitability, na.rm=TRUE)
  surface_n_niches <- rowSums(surface_presence, na.rm=TRUE)
  surface_mean_suitability[!inhull] <- NA
  surface_n_niches[!inhull] <- NA

  # No longer supporting the total area of each climate voxel in the 2D version.
  # See getNicheMetrics5D for implementation ideas -- specifically for
  # the half-step size staggering of the interval-finding code.

  # # Bring in climate raster
  # df_climate <- as.matrix(r_climate)
  # int <- cbind(
  #   findInterval(df_climate[,1], quantize(obs_env[,1], n_gradient_steps)),
  #   findInterval(df_climate[,2], quantize(obs_env[,2], n_gradient_steps)),
  #   findInterval(df_climate[,3], quantize(obs_env[,3], n_gradient_steps)),
  #   findInterval(df_climate[,4], quantize(obs_env[,4], n_gradient_steps))
  # )
  # int_vec <- apply(int, 1, paste0, collapse=".")
  # int_vec[grep("NA", int_vec)] <- NA
  #
  # # Get the total area within each climate voxel
  # area <- as.vector(area(r_climate))
  # int_area <- cbind.data.frame(int = int_vec, area) %>%
  #   group_by(int) %>%
  #   summarise(area = sum(area)) %>%
  #   mutate(int = as.character(int))
  # surface_area <- int_area$area[match(surface_ordinal_vec, int_area$int)]
  # surface_area[which(is.na(surface_area))] <- 0

  # Calculate Sorensen-like climate sensitivity
  surface_sorensen <- surface_n_limits / surface_n_niches
  # # As raster map:
  # r_sorensen <- raster(r_climate)
  # r_sorensen[] <- surface_sorensen[match(int_vec, surface_ordinal_vec)]

  # # Since metrics are only defined within intervals, discard the ends.
  # # In this case, only the upper ends remain to be discarded.
  # discard_surface <- which(rowSums(apply(surface_ordinal, 2, function(x) x == n_gradient_steps)) > 0)
  # surface_area[discard_surface] <- 0

  # Collect variables
  out_df <- cbind(surface_var,
                  # mean_suitability = surface_mean_suitability,
                  n_niches = surface_n_niches,
                  n_limits = surface_n_limits,
                  sorensen = surface_sorensen,
                  # area = surface_area,
                  inhull = inhull)
  return(out_df)
}


getNicheMetrics2D <- function(physeq, models, thresholds, x1, x2,
                              x1_range=NULL, x2_range=NULL, constant_values=NULL,
                              r_climate, n_bins_per_axis=50,
                              hull=c("2D", "4D", "none"), inflate=0, ncores=1,
                              progressbar=TRUE, dev_version=FALSE,
                              pred_vars = c("mat_celsius", "temp_seasonality", "map_mm", "mat_celsius_2", "map_mm_2", "soilInCaClpH", "organicCPercent")) {

  if(!x1 %in% pred_vars | !x2 %in% pred_vars) {
    stop("`x1` and `x2` must be among `pred_vars`")
  }
  hull <- match.arg(hull, c("2D", "4D", "none"))
  n_gradient_steps <- n_bins_per_axis + 1
  n_taxa <- length(taxa_names(physeq))
  obs_env <- as(sample_data(physeq), "data.frame")[, pred_vars]
  obs_env <- obs_env[complete.cases(obs_env),]

  env_medians <- apply(obs_env, 2, function(x) median(x, na.rm=TRUE))
  if(is.null(constant_values)) {
    env_constant <- env_medians
  } else {
    if(length(constant_values) != length(pred_vars)) {
      stop("`constant_values` must have same length as `pred_vars`.
           To use the median of observed values for a variable, enter NA for
           that variable.")
    }
    env_constant <- constant_values
    env_constant[which(is.na(constant_values))] <- env_medians[which(is.na(constant_values))]
    names(env_constant) <- pred_vars
  }
  message("Using as constant unless otherwise specified: ", paste0(pred_vars, " ", env_constant, ", "))

  if(is.null(x1_range)) {
    x1_range <- range(obs_env[,x1], na.rm=TRUE)
    message("Using as `x1_range`: ", paste0(x1_range, collapse=", "))
  }
  if(is.null(x2_range)) {
    x2_range <- range(obs_env[,x2], na.rm=TRUE)
    message("Using as `x2_range`: ", paste0(x2_range, collapse=", "))
  }

  quantize <- function(x, n) {
    seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length.out=n)
  }
  surface_var <- tidyr::expand_grid(
    # quantize(obs_env[,x1], n_gradient_steps),
    # quantize(obs_env[,x2], n_gradient_steps)
    quantize(x1_range, n_gradient_steps),
    quantize(x2_range, n_gradient_steps)
  )
  names(surface_var) <- c(x1, x2)
  for(p in pred_vars) {
    if(!p %in% colnames(surface_var)) {
      surface_var[,p] <- env_constant[p]
    }
  }

  covars <- covarNamesFromLognets(models)
  surface_var <- includeQuadraticTerms(surface_var)
  surface_var <- surface_var[,match(covars,colnames(surface_var))]
  predictors <- as(sample_data(physeq)[,covars], "data.frame")
  predictors <- predictors[,match(covars, colnames(predictors))]
  surface_var_scaled <- scaleToReference(surface_var, predictors)

  registerDoParallel(ncores)
  surface_suitability <- foreach(i = seq_along(taxa_names(physeq)), .combine=cbind) %dopar% {
    tryCatch(
      as.vector(predict(models[[i]], newx = as.matrix(surface_var_scaled), type = "response", s = "lambda.min")),
      error = function(e) return(NA))
  }
  stopImplicitCluster()
  surface_presence <- do.call(cbind, lapply(seq_along(thresholds),
                                            function(i) surface_suitability[,i] >= thresholds[[i]]))

  colnames(surface_suitability) <- colnames(surface_presence) <- taxa_names(physeq)

  # Options for constraining by a convex hull

  # # If constraining by the convex hull around observations for each pair of climate axes
  # inhull_list <- list()
  # if(identical(hull, "4D")) {
  #   for(i in 1:3) {
  #     for(j in (i+1):4) {
  #       # cvhull <- obs_env[chull(obs_env[,i], obs_env[,j]),]
  #       # inhull_list <- c(inhull_list, list(inpolygon(surface_var[,i], surface_var[,j], cvhull[,i], cvhull[,j], boundary=TRUE)))
  #       x <- obs_env[,i]
  #       y <- obs_env[,j]
  #       ch_ind <- chull(x, y)
  #       infl_x <- x[ch_ind] + inflate * (x[ch_ind] - mean(x))
  #       infl_y <- y[ch_ind] + inflate * (y[ch_ind] - mean(y))
  #       inhull_list <- c(inhull_list, list(inpolygon(surface_var[,i], surface_var[,j], infl_x, infl_y, boundary=TRUE)))
  #     }
  #   }
  #   inhull <- rowSums(do.call(cbind, inhull_list)) == 6
  # }

  # If constraining by the convex hull around observations for the two climate axes of interest
  if(identical(hull, "2D")) {
    x <- obs_env[,x1]
    y <- obs_env[,x2]
    ch_ind <- chull(x, y)
    infl_x <- x[ch_ind] + inflate * (x[ch_ind] - mean(x))
    infl_y <- y[ch_ind] + inflate * (y[ch_ind] - mean(y))
    inhull <- pracma::inpolygon(surface_var[,x1,drop=TRUE], surface_var[,x2,drop=TRUE], infl_x, infl_y, boundary=TRUE)
  }

  # Otherwise, inhull is NA
  if(identical(hull, "none")) {
    inhull <- NA
  }

  # Break here for dev version
  if(dev_version) {
    message("Returning dev version. Only suitability, n_niches, and inhull returned.")
    return(cbind(surface_var,
                 mean_suitability = rowMeans(surface_suitability, na.rm=TRUE),
                 n_niches = rowSums(surface_presence, na.rm=TRUE),
                 inhull = inhull))
  }

  surface_n_limits <- rep(NA, nrow(surface_var))
  if(progressbar) {
    message("Progress for step 2 of 2:")
    pb <- txtProgressBar(max=n_gradient_steps^2, style=3)
  }
  for(i1 in seq_len(n_gradient_steps)) {
    for(i2 in seq_len(n_gradient_steps)) {
      this_ind <- (i1-1)*n_gradient_steps + i2
      this_ord <- c(i1, i2)
      this_val <- surface_presence[this_ind,]
      if(!identical(hull, "none")) { # If a hull is used,
        if(!inhull[this_ind]) {      # and this cell isn't within the hull,
          surface_n_limits[this_ind] <- NA # then skip it.
          if(progressbar) setTxtProgressBar(pb, this_ind)
          next
        }
      }
      offsets <- rbind(
        c(-1, 0),
        c(0, -1),
        c(1, 0),
        c(0, 1)
      )
      neighbors_ord <- sweep(offsets, 2, this_ord, '+')
      neighbors_ord <- neighbors_ord[which(apply(neighbors_ord, 1, function(x) all(x %in% seq_len(n_gradient_steps)))),]
      neighbors_ind <- apply(neighbors_ord, 1, function(x) (x[1]-1)*n_gradient_steps + x[2])
      neighbors_val <- foreach(k = neighbors_ind) %do% {
        return(surface_presence[k,])
      }
      crossings <- lapply(neighbors_val, function(x) x==0 & this_val==1) # "inner" boundary
      isboundary <- rowSums(do.call(cbind, crossings)) > 0
      surface_n_limits[this_ind] <- sum(isboundary, na.rm=TRUE)
      if(progressbar) setTxtProgressBar(pb, this_ind)
    }
  }
  if(progressbar) close(pb)

  # Get the mean suitability and total number of niches within each climate voxel
  surface_mean_suitability <- rowMeans(surface_suitability, na.rm=TRUE)
  surface_n_niches <- rowSums(surface_presence, na.rm=TRUE)
  surface_mean_suitability[!inhull] <- NA
  surface_n_niches[!inhull] <- NA

  # No longer supporting the total area of each climate voxel in the 2D version.
  # See getNicheMetrics5D for implementation ideas -- specifically for
  # the half-step size staggering of the interval-finding code.

  # # Bring in climate raster
  # df_climate <- as.matrix(r_climate)
  # int <- cbind(
  #   findInterval(df_climate[,1], quantize(obs_env[,1], n_gradient_steps)),
  #   findInterval(df_climate[,2], quantize(obs_env[,2], n_gradient_steps)),
  #   findInterval(df_climate[,3], quantize(obs_env[,3], n_gradient_steps)),
  #   findInterval(df_climate[,4], quantize(obs_env[,4], n_gradient_steps))
  # )
  # int_vec <- apply(int, 1, paste0, collapse=".")
  # int_vec[grep("NA", int_vec)] <- NA
  #
  # # Get the total area within each climate voxel
  # area <- as.vector(area(r_climate))
  # int_area <- cbind.data.frame(int = int_vec, area) %>%
  #   group_by(int) %>%
  #   summarise(area = sum(area)) %>%
  #   mutate(int = as.character(int))
  # surface_area <- int_area$area[match(surface_ordinal_vec, int_area$int)]
  # surface_area[which(is.na(surface_area))] <- 0

  # Calculate Sorensen-like climate sensitivity
  surface_sorensen <- surface_n_limits / surface_n_niches
  # # As raster map:
  # r_sorensen <- raster(r_climate)
  # r_sorensen[] <- surface_sorensen[match(int_vec, surface_ordinal_vec)]

  # # Since metrics are only defined within intervals, discard the ends.
  # # In this case, only the upper ends remain to be discarded.
  # discard_surface <- which(rowSums(apply(surface_ordinal, 2, function(x) x == n_gradient_steps)) > 0)
  # surface_area[discard_surface] <- 0

  # Collect variables
  out_df <- cbind(surface_var,
                  # mean_suitability = surface_mean_suitability,
                  n_niches = surface_n_niches,
                  n_limits = surface_n_limits,
                  sorensen = surface_sorensen,
                  # area = surface_area,
                  inhull = inhull)
  return(out_df)
}



# plotNicheMetrics2D <- function(physeq, metrics_df, site_points_in_plots = c(TRUE, TRUE, TRUE),
#                                metric1="n_niches", metric2="n_limits", metric3="sorensen",
#                                label1="Niche count", label2="Niche edge\ncount",
#                                label3="Sørensen\nsensitivity",
#                                p1_breaks=NULL, p2_breaks=NULL, p3_breaks=NULL,
#                                arrange=c("vertical", "horizontal"), progressbar=TRUE) {
#   arrange <- match.arg(arrange, c("vertical", "horizontal"))
#   pred_vars <- c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality", "k")
#   obs_env_withproject <- sample_data(physeq) %>%
#     as("data.frame") %>%
#     dplyr::select(all_of(c(pred_vars, "Project"))) %>%
#     mutate(Project = factor(Project, levels=c("DoB", "Both", "NEON")))
#   common_gg_elements <- list(
#     guides(color="none"),
#     xlab("Mean annual temp. (ºC)"),
#     ylab("Annual precip. (mm)"),
#     theme_bw(),
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#           panel.background = element_blank(), legend.position = "top",
#           legend.text = element_text(angle=45, hjust=1))
#   )
#   plotting_df <- dplyr::filter(metrics_df, inhull==TRUE)
#   p1 <- ggplot(plotting_df, aes(x=mat_celsius, y=map_mm)) +
#     geom_tile(aes_string(fill=metric1)) +
#     scale_fill_viridis(option="D", label1, breaks=if(is.null(p1_breaks)) waiver() else p1_breaks)  +
#     stat_contour(aes_string(z=metric1), col = "white", size=0.1, breaks=p1_breaks) +
#     common_gg_elements
#   p2 <- ggplot(plotting_df, aes(x=mat_celsius, y=map_mm)) +
#     geom_tile(aes_string(fill=metric2)) +
#     scale_fill_viridis(option="E", label2, breaks=if(is.null(p2_breaks)) waiver() else p2_breaks) +
#     stat_contour(aes_string(z=metric2), col = "white", size=0.1, breaks=p2_breaks) +
#     common_gg_elements
#   p3 <- ggplot(plotting_df, aes(x=mat_celsius, y=map_mm)) +
#     geom_tile(aes_string(fill=metric3)) +
#     scale_fill_viridis(option="B", label3, breaks=if(is.null(p3_breaks)) waiver() else p3_breaks) +
#     stat_contour(aes_string(z=metric3), col = "white", size=0.1, breaks=p3_breaks) +
#     common_gg_elements
#   if(identical(site_points_in_plots[1], TRUE)) {
#     p1 <- p1 + geom_point(data=obs_env_withproject, aes(col=Project), size=1, pch=21, fill="white")
#   }
#   if(identical(site_points_in_plots[2], TRUE)) {
#     p2 <- p2 + geom_point(data=obs_env_withproject, aes(col=Project), size=1, pch=21, fill="white")
#   }
#   if(identical(site_points_in_plots[3], TRUE)) {
#     p3 <- p3 + geom_point(data=obs_env_withproject, aes(col=Project), size=1, pch=21, fill="white")
#   }
#   if(identical(arrange, "vertical")) {
#     return(egg::ggarrange(p1, p2, p3, ncol=1))
#   } else {
#     return(egg::ggarrange(p1, p2, p3, ncol=3))
#   }
# }


# metrics_df <- out_df
# site_points_in_plots = c(TRUE, TRUE, TRUE)
# metric1="n_niches"
# metric2="n_limits"
# metric3="sorensen"
# label1="Niche count"
# label2="Niche edge\ncount"
# label3="Sørensen\nsensitivity"
# p1_breaks=NULL
# p2_breaks=NULL
# p3_breaks=NULL
# arrange=c("vertical", "horizontal")
# progressbar=TRUE


plotNicheMetrics2D <- function(physeq, metrics_df, site_points_in_plots = c(TRUE, TRUE, TRUE),
                               metric1="n_niches", metric2="n_limits", metric3="sorensen",
                               label1="Niche count", label2="Niche edge\ncount",
                               label3="Sørensen\nsensitivity",
                               p1_breaks=NULL, p2_breaks=NULL, p3_breaks=NULL,
                               arrange=c("vertical", "horizontal"), progressbar=TRUE,
                               pred_vars = c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality")) {
  arrange <- match.arg(arrange, c("vertical", "horizontal"))
  obs_env_withproject <- sample_data(physeq) %>%
    as("data.frame") %>%
    dplyr::select(all_of(c(pred_vars, "Project"))) %>%
    mutate(Project = factor(Project, levels=c("DoB", "Both", "NEON")))
  common_gg_elements <- list(
    guides(color="none"),
    xlab("Mean annual temp. (ºC)"),
    ylab("Annual precip. (mm)"),
    theme_bw(),
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), legend.position = "top",
          legend.text = element_text(angle=45, hjust=1))
  )
  # plotting_df <- dplyr::filter(metrics_df, inhull==TRUE)
  plotting_df <- metrics_df
  p1 <- ggplot(plotting_df, aes(x=mat_celsius, y=map_mm)) +
    geom_tile(aes_string(fill=metric1)) +
    viridis::scale_fill_viridis(option="D", label1, breaks=if(is.null(p1_breaks)) waiver() else p1_breaks)  +
    stat_contour(aes_string(z=metric1), col = "white", size=0.1, breaks=p1_breaks) +
    common_gg_elements
  p2 <- ggplot(plotting_df, aes(x=mat_celsius, y=map_mm)) +
    geom_tile(aes_string(fill=metric2)) +
    viridis::scale_fill_viridis(option="E", label2, breaks=if(is.null(p2_breaks)) waiver() else p2_breaks) +
    stat_contour(aes_string(z=metric2), col = "white", size=0.1, breaks=p2_breaks) +
    common_gg_elements
  p3 <- ggplot(plotting_df, aes(x=mat_celsius, y=map_mm)) +
    geom_tile(aes_string(fill=metric3)) +
    viridis::scale_fill_viridis(option="B", label3, breaks=if(is.null(p3_breaks)) waiver() else p3_breaks) +
    stat_contour(aes_string(z=metric3), col = "white", size=0.1, breaks=p3_breaks) +
    common_gg_elements
  if(identical(site_points_in_plots[1], TRUE)) {
    p1 <- p1 + geom_point(data=obs_env_withproject, aes(col=Project), size=1, pch=21, fill="white")
  }
  if(identical(site_points_in_plots[2], TRUE)) {
    p2 <- p2 + geom_point(data=obs_env_withproject, aes(col=Project), size=1, pch=21, fill="white")
  }
  if(identical(site_points_in_plots[3], TRUE)) {
    p3 <- p3 + geom_point(data=obs_env_withproject, aes(col=Project), size=1, pch=21, fill="white")
  }
  if(identical(arrange, "vertical")) {
    return(egg::ggarrange(p1, p2, p3, ncol=1))
  } else {
    return(egg::ggarrange(p1, p2, p3, ncol=3))
  }
}



getNicheMetrics4D <- function(physeq, models, thresholds, r_climate, n_bins_per_axis=10,
                                  inflate=0.15, ncores=1, progressbar=TRUE) {
  n_gradient_steps <- n_bins_per_axis + 1
  pred_vars <- c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality", "k")
  obs_env <- as(sample_data(physeq), "data.frame")[, pred_vars]
  obs_env <- obs_env[complete.cases(obs_env),]

  n_taxa <- length(taxa_names(physeq))
  obs_env <- as(sample_data(physeq), "data.frame")[, pred_vars]
  obs_env <- obs_env[complete.cases(obs_env),]

  quantize <- function(x, n) {
    seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length.out=n)
  }
  surface_var <- expand.grid(
    quantize(obs_env[,"mat_celsius"], n_gradient_steps),
    quantize(obs_env[,"temp_seasonality"], n_gradient_steps),
    quantize(obs_env[,"map_mm"], n_gradient_steps),
    quantize(obs_env[,"prec_seasonality"], n_gradient_steps)
  )
  surface_var[,5] <- yasso_k(surface_var[,1], surface_var[,3])
  names(surface_var) <- pred_vars

  covars <- covarNamesFromLognets(models)
  surface_var <- expandTerms(surface_var)
  predictors <- as(sample_data(physeq)[,covars], "data.frame")
  keep_cols <- which(colnames(surface_var) %in% covars)
  surface_var <- surface_var[,keep_cols]
  surface_var_scaled <- scaleToReference(surface_var, predictors)

  registerDoParallel(ncores)
  surface_suitability <- foreach(i = seq_along(taxa_names(physeq)), .combine=cbind) %dopar% {
    tryCatch(
      as.vector(predict(models[[i]], newx = as.matrix(surface_var_scaled), type = "response", s = "lambda.min")),
      error = function(e) return(NA))
  }
  stopImplicitCluster()
  surface_presence <- do.call(cbind, lapply(seq_along(thresholds),
                                            function(i) surface_suitability[,i] >= thresholds[[i]]))

  colnames(surface_suitability) <- colnames(surface_presence) <- taxa_names(physeq)

  # Constrain by the convex hull around observations for each pair of climate axes
  inhull <- list()
  for(i in 1:3) {
    for(j in (i+1):4) {
      x <- obs_env[,i]
      y <- obs_env[,j]
      ch_ind <- chull(x, y)
      infl_x <- x[ch_ind] + inflate * (x[ch_ind] - mean(x))
      infl_y <- y[ch_ind] + inflate * (y[ch_ind] - mean(y))
      inhull <- c(inhull, list(inpolygon(surface_var[,i], surface_var[,j], infl_x, infl_y, boundary=TRUE)))
    }
  }
  inhyperhull <- rowSums(do.call(cbind, inhull)) == 6

  surface_ordinal <- expand.grid(
    seq_len(n_gradient_steps),
    seq_len(n_gradient_steps),
    seq_len(n_gradient_steps),
    seq_len(n_gradient_steps)
  )
  surface_ordinal_vec <- apply(surface_ordinal, 1, paste, collapse='.')

  surface_n_limits <- rep(NA, nrow(surface_var))
  if(progressbar) {
    message("Progress for step 2 of 2:")
    # pb <- txtProgressBar(max=n_gradient_steps^4, style=3)
    # pb <- txtProgressBar(max=n_gradient_steps^3, style=3)
  }
  registerDoParallel(ncores)
  # for(i1 in seq_len(n_gradient_steps)) {
  out <- foreach(i = seq_len(n_gradient_steps), .combine=rbind) %dopar% {
    n_limits_by_ind <- list()
    for(i2 in seq_len(n_gradient_steps)) {
      for(i3 in seq_len(n_gradient_steps)) {
        for(i4 in seq_len(n_gradient_steps)) {
          this_ord <- c(i, i2, i3, i4)
          this_ind <- which(surface_ordinal_vec == paste0(this_ord, collapse="."))
          this_val <- surface_presence[this_ind,]
          if(!inhyperhull[this_ind]) {
            n_limits_by_ind <- c(n_limits_by_ind, list(c(this_ind, NA)))
            next
          }
          offsets <- rbind(
            c(-1, 0, 0, 0),
            c(0, -1, 0, 0),
            c(0, 0, -1, 0),
            c(0, 0, 0, -1),
            c(1, 0, 0, 0),
            c(0, 1, 0, 0),
            c(0, 0, 1, 0),
            c(0, 0, 0, 1)
          )
          neighbors_ord <- sweep(offsets, 2, this_ord, '+')
          neighbors_ord <- neighbors_ord[which(apply(neighbors_ord, 1, function(x) all(x %in% seq_len(n_gradient_steps)))),]
          neighbors_vec <- apply(neighbors_ord, 1, paste, collapse='.')
          neighbors_ind <- which(surface_ordinal_vec %in% neighbors_vec)
          neighbors_val <- foreach(k = neighbors_ind) %do% {
            return(surface_presence[k,])
          }
          crossings <- lapply(neighbors_val, function(x) x==0 & this_val==1) # "inner" boundary
          isboundary <- rowSums(do.call(cbind, crossings)) > 0
          n_limits_by_ind <- c(n_limits_by_ind, list(c(this_ind, sum(isboundary, na.rm=TRUE))))
        }
      }
    }
    do.call(rbind, n_limits_by_ind)
  }
  stopImplicitCluster()
  surface_n_limits <- out[order(out[,1]),2]

  # Get the mean suitability and total number of niches within each climate voxel
  surface_mean_suitability <- rowMeans(surface_suitability, na.rm=TRUE)
  surface_n_niches <- rowSums(surface_presence, na.rm=TRUE)

  # Bring in climate raster
  step_sizes <- apply(surface_var[,1:4], 2, function(x) min(diff(unique(x))))
  # The purpose of adding half a step size is to find the gridcell that each location
  # is "closest" to in terms of climate space. Otherwise a location that is between
  # values 1 and 2 of an axis will be assigned to value 1, no matter how close it is
  # to value 2.
  df_climate <- as.matrix(r_climate)[,pred_vars]
  int <- cbind(
    findInterval(df_climate[,1] + step_sizes[1]/2, quantize(obs_env[,1], n_gradient_steps)),
    findInterval(df_climate[,2] + step_sizes[2]/2, quantize(obs_env[,2], n_gradient_steps)),
    findInterval(df_climate[,3] + step_sizes[3]/2, quantize(obs_env[,3], n_gradient_steps)),
    findInterval(df_climate[,4] + step_sizes[4]/2, quantize(obs_env[,4], n_gradient_steps))
  )
  int_vec <- apply(int, 1, paste0, collapse=".")
  int_vec[grep("NA", int_vec)] <- NA

  # Get the total area within each climate voxel
  area <- as.vector(area(r_climate))
  int_area <- cbind.data.frame(int = int_vec, area) %>%
    group_by(int) %>%
    summarise(area = sum(area)) %>%
    mutate(int = as.character(int))
  surface_area <- int_area$area[match(surface_ordinal_vec, int_area$int)]
  surface_area[which(is.na(surface_area))] <- 0

  # Calculate Sorensen-like climate sensitivity
  surface_sorensen <- surface_n_limits / surface_n_niches
  # As raster map:
  r_sorensen <- raster(r_climate)
  r_sorensen[] <- surface_sorensen[match(int_vec, surface_ordinal_vec)]

  # Since metrics are only defined within intervals, discard the ends.
  # In this case, only the upper ends remain to be discarded.
  discard_surface <- which(rowSums(apply(surface_ordinal, 2, function(x) x == n_gradient_steps)) > 0)
  surface_area[discard_surface] <- 0

  # Collect variables
  out_df <- cbind(surface_var, mean_suitability = surface_mean_suitability,
                  n_limits = surface_n_limits, n_niches = surface_n_niches,
                  sorensen = surface_sorensen, area = surface_area, inhull = inhyperhull)
  return(list(data = out_df, r_sor = r_sorensen))
}


recursiveList <- function(n_levels, n_elements) {
  if(n_levels <= 0 | n_elements <= 0 |
     n_levels != round(n_levels) | n_elements != round(n_elements)) {
    warning("n_levels and n_elements must be positive integers. Returning NULL.")
    return(NULL)
  }
  if(n_levels==1) {
    return(vector("list", n_elements))
  }
  if(n_levels > 1) {
    return(do.call(list, lapply(1:n_elements, function(i) recursiveList(n_levels - 1, n_elements))))
  }
}

covarNamesFromLognets <- function(models) {
  for(i in 1:length(models)) {
    if(identical(models[[i]], NA)) {
      if(i == length(models)) {
        warning("All models are NA! Returning all NAs.")
        return(rep(NA, length(models)))
      }
    } else {
      return(rownames(coef(models[[i]]))[-1])
    }
  }
}

# # Predictor variables here are hard-coded for now
# physeq <- neon_dob_prevalent
# sample_data(physeq) %>%
#   as("data.frame") %>%
#   mutate(mat_celsius_2 = mat_celsius^2,
#          map_mm_2 = map_mm^2) ->
#   sample_data(physeq)
# physeq <- prune_taxa(sample(taxa_names(physeq), 100), physeq)
#
# terms_top <- c("mat_celsius", "temp_seasonality", "map_mm", "mat_celsius_2", "map_mm_2", "soilInCaClpH", "organicCPercent")
# predictors <- as(sample_data(physeq)[,terms_top], "data.frame")
# predictors_std <- scale(predictors, center=TRUE, scale=TRUE)
# complete_records <- which(complete.cases(predictors_std))
#
# registerDoParallel(8)
# models <-  foreach(i = seq_along(taxa_names(physeq))) %dopar% {
#   y01 <- as.vector(otu_table(physeq)[,i]) > 0
#   set.seed(1010101)
#   tryCatch(cv.glmnet(as.matrix(predictors_std)[complete_records,],
#                      y01[complete_records], family = "binomial",
#                      alpha = 0, type.measure = "default", standardize=FALSE),
#            error = function(e) return(NA))
# }
# stopImplicitCluster()
# thresholds <- findThresholds(physeq, models, terms = terms_top)
# r_climate <- r_present_northam
# ncores <- 16
# var_ranges <- NULL
# n_bins_per_axis <- 9
# pred_vars = terms_top
# count_edges_on_axis = c(TRUE, FALSE, TRUE, FALSE, FALSE)

getNicheMetrics5D <- function(physeq, models, thresholds, r_climate, n_bins_per_axis=8,
                              var_ranges = NULL, inflate=0, ncores=1,
                              count_edges_on_axis = c(TRUE, TRUE, TRUE, TRUE, TRUE),
                              count_edges_on_axis2 = c(TRUE, TRUE, TRUE, TRUE, TRUE),
                              pred_vars = c("mat_celsius", "temp_seasonality", "map_mm", "mat_celsius_2", "map_mm_2", "soilInCaClpH", "organicCPercent")) {

  if(!is.logical(count_edges_on_axis)) stop("`count_edges_on_axis` should be of class 'logical' (TRUE/FALSE).")
  if(length(count_edges_on_axis) != 5) {
    stop(paste0("getNicheMetrics5D is designed to be run on a set of covariates containing
                exactly 5 (nonquadratic) terms, but `count_edges_on_axis` has length ",
                length(count_edges_on_axis)))
  }

  n_gradient_steps <- n_bins_per_axis + 1

  obs_env <- as(sample_data(physeq), "data.frame")[, pred_vars]
  obs_env <- obs_env[complete.cases(obs_env),]

  # If the square of a variable is equal to another,
  # then that other is a quadratic version. Do not use quadratic
  # versions for expand.grid.
  is_quadratic_version <- rep(FALSE, ncol(obs_env))
  for(i in 1:ncol(obs_env)) {
    is_quadratic_version <- is_quadratic_version |
      apply(obs_env, 2, function(x) all(x == obs_env[,i]^2))
  }
  pred_vars_nonquadratic <- pred_vars[!is_quadratic_version]

  if(length(pred_vars_nonquadratic) != 5) {
    warning(paste0("getNicheMetrics5D is designed to be run on a set of covariates containing exactly 5 (nonquadratic) terms, but only ",
                   length(pred_vars_nonquadratic), " detected. Interpret results with caution."))
  }

  if(is.null(var_ranges)) {
    var_ranges <- lapply(pred_vars_nonquadratic, function(x) range(obs_env[,x], na.rm=TRUE))
    message("Using as min/max of all quadratic variables to define the `var_ranges` over which to create grid.")
  }

  quantize <- function(x, n) {
    seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length.out=n)
  }
  surface_var <- do.call(tidyr::expand_grid, lapply(1:5, function(i) quantize(var_ranges[[i]], n_gradient_steps)))
  names(surface_var) <- pred_vars_nonquadratic

  step_sizes <- apply(obs_env[,pred_vars_nonquadratic], 2, function(x) min(diff(unique(x))))

  covars <- covarNamesFromLognets(models)
  surface_var <- includeQuadraticTerms(surface_var)
  surface_var <- surface_var[,match(covars,colnames(surface_var))]
  predictors <- as(sample_data(physeq)[,covars], "data.frame")
  predictors <- predictors[,match(covars, colnames(predictors))]
  surface_var_scaled <- scaleToReference(surface_var, predictors)

  # Constrain by the convex hull around observations for each pair of environmental axes.
  # Extra hull allowance for neighbors to allow n_limits to be calculated even on the edge
  # of the within-hull environmental space.
  inhull <- list()
  inhull_neighbor_allowance <- list()
  for(i in 1:4) {
    for(j in (i+1):5) {
      varname_x <- pred_vars_nonquadratic[i]
      varname_y <- pred_vars_nonquadratic[j]
      x <- obs_env[,varname_x]
      y <- obs_env[,varname_y]
      complete_ind <- which(complete.cases(cbind(x, y)))
      x <- x[complete_ind]
      y <- y[complete_ind]
      ch_ind <- chull(x, y)
      infl_x <- x[ch_ind] + inflate * (x[ch_ind] - mean(x))
      infl_y <- y[ch_ind] + inflate * (y[ch_ind] - mean(y))
      surface_x <- surface_var[,varname_x,drop=TRUE]
      surface_y <- surface_var[,varname_y,drop=TRUE]
      this_inhull <- inpolygon(
        surface_x, surface_y,
        infl_x, infl_y, boundary=TRUE
      )

      ch_ind_of_surface <- chull(surface_x[this_inhull],
                                 surface_y[this_inhull])
      surface_hull_x <- surface_x[this_inhull][ch_ind_of_surface]
      surface_hull_y <- surface_y[this_inhull][ch_ind_of_surface]
      neighbors_x <- c(rep(surface_hull_x - step_sizes[varname_x], 2), rep(surface_hull_x + step_sizes[varname_x], 2))
      neighbors_y <- rep(c(surface_hull_y - step_sizes[varname_y], surface_hull_y + step_sizes[varname_y]), 2)
      ch_ind_of_neighbors <- chull(neighbors_x, neighbors_y)
      neighbor_hull_x <- neighbors_x[ch_ind_of_neighbors]
      neighbor_hull_y <- neighbors_y[ch_ind_of_neighbors]
      complete_ind <- which(complete.cases(cbind(neighbor_hull_x, neighbor_hull_y)))
      neighbor_hull_x <- neighbor_hull_x[complete_ind]
      neighbor_hull_y <- neighbor_hull_y[complete_ind]
      this_inhull_neighbor <- inpolygon(
        surface_x, surface_y,
        neighbor_hull_x, neighbor_hull_y, boundary=TRUE
      )

      inhull <- c(inhull, list(this_inhull))
      inhull_neighbor_allowance <- c(inhull_neighbor_allowance, list(this_inhull_neighbor))
    }
  }
  inhyperhull <- rowSums(do.call(cbind, inhull)) == n_gradient_steps
  inhyperhull_neighbor_allowance <- rowSums(do.call(cbind, inhull_neighbor_allowance)) == n_gradient_steps

  surface_var_scaled <- surface_var_scaled[inhyperhull_neighbor_allowance,]

  surface_suitability <- matrix(NA, nrow = nrow(surface_var), ncol = length(taxa_names(physeq)))

  registerDoParallel(ncores)
  temp <- foreach(i = seq_along(taxa_names(physeq)), .combine=cbind) %dopar% {
    tryCatch(
      as.vector(predict(models[[i]], newx = as.matrix(surface_var_scaled), type = "response", s = "lambda.min")),
      error = function(e) return(NA))
  }
  stopImplicitCluster()
  surface_suitability[inhyperhull_neighbor_allowance,] <- temp

  surface_presence <- do.call(cbind, lapply(seq_along(thresholds),
                                            function(i) surface_suitability[,i] >= thresholds[[i]]))

  colnames(surface_suitability) <- colnames(surface_presence) <- taxa_names(physeq)

  # The offsets define edge-adjacency by rook's rule (doesn't count diagonals)
  offsets0 <- rbind(
    c(-1, 0, 0, 0, 0),
    c(0, -1, 0, 0, 0),
    c(0, 0, -1, 0, 0),
    c(0, 0, 0, -1, 0),
    c(0, 0, 0, 0, -1),
    c(1, 0, 0, 0, 0),
    c(0, 1, 0, 0, 0),
    c(0, 0, 1, 0, 0),
    c(0, 0, 0, 1, 0),
    c(0, 0, 0, 0, 1)
  )
  # The offsets can be further constrained by axes of interest
  offsets <- offsets0[rep(count_edges_on_axis, 2),]
  offsets2 <- offsets0[rep(count_edges_on_axis, 2),]

  # Can't use progress bar within a parallelized process
  # if(progressbar) {
  #   message("Progress for step 2 of 2:")
  #   pb <- txtProgressBar(max=n_gradient_steps^5, style=3)
  # }

  registerDoParallel(ncores)
  # surface_n_limits <- foreach(i = seq_len(n_gradient_steps), .combine=c) %dopar% {
  oper <- foreach(i = seq_len(n_gradient_steps),
                              .combine = nested_combine, .multicombine=TRUE,
                              .init = list(list(), list())) %dopar% {
    n_limits_by_ind <- vector("integer", length=n_gradient_steps^4)
    n_limits_by_ind2 <- vector("integer", length=n_gradient_steps^4)
    for(i2 in seq_len(n_gradient_steps)) {
      for(i3 in seq_len(n_gradient_steps)) {
        for(i4 in seq_len(n_gradient_steps)) {
          for(i5 in seq_len(n_gradient_steps)) {
            this_ind <- (i-1)*n_gradient_steps^4 + (i2-1)*n_gradient_steps^3 + (i3-1)*n_gradient_steps^2 + (i4-1)*n_gradient_steps + i5
            this_ind_within_foreach_loop <- (i2-1)*n_gradient_steps^3 + (i3-1)*n_gradient_steps^2 + (i4-1)*n_gradient_steps + i5
            this_ord <- c(i, i2, i3, i4, i5)
            if(!inhyperhull[this_ind]) {
              n_limits_by_ind[this_ind_within_foreach_loop] <- NA
              # setTxtProgressBar(pb, this_ind)
              next
            }
            this_val <- surface_presence[this_ind,]
            neighbors_ord <- sweep(offsets, 2, this_ord, '+')
            neighbors_ord <- neighbors_ord[which(apply(neighbors_ord, 1, function(x) all(x %in% seq_len(n_gradient_steps)))),]
            neighbors_ind <- apply(neighbors_ord, 1, function(x) (x[1]-1)*n_gradient_steps^4 + (x[2]-1)*n_gradient_steps^3 + (x[3]-1)*n_gradient_steps^2 + (x[4]-1)*n_gradient_steps + x[5])
            neighbors_val <- lapply(neighbors_ind, function(k) surface_presence[k,])
            crossings <- lapply(neighbors_val, function(x) x==0 & this_val==1) # "inner" boundary
            isboundary <- rowSums(do.call(cbind, crossings)) > 0 # for what species is this cell a boundary along at least one axis?
            n_limits_by_ind[this_ind_within_foreach_loop] <- sum(isboundary, na.rm=TRUE)
            # setTxtProgressBar(pb, this_ind)

            # Count limits across the optional 2nd set of offsets
            neighbors_ord <- sweep(offsets2, 2, this_ord, '+')
            neighbors_ord <- neighbors_ord[which(apply(neighbors_ord, 1, function(x) all(x %in% seq_len(n_gradient_steps)))),]
            neighbors_ind <- apply(neighbors_ord, 1, function(x) (x[1]-1)*n_gradient_steps^4 + (x[2]-1)*n_gradient_steps^3 + (x[3]-1)*n_gradient_steps^2 + (x[4]-1)*n_gradient_steps + x[5])
            neighbors_val <- lapply(neighbors_ind, function(k) surface_presence[k,])
            crossings <- lapply(neighbors_val, function(x) x==0 & this_val==1) # "inner" boundary
            isboundary <- rowSums(do.call(cbind, crossings)) > 0 # for what species is this cell a boundary along at least one axis?
            n_limits_by_ind2[this_ind_within_foreach_loop] <- sum(isboundary, na.rm=TRUE)
          }
        }
      }
    }
    list(n_limits_by_ind, n_limits_by_ind2)
  }
  stopImplicitCluster()
  surface_n_limits <- unlist(oper[[1]])
  surface_n_limits2 <- unlist(oper[[2]])

  # # Surprisingly to me, this implementation is slower:
  # t0 <- Sys.time()
  # registerDoParallel(ncores)
  # surface_n_limits <- foreach(i = seq_len(n_gradient_steps^5), .combine=c) %dopar% {
  #   if(!inhyperhull[i]) return(NA)
  #   i1 <- (i-1) %/% n_gradient_steps^4 + 1
  #   r1 <- (i-1) %% n_gradient_steps^4
  #   i2 <- r1 %/% n_gradient_steps^3 + 1
  #   r2 <- r1 %% n_gradient_steps^3
  #   i3 <- r2 %/% n_gradient_steps^2 + 1
  #   r3 <- r2 %% n_gradient_steps^2
  #   i4 <- r3 %/% n_gradient_steps + 1
  #   r4 <- r3 %% n_gradient_steps
  #   i5 <- r4 + 1
  #   this_ord <- c(i1, i2, i3, i4, i5)
  #   this_val <- surface_presence[i,]
  #   neighbors_ord <- sweep(offsets, 2, this_ord, '+')
  #   neighbors_ord <- neighbors_ord[which(apply(neighbors_ord, 1, function(x) all(x %in% seq_len(n_gradient_steps)))),]
  #   neighbors_ind <- apply(neighbors_ord, 1, function(x) (x[1]-1)*n_gradient_steps^4 + (x[2]-1)*n_gradient_steps^3 + (x[3]-1)*n_gradient_steps^2 + (x[4]-1)*n_gradient_steps + x[5])
  #   neighbors_val <- lapply(neighbors_ind, function(k) surface_presence[k,])
  #   crossings <- lapply(neighbors_val, function(x) x==0 & this_val==1) # "inner" boundary
  #   isboundary <- rowSums(do.call(cbind, crossings)) > 0 # for what species is this cell a boundary along at least one axis?
  #   return(sum(isboundary, na.rm=TRUE))
  # }
  # stopImplicitCluster()
  # t0.1 <- Sys.time()
  # t0.1 - t0

  # Get the mean suitability and total number of niches within each climate voxel
  surface_mean_suitability <- rowMeans(surface_suitability, na.rm=TRUE)
  surface_n_niches <- rowSums(surface_presence, na.rm=TRUE)

  # This next code chunk was originally for calculating the amount of raster area in each voxel. It would
  # require combining the climate and soil rasters, however.

  # # Bring in climate raster
  # step_sizes <- apply(surface_var[,1:4], 2, function(x) min(diff(unique(x))))
  # df_climate <- as.matrix(r_climate)[,intersect(pred_vars, names(r_climate))]
  # int <- as.matrix(cbind.data.frame(
  #   lapply(1:ncol(df_climate),
  #          function(i) findInterval(df_climate[,i] + step_sizes[i]/2, quantize(obs_env[,i], n_gradient_steps)))
  # ))
  # int_vec <- apply(int, 1, paste0, collapse=".")
  # int_vec[grep("NA", int_vec)] <- NA
  #
  # # Get the total area within each climate voxel
  # area <- as.vector(area(r_climate))
  # int_area <- cbind.data.frame(int = int_vec, area) %>%
  #   group_by(int) %>%
  #   summarise(area = sum(area)) %>%
  #   mutate(int = as.character(int))
  # surface_area <- int_area$area[match(surface_ordinal_vec, int_area$int)]
  # surface_area[which(is.na(surface_area))] <- 0

  # Calculate Sorensen-like climate sensitivity
  surface_sorensen <- surface_n_limits / surface_n_niches
  # # As raster map:
  # r_sorensen <- raster(r_climate)
  # r_sorensen[] <- surface_sorensen[match(int_vec, surface_ordinal_vec)]

  # # Since metrics are only defined within intervals, discard the ends.
  # # In this case, only the upper ends remain to be discarded.
  # discard_surface <- which(rowSums(apply(surface_ordinal, 2, function(x) x == n_gradient_steps)) > 0)
  # surface_area[discard_surface] <- 0

  # Collect variables
  out_df <- cbind(surface_var, mean_suitability = surface_mean_suitability,
                  n_limits = surface_n_limits, n_limits2 = surface_n_limits2,
                  n_niches = surface_n_niches,
                  sorensen = surface_sorensen, inhull = inhyperhull)
  # return(list(data = out_df, r_sor = r_sorensen))
  return(out_df)
}

getNicheMetrics5D0 <- function(physeq, models, thresholds, r_climate, n_bins_per_axis=8,
                              var_ranges = NULL, inflate=0, ncores=1,
                              count_edges_on_axis = c(TRUE, TRUE, TRUE, TRUE, TRUE),
                              pred_vars = c("mat_celsius", "temp_seasonality", "map_mm", "mat_celsius_2", "map_mm_2", "soilInCaClpH", "organicCPercent")) {

  if(!is.logical(count_edges_on_axis)) stop("`count_edges_on_axis` should be of class 'logical' (TRUE/FALSE).")
  if(length(count_edges_on_axis) != 5) {
    stop(paste0("getNicheMetrics5D is designed to be run on a set of covariates containing
                exactly 5 (nonquadratic) terms, but `count_edges_on_axis` has length ",
                length(count_edges_on_axis)))
  }

  n_gradient_steps <- n_bins_per_axis + 1

  obs_env <- as(sample_data(physeq), "data.frame")[, pred_vars]
  obs_env <- obs_env[complete.cases(obs_env),]

  # If the square of a variable is equal to another,
  # then that other is a quadratic version. Do not use quadratic
  # versions for expand.grid.
  is_quadratic_version <- rep(FALSE, ncol(obs_env))
  for(i in 1:ncol(obs_env)) {
    is_quadratic_version <- is_quadratic_version |
      apply(obs_env, 2, function(x) all(x == obs_env[,i]^2))
  }
  pred_vars_nonquadratic <- pred_vars[!is_quadratic_version]

  if(length(pred_vars_nonquadratic) != 5) {
    warning(paste0("getNicheMetrics5D is designed to be run on a set of covariates containing exactly 5 (nonquadratic) terms, but only ",
                   length(pred_vars_nonquadratic), " detected. Interpret results with caution."))
  }

  if(is.null(var_ranges)) {
    var_ranges <- lapply(pred_vars_nonquadratic, function(x) range(obs_env[,x], na.rm=TRUE))
    message("Using as min/max of all quadratic variables to define the `var_ranges` over which to create grid.")
  }

  quantize <- function(x, n) {
    seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length.out=n)
  }
  surface_var <- do.call(tidyr::expand_grid, lapply(1:5, function(i) quantize(var_ranges[[i]], n_gradient_steps)))
  names(surface_var) <- pred_vars_nonquadratic

  step_sizes <- apply(obs_env[,pred_vars_nonquadratic], 2, function(x) min(diff(unique(x))))

  covars <- covarNamesFromLognets(models)
  surface_var <- includeQuadraticTerms(surface_var)
  surface_var <- surface_var[,match(covars,colnames(surface_var))]
  predictors <- as(sample_data(physeq)[,covars], "data.frame")
  predictors <- predictors[,match(covars, colnames(predictors))]
  surface_var_scaled <- scaleToReference(surface_var, predictors)

  # Constrain by the convex hull around observations for each pair of environmental axes.
  # Extra hull allowance for neighbors to allow n_limits to be calculated even on the edge
  # of the within-hull environmental space.
  inhull <- list()
  inhull_neighbor_allowance <- list()
  for(i in 1:4) {
    for(j in (i+1):5) {
      varname_x <- pred_vars_nonquadratic[i]
      varname_y <- pred_vars_nonquadratic[j]
      x <- obs_env[,varname_x]
      y <- obs_env[,varname_y]
      complete_ind <- which(complete.cases(cbind(x, y)))
      x <- x[complete_ind]
      y <- y[complete_ind]
      ch_ind <- chull(x, y)
      infl_x <- x[ch_ind] + inflate * (x[ch_ind] - mean(x))
      infl_y <- y[ch_ind] + inflate * (y[ch_ind] - mean(y))
      surface_x <- surface_var[,varname_x,drop=TRUE]
      surface_y <- surface_var[,varname_y,drop=TRUE]
      this_inhull <- inpolygon(
        surface_x, surface_y,
        infl_x, infl_y, boundary=TRUE
      )

      ch_ind_of_surface <- chull(surface_x[this_inhull],
                                 surface_y[this_inhull])
      surface_hull_x <- surface_x[this_inhull][ch_ind_of_surface]
      surface_hull_y <- surface_y[this_inhull][ch_ind_of_surface]
      neighbors_x <- c(rep(surface_hull_x - step_sizes[varname_x], 2), rep(surface_hull_x + step_sizes[varname_x], 2))
      neighbors_y <- rep(c(surface_hull_y - step_sizes[varname_y], surface_hull_y + step_sizes[varname_y]), 2)
      ch_ind_of_neighbors <- chull(neighbors_x, neighbors_y)
      neighbor_hull_x <- neighbors_x[ch_ind_of_neighbors]
      neighbor_hull_y <- neighbors_y[ch_ind_of_neighbors]
      complete_ind <- which(complete.cases(cbind(neighbor_hull_x, neighbor_hull_y)))
      neighbor_hull_x <- neighbor_hull_x[complete_ind]
      neighbor_hull_y <- neighbor_hull_y[complete_ind]
      this_inhull_neighbor <- inpolygon(
        surface_x, surface_y,
        neighbor_hull_x, neighbor_hull_y, boundary=TRUE
      )

      inhull <- c(inhull, list(this_inhull))
      inhull_neighbor_allowance <- c(inhull_neighbor_allowance, list(this_inhull_neighbor))
    }
  }
  inhyperhull <- rowSums(do.call(cbind, inhull)) == choose(5, 2)
  inhyperhull_neighbor_allowance <- rowSums(do.call(cbind, inhull_neighbor_allowance)) == choose(5, 2)

  surface_var_scaled <- surface_var_scaled[inhyperhull_neighbor_allowance,]

  surface_suitability <- matrix(NA, nrow = nrow(surface_var), ncol = length(taxa_names(physeq)))

  registerDoParallel(ncores)
  temp <- foreach(i = seq_along(taxa_names(physeq)), .combine=cbind) %dopar% {
    tryCatch(
      as.vector(predict(models[[i]], newx = as.matrix(surface_var_scaled), type = "response", s = "lambda.min")),
      error = function(e) return(NA))
  }
  stopImplicitCluster()
  surface_suitability[inhyperhull_neighbor_allowance,] <- temp

  surface_presence <- do.call(cbind, lapply(seq_along(thresholds),
                                            function(i) surface_suitability[,i] >= thresholds[[i]]))

  colnames(surface_suitability) <- colnames(surface_presence) <- taxa_names(physeq)

  # The offsets define edge-adjacency by rook's rule (doesn't count diagonals)
  offsets <- rbind(
    c(-1, 0, 0, 0, 0),
    c(0, -1, 0, 0, 0),
    c(0, 0, -1, 0, 0),
    c(0, 0, 0, -1, 0),
    c(0, 0, 0, 0, -1),
    c(1, 0, 0, 0, 0),
    c(0, 1, 0, 0, 0),
    c(0, 0, 1, 0, 0),
    c(0, 0, 0, 1, 0),
    c(0, 0, 0, 0, 1)
  )
  # The offsets can be further constrained by axes of interest
  offsets <- offsets[rep(count_edges_on_axis, 2),]

  # Can't use progress bar within a parallelized process
  # if(progressbar) {
  #   message("Progress for step 2 of 2:")
  #   pb <- txtProgressBar(max=n_gradient_steps^5, style=3)
  # }

  registerDoParallel(ncores)
  surface_n_limits <- foreach(i = seq_len(n_gradient_steps), .combine=c) %dopar% {
    n_limits_by_ind <- vector("integer", length=n_gradient_steps^4)
    for(i2 in seq_len(n_gradient_steps)) {
      for(i3 in seq_len(n_gradient_steps)) {
        for(i4 in seq_len(n_gradient_steps)) {
          for(i5 in seq_len(n_gradient_steps)) {
            this_ind <- (i-1)*n_gradient_steps^4 + (i2-1)*n_gradient_steps^3 + (i3-1)*n_gradient_steps^2 + (i4-1)*n_gradient_steps + i5
            this_ind_within_foreach_loop <- (i2-1)*n_gradient_steps^3 + (i3-1)*n_gradient_steps^2 + (i4-1)*n_gradient_steps + i5
            this_ord <- c(i, i2, i3, i4, i5)
            if(!inhyperhull[this_ind]) {
              n_limits_by_ind[this_ind_within_foreach_loop] <- NA
              # setTxtProgressBar(pb, this_ind)
              next
            }
            this_val <- surface_presence[this_ind,]
            neighbors_ord <- sweep(offsets, 2, this_ord, '+')
            neighbors_ord <- neighbors_ord[which(apply(neighbors_ord, 1, function(x) all(x %in% seq_len(n_gradient_steps)))),]
            neighbors_ind <- apply(neighbors_ord, 1, function(x) (x[1]-1)*n_gradient_steps^4 + (x[2]-1)*n_gradient_steps^3 + (x[3]-1)*n_gradient_steps^2 + (x[4]-1)*n_gradient_steps + x[5])
            neighbors_val <- lapply(neighbors_ind, function(k) surface_presence[k,])
            crossings <- lapply(neighbors_val, function(x) x==0 & this_val==1) # "inner" boundary
            isboundary <- rowSums(do.call(cbind, crossings)) > 0 # for what species is this cell a boundary along at least one axis?
            n_limits_by_ind[this_ind_within_foreach_loop] <- sum(isboundary, na.rm=TRUE)
            # setTxtProgressBar(pb, this_ind)
          }
        }
      }
    }
    n_limits_by_ind
  }
  stopImplicitCluster()

  # # Surprisingly to me, this implementation is slower:
  # t0 <- Sys.time()
  # registerDoParallel(ncores)
  # surface_n_limits <- foreach(i = seq_len(n_gradient_steps^5), .combine=c) %dopar% {
  #   if(!inhyperhull[i]) return(NA)
  #   i1 <- (i-1) %/% n_gradient_steps^4 + 1
  #   r1 <- (i-1) %% n_gradient_steps^4
  #   i2 <- r1 %/% n_gradient_steps^3 + 1
  #   r2 <- r1 %% n_gradient_steps^3
  #   i3 <- r2 %/% n_gradient_steps^2 + 1
  #   r3 <- r2 %% n_gradient_steps^2
  #   i4 <- r3 %/% n_gradient_steps + 1
  #   r4 <- r3 %% n_gradient_steps
  #   i5 <- r4 + 1
  #   this_ord <- c(i1, i2, i3, i4, i5)
  #   this_val <- surface_presence[i,]
  #   neighbors_ord <- sweep(offsets, 2, this_ord, '+')
  #   neighbors_ord <- neighbors_ord[which(apply(neighbors_ord, 1, function(x) all(x %in% seq_len(n_gradient_steps)))),]
  #   neighbors_ind <- apply(neighbors_ord, 1, function(x) (x[1]-1)*n_gradient_steps^4 + (x[2]-1)*n_gradient_steps^3 + (x[3]-1)*n_gradient_steps^2 + (x[4]-1)*n_gradient_steps + x[5])
  #   neighbors_val <- lapply(neighbors_ind, function(k) surface_presence[k,])
  #   crossings <- lapply(neighbors_val, function(x) x==0 & this_val==1) # "inner" boundary
  #   isboundary <- rowSums(do.call(cbind, crossings)) > 0 # for what species is this cell a boundary along at least one axis?
  #   return(sum(isboundary, na.rm=TRUE))
  # }
  # stopImplicitCluster()
  # t0.1 <- Sys.time()
  # t0.1 - t0

  # Get the mean suitability and total number of niches within each climate voxel
  surface_mean_suitability <- rowMeans(surface_suitability, na.rm=TRUE)
  surface_n_niches <- rowSums(surface_presence, na.rm=TRUE)

  # This next code chunk was originally for calculating the amount of raster area in each voxel. It would
  # require combining the climate and soil rasters, however.

  # # Bring in climate raster
  # step_sizes <- apply(surface_var[,1:4], 2, function(x) min(diff(unique(x))))
  # df_climate <- as.matrix(r_climate)[,intersect(pred_vars, names(r_climate))]
  # int <- as.matrix(cbind.data.frame(
  #   lapply(1:ncol(df_climate),
  #          function(i) findInterval(df_climate[,i] + step_sizes[i]/2, quantize(obs_env[,i], n_gradient_steps)))
  # ))
  # int_vec <- apply(int, 1, paste0, collapse=".")
  # int_vec[grep("NA", int_vec)] <- NA
  #
  # # Get the total area within each climate voxel
  # area <- as.vector(area(r_climate))
  # int_area <- cbind.data.frame(int = int_vec, area) %>%
  #   group_by(int) %>%
  #   summarise(area = sum(area)) %>%
  #   mutate(int = as.character(int))
  # surface_area <- int_area$area[match(surface_ordinal_vec, int_area$int)]
  # surface_area[which(is.na(surface_area))] <- 0

  # Calculate Sorensen-like climate sensitivity
  surface_sorensen <- surface_n_limits / surface_n_niches
  # # As raster map:
  # r_sorensen <- raster(r_climate)
  # r_sorensen[] <- surface_sorensen[match(int_vec, surface_ordinal_vec)]

  # # Since metrics are only defined within intervals, discard the ends.
  # # In this case, only the upper ends remain to be discarded.
  # discard_surface <- which(rowSums(apply(surface_ordinal, 2, function(x) x == n_gradient_steps)) > 0)
  # surface_area[discard_surface] <- 0

  # Collect variables
  out_df <- cbind(surface_var, mean_suitability = surface_mean_suitability,
                  n_limits = surface_n_limits, n_niches = surface_n_niches,
                  sorensen = surface_sorensen, inhull = inhyperhull)
  # return(list(data = out_df, r_sor = r_sorensen))
  return(out_df)
}



getTransitionSpp <- function(metrics_df, physeq, models, thresholds, rows=NULL, ncores=1,
                             unlist_if_single_row=TRUE) {
  pred_vars <- c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality", "k")
  if(!identical(pred_vars[1:4], names(metrics_df)[1:4])) {
    stop("Invalid input for `metrics_df`. First 4 columns must be 'mat_celsius', 'temp_seasonality',
    'map_mm', and 'prec_seasonality'.")
  }
  allequal <- function(x) {
    length(unique(x))==1
  }
  if(!allequal(apply(metrics_df[,1:4], 2, function(x) length(unique(x))))) {
    stop("Invalid input for `metrics_df`. All 4 main climate variables must have the same number of steps.")
  }
  n_gradient_steps <- length(unique(metrics_df[,1]))
  obs_env <- as(sample_data(physeq), "data.frame")[, pred_vars]
  obs_env <- obs_env[complete.cases(obs_env),]

  n_taxa <- length(taxa_names(physeq))
  obs_env <- as(sample_data(physeq), "data.frame")[, pred_vars]
  obs_env <- obs_env[complete.cases(obs_env),]

  quantize <- function(x, n) {
    seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length.out=n)
  }
  surface_var <- expand.grid(
    quantize(obs_env[,"mat_celsius"], n_gradient_steps),
    quantize(obs_env[,"temp_seasonality"], n_gradient_steps),
    quantize(obs_env[,"map_mm"], n_gradient_steps),
    quantize(obs_env[,"prec_seasonality"], n_gradient_steps)
  )
  surface_var[,5] <- yasso_k(surface_var[,1], surface_var[,3])
  names(surface_var) <- pred_vars

  surface_var_scaled <- scaleAndExpandPredictors(surface_var, physeq, models)

  metrics_ord <- cbind(
    match(metrics_df$mat_celsius, quantize(obs_env[,"mat_celsius"], n_gradient_steps)),
    match(metrics_df$temp_seasonality, quantize(obs_env[,"temp_seasonality"], n_gradient_steps)),
    match(metrics_df$map_mm, quantize(obs_env[,"map_mm"], n_gradient_steps)),
    match(metrics_df$prec_seasonality, quantize(obs_env[,"prec_seasonality"], n_gradient_steps))
  )
  colnames(metrics_ord) <- pred_vars[1:4]
  if(any(is.na(metrics_ord))) {
    stop("Some values in `metrics_df` do not align with the climate ranges in `physeq`.
         Are you sure you used the correct input for `physeq`?")
  }
  metrics_ord_vec <- apply(metrics_ord, 1, paste, collapse='.')

  registerDoParallel(ncores)
  surface_suitability <- foreach(i = seq_along(taxa_names(physeq)), .combine=cbind) %dopar% {
    tryCatch(
      as.vector(predict(models[[i]], newx = as.matrix(surface_var_scaled), type = "response", s = "lambda.min")),
      error = function(e) return(NA))
  }
  stopImplicitCluster()
  surface_presence <- do.call(cbind, lapply(seq_along(thresholds),
                                            function(i) surface_suitability[,i] >= thresholds[[i]]))

  if(is.null(rows)) rows <- which.max(metrics_df$sorensen_std)
  transition_spp_ind <- list()
  for(j in 1:length(rows)) {
    r <- rows[j]
    i1 <- metrics_ord[r, 1]
    i2 <- metrics_ord[r, 2]
    i3 <- metrics_ord[r, 3]
    i4 <- metrics_ord[r, 4]
    this_ord <- c(i1, i2, i3, i4)
    presence_at_row <- surface_presence[r,]
    offsets <- rbind(
      c(-1, 0, 0, 0),
      c(0, -1, 0, 0),
      c(0, 0, -1, 0),
      c(0, 0, 0, -1),
      c(1, 0, 0, 0),
      c(0, 1, 0, 0),
      c(0, 0, 1, 0),
      c(0, 0, 0, 1)
    )
    neighbors_ord <- sweep(offsets, 2, this_ord, '+')
    neighbors_ord <- neighbors_ord[which(apply(neighbors_ord, 1, function(x) all(x %in% seq_len(n_gradient_steps)))),]
    neighbors_vec <- apply(neighbors_ord, 1, paste, collapse='.')
    neighbors_ind <- which(metrics_ord_vec %in% neighbors_vec)
    neighbors_val <- foreach(k = neighbors_ind) %do% {
      return(surface_presence[k,])
    }
    crossings <- lapply(neighbors_val, function(x) x==0 & presence_at_row==1) # "inner" boundary
    isboundary <- rowSums(do.call(cbind, crossings)) > 0
    transition_spp_ind[[j]] <- which(isboundary)
  }
  if(length(rows==1) & unlist_if_single_row) {
    return(unlist(lapply(transition_spp_ind, function(x) taxa_names(physeq)[x])))
  } else {
    return(lapply(transition_spp_ind, function(x) taxa_names(physeq)[x]))
  }
}


projectMetricToRaster <- function(data, metric, r_climate, pvals=NULL,
                                  override_step_warning=FALSE) {
  if(!is.null(pvals)) {
    if(length(pvals) != nrow(data)) {
      warning("The length of `pvals` must match the number of rows in `data`.")
    }
  }
  if(!metric %in% colnames(data)) {
    warning("`metric` not among column names of `data`. Returning NULL")
    return(NULL)
  }
  pred_vars <- c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality", "k")
  step_sizes <- apply(data[,1:4], 2, function(x) min(diff(unique(x))))
  df_climate <- as.matrix(r_climate)[,pred_vars]
  # int0 <- cbind(
  #   findInterval(df_climate[,1], sort(unique(data[,1]))),
  #   findInterval(df_climate[,2], sort(unique(data[,2]))),
  #   findInterval(df_climate[,3], sort(unique(data[,3]))),
  #   findInterval(df_climate[,4], sort(unique(data[,4])))
  # )
  # int_vec0 <- apply(int0, 1, paste0, collapse=".")
  # int_vec0[grep("NA", int_vec0)] <- NA
  int <- cbind(
    findInterval(df_climate[,1] + step_sizes[1]/2, sort(unique(data[,1]))),
    findInterval(df_climate[,2] + step_sizes[2]/2, sort(unique(data[,2]))),
    findInterval(df_climate[,3] + step_sizes[3]/2, sort(unique(data[,3]))),
    findInterval(df_climate[,4] + step_sizes[4]/2, sort(unique(data[,4])))
  )
  int_vec <- apply(int, 1, paste0, collapse=".")
  int_vec[grep("NA", int_vec)] <- NA

  allequal <- function(x) {
    length(unique(x))==1
  }
  if(!override_step_warning) {
    if(!allequal(apply(data[,1:4], 2, function(x) length(unique(x))))) {
      warning("First 4 columns of data, in which we expect the climate variables, do not have
         the same number of steps. Returning NULL.")
      return(NULL)
    }
  }
  n_steps <- length(unique(data[,1]))
  surface_ordinal <- expand.grid(
    seq_len(n_steps),
    seq_len(n_steps),
    seq_len(n_steps),
    seq_len(n_steps)
  )
  surface_ordinal_vec <- apply(surface_ordinal, 1, paste, collapse='.')

  # # Get the total area within each climate voxel
  # area <- as.vector(area(r_climate))
  # int_area <- cbind.data.frame(int = int_vec, area) %>%
  #   group_by(int) %>%
  #   summarise(area = sum(area)) %>%
  #   mutate(int = as.character(int))
  # surface_area <- int_area$area[match(surface_ordinal_vec, int_area$int)]
  # surface_area[which(is.na(surface_area))] <- 0

  # Make raster map(s)
  r_metric <- raster(r_climate)
  r_metric[] <- data[match(int_vec, surface_ordinal_vec), metric]
  if(is.null(pvals)) {
    return(list(r_metric))
  } else {
    r_pvals <- raster(r_climate)
    r_pvals[] <- pvals[match(int_vec, surface_ordinal_vec)]
    return(list(r_metric, r_pvals))
  }
}

projectMetricToRaster2 <- function(data, metric, r_climate,
                                   terms_nonquadratic = c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality"),
                                   pvals=NULL, override_step_warning=FALSE,
                                   constrain_to_var_ranges=TRUE) {
  if(!is.null(pvals)) {
    if(length(pvals) != nrow(data)) {
      warning("The length of `pvals` must match the number of rows in `data`.")
    }
  }
  if(!metric %in% colnames(data)) {
    warning("`metric` not among column names of `data`. Returning NULL.")
    return(NULL)
  }
  if(!all(terms_nonquadratic %in% colnames(data))) {
    warning("Not all terms are among column names of `data`. Returning NULL.")
    return(NULL)
  }

  geo_df <- as.data.frame(cbind(coordinates(r_climate), as.matrix(r_climate)[,terms_nonquadratic]))
  names(geo_df)[1:2] <- c("lon", "lat")
  geo_df <- geo_df[complete.cases(geo_df),]

  step_sizes <- apply(data[,terms_nonquadratic], 2, function(x) min(diff(unique(x))))
  if(constrain_to_var_ranges) {
    int <- do.call(cbind, lapply(terms_nonquadratic, function(x) {
    out <- findInterval(geo_df[,x] + step_sizes[x]/2, sort(unique(data[,x])))
    out[which(geo_df[,x] < min(data[,x]) | geo_df[,x] > max(data[,x]))] <- NA
      return(out)
    }))
  } else {
    int <- do.call(cbind, lapply(terms_nonquadratic, function(x) {
      findInterval(geo_df[,x] + step_sizes[x]/2, sort(unique(data[,x])))
    }))
  }
  int_vec <- apply(int, 1, paste0, collapse=".")
  int_vec[grep("NA", int_vec)] <- NA

  allequal <- function(x) {
    length(unique(x))==1
  }
  if(!override_step_warning) {
    if(!allequal(apply(data[,terms_nonquadratic], 2, function(x) length(unique(x))))) {
      warning("The climate variables in `data` do not have
         the same number of steps as eqch other. Returning NULL. Ignore this warning
              with `override_step_warning = TRUE`.")
      return(NULL)
    }
  }
  n_steps <- length(unique(data[,terms_nonquadratic[1]]))
  surface_ordinal <- do.call(tidyr::expand_grid, lapply(1:length(terms_nonquadratic), function(i) 1:n_steps))
  surface_ordinal_vec <- apply(surface_ordinal, 1, paste, collapse='.')

  # # Get the total area within each climate voxel
  # area <- as.vector(area(r_climate))
  # int_area <- cbind.data.frame(int = int_vec, area) %>%
  #   group_by(int) %>%
  #   summarise(area = sum(area)) %>%
  #   mutate(int = as.character(int))
  # surface_area <- int_area$area[match(surface_ordinal_vec, int_area$int)]
  # surface_area[which(is.na(surface_area))] <- 0

  # Make raster map(s)
  r_metric <- rasterFromXYZ(cbind(geo_df$lon, geo_df$lat, data[match(int_vec, surface_ordinal_vec), metric]))
  # r_metric <- raster(r_climate)
  # r_metric[] <- data[match(int_vec, surface_ordinal_vec), metric]
  if(is.null(pvals)) {
    return(list(r_metric))
  } else {
    r_pvals <- raster(r_climate)
    r_pvals[] <- pvals[match(int_vec, surface_ordinal_vec)]
    return(list(r_metric, r_pvals))
  }
}

projectMetricToMap <- function(data, metric, r_climate, terms_nonquadratic = c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality"),
                               borders_sp=NULL, scale_breaks=NULL, scale_text="", title_text="",
                               pvals=NULL, viridis_opt="B", constrain_to_var_ranges=TRUE,
                               mask=NULL) {
  # r_out <- projectMetricToRaster(data, metric, r_climate, pvals)
  r_out <- projectMetricToRaster2(data, metric, r_climate, terms_nonquadratic, pvals,
                                  constrain_to_var_ranges=constrain_to_var_ranges)
  r_metric <- r_out[[1]]
  if(!is.null(mask)) {
    r_metric <- raster::mask(r_metric, mask)
  }

  plotting_data <- data.frame(
    coordinates(r_metric),
    z = as.vector(r_metric)
  )
  p <- ggplot(plotting_data) +
    geom_tile(aes(x, y, fill=z)) +
    scale_fill_viridis(scale_text, option=viridis_opt, na.value="transparent",
                       breaks=if(is.null(scale_breaks)) waiver() else scale_breaks) +
    ggtitle(title_text) +
    theme_custom_map +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom")
  if(!is.null(borders_sp)) {
    borders <- st_as_sf(borders_sp)
    p <- p + geom_sf(data=st_as_sf(borders), size=0.1, col="black", fill=alpha("white", 0)) +
      coord_sf(xlim=c(-170,-55), ylim=c(18,72))
  }
  # if(!is.null(pvals)) {
  #   r_pvals <- r_out[[2]]
  #   plotting_data_sig <- data.frame(
  #     coordinates(r_pvals),
  #     pvals = as.vector(r_pvals)
  #   ) %>%
  #     mutate(keep_in_thinned = (x%%1==0.25) & (y%%1==0.25)) %>% # Thins coordinates by 1/6th across a regular grid
  #     dplyr::filter(!is.na(pvals)) %>%
  #     dplyr::filter(pvals==FALSE, keep_in_thinned==TRUE)
  #     # dplyr::filter(pvals==FALSE)
  #   p <- p + geom_point(data=plotting_data_sig, aes(x, y), col="grey80", shape=4, size=0.2)
  # }
  return(p)
}


plot_2d_nichemetric <- function(data, metric, xvar, yvar, xlab="", ylab="",
                                maxptsize = 5, inhull=TRUE, area=TRUE, viridis_opt="B") {
  dataslice <- data
  all_vars <- c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality")
  if(any(!c(xvar, yvar) %in% all_vars)) {
    stop("`xvar` and `yvar` must be among: ", paste(all_vars, collapse=", "))
  }
  const_vars <- all_vars[which(!all_vars %in% c(xvar, yvar))]
  if(inhull) {
    dataslice <- dataslice[which(dataslice$inhull==TRUE),]
  }
  range_x <- c(min(dataslice[,xvar]), max(dataslice[,xvar]))
  range_y <- c(min(dataslice[,yvar]), max(dataslice[,yvar]))
  for(v in const_vars) {
    dataslice <- dataslice[which(dataslice[,v] == quantile(data[,v], 0.49)),]
  }
  range_metric <- c(0, max(dataslice[,metric], na.rm=TRUE))
  if(area) range_area <- c(0, max(dataslice$area, na.rm=TRUE))
  g <- ggplot(dataslice, aes_string(x=xvar, y=yvar, fill=metric)) +
    geom_tile() +
    scale_fill_viridis(limits=range_metric, option=viridis_opt) +
    guides(fill="none") +
    xlab(xlab) + ylab(ylab) +
    coord_cartesian(xlim=range_x, ylim=range_y) +
    theme_custom
  if(area) {
    g <- g + geom_point(aes_string(size="area"), pch=21, col="white") +
      scale_size_continuous(range=c(0,maxptsize), limits=range_area) +
      # scale_color_viridis(limits=range_metric) +
      # guides(col="none", size="none") +
      guides(size="none")
  }
  g
}

plotAll2DNicheMetric <- function(data, metric, scale_text="", inhull=TRUE, area=TRUE, viridis_opt="B") {
  library(cowplot)
  color_legend <- cowplot::get_legend(ggplot(data, aes_string(x="mat_celsius", y="map_mm", col=metric)) + geom_point() + scale_color_viridis(scale_text, option=viridis_opt))
  g <- gridExtra::arrangeGrob(
    plot_2d_nichemetric(data, metric, "mat_celsius", "temp_seasonality", inhull=inhull, area=area, viridis_opt=viridis_opt),
    plot_2d_nichemetric(data, metric, "mat_celsius", "map_mm", inhull=inhull, area=area, viridis_opt=viridis_opt),
    plot_2d_nichemetric(data, metric, "mat_celsius", "prec_seasonality", inhull=inhull, area=area, viridis_opt=viridis_opt),
    plot_2d_nichemetric(data, metric, "temp_seasonality", "map_mm", inhull=inhull, area=area, viridis_opt=viridis_opt),
    plot_2d_nichemetric(data, metric, "temp_seasonality", "prec_seasonality", inhull=inhull, area=area, viridis_opt=viridis_opt),
    plot_2d_nichemetric(data, metric, "map_mm", "prec_seasonality", inhull=inhull, area=area, viridis_opt=viridis_opt),
    color_legend,
    layout_matrix = rbind(c(1, NA, 7),
                          c(2, 4, NA),
                          c(3, 5, 6)),
    left = textGrob("     Prec. seasonality                           MAP                               Temp. seasonality", rot=90),
    bottom = textGrob("        MAT                                    Temp. seasonality                                 MAP")
  )
  return(g)
}

permuteClimate <- function(physeq, randomseed=NULL) {
  sampledata <- as(sample_data(physeq), "data.frame")
  if(!is.null(randomseed)) set.seed(randomseed)
  permutation <- sample(seq(nrow(sampledata)), replace=TRUE)
  newsampledata <- sampledata[permutation,]
  rownames(newsampledata) <- rownames(sampledata)
  out <- physeq
  sample_data(out) <- newsampledata
  return(out)
}

permuteOccurrence <- function(physeq, randomseed=NULL, ncores=1) {
  otu <- as(otu_table(physeq), "matrix")
  if(!is.null(randomseed)) set.seed(randomseed)
  registerDoParallel(ncores)
  otu_permuted <- foreach(i = seq_along(taxa_names(physeq)), .combine = cbind) %dopar% {
    permutation <- sample(seq(nrow(otu)), replace=TRUE)
    otu[permutation, i]
  }
  stopImplicitCluster()
  colnames(otu_permuted) <- colnames(otu)
  rownames(otu_permuted) <- rownames(otu)
  out <- phyloseq(
    otu_table(otu_permuted, taxa_are_rows=FALSE),
    sample_data(physeq),
    tax_table(physeq)
  )
  return(out)
}

bootstrapSamples <- function(physeq, randomseed=NULL) {
  otu <- as(otu_table(physeq), "matrix")
  sampledata <- as(sample_data(physeq), "data.frame")
  if(!is.null(randomseed)) set.seed(randomseed)
  permutation <- sample(seq(nrow(sampledata)), replace=TRUE)
  new_otu <- otu[permutation,]
  new_sampledata <- sampledata[permutation,]
  rownames(new_otu) <- paste0("Sample_", seq(nrow(new_otu)))
  rownames(new_sampledata) <- paste0("Sample_", seq(nrow(new_sampledata)))
  out <- phyloseq(otu_table(new_otu, taxa_are_rows = FALSE),
                  sample_data(new_sampledata),
                  tax_table(physeq))
  return(out)
}

findThresholds <- function(physeq, models, terms = c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality"),
                           s=c("min", "1se"), opt=c("eq","maxtss"), ncores=1) {
  message("findThresholds: Started at ", Sys.time())
  s <- match.arg(s, c("min", "1se"))
  opt <- match.arg(opt, c("eq","maxtss"))
  s <- switch(s, "min" = "lambda.min", "1se" = "lambda.1se")

  covars <- covarNamesFromLognets(models)
  pred_data <- as(sample_data(physeq)[,terms], "data.frame")
  pred_data_std <- scale(pred_data, center=TRUE, scale=TRUE)
  ind_cols <- match(covars, colnames(pred_data_std))
  pred_data_std <- pred_data_std[,ind_cols]

  registerDoParallel(ncores)
  thres <- foreach(i = seq_along(taxa_names(physeq)), .combine = c) %dopar% {
    if(identical(models[[i]], NA)) return(NA)
    sp <- taxa_names(physeq)[i]
    y01 <- as.vector(otu_table(physeq)[,sp]) > 0
    lognet_pred <- predict(models[[i]], newx = pred_data_std, type="response", s = s)
    t_grid <- seq(min(lognet_pred, na.rm=TRUE), max(lognet_pred, na.rm=TRUE),
                  by=(max(lognet_pred, na.rm=TRUE)-min(lognet_pred, na.rm=TRUE))/100)
    eval <- dismo::evaluate(p = lognet_pred[which(y01)], a = lognet_pred[which(!y01)], tr = t_grid)
    # plot(eval@TPR - eval@TNR)

    if(identical(opt, "eq")) {
      # Threshold to minimize diff b/w sensitivity and specificity
      t <- eval@t[which.min(abs(eval@TPR - eval@TNR))]
      return(t)
    }
    if(identical(opt, "maxtss")) {
      # threshold to maximize true skill statistic
      t <- eval@t[which.max(eval@TPR + eval@TNR)]
      return(t)
    }
  }
  stopImplicitCluster()
  message("findThresholds: Ended at ", Sys.time())
  return(thres)
}

# alpha = 0: ridge regression
# alpha = 1: lasso
fitLognets <- function(physeq, alpha=0, nfolds=10, ncores=1, randomseed=NULL,
                       terms=c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality",
                               "mat_celsius_2", "map_mm_2")) {
  message("fitLognets: Started at ", Sys.time())
  if(any(!(terms %in% colnames(sample_data(physeq))))) stop("Not all terms are in physeq")

  # Reorder terms to match the order they appear in physeq
  varnames <- colnames(sample_data(physeq))
  terms <- terms[order(match(terms, varnames))]

  predictors <- as(sample_data(physeq)[,terms], "data.frame")
  predictors_std <- scale(predictors, center=TRUE, scale=TRUE)
  complete_records <- which(complete.cases(predictors_std))
  registerDoParallel(ncores)
  models <- foreach(i = seq_along(taxa_names(physeq))) %dopar% {
    y01 <- as.vector(otu_table(physeq)[,i]) > 0
    if(!is.null(randomseed)) set.seed(randomseed)
    tryCatch(
      cv.glmnet(x=as.matrix(predictors_std)[complete_records,],
                y=y01[complete_records], family = "binomial",
                alpha = alpha, type.measure = "default", standardize=FALSE,
                nfolds = nfolds),
      error = function(e) return(NA)
    )
  }
  stopImplicitCluster()
  names(models) <- taxa_names(physeq)
  message("fitLognets: Ended at ", Sys.time())
  return(models)
}

getGenusSummary <- function(physeq, n_top_genera = 20) {
  ps_genus <- tax_glom(physeq, "genus")
  id_taxa_ind <- !grepl("unidentified", as.vector(tax_table(ps_genus)[,"genus"]))
  ps_genus <- prune_taxa(id_taxa_ind, ps_genus)
  top_otus_gen <- names(sort(taxa_sums(ps_genus), TRUE))[1:min(length(taxa_sums(ps_genus)), n_top_genera)]
  ps_genus <- prune_taxa(top_otus_gen, ps_genus)
  otu_counts_by_gen <- data.frame(table(as.vector(tax_table(physeq)[,"genus"])))
  names(otu_counts_by_gen) <- c("genus", "n_otu")
  out <- cbind.data.frame(
    as(tax_table(ps_genus), "matrix"),
    mean_rel_abund = taxa_sums(ps_genus) / length(sample_sums(ps_genus))) %>%
    left_join(otu_counts_by_gen) %>%
    arrange(desc(n_otu))
  out[,-1]
}


# pred_vars <- c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality", "k")
# if(!identical(pred_vars[1:4], names(metrics_df)[1:4])) {
#   stop("Invalid input for `metrics_df`. First 4 columns must be 'mat_celsius', 'temp_seasonality',
#     'map_mm', and 'prec_seasonality'.")
# }
# allequal <- function(x) {
#   length(unique(x))==1
# }
# if(!allequal(apply(metrics_df[,1:4], 2, function(x) length(unique(x))))) {
#   stop("Invalid input for `metrics_df`. All 4 main climate variables must have the same number of steps.")
# }
# n_gradient_steps <- length(unique(metrics_df[,1]))
# obs_env <- as(sample_data(physeq), "data.frame")[, pred_vars]
# obs_env <- obs_env[complete.cases(obs_env),]
#
# n_taxa <- length(taxa_names(physeq))
# obs_env <- as(sample_data(physeq), "data.frame")[, pred_vars]
# obs_env <- obs_env[complete.cases(obs_env),]
#
# quantize <- function(x, n) {
#   seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length.out=n)
# }
# surface_var <- expand.grid(
#   quantize(obs_env[,"mat_celsius"], n_gradient_steps),
#   quantize(obs_env[,"temp_seasonality"], n_gradient_steps),
#   quantize(obs_env[,"map_mm"], n_gradient_steps),
#   quantize(obs_env[,"prec_seasonality"], n_gradient_steps)
# )
# surface_var[,5] <- yasso_k(surface_var[,1], surface_var[,3])
# names(surface_var) <- pred_vars
#
# # Bring in climate raster
# df_climate <- as.matrix(r_climate)[,pred_vars]
# int <- cbind(
#   findInterval(df_climate[,1], quantize(obs_env[,1], n_gradient_steps)),
#   findInterval(df_climate[,2], quantize(obs_env[,2], n_gradient_steps)),
#   findInterval(df_climate[,3], quantize(obs_env[,3], n_gradient_steps)),
#   findInterval(df_climate[,4], quantize(obs_env[,4], n_gradient_steps))
# )
# int_vec <- apply(int, 1, paste0, collapse=".")
# int_vec[grep("NA", int_vec)] <- NA
#
# # Get the total area within each climate voxel
# area <- as.vector(area(r_climate))
# int_area <- cbind.data.frame(int = int_vec, area) %>%
#   group_by(int) %>%
#   summarise(area = sum(area)) %>%
#   mutate(int = as.character(int))
# surface_area <- int_area$area[match(surface_ordinal_vec, int_area$int)]
# surface_area[which(is.na(surface_area))] <- 0
#
# # Calculate Sorensen-like climate sensitivity
# surface_sorensen <- surface_n_limits / surface_n_niches
# # As raster map:
# r_sorensen <- raster(r_climate)
# r_sorensen[] <- surface_sorensen[match(int_vec, surface_ordinal_vec)]
#
# # Since metrics are only defined within intervals, discard the ends.
# # In this case, only the upper ends remain to be discarded.
# discard_surface <- which(rowSums(apply(surface_ordinal, 2, function(x) x == n_gradient_steps)) > 0)
# surface_area[discard_surface] <- 0

# Split violin plot by user jan-glx on StackOverflow https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ...,
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

getNicheLimitsOnGradient <- function(
    physeq, models, thresholds, gradientvar, n_gradient_steps=100,
    constant_values=NULL, s=c("min", "1se"),
    terms = c("mat_celsius", "temp_seasonality", "map_mm",
                  "prec_seasonality", "soilInCaClpH", "nitrogenPercent", "organicCPercent"),
    progressbar = TRUE) {
  s <- match.arg(s, c("min", "1se"))
  s <- switch(s, "min" = "lambda.min", "1se" = "lambda.1se")
  if(!gradientvar %in% terms) {
    stop("`gradientvar` must be one of `terms`")
  }
  n_taxa <- length(taxa_names(physeq))
  obs_env <- as(sample_data(physeq), "data.frame")[, terms]
  obs_env <- obs_env[complete.cases(obs_env),]
  env_medians <- apply(obs_env, 2, function(x) median(x, na.rm=TRUE))
  if(is.null(constant_values)) {
    env_constant <- env_medians
  } else {
    if(length(constant_values) != length(terms)) {
      stop("`constant_values` must be same length as `terms`.
           To use the median of observed values for a variable, enter NA for
           that variable.")
    }
    env_constant <- constant_values
    env_constant[which(is.na(constant_values))] <- env_medians[which(is.na(constant_values))]
    # env_constant <- c(env_constant, yasso_k(env_constant[1], env_constant[3]))
    names(env_constant) <- terms
  }
  message("Using as constant values: ", paste0(terms, " ", env_constant, ", "))
  surface_var <- data.frame(
    seq(min(obs_env[,gradientvar], na.rm=TRUE),
        max(obs_env[,gradientvar], na.rm=TRUE),
        length.out = n_gradient_steps)
  )
  names(surface_var) <- gradientvar
  for(p in terms) {
    if(!p %in% colnames(surface_var)) {
      surface_var[,p] <- env_constant[p]
    }
  }
  # surface_var <- dplyr::select(surface_var, all_of(terms)) # reorder columns

  covars <- covarNamesFromLognets(models)
  surface_var <- includeQuadraticTerms(surface_var)
  surface_var <- surface_var[,match(covars, colnames(surface_var))]
  predictors <- as(sample_data(physeq)[,terms], "data.frame")
  predictors <- predictors[,match(covars, colnames(predictors))]
  surface_var_scaled <- scaleToReference(surface_var, predictors)

  surface_var_presence <- matrix(NA, nrow = n_gradient_steps, ncol = n_taxa,
                                 dimnames = list(NULL, taxa_names(physeq)))
  if(progressbar) pb <- txtProgressBar(min = 0, max = n_taxa, style=3)
  for(i in 1:n_taxa) {
    sp <- taxa_names(physeq)[i]
    suitability <- tryCatch(
      as.vector(predict(models[[i]], newx = as.matrix(surface_var_scaled), type = "response", s = s)),
      error = function(e) return(NA))
    surface_var_presence[,sp] <- suitability >= thresholds[[i]]
    if(progressbar) setTxtProgressBar(pb, i)
  }
  if(progressbar) close(pb)
  edge_mat_lwr <- edge_mat_upr <- rep(NA, n_taxa)
  for(i in 1:n_taxa) {
    edge_mat_lwr[i] <- suppressWarnings(min(surface_var[,gradientvar][surface_var_presence[,i]]))
    edge_mat_upr[i] <- suppressWarnings(max(surface_var[,gradientvar][surface_var_presence[,i]]))
  }
  edge_mat_lwr[is.infinite(edge_mat_lwr)] <- NA
  edge_mat_upr[is.infinite(edge_mat_upr)] <- NA

  data.frame(
    i = rep(1:n_taxa, 2),
    sp = rep(taxa_names(physeq), 2),
    edge = c(rep("lower", n_taxa),
             rep("upper", n_taxa)),
    value = c(edge_mat_lwr, edge_mat_upr)
  ) %>% arrange(i, edge)
}


# findThresholdsCV <- function(physeq, models, s=c("min", "1se"), opt=c("eq","maxacc"), ncores=1,
#                              randomseed=NULL) {
#   message("Started at ", Sys.time())
#   s <- match.arg(s, c("min", "1se"))
#   opt <- match.arg(opt, c("eq","maxacc"))
#   s <- switch(s, "min" = "lambda.min", "1se" = "lambda.1se")
#   pred_vars <- c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality", "k")
#
#   covars <- rownames(coef(models[[1]]))
#   if(any(grepl(":", covars))) inclInteractions <- TRUE else inclInteractions <- FALSE
#   if(any(grepl("_2$", covars))) inclQuadraticTerms <- TRUE else inclQuadraticTerms <- FALSE
#   message("Based on covariate names of the first model, assuming that inclQuadraticTerms = ", inclQuadraticTerms,
#           ", and inclInteractions = ", inclInteractions)
#
#   if(inclQuadraticTerms && inclInteractions) {
#     pred_data <- expandTerms(as(sample_data(physeq)[,pred_vars], "data.frame"), exclude="k")
#   } else if(inclQuadraticTerms && !inclInteractions) {
#     pred_data <- includeQuadraticTerms(as(sample_data(physeq)[,pred_vars], "data.frame"), exclude="k")
#   } else if(!inclQuadraticTerms && inclInteractions) {
#     pred_data <- includeInteractions(as(sample_data(physeq)[,pred_vars], "data.frame"), exclude="k")
#   } else {
#     pred_data <- as(sample_data(physeq)[,pred_vars], "data.frame")
#   }
#   pred_data_std <- scale(pred_data, center=TRUE, scale=TRUE)
#
#   registerDoParallel(ncores)
#   thres <- foreach(i = seq_along(taxa_names(physeq)), .combine = c) %dopar% {
#     y01 <- as.vector(otu_table(physeq)[,i]) > 0
#     t_grid <- seq(min(lognet_pred), max(lognet_pred),
#                   by=(max(lognet_pred)-min(lognet_pred))/20)
#     # ...
#     if(!is.null(randomseed)) set.seed(randomseed)
#     folds <- createFolds(seq_along(sample_names(physeq)), k = 5, list = TRUE, returnTrain = FALSE)
#     for(k in 1:5) {
#       train_ind <- unlist(folds[-k])
#       test_ind <- unlist(folds[k])
#       pred_data_std_k <- pred_data_std[train_ind,]
#       y01_k <- y01[train_ind]
#       lognet_pred <- predict(models[[i]], newx = pred_data_std_k, type="response", s = s)
#       eval <- dismo::evaluate(p = lognet_pred[which(y01)], a = lognet_pred[which(!y01)], tr = t_grid)
#     }
#     # plot(eval@TPR - eval@TNR)
#
#     if(identical(opt, "eq")) {
#       # Threshold to minimize diff b/w sensitivity and specificity
#       t <- eval@t[which.min(abs(eval@TPR - eval@TNR))]
#       return(t)
#     }
#     if(identical(opt, "maxacc")) {
#       # threshold to maximize accuracy
#       t <- eval@t[which.max(eval@TPR + eval@TNR)]
#       return(t)
#     }
#   }
#   stopImplicitCluster()
#   message("Ended at ", Sys.time())
#   return(thres)
# }
