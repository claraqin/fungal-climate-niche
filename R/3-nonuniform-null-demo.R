# R code for Qin et al. on fungal niches and climate sensitivity
# 3-nonuniform-null-demo: Justify the null model by showing its nonuniformity (for Appendix)

# Code originally from "3-nonuniform-null-demo.R"

# What does it take to achieve a uniform / random distribution of
# niche edge density and Sorensen sensitivity across climate space?

# Get subset of 500 taxa
set.seed(10101)
n_taxa <- 500
spp_subset_ind <- sort(sample(1:length(taxa_names(neon_dob_prevalent)), size=n_taxa))
neon_dob_prevalent_sppsub <- prune_taxa(taxa_names(neon_dob_prevalent)[spp_subset_ind],
                                        neon_dob_prevalent)

# Define variables
terms <- c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality",
           "mat_celsius_2", "map_mm_2", "soilInCaClpH", "nitrogenPercent", "organicCPercent")


# Method 1: Permute climate -------------------------------------------------------

# Each iteration fails to achieve a uniform null distribution; however, each
# permutation appears to create a different distribution, so the iterations
# may create a uniform distribution in aggregate.

neon_dob_permuted <- permuteClimate(neon_dob_prevalent_sppsub, randomseed=1010103)
predictors <- as(sample_data(neon_dob_permuted)[,terms], "data.frame")
predictors_std <- scale(predictors, center=TRUE, scale=TRUE)
complete_records <- which(complete.cases(predictors_std))

lognets_permuted <-  foreach(i = seq_along(spp_subset_ind)) %dopar% {
  y01 <- as.vector(otu_table(neon_dob_permuted)[,i]) > 0
  set.seed(1010101)
  tryCatch(cv.glmnet(as.matrix(predictors_std)[complete_records,],
                     y01[complete_records], family = "binomial",
                     alpha = 0, type.measure = "default", standardize=FALSE),
           error = function(e) return(NA))
}

thresholds_permuted <- findThresholds(neon_dob_permuted, lognets_permuted, pred_vars = terms)

metrics_2d_permuted <- getNicheMetrics2D(neon_dob_permuted, lognets_permuted, thresholds_permuted, "mat_celsius", "map_mm",
                                         r_present_northam, n_bins_per_axis = 10, ncores=8,
                                         pred_vars = terms)

p_permuted <- plotNicheMetrics2D(neon_dob_permuted, metrics_2d_permuted, pred_vars = terms)



# Method 2: Enforce uniform sampling intensity across climate space --------

# Each iteration appears to create a somewhat random distribution.
# It is possible that the iterations may create a uniform distribution in aggregate.

neon_dob_runif <- neon_dob_prevalent_sppsub
set.seed(1010104)
sample_data(neon_dob_runif) %>%
  as("data.frame") %>%
  mutate(mat_celsius = runif(length(mat_celsius), min(mat_celsius), max(mat_celsius)),
         map_mm = runif(length(map_mm), min(map_mm), max(map_mm))) %>%
  mutate(mat_celsius_2 = mat_celsius^2,
         map_mm_2 = map_mm^2) ->
  sample_data(neon_dob_runif)
predictors <- as(sample_data(neon_dob_runif)[,terms], "data.frame")
predictors_std <- scale(predictors, center=TRUE, scale=TRUE)
complete_records <- which(complete.cases(predictors_std))

lognets_runif <-  foreach(i = seq_along(spp_subset_ind)) %dopar% {
  y01 <- as.vector(otu_table(neon_dob_runif)[,i]) > 0
  set.seed(1010101)
  tryCatch(cv.glmnet(as.matrix(predictors_std)[complete_records,],
                     y01[complete_records], family = "binomial",
                     alpha = 0, type.measure = "default", standardize=FALSE),
           error = function(e) return(NA))
}

thresholds_runif <- findThresholds(neon_dob_runif, lognets_runif, pred_vars = terms)

metrics_2d_runif <- getNicheMetrics2D(neon_dob_runif, lognets_runif, thresholds_runif, "mat_celsius", "map_mm",
                                      r_present_northam, n_bins_per_axis = 10, ncores=8, pred_vars = terms)

p_runif <- plotNicheMetrics2D(neon_dob_runif, metrics_2d_runif, pred_vars = terms)



# Method 3: Are soil variables necessary for randomness in the null model? ----------------------------

# After trying a few iterations, it seems like including soil variables
# does help each iteration to produce a more uniform distribution.

terms_climateonly <- c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality", "mat_celsius_2", "map_mm_2")

set.seed(1010103)
sample_data(neon_dob_runif) %>%
  as("data.frame") %>%
  mutate(mat_celsius = runif(length(mat_celsius), min(mat_celsius), max(mat_celsius)),
         map_mm = runif(length(map_mm), min(map_mm), max(map_mm))) %>%
  mutate(mat_celsius_2 = mat_celsius^2,
         map_mm_2 = map_mm^2) ->
  sample_data(neon_dob_runif)
predictors <- as(sample_data(neon_dob_runif)[,terms_climateonly], "data.frame")
predictors_std <- scale(predictors, center=TRUE, scale=TRUE)
complete_records <- which(complete.cases(predictors_std))

lognets_runif2 <-  foreach(i = seq_along(spp_subset_ind)) %dopar% {
  y01 <- as.vector(otu_table(neon_dob_runif)[,i]) > 0
  set.seed(1010101)
  tryCatch(cv.glmnet(as.matrix(predictors_std)[complete_records,],
                     y01[complete_records], family = "binomial",
                     alpha = 0, type.measure = "default", standardize=FALSE),
           error = function(e) return(NA))
}

thresholds_runif2 <- findThresholds(neon_dob_runif, lognets_runif2, pred_vars = terms_climateonly)

metrics_2d_runif2 <- getNicheMetrics2D(neon_dob_runif, lognets_runif2, thresholds_runif2, "mat_celsius", "map_mm",
                                       r_present_northam, n_bins_per_axis = 10, ncores=8, pred_vars = terms_climateonly)

p_runif2 <- plotNicheMetrics2D(neon_dob_runif, metrics_2d_runif2, pred_vars = terms_climateonly)



# Method 4: In aggregate, do iterations produce a uniform null? (Soil+climate covars) ---------------------

# Short answer: no.

aggregateNullIters2D <- function(physeq, terms, n_iters=10, n_taxa=100, n_bins_per_axis=10, ncores=8) {

  n_niche_list2d <- n_limits_list2d <- sorensen_list2d <- list()
  pb <- txtProgressBar(max=n_iters, style=3)
  message("Started at ", Sys.time())

  for(i in 1:n_iters) { # Takes ~15 minutes to run 100 iterations

    # Get random subset of n_taxa taxa
    set.seed(1010100+i)
    spp_subset_ind <- sort(sample(1:length(taxa_names(physeq)), size=n_taxa))
    physeq_sppsub <- prune_taxa(taxa_names(physeq)[spp_subset_ind],
                                   physeq)

    physeq_runif <- physeq_sppsub
    set.seed(1010100+i)
    sample_data(physeq_runif) %>%
      as("data.frame") %>%
      mutate(mat_celsius = runif(length(mat_celsius), min(mat_celsius), max(mat_celsius)),
             map_mm = runif(length(map_mm), min(map_mm), max(map_mm))) %>%
      mutate(mat_celsius_2 = mat_celsius^2,
             map_mm_2 = map_mm^2) ->
      sample_data(physeq_runif)
    predictors <- as(sample_data(physeq_runif)[,terms], "data.frame")
    predictors_std <- scale(predictors, center=TRUE, scale=TRUE)
    complete_records <- which(complete.cases(predictors_std))

    registerDoParallel(8)
    lognets_runif <-  foreach(i = seq_along(spp_subset_ind)) %dopar% {
      y01 <- as.vector(otu_table(physeq_runif)[,i]) > 0
      set.seed(1010101)
      tryCatch(cv.glmnet(as.matrix(predictors_std)[complete_records,],
                         y01[complete_records], family = "binomial",
                         alpha = 0, type.measure = "default", standardize=FALSE),
               error = function(e) return(NA))
    }
    stopImplicitCluster()

    thresholds_runif <- suppressMessages(findThresholds(physeq_runif, lognets_runif, pred_vars = terms))

    metrics_2d_runif <- suppressMessages(getNicheMetrics2D(physeq_runif, lognets_runif, thresholds_runif, "mat_celsius", "map_mm",
                                                           r_present_northam, n_bins_per_axis = n_bins_per_axis, ncores=ncores, pred_vars = terms, progressbar=FALSE))

    n_niche_list2d[[i]] <- metrics_2d_runif$n_niches
    n_limits_list2d[[i]] <- metrics_2d_runif$n_limits
    sorensen_list2d[[i]] <- metrics_2d_runif$sorensen

    setTxtProgressBar(pb, i)
  }
  message("Ended at ", Sys.time())
  close(pb)

  permutation_metrics_2d <- list(
    do.call(cbind, n_niche_list2d),
    do.call(cbind, n_limits_list2d),
    do.call(cbind, sorensen_list2d))
  names(permutation_metrics_2d) <- c("n_niche", "n_limits", "sorensen")

  n_niche_mean2 <- rowMeans(do.call(cbind, n_niche_list2d), na.rm=TRUE)
  n_limits_mean2 <- rowMeans(do.call(cbind, n_limits_list2d), na.rm=TRUE)
  sorensen_mean2 <- rowMeans(do.call(cbind, sorensen_list2d), na.rm=TRUE)
  n_niche_mean2[which(is.nan(n_niche_mean2))] <- NA
  n_limits_mean2[which(is.nan(n_limits_mean2))] <- NA
  sorensen_mean2[which(is.nan(sorensen_mean2))] <- NA

  niche_metrics_2d_runif <- cbind(metrics_2d_runif[,c(1:9,14:15)],
                                  n_niche_mean2, n_limits_mean2, sorensen_mean2)
  return(niche_metrics_2d_runif)
}

# niche_metrics_2d_runif_nonprevalent <- niche_metrics_2d_runif

# niche_metrics_2d_runif_prevalent10 <- niche_metrics_2d_runif

# niche_metrics_2d_runif_prevalent20 <- niche_metrics_2d_runif

niche_metrics_2d_runif <- aggregateNullIters2D(neon_dob_prevalent, terms)

p_niche_runif_prevalent10 <- plotNicheMetrics2D(neon_dob_runif, niche_metrics_2d_runif_prevalent10,
                                                site_points_in_plots = c(TRUE, TRUE, TRUE),
                                                metric1="n_niche_mean2", metric2="n_limits_mean2", metric3="sorensen_mean2")
ggsave("plots/2d_niche_nullrunif_prevalent10.png", p_niche_runif_prevalent10, device="png", width=3, height=7, units="in")


# 5. Plot null distributions by prevalence cutoff ----------------------------

neon_dob_agg <- readRDS("data/neon_dob_agg_v2.Rds")

# Subset to OTUs present in at least 20 sites
prevalence <- apply(otu_table(neon_dob_agg), 2, function(x) sum(x>0))
neon_dob_prevalent20 <- prune_taxa(prevalence >= 20, neon_dob_agg)

niche_metrics_2d_runif_nonprevalent <- aggregateNullIters2D(neon_dob_agg, terms)

niche_metrics_2d_runif_prevalent20 <- aggregateNullIters2D(neon_dob_prevalent20, terms)

p_niche_runif_nonprevalent <- plotNicheMetrics2D(neon_dob_runif, niche_metrics_2d_runif_nonprevalent,
                                                 site_points_in_plots = c(TRUE, TRUE, TRUE),
                                                 metric1="n_niche_mean2", metric2="n_limits_mean2", metric3="sorensen_mean2")
ggsave("plots/2d_niche_nullrunif_nonprevalent.png", p_niche_runif_nonprevalent, device="png", width=3, height=7, units="in")

p_niche_runif_prevalent20 <- plotNicheMetrics2D(neon_dob_runif, niche_metrics_2d_runif_prevalent20,
                                                site_points_in_plots = c(TRUE, TRUE, TRUE),
                                                metric1="n_niche_mean2", metric2="n_limits_mean2", metric3="sorensen_mean2")
ggsave("plots/2d_niche_nullrunif_prevalent20.png", p_niche_runif_prevalent20, device="png", width=3, height=7, units="in")
