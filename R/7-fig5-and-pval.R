# R code for Qin et al. on fungal niches and climate sensitivity
# 7-fig5-and-pval.R: Make Figure 5 (niche metrics in 5D) and P-val plot for Figure 4

# Code originally from "permutations.R" and "find_uniform_null.R"

neon_dob_wellfit <- readRDS("data/neon_dob_wellfit_v1.0.Rds")

terms_important <- c("mat_celsius", "temp_seasonality", "map_mm",
                     "mat_celsius_2", "map_mm_2", "soilInCaClpH", "organicCPercent")

var_ranges <- lapply(grep("_2$", terms_important, value=TRUE, invert=TRUE),
                     function(x) range(get_variable(neon_dob_wellfit, x), na.rm=TRUE))

n_niche_boot_list5d <- n_limits_boot_list5d <- n_limits2_boot_list5d <- sorensen_boot_list5d <- list()

n_niche_permuted_list5d <- n_limits_permuted_list5d <- n_limits2_permuted_list5d <- sorensen_permuted_list5d <- list()


# Generate bootstraps on empirical data (to quantify uncertainty) --------------------

# "data/boot_metrics_5d_wellfit.Rds" currently contains the results from iterations 1-40.

n_iters <- 50

pb <- txtProgressBar(max=n_iters, style=3)
message("Started at ", Sys.time())
for(i in 1:n_iters) { # Takes 5.6 minutes per iteration with 16 cores

  # t0 <- Sys.time()
  # Bootstrap sample
  neon_dob_boot <- bootstrapSamples(neon_dob_wellfit, randomseed=1010100+i)
  # t0.1 <- Sys.time()

  lognets_boot <- suppressMessages(fitLognets(neon_dob_boot, terms = terms_important, ncores=N_CORES, randomseed = 1010101))
  # t0.2 <- Sys.time()

  thresholds_boot <- suppressMessages(findThresholds(neon_dob_boot, lognets_boot, terms = terms_important, ncores = N_CORES))
  # t0.3 <- Sys.time()

  # metrics_boot_iter <- suppressMessages(getNicheMetrics5D(
  metrics_boot_iter <- suppressMessages(getNicheMetrics5D0(
    neon_dob_boot, lognets_boot, thresholds_boot, r_present_northam, inflate = 0.15, # Can constrain by hull later
    n_bins_per_axis = 8, # Changed from 9 bins to 8 bins for 2nd revision
    var_ranges = var_ranges, ncores=N_CORES,
    count_edges_on_axis = c(TRUE, TRUE, TRUE, FALSE, FALSE),
    # count_edges_on_axis2 = c(TRUE, FALSE, TRUE, FALSE, FALSE),
    pred_vars = terms_important))
  # t0.4 <- Sys.time()
  # t0.4 - t0

  n_niche_boot_list5d[[i]] <- metrics_boot_iter$n_niches
  n_limits_boot_list5d[[i]] <- metrics_boot_iter$n_limits
  # n_limits2_boot_list5d[[i]] <- metrics_boot_iter$n_limits2
  sorensen_boot_list5d[[i]] <- metrics_boot_iter$sorensen

  # t0.5 <- Sys.time()
  # c(t0.1, t0.2, t0.3, t0.4, t0.5) - t0
  # t0.5 - t0

  setTxtProgressBar(pb, i)
}
message("Ended at ", Sys.time())
close(pb)
rm(lognets_boot)


# Generate permutated null iterations for comparison to bootstraps ---------------


# "data/permuted_metrics_5d_wellfit.Rds" contains iterations 51 through 200.
# "data/permuted_metrics_5d_wellfit_173.Rds" contains 23 iterations below 51, and then
# iterations 51 through 200.

# rerun_ind <- which(unlist(lapply(n_niche_permuted_list5d, is.null)))

n_iters2 <- 150

pb <- txtProgressBar(max=n_iters2, style=3)
message("Started at ", Sys.time())
# for(i in 1:n_iters2) {
for(i in 51:100) {

  # t2.0 <- Sys.time()
  # Permuted sample (it says "permuteClimate" but it actually permutes row numbers in
  # sample data, which can include soil variables.)
  neon_dob_permuted <- permuteClimate(neon_dob_wellfit, randomseed=1010100+i)

  lognets_permuted <- suppressMessages(fitLognets(neon_dob_permuted, terms = terms_important,
                                                  ncores=N_CORES, randomseed = 1010101))

  thresholds_permuted <- suppressMessages(findThresholds(neon_dob_permuted, lognets_permuted, terms = terms_important, ncores = N_CORES))

  metrics_permuted_iter <- suppressMessages(getNicheMetrics5D0(
    neon_dob_permuted, lognets_permuted, thresholds_permuted, r_present_northam, inflate = 0.15, # Can constrain by hull later
    n_bins_per_axis = 8, # Changed from 9 bins to 8 bins for 2nd revision
    var_ranges = var_ranges, ncores=N_CORES,
    count_edges_on_axis = c(TRUE, TRUE, TRUE, FALSE, FALSE),
    # count_edges_on_axis2 = c(TRUE, FALSE, TRUE, FALSE, FALSE),
    pred_vars = terms_important))

  n_niche_permuted_list5d[[i]] <- metrics_permuted_iter$n_niches
  n_limits_permuted_list5d[[i]] <- metrics_permuted_iter$n_limits
  # n_limits2_permuted_list5d[[i]] <- metrics_permuted_iter$n_limits2
  sorensen_permuted_list5d[[i]] <- metrics_permuted_iter$sorensen

  # t2.1 <- Sys.time()
  # t2.1 - t2.0

  setTxtProgressBar(pb, i)
}
message("Ended at ", Sys.time())
close(pb)
rm(lognets_permuted)


## Save results ----------------------

boot_metrics_5d <- list(
  do.call(cbind, n_niche_boot_list5d),
  do.call(cbind, n_limits_boot_list5d),
  do.call(cbind, sorensen_boot_list5d))
names(boot_metrics_5d) <- c("n_niche", "n_limits", "sorensen")

permuted_metrics_5d <- list(
  do.call(cbind, n_niche_permuted_list5d),
  do.call(cbind, n_limits_permuted_list5d),
  do.call(cbind, sorensen_permuted_list5d))
names(permuted_metrics_5d) <- c("n_niche", "n_limits", "sorensen")

# saveRDS(boot_metrics_5d, "data/boot_metrics_5d_wellfit_8bins.Rds")
# saveRDS(permuted_metrics_5d, "data/permuted_metrics_5d_wellfit_8bins.Rds")

# boot_metrics_5d <- readRDS("data/boot_metrics_5d_wellfit.Rds")
# permuted_metrics_5d <- readRDS("data/permuted_metrics_5d_wellfit.Rds")

# permuted_metrics_5d_51to150 <- permuted_metrics_5d
# permuted_metrics_5d_51to150$n_niche <- permuted_metrics_5d_51to150$n_niche[,51:150]
# permuted_metrics_5d_51to150$n_limits <- permuted_metrics_5d_51to150$n_limits[,51:150]
# permuted_metrics_5d_51to150$sorensen <- permuted_metrics_5d_51to150$sorensen[,51:150]
# saveRDS(permuted_metrics_5d_51to150, "data/permuted_metrics_5d_wellfit_51to150.Rds")


# Standardize and plot results --------------------------------------------

# Get hull based on full (non-permuted) dataset

terms_important_nonquadratic <- c("mat_celsius", "temp_seasonality", "map_mm", "soilInCaClpH", "organicCPercent")
obs_env <- as(sample_data(neon_dob_wellfit), "data.frame")
inflate <- 0.15

inhull_pairs <- list()
for(i in 1:4) {
  for(j in (i+1):5) {
    varname_x <- terms_important_nonquadratic[i]
    varname_y <- terms_important_nonquadratic[j]
    x <- obs_env[,varname_x]
    y <- obs_env[,varname_y]
    complete_ind <- which(complete.cases(cbind(x, y)))
    x <- x[complete_ind]
    y <- y[complete_ind]
    ch_ind <- chull(x, y)
    infl_x <- x[ch_ind] + inflate * (x[ch_ind] - mean(x))
    infl_y <- y[ch_ind] + inflate * (y[ch_ind] - mean(y))
    inhull_pairs <- c(inhull_pairs, list(inpolygon(
      metrics_boot_iter[,varname_x,drop=TRUE], metrics_boot_iter[,varname_y,drop=TRUE],
      infl_x, infl_y, boundary=TRUE
    )))
  }
}
inhull_full <- rowSums(do.call(cbind, inhull_pairs)) == 10

metrics_boot_iter[,terms_important] %>%
  mutate(inhull = inhull_full) %>%
  mutate(n_niches_mean_boot = rowMeans(boot_metrics_5d$n_niche, na.rm=TRUE),
         n_limits_mean_boot = rowMeans(boot_metrics_5d$n_limits, na.rm=TRUE),
         sorensen_mean_boot = rowMeans(boot_metrics_5d$sorensen, na.rm=TRUE),
         n_niches_sd_boot = apply(boot_metrics_5d$n_niche, 1, sd, na.rm=TRUE),
         n_limits_sd_boot = apply(boot_metrics_5d$n_limits, 1, sd, na.rm=TRUE),
         sorensen_sd_boot = apply(boot_metrics_5d$sorensen, 1, sd, na.rm=TRUE),
         n_niches_mean_permuted = rowMeans(permuted_metrics_5d$n_niche, na.rm=TRUE),
         n_limits_mean_permuted = rowMeans(permuted_metrics_5d$n_limits, na.rm=TRUE),
         sorensen_mean_permuted = rowMeans(permuted_metrics_5d$sorensen, na.rm=TRUE)) %>%
  mutate(n_niches_mean_std = n_niches_mean_boot / n_niches_mean_permuted,
         n_limits_mean_std = n_limits_mean_boot / n_limits_mean_permuted,
         sorensen_mean_std = sorensen_mean_boot / sorensen_mean_permuted,
         n_niches_sd_std = n_niches_sd_boot / n_niches_mean_permuted,
         n_limits_sd_std = n_limits_sd_boot / n_limits_mean_permuted,
         sorensen_sd_std = sorensen_sd_boot / sorensen_mean_permuted) %>%
  mutate(n_niches_cv_std = n_niches_sd_std / n_niches_mean_std,
         n_limits_cv_std = n_limits_sd_std / n_limits_mean_std,
         sorensen_cv_std = sorensen_sd_std / sorensen_mean_std) ->
  metrics_5d

# saveRDS(metrics_5d, "data/metrics_5d_wellfit_8bins.Rds")

# metrics_5d <- readRDS("data/metrics_5d_wellfit_8bins.Rds")



## Make maps --------------------------------------------------------------


metrics_5d_inhull <- metrics_5d %>%
  mutate(across(n_niches_mean_boot:sorensen_cv_std, .fns = ~ if_else(is.na(sorensen_mean_std), NA_real_, .x)))


# Mean
map_edge_mean_std <- projectMetricToMap(
  metrics_5d_inhull, metric = "n_limits_mean_std", r_climate = r_present_northam,
  terms_nonquadratic = terms_important_nonquadratic, borders_sp = north_america,
  scale_text = "Mean niche edge\ndensity, std.", viridis_opt = "D")
map_overlap_mean_std <- projectMetricToMap(
  metrics_5d_inhull, metric = "n_niches_mean_std", r_climate = r_present_northam,
  terms_nonquadratic = terms_important_nonquadratic, borders_sp = north_america,
  scale_text = "Mean niche\noverlap, std.", viridis_opt = "E")
map_sor_mean_std <- projectMetricToMap(
  metrics_5d_inhull, metric = "sorensen_mean_std", r_climate = r_present_northam,
  terms_nonquadratic = terms_important_nonquadratic, borders_sp = north_america,
  scale_text = "Mean Sørensen\nsensitivity, std.")

# CV
map_edge_cv_std <- projectMetricToMap(
  metrics_5d_inhull, metric = "n_limits_cv_std", r_climate = r_present_northam,
  terms_nonquadratic = terms_important_nonquadratic, borders_sp = north_america,
  scale_text = "CV of niche edge\ndensity, std.", viridis_opt = "D")
map_overlap_cv_std <- projectMetricToMap(
  metrics_5d_inhull, metric = "n_niches_cv_std", r_climate = r_present_northam,
  terms_nonquadratic = terms_important_nonquadratic, borders_sp = north_america,
  scale_text = "CV of niche\noverlap, std.", viridis_opt = "E",
  scale_breaks = seq(0, 1, by=0.1))
map_sor_cv_std <- projectMetricToMap(
  metrics_5d_inhull, metric = "sorensen_cv_std", r_climate = r_present_northam,
  terms_nonquadratic = terms_important_nonquadratic, borders_sp = north_america,
  scale_text = "CV of Sørensen\nsensitivity, std.")

map_sor_cv_std <- map_sor_cv_std +
  scale_fill_viridis("CV of Sørensen\nsensitivity, std.", option="B", na.value="transparent",
                     limits=c(0, 0.3))
map_sor_cv_std

map_metrics <- egg::ggarrange(
  map_edge_mean_std, map_edge_cv_std,
  map_overlap_mean_std, map_overlap_cv_std,
  map_sor_mean_std, map_sor_cv_std,
  ncol=2,
  labels=c("(a)", "(d)",
           "(b)", "(e)",
           "(c)", "(f)"),
  label.args = list(gp=gpar(fontface="bold"))
)
ggsave("plots/fig5/metrics_5d_mean_cv_8bins.png", map_metrics, device="png",
       width=5, height=8, units="in")
ggsave("plots/fig5/metrics_5d_mean_cv_8bins.pdf", map_metrics, device="pdf",
       width=5, height=8, units="in")


# Variance-to-mean ratio as an index of dispersion/nonuniformity ------------------------

vartomean <- function(x, na.rm=TRUE) {
  var(x, na.rm=na.rm) / mean(x, na.rm=na.rm)
}

rm_n_niche_boot <- which(is.na(boot_metrics_5d$n_limits & boot_metrics_5d$n_niche==0))
rm_n_niche_permuted <- which(is.na(permuted_metrics_5d$n_limits & permuted_metrics_5d$n_niche==0))
boot_metrics_5d$n_niche[rm_n_niche_boot] <- NA
permuted_metrics_5d$n_niche[rm_n_niche_permuted] <- NA

# Even though each bootstrap iteration was constrained by a convex hull, it would've been
# a convex hull based on a subset (bootstrap sample) of the full dataset. That's why there
# are some metrics that are nonmissing even though they fall outside the full-data hull.
table(!is.na(boot_metrics_5d$n_limits[,1]), inhull_full)
# To reduce model extrapolation, it's important to constrain the metrics by the full-data
# hull.

vartomean_edge_boot <- apply(boot_metrics_5d$n_limits[inhull_full,], 2, vartomean) # Bootstrapped empirical
vartomean_edge_null <- apply(permuted_metrics_5d$n_limits[inhull_full,], 2, vartomean) # Permuted null

vartomean_overlap_boot <- apply(boot_metrics_5d$n_niche[inhull_full,], 2, vartomean) # Bootstrapped empirical
vartomean_overlap_null <- apply(permuted_metrics_5d$n_niche[inhull_full,], 2, vartomean) # Permuted null

vartomean_sor_boot <- apply(boot_metrics_5d$sorensen[inhull_full,], 2, vartomean) # Bootstrapped empirical
vartomean_sor_null <- apply(permuted_metrics_5d$sorensen[inhull_full,], 2, vartomean) # Permuted null

par(mfrow=c(2,1))
hist(vartomean_edge_boot)
hist(vartomean_edge_null)
hist(vartomean_overlap_boot)
hist(vartomean_overlap_null)
hist(vartomean_sor_boot)
hist(vartomean_sor_null)
par(mfrow=c(1,1))
par(mar=c(5.1, 4.1, 4.1, 2.1))

# Permutation-based P-values:
# Each value represents the P-value for one bootstrap iteration
pval_edge <- sapply(vartomean_edge_boot, function(x) mean(x < vartomean_edge_null))
pval_overlap <- sapply(vartomean_overlap_boot, function(x) mean(x < vartomean_overlap_null))
pval_sor <- sapply(vartomean_sor_boot, function(x) mean(x < vartomean_sor_null))

ks.test(pval_edge, "punif", alternative="greater") # P < 0.001
ks.test(pval_overlap, "punif", alternative="greater") # P < 0.001
ks.test(pval_sor, "punif", alternative="greater") # P < 0.001

pvalHist <- function(pval) {
  data.frame(pval) %>%
    ggplot(aes(x=pval)) +
    geom_histogram(color=1, fill="grey80") +
    scale_x_continuous("P-value", breaks=seq(0, 1, by=0.20), limits=c(-0.05,1.05)) +
    ylab("Count") +
    theme_custom
}

p_pval_edge <- pvalHist(pval_edge)
p_pval_overlap <- pvalHist(pval_overlap)
p_pval_sor <- pvalHist(pval_sor)
p_pvals <- egg::ggarrange(p_pval_edge, p_pval_overlap, p_pval_sor,
                          labels=c("(a)", "(b)", "(c)"),
                          label.args=list(gp=gpar(fontface="bold")))

ggsave("plots/fig4/pvals_8bins.png", p_pvals, device="png", width=2.6, height=5.5, units="in")

p_2d_mean_cv_pvals <- egg::ggarrange(
  p_edge_std, p_edge_cv_std, p_pval_edge,
  p_overlap_std, p_overlap_cv_std, p_pval_overlap,
  p_sor_std, p_sor_cv_std, p_pval_sor,
  ncol=3,
  labels=c("(a)", "(d)", "(g)",
           "(b)", "(e)", "(h)",
           "(c)", "(f)", "(i)"),
  label.args = list(gp=gpar(fontface="bold"))
)

# ggsave("plots/fig4/metrics_2d_mean_cv_with_pvals.png", p_2d_mean_cv_pvals, device="png", width=8, height=5.5, units="in")
# ggsave("plots/fig4/metrics_2d_mean_cv_with_pvals.pdf", p_2d_mean_cv_pvals, device="pdf", width=8, height=5.5, units="in")


# Figure S2.4: All niche metric maps --------------------------------------


p_blank <- ggplot() + theme_bw() + theme(panel.border = element_blank())

theme_label <- list(theme_void(),
                    theme(plot.title = element_text(size=18, hjust=0.5, vjust=-10),
                          axis.title = element_text(size=18, hjust=0.5),
                          axis.title.y = element_text(angle=-90, hjust=0.5, vjust=0.5)))
p_label_obs <- ggplot() + ggtitle(~underline("Observed")) + theme_label
p_label_null <- ggplot() + ggtitle(~underline("Null")) + theme_label
p_label_std <- ggplot() + ggtitle(~underline("Standardized")) + theme_label
p_label_edge <- ggplot() + ylab(~underline("Edge density")) + theme_label
p_label_overlap <- ggplot() + ylab(~underline("Niche overlap")) + theme_label
p_label_sor <- ggplot() + ylab(~underline("Climate\nsensitivity")) + theme_label

p_combined_maps <- egg::ggarrange(
  p_label_obs, p_label_null, p_label_std, p_blank,
  projectMetricToMap(metrics_5d_inhull, metric = "n_limits_mean_boot", r_climate = r_present_northam,
                     terms_nonquadratic = terms_important_nonquadratic, borders_sp = north_america_cropped,
                     scale_text = "Mean niche edge\ndensity", viridis_opt = "D") +
    theme(legend.text = element_text(angle=45, hjust=1)),
  projectMetricToMap(metrics_5d_inhull, metric = "n_limits_mean_permuted", r_climate = r_present_northam,
                     terms_nonquadratic = terms_important_nonquadratic, borders_sp = north_america_cropped,
                     scale_text = "Mean niche edge\ndensity, null", viridis_opt = "D") +
    theme(legend.text = element_text(angle=45, hjust=1)),
  projectMetricToMap(metrics_5d_inhull, metric = "n_limits_mean_std", r_climate = r_present_northam,
                     terms_nonquadratic = terms_important_nonquadratic, borders_sp = north_america_cropped,
                     scale_text = "Mean niche edge\ndensity, std.", viridis_opt = "D") +
    theme(legend.text = element_text(angle=45, hjust=1)),
  p_label_edge,
  projectMetricToMap(metrics_5d_inhull, metric = "n_niches_mean_boot", r_climate = r_present_northam,
                     terms_nonquadratic = terms_important_nonquadratic, borders_sp = north_america_cropped,
                     scale_text = "Mean niche\noverlap", viridis_opt = "E") +
    theme(legend.text = element_text(angle=45, hjust=1)),
  projectMetricToMap(metrics_5d_inhull, metric = "n_niches_mean_permuted", r_climate = r_present_northam,
                     terms_nonquadratic = terms_important_nonquadratic, borders_sp = north_america_cropped,
                     scale_text = "Mean niche\noverlap, null", viridis_opt = "E") +
    theme(legend.text = element_text(angle=45, hjust=1)),
  projectMetricToMap(metrics_5d_inhull, metric = "n_niches_mean_std", r_climate = r_present_northam,
                     terms_nonquadratic = terms_important_nonquadratic, borders_sp = north_america_cropped,
                     scale_text = "Mean niche\noverlap, std.", viridis_opt = "E") +
    theme(legend.text = element_text(angle=45, hjust=1)),
  p_label_overlap,
  projectMetricToMap(metrics_5d_inhull, metric = "sorensen_mean_boot", r_climate = r_present_northam,
                     terms_nonquadratic = terms_important_nonquadratic, borders_sp = north_america_cropped,
                     scale_text = "Mean Sørensen\nsensitivity", viridis_opt = "B") +
    theme(legend.text = element_text(angle=45, hjust=1)),
  projectMetricToMap(metrics_5d_inhull, metric = "sorensen_mean_permuted", r_climate = r_present_northam,
                     terms_nonquadratic = terms_important_nonquadratic, borders_sp = north_america_cropped,
                     scale_text = "Mean Sørensen\nsensitivity, null", viridis_opt = "B") +
    theme(legend.text = element_text(angle=45, hjust=1)),
  projectMetricToMap(metrics_5d_inhull, metric = "sorensen_mean_std", r_climate = r_present_northam,
                     terms_nonquadratic = terms_important_nonquadratic, borders_sp = north_america_cropped,
                     scale_text = "Mean Sørensen\nsensitivity, std.", viridis_opt = "B") +
    theme(legend.text = element_text(angle=45, hjust=1)),
  p_label_sor,
  ncol=4,
  heights=c(1,2,2,2),
  widths=c(2,2,2,1),
  labels=c("", "", "", "",
           "(a)", "(d)", "(g)", "",
           "(b)", "(e)", "(h)", "",
           "(c)", "(f)", "(i)", ""),
  label.args = list(gp=gpar(fontface="bold"))
)
ggsave("plots/fig-s2.4/maps_present_combined.png", p_combined_maps, device="png",
       width=9, height=10, units="in")

