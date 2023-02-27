# R code for Qin et al. on fungal niches and climate sensitivity
# 6-fig4.R: Make Figure 4 (niche metrics in 2D)

# Code originally from "permutations.R" and "find_uniform_null.R"

neon_dob_wellfit <- readRDS("data/neon_dob_wellfit_v1.0.Rds")

terms_important <- c("mat_celsius", "temp_seasonality", "map_mm",
                     "mat_celsius_2", "map_mm_2", "soilInCaClpH", "organicCPercent")


# Generate bootstraps on empirical data (to quantify uncertainty) --------------------

n_iters <- 50

n_niche_boot_list2d <- n_limits_boot_list2d <- sorensen_boot_list2d <- list()

pb <- txtProgressBar(max=n_iters, style=3)
message("Started at ", Sys.time())
for(i in 1:n_iters) { # Takes 1.75 minutes to run 1 iteration of 20x20 bins using 16 cores

  # Bootstrap sample
  neon_dob_boot <- bootstrapSamples(neon_dob_wellfit, randomseed=1010100+i)

  lognets_boot <- suppressMessages(fitLognets(neon_dob_boot, terms = terms_important, ncores=N_CORES, randomseed = 1010101))

  thresholds_boot <- suppressMessages(findThresholds(neon_dob_boot, lognets_boot, terms = terms_important, ncores = N_CORES))

  metrics_2d_boot <- suppressMessages(getNicheMetrics2D(
    neon_dob_boot, lognets_boot, thresholds_boot,
    "mat_celsius", "map_mm",
    # "mat_celsius", "temp_seasonality",
    x1_range = range(get_variable(neon_dob_wellfit, "mat_celsius")),
    x2_range = range(get_variable(neon_dob_wellfit, "map_mm")),
    # x2_range = range(get_variable(neon_dob_wellfit, "temp_seasonality")),
    constant_values = apply(as(sample_data(neon_dob_wellfit)[,terms_important], "data.frame"), 2, median, na.rm=TRUE),
    r_climate = r_present_northam, n_bins_per_axis = 20, ncores=N_CORES, pred_vars = terms_important,
    progressbar=FALSE, inflate=0.15))

  n_niche_boot_list2d[[i]] <- metrics_2d_boot$n_niches
  n_limits_boot_list2d[[i]] <- metrics_2d_boot$n_limits
  sorensen_boot_list2d[[i]] <- metrics_2d_boot$sorensen

  setTxtProgressBar(pb, i)
}
message("Ended at ", Sys.time())
close(pb)
rm(lognets_boot)



# Generate permutated null iterations for comparison to bootstraps ---------------

n_iters2 <- 150

n_niche_permuted_list2d <- n_limits_permuted_list2d <- sorensen_permuted_list2d <- list()

pb <- txtProgressBar(max=n_iters2, style=3)
message("Started at ", Sys.time())
for(i in 51:(n_iters2+50)) {

  # Permuted sample (it says "permuteClimate" but it actually permutes row numbers in
  # sample data, which can include soil variables.)
  neon_dob_permuted <- permuteClimate(neon_dob_wellfit, randomseed=1010100+i)

  lognets_permuted <- suppressMessages(fitLognets(neon_dob_permuted, terms = terms_important,
                                                  ncores=N_CORES, randomseed = 1010101))

  thresholds_permuted <- suppressMessages(findThresholds(neon_dob_permuted, lognets_permuted, terms = terms_important, ncores = N_CORES))

  metrics_2d_permuted <- suppressMessages(getNicheMetrics2D(
    neon_dob_permuted, lognets_permuted, thresholds_permuted,
    "mat_celsius", "map_mm",
    # "mat_celsius", "temp_seasonality"
    x1_range = range(get_variable(neon_dob_wellfit, "mat_celsius")),
    x2_range = range(get_variable(neon_dob_wellfit, "map_mm")),
    # x2_range = range(get_variable(neon_dob_wellfit, "temp_seasonality")),
    constant_values = apply(as(sample_data(neon_dob_wellfit)[,terms_important], "data.frame"), 2, median, na.rm=TRUE),
    r_climate = r_present_northam, n_bins_per_axis = 20, ncores=8, pred_vars = terms_important,
    progressbar=FALSE, inflate=0.15))

  n_niche_permuted_list2d[[i]] <- metrics_2d_permuted$n_niches
  n_limits_permuted_list2d[[i]] <- metrics_2d_permuted$n_limits
  sorensen_permuted_list2d[[i]] <- metrics_2d_permuted$sorensen

  setTxtProgressBar(pb, i-50)
}
message("Ended at ", Sys.time())
close(pb)
rm(lognets_permuted)

## Save results ----------------------

boot_metrics_2d <- list(
  do.call(cbind, n_niche_boot_list2d),
  do.call(cbind, n_limits_boot_list2d),
  do.call(cbind, sorensen_boot_list2d))
names(boot_metrics_2d) <- c("n_niche", "n_limits", "sorensen")

permuted_metrics_2d <- list(
  do.call(cbind, n_niche_permuted_list2d),
  do.call(cbind, n_limits_permuted_list2d),
  do.call(cbind, sorensen_permuted_list2d))
names(permuted_metrics_2d) <- c("n_niche", "n_limits", "sorensen")

# saveRDS(boot_metrics_2d, "data/boot_metrics_2d_wellfit.Rds")
# saveRDS(permuted_metrics_2d, "data/permuted_metrics_2d_wellfit.Rds")

boot_metrics_2d <- readRDS("data/boot_metrics_2d_wellfit.Rds")
permuted_metrics_2d <- readRDS("data/permuted_metrics_2d_wellfit.Rds")



# Standardize and plot results --------------------------------------------


metrics_2d_permuted[,terms_important] %>%
  mutate(n_niches_mean_boot = rowMeans(do.call(cbind, n_niche_boot_list2d), na.rm=TRUE),
         n_limits_mean_boot = rowMeans(do.call(cbind, n_limits_boot_list2d), na.rm=TRUE),
         sorensen_mean_boot = rowMeans(do.call(cbind, sorensen_boot_list2d), na.rm=TRUE),
         n_niches_sd_boot = apply(do.call(cbind, n_niche_boot_list2d), 1, sd, na.rm=TRUE),
         n_limits_sd_boot = apply(do.call(cbind, n_limits_boot_list2d), 1, sd, na.rm=TRUE),
         sorensen_sd_boot = apply(do.call(cbind, sorensen_boot_list2d), 1, sd, na.rm=TRUE),
         n_niches_mean_permuted = rowMeans(do.call(cbind, n_niche_permuted_list2d), na.rm=TRUE),
         n_limits_mean_permuted = rowMeans(do.call(cbind, n_limits_permuted_list2d), na.rm=TRUE),
         sorensen_mean_permuted = rowMeans(do.call(cbind, sorensen_permuted_list2d), na.rm=TRUE)) %>%
  mutate(n_niches_mean_std = n_niches_mean_boot / n_niches_mean_permuted,
         n_limits_mean_std = n_limits_mean_boot / n_limits_mean_permuted,
         sorensen_mean_std = sorensen_mean_boot / sorensen_mean_permuted,
         n_niches_sd_std = n_niches_sd_boot / n_niches_mean_permuted,
         n_limits_sd_std = n_limits_sd_boot / n_limits_mean_permuted,
         sorensen_sd_std = sorensen_sd_boot / sorensen_mean_permuted) %>%
  mutate(n_niches_cv_std = n_niches_sd_std / n_niches_mean_std,
         n_limits_cv_std = n_limits_sd_std / n_limits_mean_std,
         sorensen_cv_std = sorensen_sd_std / sorensen_mean_std) ->
  metrics_2d

# Define full-data hull

terms_important_nonquadratic <- c("mat_celsius", "temp_seasonality", "map_mm", "soilInCaClpH", "organicCPercent")
obs_env <- as(sample_data(neon_dob_wellfit), "data.frame")
inflate <- 0.15

x1 <- "mat_celsius"
x2 <- "map_mm"

x <- obs_env[,x1]
y <- obs_env[,x2]
ch_ind <- chull(x, y)
infl_x <- x[ch_ind] + inflate * (x[ch_ind] - mean(x))
infl_y <- y[ch_ind] + inflate * (y[ch_ind] - mean(y))
inhull <- pracma::inpolygon(metrics_2d[,x1,drop=TRUE],
                            metrics_2d[,x2,drop=TRUE], infl_x, infl_y, boundary=TRUE)

metrics_2d$inhull <- inhull

# saveRDS(metrics_2d, "data/metrics_2d_wellfit.Rds")
metrics_2d <- readRDS("data/metrics_2d_wellfit.Rds")

common_gg_elements <- list(
  guides(color="none"),
  xlab("Mean annual temp. (ºC)"),
  ylab("Annual precip. (mm)"),
  theme_bw(),
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
  )
)

obs_env_withproject <- sample_data(neon_dob_wellfit) %>%
  as("data.frame") %>%
  dplyr::select(all_of(c(terms_important_nonquadratic, "Project"))) %>%
  mutate(Project = factor(Project, levels=c("DoB", "Both", "NEON")))

plot2DMetric <- function(metrics_df, metric, viridis_opt="D", sites_df=obs_env_withproject,
                         x1="mat_celsius", x2="map_mm") {
  dplyr::filter(metrics_df, inhull==TRUE) %>%
    ggplot(aes_string(x=x1, y=x2)) +
    geom_tile(aes_string(fill=metric)) +
    scale_fill_viridis("", option=viridis_opt)  +
    # stat_contour(aes_string(z=metric), col = "white", size=0.1) +
    geom_point(data=sites_df, aes(col=Project), alpha=0.99, size=1, pch=21, fill="white") +
    common_gg_elements
}

p_edge_obs <- plot2DMetric(metrics_2d, "n_limits_mean_boot", "D")
p_edge_null <- plot2DMetric(metrics_2d, "n_limits_mean_permuted", "D")
p_edge_std <- plot2DMetric(metrics_2d, "n_limits_mean_std", "D")
p_overlap_obs <- plot2DMetric(metrics_2d, "n_niches_mean_boot", "E")
p_overlap_null <- plot2DMetric(metrics_2d, "n_niches_mean_permuted", "E")
p_overlap_std <- plot2DMetric(metrics_2d, "n_niches_mean_std", "E")
p_sor_obs <- plot2DMetric(metrics_2d, "sorensen_mean_boot", "B")
p_sor_null <- plot2DMetric(metrics_2d, "sorensen_mean_permuted", "B")
p_sor_std <- plot2DMetric(metrics_2d, "sorensen_mean_std", "B")

p_edge_cv_std <- plot2DMetric(metrics_2d, "n_limits_cv_std", "D")
p_overlap_cv_std <- plot2DMetric(metrics_2d, "n_niches_cv_std", "E")
p_sor_cv_std <- plot2DMetric(metrics_2d, "sorensen_cv_std", "B")


p_2d_mean_cv <- egg::ggarrange(
  p_edge_std, p_edge_cv_std,
  p_overlap_std, p_overlap_cv_std,
  p_sor_std, p_sor_cv_std,
  ncol=2,
  labels=c("a", "d",
           "b", "e",
           "c", "f"),
  label.args = list(gp=gpar(fontface="bold"))
)
# plot(p_2d_mean_cv)
# ggsave("plots/fig4/metrics_2d_mean_cv_without_pvals.png", p_2d_mean_cv, device="png",
#        width=6.5, height=6, units="in")
# ggsave("plots/fig4/metrics_2d_mean_cv_without_pvals.pdf", p_2d_mean_cv, device="pdf",
#        width=6.5, height=6, units="in")



## Figure S2.3 ---------------------------------------------------------------


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

p_2d_combined <- egg::ggarrange(
  p_label_obs, p_label_null, p_label_std, p_blank,
  p_edge_obs, p_edge_null, p_edge_std, p_label_edge,
  p_overlap_obs, p_overlap_null, p_overlap_std, p_label_overlap,
  p_sor_obs, p_sor_null, p_sor_std, p_label_sor,
  ncol=4,
  heights=c(1,2,2,2),
  widths=c(2,2,2,1),
  labels=c("", "", "", "",
           "(a)", "(d)", "(g)", "",
           "(b)", "(e)", "(h)", "",
           "(c)", "(f)", "(i)", ""),
  label.args = list(gp=gpar(fontface="bold"))
)
plot(p_2d_combined)
ggsave("plots/fig-s2.3/calculation_of_std_sorensen.png", p_2d_combined, device="png",
       width=12, height=7, units="in")


# Sanity check: What does a single iteration look like? -------------------

# lognets_full <- fitLognets(neon_dob_wellfit, terms = terms_important, ncores=N_CORES, randomseed = 1010101)
#
# thresholds_full <- findThresholds(neon_dob_wellfit, lognets_full, terms = terms_important, ncores = N_CORES)
#
# metrics_2d_full <- getNicheMetrics2D(
#   neon_dob_wellfit, lognets_full, thresholds_full, "mat_celsius", "map_mm",
#   r_climate = r_present_northam, n_bins_per_axis = 20, ncores=N_CORES,
#   pred_vars = terms_important, progressbar=FALSE)
#
# plotNicheMetrics2D(neon_dob_wellfit, metrics_2d_full, label1="n_niches full",
#                    label2="n_limits full", label3="sorensen full", arrange="vertical")



# With original, climate-only predictors

# terms0 <- c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality", "mat_celsius_2", "map_mm_2")
#
# lognets_full0 <- fitLognets(neon_dob_wellfit, terms = terms0, ncores=N_CORES, randomseed = 1010101)
#
# thresholds_full0 <- findThresholds(neon_dob_wellfit, lognets_full0, terms = terms0, ncores = N_CORES)
#
# metrics_2d_full0 <- getNicheMetrics2D(
#   neon_dob_wellfit, lognets_full0, thresholds_full0, "mat_celsius", "map_mm",
#   r_climate = r_present_northam, n_bins_per_axis = 20, ncores=N_CORES,
#   pred_vars = terms0, progressbar=FALSE)
#
# plotNicheMetrics2D(neon_dob_wellfit, metrics_2d_full0, label1="n_niches clim. only",
#                    label2="n_limits clim. only", label3="sorensen clim. only", arrange="vertical")


# Variance-to-mean ratio as an index of dispersion/nonuniformity ------------------------

vartomean <- function(x, na.rm=TRUE) {
  var(x, na.rm=na.rm) / mean(x, na.rm=na.rm)
}
print_vartomean <- function(x) {
  print(c(var(x, na.rm=TRUE), mean(x, na.rm=TRUE), var(x, na.rm=TRUE) / mean(x, na.rm=TRUE)))
}

vartomean_edge_boot <- apply(boot_metrics_2d[["n_limits"]][inhull,], 2, vartomean) # Bootstrapped empirical
vartomean_edge_null <- apply(permuted_metrics_2d[["n_limits"]][inhull,], 2, vartomean) # Permuted null
vartomean_overlap_boot <- apply(boot_metrics_2d[["n_niche"]][inhull,], 2, vartomean) # Bootstrapped empirical
vartomean_overlap_null <- apply(permuted_metrics_2d[["n_niche"]][inhull,], 2, vartomean) # Permuted null
vartomean_sor_boot <- apply(boot_metrics_2d[["sorensen"]][inhull,], 2, vartomean) # Bootstrapped empirical
vartomean_sor_null <- apply(permuted_metrics_2d[["sorensen"]][inhull,], 2, vartomean) # Permuted null

par(mfrow=c(2,1))
hist(boot_metrics_2d[["n_limits"]][,1])
hist(permuted_metrics_2d[["n_limits"]][,1])
print_vartomean(boot_metrics_2d[["n_limits"]][,1])
print_vartomean(permuted_metrics_2d[["n_limits"]][,1])
hist(vartomean_edge_boot)
hist(vartomean_edge_null)
par(mfrow=c(1,1))
par(mar=c(5.1, 4.1, 4.1, 2.1))



# Permutation-based P-values
pval_edge_2d <- sapply(vartomean_edge_boot, function(x) mean(x < vartomean_edge_null))
pval_overlap_2d <- sapply(vartomean_overlap_boot, function(x) mean(x < vartomean_overlap_null))
pval_sor_2d <- sapply(vartomean_sor_boot, function(x) mean(x < vartomean_sor_null))

ks.test(pval_edge_2d, "punif", alternative="greater")
ks.test(pval_overlap_2d, "punif", alternative="greater")
ks.test(pval_sor_2d, "punif", alternative="greater")

pvalHist <- function(pval) {
  data.frame(pval) %>%
    ggplot(aes(x=pval)) +
    geom_histogram(color=1, fill="grey80") +
    scale_x_continuous("P-value", breaks=seq(0, 1, by=0.20), limits=c(-0.05,1.05)) +
    ylab("Count") +
    theme_custom
}

pvalHist(pval_edge_2d)
pvalHist(pval_overlap_2d)
pvalHist(pval_sor_2d)

# Each value represents the P-value for one bootstrap iteration

# data.frame(pval) %>%
#   ggplot(aes(x=pval)) +
#   # geom_histogram(aes(y=..density..), color=1, fill="grey80") +
#   geom_density(lwd = 1, colour = 4, fill = 4, alpha = 0.25) +
#   geom_vline(xintercept=0.05, linetype="dashed", col="red") +
#   scale_x_continuous("P-value", breaks=seq(0, 1, by=0.05)) +
#   # scale_y_continuous("Density", breaks=NULL) +
#   theme_custom
# # theme(axis.text.y = element_blank())
# ggsave("plots/fig4/pval.png", device="png", width=6.5, height=2, units="in")


# ## Plot some individual iterations to see what their distributions look like -----
#
# # Niche edge density
# p_edge_boot1 <- plot2DMetric(cbind(metrics_2d, metric=boot_metrics_2d[["n_limits"]][,1]), "metric") + ggtitle("Boot1")
# p_edge_boot2 <- plot2DMetric(cbind(metrics_2d, metric=boot_metrics_2d[["n_limits"]][,2]), "metric") + ggtitle("Boot2")
# p_edge_boot3 <- plot2DMetric(cbind(metrics_2d, metric=boot_metrics_2d[["n_limits"]][,3]), "metric") + ggtitle("Boot3")
#
# p_edge_null1 <- plot2DMetric(cbind(metrics_2d, metric=permuted_metrics_2d[["n_limits"]][,1]), "metric") + ggtitle("Null1")
# p_edge_null2 <- plot2DMetric(cbind(metrics_2d, metric=permuted_metrics_2d[["n_limits"]][,2]), "metric") + ggtitle("Null2")
# p_edge_null3 <- plot2DMetric(cbind(metrics_2d, metric=permuted_metrics_2d[["n_limits"]][,3]), "metric") + ggtitle("Null3")
#
# p_indiv_iters_2d_edge <- egg::ggarrange(
#   p_edge_boot1, p_edge_boot2, p_edge_boot3,
#   p_edge_null1, p_edge_null2, p_edge_null3, ncol=3)
#
#
# p_overlap_boot1 <- plot2DMetric(cbind(metrics_2d, metric=boot_metrics_2d[["n_niche"]][,1]), "metric") + ggtitle("Boot1")
# p_overlap_boot2 <- plot2DMetric(cbind(metrics_2d, metric=boot_metrics_2d[["n_niche"]][,2]), "metric") + ggtitle("Boot2")
# p_overlap_boot3 <- plot2DMetric(cbind(metrics_2d, metric=boot_metrics_2d[["n_niche"]][,3]), "metric") + ggtitle("Boot3")
#
# p_overlap_null1 <- plot2DMetric(cbind(metrics_2d, metric=permuted_metrics_2d[["n_niche"]][,1]), "metric") + ggtitle("Null1")
# p_overlap_null2 <- plot2DMetric(cbind(metrics_2d, metric=permuted_metrics_2d[["n_niche"]][,2]), "metric") + ggtitle("Null2")
# p_overlap_null3 <- plot2DMetric(cbind(metrics_2d, metric=permuted_metrics_2d[["n_niche"]][,3]), "metric") + ggtitle("Null3")
#
# vartomean(boot_metrics_2d[["n_niche"]][,2])
# vartomean(permuted_metrics_2d[["n_niche"]][,2])
#
# hist(boot_metrics_2d[["n_niche"]][,2])
# hist(permuted_metrics_2d[["n_niche"]][,2])
#
# p_indiv_iters_2d_overlap <- egg::ggarrange(
#   p_overlap_boot1, p_overlap_boot2, p_overlap_boot3,
#   p_overlap_null1, p_overlap_null2, p_overlap_null3, ncol=3)
#
