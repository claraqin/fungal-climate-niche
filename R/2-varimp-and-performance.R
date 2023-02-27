# R code for Qin et al. on fungal niches and climate sensitivity
# 2-varimp-and-performance.R: Justify variables used in models based on performance (for Appendix)

# Code originally from "soil.R"


# 1. Compare performance of models with soil data vs. with climate data only ----------------


neon_dob_prevalent <- readRDS("data/neon_dob_prevalent_v4.1.Rds")

# Split to train and test sets
set.seed(10101)
train_set <- sample(c(TRUE, FALSE), size=length(sample_sums(neon_dob_prevalent)),
                    prob=c(0.70, 0.30), replace=TRUE)
test_set <- !train_set
neon_dob_prevalent_train <- prune_samples(train_set, neon_dob_prevalent)
neon_dob_prevalent_test <- prune_samples(test_set, neon_dob_prevalent)
n_taxa <- 8597
spp_subset_ind <- sort(sample(1:length(taxa_names(neon_dob_prevalent_train)), size=n_taxa))
neon_dob_prevalent_train_sppsub <- prune_taxa(taxa_names(neon_dob_prevalent_train)[spp_subset_ind],
                                              neon_dob_prevalent_train)
neon_dob_prevalent_test_sppsub <- prune_taxa(taxa_names(neon_dob_prevalent_test)[spp_subset_ind],
                                             neon_dob_prevalent_test)

# Define sets of variables

# climate_only
terms1 <- c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality",
            "mat_celsius_2", "map_mm_2")
# climate_and_soil
terms2 <- c("mat_celsius", "temp_seasonality", "map_mm", "prec_seasonality",
            "mat_celsius_2", "map_mm_2", "soilInCaClpH", "nitrogenPercent", "organicCPercent")
# important_vars (a posteriori naming)
terms3 <- c("mat_celsius", "temp_seasonality", "map_mm",
            "mat_celsius_2", "map_mm_2", "soilInCaClpH", "organicCPercent")
# climate_and_ph
terms4 <- c("mat_celsius", "temp_seasonality", "map_mm",
            "mat_celsius_2", "map_mm_2", "soilInCaClpH")

predictors_train1 <- as(sample_data(neon_dob_prevalent_train_sppsub)[,terms1], "data.frame")
predictors_train_std1 <- scale(predictors_train1, center=TRUE, scale=TRUE)
complete_records1 <- which(complete.cases(predictors_train_std1))

predictors_train2 <- as(sample_data(neon_dob_prevalent_train_sppsub)[,terms2], "data.frame")
predictors_train_std2 <- scale(predictors_train2, center=TRUE, scale=TRUE)
complete_records2 <- which(complete.cases(predictors_train_std2))

predictors_train3 <- as(sample_data(neon_dob_prevalent_train_sppsub)[,terms3], "data.frame")
predictors_train_std3 <- scale(predictors_train3, center=TRUE, scale=TRUE)
complete_records3 <- which(complete.cases(predictors_train_std3))

predictors_train4 <- as(sample_data(neon_dob_prevalent_train_sppsub)[,terms4], "data.frame")
predictors_train_std4 <- scale(predictors_train4, center=TRUE, scale=TRUE)
complete_records4 <- which(complete.cases(predictors_train_std4))

registerDoParallel(N_CORES)
lognets_train1 <- foreach(i = seq_along(spp_subset_ind)) %dopar% {
  y01 <- as.vector(otu_table(neon_dob_prevalent_train_sppsub)[,i]) > 0
  set.seed(1010101)
  tryCatch(cv.glmnet(as.matrix(predictors_train_std1)[complete_records1,],
                     y01[complete_records1], family = "binomial",
                     alpha = 0, type.measure = "default", standardize=FALSE),
           error = function(e) return(NA))
}

lognets_train2 <- foreach(i = seq_along(spp_subset_ind)) %dopar% {
  y01 <- as.vector(otu_table(neon_dob_prevalent_train_sppsub)[,i]) > 0
  set.seed(1010101)
  tryCatch(cv.glmnet(as.matrix(predictors_train_std2)[complete_records2,],
                     y01[complete_records2], family = "binomial",
                     alpha = 0, type.measure = "default", standardize=FALSE),
           error = function(e) return(NA))
}

lognets_train3 <- foreach(i = seq_along(spp_subset_ind)) %dopar% {
  y01 <- as.vector(otu_table(neon_dob_prevalent_train_sppsub)[,i]) > 0
  set.seed(1010101)
  tryCatch(cv.glmnet(as.matrix(predictors_train_std3)[complete_records3,],
                     y01[complete_records3], family = "binomial",
                     alpha = 0, type.measure = "default", standardize=FALSE),
           error = function(e) return(NA))
}

lognets_train4 <- foreach(i = seq_along(spp_subset_ind)) %dopar% {
  y01 <- as.vector(otu_table(neon_dob_prevalent_train_sppsub)[,i]) > 0
  set.seed(1010101)
  tryCatch(cv.glmnet(as.matrix(predictors_train_std4)[complete_records4,],
                     y01[complete_records4], family = "binomial",
                     alpha = 0, type.measure = "default", standardize=FALSE),
           error = function(e) return(NA))
}
stopImplicitCluster()

lognet_thresholds_train1 <- findThresholds(
  neon_dob_prevalent_train_sppsub, lognets_train1,
  terms = terms1)

lognet_thresholds_train2 <- findThresholds(
  neon_dob_prevalent_train_sppsub, lognets_train2,
  terms = terms2)

lognet_thresholds_train3 <- findThresholds(
  neon_dob_prevalent_train_sppsub, lognets_train3,
  terms = terms3)

lognet_thresholds_train4 <- findThresholds(
  neon_dob_prevalent_train_sppsub, lognets_train4,
  terms = terms4)

# Predict using Lognets

newdata1 <- as(sample_data(neon_dob_prevalent_test), "data.frame")[,terms1]
newdata_std1 <- scaleToReference(newdata1, predictors_train1)

newdata2 <- as(sample_data(neon_dob_prevalent_test), "data.frame")[,terms2]
newdata_std2 <- scaleToReference(newdata2, predictors_train2)

newdata3 <- as(sample_data(neon_dob_prevalent_test), "data.frame")[,terms3]
newdata_std3 <- scaleToReference(newdata3, predictors_train3)

newdata4 <- as(sample_data(neon_dob_prevalent_test), "data.frame")[,terms4]
newdata_std4 <- scaleToReference(newdata4, predictors_train4)

registerDoParallel(N_CORES)
lognet_preds1 <- foreach(i = 1:n_taxa) %dopar% {
  if(identical(lognets_train1[[i]], NA)) return(NA)
  prediction <- predict(lognets_train1[[i]],
                        newx = as.matrix(newdata_std1),
                        type="response", s = "lambda.min")
  classification <- prediction > lognet_thresholds_train1[[i]]
  as.numeric(classification)
}

lognet_preds2 <- foreach(i = 1:n_taxa) %dopar% {
  if(identical(lognets_train2[[i]], NA)) return(NA)
  prediction <- predict(lognets_train2[[i]],
                        newx = as.matrix(newdata_std2),
                        type="response", s = "lambda.min")
  classification <- prediction > lognet_thresholds_train2[[i]]
  as.numeric(classification)
}

lognet_preds3 <- foreach(i = 1:n_taxa) %dopar% {
  if(identical(lognets_train3[[i]], NA)) return(NA)
  prediction <- predict(lognets_train3[[i]],
                        newx = as.matrix(newdata_std3),
                        type="response", s = "lambda.min")
  classification <- prediction > lognet_thresholds_train3[[i]]
  as.numeric(classification)
}

lognet_preds4 <- foreach(i = 1:n_taxa) %dopar% {
  if(identical(lognets_train4[[i]], NA)) return(NA)
  prediction <- predict(lognets_train4[[i]],
                        newx = as.matrix(newdata_std4),
                        type="response", s = "lambda.min")
  classification <- prediction > lognet_thresholds_train4[[i]]
  as.numeric(classification)
}
stopImplicitCluster()

# Validate predictions

validation_predictions <- list(
  lognet_preds1,
  lognet_preds2,
  lognet_preds3,
  lognet_preds4
)
names(validation_predictions) <- c("climate_only", "climate_and_soil", "important_vars", "climate_and_ph")

library(caret)
accuracy <- matrix(nrow=n_taxa, ncol=length(validation_predictions), dimnames=list(spp_subset_ind, names(validation_predictions)))
for(i in 1:n_taxa) {
  sp <- taxa_names(neon_dob_prevalent_train)[spp_subset_ind[i]]
  y <- as.numeric(otu_table(neon_dob_prevalent_test)[,sp] > 0)
  for(mod in colnames(accuracy)) {
    accuracy[i,mod] <- mean(y == validation_predictions[[mod]][[i]], na.rm=TRUE)
  }
}

sensitivity <- specificity <- matrix(nrow=n_taxa, ncol=length(validation_predictions), dimnames=list(spp_subset_ind, names(validation_predictions)))
for(i in 1:n_taxa) {
  sp <- taxa_names(neon_dob_prevalent_train)[spp_subset_ind[i]]
  y <- as.numeric(otu_table(neon_dob_prevalent_test)[,sp] > 0)
  if(all(c(0,1) %in% y)) {
    confusions <- lapply(
      validation_predictions,
      function(x) tryCatch(
        confusionMatrix(factor(x[[i]], levels=c(0,1)),
                        factor(y, levels=c(0,1)), positive="1"),
        error = function(e) return(NA))
    )
    sensitivity[i,] <- sapply(confusions, function(x) tryCatch(
      x$byClass["Sensitivity"], error = function(e) return(NA)))
    specificity[i,] <- sapply(confusions, function(x) tryCatch(
      x$byClass["Specificity"], error = function(e) return(NA)))
  } else {
    sensitivity[i,] <- NA
    specificity[i,] <- NA
  }
}

apply(accuracy, 2, mean, na.rm=TRUE)
apply(sensitivity, 2, mean, na.rm=TRUE)
apply(specificity, 2, mean, na.rm=TRUE)

apply(accuracy, 2, sd, na.rm=TRUE)
apply(sensitivity, 2, sd, na.rm=TRUE)
apply(specificity, 2, sd, na.rm=TRUE)

apply(sensitivity + specificity - 1, 2, mean, na.rm=TRUE) # True skill statistic
apply(sensitivity + specificity - 1, 2, sd, na.rm=TRUE) # True skill statistic

# > apply(accuracy, 2, mean, na.rm=TRUE)
#     climate_only climate_and_soil   important_vars   climate_and_ph
#        0.7237337        0.7367193        0.7378173        0.7376252
# > apply(sensitivity, 2, mean, na.rm=TRUE)
#     climate_only climate_and_soil   important_vars   climate_and_ph
#        0.6402330        0.6499021        0.6473350        0.6401724
# > apply(specificity, 2, mean, na.rm=TRUE)
#     climate_only climate_and_soil   important_vars   climate_and_ph
#        0.7382305        0.7496496        0.7508180        0.7521499
# > apply(accuracy, 2, sd, na.rm=TRUE)
#     climate_only climate_and_soil   important_vars   climate_and_ph
#        0.1044462        0.1021051        0.1005891        0.1003766
# > apply(sensitivity, 2, sd, na.rm=TRUE)
#     climate_only climate_and_soil   important_vars   climate_and_ph
#        0.3126381        0.3020326        0.3042868        0.3033765
# > apply(specificity, 2, sd, na.rm=TRUE)
#     climate_only climate_and_soil   important_vars   climate_and_ph
#        0.1239998        0.1181267        0.1183131        0.1173362
# > apply(sensitivity + specificity - 1, 2, mean, na.rm=TRUE) # True skill statistic
#     climate_only climate_and_soil   important_vars   climate_and_ph
#        0.3784635        0.3995919        0.3981530        0.3923222
# > apply(sensitivity + specificity - 1, 2, sd, na.rm=TRUE) # True skill statistic
#     climate_only climate_and_soil   important_vars   climate_and_ph
#        0.3026740        0.2993346        0.2945925        0.2958201


# 2. Summarize coefficients and variable importance ----------

coef_matrix_lognet <- matrix(NA, nrow = length(lognets_train2), ncol=10)
rownames(coef_matrix_lognet) <- names(lognets_train2)
colnames(coef_matrix_lognet) <- rownames(coef(lognets_train2[[1]], s="lambda.min"))
system.time(
  for(i in 1:nrow(coef_matrix_lognet)) { # takes 1 minute
    if(!identical(lognets_train2[[i]], NA)) {
      coef_matrix_lognet[i,] <- coef(lognets_train2[[i]], s="lambda.min")[,1]
    }
  }
)
head(coef_matrix_lognet)

colnames(coef_matrix_lognet) <- c("Intercept", "MAT", "TSEA", "MAP", "PSEA",
                                  "MAT2", "MAP2", "pH", "%N", "%C")

coef_long_lognet <- coef_matrix_lognet[,-1] %>%
  as.data.frame() %>%
  tidyr::pivot_longer(MAT:`%C`)

# Variable importance is given by taking the mean of
# the absolute values of the standardized coefficients
mean_abs_coef <- coef_long_lognet %>%
  group_by(name) %>%
  summarise(mean = mean(abs(value), na.rm=TRUE)) %>%
  arrange(desc(mean)) %>%
  mutate(name = factor(name, levels=unique(name)))
mean_abs_coef

#   name   mean
#   <fct> <dbl>
# 1 TSEA  0.408
# 2 MAT   0.392
# 3 MAP   0.334
# 4 %C    0.329
# 5 pH    0.296
# 6 MAT2  0.288
# 7 PSEA  0.254
# 8 %N    0.177
# 9 MAP2  0.175

median_abs_coef <- coef_long_lognet %>%
  group_by(name) %>%
  summarise(median = median(abs(value), na.rm=TRUE)) %>%
  arrange(desc(median)) %>%
  mutate(name = factor(name, levels=unique(name)))
median_abs_coef

#   name  median
#   <fct>  <dbl>
# 1 MAT   0.200
# 2 TSEA  0.196
# 3 pH    0.193
# 4 %C    0.176
# 5 MAT2  0.169
# 6 MAP   0.155
# 7 PSEA  0.155
# 8 %N    0.0969
# 9 MAP2  0.0899

varimp_summary <- merge(mean_abs_coef, median_abs_coef) %>%
  tidyr::pivot_longer(mean:median, names_to = "statistic")

coef_long_lognet %>%
  mutate(name = factor(name, levels=unique(median_abs_coef$name))) %>%
  ggplot(aes(y=abs(value), x=name)) +
  geom_violin(fill=NA, lwd=0.3) +
  geom_jitter(width=0.45, height=0, alpha=0.15, shape=46, col="red") +
  geom_point(data=varimp_summary, aes(shape=statistic)) +
  scale_shape_manual("", values=c(21, 24)) +
  # geom_errorbar(data = mean_abs_coef, aes(y=mean, ymin=mean, ymax=mean), width=0.6, alpha=0.5) +
  # geom_errorbar(data = median_abs_coef, aes(y=median, ymin=median, ymax=median), width=0.6, alpha=0.5, linetype="dotted") +
  ylab("Variable importance") + xlab("Variable") +
  theme_custom +
  # guides(shape="none") +
  theme(legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.position = c(0.99, 1),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        legend.box.margin = margin(0, 0, 0, 0))
ggsave("plots/fig-s1.1/coef_jitter_lognet_with_soil.png", device="png", width=5, height=4, units="in")


# 3. Plot TSS vs. #presence and min(#pres, #abs) -----------------

n_pres <- n_abs <- sens <- spec <- c()
for(i in 1:n_taxa) {
  sp <- taxa_names(neon_dob_prevalent_train)[spp_subset_ind[i]]
  y_train <- as.numeric(otu_table(neon_dob_prevalent_train)[,sp] > 0)
  n_pres[i] <- sum(y_train==1)
  n_abs[i] <- sum(y_train==0)

  y_test <- as.numeric(otu_table(neon_dob_prevalent_test)[,sp] > 0)
  if(all(c(0,1) %in% y_test)) {
    confusions <- lapply(
      validation_predictions,
      function(x) tryCatch(
        confusionMatrix(factor(x[[i]], levels=c(0,1)),
                        factor(y_test, levels=c(0,1)), positive="1"),
        error = function(e) return(NA))
    )
    sens[i] <- sapply(confusions, function(x) tryCatch(
      x$byClass["Sensitivity"],
      error = function(e) return(NA)))[3] # This is for the important_vars model
    spec[i] <- sapply(confusions, function(x) tryCatch(
      x$byClass["Specificity"],
      error = function(e) return(NA)))[3] # This is for the important_vars model
  } else {
    sens[i] <- NA
    spec[i] <- NA
  }
}
n_pres
n_abs
sens
spec

n_vs_tss <- data.frame(
  otu = taxa_names(neon_dob_prevalent_train)[spp_subset_ind],
  n_pres = n_pres,
  n_abs = n_abs,
  min_n = pmin(n_pres, n_abs),
  sens = sens,
  spec = spec,
  tss = sens + spec - 1
)

p_pres_vs_tss <- ggplot(n_vs_tss, aes(x = n_pres, y=tss)) +
  geom_point(alpha=0.1) +
  theme_bw() +
  xlab("No. presences in training set") + ylab("TSS")
p_minpa_vs_tss <- ggplot(n_vs_tss, aes(x = min_n, y=tss)) +
  geom_point(alpha=0.1) +
  theme_bw() +
  xlab("min(#Pres., #Abs.)") + ylab("TSS") +
  geom_smooth(method="lm")
p_minpa_vs_sens <- ggplot(n_vs_tss, aes(x = min_n, y=sens)) +
  geom_point(alpha=0.1) +
  theme_bw() +
  xlab("min(#Pres., #Abs.)") + ylab("Sensitivity") +
  geom_smooth(method="lm")
p_minpa_vs_spec <- ggplot(n_vs_tss, aes(x = min_n, y=spec)) +
  geom_point(alpha=0.1) +
  theme_bw() +
  xlab("min(#Pres., #Abs.)") + ylab("Specificity") +
  geom_smooth(method="lm")

p_minpa <- egg::ggarrange(
  p_minpa_vs_tss, p_minpa_vs_sens, p_minpa_vs_spec,
  ncol=3, labels=c("a", "b", "c"),
  label.args = list(gp=gpar(fontface="bold"))
)

ggsave("plots/fig-s1.2/min_pres_abs_vs_tss.png", p_minpa, device="png",
       width=8, height=2.5, units="in")

## Also plot TSS vs. guild --------

n_vs_tss2 <- n_vs_tss %>%
  left_join(cbind.data.frame(otu = taxa_names(neon_dob_prevalent_train), as(tax_table(neon_dob_prevalent_train), "matrix"))) %>%
  left_join(dplyr::select(ft, genus, guild="guild2")) %>%
  mutate(guild = as.character(guild),
         guild = if_else(guild=="EM", "EcM", guild),
         guild = if_else(guild=="saprotroph", "Sapr.", guild),
         guild = if_else(is.na(guild) | !(guild %in% c("EcM", "Sapr.")), "Other", guild),
         guild = factor(guild, levels = c("EcM", "Sapr.", "Other")))

p_tss_vs_guild <- ggplot(n_vs_tss2, aes(x=guild, y=tss)) +
  geom_boxplot(alpha=0.1) +
  theme_bw() +
  xlab("Guild") + ylab("TSS")
p_tss_vs_guild
with(n_vs_tss2, table(guild))

p_minpa_vs_guild <- ggplot(n_vs_tss2, aes(x=guild, y=min_n)) +
  geom_boxplot(alpha=0.1) +
  theme_bw() +
  xlab("Guild") + ylab("min(#Pres., #Abs.)")
p_minpa_vs_guild

p_minpa_and_guild <- egg::ggarrange(
  p_minpa_vs_tss, p_minpa_vs_sens, p_minpa_vs_spec,
  p_tss_vs_guild, p_minpa_vs_guild,
  ncol=3, labels=c("(a)", "(b)", "(c)", "(d)", "(e)"),
  label.args = list(gp=gpar(fontface="bold"))
)

ggsave("plots/fig-s1.2/tss_vs_min_pres_abs_and_guild.png", p_minpa_and_guild, device="png",
       width=8, height=5, units="in")


## Subset of taxa with TSS > 0 (neon_dob_wellfit)

spp_postss <- as.character(n_vs_tss2$otu[which(n_vs_tss2$tss>0)])

saveRDS(n_vs_tss2, "data/n_vs_tss.Rds")
saveRDS(spp_postss, "data/spp_postss.Rds")

neon_dob_wellfit <- prune_taxa(spp_postss, neon_dob_prevalent)
saveRDS(neon_dob_wellfit, "data/neon_dob_wellfit_v1.0.Rds")


# 4. Compare to maxnet -------------------------------------------------------


# Fit MaxNets

set.seed(1010101)
randomBgSites <- randomPoints(env_mask_0.15, 10000)
randomBg <- cbind.data.frame(lon=randomBgSites[,1], lat=randomBgSites[,2],
                             raster::extract(r_present_northam, randomBgSites))
df_present_northam <- as.data.frame(r_present_northam) %>%
  mutate(mat_celsius_2 = mat_celsius^2,
         map_mm_2 = map_mm^2)
temp <- raster(r_present_northam)
r_present_northam_expanded <- r_present_northam
for(i in 8:length(names(df_present_northam))) {
  temp[] <- df_present_northam[,i]
  r_present_northam_expanded <- addLayer(r_present_northam_expanded, temp)
  names(r_present_northam_expanded)[i] <- names(df_present_northam)[i]
}
names(r_present_northam_expanded)
sum(complete.cases(as.matrix(r_present_northam_expanded)))

registerDoParallel(N_CORES)
t_maxnet0 <- Sys.time()

oper <- foreach(i = spp_subset_ind,
# oper <- foreach(i = 1:500,
                .combine = nested_combine, .multicombine=TRUE,
                .init = list(list(), list())) %dopar% {
  sp <- taxa_names(neon_dob_prevalent_train)[i]

  # Get presence-absence vector for this species
  pa <- otu_table(neon_dob_prevalent_train)[,sp]
  presence_ind <- which(pa > 0)

  # ... and presence-absence vector to train model on
  trainPA <- c(rep(1, length(presence_ind)),
               rep(0, nrow(randomBg)))

  # create covariate training data by combining background with presence data
  coord_vars <- c("lon", "lat")
  trainEnv <- rbind(
    as(sample_data(neon_dob_prevalent_train), "data.frame")[presence_ind, c(terms3, coord_vars)],
    mutate(randomBg, mat_celsius_2 = mat_celsius^2, map_mm_2 = map_mm^2)[, c(terms3, coord_vars)]
  )

  # remove any records with missing observations
  keep_ind <- which(complete.cases(trainEnv))
  if(length(keep_ind) < nrow(trainEnv)) {
    if(verbose) message(paste0(sp, ": Removing ", nrow(trainEnv) - length(keep_ind), " records due to missing covariate data."))
  }
  trainEnv <- trainEnv[keep_ind,]
  trainPA <- trainPA[keep_ind]

  # tuned model
  model <- trainMaxNet(
    data = cbind(trainPA, trainEnv[,which(!names(trainEnv) %in% coord_vars)]),
    regMult = c(0.5, 1, 2, 3, 4, 5),
    classes = "l"
  )

  map <- predict(r_present_northam_expanded, model, type='cloglog')
  predictions <- raster::extract(map, trainEnv[,coord_vars])
  pred_at_pres <- predictions[trainPA == 1]
  pred_at_bg <- predictions[trainPA == 0]
  evaluation <- evaluate(p=as.vector(pred_at_pres), a=as.vector(pred_at_bg), tr=seq(0, 1, by=0.01))
  list(
    model,
    evaluation@t[which.max(evaluation@TPR + evaluation@TNR)]
  )
}

stopImplicitCluster()
t_maxnet1 <- Sys.time()
t_maxnet1 - t_maxnet0 # 50 minutes for all 8597 OTUs with 16 cores

maxnets_train <- oper[[1]]
maxnet_thresholds_train <- oper[[2]]

names(maxnets_train) <- names(maxnet_thresholds_train) <- taxa_names(neon_dob_prevalent_train)[spp_subset_ind]

# Predict using maxnets

registerDoParallel(N_CORES)
maxnet_preds <- foreach(i = 1:n_taxa) %dopar% {
# maxnet_preds <- foreach(i = 1:500) %dopar% {
  if(identical(maxnets_train[[i]], NA)) return(NA)
  prediction <- predict(maxnets_train[[i]],
                        newdata=as.matrix(newdata3),
                        type="cloglog", clamp=FALSE)
  classification <- prediction > maxnet_thresholds_train[[i]]
  as.numeric(classification)
}
stopImplicitCluster()

# Validate predictions

validation_predictions <- list(
  lognet_preds1,
  lognet_preds2,
  lognet_preds3,
  lognet_preds4,
  maxnet_preds
)
names(validation_predictions) <- c("climate_only", "climate_and_soil", "important_vars", "climate_and_ph", "maxnet_important_vars")

accuracy <- matrix(nrow=n_taxa, ncol=length(validation_predictions), dimnames=list(spp_subset_ind, names(validation_predictions)))
# for(i in 1:n_taxa) {
for(i in 1:500) {
  # sp <- taxa_names(neon_dob_prevalent_train)[spp_subset_ind[i]]
  sp <- taxa_names(neon_dob_prevalent_train)[i]
  y <- as.numeric(otu_table(neon_dob_prevalent_test)[,sp] > 0)
  for(mod in colnames(accuracy)) {
    accuracy[i,mod] <- mean(y == validation_predictions[[mod]][[i]], na.rm=TRUE)
  }
}

sensitivity <- specificity <- matrix(nrow=n_taxa, ncol=length(validation_predictions), dimnames=list(spp_subset_ind, names(validation_predictions)))
# for(i in 1:n_taxa) {
for(i in 1:500) {
  # sp <- taxa_names(neon_dob_prevalent_train)[spp_subset_ind[i]]
  sp <- taxa_names(neon_dob_prevalent_train)[i]
  y <- as.numeric(otu_table(neon_dob_prevalent_test)[,sp] > 0)
  if(all(c(0,1) %in% y)) {
    confusions <- lapply(
      validation_predictions,
      function(x) tryCatch(
        confusionMatrix(factor(x[[i]], levels=c(0,1)),
                        factor(y, levels=c(0,1)), positive="1"),
        error = function(e) return(NA))
    )
    sensitivity[i,] <- sapply(confusions, function(x) tryCatch(
      x$byClass["Sensitivity"], error = function(e) return(NA)))
    specificity[i,] <- sapply(confusions, function(x) tryCatch(
      x$byClass["Specificity"], error = function(e) return(NA)))
  } else {
    sensitivity[i,] <- NA
    specificity[i,] <- NA
  }
}

apply(accuracy, 2, mean, na.rm=TRUE)
apply(sensitivity, 2, mean, na.rm=TRUE)
apply(specificity, 2, mean, na.rm=TRUE)

apply(accuracy, 2, sd, na.rm=TRUE)
apply(sensitivity, 2, sd, na.rm=TRUE)
apply(specificity, 2, sd, na.rm=TRUE)

apply(sensitivity + specificity - 1, 2, mean, na.rm=TRUE) # True skill statistic
apply(sensitivity + specificity - 1, 2, sd, na.rm=TRUE) # True skill statistic


