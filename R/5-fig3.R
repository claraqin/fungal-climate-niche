# R code for Qin et al. on fungal niches and climate sensitivity
# 5-fig3.R: Make Figure 3 (univariate niche edges)

# Code originally from "niche_edges.R" and "esa_csee_figures.R"

neon_dob_wellfit <- readRDS("data/neon_dob_wellfit_v1.0.Rds")

terms_important <- c("mat_celsius", "temp_seasonality", "map_mm", "soilInCaClpH", "organicCPercent",
                     "mat_celsius_2", "map_mm_2")

# Fit lognet mdoels for all taxa (takes 1-2 mins with 16 cores)
lognets <- fitLognets(neon_dob_wellfit, ncores=N_CORES, randomseed=1010101,
                      terms = terms_important)

# Find thresholds for models
thresholds <- findThresholds(neon_dob_wellfit, lognets, terms = terms_important,
                             s="min", opt="maxtss", ncores=N_CORES)

# Get niche edges across gradients (MAT, MAP)

edge_mat_df_medians <- getNicheLimitsOnGradient(
  neon_dob_wellfit, lognets, thresholds, "mat_celsius",
  n_gradient_steps=100, terms = terms_important
)

edge_map_df_medians <- getNicheLimitsOnGradient(
  neon_dob_wellfit, lognets, thresholds, "map_mm",
  n_gradient_steps=100, terms = terms_important
)

# Combine with guild data

addGuildDataToNicheLimits <- function(limits_df, source_physeq, fungaltraits_df, guild_colname="guild2") {
  out <- limits_df %>%
    mutate(otu = as.character(sp), .keep="unused") %>%
    left_join(cbind.data.frame(otu = taxa_names(source_physeq), as(tax_table(source_physeq), "matrix"))) %>%
    left_join(dplyr::select(fungaltraits_df, genus, guild=one_of(guild_colname))) %>%
    mutate(guild = as.character(guild),
           guild = if_else(guild=="EM", "Ectomycorrhizal",
                           if_else(guild=="saprotroph", "Saprotrophic", guild)))
  out <- rbind(out, mutate(out, guild="All OTUs")) %>%
    dplyr::filter(guild %in% c("Ectomycorrhizal", "Saprotrophic", "All OTUs"))
}

edge_mat_byguild_medians <- addGuildDataToNicheLimits(edge_mat_df_medians, neon_dob_wellfit, ft)
edge_map_byguild_medians <- addGuildDataToNicheLimits(edge_map_df_medians, neon_dob_wellfit, ft)

# OTU counts by guild

otu_counts_medians_mat <- edge_mat_byguild_medians %>%
  dplyr::filter(!is.na(value)) %>%
  group_by(guild) %>%
  summarise(n = n_distinct(otu)) %>%
  mutate(label = paste0(guild, "\n(n = ", n, ")"))

otu_counts_medians_map <- edge_map_byguild_medians %>%
  dplyr::filter(!is.na(value)) %>%
  group_by(guild) %>%
  summarise(n = n_distinct(otu)) %>%
  mutate(label = paste0(guild, "\n(n = ", n, ")"))

# Proportion of points below min and above max

prop_truncated_mat <- edge_mat_byguild_medians %>%
  group_by(guild, edge) %>%
  summarise(n_otus = n_distinct(otu),
            n_at_bound = sum(value == if_else(edge=="lower", min(value, na.rm=TRUE),
                                              max(value, na.rm=TRUE)), na.rm=TRUE)) %>%
  mutate(p_at_bound = n_at_bound / n_otus) %>%
  dplyr::filter(guild %in% c("All OTUs", "Ectomycorrhizal", "Saprotrophic"))

prop_truncated_map <- edge_map_byguild_medians %>%
  group_by(guild, edge) %>%
  summarise(n_otus = n_distinct(otu),
            n_at_bound = sum(value == if_else(edge=="lower", min(value, na.rm=TRUE),
                                              max(value, na.rm=TRUE)), na.rm=TRUE)) %>%
  mutate(p_at_bound = n_at_bound / n_otus) %>%
  dplyr::filter(guild %in% c("All OTUs", "Ectomycorrhizal", "Saprotrophic"))

prop_truncated_mat
prop_truncated_map


# Plots -------------------------------------------------------------------

theme_set(theme_bw())

theme_dark <- theme(
  text = element_text(color="white"),
  panel.background = element_rect(fill = "transparent",
                                  colour = NA_character_), # necessary to avoid drawing panel outline
  panel.grid.major = element_blank(), # get rid of major grid
  panel.grid.minor = element_blank(), # get rid of minor grid
  plot.background = element_rect(fill = "transparent",
                                 colour = NA_character_), # necessary to avoid drawing plot outline
  legend.background = element_rect(fill = "transparent"),
  legend.box.background = element_rect(fill = "transparent"),
  legend.key = element_rect(fill = "transparent")
)

# Common plotting elements for gradient plots
common_elements_mat <- list(
  # geom_rect(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=-13, fill="grey90", col=NA),
  # geom_rect(xmin=-Inf, xmax=Inf, ymin=25.1, ymax=Inf, fill="grey90", col=NA),
  scale_y_continuous(limits=c(-13, 25.1)),
  scale_x_discrete(limits=rev, labels=setNames(otu_counts_medians_mat$label, otu_counts_medians_mat$guild)),
  scale_color_manual("Temperature\nniche edge", values=c("dodgerblue", "red"), labels=c("cold", "warm")),
  ylab("Mean annual temperature (°C)"),
  xlab(""),
  coord_flip(),
  theme_bw(),
  theme(legend.position="left", panel.grid = element_blank(), axis.title=element_text(size=10))
)

common_elements_map <- list(
  # geom_rect(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=116, fill="grey90", col=NA),
  # geom_rect(xmin=-Inf, xmax=Inf, ymin=2556, ymax=Inf, fill="grey90", col=NA),
  scale_x_discrete("", limits=rev, labels=setNames(otu_counts_medians_map$label, otu_counts_medians_map$guild)),
  scale_y_continuous("Annual precipitation (mm)", limits=c(116, 2556)),
  scale_color_manual("Precipitation\nniche edge", values=c("darkgoldenrod", "forestgreen"), labels=c("dry", "wet")),
  coord_flip(),
  theme_bw(),
  theme(legend.position="left", panel.grid = element_blank(),
        axis.title=element_text(size=10))
)

## MAT gradient plots ------------------------------------------------------

p_mat_gradient_medians_violins <-
  edge_mat_byguild_medians %>%
  dplyr::filter(value > -13 & value < 25.1) %>% # Filter within the min/max of MAT
  left_join(otu_counts_medians_mat) %>%
  ggplot(aes(y=value, x=guild, col=edge)) +
  geom_violin(width = 0.75, scale="count", fill="transparent") +
  common_elements_mat
p_gradient_mat_legend <- cowplot::get_legend(p_mat_gradient_medians_violins)
p_mat_gradient_medians_violins <- p_mat_gradient_medians_violins + guides(col="none")
p_mat_gradient_medians_violins

p_mat_gradient_medians_points <-
  edge_mat_byguild_medians %>%
  dplyr::filter(value > -13 & value < 25.1) %>% # Filter within the min/max of MAT
  left_join(otu_counts_medians_mat) %>%
  ggplot(aes(y=value, x=guild, col=edge)) +
  geom_point(size=0.1, alpha=0.5, shape=46,
             position=position_jitterdodge(jitter.width=0.7, jitter.height=0.2, dodge.width=0.75)) +
  common_elements_mat +
  guides(col = "none")
p_mat_gradient_medians_points

p_mat_truncated <- prop_truncated_mat %>%
  ggplot(aes(x=p_at_bound, y=guild, fill=edge)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values=c("dodgerblue", "red")) +
  scale_y_discrete("", limits=rev) +
  scale_x_continuous("Prop. truncated", limits=c(0, 1), labels = scales::percent) +
  guides(fill="none") +
  theme(axis.text.y = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(linetype = "dotted", color="gray"))
p_mat_truncated

p_mat_truncated_rev <- p_mat_truncated + scale_x_reverse("Prop. truncated", limits=c(1, 0), labels = scales::percent)
p_mat_truncated_rev

ggsave("plots/fig3/p_mat_violins_wellfit.pdf", p_mat_gradient_medians_violins, device="pdf", width=4, height=3, units="in")
ggsave("plots/fig3/p_mat_points_wellfit.pdf", p_mat_gradient_medians_points, device="pdf", width=4, height=3, units="in")
ggsave("plots/fig3/p_mat_points_wellfit.png", p_mat_gradient_medians_points, device="png", width=4, height=3, units="in")
ggsave("plots/fig3/p_mat_truncated_wellfit.pdf", p_mat_truncated, device="pdf", width=2, height=3, units="in")
ggsave("plots/fig3/p_mat_truncated_rev_wellfit.pdf", p_mat_truncated_rev, device="pdf", width=2, height=3, units="in")

p_mat_violins_truncated <- egg::ggarrange(
  p_mat_truncated_rev,
  p_mat_gradient_medians_violins,
  p_mat_truncated,
  ncol=3,
  widths=c(1, 2, 1)
)
ggsave("plots/fig3/p_mat_violins_and_truncated_wellfit.pdf", p_mat_violins_truncated, device="pdf", width=8, height=3, units="in")



## MAP gradient plots ------------------------------------------------------

p_map_gradient_medians_violins <-
  edge_map_byguild_medians %>%
  dplyr::filter(value > 237 & value < 2556) %>% # Filter within the min/max of MAP
  left_join(otu_counts_medians_map) %>%
  ggplot(aes(y=value, x=guild, col=edge)) +
  geom_violin(width = 0.75, scale="count") +
  common_elements_map
p_gradient_map_legend <- cowplot::get_legend(p_map_gradient_medians_violins)
p_map_gradient_medians_violins <- p_map_gradient_medians_violins + guides(col="none")
p_map_gradient_medians_violins

p_map_gradient_medians_points <-
  edge_map_byguild_medians %>%
  dplyr::filter(value > 237 & value < 2556) %>% # Filter within the min/max of MAP
  left_join(otu_counts_medians_map) %>%
  ggplot(aes(y=value, x=guild, col=edge)) +
  geom_point(size=0.1, alpha=0.5, shape=46,
             position=position_jitterdodge(jitter.width=0.7, jitter.height=12.5, dodge.width=0.75)) +
  common_elements_map +
  guides(col = "none")
p_map_gradient_medians_points

p_map_truncated <- prop_truncated_map %>%
  ggplot(aes(x=p_at_bound, y=guild, fill=edge)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values=c("darkgoldenrod", "forestgreen")) +
  scale_y_discrete("", limits=rev) +
  scale_x_continuous("Prop. truncated", limits=c(0, 1), labels = scales::percent) +
  guides(fill="none") +
  theme(axis.text.y = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(linetype = "dotted", color="gray"))
p_map_truncated

p_map_truncated_rev <- p_map_truncated + scale_x_reverse("Prop. truncated", limits=c(1, 0), labels = scales::percent)
p_map_truncated_rev

ggsave("plots/fig3/p_map_violins_wellfit.pdf", p_map_gradient_medians_violins, device="pdf", width=4, height=3, units="in")
ggsave("plots/fig3/p_map_points_wellfit.pdf", p_map_gradient_medians_points, device="pdf", width=4, height=3, units="in")
ggsave("plots/fig3/p_map_points_wellfit.png", p_map_gradient_medians_points, device="png", width=4, height=3, units="in")
ggsave("plots/fig3/p_map_truncated_wellfit.pdf", p_map_truncated, device="pdf", width=2, height=3, units="in")
ggsave("plots/fig3/p_map_truncated_rev_wellfit.pdf", p_map_truncated_rev, device="pdf", width=2, height=3, units="in")

p_map_violins_truncated <- egg::ggarrange(
  p_map_truncated_rev,
  p_map_gradient_medians_violins,
  p_map_truncated,
  ncol=3,
  widths=c(1, 2, 1)
)
ggsave("plots/fig3/p_map_violins_and_truncated_wellfit.pdf", p_map_violins_truncated, device="pdf", width=8, height=3, units="in")

## KS tests -------------

with(edge_mat_byguild_medians, ks.test(
  value[which(edge=="lower" & guild=="Ectomycorrhizal")],
  value[which(edge=="lower" & guild=="Saprotrophic")])) # 0.005

with(edge_mat_byguild_medians, ks.test(
  value[which(edge=="upper" & guild=="Ectomycorrhizal")],
  value[which(edge=="upper" & guild=="Saprotrophic")])) # 0.004

with(edge_map_byguild_medians, ks.test(
  value[which(edge=="lower" & guild=="Ectomycorrhizal")],
  value[which(edge=="lower" & guild=="Saprotrophic")])) # 0.413

with(edge_map_byguild_medians, ks.test(
  value[which(edge=="upper" & guild=="Ectomycorrhizal")],
  value[which(edge=="upper" & guild=="Saprotrophic")])) # 0.079

## Prop. absent from gradient -------------

edge_mat_byguild_medians %>%
  dplyr::filter(guild == "All OTUs") %>%
  group_by(otu) %>%
  summarise(isna = any(is.na(value))) %>%
  summarise(p_isna = mean(isna)) # 29.1%

edge_map_byguild_medians %>%
  dplyr::filter(guild == "All OTUs") %>%
  group_by(otu) %>%
  summarise(isna = any(is.na(value))) %>%
  summarise(p_isna = mean(isna)) # 30.0%



# Fig S2.1: Genus-level plots for top genera in each guild --------------------------


# Get top 10 most abundant genera in each guild
neon_dob <- readRDS("/data/ZHULAB/soil/NEON_DOB/phylo_V3.1.RDS")
neon_dob <- subset_samples(neon_dob, !is.na(lon) & !is.na(lat))
neon_dob <- subset_taxa(neon_dob, taxa_sums(neon_dob) > 0)
neon_dob <- subset_samples(neon_dob, sample_sums(neon_dob) > 0)
t0 <- Sys.time()
neon_dob_ra <- transform_sample_counts(neon_dob, function(x) x/sum(x)) # takes < 2 minutes
rm(neon_dob)
t0.1 <- Sys.time()
neon_dob_wellfit_genus <- tax_glom(prune_taxa(taxa_names(neon_dob_wellfit), neon_dob_ra), "genus") # takes 4.6 minutes
rm(neon_dob_ra)
t0.2 <- Sys.time()
diff(c(t0, t0.1, t0.2))
id_taxa_ind <- !grepl("unidentified", as.vector(tax_table(neon_dob_wellfit_genus)[,"genus"]))
neon_dob_wellfit_genus <- prune_taxa(id_taxa_ind, neon_dob_wellfit_genus)

tax_guild <- left_join(as.data.frame(as(tax_table(neon_dob_wellfit_genus), "matrix")),
                       dplyr::select(ft, genus, guild = guild2))
rownames(tax_guild) <- taxa_names(neon_dob_wellfit_genus)

neon_dob_wellfit_genus_em <- prune_taxa(rownames(tax_guild)[which(tax_guild$guild=="EM")], neon_dob_wellfit_genus)
neon_dob_wellfit_genus_sap <- prune_taxa(rownames(tax_guild)[which(tax_guild$guild=="saprotroph")], neon_dob_wellfit_genus)

abund_em <- cbind.data.frame(
  as(tax_table(neon_dob_wellfit_genus_em), "matrix"),
  mean_rel_abund = taxa_sums(neon_dob_wellfit_genus_em) / length(sample_sums(neon_dob_wellfit_genus_em))) %>%
  arrange(desc(mean_rel_abund))
abund_sap <- cbind.data.frame(
  as(tax_table(neon_dob_wellfit_genus_sap), "matrix"),
  mean_rel_abund = taxa_sums(neon_dob_wellfit_genus_sap) / length(sample_sums(neon_dob_wellfit_genus_sap))) %>%
  arrange(desc(mean_rel_abund))
head(abund_em, n=10)
head(abund_sap, n=10)

top_gen_em <- as.character(abund_em$genus[1:10])
top_gen_sap <- as.character(abund_sap$genus[1:10])

edge_mat_byguild_medians_top10 <-
  edge_mat_byguild_medians %>%
  mutate(genus = factor(genus, levels=c(top_gen_em, top_gen_sap))) %>%
  dplyr::filter(!is.na(genus), guild %in% c("Saprotrophic", "Ectomycorrhizal"))

genus_labels_df <- edge_mat_byguild_medians_top10 %>%
  dplyr::filter(!is.na(value)) %>%
  group_by(genus) %>%
  summarise(n = n_distinct(otu)) %>%
  mutate(label = paste0(genus, " (n = ", n, ")"))
genus_labels_df

edge_mat_byguild_medians_top10 %>%
  left_join(genus_labels_df) %>%
  ggplot(aes(y=value, x=genus, col=edge)) +
  facet_grid(guild~., scales="free") +
  geom_boxplot(position="dodge", fill=NA, size=0.3, outlier.color=NA) +
  geom_point(size=0.5, alpha=0.3, position=position_jitterdodge(jitter.width=0.2, jitter.height=0, dodge.width=0.75)) +
  scale_x_discrete(limits=rev, labels=setNames(genus_labels_df$label, genus_labels_df$genus)) +
  scale_color_manual("Niche edge", values=c("dodgerblue", "red"), labels=c("cold", "warm")) +
  ylab("Mean annual temperature (°C)") +
  xlab("") +
  coord_flip() +
  theme(legend.position = "top",
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major = element_line(colour="grey", linetype="dotted"))
ggsave("plots/fig-s2.1/fig-s2.1.png", device="png", width=5, height=5, units="in")



# Null-model distribution -------------------------------------------------

n_taxa <- length(taxa_names(neon_dob_wellfit))

n_iters <- 10

null_mat_edge_values <- null_map_edge_values <- matrix(nrow=n_taxa*2, ncol=n_iters)
null_mat_range <- null_map_range <- matrix(nrow=n_iters, ncol=2)

t1.0 <- Sys.time()
pb <- txtProgressBar(max=n_iters, style=3)
# for(i in 1:n_iters) {
for(i in 2:n_iters) {
  # Permute environmental data
  neon_dob_permuted <- permuteClimate(neon_dob_wellfit, randomseed=1010100+i)

  # Fit lognet mdoels for all taxa (takes 1-2 mins with 16 cores)
  lognets <- suppressMessages(fitLognets(neon_dob_permuted, ncores=N_CORES, randomseed=1010101, terms = terms_important))

  # Find thresholds for models
  thresholds <- suppressMessages(findThresholds(neon_dob_permuted, lognets, terms = terms_important, s="min", opt="maxtss", ncores=N_CORES))

  # Get niche edges across gradients (MAT, MAP)
  edge_mat_df_medians_null <- suppressMessages(getNicheLimitsOnGradient(
    neon_dob_permuted, lognets, thresholds, "mat_celsius",
    n_gradient_steps=100, terms = terms_important, progressbar=FALSE
  ))

  edge_map_df_medians_null <- suppressMessages(getNicheLimitsOnGradient(
    neon_dob_permuted, lognets, thresholds, "map_mm",
    n_gradient_steps=100, terms = terms_important, progressbar=FALSE
  ))

  # Save values
  null_mat_edge_values[,i] <- edge_mat_df_medians_null[,4]
  null_map_edge_values[,i] <- edge_map_df_medians_null[,4]
  null_mat_range[i,] <- range(get_variable(neon_dob_permuted, "mat_celsius"))
  null_map_range[i,] <- range(get_variable(neon_dob_permuted, "map_mm"))

  setTxtProgressBar(pb, i)
}
close(pb)
t1.1 <- Sys.time()
t1.1 - t1.0

null_mat_edges_l <- null_mat_edge_values[which(edge_mat_df_medians$edge=="lower"),]
null_mat_edges_u <- null_mat_edge_values[which(edge_mat_df_medians$edge=="upper"),]
null_map_edges_l <- null_map_edge_values[which(edge_map_df_medians$edge=="lower"),]
null_map_edges_u <- null_map_edge_values[which(edge_map_df_medians$edge=="upper"),]

# Count no. of edges at gradient boundary, and replace those with NA
temp <- sapply(1:n_iters, function(i) {
  if_else(null_mat_edges_l[,i]==null_mat_range[i,1], NA_real_, null_mat_edges_l[,i])
})
null_nas_mat_l <- apply(temp, 2, function(x) sum(is.na(x))) - apply(null_mat_edges_l, 2, function(x) sum(is.na(x)))
null_mat_edges_l <- temp

temp <- sapply(1:n_iters, function(i) {
  if_else(null_mat_edges_u[,i]==null_mat_range[i,2], NA_real_, null_mat_edges_u[,i])
})
null_nas_mat_u <- apply(temp, 2, function(x) sum(is.na(x))) - apply(null_mat_edges_u, 2, function(x) sum(is.na(x)))
null_mat_edges_u <- temp

temp <- sapply(1:n_iters, function(i) {
  if_else(null_map_edges_l[,i]==null_map_range[i,1], NA_real_, null_map_edges_l[,i])
})
null_nas_map_l <- apply(temp, 2, function(x) sum(is.na(x))) - apply(null_map_edges_l, 2, function(x) sum(is.na(x)))
null_map_edges_l <- temp

temp <- sapply(1:n_iters, function(i) {
  if_else(null_map_edges_u[,i]==null_map_range[i,2], NA_real_, null_map_edges_u[,i])
})
null_nas_map_u <- apply(temp, 2, function(x) sum(is.na(x))) - apply(null_map_edges_u, 2, function(x) sum(is.na(x)))
null_map_edges_u <- temp

# Number of truncated niches across null iterations (mean, SD)
null_n_truncated <- data.frame(
  gradient = c(rep("MAT", 2), rep("MAP", 2)),
  edge = c("cold", "warm", "dry", "wet"),
  mean_n_at_bound = unlist(lapply(
    list(null_nas_mat_l, null_nas_mat_u, null_nas_map_l, null_nas_map_u),
    function(x) mean(x, na.rm=TRUE))),
  sd_n_at_bound = unlist(lapply(
    list(null_nas_mat_l, null_nas_mat_u, null_nas_map_l, null_nas_map_u),
    function(x) sd(x, na.rm=TRUE)))
) %>%
  mutate(mean_p_at_bound = mean_n_at_bound / n_taxa,
         sd_p_at_bound = sd_n_at_bound / n_taxa)
null_n_truncated

# Distribution of niches across null iterations

null_edge_mat_df <- cbind(
  data.frame(
    sp = rep(taxa_names(neon_dob_wellfit), 2),
    edge = c(rep("cold", n_taxa), rep("warm", n_taxa))
  ),
  rbind(null_mat_edges_l, null_mat_edges_u)
)
null_edge_mat_long <- null_edge_mat_df %>%
  tidyr::pivot_longer(cols=-c(sp, edge), names_to="iter") %>%
  mutate(iter = factor(iter, levels=as.character(1:n_iters)))


null_edge_map_df <- cbind(
  data.frame(
    sp = rep(taxa_names(neon_dob_wellfit), 2),
    edge = c(rep("dry", n_taxa), rep("wet", n_taxa))
  ),
  rbind(null_map_edges_l, null_map_edges_u)
)
null_edge_map_long <- null_edge_map_df %>%
  tidyr::pivot_longer(cols=-c(sp, edge), names_to="iter") %>%
  mutate(iter = factor(iter, levels=as.character(1:n_iters)))


## Do single-axis niche edges have more clumping than random? ---------

vartomean(table(edge_mat_df_medians$value[which(!edge_mat_df_medians$value %in% c(-13, 25.1))]))
apply(null_edge_mat_df[,3:ncol(null_edge_mat_df)], 2, function(x) vartomean(table(x)))
mean(vartomean(table(edge_mat_df_medians$value[which(!edge_mat_df_medians$value %in% c(-13, 25.1))])) < apply(null_edge_mat_df[,3:ncol(null_edge_mat_df)], 2, function(x) vartomean(table(x))))
# Yes! MAT niche edges are more clumped than random

vartomean(table(edge_map_df_medians$value[which(!edge_map_df_medians$value %in% c(116, 2556))]))
apply(null_edge_map_df[,3:ncol(null_edge_map_df)], 2, function(x) vartomean(table(x)))
mean(vartomean(table(edge_map_df_medians$value[which(!edge_map_df_medians$value %in% c(116, 2556))])) < apply(null_edge_map_df[,3:ncol(null_edge_map_df)], 2, function(x) vartomean(table(x))), na.rm=TRUE)
# No! MAP niche edges are not significantly different from random



## Fig S2.2: Null-model distribution --------

# Common plotting elements for gradient plots
null_common_elements_mat <- list(
  scale_y_continuous(limits=c(-13, 25.1)),
  # scale_x_discrete(limits=rev, labels=setNames(otu_counts_medians_mat$label, otu_counts_medians_mat$guild)),
  scale_color_manual("Temperature\nniche edge", values=c("dodgerblue", "red"), labels=c("cold", "warm")),
  ylab("Mean annual temperature (°C)"),
  xlab(""),
  coord_flip(),
  theme_bw(),
  theme(legend.position="left", panel.grid = element_blank(), axis.title=element_text(size=10)),
  guides(col="none")
)

null_common_elements_map <- list(
  # scale_x_discrete("", limits=rev, labels=setNames(otu_counts_medians_map$label, otu_counts_medians_map$guild)),
  scale_y_continuous("Annual precipitation (mm)", limits=c(116, 2556)),
  scale_color_manual("Precipitation\nniche edge", values=c("darkgoldenrod", "forestgreen"), labels=c("dry", "wet")),
  xlab(""),
  coord_flip(),
  theme_bw(),
  theme(legend.position="left", panel.grid = element_blank(),
        axis.title=element_text(size=10)),
  guides(col="none")
)

p_null_mat_violins <- ggplot(null_edge_mat_long, aes(y=value, x=edge, col=edge, group=interaction(iter, edge))) +
  geom_violin(scale="count", fill=NA, position=position_identity(), lwd=0.05) +
  geom_violin(aes(group=NULL), scale="count", fill=NA) +
  null_common_elements_mat

p_null_map_violins <- ggplot(null_edge_map_long, aes(y=value, x=edge, col=edge, group=interaction(iter, edge))) +
  geom_violin(scale="count", fill=NA, position=position_identity(), lwd=0.05) +
  geom_violin(aes(group=NULL), scale="count", fill=NA) +
  null_common_elements_map

p_null_mat_truncated <- null_n_truncated %>%
  dplyr::filter(gradient=="MAT") %>%
  mutate(lower = mean_p_at_bound - sd_p_at_bound,
         upper = mean_p_at_bound + sd_p_at_bound) %>%
  ggplot(aes(x=mean_p_at_bound, y=edge, fill=edge)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width=0.5) +
  scale_fill_manual(values=c("dodgerblue", "red")) +
  scale_y_discrete("") +
  scale_x_continuous("Prop. truncated", limits=c(0, 1), labels = scales::percent) +
  guides(fill="none") +
  theme(axis.text.y = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(linetype = "dotted", color="gray"))

p_null_map_truncated <- null_n_truncated %>%
  dplyr::filter(gradient=="MAT") %>%
  mutate(lower = mean_p_at_bound - sd_p_at_bound,
         upper = mean_p_at_bound + sd_p_at_bound) %>%
  ggplot(aes(x=mean_p_at_bound, y=edge, fill=edge)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width=0.5) +
  scale_fill_manual(values=c("darkgoldenrod", "forestgreen")) +
  scale_y_discrete("") +
  scale_x_continuous("Prop. truncated", limits=c(0, 1), labels = scales::percent) +
  guides(fill="none") +
  theme(axis.text.y = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(linetype = "dotted", color="gray"))

p_null_violins_and_truncated <- egg::ggarrange(
  p_null_mat_violins, p_null_mat_truncated,
  p_null_map_violins, p_null_map_truncated,
  ncol = 2,
  widths = c(2,1),
  labels = c("(a)", "(b)", "(c)", "(d)"),
  label.args = list(gp=gpar(fontface="bold"))
)
ggsave("plots/fig-s2.2/null_edge_distr.png", p_null_violins_and_truncated,
       device="png", width=6, height=4, units="in")
ggsave("plots/fig-s2.2/null_edge_distr.pdf", p_null_violins_and_truncated,
       device="pdf", width=6, height=4, units="in")
