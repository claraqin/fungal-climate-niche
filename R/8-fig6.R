# R code for Qin et al. on fungal niches and climate sensitivity
# 8-fig6.R: Make Figure 6 (climate sensitivity of fungi by biome)

# Code originally from "biomes.R"

sf::sf_use_s2(FALSE)

r_sor_mean_std <- projectMetricToRaster2(metrics_5d_inhull, "sorensen_mean_std", r_present_northam,
                                         terms_nonquadratic = terms_important_nonquadratic)[[1]]

r_sor_cv_std <- projectMetricToRaster2(metrics_5d_inhull, "sorensen_cv_std", r_present_northam,
                                       terms_nonquadratic = terms_important_nonquadratic)[[1]]

# temp <- cbind.data.frame(metrics_5d_inhull[,terms_important_nonquadratic],
#                          X = boot_metrics_5d[["sorensen"]][,1] / metrics_5d_inhull$sorensen_mean_permuted)
# r_sor_mean_std <- projectMetricToRaster2(temp, "X", r_present_northam,
#                                          terms_nonquadratic = terms_important_nonquadratic)[[1]]


# Read in WWF Biomes, available at https://www.sciencebase.gov/catalog/item/508fece8e4b0a1b43c29ca22
biomes <- st_read("data/wwf_biomes/wwf_terr_ecos.shp")

biomes <- group_by(biomes, BIOME) %>%
  summarise(geometry = st_union(geometry)) %>%
  ungroup() %>%
  mutate(BIOME = as.character(BIOME))

biome_labels <- data.frame(
  BIOME = as.character(c(seq(1, 14), 98, 99)),
  label = c(
    "Tropical & Subtropical Moist Broadleaf Forests",
    "Tropical & Subtropical Dry Broadleaf Forests",
    "Tropical & Subtropical Coniferous Forests",
    "Temperate Broadleaf & Mixed Forests",
    "Temperate Conifer Forests",
    "Boreal Forests/Taiga",
    "Tropical & Subtropical Grasslands, Savannas & Shrublands",
    "Temperate Grasslands, Savannas & Shrublands",
    "Flooded Grasslands & Savannas",
    "Montane Grasslands & Shrublands",
    "Tundra",
    "Mediterranean Forests, Woodlands & Scrub",
    "Deserts & Xeric Shrublands",
    "Mangroves",
    "Undefined",
    "Undefined2"
  ),
  stringsAsFactors = FALSE
)

biomes$LABEL <- biome_labels$label[match(biomes$BIOME, biome_labels$BIOME)]
cbind(biomes$BIOME, biomes$LABEL)

biomes <- st_crop(biomes, c(xmin=-170,xmax=-55,ymin=17,ymax=72))

ggplot() +
  geom_sf(data=biomes, aes(fill=LABEL))


groupRasterByBiomes <- function(r_metric, biome_geos, biome_names=NULL) {
  out <- lapply(biome_geos, function(x) {
    temp <- as.vector(raster::mask(r_metric, as(x, "Spatial")))
    temp[which(!is.na(temp))]
  })
  if(!is.null(biome_names)) names(out) <- biome_names
  out
}

groupedRasterToLongData <- function(metric_list, name="metric") {
  temp <- do.call(rbind.data.frame, lapply(seq_along(metric_list), function(i) {
    if(is.null(names(metric_list)[i])) {
      cbind(i, unlist(metric_list[[i]]))
    } else {
      cbind(names(metric_list)[i], unlist(metric_list[[i]]))
    }
  }))
  temp[,1] <- as.character(temp[,1])
  temp[,2] <- as.numeric(as.character(temp[,2]))
  names(temp) <- c("label", name)
  temp
}

t0 <- Sys.time()

r_sor_mean_std %>%
  groupRasterByBiomes(biomes$geometry, biomes$LABEL) %>%
  groupedRasterToLongData(name="sorensen_mean") ->
  sorstd_by_biome_df

r_sor_cv_std %>%
  groupRasterByBiomes(biomes$geometry, biomes$LABEL) %>%
  groupedRasterToLongData(name="sorensen_cv") ->
  sorstd_cv_by_biome_df

t0.1 <- Sys.time()

area_list <- area(r_sor_mean_std) %>%
  groupRasterByBiomes(biomes$geometry)
area_biomes <- lapply(area_list, sum)
names(area_biomes) <- biomes$LABEL

large_biomes <- names(area_biomes)[which(area_biomes > 500000)]

t0.2 <- Sys.time()

median_sorstd_by_biome <-
  sorstd_by_biome_df %>%
  group_by(label) %>%
  summarise(median = median(sorensen_mean, na.rm=TRUE)) %>%
  arrange(median) %>%
  # mutate(BIOME = as.character(BIOME)) %>%
  left_join(biome_labels)
median_sorstd_by_biome

median_sorstd_cv_by_biome <-
  sorstd_cv_by_biome_df %>%
  group_by(label) %>%
  summarise(median = median(sorensen_cv, na.rm=TRUE)) %>%
  arrange(median) %>%
  # mutate(BIOME = as.character(BIOME)) %>%
  left_join(biome_labels)
median_sorstd_cv_by_biome

t0.3 <- Sys.time()

sorstd_by_biome_df %>%
  left_join(biome_labels) %>%
  dplyr::filter(!(BIOME %in% c("98", "99"))) %>%
  mutate(label = if_else(label %in% large_biomes, label, "Other")) %>%
  left_join(median_sorstd_by_biome) %>%
  mutate(BIOME = factor(BIOME, levels=median_sorstd_by_biome$BIOME),
         label = factor(label, levels=median_sorstd_by_biome$label[which(median_sorstd_by_biome$label %in% label)])) %>%
  dplyr::filter(!is.na(label)) ->
  sorstd_by_biome_df2

sorstd_cv_by_biome_df %>%
  left_join(biome_labels) %>%
  dplyr::filter(!(BIOME %in% c("98", "99"))) %>%
  mutate(label = if_else(label %in% large_biomes, label, "Other")) %>%
  left_join(median_sorstd_cv_by_biome) %>%
  mutate(BIOME = factor(BIOME, levels=median_sorstd_by_biome$BIOME),
         label = factor(label, levels=median_sorstd_by_biome$label[which(median_sorstd_by_biome$label %in% label)])) %>%
  dplyr::filter(!is.na(label)) ->
  sorstd_cv_by_biome_df2

t0.4 <- Sys.time()

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

sorstd_biome <- ggplot(sorstd_by_biome_df2, aes(x=sorensen_mean, fill=BIOME)) +
  facet_wrap(~label, ncol=1, scales="free_y") +
  geom_histogram() +
  geom_vline(aes(xintercept = median)) +
  scale_fill_manual("", values = rev(gg_color_hue(7))) +
  guides(fill="none") +
  ylab("") + xlab("Mean Sørensen sensitivity, std.") +
  theme_classic() +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(hjust=0, size=7))
# ggsave("plots/fig6/sorstd_mean_by_biome.png", sorstd_biome, device="png", width=2.5, height=5, units="in")

sorstd_cv_biome <- ggplot(sorstd_cv_by_biome_df2, aes(x=sorensen_cv, fill=BIOME)) +
  facet_wrap(~label, ncol=1, scales="free_y") +
  geom_histogram() +
  geom_vline(aes(xintercept = median)) +
  scale_fill_manual("", values = rev(gg_color_hue(7))) +
  guides(fill="none") +
  ylab("") + xlab("CV of Sørensen sensitivity, std.") +
  theme_classic() +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(hjust=0, size=7))
# ggsave("plots/fig6/sorstd_cv_by_biome.png", sorstd_cv_biome, device="png", width=2.5, height=5, units="in")


# Add reference map

biomes_main <- biomes %>%
  dplyr::filter(!(BIOME %in% c("98", "99"))) %>%
  mutate(LABEL = if_else(LABEL %in% large_biomes, LABEL, "Other"),
         LABEL = factor(LABEL, levels=c(median_sorstd_by_biome$label[which(median_sorstd_by_biome$label %in% LABEL)], "Other")))
# BIOME = if_else(BIOME %in% c("98", "99"), NA_character_, BIOME),
# BIOME = factor(BIOME, levels = unique(BIOME)))
cbind(biomes_main$BIOME, biomes_main$LABEL)
biome_map <- ggplot() +
  geom_sf(data=biomes_main, aes(fill=LABEL), lwd=0) +
  scale_fill_manual("", values=rev(gg_color_hue(7))) +
  guides(fill="none") +
  theme_custom_map
# ggsave("plots/fig6/biome_map.png", biome_map, device="png", width=3, height=3, units="in")
# ggsave("plots/fig6/biome_map.pdf", biome_map, device="pdf", width=3, height=3, units="in")

sorstd_mean_biome_fig <- egg::ggarrange(biome_map, sorstd_biome, ncol=1, heights = c(1,2))
ggsave("plots/fig6/sorstd_mean_biome_fig.png", sorstd_mean_biome_fig, device="png", width=2.5, height=6, units="in")
ggsave("plots/fig6/sorstd_mean_biome_fig.pdf", sorstd_mean_biome_fig, device="pdf", width=2.5, height=6, units="in")

sorstd_cv_biome_fig <- egg::ggarrange(biome_map, sorstd_cv_biome, ncol=1, heights = c(1,2))
ggsave("plots/fig6/sorstd_cv_biome_fig.png", sorstd_cv_biome_fig, device="png", width=2.5, height=6, units="in")
ggsave("plots/fig6/sorstd_cv_biome_fig.pdf", sorstd_cv_biome_fig, device="pdf", width=2.5, height=6, units="in")
