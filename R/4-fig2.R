# R code for Qin et al. on fungal niches and climate sensitivity
# 4-fig2.R: Make Figure 2 (sites in geographic space and climate space)

# Code originally from "permutations.R"

length(unique(get_variable(neon_dob, "Site")))
sites_to_project <- cbind.data.frame(
  Site = get_variable(neon_dob, "Site"),
  Project = get_variable(neon_dob, "Project")
) %>%
  group_by(Site) %>%
  summarise(n_projects = n_distinct(Project),
            Project1 = as.character(Project[1])) %>%
  arrange(desc(n_projects)) %>%
  mutate(Project = if_else(n_projects == 2, "Both", Project1))
neon_dob_site <- merge_samples(neon_dob, group="Site")
sites_df <- as(sample_data(neon_dob_site), "data.frame") %>%
  mutate(Project = sites_to_project$Project[match(sample_names(neon_dob_site), sites_to_project$Site)])

# Sites in geographic space
# p_sites_map <-
  ggplot(sites_df) +
  geom_sf(data=st_transform(st_as_sf(north_america_cropped), "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"),
  # geom_sf(data=st_as_sf(north_america_cropped),
          size=0.1, col="black", fill=alpha("white", 0)) +
  geom_point(aes(x=lon, y=lat, col=Project), shape=21, fill=alpha("white", 0.3), size=1) +
  scale_color_manual(values=c("#F8766D", "#619CFF")) +
  guides(col="none") +
  theme_custom_map

# Sites in climate space
p_sites_climate <- whittaker_base_plot() +
  geom_point(data=sites_df, aes(x=mat_celsius, y=map_mm/10, col=Project), shape=21, fill=alpha("white", 0.3), size=1) +
  xlab("Mean annual temperature (ºC)") +
  scale_y_continuous("Annual precipitation (mm)", labels=function(x) x*10) +
  scale_color_manual(values=c("#F8766D", "#619CFF")) +
  guides(col="none") +
  theme_custom +
  theme(legend.text = element_text(size=8))
  NULL
p_sites_climate

p_sites_map_climate <- egg::ggarrange(
  p_sites_map, p_sites_climate, ncol=2,
  widths = c(3, 2),
  labels=c("(a)", "(b)"),
  label.args = list(gp=gpar(fontface="bold"))
)
ggsave("plots/fig2/sites_map_and_climatespace.png", p_sites_map_climate, device="png", width=8, height=3, units="in")
ggsave("plots/fig2/sites_map_and_climatespace.pdf", p_sites_map_climate, device="pdf", width=8, height=3, units="in")
