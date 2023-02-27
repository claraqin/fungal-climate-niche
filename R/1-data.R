# R code for Qin et al. on fungal niches and climate sensitivity
# 1-data.R: Load and pre-process data for analysis


# 1. Prepare occurrence data ------------------------------------------------------------

neon_dob <- readRDS("/data/ZHULAB/soil/NEON_DOB/phylo_V3.1.RDS")

# We can't use any sites that lack coordinates
neon_dob <- subset_samples(neon_dob, !is.na(lon) & !is.na(lat))

# No need to agglomerate to a standard taxonomic level; all OTUs
# are already defined at 97% sequence similarity.

# Remove zero-abundance taxa and samples
neon_dob <- subset_taxa(neon_dob, taxa_sums(neon_dob) > 0)
neon_dob <- subset_samples(neon_dob, sample_sums(neon_dob) > 0)

# The removal of tag-swapping made no difference --
# there is no evidence of tag swapping.

# # Address potential tag-swapping
# # Take each OTU, calculate its max relative abundance in each sample.
# # Remove all occurrences of OTUs that appear in a sample below 0.1%
# # of its max relative abundance. This doesn't penalize rare species,
# # but species that appear in unusually low abundance relative to their
# # usual abundance.
# system.time(neon_dob_ra <- transform_sample_counts(neon_dob, function(x) x/sum(x))) # takes < 2 min
# range(sample_sums(neon_dob_ra))
#
# otu <- as(otu_table(neon_dob), "matrix")
# otu_ra <- as(otu_table(neon_dob_ra), "matrix")
#
# max_ra <- apply(otu_ra, 2, max)
#
# pb <- txtProgressBar(min = 0, max = length(taxa_names(neon_dob_ra)), style = 3)
# for(i in 1:ncol(otu)) {
#   rm_ind <- which(otu_ra[,i] > 0 & otu_ra[,i] < 0.001 * max_ra[i])
#   otu[rm_ind, i] <- 0
#   setTxtProgressBar(pb, i)
# }
# close(pb)
#
# neon_dob_ns <- neon_dob # ns = "no-tag-swapping"
# otu_table(neon_dob_ns) <- otu_table(otu, taxa_are_rows = FALSE)
#
# mean(sample_sums(neon_dob))
# mean(sample_sums(neon_dob_ns))
# mean(taxa_sums(neon_dob))
# mean(taxa_sums(neon_dob_ns))
#
# rm(neon_dob_ra)
# rm(neon_dob_ns)
# rm(otu)
# rm(max_ra)

# Aggregate sites within 10-minute of a degree; this is the resolution of the
# climate data
xgrid <- 1/6*seq(1,600)-160 # Based on longitudinal extent of sample data
ygrid <- 1/6*seq(1,360)+15 # Based on latitudinal extent of sample data
sample_data(neon_dob) %>%
  as("data.frame") %>%
  mutate(coordsSite = paste(round(xgrid[findInterval(lon, xgrid)], 2),
                            round(ygrid[findInterval(lat, ygrid)], 2), sep="_")) ->
  sample_data(neon_dob)
cbind.data.frame(
  coordsSite = get_variable(neon_dob, "coordsSite"),
  Site = get_variable(neon_dob, "Site"),
  Project = get_variable(neon_dob, "Project")
) %>%
  dplyr::filter(Project == "NEON") %>%
  group_by(Site) %>%
  summarise(n_coords = n_distinct(coordsSite)) %>%
  arrange(desc(n_coords)) # Number of operational sites that each NEON site gets split into
sites_to_project <- cbind.data.frame(
  coordsSite = get_variable(neon_dob, "coordsSite"),
  Project = get_variable(neon_dob, "Project")
) %>%
  group_by(coordsSite) %>%
  summarise(n_projects = n_distinct(Project),
            Project1 = as.character(Project[1])) %>%
  arrange(desc(n_projects)) %>%
  mutate(Project = if_else(n_projects == 2, "Both", Project1))
neon_dob_agg <- merge_samples(neon_dob, group="coordsSite")
sample_data(neon_dob_agg) %>%
  as("data.frame") %>%
  mutate(Project = sites_to_project$Project[match(sample_names(neon_dob_agg), sites_to_project$coordsSite)],
         lon = as.numeric(sapply(strsplit(sample_names(neon_dob_agg), split="_"), function(x) x[1])),
         lat = as.numeric(sapply(strsplit(sample_names(neon_dob_agg), split="_"), function(x) x[2])),
         coordsSite = rownames(.)) ->
  sample_data(neon_dob_agg)
# Will need to extract climate values from climate raster again.
# Do this in next section.

# 0-1 transform
neon_dob_agg <- transform_sample_counts(neon_dob_agg, function(x) ifelse(x>0, 1, 0))

# No. of sites that each OTU is present in
prevalence <- apply(otu_table(neon_dob_agg), 2, function(x) sum(x>0))
hist(prevalence)
sum(prevalence >= 20); mean(prevalence >= 20)
sum(prevalence >= 10); mean(prevalence >= 10)
sum(prevalence >= 5); mean(prevalence >= 5)

# Subset to OTUs present in at least 10 sites
neon_dob_prevalent <- prune_taxa(prevalence >= 10, neon_dob_agg)
# neon_dob_prevalent20 <- prune_taxa(prevalence >= 20, neon_dob_agg)



# 2. Get soil data for NEON sites -----------------------------------


# Code originally from "soil.R"

# Get soil physical and chemical properties, initial site characterization
soil <- loadByProduct(dpID="DP1.10047.001")

# Aggregate to coordinate-based sites
soil$spc_biogeochem %>%
  # Get soil horizons that are completely less than 30 cm in depth
  dplyr::filter(horizonID %in% soil$spc_perhorizon$horizonID[which(soil$spc_perhorizon$horizonBottomDepth < 30)]) %>%
  left_join(dplyr::select(soil$spc_perplot, plotID, decimalLongitude, decimalLatitude), by="plotID") %>%
  mutate(coordsSite = paste(round(xgrid[findInterval(decimalLongitude, xgrid)], 2),
                            round(ygrid[findInterval(decimalLatitude, ygrid)], 2), sep="_")) %>%
  group_by(coordsSite) %>%
  summarise(across(c(phCacl2, phH2o, nitrogenTot, carbonTot, ctonRatio), .fns = ~ mean(.x, na.rm=FALSE))) %>%
  mutate(across(c(nitrogenTot, carbonTot), .fns = ~ .x/10)) %>%
  # NOTE: Although we call it "soilInCaClpH", we actually use pH base in a water solution.
  # Similarly, although we call it "organicCPercent," we actually use the percent TOTAL carbon.
  dplyr::select(coordsSite, soilInCaClpH = phH2o, nitrogenPercent = nitrogenTot, organicCPercent = carbonTot) ->
  soilchem_coordsSites

sample_data(neon_dob_prevalent) %>%
  as("data.frame") %>%
  dplyr::filter(Project=="DoB") %>%
  group_by(coordsSite) %>%
  summarise(across(c(soilInCaClpH, nitrogenPercent, organicCPercent), .fns = ~ mean(.x, na.rm=TRUE))) -> soilchem_coordsSites_dob

soilchem_coordsSites <- rbind(soilchem_coordsSites,
                              soilchem_coordsSites_dob[which(!soilchem_coordsSites_dob$coordsSite %in% soilchem_coordsSites$coordsSite),])

# Add soil data to neon_dob_prevalent
sample_data(neon_dob_prevalent) %>% as("data.frame") %>%
  mutate(soilInCaClpH = soilchem_coordsSites$soilInCaClpH[match(.$coordsSite, soilchem_coordsSites$coordsSite)],
         nitrogenPercent = soilchem_coordsSites$nitrogenPercent[match(.$coordsSite, soilchem_coordsSites$coordsSite)],
         organicCPercent = soilchem_coordsSites$organicCPercent[match(.$coordsSite, soilchem_coordsSites$coordsSite)]) ->
  temp
apply(dplyr::select(temp, soilInCaClpH, nitrogenPercent, organicCPercent), 2, function(x) mean(!is.na(x)))
# soilInCaClpH nitrogenPercent organicCPercent
#    0.9375000       0.9609375       0.9609375
sample_data(neon_dob_prevalent) <- temp
sample_data(neon_dob_agg) <- temp


# 3. Prepare guild data ------------------------------------------------------


# First see if genus is unique, or if higher taxonomic ranks are necessary for matching
tax_table(neon_dob_prevalent) %>%
  as("matrix") %>% as.data.frame() %>%
  group_by(genus) %>%
  summarise(in_n_families = n_distinct(family)) %>%
  arrange(desc(in_n_families))
# Genus is unambiguous!

# Load FungalTraits
ft <- read.csv("data/FungalTraits_1.2_ver_16Dec_2020.csv", na.strings="")
ft <- ft[which(!duplicated(ft)),]
names(ft) <- tolower(names(ft))
names(ft)

# Is genus a unique identifier in FungalTraits?
dupes_genus <- ft$genus[which(duplicated(ft$genus))]
dupes_genus
ft[which(ft$genus %in% dupes_genus),c("class", "order", "family", "genus", "primary_lifestyle")]
as.character(dupes_genus) %in% tax_table(neon_dob_prevalent)[,"genus"]

# Yes! Thus, safe to remove duplicates.
ft <- ft[-which(duplicated(ft$genus)),]

# Create new column(s), representing broader groupings
guild2 <- rep(NA, nrow(ft))
guild2[grep("ectomycorrhizal", ft$primary_lifestyle)] <- "EM"
guild2[grep("arbuscular_mycorrhizal", ft$primary_lifestyle)] <- "AM"
guild2[grep("saprotroph", ft$primary_lifestyle)] <- "saprotroph"
guild2[grep("pathogen", ft$primary_lifestyle)] <- "pathogen"
guild2[grep("symbiotroph", ft$primary_lifestyle)] <- "symbiotroph"
guild2[grep("symbiont", ft$primary_lifestyle)] <- "symbiont"
guild2[grep("parasite", ft$primary_lifestyle)] <- "parasite"
guild2[grep("lichenized", ft$primary_lifestyle)] <- "lichenized"
guild2[grep("epiphyte", ft$primary_lifestyle)] <- "epiphyte"
guild2[grep("endophyte", ft$primary_lifestyle)] <- "endophyte"
guild2[is.na(guild2) & ft$primary_lifestyle!="unspecified"] <- "other"
ft$guild2 <- as.factor(guild2)
table(ft$guild2, useNA="ifany")

# Use a subset of trait data for this analysis
ft <- ft %>%
  dplyr::select(
    class:genus, primary_lifestyle, secondary_lifestyle,
    ectomycorrhiza_exploration_type_template,
    growth_form_template, guild2) %>%
  dplyr::filter(!is.na(primary_lifestyle))

# What are the consequences of subsetting the NEON-DoB dataset, for our
# ability to conduct trait analysis?
sum(tax_table(neon_dob_prevalent)[,"genus"] %in% ft$genus); mean(tax_table(neon_dob_prevalent)[,"genus"] %in% ft$genus)
# 5006 (58%) OTUs have genus-level guild annotation.


# 4. Load and pre-process climate data ------------------------------------


## 4.1. Load climate data ------------------------------------------------


# Load present (actually 1970-2000 average) climate data
r_present <- getData("worldclim",var="bio",res=10)
r_present <- r_present[[c(1,4,12,15)]]
names(r_present) <- c("mat_celsius","temp_seasonality","map_mm","prec_seasonality")

# Run necessary transformations on wordclim-provided temperature data
r_present$mat_celsius <- r_present$mat_celsius/10
r_present$temp_seasonality <- r_present$temp_seasonality/1000

# Crop climate data to study region
# Let's use North America between 18 and 72 degrees North, excluding Greenland
north_america <- ne_countries(continent="North America")
plot(north_america)
b <- as(extent(-170,-55,18,72), "SpatialPolygons")
crs(b) <- crs(north_america)
north_america_cropped <- raster::crop(north_america, b)
plot(north_america_cropped)

# Crop to region
r_present_northam <- raster::mask(raster::crop(r_present, north_america_cropped), north_america_cropped)

# Crop to remove Great Lakes
greatlakes <- rnaturalearth::ne_download(
  scale = 110, type = 'lakes', category = 'physical'
) %>%
  sf::st_as_sf(lakes110, crs = 4269) %>%
  dplyr::filter(name_en %in% c("Superior", "Michigan", "Huron", "Erie", "Ontario"))
clipOutPoly <- function(r, poly) {
  r_poly <- raster::mask(r, poly)
  r[which(!is.na(as.matrix(r_poly)))] <- NA
  r
}
r_present_northam <- clipOutPoly(r_present_northam, greatlakes)

# Load climate data into Phyloseq sample data
# Calculate decomposition coefficient (k) in sample data
addClimateToPhyseq <- function(physeq, r_climate) {
  data <- sample_data(physeq) %>% as("data.frame") %>%
    mutate(mat_celsius = raster::extract(r_climate[["mat_celsius"]], cbind(lon, lat)),
           map_mm = raster::extract(r_climate[["map_mm"]], cbind(lon, lat)),
           temp_seasonality = raster::extract(r_climate[["temp_seasonality"]], cbind(lon, lat)),
           prec_seasonality = raster::extract(r_climate[["prec_seasonality"]], cbind(lon, lat)))
  # If climate is NA for any sites, it's most likely because it falls just outside the raster cells.
  # Fill in with nearest cell.
  na_climate_ind <- which(is.na(data$mat_celsius))
  for(i in na_climate_ind) {
    xy <- cbind(data$lon[i], data$lat[i])
    nearest_ind <- which.min(replace(distanceFromPoints(r_climate[[1]], xy), is.na(r_climate[[1]]), NA))
    values <- r_climate@data@values[nearest_ind,]
    data$mat_celsius[i] <- values["mat_celsius"]
    data$map_mm[i] <- values["map_mm"]
    data$temp_seasonality[i] <- values["temp_seasonality"]
    data$prec_seasonality[i] <- values["prec_seasonality"]
  }
  # which(is.na(data$mat_celsius)) # All filled in!
  # data$k <- yasso_k(data$mat_celsius, data$map_mm)
  sample_data(physeq) <- data
  return(physeq)
}

neon_dob_prevalent <- addClimateToPhyseq(neon_dob_prevalent, r_present_northam)

# Get quadratic terms
sample_data(neon_dob_prevalent) %>%
  as("data.frame") %>%
  mutate(mat_celsius_2 = mat_celsius^2,
         map_mm_2 = map_mm^2) ->
  temp
sample_data(neon_dob_prevalent) <- temp
sample_data(neon_dob_agg) <- temp


# 5. Get soil rasters for spatial predictions and filling in missing data ----------------------------------------------------


# Data downloaded from ISRIC https://files.isric.org/soilgrids/latest/data_aggregated/5000m/phh2o/
# and https://files.isric.org/soilgrids/latest/data_aggregated/5000m/soc/

# pH in H2O

r_ph <- raster("phh2o_5-15cm_mean_5000.tif")
names(r_ph) <- "soilInCaClpH"
r_ph_reproj <- projectRaster(r_ph, crs = crs(r_present_northam))
r_ph_northam <- raster::mask(raster::crop(r_ph_reproj, north_america_cropped), north_america_cropped)
r_ph_northam_resample <- raster::resample(r_ph_northam, r_present_northam)
r_present_northam <- addLayer(r_present_northam, r_ph_northam_resample / 10)

rm(r_ph)
rm(r_ph_reproj)
rm(r_ph_northam)
rm(r_ph_northam_resample)

# Soil organic carbon content, a proxy for total percent carbon

r_soc <- raster("soc_5-15cm_mean_5000.tif")
names(r_soc) <- "organicCPercent"
r_soc_reproj <- projectRaster(r_soc, crs = crs(r_present_northam))
r_soc_northam <- raster::mask(raster::crop(r_soc_reproj, north_america_cropped), north_america_cropped)
r_soc_northam_resample <- raster::resample(r_soc_northam, r_present_northam)
plot(r_soc_northam_resample)
r_present_northam <- addLayer(r_present_northam, r_soc_northam_resample / 100)

rm(r_soc)
rm(r_soc_reproj)
rm(r_soc_northam)
rm(r_soc_northam_resample)

# Check for alignment between raster values and site-level values

soil_comparison_df <- data.frame(
  variable = c(rep("pH", nrow(sample_data(neon_dob_prevalent))), rep("SOC", nrow(sample_data(neon_dob_prevalent)))),
  val_neon_dob = c(get_variable(neon_dob_prevalent, "soilInCaClpH"), get_variable(neon_dob_prevalent, "organicCPercent")),
  val_raster = c(
    raster::extract(r_present_northam[["soilInCaClpH"]], cbind(get_variable(neon_dob_prevalent, "lon"), get_variable(neon_dob_prevalent, "lat"))),
    raster::extract(r_present_northam[["organicCPercent"]], cbind(get_variable(neon_dob_prevalent, "lon"), get_variable(neon_dob_prevalent, "lat"))))
)

ggplot(soil_comparison_df, aes(x = val_raster, y = val_neon_dob)) +
  facet_grid(variable~.) +
  geom_point() +
  geom_abline(slope=1, intercept=0)
# Fairly good match for soil pH; not a fantastic match for SOC.


# Fill in missing data

sum(is.na(get_variable(neon_dob_prevalent, "soilInCaClpH"))) # n=8
sum(is.na(get_variable(neon_dob_prevalent, "organicCPercent"))) # n=5
sum(is.na(get_variable(neon_dob_prevalent, "nitrogenPercent"))) # n=5
sum(is.na(get_variable(neon_dob_prevalent, "soilInCaClpH")) |
      is.na(get_variable(neon_dob_prevalent, "organicCPercent"))) # n=10

sample_data(neon_dob_prevalent) %>%
  as("data.frame") %>%
  # mutate(soilInCaClpH = raster::extract(r_present_northam[["soilInCaClpH"]], cbind(lon, lat)),
  #        organicCPercent = raster::extract(r_present_northam[["organicCPercent"]], cbind(lon, lat)), organicCPercent) ->
  mutate(soilInCaClpH = if_else(is.na(soilInCaClpH), raster::extract(r_present_northam[["soilInCaClpH"]], cbind(lon, lat)), soilInCaClpH),
         organicCPercent = if_else(is.na(organicCPercent), raster::extract(r_present_northam[["organicCPercent"]], cbind(lon, lat)), organicCPercent)) ->
  temp


sample_data(neon_dob_prevalent)[which(is.na(get_variable(neon_dob_prevalent, "soilInCaClpH"))),]
which(is.na(temp$organicCPercent))

# Some pH values are still missing because they're just slightly out of the bounds of the raster.
# Get the value of the nearest raster cell.
na_ind <- which(is.na(temp$soilInCaClpH))
for(i in na_ind) {
  xy <- cbind(temp$lon[i], temp$lat[i])
  nearest_ind <- which.min(replace(distanceFromPoints(r_present_northam[["soilInCaClpH"]], xy), is.na(r_present_northam[["soilInCaClpH"]]), NA))
  # values <- r_climate@data@values[nearest_ind,]
  values <- as.matrix(r_present_northam)[nearest_ind,]
  temp$soilInCaClpH[i] <- values["soilInCaClpH"]
}

sum(is.na(temp$soilInCaClpH))
sum(is.na(temp$organicCPercent))

sample_data(neon_dob_prevalent) <- temp
sample_data(neon_dob_agg) <- temp


### *Flag*: neon_dob_prevalent and neon_dob_agg objects are saved here --------

# saveRDS(neon_dob_prevalent, "data/neon_dob_prevalent_v4.1.Rds")
# saveRDS(neon_dob_agg, "data/neon_dob_agg_v2.1.Rds")

neon_dob_prevalent <- readRDS("data/neon_dob_prevalent_v4.1.Rds")
neon_dob_agg <- readRDS("data/neon_dob_agg_v2.1.Rds")


## 4.2. Get (inflated) convex hulls in climate space -------------------------------------------------


# Removed a subsection of code here that was used to compare present vs.
# future climate spaces, because we decided to use present climates only.
# This is still available in the original "main.R" file.

# Convex hull
# From Brian Steidinger
in.site.poly.inf <- function(sitedata, geodf, xvar, yvar, inflate = 0.25) {
  site.data = sitedata[, c(xvar, yvar)]
  x = site.data[,1]
  y = site.data[,2]
  Mx = mean(x)
  My = mean(y)
  CH=chull(site.data)
  BumpX = x[CH] + inflate*(x[CH]-Mx)
  BumpY = y[CH] + inflate*(y[CH]-My)
  geo.df = geodf[, c(xvar, yvar)]
  inpolygon(geo.df[,1], geo.df[,2], BumpX, BumpY, boundary = T)
}

# Wrapper function
getMaskFromClimateHull <- function(sitedata, geo_r, inflate = 0.15) {
  geo_df <- as.data.frame(cbind(coordinates(geo_r), as.matrix(geo_r)))
  names(geo_df)[1:2] <- c("lon", "lat")
  geo_df <- geo_df[complete.cases(geo_df),]

  # matrix for logical test of variable combinations
  var_comb <- c("MAT.TMPSEA", "MAT.MAP", "MAT.PRSEA", "MAP.TMPSEA", "MAP.PRSEA", "TMPSEA.PRSEA")
  var_mat <- matrix(nrow=nrow(geo_df), ncol=6)
  colnames(var_mat) <- var_comb
  var_mat <- data.frame(var_mat)

  # logical tests
  var_mat$MAT.TMPSEA <- in.site.poly.inf(sitedata, geo_df, "mat_celsius", "temp_seasonality", inflate)
  var_mat$MAT.MAP <- in.site.poly.inf(sitedata, geo_df, "mat_celsius", "map_mm", inflate)
  var_mat$MAT.PRSEA <- in.site.poly.inf(sitedata, geo_df, "mat_celsius", "prec_seasonality", inflate)
  var_mat$MAP.TMPSEA <- in.site.poly.inf(sitedata, geo_df, "map_mm", "temp_seasonality", inflate)
  var_mat$MAP.PRSEA <- in.site.poly.inf(sitedata, geo_df, "map_mm", "prec_seasonality", inflate)
  var_mat$TMPSEA.PRSEA <- in.site.poly.inf(sitedata, geo_df, "temp_seasonality", "prec_seasonality", inflate)

  # rowSums to see total number of polygons in
  var_mat$in.or.out <- rowSums(var_mat)

  # rasterize the points
  r_within <- rasterFromXYZ(cbind(geo_df$lon, geo_df$lat, var_mat$in.or.out))

  r_within_01 <- r_within
  r_within_01[r_within==6] <- 1
  r_within_01[r_within!=6] <- 0

  # how much area inclusive
  within_area_df <- data.frame(rasterToPoints(raster::area(r_within_01, na.rm=T)))
  within_area <- sum(within_area_df$layer)

  # Define mask based on hull
  climatehull_mask <- r_within_01
  climatehull_mask[climatehull_mask==0] <- NA

  return(list(mask = climatehull_mask, area = within_area, map = r_within_01))
}

# Wrapper function
getMaskFromEnvHull <- function(sitedata, geo_r,
                               terms_nonquadratic = c("mat_celsius", "map_mm", "temp_seasonality", "soilInCaClpH", "organicCPercent"),
                               inflate = 0.15) {
  geo_df <- as.data.frame(cbind(coordinates(geo_r), as.matrix(geo_r)))
  names(geo_df)[1:2] <- c("lon", "lat")
  geo_df <- geo_df[complete.cases(geo_df),]

  # matrix for logical test of variable combinations
  # var_comb <- c("MAT.TMPSEA", "MAT.MAP", "MAT.PRSEA", "MAP.TMPSEA", "MAP.PRSEA", "TMPSEA.PRSEA")
  var_mat <- matrix(nrow=nrow(geo_df), ncol=10)
  # colnames(var_mat) <- var_comb
  var_mat <- data.frame(var_mat)

  # logical tests
  inhull_list <- list()
  for(i in 1:(length(terms_nonquadratic)-1)) {
    for(j in (i+1):length(terms_nonquadratic)) {
      xvar <- terms_nonquadratic[i]
      yvar <- terms_nonquadratic[j]
      x <- sitedata[,xvar]
      y <- sitedata[,yvar]
      complete_ind <- which(complete.cases(cbind(x, y)))
      x <- x[complete_ind]
      y <- y[complete_ind]
      ch_ind <- chull(x, y)
      infl_x <- x[ch_ind] + inflate * (x[ch_ind] - mean(x))
      infl_y <- y[ch_ind] + inflate * (y[ch_ind] - mean(y))
      inhull_list <- c(inhull_list, list(inpolygon(geo_df[,xvar], geo_df[,yvar], infl_x, infl_y, boundary=TRUE)))
    }
  }

  # rowSums to see total number of polygons in
  inhull_bool <- rowSums(do.call(cbind, inhull_list)) == 10

  # rasterize the points
  r_within_01 <- rasterFromXYZ(cbind(geo_df$lon, geo_df$lat, inhull_bool))

  # how much area inclusive
  within_area_df <- data.frame(rasterToPoints(raster::area(r_within_01, na.rm=T)))
  within_area <- sum(within_area_df$layer)

  # Define mask based on hull
  climatehull_mask <- r_within_01
  climatehull_mask[climatehull_mask==0] <- NA

  return(list(mask = climatehull_mask, area = within_area, map = r_within_01))
}

# Convex hulls for present climates
obsEnv <- as(sample_data(neon_dob_prevalent), "data.frame")

# out_present1 <- getMaskFromClimateHull(obsEnv, r_present_northam, inflate = 0)
# out_present2 <- getMaskFromClimateHull(obsEnv, r_present_northam, inflate = 0.15)
# out_present3 <- getMaskFromClimateHull(obsEnv, r_present_northam, inflate = 0.25)
# out_present4 <- getMaskFromClimateHull(obsEnv, r_present_northam, inflate = 0.50)
# out_present5 <- getMaskFromClimateHull(obsEnv, r_present_northam, inflate = 0.75)
# out_present6 <- getMaskFromClimateHull(obsEnv, r_present_northam, inflate = 1.00)

out_present1 <- getMaskFromEnvHull(obsEnv, r_present_northam, inflate = 0)
out_present1.0 <- getMaskFromEnvHull(obsEnv, r_present_northam, inflate = 0.10)
out_present2 <- getMaskFromEnvHull(obsEnv, r_present_northam, inflate = 0.15)
out_present3 <- getMaskFromEnvHull(obsEnv, r_present_northam, inflate = 0.25)
out_present4 <- getMaskFromEnvHull(obsEnv, r_present_northam, inflate = 0.50)
out_present5 <- getMaskFromEnvHull(obsEnv, r_present_northam, inflate = 0.75)
out_present6 <- getMaskFromEnvHull(obsEnv, r_present_northam, inflate = 1.00)

# png(file="plots/inflation_hull_map_present.png", width=600, height=800)
# par(mfrow=c(3,2))
# plot(out_present1$map, main="Present analogous climates,\ninflation = 0", axes=FALSE, legend=FALSE)
# plot(out_present2$map, main="Present analogous climates,\ninflation = 0.15", axes=FALSE, legend=FALSE)
# plot(out_present3$map, main="Present analogous climates,\ninflation = 0.25", axes=FALSE, legend=FALSE)
# plot(out_present4$map, main="Present analogous climates,\ninflation = 0.50", axes=FALSE, legend=FALSE)
# plot(out_present5$map, main="Present analogous climates,\ninflation = 0.75", axes=FALSE, legend=FALSE)
# plot(out_present6$map, main="Present analogous climates,\ninflation = 1.00", axes=FALSE, legend=FALSE)
# par(mfrow=c(1,1))
# dev.off()

par(mfrow=c(3,1))
plot(out_present1$map, main="Present analogous environments,\ninflation = 0", axes=FALSE, legend=FALSE)
# plot(out_present1.0$map, main="Present analogous environments,\ninflation = 0.10", axes=FALSE, legend=FALSE)
plot(out_present2$map, main="Present analogous environments,\ninflation = 0.15", axes=FALSE, legend=FALSE)
plot(out_present3$map, main="Present analogous environments,\ninflation = 0.25", axes=FALSE, legend=FALSE)
par(mfrow=c(1,1))

env_mask_0.15 <- out_present2$mask
