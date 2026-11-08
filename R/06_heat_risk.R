# =====================================================================
# 06_heat_risk.R  Heat risk assessment and Critical Heat Risk Zones
# Chapter 3, Section 3.10 (Objective 3). Follows Obe et al. (2023).
# =====================================================================
source("R/functions.R"); library(gstat)
terraOptions(memfrac = 0.5)   # let terra work in chunks on disk rather than all in memory
gc()
tmpl <- template(100)
BAND <- 300                                   # Gi* distance band (m); every cell has neighbours
SENS_BANDS <- c(200, 500)                     # sensitivity bands (500 m needs the most memory; use c(200) if it stalls)

# ---- Hazard: LST 2020-2023 averaged to 100 m ----
H_raw <- resample(rast(file.path(GEE, "LST_2020_2023.tif")), tmpl, method = "average")

# ---- WorldPop 2020: downloaded files in data/raw if present, otherwise the GEE exports in data/gee ----
ppp_file <- list.files(RAW, pattern = "gha_ppp_2020.*\\.tif$", full.names = TRUE)
age_files <- list.files(RAW, pattern = "^gha_[fm]_[0-9]+_2020.*\\.tif$", full.names = TRUE)
aoi_ll <- buffer(gar_v, 5000)
if (length(ppp_file) > 0 && length(age_files) >= 36) {
  message("WorldPop: using downloaded files in data/raw")
  pop <- rast(ppp_file[1]); aoi <- project(aoi_ll, crs(pop)); pop <- crop(pop, aoi)
  age <- as.integer(sub("^gha_[fm]_([0-9]+)_2020.*$", "\\1", basename(age_files)))
  # Add the age rasters one at a time (loading all 36 together can exhaust memory)
  vul <- NULL; tot <- NULL
  for (i in seq_along(age_files)) {
    r <- crop(rast(age_files[i]), aoi); r <- ifel(is.na(r), 0, r)
    tot <- if (is.null(tot)) r else tot + r
    if (age[i] %in% c(0, 1, 5, 65, 70, 75, 80)) vul <- if (is.null(vul)) r else vul + r
    rm(r); if (i %% 6 == 0) gc()
  }
} else {
  message("WorldPop: using GEE exports in data/gee")
  pop <- rast(file.path(GEE, "worldpop_total_2020.tif")); aoi <- project(aoi_ll, crs(pop)); pop <- crop(pop, aoi)
  vul <- crop(rast(file.path(GEE, "worldpop_vulnerable_2020.tif")), aoi)
  tot <- pop
}

# ---- Exposure: population density (people per km2) ----
E_raw <- project(pop / cellSize(pop, unit = "km"), tmpl, method = "bilinear")

# ---- Vulnerable age group: under 10 (0, 1, 5) and 65+ (65, 70, 75, 80) ----
v_share_raw <- project(ifel(tot > 0, vul / tot, NA), tmpl, method = "bilinear")      # main: share
v_count_raw <- project(vul / cellSize(vul, unit = "km"), tmpl, method = "bilinear")  # sensitivity: count

# ---- Land mask: inside the region and not water (LCZ G = 17) ----
lcz <- resample(rast(file.path(GEE, "LCZ_100m.tif")), tmpl, method = "near")
land <- mask(ifel(lcz != 17, 1, NA), gar_v)

# ---- Relative Wealth Index: ordinary kriging to the 100 m grid ----
# Plain coordinate tables (no sf/stars objects inside gstat), duplicate locations removed,
# progress shown in the console. Kriged on a 500 m grid (RWI itself is ~2.4 km), then bilinear to 100 m.
rwi_sf <- read_csv(list.files(RAW, pattern = "relative_wealth.*\\.csv$", full.names = TRUE)[1], show_col_types = FALSE) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |> st_transform(CRS_UTM) |>
  st_filter(st_buffer(st_union(gar), 10000))
rwi_df <- data.frame(st_coordinates(rwi_sf), rwi = rwi_sf$rwi)
rwi_df <- rwi_df[!duplicated(rwi_df[, c("X", "Y")]), ]
message("RWI points used for kriging: ", nrow(rwi_df))
vg <- variogram(rwi ~ 1, locations = ~ X + Y, data = rwi_df)
vfit <- try(fit.variogram(vg, vgm(c("Exp", "Sph", "Gau"))), silent = TRUE)
if (inherits(vfit, "try-error") || any(vfit$range < 0) || isTRUE(attr(vfit, "singular"))) {
  message("Automatic variogram fit failed: using an exponential model with starting values from the data")
  vfit <- fit.variogram(vg, vgm(psill = var(rwi_df$rwi), model = "Exp", range = max(vg$dist) / 3, nugget = 0))
}
print(vfit)
KRIGE_RES <- 500
tmpl_k <- template(KRIGE_RES)
land_k <- mask(setValues(tmpl_k, 1), buffer(gar_v, 1000))
grid_df <- as.data.frame(land_k, xy = TRUE, na.rm = TRUE)[, c("x", "y")]; names(grid_df) <- c("X", "Y")
message("Kriging RWI at ", nrow(grid_df), " points (", KRIGE_RES, " m grid) - progress is shown below")
kr <- krige(rwi ~ 1, locations = ~ X + Y, data = rwi_df, newdata = grid_df, model = vfit, nmax = 30, debug.level = -1)
RWI_k <- rast(data.frame(x = kr$X, y = kr$Y, z = kr$var1.pred), type = "xyz", crs = CRS_UTM)
RWI_raw <- resample(RWI_k, tmpl, method = "bilinear")
png(file.path(MAP, "rwi_variogram.png"), width = 1600, height = 1200, res = 250); print(plot(vg, vfit)); dev.off()

# ---- Vegetation proportion (PV) from Step 3 ----
PV_raw <- resample(rast(file.path(GEE, "PV_2020_2023.tif")), tmpl, method = "average")

raw <- mask(c(H_raw, E_raw, v_share_raw, v_count_raw, RWI_raw, PV_raw, lcz), land)
names(raw) <- c("H", "E", "v", "v_count", "RWI", "PV", "LCZ")
writeRaster(raw, file.path(MAP, "risk_inputs_raw.tif"), overwrite = TRUE)
suppressWarnings(rm(list = intersect(ls(), c("ages", "vul", "tot", "pop", "kr", "grid_pts", "RWI_k", "land_k")))); invisible(gc())

# ---- Normalise (Equation 3.4), vulnerability (3.3) and heat risk (3.5) ----
n <- rast(lapply(c("H", "E", "v", "v_count", "RWI", "PV"), function(k) minmax(raw[[k]])))
names(n) <- c("H", "E", "v", "v_count", "RWI", "PV")
V   <- (n$v + (1 - n$RWI) + (1 - n$PV)) / 3
HRL <- n$H * n$E * V
out <- c(n$H, n$E, V, HRL); names(out) <- c("H_norm", "E_norm", "V_norm", "HRL")
writeRaster(out, file.path(MAP, "heat_risk_layers.tif"), overwrite = TRUE)

# ---- Table of all cells ----
df <- as.data.frame(c(out, n, raw$LCZ), xy = TRUE, na.rm = TRUE)
df <- df[, !duplicated(names(df))]
pts <- as.matrix(df[, c("x", "y")])
cells <- cellFromXY(tmpl, pts)                     # grid cell of each row, for the raster-based Gi*

# ---- Components and HRL by LCZ type, with Welch tests on a sample ----
by_lcz <- df |> group_by(LCZ) |> summarise(cells = n(), area_km2 = n() * 0.01,
  across(c(H_norm, E_norm, V_norm, HRL), list(mean = mean, sd = sd)), .groups = "drop")
write_csv(by_lcz, file.path(TAB, "risk_components_by_lcz.csv"))
smp <- sample_by_class(df |> mutate(LCZ = factor(LCZ)), "LCZ", n = 1000)
welch <- lapply(c("H_norm", "E_norm", "V_norm", "HRL"), function(k) {
  w <- oneway.test(as.formula(paste(k, "~ LCZ")), data = smp, var.equal = FALSE)
  tibble(layer = k, F = unname(w$statistic), p_value = w$p.value)
}) |> bind_rows()
write_csv(welch, file.path(TAB, "risk_welch_by_lcz.csv"))

# ---- Global Moran's I, then Gi* and Critical Heat Risk Zones ----
mi <- moran.test(df$HRL, nb2listw(dnearneigh(pts, 0, 150), style = "W", zero.policy = TRUE), zero.policy = TRUE)
df$gi <- gi_cells(df$HRL, cells, tmpl, BAND)        # Gi*, inverse distance within BAND, self included
df$gi_class <- gi_class(df$gi)
df$CHRZ <- df$gi >= 2.576
# False Discovery Rate (Benjamini-Hochberg) correction, as in Obe et al. (2023): reported as a sensitivity variant
df$p_fdr <- p.adjust(2 * pnorm(-abs(df$gi)), method = "BH")
df$CHRZ_fdr <- df$gi > 0 & df$p_fdr < 0.01

pts_sf <- st_as_sf(df[, c("x", "y")], coords = c("x", "y"), crs = CRS_UTM)
dist <- st_join(pts_sf, gar["admin2Name"])
df$district <- dist$admin2Name
miss <- which(is.na(df$district))                    # coastal cells just outside the district outlines
if (length(miss)) df$district[miss] <- gar$admin2Name[st_nearest_feature(pts_sf[miss, ], gar)]

write_csv(tibble(moran_I = unname(mi$estimate[1]), z = unname(mi$statistic), p_value = mi$p.value),
          file.path(TAB, "hrl_global_moran.csv"))
write_csv(df |> count(gi_class) |> mutate(area_km2 = n * 0.01), file.path(TAB, "hrl_gi_classes.csv"))
write_csv(df |> filter(CHRZ) |> count(LCZ) |> mutate(area_km2 = n * 0.01, share = n / sum(n)), file.path(TAB, "chrz_by_lcz.csv"))
write_csv(df |> filter(CHRZ) |> count(district) |> mutate(area_km2 = n * 0.01, share = n / sum(n)) |> arrange(desc(n)),
          file.path(TAB, "chrz_by_district.csv"))
cat("CHRZ area:", sum(df$CHRZ) * 0.01, "km2\n")

chrz_r <- rasterize(vect(df[df$CHRZ, ], geom = c("x", "y"), crs = CRS_UTM), tmpl, field = 1)
gi_r   <- rasterize(vect(df, geom = c("x", "y"), crs = CRS_UTM), tmpl, field = "gi")
writeRaster(c(gi_r, chrz_r), file.path(MAP, "hrl_gi_and_chrz.tif"), overwrite = TRUE)
saveRDS(df, file.path(MAP, "risk_cells.rds"))       # checkpoint: Step 7 can run even if the sensitivity part fails

# ---- Sensitivity analysis (Section 3.10.2) ----
chrz_of <- function(x, band = BAND) gi_cells(x, cells, tmpl, band) >= 2.576
main <- df$CHRZ
alt <- list(
  additive        = (df$H_norm + df$E_norm + df$V_norm) / 3,
  without_PV      = df$H_norm * df$E_norm * (df$v + (1 - df$RWI)) / 2,
  age_as_count    = df$H_norm * df$E_norm * (df$v_count + (1 - df$RWI) + (1 - df$PV)) / 3,
  without_age     = df$H_norm * df$E_norm * ((1 - df$RWI) + (1 - df$PV)) / 2,
  without_RWI     = df$H_norm * df$E_norm * (df$v + (1 - df$PV)) / 2)
sens <- lapply(names(alt), function(k) {
  ch <- chrz_of(alt[[k]])
  tibble(variant = k, chrz_km2 = sum(ch) * 0.01, share_of_main_retained = sum(ch & main) / sum(main))
}) |> bind_rows()
for (b in SENS_BANDS) {
  ch <- chrz_of(df$HRL, b)
  sens <- bind_rows(sens, tibble(variant = paste0("band_", b, "m"), chrz_km2 = sum(ch) * 0.01,
                                 share_of_main_retained = sum(ch & main) / sum(main)))
}
sens <- bind_rows(sens, tibble(variant = "fdr_corrected_99", chrz_km2 = sum(df$CHRZ_fdr) * 0.01,
                               share_of_main_retained = sum(df$CHRZ_fdr & main) / sum(main)))
write_csv(sens, file.path(TAB, "hrl_sensitivity.csv")); print(sens)
saveRDS(df, file.path(MAP, "risk_cells.rds"))
