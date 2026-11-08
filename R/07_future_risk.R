# =====================================================================
# 07_future_risk.R  Future surface heat hotspots and future heat risk
# Chapter 3, Section 3.11 (Objective 4)
# Needs: delta_*.tif from gee/04_cmip6_delta.js, and Step 6 outputs.
# =====================================================================
source("R/functions.R")
lst <- rast(file.path(GEE, "LST_2020_2023.tif"))
keys <- c("ssp245_2041_2060", "ssp245_2081_2100", "ssp585_2041_2060", "ssp585_2081_2100")

# ---- 1. Projected LST = observed LST + ensemble mean delta (bilinear to 30 m) ----
proj <- list(); dsum <- list()
for (k in keys) {
  d <- rast(file.path(GEE, paste0("delta_", k, ".tif")))           # band 1 = mean, band 2 = sd across models
  d30 <- project(d, lst, method = "bilinear")
  proj[[k]] <- lst + d30[[1]]
  writeRaster(proj[[k]], file.path(MAP, paste0("LST_projected_", k, ".tif")), overwrite = TRUE)
  dsum[[k]] <- tibble(projection = k,
                      mean_dT = global(mask(d30[[1]], gar_v), "mean", na.rm = TRUE)[1, 1],
                      mean_model_sd = global(mask(d30[[2]], gar_v), "mean", na.rm = TRUE)[1, 1],
                      mean_LST = global(mask(proj[[k]], gar_v), "mean", na.rm = TRUE)[1, 1])
}
write_csv(bind_rows(dsum), file.path(TAB, "future_delta_summary.csv"))

# ---- 2. Future surface heat hotspots (main Gi* specification, 200 m grid) ----
base <- run_hotspots(lst, 200, alternatives = FALSE)
fut  <- lapply(proj, run_hotspots, res = 200, alternatives = FALSE)
stopifnot(all(sapply(fut, nrow) == nrow(base)))                   # same grid, same cell order
hot <- sapply(fut, function(g) g$gi_main >= 1.96)
base$persistent <- rowSums(hot) == 4
base$scenario_contingent <- rowSums(hot[, grep("ssp585", keys)]) > 0 & rowSums(hot[, grep("ssp245", keys)]) == 0
p95 <- quantile(base$LST, 0.95)
hs_tab <- bind_rows(
  tibble(projection = "baseline_2020_2023", hot_km2 = sum(base$gi_main >= 1.96) * 0.04, cells_above_base_p95 = mean(base$LST > p95)),
  lapply(keys, function(k) tibble(projection = k, hot_km2 = sum(hot[, k]) * 0.04,
                                  cells_above_base_p95 = mean(fut[[k]]$LST > p95))) |> bind_rows(),
  tibble(projection = c("persistent", "scenario_contingent"),
         hot_km2 = c(sum(base$persistent), sum(base$scenario_contingent)) * 0.04))
write_csv(hs_tab, file.path(TAB, "future_hotspots_summary.csv"))
st_write(base, file.path(MAP, "future_hotspots_200m.gpkg"), delete_dsn = TRUE, quiet = TRUE)

# ---- 3. Future heat risk: projected hazard, baseline exposure and vulnerability ----
df  <- readRDS(file.path(MAP, "risk_cells.rds"))
xy  <- as.matrix(df[, c("x", "y")])
tmpl100 <- template(100)
cell_vals <- function(r) values(r)[cellFromXY(r, xy)]
H_base <- cell_vals(rast(file.path(MAP, "risk_inputs_raw.tif"))[["H"]])
H_fut  <- lapply(proj, function(r) cell_vals(resample(r, tmpl100, method = "average")))

# Joint normalisation across baseline and all four projections, so projected warming is kept
lo <- min(c(H_base, unlist(H_fut)), na.rm = TRUE); hi <- max(c(H_base, unlist(H_fut)), na.rm = TRUE)
EV <- df$E_norm * df$V_norm
hrl <- c(list(baseline = (H_base - lo) / (hi - lo) * EV), lapply(H_fut, function(h) (h - lo) / (hi - lo) * EV))

cells <- cellFromXY(tmpl100, xy)
chrz <- sapply(hrl, function(x) gi_cells(x, cells, tmpl100, 300) >= 2.576)   # raster-based Gi*, as in Step 6
df$chrz_persistent <- rowSums(chrz[, keys]) == 4
df$chrz_contingent <- rowSums(chrz[, grep("ssp585", keys)]) > 0 & rowSums(chrz[, grep("ssp245", keys)]) == 0

risk_tab <- tibble(projection = colnames(chrz), chrz_km2 = colSums(chrz) * 0.01,
                   mean_HRL = sapply(hrl, mean, na.rm = TRUE)) |>
  bind_rows(tibble(projection = c("persistent", "scenario_contingent"),
                   chrz_km2 = c(sum(df$chrz_persistent), sum(df$chrz_contingent)) * 0.01))
write_csv(risk_tab, file.path(TAB, "future_chrz_summary.csv")); print(risk_tab)
write_csv(df |> filter(chrz_persistent) |> count(district) |> mutate(area_km2 = n * 0.01) |> arrange(desc(n)),
          file.path(TAB, "future_chrz_persistent_by_district.csv"))
write_csv(df |> filter(chrz_persistent) |> count(LCZ) |> mutate(area_km2 = n * 0.01),
          file.path(TAB, "future_chrz_persistent_by_lcz.csv"))

for (k in colnames(chrz)) {
  r <- rasterize(vect(cbind(df[, c("x", "y")], v = as.integer(chrz[, k])), geom = c("x", "y"), crs = CRS_UTM), tmpl100, field = "v")
  writeRaster(r, file.path(MAP, paste0("CHRZ_", k, ".tif")), overwrite = TRUE)
}
