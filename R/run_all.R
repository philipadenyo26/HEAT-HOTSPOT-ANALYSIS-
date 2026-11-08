# =====================================================================
# run_all.R  Runs the whole R analysis in order (Steps 1 to 7).
# Open RStudio with the thesis_analysis folder as the working directory, then:
#     source("R/run_all.R")
# Before each step it checks the input files. Steps whose inputs are missing are
# skipped (the missing files are listed); a step that fails does not stop the others,
# except that Step 7 needs Step 6. A summary is printed at the end.
# To run only some steps:  STEPS <- c(3, 5); source("R/run_all.R")
# =====================================================================
if (!exists("STEPS")) STEPS <- 0:8

g <- function(f) file.path("data/gee", f)
r <- function(f) file.path("data/raw", f)
has_worldpop <- function() length(list.files("data/raw", pattern = "gha_ppp_2020.*\\.tif$")) > 0 ||
  file.exists(g("worldpop_total_2020.tif"))
has_rwi <- function() length(list.files("data/raw", pattern = "relative_wealth.*\\.csv$")) > 0

steps <- list(
  list(n = 0, file = "R/00_setup.R",          title = "Setup (packages and folders)", needs = character(0)),
  list(n = 1, file = "R/01_eda.R",            title = "Exploratory data analysis",
       needs = c(r("GAR_districts.shp"), g("NDVI_2020_2023.tif"), g("NDBI_2020_2023.tif"))),
  list(n = 2, file = "R/02_uhi_validation.R", title = "UHI intensity and ERA5 validation",
       needs = c(g("LST_2015_2019.tif"), g("LST_2020_2023.tif"), g("NDVI_2015_2019.tif"), g("NDBI_2015_2019.tif"),
                 g("NDVI_2020_2023.tif"), g("NDBI_2020_2023.tif"), r("era5_tmax_dry_2020_2023.tif"))),   # USGS check runs too if its file is present
  list(n = 3, file = "R/03_lcz_profiles.R",   title = "LCZ thermal profiles and global map agreement",
       needs = c(g("LCZ_100m.tif"), g("LST_2020_2023.tif"))),
  list(n = 4, file = "R/04_hotspots.R",       title = "Getis-Ord Gi* and LISA hotspots",
       needs = c(g("LST_2020_2023.tif"))),
  list(n = 5, file = "R/05_gwr.R",            title = "Geographically Weighted Regression",
       needs = c(g("LST_2020_2023.tif"), g("NDVI_2020_2023.tif"), g("NDBI_2020_2023.tif"), g("LCZ_100m.tif"))),
  list(n = 6, file = "R/06_heat_risk.R",      title = "Heat risk and Critical Heat Risk Zones",
       needs = c(g("LST_2020_2023.tif"), g("PV_2020_2023.tif"), g("LCZ_100m.tif")), extra = "worldpop_rwi"),
  list(n = 7, file = "R/07_future_risk.R",    title = "Future heat risk (CMIP6)",
       needs = c(g("delta_ssp245_2041_2060.tif"), g("delta_ssp245_2081_2100.tif"),
                 g("delta_ssp585_2041_2060.tif"), g("delta_ssp585_2081_2100.tif"),
                 "outputs/maps/risk_cells.rds", "outputs/maps/risk_inputs_raw.tif")),
  list(n = 8, file = "R/08_figures.R",        title = "Story figures (each made only if its inputs exist)",
       needs = c(g("LST_2020_2023.tif")))
)

log <- data.frame(step = integer(), title = character(), status = character(), minutes = numeric(),
                  note = character(), stringsAsFactors = FALSE)

for (s in steps) {
  if (!(s$n %in% STEPS)) next
  missing <- s$needs[!file.exists(s$needs)]
  if (!is.null(s$extra)) {
    if (!has_worldpop()) missing <- c(missing, "WorldPop (data/raw/gha_ppp_2020*.tif or data/gee/worldpop_total_2020.tif)")
    if (!has_rwi()) missing <- c(missing, "Relative Wealth Index CSV (data/raw/*relative_wealth*.csv)")
  }
  if (s$n == 3 && !file.exists(g("lcz_local_vs_global.csv")) && !file.exists(g("global_LCZ_100m.tif")))
    missing <- c(missing, g("global_LCZ_100m.tif"))
  cat("\n==============================================================\n")
  cat(sprintf("STEP %d: %s\n", s$n, s$title))
  if (length(missing)) {
    cat("  SKIPPED - missing input files:\n"); cat(paste0("   - ", missing, "\n"), sep = "")
    log[nrow(log) + 1, ] <- list(s$n, s$title, "SKIPPED", NA, paste(basename(missing), collapse = "; "))
    next
  }
  t0 <- Sys.time()
  ok <- tryCatch({ source(s$file, local = new.env()); TRUE },
                 error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); conditionMessage(e) })
  mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  invisible(gc())                                  # free memory before the next step
  log[nrow(log) + 1, ] <- list(s$n, s$title, if (isTRUE(ok)) "DONE" else "ERROR", mins, if (isTRUE(ok)) "" else ok)
  cat(sprintf("  %s in %.1f min\n", if (isTRUE(ok)) "DONE" else "ERROR", mins))
}

cat("\n==================== SUMMARY ====================\n")
print(log, row.names = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(log, "outputs/tables/run_log.csv", row.names = FALSE)
cat("\nTables are in outputs/tables, maps in outputs/maps. Run log: outputs/tables/run_log.csv\n")
