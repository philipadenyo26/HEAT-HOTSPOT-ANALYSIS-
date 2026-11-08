# =====================================================================
# 05_gwr.R  Geographically Weighted Regression of LST drivers
# Chapter 3, Section 3.9.2 (Objective 2)
# =====================================================================
source("R/functions.R"); library(GWmodel); library(car)

tmpl <- template(200)
lst  <- resample(rast(file.path(GEE, "LST_2020_2023.tif")),  tmpl, method = "average")
ndvi <- resample(rast(file.path(GEE, "NDVI_2020_2023.tif")), tmpl, method = "average")
ndbi <- resample(rast(file.path(GEE, "NDBI_2020_2023.tif")), tmpl, method = "average")
lcz  <- resample(rast(file.path(GEE, "LCZ_100m.tif")), tmpl, method = "mode")    # majority LCZ per cell
st <- mask(c(lst, ndvi, ndbi, lcz), gar_v); names(st) <- c("LST", "NDVI", "NDBI", "LCZ")
dat <- as.data.frame(st, xy = TRUE, na.rm = TRUE) |> filter(LCZ != 17)

# LCZ dummies with LCZ D (code 14) as reference. 'grouped' merges classes if local collinearity appears.
make_design <- function(d, grouped = FALSE) {
  cls <- as.character(d$LCZ)
  if (grouped) cls <- recode(cls, `2` = "compact", `3` = "compact", `5` = "open_mid", `6` = "open", `8` = "large_ind",
                             `10` = "large_ind", `9` = "sparse", `11` = "trees", `12` = "trees",
                             `16` = "bare", `14` = "14")
  f <- relevel(factor(cls), ref = "14")
  X <- model.matrix(~ f)[, -1, drop = FALSE]
  colnames(X) <- paste0("LCZ_", sub("^f", "", colnames(X)))
  bind_cols(d, as.data.frame(X))
}

fit_gwr <- function(dat, grouped = FALSE, drop = NULL, kernel = "bisquare") {
  dd <- make_design(dat, grouped)
  preds <- setdiff(c("NDVI", "NDBI", grep("^LCZ_", names(dd), value = TRUE)), drop)
  form <- as.formula(paste("LST ~", paste(preds, collapse = " + ")))
  set.seed(42)
  samp <- dd[sample(nrow(dd), min(5000, nrow(dd))), ]          # sample for bandwidth and model fit
  sp <- as(st_as_sf(samp, coords = c("x", "y"), crs = CRS_UTM), "Spatial")
  bw <- bw.gwr(form, data = sp, approach = "AICc", kernel = kernel, adaptive = TRUE)
  gw <- gwr.basic(form, data = sp, bw = bw, kernel = kernel, adaptive = TRUE)
  list(form = form, samp = samp, sp = sp, bw = bw, gw = gw, dd = dd, kernel = kernel)
}

# 1. Collinearity check on the full grid (Section 3.6)
ols_all <- lm(LST ~ ., data = make_design(dat) |> select(LST, NDVI, NDBI, starts_with("LCZ_")))
vifs <- vif(ols_all); write.csv(as.data.frame(vifs), file.path(TAB, "gwr_vif.csv")); print(vifs)
high_vif <- any(vifs[c("NDVI", "NDBI")] > 10)

# 2. GWR (separate NDVI and NDBI models if VIF > 10)
# Robust fit: bisquare kernel first; if a local matrix is singular (an LCZ class absent from some
# neighbourhoods), refit with a Gaussian kernel (weights never fall to zero); if that still fails,
# group related LCZ classes and use the Gaussian kernel.
robust_fit <- function(dat, drop = NULL) {
  for (opt in list(list(grouped = FALSE, kernel = "bisquare"), list(grouped = FALSE, kernel = "gaussian"),
                   list(grouped = TRUE, kernel = "gaussian"))) {
    fit <- try(fit_gwr(dat, grouped = opt$grouped, drop = drop, kernel = opt$kernel), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      message("GWR fitted with kernel = ", opt$kernel, ", grouped LCZ = ", opt$grouped)
      fit$grouped <- opt$grouped
      return(fit)
    }
    message("GWR attempt failed (kernel = ", opt$kernel, ", grouped = ", opt$grouped, "): trying the next option")
  }
  stop("GWR could not be fitted with any option")
}
runs <- if (high_vif) list(ndvi = robust_fit(dat, drop = "NDBI"), ndbi = robust_fit(dat, drop = "NDVI")) else list(main = robust_fit(dat))

for (k in names(runs)) {
  r <- runs[[k]]
  # OLS on the same sample for comparison by AICc
  ols <- lm(r$form, data = r$samp); n <- nrow(r$samp); p <- length(coef(ols)) + 1
  aicc_ols <- AIC(ols) + 2 * p * (p + 1) / (n - p - 1)
  # Residual spatial autocorrelation (8 nearest neighbours)
  nb <- knn2nb(knearneigh(sp::coordinates(r$sp), k = 8))
  mt <- moran.test(r$gw$SDF$residual, nb2listw(nb, style = "W"))
  mt_ols <- moran.test(residuals(ols), nb2listw(nb, style = "W"))   # same neighbours, OLS residuals
  write_csv(tibble(model = k, kernel = r$kernel, grouped_lcz = r$grouped, bandwidth_neighbours = r$bw,
                   AICc_OLS = aicc_ols, AICc_GWR = r$gw$GW.diagnostic$AICc,
                   R2_OLS = summary(ols)$r.squared, R2_GWR = r$gw$GW.diagnostic$gw.R2,
                   residual_MoranI = unname(mt$estimate[1]), residual_MoranI_p = mt$p.value,
                   OLS_residual_MoranI = unname(mt_ols$estimate[1]),
                   autocorrelation_removed_pct = 100 * (1 - unname(mt$estimate[1]) / unname(mt_ols$estimate[1]))),
            file.path(TAB, paste0("gwr_diagnostics_", k, ".csv")))

  # Coefficient surfaces at every grid cell, using the fitted bandwidth
  all_sp <- as(st_as_sf(r$dd, coords = c("x", "y"), crs = CRS_UTM), "Spatial")
  surf <- try(gwr.basic(r$form, data = r$sp, regression.points = all_sp, bw = r$bw, kernel = r$kernel, adaptive = TRUE), silent = TRUE)
  if (inherits(surf, "try-error")) {                     # singular at some map cells: Gaussian weights never reach zero
    message("Coefficient surface: switching to the Gaussian kernel")
    surf <- gwr.basic(r$form, data = r$sp, regression.points = all_sp, bw = r$bw, kernel = "gaussian", adaptive = TRUE)
  }
  sdf <- st_as_sf(surf$SDF)
  for (v in setdiff(names(sdf), "geometry")) {
    rr <- rasterize(vect(sdf), tmpl, field = v)
    writeRaster(rr, file.path(MAP, paste0("gwr_", k, "_", v, ".tif")), overwrite = TRUE)
  }
  st_write(st_as_sf(r$gw$SDF), file.path(MAP, paste0("gwr_", k, "_sample_fit.gpkg")), delete_dsn = TRUE, quiet = TRUE)
}
