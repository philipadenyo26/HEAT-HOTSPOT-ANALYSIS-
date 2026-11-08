# =====================================================================
# 02_uhi_validation.R  UHI intensity and ERA5 validation (Phase 3, Objective 1)
# =====================================================================
source("R/functions.R")

rural_districts <- c("Ada East", "Ada West", "Shai Osudoku", "Ningo/Prampram")
rural_v <- vect(gar[gar$admin2Name %in% rural_districts, ])

uhi_period <- function(period) {
  lst  <- rast(file.path(GEE, paste0("LST_", period, ".tif")))
  ndvi <- rast(file.path(GEE, paste0("NDVI_", period, ".tif")))
  ndbi <- rast(file.path(GEE, paste0("NDBI_", period, ".tif")))
  urban <- mask(lst, ndbi > 0, maskvalues = c(0, NA))
  rural <- mask(mask(lst, ndvi > 0.4, maskvalues = c(0, NA)), rural_v)
  t_urban <- global(urban, "mean", na.rm = TRUE)[1, 1]
  t_rural <- global(rural, "mean", na.rm = TRUE)[1, 1]
  uhii_map <- lst - t_rural                      # pixel-level intensity relative to the rural reference
  writeRaster(mask(uhii_map, gar_v), file.path(MAP, paste0("UHII_", period, ".tif")), overwrite = TRUE)
  tibble(period, mean_urban_lst = t_urban, mean_rural_lst = t_rural, UHII = t_urban - t_rural)
}
uhii <- bind_rows(uhi_period("2015_2019"), uhi_period("2020_2023"))

# Second definition: "urban" = built LCZ classes (2, 3, 5, 6, 8, 10), so dry bare plains (NDBI > 0) are excluded
BUILT <- c(2, 3, 5, 6, 8, 10)
if (file.exists(file.path(GEE, "LCZ_100m.tif"))) {
  lczr <- rast(file.path(GEE, "LCZ_100m.tif"))
  uhi_lcz <- function(period) {
    lst  <- rast(file.path(GEE, paste0("LST_", period, ".tif")))
    ndvi <- rast(file.path(GEE, paste0("NDVI_", period, ".tif")))
    l30 <- resample(lczr, lst, method = "near")
    urban <- mask(mask(lst, l30 %in% BUILT, maskvalues = c(0, NA)), gar_v)
    rural <- mask(mask(lst, ndvi > 0.4, maskvalues = c(0, NA)), rural_v)
    t_u <- global(urban, "mean", na.rm = TRUE)[1, 1]; t_r <- global(rural, "mean", na.rm = TRUE)[1, 1]
    tibble(period, mean_urban_lst_lcz = t_u, UHII_lcz = t_u - t_r)
  }
  uhii <- uhii |> left_join(bind_rows(uhi_lcz("2015_2019"), uhi_lcz("2020_2023")), by = "period")
}
write_csv(uhii, file.path(TAB, "uhii_by_period.csv")); print(uhii)

# ---- ERA5 validation: spatial consistency of anomaly patterns ----
lst  <- rast(file.path(GEE, "LST_2020_2023.tif"))
era5 <- rast(file.path(RAW, "era5_tmax_dry_2020_2023.tif"))
era5_30 <- project(era5, lst, method = "bilinear")         # bilinear interpolation to the 30 m grid
st <- mask(c(lst, era5_30), gar_v); names(st) <- c("lst", "era5")
d <- as.data.frame(st, na.rm = TRUE)
d <- d |> mutate(lst_anom = lst - mean(lst), era5_anom = era5 - mean(era5))
set.seed(7); ds <- slice_sample(d, n = min(100000, nrow(d)))
val <- tibble(pearson_r = cor(ds$lst_anom, ds$era5_anom),
              rmse_K = sqrt(mean((ds$lst_anom - ds$era5_anom)^2)),
              accept = pearson_r > 0.70 & rmse_K < 3.0)
write_csv(val, file.path(TAB, "lst_era5_validation.csv")); print(val)

# ---- Independent check against the USGS Landsat Level-2 surface temperature product (30 m) ----
if (file.exists(file.path(GEE, "USGS_ST_2020_2023.tif"))) {
  usgs <- resample(rast(file.path(GEE, "USGS_ST_2020_2023.tif")), lst, method = "near")
  pr <- mask(c(lst, usgs), gar_v); names(pr) <- c("lst", "usgs")
  dd <- as.data.frame(pr, na.rm = TRUE)
  set.seed(11); dd <- dd[sample(nrow(dd), min(200000, nrow(dd))), ]
  cross <- tibble(pearson_r = cor(dd$lst, dd$usgs),
                  rmse_C = sqrt(mean((dd$lst - dd$usgs)^2)),
                  mean_bias_C = mean(dd$lst - dd$usgs),
                  anomaly_r = cor(dd$lst - mean(dd$lst), dd$usgs - mean(dd$usgs)),
                  pixels = nrow(dd))
  write_csv(cross, file.path(TAB, "lst_usgs_crosscheck.csv")); print(cross)
  set.seed(12); ds <- dd[sample(nrow(dd), min(20000, nrow(dd))), ]
  gx <- ggplot(ds, aes(usgs, lst)) + geom_bin2d(bins = 80) + geom_abline(linetype = 2) +
    labs(x = "USGS Level-2 surface temperature (°C)", y = "Single-Channel LST (°C)") + theme_minimal()
  ggsave(file.path(MAP, "lst_vs_usgs.png"), gx, width = 6, height = 5, dpi = 300)
}

# ---- Add-on: linear calibration of the Single-Channel LST to the USGS product ----
# LST_cal = a + b * LST_SC. Fitted on half of the districts, tested on the other half (spatial hold-out),
# then refitted on all pixels. The Single-Channel LST remains the primary product; calibrated values are
# reported alongside. Hotspots, Gi*, normalised hazard and CHRZ are unchanged by a linear transform.
if (file.exists(file.path(GEE, "USGS_ST_2020_2023.tif"))) {
  dist_r <- rasterize(gar_v, lst, field = "admin2Name")
  pr2 <- mask(c(lst, usgs, dist_r), gar_v); names(pr2) <- c("lst", "usgs", "district")
  dd2 <- as.data.frame(pr2, na.rm = TRUE)
  set.seed(21); dd2 <- dd2[sample(nrow(dd2), min(300000, nrow(dd2))), ]
  set.seed(22); fit_d <- sample(unique(dd2$district), floor(length(unique(dd2$district)) / 2))
  fit_set <- dd2[dd2$district %in% fit_d, ]; test_set <- dd2[!dd2$district %in% fit_d, ]
  m_half <- lm(usgs ~ lst, data = fit_set)
  pred <- predict(m_half, newdata = test_set)
  holdout <- tibble(fit_districts = length(fit_d), test_pixels = nrow(test_set),
                    rmse_before_C = sqrt(mean((test_set$lst - test_set$usgs)^2)),
                    bias_before_C = mean(test_set$lst - test_set$usgs),
                    rmse_after_C = sqrt(mean((pred - test_set$usgs)^2)),
                    bias_after_C = mean(pred - test_set$usgs))
  m_all <- lm(usgs ~ lst, data = dd2)
  cal <- tibble(a = unname(coef(m_all)[1]), b = unname(coef(m_all)[2]), r2 = summary(m_all)$r.squared)
  write_csv(bind_cols(cal, holdout), file.path(TAB, "lst_calibration.csv")); print(bind_cols(cal, holdout))

  # Calibrated rasters (both periods; the 2015-2019 calibration assumes the same relationship holds)
  for (pp in c("2015_2019", "2020_2023")) {
    r <- rast(file.path(GEE, paste0("LST_", pp, ".tif")))
    writeRaster(cal$a + cal$b * r, file.path(MAP, paste0("LST_", pp, "_calibrated.tif")), overwrite = TRUE)
  }
  # UHII in calibrated units: a difference scales by b
  write_csv(uhii |> mutate(mean_urban_lst_cal = cal$a + cal$b * mean_urban_lst,
                           mean_rural_lst_cal = cal$a + cal$b * mean_rural_lst,
                           UHII_cal = cal$b * UHII) |>
              mutate(across(any_of("UHII_lcz"), ~ cal$b * .x, .names = "{.col}_cal")),
            file.path(TAB, "uhii_by_period_calibrated.csv"))
}
