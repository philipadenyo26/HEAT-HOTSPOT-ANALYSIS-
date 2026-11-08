# =====================================================================
# 09_gwr_summary.R  Summary table and maps of the GWR coefficients (Objective 2)
# Reads the coefficient rasters written by Step 5 (outputs/maps/gwr_main_*.tif).
# Writes outputs/tables/gwr_coefficient_summary.csv and outputs/figures/F5b_gwr_coefficients.png
# =====================================================================
source("R/functions.R")
suppressPackageStartupMessages(library(patchwork))
files <- list.files(MAP, pattern = "^gwr_main_.*\\.tif$", full.names = TRUE)
files <- files[!grepl("residual|Local_R2|yhat|_SE\\.tif|_TV\\.tif|_y\\.tif|sample_fit", files)]
summ <- lapply(files, function(f) {
  v <- values(mask(rast(f), gar_v), na.rm = TRUE)
  tibble(term = sub("^gwr_main_(.*)\\.tif$", "\\1", basename(f)),
         median = median(v), p25 = quantile(v, 0.25), p75 = quantile(v, 0.75),
         min = min(v), max = max(v), share_negative = mean(v < 0), cells = length(v))
}) |> bind_rows()
write_csv(summ, file.path(TAB, "gwr_coefficient_summary.csv")); print(summ, n = 30)

r2f <- list.files(MAP, pattern = "^gwr_main_Local_R2\\.tif$", full.names = TRUE)
if (length(r2f)) {
  v <- values(mask(rast(r2f), gar_v), na.rm = TRUE)
  write_csv(tibble(median_local_R2 = median(v), p10 = quantile(v, 0.1), p90 = quantile(v, 0.9)),
            file.path(TAB, "gwr_local_r2_summary.csv"))
}
mapit <- function(term, title) {
  f <- file.path(MAP, paste0("gwr_main_", term, ".tif")); if (!file.exists(f)) return(NULL)
  d <- as.data.frame(mask(rast(f), gar_v), xy = TRUE, na.rm = TRUE); names(d)[3] <- "b"
  lim <- max(abs(quantile(d$b, c(0.02, 0.98))))
  ggplot() + geom_raster(data = d, aes(x, y, fill = b)) + geom_sf(data = gar, fill = NA, colour = "grey25", linewidth = 0.1) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-lim, lim), oob = scales::squish, name = "Coefficient") +
    labs(title = title) + coord_sf(crs = CRS_UTM) + theme_minimal(base_size = 10) +
    theme(axis.title = element_blank(), axis.text = element_blank(), panel.grid = element_blank())
}
p <- Filter(Negate(is.null), list(mapit("NDVI", "NDVI"), mapit("NDBI", "NDBI"), mapit("LCZ_3", "LCZ 3 compact low rise (vs low plants)"),
                                   mapit("LCZ_6", "LCZ 6 open low rise (vs low plants)")))
if (length(p)) ggsave(file.path("outputs/figures", "F5b_gwr_coefficients.png"), wrap_plots(p, ncol = 2), width = 11, height = 8, dpi = 300, bg = "white")
message("GWR summary done")
