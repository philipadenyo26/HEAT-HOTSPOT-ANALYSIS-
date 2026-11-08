# =====================================================================
# 01_eda.R  Exploratory data analysis (Chapter 3, Phase 1)
# Part A runs after Step 1; Part B runs after Step 3 (needs LST).
# =====================================================================
source("R/functions.R"); library(skimr)

# ---------- Part A: indices, population and wealth layers ----------
ndvi <- rast(file.path(GEE, "NDVI_2020_2023.tif"))
ndbi <- rast(file.path(GEE, "NDBI_2020_2023.tif"))
print(global(c(ndvi, ndbi), c("min", "max", "mean", "sd"), na.rm = TRUE))

# Plausible ranges: both indices must lie between -1 and 1
stopifnot(global(ndvi, "max", na.rm = TRUE)[1, 1] <= 1, global(ndbi, "min", na.rm = TRUE)[1, 1] >= -1)

ppp_file <- list.files(RAW, pattern = "gha_ppp_2020.*\\.tif$", full.names = TRUE)
pop <- if (length(ppp_file) > 0) rast(ppp_file[1]) else rast(file.path(GEE, "worldpop_total_2020.tif"))
pop <- crop(pop, project(gar_v, crs(pop)))
cat("WorldPop cells:", ncell(pop), " missing:", global(is.na(pop), "sum")[1, 1], "\n")
cat("Total population in extent:", round(global(pop, "sum", na.rm = TRUE)[1, 1]), "\n")

rwi <- read_csv(list.files(RAW, pattern = "relative_wealth.*\\.csv$", full.names = TRUE)[1], show_col_types = FALSE)
rwi_gar <- st_as_sf(rwi, coords = c("longitude", "latitude"), crs = 4326) |> st_transform(CRS_UTM) |> st_filter(gar)
print(skim(st_drop_geometry(rwi_gar)))

# ---------- Part B: LST distributions (run after Step 3) ----------
if (file.exists(file.path(GEE, "LST_2020_2023.tif"))) {
  lst <- rast(file.path(GEE, "LST_2020_2023.tif"))
  st <- c(lst, ndvi, ndbi); names(st) <- c("LST", "NDVI", "NDBI")
  df <- as.data.frame(st, na.rm = TRUE) |>
    mutate(cover = case_when(NDVI < 0 ~ "Water", NDVI > 0.4 ~ "Vegetated", NDBI > 0 ~ "Built-up", TRUE ~ "Mixed"))
  smry <- df |> group_by(cover) |>
    summarise(n = n(), mean = mean(LST), median = median(LST), sd = sd(LST),
              iqr = IQR(LST), p05 = quantile(LST, 0.05), p95 = quantile(LST, 0.95))
  write_csv(smry, file.path(TAB, "eda_lst_by_cover.csv")); print(smry)

  set.seed(1); dfs <- slice_sample(df, n = min(200000, nrow(df)))
  g1 <- ggplot(dfs, aes(LST, fill = cover)) + geom_histogram(bins = 60, alpha = 0.7, position = "identity") +
    labs(x = "LST (°C)", y = "Pixels") + theme_minimal()
  g2 <- ggplot(dfs, aes(cover, LST)) + geom_boxplot(outlier.size = 0.3) + labs(x = NULL, y = "LST (°C)") + theme_minimal()
  ggsave(file.path(MAP, "eda_lst_histogram.png"), g1, width = 7, height = 4, dpi = 300)
  ggsave(file.path(MAP, "eda_lst_boxplot.png"), g2, width = 6, height = 4, dpi = 300)

  ps <- read_csv(file.path(GEE, "lst_per_scene.csv"), show_col_types = FALSE) |> mutate(date = as.Date(date))
  g3 <- ggplot(ps, aes(date, lst_median)) + geom_point() + geom_smooth(se = FALSE, method = "loess") +
    labs(x = NULL, y = "Median LST over the region (°C)") + theme_minimal()
  ggsave(file.path(MAP, "eda_lst_per_scene.png"), g3, width = 7, height = 4, dpi = 300)
}
