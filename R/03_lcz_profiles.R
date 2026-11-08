# =====================================================================
# 03_lcz_profiles.R  LCZ thermal profiles and global map comparison
# Chapter 3, Phase 4 (Objective 1)
# =====================================================================
source("R/functions.R")

lcz_names <- c(`2` = "LCZ 2 Compact mid-rise", `3` = "LCZ 3 Compact low-rise", `5` = "LCZ 5 Open mid-rise", `6` = "LCZ 6 Open low-rise",
               `8` = "LCZ 8 Large low-rise", `9` = "LCZ 9 Sparsely built", `10` = "LCZ 10 Heavy industry",
               `11` = "LCZ A Dense trees", `12` = "LCZ B Scattered trees", `14` = "LCZ D Low plants",
               `16` = "LCZ F Bare soil", `17` = "LCZ G Water")

lcz <- rast(file.path(GEE, "LCZ_100m.tif"))
lst <- rast(file.path(GEE, "LST_2020_2023.tif"))
lst100 <- resample(lst, lcz, method = "average")            # mean LST of the 30 m pixels in each 100 m cell
st <- mask(c(lcz, lst100), gar_v); names(st) <- c("LCZ", "LST")

df <- as.data.frame(st, na.rm = TRUE) |> filter(LCZ != 17) |>
  mutate(LCZ = factor(lcz_names[as.character(LCZ)], levels = lcz_names))

# Area and mean LST per class (all cells)
profile <- df |> group_by(LCZ) |>
  summarise(cells = n(), area_km2 = n() * 0.01, mean_LST = mean(LST), sd_LST = sd(LST),
            median_LST = median(LST), .groups = "drop")
if (file.exists(file.path(TAB, "lst_calibration.csv"))) {      # add calibrated values (add-on; SC remains primary)
  cal <- read_csv(file.path(TAB, "lst_calibration.csv"), show_col_types = FALSE)
  profile <- profile |> mutate(mean_LST_cal = cal$a + cal$b * mean_LST, sd_LST_cal = cal$b * sd_LST,
                               median_LST_cal = cal$a + cal$b * median_LST)
}
write_csv(profile, file.path(TAB, "lcz_thermal_profile.csv")); print(profile)

# Welch tests on a random sample per class (up to 1,000 cells each)
smp <- sample_by_class(df, "LCZ", n = 1000) |> droplevels()
welch <- oneway.test(LST ~ LCZ, data = smp, var.equal = FALSE)
print(welch)
pw <- pairwise.t.test(smp$LST, smp$LCZ, pool.sd = FALSE, p.adjust.method = "none")
write_csv(tibble(statistic = "Welch F", value = unname(welch$statistic), df1 = welch$parameter[1],
                 df2 = welch$parameter[2], p_value = welch$p.value), file.path(TAB, "lcz_welch_anova.csv"))
write.csv(pw$p.value, file.path(TAB, "lcz_pairwise_welch_p.csv"))

g <- ggplot(smp, aes(LCZ, LST)) + geom_boxplot(outlier.size = 0.3, fill = "grey90") +
  coord_flip() + labs(x = NULL, y = "LST (°C), dry season 2020 to 2023") + theme_minimal()
ggsave(file.path(MAP, "lcz_lst_boxplot.png"), g, width = 7, height = 5, dpi = 300)

# Agreement with the global LCZ map (Demuzere et al., 2022)
# Uses lcz_local_vs_global.csv if it exists; otherwise compares LCZ_100m with global_LCZ_100m.tif directly.
if (file.exists(file.path(GEE, "lcz_local_vs_global.csv"))) {
  cmp <- read_csv(file.path(GEE, "lcz_local_vs_global.csv"), show_col_types = FALSE)
} else {
  glob <- resample(rast(file.path(GEE, "global_LCZ_100m.tif")), lcz, method = "near")
  pair <- mask(c(lcz, glob), gar_v); names(pair) <- c("local", "global")
  cmp <- as.data.frame(pair, na.rm = TRUE)
  set.seed(1); cmp <- cmp[sample(nrow(cmp), min(20000, nrow(cmp))), ]
}
agree <- mean(cmp$local == cmp$global, na.rm = TRUE)
xt <- table(local = lcz_names[as.character(cmp$local)], global = cmp$global)
write.csv(as.data.frame.matrix(xt), file.path(TAB, "lcz_local_vs_global_crosstab.csv"))
write_csv(tibble(points = nrow(cmp), overall_agreement = agree), file.path(TAB, "lcz_local_vs_global_agreement.csv"))
cat("Overall agreement with the global LCZ map:", round(100 * agree, 1), "%\n")

# ---- LST bias against the USGS product by LCZ (is the offset uniform across land cover?) ----
if (file.exists(file.path(GEE, "USGS_ST_2020_2023.tif"))) {
  lst30 <- rast(file.path(GEE, "LST_2020_2023.tif"))
  usgs30 <- resample(rast(file.path(GEE, "USGS_ST_2020_2023.tif")), lst30, method = "near")
  b100 <- resample(lst30 - usgs30, lcz, method = "average")
  db <- as.data.frame(mask(c(lcz, b100), gar_v), na.rm = TRUE); names(db) <- c("LCZ", "bias")
  bias_by_lcz <- db |> mutate(LCZ = lcz_names[as.character(LCZ)]) |> group_by(LCZ) |>
    summarise(mean_bias_C = mean(bias), sd_bias_C = sd(bias), cells = n(), .groups = "drop")
  write_csv(bias_by_lcz, file.path(TAB, "lst_bias_by_lcz.csv")); print(bias_by_lcz)
}
