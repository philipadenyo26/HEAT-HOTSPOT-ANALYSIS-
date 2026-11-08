# =====================================================================
# 08_figures.R  Figures that tell the story of heat in Greater Accra
#   1 Where it is hot      F1 LST map, F2 change between periods, F3 LCZ map
#   2 Why it is hot        F4 LST by LCZ, F5 vegetation cools
#   3 Where heat clusters  F6 surface heat hot spots
#   4 Who is at risk       F7 Critical Heat Risk Zones, F8 districts, F9 people, F10 Accra vs Lagos
#   5 What is coming       F11 future warming and risk, F12 projected LST maps
# Each figure is made only if its input files exist. Output: outputs/figures/*.png (300 dpi)
# =====================================================================
source("R/functions.R")
suppressPackageStartupMessages({ library(scales); library(patchwork); library(ggrepel) })
FIG <- "outputs/figures"; dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
has <- function(...) all(file.exists(c(...)))
save_fig <- function(p, name, w = 8, h = 6) {
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  message("  saved ", name)
}
try_fig <- function(name, expr) tryCatch(expr, error = function(e) message("  ", name, " skipped: ", conditionMessage(e)))

# ---- Common look ----
theme_story <- theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(colour = "grey30"),
        plot.caption = element_text(colour = "grey45", size = 8, hjust = 0), panel.grid.minor = element_blank(),
        legend.position = "right", plot.title.position = "plot", plot.caption.position = "plot")
theme_map <- theme_story + theme(axis.title = element_blank(), panel.grid = element_line(colour = "grey92"))
heat_pal <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c", "#67001f")
src_ls <- "Source: Landsat 8/9 dry-season composites (Nov-Mar), Single-Channel LST; boundaries: OCHA."
lcz_col <- c(`2` = "#d10000", `3` = "#ff0000", `5` = "#ff6600", `6` = "#ff9955", `8` = "#bcbcbc", `9` = "#ffccaa",
             `10` = "#555555", `11` = "#006a00", `12` = "#00aa00", `14` = "#b9db79", `16` = "#fbf7ae", `17` = "#6a6aff")
lcz_lab <- c(`2` = "2 Compact mid-rise", `3` = "3 Compact low-rise", `5` = "5 Open mid-rise", `6` = "6 Open low-rise",
             `8` = "8 Large low-rise", `9` = "9 Sparsely built", `10` = "10 Heavy industry", `11` = "A Dense trees",
             `12` = "B Scattered trees", `14` = "D Low plants", `16` = "F Bare soil", `17` = "G Water")
to_df <- function(r, fact = 1) {
  if (fact > 1) r <- aggregate(r, fact, fun = "mean", na.rm = TRUE)
  d <- as.data.frame(mask(r, gar_v), xy = TRUE, na.rm = TRUE); names(d)[3] <- "value"
  rs <- res(r); d$x <- round(d$x / rs[1], 1) * rs[1]; d$y <- round(d$y / rs[2], 1) * rs[2]; d
}
district_mean <- function(r) {
  e <- terra::extract(r, gar_v, fun = mean, na.rm = TRUE)
  tibble(district = gar$admin2Name, value = e[[2]])
}
labs_pts <- st_point_on_surface(gar) |> mutate(X = st_coordinates(geometry)[, 1], Y = st_coordinates(geometry)[, 2])

# =====================================================================
# 1  WHERE IT IS HOT
# =====================================================================
if (has(file.path(GEE, "LST_2020_2023.tif"))) try_fig("F1", {
  lst <- rast(file.path(GEE, "LST_2020_2023.tif"))
  d <- to_df(lst, 4)
  dm <- district_mean(lst) |> arrange(desc(value))
  top <- labs_pts |> inner_join(head(dm, 6), by = c("admin2Name" = "district"))
  p <- ggplot() + geom_raster(data = d, aes(x, y, fill = value)) +
    geom_sf(data = gar, fill = NA, colour = "grey20", linewidth = 0.15) +
    geom_label_repel(data = top, aes(X, Y, label = paste0(admin2Name, "\n", round(value, 1), " °C")),
                     size = 2.5, label.size = 0.1, alpha = 0.9, min.segment.length = 0,
                     nudge_y = -12000, direction = "x", segment.colour = "grey30") +   # labels moved out to sea
    scale_fill_gradientn(colours = heat_pal, name = "LST (°C)") +
    labs(title = "The hottest surfaces sit in Accra's dense core and coastal belt",
         subtitle = sprintf("Mean dry-season land surface temperature, 2020-2023. Labelled: the six hottest districts (region mean %.1f °C)",
                            global(mask(lst, gar_v), "mean", na.rm = TRUE)[1, 1]),
         caption = src_ls) + coord_sf(crs = CRS_UTM) + theme_map
  save_fig(p, "F1_LST_2020_2023_map", 9, 6.5)
})

if (has(file.path(GEE, "LST_2020_2023.tif"), file.path(GEE, "LST_2015_2019.tif"))) try_fig("F2", {
  # Relative change: each pixel's LST minus that period's regional mean, then the difference between periods.
  # This removes shifts common to the whole region (different scenes, weather), leaving where surfaces
  # warmed faster or slower than the region as a whole.
  a <- mask(rast(file.path(GEE, "LST_2015_2019.tif")), gar_v); b <- mask(resample(rast(file.path(GEE, "LST_2020_2023.tif")), a), gar_v)
  ma <- global(a, "mean", na.rm = TRUE)[1, 1]; mb <- global(b, "mean", na.rm = TRUE)[1, 1]
  raw_d <- b - a; rel <- (b - mb) - (a - ma)
  dist_tab <- district_mean(raw_d) |> rename(raw_change_C = value) |>
    left_join(district_mean(rel) |> rename(relative_change_C = value), by = "district") |> arrange(desc(relative_change_C))
  write_csv(dist_tab |> mutate(region_mean_change_C = mb - ma), file.path(TAB, "lst_change_by_district.csv"))
  d <- to_df(rel, 4)
  lim <- max(abs(quantile(d$value, c(0.02, 0.98))))
  p1 <- ggplot() + geom_raster(data = d, aes(x, y, fill = value)) +
    geom_sf(data = gar, fill = NA, colour = "grey20", linewidth = 0.15) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-lim, lim), oob = squish,
                         name = "Relative\nchange (°C)") + coord_sf(crs = CRS_UTM) + theme_map
  p2 <- ggplot(head(dist_tab, 10), aes(reorder(district, relative_change_C), relative_change_C)) + geom_col(fill = "#b2182b") +
    coord_flip() + labs(x = NULL, y = "Warming beyond the regional average (°C)", title = "Districts warming fastest") + theme_story
  p <- (p1 | p2) + plot_layout(widths = c(2.2, 1)) +
    plot_annotation(title = sprintf("The urban core is warming faster than the region: %s and %s gained %.1f and %.1f °C more than the regional average",
                                    dist_tab$district[1], dist_tab$district[2], dist_tab$relative_change_C[1], dist_tab$relative_change_C[2]),
                    subtitle = sprintf("Change in dry-season LST from 2015-2019 to 2020-2023, relative to the regional mean (regional mean change %+.1f °C)", mb - ma),
                    caption = paste(src_ls, "Relative change removes region-wide differences between the two sets of scenes."), theme = theme_story)
  save_fig(p, "F2_LST_change_map", 13, 6)
})

if (has(file.path(GEE, "LCZ_100m.tif"))) try_fig("F3", {
  lcz <- rast(file.path(GEE, "LCZ_100m.tif"))
  d <- as.data.frame(mask(lcz, gar_v), xy = TRUE, na.rm = TRUE); names(d)[3] <- "code"
  d$code <- factor(as.character(d$code), levels = names(lcz_col))
  area <- d |> count(code) |> mutate(km2 = n * 0.01)
  p1 <- ggplot() + geom_raster(data = d, aes(x, y, fill = code)) +
    geom_sf(data = gar, fill = NA, colour = "grey15", linewidth = 0.15) +
    scale_fill_manual(values = lcz_col, labels = lcz_lab, name = "Local Climate Zone", drop = TRUE) +
    coord_sf(crs = CRS_UTM) + theme_map
  p2 <- ggplot(area, aes(reorder(lcz_lab[as.character(code)], km2), km2, fill = code)) + geom_col(show.legend = FALSE) +
    scale_fill_manual(values = lcz_col) + coord_flip() + labs(x = NULL, y = "Area (km²)") + theme_story
  p <- (p1 | p2) + plot_layout(widths = c(2.2, 1)) +
    plot_annotation(title = "Local Climate Zones of the Greater Accra Region",
                    subtitle = "Random Forest classification, 100 m (overall accuracy 0.64, kappa 0.61 on independent validation squares)",
                    caption = "Source: Sentinel-2 (2020-2022 dry season) and Google Open Buildings 2.5D; Random Forest in Google Earth Engine.",
                    theme = theme_story)
  save_fig(p, "F3_LCZ_map", 13, 6.5)
})

# =====================================================================
# 2  WHY IT IS HOT
# =====================================================================
if (has(file.path(GEE, "LCZ_100m.tif"), file.path(GEE, "LST_2020_2023.tif"))) try_fig("F4", {
  lcz <- rast(file.path(GEE, "LCZ_100m.tif"))
  lst100 <- resample(rast(file.path(GEE, "LST_2020_2023.tif")), lcz, method = "average")
  d <- as.data.frame(mask(c(lcz, lst100), gar_v), na.rm = TRUE); names(d) <- c("code", "LST")
  d <- d |> filter(code != 17) |> mutate(code = as.character(code))
  ord <- d |> group_by(code) |> summarise(m = median(LST)) |> arrange(m)
  d$code <- factor(d$code, levels = ord$code)
  trees <- median(d$LST[d$code == "11"]); hottest <- max(ord$m)
  p <- ggplot(d, aes(code, LST, fill = code)) +
    geom_hline(yintercept = trees, linetype = 2, colour = "#006a00") +
    geom_boxplot(outlier.shape = NA, width = 0.6, show.legend = FALSE) +
    annotate("text", x = length(levels(d$code)) + 0.45, y = trees, label = "Dense trees (median)", colour = "#006a00",
             hjust = -0.05, size = 3) +
    scale_fill_manual(values = lcz_col) + scale_x_discrete(labels = lcz_lab) +
    coord_flip(ylim = quantile(d$LST, c(0.01, 0.99))) +
    labs(x = NULL, y = "Land surface temperature (°C), dry season 2020-2023",
         title = sprintf("Dense neighbourhoods run up to %.1f °C hotter than tree-covered land", hottest - trees),
         subtitle = "LST by Local Climate Zone (100 m cells). Differences between zones are significant (Welch's ANOVA, p < 0.001)",
         caption = src_ls) + theme_story
  save_fig(p, "F4_LST_by_LCZ", 9, 6)
})

if (has(file.path(GEE, "LST_2020_2023.tif"), file.path(GEE, "NDVI_2020_2023.tif"))) try_fig("F5", {
  st <- mask(c(rast(file.path(GEE, "NDVI_2020_2023.tif")), rast(file.path(GEE, "LST_2020_2023.tif"))), gar_v)
  d <- as.data.frame(st, na.rm = TRUE); names(d) <- c("NDVI", "LST")
  set.seed(5); d <- d[sample(nrow(d), min(150000, nrow(d))), ] |> filter(NDVI > 0)
  b <- coef(lm(LST ~ NDVI, data = d))[2]
  p <- ggplot(d, aes(NDVI, LST)) + geom_hex(bins = 70) +
    geom_smooth(method = "lm", colour = "black", linewidth = 0.8, se = FALSE) +
    scale_fill_viridis_c(option = "magma", trans = "log10", name = "Pixels") +
    labs(x = "Vegetation (NDVI)", y = "Land surface temperature (°C)",
         title = sprintf("Every 0.1 increase in vegetation cover is linked to about %.1f °C cooler surfaces", -b / 10),
         subtitle = "Pixel-level relationship, dry season 2020-2023 (random sample of 150,000 land pixels)",
         caption = src_ls) + theme_story
  save_fig(p, "F5_vegetation_cools", 8, 6)
})

# =====================================================================
# 3  WHERE HEAT CLUSTERS
# =====================================================================
if (has(file.path(MAP, "hotspots_200m.gpkg"))) try_fig("F6", {
  hs <- st_read(file.path(MAP, "hotspots_200m.gpkg"), quiet = TRUE)
  cols <- c("Cold 99" = "#2166ac", "Cold 95" = "#67a9cf", "Cold 90" = "#d1e5f0", "Not significant" = "grey93",
            "Hot 90" = "#fddbc7", "Hot 95" = "#ef8a62", "Hot 99" = "#b2182b")
  hs$gi_class <- factor(hs$gi_class, levels = names(cols))
  share <- st_drop_geometry(hs) |> count(gi_class) |> mutate(pct = 100 * n / sum(n))
  hot99 <- share$pct[share$gi_class == "Hot 99"]
  p <- ggplot(hs) + geom_sf(aes(fill = gi_class), colour = NA) +
    geom_sf(data = gar, fill = NA, colour = "grey20", linewidth = 0.15) +
    scale_fill_manual(values = cols, name = "Getis-Ord Gi*", drop = FALSE) +
    labs(title = sprintf("Heat is not random: %.0f%% of the region forms statistically significant hot spots (99%%)", hot99),
         subtitle = "Clusters of high land surface temperature, 200 m grid, first- and second-order queen contiguity",
         caption = src_ls) + coord_sf(crs = CRS_UTM) + theme_map
  save_fig(p, "F6_surface_hotspots", 9, 6.5)
})

if (has(file.path(MAP, "hotspots_urban_200m.gpkg"))) try_fig("F6b", {
  hu <- st_read(file.path(MAP, "hotspots_urban_200m.gpkg"), quiet = TRUE)
  cols <- c("Cold 99" = "#2166ac", "Cold 95" = "#67a9cf", "Cold 90" = "#d1e5f0", "Not significant" = "grey93",
            "Hot 90" = "#fddbc7", "Hot 95" = "#ef8a62", "Hot 99" = "#b2182b")
  hu$gi_class <- factor(hu$gi_class, levels = names(cols))
  pct <- 100 * mean(hu$gi_class == "Hot 99")
  top <- if (file.exists(file.path(TAB, "gi_urban_hot99_by_district.csv")))
    head(read_csv(file.path(TAB, "gi_urban_hot99_by_district.csv"), show_col_types = FALSE)$admin2Name, 3) else character(0)
  p <- ggplot(hu) + geom_sf(aes(fill = gi_class), colour = NA) +
    geom_sf(data = gar, fill = NA, colour = "grey20", linewidth = 0.15) +
    scale_fill_manual(values = cols, name = "Getis-Ord Gi*", drop = FALSE) +
    labs(title = sprintf("Within the built-up city, %.0f%% of urban land forms 99%% hot spots%s", pct,
                         if (length(top)) paste0(", led by ", paste(top, collapse = ", ")) else ""),
         subtitle = "Gi* computed among built-up cells only (LCZ 2, 3, 5, 6, 8, 10), 200 m grid",
         caption = src_ls) + coord_sf(crs = CRS_UTM) + theme_map
  save_fig(p, "F6b_urban_hotspots", 9, 6.5)
})

# =====================================================================
# 4  WHO IS AT RISK
# =====================================================================
if (has(file.path(MAP, "risk_cells.rds"), file.path(MAP, "risk_inputs_raw.tif"))) {
  rc <- readRDS(file.path(MAP, "risk_cells.rds"))
  raw <- rast(file.path(MAP, "risk_inputs_raw.tif"))
  cells <- cellFromXY(raw, as.matrix(rc[, c("x", "y")]))
  rc$people <- values(raw[["E"]])[cells] * 0.01          # residents per 100 m cell (density x 0.01 km2)
  rc$vulnerable <- values(raw[["v_count"]])[cells] * 0.01  # residents under 10 or 65+
  rc$people[is.na(rc$people)] <- 0; rc$vulnerable[is.na(rc$vulnerable)] <- 0

  try_fig("F7", {
    chrz_km2 <- sum(rc$CHRZ) * 0.01
    p <- ggplot() + geom_raster(data = rc, aes(x, y, fill = HRL)) +
      geom_raster(data = rc[rc$CHRZ, ], aes(x, y), fill = "#67001f") +
      geom_sf(data = gar, fill = NA, colour = "grey20", linewidth = 0.15) +
      scale_fill_gradient(low = "grey95", high = "#f4a582", name = "Heat risk\n(H × E × V)") +
      labs(title = sprintf("%.0f km² of Greater Accra are Critical Heat Risk Zones", chrz_km2),
           subtitle = "Dark red: 99% hot spots of heat risk (hazard × exposure × vulnerability), 100 m grid",
           caption = "Hazard: Landsat LST; exposure: WorldPop 2020; vulnerability: age (<10, 65+), Relative Wealth Index, vegetation.") +
      coord_sf(crs = CRS_UTM) + theme_map
    save_fig(p, "F7_critical_heat_risk_zones", 9, 6.5)
  })

  try_fig("F8", {
    by_d <- rc |> filter(CHRZ) |> group_by(district) |>
      summarise(km2 = n() * 0.01, people = sum(people), vulnerable = sum(vulnerable), .groups = "drop") |>
      arrange(desc(people)) |> head(10)
    p1 <- ggplot(by_d, aes(reorder(district, km2), km2)) + geom_col(fill = "#b2182b") + coord_flip() +
      labs(x = NULL, y = "Area in CHRZ (km²)") + theme_story
    p2 <- ggplot(by_d, aes(reorder(district, km2), people / 1000)) + geom_col(fill = "#67001f") +
      geom_col(aes(y = vulnerable / 1000), fill = "#fdae61") + coord_flip() +
      labs(x = NULL, y = "Residents in CHRZ (thousands); orange = under 10 or 65+") + theme_story
    p <- (p1 | p2) + plot_annotation(
      title = "Where to act first: the districts with the most people in Critical Heat Risk Zones",
      subtitle = "Top 10 districts by residents living in CHRZ", caption = "Population: WorldPop 2020.", theme = theme_story)
    save_fig(p, "F8_CHRZ_by_district", 12, 5.5)
  })

  try_fig("F9", {
    tot <- sum(rc$people); totv <- sum(rc$vulnerable)
    inz <- sum(rc$people[rc$CHRZ]); inzv <- sum(rc$vulnerable[rc$CHRZ])
    area_share <- 100 * mean(rc$CHRZ)
    d <- tibble(what = factor(c("Land area", "All residents", "Children <10 and adults 65+"),
                              levels = c("Land area", "All residents", "Children <10 and adults 65+")),
                pct = c(area_share, 100 * inz / tot, 100 * inzv / totv))
    p <- ggplot(d, aes(what, pct, fill = what)) + geom_col(width = 0.6, show.legend = FALSE) +
      geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.4, size = 4) +
      scale_fill_manual(values = c("grey60", "#b2182b", "#67001f")) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(x = NULL, y = "Share inside Critical Heat Risk Zones (%)",
           title = sprintf("About %s people, including %s young children and older adults, live in Critical Heat Risk Zones",
                           format(round(inz, -3), big.mark = ","), format(round(inzv, -3), big.mark = ",")),
           subtitle = "Risk is concentrated where people live: the zones hold a far larger share of residents than of land",
           caption = "Population: WorldPop 2020 (100 m). Shares are of the population and land within the study area.") +
      theme_story
    save_fig(p, "F9_people_in_CHRZ", 9, 5.5)
  })

  try_fig("F10", {
    acc <- rc |> filter(CHRZ) |> count(LCZ) |> mutate(pct = 100 * n / sum(n), city = "Greater Accra (this study)",
                                                      code = as.character(LCZ))
    lagos <- tibble(code = c("3", "6", "7", "5", "2/8"), pct = c(59, 21, 17, 2, 1), city = "Lagos (Obe et al., 2023)")
    acc_top <- acc |> select(code, pct, city)
    d <- bind_rows(acc_top, lagos)
    lab2 <- c(lcz_lab, `7` = "7 Lightweight low-rise", `2/8` = "2 and 8")
    col2 <- c(lcz_col, `7` = "#faee05", `2/8` = "#999999")
    p <- ggplot(d, aes(city, pct, fill = code)) + geom_col(width = 0.6, colour = "white") +
      geom_text(aes(label = ifelse(pct >= 4, sprintf("%.0f%%", pct), "")), position = position_stack(vjust = 0.5), size = 3) +
      scale_fill_manual(values = col2, labels = lab2, name = "LCZ") + coord_flip() +
      labs(x = NULL, y = "Share of Critical Heat Risk Zone area (%)",
           title = "In both West African cities, dense low-rise neighbourhoods carry most of the critical heat risk",
           subtitle = "Composition of Critical Heat Risk Zones by Local Climate Zone",
           caption = "Lagos values from Obe et al. (2023), City and Environment Interactions 20, 100128.") + theme_story
    save_fig(p, "F10_CHRZ_Accra_vs_Lagos", 10, 4.5)
  })
}

# =====================================================================
# 5  WHAT IS COMING
# =====================================================================
if (has(file.path(TAB, "future_delta_summary.csv"), file.path(TAB, "future_chrz_summary.csv"))) try_fig("F11", {
  dl <- read_csv(file.path(TAB, "future_delta_summary.csv"), show_col_types = FALSE)
  ch <- read_csv(file.path(TAB, "future_chrz_summary.csv"), show_col_types = FALSE)
  nice <- function(k) sub("ssp245", "SSP2-4.5", sub("ssp585", "SSP5-8.5", sub("_(\\d{4})_(\\d{4})", "\n\\1-\\2", k)))
  dl <- dl |> mutate(scen = ifelse(grepl("585", projection), "SSP5-8.5", "SSP2-4.5"), lab = nice(projection))
  p1 <- ggplot(dl, aes(lab, mean_dT, colour = scen)) +
    geom_pointrange(aes(ymin = mean_dT - mean_model_sd, ymax = mean_dT + mean_model_sd), size = 0.6) +
    scale_colour_manual(values = c("SSP2-4.5" = "#fdae61", "SSP5-8.5" = "#b2182b"), name = NULL) +
    labs(x = NULL, y = "Dry-season warming vs 1981-2010 (°C)", title = "Projected warming") + theme_story
  ch2 <- ch |> filter(!projection %in% c("persistent", "scenario_contingent")) |>
    mutate(lab = ifelse(projection == "baseline", "Baseline\n2020-2023", nice(projection)),
           scen = case_when(grepl("585", projection) ~ "SSP5-8.5", grepl("245", projection) ~ "SSP2-4.5", TRUE ~ "Baseline"))
  base_hrl <- ch2$mean_HRL[ch2$scen == "Baseline"]; base_km2 <- ch2$chrz_km2[ch2$scen == "Baseline"]
  fut <- ch2 |> filter(scen != "Baseline") |> mutate(rise = 100 * (mean_HRL / base_hrl - 1))
  persist <- ch$chrz_km2[ch$projection == "persistent"]
  p2 <- ggplot(fut, aes(lab, rise, fill = scen)) + geom_col(width = 0.65) +
    geom_text(aes(label = sprintf("+%.0f%%", rise)), vjust = -0.4, size = 3.5) +
    scale_fill_manual(values = c("SSP2-4.5" = "#fdae61", "SSP5-8.5" = "#b2182b"), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Rise in mean heat risk vs today (%)", title = "Heat risk intensity") + theme_story
  p <- (p1 | p2) + plot_annotation(
    title = sprintf("Tomorrow's danger zones are today's, only worse: mean heat risk rises by up to %.0f%%, while %.0f km² stay critical under every pathway",
                    max(fut$rise), persist),
    subtitle = sprintf("CMIP6 ensemble (NEX-GDDP-CMIP6, 6 models); bars on the left show ±1 SD across models. Critical Heat Risk Zones: %.0f km² today, %.0f-%.0f km² in the projections. Exposure and vulnerability held at today's levels",
                       base_km2, min(fut$chrz_km2), max(fut$chrz_km2)),
    caption = "Delta-change applied to observed LST; risk re-assessed with Gi* (99%) after joint normalisation of hazard.", theme = theme_story)
  save_fig(p, "F11_future_warming_and_risk", 12, 5.5)
})

proj_files <- file.path(MAP, paste0("LST_projected_", c("ssp245_2041_2060", "ssp245_2081_2100", "ssp585_2041_2060", "ssp585_2081_2100"), ".tif"))
if (has(proj_files)) try_fig("F12", {
  labs4 <- c("SSP2-4.5, 2041-2060", "SSP2-4.5, 2081-2100", "SSP5-8.5, 2041-2060", "SSP5-8.5, 2081-2100")
  d <- bind_rows(lapply(seq_along(proj_files), function(i) to_df(rast(proj_files[i]), 7) |> mutate(panel = labs4[i])))
  d$panel <- factor(d$panel, levels = labs4)
  p <- ggplot() + geom_raster(data = d, aes(x, y, fill = value)) +
    geom_sf(data = gar, fill = NA, colour = "grey25", linewidth = 0.1) + facet_wrap(~panel, ncol = 2) +
    scale_fill_gradientn(colours = heat_pal, name = "LST (°C)") +
    labs(title = "Under high emissions, today's hottest neighbourhoods become hotter still by the end of the century",
         subtitle = "Projected dry-season land surface temperature (observed 2020-2023 LST + CMIP6 ensemble warming)",
         caption = src_ls) + coord_sf(crs = CRS_UTM) + theme_map + theme(axis.text = element_blank())
  save_fig(p, "F12_projected_LST_maps", 11, 9)
})

message("Figures are in ", FIG)
