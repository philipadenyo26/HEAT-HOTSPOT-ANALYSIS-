# =====================================================================
# 04_hotspots.R  Getis-Ord Gi* and Anselin LISA on surface heat
# Chapter 3, Section 3.9.1 (Objective 2). Uses functions.R helpers.
# Takes a few minutes: the 200 m grid has roughly 90,000 cells.
# =====================================================================
source("R/functions.R")
lst30 <- rast(file.path(GEE, "LST_2020_2023.tif"))

# run_hotspots() is defined in functions.R

g200 <- run_hotspots(lst30, 200)
g400 <- run_hotspots(lst30, 400)
st_write(g200, file.path(MAP, "hotspots_200m.gpkg"), delete_dsn = TRUE, quiet = TRUE)
st_write(g400, file.path(MAP, "hotspots_400m.gpkg"), delete_dsn = TRUE, quiet = TRUE)

# Hotspots within the built-up city only (built LCZ classes 2, 3, 5, 6, 8, 10): where heat clusters
# among neighbourhoods, separate from the hot dry plains that dominate the region-wide map
if (file.exists(file.path(GEE, "LCZ_100m.tif"))) {
  l30 <- resample(rast(file.path(GEE, "LCZ_100m.tif")), lst30, method = "near")
  lst_urb <- mask(lst30, l30 %in% c(2, 3, 5, 6, 8, 10), maskvalues = c(0, NA))
  gu <- run_hotspots(lst_urb, 200, alternatives = FALSE)
  st_write(gu, file.path(MAP, "hotspots_urban_200m.gpkg"), delete_dsn = TRUE, quiet = TRUE)
  write_csv(st_drop_geometry(gu) |> count(gi_class) |> mutate(area_km2 = n * 0.04, share = n / sum(n)),
            file.path(TAB, "gi_classes_urban_200m.csv"))
  hot <- st_join(st_centroid(gu[gu$gi_main >= 2.576, ]), gar["admin2Name"])
  write_csv(st_drop_geometry(hot) |> count(admin2Name, name = "cells") |> mutate(area_km2 = cells * 0.04) |>
              arrange(desc(cells)), file.path(TAB, "gi_urban_hot99_by_district.csv"))
}

# Summary tables (area in km2: 0.04 per 200 m cell)
gi_tab <- st_drop_geometry(g200) |> count(gi_class) |> mutate(area_km2 = n * 0.04, share = n / sum(n))
lisa_tab <- st_drop_geometry(g200) |> count(lisa) |> mutate(area_km2 = n * 0.04)
rob_tab <- st_drop_geometry(g200) |> count(robust) |> mutate(area_km2 = n * 0.04)
write_csv(gi_tab, file.path(TAB, "gi_classes_200m.csv"))
write_csv(lisa_tab, file.path(TAB, "lisa_classes_200m.csv"))
write_csv(rob_tab, file.path(TAB, "gi_robustness_200m.csv"))

# Resolution check: share of robust 200 m hot spots that fall in 400 m hot spots (95%)
hot200 <- st_centroid(g200[g200$robust == "Robust hot spot", ])
hot400 <- g400[g400$gi_main >= 1.96, ]
inside <- lengths(st_intersects(hot200, hot400)) > 0
write_csv(tibble(robust_hot_200m = nrow(hot200), share_in_400m_hot = mean(inside)),
          file.path(TAB, "gi_resolution_check.csv"))

# Maps
cols <- c("Cold 99" = "#2166ac", "Cold 95" = "#67a9cf", "Cold 90" = "#d1e5f0", "Not significant" = "grey92",
          "Hot 90" = "#fddbc7", "Hot 95" = "#ef8a62", "Hot 99" = "#b2182b")
m1 <- ggplot(g200) + geom_sf(aes(fill = gi_class), colour = NA) + scale_fill_manual(values = cols, name = "Gi*") +
  geom_sf(data = gar, fill = NA, colour = "grey30", linewidth = 0.2) + theme_void()
m2 <- ggplot(g200) + geom_sf(aes(fill = lisa), colour = NA) +
  scale_fill_manual(values = c("High-High" = "#b2182b", "Low-Low" = "#2166ac", "High-Low outlier" = "#f4a582",
                               "Low-High outlier" = "#92c5de", "Not significant" = "grey92"), name = "LISA") +
  geom_sf(data = gar, fill = NA, colour = "grey30", linewidth = 0.2) + theme_void()
ggsave(file.path(MAP, "map_gi_star_200m.png"), m1, width = 8, height = 5, dpi = 400)
ggsave(file.path(MAP, "map_lisa_200m.png"), m2, width = 8, height = 5, dpi = 400)
