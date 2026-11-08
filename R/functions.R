# =====================================================================
# functions.R  Shared settings and helper functions (sourced by every script)
# =====================================================================
suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(tidyr)
  library(readr); library(ggplot2); library(spdep)
})
CRS_UTM <- "EPSG:32630"
GEE <- "data/gee"; RAW <- "data/raw"; TAB <- "outputs/tables"; MAP <- "outputs/maps"

# Study area: the 29 districts, projected to UTM 30N
gar <- st_read(file.path(RAW, "GAR_districts.shp"), quiet = TRUE) |> st_transform(CRS_UTM)
gar_v <- vect(gar)

# Min-max normalisation of a raster (Equation 3.4); optional fixed limits for joint normalisation
minmax <- function(r, lo = NULL, hi = NULL) {
  if (is.null(lo)) { mm <- global(r, c("min", "max"), na.rm = TRUE); lo <- mm[1, 1]; hi <- mm[1, 2] }
  (r - lo) / (hi - lo)
}

# Empty template grid over the study area at a given resolution (metres)
template <- function(res) {
  rast(ext(gar_v), resolution = res, crs = CRS_UTM)
}

# Raster cells to sf polygons (one polygon per non-NA cell)
cells_to_sf <- function(r) st_as_sf(as.polygons(r, aggregate = FALSE, na.rm = TRUE))

# Getis-Ord Gi* classes at 90, 95 and 99% confidence
gi_class <- function(z) {
  cut(z, breaks = c(-Inf, -2.576, -1.96, -1.645, 1.645, 1.96, 2.576, Inf),
      labels = c("Cold 99", "Cold 95", "Cold 90", "Not significant", "Hot 90", "Hot 95", "Hot 99"),
      right = FALSE)
}

# Gi* z-scores for a vector x with a neighbour list (self is added for the star statistic)
gistar <- function(x, nb, glist = NULL, style = "W") {
  nbs <- include.self(nb)
  if (!is.null(glist)) {
    lw <- nb2listw(nbs, glist = glist, style = style, zero.policy = TRUE)
  } else {
    lw <- nb2listw(nbs, style = style, zero.policy = TRUE)
  }
  as.numeric(localG(x, lw, zero.policy = TRUE))
}

# Random sample of raster cells per class, to keep tests from being dominated by
# millions of spatially dependent pixels
sample_by_class <- function(df, class_col, n = 1000, seed = 42) {
  set.seed(seed)
  df |> group_by(.data[[class_col]]) |> slice_sample(n = n) |> ungroup()
}

# Gi* and LISA on LST aggregated to a grid of 'res' metres (Section 3.9.1)
run_hotspots <- function(lst30, res, alternatives = TRUE) {
  r <- mask(resample(lst30, template(res), method = "average"), gar_v)
  names(r) <- "LST"
  g <- cells_to_sf(r)
  x <- g$LST
  pts <- st_coordinates(st_centroid(st_geometry(g)))

  # Weights: main = first- and second-order queen; alternatives = first-order queen, 400 m, 800 m
  q1  <- poly2nb(g, queen = TRUE)
  q12 <- nblag_cumul(nblag(q1, 2))
  if (alternatives) { d400 <- dnearneigh(pts, 0, 400); d800 <- dnearneigh(pts, 0, 800) }

  g$gi_main <- gistar(x, q12)
  g$gi_class <- gi_class(g$gi_main)
  if (!alternatives) return(g)
  g$gi_q1   <- gistar(x, q1)
  g$gi_d400 <- gistar(x, d400)
  g$gi_d800 <- gistar(x, d800)

  # Robust hot and cold spots: significant (95%) under at least two of the three alternatives
  hot  <- (g$gi_q1 >= 1.96) + (g$gi_d400 >= 1.96) + (g$gi_d800 >= 1.96)
  cold <- (g$gi_q1 <= -1.96) + (g$gi_d400 <= -1.96) + (g$gi_d800 <= -1.96)
  g$robust <- case_when(hot >= 2 ~ "Robust hot spot", cold >= 2 ~ "Robust cold spot", TRUE ~ "Not robust")

  # Anselin Local Moran's I (main weights)
  lw  <- nb2listw(q12, style = "W", zero.policy = TRUE)
  li  <- localmoran(x, lw, zero.policy = TRUE)
  z   <- as.numeric(scale(x)); lag <- lag.listw(lw, z, zero.policy = TRUE)
  p   <- li[, "Pr(z != E(Ii))"]
  g$lisa <- case_when(p >= 0.05 ~ "Not significant",
                      z > 0 & lag > 0 ~ "High-High", z < 0 & lag < 0 ~ "Low-Low",
                      z > 0 & lag < 0 ~ "High-Low outlier", TRUE ~ "Low-High outlier")
  g
}


# Inverse-distance spatial weights with the cell itself included (Gi*), for 100 m risk grids
idw_listw <- function(pts, band) {
  nb <- include.self(dnearneigh(pts, 0, band))
  d  <- nbdists(nb, pts)
  nb2listw(nb, glist = lapply(d, function(x) 1 / pmax(x, 100)), style = "B", zero.policy = TRUE)
}
gi_z <- function(x, lw) as.numeric(localG(x, lw, zero.policy = TRUE))

# ---------------------------------------------------------------------
# Getis-Ord Gi* on a regular raster grid with inverse-distance weights (self included),
# computed with moving-window sums. Same statistic as localG() with the IDW listw above,
# but uses a fraction of the memory and runs in seconds on hundreds of thousands of cells.
#   r: SpatRaster (NA outside the analysis area); band: distance band in metres
# ---------------------------------------------------------------------
gi_focal <- function(r, band, res = 100) {
  k <- floor(band / res)
  off <- expand.grid(j = -k:k, i = -k:k)
  d <- res * sqrt(off$i^2 + off$j^2)
  W <- matrix(ifelse(d <= band + 1e-6, 1 / pmax(d, res), 0), nrow = 2 * k + 1, byrow = TRUE)
  valid <- !is.na(r)
  x0 <- ifel(valid, r, 0); v1 <- ifel(valid, 1, 0)
  SW  <- focal(v1, W,   fun = "sum", fillvalue = 0)
  SWX <- focal(x0, W,   fun = "sum", fillvalue = 0)
  SW2 <- focal(v1, W^2, fun = "sum", fillvalue = 0)
  vals <- values(r, na.rm = TRUE); n <- length(vals)
  xbar <- mean(vals); S <- sqrt(mean(vals^2) - xbar^2)
  G <- (SWX - xbar * SW) / (S * sqrt((n * SW2 - SW^2) / (n - 1)))
  mask(G, r)
}
# Gi* for a vector of cell values: v[i] belongs to grid cell cells[i] of template tmpl
gi_cells <- function(v, cells, tmpl, band) {
  r <- rast(tmpl); r[cells] <- v
  values(gi_focal(r, band, res = res(tmpl)[1]))[cells]
}
