# =====================================================================
# 00_setup.R  Run once. Installs packages and creates folders.
# Open RStudio with the thesis_analysis folder as the working directory.
# =====================================================================
pkgs <- c("terra", "sf", "dplyr", "tidyr", "readr", "ggplot2", "skimr",
          "spdep", "GWmodel", "gstat", "car", "scales", "patchwork", "ggrepel", "hexbin")
new <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new)) install.packages(new)
for (d in c("data/raw", "data/gee", "outputs/tables", "outputs/maps")) dir.create(d, recursive = TRUE, showWarnings = FALSE)
message("Setup complete")
