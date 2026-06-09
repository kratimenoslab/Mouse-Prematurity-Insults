suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(here)
})
set.seed(42)

# ---- Project layout ---------------------------------------------------------
project_dir <- "."
data_dir    <- file.path(project_dir, "data")
out_dir     <- file.path(project_dir, "derived")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

group_levels <- c("Sal-N","Sal-Hx","LPS-N","LPS-Hx")

# Normalize any synonym for the four canonical experimental groups.
normalize_group <- function(g) {
  case_when(
    g %in% c("Sal-Nx","NX","Nx","Control")           ~ "Sal-N",
    g %in% c("Sal-Hx","HX","Hx","Hypoxia")           ~ "Sal-Hx",
    g %in% c("LPS-Nx","LPS")                         ~ "LPS-N",
    g %in% c("LPS-Hx","DH","Double Hit","Double-Hit") ~ "LPS-Hx",
    TRUE                                             ~ g
  )
}

# ---- Primary measurement: total puncta density per animal -------------------
# One value per cerebellum (per-animal Total VGLUT2 puncta / Total ML area).
total_density <- read_excel(file.path(data_dir, "total_density_per_animal.xlsx"),
                            .name_repair = "minimal") |>
  pivot_longer(everything(), names_to = "group_raw",
               values_to = "total_density_per_um2") |>
  filter(!is.na(total_density_per_um2)) |>
  mutate(group_label = factor(normalize_group(group_raw), levels = group_levels),
         animal_id   = paste0(group_raw, "_",
                              ave(group_raw, group_raw, FUN = seq_along))) |>
  select(animal_id, group_label, total_density_per_um2)

write_csv(total_density, file.path(out_dir, "per_animal_total_density.csv"))

# ---- Secondary measurement: per-region (Inner / Outer ML) puncta data -------
# Each per-animal summary spreadsheet carries the spot counts and areas for
# the Inner (blue) and Outer (yellow) Molecular Layer regions.
summary_files <- list.files(file.path(data_dir, "per_animal_summaries"),
                            pattern = "_summary\\.xlsx$",
                            recursive = TRUE, full.names = TRUE)

read_summary <- function(p) {
  d <- suppressWarnings(read_excel(p, .name_repair = "minimal"))
  v <- setNames(d$value, d$metric)
  tibble(
    animal_id    = sub("_summary\\.xlsx$", "", basename(p)),
    group_folder = basename(dirname(p)),
    outer_count  = as.numeric(v["yellow_total"]),     # yellow = Outer ML
    inner_count  = as.numeric(v["blue_total"]),       # blue   = Inner ML
    outer_area   = as.numeric(v["yellow_area_um2"]),
    inner_area   = as.numeric(v["blue_area_um2"])
  )
}

per_band <- bind_rows(lapply(summary_files, read_summary)) |>
  mutate(
    group_label              = factor(normalize_group(group_folder),
                                      levels = group_levels),
    outer_density_per_um2    = outer_count / outer_area,
    inner_density_per_um2    = inner_count / inner_area,
    outer_inner_ratio        = outer_density_per_um2 / inner_density_per_um2,
    total_punctae            = outer_count + inner_count,
    total_ml_area_um2        = outer_area  + inner_area,
    total_density_per_um2    = total_punctae / total_ml_area_um2
  ) |>
  select(animal_id, group_label,
         inner_count, outer_count, inner_area, outer_area,
         inner_density_per_um2, outer_density_per_um2, outer_inner_ratio,
         total_punctae, total_ml_area_um2, total_density_per_um2)

write_csv(per_band, file.path(out_dir, "per_animal_per_band.csv"))

cat("== Primary measurement (per-animal total VGLUT2 density) ==\n")
print(total_density |> group_by(group_label) |>
  summarise(n    = n(),
            mean = signif(mean(total_density_per_um2), 3),
            sem  = signif(sd(total_density_per_um2) / sqrt(n()), 3),
            .groups = "drop"))

cat("\n== Secondary measurement (per-band Inner / Outer ML) ==\n")
print(per_band |> group_by(group_label) |>
  summarise(n = n(),
            inner_density     = signif(mean(inner_density_per_um2), 3),
            outer_density     = signif(mean(outer_density_per_um2), 3),
            outer_inner_ratio = signif(mean(outer_inner_ratio),     3),
            .groups = "drop"))

cat("\nWrote:\n  ", file.path(out_dir, "per_animal_total_density.csv"),
    "\n  ", file.path(out_dir, "per_animal_per_band.csv"), "\n")
