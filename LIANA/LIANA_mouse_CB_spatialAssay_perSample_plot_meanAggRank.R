library(spatstat.explore)
library(Seurat)
library(liana)
library(dplyr)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(liana)

## LOAD DATA
#load('/mnt/morbo/Data/Users/kwoyshner/cerebellum/data/mouse_msobj_transformed_analysis.rda')
outputDir <- 'liana/'
prefix <- 'mouseSpatial_PCGC_perSample_'
setwd(outputDir)
# Read in to avoid regenerating liana results
alllianas_agg <- read.csv(paste0(outputDir, prefix, "alllianasagg.csv"), row.names = 1)

alllianas_agg <- alllianas_agg %>%
  mutate(
    condition = sub(":.*$", "", uniquepanel),
    sample_id = sub("^.*:", "", uniquepanel),
    interaction_id = paste(source, target, ligand.complex, receptor.complex, sep = "|"),
    lr_pair_id = paste(ligand.complex, receptor.complex, sep = "|")
  )

cond_order <- c("sal-Nx", "sal-Hx", "LPS-Nx", "DH")
cond_order <- cond_order[cond_order %in% unique(alllianas_agg$condition)]


# Shared plot theme
liana_theme <- theme(
  axis.text.x      = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
  axis.text.y      = element_text(size = 12),
  axis.title       = element_text(size = 12),
  axis.title.x.top = element_text(size = 12),
  strip.text.x     = element_text(size = 12),
  plot.title       = element_text(size = 14, face = "bold", hjust = 0.5)
)



#####################################################
######## FUNCTION TO PREPARE TOP INTERACTIONS ########
#####################################################

prepare_liana_condition_plot_data <- function(
    df,
    condition_name,
    selection_method = c("top20_total", "top10_each_direction"),
    n_total = 20,
    n_each_direction = 10,
    required_n_samples = 2
) {
  
  selection_method <- match.arg(selection_method)
  
  # Filter to one condition and PC/GC only
  cond_df <- df %>%
    filter(
      condition == condition_name
    )
  
  # Count samples per unique directional interaction.
  # Unique match = source-target-ligand-receptor.
  interaction_sample_presence <- cond_df %>%
    distinct(source, target, ligand.complex, receptor.complex, interaction_id, sample_id) %>%
    count(source, target, ligand.complex, receptor.complex, interaction_id, name = "n_samples_present")
  
  # Keep only source-target-ligand-receptor interactions PRESENT in both samples
  present_both <- interaction_sample_presence %>%
    filter(n_samples_present >= required_n_samples)
  
  cond_present_both <- cond_df %>%
    semi_join(
      present_both,
      by = c("source", "target", "ligand.complex", "receptor.complex", "interaction_id")
    )
  
  # Rank each directional interaction min rank per sample.
  # Lower aggregate_rank = better.
  interaction_ranked <- cond_present_both %>%
    group_by(source, target, ligand.complex, receptor.complex, interaction_id, lr_pair_id) %>%
    summarise(
      mean_aggregate_rank = mean(aggregate_rank, na.rm = TRUE),
      mean_sca_LRscore = mean(sca.LRscore, na.rm = TRUE),
      min_aggregate_rank = min(aggregate_rank, na.rm = TRUE),
      n_samples_present = n_distinct(sample_id),
      .groups = "drop"
    ) %>%
    arrange(mean_aggregate_rank, desc(mean_sca_LRscore))

  #####################################################
  ######### SELECT TOP INTERACTIONS ###################
  #####################################################
  
  if (selection_method == "top20_total") {
    
    # Top 20 total directional source-target-ligand-receptor interactions
    selected_interactions <- interaction_ranked %>%
      arrange(mean_aggregate_rank, desc(mean_sca_LRscore)) %>%
      slice_head(n = n_total)
    
  } else if (selection_method == "top10_each_direction") {
    
    # Top 10 PC->GC and top 10 GC->PC
    selected_interactions <- interaction_ranked %>%
      filter(
        (source == "PC" & target == "GC") |
          (source == "GC" & target == "PC")
      ) %>%
      group_by(source, target) %>%
      arrange(mean_aggregate_rank, desc(mean_sca_LRscore), .by_group = TRUE) %>%
      slice_head(n = n_each_direction) %>%
      ungroup()
  }
  
  #####################################################
  ######### SUBSET FULL CONDITION DATA ################
  #####################################################
  
  # once top pairs are chosen, subset the whole dataset to every row that has this ligand-receptor pair.
  #
  # so if Tnc-Ptprz1 is selected in PC->GC, you keep all rows in that
  # condition with ligand.complex == Tnc and receptor.complex == Ptprz1,
  # including PC->PC, GC->GC, etc., if present.
  #
  
  selected_lr_pairs <- selected_interactions %>%
    distinct(ligand.complex, receptor.complex, lr_pair_id)
  
  cond_selected_full <- cond_df %>%
    semi_join(
      selected_lr_pairs,
      by = c("ligand.complex", "receptor.complex", "lr_pair_id")
    )
  
  #####################################################
  ######### AVERAGE ACROSS SAMPLES ####################
  #####################################################
  
  # LIANA dotplot defaults use:
  # color  = natmi.edge_specificity
  # size   = sca.LRscore
  #
  # We average those across samples for each condition/source/target/L/R.
  
  cond_avg <- cond_selected_full %>%
    group_by(
      condition,
      source,
      target,
      ligand.complex,
      receptor.complex,
      lr_pair_id
    ) %>%
    summarise(
      aggregate_rank = mean(aggregate_rank, na.rm = TRUE),
      mean_rank = mean(mean_rank, na.rm = TRUE),
      natmi.edge_specificity = mean(natmi.edge_specificity, na.rm = TRUE),
      sca.LRscore = mean(sca.LRscore, na.rm = TRUE),
      n_samples = n_distinct(sample_id),
      .groups = "drop"
    )
  
  # Alphabetical plotting order by ligand, then receptor
  cond_avg <- cond_avg %>%
    arrange(ligand.complex, receptor.complex) %>%
    mutate(
      ligand.complex = factor(ligand.complex, levels = unique(ligand.complex)),
      receptor.complex = factor(receptor.complex, levels = unique(receptor.complex))
    )
  
  return(list(
    plot_data = cond_avg,
    selected_interactions = selected_interactions,
    selected_lr_pairs = selected_lr_pairs
  ))
}


#####################################################
######## FUNCTION TO MAKE LIANA DOTPLOT #############
#####################################################

make_liana_plot <- function(plot_data, condition_name, title_suffix = NULL,
                            spec_limits = NULL, expr_limits = NULL) {
  
  title <- condition_name
  if (!is.null(title_suffix)) {
    title <- paste0(condition_name, " - ", title_suffix)
  }
  
  p <- plot_data %>%
    liana_dotplot(
      source_groups = c("PC", "GC"),
      target_groups = c("PC", "GC")
    ) +
    #scale_size_continuous(
    #  limits = spec_limits,
    #  #range = c(0.25, 2.5),
    #  name = "Interaction\nSpecificity"
    #) +
    scale_color_viridis_c(
      limits = expr_limits,
      name = "Expression\nMagnitude"
    ) +
    liana_theme +
    ggtitle(title)
  
  return(p)
}

#####################################################
############# METHOD A: TOP 20 TOTAL ################
#####################################################

plots_top20_total <- list()
plotdata_top20_total <- list()
selected_top20_total <- list()

#####################################################
### First pass: prepare data, but do not plot yet ####
#####################################################

for (cond in cond_order) {
  
  message("Preparing condition: ", cond, " using top20_any")
  
  res <- prepare_liana_condition_plot_data(
    df = alllianas_agg,
    condition_name = cond,
    selection_method = "top20_total",
    n_total = 20,
    required_n_samples = 2
  )
  
  plotdata_top20_total[[cond]] <- res$plot_data
  selected_top20_total[[cond]] <- res$selected_interactions
}

#####################################################
### Calculate shared limits across all conditions ####
#####################################################

plot_df_top20_total_all <- bind_rows(plotdata_top20_total)

spec_limits_top20_total <- range(
  plot_df_top20_total_all$natmi.edge_specificity,
  na.rm = TRUE
)

expr_limits_top20_total <- range(
  plot_df_top20_total_all$sca.LRscore,
  na.rm = TRUE
)

print(spec_limits_top20_total)
print(expr_limits_top20_total)

#####################################################
### Second pass: make plots using shared limits ######
#####################################################

for (cond in cond_order) {
  
  message("Plotting condition: ", cond, " using top20_any")
  
  p <- make_liana_plot(
    plot_data = plotdata_top20_total[[cond]],
    condition_name = cond,
    title_suffix = "Top 20 total",
    spec_limits = spec_limits_top20_total,
    expr_limits = expr_limits_top20_total
  )
  
  plots_top20_total[[cond]] <- p
  
  ggsave(
    filename = paste0(outputDir, prefix, "liana_dotplot_", cond, "_top20_anyPCGC_presentBothSamples_meanRank_avgMagSpecSamples.pdf"),
    plot = p,
    width = 14,
    height = 10,
    units = "in"
  )
}

#####################################################
### Save selected interactions ######################
#####################################################

selected_top20_total_df <- bind_rows(selected_top20_total, .id = "condition")

write.csv(
  selected_top20_total_df,
  paste0(outputDir, prefix, "selected_interactions_top20_anyPCGC_presentBothSamples_meanRank_avgMagSpecSamples.csv"),
  row.names = FALSE
)

#####################################################
### Combined plot, method A #########################
#####################################################

combined_top20_total <- wrap_plots(plots_top20_total[cond_order], ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  filename = paste0(outputDir, prefix, "liana_dotplot_combined_top20_anyPCGC_presentBothSamples_meanRank_avgMagSpecSamples.pdf"),
  plot = combined_top20_total,
  width = 18,
  height = 14,
  units = "in"
)


#####################################################
######## METHOD B: TOP 10 PC->GC AND GC->PC #########
#####################################################

plots_top10_direction <- list()
plotdata_top10_direction <- list()
selected_top10_direction <- list()

#####################################################
### First pass: prepare data, but do not plot yet ####
#####################################################

for (cond in cond_order) {
  
  message("Preparing condition: ", cond, " using top10_each_direction")
  
  res <- prepare_liana_condition_plot_data(
    df = alllianas_agg,
    condition_name = cond,
    selection_method = "top10_each_direction",
    n_each_direction = 10,
    required_n_samples = 2
  )
  
  plotdata_top10_direction[[cond]] <- res$plot_data
  selected_top10_direction[[cond]] <- res$selected_interactions
}

#####################################################
### Calculate shared limits across all conditions ####
#####################################################

plot_df_top10_direction_all <- bind_rows(plotdata_top10_direction)

spec_limits_top10_direction <- range(
  plot_df_top10_direction_all$natmi.edge_specificity,
  na.rm = TRUE
)

expr_limits_top10_direction <- range(
  plot_df_top10_direction_all$sca.LRscore,
  na.rm = TRUE
)

print(spec_limits_top10_direction)
print(expr_limits_top10_direction)

#####################################################
### Second pass: make plots using shared limits ######
#####################################################

for (cond in cond_order) {
  
  message("Plotting condition: ", cond, " using top10_each_direction")
  
  p <- make_liana_plot(
    plot_data = plotdata_top10_direction[[cond]],
    condition_name = cond,
    title_suffix = "Top 10 PC→GC and GC→PC",
    spec_limits = spec_limits_top10_direction,
    expr_limits = expr_limits_top10_direction
  )
  
  plots_top10_direction[[cond]] <- p
  
  ggsave(
    filename = paste0(outputDir, prefix, "liana_dotplot_", cond, "_top10_eachDirection_presentBothSamples_meanRank_avgMagSpecSamples.pdf"),
    plot = p,
    width = 14,
    height = 10,
    units = "in"
  )
}

#####################################################
### Save selected interactions ######################
#####################################################

selected_top10_direction_df <- bind_rows(selected_top10_direction, .id = "condition")

write.csv(
  selected_top10_direction_df,
  paste0(outputDir, prefix, "selected_interactions_top10_eachDirection_presentBothSamples_meanRank_avgMagSpecSamples.csv"),
  row.names = FALSE
)

#####################################################
### Combined plot, method B #########################
#####################################################

combined_top10_direction <- wrap_plots(plots_top10_direction[cond_order], ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  filename = paste0(outputDir, prefix, "liana_dotplot_combined_top10_eachDirection_presentBothSamples_meanRank_avgMagSpecSamples.pdf"),
  plot = combined_top10_direction,
  width = 18,
  height = 14,
  units = "in"
)