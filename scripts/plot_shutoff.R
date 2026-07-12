#!/usr/bin/env Rscript
# Plot host vs phage read-fraction timecourse from results/shutoff_timecourse.tsv
#
# Usage: Rscript scripts/plot_shutoff.R results/shutoff_timecourse.tsv results/plots/

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_shutoff.R <shutoff_timecourse.tsv> <output_plot_dir>")
}
in_path <- args[1]
out_dir <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

df <- read_tsv(in_path, show_col_types = FALSE)

if (nrow(df) == 0) {
  stop("Input table is empty - nothing to plot.")
}

# --- Averaged stacked-area plot across replicates ---
avg_df <- df %>%
  group_by(timepoint_min) %>%
  summarise(
    host_fraction = mean(host_fraction, na.rm = TRUE),
    phage_fraction = mean(phage_fraction, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(host_fraction, phage_fraction),
               names_to = "origin", values_to = "fraction") %>%
  mutate(origin = recode(origin, host_fraction = "Host", phage_fraction = "Phage"))

p_area <- ggplot(avg_df, aes(x = timepoint_min, y = fraction, fill = origin)) +
  geom_area(position = "stack", alpha = 0.9) +
  scale_fill_manual(values = c("Host" = "#3B6FA0", "Phage" = "#D9622B")) +
  scale_x_continuous(breaks = sort(unique(avg_df$timepoint_min))) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Host transcriptional shutoff over phage infection",
    subtitle = "Mean read fraction by origin, averaged across replicates",
    x = "Minutes post-infection",
    y = "Fraction of aligned reads",
    fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "shutoff_stacked_area.png"), p_area, width = 8, height = 5, dpi = 300)

# --- Per-replicate line plot, to check reproducibility before trusting the average ---
rep_df <- df %>%
  mutate(replicate = factor(replicate)) %>%
  select(timepoint_min, replicate, host_fraction, phage_fraction) %>%
  pivot_longer(cols = c(host_fraction, phage_fraction),
               names_to = "origin", values_to = "fraction") %>%
  mutate(origin = recode(origin, host_fraction = "Host", phage_fraction = "Phage"))

p_lines <- ggplot(rep_df, aes(x = timepoint_min, y = fraction, color = origin,
                               linetype = replicate, group = interaction(origin, replicate))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Host" = "#3B6FA0", "Phage" = "#D9622B")) +
  scale_x_continuous(breaks = sort(unique(rep_df$timepoint_min))) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Host vs phage read fraction, by replicate",
    x = "Minutes post-infection",
    y = "Fraction of aligned reads",
    color = NULL,
    linetype = "Replicate"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "shutoff_by_replicate.png"), p_lines, width = 8, height = 5, dpi = 300)

message("Wrote plots to ", out_dir)
