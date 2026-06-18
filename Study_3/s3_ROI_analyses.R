  # BETA ANALYSES
# ============================================================
library(tidyverse)
library(rstatix)
library(ggpubr)
library(openxlsx)
library(ggplot2)
library(cowplot)
library(png)

# ============================================================
# FILE PATHS
# ============================================================
roi_files <- list(
  "lIFG [-44, 22, -8]"   = "Betas/Betas_leftIFG.txt",
  "rIFG [56, 24, 18]"    = "Betas/Betas_rightIFG.txt",
  "lMTG [-62, -34, -6]"  = "Betas/Betas_leftTemp.txt",
  "rMTG [62, -30, -6]"   = "Betas/Betas_rightTemp.txt",
  "lTPJ [-44, -56, 40]"  = "Betas/Betas_leftTPJ.txt",
  "rTPJ [50, -56, 30]"   = "Betas/Betas_rightTPJ.txt",
  "lmPFC [-4, 42, 28]"   = "Betas/Betas_leftmPFC.txt",
  "rmPFC [6, 42, 36]"    = "Betas/Betas_rightmPFC.txt",
  "lPCun [-4, -66, 32]"  = "Betas/Betas_leftPCun.txt"
)

roi_order_shared <- c(
  "lIFG [-44, 22, -8]",
  "rIFG [56, 24, 18]",
  "lMTG [-62, -34, -6]",
  "rMTG [62, -30, -6]",
  "lTPJ [-44, -56, 40]",
  "rTPJ [50, -56, 30]",
  "lmPFC [-4, 42, 28]",
  "rmPFC [6, 42, 36]",
  "lPCun [-4, -66, 32]"
)

task_levels       <- c("prosody", "semantic", "irony", "sarcasm", "tom")
task_order_shared <- task_levels

task_labels <- setNames(
  c("Prosody", "Semantic", "Irony", "Sarcasm", "ToM"),
  task_levels
)

results_dir <- "plots/s3_betas"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

exclude_ids <- c(11, 13, 15, 17, 22)

# ============================================================
# COLORBLIND-SAFE PALETTE  (Okabe–Ito)
# ============================================================
cb_palette <- c("positive" = "#0072B2",
                "negative" = "#E69F00")

# ============================================================
# BUILD 4-CONDITION DATASET
# ============================================================
process_4cond <- function(file, roi_name) {
  df <- read.table(file, header = TRUE, sep = "\t") %>%
    mutate(id = row_number()) %>%
    filter(!id %in% exclude_ids)
  
  map_dfr(task_levels, function(t) {
    tibble(
      id    = df$id,
      ROI   = roi_name,
      task  = t,
      CP_SP = rowMeans(select(df, contains(paste0(t, ".CP_SPpos"))), na.rm = TRUE),
      CP_SN = rowMeans(select(df, contains(paste0(t, ".CP_SNpos"))), na.rm = TRUE),
      CN_SP = rowMeans(select(df, contains(paste0(t, ".CN_SPneg"))), na.rm = TRUE),
      CN_SN = rowMeans(select(df, contains(paste0(t, ".CN_SNneg"))), na.rm = TRUE)
    )
  })
}

all_4cond <- imap_dfr(roi_files, process_4cond) %>%
  pivot_longer(cols = c(CP_SP, CP_SN, CN_SP, CN_SN),
               names_to = "condition", values_to = "beta") %>%
  mutate(
    prosody_val  = ifelse(str_starts(condition, "CP"), "positive", "negative"),
    semantic_val = ifelse(str_ends(condition,   "SP"), "positive", "negative"),
    task         = factor(task, levels = task_levels)
  )

# ============================================================
# 2×2 REPEATED-MEASURES ANOVA
# ============================================================
anova_2x2 <- all_4cond %>%
  group_by(ROI, task) %>%
  anova_test(
    dv     = beta,
    wid    = id,
    within = c(prosody_val, semantic_val)
  )

anova_table <- get_anova_table(anova_2x2) %>%
  as_tibble() %>%
  ungroup()

# ============================================================
# COMPUTE FDR CORRECTION — per effect family
# Applied once here, used everywhere below
# ============================================================
anova_table <- anova_table %>%
  mutate(p_num = as.numeric(gsub("[^0-9eE.+-]", "", p))) %>%
  group_by(Effect) %>%
  mutate(p_fdr = p.adjust(p_num, method = "fdr")) %>%
  ungroup()

# ============================================================
# EXPORT FULL ANOVA RESULTS TO EXCEL
# ============================================================
anova_clean <- anova_table %>%
  mutate(
    ROI  = factor(ROI,  levels = roi_order_shared),
    task = factor(task, levels = task_order_shared),
    Effect = case_when(
      Effect == "prosody_val"              ~ "Prosody",
      Effect == "semantic_val"             ~ "Semantics",
      Effect == "prosody_val:semantic_val" ~ "Prosody × Semantics",
      TRUE ~ Effect
    ),
    Effect = factor(Effect, levels = c("Prosody", "Semantics", "Prosody × Semantics")),
    F      = round(F, 3),
    ges    = round(as.numeric(ges), 3),
    p_fdr_fmt = case_when(
      p_fdr < .001 ~ "< .001 ***",
      p_fdr < .01  ~ paste0(round(p_fdr, 3), " **"),
      p_fdr < .05  ~ paste0(round(p_fdr, 3), " *"),
      TRUE         ~ as.character(round(p_fdr, 3))
    )
  ) %>%
  arrange(ROI, task, Effect) %>%
  select(ROI, task, Effect, DFn, DFd, F, p_fdr_fmt, ges) %>%
  rename("p (FDR-corrected)" = p_fdr_fmt)

# ── Interaction summary wide table ───────────────────────────
summary_wide <- anova_table %>%
  filter(Effect == "prosody_val:semantic_val") %>%
  mutate(
    ROI  = factor(ROI,  levels = roi_order_shared),
    task = factor(task, levels = task_order_shared),
    cell = case_when(
      p_fdr < .001 ~ paste0("F(", DFn, ",", DFd, ")=", round(F, 2), ", p < .001 ***"),
      p_fdr < .01  ~ paste0("F(", DFn, ",", DFd, ")=", round(F, 2), ", p=", round(p_fdr, 3), " **"),
      p_fdr < .05  ~ paste0("F(", DFn, ",", DFd, ")=", round(F, 2), ", p=", round(p_fdr, 3), " *"),
      TRUE         ~ "ns"
    )
  ) %>%
  select(ROI, task, cell) %>%
  pivot_wider(names_from = task, values_from = cell) %>%
  arrange(ROI) %>%
  select(ROI, any_of(task_order_shared))

# ── Write to Excel ────────────────────────────────────────────
wb <- createWorkbook()

addWorksheet(wb, "Full_ANOVA")
writeData(wb, "Full_ANOVA", anova_clean)

addWorksheet(wb, "Interaction_Summary")
writeData(wb, "Interaction_Summary", summary_wide)

for (roi in roi_order_shared) {
  sheet_name <- gsub("[\\[\\], -]", "", roi)
  sheet_name <- substr(sheet_name, 1, 31)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name,
            anova_clean %>% filter(ROI == roi) %>% select(-ROI))
}

saveWorkbook(wb,
             file.path(results_dir, "ANOVA_results_full.xlsx"),
             overwrite = TRUE)
cat("\nFull ANOVA results saved.\n")

# ============================================================
# HELPER: significance label (FDR-corrected)
# ============================================================
sig_label <- function(p) {
  case_when(
    p < .001 ~ "***",
    p < .01  ~ "**",
    p < .05  ~ "*",
    TRUE     ~ ""
  )
}

# ============================================================
# HELPER: bracket annotations using FDR-corrected p
# ============================================================
get_interaction_annot <- function(anova_tbl, roi, tsk) {
  row <- anova_tbl %>%
    filter(ROI == roi, task == tsk,
           Effect == "prosody_val:semantic_val")
  if (nrow(row) == 0) return(NULL)
  list(
    F_val = round(row$F, 2),
    df1   = row$DFn,
    df2   = row$DFd,
    p     = row$p_fdr,      # ← FDR-corrected p
    label = sig_label(row$p_fdr)
  )
}

# ============================================================
# SUMMARY DATA  (mean ± SE per cell)
# ============================================================
cond_summary <- all_4cond %>%
  group_by(ROI, task, prosody_val, semantic_val) %>%
  summarise(
    mean = mean(beta, na.rm = TRUE),
    se   = sd(beta,  na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# ============================================================
# PLOT FUNCTION
# ============================================================
make_2x2_plot <- function(data_raw,
                          data_summary,
                          anova_tbl,
                          roi_order,
                          task_order,
                          plot_title = "2x2: Prosody x Semantic valence") {
  
  d_raw <- data_raw %>%
    filter(ROI %in% roi_order, task %in% task_order) %>%
    mutate(
      ROI          = factor(ROI,          levels = roi_order),
      task         = factor(task,         levels = task_order),
      prosody_val  = factor(prosody_val,  levels = c("positive", "negative")),
      semantic_val = factor(semantic_val, levels = c("positive", "negative"))
    )
  
  d_sum <- data_summary %>%
    filter(ROI %in% roi_order, task %in% task_order) %>%
    mutate(
      ROI          = factor(ROI,          levels = roi_order),
      task         = factor(task,         levels = task_order),
      prosody_val  = factor(prosody_val,  levels = c("positive", "negative")),
      semantic_val = factor(semantic_val, levels = c("positive", "negative"))
    )
  
  # Y-axis is driven by the summary statistics (means ± SE) plus the
  # significance bracket — NOT by the spread of individual betas. This keeps the
  # bars filling each panel so the Prosody × Semantics interaction is legible;
  # a few extreme raw dots fall outside the window and are cropped (see the
  # d_raw clipping below). `sum_range` sets all vertical spacing so brackets and
  # labels scale with the means rather than the (much larger) raw spread.
  y_ceil <- d_sum %>%
    group_by(ROI, task) %>%
    summarise(
      y_top_sum = max(mean + se, na.rm = TRUE),
      y_bot_sum = min(mean - se, na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    mutate(
      y_top_sum = pmax(y_top_sum, 0),          # always keep the 0 baseline in view
      y_bot_sum = pmin(y_bot_sum, 0),
      sum_range = y_top_sum - y_bot_sum,
      y_top         = y_top_sum,
      y_dot_top     = y_top + 0.18 * sum_range,        # crop dots just above the error bars
      y_bracket_bot = y_dot_top + 0.06 * sum_range,    # bracket sits in clear space above dots
      y_bracket_top = y_dot_top + 0.14 * sum_range,
      y_label       = y_dot_top + 0.24 * sum_range,
      y_axis_max    = y_dot_top + 0.36 * sum_range,
      y_axis_min    = y_bot_sum - 0.18 * sum_range
    )

  y_expand_df <- y_ceil %>%
    select(ROI, task, y_axis_max, y_axis_min) %>%
    mutate(
      ROI  = factor(ROI,  levels = roi_order),
      task = factor(task, levels = task_order)
    ) %>%
    pivot_longer(cols = c(y_axis_max, y_axis_min),
                 names_to = "bound", values_to = "y_expand")

  # Crop raw dots to the summary-driven window so they no longer stretch the
  # axis (free_y autoscales to plotted data, so out-of-range points must be
  # dropped rather than relying on coord_cartesian, which can't clip per facet).
  d_raw <- d_raw %>%
    left_join(select(y_ceil, ROI, task, y_axis_min, y_dot_top),
              by = c("ROI", "task")) %>%
    filter(beta >= y_axis_min, beta <= y_dot_top) %>%
    select(-y_axis_min, -y_dot_top)
  
  annot_df <- map_dfr(roi_order, function(r) {
    map_dfr(task_order, function(t) {
      info <- get_interaction_annot(anova_tbl, r, t)
      if (is.null(info) || info$label == "") return(NULL)
      yc <- filter(y_ceil, ROI == r, task == t)
      tibble(
        ROI           = r,
        task          = t,
        label         = info$label,
        y_bracket_bot = yc$y_bracket_bot,
        y_bracket_top = yc$y_bracket_top,
        y_label       = yc$y_label
      )
    })
  }) %>%
    mutate(
      ROI  = factor(ROI,  levels = roi_order),
      task = factor(task, levels = task_order)
    )
  
  dodge_width <- 0.65
  pd_jitter   <- position_jitterdodge(
    jitter.width  = 0.07,
    jitter.height = 0,
    dodge.width   = dodge_width,
    seed          = 42
  )
  pd_clean <- position_dodge(width = dodge_width)
  
  p <- ggplot(mapping = aes(x = prosody_val, colour = semantic_val,
                            fill = semantic_val)) +
    
    geom_blank(
      data        = y_expand_df,
      aes(x = "positive", y = y_expand),
      inherit.aes = FALSE
    ) +
    
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey65", linewidth = 0.35) +
    
    geom_bar(
      data     = d_sum,
      aes(y = mean),
      position = pd_clean,
      stat     = "identity",
      width    = 0.55,
      alpha    = 0.50
    ) +
    
    geom_point(
      data     = d_raw,
      aes(y = beta),
      position = pd_jitter,
      shape    = 16,
      size     = 1.1,
      alpha    = 0.45,
      stroke   = 0
    ) +
    
    geom_errorbar(
      data      = d_sum,
      aes(y = mean, ymin = mean - se, ymax = mean + se),
      position  = pd_clean,
      width     = 0.16,
      linewidth = 0.55,
      colour    = "grey20"
    ) +
    
    geom_point(
      data     = d_sum,
      aes(y = mean),
      position = pd_clean,
      shape    = 18,
      size     = 3.0
    ) +
    
    geom_segment(
      data = annot_df,
      aes(x = 1.08, xend = 1.08,
          y = y_bracket_bot, yend = y_bracket_top),
      colour = "grey30", linewidth = 0.45,
      inherit.aes = FALSE
    ) +
    
    geom_segment(
      data = annot_df,
      aes(x = 1.92, xend = 1.92,
          y = y_bracket_bot, yend = y_bracket_top),
      colour = "grey30", linewidth = 0.45,
      inherit.aes = FALSE
    ) +
    
    geom_segment(
      data = annot_df,
      aes(x = 1.08, xend = 1.92,
          y = y_bracket_top, yend = y_bracket_top),
      colour = "grey30", linewidth = 0.45,
      inherit.aes = FALSE
    ) +
    
    geom_text(
      data = annot_df,
      aes(x = 1.5, y = y_label, label = label),
      size      = 3.5,
      colour    = "grey20",
      hjust     = 0.5,
      fontface  = "bold",
      inherit.aes = FALSE
    ) +
    
    facet_grid(
      ROI ~ task,
      scales   = "free_y",
      labeller = labeller(task = task_labels)
    ) +
    
    scale_colour_manual(values = cb_palette, name = "Semantic valence") +
    scale_fill_manual(  values = cb_palette, name = "Semantic valence") +
    
    theme_bw(base_size = 8) +
    theme(
      plot.background    = element_rect(fill = "white", colour = NA),
      panel.background   = element_rect(fill = "white", colour = NA),
      strip.text.x       = element_text(size = 9, face = "bold"),
      strip.text.y       = element_text(size = 9, face = "bold", angle = 0),
      strip.background   = element_rect(fill = "grey92", colour = NA),
      legend.position    = "bottom",
      legend.key.size    = unit(4, "mm"),
      legend.text        = element_text(size = 9),
      legend.title       = element_text(size = 9, face = "bold"),
      legend.background  = element_rect(fill = "white", colour = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x        = element_text(angle = 30, hjust = 1, size = 9),
      axis.text.y        = element_text(size = 9),
      axis.title         = element_text(size = 10),
      plot.caption       = element_text(size = 6, colour = "grey40", hjust = 0)
    ) +
    
    labs(
      title   = NULL,
      x       = "Prosody valence",
      y       = expression("Signal change (%) ± SE"),
    )
  
  p
}

# ============================================================
# ALL-ROI PLOT
# ============================================================
p_all <- make_2x2_plot(
  data_raw     = all_4cond,
  data_summary = cond_summary,
  anova_tbl    = anova_table,
  roi_order    = names(roi_files),
  task_order   = task_levels
)

ggsave(
  file.path(results_dir, "2x2_interaction_plot_all_ROIs.png"),
  p_all,
  width = 180, height = 200, units = "mm", dpi = 300
)
cat("All-ROI plot saved.\n")

# ============================================================
# LEFT IFG ONLY PLOT  (Figure 4A — manuscript figure)
# ------------------------------------------------------------
# This standalone left-IFG 2x2 panel IS the manuscript figure: the title and
# enlarged export below are folded straight in, so there is no separate
# "Figure 4A" object downstream. Satisfies R5.10 (Panel A enlarged on its own).
# ============================================================
p_leftIFG <- make_2x2_plot(
  data_raw     = all_4cond,
  data_summary = cond_summary,
  anova_tbl    = anova_table,
  roi_order    = "lIFG [-44, 22, -8]",
  task_order   = task_levels
) +
  ggtitle("Beta estimates of the Prosody × Semantics interaction in the left IFG") +
  theme(
    plot.margin = margin(6, 6, 4, 6),
    plot.title  = element_text(size = 11, colour = "black", face = "plain",
                               hjust = 0, margin = margin(b = 6))
  )

ggsave(
  file.path(results_dir, "Figure4A_lIFG_betas.jpg"),
  p_leftIFG,
  width = 180, height = 60, units = "mm", dpi = 600   # good 180x55 proportions + room for the title
)
cat("Figure 4A (lIFG betas) saved.\n")

# ============================================================
# HEATMAP — FDR-corrected p-values
# ============================================================
heatmap_df <- anova_table %>%
  filter(Effect == "prosody_val:semantic_val") %>%
  mutate(
    sig_bool   = p_fdr < .05,           # ← FDR threshold
    F_plot     = ifelse(sig_bool, F, NA_real_),
    ROI        = factor(ROI, levels = roi_order_shared),
    task       = factor(task, levels = task_order_shared),
    stars      = sig_label(p_fdr),      # ← FDR-based stars
    tile_label = ifelse(sig_bool,
                        paste0(stars, "\nF=", round(F, 1)),
                        "ns")
  )

f_max <- max(heatmap_df$F_plot, na.rm = TRUE)

p_heatmap <- ggplot(heatmap_df, aes(x = task, y = ROI)) +
  
  geom_tile(
    data   = filter(heatmap_df, !sig_bool),
    fill   = "grey92", colour = "white", linewidth = 0.6
  ) +
  
  scale_y_discrete(limits = rev(roi_order_shared)) +
  
  geom_tile(
    data   = filter(heatmap_df, sig_bool),
    aes(fill = F_plot),
    colour = "white", linewidth = 0.6
  ) +
  
  geom_tile(
    data      = filter(heatmap_df, ROI == "lIFG [-44, 22, -8]"),
    fill      = NA,
    colour    = "#0072B2",
    linewidth = 1.2
  ) +
  
  geom_text(
    aes(label  = tile_label,
        colour = ifelse(sig_bool & !is.na(F_plot) & F_plot > f_max * 0.55,
                        "white", "grey20")),
    size       = 3.5,
    lineheight = 0.85,
    fontface   = "bold"
  ) +
  scale_colour_identity() +
  
  scale_fill_gradientn(
    colours  = c("#D6EAF8", "#2E86C1", "#1A5276"),
    na.value = "grey92",
    limits   = c(0, f_max),
    name     = "F-value\n(interaction)"
  ) +
  
  scale_x_discrete(
    position = "top",
    limits   = task_order_shared,
    labels   = task_labels
  ) +
  
  theme_minimal(base_size = 9) +
  theme(
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.grid        = element_blank(),
    axis.text.x       = element_text(size = 9,  face = "bold", colour = "grey20"),
    axis.text.y       = element_text(size = 9,  colour = "grey20"),
    axis.title        = element_blank(),
    legend.position   = "right",
    legend.key.height = unit(12, "mm"),
    legend.key.width  = unit(3,  "mm"),
    legend.text       = element_text(size = 9),
    legend.title      = element_text(size = 9,  face = "bold"),
    plot.caption      = element_text(size = 6, colour = "grey40", hjust = 0),
    plot.margin       = margin(4, 4, 4, 4)
  )

ggsave(
  file.path(results_dir, "heatmap_interaction_allROIs.png"),
  p_heatmap,
  width = 130, height = 100, units = "mm", dpi = 300
)
cat("Heatmap saved.\n")

# ============================================================
# Glass-brain ROI localisation (horizontal) — loaded here, used by the magnitude
# figure below. The former combined 3-panel figure was removed; the manuscript
# now has Figure 4 (magnitude + ROI localisation), Figure 5 (lIFG betas), and
# Figure 6 (heatmap).
# ============================================================
glass_brain_img <- readPNG("plots/s3_Betas/roi_glass_brain_horizontal.png")

# ============================================================
# SPLIT FIGURE 4 FOR LEGIBILITY (Reviewer R5.10)
# ------------------------------------------------------------
# The Reviewer asked to separate Panel A from Panel B so each can be
# enlarged, and noted Panel C may be combined with Panel B or stand alone.
# We therefore export (4A) the left-IFG beta panel on its own and
# (4B) the F-value heatmap combined with the glass-brain localisation.
# ============================================================

# ── Figure 4A is now built directly in the LEFT IFG ONLY PLOT section above ──

# ── Figure 4B — F-value heatmap, standalone & enlarged ─────────────────────
# Glass-brain ROI localisation now lives with the magnitude figure (it indexes
# the same 9 ROIs), so it is not duplicated here. This satisfies R5.10:
# Panel A (left-IFG betas) and Panel B (heatmap) are now separate figures.
fig4B <- p_heatmap +
  theme(plot.margin = margin(14, 4, 2, 4)) +
  ggtitle("F-values of the Prosody × Semantics interaction") +
  theme(plot.title = element_text(size = 11, colour = "black", face = "plain",
                                  hjust = 0, margin = margin(b = 4)))

ggsave(
  file.path(results_dir, "Figure4B_heatmap.jpg"),
  fig4B,
  width = 180, height = 130, units = "mm", dpi = 600
)
cat("Figure 4B (heatmap) saved.\n")

# ============================================================
# LITERAL vs NON-LITERAL — absolute activation, collapsed over tasks
# ------------------------------------------------------------
# Reviewer R3.4(a): the prosody x semantics interaction tells us where the
# four cells diverge, but a region can show a strong interaction while BOTH
# the literal and non-literal signals remain weak (or both strong). To show
# *for which regions there is more activation for non-literal than literal,
# independent of the interaction*, we collapse over the five tasks and over
# the two cells within each congruency level, then compare the congruency
# main effect (Literal = congruent CP_SP/CN_SN; Non-literal = incongruent
# CP_SN/CN_SP) per ROI.
# ============================================================
cong_levels  <- c("Literal", "Non-literal")
cong_palette <- c("Literal" = "#009E73", "Non-literal" = "#CC79A7")

# ── Per-subject value: mean over the 5 tasks + the 2 cells per congruency ──
congruency_subj <- all_4cond %>%
  mutate(congruency = ifelse(prosody_val == semantic_val,
                             "Literal", "Non-literal")) %>%
  group_by(ROI, id, congruency) %>%
  summarise(beta = mean(beta, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    ROI        = factor(ROI, levels = roi_order_shared),
    congruency = factor(congruency, levels = cong_levels)
  )

# ── Paired t-test per ROI (Non-literal − Literal) + FDR across ROIs ────────
congruency_wide <- congruency_subj %>%
  pivot_wider(names_from = congruency, values_from = beta)

congruency_stats <- congruency_wide %>%
  group_by(ROI) %>%
  summarise(
    mean_literal     = mean(Literal,        na.rm = TRUE),
    mean_nonliteral  = mean(`Non-literal`,  na.rm = TRUE),
    diff             = mean(`Non-literal` - Literal, na.rm = TRUE),
    t                = t.test(`Non-literal`, Literal, paired = TRUE)$statistic,
    df               = t.test(`Non-literal`, Literal, paired = TRUE)$parameter,
    p                = t.test(`Non-literal`, Literal, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_fdr = p.adjust(p, method = "fdr"),
    stars = sig_label(p_fdr),
    ROI   = factor(ROI, levels = roi_order_shared)
  ) %>%
  arrange(ROI)

# ── Mean ± SE per cell for plotting ───────────────────────────────────────
congruency_summary <- congruency_subj %>%
  group_by(ROI, congruency) %>%
  summarise(
    mean = mean(beta, na.rm = TRUE),
    se   = sd(beta,  na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# ── Export the table — DESCRIPTIVE ONLY ───────────────────────────────────
# Inferential columns (t, p, p_fdr) are intentionally dropped: the ROIs were
# selected from the Non-literal > Literal contrast, so per-ROI NL>L tests are
# circular (double-dipping) and must not be reported. Only the absolute means
# and their difference are descriptive and safe to report (R3.4-1).
write.csv(
  congruency_stats %>%
    select(ROI, mean_literal, mean_nonliteral, diff) %>%
    mutate(across(c(mean_literal, mean_nonliteral, diff), ~ round(.x, 4))),
  file.path(results_dir, "literal_vs_nonliteral_byROI.csv"),
  row.names = FALSE
)
cat("Literal vs Non-literal table saved (descriptive only).\n")

# ── Plotting positions ────────────────────────────────────────────────────
# NOTE: no significance stars are drawn here. The 9 ROIs were DEFINED from the
# Non-literal > Literal contrast, so a per-ROI NL>L test is significant by
# construction (double-dipping / circular) and is NOT reported as an
# independent effect (see manuscript l.156, response R3.4a). This figure is
# purely DESCRIPTIVE: it shows the inter-region magnitude gradient (R3.4-1).
dodge_w <- 0.62
pd_pts  <- position_jitterdodge(jitter.width = 0.10, jitter.height = 0,
                                dodge.width = dodge_w, seed = 42)
pd_bar  <- position_dodge(width = dodge_w)

p_congruency <- ggplot(mapping = aes(x = ROI, fill = congruency,
                                     colour = congruency)) +

  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey65", linewidth = 0.35) +

  geom_bar(
    data     = congruency_summary,
    aes(y = mean),
    stat     = "identity",
    position = pd_bar,
    width    = 0.55,
    alpha    = 0.55
  ) +

  geom_point(
    data     = congruency_subj,
    aes(y = beta),
    position = pd_pts,
    shape    = 16, size = 0.8, alpha = 0.30, stroke = 0
  ) +

  geom_errorbar(
    data      = congruency_summary,
    aes(y = mean, ymin = mean - se, ymax = mean + se),
    position  = pd_bar,
    width     = 0.18, linewidth = 0.5, colour = "grey20"
  ) +

  scale_fill_manual(  values = cong_palette, name = NULL) +
  scale_colour_manual(values = cong_palette, name = NULL) +

  theme_bw(base_size = 9) +
  theme(
    plot.background    = element_rect(fill = "white", colour = NA),
    panel.background   = element_rect(fill = "white", colour = NA),
    legend.position    = "bottom",
    legend.key.size    = unit(4, "mm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(angle = 30, hjust = 1, size = 8),
    axis.text.y        = element_text(size = 9),
    axis.title         = element_text(size = 10)
  ) +

  labs(
    x = NULL,
    y = expression("Signal change (%) ± SE")
  )

ggsave(
  file.path(results_dir, "literal_vs_nonliteral_allROIs.png"),
  p_congruency,
  width = 180, height = 95, units = "mm", dpi = 300
)
cat("Literal vs Non-literal plot saved.\n")

# ============================================================
# MAGNITUDE FIGURE — horizontal glass brain (top) + NL vs L bars (bottom)
# ------------------------------------------------------------
# Answers R3.4-1 (which regions are more / less activated in non-literal than
# literal speech, and where they sit) by pairing the anatomical localisation
# with the inter-ROI magnitude gradient. Stacked VERTICALLY with the horizontal
# glass brain (A) as a slim strip on top and the magnitude bars (B) keeping the
# FULL figure width beneath — they then read large once the figure is set to
# \textwidth in the manuscript. Purely DESCRIPTIVE — no significance stars, since
# the ROIs were selected from the Non-literal > Literal contrast (see l.156).
# Uses the horizontal glass brain (glass_brain_img) from ROI_plot.py.
# ============================================================
panel_loc <- ggdraw() +
  draw_image(glass_brain_img, x = 0, y = 0, width = 1, height = 1)

magnitude_fig <- plot_grid(
  panel_loc,
  p_congruency,
  nrow           = 2,
  rel_heights    = c(0.5, 1),   # brain a slim strip on top; bars take full width and most of the height below
  labels         = c("A", "B"),
  label_size     = 12,
  label_fontface = "bold",
  label_x        = 0,
  label_y        = 1
)

ggsave(
  file.path(results_dir, "Figure_magnitude_NLvsL_glassbrain.jpg"),
  magnitude_fig,
  width = 190, height = 150, units = "mm", dpi = 600
)
cat("Magnitude figure (NL vs L + glass brain) saved.\n")