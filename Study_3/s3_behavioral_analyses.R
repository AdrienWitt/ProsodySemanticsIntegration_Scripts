#Study 3
library(dplyr)
library(ggplot2)
library(car)
library(lme4)
library(ggpubr)
library(emmeans)
library(lmerTest)
library(performance)
library(tidyverse)


# ============================================================
# DATA LOADING
# ============================================================

files <- Sys.glob("data/*/Resultfile*.txt")
file_list <- as.list(files)

read <- function(x) { read.table(x, header = TRUE) }

df_list <- lapply(file_list, read)
s3 <- do.call(rbind, df_list)
s3 <- s3 %>% mutate(across(where(is.character), as.factor))

# Remove left handed
s3 <- subset(s3, !(subID %in% c(11, 13, 15, 17, 22)))

# Individual variables
mean(s3$age)
sd(s3$age)
table(s3$genre)

# ============================================================
# COLORBLIND-FRIENDLY PALETTES (Okabe-Ito)
# ============================================================

okabe_sem   <- c("Negative" = "#E69F00", "Positive" = "#56B4E9")
okabe_valid <- c("Valid"     = "#009E73",
                 "Non-valid" = "#E69F00",
                 "No key"    = "#CC79A7")

# ============================================================
# RT ANALYSES
# ============================================================

rt_run <- s3 %>%
  group_by(run) %>%
  summarise(
    mean_RT = mean(RT, na.rm = TRUE),
    se_RT   = sd(RT, na.rm = TRUE) / sqrt(sum(!is.na(RT)))
  )

rt_task <- s3 %>%
  filter(!is.na(task)) %>%
  group_by(task) %>%
  summarise(
    mean_RT = mean(RT, na.rm = TRUE),
    se_RT   = sd(RT, na.rm = TRUE) / sqrt(sum(!is.na(RT)))
  )

rt_task$task <- factor(rt_task$task,
                       levels = c("prosody", "semantic", "sarcasm", "irony", "tom"))

# Linear mixed models
s3$run_numeric <- as.numeric(as.factor(s3$run))
lm_run  <- lmer(RT ~ run_numeric + (1|subID), data = s3)
summary(lm_run)

lm_task <- lmer(RT ~ task + (1|subID), data = s3)
summary(lm_task)

# ============================================================
# RT PLOTS
# ============================================================

p_run <- ggplot(rt_run, aes(x = run, y = mean_RT)) +
  geom_point(size = 3, color = "#0072B2") +
  geom_errorbar(aes(ymin = mean_RT - se_RT, ymax = mean_RT + se_RT),
                width = 0.2, color = "#0072B2") +
  geom_smooth(aes(group = 1), method = "lm",
              color = "#0072B2", fill = "#56B4E9", alpha = 0.2) +
  theme_classic() +
  labs(x = "Run", y = "Mean Reaction Time (s)",
       title = "Reaction Time across Runs")

p_task <- ggplot(rt_task, aes(x = task, y = mean_RT)) +
  geom_bar(stat = "identity", fill = "#0072B2", alpha = 0.8, width = 0.6) +
  geom_errorbar(aes(ymin = mean_RT - se_RT, ymax = mean_RT + se_RT),
                width = 0.2, color = "#0072B2") +
  theme_classic() +
  labs(x = "Task", y = "Mean Reaction Time (s)",
       title = "Reaction Time across Tasks")

# ============================================================
# VALIDATION ANALYSIS
# ============================================================

valid_run <- s3 %>%
  group_by(run) %>%
  summarise(
    n_total    = n(),
    n_valid    = sum(Validation_Type == "valid",    na.rm = TRUE),
    n_nonvalid = sum(Validation_Type == "nonvalid", na.rm = TRUE),
    n_nokey    = sum(Validation_Type == "nokey",    na.rm = TRUE),
    prop_valid    = n_valid    / n_total,
    prop_nonvalid = n_nonvalid / n_total,
    prop_nokey    = n_nokey    / n_total
  )

valid_long <- valid_run %>%
  select(run,
         "Valid"     = prop_valid,
         "Non-valid" = prop_nonvalid,
         "No key"    = prop_nokey) %>%
  pivot_longer(cols = c("Valid", "Non-valid", "No key"),
               names_to  = "type",
               values_to = "proportion")

p_valid <- ggplot(valid_long, aes(x = run, y = proportion,
                                  color = type, group = type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  theme_classic() +
  scale_color_manual(values = okabe_valid) +
  labs(x = "Run", y = "Proportion", color = "Response type",
       title = "Response validity across runs")

ggarrange(p_task, p_run, p_valid,
          ncol = 2, nrow = 2,
          labels = c("A", "B", "C"))

ggsave("plots/supplementary_RT_validation.jpg", width = 12, height = 10)

# ============================================================
# FILTER VALID OBSERVATIONS AND ADD VARIABLES
# ============================================================

s3 <- s3 %>% filter(Validation_Type %in% c("valid", "directvalide"))

s3$cond_sit <- paste0(s3$Condition_name, s3$Situation)
s3$con  <- as.factor(substr(s3$Condition_name, start = 1, stop = 2))
s3$sem  <- as.factor(substr(s3$Condition_name, start = 4, stop = 5))
s3$pros <- as.factor(substr(s3$Condition_name, start = 6, stop = 8))

# ============================================================
# MODELS
# ============================================================

run_model <- function(data) {
  lmer(Evaluation_Score ~ sem*pros +
         (1 + pros + sem|subID) + (1|Context) + (1|Statement),
       data = data)
}

data_pros <- s3 %>% filter(task == "prosody")
data_sem  <- s3 %>% filter(task == "semantic")
data_iro  <- s3 %>% filter(task == "irony")
data_sar  <- s3 %>% filter(task == "sarcasm")
data_tom  <- s3 %>% filter(task == "tom")

pros_model <- run_model(data_pros)
sem_model  <- run_model(data_sem)
iro_model  <- run_model(data_iro)
sar_model  <- run_model(data_sar)
tom_model  <- run_model(data_tom)

get_stats <- function(model) {
  anova_res <- data.frame(anova(model, type = 2)) %>%
    mutate(across(everything(), ~ round(., 3)))
  r2_vals <- r2(model)
  list(anova       = anova_res,
       marginal    = r2_vals$R2_marginal,
       conditional = r2_vals$R2_conditional)
}

stats_pros <- get_stats(pros_model)
stats_sem  <- get_stats(sem_model)
stats_iro  <- get_stats(iro_model)
stats_sar  <- get_stats(sar_model)
stats_tom  <- get_stats(tom_model)

# Estimated marginal means for each model
# Uncomment to compute — can be slow and may produce singularity warnings

# mm_pros <- data.frame(emmeans(pros_model, ~ sem*pros)) %>%
#   mutate(across(-c(1:2), ~ round(., 3))) %>%
#   mutate(ic = paste0("[", pull(., 6), " ; ", pull(., 7), "]"))

# mm_sem <- data.frame(emmeans(sem_model, ~ sem*pros)) %>%
#   mutate(across(-c(1:2), ~ round(., 3))) %>%
#   mutate(ic = paste0("[", pull(., 6), " ; ", pull(., 7), "]"))

# mm_iro <- data.frame(emmeans(iro_model, ~ sem*pros)) %>%
#   mutate(across(-c(1:2), ~ round(., 3))) %>%
#   mutate(ic = paste0("[", pull(., 6), " ; ", pull(., 7), "]"))

# mm_sar <- data.frame(emmeans(sar_model, ~ sem*pros)) %>%
#   mutate(across(-c(1:2), ~ round(., 3))) %>%
#   mutate(ic = paste0("[", pull(., 6), " ; ", pull(., 7), "]"))

# mm_tom <- data.frame(emmeans(tom_model, ~ sem*pros)) %>%
#   mutate(across(-c(1:2), ~ round(., 3))) %>%
#   mutate(ic = paste0("[", pull(., 6), " ; ", pull(., 7), "]"))

# ============================================================
# HELPER FUNCTIONS FOR PLOTS
# ============================================================

recode_vars <- function(data) {
  data %>% mutate(
    con      = fct_recode(con,  "Positive Context" = "CP", "Negative Context" = "CN"),
    Semantic = fct_recode(sem,  "Negative" = "SN", "Positive" = "SP"),
    Prosody  = fct_recode(pros, "Negative" = "neg", "Positive" = "pos")
  )
}

make_plot <- function(data, y_var, title, y_label) {
  ggplot(data, aes(x = Prosody, y = .data[[y_var]], fill = Semantic)) +
    geom_boxplot() +
    theme(
      axis.text.x      = element_text(angle = 30, vjust = 0.6, size = 12),
      axis.text.y      = element_text(size = 12),
      axis.title       = element_text(size = 13),
      plot.title       = element_text(hjust = 0.5, size = 14),
      strip.text       = element_text(size = 12),          # facet panel labels
      legend.title     = element_text(size = 12, face = "bold"),
      legend.text      = element_text(size = 11),
      panel.background = element_rect(fill = "white")
    ) +
    labs(title = title, x = "Prosody", y = y_label) +
    facet_wrap(~ con, scales = "free_x", drop = TRUE) +
    scale_fill_manual(values = okabe_sem)
}

# ============================================================
# PREPARE DATA FOR PLOTS
# ============================================================

data_pros <- data_pros[!is.na(data_pros$Evaluation_Score), ]
data_pros$pred <- predict(pros_model, type = "response")
data_pros <- recode_vars(data_pros)

data_sem <- data_sem[!is.na(data_sem$Evaluation_Score), ]
data_sem$pred <- predict(sem_model, type = "response")
data_sem <- recode_vars(data_sem)

data_iro <- data_iro[!is.na(data_iro$Evaluation_Score), ]
data_iro$pred <- predict(iro_model, type = "response")
data_iro <- recode_vars(data_iro)

data_sar <- data_sar[!is.na(data_sar$Evaluation_Score), ]
data_sar$pred <- predict(sar_model, type = "response")
data_sar <- recode_vars(data_sar)

data_tom <- data_tom[!is.na(data_tom$Evaluation_Score), ]
data_tom$pred <- predict(tom_model, type = "response")
data_tom <- recode_vars(data_tom)

# ============================================================
# SAVE ALL PLOTS
# ============================================================

tasks <- list(
  list(data = data_pros, title = "Evaluation of the prosody valence",
       panel = "A)", stats = stats_pros,
       path_m = "plots/s3_model/prosody_model.jpg",
       path_d = "plots/s3_data/prosody_data.jpg"),
  list(data = data_sem,  title = "Evaluation of the semantic valence",
       panel = "B)", stats = stats_sem,
       path_m = "plots/s3_model/semantic_model.jpg",
       path_d = "plots/s3_data/semantic_data.jpg"),
  list(data = data_iro,  title = "Evaluation of the degree of perceived irony",
       panel = "C)", stats = stats_iro,
       path_m = "plots/s3_model/irony_model.jpg",
       path_d = "plots/s3_data/irony_data.jpg"),
  list(data = data_sar,  title = "Evaluation of the degree of perceived sarcasm",
       panel = "D)", stats = stats_sar,
       path_m = "plots/s3_model/sarcasm_model.jpg",
       path_d = "plots/s3_data/sarcasm_data.jpg"),
  list(data = data_tom,  title = "Evaluation of the degree of uncomfortableness",
       panel = "E)", stats = stats_tom,
       path_m = "plots/s3_model/tom_mod.jpg",
       path_d = "plots/s3_data/tom_data.jpg")
)

# ============================================================
# SAVE INDIVIDUAL PLOTS + COLLECT FOR COMBINED FIGURES
# ============================================================

plots_model <- list()
plots_data  <- list()

for (t in tasks) {
  
  panel_title <- t$title
  
  # Individual clean plots
  p_model <- make_plot(t$data, "pred",             t$title, "Evaluation score")
  p_data  <- make_plot(t$data, "Evaluation_Score", t$title, "Evaluation score")
  ggsave(t$path_m, plot = p_model, width = 7, height = 7, dpi = 300)
  ggsave(t$path_d, plot = p_data,  width = 7, height = 7, dpi = 300)
  
  # Store for combined figures
  plots_model[[t$panel]] <- make_plot(t$data, "pred", panel_title, "Evaluation score")
  plots_data[[t$panel]]  <- make_plot(t$data, "Evaluation_Score", panel_title, "Evaluation score")
}

# ============================================================
# COMBINED FIGURES -
# ============================================================

shared_legend <- get_legend(
  plots_model[["A)"]] + theme(legend.position = "bottom",
                              legend.title = element_text(size = 10),
                              legend.text  = element_text(size = 9))
)

make_combined <- function(plots, save_path) {
  
  row1 <- ggarrange(
    plots[["A)"]] + theme(legend.position = "none") + 
      labs(tag = "A"),
    plots[["B)"]] + theme(legend.position = "none") + 
      labs(tag = "B"),
    ncol = 2
  )
  
  row2 <- ggarrange(
    plots[["C)"]] + theme(legend.position = "none") + 
      labs(tag = "C"),
    plots[["D)"]] + theme(legend.position = "none") + 
      labs(tag = "D"),
    ncol = 2
  )
  
  row3 <- ggarrange(
    ggplot() + theme_void(),
    plots[["E)"]] + theme(legend.position = "none") + 
      labs(tag = "E"),
    ggplot() + theme_void(),
    ncol = 3,
    widths = c(0.5, 1, 0.5)
  )
  
  fig <- ggarrange(
    row1, row2, row3,
    shared_legend,
    nrow = 4,
    heights = c(1, 1, 1, 0.1)
  )
  
  ggsave(save_path, width = 14, height = 14, dpi = 300)
}

make_combined(plots_model, save_path = "plots/s3_model/combined_model.jpg")

make_combined(plots_data, save_path = "plots/s3_data/combined_data.jpg")