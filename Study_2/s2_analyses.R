## Study 2

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

files <- Sys.glob("data/*/*.txt")
file_list <- as.list(files)

read <- function(x) { read.table(x, header = TRUE) }

df_list <- lapply(file_list, read)
s2 <- do.call(rbind, df_list)
s2 <- s2 %>% mutate(across(where(is.character), as.factor))

# Individual variables
mean(s2$age)
sd(s2$age)
mean(s2$Time_Context)
sd(s2$Time_Context)
mean(s2$Time_Statement)
sd(s2$Time_Statement)

# Exclusion of a bugged stimulus
s2 <- s2[!s2$Statement == "SNnegf1_1.wav", ]

# Only valid observations
s2 <- s2 %>% filter(Validation_Type %in% c("valid", "directvalide") |
                      Validation_Type_2 %in% c("valid", "directvalide"))

# Add useful variables
s2$cond_sit <- paste0(s2$Condition_name, s2$Situation)
s2$con  <- as.factor(substr(s2$Condition_name, start = 1, stop = 2))
s2$sem  <- as.factor(substr(s2$Condition_name, start = 4, stop = 5))
s2$pros <- as.factor(substr(s2$Condition_name, start = 6, stop = 8))

# ============================================================
# OUTLIER REMOVAL
# ============================================================

identify_outliers <- function(x, coef = 1.5) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  which(x > q3 + coef * iqr | x < q1 - coef * iqr)
}

remove_outliers <- function(data, score_var) {
  total_outliers <- 0
  for (task in unique(data$task)) {
    for (cond in unique(data$cond_sit)) {
      idx <- data$task == task & data$cond_sit == cond
      outliers <- identify_outliers(data[[score_var]][idx])
      data[[score_var]][idx][outliers] <- NA
      total_outliers <- total_outliers + length(outliers)
    }
  }
  message(paste("Removed a total of", total_outliers, "outliers from", score_var))
  data
}

s2 <- remove_outliers(s2, "Evaluation_Score")
s2 <- remove_outliers(s2, "Evaluation_Score_2")

# ============================================================
# COLORBLIND-FRIENDLY PALETTES (Okabe-Ito)
# ============================================================

okabe_sem <- c("Negative" = "#E69F00", "Positive" = "#56B4E9")

# ============================================================
# MODELS
# ============================================================

run_model <- function(data, score_var) {
  lmer(as.formula(paste0(score_var, " ~ con*sem*pros +
       (1 + con + sem + pros|subID) + (1|Context) + (1|Statement)")),
       data = data)
}

get_stats <- function(model) {
  anova_res <- data.frame(anova(model, type = 2)) %>%
    mutate(across(everything(), ~ round(., 3)))
  r2_vals <- r2(model)
  list(anova       = anova_res,
       marginal    = r2_vals$R2_marginal,
       conditional = r2_vals$R2_conditional)
}

data_pros <- s2 %>% filter(task == "prosody")
data_sem  <- s2 %>% filter(task == "semantic")
data_iro  <- s2 %>% filter(task == "sarcasm")
data_sar  <- s2 %>% filter(task == "sarcasm")
data_tom  <- s2 %>% filter(task == "tom")

pros_model <- run_model(data_pros, "Evaluation_Score")
sem_model  <- run_model(data_sem,  "Evaluation_Score")
iro_model  <- run_model(data_iro,  "Evaluation_Score")
sar_model  <- run_model(data_sar,  "Evaluation_Score_2")
tom_model  <- run_model(data_tom,  "Evaluation_Score")

stats_pros <- get_stats(pros_model)
stats_sem  <- get_stats(sem_model)
stats_iro  <- get_stats(iro_model)
stats_sar  <- get_stats(sar_model)
stats_tom  <- get_stats(tom_model)

# Estimated marginal means - uncomment to compute
# mm_pros <- data.frame(emmeans(pros_model, ~ con*sem*pros)) %>%
#   mutate(across(-c(1:3), ~ round(., 3))) %>%
#   mutate(ic = paste0("[", pull(., 7), " ; ", pull(., 8), "]"))

# mm_sem <- data.frame(emmeans(sem_model, ~ con*sem*pros)) %>%
#   mutate(across(-c(1:3), ~ round(., 3))) %>%
#   mutate(ic = paste0("[", pull(., 7), " ; ", pull(., 8), "]"))

# mm_iro <- data.frame(emmeans(iro_model, ~ con*sem*pros)) %>%
#   mutate(across(-c(1:3), ~ round(., 3))) %>%
#   mutate(ic = paste0("[", pull(., 7), " ; ", pull(., 8), "]"))

# mm_sar <- data.frame(emmeans(sar_model, ~ con*sem*pros)) %>%
#   mutate(across(-c(1:3), ~ round(., 3))) %>%
#   mutate(ic = paste0("[", pull(., 7), " ; ", pull(., 8), "]"))

# mm_tom <- data.frame(emmeans(tom_model, ~ con*sem*pros)) %>%
#   mutate(across(-c(1:3), ~ round(., 3))) %>%
#   mutate(ic = paste0("[", pull(., 7), " ; ", pull(., 8), "]"))

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
    facet_wrap(~ con, scales = "fixed") +
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

data_sar <- data_sar[!is.na(data_sar$Evaluation_Score_2), ]
data_sar$pred <- predict(sar_model, type = "response")
data_sar <- recode_vars(data_sar)

data_tom <- data_tom[!is.na(data_tom$Evaluation_Score), ]
data_tom$pred <- predict(tom_model, type = "response")
data_tom <- recode_vars(data_tom)

# ============================================================
# TASK LIST FOR PLOT LOOP
# Note: sarcasm uses Evaluation_Score_2 for data plot
# ============================================================

tasks <- list(
  list(data = data_pros, panel = "A)", title = "Evaluation of the prosody valence",
       score_var = "Evaluation_Score",
       path_m = "plots/s2_model/prosody_model.jpg",
       path_d = "plots/s2_data/prosody_data.jpg"),
  list(data = data_sem,  panel = "B)", title = "Evaluation of the semantic valence",
       score_var = "Evaluation_Score",
       path_m = "plots/s2_model/semantic_model.jpg",
       path_d = "plots/s2_data/semantic_data.jpg"),
  list(data = data_iro,  panel = "C)", title = "Evaluation of the degree of perceived irony",
       score_var = "Evaluation_Score",
       path_m = "plots/s2_model/irony_model.jpg",
       path_d = "plots/s2_data/irony_data.jpg"),
  list(data = data_sar,  panel = "D)", title = "Evaluation of the degree of perceived sarcasm",
       score_var = "Evaluation_Score_2",
       path_m = "plots/s2_model/sarcasm_model.jpg",
       path_d = "plots/s2_data/sarcasm_data.jpg"),
  list(data = data_tom,  panel = "E)", title = "Evaluation of the degree of uncomfortableness",
       score_var = "Evaluation_Score",
       path_m = "plots/s2_model/tom_mod.jpg",
       path_d = "plots/s2_data/tom_data.jpg")
)

# ============================================================
# SAVE INDIVIDUAL PLOTS + COLLECT FOR COMBINED FIGURES
# ============================================================

plots_model <- list()
plots_data  <- list()

for (t in tasks) {
  
  # Individual clean plots
  p_model <- make_plot(t$data, "pred",      t$title, "Evaluation score")
  p_data  <- make_plot(t$data, t$score_var, t$title, "Evaluation score")
  ggsave(t$path_m, plot = p_model, width = 7, height = 7, dpi = 300)
  ggsave(t$path_d, plot = p_data,  width = 7, height = 7, dpi = 300)
  
  # Store for combined figures
  plots_model[[t$panel]] <- make_plot(t$data, "pred",      t$title, "Evaluation score")
  plots_data[[t$panel]]  <- make_plot(t$data, t$score_var, t$title, "Evaluation score")
}

# ============================================================
# COMBINED FIGURES - 2 plots per row, last row centered
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

make_combined(plots_model, save_path = "plots/s2_model/combined_model.jpg")

make_combined(plots_data, save_path = "plots/s2_data/combined_data.jpg")