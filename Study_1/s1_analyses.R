## Study 1

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

ps1 <- read.csv("data/matriceNewP.csv", sep = ";", row.names = NULL)
ps1 <- ps1 %>% mutate(across(where(is.character), as.factor))
ps1$sem  <- factor(ps1$sem,  levels = c("n", "pseudo", "p"))
ps1$pros <- factor(ps1$pros, levels = c("n", "m", "p"))

# ============================================================
# COLORBLIND-FRIENDLY PALETTE (Okabe-Ito)
# ============================================================

okabe_pros <- c("Negative" = "#E69F00",
                "Monotone" = "#999999",
                "Positive" = "#56B4E9")

# ============================================================
# ANALYSES
# ============================================================

# Context
data_context <- ps1 %>% filter(sentenceType == "context")

## Valence
vc       <- lmer(scoreV ~ sem + (1 + sem|ID) + (1|stimuli), data = data_context)
vc_anova <- round(data.frame(anova(vc, type = 2) %>% mutate(across(-6, ~ round(., 3)))), 3)
vc_emm   <- emmeans(vc, ~ sem)
vc_results <- pairs(vc_emm)
r2_values  <- r2(vc)
marginal_R2    <- r2_values$R2_marginal
conditional_R2 <- r2_values$R2_conditional

## Arousal
ac       <- lmer(scoreA ~ sem + (1 + sem|ID) + (1|stimuli), data = data_context)
ac_anova <- round(data.frame(anova(ac, type = 2) %>% mutate(across(-6, ~ round(., 3)))), 3)
ac_emm   <- emmeans(ac, ~ sem)
ac_results <- pairs(ac_emm)
r2_values  <- r2(ac)
marginal_R2    <- r2_values$R2_marginal
conditional_R2 <- r2_values$R2_conditional

# Statement
data_statement <- ps1 %>% filter(sentenceType == "statement")

## Valence
vs       <- lmer(scoreV ~ sem*pros + (1 + sem + pros|ID) + (1|stimuli), data = data_statement)
vs_anova <- data.frame(anova(vs, type = 2)) %>% mutate(across(-6, ~ round(., 3)))
r2_values <- r2(vs)
marginal_R2    <- r2_values$R2_marginal
conditional_R2 <- r2_values$R2_conditional

vs_pc      <- emmeans(vs, ~ pros)
vs_results <- data.frame(pairs(vs_pc)) %>% mutate(across(-1, ~ round(., 3)))

## Arousal
as_model <- lmer(scoreA ~ sem*pros + (1 + sem + pros|ID) + (1|stimuli), data = data_statement)
as_anova <- data.frame(anova(as_model, type = 2)) %>% mutate(across(-6, ~ round(., 3)))
r2_values <- r2(as_model)
marginal_R2    <- r2_values$R2_marginal
conditional_R2 <- r2_values$R2_conditional

as_pc      <- emmeans(as_model, ~ pros)
as_results <- data.frame(pairs(as_pc)) %>% mutate(across(-1, ~ round(., 3)))

# ============================================================
# PREPARE DATA FOR PLOTS
# ============================================================

ps1 <- ps1 %>% filter(!is.na(cond))
ps1$pros <- fct_recode(ps1$pros,
                       "Negative" = "n",
                       "Monotone" = "m",
                       "Positive" = "p")

x_labels <- c("Negative Context", "Positive Context", "Pseudo Context",
              "Negative Statement", "Positive Statement", "Pseudo Statement")

# ============================================================
# HELPER FUNCTION FOR PLOTS
# ============================================================

make_plot_s1 <- function(data, y_var, title, y_label) {
  ggplot(data, aes(x = cond, y = .data[[y_var]], fill = pros)) +
    geom_boxplot() +
    scale_fill_manual(values = okabe_pros, name = "Prosody") +
    theme(
      axis.text.x      = element_text(angle = 65, vjust = 0.6),
      panel.background = element_rect(fill = "white"),
      plot.title       = element_text(hjust = 0.5)
    ) +
    labs(title = title, x = "Conditions", y = y_label) +
    scale_x_discrete(labels = x_labels)
}

# ============================================================
# BUILD PLOTS
# ============================================================

p_valence <- make_plot_s1(ps1, "scoreV",
                          "Evaluation of the prosody valence of the stimuli",
                          "Score valence")

p_arousal <- make_plot_s1(ps1, "scoreA",
                          "Evaluation of the emotional intensity of the stimuli",
                          "Score arousal")

# ============================================================
# SAVE INDIVIDUAL PLOTS
# ============================================================

ggsave("plots/valence.png", plot = p_valence, width = 7, height = 7, dpi = 300)
ggsave("plots/arousal.png", plot = p_arousal, width = 7, height = 7, dpi = 300)

# ============================================================
# COMBINED FIGURE
# ============================================================

shared_legend <- get_legend(
  p_valence + theme(legend.position = "bottom",
                    legend.title = element_text(size = 10),
                    legend.text  = element_text(size = 9))
)

fig_combined <- ggarrange(
  p_valence + theme(legend.position = "none") + labs(tag = "A"),
  p_arousal + theme(legend.position = "none") + labs(tag = "B"),
  shared_legend,
  nrow = 3,
  heights = c(1, 1, 0.1)
)

ggsave("plots/pilot.png", width = 7, height = 7, dpi = 300)