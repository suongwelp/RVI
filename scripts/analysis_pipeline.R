################################################################################
# Analysis Pipeline for Rubber Voice Illusion Study
# Author: Suong Welp
# Date: July 21, 2025
# 
# This script analyzes behavioral ratings, fundamental frequency (F0) changes,
# and N1 ERP responses in a voice feedback manipulation experiment
################################################################################

# Setup -------------------------------------------------------------------
rm(list = ls())

# Load required packages
packages <- c("dplyr", "tidyverse", "ggplot2", "ggsignif", "rstatix", 
              "ggpubr", "ez", "emmeans", "broom", "viridis", "lme4", 
              "lmerTest", "ggpattern", "gridExtra","grid")

invisible(lapply(packages, library, character.only = TRUE))

# Set working directory and load data
#setwd(getwd())
setwd("C:/phD/nature_manuscript/Supplementary/data")

# Load datasets
demo <- read.csv("Dat0_Demo_SPQ_PDI.csv", header = TRUE)
questionnaire <- read.csv("Dat1_illusion_questionnaire.csv", header = TRUE)
F0 <- read.csv("Dat2_F0.csv", header = TRUE)
N1raw <- read.csv("Dat3_N1.csv", header = TRUE)

# Helper Functions --------------------------------------------------------

# Geometric mean function
geo_mean <- function(x) exp(mean(log(x), na.rm = TRUE))

# Create group classifications
create_groups <- function(data) {data %>% mutate(
      sample_SPQ = ifelse(SPQ >= 22, 'high-schizotypy', 'control'),
      sample_PDI = ifelse(PDI >= 5, 'high-schizotypy', 'control'))}

# Format p-values for plotting
format_pvalue <- function(p) {case_when(
    p < 0.001 ~ "***",p < 0.01  ~ "**", 
    p < 0.05  ~ "*",TRUE      ~ "n.s.")}

# Data Preprocessing ------------------------------------------------------
# Prepare questionnaire data
questionnaire <- questionnaire %>%
  pivot_longer(cols = SA_i_speak:SA_i_control, 
               names_to = "question", 
               values_to = "value") %>%
  mutate(value = as.numeric(value))

# Calculate condition means for questionnaire
conditionMeanQ <- questionnaire %>%
  group_by(sub, condition, question) %>%
  summarise(value = mean(value), .groups = "drop")

# Define samples based on cutoffs
demo <- create_groups(demo)

#===============================================================================
# QUESTIONNAIRE ANALYSIS
#===============================================================================

# Table 1: RVI experience (SO/SA ratings > 0) ----
rvi_table <- conditionMeanQ %>%
  # Create RVI Yes/No based on value > 0
  mutate(
    RVI = ifelse(value > 0, "Yes", "No"),
    # Create condition labels to match the table
    condition_label = case_when(
      condition == "veridical" ~ "Veridical",
      condition == "match" ~ "Stranger-match", 
      condition == "mismatch" ~ "Stranger-mismatch",
      TRUE ~ condition)) %>%
  # Count unique participants by condition, question, and RVI
  group_by(condition_label, question, RVI) %>%
  summarise(count = n_distinct(sub), .groups = "drop") %>%
  # Reshape to wide format
  pivot_wider(names_from = RVI, values_from = count, values_fill = 0) %>%
  # Arrange to match the table order
  arrange(factor(condition_label, levels = c("Veridical", "Stranger-match", "Stranger-mismatch")),
          question)

# Print the table
print(rvi_table)

# Figure 1: SO and SA ratings by condition ----
medians <- aggregate(value ~  condition+question, conditionMeanQ, median)
medians$value <- round(medians$value,digits = 2)

p1 <- ggplot(conditionMeanQ, aes(x = question, y = value, fill = question)) +
  geom_boxplot(aes(linetype = question), show.legend = FALSE, alpha = 0.7, size = 1,linewidth = 0.7) +
  scale_fill_viridis_d()+
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash")) +  # Different line types
  facet_wrap(~ condition, ncol = 1,labeller = labeller(condition = c("match" = "stranger-match","mismatch" = "stranger-mismatch"))) + 
  theme(panel.spacing = unit(1, "lines"),legend.position = "none") +
  #ggtitle("Sense of Ownership & Sense of Agency from all participants") +
  scale_y_continuous(breaks = c(-3, 0, 3), labels = c("disagree", "indifferent", "agree")) + 
  coord_flip() +   labs(y = "rating", x = "item") +
  stat_summary(fun = median, geom = "point", shape = 1, size = 1) +
  geom_text(data = medians, aes(label = value, y = value + 0.5))


# Table 2: Ratings by block and order ----
means <- aggregate(value ~  condition+question, conditionMeanQ, mean)
means$value <- round(means$value,digits = 2)

# Block analysis
block <- questionnaire %>% 
  dplyr::group_by(condition,question,block) %>% 
  summarise(median= round(median(value),2),mean = round(mean(value),2))

questionnaire <- questionnaire %>%mutate(across(c(sub, block, condition, question), as.factor))

block_anova <- ezANOVA(data = questionnaire,dv = value,wid = sub,
    within = .(block),between = .(condition, question), detailed = TRUE)
print(block_anova)

# Order analysis
order <- questionnaire %>%
  filter(order != 0) %>%
  group_by(sub, condition, question, order) %>%
  summarise(
    median = round(median(value), 2),
    mean = round(mean(value), 2),
    .groups = "drop") %>%
  mutate(across(c(sub, order, condition, question), as.factor))


# Random intercepts for subjects; add random slopes if convergence allows
m <- lmer(mean ~ order * condition * question + (1|sub), data = order)
anova(m, type = 3)   # omnibus F-tests
emm <- emmeans(m, ~ order | condition * question)
# Compare order 1 vs 2 within each condition×question
pairs(emm, adjust = "holm")

## ANOVA all participants with MEAN
# RM-ANOVA with condition and question as within-subject factors
rm_anova <- ezANOVA(
  data = conditionMeanQ,
  dv = value,
  wid = sub,  # participant ID
  within = .(condition, question),  # both are within-subject factors
  detailed = TRUE)
print(rm_anova)

# Post-hoc tests by condition and sample
pairwise.t.test(conditionMeanQ$value, conditionMeanQ$condition, 
                paired = TRUE, p.adjust.method = "bonferroni")
# For questions
run_pairwise <- function(data) {pairwise.t.test(data$value, data$condition, 
                  paired = TRUE, p.adjust.method = "bonferroni")}

# Apply to each question
results_by_question <- conditionMeanQ %>%
  group_by(question) %>%  group_modify(~ {result <- run_pairwise(.x)
    cat("\nQuestion:", unique(.x$question), "\n")
    print(result)
    data.frame(question = unique(.x$question))})

# descriptive table Appendix 3 - table 8
conditionMeanQ <- merge(conditionMeanQ,demo,by="sub")
meansampleSPQ <- conditionMeanQ %>% 
  dplyr::group_by(condition,question,sample_SPQ) %>% 
  summarise(median= round(median(value),2),mean = round(mean(value),2))

meansamplePDI <- conditionMeanQ %>% 
  dplyr::group_by(condition,question,sample_PDI) %>% 
  summarise(median= round(median(value),2),mean = round(mean(value),2))

# Group comparisons (SPQ and PDI)
ques_SPQ_anova1 <- ezANOVA(
  data = conditionMeanQ,
  dv = value,
  wid = sub,
  within = .(condition, question),
  between = sample_SPQ,  # sample as between-subjects factor
  detailed = TRUE)
print(ques_SPQ_anova1)

ques_PDI_anova1 <- ezANOVA(
  data = conditionMeanQ,
  dv = value,
  wid = sub,
  within = .(condition, question),
  between = sample_PDI,  # sample as between-subjects factor
  detailed = TRUE)
print(ques_PDI_anova1)

# Posthoc for each question in each sample (SPQ)
questions <- unique(conditionMeanQ$question)
samples <- unique(conditionMeanQ$sample_SPQ)
for(quest in questions) {
  for(samp in samples) {
    cat("\n=== Question:", quest, "| Sample:", samp, "===\n")
    subset_data <- conditionMeanQ[conditionMeanQ$question == quest & 
                                    conditionMeanQ$sample_SPQ == samp, ]
    if(nrow(subset_data) > 0) {
      result <- pairwise.t.test(subset_data$value, subset_data$condition, 
                                paired = TRUE, p.adjust.method = "bonferroni")
      print(result)}cat("\n")}}

# Posthoc for each question in each sample (PDI)
samples <- unique(conditionMeanQ$sample_PDI)
for(quest in questions) {
  for(samp in samples) {
    cat("\n=== Question:", quest, "| Sample:", samp, "===\n")
    subset_data <- conditionMeanQ[conditionMeanQ$question == quest & 
                                    conditionMeanQ$sample_PDI == samp, ]
    if(nrow(subset_data) > 0) {
      result <- pairwise.t.test(subset_data$value, subset_data$condition, 
                                paired = TRUE, p.adjust.method = "bonferroni")
      print(result)} 
    cat("\n")}
  }

#===============================================================================
# F0 ANALYSIS
#===============================================================================

# Prepare F0 data with geometric means
F0_all <- merge(F0,demo,by="sub")
geo_data <- F0_all %>%
  group_by(sub, condition, time) %>%
  summarise(gm_F0 = geo_mean(Mean.Pitch), .groups = "drop")
geo_wide <- geo_data %>%
  pivot_wider(names_from = time, values_from = gm_F0)
geo_wide <- geo_wide %>%
  mutate(semitone_shift = 12 * log2(post / pre))

# Table 3: F0 descriptives by condition and sample ----
condition_mean_F0 <- geo_wide %>%
  group_by(condition) %>%
  summarise(mean_pre_F0 = round(mean(pre, na.rm = TRUE),2),
    sd_pre_F0 = round(sd(pre, na.rm = TRUE),2),
    mean_post_F0 = round(mean(post, na.rm = TRUE),2),
    sd_post_F0 = round(sd(post, na.rm = TRUE),2),
    mean_semitone_shift = round(mean(semitone_shift, na.rm = TRUE),2),
    sd_semitone_shift = round(sd(semitone_shift, na.rm = TRUE),2),
    n_participants = n(),.groups = "drop")

# per sample: PDI
sample_means_F0_PDI <- F0_all %>%
  group_by(sub, condition, time, sample_PDI) %>%  # add 'sample' if it exists
  summarise(gm_F0 = geo_mean(Mean.Pitch), .groups = "drop") %>%
  pivot_wider(names_from = time, values_from = gm_F0) %>%
  mutate(semitone_shift = 12 * log2(post / pre)) %>%
  group_by(condition, sample_PDI) %>%
  summarise(mean_pre_F0 = round(mean(pre, na.rm = TRUE),2),
    sd_pre_F0 = round(sd(pre, na.rm = TRUE),2),
    mean_post_F0 = round(mean(post, na.rm = TRUE),2),
    sd_post_F0 = round(sd(post, na.rm = TRUE),2),
    mean_semitone_shift = round(mean(semitone_shift, na.rm = TRUE),2),
    sd_semitone_shift = round(sd(semitone_shift, na.rm = TRUE),2),
    n_participants = n(),.groups = "drop")

# per sample: SPQ
sample_means_F0_SPQ <- F0_all %>%
  group_by(sub, condition, time, sample_SPQ) %>%  # add 'sample' if it exists
  summarise(gm_F0 = geo_mean(Mean.Pitch), .groups = "drop") %>%
  pivot_wider(names_from = time, values_from = gm_F0) %>%
  mutate(semitone_shift = 12 * log2(post / pre)) %>%
  group_by(condition, sample_SPQ) %>%
  summarise(mean_pre_F0 = round(mean(pre, na.rm = TRUE),2),
    sd_pre_F0 = round(sd(pre, na.rm = TRUE),2),
    mean_post_F0 = round(mean(post, na.rm = TRUE),2),
    sd_post_F0 = round(sd(post, na.rm = TRUE),2),
    mean_semitone_shift = round(mean(semitone_shift, na.rm = TRUE),2),
    sd_semitone_shift = round(sd(semitone_shift, na.rm = TRUE),2),
    n_participants = n(),.groups = "drop")


geo_wide <- merge(geo_wide,demo,by="sub")

# anova SPQ group differences
rm_anova_F0_SPQ <- ezANOVA(
  data = geo_wide,
  dv = semitone_shift,
  wid = sub,
  within = condition,      # within-subjects factor
  between = sample_SPQ,        # between-subjects factor
  detailed = TRUE)
print(rm_anova_F0_SPQ)

# anova PDI group differences
rm_anova_F0_PDI <- ezANOVA(
  data = geo_wide,
  dv = semitone_shift,
  wid = sub,
  within = condition,      # within-subjects factor
  between = sample_PDI,        # between-subjects factor
  detailed = TRUE)
print(rm_anova_F0_PDI)

# Ensure time is ordered: pre before post
geo_long <- geo_wide %>%
  select(sub, condition, pre, post) %>%
  pivot_longer(cols = c(pre, post), names_to = "time", values_to = "F0") %>%
  mutate(time = factor(time, levels = c("pre", "post")))

# Figure 2: F0 pre vs post ----
# Paired t-tests and formatted p-values
sig_data <- geo_wide %>%
  group_by(condition) %>%
  summarise(p = t.test(pre, post, paired = TRUE)$p.value, .groups = "drop") %>%
  mutate(p_adj = p.adjust(p, method = "bonferroni"),
         label = format_pvalue(p_adj))

# Plot
p2 <- ggplot(geo_long, aes(x = time, y = F0, group = sub, color = condition, shape = condition)) +
  geom_line(alpha = 0.4) +
  geom_point(size = 2) +
  geom_hline(yintercept = 147, linetype = "dashed", color = "#440154", alpha = 0.7) +
  geom_hline(yintercept = 294, linetype = "dashed", color = "#35b779", alpha = 0.7) +
  scale_color_manual(values = c("match" = "#440154", "mismatch" = "#31688e")) +
  scale_shape_manual(values = c("match" = 16, "mismatch" = 17)) +
  facet_wrap(~condition, labeller = labeller(condition = c("match" = "stranger-match", "mismatch" = "stranger-mismatch"))) +
  labs(x = "Time", y = "Geometric Mean F0 (Hz)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  geom_text(data = sig_data, aes(x = 1.5, y = max(geo_long$F0, na.rm = TRUE) * 1.05, label = label),
            inherit.aes = FALSE, size = 6)
# Print the base plot
grid.newpage()
g <- ggplotGrob(p2)
grid.draw(g)

# Add annotation between facets
grid.text("female_stranger: 294 Hz",x = 0.545, y = 0.825,  
  gp = gpar(col = "#35b779", fontsize = 9, fontface = "italic"))
grid.text("male_stranger: 147 Hz",x = 0.55, y = 0.355, 
  gp = gpar(col = "#440154", fontsize = 9, fontface = "italic"))


# Figure 3: F0-Questionnaire correlations ----
F0_Qdat <- merge(geo_wide,conditionMeanQ[,c(1:4)],by = c("sub","condition"))

# SPQ cutoff correlations
cor_stats <- F0_Qdat %>% group_by(question, sample_SPQ) %>%
  summarise(
    r = cor(value, semitone_shift, use = "complete.obs"),
    p = cor.test(value, semitone_shift)$p.value,
    .groups = "drop") %>%
  mutate(label = paste0(sample_SPQ, ": r = ", round(r, 2), " ", format_pvalue(p)),
    y_pos = ifelse(sample_SPQ == "control", -3.5, -4.0),x_pos = -Inf)

p3a <- ggplot(F0_Qdat, aes(x = value, y = semitone_shift, color = sample_SPQ,
                        shape = sample_SPQ, linetype = sample_SPQ)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, size = 1, alpha = 0.2) +
  facet_wrap(~question, scales = "free_x") +
  geom_text(data = cor_stats, aes(x = x_pos, y = y_pos, label = label),
            inherit.aes = FALSE,hjust = -0.05, size = 3.5, color = "black") +
  scale_color_manual(values = c("control" = "#5aae61", "high-schizotypy" = "#440154"), name = "Group") +
  scale_shape_manual(values = c("control" = 16, "high-schizotypy" = 17), name = "Group") +
  scale_linetype_manual(values = c("control" = "solid", "high-schizotypy" = "dashed"), name = "Group") +
  labs(title = "a) SPQ Cutoff", x = "Questionnaire Rating", y = "Semitone Shift") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(linetype = c("solid", "dashed"), shape = c(16, 17))),
         shape = "none", linetype = "none")

# PDI cutoff correlations
cor_stats <- F0_Qdat %>%
  group_by(question, sample_PDI) %>%
  summarise(
    r = cor(value, semitone_shift, use = "complete.obs"),
    p = cor.test(value, semitone_shift)$p.value,
    .groups = "drop") %>%
  mutate(
    label = paste0(sample_PDI, ": r = ", round(r, 2), " ", format_pvalue(p)),
    y_pos = ifelse(sample_PDI == "control", -3.5, -4.0),x_pos = -Inf)

p3b <- ggplot(F0_Qdat, aes(x = value, y = semitone_shift, color = sample_PDI,
                        shape = sample_PDI, linetype = sample_PDI)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, size = 1, alpha = 0.2) +
  facet_wrap(~question, scales = "free_x") +
  geom_text(data = cor_stats,
            aes(x = x_pos, y = y_pos, label = label),
            inherit.aes = FALSE,
            hjust = -0.05, size = 3.5, color = "black") +
  scale_color_manual(values = c("control" = "#5aae61", "high-schizotypy" = "#440154"), name = "Group") +
  scale_shape_manual(values = c("control" = 16, "high-schizotypy" = 17), name = "Group") +
  scale_linetype_manual(values = c("control" = "solid", "high-schizotypy" = "dashed"), name = "Group") +
  labs(title = "b) PDI Cutoff", x = "Questionnaire Rating", y = "Semitone Shift") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(linetype = c("solid", "dashed"), shape = c(16, 17))),
         shape = "none", linetype = "none")

# Combine plots
gridExtra::grid.arrange(p3a, grid::nullGrob(), p3b, ncol = 3, widths = c(1, 0.1, 1))

#===============================================================================
# N1 AMPLITUDE ANALYSIS
#===============================================================================

# Prepare N1 data
N1 <- N1raw %>%
  pivot_longer(
    cols = c(veridical_speak, stranger_match, stranger_listen, stranger_mis, veridical_listen),
    names_to = "condition",
    values_to = "amp")
N1$condition <- factor(N1$condition , levels=c('veridical_speak','veridical_listen', 'stranger_match', 
                                                     'stranger_mis', 'stranger_listen' ))
N1_all <- merge(N1, demo, by = "sub")

# Table 4: N1 amplitude descriptives ----
table4_all <- N1_all %>%
  group_by(condition) %>%
  summarise(mean_amp = round(mean(amp), 2), sd_amp = round(sd(amp), 2), .groups = "drop")

table4_SPQ <- N1_all %>%
  group_by(condition, sample_SPQ) %>%
  summarise(mean_amp = round(mean(amp), 2), sd_amp = round(sd(amp), 2), .groups = "drop")

table4_PDI <- N1_all %>%
  group_by(condition, sample_PDI) %>%
  summarise(mean_amp = round(mean(amp), 2), sd_amp = round(sd(amp), 2), .groups = "drop")

# Figure 5: N1 amplitude by condition ----
means <- aggregate(amp ~ condition, N1_all, mean)
means$amp <- round(means$amp, digits = 2)

p5 <- ggplot(N1_all, aes(x=condition, y=amp, group=condition, fill=condition)) + 
  geom_boxplot(alpha = 0.7) +
  stat_compare_means(comparisons = list(
    c('stranger_mis', 'stranger_listen'),
    c('stranger_match', 'stranger_listen'), 
    c('veridical_speak', 'veridical_listen'),
    c('stranger_match', 'stranger_mis')),
  method = "t.test",
  p.adjust.method = "bonferroni",
  label = "p.signif",
  step.increase = 0.1,size=5) +
  scale_fill_viridis_d() +
  stat_summary(fun = mean, geom = "point", shape = 1, size = 1) +
  geom_text(data = means, aes(label = amp, y = amp + 1)) +
  labs(title = 'b) N1 amplitudes in different conditions for all participants', y = "Amplitude [uV]") + 
  theme_minimal() + theme(legend.position = "none",
    plot.title = element_text(size = 12))

# # N1 ANOVAs for main condition
rm_anova_n1_all <- ezANOVA(
  data = N1_all,
  dv = amp,
  wid = sub,
  within = condition,
  detailed = TRUE)
print(rm_anova_n1_all)

# Pairwise t-tests with Bonferroni correction
pairwise.t.test(N1_all$amp, N1_all$condition, 
                paired = TRUE, p.adjust.method = "bonferroni")


# anova for each group
N1_all$sample <- ifelse(N1_all$PDI >= 5, 'high-schizotypy', 'control')
rm_anova_n1_group <- ezANOVA(
  data = N1_all,
  dv = amp,
  wid = sub,
  within = condition,
  between = sample,  # if you want to include sample as between-subjects factor
  detailed = TRUE)
print(rm_anova_n1_group)

# Post-hoc tests
sample_diff <- N1_all %>%
  group_by(condition) %>%
  summarise(p_raw = t.test(amp[sample == "control"], 
                   amp[sample == "high-schizotypy"])$p.value,
    control_mean = mean(amp[sample == "control"], na.rm = TRUE),
    schizo_mean = mean(amp[sample == "high-schizotypy"], na.rm = TRUE),
    .groups = "drop") %>%
  mutate(p_adjusted = p.adjust(p_raw, method = "bonferroni"))

print(sample_diff)

# Figure 6: N1 amplitude by group ----
# Group comparison statistics
sample_diff <- N1_all %>%
  group_by(condition) %>%
  summarise(
    p_raw = t.test(amp[sample_PDI == "control"], amp[sample_PDI == "high-schizotypy"])$p.value,
    .groups = "drop") %>%
  mutate(p_adjusted = p.adjust(p_raw, method = "bonferroni"))

# for SPQ cutoff
meanSPQ <- aggregate(amp ~ condition + sample_SPQ, N1_all, mean)
meanSPQ$amp <- round(meanSPQ$amp, digits = 2)

# Calculate both raw and corrected p-values
p_stats <- N1_all %>%  group_by(condition) %>%
  summarise(p_raw = t.test(amp[sample_SPQ == "control"], 
                           amp[sample_SPQ == "high-schizotypy"])$p.value,.groups = "drop") %>%
  mutate(p_adj = p.adjust(p_raw, method = "bonferroni"),
         # Add stars for significant p-values
         p_raw_formatted = paste0("p=", sub("^0+", "", round(p_raw, 3)), ifelse(p_raw < 0.05, "*", "")),
         p_adj_formatted = paste0("p.adj=", sub("^0+", "", round(p_adj, 3)),ifelse(p_adj < 0.05, "*", "")),
         label = paste0(p_raw_formatted, "\n", p_adj_formatted))

# Plot Figure 6 A)
p6a <- ggplot(N1_all, aes(x = sample_SPQ, y = amp, fill = condition)) +
  geom_boxplot(color = alpha("black", 0.3)) +
  scale_fill_viridis_d(alpha = 0.5) +
  facet_wrap(~ factor(condition, levels = c('veridical_speak','veridical_listen', 'stranger_match', 'stranger_mis', 'stranger_listen')),nrow = 1)+
  ggtitle("A) SPQ Cutoff") +
  stat_summary(fun = mean, geom = "point", shape = 1, size = 1) +
  geom_text(data = meanSPQ, aes(label = amp, y = amp + 1)) +
  # Add manual significance lines
  geom_segment(aes(x = 1, y = 5, xend = 2, yend = 5), 
               data = data.frame(condition = levels(N1_all$condition))) +
  # Add both raw and adjusted p-values with larger text
  geom_text(data = p_stats, aes(x = 1.5, y = 7.5, label = label), 
            inherit.aes = FALSE, size = 3.5) +  ylim(-20, 10) + 
  labs(y = "N1 Amplitude [μV]", x = NULL) +  
  theme_minimal() +
  theme(legend.position = "none")

# for PDI cutoff
meanPDI <- aggregate(amp ~ condition + sample_PDI, N1_all, mean)
meanPDI$amp <- round(meanPDI$amp, digits = 2)
# Calculate both raw and corrected p-values
p_stats <- N1_all %>%  group_by(condition) %>%
  summarise(p_raw = t.test(amp[sample_PDI == "control"], 
                           amp[sample_PDI == "high-schizotypy"])$p.value,.groups = "drop") %>%
  mutate(p_adj = p.adjust(p_raw, method = "bonferroni"),
         # Add stars for significant p-values
         p_raw_formatted = paste0("p=", sub("^0+", "", round(p_raw, 3)), ifelse(p_raw < 0.05, "*", "")),
         p_adj_formatted = paste0("p.adj=", sub("^0+", "", round(p_adj, 3)),ifelse(p_adj < 0.05, "*", "")),
         label = paste0(p_raw_formatted, "\n", p_adj_formatted))

# Plot Figure 6 B)
p6b <- ggplot(N1_all, aes(x = sample_PDI, y = amp, fill = condition)) +
  geom_boxplot(color = alpha("black", 0.3)) +
  scale_fill_viridis_d(alpha = 0.5) +
  facet_wrap(~ factor(condition, levels = c('veridical_speak','veridical_listen', 'stranger_match', 'stranger_mis', 'stranger_listen')), 
             nrow = 1)+
  ggtitle("B) PDI Cutoff") +
  stat_summary(fun = mean, geom = "point", shape = 1, size = 1) +
  geom_text(data = meanPDI, aes(label = amp, y = amp + 1)) +
  # Add manual significance lines
  geom_segment(aes(x = 1, y = 5, xend = 2, yend = 5), data = data.frame(condition = levels(N1_all$condition))) +
  # Add both raw and adjusted p-values with larger text
  geom_text(data = p_stats, aes(x = 1.5, y = 7.5, label = label), inherit.aes = FALSE, size = 3.5) +  
  ylim(-20, 10) + labs(y = "N1 Amplitude [μV]", x = NULL) +  
  theme_minimal() + theme(legend.position = "none")

#===============================================================================
# N1 SUPPRESSION ANALYSIS
#===============================================================================

# N1 SIS in absolute value, uV
N1SIS_uV <- N1raw %>% mutate(
    veridical = veridical_speak - veridical_listen,
    match = stranger_match - stranger_listen,
    mismatch = stranger_mis - stranger_listen) %>%  
   select(sub, veridical, match, mismatch) %>%
   pivot_longer(cols = c(veridical, match, mismatch),
    names_to = "condition",
    values_to = "amp_uV")

N1SIS_p <- N1raw %>%mutate(
    veridical = 1-veridical_speak/veridical_listen,
    match = 1-stranger_match/stranger_listen,
    mismatch = 1-stranger_mis/stranger_listen) %>%
  select(sub, veridical, match, mismatch) %>%
  pivot_longer(
    cols = c(veridical, match, mismatch),
    names_to = "condition",
    values_to = "amp_p")

N1SIS <- merge(N1SIS_uV,N1SIS_p , by = c("sub","condition"))
N1SIS <- merge(N1SIS,demo , by = "sub")
N1SIS$condition <- factor(N1SIS$condition,
                          levels = c("veridical","match", "mismatch" ),ordered = TRUE)

# Table 5: N1 suppression descriptives ----
table5_all <- N1SIS %>% 
  dplyr::group_by(condition) %>% 
  summarise(mean_uV = round(mean(amp_uV),2),sd_uV = round(sd(amp_uV),2),
            mean_p = round(mean(amp_p)*100,2),sd_p = round(sd(amp_p)*100,2))

table5_SPQ <- N1SIS %>% 
  dplyr::group_by(condition,sample_SPQ) %>% 
  summarise(mean_uV = round(mean(amp_uV),2),sd_uV = round(sd(amp_uV),2),
            mean_p = round(mean(amp_p)*100,2),sd_p = round(sd(amp_p)*100,2))
table5_PDI <- N1SIS %>% 
  dplyr::group_by(condition,sample_PDI) %>% 
  summarise(mean_uV = round(mean(amp_uV),2),sd_uV = round(sd(amp_uV),2),
            mean_p = round(mean(amp_p)*100,2),sd_p = round(sd(amp_p)*100,2))

# N1 suppression ANOVAs
# SPQ cutoff
anova_n1_SPQ <- ezANOVA(
  data = N1SIS,
  dv = amp_p,
  wid = sub,
  within = condition,
  between = sample_SPQ,  
  detailed = TRUE)
print(anova_n1_SPQ)

# PDI cutoff
anova_n1_PDI <- ezANOVA(
  data = N1SIS,
  dv = amp_p,
  wid = sub,
  within = condition,
  between = sample_PDI,  
  detailed = TRUE)
print(anova_n1_PDI)

# Figure 7: N1 suppression by group ----
# figure 7 A) SPQ
stat.test <- N1SIS %>%
  group_by(condition) %>%
  wilcox_test(amp_p ~ sample_SPQ) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

print(stat.test)

# Add p-values onto the bar plots
stat.test <- stat.test %>%
  add_xy_position(fun = "mean", x = "condition", dodge = 0.8) 
stat.test$y.position <- 75

p7_SPQ <- ggbarplot(table5_SPQ, x = "condition", y = "mean_p", 
          fill = "sample_SPQ",position = position_dodge(0.8),
          palette = c("control" = "#21908C", "high-schizotypy" = "#5DC863"),
          legend = "bottom") +  ggtitle("a1) SPQ Cutoff") +
  stat_summary(fun = mean, geom = "point", shape = 1, size = 1) +
  geom_text(data = table5_SPQ, aes(x = condition, y = mean_p + 0.05, label = mean_p, group = sample_SPQ),
            position = position_dodge(width = 1), vjust = -0.5, size = 4) +
  labs(x = "condition", y = "N1 Suppression in %",fill= "Group") +
  scale_x_discrete(
    limits = c("veridical", "match", "mismatch"),
    labels = c("veridical", "stranger-match", "stranger-mismatch")) +
  stat_pvalue_manual(stat.test,label = "p.adj.signif",tip.length = 0,bracket.nudge.y = 1) +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "bottom")

# figure 7 B) PDI
stat.test <- N1SIS %>%
  group_by(condition) %>%
  wilcox_test(amp_p ~ sample_PDI) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

print(stat.test)

# Add p-values onto the bar plots
stat.test <- stat.test %>%
  add_xy_position(fun = "mean", x = "condition", dodge = 0.8) 
stat.test$y.position <- 75

p7_PDI <- ggbarplot(table5_PDI, x = "condition", y = "mean_p", 
  fill = "sample_PDI", 
  position = position_dodge(0.8),
  palette = c("control" = "#21908C", "high-schizotypy" = "#5DC863"),
  legend = "bottom") +  ggtitle("a2) PDI Cutoff") +
  stat_summary(fun = mean, geom = "point", shape = 1, size = 1) +
  geom_text(data = table5_PDI, aes(x = condition, y = mean_p + 0.05, label = mean_p, group = sample_PDI),
            position = position_dodge(width = 1), vjust = -0.5, size = 4) +
  labs(x = "condition", y = "N1 Suppression in %",fill="Group") +
  scale_x_discrete(
    limits = c("veridical", "match", "mismatch"),
    labels = c("veridical", "stranger-match", "stranger-mismatch")) +
  stat_pvalue_manual(stat.test,label = "p.adj.signif",tip.length = 0,bracket.nudge.y = 1) +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(face = "bold"),
    legend.position = "bottom")


#===============================================================================
# CORRELATION AND REGRESSION ANALYSIS
#===============================================================================

SOSA <- conditionMeanQ %>%
  pivot_wider(names_from = question, values_from = value) %>%
  mutate(
    ownership = rowMeans(select(., SO_my_voice, SO_my_voice_modified), na.rm = TRUE),
    agency = rowMeans(select(., SA_i_control, SA_i_speak), na.rm = TRUE))

all_dat <- merge(N1SIS,geo_wide[,c("sub","condition","semitone_shift")],by=c("sub","condition"))
all_dat <- merge(all_dat,SOSA[,c("sub","condition","ownership","agency")],by=c("sub","condition"))
match_data <- subset(all_dat, condition == "match")
match_data <- match_data %>%rename(N1 = amp_p,F0shift = semitone_shift)

control_model_PDI <- lm(N1 ~ F0shift, data = subset(match_data, sample_PDI == "control"))
high_model_PDI <- lm(N1 ~ F0shift, data = subset(match_data, sample_PDI == "high-schizotypy"))
summary(control_model_PDI)
summary(high_model_PDI)


control_model_SPQ <- lm(N1 ~ F0shift, data = subset(match_data, sample_SPQ == "control"))
high_model_SPQ <- lm(N1 ~ F0shift, data = subset(match_data, sample_SPQ == "high-schizotypy"))
summary(control_model_SPQ)
summary(high_model_SPQ)

# Combine model summaries (PDI)
beta_stats_PDI <- tibble(
  sample_PDI = c("control", "high-schizotypy"),
  model = list(control_model_PDI, high_model_PDI)) %>%
  mutate(coef_summary = map(model, ~ summary(.x)$coefficients),
    beta = map_dbl(coef_summary, ~ .x["F0shift", "Estimate"]),
    p = map_dbl(coef_summary, ~ .x["F0shift", "Pr(>|t|)"]),
    label = paste0(sample_PDI, ": β = ", round(beta, 2), " (", format_pvalue(p), ")"),
    y_pos = ifelse(sample_PDI == "control", 1.3, 1.2),x_pos = 0)

# SPQ
beta_stats_SPQ <- tibble(
  sample_SPQ = c("control", "high-schizotypy"),
  model = list(control_model_SPQ, high_model_SPQ)) %>%
  mutate(coef_summary = map(model, ~ summary(.x)$coefficients),
         beta = map_dbl(coef_summary, ~ .x["F0shift", "Estimate"]),
         p = map_dbl(coef_summary, ~ .x["F0shift", "Pr(>|t|)"]),
         label = paste0(sample_SPQ, ": β = ", round(beta, 2), " (", format_pvalue(p), ")"),
         y_pos = ifelse(sample_SPQ == "control", 1.3, 1.2),x_pos = 0)

# Figure 8: N1-F0 correlations by group ----
custom_colors <- c("control" = "#5aae61", "high-schizotypy" = "#440154")
custom_shapes <- c("control" = 16, "high-schizotypy" = 17)
custom_linetypes <- c("control" = "solid", "high-schizotypy" = "dashed")

# Fig 8A) SPQ 
p8_SPQ <- ggplot(match_data, aes(x = F0shift, y = N1, color = sample_SPQ, shape = sample_SPQ, linetype = sample_SPQ)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1, alpha = 0.2) +
  geom_text(data = beta_stats_SPQ, aes(x = x_pos, y = y_pos, label = label),
            inherit.aes = FALSE, hjust = -0.05, size = 3.5, color = "black") +
  scale_color_manual(values = custom_colors, name = "Group") +
  scale_shape_manual(values = custom_shapes, name = "Group") +
  scale_linetype_manual(values = custom_linetypes, name = "Group") +
  labs(title = "b1) SPQ Cutoff", x = "F0 Shift (semitone)", y = "N1 Suppression (%)") +
  theme_minimal(base_size = 10) + theme(plot.title = element_text(hjust = 0), legend.position = "bottom")

# Fig 8B) PDI 
p8_PDI <- ggplot(match_data, aes(x = F0shift, y = N1, color = sample_PDI, shape = sample_PDI, linetype = sample_PDI)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1, alpha = 0.2) +
  geom_text(data = beta_stats_PDI, aes(x = x_pos, y = y_pos, label = label),
            inherit.aes = FALSE, hjust = -0.05, size = 3.5, color = "black") +
  scale_color_manual(values = custom_colors, name = "Group") +
  scale_shape_manual(values = custom_shapes, name = "Group") +
  scale_linetype_manual(values = custom_linetypes, name = "Group") +
  labs(title = "b2) PDI Cutoff", x = "F0 Shift (semitone)", y = "N1 Suppression (%)") +
  theme_minimal(base_size = 10) + theme(plot.title = element_text(hjust = 0), legend.position = "bottom")

# Combine plots
gridExtra::grid.arrange(
  textGrob("a) N1 Suppression Group Difference",  x = 0.1, hjust = 0,gp = gpar(fontsize = 12)),nullGrob(), nullGrob(),
  p7_SPQ, grid::nullGrob(), p7_PDI,
  textGrob("b) Correlation between N1 Suppression & F0 Shift in each Group",  x = 0.1, hjust = 0,gp = gpar(fontsize = 12)),nullGrob(), nullGrob(),
  p8_SPQ, grid::nullGrob(), p8_PDI,
  ncol = 3,widths = c(1, 0.1, 1),heights = c(0.2,1.5, 0.2, 1.5))

# Table 6: Hierarchical regression models ----

m1 <- lm(N1 ~ F0shift, data = match_data)

# PDI hierarchy
m2_PDI <- lm(N1 ~ F0shift + sample_PDI,data = match_data)
m3_PDI <- lm(N1 ~ F0shift * sample_PDI,data = match_data)
m4_PDI <- lm(N1 ~ F0shift * sample_PDI + ownership + agency,data = match_data)
m5_PDI <- lm(N1 ~ F0shift * sample_PDI + ownership * sample_PDI + agency * sample_PDI,data = match_data)

# SPQ hierarchy
m2_SPQ <- lm(N1 ~ F0shift + sample_SPQ, data = match_data)
m3_SPQ <- lm(N1 ~ F0shift * sample_SPQ, data = match_data)
m4_SPQ <- lm(N1 ~ F0shift * sample_SPQ + ownership + agency,data = match_data)
m5_SPQ <- lm(N1 ~ F0shift * sample_SPQ + ownership * sample_SPQ + agency * sample_SPQ,data = match_data)

make_full_tbl <- function(model_list, added_terms) {
  
  ## ΔDF, F, p from anova
  aov_df <- as.data.frame(do.call(anova, unname(model_list)))
  aov_tbl <- tibble(
    Model = names(model_list),
    ΔDF   = c(NA, aov_df$Df[-1]),
    F     = c(NA, round(aov_df$F[-1], 2)),
    p_raw = c(NA, aov_df$`Pr(>F)`[-1])) %>%
    mutate(Sig = case_when(
        is.na(p_raw)         ~ "",
        p_raw < .001         ~ "***",
        p_raw < .01          ~ "**",
        p_raw < .05          ~ "*",
        TRUE                 ~ ""),
      p = ifelse(is.na(p_raw), "", formatC(p_raw, format = "e", digits = 2))) %>%
    select(-p_raw)
  ## adj R², AIC
  fit_tbl <- imap_dfr(model_list, function(mod, nm) {g <- glance(mod)
    tibble(Model = nm,
           adj_R2 = round(g$adj.r.squared, 3),
           AIC    = round(AIC(mod), 2))}) %>%
    mutate(ΔAIC = AIC - min(AIC))
  ## merge & add Added-terms column
  aov_tbl %>%
    left_join(fit_tbl, by = "Model") %>%
    mutate(Added_terms = added_terms) %>%
    select(Model, Added_terms, ΔDF, F, p, Sig, adj_R2, AIC, ΔAIC)}

pdi_models <- list(
  m1 = m1,
  m2 = m2_PDI,
  m3 = m3_PDI,
  m4 = m4_PDI,
  m5 = m5_PDI
)
pdi_added <- c(
  "—",
  "+ PDI main effect",
  "+ F0 × PDI interaction",
  "+ SO & SA (additive)",
  "+ SO × PDI & SA × PDI"
)

spq_models <- list(
  m1 = m1,
  m2 = m2_SPQ,
  m3 = m3_SPQ,
  m4 = m4_SPQ,
  m5 = m5_SPQ
)
spq_added <- c(
  "—",
  "+ SPQ main effect",
  "+ F0 × SPQ interaction",
  "+ SO & SA (additive)",
  "+ SO × SPQ & SA × SPQ"
)

pdi_table <- make_full_tbl(pdi_models, pdi_added)
spq_table <- make_full_tbl(spq_models, spq_added)

cat("\n######  SPQ hierarchy  ######\n")
print(spq_table, n = Inf)

cat("\n######  PDI hierarchy  ######\n")
print(pdi_table, n = Inf)



