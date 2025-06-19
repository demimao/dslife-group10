
library(broom)
library(readxl)
library(dplyr)
library(tidyr)
library(tidyverse)
library(tibble)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggrepel)
# read data
mut_df <- read_excel("41586_2021_4278_MOESM4_ESM.xlsx", sheet = 2)
clin_df <- read_excel("41586_2021_4278_MOESM4_ESM.xlsx", sheet = 1)

# Count the number of mutations of each gene
gene_freq <- table(mut_df$Hugo_Symbol)

# Calculate the total number of mutations
total_mutations <- sum(gene_freq)

# Calculate the frequency percentage of each gene
gene_freq_percent <- (gene_freq / total_mutations) * 100

# Sort in descending order of frequency
gene_freq_percent_sorted <- sort(gene_freq_percent, decreasing = TRUE)

# Print the percentages of all gene frequencies
print(gene_freq_percent_sorted)

# Print the percentages of top 5 genes frequencies
top5_genes_percent <- head(gene_freq_percent_sorted, 5)
print(top5_genes_percent)
# Gene list
genes_of_interest <- c("TP53", "TTN", "PIK3CA", "SYNE1", "HMCN1")
# Extract whether each sample contains TP53 TTN PIK3CA SYNE1 HMCN1

mut_status <- mut_df %>%
  filter(Hugo_Symbol %in% c("TP53", "TTN", "PIK3CA", "SYNE1", "HMCN1")) %>%
  select(Donor.ID, Hugo_Symbol) %>%
  distinct() %>%
  mutate(has_mut = 1) %>%
  pivot_wider(names_from = Hugo_Symbol, values_from = has_mut, values_fill = 0)

# View the sample
head(mut_status)
# Check which gene sequences are missing
missing_genes <- setdiff(genes_of_interest, colnames(mut_status))
# Add a value of 0 to the missing gene sequence
for (gene in missing_genes) {
  mut_status[[gene]] <- 0
}
# Combine mutation information with clinical information
merged_df <- clin_df %>%
  left_join(mut_status, by = "Donor.ID") %>%
  mutate(
    TP53 = ifelse(is.na(TP53), 0, TP53),
    TTN = ifelse(is.na(TTN), 0, TTN),
    PIK3CA = ifelse(is.na(PIK3CA), 0, PIK3CA),
    SYNE1 = ifelse(is.na(SYNE1), 0, SYNE1),
    HMCN1 = ifelse(is.na(HMCN1), 0, HMCN1),
    pCR = ifelse(pCR.RD == "pCR", 1, 0)  # Convert to a 0/1 variable
  )
# TP53 and pCR
model_tp53 <- glm(pCR ~ TP53, data = merged_df, family = binomial)
summary(model_tp53)
exp(cbind(OR = coef(model_tp53), confint(model_tp53)))  # OR + CI

# TTN and pCR
model_ttn <- glm(pCR ~ TTN, data = merged_df, family = binomial)
summary(model_ttn)
exp(cbind(OR = coef(model_ttn), confint(model_ttn)))  # OR + CI

# PIK3CA and pCR
model_pik3ca <- glm(pCR ~ PIK3CA, data = merged_df, family = binomial)
summary(model_pik3ca)
exp(cbind(OR = coef(model_pik3ca), confint(model_pik3ca)))  # OR + CI

# SYNE1 and pCR
model_syne1 <- glm(pCR ~ SYNE1, data = merged_df, family = binomial)
summary(model_syne1)
exp(cbind(OR = coef(model_syne1), confint(model_syne1)))  # OR + CI

# HMCN1 and pCR
model_hmcn1 <- glm(pCR ~ HMCN1, data = merged_df, family = binomial)
summary(model_hmcn1)
exp(cbind(OR = coef(model_hmcn1), confint(model_hmcn1)))  # OR + CI

# Extract the model information function
extract_or_ci <- function(model, gene) {
  est <- summary(model)$coefficients[2, ]
  ci <- confint(model)[2, ]
  data.frame(
    Gene = gene,
    OR = exp(est["Estimate"]),
    CI_lower = exp(ci[1]),
    CI_upper = exp(ci[2]),
    p_value = est["Pr(>|z|)"]
  )
}

# Construct the statistical result table
results_df <- bind_rows(
  extract_or_ci(model_tp53, "TP53"),
  extract_or_ci(model_ttn, "TTN"),
  extract_or_ci(model_pik3ca, "PIK3CA"),
  extract_or_ci(model_syne1, "SYNE1"),
  extract_or_ci(model_hmcn1, "HMCN1")
)

# View
print(results_df)


# Drawing
# Add group colors (red for >1, blue for <1)
results_df <- results_df %>%
  mutate(
    label = paste0("OR = ", round(OR, 2), 
                   "\nCI: [", round(CI_lower, 2), ", ", round(CI_upper, 2), "]",
                   "\nP = ", signif(p_value, 2)),
    color = case_when(
      p_value > 0.05 ~ "nonsignificant",
      OR >= 1 ~ "risk",
      OR < 1 ~ "protective"
    )
  )


ggplot(results_df, aes(x = reorder(Gene, OR), y = OR)) +
  geom_point(aes(color = color), size = 4) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, color = color), width = 0.2, size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  coord_flip() +
  scale_color_manual(
    values = c(
      "risk" = "#E64B35FF",         # 红色
      "protective" = "#4DBBD5FF",   # 蓝色
      "nonsignificant" = "gray60"   # 灰色
    )
  ) +
  scale_y_log10() +
  geom_text(aes(label = label), hjust = -0.1, size = 4.2) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +
  labs(
    y = "Odds Ratio (log scale)",
    title = "Forest Plot: Gene Mutations and pCR Response"
  )




# Mutation number per sample (can be used as a simplified version of TMB)
tmb_df <- mut_df %>%
  group_by(Donor.ID) %>%
  summarise(TMB = n())  # If each variation is one locus, then n() is approximately the total number of mutations

head(tmb_df)
# Combine TMB with clinical data
merged_df <- merged_df %>%
  mutate(
    HER2.status = toupper(HER2.status),  # 全部转大写
    HER2_group = case_when(
      HER2.status == "POS" ~ "HER2+",
      HER2.status == "NEG" ~ "HER2-",
      TRUE ~ NA_character_
    )
  )
head(merged_df)
# Overall comparison of TMB (median + test)
library(ggpubr)

# View the median
merged_df %>%
  group_by(pCR_binary) %>%
  summarise(median_TMB = median(TMB, na.rm = TRUE))

# Wilcoxon test
wilcox.test(TMB ~ pCR_binary, data = merged_df)
# HER2- group
wilcox.test(TMB ~ pCR_binary, data = merged_df %>% filter(HER2_group == "HER2-"))

# HER2+ group
wilcox.test(TMB ~ pCR_binary, data = merged_df %>% filter(HER2_group == "HER2+"))

# Drawing
ggplot(merged_df, aes(x = pCR_binary, y = TMB, fill = pCR_binary)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
  stat_compare_means(method = "wilcox.test", label.y = max(merged_df$TMB, na.rm = TRUE) + 5) +
  labs(title = "TMB distribution by pCR status", y = "TMB (mutation count)", x = NULL) +
  theme_minimal() +
  scale_fill_manual(values = c("pCR" = "#E64B35", "RD" = "#4DBBD5"))
ggplot(merged_df, aes(x = pCR_binary, y = TMB, fill = pCR_binary)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
  stat_compare_means(method = "wilcox.test") +
  facet_wrap(~ HER2_group) +
  labs(title = "TMB by pCR and HER2 status", y = "TMB (mutation count)", x = NULL) +
  theme_minimal() +
  scale_fill_manual(values = c("pCR" = "#E64B35", "RD" = "#4DBBD5"))





# Import HRD data
hrd_df <- read_tsv("transneo-diagnosis-HRD.tsv")


# Create pCR binarization variables & HER2 groups (HER2- / HER2+)
clin_df <- clin_df %>%
  mutate(
    pCR_binary = ifelse(pCR.RD == "pCR", 1, 0),
    HER2_group = ifelse(HER2.status == "POS", "HER2+", "HER2-")
  )

# Combine HRD with clinical data
merged <- clin_df %>%
  inner_join(hrd_df, by = c("Donor.ID" = "Trial.ID"))

# Association analysis of HRD and pCR in the whole cohort: Logistic regression
model_hrd <- glm(pCR_binary ~ HRD.sum, data = merged, family = binomial)
summary(model_hrd)
exp(cbind(OR = coef(model_hrd), confint(model_hrd)))  # OR + CI


# HRD.sum non-parametric test in pCR vs non-PCR
wilcox.test(HRD.sum ~ pCR_binary, data = merged)


# HER2 subtype stratified analysis

# HER2- group
wilcox.test(HRD.sum ~ pCR_binary, data = merged %>% filter(HER2_group == "HER2-"))

# HER2+ group
wilcox.test(HRD.sum ~ pCR_binary, data = merged %>% filter(HER2_group == "HER2+"))

# -------------------------
# Visualization of HRD vs pCR (full cohort)
# -------------------------
ggplot(merged, aes(x = as.factor(pCR_binary), y = HRD.sum, fill = as.factor(pCR_binary))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, shape = 21, size = 2, alpha = 0.7, color = "black") +
  scale_fill_manual(values = c("red", "steelblue"), labels = c("RD", "pCR")) +
  scale_x_discrete(labels = c("RD", "pCR")) +
  labs(
    x = "Response (pCR vs RD)",
    y = "HRD Sum Score",
    title = "HRD score by pCR status"
  ) +
  theme_minimal(base_size = 14)



# read CIN data
cin_df <- read_tsv("transneo-diagnosis-ASCAT-CIN.tsv")
colnames(cin_df) <- c("Trial.ID", "CIN")


#  Data cleaning and merging
merged_df <- clin_df %>%
  rename(Trial.ID = Donor.ID) %>%
  inner_join(cin_df, by = "Trial.ID") %>%
  mutate(
    pCR_binary = ifelse(pCR.RD == "pCR", 1, 0),
    RCB_class = factor(RCB.category, levels = c("RCB-I", "RCB-II", "RCB-III", "pCR"))
  )

# Test whether CIN is correlated with the monotony of the RCB class
kruskal_test <- kruskal.test(CIN ~ RCB.category, data = merged_df)
print(kruskal_test)

# Compare the differences in CIN between pCR and RD
wilcox_test <- wilcox.test(CIN ~ pCR_binary, data = merged_df)
print(wilcox_test)

# Visualize CIN and RCB classes
p1 <- ggplot(merged_df, aes(x = RCB.category, y = CIN)) +
  geom_boxplot(aes(fill = RCB.category), outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  labs(title = "CIN vs RCB category", x = "RCB category", y = "CIN") +
  stat_compare_means(method = "kruskal.test", label.y = max(merged_df$CIN) * 1.05)

# Visualize the distribution of CIN in pCR vs RD
p2 <- ggplot(merged_df, aes(x = factor(pCR_binary, labels = c("RD", "pCR")), y = CIN)) +
  geom_boxplot(aes(fill = factor(pCR_binary)), outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_fill_manual(values = c("red", "skyblue")) +
  theme_minimal() +
  labs(title = "CIN vs pCR status", x = "pCR status", y = "CIN") +
  stat_compare_means(method = "wilcox.test", label.y = max(merged_df$CIN) * 1.05)

# Display graphics
ggarrange(p1, p2, ncol = 2, labels = c("A", "B"))


# cleaning data
clin_clean <- clin_df %>%
  mutate(
    iC10 = as.factor(iC10),
    pCR_binary = ifelse(pCR.RD == "pCR", 1, 0)  # 1 represents pCR，0 represents RD
  )

# Descriptive statistics: pCR rates were calculated by iC10 grouping
pcr_summary <- clin_clean %>%
  group_by(iC10) %>%
  summarise(
    n = n(),
    pCR_count = sum(pCR_binary, na.rm = TRUE),
    pCR_rate = round(mean(pCR_binary, na.rm = TRUE), 2)
  ) %>%
  arrange(desc(pCR_rate))

print(pcr_summary)

# Chi-square test: Whether iC10 and pCR.RD are independent
table_ic_pcr <- table(clin_clean$iC10, clin_clean$pCR.RD)
chisq_result <- chisq.test(table_ic_pcr)
print(chisq_result)

# Logistic regression model: pCR ~ iC10 (iC1 for the reference group)
model <- glm(pCR_binary ~ iC10, data = clin_clean, family = binomial)
summary(model)

# Show the OR (Odds Ratio) and confidence interval of the regression results
exp_coef <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE)
print(exp_coef)

# Visualization: Bar chart of pCR rates for different iC10 subtypes
ggplot(pcr_summary, aes(x = iC10, y = pCR_rate)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = paste0(pCR_count, "/", n)), vjust = -0.5) +
  ylim(0, 1) +
  labs(
    title = "pCR Rate by iC10 Subtype",
    x = "iC10 Subtype",
    y = "pCR Rate"
  ) +
  theme_minimal()

