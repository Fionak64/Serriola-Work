### 04_height_env_regression_stepwise_Lser_K9.R

### Author: Fiona Krammer

### Description:
### This script investigates plant height variation in Lactuca serriola across genetic groups (K = 9)
### and environmental gradients. It performs ANOVA to test for significant differences in height between
### genetic groups, and visualizes group-wise distributions and spatial trends (latitude, longitude).
### The script includes correlation analysis between height and environmental variables, followed by
### a stepwise regression to identify key environmental predictors of height using standardized coefficients.

# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------
library(openxlsx)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggpubr)
library(forcats)
library(effects)

# ------------------------------------------------------------------------------
# Load Data
# ------------------------------------------------------------------------------
phenotype <- read.xlsx("./Data/Phenotypic_data.xlsx", colNames = TRUE, startRow = 2)
geneticgroups <- read.xlsx("./Data/Lser200_pop_info_apr2025_long_form_changedregions.xlsx", sheet = "K9")
meta.1 <- read.xlsx("./Data/More_detailed_Serriola_meta_envser_July2021_with_Feb2025_inspection.xlsx", sheet = "Use")

# ------------------------------------------------------------------------------
# Prepare Phenotype Data
# ------------------------------------------------------------------------------
filtered_phenotype <- phenotype[phenotype$Accession %in% geneticgroups$LKID,]

# Subset by replicate and day
get_pheno <- function(replicate, day) {
  filtered_phenotype %>%
    filter(Replicate == replicate, Day == day) %>%
    select(4:ncol(.))
}

pheno1.78 <- get_pheno("1", "Day-78")
pheno2.78 <- get_pheno("2", "Day-78")
pheno1.93 <- get_pheno("1", "Day-93")
pheno2.93 <- get_pheno("2", "Day-93")

# Remove NA-only columns (identify from any of the matrices)
remove_na_cols <- function(pheno1, pheno2) {
  empty_cols <- names(pheno1)[colSums(is.na(pheno1)) == nrow(pheno1)]
  pheno1 <- pheno1[, !(names(pheno1) %in% empty_cols)]
  pheno2 <- pheno2[, !(names(pheno2) %in% empty_cols)]
  list(pheno1 = pheno1, pheno2 = pheno2)
}

ph78 <- remove_na_cols(pheno1.78, pheno2.78)
ph93 <- remove_na_cols(pheno1.93, pheno2.93)

# Calculate means
average.78 <- (ph78$pheno1 + ph78$pheno2) / 2
average.93 <- (ph93$pheno1 + ph93$pheno2) / 2

# Keep only ".mean" columns
average.78 <- average.78[, grepl("\\.mean$", names(average.78))]
average.93 <- average.93[, grepl("\\.mean$", names(average.93))]

# Remove all-NA columns
average.78 <- average.78[, colSums(is.na(average.78)) < nrow(average.78)]

# Add accession IDs
average.78$LKID <- filtered_phenotype$Accession[filtered_phenotype$Replicate == "1" & filtered_phenotype$Day == "Day-78"]
average.93$LKID <- filtered_phenotype$Accession[filtered_phenotype$Replicate == "1" & filtered_phenotype$Day == "Day-93"]

# ------------------------------------------------------------------------------
# Prepare Genetic Group Info
# ------------------------------------------------------------------------------
geneticgroups <- geneticgroups %>%
  mutate(LKID = trimws(LKID)) %>%
  filter(!is.na(region), ggroups != 9) %>%
  mutate(ggroups = ifelse(ggroups == 5, 2, ggroups))

geneticgroups <- geneticgroups[!is.na(geneticgroups$region), ]

average.78 <- left_join(average.78, geneticgroups[, c("LKID", "ggroups", "region")], by = "LKID")
average.93 <- left_join(average.93, geneticgroups[, c("LKID", "ggroups")], by = "LKID")

average.78 <- average.78[!is.na(average.78$ggroups), ]
average.93 <- average.93[!is.na(average.93$ggroups), ]

average.78$Group <- factor(average.78$ggroups)
average.93$Group <- factor(average.93$ggroups)

# ------------------------------------------------------------------------------
# ANOVA
# ------------------------------------------------------------------------------
getp <- function(expression, groups) {
  aout <- anova(aov(expression ~ groups))
  aout$`Pr(>F)`[1]
}

pheno_vars.78 <- average.78[, grepl("\\.mean$", names(average.78))]
pheno_vars.93 <- average.93[, grepl("\\.mean$", names(average.93))]

pvalues.78 <- apply(pheno_vars.78, 2, function(x) getp(x, average.78$Group))
pvalues.93 <- apply(pheno_vars.93, 2, function(x) getp(x, average.93$Group))

log10pvalues.78 <- -log10(pvalues.78)
log10pvalues.93 <- -log10(pvalues.93)

# ------------------------------------------------------------------------------
# Plot: Height at Day 78 by Genetic Group
# ------------------------------------------------------------------------------
create_height_plot <- function(df, day_label, pvalues, value_col = "height.mean") {
  max_height <- max(df[[value_col]], na.rm = TRUE)
  num_groups <- length(unique(df$Group))
  
  ggplot(df, aes(x = Group, y = .data[[value_col]])) +
    geom_boxplot(outlier.shape = NA, fill = "gray90", color = "black", width = 0.6, size = 0.3) +
    geom_jitter(aes(color = region), width = 0.15, size = 3, alpha = 0.7) +
    scale_color_manual(values = origin_colors) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste("Mean Height by Genetic Group (Day", day_label, ")"),
      x = "Genetic Group", y = "Mean Height (m)", color = "Region"
    ) +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.4)
    ) +
    annotate("text", x = num_groups - 1, y = max_height,
             label = paste0("p = ", signif(pvalues[value_col], 3)),
             hjust = 0, vjust = 1, size = 4, fontface = "italic")
}

plot_day78 <- create_height_plot(average.78, "78", pvalues.78)
ggsave("Boxplot_MeanHeight_Day78.png", plot_day78, width = 12, height = 9, dpi = 300)

# ------------------------------------------------------------------------------
# Latitude / Longitude Regressions
# ------------------------------------------------------------------------------
meta.1$LKID <- trimws(as.character(meta.1$LKID))

df_height <- data.frame(Accession = average.78$LKID, Height = average.78$height.mean)
df_merged <- merge(df_height, meta.1[, c("LKID", "LATITUDE", "LONGITUDE")],
                   by.x = "Accession", by.y = "LKID", all.x = TRUE) %>%
  filter(!is.na(Height), !is.na(LATITUDE), !is.na(LONGITUDE))

lm_fit.lat <- lm(Height ~ LATITUDE, data = df_merged)
lm_fit.long <- lm(Height ~ LONGITUDE, data = df_merged)
print(summary(lm_fit.lat))
print(summary(lm_fit.long))

# Merge coordinates into phenotype
average.78 <- left_join(average.78, meta.1[, c("LKID", "LATITUDE", "LONGITUDE")], by = "LKID") %>%
  filter(!is.na(height.mean), !is.na(LATITUDE), !is.na(LONGITUDE))

# ------------------------------------------------------------------------------
# Plot Height vs Latitude/Longitude
# ------------------------------------------------------------------------------
plot_height_vs_coord <- function(df, coord_col, subplot_label) {
  ggplot(df, aes(x = .data[[coord_col]], y = height.mean)) +
    geom_point(size = 2, alpha = 0.8, color = "#1f78b4") +
    geom_smooth(method = "lm", se = TRUE, color = "#e31a1c", size = 1.2) +
    stat_regline_equation(aes(label = ..eq.label..), label.x.npc = 0.02, label.y.npc = 0.95) +
    stat_regline_equation(aes(label = ..rr.label..), label.x.npc = 0.02, label.y.npc = 0.87) +
    labs(x = paste0(coord_col, " (°)"), y = "Mean Plant Height (m)") +
    annotate("text", x = -Inf, y = Inf, label = subplot_label,
             hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold") +
    theme_minimal(base_size = 14)
}

p_lat_78 <- plot_height_vs_coord(average.78, "LATITUDE", "A")
p_long_78 <- plot_height_vs_coord(average.78, "LONGITUDE", "B")
combined_plot <- plot_grid(p_lat_78, p_long_78, nrow = 1)
ggsave("height_vs_lat_long_day78.png", combined_plot, width = 12, height = 5, dpi = 300)

# ------------------------------------------------------------------------------
# Correlation Heatmap
# ------------------------------------------------------------------------------
env_data <- meta.1[, 155:173]
env_data$altitude <- meta.1$Altitude
rownames(env_data) <- meta.1$LKID
env_data <- env_data[rowSums(is.na(env_data)) == 0,]

rownames(average.78) <- average.78$LKID
sig.pheno.78 <- average.78[, log10pvalues.78 > 1.3]
sig.pheno.78 <- sig.pheno.78[, 1:(ncol(sig.pheno.78) - 4)]

common_ids <- intersect(rownames(sig.pheno.78), rownames(env_data))
pheno.sub <- sig.pheno.78[common_ids, , drop = FALSE]
env.sub <- env_data[common_ids, ]

get_full_correlation_df <- function(env_data, pheno_data) {
  expand.grid(EnvVar = colnames(env_data), PhenoVar = colnames(pheno_data)) %>%
    rowwise() %>%
    do({
      x <- env_data[[.$EnvVar]]
      y <- pheno_data[[.$PhenoVar]]
      if (sum(complete.cases(x, y)) > 2) {
        test <- cor.test(x, y, method = "pearson")
        data.frame(EnvVar = .$EnvVar, PhenoVar = .$PhenoVar,
                   Correlation = test$estimate, Pvalue = test$p.value)
      } else {
        data.frame(EnvVar = .$EnvVar, PhenoVar = .$PhenoVar,
                   Correlation = NA, Pvalue = NA)
      }
    }) %>%
    bind_rows() %>%
    mutate(Significance = case_when(
      Pvalue < 0.001 ~ "***",
      Pvalue < 0.01 ~ "**",
      Pvalue < 0.05 ~ "*",
      TRUE ~ ""
    ))
}

cor_df <- get_full_correlation_df(env.sub, pheno.sub)

gg_cor <- ggplot(cor_df, aes(x = PhenoVar, y = fct_rev(EnvVar), fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste0(round(Correlation, 2), Significance)), size = 3) +
  scale_fill_gradient2(low = "#d7301f", mid = "white", high = "#1a9850",
                       midpoint = 0, limits = c(-1, 1), name = "Pearson\nr") +
  theme_minimal(base_size = 13) +
  labs(title = "Pearson Correlation: Phenotypes (Day 78) vs Environmental Variables") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

ggsave("pheno78_env_correlation.png", gg_cor, width = 10, height = 6.5, dpi = 300)

# Combine figures for Figure S5
final_plot <- plot_grid(combined_plot, gg_cor, labels = c("", "C"), ncol = 1,
                        rel_heights = c(1, 1.5), label_size = 14)
ggsave("Figure_S5.png", final_plot, width = 9, height = 10, dpi = 300)


# ------------------------------------------------------------------------------
# Stepwise Model Height
# ------------------------------------------------------------------------------

# Create dataset
height.78 <- pheno.sub[, "height.mean", drop = FALSE]
data.78 <- cbind(height.78, env.sub)
data.78 <- data.78[, !(names(data.78) %in% c("Group", "Region"))]

data.std <- data.78
num_vars <- sapply(data.std, is.numeric)
data.std[num_vars] <- scale(data.std[num_vars])

model.std <- lm(height.mean ~ ., data = data.std)
step.model.std <- step(model.std, direction = "both", trace = FALSE)
# Plot partial effects
plot(allEffects(step.model.std), main = "Effects on Height")

coefs <- summary(step.model.std)$coefficients
coef_df <- data.frame(
  Predictor = rownames(coefs)[-1],  # exclude intercept
  Coefficient = coefs[-1, "Estimate"],
  P_value = coefs[-1, "Pr(>|t|)"],
  stringsAsFactors = FALSE
)


# Use bio_labels for better labeling (ensure this vector exists)
if (!exists("bio_labels")) {
  bio_labels <- setNames(coef_df$Predictor, coef_df$Predictor)  # fallback: use raw names
}

coef_df$PredictorLabel <- bio_labels[coef_df$Predictor]
coef_df$PredictorLabel[is.na(coef_df$PredictorLabel)] <- coef_df$Predictor[is.na(coef_df$PredictorLabel)]

# Determine direction
coef_df$Direction <- ifelse(coef_df$Coefficient > 0, "Positive", "Negative")

# Plot
stepwise_height <- ggplot(coef_df, aes(x = fct_reorder(PredictorLabel, Coefficient), y = Coefficient, fill = Direction)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  coord_flip() +
  geom_text(aes(label = paste0("p = ", signif(P_value, 2))),
            hjust = 0.15,
            nudge_y = max(abs(coef_df$Coefficient), na.rm = TRUE) * 0.05,
            size = 4, fontface = "italic", color = "black") +
  scale_fill_manual(values = c("Positive" = "#2c7fb8", "Negative" = "indianred1")) +
  labs(
    title = "Standardized Coefficients for Height",
    x = "Predictor",
    y = "Standardized Coefficient",
    fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 13),
    plot.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.4),
    legend.position = "none"
  ) +
  expand_limits(y = max(abs(coef_df$Coefficient), na.rm = TRUE) * 1.2)


