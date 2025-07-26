### 02_env_analysis_ANOVA_visualization_Lser_K9.R

### Author: Fiona Krammer
### Description:
### This script analyzes environmental variation across genetic groups (K = 9) in Lactuca serriola.
### It performs ANOVA to identify significant differences in environmental variables between groups,
### generates correlation heatmaps of environmental traits, and creates visualizations (boxplots,
### density plots) for key variables including temperature, precipitation, and altitude.


library(openxlsx)
library(dplyr)
library(gplots)
library(GGally)
library(cowplot)

# ------------------------------------------------------------------------------
# PREPARE ENVIRONMENTAL DATA FILES
# ------------------------------------------------------------------------------

meta.1 <- read.xlsx("./Data/More_detailed_Serriola_meta_envser_July2021_with_Feb2025_inspection.xlsx", sheet = "Use")
geneticgroups <- read.xlsx("./Data/Lser200_pop_info_apr2025_long_form_changedregions.xlsx", sheet = "K9")
meta <- meta.1[meta.1$LKID %in% geneticgroups$LKID]

# ------------------------------------------------------------------------------
# PREPARE ENVIRONMENTAL DATA FILES
# ------------------------------------------------------------------------------

# Only environmental data
rownames(meta.1) <- meta.1$LKID
meta <- subset(meta.1, select = -c(1:154))
meta <- subset(meta, select = -c(20:ncol(meta)))
meta$altitude <- meta.1$Altitude

# Filter phenotype file so only the LKIDs from the genetic groups file are left
filtered.env <- meta[rownames(meta) %in% geneticgroups$LKID,]


# ------------------------------------------------------------------------------
# PREPARE GENETIC GROUPS
# ------------------------------------------------------------------------------

# Trim LKID and map to groups
geneticgroups <- geneticgroups[!is.na(geneticgroups$region), ]
geneticgroups$LKID <- trimws(geneticgroups$LKID)
LKID_to_group <- setNames(geneticgroups$ggroups, geneticgroups$LKID)

# Assign groups to meta rownames
groups <- LKID_to_group[rownames(filtered.env)]

# Move any group 5 samples to group 2 -> most similar, we don't want a group with only 1 sample
groups[groups == 5] <- 2

# Convert to factor (drop unused levels automatically after filtering)
groups <- factor(groups)

# Check for unmatched samples
unmatched <- rownames(filtered.env)[!rownames(filtered.env) %in% names(groups)]
length(groups) == nrow(filtered.env)  # Should be TRUE
print(table(groups))

# Remove group 9 samples from all datasets
keep_samples <- groups != 9
meta <- filtered.env[keep_samples, , drop = FALSE]
groups <- droplevels(groups[keep_samples])


# ------------------------------------------------------------------------------
# REMOVING LKIDS DUE TO NA ROWS
# ------------------------------------------------------------------------------

# There are some lines where there is no environmental data available. I'll remove these rows
complete_rows <- complete.cases(meta)
meta_clean <- meta[complete_rows, ]
groups_clean <- groups[complete_rows]
length(groups_clean) == nrow(meta_clean)


# ------------------------------------------------------------------------------
# ANOVA
# ------------------------------------------------------------------------------

# ANOVA function
getp <- function(expression, groups) {
  aout <- anova(aov(expression ~ groups))
  return(aout$`Pr(>F)`[1])  # Extract first p-value
}

pvalues.env <- apply(meta_clean, 2, function(x) getp(x, groups_clean))

# Transform pvalues: -log10
log10pvalues.env <- -log10(pvalues.env)
log10pvalues.env <- as.matrix(log10pvalues.env)
rownames(log10pvalues.env) <- colnames(meta)
log10pvalues.env_sorted <- apply(log10pvalues.env, 2, function(x) sort(x, decreasing = TRUE))

length(log10pvalues.env[log10pvalues.env>1.3,]) # 20
length(log10pvalues.env[log10pvalues.env>10,]) # 11


# ------------------------------------------------------------------------------
# ENVIRONMENTAL DATA BOXPLOTS
# ------------------------------------------------------------------------------

region_map <- setNames(geneticgroups$region, geneticgroups$LKID)
meta_clean$Region <- region_map[rownames(meta_clean)]
meta_clean$Group <- groups_clean

meta_envplot <- meta_clean %>%
  filter(Group %in% 1:8, !is.na(Region))

origin_colors <- c(
  "Turkey"        = "#E69F00",
  "Europe"        = "#56B4E9",
  "Caspian Sea"   = "#009E73",
  "Middle East"   = "#F0E442",
  "Central Asia"  = "#CC79A7"
)

# FIGURE 3A: Bio01 (Annual Mean Temperature) - converted from tenths of °C to °C
plot_bio01 <- ggplot(meta_envplot, aes(x = as.factor(Group), y = .data[["bio01"]] / 10)) +
  geom_boxplot(outlier.shape = NA, fill = "gray90", color = "black", width = 0.6, lwd = 0.3) +
  geom_jitter(aes(color = Region), width = 0.15, size = 1.4, alpha = 0.7) +
  scale_color_manual(values = origin_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Annual Mean Temperature",
    x = "Genetic Group",
    y = "Annual Mean Temperature (°C)"
  ) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
  ) +
  annotate(
    "text",
    x = 6,  # Adjust based on your number of groups
    y = max(meta_envplot$bio01, na.rm = TRUE) / 10 + 0.2,
    label = paste0("p = ", signif(pvalues.env["bio01"], 3)),
    hjust = 0,
    vjust = 0,
    size = 4,
    fontface = "italic"
  )

# FIGURE 3B: Bio12
plot_bio12 <- ggplot(meta_envplot, aes(x = as.factor(Group), y = .data[["bio12"]])) +
  geom_boxplot(outlier.shape = NA, fill = "gray90", color = "black", width = 0.6, lwd = 0.3) +
  geom_jitter(aes(color = Region), width = 0.15, size = 1.4, alpha = 0.7) +
  scale_color_manual(values = origin_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = paste0("Annual Precipitation"),
    x = "Genetic Group",
    y = "Annual Precipitation (mm)",
    color = "Region"
  ) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
  ) +
  annotate(
    "text",
    x = 6,  # Adjust this value as needed depending on how many groups you have (e.g., if you have 8 groups)
    y = max(meta_envplot$bio12, na.rm = TRUE),
    label = paste0("p = ", signif(pvalues.env["bio12"], 3)),
    hjust = 0,
    vjust = 0,
    size = 4,
    fontface = "italic"
  )

# FIGURE 3C: Distribution of Altitude Across Regions
plot_box <- ggplot(meta_envplot, aes(x = Region, y = altitude, fill = Region)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6, lwd = 0.3) +
  geom_jitter(aes(color = Region), width = 0.15, size = 1.4, alpha = 0.7) +
  scale_fill_manual(values = origin_colors) +
  scale_color_manual(values = origin_colors) +
  theme_minimal(base_size = 14) +
  labs(title = "Altitude Distribution by Region", y = "Altitude (m)", x = NULL) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 9.5, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 13),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
  )


# FIGURE 3D: Altitude
plot_altitude2 <- ggplot(meta_envplot, aes(x = as.factor(Group), y = .data[["altitude"]])) +
  geom_boxplot(outlier.shape = NA, fill = "gray90", color = "black", width = 0.6, lwd = 0.3) +
  geom_jitter(aes(color = Region), width = 0.15, size = 1.4, alpha = 0.7) +
  scale_color_manual(values = origin_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = paste0("Altitude"),
    x = "Genetic Group",
    y = "Altitude (m)",
    color = "Region"
  ) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
  ) +
  annotate(
    "text",
    x = 6,  # Adjust this value as needed depending on how many groups you have (e.g., if you have 8 groups)
    y = max(meta_envplot$altitude, na.rm = TRUE),
    label = paste0("p = ", signif(pvalues.env["altitude"], 3)),
    hjust = 0,
    vjust = 0,
    size = 4,
    fontface = "italic"
  )

# FIGURE 3: merge plots with cowplot
combined_plot <- plot_grid(
  plot_bio01,
  plot_bio12,
  plot_box,
  plot_altitude2,
  labels = c("A", "B", "C", "D"),
  label_size = 14,
  ncol = 2,
  align = "hv",
  axis = "tblr"
)

# Save FIGURE 3
ggsave("combined_plots_figure3.png", combined_plot, width = 10, height = 8, dpi = 300)

# FIGURE S2B: Density plot altitude
plot.density <- ggplot(meta_envplot, aes(x = altitude, fill = Region)) +
  geom_density(alpha = 0.8) +
  scale_fill_manual(values = origin_colors) +
  theme_minimal() +
  labs(title = "Altitude Density by Region", x = "Altitude (m)", y = "Density") +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 9.5, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 13, face = "bold"),
    plot.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
  )


# ------------------------------------------------------------------------------
# CORRELATION HEAT MAPS
# ------------------------------------------------------------------------------

# Remove "Region" and "Group" columns
meta_cor <- meta_clean[, !(names(meta_clean) %in% c("Region", "Group"))]

# Z-score normalization
meta_scaled <- scale(meta_cor)

# Convert to data frame to ensure compatibility
meta_scaled_df <- as.data.frame(meta_scaled)

# FIGURE S2A: Correlation
p <- ggcorr(meta_scaled_df,
            method = c("everything", "pearson"),
            label = TRUE,
            label_round = 2,
            label_size = 3.2,      # Larger label font for readability
            low = "blue", mid = "white", high = "red",
            midpoint = 0,
            size = 3.8) +          # Slightly thicker tiles for better visibility
  ggtitle("Correlation Between Environmental Parameters") +
  theme_minimal(base_size = 20) +  # Bigger overall font size for axis, legend, etc.
  theme(
    plot.title = element_blank(),   # Bigger, bold, centered title
    axis.text = element_text(size = 14, color = "black"),              # Bigger axis text
    axis.title = element_text(size = 16, face = "bold"),               # Bold axis titles
    legend.title = element_text(size = 16, face = "bold"),             # Bold legend title
    legend.text = element_text(size = 14)                              # Larger legend text
  )

# FIGURE S2: Merge figures with cowplot
figureS2 <- plot_grid(p, plot.density,
          nrow = 2, labels = c("A", "B"), rel_heights = c(1.5, 1), label_size = 16)
# Save high-resolution figure
ggsave("Figure_S2.png", figureS2, width = 8, height = 12, dpi = 300)


