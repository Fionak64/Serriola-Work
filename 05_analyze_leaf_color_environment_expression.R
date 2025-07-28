### 05_analyze_leaf_color_environment_expression.R
### Author: Fiona Krammer
### Description: This script analyzes variation in leaf greenness of Lactuca serriola using RGB values from drone images.
###              It assesses differences across genetic groups, performs correlation analyses with environmental variables
###              and gene expression, and visualizes representative drone images and key associations.


library(dplyr)
library(ggplot2)
library(cowplot)
library(openxlsx)
library(Hmisc)
library(broom)
library(sf)
library(rnaturalearth)
library(scales)

#-------------------------------------------------------------------------------
# LOAD FILES
#-------------------------------------------------------------------------------

load("./Data/leaf_shape_vector.out")
load(file = "./Data/obj_all.pl_ser_rep1_2506_rgb_dsm_msp_red_nd.out.crdownload")
geneticgroups <- read.xlsx("./Data/Lser200_pop_info_apr2025_long_form_changedregions.xlsx", sheet = "K9")
meta.1 <- read.xlsx("./Data/More_detailed_Serriola_meta_envser_July2021_with_Feb2025_inspection.xlsx", sheet = "Use")
cnt.cpm <- read.csv("./Data/CPM_SER_Tissues_Merged.csv", row.names = 1)
load("./Data/obj_gene.meta.out")
GO.file <- read.delim("./Data/20231004_Lactuca_sativa.annotation_overview_with_V8_positions.tsv")

#-------------------------------------------------------------------------------
# GENERATE DATA FRAME WITH MEDIAN COLOR VALUE PER LINE
#-------------------------------------------------------------------------------

line_rgb <- all.pl %>%
  group_by(use.lk) %>%
  summarise(
    mean_R = mean(red, na.rm = TRUE),
    mean_G = mean(green, na.rm = TRUE),
    mean_B = mean(blue, na.rm = TRUE)
  ) %>%
  rename(LKID = use.lk)

#-------------------------------------------------------------------------------
# GENETIC GROUPS VECTOR
#-------------------------------------------------------------------------------

geneticgroups <- geneticgroups[!is.na(geneticgroups$region), ]
geneticgroups <- geneticgroups[!geneticgroups$ggroups == 9, ]
geneticgroups$ggroups[geneticgroups$ggroups == 5] <- 2
geneticgroups$LKID <- trimws(geneticgroups$LKID)
LKID_to_group <- setNames(geneticgroups$ggroups, geneticgroups$LKID)

# Assign groups to LKIDs in line_rgb
groups <- LKID_to_group[line_rgb$LKID]

# Convert to factor (drop unused levels automatically after filtering)
groups <- factor(groups)

# Identify LKIDs that do not have a group assigned
unmatched <- line_rgb$LKID[!line_rgb$LKID %in% names(groups)]
line_rgb <- line_rgb[!line_rgb$LKID %in% unmatched, ]
groups <- groups[names(groups) %in% line_rgb$LKID]

length(groups) == nrow(line_rgb)  # Should be TRUE
print(table(groups))

# Add the groups to line_rgb
line_rgb$Group <- groups[line_rgb$LKID]
line_rgb <- line_rgb %>% left_join(geneticgroups, by = "LKID")

# Add Leaf info
line_rgb$Leaf <- leafshape_vector[match(line_rgb$LKID, names(leafshape_vector))]

#-------------------------------------------------------------------------------
# PLOT MEAN VALUES PER LINE
#-------------------------------------------------------------------------------

# Define region map
region_map <- setNames(geneticgroups$region, geneticgroups$LKID)

# Define color palette with NA group
origin_colors <- c(
  "Turkey"       = "#E69F00",
  "Europe"       = "#56B4E9",
  "Caspian Sea"  = "#009E73",
  "Middle East"  = "#F0E442",
  "Central Asia" = "#CC79A7"
)

# Plot function
plot_lab_greenleaf <- function(df, value_col, label_name) {
  
  df <- df %>%
    mutate(
      Group = as.factor(Group),
      region = region_map[LKID],
      region = ifelse(is.na(region), "NA", region)
    ) %>%
    filter(Group %in% as.character(1:8), !is.na(.data[[value_col]]))
  
  # ANOVA and p-value
  formula <- as.formula(paste0(value_col, " ~ Group"))
  anova_result <- aov(formula, data = df)
  p_value <- summary(anova_result)[[1]][["Pr(>F)"]][1]
  
  # Determine y-position for p-value annotation
  y_max <- max(df[[value_col]], na.rm = TRUE)
  
  plot_main <- ggplot(df, aes(x = Group, y = .data[[value_col]])) +
    geom_boxplot(outlier.shape = NA, fill = "gray90", color = "black", width = 0.6, lwd = 0.3) +
    geom_jitter(aes(color = region), width = 0.15, size = 1.4, alpha = 0.7) +
    scale_color_manual(values = origin_colors) +
    theme_minimal(base_size = 14) +
    labs(
      title = "",
      x = "Genetic Group",
      y = label_name,
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
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
    ) +
    annotate(
      "text",
      x = 6,
      y = y_max,
      label = paste0("p = ", signif(p_value, 3)),
      hjust = 0,
      vjust = 0,
      size = 4,
      fontface = "italic"
    )
  
  return(plot_main)
}

geneticgroups_leaf_greenness <- plot_lab_greenleaf(df = line_rgb, value_col = "mean_G", label_name = "Leaf Greenness")
ggsave("color_G_by_group.png", geneticgroups_leaf_greenness, width = 9, height = 6, dpi = 300)

#-------------------------------------------------------------------------------
# LEAF GREENNESS ~ ENVIRONMENT CORRELATION ANALYSIS
#-------------------------------------------------------------------------------

# Prepare environmental data
rownames(meta.1) <- meta.1$LKID
meta.0 <- meta.1[meta.1$LKID %in% names(groups), ]
meta <- meta.0[, -(1:154)]
meta <- meta[, 1:19]  # Keep only environmental BIO variables
meta$altitude <- meta.0$Altitude


# Merge leaf greenness with environmental metadata
df_corr <- line_rgb %>%
  select(LKID, mean_G) %>%
  filter(LKID %in% rownames(meta)) %>%
  left_join(meta %>% mutate(LKID = rownames(meta)), by = "LKID")

# Calculate Spearman correlation
rcorr_result <- rcorr(as.matrix(df_corr %>% select(-LKID)), type = "spearman")
corr_table <- data.frame(
  Variable = rownames(rcorr_result$r),
  rho = rcorr_result$r["mean_G", ],
  p = rcorr_result$P["mean_G", ],
  stringsAsFactors = FALSE
)

# Remove self-correlation
corr_table <- corr_table %>%
  filter(Variable != "mean_G")

# Set BIO labels
bio_labels <- c(
  bio01 = "Annual Mean Temperature",
  bio02 = "Mean Diurnal Range",
  bio03 = "Isothermality",
  bio04 = "Temperature Seasonality",
  bio05 = "Max Temp. of Warmest Month",
  bio06 = "Min Temp. of Coldest Month",
  bio07 = "Temperature Annual Range",
  bio08 = "Mean Temp. of Wettest Quarter",
  bio09 = "Mean Temp. of Driest Quarter",
  bio10 = "Mean Temp. of Warmest Quarter",
  bio11 = "Mean Temp. of Coldest Quarter",
  bio12 = "Annual Precipitation",
  bio13 = "Precip. of Wettest Month",
  bio14 = "Precip. of Driest Month",
  bio15 = "Precip. Seasonality",
  bio16 = "Precip. of Wettest Quarter",
  bio17 = "Precip. of Driest Quarter",
  bio18 = "Precip. of Warmest Quarter",
  bio19 = "Precip. of Coldest Quarter",
  Altitude = "Altitude"
)

# Annotate and format
corr_table <- corr_table %>%
  mutate(
    label = ifelse(Variable %in% names(bio_labels), bio_labels[Variable], Variable),
    signif = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE ~ ""
    ),
    label_plot = paste0(label, " ", signif)
  ) %>%
  arrange(desc(rho))

# Plot: Figure S7D
leaf_env_correlation <- ggplot(corr_table, aes(x = reorder(label_plot, rho), y = rho, fill = rho > 0)) +
  geom_col(show.legend = FALSE) +
  geom_text(
    aes(label = paste0("p = ", formatC(p, format = "e", digits = 2))),
    hjust = ifelse(corr_table$rho > 0, -0.1, 1.1),
    size = 3.5
  ) +
  scale_fill_manual(values = c("#D55E00", "#56B4E9")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Environmental Correlates of Leaf Greenness",
    subtitle = "Spearman's rank correlation (ρ) with significance",
    x = NULL, y = expression(Spearman~rho)
  ) +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    axis.text = element_text(size = 11),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  coord_flip() +
  expand_limits(y = c(-1, 1))

ggsave("./Results/LeafGreenness_Environment_Correlation_Plot.png", leaf_env_correlation, width = 8, height = 6, dpi = 300)
print(leaf_env_correlation)

#-------------------------------------------------------------------------------
# STEPWISE APPROACH
#-------------------------------------------------------------------------------

# Drop LKID and remove rows with missing data
df_model <- df_corr %>% dplyr::select(-LKID) %>% na.omit()

# Full model with all variables
full_model <- lm(mean_G ~ ., data = df_model)

# Null model (intercept only)
null_model <- lm(mean_G ~ 1, data = df_model)

# Stepwise regression (both directions)
step_model <- step(null_model,
                   scope = list(lower = null_model, upper = full_model),
                   direction = "both",
                   trace = TRUE)

# Output model summary
summary(step_model)

# Extract coefficients and remove intercept
stepwise_coefs <- tidy(step_model) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    direction = ifelse(estimate > 0, "Positive", "Negative"),
    term_clean = gsub("`", "", term),  # Clean any backticks
    label = bio_labels[term_clean],
    label = ifelse(is.na(label), term_clean, label)  # fallback to term if not found
  )

# Sort terms by absolute effect size
stepwise_coefs <- stepwise_coefs %>%
  mutate(label = reorder(label, estimate))

# Plot
stepwise_barplot <- ggplot(stepwise_coefs, aes(x = label, y = estimate, fill = direction)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.2f", estimate)),
            hjust = ifelse(stepwise_coefs$estimate > 0, -0.15, 1.1),
            color = "black", size = 4) +
  scale_fill_manual(values = c("Positive" = "#56B4E9", "Negative" = "#D55E00")) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    title = "",
    x = NULL,
    y = "Effect size (β coefficient)"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    panel.grid.major.y = element_blank()
  ) +
  expand_limits(y = c(min(stepwise_coefs$estimate) * 1.2, max(stepwise_coefs$estimate) * 1.2))

#-------------------------------------------------------------------------------
# CORRELATION LEAF GREENNESS WITH LONGITUDE & LATITUDE
#-------------------------------------------------------------------------------

ser.lonlat <- left_join(meta.0, line_rgb[, c("LKID", "mean_G")], by = "LKID")

# Linear model: mean_G vs Longitude
lm_long <- lm(mean_G ~ LONGITUDE, data = ser.lonlat)

# Linear model: mean_G vs Latitude
lm_lat <- lm(mean_G ~ LATITUDE, data = ser.lonlat)

# Extract stats for Longitude
summary_long <- summary(lm_long)
r2_long <- round(summary_long$r.squared, 3)
pval_long <- summary_long$coefficients[2, 4]
pval_long_text <- ifelse(pval_long < 0.001, "p < 0.001", paste0("p = ", signif(pval_long, 3)))

# Extract stats for Latitude
summary_lat <- summary(lm_lat)
r2_lat <- round(summary_lat$r.squared, 3)
pval_lat <- summary_lat$coefficients[2, 4]
pval_lat_text <- ifelse(pval_lat < 0.001, "p < 0.001", paste0("p = ", signif(pval_lat, 3)))

# Plot mean_G vs LONGITUDE with stats annotation
p_long <- ggplot(ser.lonlat, aes(x = LONGITUDE, y = mean_G)) +
  geom_point(shape = 21, fill = "#41b6c4", color = "black", size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "darkblue", size = 1) +
  labs(
    x = "Longitude",
    y = "Leaf Greenness"
    #title = "Leaf Greenness vs Longitude"
  ) +
  annotate("text", x = min(ser.lonlat$LONGITUDE, na.rm=TRUE), y = max(ser.lonlat$mean_G, na.rm=TRUE), 
           label = paste0("R² = ", r2_long, "\n", pval_long_text), 
           hjust = 0, vjust = 1, size = 5, color = "black") +
  theme_minimal(base_size = 14)

# Plot mean_G vs LATITUDE with stats annotation
p_lat <- ggplot(ser.lonlat, aes(x = LATITUDE, y = mean_G)) +
  geom_point(shape = 21, fill = "#41b6c4", color = "black", size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "darkblue", size = 1) +
  labs(
    x = "Latitude",
    y = "Leaf Greenness"
    #title = "Leaf Greenness vs Latitude"
  ) +
  annotate("text", x = min(ser.lonlat$LATITUDE, na.rm=TRUE), y = max(ser.lonlat$mean_G, na.rm=TRUE), 
           label = paste0("R² = ", r2_lat, "\n", pval_lat_text), 
           hjust = 0, vjust = 1, size = 5, color = "black") +
  theme_minimal(base_size = 14)

# Combine B and C side by side
row2 <- plot_grid(
  p_long, p_lat,
  labels = c("B", "C"),
  label_size = 14,
  nrow = 1,
  rel_widths = c(1, 1)
)

# ------------------------------------------------------------------------------
# MAP
# ------------------------------------------------------------------------------

# Ensure proper coordinate handling
ser.lonlat <- ser.lonlat[!ser.lonlat$Country == "KEN", ]
sf::sf_use_s2(FALSE)
st_crs(ser.lonlat)
ser.lonlat$mean

# Load world map and crop to Europe
world_map <- ne_countries(scale = 10, returnclass = 'sf')
europe_map <- st_crop(world_map, xmin = -10, xmax = 80, ymin = 20, ymax = 60)

p <- ggplot() +
  # Base map
  geom_sf(data = europe_map, fill = 'gray98', color = "gray70", size = 0.3) +
  
  # Greenness points with black outline
  geom_point(
    data = ser.lonlat,
    aes(x = LONGITUDE, y = LATITUDE, fill = mean_G),
    shape = 21,
    size = 4,
    color = "black",
    stroke = 0.1  # Adjust thickness of the black outline
  ) +
  
  # Improved color scale (now for fill)
  scale_fill_gradientn(
    colors = c("#fefcc7", "#a1dab4", "#41b6c4", "#2c7fb8", "#253494"),
    name = "Leaf Greenness",
    guide = guide_colorbar(barwidth = 1.2, barheight = 10)
  ) +
  
  # Titles
  labs(
    title = "Geographical Variation in Leaf Greenness of L. serriola",
    subtitle = "Mean green pixel intensity from RGB drone imagery",
    x = "Longitude", y = "Latitude"
  ) +
  
  coord_sf(expand = FALSE) +
  
  # Theme
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    axis.text = element_text(color = "black", size = 12),
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_line(color = "gray90", size = 0.2),
    plot.margin = margin(12, 12, 12, 12)
  )

# Final combined figure
combined_plot_figureS7 <- plot_grid(
  p,      # A: map
  row2,   # B and C: regressions
  leaf_env_correlation,  # D: leaf-env correlation
  labels = c("A", "", "D"),
  label_size = 14,
  ncol = 1,
  rel_heights = c(2.5, 2, 3)
)
combined_plot_figureS7 # Figure S7

#-------------------------------------------------------------------------------
# CORRELATION LEAF GREENNESS WITH GENE EXPRESSION
#-------------------------------------------------------------------------------

# Prepare cnt.cpm file
leaf.selection <- cnt.cpm[, grepl("leaf_", colnames(cnt.cpm))]
colnames(leaf.selection) <- sub("^leaf_", "", colnames(leaf.selection))

# Leaf counts: get log2 counts per million for leaf samples
occ.leaf1 <- apply(leaf.selection, 1, function(x) sum(x > 0)) # Count genes expressed in more than 0 samples
upd.leaf.selection <- leaf.selection[occ.leaf1 > 20, ] # Keep genes expressed in more than 20 samples
log2.cpm.leaf <- log2((upd.leaf.selection + 1) / apply(upd.leaf.selection, 1, mean, na.rm = TRUE)) # Row-wise transformation

# match files
shared_LKIDs <- intersect(line_rgb$LKID, colnames(log2.cpm.leaf))
green_leaf <- line_rgb[line_rgb$LKID %in% shared_LKIDs, ]
log2.cpm.leaf <- log2.cpm.leaf[, shared_LKIDs]

# Spearman correlation analysis for each gene
cor_results <- apply(log2.cpm.leaf, 1, function(expr) {
  cor(expr, green_leaf$mean_G, method = "spearman")
})

# Compute p-values as well
p_values <- apply(log2.cpm.leaf, 1, function(expr) {
  cor.test(expr, green_leaf$mean_G, method = "spearman")$p.value
})

# Combine results into a dataframe
correlation_df <- data.frame(
  Gene = rownames(log2.cpm.leaf),
  Spearman_rho = cor_results,
  P_value = p_values
)

# Adjust for multiple testing
correlation_df$FDR <- p.adjust(correlation_df$P_value, method = "fdr")

# Sort by correlation strength
correlation_df <- correlation_df[order(correlation_df$Spearman_rho, decreasing = TRUE), ]

# View top correlated genes
head(correlation_df)

top.genes <- correlation_df$Gene[1:200]

homologs.top.genes <- GO.file$`best.hit.ARAPORT11..mmseqs2.`[GO.file$gene.ID %in% top.genes]
homologs.top.genes <- homologs.top.genes[homologs.top.genes != ""]
homologs.top.genes <- gsub("\\.\\d+$", "", trimws(homologs.top.genes))
homologs.top.genes <- homologs.top.genes[!duplicated(homologs.top.genes)]
sum(duplicated(homologs.top.genes)) # 0

top.100 <- homologs.top.genes[1:100]
writeLines(top.100, "homologs_top100_greenleaves.txt")

#-------------------------------------------------------------------------------
# DRONE IMAGES FIGURE 6
#-------------------------------------------------------------------------------

plot_cropped_leaf <- function(data, lkid, crop_fraction = 0.6) {
  # 1. Filter for the selected LKID
  df <- data[data$use.lk == lkid, ]
  
  if (nrow(df) == 0) {
    stop(paste("No data found for", lkid))
  }
  
  # 2. Compute square crop box
  xmin <- min(df$x)
  xmax <- max(df$x)
  ymin <- min(df$y)
  ymax <- max(df$y)
  
  side <- crop_fraction * min(xmax - xmin, ymax - ymin)
  xmid <- (xmin + xmax) / 2
  ymid <- (ymin + ymax) / 2
  
  xlim <- c(xmid - side / 2, xmid + side / 2)
  ylim <- c(ymid - side / 2, ymid + side / 2)
  
  # 3. Crop
  square_crop <- subset(df,
                        x >= xlim[1] & x <= xlim[2] &
                          y >= ylim[1] & y <= ylim[2])
  
  if (nrow(square_crop) == 0) {
    stop(paste("Crop region for", lkid, "is empty. Try increasing crop_fraction."))
  }
  
  # 4. Convert to RGB
  use_color <- rgb(square_crop$red, square_crop$green, square_crop$blue, maxColorValue = 255)
  
  # 5. Plot
  plot_out <- ggplot(square_crop) +
    geom_tile(aes(x, y), fill = use_color) +
    coord_fixed() +
    theme_void() +
    theme(
      plot.margin = margin(0, 0, 0, 0),
      plot.title = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  
  return(plot_out)
}

plot.lk329 <- plot_cropped_leaf(all.pl, "LK329")
plot.lk235 <- plot_cropped_leaf(all.pl, "LK235")
plot.lk381 <- plot_cropped_leaf(all.pl, "LK381")
plot.lk298 <- plot_cropped_leaf(all.pl, "LK298")
plot.lk388 <- plot_cropped_leaf(all.pl, "LK388")

drone.images <- plot_grid(plot.lk329, plot.lk235, plot.lk381, plot.lk298, 
                          plot.lk388, labels = c("A", "B", "C", "D", "E"), 
                          label_size = 16, nrow = 1)

lower.plot <- plot_grid(geneticgroups_leaf_greenness, stepwise_barplot, labels = c("F", "G"), label_size = 16, nrow = 1)

# Figure 6
combined_figure <- plot_grid(drone.images, lower.plot, nrow = 2)




