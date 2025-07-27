### 03_pcoa_gene_expression_Lser_K9.R

### Author: Fiona Krammer  
### Description:  
### This script analyzes gene expression patterns across genetic groups (K = 9) in Lactuca serriola.  
### It processes log2 CPM expression data for leaf and root tissues, performs ANOVA to identify genes  
### with group-specific expression, and uses PCoA to visualize expression differences by region and group.  
### Log2 mean expression ratios are calculated and visualized in heatmaps, and PC-associated genes are  
### extracted and matched to Arabidopsis homologs for downstream GO enrichment analysis.  
### Additionally, genes with consistent group-specific up- or down-regulation are identified per tissue  
### and linked to Arabidopsis homologs.

# Load libraries
library(ggplot2)
library(readxl)
library(dplyr)
library(cowplot)
library(gplots)
library(RColorBrewer)
library(magick)

#-------------------------------------------------------------------------------
# LOAD FILES
#-------------------------------------------------------------------------------

cnt.cpm <- read.csv("./Data/CPM_SER_Tissues_Merged.csv", row.names = 1)
load("./Data/obj_gene.meta.out")
genetic.groups.9 <- read_excel("./Data/Lser200_pop_info_apr2025_long_form_changedregions.xlsx", sheet = "K9")
GO.file <- read.delim("./Data/20231004_Lactuca_sativa.annotation_overview_with_V8_positions.tsv")
load("./Data/obj_geneGO_ATH_GO_GOSLIM_20200501 1.out")
load("./Data/obj_GOnames_20200529 1.out")

#-------------------------------------------------------------------------------
# DATA PREPROCESSING
#-------------------------------------------------------------------------------

# Filter and clean metadata
genetic.groups.9 <- genetic.groups.9[, 1:(ncol(genetic.groups.9)-3)]
genetic.groups.9 <- genetic.groups.9[, !names(genetic.groups.9) %in% "noofg"]
genetic.groups.9 <- genetic.groups.9[!is.na(genetic.groups.9$region), ]

# Separate leaf and root samples
leaf.selection <- cnt.cpm[, grepl("leaf_", colnames(cnt.cpm))]
root.selection <- cnt.cpm[, grepl("root_", colnames(cnt.cpm))]

# Clean column names
colnames(leaf.selection) <- sub("^leaf_", "", colnames(leaf.selection))
colnames(root.selection) <- sub("^root_", "", colnames(root.selection))

# Log2 CPM transformation
log2_cpm_transform <- function(data, threshold = 20) {
  expressed_genes <- apply(data, 1, function(x) sum(x > 0)) > threshold
  filtered <- data[expressed_genes, ]
  log2((filtered + 1) / rowMeans(filtered, na.rm = TRUE))
}

log2.cpm.leaf <- log2_cpm_transform(leaf.selection)
log2.cpm.root <- log2_cpm_transform(root.selection)

#-------------------------------------------------------------------------------
# GROUP ASSIGNMENT
#-------------------------------------------------------------------------------

assign_groups <- function(data, metadata, expression_matrix) {
  metadata$LKID <- trimws(metadata$LKID)
  filtered_meta <- metadata[metadata$LKID %in% colnames(expression_matrix), ]
  group_map <- setNames(filtered_meta$ggroups, filtered_meta$LKID)
  groups <- group_map[colnames(expression_matrix)]
  groups[groups == 5] <- 2
  groups <- factor(groups)
  expression_matrix <- expression_matrix[, !is.na(groups)]
  groups <- groups[!is.na(groups)]
  stopifnot(length(groups) == ncol(expression_matrix))
  list(groups = groups, expression = expression_matrix)
}

leaf_info <- assign_groups(genetic.groups.9, genetic.groups.9, log2.cpm.leaf)
groups.leaf <- leaf_info$groups
log2.cpm.leaf <- leaf_info$expression

root_info <- assign_groups(genetic.groups.9, genetic.groups.9, log2.cpm.root)
groups.root <- root_info$groups
log2.cpm.root <- root_info$expression


#-------------------------------------------------------------------------------
# ANOVA ANALYSIS
#-------------------------------------------------------------------------------

getp <- function(expr, group) {
  aout <- anova(aov(expr ~ group))
  aout$`Pr(>F)`[1]
}

pvalues.leaf <- apply(log2.cpm.leaf, 1, function(x) getp(x, groups.leaf))
pvalues.root <- apply(log2.cpm.root, 1, function(x) getp(x, groups.root))

log10pvalues.leaf <- -log10(pvalues.leaf)
log10pvalues.root <- -log10(pvalues.root)

#-------------------------------------------------------------------------------
# PCOA PLOTTING
#-------------------------------------------------------------------------------

process_pcoa_logCPM <- function(cnt_cpm_lg, log10pvalues, groups, genetic_groups, title_suffix, output_file, threshold) {
  # Filter significant genes
  filtered_genes <- rownames(cnt_cpm_lg)[log10pvalues > threshold]
  filtered_expression <- cnt_cpm_lg[filtered_genes, ]
  
  # PCoA
  pcoa_result <- prcomp(t(filtered_expression), scale. = TRUE)
  LKIDs <- colnames(filtered_expression)
  
  # Metadata
  pcoa_df <- data.frame(
    PCo1 = pcoa_result$x[, 1],
    PCo2 = pcoa_result$x[, 2],
    LKID = LKIDs,
    Group = as.factor(groups[match(LKIDs, names(groups))]),
    Region = genetic_groups$region[match(LKIDs, genetic_groups$LKID)]
  )
  
  # Set Region colors
  region_colors <- c(
    "Turkey"       = "#E69F00",
    "Europe"       = "#56B4E9",
    "Caspian Sea"  = "#009E73",
    "Middle East"  = "#F0E442",
    "Central Asia" = "#CC79A7"
  )
  
  # Define shapes for 8 genetic groups
  group_shapes <- c(21, 22, 23, 24, 25, 7, 8, 21) # last shape repeated for 8th group (adjust if needed)
  group_levels <- levels(pcoa_df$Group)
  shape_mapping <- setNames(group_shapes[seq_along(group_levels)], group_levels)
  
  # Split data by filled vs hollow shapes
  filled_shapes <- c(21, 22, 23, 24, 25)
  pcoa_df$ShapeType <- ifelse(shape_mapping[as.character(pcoa_df$Group)] %in% filled_shapes, "filled", "hollow")
  pcoa_df$Shape <- shape_mapping[as.character(pcoa_df$Group)]
  
  # Plot
  library(ggplot2)
  
  pcoa_plot <- ggplot() +
    # Filled shapes layer
    geom_point(
      data = subset(pcoa_df, ShapeType == "filled"),
      aes(x = PCo1, y = PCo2, shape = Group, fill = Region),
      size = 5, color = "black", stroke = 0.6, alpha = 0.9
    ) +
    # Hollow shapes layer
    geom_point(
      data = subset(pcoa_df, ShapeType == "hollow"),
      aes(x = PCo1, y = PCo2, shape = Group, color = Region),
      size = 5, fill = NA, stroke = 1, alpha = 0.9
    ) +
    scale_shape_manual(values = shape_mapping) +
    scale_fill_manual(values = region_colors, na.translate = FALSE,
                      guide = guide_legend(override.aes = list(shape = 22, size = 5, color = NA))) +
    scale_color_manual(values = region_colors, na.translate = FALSE,
                       guide = "none") +
    theme_minimal(base_size = 14) +
    labs(
      x = "PC 1",
      y = "PC 2",
      fill = "Region",
      color = "Region",
      shape = "Genetic Group"
    ) +
    theme(
      legend.position = "right",
      axis.text = element_text(size = 13),
      axis.title = element_text(size = 16, face = "bold"),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 13),
      panel.grid = element_line(color = "gray90"),
      plot.title = element_blank()
    )
  
  # Save or return plot
  if (!is.null(output_file)) {
    ggsave(output_file, plot = pcoa_plot, width = 10, height = 8, dpi = 300)
  } else {
    return(pcoa_plot)
  }
}


plot_leaf <- process_pcoa_logCPM(log2.cpm.leaf, log10pvalues.leaf, groups.leaf, genetic.groups.9, NULL, NULL, 7)
plot_root <- process_pcoa_logCPM(log2.cpm.root, log10pvalues.root, groups.root, genetic.groups.9, NULL, NULL, 7)

combined_plot <- plot_grid(plot_leaf, plot_root, labels = c("A", "B"), label_size = 20, ncol = 2, align = "hv")
ggsave("Combined_PCoA_Leaf_Root.png", combined_plot, width = 16, height = 8, dpi = 300)

#-------------------------------------------------------------------------------
# LOG2 MEAN RATIO HEATMAPS
#-------------------------------------------------------------------------------

calculate_log_mean_ratio <- function(data, groups) {
  group_means <- sapply(levels(groups), function(g) {
    rowMeans(data[, groups == g, drop = FALSE], na.rm = TRUE)
  })
  log2((group_means + 1e-6) / (rowMeans(group_means, na.rm = TRUE) + 1e-6))
}

leaf.selection.lmr <- leaf.selection[rownames(leaf.selection) %in% rownames(log2.cpm.leaf), colnames(leaf.selection) %in% colnames(log2.cpm.leaf)]
root.selection.lmr <- root.selection[rownames(root.selection) %in% rownames(log2.cpm.root), colnames(root.selection) %in% colnames(log2.cpm.root)]

log_mean_ratio.leaf <- calculate_log_mean_ratio(leaf.selection.lmr, groups.leaf)
log_mean_ratio.root <- calculate_log_mean_ratio(root.selection.lmr, groups.root)

plot_heatmap <- function(data, log10p, groups, filename) {
  groups <- factor(groups)
  groups <- groups[colnames(data)]
  col_colors <- brewer.pal(length(levels(groups)), "Set2")[as.numeric(groups)]
  png(filename, width = 3000, height = 2000, res = 300)
  heatmap.2(data[log10p > 7, ],
            ColSideColors = col_colors,
            trace = "none",
            scale = "row",
            col = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
            cexCol = 1.2,
            cexRow = 1.0,
            key.title = "Log2 Ratio",
            key.xlab = "Z-score",
            margins = c(8, 10),
            density.info = "none")
  dev.off()
}

plot_heatmap(log_mean_ratio.leaf, log10pvalues.leaf, groups.leaf, "Heatmap_Leaf_Log2MeanRatio_SigPval7_K9Groups.png")
plot_heatmap(log_mean_ratio.root, log10pvalues.root, groups.root, "Heatmap_Root_Log2MeanRatio_SigPval7_K9Groups.png")

# Combine heatmaps
leaf_img <- image_read("Heatmap_Leaf_Log2MeanRatio_SigPval7_K9Groups.png")
root_img <- image_read("Heatmap_Root_Log2MeanRatio_SigPval7_K9Groups.png")
leaf_labeled <- image_annotate(leaf_img, "A", size = 60, gravity = "northwest", location = "+20+20")
root_labeled <- image_annotate(root_img, "B", size = 60, gravity = "northwest", location = "+20+20")
combined <- image_append(c(leaf_labeled, root_labeled), stack = TRUE)
image_write(combined, "combined_heatmaps_FigureS4_AB.png")

#-------------------------------------------------------------------------------
# ANALYSIS OF PC-ASSOCIATED GENES
#-------------------------------------------------------------------------------

extract_pc_genes <- function(pcoa_result, pc_index = 1, cutoff = 0.01) {
  loadings <- pcoa_result$rotation[, pc_index]
  names(loadings[abs(loadings) > cutoff])
}

pcoa_leaf <- prcomp(t(log2.cpm.leaf), scale. = TRUE)
pc1.leaf.genes <- extract_pc_genes(pcoa_leaf, 1)
pc2.leaf.genes <- extract_pc_genes(pcoa_leaf, 2)

pcoa_root <- prcomp(t(log2.cpm.root), scale. = TRUE)
pc1.root.genes <- extract_pc_genes(pcoa_root, 1)
pc2.root.genes <- extract_pc_genes(pcoa_root, 2)

extract_arabidopsis_homologs <- function(genes, GO.file) {
  homologs <- GO.file$best.hit.ARAPORT11..mmseqs2.[GO.file$gene.ID %in% genes]
  homologs <- homologs[homologs != ""]
  sub("\\.\\d+$", "", homologs)
}

PC1.leaf.arab <- extract_arabidopsis_homologs(pc1.leaf.genes, GO.file)
PC2.leaf.arab <- extract_arabidopsis_homologs(pc2.leaf.genes, GO.file)
PC1.root.arab <- extract_arabidopsis_homologs(pc1.root.genes, GO.file)
PC2.root.arab <- extract_arabidopsis_homologs(pc2.root.genes, GO.file)

# Arabidopsis homologs were uploaded on StringDB for GO enrichment analysis

#-------------------------------------------------------------------------------
# GROUP GENE EXPRESSION PATTERNS
#-------------------------------------------------------------------------------

logFC_threshold <- 0.5

group_genes <- function(data, pvals, threshold) {
  up <- apply(data[pvals > 5, ], 2, function(x) names(x[x > threshold]))
  down <- apply(data[pvals > 5, ], 2, function(x) names(x[x < -threshold]))
  list(up = up, down = down)
}

leaf_patterns <- group_genes(log_mean_ratio.leaf, log10pvalues.leaf, logFC_threshold)
root_patterns <- group_genes(log_mean_ratio.root, log10pvalues.root, logFC_threshold)

extract_homologs_per_group <- function(GO.file, up_genes, down_genes) {
  result <- list()
  for (i in 1:8) {
    gid <- as.character(i)
    get_homologs <- function(ids) {
      hits <- GO.file$best.hit.ARAPORT11..mmseqs2.[GO.file$gene.ID %in% ids[[gid]]]
      hits <- hits[hits != ""]
      unique(sub("\\.\\d+$", "", hits))
    }
    result[[paste0("g", gid, ".up")]] <- get_homologs(up_genes)
    result[[paste0("g", gid, ".down")]] <- get_homologs(down_genes)
  }
  result
}

shoot_gene_lists <- extract_homologs_per_group(GO.file, leaf_patterns$up, leaf_patterns$down)
root_gene_lists <- extract_homologs_per_group(GO.file, root_patterns$up, root_patterns$down)

# Merge gene lists
combine_gene_lists <- function(gene_list) {
  lapply(1:8, function(i) {
    gid <- as.character(i)
    c(gene_list[[paste0("g", gid, ".up")]], gene_list[[paste0("g", gid, ".down")]])
  }) |> setNames(paste0("g", 1:8))
}

shoot_merged <- combine_gene_lists(shoot_gene_lists)
root_merged <- combine_gene_lists(root_gene_lists)

# Example: root group 2 genes
root_merged$g2

# Merged Arabidopsis homologs were uploaded on StringDB for GO enrichment analysis