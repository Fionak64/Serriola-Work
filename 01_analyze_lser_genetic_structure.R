### 01_analyze_lser_genetic_structure.R
### Author: Basten Snoek (adapted by Fiona Krammer)
### Description: This script performs hierarchical clustering on SNP data of Lactuca serriola genotypes
###              and visualizes genetic groupings with annotated region and origin.

# Load packages
library(ggplot2)
library(viridis)
library(openxlsx)
library(cowplot)
library(ggrepel)
library(ggdendro)

# ---------------------- Load and prepare data ----------------------

# Load SNP genotype object
load("./Data/obj_ldb.all_Lsatser200_no_LK283_20250120.out")

# Load sample metadata
satser.meta <- read.xlsx("./Data/Satser_line_meta_for_ldb.xlsx")

# Convert genotype strings to numeric matrix
ll.mat <- sapply(strsplit(ldb.all$use.id, ""), as.numeric)

# Select *L. serriola* samples (exclude suspected hybrids and errors)
ser.mat <- t(ll.mat[c(199:280, 282:396), ])
colnames(ser.mat) <- satser.meta$LKID[c(199:280, 282:396)]

not.use <- c("LK201","LK202","LK203","LK204","LK206","LK207","LK208","LK212",
             "LK223","LK228","LK251","LK256","LK258","LK306","LK373","LK377",
             "LK378","LK379","LK383","LK205","LK220","LK263","LK397","LK398")
ser.mat <- ser.mat[, !colnames(ser.mat) %in% not.use]
use.lk <- colnames(ser.mat)

# ---------------------- Region grouping ----------------------

# Define countries into regions
Central_Asia <- c("UZB","AFG","KGZ")
Kasp_Sea     <- c("ARM","AZE","GEO","RUS","IRN")
NW_Europa    <- c("BEL","DEU","DNK","GBR","NLD","NOR","SWE")
East_Europa  <- c("CZE","HUN","POL","SVK","ROU","UKR","GRC","BGR")
Middle_Eastern <- c("EGY","IRQ","ISR","KEN","SYR")
SW_Europa    <- c("ESP","FRA","ITA","MLT","CHE","PRT","HRV","SVN")
Black_Sea    <- c("TUR")

# Assign region labels
selc <- satser.meta$LKID %in% use.lk
country <- satser.meta$Country[selc]
region <- rep(NA, length(country))
region[country %in% Central_Asia]   <- 'Central Asia'
region[country %in% Kasp_Sea]       <- 'Caspian Sea'
region[country %in% NW_Europa]      <- 'Europe'
region[country %in% East_Europa]    <- 'Europe'
region[country %in% Middle_Eastern] <- 'Middle East'
region[country %in% SW_Europa]      <- 'Europe'
region[country %in% Black_Sea]      <- 'Turkey'
names(region) <- satser.meta$LKID[selc]

# ---------------------- Clustering ----------------------

ser.cov <- cov(ser.mat)
hc.ser <- hclust(dist(ser.cov))
ctr.ser <- cutree(hc.ser, 2:20)

# ---------------------- Dendrogram with regions ----------------------

ddata <- dendro_data(hc.ser, type = "rectangle")
ok <- label(ddata)
region <- region[ok$label]

ggd <- ggplot(segment(ddata)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.3) +
  annotate("text", x = c(182, 182), y = c(1, -4), label = c("Group", "Region")) +
  annotate("label", x = c(170,155,130,75,10), y = c(-5,-7,-5,-5,-5),
           label = c("Central Asia", "Caspian Sea", "Middle East", "Europe", "Turkey"),
           angle = 90, col = c("green3","cyan","red4","black","gold3")) +
  geom_text(data = label(ddata),
            aes(x = x, y = y - 2, label = satser.meta[ok$label,]$Country),
            hjust = 1, vjust = 0.5, size = 1.5) +
  geom_point(data = ok,
             aes(x = x, y = y, color = as.factor(ctr.ser[ok$label, 8])),
             size = 2, shape = 15) +
  geom_point(data = ok,
             aes(x = x, y = y - 3, fill = region),
             size = 2, shape = 22) +
  scale_color_manual(values = c("black", "red", "red4", "gold", "firebrick", "cyan3", "green3", "green1", "green4")) +
  scale_fill_manual(values = c("cyan", "green3", "black", "red4", "gold2")) +
  guides(color = guide_legend(title = "Group", override.aes = list(size = 6)),
         fill = guide_legend(title = "Region", override.aes = list(size = 6))) +
  coord_flip() +
  scale_y_reverse(expand = c(0.2, 0)) +
  theme_dendro() +
  theme(legend.background = element_rect(color = "grey", linewidth = 0.2),
        legend.text = element_text(size = 12))

# ---------------------- Reranking and color mapping ----------------------

rerank <- function(x) {
  as.numeric(rank(-table(x), ties.method = "random")[as.character(x)])
}
rr.ser <- as.numeric(apply(ctr.ser, 2, rerank))

tree.mat <- matrix(NA, ncol = 19, nrow = nrow(ctr.ser))
tree.mat[, 19] <- ctr.ser[, 19]

for (i in 1:18) {
  tab <- table(ctr.ser[, i], ctr.ser[, 19])
  tree.mat[, i] <- apply(tab, 1, which.max)[ctr.ser[, i]]
}
tree.ser <- as.numeric(tree.mat)

# Define color palette
selc.col <- order(-table(tree.ser))
use.col <- c("black","grey","seagreen1","red1","red4","green1","green4","blue1",
             "orange4","gold","purple","cyan1","cyan4",turbo(8)[2:8])
use.col[selc.col] <- use.col
use.col[1] <- "magenta"

# ---------------------- Final heatmap figure ----------------------

LKID <- colnames(ser.mat)
country <- satser.meta$Country[satser.meta$LK_ID %in% LKID]
LK.country <- apply(cbind(LKID, country), 1, paste, collapse = "_")
region <- region[LKID]

to.pl <- data.frame(
  noofg = rep(2:20, each = length(LKID)),
  ggroups = as.numeric(ctr.ser),
  LKID,
  region,
  rr.ser,
  tree.ser,
  LK.country
)

to.pl$LKID <- factor(to.pl$LKID, levels = LKID[hc.ser$order])
to.pl$LK.country <- factor(to.pl$LK.country, levels = LK.country[hc.ser$order])

ser.gen.pop.fig <- ggplot(to.pl) +
  geom_tile(aes(LK.country, noofg, fill = as.factor(tree.ser)), col = "white") +
  scale_fill_manual(values = use.col) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(. ~ region, space = "free_x", scale = "free_x") +
  xlab("Genotype") + ylab("Number of clusters") +
  labs(fill = "Genetic Groups") +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16),
    strip.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.3, "cm")
  )

# Save figure
ggsave(
  filename = "./Figures/ser_gen_pop_fig_highres.png",
  plot = ser.gen.pop.fig,
  width = 14,
  height = 6,
  dpi = 600
)

# Export long-format cluster table
write.xlsx(to.pl, file = "./Outputs/Lser200_pop_info_apr2025_long_form.xlsx")
