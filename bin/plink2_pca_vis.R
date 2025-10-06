#!/usr/bin/Rscript

# Author: Benben Miao
# Email: benben.miao@outlook.com
# GitHub: https://github.com/benben-miao
# R: v4.3.3

# 1. Library Packages
suppressPackageStartupMessages({
	library(argparse)
	library(ggplot2)
	library(dplyr)
})

# 2. Argument Parser
parser <- ArgumentParser(description = "PCA Visualization for PLINK PCA Eigenvec")

# 3. Add Argument
parser$add_argument(
	"--eigenvec",
	type = "character",
	default = "./pca.eigenvec",
	required = TRUE,
	help = "PLINK PCA results [pca.eigenvec]."
)
parser$add_argument(
	"--eigenval",
	type = "character",
	default = "./pca.eigenval",
	required = TRUE,
	help = "PLINK PCA results [pca.eigenval]."
)
parser$add_argument(
	"--sample_pop",
	type = "character",
	default = "./sample_pop.txt",
	required = TRUE,
	help = "Table with samples in col1, pops in col2 [sample_pop.txt]."
)
parser$add_argument(
	"--select_pops",
	type = "character",
	default = "GroupA,GroupB",
	required = TRUE,
	help = "Select populations ['GroupA,GroupB']."
)
parser$add_argument(
	"--x_pc",
	type = "integer",
	default = 1,
	required = FALSE,
	help = "Which PC for X Axis [1]."
)
parser$add_argument(
	"--y_pc",
	type = "integer",
	default = 2,
	required = FALSE,
	help = "Which PC for Y Axis [2]."
)
parser$add_argument(
	"--point_size",
	type = "numeric",
	default = 4,
	required = FALSE,
	help = "Points size [4]."
)
parser$add_argument(
	"--point_alpha",
	type = "numeric",
	default = 0.8,
	required = FALSE,
	help = "Points alpha [0.8]."
)

# 4. Parse Arguments
args <- parser$parse_args()
x_pc <- args$x_pc
y_pc <- args$y_pc

## 1. Read
vec <- read.table(
	args$eigenvec,
	comment.char = "",
	header = TRUE,
	row.names = 1,
	sep = "\t",
	stringsAsFactors = FALSE
)
val <- read.table(
	args$eigenval,
	header = FALSE,
	col.names = "Eigenvalue"
)
pop <- read.table(
	args$sample_pop,
	header = FALSE,
	sep = "\t",
	stringsAsFactors = FALSE,
	col.names = c("IID", "Population")
)

## 2. Data
vec_pop <- merge(vec, pop, by = "IID")

pop_select <- trimws(unlist(strsplit(args$select_pops, ",")))
vec_pop_select <- vec_pop %>%
	filter(Population %in% pop_select)

# Proportion of Variance Explained
pve <- round((val$Eigenvalue / sum(val$Eigenvalue)) * 100, 2)

## 3. Plot
p <- ggplot(
	data = vec_pop_select,
	aes(x = vec_pop_select[, x_pc + 1],
		y = vec_pop_select[, y_pc + 1],
		color = Population,
		fill = Population
	)) +
	geom_point(
		size = args$point_size,
		alpha = args$point_alpha) +
	stat_ellipse(
		geom = "polygon",
		position = "identity",
		type = "t", # c("t", "norm", "euclid")
		level = 0.95,
		segments = 100,
		na.rm = FALSE,
		show.legend = NA,
		inherit.aes = TRUE,
		alpha = 0.2
	) +
	labs(
		x = paste0("PC", x_pc, " (", pve[x_pc], "%)"),
		y = paste0("PC", y_pc, " (", pve[y_pc], "%)"),
	) +
	theme_light() +
	theme(
		panel.border = element_rect(colour = "black", linewidth = 1, fill = NA),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.ticks = element_line(colour = "black", linewidth = 1),
		axis.ticks.length = unit(0.2, "cm"),
		text = element_text(color = "black", size = 12),
		axis.text.x = element_text(color = "black", angle = 0, hjust = 0.5),
		axis.text.y = element_text(color = "black")
	)

## 4. Save
ggsave(
	paste(dirname(args$eigenvec), "/", paste(pop_select, collapse = "_"), "_PC_", x_pc, y_pc, ".pdf", sep = ""),
	plot = p,
	device = "pdf",
	units = "in",
	width = 10,
	height = 8,
	dpi = 300
)

ggsave(
	paste(dirname(args$eigenvec), "/", paste(pop_select, collapse = "_"), "_PC_", x_pc, y_pc, ".jpeg", sep = ""),
	plot = p,
	device = "jpeg",
	units = "in",
	width = 10,
	height = 8,
	dpi = 300
)
