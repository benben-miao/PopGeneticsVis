#!/usr/bin/Rscript

# Author: Benben Miao
# Email: benben.miao@outlook.com
# GitHub: https://github.com/benben-miao
# R: v4.3.3

# 1. Library
suppressPackageStartupMessages({
	library(argparse)
	library(dplyr)
	library(ggplot2)
})

# 2. Argument Parser
parser <- ArgumentParser(description = "Meth.Diff Hyper & Hypo Violin Plot")

# 3. Add Argumnet
parser$add_argument(
	"--methylkit_dmr",
	type = "character",
	default = "./methylkit_dmr.txt",
	required = TRUE,
	help = "MethylKit DMR results (chr, start, end, ..., meth.diff) [./methylkit_dmr.txt]."
)
parser$add_argument(
	"--violin_scale",
	type = "character",
	default = "count",
	choices = c(
		"area",
		"count",
		"width"
	),
	required = FALSE,
	help = "Violin area and width [count]."
)
parser$add_argument(
	"--violin_width",
	type = "numeric",
	default = 1,
	required = FALSE,
	help = "Violin width percent [1]."
)
parser$add_argument(
	"--point_size",
	type = "numeric",
	default = 1.5,
	required = FALSE,
	help = "Point size [1.5]."
)

# 4. Parse Args
args <- parser$parse_args()

## 1. Read
dmr <- read.table(
	args$methylkit_dmr,
	header = TRUE,
	sep = "\t",
	stringsAsFactors = FALSE
)

## 2. Data
chrs <- paste0("chr", 1:18)
data <- dmr %>%
	filter(chr %in% chrs) %>%
	mutate(chr = factor(chr, levels = chrs))

hyper_hypo <- data %>%
	group_by(chr) %>%
	summarise(
		n_positive = sum(meth.diff > 0, na.rm = TRUE),
		n_negative = sum(meth.diff < 0, na.rm = TRUE),
		max_y = max(meth.diff, na.rm = TRUE),
		.groups = 'drop'
	) %>%
	mutate(
		chr = factor(chr, levels = chrs),
		y_pos = max_y + (max(data$meth.diff) - min(data$meth.diff)) * 0.05,
		label = sprintf("+%d / -%d", n_positive, n_negative)
	)

## 3. Plot
p <- ggplot(data, aes(x = chr, y = meth.diff)) +
	geom_violin(
		stat = "ydensity",
		position = "dodge",
		draw_quantiles = NULL,
		trim = TRUE,
		bounds = c(-Inf, Inf),
		scale = args$violin_scale,
		na.rm = FALSE,
		orientation = NA,
		show.legend = NA,
		inherit.aes = TRUE,
		width = args$violin_width,
		fill = "white",
		color = "black") +
	geom_point(
		aes(color = meth.diff > 0),
		position = position_jitter(width = 0.2, height = 0),
		size = args$point_size) +
	scale_color_manual(
		values = c(`TRUE` = "#ff000088", `FALSE` = "#00880088"),
		name = "Diff",
		labels = c("Hypo", "Hyper")) +
	geom_text(
		data = hyper_hypo,
		aes(x = chr, y = 0, label = label),
		size = 3,
		hjust = 0.5,
		vjust = -1,
		angle = -90,
		color = "black") +
	scale_x_discrete(drop = FALSE) +
	labs(x = "Chromosome", y = "Meth Diff (Hyper/Hypo)") +
	theme_light() +
	theme(
		panel.border = element_rect(colour = "black", linewidth = 1, fill = NA),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.ticks = element_line(colour = "black", linewidth = 1),
		axis.ticks.length = unit(0.2, "cm"),
		text = element_text(color = "black", size = 12),
		axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
		axis.text.y = element_text(color = "black"),
	)

## 4. Save
ggsave(
	paste0(args$methylkit_dmr, ".violin.pdf"),
	plot = p,
	device = "pdf",
	units = "in",
	width = 10,
	height = 6,
	dpi = 300
)

ggsave(
	paste0(args$methylkit_dmr, ".violin.jpeg"),
	plot = p,
	device = "jpeg",
	units = "in",
	width = 10,
	height = 6,
	dpi = 300
)
