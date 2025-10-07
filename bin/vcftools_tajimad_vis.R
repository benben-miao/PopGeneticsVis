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
	library(ggpubr)
})

# 2. Argument Parser
parser <- ArgumentParser(description = "VCFtools TajimaD Violin Box Vis")

# 3. Add Argument
parser$add_argument(
	"--tajimad_win",
	type = "character",
	default = "./all.Tajima.D",
	required = TRUE,
	help = "VCFtools TajimaD with ordered POP [all.Tajima.D]."
)
parser$add_argument(
	"--select_pops",
	type = "character",
	default = "GroupA,GroupB",
	required = TRUE,
	help = "Select populations ['GroupA,GroupB']."
)
parser$add_argument(
	"--box_outliers",
	type = "logical",
	default = FALSE,
	required = FALSE,
	help = "Show box outliers uppercase [FALSE]."
)

# 4. Parse Arguments
args <- parser$parse_args()

## 1. Sample Pop
pop_select <- trimws(unlist(strsplit(args$select_pops, ",")))

# 4. Read PI
tajimad_win <- read.table(
	args$tajimad_win,
	header = TRUE,
	comment.char = "",
	sep = "\t",
	stringsAsFactors = FALSE) %>%
	filter(POP %in% pop_select) %>%
	mutate(
		TajimaD = as.numeric(TajimaD),
		POP = factor(POP, levels = unique(POP))
	) %>%
	filter(!is.na(TajimaD) & is.finite(TajimaD))

# 5. Plot
p <- ggplot(tajimad_win, aes(x = POP, y = TajimaD, fill = POP, color = POP)) +
	geom_violin(
		stat = "ydensity",
		position = "dodge",
		# draw_quantiles = NULL,
		trim = FALSE,
		bounds = c(-Inf, Inf),
		scale = "count",
		na.rm = FALSE,
		orientation = NA,
		show.legend = TRUE,
		inherit.aes = TRUE,
		width = 0.8,
		alpha = 0.3,
		linetype = 1,
		linewidth = 1,
		colour = "#00000000") +
	geom_boxplot(
		stat = "boxplot",
		position = "dodge2",
		outliers = args$box_outliers,
		outlier.colour = NULL,
		outlier.color = "#00000088",
		outlier.fill = NULL,
		outlier.shape = 16,
		outlier.size = 0.5,
		outlier.stroke = 0.5,
		outlier.alpha = 0.1,
		notch = FALSE,
		notchwidth = 0.5,
		staplewidth = 0,
		varwidth = FALSE,
		na.rm = FALSE,
		orientation = NA,
		show.legend = FALSE,
		inherit.aes = TRUE,
		width = 0.2,
		alpha = 0.5,
		linewidth = 1) +
	# geom_jitter(
	# 	stat = "identity",
	# 	# position = "jitter",
	# 	width = 0.1,
	# 	height = NULL,
	# 	na.rm = FALSE,
	# 	show.legend = NA,
	# 	inherit.aes = TRUE,
	# 	size = args$point_size,
	# 	alpha = 0.5) +
	labs(x = "Population", y = "TajimaD") +
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

n_pops <- length(unique(tajimad_win$POP))
if (n_pops == 2) {
	p <- p +
		stat_compare_means(
			method = "wilcox.test", # c("t.test", "wilcox.test")
			paired = FALSE,
			method.args = list(),
			ref.group = NULL,
			comparisons = NULL,
			hide.ns = FALSE,
			label.sep = ", ",
			label = "p.format", # c("p.signif", "p.format")
			label.x.npc = "left",
			label.y.npc = "top",
			label.x = NULL,
			label.y = NULL,
			vjust = 0,
			tip.length = 0.03,
			bracket.size = 0.3,
			step.increase = 0.1,
			symnum.args = list(),
			geom = "text",
			position = "identity",
			na.rm = FALSE,
			show.legend = NA,
			inherit.aes = TRUE
		)
} else if (n_pops > 2) {
	p <- p +
		stat_compare_means(
			method = "wilcox.test", # c("t.test", "wilcox.test")
			paired = FALSE,
			method.args = list(),
			ref.group = NULL,
			comparisons = combn(unique(as.character(tajimad_win$POP)), 2, simplify = FALSE),
			hide.ns = FALSE,
			label.sep = ", ",
			label = "p.format", # c("p.signif", "p.format")
			label.x.npc = "left",
			label.y.npc = "top",
			label.x = NULL,
			label.y = NULL,
			vjust = 0,
			tip.length = 0.03,
			bracket.size = 0.3,
			step.increase = 0.1,
			symnum.args = list(),
			geom = "text",
			position = "identity",
			na.rm = FALSE,
			show.legend = NA,
			inherit.aes = TRUE
		)
}

## 4. Save
ggsave(
	paste(dirname(args$tajimad_win), "/",
		  paste(pop_select, collapse = "_"),
		  "_TajimaD", ".pdf",
		  sep = ""),
	plot = p,
	device = "pdf",
	units = "in",
	width = 10,
	height = 5,
	dpi = 300
)

ggsave(
	paste(
		dirname(args$tajimad_win), "/",
		paste(pop_select, collapse = "_"),
		"_TajimaD", ".jpeg",
		sep = ""),
	plot = p,
	device = "jpeg",
	units = "in",
	width = 10,
	height = 5,
	dpi = 300
)
