#!/usr/bin/Rscript

# Author: Benben Miao
# Email: benben.miao@outlook.com
# GitHub: https://github.com/benben-miao
# R: v4.3.3

# 1. Library
suppressPackageStartupMessages({
	library(argparse)
	library(data.table)
	library(dplyr)
	library(ggplot2)
})

# 2. Argument Parser
parser <- ArgumentParser(description = "Methyl Ratio Vis for SytosineReport")

# 3. Add Argumnet
parser$add_argument(
	"--cytosine_report",
	type = "character",
	default = "./cytosine_report.txt",
	required = TRUE,
	help = "Cytosine report from Bismark/MethylDackel [./cytosine_report.txt]."
)
parser$add_argument(
	"--sample_num",
	type = "integer",
	default = 1000000,
	required = FALSE,
	help = "Methyl site number for visualization [1000000]."
)

# 4. Parse Args
args <- parser$parse_args()

## 1. data.table -> fread
cytosine_report <- fread(
	args$cytosine_report,
	header = FALSE,
	sep = "\t",
	col.names = c("chr", "start", "strand", "c_count", "t_count", "context", "context3"),
	stringsAsFactors = FALSE
)

## 2. Data
#chrs <- paste0("chr", 1:18)
methyl_ratio <- cytosine_report %>%
# filter(chr %in% chrs) %>%
filter(c_count + t_count > 0) %>%
select(chr, start, c_count, t_count, context) %>%
mutate(ratio = c_count / (c_count + t_count)) %>%
mutate(context = if_else(context == "CG", "CpG", context)) %>%
mutate(context = factor(context, levels = c("CpG", "CHG", "CHH")))

## 3. Plot
p <- methyl_ratio %>%
filter(ratio > 0) %>%
sample_n(args$sample_num) %>%
ggplot(aes(x = ratio, color = context)) +
geom_density() +
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
	paste0(args$cytosine_report, ".methratio.pdf"),
	plot = p,
	device = "pdf",
	units = "in",
	width = 10,
	height = 6,
	dpi = 300
)

ggsave(
	paste0(args$cytosine_report, ".methratio.jpeg"),
	plot = p,
	device = "jpeg",
	units = "in",
	width = 10,
	height = 6,
	dpi = 300
)
