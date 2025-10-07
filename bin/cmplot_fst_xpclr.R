#!/usr/bin/Rscript

# Author: Benben Miao
# Email: benben.miao@outlook.com
# GitHub: https://github.com/benben-miao
# R: v4.3.3

# 1. Library
suppressPackageStartupMessages({
	library(argparse)
	library(CMplot)
	library(dplyr)
	library(tidyr)
	library(stringr)
})

# 2. Argument Parser
parser <- ArgumentParser(description = "CMplot for FST or XPCLR results annotation")

# 3. Add Argumnet
parser$add_argument(
	"--method",
	type = "character",
	default = "FST",
	choices = c(
		"FST",
		"XPCLR"
	),
	required = FALSE,
	help = "SNP selection methods [FST]."
)
parser$add_argument(
	"--fst_xpclr",
	type = "character",
	default = "fst.anno",
	required = TRUE,
	help = "FST or XPCLR results with anno [fst.anno]."
)
parser$add_argument(
	"--plot_type",
	type = "character",
	default = "m",
	choices = c(
		"m",
		"c",
		"d",
		"q"
	),
	required = FALSE,
	help = "Plot type m(Manhattan), c(Circle), d(Density), q(QQ-plot) [m]."
)

# 4. Parse Args
args <- parser$parse_args()

fst_xpclr <- read.table(
	args$fst_xpclr,
	header = TRUE,
	sep = "\t",
	fill = TRUE,
	na.strings = "NA",
	stringsAsFactors = FALSE)

if (args$method == "FST") {
	top5_top1 <- quantile(
		fst_xpclr$WEIGHTED_FST,
		probs = c(0.95, 0.99),
		na.rm = TRUE
	)
	top5_value <- top5_top1[1]
	top1_value <- top5_top1[2]

	fst_top5 <- fst_xpclr %>%
		arrange(desc(WEIGHTED_FST)) %>%
		filter(WEIGHTED_FST >= top5_value) %>%
		mutate(id = paste0(CHROM, "_", BIN_START))
	# write.table(
	# 	fst_top5,
	# 	paste0(args$fst_xpclr, ".top5"),
	# 	append = FALSE,
	# 	quote = FALSE,
	# 	sep = "\t",
	# 	eol = "\n",
	# 	na = "NA",
	# 	dec = ".",
	# 	row.names = FALSE,
	# 	col.names = TRUE,
	# 	qmethod = "escape",
	# 	fileEncoding = ""
	# )

	fst_top5_gene <- fst_top5 %>%
		select(gene_id) %>%
		filter(!is.na(gene_id), str_trim(gene_id) != "") %>%
		mutate(gene_id = str_replace_all(gene_id, "\\(p\\)|\\(g\\)", "")) %>%
		separate_rows(gene_id, sep = ",") %>%
		mutate(gene_id = str_trim(gene_id)) %>%
		arrange(gene_id) %>%
		distinct(gene_id)
	# write.table(
	# 	fst_top5_gene,
	# 	paste0(args$fst_xpclr, ".top5.gene"),
	# 	append = FALSE,
	# 	quote = FALSE,
	# 	sep = "\t",
	# 	eol = "\n",
	# 	na = "NA",
	# 	dec = ".",
	# 	row.names = FALSE,
	# 	col.names = TRUE,
	# 	qmethod = "escape",
	# 	fileEncoding = ""
	# )

	data <- fst_top5 %>%
		mutate(CHROM = as.numeric(gsub("chr", "", CHROM))) %>%
		select(id, CHROM, BIN_START, WEIGHTED_FST)

	highlight <- fst_top5 %>%
		rowwise() %>%
		mutate(
			gene_id = gene_id %>%
				str_split(",\\s*") %>%
				unlist() %>%
				str_trim() %>%
				str_remove_all("\\(p\\)|\\(g\\)") %>%
				unique() %>%
				paste(collapse = ",")
		) %>%
		ungroup() %>%
		filter(!is.na(gene_id), gene_id != "") %>%
		slice_max(order_by = WEIGHTED_FST, n = 50) %>%
		select(id, gene_id)
} else if (args$method == "XPCLR") {
	top5_top1 <- quantile(
		fst_xpclr$xpclr_norm,
		probs = c(0.95, 0.99),
		na.rm = TRUE
	)
	top5_value <- top5_top1[1]
	top1_value <- top5_top1[2]

	xpclr_top5 <- fst_xpclr %>%
		arrange(desc(xpclr_norm)) %>%
		filter(xpclr_norm >= top5_value) %>%
		mutate(id = paste0(chrom, "_", start))
	# write.table(
	# 	xpclr_top5,
	# 	paste0(args$fst_xpclr, ".top5"),
	# 	append = FALSE,
	# 	quote = FALSE,
	# 	sep = "\t",
	# 	eol = "\n",
	# 	na = "NA",
	# 	dec = ".",
	# 	row.names = FALSE,
	# 	col.names = TRUE,
	# 	qmethod = "escape",
	# 	fileEncoding = ""
	# )

	xpclr_top5_gene <- xpclr_top5 %>%
		select(gene_id) %>%
		filter(!is.na(gene_id), str_trim(gene_id) != "") %>%
		mutate(gene_id = str_replace_all(gene_id, "\\(p\\)|\\(g\\)", "")) %>%
		separate_rows(gene_id, sep = ",") %>%
		mutate(gene_id = str_trim(gene_id)) %>%
		arrange(gene_id) %>%
		distinct(gene_id)
	# write.table(
	# 	xpclr_top5_gene,
	# 	paste0(args$fst_xpclr, ".top5.gene"),
	# 	append = FALSE,
	# 	quote = FALSE,
	# 	sep = "\t",
	# 	eol = "\n",
	# 	na = "NA",
	# 	dec = ".",
	# 	row.names = FALSE,
	# 	col.names = TRUE,
	# 	qmethod = "escape",
	# 	fileEncoding = ""
	# )

	data <- xpclr_top5 %>%
		mutate(chrom = as.numeric(gsub("chr", "", chrom))) %>%
		select(id, chrom, start, xpclr_norm)

	highlight <- xpclr_top5 %>%
		rowwise() %>%
		mutate(
			gene_id = gene_id %>%
				str_split(",\\s*") %>%
				unlist() %>%
				str_trim() %>%
				str_remove_all("\\(p\\)|\\(g\\)") %>%
				unique() %>%
				paste(collapse = ",")
		) %>%
		ungroup() %>%
		filter(!is.na(gene_id), gene_id != "") %>%
		slice_max(order_by = xpclr_norm, n = 50) %>%
		select(id, gene_id)
}

if (args$plot_type == "m") {
	width <- 10
	height <- 5
	threshold <- c(top5_value, top1_value)
} else if (args$plot_type == "c") {
	width <- 7
	height <- 7
	threshold <- NULL
} else if (args$plot_type == "q") {
	width <- 7
	height <- 7
	threshold <- NULL
} else if (args$plot_type == "d") {
	width <- 10
	height <- 7
	threshold <- NULL
}

pdf(
	paste(args$fst_xpclr, args$plot_type, "pdf", sep = "."),
	width = width,
	height = height
)
CMplot(
	data,
	LOG10 = FALSE,
	### chromosome
	bin.size = 1e6, # 1Mb
	band = 1,
	chr.labels = NULL,
	chr.border = FALSE,
	chr.pos.max = FALSE,
	chr.labels.angle = 0,
	chr.den.col = c("#0088ff", "#ff8800", "#ff0000"),
	cir.band = 0.1,
	cir.chr = TRUE,
	cir.chr.h = 0.5,
	cir.axis = TRUE,
	cir.axis.col = "black",
	cir.axis.grid = TRUE,
	### multitrack
	multraits = FALSE,
	multracks = FALSE,
	multracks.xaxis = FALSE,
	H = 1,
	# ylim = c(0.1, 0.8),
	### points
	pch = 16,
	type = "p", # "p" (point), "l" (cross line), "h" (vertical lines)
	col = rep(c("#aaaaaa55", "#33333355"), 9),
	# points.alpha = 100L,
	### circle
	plot.type = args$plot_type, # "m","c","q","d"
	cex = c(0.7, 1, 1), # circle, manhattan, qqplot
	r = 1, # circle
	outward = TRUE, # circle
	axis.cex = 1, # circle
	axis.lwd = 1, # circle
	lab.cex = 1.5, # circle
	lab.font = 1, # circle
	### manhattan
	# ylab = expression(-log[10](italic(p))), # manhattan and qqplot
	ylab = "WEIGHTED FST (Top 5%)",
	ylab.pos = 3, # manhattan
	xticks.pos = 1, # manhattan
	mar = c(0, 0, 0, 0), # manhattan, bottom, left, up, and right
	mar.between = 1, # manhattan
	### threshold
	threshold = threshold,
	threshold.col = c("#0000ff", "#ff0000"),
	threshold.lty = c(2, 2),
	threshold.lwd = 2,
	### signal
	amplify = FALSE,
	signal.cex = 2,
	signal.pch = 16,
	signal.line = 2,
	signal.col = NULL,
	### highlight
	highlight = highlight$id,
	highlight.text = as.vector(highlight$gene_id),
	highlight.cex = 1,
	highlight.pch = 16,
	highlight.type = "p",
	highlight.col = "#ff0000",
	highlight.text.font = 1, # 1: normal, 2: bold, 3: italic
	highlight.text.cex = 0.7,
	highlight.text.col = "#000000",
	### save
	conf.int = TRUE,
	conf.int.col = "white",
	file.output = FALSE,
	file.name = "file_name",
	file = "pdf", # "jpg", "pdf", "tiff", "png"
	box = FALSE,
	main = "",
	main.cex = 1,
	main.font = 1,
	legend.ncol = NULL,
	legend.cex = 1,
	legend.pos = "right", # "left","middle","right"
	dpi = 300,
	width = NULL,
	height = NULL,
	verbose = TRUE
)
dev.off()
