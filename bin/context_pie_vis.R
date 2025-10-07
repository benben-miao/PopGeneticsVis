#!/usr/bin/Rscript

# Author: Benben Miao
# Email: benben.miao@outlook.com
# GitHub: https://github.com/benben-miao
# R: v4.3.3

# 1. Library
suppressPackageStartupMessages({
	library(argparse)
})

# 2. Argument Parser
parser <- ArgumentParser(description = "Multi-Pie row: subjects; col: parts")

# 3. Add Argumnet
parser$add_argument(
	"--subjects_parts",
	type = "character",
	default = "./subjects_parts.txt",
	required = TRUE,
	help = "Row: subjects; Col: parts [./subjects_parts.txt]."
)
parser$add_argument(
	"--col_num",
	type = "integer",
	default = 3,
	required = FALSE,
	help = "Columns/Pies in row [3]."
)

# 4. Parse Args
args <- parser$parse_args()

subjects_parts <- args$subjects_parts
col_num <- args$col_num

## 1. Data
data <- read.table(
	subjects_parts,
	header = TRUE,
	sep = "\t",
	stringsAsFactors = FALSE
)

pie_num <- nrow(data)
pie_num_col <- col_num
pie_num_row <- ceiling(pie_num / pie_num_col)


parts <- colnames(data)[-1]
parts_num <- length(parts)

if (parts_num <= 8) {
	colors <- scales::brewer_pal(type = "qual", palette = "Set1")(parts_num)
} else {
	colors <- rainbow(parts_num)
}
names(colors) <- parts

## 2. Plot
plot <- function() {
	par(mfrow = c(pie_num_row, pie_num_col),
			mar = c(0, 2, 2, 0),
			oma = c(0, 0, 2, 0))
	# mtext("Multiple Pies Visualization", outer = TRUE, cex = 1.5, font = 2)

	for (i in 1:pie_num) {
		row_data <- as.numeric(data[i, -1])
		names(row_data) <- parts

		labels <- paste0(names(row_data), "\n", format(row_data, big.mark = ","))
		pie(row_data,
				labels = labels,
				main = data[i, 1],
				col = colors,
				cex = 1.0)
	}
}

## 3. Save
pdf(paste0(subjects_parts,".pdf"), width = 10, height = 6)
plot()
dev.off()

jpeg(paste0(subjects_parts,".jpeg"), units = "in", width = 10, height = 6, res = 300)
plot()
dev.off()
