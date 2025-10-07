#!/usr/bin/Rscript

# Author: Benben Miao
# Email: benben.miao@outlook.com
# GitHub: https://github.com/benben-miao
# R: v4.3.3

# 1. Library Packages
suppressPackageStartupMessages({
	library(argparse)
	library(methylKit)
	library(GenomicFeatures)
	library(rtracklayer)
	library(GenomicRanges)
	library(dplyr)
})

# 2. ArgumentParser
parser <- ArgumentParser(description = "MethylKit Analysis Workflow")

# 3. Add Argumnet
parser$add_argument(
	"--control",
	type = "character",
	required = TRUE,
	help = "Control group name."
)
parser$add_argument(
	"--treated",
	type = "character",
	required = TRUE,
	help = "Treated group name."
)
parser$add_argument(
	"--in_dir",
	type = "character",
	default = "./",
	required = TRUE,
	help = "Input directory [./]."
)
parser$add_argument(
	"--out_dir",
	type = "character",
	default = "./",
	required = TRUE,
	help = "Output directory [./]."
)
parser$add_argument(
	"--threads",
	type = "integer",
	default = 1,
	help = "Threads, more CPUs more Mems need [1]."
)
parser$add_argument(
	"--nrep",
	type = "integer",
	default = 3,
	help = "Samples replicates each group [3]."
)
parser$add_argument(
	"--level",
	type = "character",
	default = "DMR",
	choices = c(
		"DMC", 
		"DMR", 
		"DMP"
	),
	help = "Methylation level [DMR]."
)
parser$add_argument(
	"--context",
	type = "character",
	default = "CpG",
	choices = c(
		"CpG", 
		"CHG", 
		"CHH"
	),
	help = "Methylation context [CpG]."
)
parser$add_argument(
	"--mincov",
	type = "integer",
	default = 10,
	help = "Minimum coverage [10]."
)
parser$add_argument(
	"--win_size",
	type = "integer",
	default = 2000,
	help = "Window size for DMR (bp) [2000]."
)
parser$add_argument(
	"--step_size",
	type = "integer",
	default = 2000,
	help = "Step size for DMR (bp) [2000]."
)
parser$add_argument(
	"--corr_method",
	type = "character",
	default = "pearson",
	choices = c(
		"pearson", 
		"kendall", 
		"spearman"
	),
	help = "Correlation method [pearson]."
)
parser$add_argument(
	"--adjust_method",
	type = "character",
	default = "SLIM",
	choices = c(
		"SLIM",
		"holm",
		"hochberg",
		"hommel",
		"bonferroni",
		"BH",
		"BY",
		"fdr",
		"none",
		"qvalue"
	),
	help = "P-value adjustment method [SLIM]."
)
parser$add_argument(
	"--test_method",
	type = "character",
	default = "Chisq",
	choices = c(
		"F", 
		"Chisq", 
		"fast.fisher", 
		"midPval"
	),
	help = "Statistical test method [Chisq]."
)
parser$add_argument(
	"--qvalue_cutoff",
	type = "double",
	default = 0.05,
	help = "Q-value cutoff [0.05]."
)
parser$add_argument(
	"--meth_cutoff",
	type = "double",
	default = 25,
	help = "Methylation difference cutoff [25]."
)
parser$add_argument(
	"--gtf",
	type = "character",
	default = NULL,
	help = "GTF/GFF of genome [NULL]."
)
parser$add_argument(
	"--promoter_up",
	type = "integer",
	default = 2000,
	help = "Promoter upstream distance (bp) [2000]"
)
parser$add_argument(
	"--promoter_down",
	type = "integer",
	default = 200,
	help = "Promoter downstream distance (bp) [200]"
)
parser$add_argument(
	"--plot_width",
	type = "double",
	default = 10.00,
	help = "Plot width (inch) [10.00]."
)
parser$add_argument(
	"--plot_height",
	type = "double",
	default = 6.18,
	help = "Plot height (inch) [10.00]."
)

# 4. Args
args <- parser$parse_args()

if (!dir.exists(args$out_dir))
	dir.create(args$out_dir, recursive = TRUE)

## 1. Files
cytosine_files <- file.path(
	args$in_dir, 
	paste0(
		c(rep(args$control, args$nrep), rep(args$treated, args$nrep)), 
		"_", 
		rep(1:args$nrep, 2), 
		".cytosine_report.txt")
)
stopifnot(all(file.exists(cytosine_files)))
cytosine_list <- as.list(cytosine_files)

## 2. methRead
meth_obj <- methRead(
	location = cytosine_list,
	sample.id = as.list(c(
		paste0(args$control, "_", 1:args$nrep), paste0(args$treated, "_", 1:args$nrep)
	)),
	assembly = "HALdis.HiFi.Genome.hapX_Mitochondrion.Chrom",
	dbtype = NA,
	pipeline = "bismarkCytosineReport",
	header = TRUE,
	sep = "\t",
	context = args$context,
	resolution = "base",
	treatment = c(rep(0, args$nrep), rep(1, args$nrep)),
	mincov = args$mincov
)

## 3. filterByCoverage
meth_filt <- filterByCoverage(
	meth_obj, 
	lo.count = args$mincov,
	lo.perc = NULL,
	hi.count = NULL,
	hi.perc = 99.9,
	chunk.size = 1e+06,
	save.db = FALSE
)

## 4. normalizeCoverage
meth_norm <- normalizeCoverage(
	meth_filt, 
	method = "median",
	chunk.size = 1e+06,
	save.db = FALSE
)

## 5. txdb
txdb <- NULL
if (!is.null(args$gtf)){
	if (!file.exists(args$gtf)) stop("GTF not found: ", args$gtf)
	txdb <- makeTxDbFromGFF(args$gtf, format = "gtf")
}

annotate_regions <- function(
	df,
	txdb,
	promoter_up = 2000,
	promoter_down = 200){
	
	if (is.null(txdb)){
		df$anno_genes <- NA
		df$anno_promoters <- NA
		return(df)
	}

	if (inherits(df, "methylDiff")) {
		df <- as.data.frame(df)
	}
	
	genes_gr <- genes(txdb)
	promoters_gr <- promoters(
		txdb,
		upstream = promoter_up,
		downstream = promoter_down
	)

	gr_df <- GRanges(
		seqnames = df$chr,
		ranges = IRanges(start = as.integer(df$start), 
		end = as.integer(df$end)),
		strand = "*"
	)

	ol_genes <- findOverlaps(gr_df, genes_gr, ignore.strand = TRUE)
	genes_col <- rep(NA, length(gr_df))
	genes_col[queryHits(ol_genes)] <- mcols(genes_gr)$gene_id[subjectHits(ol_genes)]

  ol_promoters <- findOverlaps(gr_df, promoters_gr)
  promoters_col <- rep(NA, length(gr_df))
  promoters_col[queryHits(ol_promoters)] <- mcols(promoters_gr)$tx_name[subjectHits(ol_promoters)]
  
  df$anno_genes <- genes_col
  df$anno_promoters <- promoters_col
  
  return(df)
}

## 6. Methyl level
unite_plot <- NULL
diff_plot <- NULL

if (args$level == "DMC"){
	dmc_unit <- unite(
		meth_norm,
		destrand = FALSE,
		min.per.group = 1L,
		chunk.size = 1e+06,
		mc.cores = args$threads,
		save.db = FALSE
	)
	unite_plot <- dmc_unit

	dmc_diff <- calculateDiffMeth(
		dmc_unit,
		overdispersion = "none",
		adjust = args$adjust_method,
		effect = "wmean",
		test = args$test_method,
		mc.cores = args$threads,
		slim = TRUE,
		weighted.mean = TRUE,
		chunk.size = 1e+06,
		save.db = FALSE
	)
	diff_plot <- dmc_diff
	
	dmc_diff_df <- as.data.frame(dmc_diff)
	write.table(
		dmc_diff_df,
		file.path(
			args$out_dir,
			paste(args$control, "vs", args$treated, args$level, args$context, "base.txt", sep = "_")),
		sep = "\t",
		quote = FALSE,
		row.names = FALSE
	)

	dmc_diff_sig <- getMethylDiff(
		dmc_diff,
		difference = args$meth_cutoff,
		qvalue = args$qvalue_cutoff,
		type = "all",
		chunk.size = 1e+06,
		save.db = FALSE
	)
	dmc_diff_sig_df <- as.data.frame(dmc_diff_sig)
	write.table(
		dmc_diff_sig_df,
		file.path(
			args$out_dir,
			paste(args$control, "vs", args$treated, args$level, args$context, "base_Sig.txt", sep = "_")),
		sep = "\t",
		quote = FALSE,
		row.names = FALSE
	)

	dmc_diff_sig_hyper <- getMethylDiff(
		dmc_diff,
		difference = args$meth_cutoff,
		qvalue = args$qvalue_cutoff,
		type = "hyper",
		chunk.size = 1e+06,
		save.db = FALSE
	)
	write.table(
		dmc_diff_sig_hyper,
		file.path(
			args$out_dir,
			paste(args$control, "vs", args$treated, args$level, args$context, "base_Sig_Hyper.txt", sep = "_")),
		sep = "\t",
		quote = FALSE,
		row.names = FALSE
	)
	
	dmc_diff_sig_hypo <- getMethylDiff(
		dmc_diff,
		difference = args$meth_cutoff,
		qvalue = args$qvalue_cutoff,
		type = "hypo",
		chunk.size = 1e+06,
		save.db = FALSE
	)
	write.table(
		dmc_diff_sig_hypo,
		file.path(
			args$out_dir,
			paste(args$control, "vs", args$treated, args$level, args$context, "base_Sig_Hypo.txt", sep = "_")),
		sep = "\t",
		quote = FALSE,
		row.names = FALSE
	)

	if (!is.null(txdb)){
		dmc_sig_df_anno <- annotate_regions(
			dmc_diff_sig_df,
			txdb,
			promoter_up = args$promoter_up,
			promoter_down = args$promoter_down
		)
		write.table(
			dmc_sig_df_anno,
			file.path(
				args$out_dir,
				paste(args$control, "vs", args$treated, args$level, args$context, "base_Sig_Anno.txt", sep = "_")),
			sep = "\t",
			quote = FALSE,
			row.names = FALSE
		)
	}
}

if (args$level == "DMR"){	
	dmr_tile <- tileMethylCounts(
		meth_norm,
		win.size = args$win_size,
		step.size = args$step_size,
		cov.bases = 0,
		mc.cores = args$threads,
		save.db = FALSE
	)
	
	dmr_unite <- unite(
		dmr_tile,
		destrand = FALSE,
		min.per.group = 1L,
		chunk.size = 1e+06,
		mc.cores = args$threads,
		save.db = FALSE
	)
	unite_plot <- dmr_unite

	dmr_diff <- calculateDiffMeth(
		dmr_unite,
		overdispersion = "none",
		adjust = args$adjust_method,
		effect = "wmean",
		test = args$test_method,
		mc.cores = args$threads,
		slim = TRUE,
		weighted.mean = TRUE,
		chunk.size = 1e+06,
		save.db = FALSE
	)
	diff_plot <- dmr_diff
	
	dmr_diff_df <- as.data.frame(dmr_diff)
	write.table(
		dmr_diff_df,
		file.path(
			args$out_dir,
			paste(args$control, "vs", args$treated, args$level, args$context, "2k.txt", sep = "_")),
		sep = "\t",
		quote = FALSE,
		row.names = FALSE
	)

	dmr_diff_sig <- getMethylDiff(
		dmr_diff,
		difference = args$meth_cutoff,
		qvalue = args$qvalue_cutoff,
		type = "all",
		chunk.size = 1e+06,
		save.db = FALSE
	)
	dmr_diff_sig_df <- as.data.frame(dmr_diff_sig)
	write.table(
		dmr_diff_sig_df,
		file.path(
			args$out_dir,
			paste(args$control, "vs", args$treated, args$level, args$context, "2k_Sig.txt", sep = "_")),
		sep = "\t",
		quote = FALSE,
		row.names = FALSE
	)

	dmr_diff_sig_hyper <- getMethylDiff(
		dmr_diff,
		difference = args$meth_cutoff,
		qvalue = args$qvalue_cutoff,
		type = "hyper",
		chunk.size = 1e+06,
		save.db = FALSE
	)
	write.table(
		dmr_diff_sig_hyper,
		file.path(
			args$out_dir,
			paste(args$control, "vs", args$treated, args$level, args$context, "2k_Sig_Hyper.txt", sep = "_")),
		sep = "\t",
		quote = FALSE,
		row.names = FALSE
	)
	
	dmr_diff_sig_hypo <- getMethylDiff(
		dmr_diff,
		difference = args$meth_cutoff,
		qvalue = args$qvalue_cutoff,
		type = "hypo",
		chunk.size = 1e+06,
		save.db = FALSE
	)
	write.table(
		dmr_diff_sig_hypo,
		file.path(
			args$out_dir,
			paste(args$control, "vs", args$treated, args$level, args$context, "2k_Sig_Hypo.txt", sep = "_")),
		sep = "\t",
		quote = FALSE,
		row.names = FALSE
	)

	if (!is.null(txdb)){
		dmr_sig_df_anno <- annotate_regions(
			dmr_diff_sig_df,
			txdb,
			promoter_up = args$promoter_up,
			promoter_down = args$promoter_down
		)
		write.table(
			dmr_sig_df_anno,
			file.path(
				args$out_dir,
				paste(args$control, "vs", args$treated, args$level, args$context, "2k_Sig_Anno.txt", sep = "_")),
			sep = "\t",
			quote = FALSE,
			row.names = FALSE
		)
	}
}

if (args$level == "DMP"){
	if (is.null(txdb)) stop("--gtf must be provided for dmp analysis")
	
	promoters_gr <- promoters(
		txdb,
		upstream = args$promoter_up,
		downstream = args$promoter_down
	)

	region_counts <- regionCounts(
		meth_norm,
		regions = promoters_gr,
		save.db = FALSE
	)
	
	dmp_unite <- unite(
		region_counts,
		destrand = FALSE,
		min.per.group = 1L,
		chunk.size = 1e6,
		mc.cores = args$threads,
		save.db = FALSE
	)
	unite_plot <- dmp_unite

	dmp_diff <- calculateDiffMeth(
		dmp_unite,
		overdispersion = "none",
		adjust = args$adjust_method,
		effect = "wmean",
		test = args$test_method,
		mc.cores = args$threads,
		slim = TRUE,
		weighted.mean = TRUE,
		chunk.size = 1e+06,
		save.db = FALSE
	)
	diff_plot <- dmp_diff
	
	dmp_diff_df <- as.data.frame(dmp_diff)
	write.table(
		dmp_diff_df,
		file.path(
			args$out_dir,
			paste(args$control, "vs", args$treated, args$level, args$context, "2k.txt", sep = "_")),
		sep = "\t",
		quote = FALSE,
		row.names = FALSE
	)

	dmp_diff_sig <- getMethylDiff(
		dmp_diff,
		difference = args$meth_cutoff,
		qvalue = args$qvalue_cutoff,
		type = "all",
		chunk.size = 1e+06,
		save.db = FALSE
	)
	dmp_diff_sig_df <- as.data.frame(dmp_diff_sig)
	write.table(
		dmp_diff_sig_df,
		file.path(
			args$out_dir,
			paste(args$control, "vs", args$treated, args$level, args$context, "2k_Sig.txt", sep = "_")),
		sep = "\t",
		quote = FALSE,
		row.names = FALSE
	)

	dmp_diff_sig_hyper <- getMethylDiff(
		dmp_diff,
		difference = args$meth_cutoff,
		qvalue = args$qvalue_cutoff,
		type = "hyper",
		chunk.size = 1e+06,
		save.db = FALSE
	)
	write.table(
		dmp_diff_sig_hyper,
		file.path(
			args$out_dir,
			paste(args$control, "vs", args$treated, args$level, args$context, "2k_Sig_Hyper.txt", sep = "_")),
		sep = "\t",
		quote = FALSE,
		row.names = FALSE
	)
	
	dmp_diff_sig_hypo <- getMethylDiff(
		dmp_diff,
		difference = args$meth_cutoff,
		qvalue = args$qvalue_cutoff,
		type = "hypo",
		chunk.size = 1e+06,
		save.db = FALSE
	)
	write.table(
		dmp_diff_sig_hypo,
		file.path(
			args$out_dir,
			paste(args$control, "vs", args$treated, args$level, args$context, "2k_Sig_Hypo.txt", sep = "_")),
		sep = "\t",
		quote = FALSE,
		row.names = FALSE
	)

	if (!is.null(txdb)){
		dmp_sig_df_anno <- annotate_regions(
			dmp_diff_sig_df,
			txdb,
			promoter_up = args$promoter_up,
			promoter_down = args$promoter_down
		)
		write.table(
			dmp_sig_df_anno,
			file.path(
				args$out_dir,
				paste(args$control, "vs", args$treated, args$level, args$context, "2k_Sig_Anno.txt", sep = "_")),
			sep = "\t",
			quote = FALSE,
			row.names = FALSE
		)
	}
}

## 7. getCorrelation
pdf(file.path(
	args$out_dir, 
	paste(args$control, "vs", args$treated, args$level, args$context, "Corr.pdf", sep = "_")),
	width = args$plot_width,
	height = args$plot_height)
getCorrelation(
	unite_plot,
	method = args$corr_method,
	plot = TRUE
)
dev.off()

## 8. PCASamples
pdf(file.path(
	args$out_dir,
	paste(args$control, "vs", args$treated, args$level, args$context, "PCA.pdf", sep = "_")),
	width = args$plot_width,
	height = args$plot_height)
PCASamples(
	unite_plot,
	screeplot = FALSE,
	adj.lim = c(0.0004,0.1),
	scale = TRUE,
	center = TRUE,
	comp = c(1, 2),
	transpose = TRUE,
	sd.filter = TRUE,
	sd.threshold = 0.5,
	filterByQuantile = TRUE
)
dev.off()

pdf(file.path(
	args$out_dir,
	paste(args$control, "vs", args$treated, args$level, args$context, "PCA_Scree.pdf", sep = "_")),
	width = args$plot_width,
	height = args$plot_height)
PCASamples(
	unite_plot,
	screeplot = TRUE,
	adj.lim = c(0.0004,0.1),
	scale = TRUE,
	center = TRUE,
	comp = c(1, 2),
	transpose = TRUE,
	sd.filter = TRUE,
	sd.threshold = 0.5,
	filterByQuantile = TRUE
)
dev.off()

## 11. diffMethPerChr
pdf(file.path(
	args$out_dir, 
	paste(args$control, "vs", args$treated, args$level, args$context, "DiffPerChr.pdf", sep = "_")),
	width = args$plot_width,
	height = args$plot_height)
diffMethPerChr(
	diff_plot,
	plot = TRUE,
	qvalue.cutoff = args$qvalue_cutoff,
	meth.cutoff = args$meth_cutoff,
	keep.empty.chrom = FALSE
)
dev.off()
