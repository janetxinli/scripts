#!/usr/bin/env Rscript

library(tidyverse)
library(gprofiler2)
library(argparser)

# Parse arguments
p <- arg_parser("Identify enriched GO terms in orthogroup genes")
p <- add_argument(
    p,
    "genes",
    help="Gene list file (one gene per line)"
)
p <- add_argument(
    p,
    "--output",
    help="Output file name [orthogroup_genes.gsea.txt]",
    default="orthogroup_genes.gsea.txt"
)
p <- add_argument(
    p,
    "--cutoff",
    help="Cutoff for false discovery rate [0.01]",
    default=0.01
)
p <- add_argument(
    p,
    "--gmt",
    help="Create a custom GMT file",
    default=NA
)
p <- add_argument(
    p,
    "--organism",
    help="Organism for GSEA",
    default=NA
)
p <- add_argument(
   p, 
   "--all",
   flag=TRUE,
   help="Print results for all terms (sets significant=FALSE in gost call)"
)

argv <- parse_args(p)
stopifnot(xor(is.na(argv$organism), is.na(argv$gmt)))

# Load in custom GMT file
if (!is.na(argv$gmt)) {
    org <- upload_GMT_file("gsea/orthogroup_genes.gmt")
} else {
    org <- argv$organism
}

# Load in genes
gene_list <- read_csv(argv$genes, col_names=F)
gene_list <- gene_list$X1

# Run GOSt
res <- gost(
    gene_list,
    organism=org,
    correction_method="fdr",
    user_threshold=argv$cutoff,
    significant=!argv$all
)

res <- res$result %>% select(-parents)

# Save results
write_delim(res, file=argv$output, delim="\t")
