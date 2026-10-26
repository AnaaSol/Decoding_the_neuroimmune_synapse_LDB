#!/usr/bin/env Rscript
################################################################################
# Phase 5 / step 5a: recover atlas effector genes that Phase 3 QC's min-cells
# filter removed from the integrated snRNA object (e.g. the lncRNA SNCA-AS1),
# by reading their counts straight from the raw 10x matrices for exactly the
# nuclei retained in GSE178146_integrated.rds. Output is a small sparse
# genes x nuclei matrix that 05_effector_genes_in_lbd_cortex.R appends.
################################################################################
suppressPackageStartupMessages({ library(Matrix); library(Seurat) })
project_dir <- "/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB"
genes <- c("SNCA-AS1", "SNCA", "MMRN1", "IDUA")
s <- readRDS(file.path(project_dir, "data/GSE178146_processed/GSE178146_integrated.rds"))
keep_cells <- colnames(s); rm(s); invisible(gc())
out <- list()
for (smp in unique(sub("_.*", "", keep_cells))) {
  d <- file.path(project_dir, "data/GSE178146_extracted", paste0(smp, "_raw_feature_bc_matrix"))
  d <- list.files(d, pattern = "matrix.mtx.gz$", recursive = TRUE, full.names = TRUE)[1]
  d <- dirname(d)
  feat <- read.delim(gzfile(file.path(d, "features.tsv.gz")), header = FALSE)
  bc <- paste0(smp, "_", readLines(gzfile(file.path(d, "barcodes.tsv.gz"))))
  rows <- match(genes, feat$V2)
  # stream the MatrixMarket triplets and keep only the target gene rows
  # (readMM on a raw droplet matrix is prohibitively slow)
  trip <- read.table(pipe(sprintf("zcat %s | awk 'NR>3 && (%s)'", shQuote(file.path(d, "matrix.mtx.gz")),
                                  paste0("$1==", rows, collapse = " || "))),
                     col.names = c("i", "j", "x"))
  cols <- which(bc %in% keep_cells)
  trip <- trip[trip$j %in% cols, ]
  x <- sparseMatrix(i = match(trip$i, rows), j = match(trip$j, cols), x = trip$x,
                    dims = c(length(genes), length(cols)), dimnames = list(genes, bc[cols]))
  out[[smp]] <- x
  cat(sprintf("%s: %d nuclei, SNCA-AS1 UMIs=%d\n", smp, length(cols), sum(x["SNCA-AS1", ])))
}
x <- do.call(cbind, out)
stopifnot(all(colnames(x) %in% keep_cells), ncol(x) == length(keep_cells))
saveRDS(x, file.path(project_dir, "data/coloc_atlas/effector_expression/raw_recovered_gene_counts.rds"))
cat("Done:", nrow(x), "genes x", ncol(x), "nuclei\n")
