#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rhdf5)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
source_file <- if (length(args) >= 1) args[[1]] else "data/e12_celltypes_spatialclustering.h5ad"
out_file <- if (length(args) >= 2) args[[2]] else "assets/data/site-data.json"

if (!file.exists(source_file)) {
  stop("Source H5AD file was not found: ", source_file)
}

read_cat <- function(path, field) {
  categories <- h5read(path, paste0("/obs/", field, "/categories"))
  codes <- as.integer(h5read(path, paste0("/obs/", field, "/codes")))
  values <- rep(NA_character_, length(codes))
  ok <- !is.na(codes) & codes >= 0
  values[ok] <- categories[codes[ok] + 1]
  list(categories = as.character(categories), codes = codes, values = values)
}

read_numeric <- function(path, key) {
  as.numeric(h5read(path, paste0("/obs/", key)))
}

orient_embedding <- function(mat, n_obs) {
  mat <- as.matrix(mat)
  if (ncol(mat) == n_obs) {
    mat <- t(mat)
  }
  if (nrow(mat) != n_obs) {
    stop("Embedding shape does not match observation count.")
  }
  mat
}

round_vec <- function(x, digits = 3) {
  as.numeric(round(x, digits))
}

clean_scalar <- function(x, digits = 3) {
  if (!is.finite(x)) return(NA_real_)
  as.numeric(round(x, digits))
}

safe_colors <- function(path, name, fallback_n) {
  dataset <- paste0("/uns/", name, "_colors")
  out <- tryCatch(as.character(h5read(path, dataset)), error = function(e) character())
  if (length(out) == fallback_n) out else rep(NA_character_, fallback_n)
}

fallback_palette <- c(
  "#2F80ED", "#F2994A", "#27AE60", "#EB5757", "#9B51E0",
  "#00A6A6", "#B7791F", "#D946EF", "#64748B", "#16A34A",
  "#E11D48", "#0891B2", "#7C3AED", "#CA8A04", "#475569"
)

category_fields <- c(
  "cell_types", "finer_cell_types", "leiden", "hybrid_leiden",
  "library_id", "timepoint", "excitatory", "inhibitory"
)

cats <- setNames(lapply(category_fields, function(field) read_cat(source_file, field)), category_fields)
n_obs <- length(cats$cell_types$codes)

spatial <- orient_embedding(h5read(source_file, "/obsm/spatial"), n_obs)
umap <- orient_embedding(h5read(source_file, "/obsm/X_umap"), n_obs)
generic <- orient_embedding(h5read(source_file, "/obsm/generic"), n_obs)
spatial3d <- orient_embedding(h5read(source_file, "/obsm/spatial3d"), n_obs)

n_counts <- read_numeric(source_file, "n_counts")
n_genes <- read_numeric(source_file, "n_genes_by_counts")
total_counts <- read_numeric(source_file, "total_counts")

cells <- vector("list", n_obs)
for (i in seq_len(n_obs)) {
  cells[[i]] <- c(
    clean_scalar(spatial[i, 1]), clean_scalar(spatial[i, 2]),
    clean_scalar(umap[i, 1]), clean_scalar(umap[i, 2]),
    clean_scalar(generic[i, 1]), clean_scalar(generic[i, 2]),
    clean_scalar(spatial3d[i, 1]), clean_scalar(spatial3d[i, 2]), clean_scalar(spatial3d[i, 3]),
    cats$cell_types$codes[i], cats$finer_cell_types$codes[i],
    cats$leiden$codes[i], cats$hybrid_leiden$codes[i],
    cats$library_id$codes[i], cats$timepoint$codes[i],
    cats$excitatory$codes[i], cats$inhibitory$codes[i],
    round(n_counts[i], 2), round(n_genes[i], 2), round(total_counts[i], 2)
  )
}

annotation_summary <- lapply(names(cats), function(field) {
  info <- cats[[field]]
  colors <- safe_colors(source_file, field, length(info$categories))
  if (all(is.na(colors))) {
    colors <- rep(fallback_palette, length.out = length(info$categories))
  }
  counts <- tabulate(info$codes[info$codes >= 0] + 1, nbins = length(info$categories))
  data.frame(
    label = info$categories,
    count = as.integer(counts),
    color = colors,
    stringsAsFactors = FALSE
  )
})
names(annotation_summary) <- names(cats)

var_names <- as.character(h5read(source_file, "/var/_index"))
gene_stats <- data.frame(
  gene = var_names,
  mean = round_vec(h5read(source_file, "/var/mean"), 5),
  mean_counts = round_vec(h5read(source_file, "/var/mean_counts"), 5),
  n_cells = as.integer(h5read(source_file, "/var/n_cells")),
  total_counts = as.integer(h5read(source_file, "/var/total_counts")),
  pct_dropout_by_counts = round_vec(h5read(source_file, "/var/pct_dropout_by_counts"), 3),
  stringsAsFactors = FALSE
)
gene_stats <- gene_stats[order(gene_stats$n_cells, decreasing = TRUE), ]

h5_tree <- h5ls(source_file, recursive = TRUE)
layers <- h5_tree$name[h5_tree$group == "/layers"]
embeddings <- h5_tree$name[h5_tree$group == "/obsm"]

data <- list(
  metadata = list(
    title = "E12 cerebellar spatial transcriptomics + scRNA-seq",
    source_file = basename(source_file),
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    n_cells = n_obs,
    n_genes = length(var_names),
    layers = as.character(layers),
    embeddings = as.character(embeddings),
    duplicate_note = "The sibling file named 'e12_celltypes_spatialclustering.h5ad 2' has the same checksum and is treated as a duplicate."
  ),
  schema = c(
    "spatial_x", "spatial_y", "umap_x", "umap_y", "generic_x", "generic_y",
    "spatial3d_x", "spatial3d_y", "spatial3d_z",
    "cell_types", "finer_cell_types", "leiden", "hybrid_leiden",
    "library_id", "timepoint", "excitatory", "inhibitory",
    "n_counts", "n_genes_by_counts", "total_counts"
  ),
  annotations = annotation_summary,
  genes = gene_stats,
  qc = list(
    n_counts = as.list(summary(n_counts)),
    n_genes_by_counts = as.list(summary(n_genes)),
    total_counts = as.list(summary(total_counts))
  ),
  cells = cells
)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
write_json(data, out_file, auto_unbox = TRUE, digits = NA, dataframe = "rows", pretty = FALSE, na = "null")
cat("Wrote", out_file, "with", n_obs, "cells and", length(var_names), "genes\n")
