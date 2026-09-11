#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
source_dir <- if (length(args) >= 1) args[[1]] else "data/scrna-source"
out_dir <- if (length(args) >= 2) args[[2]] else "assets/data/scrna-counts"
dataset_filter <- if (length(args) >= 3) args[[3]] else NULL
genes_per_shard <- if (length(args) >= 4) as.integer(args[[4]]) else 64L

datasets <- list(
  list(file = "URL_082026.rds", slug = "url"),
  list(file = "RLprog.rds", slug = "rl-prog"),
  list(file = "chicken_excitatory.rds", slug = "chicken-excitatory"),
  list(file = "chicken_inhib.rds", slug = "chicken-inhibitory"),
  list(file = "all_inhib.rds", slug = "all-inhib"),
  list(file = "combined_alltp08_2026.rds", slug = "combined-alltp08-2026")
)
if (!is.null(dataset_filter)) datasets <- Filter(function(x) x$slug == dataset_filter, datasets)
if (!length(datasets)) stop("No matching dataset")

json_string <- function(x) encodeString(as.character(x), quote = '"')

add_assay5_names <- function(matrix, assay_slots) {
  if (is.null(rownames(matrix)) && !is.null(assay_slots[["features"]])) {
    names <- attr(assay_slots[["features"]], "dimnames")[[1]]
    if (length(names) == nrow(matrix)) rownames(matrix) <- names
  }
  if (is.null(colnames(matrix)) && !is.null(assay_slots[["cells"]])) {
    names <- attr(assay_slots[["cells"]], "dimnames")[[1]]
    if (length(names) == ncol(matrix)) colnames(matrix) <- names
  }
  matrix
}

plain_logical_map <- function(map) {
  values <- map
  attr(values, "class") <- NULL
  values
}

join_assay5_count_layers <- function(layers, assay_slots) {
  count_layer_names <- names(layers)[startsWith(names(layers), "counts")]
  if (!length(count_layer_names)) return(NULL)

  feature_map <- assay_slots[["features"]]
  cell_map <- assay_slots[["cells"]]
  if (is.null(feature_map) || is.null(cell_map)) return(NULL)
  feature_values <- plain_logical_map(feature_map)
  cell_values <- plain_logical_map(cell_map)
  feature_names <- attr(feature_map, "dimnames")[[1]]
  cell_names <- attr(cell_map, "dimnames")[[1]]
  map_layers <- attr(feature_map, "dimnames")[[2]]

  matrices <- lapply(count_layer_names, function(layer_name) {
    map_index <- match(layer_name, map_layers)
    if (is.na(map_index)) stop("No feature/cell map found for layer ", layer_name)
    matrix <- layers[[layer_name]]
    row_names <- feature_names[feature_values[, map_index]]
    column_names <- cell_names[cell_values[, map_index]]
    if (nrow(matrix) != length(row_names) || ncol(matrix) != length(column_names)) {
      stop("Layer dimensions do not match feature/cell maps for ", layer_name)
    }
    dimnames(matrix) <- list(row_names, column_names)
    matrix
  })

  reference_genes <- rownames(matrices[[1]])
  if (any(vapply(matrices, function(matrix) !identical(rownames(matrix), reference_genes), logical(1)))) {
    stop("Split count layers do not share the same ordered gene set")
  }
  do.call(cbind, matrices)
}

count_matrix <- function(slots, cell_count) {
  assays <- slots[["assays"]]
  for (assay_name in c("RNA", "SCT")) {
    if (!assay_name %in% names(assays)) next
    assay_slots <- attributes(assays[[assay_name]])
    candidate <- assay_slots[["counts"]]
    if (!is.null(candidate) && ncol(candidate) == cell_count) {
      return(add_assay5_names(candidate, assay_slots))
    }
    layers <- assay_slots[["layers"]]
    if (!is.null(layers) && "counts" %in% names(layers) && ncol(layers[["counts"]]) == cell_count) {
      return(add_assay5_names(layers[["counts"]], assay_slots))
    }
    if (!is.null(layers)) {
      joined_counts <- join_assay5_count_layers(layers, assay_slots)
      if (!is.null(joined_counts) && ncol(joined_counts) == cell_count) return(joined_counts)
    }
  }
  stop("No complete sparse count matrix was found")
}

write_shard <- function(matrix_t, start_gene, end_gene, path) {
  connection <- gzfile(path, open = "wb", compression = 9)
  on.exit(close(connection), add = TRUE)
  writeBin(as.integer(end_gene - start_gene + 1L), connection, size = 4, endian = "little")
  pointers <- attr(matrix_t, "p")
  indices <- attr(matrix_t, "i")
  values <- attr(matrix_t, "x")
  for (gene_index in start_gene:end_gene) {
    from <- pointers[[gene_index]] + 1L
    to <- pointers[[gene_index + 1L]]
    nnz <- max(0L, to - from + 1L)
    writeBin(as.integer(nnz), connection, size = 4, endian = "little")
    if (nnz > 0L) {
      positions <- from:to
      writeBin(as.integer(indices[positions]), connection, size = 4, endian = "little")
      writeBin(as.numeric(values[positions]), connection, size = 4, endian = "little")
    }
  }
}

export_dataset <- function(spec) {
  source_file <- file.path(source_dir, spec$file)
  message("Reading ", source_file)
  object <- readRDS(source_file)
  slots <- attributes(object)
  metadata <- slots[["meta.data"]]
  counts <- count_matrix(slots, nrow(metadata))
  if (!inherits(counts, "sparseMatrix")) counts <- as(counts, "sparseMatrix")

  matrix_cells <- colnames(counts)
  metadata_cells <- rownames(metadata)
  if (!is.null(matrix_cells)) {
    order_index <- match(metadata_cells, matrix_cells)
    if (anyNA(order_index)) stop("Count matrix cells do not align with metadata")
    counts <- counts[, order_index, drop = FALSE]
  }

  matrix_t <- as(t(counts), "dgCMatrix")
  genes <- colnames(matrix_t)
  if (is.null(genes)) genes <- rownames(counts)
  dataset_dir <- file.path(out_dir, spec$slug)
  dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)
  shard_count <- ceiling(length(genes) / genes_per_shard)
  shards <- vector("list", shard_count)
  for (shard_index in seq_len(shard_count)) {
    start_gene <- (shard_index - 1L) * genes_per_shard + 1L
    end_gene <- min(length(genes), shard_index * genes_per_shard)
    filename <- sprintf("counts-%04d.bin.gz", shard_index - 1L)
    message("Writing ", spec$slug, " shard ", shard_index, "/", shard_count)
    write_shard(matrix_t, start_gene, end_gene, file.path(dataset_dir, filename))
    shards[[shard_index]] <- list(file = filename, start = start_gene - 1L, count = end_gene - start_gene + 1L)
  }

  genes_json <- paste(json_string(genes), collapse = ",")
  shards_json <- paste(vapply(shards, function(shard) {
    paste0('{"file":', json_string(shard$file), ',"start":', shard$start, ',"count":', shard$count, "}")
  }, character(1)), collapse = ",")
  index_json <- paste0(
    '{"dataset":', json_string(spec$slug),
    ',"source_file":', json_string(spec$file),
    ',"encoding":"gzip-dense-records-v1","value_type":"float32","index_type":"uint32"',
    ',"n_cells":', nrow(metadata), ',"n_genes":', length(genes),
    ',"genes":[', genes_json, '],"shards":[', shards_json, "]}"
  )
  writeLines(index_json, file.path(dataset_dir, "index.json"), useBytes = TRUE)
  total_size <- sum(file.info(list.files(dataset_dir, full.names = TRUE))$size)
  message("Wrote ", spec$slug, ": ", length(genes), " genes, ", format(total_size, big.mark = ","), " bytes")
  rm(object, slots, metadata, counts, matrix_t)
  invisible(gc())
}

invisible(lapply(datasets, export_dataset))
