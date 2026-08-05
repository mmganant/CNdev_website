#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
source_dir <- if (length(args) >= 1) args[[1]] else "data/scrna-source"
out_dir <- if (length(args) >= 2) args[[2]] else "assets/data/scrna"
dataset_filter <- if (length(args) >= 3) args[[3]] else NULL

datasets <- list(
  list(
    file = "URL.rds",
    slug = "url",
    title = "URL integrated atlas",
    fields = c("fine_cell_type", "cell_type", "cell_type_vague", "TimePoint", "Sample", "branch", "DCN_status")
  ),
  list(
    file = "all_inhib.rds",
    slug = "all-inhib",
    title = "Integrated inhibitory atlas",
    fields = c("cellType", "precisest_label", "broad_lineage", "dev_state", "species", "stage", "timepoint", "dataset", "Age", "cluster")
  ),
  list(
    file = "combined_alltp08_2026.rds",
    slug = "combined-alltp08-2026",
    title = "Combined developmental atlas",
    fields = c("final.clusters2", "classes", "cell_type", "fine_cell_type", "TimePoint", "med_2p_branch_final", "RL_subcluster_named")
  )
)
if (!is.null(dataset_filter)) {
  datasets <- Filter(function(spec) identical(spec$slug, dataset_filter), datasets)
  if (!length(datasets)) stop("Unknown dataset slug: ", dataset_filter)
}

palette <- c(
  "#2f80ed", "#f2994a", "#27ae60", "#eb5757", "#9b51e0", "#00a6a6",
  "#b7791f", "#d946ef", "#64748b", "#16a34a", "#e11d48", "#0891b2",
  "#7c3aed", "#ca8a04", "#475569", "#f97316", "#14b8a6", "#8b5cf6"
)

gene_panel <- c(
  "Atoh1", "Ptf1a", "Pax6", "Pax2", "Sox2", "Mki67", "Top2a", "Neurod1",
  "Neurod2", "Barhl1", "Eomes", "Meis2", "Lhx9", "Tbr1", "Foxp2", "Tle4",
  "Pou4f1", "Slc17a6", "Slc17a7", "Gad1", "Gad2", "Slc32a1", "Reln", "Calb1",
  "Calb2", "Pcp2", "Gabra6", "Aldoc", "Gfap", "Aqp4", "Olig2", "Pdgfra"
)

clean_number <- function(x, digits = 3) {
  out <- round(as.numeric(x), digits)
  out[!is.finite(out)] <- NA_real_
  out
}

json_string <- function(x) {
  ifelse(is.na(x), "null", encodeString(as.character(x), quote = '"', na.encode = FALSE))
}

json_number <- function(x) {
  ifelse(is.na(x) | !is.finite(x), "null", format(x, trim = TRUE, scientific = FALSE, digits = 15))
}

write_dataset_json <- function(path, result) {
  connection <- file(path, open = "wb")
  on.exit(close(connection), add = TRUE)
  writeLines("{", connection, useBytes = TRUE)
  metadata <- result$metadata
  metadata_json <- paste0(
    '"title":', json_string(metadata$title),
    ',"source_file":', json_string(metadata$source_file),
    ',"generated_at":', json_string(metadata$generated_at),
    ',"n_cells":', metadata$n_cells,
    ',"assays":[', paste(json_string(metadata$assays), collapse = ","), "]"
  )
  writeLines(paste0('"metadata":{', metadata_json, "},"), connection, useBytes = TRUE)
  writeLines(paste0('"schema":[', paste(json_string(result$schema), collapse = ","), "],"), connection, useBytes = TRUE)
  writeLines('"annotations":{', connection, useBytes = TRUE)
  annotation_names <- names(result$annotations)
  for (field_index in seq_along(annotation_names)) {
    field <- annotation_names[[field_index]]
    rows <- result$annotations[[field]]
    row_json <- vapply(seq_len(nrow(rows)), function(i) {
      paste0('{"label":', json_string(rows$label[[i]]), ',"count":', rows$count[[i]], ',"color":', json_string(rows$color[[i]]), "}")
    }, character(1))
    suffix <- if (field_index < length(annotation_names)) "," else ""
    writeLines(paste0(json_string(field), ":[", paste(row_json, collapse = ","), "]", suffix), connection, useBytes = TRUE)
  }
  writeLines('},"genes":[', connection, useBytes = TRUE)
  for (gene_index in seq_along(result$genes)) {
    gene <- result$genes[[gene_index]]
    suffix <- if (gene_index < length(result$genes)) "," else ""
    line <- paste0(
      '{"gene":', json_string(gene$gene),
      ',"max":', json_number(gene$max),
      ',"values":[', paste(json_number(gene$values), collapse = ","), "]}", suffix
    )
    writeLines(line, connection, useBytes = TRUE)
  }
  writeLines('],"cells":[', connection, useBytes = TRUE)
  for (i in seq_along(result$cells)) {
    suffix <- if (i < length(result$cells)) "," else ""
    writeLines(paste0("[", paste(json_number(result$cells[[i]]), collapse = ","), "]", suffix), connection, useBytes = TRUE)
  }
  writeLines("]}", connection, useBytes = TRUE)
}

expression_matrix <- function(slots) {
  assays <- slots[["assays"]]
  preferred <- c("SCT", "RNA")
  for (assay_name in preferred[preferred %in% names(assays)]) {
    assay_slots <- attributes(assays[[assay_name]])
    add_assay5_names <- function(matrix) {
      if (is.null(rownames(matrix)) && !is.null(assay_slots[["features"]])) {
        feature_names <- attr(assay_slots[["features"]], "dimnames")[[1]]
        if (length(feature_names) == nrow(matrix)) rownames(matrix) <- feature_names
      }
      if (is.null(colnames(matrix)) && !is.null(assay_slots[["cells"]])) {
        cell_names <- attr(assay_slots[["cells"]], "dimnames")[[1]]
        if (length(cell_names) == ncol(matrix)) colnames(matrix) <- cell_names
      }
      matrix
    }
    if (!is.null(assay_slots[["data"]])) return(add_assay5_names(assay_slots[["data"]]))
    layers <- assay_slots[["layers"]]
    if (!is.null(layers)) {
      data_layers <- names(layers)[startsWith(names(layers), "data")]
      if ("data" %in% data_layers) return(add_assay5_names(layers[["data"]]))
      if (length(data_layers) == 1) return(add_assay5_names(layers[[data_layers[[1]]]]))
    }
  }
  NULL
}

extract_dataset <- function(spec) {
  source_file <- file.path(source_dir, spec$file)
  if (!file.exists(source_file)) stop("Source RDS file was not found: ", source_file)

  message("Reading ", source_file)
  object <- readRDS(source_file)
  slots <- attributes(object)
  metadata <- slots[["meta.data"]]
  reductions <- slots[["reductions"]]
  if (is.null(metadata) || is.null(reductions[["umap"]])) {
    stop(spec$file, " does not contain Seurat metadata and a UMAP reduction")
  }

  embedding <- attributes(reductions[["umap"]])[["cell.embeddings"]]
  if (is.null(embedding) || ncol(embedding) < 2 || nrow(embedding) != nrow(metadata)) {
    stop("UMAP dimensions do not match metadata in ", spec$file)
  }

  fields <- spec$fields[spec$fields %in% names(metadata)]
  annotations <- list()
  codes <- list()
  for (field in fields) {
    values <- as.character(metadata[[field]])
    values[is.na(values) | !nzchar(values)] <- "Not available"
    levels <- unique(values)
    field_codes <- match(values, levels) - 1L
    counts <- tabulate(field_codes + 1L, nbins = length(levels))
    annotations[[field]] <- data.frame(
      label = levels,
      count = as.integer(counts),
      color = rep(palette, length.out = length(levels)),
      stringsAsFactors = FALSE
    )
    codes[[field]] <- field_codes
  }

  qc_fields <- c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent.mt")
  qc_fields <- qc_fields[qc_fields %in% names(metadata)]
  schema <- c("umap_x", "umap_y", fields, qc_fields)
  columns <- c(
    list(clean_number(embedding[, 1]), clean_number(embedding[, 2])),
    unname(codes),
    lapply(qc_fields, function(field) clean_number(metadata[[field]], 2))
  )
  cells <- lapply(seq_len(nrow(metadata)), function(i) {
    unname(vapply(columns, function(column) column[[i]], numeric(1)))
  })

  expression <- expression_matrix(slots)
  genes <- list()
  if (!is.null(expression)) {
    expression_names <- rownames(expression)
    present_genes <- gene_panel[gene_panel %in% expression_names]
    expression_cells <- colnames(expression)
    metadata_cells <- rownames(metadata)
    order_index <- if (is.null(expression_cells) && ncol(expression) == nrow(metadata)) {
      seq_len(nrow(metadata))
    } else {
      match(metadata_cells, expression_cells)
    }
    if (anyNA(order_index)) stop("Expression matrix cells do not align with metadata in ", spec$file)
    genes <- lapply(present_genes, function(gene) {
      values <- clean_number(as.numeric(expression[gene, order_index]), 3)
      list(
        gene = gene,
        max = clean_number(as.numeric(quantile(values, 0.99, na.rm = TRUE, names = FALSE)), 3),
        values = values
      )
    })
  }

  assays <- names(slots[["assays"]])
  result <- list(
    metadata = list(
      title = spec$title,
      source_file = spec$file,
      generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      n_cells = nrow(metadata),
      assays = assays
    ),
    schema = schema,
    annotations = annotations,
    genes = genes,
    cells = cells
  )

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(out_dir, paste0(spec$slug, ".json"))
  write_dataset_json(output_file, result)
  message("Wrote ", output_file, " with ", nrow(metadata), " cells")
  rm(object, slots, metadata, reductions, embedding, result, cells, columns, expression, genes)
  invisible(gc())

  list(
    id = spec$slug,
    title = spec$title,
    source_file = spec$file,
    data_url = paste0("assets/data/scrna/", spec$slug, ".json")
  )
}

manifest <- lapply(datasets, extract_dataset)
if (is.null(dataset_filter)) {
  manifest_rows <- vapply(manifest, function(item) {
    paste0(
      '{"id":', json_string(item$id),
      ',"title":', json_string(item$title),
      ',"source_file":', json_string(item$source_file),
      ',"data_url":', json_string(item$data_url), "}"
    )
  }, character(1))
  writeLines(paste0('{"datasets":[', paste(manifest_rows, collapse = ","), "]}"), file.path(out_dir, "manifest.json"), useBytes = TRUE)
  message("Wrote scRNA-seq manifest")
}
