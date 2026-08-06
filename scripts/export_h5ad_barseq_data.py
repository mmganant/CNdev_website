#!/usr/bin/env python3
"""Export a BARseq H5AD file to the compact JSON schema used by the website."""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np

PALETTE = [
    "#2F80ED", "#F2994A", "#27AE60", "#EB5757", "#9B51E0", "#00A6A6",
    "#B7791F", "#D946EF", "#64748B", "#16A34A", "#E11D48", "#0891B2",
    "#7C3AED", "#CA8A04", "#475569",
]
CATEGORY_FIELDS = [
    "cell_types", "finer_cell_types", "leiden", "hybrid_leiden",
    "library_id", "timepoint", "excitatory", "inhibitory",
]


def decode(values):
    return [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values]


def categorical(handle, field, count):
    path = f"obs/{field}"
    if path not in handle:
        return [], np.full(count, -1, dtype=np.int32)
    node = handle[path]
    if isinstance(node, h5py.Group) and "categories" in node and "codes" in node:
        return decode(node["categories"][:]), node["codes"][:].astype(np.int32)
    values = decode(node[:])
    levels = list(dict.fromkeys(values))
    lookup = {value: index for index, value in enumerate(levels)}
    return levels, np.asarray([lookup[value] for value in values], dtype=np.int32)


def numeric(handle, field, count):
    path = f"obs/{field}"
    return handle[path][:].astype(float) if path in handle else np.full(count, np.nan)


def rounded(value, digits=3):
    value = float(value)
    return round(value, digits) if np.isfinite(value) else None


def export(source, output, title):
    with h5py.File(source, "r") as handle:
        umap = handle["obsm/X_umap"][:]
        count = umap.shape[0]
        spatial = handle["obsm/spatial"][:]
        generic = handle["obsm/generic"][:] if "obsm/generic" in handle else spatial
        spatial3d = handle["obsm/spatial3d"][:] if "obsm/spatial3d" in handle else np.column_stack((spatial, np.zeros(count)))
        categories = {}
        codes = {}
        for field in CATEGORY_FIELDS:
            levels, field_codes = categorical(handle, field, count)
            codes[field] = field_codes
            colors_path = f"uns/{field}_colors"
            colors = decode(handle[colors_path][:]) if colors_path in handle else []
            if len(colors) != len(levels):
                colors = [PALETTE[index % len(PALETTE)] for index in range(len(levels))]
            counts = np.bincount(field_codes[field_codes >= 0], minlength=len(levels)) if levels else []
            categories[field] = [
                {"label": label, "count": int(counts[index]), "color": colors[index]}
                for index, label in enumerate(levels)
            ]

        n_counts = numeric(handle, "n_counts", count)
        n_genes = numeric(handle, "n_genes_by_counts", count)
        total_counts = numeric(handle, "total_counts", count)
        var_names = decode(handle["var/_index"][:])
        genes = []
        for index, gene in enumerate(var_names):
            genes.append({
                "gene": gene,
                "mean": rounded(handle["var/mean"][index], 5),
                "mean_counts": rounded(handle["var/mean_counts"][index], 5),
                "n_cells": int(handle["var/n_cells"][index]),
                "total_counts": int(handle["var/total_counts"][index]),
                "pct_dropout_by_counts": rounded(handle["var/pct_dropout_by_counts"][index], 3),
            })
        genes.sort(key=lambda item: item["n_cells"], reverse=True)

        payload = {
            "metadata": {
                "title": title,
                "source_file": source.name,
                "generated_at": datetime.now(timezone.utc).isoformat(),
                "n_cells": count,
                "n_genes": len(var_names),
                "layers": list(handle.get("layers", {}).keys()),
                "embeddings": list(handle["obsm"].keys()),
            },
            "schema": [
                "spatial_x", "spatial_y", "umap_x", "umap_y", "generic_x", "generic_y",
                "spatial3d_x", "spatial3d_y", "spatial3d_z", *CATEGORY_FIELDS,
                "n_counts", "n_genes_by_counts", "total_counts",
            ],
            "annotations": categories,
            "genes": genes,
            "qc": {
                "n_counts": {"Min.": rounded(np.nanmin(n_counts)), "Median": rounded(np.nanmedian(n_counts)), "Mean": rounded(np.nanmean(n_counts)), "Max.": rounded(np.nanmax(n_counts))},
                "n_genes_by_counts": {"Min.": rounded(np.nanmin(n_genes)), "Median": rounded(np.nanmedian(n_genes)), "Mean": rounded(np.nanmean(n_genes)), "Max.": rounded(np.nanmax(n_genes))},
                "total_counts": {"Min.": rounded(np.nanmin(total_counts)), "Median": rounded(np.nanmedian(total_counts)), "Mean": rounded(np.nanmean(total_counts)), "Max.": rounded(np.nanmax(total_counts))},
            },
            "cells": [],
        }
        for index in range(count):
            payload["cells"].append([
                rounded(spatial[index, 0]), rounded(spatial[index, 1]),
                rounded(umap[index, 0]), rounded(umap[index, 1]),
                rounded(generic[index, 0]), rounded(generic[index, 1]),
                rounded(spatial3d[index, 0]), rounded(spatial3d[index, 1]), rounded(spatial3d[index, 2]),
                *[int(codes[field][index]) for field in CATEGORY_FIELDS],
                rounded(n_counts[index], 2), rounded(n_genes[index], 2), rounded(total_counts[index], 2),
            ])

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")
    print(f"Wrote {output} with {count} cells and {len(var_names)} genes")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise SystemExit("Usage: export_h5ad_barseq_data.py SOURCE.h5ad OUTPUT.json [TITLE]")
    export(Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3] if len(sys.argv) > 3 else Path(sys.argv[1]).stem)
