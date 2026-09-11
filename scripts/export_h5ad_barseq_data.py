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
    "#7C3AED", "#CA8A04", "#475569", "#F97316", "#14B8A6", "#8B5CF6",
]

FIELD_PALETTES = {
    "integrated_cell_type": {
        "purkinje cells": "#0000ff", "midbrain-derived + int/lat dcn": "#ff0000",
        "medial dcn": "#ffa500", "dcn": "#ffa500", "astroglia": "#808080",
        "": "#ffffff", "choroid plexus": "#000000", "granule cells": "#008000",
        "midbrain-fated cells": "#a52a2a", "molecular layer interneurons": "#ffc0cb",
        "outside cb": "#800080", "glia/oligodendrocytes": "#add8e6",
        "i1": "#ffff00", "unknown": "#bbbbbb",
    },
    "CN_exc_inhib": {
        "other": "#d3d3d3", "dcn": "#800080", "i1": "#ff0000",
        "mli": "#008000", "i2/3": "#008000",
    },
}

FIELD_ORDERS = {
    "integrated_cell_type": [
        "outside cb", "choroid plexus", "astroglia", "glia/oligodendrocytes",
        "molecular layer interneurons", "granule cells", "dcn", "purkinje cells", "",
        "midbrain-derived + int/lat dcn", "medial dcn", "midbrain-fated cells", "i1", "unknown",
    ],
    "CN_exc_inhib": ["other", "dcn", "i1", "mli"],
}


def decode(values):
    return [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values]


def label_color(label, field=None):
    key = " ".join(str(label).strip().lower().replace("_", " ").split())
    if key in FIELD_PALETTES.get(field, {}):
        return FIELD_PALETTES[field][key]
    if key in {"na", "n/a", "nan", "none", "null", "not available", "missing", "unknown", "unassigned"}:
        return "#9aa39f"
    overrides = {
        "extracerebellar-fated": "#111111", "intp": "#2ca25f", "inta/lat": "#f28e2b",
        "inta": "#e76f9a", "lat": "#d62728", "medearly": "#8c564b", "med early": "#8c564b",
        "early medial": "#8c564b", "medlate": "#377eb8", "med late": "#377eb8",
        "late medial": "#377eb8", "rl": "#78cbe6", "vz": "#78cbe6",
        "int/lat prog": "#f2c94c", "int+latprog": "#f2c94c", "i1": "#f28e2b",
        "i2/3": "#2ca25f", "i2": "#2ca25f", "i3": "#2ca25f",
        "other": "#d9dedb", "others": "#d9dedb",
    }
    if key in overrides:
        return overrides[key]
    value = 0
    for byte in str(label).encode("utf-8"):
        value = (value * 33 + byte) % 2147483647
    return PALETTE[value % len(PALETTE)].lower()


def reorder_categories(field, levels, codes):
    desired = FIELD_ORDERS.get(field)
    if not desired:
        return levels, codes
    rank = {label: index for index, label in enumerate(desired)}
    old_order = sorted(
        range(len(levels)),
        key=lambda index: (
            rank.get(" ".join(levels[index].strip().lower().replace("_", " ").split()), len(desired)),
            index,
        ),
    )
    remap = np.empty(len(levels), dtype=np.int32)
    for new_index, old_index in enumerate(old_order):
        remap[old_index] = new_index
    reordered_codes = np.asarray([remap[code] if code >= 0 else code for code in codes], dtype=np.int32)
    return [levels[index] for index in old_order], reordered_codes


def metadata_column(node):
    if isinstance(node, h5py.Group) and "categories" in node and "codes" in node:
        levels = decode(node["categories"][:])
        codes = node["codes"][:].astype(np.int32)
        if len(levels) > 250:
            try:
                numeric_levels = np.asarray(levels, dtype=float)
                return "numeric", None, np.asarray([numeric_levels[code] if code >= 0 else np.nan for code in codes])
            except ValueError:
                pass
        return "annotation", levels, codes
    values = node[:]
    if values.dtype.kind in "SUO":
        text = decode(values)
        levels = list(dict.fromkeys(text))
        lookup = {value: index for index, value in enumerate(levels)}
        return "annotation", levels, np.asarray([lookup[value] for value in text], dtype=np.int32)
    return "numeric", None, values.astype(float)


def numeric(handle, field, count):
    path = f"obs/{field}"
    return handle[path][:].astype(float) if path in handle else np.full(count, np.nan)


def rounded(value, digits=3):
    value = float(value)
    return round(value, digits) if np.isfinite(value) else None


def safe_int(value):
    return int(value) if np.isfinite(value) else 0


def read_var_names(handle):
    for field in ("_index", "Names"):
        if field in handle["var"]:
            return decode(handle[f"var/{field}"][:])
    raise KeyError("No gene-name field was found in var")


def export(source, output, title):
    with h5py.File(source, "r") as handle:
        spatial = handle["obsm/spatial"][:]
        count = spatial.shape[0]
        umap = handle["obsm/X_umap"][:] if "obsm/X_umap" in handle else spatial
        generic = handle["obsm/generic"][:] if "obsm/generic" in handle else spatial
        spatial3d = handle["obsm/spatial3d"][:] if "obsm/spatial3d" in handle else np.column_stack((spatial, np.zeros(count)))
        categories, codes, numeric_metadata = {}, {}, {}
        for field in handle["obs"].keys():
            if field == "_index":
                continue
            kind, levels, values = metadata_column(handle[f"obs/{field}"])
            if kind == "numeric":
                numeric_metadata[field] = values
                continue
            levels, values = reorder_categories(field, levels, values)
            codes[field] = values
            counts = np.bincount(values[values >= 0], minlength=len(levels)) if levels else []
            categories[field] = [
                {"label": label, "count": int(counts[index]), "color": label_color(label, field)}
                for index, label in enumerate(levels)
            ]

        n_counts = numeric(handle, "n_counts", count)
        n_genes = numeric(handle, "n_genes_by_counts", count)
        total_counts = numeric(handle, "total_counts", count)
        var_names = read_var_names(handle)
        genes = []
        for index, gene in enumerate(var_names):
            mean = handle["var/mean"][index] if "var/mean" in handle else 0
            mean_counts = handle["var/mean_counts"][index] if "var/mean_counts" in handle else mean
            n_cells = handle["var/n_cells"][index] if "var/n_cells" in handle else 0
            total = handle["var/total_counts"][index] if "var/total_counts" in handle else 0
            dropout = handle["var/pct_dropout_by_counts"][index] if "var/pct_dropout_by_counts" in handle else 0
            genes.append({
                "gene": gene,
                "mean": rounded(mean, 5),
                "mean_counts": rounded(mean_counts, 5),
                "n_cells": safe_int(n_cells),
                "total_counts": safe_int(total),
                "pct_dropout_by_counts": rounded(dropout, 3),
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
            "schema": ["spatial_x", "spatial_y", *codes.keys(), *numeric_metadata.keys()],
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
                *[int(values[index]) for values in codes.values()],
                *[rounded(values[index], 4) for values in numeric_metadata.values()],
            ])

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")
    print(f"Wrote {output} with {count} cells and {len(var_names)} genes")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise SystemExit("Usage: export_h5ad_barseq_data.py SOURCE.h5ad OUTPUT.json [TITLE]")
    export(Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3] if len(sys.argv) > 3 else Path(sys.argv[1]).stem)
