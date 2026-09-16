#!/usr/bin/env python3
"""Apply curated BARseq label harmonization to browser JSON exports."""

import json
import re
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "assets" / "data"
MANIFEST_PATH = DATA_DIR / "barseq-manifest.json"

PRECISE_COLORS = {
    "Outside Cb": "#800080",
    "Choroid Plexus": "#000000",
    "Astroglia": "#808080",
    "Glia/Oligodendrocytes": "#add8e6",
    "Molecular Layer Interneurons": "#ffc0cb",
    "Granule cells": "#008000",
    "Purkinje Cells": "#0000ff",
    "excCN": "#ff0000",
    "i1": "#ffff00",
    "Interneurons": "#ffb6c1",
    "Unknown": "#bbbbbb",
}

INTERNEURON_SUBTYPE_COLORS = [
    "#ffb6c1",  # light pink
    "#1f77b4",  # blue
    "#ff7f0e",  # orange
    "#9467bd",  # purple
    "#e377c2",  # pink
    "#17becf",  # cyan
    "#8c564b",  # brown
    "#bcbd22",  # olive
]

P0_PRECISE_COLORS = {
    "Outside Cb": "#800080",
    "Choroid Plexus": "#000000",
    "Astroglia": "#808080",
    "Glia/Oligodendrocytes": "#9ecae1",
    "Oligodendrocytes": "#3182bd",
    "Cb_prog1": "#17becf",
    "Interneurons": "#ffb6c1",
    "External Granule Layer": "#98df8a",
    "Inner Granule Layer": "#2ca02c",
    "Cb_prog2": "#006d2c",
    "Cb": "#00441b",
    "Purkinje Cells": "#0000ff",
    "excCN": "#ff0000",
    "Midbrain": "#a52a2a",
    "i1 Neurons": "#ffff00",
    "Unknown": "#bbbbbb",
}

P4_PRECISE_COLORS = {
    "Interneurons": "#ffb6c1",
    "Glia": "#9ecae1",
    "Glia1": "#6baed6",
    "Oligodendrocytes": "#3182bd",
}

CN_COLOR_OVERRIDES = {"Interneurons": "#008000"}

E17_PRECISE_ORDER = list(PRECISE_COLORS)
P0_PRECISE_ORDER = [
    "Outside Cb",
    "Choroid Plexus",
    "Astroglia",
    "Glia/Oligodendrocytes",
    "Oligodendrocytes",
    "Cb_prog1",
    "Interneurons",
    "External Granule Layer",
    "Inner Granule Layer",
    "Cb_prog2",
    "Cb",
    "Purkinje Cells",
    "excCN",
    "Midbrain",
    "i1 Neurons",
    "Unknown",
]
E17_INTEGRATED_TO_PRECISE = {
    "Outside Cb": "Outside Cb",
    "Choroid Plexus": "Choroid Plexus",
    "Astroglia": "Astroglia",
    "Glia/Oligodendrocytes": "Glia/Oligodendrocytes",
    "Molecular Layer Interneurons": "Molecular Layer Interneurons",
    "Granule cells": "Granule cells",
    "Purkinje Cells": "Purkinje Cells",
    "Midbrain-derived + Int/Lat DCN": "excCN",
    "Medial DCN": "excCN",
    "Midbrain-fated cells": "excCN",
    "i1": "i1",
    "i2/3": "Interneurons",
    "Unknown": "Unknown",
}


def annotation_labels(payload, field):
    return [row["label"] for row in payload["annotations"][field]]


def labels_for_cells(payload, field):
    column = payload["schema"].index(field)
    labels = annotation_labels(payload, field)
    return [labels[cell[column]] if cell[column] >= 0 else "Unknown" for cell in payload["cells"]]


def integrated_palette_color(label):
    """Map precise-label aliases onto the shared integrated-atlas palette."""
    key = re.sub(r"\s+", " ", str(label).strip().lower().replace("_", " "))
    if "outside" in key:
        return PRECISE_COLORS["Outside Cb"]
    if key == "cp" or "choroid" in key:
        return PRECISE_COLORS["Choroid Plexus"]
    if "astro" in key:
        return PRECISE_COLORS["Astroglia"]
    if "molecular layer" in key:
        return PRECISE_COLORS["Molecular Layer Interneurons"]
    if "glia" in key or "oligodendro" in key:
        return PRECISE_COLORS["Glia/Oligodendrocytes"]
    if "granule" in key or key == "rl" or key == "cb" or key.startswith("cb prog"):
        return PRECISE_COLORS["Granule cells"]
    if "purkinje" in key:
        return PRECISE_COLORS["Purkinje Cells"]
    if key == "i1" or key.startswith("i1 "):
        return PRECISE_COLORS["i1"]
    if any(term in key for term in ("int/lat", "int+lat", "inta", "intp", "lat prog")) or key in {"lat", "lat vz", "int vz"}:
        return "#ff0000"
    if "midbrain" in key or key == "isthmus":
        return "#a52a2a"
    if key == "exccn" or re.search(r"(^| )(dcn|cn)\d*( |$)", key) or key.startswith("med"):
        return PRECISE_COLORS["excCN"]
    if "interneuron" in key or "inhib" in key:
        return PRECISE_COLORS["Interneurons"]
    return PRECISE_COLORS["Unknown"]


def is_interneuron_subtype(label):
    key = re.sub(r"\s+", " ", str(label).strip().lower().replace("_", " "))
    return bool(
        re.fullmatch(r"interneurons?\s*\d*", key)
        or re.fullmatch(r"inhib(?:\s*prog)?\s*\d*", key)
    )


def recode(payload, field, labels, preferred_order=None, color_overrides=None):
    column = payload["schema"].index(field)
    previous = {row["label"]: row.get("color", "#bbbbbb") for row in payload["annotations"][field]}
    counts = Counter(labels)
    order = []
    if preferred_order:
        order.extend(label for label in preferred_order if counts[label])
    order.extend(label for label in dict.fromkeys(labels) if label not in order)
    lookup = {label: code for code, label in enumerate(order)}
    for cell, label in zip(payload["cells"], labels):
        cell[column] = lookup[label]
    color_overrides = color_overrides or PRECISE_COLORS
    rows = [
        {
            "label": label,
            "count": counts[label],
            "color": (
                color_overrides.get(label, integrated_palette_color(label))
                if field == "finer_cell_types"
                else color_overrides.get(label, previous.get(label, "#bbbbbb"))
            ),
        }
        for label in order
    ]
    if field == "finer_cell_types":
        interneuron_rows = [row for row in rows if is_interneuron_subtype(row["label"])]
        if len(interneuron_rows) > 1:
            for index, row in enumerate(interneuron_rows):
                row["color"] = INTERNEURON_SUBTYPE_COLORS[index % len(INTERNEURON_SUBTYPE_COLORS)]
    payload["annotations"][field] = rows


def harmonize_dataset(dataset_id, payload):
    cn_labels = labels_for_cells(payload, "CN_exc_inhib")
    cn_labels = ["Interneurons" if label == "i2/3" else label for label in cn_labels]
    recode(
        payload,
        "CN_exc_inhib",
        cn_labels,
        ["Other", "DCN", "i1", "Interneurons"],
        CN_COLOR_OVERRIDES,
    )

    if dataset_id == "E17":
        integrated = labels_for_cells(payload, "integrated_cell_type")
        missing = sorted(set(integrated) - set(E17_INTEGRATED_TO_PRECISE))
        if missing:
            raise ValueError(f"Unmapped E17 integrated labels: {missing}")
        precise = [E17_INTEGRATED_TO_PRECISE[label] for label in integrated]
        recode(payload, "finer_cell_types", precise, E17_PRECISE_ORDER)
        return

    precise = labels_for_cells(payload, "finer_cell_types")
    precise = [
        "Interneurons" if cn_label == "Interneurons" else label
        for label, cn_label in zip(precise, cn_labels)
    ]
    if "Cb" in precise:
        integrated = labels_for_cells(payload, "integrated_cell_type")
        precise = [
            E17_INTEGRATED_TO_PRECISE[integrated_label] if label == "Cb" else label
            for label, integrated_label in zip(precise, integrated)
        ]
    if dataset_id == "P0":
        precise = ["Interneurons" if label == "Glia" else label for label in precise]
        integrated = labels_for_cells(payload, "integrated_cell_type")
        precise = [
            E17_INTEGRATED_TO_PRECISE[integrated_label]
            if label in {"CP", "CN"}
            else label
            for label, integrated_label in zip(precise, integrated)
        ]
        precise = [
            {
                "Molecular Layer Interneurons": "Interneurons",
                "Interneurons+Glia": "Glia/Oligodendrocytes",
                "Granule cells": "Inner Granule Layer",
                "i1": "i1 Neurons",
            }.get(label, label)
            for label in precise
        ]
        recode(payload, "finer_cell_types", precise, P0_PRECISE_ORDER, P0_PRECISE_COLORS)
        cn_labels = [
            "DCN" if precise_label == "excCN" else "Other" if cn_label == "DCN" else cn_label
            for precise_label, cn_label in zip(precise, cn_labels)
        ]
        recode(
            payload,
            "CN_exc_inhib",
            cn_labels,
            ["Other", "DCN", "i1", "Interneurons"],
            CN_COLOR_OVERRIDES,
        )
        return
    if dataset_id == "P4":
        recode(payload, "finer_cell_types", precise, color_overrides=P4_PRECISE_COLORS)
        return
    recode(payload, "finer_cell_types", precise)


def main():
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    for dataset in manifest["datasets"]:
        data_path = ROOT / dataset["data_url"].split("?", 1)[0]
        payload = json.loads(data_path.read_text(encoding="utf-8"))
        harmonize_dataset(dataset["id"], payload)
        data_path.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")
        dataset["data_url"] = f"{dataset['data_url'].split('?', 1)[0]}?v=20260916-7"
        print(f"Updated {dataset['id']}: {data_path.name}")
    MANIFEST_PATH.write_text(json.dumps(manifest, separators=(",", ":")), encoding="utf-8")


if __name__ == "__main__":
    main()
