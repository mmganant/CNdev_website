#!/usr/bin/env python3
"""Apply curated BARseq label harmonization to browser JSON exports."""

import json
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
    "excCN": "#ffa500",
    "i1": "#ffff00",
    "Interneurons": "#008000",
    "Unknown": "#bbbbbb",
}

E17_PRECISE_ORDER = list(PRECISE_COLORS)
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


def recode(payload, field, labels, preferred_order=None):
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
    payload["annotations"][field] = [
        {
            "label": label,
            "count": counts[label],
            "color": PRECISE_COLORS.get(label, previous.get(label, "#bbbbbb")),
        }
        for label in order
    ]


def harmonize_dataset(dataset_id, payload):
    cn_labels = labels_for_cells(payload, "CN_exc_inhib")
    cn_labels = ["Interneurons" if label == "i2/3" else label for label in cn_labels]
    recode(payload, "CN_exc_inhib", cn_labels, ["Other", "DCN", "i1", "Interneurons"])

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
    if dataset_id == "P0":
        precise = ["Interneurons" if label == "Glia" else label for label in precise]
    recode(payload, "finer_cell_types", precise)


def main():
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    for dataset in manifest["datasets"]:
        data_path = ROOT / dataset["data_url"].split("?", 1)[0]
        payload = json.loads(data_path.read_text(encoding="utf-8"))
        harmonize_dataset(dataset["id"], payload)
        data_path.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")
        dataset["data_url"] = f"{dataset['data_url'].split('?', 1)[0]}?v=20260916-1"
        print(f"Updated {dataset['id']}: {data_path.name}")
    MANIFEST_PATH.write_text(json.dumps(manifest, separators=(",", ":")), encoding="utf-8")


if __name__ == "__main__":
    main()
