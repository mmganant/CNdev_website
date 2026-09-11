#!/usr/bin/env python3
"""Export an H5AD matrix as browser-loadable, gene-sharded sparse counts."""

import gzip
import json
import struct
import sys
from pathlib import Path

import h5py
import numpy as np


def decode(values):
    return [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values]


def read_genes(handle):
    for field in ("_index", "Names"):
        if field in handle["var"]:
            return decode(handle[f"var/{field}"][:])
    raise KeyError("No gene-name field was found in var")


def export(source: Path, output: Path, matrix_path="layers/counts", genes_per_shard=64):
    output.mkdir(parents=True, exist_ok=True)
    with h5py.File(source, "r") as handle:
        matrix = handle[matrix_path]
        genes = read_genes(handle)
        n_cells, n_genes = matrix.shape
        shards = []
        for start in range(0, n_genes, genes_per_shard):
            end = min(start + genes_per_shard, n_genes)
            filename = f"counts-{start // genes_per_shard:04d}.bin.gz"
            with gzip.open(output / filename, "wb", compresslevel=9) as stream:
                stream.write(struct.pack("<I", end - start))
                block = np.asarray(matrix[:, start:end], dtype=np.float32)
                for column in range(block.shape[1]):
                    values = block[:, column]
                    indices = np.flatnonzero(values).astype("<u4", copy=False)
                    stream.write(struct.pack("<I", len(indices)))
                    stream.write(indices.tobytes())
                    stream.write(values[indices].astype("<f4", copy=False).tobytes())
            shards.append({"file": filename, "start": start, "count": end - start})
            print(f"Wrote {filename} ({end}/{n_genes} genes)")

    index = {
        "dataset": source.stem.replace("_", "-"),
        "source_file": source.name,
        "matrix": matrix_path,
        "encoding": "gzip-dense-records-v1",
        "value_type": "float32",
        "index_type": "uint32",
        "n_cells": n_cells,
        "n_genes": n_genes,
        "genes": genes,
        "shards": shards,
    }
    (output / "index.json").write_text(json.dumps(index, separators=(",", ":")), encoding="utf-8")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise SystemExit("Usage: export_h5ad_sparse_counts.py SOURCE.h5ad OUTPUT_DIR [MATRIX_PATH]")
    export(Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3] if len(sys.argv) > 3 else "layers/counts")
