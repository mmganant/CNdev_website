#!/usr/bin/env python3
"""Build lightweight spatial previews from exported BARseq site JSON files."""

import json
import sys
from pathlib import Path


def main(output, specs):
    previews = {}
    for dataset_id, source in specs:
        data = json.loads(Path(source).read_text(encoding="utf-8"))
        field = "finer_cell_types"
        field_index = data["schema"].index(field)
        step = max(1, len(data["cells"]) // 6000)
        previews[dataset_id] = {
            "categories": data["annotations"].get(field, []),
            "cells": [[cell[0], cell[1], cell[field_index]] for cell in data["cells"][::step]],
        }
    Path(output).write_text(json.dumps({"datasets": previews}, separators=(",", ":")), encoding="utf-8")


if __name__ == "__main__":
    if len(sys.argv) < 4 or len(sys.argv[2:]) % 2:
        raise SystemExit("Usage: export_barseq_previews.py OUTPUT.json ID FILE [ID FILE ...]")
    main(sys.argv[1], list(zip(sys.argv[2::2], sys.argv[3::2])))
