# Cerebellar Development Transcriptomic Atlas

A static, browser-based atlas for exploring BARseq3 spatial transcriptomics and
scRNA-seq datasets across cerebellar development. The site is designed for
GitHub Pages and does not require a build step or server-side application.

## Project structure

```text
.
|-- index.html                 Page structure and accessible UI controls
|-- assets/
|   |-- styles.css             Shared visual system and responsive layouts
|   |-- app.js                 BARseq3 spatial explorer
|   |-- scrna.js               scRNA-seq explorer
|   |-- data/                  Browser-ready datasets, manifests, and counts
|   `-- images/                Site, schematic, and timepoint imagery
|-- data/                      Local raw-data staging area (not published)
`-- scripts/                   Export and validation utilities
```

## Run locally

Serve the repository root with any static HTTP server. Python is sufficient:

```sh
python -m http.server 4173
```

Then open `http://127.0.0.1:4173/`. Opening `index.html` directly is not
supported because the explorers load dataset files with `fetch()`.

## Validate before publishing

The project has no runtime dependencies. With Node.js 20 or newer, run:

```sh
npm run check
```

This checks the local assets referenced by `index.html`, both dataset manifests,
the browser-data files, and the sparse gene-count indexes. The original
projection-level data check remains available as:

```sh
npm run check:projection
```

## Data flow

- `assets/data/barseq-manifest.json` lists the published BARseq3 timepoints.
- `assets/data/scrna/manifest.json` lists the published scRNA-seq datasets.
- Browser-data JSON files contain coordinates and metadata used to draw plots.
- Sparse count indexes and shards provide on-demand all-gene visualization
  without shipping dense RDS or H5AD objects to every visitor.
- Raw source objects stay under `data/` and are excluded from Git.

See [data/README.md](data/README.md) for export commands and instructions for
adding datasets.

## Publishing

GitHub Pages serves the files in this repository directly. After changing data
or front-end files, run the validation command, commit the generated browser
assets, and push the branch used by the Pages workflow.
