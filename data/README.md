# Local source data

Place the raw AnnData files here when regenerating the browser export.

The website does not commit the raw `.h5ad` files because they are too large for a normal GitHub repository. Run:

```sh
Rscript scripts/export_h5ad_site_data.R
```

to regenerate `assets/data/site-data.json` from `data/e12_celltypes_spatialclustering.h5ad`.

Additional BARseq H5AD datasets can be exported with:

```sh
python scripts/export_h5ad_barseq_data.py /path/to/input.h5ad assets/data/site-data-stage.json "Stage BARseq3"
```

Then add the generated file to `assets/data/barseq-manifest.json` so it appears in the stage selector.

The scRNA-seq browser exports are generated from local Seurat `.rds` objects. They are intentionally
not committed because the raw objects are multi-gigabyte files. Run:

```sh
Rscript scripts/export_rds_scrna_data.R /path/to/rds/folder assets/data/scrna
```

to extract only UMAP coordinates, selected annotations, and lightweight QC values into browser-ready JSON.
