# Local source data

Place the raw AnnData files here when regenerating the browser export.

The website does not commit the raw `.h5ad` files because they are too large for a normal GitHub repository. Run:

```sh
Rscript scripts/export_h5ad_site_data.R
```

to regenerate `assets/data/site-data.json` from `data/e12_celltypes_spatialclustering.h5ad`.
