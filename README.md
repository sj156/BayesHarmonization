# BayesianHarmonization-

## SuLY_EDA

Exploratory visualization of the SuLY-MPRAGE traveling-subject data
(276 skull-stripped, co-registered MPRAGE scans; 17 subjects; 5 scanners;
protocol settings Param01–05). Raw `.nii.gz` data is not included.

| File / folder | Content |
|---|---|
| `00_data_overview.ipynb` | Data overview notes |
| `suly_visualization.R` | Raw-intensity plots (magma colormap) → `viz_png/` |
| `suly_visualization_standardized.R` | Per-scan standardized plots, `x / median(brain) - 1` or `x / mean(brain) - 1`, blue–white–red → `viz_png_standardized_median/`, `viz_png_standardized_mean/` |

Each output folder is organized as `<subject>/<scanner>/<param>/` and holds three
plots per scan: `_three_views.png` (mid-slice sagittal/coronal/axial),
`_axial_montage.png` (25 axial slices), `_mip_axis3.png` (maximum intensity projection).

Scripts expect the data at `../SuLY-MPRAGE-Handoff/` relative to `SuLY_EDA/`:

```
Rscript SuLY_EDA/suly_visualization.R [sub-XXXX]
Rscript SuLY_EDA/suly_visualization_standardized.R [sub-XXXX] [median|mean]
```
