# **Immune hotspot analysis**

Scripts for data analysis conducted for paper “Spatial Positioning of Immune Hotspots Reflects the Interplay between B and T Cells in Lung Squamous Cell Carcinoma”.

## Implementation

- Clone the repository

  ```bash
  git clone https://github.com/yuerua/t_b_spatial_IH.git
  ```

- Download data to the same directory. 

  ```bash
  data
  ├── TCGA
  │   ├── Bcell_score_LUSC.RData
  │   ├── TCGA_LUAD-LUSC_path_EST_ABS_BoLi_DAV.RData
  │   ├── TCGA_LUAD-biolinksGeneExp_immuneHotspot.RData
  │   ├── TCGA_LUSC-biolinksGeneExp_immuneHotspot.RData
  │   ├── lusc_mutload.RData
  │   └── survival_LUSC_LUAD.RData
  └── validation
      └── validation.RData
  ```

- Run `main.rmd` to reproduce the results and figures

- Find scripts for calculating immune hotspots, tumour mutation burden (TMB), and classification accuracy in corresponding folders.

  



