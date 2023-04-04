# Meta-analysis of gene activity contributions and correlation with gene expression, through GAGAM.

This repo includes all the material and methods supporting the unpublished work "Meta-analysis of gene activity contributions and correlation with gene expression, through GAGAM"

One can find the Preprint at ...

## Repo structure and organization

~~~~
├── README.md
├── DATA                                                   //Data files supporting GAGAM construction and meta-analysis
│   ├── 10x_PBMC_Multiome_Controller                       //Starting data files to support GAGAM construction
│   │   ├── conns                                          //Housekeeping genes lists for human and mouse for RAGI computation
│   │   └── Labeled peaks                                  //labeled peak file for GAGAM construction
│   └── Genomes                                            //Genomic assemblies data for GAGAM construction
├── RESULTS                                                //IMAGES and Tables output of the meta-analysis
│   ├── IMAGES                                             
│   │   ├── ENHANCER                                       //A-E plots for the Enhancer contribution
│   │   ├── EXON                                           //A-E plots for the Exon contribution
│   │   └── PROM                                           //A-E plots for the Promoter contribution                         
│   └── Tables                                             //Correlation tables, coherence tables, peak information table           
└── Script
    └── MAGA.R                                         //R scripts for GAGAM construction and meta-analysis
~~~~

## Script Usage

Only one script is necessary to obtain the results, the IWBBIO23.r file.
It requires to download the data from the source at:
- 10X Multiome Controller PBMC: https://www.10xgenomics.com/resources/datasets/10-k-human-pbm-cs-multiome-v-1-0-chromium-controller-1-standard-2-0-0

Organize the downloaded files in the folder:
- 10x_PBMC_Multiome_Controller
  - filtered_feature_bc_matrix
    - barcodes.tsv
    - matrix.mtx
    - features.bed

Then from inside the Script folder, execute the script with the syntax:
```
Rscript MAGA.R
```

It will output all the Figures and Tables in the TMPResults which has the same structure as the Results folder in this repo.



