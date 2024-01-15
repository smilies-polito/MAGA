# Meta-analysis of gene activity contributions and correlation with gene expression, through GAGAM.

This repo includes all the material and methods supporting the published work "Meta-analysis of gene activity contributions and correlation with gene expression, through GAGAM", and the unpublished work "Cross-omic Transcription Factors meta-analysis: an insight on TFs accessibility and expression correlation". 
One can find the first at https://link.springer.com/chapter/10.1007/978-3-031-34960-7_14, and the second at .

## Repo structure and organization

~~~~
├── README.md
├── DATA                                                   //Data files supporting GAGAM construction and meta-analysis
│   ├── 10x_PBMC_Multiome_Controller                       //Starting data files to support GAGAM construction
│   │   ├── conns                                          //Housekeeping genes lists for human and mouse for RAGI computation
│   │   └── Labeled peaks                                  //labeled peak file for GAGAM construction
│   └── Genomes                                            //Genomic assemblies data for GAGAM construction
├── RESULTS
│   ├── MAGA                                                //IMAGES and Tables output of the meta-analysis
│   │   ├── IMAGES                                             
│   │   │   ├── ENHANCER                                    //A-E plots for the Enhancer contribution
│   │   │   ├── EXON                                        //A-E plots for the Exon contribution
│   │   │   └── PROM                                        //A-E plots for the Promoter contribution                         
│   │   └── Tables                                          //Correlation tables, coherence tables, peak information table     
│   └── MAGA TF                                             //IMAGES of the TFs meta-analysis
│       ├── IMAGES                                          //IMAGES and Tables output of the meta-analysis 
│       │   ├── FOOT_ENHD                                   //FTP-E plots for the Enhancer contribution
│       │   ├── FOOT_ENHD_GENE                              //FTP-E plots for the Enhancer contribution (cell-types in genes)
│       │   ├── FOOT_PROM                                   //FTP-E plots for the Promoter contribution
│       │   ├── FOOT_PROM_GENE                              //FTP-E plots for the Promoter contribution (cell-types in genes)
│       │   ├── MOTIF_ENHD                                  //ME-E plots for the Enhancer contribution
│       │   ├── MOTIF_PROM                                  //ME-E plots for the Promoter contribution
│       │   └── MISC                                        //Miscellaneous images produced and described in the article                        
│       ├── motif_activity.RDS                              //R object containing all the processed data
│       ├── TF_footprint_prom                               //TFs Footprint scores for all the TFs in promoters
│       └── TF_footprint                                    //TFs Footprint scores for all the TFs in enhancers
└── Script
    ├── MAGA_TF.R                                           //R scripts for TF meta-analysis
    ├── chromvar.R                                          //R scripts for running the necessary chromvar analyses
    └── MAGA.R                                              //R scripts for GAGAM construction and meta-analysis
~~~~

## Script Usage

### MAGA
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

### MAGA TF
Two scripts are necessary to obtain the results, MAGA_TF.r file.
It requires to download the data from the source at:
- 10X Multiome Controller PBMC: https://www.10xgenomics.com/resources/datasets/10-k-human-pbm-cs-multiome-v-1-0-chromium-controller-1-standard-2-0-0

Organize the downloaded files in the folder:
- 10x_PBMC_Multiome_Controller
  - filtered_feature_bc_matrix
    - barcodes.tsv
    - matrix.mtx
    - features.bed

Then from inside the main folder, execute the script with the syntax:
```
Rscript Script/MAGA_TF.R
Rscript Script/chromvar.R
```

It will output all the Figures and Tables in the TMPResults which has the same structure as the Results folder in this repo.



