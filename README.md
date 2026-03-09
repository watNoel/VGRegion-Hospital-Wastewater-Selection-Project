# VGRegion-Hospital-Wastewater-Selection-Project
On selection toward antibiotic-resistant bacteria from hospital and municipal wastewater in the Västra Götaland Region of west Sweden

##  Statistical analysis of resistance rates in *E. coli* and ARG carriage rates 
To run the scripts for the statistical analysis of resistance rates in e.coli and ARG carriage rates,*ecoli_resistance_rates_models.R* and *metagenomic_arg_models.R* , the following R packages are required along with an R installation ( R version 4.4.0 tested on windows x64-based laptop). Typical installation time of the packages is around 5-10 minutes. 

| Package   | Version   | Purpose                                         |
|-----------|-----------|-------------------------------------------------|
| tidyverse | 2.0.0     | Data wrangling                                  |
| glmmTMB   | 1.1.14    | Negative binomial mixed effects modelling       |
| emmeans   | 2.0.1    | Estimated marginal means and contrasts          |
| DHARMa    | 0.4.7    | Checking residuals and model fit          |

For the script *ecoli_resistance_rates_models.R*, the input dataset is the *Merged Counts E. coli* sheet of the source data excel file in the submitted paper. The file is also present in the /inputs folder and named *merged_counts_ecoli.xlsx*:
The expected output files,, providing comprehensive statistical outputs from the modelling, are found under *outputs/ecoli_stats.xlsx*

For the script *metagenomic_arg_models.R*, metadata for the samples are needed, and a file with the ARG counts in the different samples.
These files are found under *inputs/metadata.xlsx* and *inputs/ResFinder_DB_counts_by_group.xlsx*. 
The expected output files, providing comprehensive statistical outputs from the modelling, are found under *outputs/metagenomics_stats.xlsx*


For the scripts to function, take care to update the paths to the different datasets in the scripts.
The scripts assume there is a folder called output in the directory where you run the code. If not, create such a folder or change the output destination of the created xlsx files.
NOTE: Prior to publication, all inputs and outtputs have been provided as encrypted files and a passkey provided to reviewers. Upon publication, this encryption will be removed.


## Metagenomic ARG Detection Pipeline

The script `arg_counter.py` processes Illumina paired-end metagenomic reads to identify and quantify antibiotic resistance genes (ARGs). The workflow includes read quality control, trimming, database preparation, sequence alignment against resistance gene databases, result filtering, ARG counting, and taxonomic profiling.

The pipeline was tested on a server running Red Hat Enterprise Linux (v7.9) consisting of 80 Intel(R) Xeon(R) Gold 6230 CPU @ 2.10GHz CPUs and 754 G RAM. The base programming language was **Python 3.9.15** and requires several Python libraries and external bioinformatics tools installed via **Conda 22.9.0**.

### Python Dependencies

| Library | Version | Purpose |
|---|---|---|
| pandas | 1.5.2 | DataFrame manipulation and processing of alignment results |
| numpy | 1.23.5 | Numerical operations and data processing |
| tqdm | 4.67.1 | Progress bars for long-running steps |

All other modules used in the script (e.g. `os`, `subprocess`, `re`, `logging`, `multiprocessing`, `gzip`, `tempfile`, `shutil`, `json`, `glob`) are part of the Python standard library.

### External Bioinformatics Tools

| Tool | Version | Purpose |
|---|---|---|
| FastQC | 0.11.9 | Quality control of raw and trimmed sequencing reads |
| MultiQC | 1.14 | Aggregation of QC reports across samples |
| BBMap (bbduk) | 39.01 | Adapter removal and quality trimming |
| VSEARCH | 2.22.1 | Sequence clustering and duplicate removal |
| EMBOSS (transeq) | 6.6.0 | Translation of nucleotide sequences to protein sequences |
| DIAMOND | 0.9.14 | Fast alignment of reads against protein reference databases |
| Sylph | 0.8.1 | Metagenomic taxonomic profiling |
| Sylph-tax | 1.5.1 | Taxonomic annotation and abundance summarization |

### Example Conda Environment Installation

The required Python libraries and bioinformatics tools can be installed using **Conda**. The following example creates an environment that contains the main dependencies required to run the pipeline for steps 1-6 (excluding taxanomic analysis). Specific steps in the pipeline can be specified by binary switching (i.e. `step_1=True`, `step_2=False`).

```bash
conda create -n arg_counter \
  python=3.9.15 \
  pandas=1.5.2 \
  numpy=1.23.5 \
  tqdm=4.67.1 \
  fastqc=0.11.9 \
  multiqc=1.14 \
  bbmap=39.01 \
  vsearch=2.22.1 \
  emboss=6.6.0 \
  diamond=0.9.14 \
  -c bioconda -c conda-forge
  ```


It is important to note that in order to run the taxononmic analysis, it is highly recommended to run steps 1-6 first. This ensures that all nessesary files and directories are generated. 

When running the taxonomic analysis `step_7` it is nessesary to create an independent conda environment.

```bash
conda create -n sylph \
  sylph=0.8.1 \
  sylph-tax=1.5.1 \
  -c bioconda -c conda-forge
  ```

### Pipeline Overview

The pipeline consists of several steps that can be enabled or disabled in the script:

1. **Initial quality control**  
   Raw sequencing reads are assessed using FastQC and summarized with MultiQC.

2. **Read trimming**  
   Adapter sequences and low-quality bases are removed using BBduk.

3. **Post-trimming quality control**  
   FastQC and MultiQC are run again to verify trimming performance.

4. **Database preparation**  
   Reference resistance gene sequences are clustered with VSEARCH, translated to protein sequences with EMBOSS Transeq, filtered for valid open reading frames, and formatted as DIAMOND databases.

5. **ARG detection**  
   Trimmed reads are aligned against the reference ARG database using DIAMOND BLASTX.

6. **ARG counting and phenotype annotation**  
   Alignment results are filtered, merged with resistance phenotype information, and summarized into counts per resistance class.

7. **Taxonomic profiling**  
   Metagenomic composition can be estimated using Sylph and Sylph-tax with the **GTDB_r226** reference database.

### Usage

The pipeline is configured directly within the script using several variables:

- `wdir` – project working directory  
- `data_dir` – directory containing raw FASTQ files  
- `out_dir` – output directory for results  
- `db_dir` – directory containing reference databases  

The python metagenomics pipeline file `arg_counter.py` should be placed in the working directory (`wdir`).

Ensure that within the `wdir` there exists a `data/` directory that contains `output/` where the results will be generated and the `dataDirectory/` folder that contains all the Illumina paired-end FASTQ files. Moreover, ensure that `databases/` also exists where the `ResFinder_DB/` is located that contains all the files downloadable from **ResFinder** `git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git`.

The pipeline expects Illumina paired-end FASTQ files (`*.fastq.gz`) with filenames following the pattern:

```bash
sampleID_L00X_R1_001.fastq.gz
sampleID_L00X_R2_001.fastq.gz
```

where `X` corresponds to sequencing lanes.

Pipeline steps can be toggled by modifying the boolean variables at the beginning of the script:

```python
#Quality check
step_1 = True
#Trim reads
step_2 = True
#Quality check after trimming
step_3 = True
#Create Databases
step_4 = False
#Create Resfinder DB
step_4_1 = False
#Run Diamond
step_5 = False
#Perform analys on count data
step_6 = False
# Taxonomic Profiling
step_7 = False
```
In order for successful completion of the pipeline the steps should be run in sequeuence.

**Example execution**:

Activate the environment before running the pipeline that excludes taxonomic analysis (steps 1-6):

```bash
conda activate arg_counter
```

Activate the environment before running the pipeline for only taxanomic analysis `step_7`:

```bash
conda activate sylph
```
**Steps 1-6 and 7 cannot be run in the same session due to the requirement of different active conda environments as shown above!**

And in the bash terminal:

```bash
python arg_counter.py
```

### Output

The pipeline produces multiple output directories containing:

* Quality control reports (FastQC, MultiQC)

* Trimmed sequencing reads

* DIAMOND alignment results

* Filtered ARG detection tables

* ARG counts summarized by resistance class

* Relative taxonomic abundance tables

The following files are produced and used for downstream analysis in R:

```
ResFinder_DB_counts_by_group.tsv
merged_taxonomy_abundance.tsv
```
---

## R Scripts: ARG Counts, Beta Diversity, and Taxonomy Analysis

This section describes the R scripts used for ARG counts, beta diversity, and taxonomic composition analyses. The tables list the required R packages, versions, and their purposes.

### **ARG_Counts_Script.R**

| Package        | Version | Purpose                                           |
|----------------|---------|-------------------------------------------------|
| ggplot2        | 3.5.2   | Data visualization, stacked bar plots          |
| dplyr          | x.x.x   | Data wrangling and transformation              |
| readr          | x.x.x   | File import and TSV/CSV reading                |
| tidyverse      | x.x.x   | Pipes and helper functions                      |
| stringi        | x.x.x   | String cleaning and manipulation               |
| RColorBrewer   | x.x.x   | Color palettes for plots                        |
| showtext       | x.x.x   | Embedding fonts (Arial) in plots               |
| svglite        | x.x.x   | SVG export of plots                             |

**Purpose:** Load metadata and ResFinder ARG counts, clean gene/class names, normalize counts to CPM, aggregate top ARG classes per sample group, and generate stacked bar plots.

### **Beta_Diversity_Script.R**

| Package        | Version | Purpose                                           |
|----------------|---------|-------------------------------------------------|
| vegan          | 2.7.1   | Bray–Curtis distance and diversity metrics     |
| ggplot2        | 3.5.2   | PCoA visualization                              |
| dplyr          | x.x.x   | Data manipulation                               |
| tibble         | x.x.x   | Data frame handling                             |
| pairwiseAdonis | x.x.x   | Pairwise PERMANOVA                              |
| tidyverse      | x.x.x   | General data wrangling                          |
| ggpubr         | x.x.x   | Plot annotation and formatting                  |
| grid           | x.x.x   | Plot layout adjustments                         |
| svglite        | x.x.x   | Export PCoA figures in SVG                      |

**Purpose:** Compute Bray–Curtis distances, perform PCoA, PERMANOVA, beta-dispersion, pairwise tests, and generate PCoA plots with convex hulls representing sample groups.

### **Taxonomy_Script.R**

| Package        | Version | Purpose                                           |
|----------------|---------|-------------------------------------------------|
| tidyverse      | x.x.x   | Data wrangling, pivoting, plotting             |
| ggsci          | x.x.x   | Scientific color palettes                       |
| randomcoloR    | x.x.x   | Generate distinct colors for taxa              |
| scales         | x.x.x   | Axis formatting for plots                       |

**Purpose:** Load taxonomy and metadata, extract taxa at a specified level (e.g., phylum), compute relative abundance, select top N taxa per group, generate stacked bar plots with consistent color and factor ordering.


## Usage

All three scripts are designed to be run in R (tested on **R 4.5.1** within the **RStudio IDE 2025.05.1, Build 513**). Typical usage steps:

1. Set the working directory in each script to where your metadata and input files are stored:
   ```r
   setwd("path/to/files")
    ```
2. Ensure the metadata CSV and taxonomy/ARG input files exist and are named correctly within the same working directory:

* Metadata: `metadata.csv` - Found as part of the supplementry files.

* ARG counts: `ResFinder_DB_counts_by_group.tsv` - Output of the metagenomics pipeline.

* Taxonomy: `merged_taxonomy_abundance.tsv` - Output of the metagenomics pipeline.

* Phenotypes mapping: `phenotypes.txt` - File from the ResFinder repository.

Optional: Adjust parameters within the scripts:

* `top_n` (`Taxonomy_Script.R`) – number of top taxa to visualize

* `tax_level` (`Taxonomy_Script.R`) – e.g., `"p__"` for phylum

* Filtering by ARG class/gene (`ARG_Counts_Script.R`)

4. Run the scripts in R or RStudio:

```R
source("ARG_Counts_Script.R")
source("Beta_Diversity_Script.R")
source("Taxonomy_Script.R")
```

### Output
Each script generates plots and summary tables:

**ARG_Counts_Script.R**

* Stacked bar plots of top ARG classes per sample group (`Plots/arg_counts.tiff`, `Plots/arg_counts.svg`)

* Normalized data and summary tables embedded in R objects (`normalized_data`, `summary_data`)

**Beta_Diversity_Script.R**

* PCoA plots with convex hulls (`Plots/beta_div.tiff`, `Plots/beta_div.svg`)

* PERMANOVA and beta-dispersion results in R objects

* Optional pairwise PERMANOVA tables

**Taxonomy_Script.R**

* Stacked bar plots of relative abundances at the selected taxonomic level (`Plots/tax_phyla.tiff`, `Plots/tax_phyla.svg`)

* Summarized relative abundance tables per group embedded in R objects (`complete_df_summarized`)

* All figures are exported as both TIFF (high-res) and SVG (vector) formats.
