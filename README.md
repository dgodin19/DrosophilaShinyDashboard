# Thermal Plasticity RNA-seq Explorer

An R Shiny dashboard for exploring RNA-seq gene expression data, running differential expression analysis, and visualizing gene set enrichment results from a study of thermal acclimation in *Drosophila melanogaster* from Salachan et al. (2022).

## Background

This dashboard explores an RNA-seq dataset generated to study thermal acclimation in *Drosophila melanogaster*. The experiment compared flies acclimated to constant (19 ± 0°C) versus fluctuating (19 ± 8°C) temperatures, with genome-wide transcriptomic profiling carried out across three biological groups: adult males, adult females, and larvae.

The dataset consists of RNA-seq gene expression counts across 80 samples, with metadata capturing:

- **treatment** — constant vs. fluctuating temperature acclimation
- **sex** — male vs. female (adults)
- **timepoint** — sampling time
- **lifestage** — adult vs. larva
- **replicate** — biological replicate ID

Two counts matrices are used: a **CPM-normalized matrix** for exploratory analysis and a **raw counts matrix** for differential expression with `limma voom`. Because the number of factor levels varies across treatment/sex/lifestage/timepoint combinations, the differential expression model design is built dynamically from whichever metadata variables the user selects, rather than assuming a fixed set of contrasts. Gene set enrichment results were generated separately with `fgsea`, ranking genes by their differential expression statistics, and are loaded into the app alongside the expression data.

This dashboard was built to interactively explore that dataset: the sample metadata, the normalized/raw count matrices, the differential expression results, and the fgsea gene set enrichment results.

## What the App Does

The app is organized into tabs covering data exploration, differential expression analysis, and enrichment visualization for 80 RNA-seq samples with metadata on **treatment, sex, timepoint, lifestage,** and **replicate**.

### Sample Information Exploration
Explore the metadata before running any analysis.
- **Summary** — row/column counts, data types, and distinct values per column
- **Sample Data Table** — sortable, searchable table (DT)
- **Variable Visualization** — bar plots of metadata variable distributions, single variable or all at once

### Counts Matrix Exploration
Explore the structure of the normalized (CPM) counts matrix and test filtering strategies.
- **Counts Table** — preview of the filtered normalized counts
- **Filtering Summary** — sample count, total genes, genes passing filters, genes removed
- **Filtering Diagnostics** — scatter plots of median expression vs. variance and median expression vs. number of zero counts, with genes passing filters highlighted
- **Clustered Heatmap** — top 30 most variable genes, optional log transform, clustered by gene and sample

### PCA
Principal component analysis on the CPM-normalized, filtered counts.
- Choose X and Y principal components
- Color samples by any metadata variable
- Axes labeled with percent variance explained

### Differential Expression
Runs `limma voom` on the raw counts matrix.
- **Results Table** — log fold change, average expression, p value, adjusted p value
- **P value Histogram** — distribution of p values across genes
- **Log2 Fold Change Plot** — distribution of effect sizes
- **MA Plot** — log fold change vs. average expression
- **Volcano Plot** — highlights genes with strong significance and fold change

### Gene Set Enrichment Analysis
Uses `fgsea` results generated separately, ranked by differential expression statistics.
- **Top Pathways** — bar plot of most significant pathways by adjusted p value, slider to control how many are shown, click a pathway to see its table entry
- **Results Table** — sortable, filterable by adjusted p value threshold and enrichment direction (positive/negative), with download support
- **NES vs. Significance Plot** — normalized enrichment score vs. −log10 adjusted p value, significant pathways highlighted

## Analysis Pipeline

1. Raw and CPM-normalized counts matrices, sample metadata, and fgsea results are uploaded to the app.
2. Metadata variables (treatment, sex, timepoint, lifestage) are used to build the differential expression model design.
3. Differential expression is computed with `limma voom`.
4. Gene set enrichment is computed separately with `fgsea`, ranking genes by the differential expression statistic.

## Shiny Functionalities Used

- File uploads for counts matrix, metadata, raw counts, and fgsea results
- Reactive expressions for filtering and analysis updates
- Sliders for filtering thresholds and pathway display counts
- Select inputs for PCA components and metadata grouping
- Interactive tables via DT
- Dynamic plots that update as inputs change

## Known Limitations / Design Notes

- **Design contrasts:** Attempting to specify explicit contrasts between groups ran into `"contrasts can be applied only to factors with 2 or more levels"`, since the number of factor levels varies by dataset. The current workaround lets the user pick a design formula, and the analysis uses the last coefficient from the model output rather than a user-specified contrast.
- **App structure:** Nested tab panels and reactive expressions were used to keep multiple analyses organized and outputs updating dynamically without cluttering the interface.

## Data Requirements

To run the app, you'll need:
- A sample metadata table (treatment, sex, timepoint, lifestage, replicate)
- A CPM-normalized counts matrix
- A raw counts matrix (for differential expression)
- fgsea gene set enrichment results

An example set of data is in the data folder. 

## References

Salachan, Paul Vinu, and Jesper Givskov Sørensen. 2022. “Molecular Mechanisms Underlying Plasticity in a Thermally Varying Environment.” Molecular Ecology, April. https://doi.org/10.1111/mec.16463.

An example set of data is in the data folder. 
