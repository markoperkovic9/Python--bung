# Explanation of some points in ArchR-Guide


Quality Control (QC) in ArchR

In scATAC-seq, we don't just look at "how many reads" a cell has. We look at Signal vs. Noise.
Two Main Metrics:

    TSS Enrichment Score: This measures how many fragments fall near Transcription Start Sites (TSS) compared to the "background" genomic noise. High TSS enrichment = high-quality data.

    Fragment Size Distribution: Because DNA is wrapped around nucleosomes, a good ATAC-seq library should show a distinct "pattern" (peaks representing fragments that are sub-nucleosomal, mono-nucleosomal, etc.).


What exactly are Arrow files?

![Erstes Bild](image.png)

1. Left Side: Parallel Addition of Data (The "Writing" Process)
When you run a command like addTileMatrix() (which counts fragments in 500-bp windows across the genome), ArchR doesn't try to calculate everything at once.

Chromosome-wise Processing: ArchR looks at the Accessible Fragments already stored in the Arrow file. These fragments are physically organized by chromosome (chr1, chr2, etc.).

Multithreading: If you have 8 CPU cores, ArchR can process 8 chromosomes simultaneously. Each core calculates the matrix for its specific chromosome in isolation.

Direct-to-Disk: Once a chromosome's matrix is calculated, it is written directly back into the HDF5-based Arrow file on your hard drive.

The Benefit: You never have the full genome-wide matrix of all cells in your RAM at once. You only ever hold one or two chromosomes at a time, keeping your "memory footprint" very low.

2. Right Side: Parallel Retrieval of Data (The "Reading" Process)
Later in your analysis, you might want to extract specific data (e.g., marker genes or specific peaks) to create a plot or perform a statistical test.

Parallel Feature Filtration: ArchR reads the Feature Matrix from the Arrow file. Again, it does this in parallel across chromosomes.

The "Funnel" (Filtering): Not all data is needed for every step. ArchR filters the data at the chromosome level (e.g., "give me only the cells in Cluster 1" or "give me only these 500 specific genes").

Matrix Blocks: It extracts these small "blocks" of data from the disk.

Consolidation: These small, filtered blocks are finally merged into a single Filtered Matrix in your R session.

The Benefit: Because you've filtered the data before bringing it fully into R, the resulting object is small enough to work with comfortably for plotting and downstream analysis.

## Input File Types in ArchR

To build Arrow files, ArchR requires specific input formats. While it can handle several types, some are much more efficient than others.

1. Fragment Files (Highly Recommended)
The fragment file is the standard input for ArchR. It is typically produced by pipelines like CellRanger ATAC or SnapATAC.

Format: .tsv.gz (tab-separated, compressed) accompanied by a .tbi (index) file.

Content: A table where each row represents a single DNA fragment captured by the assay.

Columns: Usually Chromosome, Start, End, CellBarcode, and DuplicateCount.

Why it's preferred: It is much smaller than a BAM file and contains exactly what ArchR needs: the start and end positions of the Tn5 transposase insertion.

2. BAM Files
ArchR can also process aligned BAM files.

Format: .bam

Drawback: BAM files are much larger and contain more information than necessary for ATAC-seq (like the full sequence of the read).

Process: If you provide BAM files, ArchR will internally convert them into fragment data during the Arrow file creation. This takes more time and disk space.

3. Genome Annotations
Since ArchR needs to know where genes and transcription start sites (TSS) are, you must provide genomic coordinates.

BSgenome: ArchR uses Bioconductor BSgenome packages (e.g., BSgenome.Hsapiens.UCSC.hg38 for humans).

GeneAnnotation: Information about gene boundaries and exons.

GenomeAnnotation: Information about chromosome sizes.

4. Blacklist Regions
These are "problematic" areas of the genome where the signal is often artificially high (due to repetitive sequences or mapping errors).

Function: ArchR uses these coordinates to filter out "noise" that could lead to false-positive peaks.

Source: Usually provided by the ENCODE project.

### Summary of Input File Types

| File Type | Extension | Role |
| :--- | :--- | :--- |
| **Fragment File** | `.tsv.gz` | Primary data; contains Tn5 insertion sites per cell. |
| **Index File** | `.tbi` | Required for fast random access to the fragment file. |
| **BAM File** | `.bam` | Alternative primary data (aligned reads); slower to process. |
| **Genome Reference** | `BSgenome` | Maps coordinates to a specific species (e.g., hg38, mm10). |
| **Blacklist** | `.bed` / `GRanges` | Regions to ignore to reduce background noise. |


## 3.5 Setting the Seed (Reproducibility)

ArchR uses randomized algorithms for tasks like **UMAP**, **Clustering**, and **Doublet Inference**. To ensure results are reproducible:

* **The Command:** `set.seed(1)`
* **Purpose:** Ensures that the "random" parts of the code always produce the same output.
* **Best Practice:** Always set the seed at the beginning of your script and keep it consistent throughout your project.

> **Thesis Tip:** In your "Materials and Methods" section, it is professional to mention: *"Analyses were performed using ArchR (vX.X) with a fixed random seed of 1 to ensure reproducibility."*

## 3.5.2 HDF5 File Locking (HPC Specifics)

When working on a University HPC, the file system (Lustre, GPFS) often conflicts with HDF5's default locking mechanism.

* **Problem:** Network-based storage may not support "File Locking," causing ArchR to crash during Arrow file creation.
* **Solution:** Use `addArchRFileLocking(locking = FALSE)` if the HPC throws HDF5-related errors.
* **Consequence:** Disabling locking allows for **subthreading** (increased performance) but requires caution to ensure only one process writes to an Arrow file at a time.

> **HPC Checklist:** > 1. Check if the filesystem is Lustre or GPFS.
> 2. Always include `addArchRFileLocking()` in the setup block of the R script.

## 3.5.3 Setting the Reference Genome

ArchR requires a biological "map" to interpret genomic coordinates. 

* **Command:** `addArchRGenome("hg38")` (or "mm10" for mouse).
* **What it provides:** 1. **Chromosome sizes:** Essential for tiling the genome.
    2. **Gene locations:** Required to link ATAC-seq peaks to specific genes.
    3. **TSS locations:** Critical for Quality Control (QC) metrics.

> **HPC Tip:** If your cluster has no internet access, ensure the corresponding `BSgenome` R-package is pre-installed in your library before running this command.

## 3.5.4 Organizing Input Files

To process your own data, you need to create a named character vector in R. This vector tells ArchR two things:

* Where the data is located (the file path).

* What each sample should be called (the sample name).

Before creating Arrow files, you must create a named vector in R that points to your raw fragment files.

* **Variable:** `inputFiles`
* **Structure:** A character vector where the **value** is the path and the **name** is the Sample ID.
* **Impact:** The names used here will become the filenames for your `.arrow` files (e.g., `Sample1.arrow`).

### Example Code for Thesis Data:
```r
# Point to your actual fragment files on the HPC storage
inputFiles <- c(
  "Control_Rep1" = "/path/to/control_fragments.tsv.gz",
  "Stimulated_Rep1" = "/path/to/stimulated_fragments.tsv.gz"
) 
``` 

HPC Best Practice: Always use Absolute Paths (starting with /) to ensure your script finds the data when submitted as a batch job (SLURM/LSF).

## 3.5.5 Environment Setup: Genome & Parallelization

Every ArchR session must begin by defining the reference genome and the number of threads for parallel processing.

### 1. Reference Genome
* **Command:** `addArchRGenome("hg19")`
* **Purpose:** Provides gene coordinates and TSS locations. 
* **Note:** Ensure this matches the genome used for initial read alignment.

### 2. Parallel Threads
* **Command:** `addArchRThreads(threads = 16)`
* **Purpose:** Speeds up "per-chromosome" operations (e.g., TileMatrix creation).
* **HPC Rule:** The number of threads in R **must match** the number of cores requested in your HPC batch script (e.g., `--cpus-per-task`).

```r
# Example Setup for HPC
addArchRGenome("hg19")           # Matches tutorial data
addArchRThreads(threads = 16)    # Matches requested HPC cores
```

Warning: If you set threads = 32 in R but only requested cpus-per-task=1 from the HPC, your job will be extremely slow (due to "context switching") or might even be killed by the HPC for resource violation.

Thesis Tip: In your Methods section, always report the genome version (e.g., GRCh37/hg19) and mention that parallel processing was used to handle the large data scale. 
For your own actual thesis data, you will likely use hg38 (the current standard). Always ensure your reference genome matches the one used during the alignment step (e.g., by CellRanger).

## 3.6 Processing Raw Data: createArrowFiles()

This step converts raw fragment files into `.arrow` files. It is the "Quality Control" gatekeeper of the pipeline.

### Function Parameters:
| Parameter | Value | Description |
| :--- | :--- | :--- |
| `minTSS` | 4 | Min. Transcription Start Site enrichment. Ensures high signal-to-noise ratio. |
| `minFrags` | 1000 | Min. unique fragments per cell. Filters out "empty" or poorly sequenced cells. |
| `addTileMat` | TRUE | Generates a 500-bp bin matrix for genome-wide accessibility. |
| `addGeneScoreMat` | TRUE | Predicts gene expression based on local chromatin accessibility. |

### Why 500-bp bins?
ArchR uses a **"Base-pair Resolution"** approach but stores data in 500-bp tiles for computational efficiency. This provides a good balance between detail and speed for dimensionality reduction.

> **HPC Performance Note:** This step will take longer on your HPC if you have many samples. Ensure you have requested enough **time** in your SLURM script (e.g., 2–4 hours for a large dataset) even though the tutorial only takes 15 minutes.

## 3.6.1 The ArrowFiles Object (The Pointer)

A key concept in ArchR is that the data stays on the disk, not in the RAM.

* **Definition:** The `ArrowFiles` object is a **named character vector** containing the file paths to the physical `.arrow` files.
* **Mechanism:** ArchR uses this vector to know where to read/write data during downstream steps.
* **Memory Efficiency:** This allows the user to handle massive datasets because R only "points" to the data rather than "loading" it.

> **Note:** If you move your `.arrow` files to a different folder, this vector will break. You would need to update the paths in the vector for ArchR to find them again.

## 3.7 Per-cell Quality Control (QC)

Quality control is the most critical step in scATAC-seq analysis. It ensures that the biological signals we observe are not driven by technical noise, dead cells, or sequencing artifacts. ArchR evaluates every cell based on three main pillars.

---

### 3.7.1 Unique Nuclear Fragments
Instead of looking at total sequencing reads, ArchR focuses on **Unique Nuclear Fragments** to measure the usable data per cell.

* **Mitochondrial Filtering:** Fragments mapping to the mitochondrial genome (mtDNA) are removed. Since mtDNA is not wrapped in histones, it is highly accessible and can overwhelm the signal without providing info on nuclear regulation.
* **Deduplication:** Identical fragments (starting and ending at the same position) are collapsed into a single count to remove PCR bias.
* **Thresholds:** A standard minimum is **> 1,000 unique fragments**. Cells below this are too "sparse" to be accurately assigned to a cell type.



---

### 3.7.2 TSS Enrichment Score (Signal-to-Background)
The Transcription Start Site (TSS) Enrichment Score is the gold-standard metric for scATAC-seq quality. It measures how "concentrated" the Tn5 cuts are at known promoters compared to the genomic background.

* **The Biology:** Healthy cells have open, accessible chromatin at TSS regions. In dead or dying cells, chromatin collapses (de-chromatinization), leading to random Tn5 cutting across the entire genome.
* **The Calculation:**
  ArchR calculates the ratio of fragments at the center of the TSS vs. the "flanks" (background noise).
  
  $$\text{TSS Enrichment Score} = \frac{\text{Mean signal in the TSS center (+/- 50 bp)}}{\text{Mean signal in the flanking background regions}}$$

* **Interpretation:**
  * **Score > 10:** Excellent signal.
  * **Score 4 - 10:** Typical/Acceptable quality.
  * **Score < 4:** High background noise (likely dead cells).



---

### 3.7.3 Fragment Size Distribution (Nucleosomal Periodicity)
Because genomic DNA is wrapped around **nucleosomes** (~147 bp of DNA per wrap), the lengths of the fragments generated by Tn5 follow a specific pattern.

* **The Pattern:** You should see clear "peaks" in your distribution plot:
    * **Sub-nucleosomal (< 100 bp):** Cuts made between nucleosomes.
    * **Mono-nucleosomal (~200 bp):** Fragments spanning one nucleosome.
    * **Multi-nucleosomal (> 400 bp):** Fragments spanning two or more nucleosomes.
* **Why it matters:** This "staircase" pattern (periodicity) proves that the Tn5 enzyme successfully captured the natural chromatin architecture of the cell.



---

### 3.7.4 QC Summary Table 

| Metric | Threshold | Indication of Quality |
| :--- | :--- | :--- |
| **Unique Fragments** | > 1,000 - 2,500 | Sufficient data depth for dimensionality reduction. |
| **TSS Score** | > 4 (Ideal > 7) | High signal-to-noise ratio; healthy cells. |
| **Nucleosomal Periodicity** | Visible Peaks | Integrity of chromatin structure is preserved. |
| **Mitochondrial Rate** | < 10% | Effective cell lysis and high nuclear purity. |

> **Thesis Tip:** In your Methods section, explicitly state the thresholds used for `minTSS` and `minFrags`. These are often the first parameters a reviewer will check to validate your cell-filtering strategy.

### 3.7.5 Visualizing Quality Control (QC) Results

ArchR automatically generates a `QualityControl` folder. Each sample is evaluated by two primary visualizations to ensure biological signal exceeds technical noise.

View plots here: https://www.archrproject.com/bookdown/per-cell-quality-control.html

#### 1. TSS Enrichment vs. Fragments
This plot identifies high-quality cells. We look for a high density of cells above our TSS threshold (Y-axis) and fragment threshold (X-axis).

For the purpose of this guide, we will look at the data for `PBMC` 

* **TSS Threshold:** Indicated by the horizontal dashed line.
* **Frag Threshold:** Indicated by the vertical dashed line.
* **Observation:** Our tutorial samples show exceptional quality, with median TSS scores ranging from ~15 to ~20.

#### 2. Fragment Size Distribution
This plot confirms that the library preparation successfully captured the nucleosomal structure of the DNA.


* **Sub-nucleosomal Peak:** Highest signal at <100 bp.
* **Mono-nucleosomal Peak:** Clear shoulder/peak at ~200 bp.
* **Interpretation:** The presence of these discrete peaks (periodicity) is a hallmark of a successful scATAC-seq experiment.

> **Thesis Note:** These plots are the "ID card" for your data. Always include the median TSS enrichment and the total number of cells passing filter for every sample in your manuscript.

## 3.8 Using BAM Files for Arrow File Creation

Although **Fragment Files** are the gold standard, ArchR supports **BAM files** if necessary.

### 1. Identifying the Cell Barcode
Unlike fragment files, where the barcode is a plain column, BAM files store barcodes in metadata tags.
* **10x Genomics:** Usually the `CB` tag (e.g., `CB:Z:GCGGGTTAGAACGTCG-1`).
* **Command Parameter:** Use `bcTag = "CB"` in `createArrowFiles()`.

### 2. Filtering with bamFlag
You must define which alignments are "high quality" using the `bamFlag` list.
```r
# Example for high-quality paired-end reads
bamFlag <- list(
  isMinusStrand = FALSE, 
  isProperPair = TRUE, 
  isDuplicate = FALSE
)
```
### Tip for Thesis:
If you are forced to use BAM files, always specify in your **Methods** which BAM tag you used for cell identification (e.g., "Cell barcodes were extracted from the CB tag of the position-sorted BAM files"). This ensures your analysis can be perfectly replicated by others. 

HPC Note: Converting BAM files to Arrow format is much more CPU and I/O intensive than using fragment files. If using BAMs on the HPC, ensure you request sufficient wall-time for your job.


### 3.8.1 Anatomy of a scATAC-seq BAM Record

Before running ArchR on BAM files, you should inspect the first few lines using `samtools view` to identify the barcode tag.

**Example Record:**
`A00519... chr1 9997 ... CB:Z:GCGGGTTAGAACGTCG-1 ...`

* **Target Tag:** `CB:Z:GCGGGTTAGAACGTCG-1` 
* **ArchR Setting:** `bcTag = "CB"` 

#### Key Troubleshooting Parameters:
* **`bcTag`**: Tells ArchR which tag contains the cell barcode (e.g., "CB" or "CR").
* **`gsubExpression`**: If your barcodes have extra characters you don't want (like the `-1` at the end), you can use this parameter to "clean" the string as it is read.
* **`bamFlag`**: Filters the reads based on their bitwise flags (e.g., keeping only "properly paired" reads or removing duplicates).



> **Thesis Tip:** If your data comes from a non-standard pipeline (anything other than CellRanger), always verify the `bcTag` manually. Using the wrong tag will result in 0 cells passing the filter.

# 4. Doublet Inference

Doublets are a primary source of technical noise where one droplet contains multiple nuclei, leading to "average" cellular profiles. ArchR identifies and removes these computationally to prevent the formation of artifactual clusters.

## 4.1 The Mechanism of Doublet Identification

https://www.archrproject.com/bookdown/how-does-doublet-identification-work-in-archr.html

ArchR utilizes a "Synthetic Doublet" approach to identify artifacts without requiring external genotype data.

### The Core Process:


1. **Simulation:** ArchR creates **synthetic doublets** by mixing reads from thousands of combinations of your real cells. 
If Cell A has open peaks at Region 1, and Cell B has open peaks at Region 2, the synthetic doublet will have signal at both regions. This is done for thousands of combinations to ensure that every possible pairing of cell types (e.g., T-cell + B-cell, Monocyte + T-cell) is represented in the simulation.

\
2. **Projection:** These simulated cells are projected into the project's embedding (e.g., UMAP). So this is the same UMAP, where my real cells are located.

Predictable Location: Because a synthetic doublet is a 50/50 mix of two cells, the dimensionality reduction algorithm (UMAP) will naturally place that synthetic point in the "no man's land" between the two parent cell clusters.

\
3. **Neighbor Analysis:** Real cells that cluster closely with synthetic doublets are identified.  For every real cell, ArchR identifies its $k$ nearest neighbors (usually $k=10$ or $15$).
It counts how many of those neighbors are synthetic doublets vs. real cells.

\
4. **Scoring:** To make the results statistically robust, ArchR doesn't just do this once.

Iteration: By repeating this "synthesis and neighbor-check" thousands of times, ArchR builds a distribution.

Enrichment Score: A cell that is consistently surrounded by synthetic doublets across many iterations receives a high Doublet Enrichment Score. This indicates that the cell's accessibility profile is indistinguishable from a known hybrid simulation.

\
![alt text](image-1.png)

#### How to read a Doublet UMAP:
* **High Enrichment (Purple/Dark):** These cells are likely doublets. They often form "bridges" or "tendrils" connecting two distinct clusters.


\
![alt text](image-2.png)

#### Validation of Doublet Detection using Genotyped Cell Lines

* Colored Clusters: Each color represents a distinct, known cell line (e.g., HeLa, Jurkat, K562) identified by its unique genetic SNPs (genotype).

* Black Dots (Demuxlet Doublets): These are "true" doublets identified by demuxlet, a tool that uses genotype information to find droplets containing more than one cell type.

* The "Bridge" Phenomenon: The black doublet points are not randomly scattered; instead, they form clear bridges or "strings" that connect the distinct colored clusters.

\
Predictability: The plot demonstrates that doublets have a predictable mathematical profile—they always appear as a mathematical average of two parent cells, placing them between clusters in low-dimensional space.

Ground Truth: Because we can see the black dots exactly where the "bridges" are, we have high confidence that ArchR's simulation method (which creates synthetic bridges) is looking in the right place.

### 4.1.3 Quantifying Accuracy: ROC Curves and AUC

To prove the reliability of the doublet scoring, ArchR's performance was measured using an ROC curve against known "ground truth" doublets.

* **ArchR Prediction Accuracy:** Reaches an **AUC of ~0.92**. In statistics, an AUC above 0.9 is considered "excellent" for a classification task.
* **ArchR vs. Fragment Counts:** A common mistake is assuming that doublets are simply cells with the highest number of fragments. This plot shows that using fragment counts is a poor predictor (**AUC ~0.64**).
* **Conclusion:** ArchR’s simulation method effectively distinguishes between "deeply sequenced single cells" and "true doublets," which a simple fragment cutoff would fail to do.

\
![alt text](image-3.png)
* **True Positive Rate (Y-axis)**: The percentage of "real" doublets correctly identified.
* **False Positive Rate (X-axis)**: The percentage of single cells mistakenly flagged as doublets.
* **AUC (Area Under the Curve)**: A single number representing total accuracy. A score of 1.0 is a perfect test, while 0.5 (the dashed diagonal line) is no better than a random guess.

### 4.1.4 The "Cleaned" Dataset: Post-Doublet Removal

Following the computational removal of doublets, the data structure shifts from a connected "web" to distinct, isolated clusters.

* **Visual Impact:** The "bridges" of black dots that previously connected unrelated cell types have been eliminated.
* **Biological Accuracy:** The resulting UMAP now accurately reflects the expected biological composition: 10 distinct, non-overlapping cell lines.
* **Significance for Analysis:** Removing these "hybrid" cells is essential before performing clustering or trajectory analysis. Without this step, doublets could be misinterpreted as transitional cell states or rare progenitor populations.

* **Residual Black Dots**: 
These are likely homotypic doublets (e.g., two Jurkat cells in one droplet). Because they have the same profile as a single Jurkat cell, they do not form "bridges" and are nearly impossible to remove computationally—but they also do not cause "phantom clusters," so they are less harmful to your analysis.

\
![alt text](image-4.png)
> **Thesis Note:** Use this "Before and After" comparison to demonstrate the effectiveness of your QC pipeline. It proves that ArchR effectively identifies heterotypic doublets (those between different cell types) which are the most damaging to cluster identity.

## 4.2 Inferring scATAC-seq Doublets

After creating the `ArchRProject`, we must identify potential doublets. This process uses the "Synthetic Doublet" simulation to assign a probability score to every cell.

### 4.2.1 The Code
```r
proj <- addDoubletScores(
  input = proj,
  k = 10, 
  knnMethod = "UMAP", 
  LSIMethod = 1
)
```


#### Configuration Parameters
| Parameter | Setting | Description |
| :--- | :--- | :--- |
| `k` | 10 | Number of neighbors for score calculation. |
| `knnMethod` | "UMAP" | Embedding space used for neighbor search. |
| `LSIMethod` | 1 | LSI projection version. |

#### HPC Log Interpretation
* **Reading Fragments:** ArchR is pulling data from disk.
* **In Silico Synthesis:** Creating thousands of synthetic doublet profiles in RAM.
* **Success Message:** Metadata columns (`DoubletScore` and `DoubletEnrichment`) successfully added to Arrow files.

#### Validation Metrics
* **ArchR Accuracy:** AUC ~0.92 (Excellent).
* **Naive Frag Count Accuracy:** AUC ~0.64 (Poor).

> **Thesis Note:** Emphasize that ArchR's simulation approach outperforms simple fragment cutoffs by nearly 30% in accuracy (AUC 0.92 vs 0.64), ensuring that high-depth single cells are not accidentally discarded.

### 4.2.3 Troubleshooting: Low $R^2$ and Homotypic Doublets

During processing, ArchR reports an $R^2$ value for each sample's projection. This metric is a "health check" for the doublet simulation.

#### 1. Understanding $R^2$ Thresholds
* **$R^2 > 0.9$:** Reliable projection; the sample has enough heterogeneity to identify "bridge" doublets.
* **$R^2 < 0.9$:** Low heterogeneity; cells are too similar, meaning most doublets are **homotypic** (two similar cells in one droplet).

#### 2. Why ArchR Skips Samples
If $R^2$ is low, ArchR will skip doublet prediction because the results would likely be inaccurate. In a uniform population, a doublet's signal is indistinguishable from a single cell's signal.

#### 3. The LSI Workaround
If you must run doublet detection on a low-heterogeneity sample, you can force the calculation in the LSI subspace:

```r
proj <- addDoubletScores(
  input = proj,
  k = 10, 
  knnMethod = "LSI", 
  force = TRUE 
)
```

### 4.2.4 Interpreting Doublet QC Plots

After running `addDoubletScores()`, ArchR populates the `QualityControl` folder with three key visualizations per sample. These plots confirm whether the "Synthetic Doublet" bridges align with the intermediate cells in your actual data.

| Plot Type | Metric | Primary Use |
| :--- | :--- | :--- |
| **Doublet Enrichment** | Relative density vs. expected | **Primary metric** used for doublet identification and filtering. |
| **Doublet Scores** | $-\log_{10}(\text{p-adj})$ | Statistical significance; used as a secondary validation. |
| **Doublet Density** | Projection density | Visualizes where the "synthetic" artifacts are located. |

#### Representative Visualizations (BMMC Sample)

**A. Doublet Enrichment**
![alt text](image-5.png)
*High enrichment (darker colors) indicates cells that likely represent a mixture of two distinct cell types.*

**B. Doublet Scores**
![alt text](image-6.png)
*Shows the statistical probability of a cell being a doublet; used to ensure filtering is not arbitrary.*

**C. Doublet Density**
![alt text](image-7.png)
*This "underlay" plot shows where the simulation placed the synthetic doublets. Note how the density is highest in the "bridges" between biological clusters.*

> **Thesis Tip:** Always check the **Doublet Density** plot first. If the synthetic density is overlapping heavily with your main clusters (rather than the bridges), it may indicate that your sample has low heterogeneity, and you should be careful not to over-filter your data.

## 4.3 Using demuxlet with ArchR (Genotype Validation)

If your dataset contains a pool of multiple donors, you can use **demuxlet** to provide a "Ground Truth" for doublet detection.

### 4.3.1 What is demuxlet?
Demuxlet is a tool that utilizes genetic polymorphisms (SNPs) to assign each barcode to a specific donor. It identifies:
* **Singlets:** Barcodes matching a single donor.
* **Doublets:** Barcodes containing a mixture of SNPs from two different donors.
* **Ambiguous:** Barcodes with insufficient SNP coverage.

### 4.3.2 Why integrate with ArchR?
Integrating `demuxlet` allows you to validate ArchR's computational `DoubletEnrichment` scores. In the ArchR "Mixology" study, computational scores matched demuxlet's physical assignments with an **AUC > 0.90**.

### 4.3.3 Integration Workflow
1. **Run demuxlet:** Perform genotype-based demultiplexing on your alignment (BAM) files outside of ArchR.
2. **Import results:** Use `addCellColData()` to add the assignments to your `ArchRProject`.
3. **Compare:** Plot the `DoubletEnrichment` scores and color them by the `Demuxlet` assignment to see if the "bridges" align with the known doublets.

| Feature | ArchR Simulation | demuxlet |
| :--- | :--- | :--- |
| **Requirements** | scATAC-seq data only | Genotype (VCF) files + Mixed Donors |
| **Primary Use** | Standard QC for any sample | Gold standard validation |
| **Detects...** | Heterotypic doublets (bridges) | Doublets between different genotypes |

> **Thesis Tip:** If your research uses donor pooling, showing a correlation between ArchR’s Doublet Scores and demuxlet’s assignments is a powerful way to prove your data's technical rigor.



# 5. Creating an ArchRProject

The `ArchRProject` is the central object used for all downstream analysis. It links multiple Arrow files into a single project, allowing for cross-sample comparisons and unified clustering.

## 5.0 Overview
The `ArchRProject` is unique because it is a **small, memory-resident object** that manages **large, disk-resident data** (Arrow files). By interacting with this object, you can rapidly perform complex calculations without crashing your R session due to memory limits.

### Why use an ArchRProject?
1. **Unified Workflow:** It is the foundation for almost every ArchR function.
2. **Reproducibility:** Projects can be saved and re-loaded later to maintain analysis continuity.
3. **Collaboration:** Makes it easy to zip and share an entire project with a collaborator.

## 5.1 Initialization

```r
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "HPC_ArchR_Project",
  copyArrows = TRUE
)
```

#### Parameter Deep Dive

**ArrowFiles**: A character vector of paths to the .arrow files you want to include in the analysis.

**outputDirectory**: The name of the folder where ArchR will save your project object, plots, and intermediate data.

**copyArrows**: If TRUE (Recommended), ArchR creates a copy of your Arrow files inside the outputDirectory. This makes the project self-contained.

If FALSE, it only creates links to the original files. If you move or delete those original files, your project will break.


We can examine the contents of our ArchRProject using `proj(name)`


#### What the Project Object Contains
Once created, the `proj object` is a lightweight container for the following:

**Cell Metadata** (cellColData): Stores information for every cell (e.g., sample name, TSS scores, total fragments).

**Sample Statistics**: Summary metrics for each sample.

**Analysis State**: Tracks which matrices have been added (Tile, GeneScore, Peak) and what clustering/UMAP parameters were used.

### Saving an ArchRProject

```r
### Corrected save function with quotes and proper parameter naming
saveArchRProject(
    ArchRProj = projHeme1,
    outputDirectory = "/usr/people/BZEDVZ/18perkov/ArchR_Tutorial", # Must be in quotes!
    overwrite = TRUE,
    load = TRUE,
    dropCells = FALSE,
    logFile = createLogFile("saveArchRProject"),
    threads = 16 # Usually set as an integer
)
```
__Very important!__ The load parameter determines whether or not the saveArchRProject() function will return the saved ArchRProject object which you would assign to overwrite the original ArchRProject object or provide a new ArchRProject name using <-. This effectively saves and loads the ArchRProject from its new location. If load = FALSE, then this process does NOT update the ArchRProject object that is active in your current R session. Specifically, the object named projHeme1 in the current R session will still point to the original location of the Arrow files, not the copied Arrow files that reside in the specified outputDirectory. You might use this behavior if you were saving your ArchRProject and shutting down your R session so that the project can be reloaded at a later time.

__Also important!__ The term ArchRProject can be confusing to some users because we use this both to refer to the object that is actively loaded within the R environment and to the file that has been saved on disk as part of the saveArchRProject() process. However, it is extremely important to understand that these are not equivalent. Manipulations that you perform on an ArchRProject that is actively loaded within the R environment do not automatically propagate to the on-disk file. You must run the saveArchRProject() function to store those changes.
### Loading an existing ArchRProject
Once a project has been saved to the HPC disk, you do not need to re-run the Arrow file creation or initial QC. You can simply "point" R to the directory where the project lives.

#### The Implementation
To load your project, provide the path to the **folder** containing your data. You do not need to point to a specific file; ArchR will automatically look for the `.rds` state file within that folder.

```r
### Load the project from your specific HPC directory
projHeme1 <- loadArchRProject(path = "/usr/people/BZEDVZ/18perkov/ArchR_Tutorial")

### Verify the project loaded correctly
projHeme1
### Loading an ArchRProject
```


## 5.2 Inspecting Cell Metadata (cellColData)

The `cellColData` stores all metadata associated with individual cells. This includes the QC metrics generated during Arrow file creation and the scores from doublet inference.

### 5.2.1 How to Access Metadata in R
You can interact with this data like a standard R DataFrame or by using the `$` accessor on the project object:

```r
# View the first few rows
head(proj@cellColData)

# Access a specific metric for all cells
tss_scores <- proj$TSSEnrichment

# Summary of fragment counts
summary(proj$nFrags)
```



### 5.2.4 Detailed cellColData Definitions

The `cellColData` object stores critical metadata for every individual cell in the project. The following table defines the columns calculated during Arrow file creation:

| Column Name | Description |
| :--- | :--- |
| **TSSEnrichment** | The per-cell Transcription Start Site (TSS) enrichment score. |
| **ReadsInTSS** | The number of reads that fall within TSS regions (default is 100 bp around TSS). |
| **ReadsInPromoter** | The number of reads that fall in promoter regions (default is -2000 to +100 from the TSS). |
| **PromoterRatio** | The ratio of reads in promoters to reads outside of promoters. |
| **ReadsInBlacklist** | The number of reads that fall in defined genomic blacklist regions. |
| **BlacklistRatio** | The ratio of reads in blacklist regions to reads outside of blacklist regions. |
| **NucleosomeRatio** | Represents the ratio of reads mapping to nucleosome-sized fragments, calculated as: $(nDiFrags + nMultiFrags) / nMonoFrags$. |
| **nFrags** | The total number of unique nuclear fragments recovered per cell. |
| **nMonoFrags** | The number of fragments with a length less than $2 \times \text{nucLength}$ (where `nucLength` is 147 bp by default). |
| **nDiFrags** | The number of fragments with a length $\geq 2 \times \text{nucLength}$ but $< 3 \times \text{nucLength}$. |
| **nMultiFrags** | The number of fragments with a length $\geq 3 \times \text{nucLength}$. |
| **PassQC** | Equal to $1$ if the cell passed initial QC filters or $0$ if it did not. |

---

### Technical Reference
In ArchR, fragment sizes are categorized based on **nucLength**, which is the length of DNA wrapped around a single nucleosome (147 bp by default). This allows the calculation of the **NucleosomeRatio** to assess chromatin state:

$$\text{NucleosomeRatio} = \frac{nDiFrags + nMultiFrags}{nMonoFrags}$$

### Deep Dive: Understanding cellColData Metrics

The `cellColData` object is the "biopsy report" for every individual cell in your project. These metrics, calculated during Arrow file creation, allow us to distinguish high-quality biological signals from technical noise.

---

### 1. Signal-to-Noise & Purity Metrics
These metrics measure how well the Tn5 enzyme targeted open, active regions of the genome compared to random background noise.

* **TSSEnrichment**
    * **Definition:** The ratio of fragments centered at Transcription Start Sites (TSS) vs. the surrounding genomic background.
    * **Biological Context:** Active promoters are depleted of nucleosomes and are "wide open." A high score (typically > 4) indicates a high-quality cell with a clear signal.
* **ReadsInTSS**
    * **Definition:** The raw number of reads falling within a 100 bp window around known TSS coordinates.
* **ReadsInPromoter & PromoterRatio**
    * **Definition:** Measures fragments in the broader promoter neighborhood (default: -2000 to +100 bp from TSS).
    * **Importance:** A high **PromoterRatio** indicates an efficient library prep where sequencing depth was spent on regulatory regions rather than "genomic dark matter."



---

### 2. Artifact & Background Metrics
These metrics identify fragments coming from "sticky" or problematic areas of the genome.

* **ReadsInBlacklist & BlacklistRatio**
    * **Definition:** The count and proportion of fragments falling into the **ENCODE Blacklist** (genomic regions that consistently show high, non-specific background noise).
    * **Importance:** High ratios suggest the cell may have been dying or the DNA was degraded, leading to "sticky" non-specific noise that can interfere with clustering.

---

### 3. Structural & Fragment Size Metrics
ATAC-seq leverages the fact that DNA wraps around histones (nucleosomes) in 147 bp increments. Tn5 cuts only in the "linker" DNA between these beads, creating a characteristic "staircase" pattern.



* **nFrags**
    * **Definition:** The total unique nuclear fragments recovered. This is your primary measure of **sequencing depth** per cell.
* **nMonoFrags**
    * **Definition:** Fragments shorter than two nucleosomes (~147-294 bp). These are usually the most informative "open chromatin" fragments.
* **nDiFrags & nMultiFrags**
    * **Definition:** Fragments spanning two or three+ nucleosomes. These represent more compact or "closed" chromatin states.
* **NucleosomeRatio**
    * **Formula:** 
  $$\text{NucleosomeRatio} = \frac{nDiFrags + nMultiFrags}{nMonoFrags}$$
    * **Importance:** A high ratio indicates poor Tn5 penetration or DNA clumping, while a low ratio indicates high-quality, accessible chromatin.



---

### 4. The Analysis Gatekeeper

* **PassQC**
    * **Definition:** A binary flag (**1** for Pass, **0** for Fail).
    * **Importance:** ArchR labels low-quality cells but doesn't delete them immediately. When you run clustering or UMAP, ArchR automatically filters the project to include only cells where `PassQC = 1`.

---

## Technical Summary Table

| Metric Category        | Key Metric        | High Value Interpretation                     |
| :--------------------- | :---------------- | :-------------------------------------------- |
| **Purity**             | `TSSEnrichment`   | Strong biological signal; real cell.          |
| **Efficiency**         | `PromoterRatio`   | Targeted sequencing; high-quality library.    |
| **Technical Noise**    | `BlacklistRatio`  | Potential artifact or dying cell.             |
| **Physical Integrity** | `NucleosomeRatio` | Poor enzyme access or DNA degradation.        |
| **Depth**              | `nFrags`          | Highly sequenced cell with rich data density. |

> **Note:** ArchR uses a default **nucLength** of **147 bp** (the length of DNA wrapped around a human nucleosome) to categorize all fragment sizes.

## 5.3 Manipulating an ArchR Project

### 5.3.1 Interacting with Metadata using the `$` Accessor

The `$` operator is a convenient shorthand in ArchR that allows you to bypass the `@cellColData` slot and access metadata columns directly from the project object. This is identical to how you would access columns in a standard R data frame.

#### Example 1: Accessing Cell Names
Every cell in ArchR is uniquely identified by its sample name and barcode. You can retrieve these using `$cellNames`:

```r
head(proj$cellNames)
``` 

We can access the sample names associated with each cell:

```r
head(proj$Sample)
```
We can access the TSS Enrichment Scores for each cell:
```r
quantile(proj$TSSEnrichment)
##       0%      25%      50%      75%     100% 
##  4.10900 13.92550 16.81500 19.93025 41.98000
``` 


When we run `quantile(proj$TSSEnrichment)`, we are looking at the statistical "spread" of our data quality. This tells us how consistent the signal is across our entire cell population.

| Quantile         | Score      | Interpretation                                                                                                                                          |
| :--------------- | :--------- | :------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **0% (Min)**     | **4.109**  | **Filter Check:** Confirms every cell in the project successfully passed the minimum QC threshold (TSS > 4).                                            |
| **25%**          | **13.925** | **Lower Quartile:** Even the lower-performing cells in this dataset have excellent signal-to-noise ratios.                                              |
| **50% (Median)** | **16.815** | **Core Quality:** This is your primary "Health Metric." A median of ~17 indicates an exceptionally clean library with very high promoter accessibility. |
| **75%**          | **19.930** | **Upper Quartile:** 25% of your cells have a signal-to-noise ratio of 20:1 or better.                                                                   |
| **100% (Max)**   | **41.980** | **Upper Limit:** The highest quality single-cell capture in the sample.                                                                                 |

---

#### Why these numbers matter for the Thesis:
In scATAC-seq, the median TSS enrichment score is the most cited metric for data quality.
* **TSS < 4:** Generally considered "failing" or background noise.
* **TSS 4–10:** Standard/Acceptable quality.
* **TSS > 10:** High-quality data.
* **TSS ~17 (Our Result):** Represents state-of-the-art library prep with very high biological signal.
  
### 5.3.2 Example 2: Subsetting an ArchRProject by Cells

Subsetting is a powerful feature that allows you to create a new, smaller `ArchRProject` containing only a specific group of cells. This is essential if you want to analyze a single sample in isolation or remove specific clusters.

#### 1. The Subsetting Syntax
The `ArchRProject` behaves like a standard R matrix. You can subset it using the bracket notation: `proj[rows, columns]`.

* **Rows:** Represent individual cells.
* **Columns:** Represent different data types (rarely used in basic subsetting).

#### 2. Practical Examples

**A. Subsetting by Cell Names**
If you have a specific list of cell names (barcodes), you can filter the project to keep only those cells:

```r
# Create a list of the first 100 cell names
cells_to_keep <- head(proj$cellNames, 100)

# Subset the project
subset_proj <- proj[cells_to_keep, ]
```

**B. Subsetting by Metadata Condition**

Keep only cells where the 'Sample' metadata matches BMMC
```r
projBMMC <- proj[proj$Sample == "scATAC_BMMC_R1", ]
``` 

#### 3. What happens to the data?

When you subset an ArchRProject, ArchR does not create new Arrow files on your disk. Instead, it creates a new, lightweight project object in memory that "points" only to the subset of cells in the original Arrow files.

Why this is useful for your Thesis:
Speed: You can quickly test a hypothesis on a single sample without re-processing the entire dataset.

Clarity: You can remove "junk" clusters or doublets from your project to make your final UMAP plots cleaner.

Comparison: You can create two separate projects (e.g., proj_Control and proj_Treated) to run specific differential accessibility tests.

Technical Note: After subsetting, your cell counts and metadata distributions (like those seen in quantile()) will automatically update to reflect only the remaining cells.

It is important to note that some operations will need to be re-run after project subsetting. Much of this is dictated by common sense. For example, if you want to look at only the subsetted cells on a new UMAP embedding, you need to create that new embedding first. Similarly, if you want to perform sub-clustering, you should almost certainly re-run dimensionality reduction (LSI) first.

The primary disadvantage of subsetArchRProject() is that it makes copies of the Arrow files which can be quite large for bigger data sets. Nevertheless, this is the absolute most stable way to subset a project and is the only way that we recommend.

### Physical Subsetting with `subsetArchRProject()`
---

#### R Implementation Example
In this example, we identify all cells belonging to the "BMMC" sample and save them as a separate project.

```r
# Step 1: Identify the cell indices for the BMMC sample
idxSample <- BiocGenerics::which(projHeme1$Sample %in% "scATAC_BMMC_R1")

# Step 2: Create the physical subset
projSubset <- subsetArchRProject(
  ArchRProj = projHeme1,               # The original project
  cells = projHeme1$cellNames[idxSample], # Barcodes to keep
  outputDirectory = "BMMC_Only_Subset",   # New directory on the HPC
  dropCells = TRUE,                    # Physically remove non-selected cells
  force = TRUE                         # Overwrite if the folder exists
)

``` 
For further options, visit: https://www.archrproject.com/bookdown/manipulating-an-archrproject.html




### Example 3: Adding Data to an ArchRProject

**1) Creating Custom Metadata**

Often, the original sample names (e.g., scATAC_BMMC_R1) are too long for plotting. You can use R functions like `gsub()` to clean them up:
```r
### Create "bioNames" by removing "scATAC_" and "_R1"
bioNames <- gsub("_R2|_R1|scATAC_","", proj$Sample)

### Look at the first few cleaned names
head(bioNames)
### Output: [1] "BMMC" "BMMC" "BMMC" ...
```
**2) Method A: The $ Accessor (Quick & Dirty)**

If you have a vector that is exactly the same length as the number of cells in your project, you can assign it directly:
```r
proj$bioNames <- bioNames
``` 

**3) Method B: `addCellColData()` (Precise & Robust):**
   
This function is more powerful because it allows you to add data to only a subset of cells. ArchR will automatically fill in the missing entries with NA.
```r
# Example: Adding data only to the first 10 cells
proj <- addCellColData(
  ArchRProj = proj, 
  data = bioNames[1:10], 
  cells = proj$cellNames[1:10], 
  name = "bioNames_Subset"
)
```

**4) Verifying the Addition:**
You can use `getCellColData()` to retrieve specific columns and compare them:
```r
# Retrieve both the full and subsetted columns
metadata_check <- getCellColData(proj, select = c("bioNames", "bioNames_Subset"))
head(metadata_check)
```
**Thesis Tip:** Always create a "Clean_Sample_Name" column early in your workflow. Using short names like "Control" and "Mutant" instead of long file names will make your downstream UMAP legends and heatmaps much more professional and readable.

### Example 4: Obtaining Columns from cellColData

While the `$` accessor is convenient for quick checks, the `getCellColData()` function is the specialized tool for retrieving metadata. It is more flexible because it allows you to select multiple columns at once, perform mathematical operations during retrieval, and always returns a clean `S4 Vectors DataFrame`.

**1) Retrieving a Single Column by Name**

You can easily pull out specific metrics, such as the total unique fragments (`nFrags`), to use in external plotting or statistical tests:
```r
### Retrieve a single column
df <- getCellColData(projHeme1, select = "nFrags")

### View the resulting DataFrame
head(df)
### Output: DataFrame with barcodes as row names and "nFrags" as the column.
```
**2) Performing Operations "On the Fly"**

One of the most powerful features of `getCellColData()` is the ability to perform calculations directly within the `select` parameter. This avoids creating unnecessary permanent columns in your project:
```r
### Perform log10 transformation and basic arithmetic during retrieval
df <- getCellColData(projHeme1, select = c("log10(nFrags)", "nFrags - 1"))

### View the transformed data
head(df)
### Output: A DataFrame with two columns: "log10(nFrags)" and "nFrags - 1"
```

**3) Selecting Multiple Metrics**

If you need to compare two different quality control metrics (like Sequencing Depth vs. TSS Enrichment), you can request them together in a single vector:
```r
### Retrieve multiple QC metrics simultaneously
qc_metrics <- getCellColData(projHeme1, select = c("nFrags", "TSSEnrichment"))

### Check the summary stats for the retrieved data
summary(qc_metrics)
```

**4) Why use `getCellColData()` instead of `$`?**
   
Calculations: It supports inline math (e.g., log10, sqrt, addition/subtraction).

Format: It consistently returns a DataFrame object which preserves cell barcodes as row names, ensuring your data never gets un-synced.

Multi-selection: It can handle a list of many columns at once, which is much cleaner than calling proj$Column five separate times.

**Thesis Tip**: When reporting Sequencing Depth in your thesis, it is standard practice to use the log10-transformed fragment count. Instead of permanently cluttering your project with a "log10_nFrags" column, use getCellColData(proj, select = "log10(nFrags)") to generate your plots and tables on the fly. This keeps your ArchRProject file size smaller and your workspace cleaner.

### Example 5: Plotting QC metrics - log10(Unique Fragments) vs TSS enrichment score

**Selecting the Two Most Robust Metrics** 

While `cellColData` contains 15+ columns, these two are considered the most reliable indicators of technical success:

  * **log10(nFrags)**: Represents the number of unique nuclear fragments (sequencing depth). We use the log10 transformation because fragment counts can vary from 1,000 to 100,000; the log scale makes this range linear and easier to visualize.

  * __TSSEnrichment:__ Measures the signal-to-background ratio. Cells with high signal at Transcription Start Sites are "true" cells, while low scores represent random genomic noise.

\
We use getCellColData() to pull both metrics at once. Notice that we perform the log10 operation directly inside the function call:
```r
### Extract depth and signal quality simultaneously
df <- getCellColData(projHeme1, select = c("log10(nFrags)", "TSSEnrichment"))

### View the first few rows of your new QC dataframe
head(df)
```
__What does the Output tell us?__

The resulting df (DataFrame) links every cell barcode to its specific QC values. This is the raw data used to generate the "cloud" plot seen in Section 5.2.8.

  * __High log10(nFrags) + High TSSEnrichment:__ These are your "Healthy Cells" (Top-Right of the cloud).

  * __Low log10(nFrags):__ These are "Under-sequenced" cells that lack enough data to be analyzed reliably.

  * __Low TSSEnrichment:__ These are "Noisy" cells (background noise or dying cells) that will likely fail to cluster correctly.

  * __Why we do this before Clustering:__ By looking at these metrics together, we can verify if the QC Cutoffs we set during Arrow file creation (e.g., filterTSS = 4 and filterFrags = 1000) were appropriate. If the "cloud" of cells is too close to the dashed lines, it may indicate that we need to be more stringent with our filters to ensure high-quality downstream results.

### Plotting
```r
p <- ggPoint(
    x = df[,1], # Log10 Unique Fragments
    y = df[,2], # TSS Enrichment
    colorDensity = TRUE, # Colors points by how crowded the area is
    continuousSet = "sambaNight",  # High-contrast color palette
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

p
``` 

#### Parameter Breakdown: Customizing the QC Plot

The `ggPoint()` function is a highly optimized plotting utility. Understanding its parameters allows you to fine-tune your Quality Control visualizations to ensure they are both mathematically accurate and publication-ready.

#### Detailed Parameter Reference

| Parameter           | Value / Function     | Purpose                                                                                                                                                                       |
| :------------------ | :------------------- | :---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **`x`**             | `df[,1]`             | Maps the **Log10(Unique Fragments)** to the horizontal axis.                                                                                                                  |
| **`y`**             | `df[,2]`             | Maps the **TSS Enrichment Score** to the vertical axis.                                                                                                                       |
| **`colorDensity`**  | `TRUE`               | Calculates local point density. This highlights the "core" of your cell population by coloring the most crowded areas differently.                                            |
| **`continuousSet`** | `"sambaNight"`       | Specifies the color palette. "sambaNight" is a high-contrast theme (dark blue to bright yellow) that makes density variations easy to see.                                    |
| **`xlim`**          | `c(log10(500), ...)` | Sets the x-axis range. Starting at 500 fragments provides visual context for the 1,000 fragment cutoff.                                                                       |
| **`ylim`**          | `c(0, ...)`          | Sets the y-axis range. Starting at 0 allows you to see the full "floor" of the signal-to-noise ratio.                                                                         |
| **`quantile`**      | `probs = 0.99`       | **Outlier Control:** By setting the upper limit to the 99th percentile, ArchR prevents a few extreme "super-cells" from compressing the rest of your data into a tiny corner. |

---

#### The Significance of the "Dashed Lines"

To finalize the plot, we add standard `ggplot2` layers to represent our biological and technical filters:

1. **`geom_hline(yintercept = 4)`**:
   * **The TSS Filter:** This horizontal dashed line represents the minimum acceptable signal-to-background ratio. Cells below this line are likely background noise (genomic "soup").
2. **`geom_vline(xintercept = 3)`**:
   * **The Depth Filter:** Since we are on a log10 scale, $10^3 = 1,000$. This vertical line marks the cutoff for sequencing depth. Cells to the left of this line do not have enough data for confident analysis.

---

### Why the 99th Percentile Matters
In almost every scATAC-seq run, you will have a few "outlier" droplets with massive fragment counts (e.g., 500,000+ fragments). If your plot scales to accommodate those few dots, your actual cell population (usually between 2,000 and 10,000 fragments) will appear as a microscopic smudge on the far left. 

Using `quantile(df[,1], probs = 0.99)` tells ArchR: *"Ignore the top 1% of extreme values when deciding how wide to make the plot."* This ensures the "orange core" of your data is always front and center.

* **Thesis Tip:** If you notice that your "density core" (the brightest yellow/orange part) is very close to or overlapping with the dashed lines, it is a sign that your filters are too lenient. You may need to increase your TSS threshold to 6 or 8 for your final analysis to ensure you are only looking at the highest-quality nuclei.

![alt text](image-9.png)

To save an editable vectorized version of this plot, we use `plotPDF()`. This saves the plot within the “Plots” sub-directory of our ArchRProject directory (defined by `getOutputDirectory(projHeme1)`).

```r
plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = projHeme1, addDOC = FALSE)
## Plotting Ggplot!
```

# 5.3 Plotting Sample Statistics

Once multiple samples are integrated into an `ArchRProject`, it is vital to compare their quality metrics side-by-side. This ensures that downstream results (like clustering) are driven by biology rather than technical differences between samples.

### 5.3.1 Introduction to `plotGroups()`
The `plotGroups()` function is the universal tool in ArchR for visualizing distributions across categories. It primarily produces:
* **Ridge Plots (`ridges`)**: Ideal for comparing the "shape" of distributions across many groups.
* **Violin Plots (`violin`)**: Ideal for seeing the density and specific percentiles (when paired with boxplots).

---

### 5.3.2 Example 1: Ridge Plots for TSS Enrichment

In this example, we visualize the **TSS Enrichment** distribution for every sample in the project. This allows us to see if the "signal-to-noise" is consistent across our experiment.

#### Code Breakdown

```r
p1 <- plotGroups(
    ArchRProj = proj,         # Your project object
    groupBy = "Sample",       # The metadata column to group by
    colorBy = "cellColData",  # Tells ArchR the data is in the metadata table
    name = "TSSEnrichment",   # The specific metric to plot
    plotAs = "ridges",        # The type of plot
    baseSize = 10             # Adjusts the global font size for the plot
)
```
The `plotGroups()` function is highly flexible. Below is a breakdown of the key parameters used to customize your sample statistic plots:

| Parameter        | Definition                         | Purpose                                                                                                                 |
| :--------------- | :--------------------------------- | :---------------------------------------------------------------------------------------------------------------------- |
| **`ArchRProj`**  | The ArchRProject object.           | Tells the function which project to pull metadata and data from.                                                        |
| **`groupBy`**    | A string (e.g., `"Sample"`)        | Defines how to categorize the cells on the Y-axis. You can group by samples, clusters, or any custom metadata column.   |
| **`colorBy`**    | A string (usually `"cellColData"`) | Specifies the "slot" in the ArchRProject where the data is stored. For QC metrics, this is always `cellColData`.        |
| **`name`**       | A string (e.g., `"TSSEnrichment"`) | The exact name of the column you want to plot. This can also include math like `"log10(nFrags)"`.                       |
| **`plotAs`**     | `"ridges"` or `"violin"`           | Sets the visual style. **Ridges** show overlapping density curves; **Violin** shows the density "envelope" of the data. |
| **`baseSize`**   | Numeric (default is often 10-12)   | A global multiplier for text size. Increasing this makes axis labels and titles larger for presentations or posters.    |
| **`alpha`**      | Numeric (0 to 1)                   | Sets the transparency of the plot. Especially useful for violin plots to see the internal boxplots more clearly.        |
| **`addBoxPlot`** | Logical (`TRUE`/`FALSE`)           | (Violin only) If `TRUE`, adds a box-and-whisker plot inside the violin to show the median and quartiles.                |

---

### Understanding the Logic
When you call `plotGroups()`, you are essentially telling ArchR:
> "Go into my **`ArchRProj`**, look into the **`colorBy`** table, find the column named **`name`**, and group that data by the categories found in **`groupBy`**. Finally, render the results as **`plotAs`**."

* **Thesis Tip:** Consistency is key for your final figures. Once you find a `baseSize` and `alpha` value that looks good, use those same values for every ridge and violin plot in your thesis to give your document a cohesive, professional look.

### Saving Plots with `plotPDF()`

Unlike standard R saving functions, `plotPDF()` is "Project Aware." It automatically handles the file pathing to keep your HPC workspace organized.

#### 1. How the Path is Determined
When you run `plotPDF()`, ArchR looks at the `outputDirectory` stored inside your `ArchRProject` object. It then:
1. Navigates to that directory.
2. Creates a sub-folder called **"Plots"** (if it doesn't already exist).
3. Saves your PDF there using the provided `name`.

```r
# Save multiple plots (p1, p2, p3, p4) into one vectorized PDF
plotPDF(
    p1, p2, p3, p4, 
    name = "QC-Sample-Statistics.pdf", 
    ArchRProj = projHeme1, 
    addDOC = FALSE, 
    width = 4, 
    height = 4
)  
```

## 5.4 Plotting Sample Fragment Size Distributions

The fragment size distribution is a direct reflection of the "beads on a string" structure of your chromatin. High-quality scATAC-seq data should always show a clear periodicity (peaks and valleys) corresponding to nucleosome spacing. Fragment size distributions in ATAC-seq can be quite variable across samples, cell types, and batches. Slight differences like those shown below are common and do not necessarily correlate with differences in data quality.

#### 1. The `plotFragmentSizes()` Function
ArchR accesses the pre-computed fragment lengths in the Arrow files to visualize these distributions across all samples simultaneously.

```r
# Generate the fragment size distribution plot
p1 <- plotFragmentSizes(ArchRProj = projHeme1)

# Display the plot
p1
```

#### 2. Parameter Definitions

`ArchRProj` = The project object containing the links to your Arrow files.

`maxSize`	(Optional) = The maximum fragment size to plot. Default is usually 750 bp.

`groupBy`	(Optional) = Allows you to group the distribution by Sample (default) or Cluster.

![alt text](image-11.png)

#### 3. How to Interpret the Plot

  Sub-nucleosomal Peak: The largest peak on the left (<150 bp). Represents open, accessible regions.

  Nucleosomal Periodicity: You should see smaller, diminishing peaks at ~200 bp, ~400 bp, and ~600 bp.

  Quality Check: If your plot is just a flat line or a single smooth "hump" without distinct peaks/valleys, your library likely has poor chromatin structure (potentially due to over-fragmentation or cell death).

    Thesis Tip: Include this plot to demonstrate that your library has the expected nucleosomal periodicity. Even if the heights of the peaks vary slightly between your experimental groups, the presence of the "staircase" proves that the Tn5 enzyme was acting on healthy, protein-bound chromatin rather than naked, degraded DNA.

### 5.4.2 Plotting TSS Enrichment Profiles

The TSS enrichment profile provides a visual "fingerprint" of your library's signal-to-noise ratio. A high-quality profile should show a sharp, symmetric peak at the center of the TSS.

#### 1. The `plotTSSEnrichment()` Function
This function aggregates the signal across thousands of known transcription start sites to show the average accessibility "landscape" of your samples.

```r
# Generate the TSS enrichment profile plot
p2 <- plotTSSEnrichment(ArchRProj = projHeme1)

# Display the plot
p2
``` 
![alt text](image-12.png)

#### 2. Visual Breakdown of the Plot:

* **The Center Peak (Distance = 0):** This represents the **Nucleosome-Free Region (NFR)**. High-quality libraries show a sharp, tall spike here, indicating that the Tn5 enzyme successfully targeted open promoters.
* **The +1 Nucleosome Shoulder:** The distinct "bump" observed to the right of the central peak (roughly at +150 to +200 bp) is caused by the **+1 nucleosome**. The presence of this shoulder is a strong indicator of biological integrity and high-resolution data.
* **The Enrichment Score:** The Y-axis represents the fold-enrichment over the background. In this plot, the PBMC sample (green) shows a maximum enrichment of >20, which is considered exceptional.

#### Comparison Across Samples:
In the provided plot, we see three distinct lines:
1.  **scATAC_PBMC_R1 (Green)**: Highest signal-to-noise.
2.  **scATAC_CD34_BMMC_R1 (Blue)**: Intermediate signal.
3.  **scATAC_BMMC_R1 (Red)**: Lowest signal of the three, though still far above the standard QC threshold of 4.

> **Thesis Tip:** When discussing this plot, emphasize the **symmetry** and the **+1 nucleosome shoulder**. A messy or "jagged" plot usually indicates low cell numbers or high background noise. The smooth, distinct peaks shown here justify the use of these samples for high-resolution regulatory analysis.

# 5.6 Filtering Doublets from an ArchRProject

Doublets occur when two nuclei are captured in a single droplet. In scATAC-seq, these artifacts can lead to "ghost clusters" or false transitions in trajectory analysis. `filterDoublets()` removes these cells based on the scores calculated during the initial QC phase.

More on this topic here: https://www.archrproject.com/bookdown/filtering-doublets-from-an-archrproject.html


# 6. Dimensionality reduction with ArchR


![alt text](image-14.png)

![alt text](image-15.png)

## 6.2. Iterative Latent Semantic Indexing (LSI)


In single-cell analysis, "Dimensionality Reduction" is the process of compressing thousands of genomic features into a few coordinates (like a UMAP) that we can actually visualize and cluster.

### 1. The Comparison: scRNA-seq vs. scATAC-seq

| Feature                 | scRNA-seq (Gene Expression)                                                  | scATAC-seq (Chromatin)                                                                                                          |
| :---------------------- | :--------------------------------------------------------------------------- | :------------------------------------------------------------------------------------------------------------------------------ |
| **Data Type**           | **Continuous/Counts:** You can have 0, 1, 10, or 100 transcripts.            | **Binary:** A region is either open (1) or closed (0).                                                                          |
| **Selection Strategy**  | **Highly Variable Genes (HVGs):** Easy to find genes that vary across cells. | **Highly Accessible Peaks:** Selecting the "most open" regions fails.                                                           |
| **The "Noise" Problem** | Noise is usually technical (dropout).                                        | "Most accessible" regions are often "housekeeping" sites (open in all cells), which adds noise and masks cell-type differences. |

---

## 2. The Solution: The "Iterative LSI" Approach
Because we cannot identify "variable" peaks in a binary matrix effectively, ArchR uses a multi-step refinement process.

![alt text](image-16.png)


### Phase 1: The "Rough Sketch"
ArchR starts by looking at the **most accessible tiles** (usually 500bp windows). 
* **Action:** It runs an initial Latent Semantic Indexing (LSI) transformation.
* **Goal:** To get a "low-resolution" look at the data. 
* **Result:** This identifies major cell lineages (e.g., separating all T-cells from all B-cells) without getting bogged down in technical "batch" differences.

### Phase 2: Feature Refinement
Once these "rough" clusters are identified, ArchR calculates the average accessibility for each group.
* **Action:** It asks, *"Which peaks are actually different between Cluster A and Cluster B?"*
* **Goal:** To find the **Variable Peaks**.
* **Result:** This filters out the "housekeeping" noise and focuses on biologically informative regions.

### Phase 3: The "High-Definition" Map
ArchR runs LSI **a second time** (or more), but this time it only uses the variable features found in Phase 2.
* **Result:** A much cleaner dimensionality reduction that minimizes batch effects and provides a more accurate biological "fingerprint" of each cell.

---

## 3. Implementation in ArchR
The `addIterativeLSI()` function automates this entire process.

```R
projHeme2 <- addIterativeLSI(
    ArchRProj = projHeme2,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30
)
```

## Parameter Breakdown: `addIterativeLSI()`

The `addIterativeLSI()` function is the primary tool in ArchR for dimensionality reduction. It uses a multi-pass approach to distinguish biological signal from technical noise.

### Configuration Table

| Parameter           | Value            | Description                                                                                                                                                                 |
| :------------------ | :--------------- | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **`ArchRProj`**     | `projHeme2`      | The **ArchRProject** object to which the dimensionality reduction will be added.                                                                                            |
| **`useMatrix`**     | `"TileMatrix"`   | The input data matrix. Using the `TileMatrix` (500bp windows) allows for an unbiased initial pass before peaks are even called.                                             |
| **`name`**          | `"IterativeLSI"` | The name given to this specific reduction. This allows you to store multiple runs (e.g., with different parameters) in the same project.                                    |
| **`iterations`**    | `2`              | The number of times the LSI process is repeated. The first pass finds broad clusters; subsequent passes use features variable across those clusters to refine the results.  |
| **`clusterParams`** | `list(...)`      | A list of parameters passed to the clustering algorithm (uses `Seurat::FindClusters`). These "internal" clusters are used to identify variable features between LSI rounds. |
| **`resolution`**    | `0.2`            | The granularity of the internal clustering. Higher values lead to more clusters. For the first iteration, a low resolution is preferred to capture major cell lineages.     |
| **`sampleCells`**   | `10000`          | The number of cells sampled to perform the "Estimated LSI" procedure. This allows the function to scale to millions of cells without crashing your RAM.                     |
| **`n.start`**       | `10`             | The number of random starting points for the K-means clustering step, ensuring the identified clusters are stable.                                                          |
| **`varFeatures`**   | `25000`          | The number of top "variable features" (tiles or peaks) used for the final LSI. These features drive the separation in your UMAP.                                            |
| **`dimsToUse`**     | `1:30`           | The LSI dimensions (components) to be retained. Typically, the first 30 dimensions capture the majority of biological variance. __Biological Intuition:__ If your clusters don't make sense (e.g., they don't match known markers), you can manually override the algorithm. __Action:__ Change dimsToUse = 1:30 to dimsToUse = 2:30. This completely ignores the first dimension and forces the UMAP/Clustering to rely on the remaining, usually cleaner, dimensions.                                             |

---



![alt text](image-17.png) ![alt text](image-18.png)

## Version 1 Interpretation: Broad Lineage Identification

Version 1 represents the baseline "QC pass" for dimensionality reduction. Using **2 iterations**, a fixed **resolution of 0.2**, and **25,000 variable features**, this run provides a broad, low-resolution view of the major cell lineages in the dataset.

---

### 1. Global Topology and Structure
The UMAP in Version 1 shows cells forming a largely **continuous and singular mass**.
* **Connectivity:** The different biological groups are not yet clearly separated into distinct "islands".
* **Feature Impact:** By using 25,000 features, the model retains a high amount of genomic information, including "housekeeping" accessibility common to many cells.
* **Biological Signal:** This results in a "cloud-like" structure where broad differences are visible, but subtle transition states between closely related cell types remain blurred.



### 2. Cluster Granularity
The algorithm identifies **7 distinct clusters** (C1 through C7).
* **Major Lineages:** These clusters capture the primary biological "continents" of the data, such as mature vs. progenitor populations.
* **Cluster 2 (Dark Blue):** This group is positioned at the bottom of the mass, clearly separating a specific lineage (CD34+ progenitors) from the rest of the cells.
* **Resolution Limits:** Because only two iterations were performed with a 0.2 resolution, the algorithm misses smaller sub-populations that only become distinct when background noise is further filtered.

### 3. Sample Integration and Batch Effects
The `SampleName` plot confirms that the dimensionality reduction successfully handled technical variation between samples.
* **Successful Mixing:** Samples **1-scATAC_BMMC_R1** (red) and **3-scATAC_PBMC_R1** (green) are well-integrated and overlap extensively.
* **Biological Isolation:** Sample **2-scATAC_CD34_BMMC_R1** (blue) remains relatively isolated in Cluster 2.
* **Interpretation:** This confirms that the separation seen is driven by **biology** (progenitor status) rather than **technical batch**, which would have caused the samples to group separately by name.

---

### Summary Table for Version 1

| Metric                 | Observation          | Interpretation                                                                            |
| :--------------------- | :------------------- | :---------------------------------------------------------------------------------------- |
| **Cluster Count**      | 7 Clusters           | Identifies major cell populations; lacks sub-type resolution.               |
| **Topology**           | Connected "Mass"     | High feature count (25k) preserves common signals, keeping clusters close.    |
| **Integration**        | Well-mixed Red/Green | Strong batch correction between BMMC and PBMC samples[cite: 6, 8, 21].                    |
| **Lineage Definition** | Broad Lineages       | Ideal for a first-pass analysis to verify data quality and major cell groups. |



**Conclusion:** Version 1 serves as an essential check to ensure that different samples are integrating properly before moving to higher-resolution runs like Version 4.

## Version 4

We used this function for another dimensionality reduction, just to visualize differences:

```r
projHeme2 <- addIterativeLSI(
    ArchRProj = projHeme2,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI2", 
    iterations = 4, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.1, 0.2, 0.4), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 15000, 
    dimsToUse = 1:30
)
```

The output of the 4th iteration is as follows:

![alt text](image-19.png), ![alt text](image-20.png)

---

### 1. Global Topology: The "Two-Island" Structure
The most significant change in Version 4 is the shift in topology; the continuous mass from Version 1 has split into two distinct "islands".
* **Lineage Separation:** This clear physical separation on the UMAP indicates that the algorithm has successfully isolated the primary biological lineages (likely Lymphoid vs. Myeloid/Progenitor).
* **Reduced Background:** Lowering the `varFeatures` to 15,000 removed thousands of "housekeeping" peaks that previously blurred these groups together, allowing for a cleaner mathematical separation.
* **Distinct States:** The distance between these islands suggests that the chromatin accessibility profiles of these groups are fundamentally different, providing a superior starting point for cell-type identification.



### 2. Cluster Granularity (11 Clusters)
Version 4 identifies **11 clusters** (C1 through C11), a significant increase from the 7 clusters found in Version 1.
* **Sub-type Discovery:** The high final resolution (0.4) allows for the identification of subtle sub-populations. For example, the smaller island on the right is now divided into three distinct groups (Clusters 5, 6, and 7).
* **Stable Grouping:** On the larger left island, clusters 8, 9, 10, and 11 show a clear "flow" or transition, which likely represents different stages of cell differentiation or maturation.

### 3. Advanced Integration & Batch Correction
The `SampleName` plot for Version 4 shows excellent integration, confirming that the refinement is biological, not technical.
* **Seamless Mixing:** Samples **1-scATAC_BMMC_R1** (red) and **3-scATAC_PBMC_R1** (green) are perfectly blended within the same clusters on both islands.
* **The "Anchor" Effect:** By starting at a very low resolution (0.1), ArchR first identified the major biological groups across all samples. This ensured that the subsequent higher-resolution passes refined the biological signal rather than separating cells by their original sample batch.
* **Progenitor Status:** The **2-scATAC_CD34_BMMC_R1** sample (blue) remains a distinct population (concentrated in Clusters 2, 3, and 4), confirming that these progenitor cells have a unique and recognizable chromatin landscape.

---

### Comparison Summary: Version 1 vs. Version 4

| Metric       | Version 1 (The Rough Draft)                        | Version 4 (The Final Map)                                         |
| :----------- | :------------------------------------------------- | :---------------------------------------------------------------- |
| **Topology** | Single, connected "cloud". | Two distinct, lineage-specific islands. |
| **Clusters** | 7 Clusters (Broad) .      | 11 Clusters (High-resolution)].          |
| **Noise**    | Higher (25k features).     | Lower (15k features); more specific.    |
| **Use Case** | Quick verification of sample mixing.               | Deep discovery of sub-types and lineages.                         |




### Understanding the "Iterative" Logic

ArchR doesn't just run a single PCA-like transformation. It follows a "Zoom-In" logic:

1. **Iteration 1 (The Wide Lens):** Uses the most accessible tiles to find the "big" differences (e.g., Lymphocytes vs. Myeloid cells).
2. **Feature Selection:** It identifies which peaks are actually different *between* those big groups.
3. **Iteration 2 (The Macro Lens):** Re-runs the LSI using only those specific, biologically informative peaks. 



### Important Note: Dimension 1 Correlation
In scATAC-seq, **LSI Dimension 1** is often highly correlated with **Sequencing Depth** (how many fragments a cell has). 
* **Action:** After running this command, check the correlation. 
* **Adjustment:** If LSI1 is purely technical noise, you should change `dimsToUse` to `2:30` in your downstream UMAP and Clustering functions.




## 6.3 Estimated LSI: Scaling to Massive Datasets

When working with "atlas-scale" datasets (hundreds of thousands to millions of cells), computing a full LSI matrix is often impossible due to RAM (memory) limitations. **Estimated LSI** is ArchR's optimization for these massive projects.

---

### 1. The Core Concept: "Landmarks & Projection"
Instead of calculating the mathematical relationship between every single cell at once, ArchR uses a representative subset to build the coordinate system.

* **Landmark Cells:** ArchR selects a random subset of cells to define the "rules" of the LSI space.
* **Projection:** The remaining cells are then "projected" into that pre-defined space.
* **Disk-Backed Processing:** Because ArchR reads the non-landmark cells from the Arrow files on your disk one by one, it never needs to load the entire dataset into your RAM at once.



---

### 2. The Step-by-Step Workflow
1.  **Selection:** A specific number of "landmark" cells is randomly sampled from the project.
2.  **LSI on Landmarks:** ArchR performs standard LSI (TF-IDF and SVD) on only these landmarks.
3.  **IDF Transfer:** The **Inverse Document Frequency (IDF)** values—the weights that tell us which peaks are most informative—are saved from the landmark set.
4.  **Projection:** All other cells are normalized using those landmark weights and moved into the landmark-defined space.

---

### 3. Key Parameters in `addIterativeLSI()`
You trigger Estimated LSI by adding these two parameters to your function call:

| Parameter              | Function                                                                               |
| :--------------------- | :------------------------------------------------------------------------------------- |
| **`sampleCellsFinal`** | The number of cells to be used as the "landmark" set (e.g., 50,000).                   |
| **`projectCellsPre`**  | A logical (`TRUE`/`FALSE`) that tells ArchR to use the landmark subset for projection. |

---

### 4. Important Considerations
* **Landmark Diversity:** The landmark set must be large enough to capture rare cell types. If a cell type makes up only 0.1% of your data and you only pick 1,000 landmarks, you might miss that population entirely in your initial map.
* **Memory vs. Accuracy:** While estimated LSI is slightly less precise than a full LSI calculation, the difference is negligible for large datasets and is often the only way to process them on standard hardware.

Estimated LSI allows you to analyze **millions of cells** without needing a supercomputer. It builds a high-quality "anchor" map using a subset of your data and then fits the rest of the cells into that map efficiently.

## 6.4 Batch Effect Correction with Harmony

While ArchR's **Iterative LSI** is designed to minimize technical noise, some datasets exhibit "strong" batch effects where cells group by sample or processing date rather than biology. In these cases, ArchR uses **Harmony**, a popular batch-correction algorithm originally developed for scRNA-seq.

---

### 1. What is Harmony?
Harmony is an algorithm that "aligns" clusters across different batches. It projects cells into a shared space and then iteratively moves them until the clusters from different samples overlap, provided they share the same biological identity.

* **Input:** An existing dimensionality reduction (e.g., your `IterativeLSI` object).
* **Output:** A new, corrected dimensionality reduction object (e.g., `Harmony`).
* **The Goal:** To ensure that a "T-cell" from Sample A and a "T-cell" from Sample B are located at the same coordinates in your UMAP.



---

### 2. Parameter Breakdown: `addHarmony()`

| Parameter         | Value            | Description                                                                                   |
| :---------------- | :--------------- | :-------------------------------------------------------------------------------------------- |
| **`ArchRProj`**   | `projHeme2`      | Your active ArchRProject.                                                                     |
| **`reducedDims`** | `"IterativeLSI"` | The name of the input dimensionality reduction you want to correct.                           |
| **`name`**        | `"Harmony"`      | The name you want to give to the new, corrected reducedDims object.                           |
| **`groupBy`**     | `"Sample"`       | The column in your metadata that defines the batches (e.g., "Sample", "Replicate", or "Day"). |

---

### 3. Implementation Code

```R
# Correcting for batch effects using Harmony
projHeme2 <- addHarmony(
    ArchRProj = projHeme2,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
)

# Note: Harmony typically converges quickly (e.g., 3-10 iterations)
```

# Chapter 7: Clustering with ArchR - Basics & Biological Background

Clustering is a fundamental step in single-cell analysis that allows us to organize thousands of individual cells into meaningful groups based on their molecular profiles.

### 1. Biological Significance: Why Cluster?
In a complex tissue sample, different cells perform different functions. These functions are driven by which parts of the DNA are "open" (accessible) and available for transcription. 
* **Identifying Cell Identity:** By grouping cells with similar chromatin accessibility patterns, we can identify distinct cell types (e.g., T-cells vs. B-cells) and cell states (e.g., resting vs. activated).
* **Unbiased Discovery:** Clustering is an "unsupervised" process. We don't tell the computer what cell types exist; instead, the computer finds groups that naturally exist in the data based on their shared epigenetic fingerprints.


---

### 2. The Logic: From LSI to Clusters
Clustering is rarely performed on the raw, sparse data matrix. Instead, it happens in the "Reduced Dimension" space created by LSI or Harmony.

1.  **Reduced Dimensions:** We use the coordinates from LSI (e.g., 30 dimensions) as the input.
2.  **Nearest Neighbor Graph:** The algorithm looks at every cell and identifies its closest "neighbors" in that 30D space.
3.  **Community Detection:** The computer then finds "communities" (clusters) of cells that are more densely connected to each other than to the rest of the graph.

---

### 3. Core Algorithmic Concepts
ArchR leverages state-of-the-art tools from the scRNA-seq world (like **Seurat**) to perform these calculations because they are robust and highly scaleable.

* **Graph-Based Clustering:** This is the industry standard. It treats cells like nodes in a social network and finds the "friend groups".
* **Louvain/Leiden Algorithms:** These are the specific mathematical formulas used to find these communities. In ArchR, this process is **deterministic**, meaning that if you run the same data twice with the same settings, you will get the exact same clusters.

---



| Concept       | Explanation                                                          |
| :------------ | :------------------------------------------------------------------- |
| **Input**     | Reduced dimensions (e.g., `IterativeLSI` or `Harmony`).              |
| **Goal**      | Discover distinct cell types/states without prior labeling.          |
| **Technique** | Nearest-Neighbor graph construction followed by community detection. |
| **Tools**     | ArchR wraps standard scRNA-seq methods (Seurat, scran).              |


**Next Steps:** Once clusters are identified, we can calculate **Marker Peaks** and **Gene Scores** to determine what those clusters actually represent (e.g., Cluster 1 = Monocytes).

## 7.1 Clustering using Seurat's FindClusters() function

In ArchR, the most common and successful method for identifying cell groups is the **graph-based clustering** approach implemented by the `Seurat` package. This method is deterministic and highly effective for high-dimensional single-cell data.

---

### 1. The Method: Graph-Based Community Detection
Instead of calculating the distance between every pair of cells (which is slow), this method builds a "web" of connections:
1.  **Nearest Neighbor Graph:** It identifies the $k$ most similar cells (neighbors) for every cell in the LSI/Harmony space.
2.  **SNN (Shared Nearest Neighbor):** It weights the edges between cells based on how many neighbors they share.
3.  **Louvain Algorithm:** It uses the Louvain community detection algorithm to find groups of cells that are more tightly knit than others.



---

### 2. Implementation: `addClusters()`
The `addClusters()` function acts as a wrapper that passes your dimensionality reduction coordinates directly to Seurat's clustering engine.

```R
projHeme2 <- addClusters(
    input = projHeme2,
    reducedDims = "IterativeLSI", # Or "Harmony"
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8              # Higher = more clusters; Lower = fewer
)
```

To access these clusters we can use the `$` accessor which shows the cluster ID for each single cell.
```r
head(projHeme2$Clusters)
## [1] "C9"  "C11" "C4"  "C4"  "C4"  "C7"
```

We can tabulate the number of cells present in each cluster:
```r
table(projHeme2$Clusters)
## 
##   C1  C10  C11  C12   C2   C3   C4   C5   C6   C7   C8   C9 
## 1532  903 1250  633 1120  314  351  386  702 1261 1377  421
```

To better understand which samples reside in which clusters, we can create a cluster confusion matrix across each sample using the confusionMatrix() function.
```r
cM <- confusionMatrix(paste0(projHeme2$Clusters), paste0(projHeme2$Sample))
cM
## 12 x 3 sparse Matrix of class "dgCMatrix"
##     scATAC_BMMC_R1 scATAC_CD34_BMMC_R1 scATAC_PBMC_R1
## C9             254                   5            162
## C11           1202                   .             48
## C4             351                   .              .
## C7             310                 940             11
## C1            1489                  10             33
## C6             171                 531              .
## C8             139                1238              .
## C12             86                   .            547
## C3             160                 144             10
## C10            322                   .            581
## C2             117                   2           1001
## C5              88                 298              .
```

To plot this confusion matrix as a heatmap, we use the pheatmap package:
```r 
library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM), 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black"
)
p
```
\
__As a result, we get this heatmap:__
![alt text](image-21.png)

There are times where the relative location of cells within the 2-dimensional embedding does not agree perfectly with the identified clusters. More explicitly, cells from a single cluster may appear in multiple different areas of the embedding. In these contexts, it may be appropriate to adjust the clustering parameters or embedding parameters until there is agreement between the two.



### 1. Understanding the Heatmap (The Confusion Matrix)
It shows a **normalized Confusion Matrix** visualized as a heatmap. It tells us exactly how our biological samples are distributed across the identified clusters.

* **Normalization:** The code `cM / Matrix::rowSums(cM)` is crucial. It ensures that the color intensity represents the **proportion** of each cluster's cells found in a given sample, rather than the raw count.
* **Reading the Blocks:**
    * **The CD34+ Progenitor Group:** Clusters **C3, C8, C5, C7, and C6** are almost exclusively found in the `scATAC_CD34_BMMC_R1` sample (the dark blue blocks on the left). This confirms these are progenitor-specific cell types.
    * **The BMMC Group:** Clusters **C4, C11, and C1** are specific to the `scATAC_BMMC_R1` sample.
    * **The PBMC Group:** Clusters **C10, C12, and C2** are dominant in the `scATAC_PBMC_R1` sample.
* **The "Shared" Clusters:** Notice **C9**. It has lighter blue across multiple columns, indicating it contains a mixture of cells from different samples—likely a common mature cell type present in both bone marrow and peripheral blood.


---

### 2. The Logic: Why Use Graph-Based Clustering?
ArchR uses the **Seurat** engine because graph-based clustering is "community-driven" rather than distance-driven.
1.  **SNN Graph:** It creates a web of connections where "friendships" (edges) are stronger if cells share many of the same neighbors.
2.  **Louvain Algorithm:** It identifies "neighborhoods" where cells are more connected to each other than to the rest of the map.
3.  **Determinism:** This process is **deterministic**—unlike some older methods, running this with the same settings will always yield the exact same cluster IDs.

---

### 3. Resolving Discrepancies: Clusters vs. Embeddings
The final paragraph of your text addresses a common frustration: **"Why does my UMAP look different than my Clusters?"**

* **30D vs. 2D:** Clustering happens in the high-dimensional LSI space (e.g., 30 dimensions). A UMAP is a **projection** into 2D. In 30 dimensions, two groups might be clearly separate, but the 2D "shadow" (UMAP) might make them look like they overlap.
* **Disconnected Clusters:** Sometimes, cells in Cluster 1 might appear in two different "islands" on your UMAP. 
    * **The Interpretation:** Usually, the high-dimensional clustering is more mathematically accurate than the 2D visualization.
* **How to Fix It:** If the disagreement is severe, you have two options:
    1.  **Adjust Clustering:** Change the `resolution` (lower it to merge split groups, raise it to separate them).
    2.  **Adjust Embedding:** Tweak UMAP parameters like `nNeighbors` or `minDist` to better represent the high-dimensional structure.

---

### Summary 

| Feature                     | Observation                                           | Action                                            |
| :-------------------------- | :---------------------------------------------------- | :------------------------------------------------ |
| **High Resolution (0.8)**   | Leads to more, smaller clusters.                      | Lower to 0.4 if clusters look too fragmented.     |
| **Sample-Specific Cluster** | Cluster appears in only one sample in the heatmap.    | Verify if this is real biology or a batch effect. |
| **Split Clusters on UMAP**  | One cluster ID is found in two separate UMAP islands. | Re-evaluate LSI dimensions or UMAP settings.      |


**Pro-Tip:** Always trust your Confusion Matrix over your UMAP. The heatmap tells you the truth about the data's structure, while the UMAP is just a "pretty picture" summary that can sometimes be misleading.

__IMPORTANT__: Clustering and Visualization are two different things! They use the same input but use completely different algorithms. 

# 8.0 Single-cell Embeddings

Embeddings take the high-dimensional results from IterativeLSI and project them onto a 2D plane. This is the stage where your data finally "looks" like a single-cell experiment. We call these “embeddings” because they are strictly used to visualize the clusters and are not used to identify clusters which is done in an LSI sub-space as mentioned in previous chapters.

### 8.1 Uniform Manifold Approximation and Projection (UMAP)

UMAP is the preferred embedding in ArchR because of its speed and ability to preserve both local and global relationships between cells.

#### Implementation Code

```r
# Calculate UMAP based on the IterativeLSI results
projHeme2 <- addUMAP(
  ArchRProj = projHeme2, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, # Local vs Global: Higher values create a more "global" view; lower values focus on very local similarities.
  minDist = 0.5, #Tightness: Controls ho wtightly packed the points are. Lower values (e.g., 0.1) result in tighter, denser clusters.
  metric = "cosine", #The distance math used. "cosine" is typically better for sparse scATAC-seq data than "euclidean".        
  force = TRUE
)
```

You can list the available embeddings objects in an ArchRProject using the slot extraction opperator @:

`projHeme2@embeddings`

To plot the UMAP results, we use the `plotEmbedding()` function and pass the name of the UMAP embedding we just generated (“UMAP”). We can tell ArchR how to color the cells by using a combination of `colorBy` which tells ArchR which matrix to use to find the specified metadata column provided to `name`.
```r
p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
```
![alt text](image-22.png)
---
* Instead of coloring by “Sample” as above, we can color by “Clusters” which were identified in a previous chapter.
```r 
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
```
![alt text](image-23.png)
---
* We can visualize these two plots side by side using the ggAlignPlots() function, specifying a horizontal orientation using type = "h"
```r
ggAlignPlots(p1, p2, type = "h")
```
---
To save an editable vectorized version of this plot, we use `plotPDF()`
```r
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
```
### Interpreting Single-Cell UMAPs

The UMAP (Uniform Manifold Approximation and Projection) is a 2D map where each dot represents a single cell. Cells that are plotted close together have highly similar chromatin accessibility profiles (and are therefore likely the same cell type), while cells located far apart are biologically distinct. 

Visualizing the UMAP by different metadata columns is crucial for validating both the technical quality and the biological reality of your dataset.

#### 1. Plot 1: UMAP Colored by Sample (Validating Integration)

This plot serves as your primary **Batch Effect Check**. It shows the contribution of each individual sample to the overall map:
* **Red (`scATAC_BMMC_R1`):** Bone Marrow Mononuclear Cells.
* **Blue (`scATAC_CD34_BMMC_R1`):** CD34+ enriched Bone Marrow Cells (Stem/Progenitor cells).
* **Green (`scATAC_PBMC_R1`):** Peripheral Blood Mononuclear Cells.

**What does this tell us?**
* **Biological Integration:** The Red (BMMC) and Green (PBMC) cells overlap heavily in the bottom-left "island" (Clusters 9, 10, 11) and the bottom-right tail (Cluster 3). This is an excellent result. It demonstrates that mature immune cells found in both the bone marrow and the blood (such as T-cells or B-cells) have nearly identical chromatin profiles, regardless of their extraction site.
* **Biological Distinctness:** The top right section of the UMAP (Clusters 6, 7, 14) is dominated by the Blue (CD34+) and Red (BMMC) samples, with almost no Green (PBMC) cells. This aligns perfectly with known biology: CD34+ progenitor cells reside in the bone marrow and are rarely found in circulating peripheral blood.
* **Conclusion:** Because the overlap aligns with biological expectations rather than segregating purely by sample origin, **we do not have a severe batch effect.** The integration is successful, and we can proceed without running Harmony batch correction.

#### 2. Plot 2: UMAP Colored by Clusters (Validating Structure)

This plot displays the result of the Seurat-based clustering algorithm. ArchR has mathematically partitioned the cells into 14 distinct populations based on their LSI coordinates.

**What does this tell us?**
* **Discrete Populations (The "Islands"):** The large, detached group on the left side of the UMAP (Clusters 9, 10, 11) represents a cell lineage that is highly distinct from the rest of the cells. In blood and bone marrow datasets, an isolated island like this is very often the Lymphoid lineage (T-cells, B-cells, NK cells).
* **Continuous Populations (The "Continuum"):** The large, interconnected mass on the right side of the plot (Clusters 1 through 7, 12, 14) represents a developmental continuum. Notice the lack of hard visual breaks between these clusters. This visualizes a biological trajectory—likely hematopoietic stem cells (top right) slowly differentiating into mature myeloid cells (bottom right) through various intermediate, transitional states.
* **Cluster Density:** Cluster 2 (dark blue, bottom right) is very dense and tightly packed, meaning those cells are highly homogenous. In contrast, Cluster 13 (light purple, top center) is much more diffuse, indicating greater biological variability and a broader gradient of cell states within that specific grouping.

> **Thesis Tip:** When presenting these UMAPs in your results chapter, always show the "Sample" and "Cluster" plots side-by-side. This immediately answers the two most common questions from reviewers: "Did the samples integrate properly across different biological replicates?" and "How many biologically distinct populations did the algorithm successfully identify?"
>


## 8.2 t-Stocastic Neighbor Embedding (t-SNE)

https://www.archrproject.com/bookdown/t-stocastic-neighbor-embedding-t-sne.html

## 8.3 Dimensionality Reduction After Harmony

If your initial UMAP colored by `Sample` showed significant technical separation—meaning cells of the same biological type were segregated purely by which sequencing run they belonged to—you need to apply **Batch Correction**. 

ArchR natively integrates with **Harmony**, an algorithm designed to align single-cell data across different batches while preserving true biological variation.

### Running Harmony
Harmony works directly on the LSI coordinates. It pulls similar cells from different batches together in high-dimensional space. 

```r
# Apply Harmony batch correction to the LSI dimensions
projHeme2 <- addHarmony(
    ArchRProj = projHeme2,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample", # The metadata column causing the batch effect
    force = TRUE
)
```

### Calculate a new UMAP using the Harmony-corrected dimensions
```r
projHeme2 <- addUMAP(
    ArchRProj = projHeme2, 
    reducedDims = "Harmony", # CRITICAL: Point to the corrected data!
    name = "UMAPHarmony",    # Give it a unique name so you don't overwrite the old UMAP
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = TRUE
)
```

### Plot the new Harmony UMAP colored by Sample
```r
p3 <- plotEmbedding(
    ArchRProj = projHeme2, 
    colorBy = "cellColData", 
    name = "Sample", 
    embedding = "UMAPHarmony"
)
```
### Plot the new Harmony UMAP colored by Clusters
```r 
p4 <- plotEmbedding(
    ArchRProj = projHeme2, 
    colorBy = "cellColData", 
    name = "Clusters", 
    embedding = "UMAPHarmony"
)
```

### Display or save
```r
plotPDF(p3, p4, name = "Plot-UMAP-Harmony-Sample-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE)
```

![alt text](image-24.png)

![alt text](image-25.png)

### 8.5 Visual Comparison: Uncorrected LSI vs. Harmony Batch Correction

When integrating multiple single-cell samples, it is critical to determine whether the cells are grouping by **true biological state** or by **technical artifacts** (batch effects). The side-by-side comparison of your UMAPs before and after Harmony provides a perfect textbook example of why batch correction is often necessary.

#### 1. The Pre-Harmony (Uncorrected) UMAPs
*Referencing Pages 1 & 2 (UMAP of IterativeLSI)*

* **The Sample Separation:** In the uncorrected UMAP colored by `Sample`, there is distinct spatial segregation between the biological replicates. For instance, `scATAC_PBMC_R1` (Green) and `scATAC_BMMC_R1` (Red) form adjacent but largely separated territories.
* **The Problem:** Biologically, we know that both Bone Marrow (BMMC) and Peripheral Blood (PBMC) contain identical mature immune cells (like circulating T-cells). In a perfect world, the red and green dots for these specific cell types should be perfectly stacked on top of each other. Their separation indicates that the LSI algorithm is picking up on the "technical signature" of the sequencing run rather than just the biology.
* **Cluster Artifacts:** Because of this separation, the 14 clusters identified are likely biased. Some clusters might simply be "PBMC T-cells" and "BMMC T-cells" split into two groups, which artificially inflates your cluster count.

#### 2. The Post-Harmony (Corrected) UMAPs
*Referencing Pages 3 & 4 (UMAPHarmony of Harmony)*

* **The Sample Integration:** The `UMAPHarmony` colored by `Sample` represents the data after the Harmony algorithm has aligned the shared biological states. 
* **The Fix:** Harmony effectively forces cells with similar biological signatures to overlap in high-dimensional space, regardless of their origin batch. You will typically see the Red and Green samples thoroughly mixed within the mature immune cell clusters. Crucially, Harmony is "biology-aware"—it should leave the unique `scATAC_CD34_BMMC_R1` (Blue) progenitor cells distinct, as they do not have a biological equivalent in the peripheral blood sample.
* **Refined Biological Structure:** The topology of the `UMAPHarmony` colored by `Clusters` shifts significantly. We still see the detached "lymphoid" island (Clusters 9, 10, 11) and the main developmental continuum (Clusters 1 through 7), but the cells within them are now organized by their true biological identity. 

#### 3. Discussion & Conclusion

**Why does this matter?**
If you proceed with the uncorrected IterativeLSI dimensions, your downstream analysis will be flawed. For example, if you look for differentially accessible peaks between Cluster X and Cluster Y, you might just be finding the technical noise between the PBMC and BMMC sequencing runs. 

By using the **Harmony-corrected UMAP**, you ensure that your clusters represent genuine cell types and states. This makes your subsequent steps—identifying marker genes, assigning cell type labels, and building developmental trajectories—biologically accurate and highly robust.

> **Thesis Tip:** Include these four plots as a 2x2 multi-panel figure in your methodology or supplementary section. Use it to explicitly justify your use of Harmony. A statement like, *"Initial dimensionality reduction via IterativeLSI revealed prominent batch effects driven by sample origin (Fig A). Application of Harmony successfully integrated shared mature lineages across BMMC and PBMC samples while preserving the unique CD34+ progenitor populations (Fig C),"* demonstrates deep analytical rigor to your committee.
>

### 8.4 Highlighting specific cells on an embedding

#### Basics
https://www.archrproject.com/bookdown/highlighting-specific-cells-on-an-embedding.html

#### Identifying Cell Types with Gene Scores

In scATAC-seq, we do not directly measure RNA transcripts. Instead, ArchR calculates a **Gene Score** by summing the chromatin accessibility signal across the gene body and its surrounding regulatory elements. This score serves as a highly accurate proxy for gene expression, allowing us to identify cell clusters using known biological markers.

#### 1. Highlighting Specific Clusters
Sometimes, projecting a marker gene across the entire UMAP can be noisy. To clean up the visualization, you can use the `highlightCells` parameter to isolate specific clusters of interest. All other cells will be rendered as a grey background.

```r
# Plot CD14 Gene Scores, highlighting ONLY Clusters 1 through 5
plotEmbedding(
  ArchRProj = projHeme2,
  embedding = "UMAP",
  colorBy = "GeneScoreMatrix",
  name = "CD14",
  size = 1,
  sampleCells = NULL,
  highlightCells = getCellNames(ArchRProj = projHeme2)[which(projHeme2@cellColData$Clusters %in% c("C1","C2","C3","C4","C5"))],
  baseSize = 10,
  plotAs = "points"
)
```

![alt text](image-26.png)
#### 2. __Interpreting the Marker Plot__

  The Marker: CD14 is a standard marker for the monocyte lineage.  

  The Visual: The plot utilizes a heatmap gradient (Log2(NormCounts + 1)).  

  Grey Cells: Cells outside of our target clusters (C1-C5).  

  Dark Blue Cells: Target clusters with closed, inaccessible chromatin at the CD14 locus.

  Yellow/Pink Cells: Target clusters with highly accessible CD14 chromatin, confirming their identity as monocytes.

  Biological Takeaway: The presence of a gradient within this highlighted continuum indicates a developmental trajectory. We can visually trace the maturation of these cells as the CD14 gene becomes progressively more accessible.

  * Thesis Tip: When assigning identities to your clusters, single-marker UMAPs are highly persuasive. If you label a cluster as "Monocytes" in your text, referencing a supplemental figure showing bright yellow CD14 gene scores precisely over that cluster provides the necessary biological proof for your computational claims.

## 8.5  Importing an embedding from external software

https://www.archrproject.com/bookdown/importing-an-embedding-from-external-software.html

# 9.0 Gene Scores and Marker Genes: The Biological Background

To effectively analyze your Atrial Fibrillation (AF) scATAC-seq dataset, you must translate raw chromatin accessibility peaks into biologically meaningful units: **genes**. Because you are working purely with DNA accessibility (not RNA transcripts), understanding how ArchR bridges this gap is fundamental to interpreting your cardiomyocyte clusters and mapping variant effects.

### 1. The scATAC-seq Challenge: We Don't Measure RNA
In a standard scRNA-seq experiment, you directly count the mRNA transcripts produced by a cell to determine which genes are turned "on." In scATAC-seq, you do not measure RNA. Instead, you measure the physical "openness" of the chromatin. 

While open chromatin is a prerequisite for transcription, it is not a direct 1:1 measurement of it. A region might be open because a transcription factor is bound, but the gene isn't actively firing yet (a "poised" state). Therefore, we must *infer* gene expression based on the epigenomic landscape.

### 2. What is a "Gene Score"?
A **Gene Score** (or Gene Activity Score) is ArchR’s mathematical prediction of how highly expressed a gene is, based entirely on the surrounding open chromatin. 



ArchR calculates this score by looking at a massive genomic window around every gene and summing the ATAC-seq signal based on specific biological rules:
* **The Promoter:** The region immediately upstream of the Transcription Start Site (TSS) is given the highest weight. If the promoter is closed, the gene is almost certainly off.
* **The Gene Body:** Accessibility across the actual coding region of the gene is also heavily weighted, as open chromatin here indicates active transcription machinery moving through the DNA.
* **Distal Enhancers (Distance-Weighted):** ArchR looks at open peaks far away from the gene (up to hundreds of kilobases). Because enhancers loop over to touch promoters in 3D space, open enhancers contribute to the Gene Score. ArchR uses a "distance decay" model, meaning a peak 10kb away contributes more to the score than a peak 100kb away.

### 3. What is a "Marker Gene"?
A **Marker Gene** is a gene whose high expression (and therefore, high Gene Score) uniquely defines a specific cell identity, anatomical region, or disease state. In your thesis, you will use marker genes to annotate your UMAP clusters.

Since you are comparing different regions of the heart, you must rely on region-specific cardiomyocyte (CM) markers:
* **Pan-Cardiomyocyte Markers:** *TNNT2*, *MYH6*. If a cluster has high gene scores for these, it is a CM (not a fibroblast or endothelial cell).
* **Atrial Markers:** *NPPA*, *MYL4*. 
* **Left-Atrium Specific:** *PITX2* (This is the most critical marker for your AF thesis, as *PITX2* defines left-atrial identity).
* **Ventricular Markers:** *MYL2*, *MYH7*. 

When you plot the Gene Score of *MYL2* on your UMAP, it will "light up" the ventricular clusters and remain dark in the atrial clusters.

### 4. Relevance to Atrial Fibrillation (AF) Pathophysiology
During Atrial Fibrillation, the atria undergo massive electrical and structural remodeling. This means the epigenetic landscape changes, and consequently, the Gene Scores will shift.

* **Fetal Gene Program:** Stressed cardiomyocytes often revert to a fetal state. You may see AF clusters showing high Gene Scores for fetal markers (*NPPA*, *NPPB*) compared to healthy tissue.
* **Ion Channel Remodeling:** You can use Gene Scores to investigate whether the chromatin around key potassium and sodium channels (like *KCNQ1* or *SCN5A*) closes during AF, leading to the electrical chaos characteristic of the disease.

### 5. Connecting to Variant Effect Prediction (VEP)
This is the core of your thesis. Why do we care about Gene Scores when looking at genetic variants?

GWAS studies have identified hundreds of single nucleotide polymorphisms (SNPs) associated with AF. Over 90% of these sit in non-coding enhancers, not in the genes themselves. 
1. **The Problem:** If a SNP sits in an enhancer in the middle of nowhere, how do you know which gene it causes to malfunction?
2. **The ArchR Solution:** Because ArchR's Gene Scores incorporate distal enhancers into their calculations, you can mathematically link an enhancer peak to a specific gene's promoter. 
3. **Variant Effect Prediction:** If an AF patient has a mutation in an enhancer, your predictive model will try to determine if that mutation destroys a Transcription Factor binding site. If it does, the enhancer closes. Because the enhancer closes, the **Gene Score** for the linked target gene (e.g., *PITX2*) drops. This mechanical chain of events is what you are predicting.

> **Thesis Application:** In your analysis pipeline, you will first use well-established Marker Genes to confidently label your LA, RA, LV, and RV cardiomyocyte clusters. Once the identities are locked in, you will calculate Gene Scores across the entire genome to find out which specific genes are epigenetically silenced or activated during Atrial Fibrillation.


## 9.1 How ArchR Calculates Gene Scores: The Default Model

In scATAC-seq, predicting RNA expression from chromatin accessibility is a complex mathematical challenge. The paragraph from the ArchR tutorial highlights the winning algorithm they developed after testing 50 different variations. 

To make your Variant Effect Prediction (VEP) accurate, you need to understand how ArchR decides which ATAC-seq peaks belong to which genes. The model relies on three fundamental biological rules:

#### 1. Accessibility within the entire gene body
* **The Concept:** Traditional models only looked at the promoter (the region right before the gene starts). ArchR looks at the promoter *and* the entire coding sequence of the gene itself.
* **The Biology:** When a gene is actively being transcribed, RNA Polymerase physically moves through the gene body, forcing the chromatin to open up. Therefore, a high concentration of ATAC-seq reads spanning the entire length of a gene is a massive biological indicator that the gene is turned "on."

#### 2. Exponential weighting of distal regulatory elements
* **The Concept:** ArchR includes peaks that are located far away from the gene (distal elements), but it applies a "distance penalty." 
* **The Biology:** Distal elements are **enhancers**. Because DNA exists in 3D space, an enhancer can fold over to touch a promoter. However, the further away an enhancer is on the linear DNA strand, the less likely it is to interact with that specific gene. 
* **The Math:** ArchR uses an *exponential decay* function. An open peak 10 kilobases (kb) away from the *PITX2* promoter will add a significant amount to the *PITX2* Gene Score. A peak 100 kb away will add a much smaller fraction. A peak 500 kb away will add almost nothing.



#### 3. Imposed gene boundaries
* **The Concept:** ArchR prevents enhancers from "jumping" over neighboring genes to artificially inflate a Gene Score.
* **The Biology:** The genome is organized into insulated neighborhoods (Topologically Associating Domains, or TADs). Enhancers are generally restricted to acting upon the genes within their own neighborhood. If ArchR sees an enhancer, it looks for the closest gene. If there is another gene sitting *between* the enhancer and your target gene, ArchR assumes the enhancer belongs to the neighbor, establishing a "boundary" that blocks the signal.
* **Why it matters:** This drastically reduces false positives. Without these boundaries, a massive enhancer for a cardiac muscle gene might accidentally inflate the Gene Score of an unrelated neighboring gene.

> **Thesis Tip:** When discussing your Variant Effect Prediction methodology, explicitly mention that you utilized ArchR's distance-weighted, boundary-imposed Gene Score model. This demonstrates to your committee that your mapping of non-coding AF variants to target genes accounts for the complex 3D folding and boundary insulation of the human genome, rather than just simple linear proximity.

## 9.2 Identification of Marker Feautures

![alt text](image-27.png)



After calculating UMAPs and establishing your clusters, you are left with mathematical groupings of cells (e.g., Cluster 1, Cluster 2). To turn these abstract numbers into biological cell types (e.g., Left Atrial Cardiomyocytes, Fibroblasts, Macrophages), you must identify **Marker Features**.

A marker feature is any biological metric (a Gene Score, a specific ATAC-seq peak, or a Transcription Factor motif) that is uniquely accessible in one specific cluster compared to the rest of the dataset.

#### 1. The scATAC-seq Bias Problem
Finding markers in scATAC-seq is notoriously tricky due to technical noise. If Cluster A simply has a higher sequencing depth (more fragments per cell) than Cluster B, a naive statistical test will tell you that *every* gene is a marker for Cluster A, purely because there is more data there.

#### 2. ArchR's Solution: Bias-Matched Backgrounds
ArchR overcomes this using a highly robust algorithm within the `getMarkerFeatures()` function. When testing if a gene is a marker for Cluster 1, ArchR does not just compare Cluster 1 to all other cells. Instead, it carefully selects a "background" group of cells from the other clusters that perfectly match the cells in Cluster 1 based on two critical quality metrics:
* **TSS Enrichment** (Signal-to-noise ratio)
* **log10(nFrags)** (Sequencing depth)

By matching these biases, ArchR ensures that the markers it finds are driven by true biological differences, not technical artifacts.

#### 3. Executing the Marker Search
For your initial cell type annotation, you will want to find marker **Gene Scores**. This function performs a Wilcoxon rank-sum test to identify genes that are significantly more accessible in each cluster.

```r
# Identify marker Gene Scores for all clusters
markersGS <- getMarkerFeatures(
    ArchRProj = projHeme2, 
    useMatrix = "GeneScoreMatrix", #We use "GeneScoreMatrix" first to find marker genes. Later, you can change this to "PeakMatrix" to find marker enhancers.
    groupBy = "Clusters", #The metadata column containing your groups. Usually "Clusters".
    bias = c("TSSEnrichment", "log10(nFrags)"), #The technical metrics ArchR must control for to prevent false positives.
    testMethod = "wilcoxon"
)
```

#### 4. Visualitzing the Markers

* __Heatmaps__

The output of `getMarkerFeatures()` is a massive matrix of p-values and Fold Changes. The best way to view the top markers across all clusters simultaneously is a __Marker Heatmap__. 
```r 
# Extract the top 40 marker genes per cluster
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

# Generate the Heatmap
heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  nLabel = 3, # Labels the top 3 genes per cluster on the y-axis
  transpose = TRUE
)

# Draw the plot
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
```

![alt text](image-29.png)
### 10.1 Interpreting the Marker Gene Heatmap

The marker gene heatmap is essentially the "decoder ring" for your scATAC-seq dataset. It takes the abstract mathematical clusters from your UMAP and assigns them concrete biological identities based on the accessibility of known regulatory genes.

Here is a detailed breakdown of how to read this specific plot and what it reveals about your blood and bone marrow cells.

#### 1. Understanding the Plot Anatomy
* **The Rows (C1 - C14):** Each row represents one of the 14 clusters ArchR identified in your dataset.
* **The Columns (2645 features):** Each thin vertical line is a specific gene that passed your statistical threshold (`FDR <= 0.01 & Log2FC >= 1.25`). 
* **The Color Scale (Z-Scores):** The heatmap does not show absolute expression; it shows *relative* accessibility across clusters. 
    * **Yellow (+2):** This gene's chromatin is highly open/accessible in this specific cluster compared to the average.
    * **Blue (-2):** This gene's chromatin is closed/inaccessible in this cluster.
* **The Labels (Top):** ArchR looks at your custom `markerGenes` list and places a text label above the specific column where that gene is plotted.

#### 2. Biological Annotation of Your Clusters
By matching the bright yellow "blocks" to the labels at the top, we can confidently assign cell types to your clusters:

* **The B-Cell Lineage (Clusters 12, 13, 14):** Look at the labels for **EBF1, MS4A1, PAX5, and MME**. Directly below these labels, you see a massive block of bright yellow exclusively in rows C12, C13, and C14. This perfectly confirms that these clusters represent your B-cell populations.
* **The T-Cell Lineage (Clusters 8, 9, 10, 11):** The T-cell markers are split, which reveals sub-types! **CD3D and IL7R** (often associated with naive or helper T-cells) are highly accessible in C8, C9, and C10. However, **TBX21 and CD8A** (markers for cytotoxic CD8+ T-cells) light up brightly in C11.
* **The Erythroid Lineage (Clusters 6, 7):**
  **GATA1**, a master transcription factor for red blood cell development, shows strong accessibility in C6 and C7.
* **The Myeloid/Monocyte Lineage (Clusters 4, 5):**
  **MPO** (Myeloperoxidase), a classic myeloid/granulocyte marker, is highly enriched in C4 and C5.

#### 3. Troubleshooting Missing Labels
You might notice that some genes from your R code (like `CD34`, `CD14`, `IRF8`) are missing from the top of the heatmap. 

**Why does this happen?**
1. **Statistical Cutoff:** A gene is only plotted if it passes the `cutOff = "FDR <= 0.01 & Log2FC >= 1.25"` threshold in *at least one* cluster. If `CD14` didn't meet this strict fold-change requirement, ArchR drops it from this specific plot to prevent visualizing statistical noise.
2. **Visual Overlap:** Sometimes, if two marker genes are located right next to each other in the matrix, ArchR drops one label to prevent the text from overlapping and becoming unreadable.

**The Solution:** Notice the large, unlabelled bright yellow blocks in **C1, C2, and C3**. Given your input list, these are highly likely your `CD34+` Early Progenitors. To prove this, you can generate a specific UMAP colored by the `CD34` Gene Score (just like we did earlier for CD14) to visually confirm its presence, bypassing the strict cutoff of the heatmap.

> **Thesis Tip:** In your final manuscript, this heatmap is the definitive proof of your cell typing. You will state: "Clusters were annotated based on the differential accessibility of canonical lineage markers (Fig X). For example, clusters 12-14 were annotated as B-cells due to significant enrichment of *MS4A1* and *PAX5* gene scores."

### __Volcano Plots__

If you want to focus deeply on a single cluster (for example, proving that Cluster 5 is your Left Atrial Cardiomyocyte population), you use a Volcano Plot. This plots every gene based on its Fold Change (x-axis) and its statistical significance (y-axis).
```r
# Plot a Volcano plot specifically for Cluster 5
p <- plotMarkers(
    seMarker = markersGS, 
    name = "C5", 
    cutOff = "FDR <= 0.01 & Log2FC >= 1", 
    plotAs = "MA" # Can be "MA" or "Volcano"
)
p
```
![alt text](image-28.png)

> **Thesis Tip:** In your AFib project, you will heavily rely on getMarkerFeatures(). First, you will use it with the GeneScoreMatrix to confirm which clusters are LA, RA, LV, and RV cardiomyocytes. Later, you will run this exact same function but group your cells by Disease State (e.g., AFib vs. Healthy) and use the PeakMatrix to find the specific enhancer peaks that structurally remodel during the disease!

### 9.4 Visualizing Marker Genes on an Embedding

While the Marker Gene Heatmap provides a fantastic global overview of all clusters at once, it abstracts away the relationships between the cells. To fully validate your cluster annotations and observe developmental trajectories, you must project these marker genes directly back onto your 2D UMAP.

This allows you to see not just *if* a gene is active in a cluster, but *how* it is active. Does it turn on suddenly, or does it gradually increase in accessibility along a continuous biological trajectory?

#### 1. Generating Multiple Marker UMAPs
Instead of plotting one gene at a time, ArchR allows you to pass your entire list of marker genes into the `plotEmbedding()` function. This will generate a list of UMAP plots—one for each gene—colored by its Gene Score.

```r
# Define the specific marker genes of interest
markerGenes <- c(
  "CD34", # Early Progenitor
  "GATA1", # Erythroid
  "PAX5", "MS4A1", "EBF1", "MME", # B-Cell Trajectory
  "CD14", "CEBPB", "MPO", # Monocytes
  "IRF8", 
  "CD3D", "CD8A", "TBX21", "IL7R" # TCells
)

# Generate a list of UMAP plots colored by GeneScoreMatrix
p <- plotEmbedding(
    ArchRProj = projHeme2, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    quantCut = c(0.01, 0.95) # Scales the color gradient
    imputeWeights = NULL
)

# To plot a specific gene, we can subset this plot list:
p$GATA1

# To plot all genes we can use `cowplot`to arrange the various marker genes into a single plot.
p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
```

#### Interpreting the CowPlot
  
  "C:\Users\marko\OneDrive - Universität Graz\Dokumente\Uni\Masterarbeit\Plots ArchR Tutorial\Plot_UMAP_Marker_Genes_WO_Imputation.pdf"

By overlaying the Gene Scores of specific markers directly onto your UMAP, you effectively create a biological map of your data. The cowplot grid you generated allows us to visually trace the exact developmental lineages and validate the identities of the discrete "islands" and "continuums" we saw in the clustering phase.

Here is the detailed interpretation of your multi-panel UMAP to include in your thesis results.

#### 1. The Progenitor Root (CD34)
* **Visual Signature:** *CD34* is heavily enriched at the very top of the large, interconnected continuum on the right side of the UMAP.
* **Biological Meaning:** CD34 is a classic marker for Hematopoietic Stem Cells (HSCs) and multipotent progenitors. Because it lights up exactly at the "apex" of the interconnected mass, we can definitively establish this top-right region as the root of your developmental trajectories. All other myeloid and erythroid cells branch out from this point.

#### 2. The Erythroid Branch (GATA1)
* **Visual Signature:** Moving slightly down and to the left from the CD34+ root, *GATA1* lights up in a distinct offshoot branch of the main continuum.
* **Biological Meaning:** GATA1 is the master transcription factor for red blood cell (erythrocyte) and megakaryocyte development. This shows a clear developmental fork: some progenitor cells are moving away from the main trunk to commit to the erythroid lineage.

#### 3. The Myeloid / Monocyte Trajectory (MPO & CD14)
* **Visual Signature:** If you follow the main trunk down from the CD34+ root, it forms a long, vertical tail. *MPO* (Myeloperoxidase) lights up brightly in the middle of this tail. *CD14* lights up intensely at the very bottom tip of this tail.
* **Biological Meaning:** This beautifully visualizes continuous cell differentiation! The cells flow from stem cells (CD34+), mature into early myeloid/granulocyte progenitors (MPO+ in the middle), and finally terminally differentiate into mature monocytes (CD14+) at the bottom.

#### 4. The B-Cell Island (PAX5, MS4A1, MME)
* **Visual Signature:** These three markers completely illuminate the detached island located at the top of the UMAP space. 
* **Biological Meaning:** * *PAX5* is a master B-cell lineage transcription factor.
    * *MS4A1* encodes the CD20 protein (the target of the drug Rituximab), marking mature B-cells.
    * *MME* encodes CD10, often marking early/pre-B cells.
    * Because this island is physically separated from the main CD34+ continuum, it indicates these are fully committed, mature cells circulating in the blood or bone marrow, distinct from the active myeloid differentiation happening on the right.

#### 5. The T-Cell Island (CD3D & CD8A)
* **Visual Signature:** *CD3D* lights up the entirety of the large, detached island on the far left. *CD8A* only lights up the bottom-right portion of that exact same island.
* **Biological Meaning:** This demonstrates ArchR's ability to capture both broad lineages and specific sub-types. *CD3D* is a pan-T-cell marker, confirming the entire left island consists of T-lymphocytes. *CD8A* highlights the specific spatial territory within that island occupied by Cytotoxic (CD8+) T-cells, leaving the dark regions to likely represent CD4+ Helper T-cells.

> **Thesis Tip:** When writing your figure legend for this cowplot, emphasize the spatial dynamics. Use phrasing such as: *"Feature plots of inferred gene activity demonstrate clear developmental trajectories. Stem cell marker CD34 localizes to the apex of the central continuum, which bifurcates into a GATA1+ erythroid branch and an MPO+/CD14+ myeloid differentiation axis. Lymphoid lineages form distinct, terminally differentiated clusters expressing canonical T-cell (CD3D+) and B-cell (MS4A1+) markers."* This clearly demonstrates your mastery of how spatial UMAP positioning reflects underlying cardiovascular and hematopoietic biology.

### 9.5 Marker Genes Imputation with MAGIC

In the previous section, you may have noticed that the UMAPs colored by Gene Scores look a bit "grainy." Some cells in the middle of a clear cluster might be dark (showing a score of 0), even though we know biologically they should be expressing that marker. This brings us to the biggest technical hurdle in scATAC-seq: **Data Sparsity**.

#### 1. The Sparsity Problem
scATAC-seq data is notoriously sparse, suffering heavily from "dropouts." A dropout occurs when a region of chromatin is genuinely open in a cell, but the sequencing machine simply failed to capture that specific DNA fragment by chance. Because we only have two copies of DNA per cell, missing just one or two fragments can completely wipe out the Gene Score for that cell.

#### 2. The Solution: MAGIC Imputation
To fix this visual noise, ArchR integrates **MAGIC** (Markov Affinity-based Graph Imputation of Cells). 



MAGIC works by looking at the cell's "neighborhood." If a single cell has a Gene Score of 0 for *CD14*, but the algorithm sees that its 30 closest neighbors (based on the overall LSI dimensions) all have extremely high *CD14* scores, MAGIC assumes the 0 is a technical dropout. It "borrows" the signal from the neighbors to smooth out and impute the missing value.

#### 3. Implementing MAGIC in ArchR
Adding MAGIC weights is computationally straightforward. ArchR calculates the diffusion matrix and stores the imputation weights directly in the project object.

```r
# Calculate and add MAGIC imputation weights to your project
projHeme2 <- addImputeWeights(ArchRProj = projHeme2)
```

#### 4. Visualizing Imputed Gene Scores
```r
#Plot the marker genes WITH MAGIC Imputation
p_imputed <- plotEmbedding(
    ArchRProj = projHeme2, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projHeme2) # This applies the smoothing!
)

# Save the smoothed plots to a new PDF
plotPDF(plotList = p_imputed, 
        name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
        ArchRProj = projHeme2, 
        addDOC = FALSE, 
        width = 5, 
        height = 5)
```
#### Important Caveats

MAGIC is a Visualization Tool, Not Ground Truth!
While MAGIC makes UMAPs look beautiful and makes developmental trajectories incredibly clear, you must be careful how you use it in your Atrial Fibrillation (AF) thesis.

* Use it for: Generating clean figures for your manuscript, confirming cluster annotations, and visually demonstrating continuous lineage trajectories.

* Do NOT use it for: Differential accessibility testing or your Variant Effect Prediction models. Because MAGIC forces neighbors to look similar, it artificially destroys biological variance and will create massive false-positive p-values if you run statistics on the imputed numbers. Always perform your core statistical tests (like `getMarkerFeatures()`) on the raw, un-imputed data.

* Thesis Tip: A great supplementary figure for your thesis would be a side-by-side comparison. Show a specific marker gene (like PITX2 for the left atrium) plotted without imputation (showing the raw, sparse reality) right next to the plot with MAGIC imputation (showing the smoothed, biological consensus). This demonstrates transparency in your bioinformatics methodology.

### 9.5.1 Analyzing the MAGIC-Imputed Marker Gene UMAPs

By comparing these new imputed UMAPs directly to your previous un-imputed plots, the power of MAGIC (Markov Affinity-based Graph Imputation of Cells) becomes immediately obvious. The "salt-and-pepper" graininess caused by technical scATAC-seq dropouts has been completely completely smoothed out. 

Instead of isolated dots of expression, you now have clear, continuous topological maps of gene activity. Here is how to interpret these imputed visualizations for your analysis.

#### 1. Visualizing Continuous Differentiation (The Right Continuum)
In the un-imputed plots, the developmental trajectories required a bit of imagination to connect the dots. With MAGIC imputation, the flow of differentiation is mathematically smoothed and visually undeniable:
* **The Root:** Look at **CD34** (Page 1). The entire top portion of the right-hand continuum is a solid, glowing mass. This clearly defines the progenitor pool.
* **The Branches:** As you move down from the CD34+ region, the continuum splits. 
    * Moving to the left branch, **GATA1** (Page 2) smoothly lights up, definitively marking the transition into the Erythroid lineage.
    * Moving straight down the main trunk, **MPO** (Page 7) lights up the middle section (early myeloid/granulocyte differentiation).
    * Finally, at the very bottom tip of the trunk, **CD14** (Page 6) dominates. 
* **The Takeaway:** Because MAGIC shares information between nearest neighbors, the transition from CD34 $\rightarrow$ MPO $\rightarrow$ CD14 is rendered as a perfect, smooth color gradient. This visually proves that these cells are not discrete clusters, but rather a single continuous population caught in various stages of biological maturation.

#### 2. Solidifying Discrete Populations (The "Islands")
For the fully mature, circulating immune cells, MAGIC imputation confirms their high degree of homogeneity.
* **The B-Cell Island:** Look at **PAX5**, **MS4A1**, and **MME** (Pages 3, 4, and 5). Previously, these markers showed patchy expression in the top detached island. Now, the entire island is uniformly illuminated. This confirms that practically every cell in that cluster shares a unified B-cell regulatory program.
* **The T-Cell Island:** Look at **CD3D** (Page 8). The entire left island is deeply, uniformly enriched for this pan-T-cell marker. 

#### 3. Revealing Sub-Cluster Structure
Even within a solid island, MAGIC helps define sub-populations without the distraction of dropout noise.
* Look at **CD8A** (Page 9) compared to **CD3D** (Page 8). While the whole left island is CD3D+, only the bottom-right "lobe" of that island is warmly colored for CD8A. Imputation makes the boundary between the CD8+ (cytotoxic) and CD8- (likely CD4+ helper) T-cells remarkably crisp.

> **Thesis Tip:** In your final manuscript, these MAGIC-imputed feature plots are exactly what you want to use for your primary figures. They are aesthetically superior and instantly convey your biological narrative to the reader. You can write: *"To account for the inherent sparsity of single-cell chromatin accessibility data, Gene Scores were visually smoothed using MAGIC imputation. Feature plots of the imputed scores resolved highly continuous developmental gradients—such as the gradual acquisition of CD14 along the myeloid axis—and confirmed the homogeneous identity of discrete lymphoid clusters."*

## 9.6: Module Scores

![alt text](image-30.png)

### Interpreting Gene Module Scores

While looking at individual marker genes (like *CD34* or *CD3D*) is highly informative, relying on a single gene can sometimes be risky due to technical dropouts or biological noise. **Module Scoring** provides a highly robust alternative by calculating the aggregate accessibility of an *entire list* of genes associated with a specific biological state or cell type.

By evaluating a whole "module" or "signature" of genes simultaneously, ArchR mathematically smooths out the noise, providing a much higher-confidence prediction of cell identity.

#### 1. The B-Cell Module Score (`Module.BScore`)
The first plot visualizes the `Module.BScore` plotted across your `UMAP Dimension 1` and `UMAP Dimension 2`. 
* **The Visual Signature:** The bright red and yellow cells (indicating high values) are perfectly restricted to the small, detached island at the top-center of the UMAP. The entire rest of the map is dark blue (low values).
* **The Interpretation:** Because this score aggregates multiple B-cell-specific genes into a single metric, this plot serves as absolute, undeniable confirmation that this specific island represents the B-cell lineage. The lack of "background noise" in the other clusters demonstrates the high specificity of module scoring.

#### 2. The T-Cell Module Score (`Module.TScore`)
The second plot visualizes the `Module.TScore` across the exact same spatial embedding.
* **The Visual Signature:** The bright red and yellow cells have shifted entirely to the large, detached island on the far left side of the UMAP. 
* **The Interpretation:** This confirms the left island is your T-cell population. Notice how the coloration is fairly uniform across that entire left island. Earlier, when you plotted *CD8A* alone, it only lit up half of this island. By using a broader *pan-T-cell* module, you capture the entire lineage (both CD4+ and CD8+ cells) at once.

#### 3. Application to Your Atrial Fibrillation Thesis
Module scoring will be one of the most powerful tools in your cardiovascular analysis pipeline. You will not be looking at blood cells; you will be looking at millions of cardiomyocytes under stress.

* **Region-Specific Identity:** Instead of relying just on *PITX2* to find your Left Atrial cells, you can curate a "Left Atrium Module" containing 50 known LA-specific genes. Plotting this module score will definitively isolate your LA cells from your RV or LV cells.
* **Disease State Scoring:** You can create an "Atrial Fibrillation Stress Module" using known disease-associated genes (e.g., fibrosis markers, fetal gene program markers like *NPPA/NPPB*, or altered ion channels). By plotting this module score on your UMAP, you can visually identify which specific sub-clusters of cardiomyocytes are experiencing the most severe epigenetic remodeling during the disease.

## 9.7 Track Plotting with ArchR Browser

These track plots are the bread and butter of scATAC-seq analysis in ArchR. They visualize the chromatin accessibility landscape (peaks) across your different cell clusters (C1–C14) at specific genomic loci. 

To interpret these tracks **unbiasedly**, you need to separate the visual data from your biological expectations. Instead of looking at a cluster and trying to force a cell-type label onto it, you should build a matrix of accessibility signatures and let the data tell you what each cluster represents.


### 1. Orient Yourself to the Plot Anatomy
Before making any biological assumptions, understand exactly what the plots are showing:
* **Y-Axis (Signal):** This represents the normalized ATAC-seq signal. Higher peaks mean the chromatin is more "open" or accessible in that specific cluster. 
* **X-Axis (Coordinates):** These are the genomic coordinates. 
* **Bottom Panel (Gene Models):** This shows the genes in that region. The thickest blue/red bars represent exons, the thin lines are introns, and the flat end of the gene model indicates the Transcription Start Site (TSS) / promoter region.

### 2. Isolate the Promoter/TSS
Chromatin accessibility at the promoter is the strongest indicator of potential gene expression. 
* Locate your marker gene in the bottom panel (e.g., **CD14**, **IL7R**).
* Draw an imaginary vertical line up from the TSS of that gene through all the cluster tracks (C1-C14). 
* **Rule of thumb:** Only consider a cluster "positive" for that marker if there is a distinct, sharp peak directly over the TSS or immediately upstream. Ignore background noise.

### 3. Create a "Blind" Matrix
To remain unbiased, evaluate the tracks without looking at your comments (e.g., `#B-Cell Trajectory`). 
Create a simple table or spreadsheet. On the X-axis, list your clusters (C1-C14). On the Y-axis, list your marker genes (CD34, GATA1, PAX5, MS4A1, etc.).
* Go through the PDF page by page.
* Score each cluster for each gene as `+` (strong peak at TSS), `+/-` (weak/ambiguous peak), or `-` (no peak/flat line).
* *Example using your T-Cell markers:* Evaluate **CD3D**, **CD8A**, **TBX21**, and **IL7R**. Note exactly which clusters show open chromatin at these promoters.

### 4. Look for Concordance (The Unbiased Check)
This is where the unbiased approach pays off. Biological markers rarely act alone. 
* **Lineage Agreement:** If a cluster is truly a B-Cell, it shouldn't just have a peak at **PAX5**; it should ideally also show accessibility at **MS4A1**. If C8 has a huge peak at PAX5 but absolutely nothing at MS4A1, treat C8's identity with skepticism—it might be an intermediate state or a different lineage entirely.
* **Mutual Exclusivity:** Look for expected divergence. A cluster strongly positive for the Monocyte marker **CD14** should generally lack peaks at the T-cell markers like **CD8A**. If a cluster shows strong peaks for *both*, it could indicate a double-cell (multiplet) artifact in your sequencing data rather than a novel biological state.

### 5. Assign Identities
Only after your matrix is fully populated should you reveal the biological categories to yourself. Match your matrix clusters to the expected profiles:
* **Early Progenitor:** Clusters positive only for CD34.
* **T-Cells:** Clusters positive for CD3D, IL7R, etc.
* **Sub-clustering:** Notice differences within a lineage. For example, if Clusters 4, 5, and 6 all have **CD3D** peaks, but only Cluster 6 has a **CD8A** peak, you have just unbiasedly identified your CD8+ T-cell subpopulation versus other T-cells.


https://www.archrproject.com/reference/plotBrowserTrack.html