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

* **Thesis Tip:** Consistency is key for your final figures. Once you find a `baseSize` and `alpha` value that looks good in your VS Code preview, use those same values for every ridge and violin plot in your thesis to give your document a cohesive, professional look.