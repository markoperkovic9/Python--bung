# Structured Learning Plan: Biological Foundations for Cardiac scATAC-seq Analysis

This plan is designed to bridge your computational skills in ArchR with the deep biological context required for your thesis. It is broken down into four sequential phases, moving from macro-anatomy down to molecular genetics.

## Phase 1: The Healthy Heart & Region-Specific Cardiomyocyte Biology
*Goal: Understand why a left atrial cardiomyocyte has a fundamentally different epigenome than a left ventricular cardiomyocyte.*

### 1. Developmental Origins
* **Concepts to Master:** First Heart Field (FHF) vs. Second Heart Field (SHF). 
* **Why it matters:** The Left Ventricle is primarily FHF, while the Right Ventricle and Atria are heavily SHF. These embryological origins leave permanent "epigenetic scars" (accessible regions) that you will see in your data.
* **Action Items:**
    * [ ] Read a review on mammalian heart development focusing on cardiac progenitor fields.
    * [ ] Map out which anatomical structures derive from the FHF vs. SHF.

### 2. Atrial vs. Ventricular Physiology
* **Concepts to Master:** Electrophysiology, contractile force, and metabolism.
* **Why it matters:** Atrial cells are specialized for fast conduction and low pressure; ventricular cells are specialized for high pressure and slow fatigue. Their open chromatin will reflect these functional differences (e.g., different ion channels).
* **Action Items:**
    * [ ] Review the cardiac action potential differences between atrial and ventricular myocytes.
    * [ ] Identify the classic transcriptomic markers for these regions (e.g., *MYL4*, *NPPA*, *PITX2* for Atria; *MYL2*, *MYH7* for Ventricles).

## Phase 2: The Disease Context - Atrial Fibrillation (AF)
*Goal: Understand the pathology of AF to know what abnormal signals to look for in your atrial scATAC-seq data.*

### 1. Pathophysiology of AF
* **Concepts to Master:** Electrical remodeling (ion channel dysfunction) and Structural remodeling (fibrosis, hypertrophy, atrial enlargement).
* **Why it matters:** AF causes massive transcriptomic and epigenetic shifts. You need to know if the changes you observe in your data are driving the disease or are secondary responses to the disease.
* **Action Items:**
    * [ ] Read clinical reviews on the mechanisms of Atrial Fibrillation.
    * [ ] Study the role of fibroblasts vs. cardiomyocytes in AF-induced structural remodeling (ensure your CM clusters are highly purified).

### 2. The Genetics of AF
* **Concepts to Master:** The 4q25 locus and *PITX2*.
* **Why it matters:** *PITX2* is the master regulator of left-atrial identity and the most famous risk locus for AF. 
* **Action Items:**
    * [ ] Read the landmark papers linking *PITX2* non-coding mutations to AFib risk.
    * [ ] Explore the GWAS Catalog for "Atrial Fibrillation" to see the top associated genetic loci.

## Phase 3: Cardiac Epigenetics & Gene Regulation
*Goal: Translate biological changes into chromatin accessibility concepts.*

### 1. Cis-Regulatory Elements (CREs)
* **Concepts to Master:** Promoters vs. Enhancers.
* **Why it matters:** Promoters are usually open across all cell types, but Enhancers are highly specific to cell state, region, and disease. Most of your region-specific and AF-specific discoveries will happen at Enhancers.
* **Action Items:**
    * [ ] Review the mechanisms of enhancer-promoter looping.
    * [ ] Learn how scATAC-seq defines a "peak" and how that relates to an active enhancer.

### 2. Core Cardiac Transcription Factors (TFs)
* **Concepts to Master:** The cardiac TF network (GATA4, MEF2C, NKX2-5, TBX5, HAND2, PITX2).
* **Why it matters:** TFs bind to the open DNA you are measuring. If a region of chromatin closes during AF, it is usually because a TF stopped binding there.
* **Action Items:**
    * [ ] Look up the DNA binding motifs for these core cardiac TFs (using databases like JASPAR).
    * [ ] Understand the concept of "Pioneer Factors" that can open closed chromatin.

## Phase 4: Statistical Genetics & Variant Effect Prediction (VEP)
*Goal: Connect the GWAS risk variants to your scATAC-seq peaks.*

### 1. From GWAS to Fine-Mapping
* **Concepts to Master:** Linkage Disequilibrium (LD) and causal variants.
* **Why it matters:** GWAS provides a "block" of correlated SNPs. Your job is to overlay these SNPs onto your scATAC-seq peaks to find the *one* SNP that actually sits inside an open enhancer in an atrial cardiomyocyte.
* **Action Items:**
    * [ ] Learn the basics of Linkage Disequilibrium.
    * [ ] Read papers on "epigenomic fine-mapping" of complex traits.

### 2. Variant Effect Prediction Mechanisms
* **Concepts to Master:** chromatin accessibility Quantitative Trait Loci (caQTLs) and motif disruption.
* **Why it matters:** The core hypothesis of your thesis: A genetic mutation (SNP) changes a DNA letter -> The TF can no longer bind -> The enhancer closes -> The target gene turns off -> AF develops.
* **Action Items:**
    * [ ] Read about how single nucleotide variants disrupt TF binding motifs.
    * [ ] Familiarize yourself with the concepts behind deep learning models used for VEP (like BPNet, Enformer, or ChromBPNet) which predict how a specific DNA sequence change will alter the ATAC-seq signal.