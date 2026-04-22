# Explanation of some points in ArchR-Guide


Quality Control (QC) in ArchR

In scATAC-seq, we don't just look at "how many reads" a cell has. We look at Signal vs. Noise.
Two Main Metrics:

    TSS Enrichment Score: This measures how many fragments fall near Transcription Start Sites (TSS) compared to the "background" genomic noise. High TSS enrichment = high-quality data.

    Fragment Size Distribution: Because DNA is wrapped around nucleosomes, a good ATAC-seq library should show a distinct "pattern" (peaks representing fragments that are sub-nucleosomal, mono-nucleosomal, etc.).


What exactly are Arrow files?

![Erstes Bild](image.png)

