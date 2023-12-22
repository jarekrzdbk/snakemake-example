# Example of Snakemake Workflow Implementation for Transcriptomics Tasks

## Github Pipeline:
**Build Status**  
![Snakemake CI](https://github.com/jarekrzdbk/snakemake-homework/actions/workflows/main.yml/badge.svg)

## Implemented Tasks:

1. **Performed FastQC**

2. **Performed MultiQC**

3. **Trimmed Barcodes**

4. **Performed FastQC on Trimmed Sequences**

5. **Performed MultiQC on Trimmed Sequences**

The biggest indicator of successful removal of adapters is the Adapter Content plot:

| Before Trimming                                  | After Trimming                                   |
|:------------------------------------------------:|:------------------------------------------------:|
| ![notTrimmed](qc_plots/fastqc_adapter_content_plot.png)  |  ![trimmed](qc_plots/fastqc_adapter_content_plot_trimmed.png) |

6. **Used STAR Aligner**

7. **Indexed BAM with Samtools Index, and FeatureCounts**

Here we specify strandness:  
0: Not stranded  
1: Stranded  
2: Reversely stranded

For Collibri, use 1 as it is stranded.  
For KAPA, use 2 as it is reversely stranded.

8. **Used DESeq2 to Perform Differential Expression (DE) Analysis Comparing UHRR vs HBR**

As a result, obtained DE genes with p-adjusted values:

For Collibri:  
![collibri](/condition_treated_results.csv)

And KAPA:  
![KAPA](/condition_treated_results_kapa.csv)

Also, volcano plots:

For Collibri:  
![collibri-volcano](/deseq2_files/figure-gfm/collibri-volcano-1.png)

And KAPA:  
![kapa-volcano](/deseq2_files/figure-gfm/kapa-volcano-1.png)

It has a Snakemake file, which will create the plot.

9. **Performed PCA Using DE Genes**

Plot for Collibri:  
![collibri-pca](/deseq2_files/figure-gfm/collibri-pca-1.png)

And KAPA:  
![kapa-pca](/deseq2_files/figure-gfm/kapa-pca-1.png)

It is clear that HBR and UHRR samples are separated into two distinct clusters.

More details and figures in separate `deseq2.md`.  
![deseq2.md](/deseq2.md)

10. **Performed GSEA Using fgsea**

We should use shrunken fold changes.

Performing on Reactome pathways did not give promising results, so performed on KEGG pathways.

Collibri showed slightly better results (with smaller p-adjusted values) than KAPA but found similar pathways. Some pathways were related to cancer:

- KEGG_BASAL_CELL_CARCINOMA
- KEGG_COLORECTAL_CANCER
- KEGG_ENDOMETRIAL_CANCER
- KEGG_WNT_SIGNALING_PATHWAY
- KEGG_PATHWAYS_IN_CANCER
- KEGG_REGULATION_OF_ACTIN_CYTOSKELETON

For Collibri:  
![collibri-pathways](/gsea/collibri_pathways.png)

For KAPA:  
![kapa-pathways](/gsea/kapa-pathways.png)

## DAG of Snakemake Workflow:

![dag](/dag.svg)

