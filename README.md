
# Comparative lncRNA-seq Analysis: Diploid vs. Tetraploid *Brassica oleracea*

This repository provides a detailed pipeline to compare long non-coding RNA (lncRNA) sequencing data between diploid and tetraploid samples of *Brassica oleracea*. 

## Prerequisites

Install the following tools:

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (quality control)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) or [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) (read trimming)
- [HISAT2](https://daehwankimlab.github.io/hisat2/) (alignment)
- [StringTie](https://ccb.jhu.edu/software/stringtie/) (optional: transcript assembly)
- [featureCounts](http://subread.sourceforge.net/) (read counting)
- [Samtools](http://www.htslib.org/) (BAM file handling)
- [R](https://www.r-project.org/) with [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) (differential expression analysis)
- [MultiQC](https://multiqc.info/) (optional: aggregate QC reports)

## Input Files

- Diploid FASTQ: `pC47-2_good_1.fq.gz`, `pC47-2_good_2.fq.gz`
- Tetraploid FASTQ: `LncRNA_TC47-2_good_1.fq.gz`, `LncRNA_TC47-2_good_2.fq.gz`
- Reference genome: `Brassica_oleracea_JZS_v2.fasta`
- lncRNA annotation: `Brassica_oleracea_lncRNAs.gtf`
- (Optional) lncRNA sequences for validation: `Brassica_oleracea_lncRNAs.fasta`

## Pipeline

### 1. Quality Control

```bash
mkdir fastqc_results

fastqc -o fastqc_results pC47-2_good_1.fq.gz pC47-2_good_2.fq.gz   LncRNA_TC47-2_good_1.fq.gz LncRNA_TC47-2_good_2.fq.gz

# (Optional) Aggregate FastQC reports
multiqc fastqc_results -o multiqc_results
```

**Output:** HTML reports in `fastqc_results/` and `multiqc_results/`.

---

### 2. Read Trimming

```bash
mkdir trimmed_reads

for sample in pC47-2 LncRNA_TC47-2; do
  trimmomatic PE -phred33     ${sample}_good_1.fq.gz ${sample}_good_2.fq.gz     trimmed_reads/${sample}_1_paired.fq.gz trimmed_reads/${sample}_1_unpaired.fq.gz     trimmed_reads/${sample}_2_paired.fq.gz trimmed_reads/${sample}_2_unpaired.fq.gz     ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
```

Use the `_paired.fq.gz` files for downstream analysis.

---

### 3. Index the Reference Genome

```bash
mkdir genome_index

hisat2-build Brassica_oleracea_JZS_v2.fasta genome_index/Brassica_oleracea_JZS_v2
```

---

### 4. Align Reads

```bash
mkdir alignments

for sample in pC47-2 LncRNA_TC47-2; do
  hisat2 -p 8 --dta -x genome_index/Brassica_oleracea_JZS_v2     -1 trimmed_reads/${sample}_1_paired.fq.gz     -2 trimmed_reads/${sample}_2_paired.fq.gz     -S alignments/${sample}.sam
  
  samtools view -bS alignments/${sample}.sam | samtools sort -o alignments/${sample}.sorted.bam
  samtools index alignments/${sample}.sorted.bam
  rm alignments/${sample}.sam
done
```

---

### 5. Quantify lncRNA Expression

```bash
mkdir quantification

featureCounts -p -t exon -g gene_id   -a Brassica_oleracea_lncRNAs.gtf   -o quantification/lncRNA_counts.txt   alignments/pC47-2.sorted.bam alignments/LncRNA_TC47-2.sorted.bam
```

---

### 6. Differential Expression Analysis (DESeq2)

```r
# Install DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

library(DESeq2)

# Read count matrix
counts <- read.table("quantification/lncRNA_counts.txt", header=TRUE, row.names=1, skip=1)
counts <- counts[, 6:ncol(counts)]
colnames(counts) <- c("pC47-2", "LncRNA_TC47-2")

# Metadata
colData <- data.frame(
  sample = c("pC47-2", "LncRNA_TC47-2"),
  ploidy = c("Diploid", "Tetraploid"),
  row.names = c("pC47-2", "LncRNA_TC47-2")
)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ ploidy)

dds <- DESeq(dds)
res <- results(dds, contrast=c("ploidy", "Tetraploid", "Diploid"))

write.csv(as.data.frame(res), "differential_expression_results.csv")

sig_res <- subset(res, padj < 0.05)
write.csv(as.data.frame(sig_res), "significant_differential_expression.csv")
```

---

### 7. Visualization (Optional)

**Volcano Plot:**

```r
library(ggplot2)
res_df <- as.data.frame(res)
res_df$significant <- res_df$padj < 0.05 & !is.na(res_df$padj)

ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
  geom_point() +
  theme_minimal() +
  labs(title="Volcano Plot: Tetraploid vs Diploid", x="Log2 Fold Change", y="-Log10 Adjusted P-value")

ggsave("volcano_plot.png")
```

**Heatmap:**

```r
library(pheatmap)
sig_counts <- counts[rownames(sig_res), ]
pheatmap(log2(sig_counts + 1), scale="row", main="Heatmap of Significant lncRNAs")
```

---

### 8. Functional Annotation (Optional)

- BLAST lncRNA sequences against known databases (e.g., NONCODE, LncBook).
- Perform co-expression analysis (e.g., using WGCNA) to find related protein-coding genes.
- Explore nearby genes based on lncRNA locations from the GTF file.

---

## Example Directory Structure

```
project/
├── fastqc_results/
├── trimmed_reads/
├── genome_index/
├── alignments/
├── quantification/
├── differential_expression_results.csv
├── significant_differential_expression.csv
├── volcano_plot.png
├── Brassica_oleracea_JZS_v2.fasta
├── Brassica_oleracea_lncRNAs.gtf
├── pC47-2_good_1.fq.gz
├── pC47-2_good_2.fq.gz
├── LncRNA_TC47-2_good_1.fq.gz
├── LncRNA_TC47-2_good_2.fq.gz
```

---

## Notes

- **Resources:** RNA-seq analyses are computationally intensive. Use high-performance computing if needed.
- **Replicates:** This pipeline assumes one diploid and one tetraploid sample. Modify the DESeq2 input if more replicates are available.
- **Adjust Parameters:** You may need to tune trimming and alignment parameters based on your data.
- **Validation:** To discover novel lncRNAs, use StringTie and compare assembled transcripts to the lncRNA annotation.

---

## Contact

Feel free to open an issue or pull request if you have questions or suggestions!
