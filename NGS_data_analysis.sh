#!/bin/bash
# COMPLETE NGS ANALYSIS IN ONE SCRIPT (2025)
# WGS + RNA-Seq + miRNA-Seq → Final Figures & Tables
# Just run: bash NGS_data_analysis.sh

set -e  # Stop on any error

# Activate environment
source ~/miniconda/bin/activate ngs_master

# INPUT (change these)
RAW_DIR="raw_fastq"           # Your .fastq.gz files here
SAMPLE_SHEET="samples.csv"    # columns: sample,condition,replicate,type (dna/rna/mirna)
GENOME="genome/hg38.fa"       # Reference genome
GTF="genome/gencode.v44.gtf"  # For RNA-Seq
MIRBASE="mirbase/hsa.mature.fa"  # For miRNA

mkdir -p results/{qc,trimmed,aligned,counts,variants,mirna,de,figures}

echo "Step 1: Quality Check"
fastqc $RAW_DIR/*.fastq.gz -o results/qc/ -t 20
multiqc results/qc/ -o results/qc/

echo "Step 2: Trimming"
for f in $RAW_DIR/*_R1*.fastq.gz; do
    base=$(basename $f _R1.fastq.gz)
    trimmomatic PE -threads 16 \
        $RAW_DIR/${base}_R1.fastq.gz $RAW_DIR/${base}_R2.fastq.gz \
        results/trimmed/${base}_R1_trim.fastq.gz results/trimmed/${base}_R1_unpaired.fastq.gz \
        results/trimmed/${base}_R2_trim.fastq.gz results/trimmed/${base}_R2_unpaired.fastq.gz \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

echo "Step 3: Alignment & Processing"
while IFS=, read -r sample condition rep type; do
    if [[ $type == "dna" ]]; then
        # WGS Pipeline
        bwa mem -t 16 $GENOME results/trimmed/${sample}_R1_trim.fastq.gz results/trimmed/${sample}_R2_trim.fastq.gz | \
            samtools sort -@ 16 -o results/aligned/${sample}.bam -
        samtools index results/aligned/${sample}.bam
        gatk MarkDuplicates -I results/aligned/${sample}.bam -O results/aligned/${sample}.dedup.bam -M metrics.txt
        gatk HaplotypeCaller -R $GENOME -I results/aligned/${sample}.dedup.bam -O results/variants/${sample}.vcf.gz

    elif [[ $type == "rna" ]]; then
        # RNA-Seq Pipeline
        hisat2 -p 16 --dta -x $GENOME -1 results/trimmed/${sample}_R1_trim.fastq.gz -2 results/trimmed/${sample}_R2_trim.fastq.gz | \
            samtools sort -@ 16 -o results/aligned/${sample}.bam -
        featureCounts -T 16 -p -t exon -g gene_id -a $GTF -o results/counts/${sample}_counts.txt results/aligned/${sample}.bam

    elif [[ $type == "mirna" ]]; then
        # miRNA-Seq Pipeline
        cutadapt -a TGGAATTCTCGGGTGCCAAGG -O 6 --minimum-length 18 --maximum-length 30 \
            -o results/trimmed/${sample}_trimmed.fastq.gz $RAW_DIR/${sample}.fastq.gz
        mirdeep2.pl results/trimmed/${sample}_trimmed.fastq.gz $GENOME results/trimmed/*.fastq.gz \
            $MIRBASE none none -t Human -v -P results/mirna/${sample}
    fi
done < <(tail -n +2 $SAMPLE_SHEET)

echo "Step 4: Differential Expression (RNA-Seq)"
Rscript -e "
library(DESeq2)
countFiles <- list.files('results/counts', pattern='_counts.txt', full.names=T)
counts <- sapply(countFiles, function(x) read.table(x, header=T, row.names=1)[,6])
colnames(counts) <- basename(countFiles) %>% gsub('_counts.txt','',.)
coldata <- read.csv('samples.csv')
coldata <- coldata[coldata\$type=='rna',]
dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)
dds <- DESeq(dds)
res <- results(dds, contrast=c('condition','Tumor','Normal'))
write.csv(as.data.frame(res), 'results/de/deseq2_results.csv')
png('results/figures/volcano.png'); plotMA(res); dev.off()
"

echo "Step 5: Final MultiQC Report"
multiqc results/ -o results/figures/

echo "ANALYSIS COMPLETE!"
echo "Results → results/ folder"
echo "Figures → results/figures/"
echo "Final Report → results/figures/multiqc_report.html"
