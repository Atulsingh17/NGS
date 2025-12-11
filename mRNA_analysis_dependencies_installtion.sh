#!/bin/bash
# Installation script for mRNA sequencing data analysis dependencies
# Author: Atul_singh
# Run as: sudo bash mRNA_analysis_dependencies_installtion.sh

echo "Updating system..."
sudo apt update && sudo apt upgrade -y

echo "Installing system dependencies..."
sudo apt install -y openjdk-11-jre-headless wget curl git build-essential zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev libncursesw5-dev r-base r-base-dev fastqc trimmomatic hisat2

echo "Installing Miniconda..."
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda
source $HOME/miniconda/bin/activate
conda init bash
source ~/.bashrc

echo "Creating mRNA environment..."
conda create -n mrna_env python=3.10 -y
conda activate mrna_env

echo "Installing Python packages via conda/pip..."
conda install -y -c bioconda star hisat2 featurecounts htseq subread multiqc bedtools
pip install pysam rpy2

echo "Installing R packages for DE analysis (DESeq2, edgeR, etc.)..."
R -e "install.packages('BiocManager', repos='http://cran.us.r-project.org')"
R -e "BiocManager::install(c('DESeq2', 'edgeR', 'limma', 'ggplot2', 'pheatmap', 'clusterProfiler'))"

echo "Installing STAR (if not via conda)..."
wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.11b.tar.gz
tar xzf 2.7.11b.tar.gz
cd STAR-2.7.11b/source
make STAR
sudo mv STAR /usr/local/bin/

echo "Adding paths to .bashrc..."
echo 'export PATH=$HOME/miniconda/bin:/usr/local/bin:$PATH' >> ~/.bashrc
source ~/.bashrc

echo "Testing installations..."
fastqc --version
hisat2 --version
STAR --version
featureCounts -v
R -e "library(DESeq2); print('DESeq2 installed!')"

echo "mRNA dependencies installed! Activate with: conda activate mrna_env"
