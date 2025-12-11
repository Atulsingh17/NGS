#!/bin/bash
# Installation script for WGS data analysis dependencies
# Author: Atul_singh
# Run as: sudo bash WGS_analysis_dependencies_installation.sh

echo "Updating system..."
sudo apt update && sudo apt upgrade -y

echo "Installing system dependencies..."
sudo apt install -y openjdk-11-jre-headless wget curl git build-essential zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev libncursesw5-dev samtools bcftools bwa fastqc trimmomatic

echo "Installing Miniconda..."
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda
source $HOME/miniconda/bin/activate
conda init bash
source ~/.bashrc

echo "Creating WGS environment..."
conda create -n wgs_env python=3.10 -y
conda activate wgs_env

echo "Installing Python packages via conda/pip..."
conda install -y -c bioconda gatk4 picard bwa samtools bcftools fastqc multiqc bedtools vcftools
pip install pysam pyvcf

echo "Installing GATK (if not via conda)..."
wget https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip
unzip gatk-4.5.0.0.zip
sudo mv gatk-4.5.0.0 /usr/local/gatk

echo "Adding paths to .bashrc..."
echo 'export PATH=$HOME/miniconda/bin:/usr/local/gatk:$PATH' >> ~/.bashrc
source ~/.bashrc

echo "Testing installations..."
fastqc --version
bwa
samtools --version
gatk --version

echo "WGS dependencies installed! Activate with: conda activate wgs_env"
