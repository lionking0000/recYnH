rec-YnH
==========

Table of Contents:

- [Summary](#summary)
- [Requirements](#requirements)
- [Installation](#installation)
	- [rec-YnH](#vast-tools-1)
- [Usage](#usage)
	- [Help](#help)
	- [Quick Usage](#quick-usage)
	- [Alignment](#alignment)
	- [Merging Outputs](#merging-outputs)
	- [Strand-specific RNAseq data](#strand-specific-rnaseq-data)
	- [Combining Results](#combining-results)
	- [Comparing PSIs Between Samples](#comparing-psis-between-samples)
	- [Differential Splicing Analysis](#differential-splicing-analysis)
	- [Plotting](#plotting)
	- [Simplifying Combine Table](#simplifying-combine-table)
- [Issues](#issues)
- [Contributions](#contributions)
- [Citation](#citation)
- [References](#references)
	
Summary
-------
This program is taking recYnH sequencing files and generating recYnH interaction score matrix to correponding genes.

Requirements
------------

rec-YnH requires the following software:
 * python 2.7 or higher, with the following packages installed (see Installation Section):
   * scipy, https://pypi.org/project/scipy/
 * cutadapt, http://cutadapt.readthedocs.io/en/stable/installation.html
 * blast+ 2.6.0 or higher, http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+
 
 * R 3.1 or higher, with the following packages installed (see Installation Section):
   * optparse
   * RColorBrewer
   * reshape2
   * ggplot2 >= v2.0
   * MASS
   * devtools
   * [psiplot](https://github.com/kcha/psiplot)
 
 * bowtie 1.0.0 (Langmead et al., 2009), http://bowtie-bio.sourceforge.net/index.shtml
 
Installation
------------

~~~~
> git clone https://github.com/lionking0000/recYnH.git
~~~~

R packages

RUN Rscript -e "install.packages(c('pheatmap','RColorBrewer'))"
RUN Rscript -e "install.packages(c('gtools','reshape2'))"
RUN Rscript -e "install.packages(c('psych','clipr'))"
RUN Rscript -e "install.packages(c('swfscMisc','PerformanceAnalytics'))"
RUN Rscript -e "install.packages(c('mixtools','pROC'))"
RUN Rscript -e "install.packages(c('outliers','readxl'))"
RUN Rscript -e "install.packages(c('d3heatmap','matrixStats'))"

Alignment
------------

Here is the example of how to do alignment (mapping) FASTQ files into reference sequences (FASTA).

Example 1
~~~~
python recYnH.py align -i1 ./example/db/A463-MGj69.RBP-MAP.-150.fa -f1 ./example/fastq/S1_WD_R1.300000.fastq.gz -f2 ./example/fastq/S1_WD_R2.300000.fastq.gz -o ./example/output/exp1 -n S1_WD
~~~~
~~~~
python recYnH.py align -i1 ./example/db/A463-MGj69.RBP-MAP.-150.fa -f1 ./example/fastq/S2_QD_R1.300000.fastq.gz -f2 ./example/fastq/S2_QD_R2.300000.fastq.gz -o ./example/output/exp1 -n S2_QD
~~~~

Example 2
~~~~
python recYnH.py align -i1 ./example/db/roth2016_control_set_plus_control.-150.fa -f1 ./example/fastq/EXP1_W_R1.fastq -f2 ./example/fastq/EXP1_W_R2.fastq -o ./example/output/ -n EXP1_W
~~~~

~~~~
python recYnH.py align -i1 ./example/db/roth2016_control_set_plus_control.-150.fa -f1 ./example/fastq/EXP1_Q_R1.fastq -f2 ./example/fastq/EXP1_Q_R2.fastq -o ./example/output/ -n EXP1_Q
~~~~


Merge two matrix files to generate Interaction Score matrix
------------

~~~~
python recYnH.py merge -m1 ./example/output/EXP1_W -m2 ./example/output/EXP1_Q -o ./example/output/ -n EXP1
~~~~

zcat ../../../2017-11-03_MiSeq/S1_W_R1.fastq.gz  | head -n 400000 > EXP1_W_R1.fastq