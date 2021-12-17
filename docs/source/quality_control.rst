***************
Quality Control
***************

Next generation sequencing inherently brings about techinal bias on sequenced reads. It is essential to check the quality of sequencing dataset in some aspects.

Illumina
########

For the each fastq of paired-end dataset, `FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_  is performed. Two fastqc results are aggregated into one report with `Multiqc <https://multiqc.info/>`_.

Nanopore
########

`NanoPlot <https://github.com/wdecoster/NanoPlot>`_ is used.