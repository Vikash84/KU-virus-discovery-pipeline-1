**************
Read filtering
**************

By filtering some noisy reads, following analysis results in stronger signal of truly existing viruses as well as it becomes faster.

Illumina
########

Along the PCR step, some genomic sequences produce multiple copies of reads. This can lead to wrong interpreation. The pipeline gets rid of duplicated reads with `Prinseq <http://prinseq.sourceforge.net/>`_.

Too short and low quality reads and primer sequence give noise to true viral sequences. We can remove them with `Trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_.

Default: `ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`

The reads originating from host genome also should be removed. For Illumina reads, `bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_ is employed.

Nanopore
########

`Nanofilt <https://github.com/wdecoster/nanofilt>`_ is used for the length filtration. You can give read qulity constraint by modifying `nextflow.config`.

Default: `--length 300 --readtype 1D`

`Minimap2 <https://github.com/lh3/minimap2>`_ aligns Nanopore reads onto host genome.