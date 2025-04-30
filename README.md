# WASP-Indel: Indel sensitive pipeline for unbiased read mapping and molecular QTL discovery

## Introduction

WASP-Indel is based on the software tool WASP suite for unbiased allele-specific read mapping and
[the paper](http://biorxiv.org/content/early/2014/11/07/011221): van de Geijn B\*, McVicker G\*, Gilad Y, Pritchard JK. "WASP: allele-specific software for robust discovery of molecular quantitative trait loci".

The main practical advance of WASP-Indel, over the original WASP, is the ability to determine the impact of indels, as well as SNPs, on molecular traits. WASP-Indel also benefits from a simplification of the required input files by eliminating the requirement of converting the input genome FASTA files and VCF to the HDF5 format. There are additional advances in identifying potentially mis-genotyped variants and in the processing of large number of phenotype/genotype combinations.

WASP-Indel has two parts, which can be used independently of each
other:

1. Read filtering tools that correct for biases in allele-specific
   mapping. 

2. A Combined Haplotype Test (CHT) that tests for genetic association
   with a molecular trait using counts of mapped and allele-specific
   reads.

The following directories and files are included with WASP.
Each directory contains its own README file:

* [mapping](./mapping) - Mappability filtering pipeline for correcting allelic mapping biases

* [CHT_input](./CHT_input) - Generate the input for the Combined Haplotype Test from the unbiased allele specific mapping.

* [CHT](./CHT) - Code for running the Combined Haplotype Test


This is a new repo to adjust WASP_indel software to ASAP project
Origual WASP-Indel repository: https://github.com/adam-rabinowitz/wasp_indel_v1.0

Citations:
Sigalova, O. M., Forneris, M., Stojanovska, F., Zhao, B., Viales, R. R., Rabinowitz, A., Hammal, F., Ballester, B., Zaugg, J. B., & Furlong, E. E. M. (2025). Integrating genetic variation with deep learning provides context for variants impacting transcription factor binding during embryogenesis. Genome Research, 35, 1–16. https://doi.org/10.1101/GR.279652.124
Van De Geijn, B., Mcvicker, G., Gilad, Y., & Pritchard, J. K. (2015). WASP: Allele-specific software for robust molecular quantitative trait locus discovery. Nature Methods, 12(11), 1061–1063. https://doi.org/10.1038/nmeth.3582


