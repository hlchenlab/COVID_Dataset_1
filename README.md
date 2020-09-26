# COVID_Dataset_1

## Overview: 
The devastating impacts of the 2020 COVID-19 pandemic has demonstrated the power of viral infection in susceptible hosts. As such, there is a growing importance and urgency to understand the dynamics of SARS-CoV-2 infection induced changes in host gene expression.
Empirical evidence has shown that Golden Hamster (Mesocricetus auratus) has been proven to be a suitable host model to assess SARS-CoV-2 kinetics and as such was choosen as the live model for our study. The experiment design and details of sequencing can be found online along with the raw fastq sequencing data on the Gene Expression Omnibus database (https://www.ncbi.nlm.nih.gov/geo/) with the GEO accession number: GSE156005 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156005).

## Description of Contents:

Each folder possesses the relevant bash / or Rscript and the relevant metadata used in the original analysis of the hamster data. 

Alignment - Performed using STAR (V2.7.2a) (Dobin et al., 2013). Index generation was buit using the combined gemeones of Golden Hamster (NCBI ID:GCA_000215625.1) and the SARS-CoV-2 genomes (RefSeq ID; NC_045512.2). A combined GTF and FASTA file are in this directory for this purpose. 

Differential Expression - DESeq2 (Love et al., 2014) was used. R code for volcano plots and PCA plots is also included.

Gene Enrichment Analysis- The recent version of Gene Ontology (download on January 2020) was used for all gene enrichment analysis.

## References:

DOBIN, A., DAVIS, C. A., SCHLESINGER, F., DRENKOW, J., ZALESKI, C., JHA, S., BATUT, P., CHAISSON, M. & GINGERAS, T. R. 2013. STAR: ultrafast universal RNA-seq aligner. Bioinformatics (Oxford, England), 29, 15-21.

LOVE, M. I., HUBER, W. & ANDERS, S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15, 550.
