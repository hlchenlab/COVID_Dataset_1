# COVID_Dataset_1

## Overview: 
The devastating impacts of the 2020 COVID-19 pandemic has demonstrated the power of viral pathogencity in susceptible hosts. As such, there is a growing importance and urgency to understand the dynamics of SARS-CoV-2 infection induced changes in host gene expression.
Empirical evidence has shown that Golden Hamster (Mesocricetus auratus) has been proven to be a suitable host model to assess SARS-CoV-2 kinetics and as such was choosen as the live model for our study. The raw fastq sequencing data is available on the Gene Expression Omnibus database (https://www.ncbi.nlm.nih.gov/geo/) with the GEO accession number: GSE156005 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156005).

### Whole genome sequencing analysis of SARS-CoV-2
Consensus genomes were constructed based on the modified Artic Network nCoV-2019 novel coronavirus bioinformatics protocol [Version 1] (https://artic.network/ncov-2019). Basecalled sequencing reads were demultiplexed with Porechop v0.2.4, followed by the mapping to the SARS-CoV-2 reference genome Wuhan-Hu-1 reference genome (Accession No.: NC_045512) using BWA v0.7.17 and Samtools v1.7. The primer sequences were then removed with the align-trim module. The BAM file was used for variant calling by Medaka v0.11.5 with threshold value set at 1 to ensure haploid decoding. The consensus genome was assembled based on the VCF from Medaka v0.11.5 and BAM file. Depth per base required for variant calling was set at >20X with a minimum required Phred score of 10.

## Description of Contents:

Each folder possesses the relevant bash / or Rscript and the relevant metadata/accessory files used in the original analysis of this dataset. 

### Alignment 
Performed using STAR (V2.7.2a) (Dobin et al., 2013). Index generation was buit using the combined geneomnes of Golden Hamster (NCBI ID:GCA_000215625.1) and the SARS-CoV-2 genomes (NCBI ID: NC_045512.2)(https://www.ncbi.nlm.nih.gov/nuccore/1798174254 for fasta file). Hamster (Mesocricetus auratus) genome and annotation were downloaded from Ensembl (version 100) (ftp://ftp.ensembl.org/pub/release-100/). A gtf file based on the NC_045512.2 entry has been included in this directory. Annotations and FASTA files for both were subsequently combined and used to generate the STAR index. 

### Differential Expression 
DESeq2 (Love et al., 2014) was used. R code for differential analysis, volcano plots and PCA plots is also included. A metadata file has been provided to give an accurate description of samples with a raw-count dataframe derived from STAR output. We also provide a table of SARS-CoV-2 gene names to enable conversion of gene-ids for analysis if needed.

### Gene Enrichment Analysis
A recent version of Gene Ontology (http://geneontology.org/docs/download-ontology, downloaded January 2020) was used for all gene enrichment analysis. We have included an R-compatible GO annotation file to fit the use for our R code if needed. An Entrez-to-Ensembl gene id conversion table (R-compatible) was built using GENCODE (v29) gene annotation (https://www.gencodegenes.org/human/release_29.html) is also available in this directory.

## Other analysis relating to this dataset:

### Phylogenetic analysis
 
Consensus genomes were aligned by ClustalO (v1.2.4) (http://www.clustal.org/) with output format fixed at .phylip. The output phylip file was then uploaded onto PhyML 3.0 [http://www.atgc-montpellier.fr/phyml/] to calculate the optimal substitution model using Akaike Information Criterion. The selected model was then imported into PhyML (v3.3.20200621) in mpi format to allow multithread analysis. Bootstrap replicates was set at 1000×, and maximum-likelihood phylogenetic tree was rooted on the earliest published genome (accession no.: NC_045512.2). The phylogenetic tree was visualized by META-X (v10.1.7) (https://www.megasoftware.net/).

## References:

DOBIN, A., DAVIS, C. A., SCHLESINGER, F., DRENKOW, J., ZALESKI, C., JHA, S., BATUT, P., CHAISSON, M. & GINGERAS, T. R. 2013. STAR: ultrafast universal RNA-seq aligner. Bioinformatics (Oxford, England), 29, 15-21.

LOVE, M. I., HUBER, W. & ANDERS, S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15, 550.
