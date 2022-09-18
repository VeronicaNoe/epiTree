# Environmental and genealogical signals on DNA methylation in a widespread apomictic dandelion lineage

This repository reproduces the results presented in [*"Environmental and genealogical signals on DNA methylation in a widespread apomictic dandelion lineage"*](link2journal) by V.N. Ibañez, M. van Antro, C. Peña Ponton, S. Ivanovic, C.A.M. Wagemaker, F. Gawehns, K.J.F. Verhoeven. [download preprint version](link here).

DNA methylation that occurs in CG sequence context shows transgenerational stability and high epimutation rate, and can thus provide genealogical information at short time scales. Here, [epiGBS2 protocol](https://github.com/nioo-knaw/epiGBS2) is used to analyze DNA methylation between accessions from a geographically widespread, apomictic common dandelion *(Taraxacum officinale)* lineage grown experimentally under different light conditions.

We show that the light treatment induced  differentially methylated cytosines (DMCs) in all sequence contexts, with a bias toward transposable elements. Accession differences were associated with DMCs in CG context. Hierarchical clustering of samples based on total mCG profiles revealed a perfect clustering of samples by accession identity, irrespective of light conditions. 

Using microsatellite information as a benchmark of genetic divergence within the clonal lineage, we show that genetic divergence between accessions correlates strongly with overall mCG profiles. However, our results suggest that environmental effects that do occur in CG context might produce a heritable signal that partly dilutes the genealogical signal. 

Our methodology  can be used as a tool for reconstructing micro-evolutionary genealogy, particularly for systems lacking genetic variation, such as clonal and vegetatively propagated plants.

## Processing of raw multiplexed read sequences

The 40 samples of this study were multiplexed together with 40 additional epiGBS samples in the same sequencing library. 
These additional samples consist of the same 40 experimental plants after digestion with the restriction enzymes AseI and NsiI. 
In this manuscript we present results from the *Csp6*I - *Nsi*I digested epiGBS samples, as these yielded higher sequencing output than the *Ase*I - *Nsi*I based samples. The results of the latter are presented as supplementary information.

In order to proceed, you will need the [raw multiplexed read sequences](link2ENA), [*Csp*6-*Nsi*I_barcode.tsv](link2zenodo) and [*Ase*6-*Nsi*I_barcode.tsv](link2zenodo) files. In addition, the [*Csp*6-*Nsi*I_config.yaml](link2zenodo) and [*Ase*6-*Nsi*I_config.yaml](link2zenodo) configuration files are provided.
Raw read data should be processed following the steps in [*Preparation to run the pipeline*](https://github.com/nioo-knaw/epiGBS2#preparation-to-run-the-pipeline). 

The outputs from epiGBS pipeline used in the following section are [methylation.bed](link2zenodo), [consensus_cluster.renamed_csp6.fa](link2zenodo) and [consensus_cluster.renamed_aseI.fa](link2zenodo).

To compare, the reports for [*Csp*6-*Nsi*I](link2zenodo) and [*Ase*I-*Nsi*I](link2zenodo) obtained after run the epiGBS pipeline are also provided.

## Downstream analysis of methylation data

Below there are notebooks demonstrating how data was processing in the article using a subset of the sequenced cytosines.

|# |Script|Description| Notebook|
|:-:|----|:------:|:---:|
|1|[01_filterMethylation.R](https://github.com/VeronicaNoe/epiTree/blob/main/Rscripts/01_filterMethylation.R)| Generate a filtered methylation file |[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/VeronicaNoe/epiTree/blob/main/notebooks/01_filterMethylation.ipynb)|
|2|[02_characterizeOverallMethylation.R](https://github.com/VeronicaNoe/epiTree/blob/main/Rscripts/02_characterizeOverallMethylation.R)| Generate a filtered methylation file |[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/VeronicaNoe/epiTree/blob/main/notebooks/02_characterizeOverallMethylation.ipynb)|
|3|[03_distances.R](https://github.com/VeronicaNoe/epiTree/blob/main/Rscripts/03_distances.R)| Generate a filtered methylation file |[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/VeronicaNoe/epiTree/blob/main/notebooks/03_distances.ipynb)|
|4|[04_differentialCytosineMethylationWithDSS.R](https://github.com/VeronicaNoe/epiTree/blob/main/Rscripts/04_differentialCytosineMethylationWithDSS.R)| Generate a filtered methylation file |[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/VeronicaNoe/epiTree/blob/main/notebooks/04_differentialCytosineMethylationWithDSS.ipynb)|
|5|[05_manhattanPlot.R](https://github.com/VeronicaNoe/epiTree/blob/main/Rscripts/05_manhattanPlot.R)| Generate a filtered methylation file |[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/VeronicaNoe/epiTree/blob/main/notebooks/05_manhattanPlot.ipynb)|
|6|[06_filterMethylation.R](https://github.com/VeronicaNoe/epiTree/blob/main/Rscripts/01_filterMethylation.R)| Generate a filtered methylation file |[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/VeronicaNoe/epiTree/blob/main/notebooks/01_filterMethylation.ipynb)|
