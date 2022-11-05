# Environmental and genealogical effects on DNA methylation in a widespread apomictic dandelion lineage

This repository reproduces the results presented in [*"Environmental and genealogical effects on DNA methylation in a widespread apomictic dandelion lineage"*](link2journal) by V.N. Ibañez, M. van Antro, C. Peña Ponton, S. Ivanovic, C.A.M. Wagemaker, F. Gawehns, K.J.F. Verhoeven.

DNA methylation that occurs in CG sequence context shows transgenerational stability and high epimutation rate, and can thus provide genealogical information at short time scales. Here, [epiGBS2 protocol](https://github.com/nioo-knaw/epiGBS2) is used to analyze DNA methylation between accessions from a geographically widespread, apomictic common dandelion *(Taraxacum officinale)* lineage grown experimentally under different light conditions.

We show that the light treatment induced  differentially methylated cytosines (DMCs) in all sequence contexts, with a bias toward transposable elements. Accession differences were associated with DMCs in CG context. Hierarchical clustering of samples based on total mCG profiles revealed a perfect clustering of samples by accession identity, irrespective of light conditions. 

Using microsatellite information as a benchmark of genetic divergence within the clonal lineage, we show that genetic divergence between accessions correlates strongly with overall mCG profiles. However, our results suggest that environmental effects that do occur in CG context might produce a heritable signal that partly dilutes the genealogical signal. 

Our methodology  can be used as a tool for reconstructing micro-evolutionary genealogy, particularly for systems lacking genetic variation, such as clonal and vegetatively propagated plants.

## Processing of raw multiplexed read sequences

A total of 80 samples were multiplexed together in the same sequence library. Half of them, were digested with the restriction enzymes *AseI* or  *Csp6* and NsiI. 
In this manuscript we present results from the *Csp6*I - *Nsi*I digested epiGBS samples, as these yielded higher sequencing output than the *Ase*I - *Nsi*I based samples. However, the following scripts can handle both sets of data.

In order to proceed, you will need the following files:
  - [raw multiplexed read sequences](https://doi.org/10.5281/zenodo.6793166), 
  - [*Csp*6-*Nsi*I_barcode.tsv](https://doi.org/10.5281/zenodo.6793166),
  - [*Ase*6-*Nsi*I_barcode.tsv](https://doi.org/10.5281/zenodo.6793166),
  - [methylation.bed](https://doi.org/10.5281/zenodo.6793166), 
  - [consensus_cluster.renamed_csp6.fa](https://doi.org/10.5281/zenodo.6793166) and 
  - [consensus_cluster.renamed_aseI.fa](https://doi.org/10.5281/zenodo.6793166).
  
Also, the configuration files used with the [epiGBS pipeline](https://github.com/nioo-knaw/epiGBS2) are provided.
  - [*Csp*6-*Nsi*I_config.yaml](https://doi.org/10.5281/zenodo.6793166) and 
  - [*Ase*6-*Nsi*I_config.yaml](https://doi.org/10.5281/zenodo.6793166).

Raw read data should be processed with [epiGBS pipeline](https://github.com/nioo-knaw/epiGBS2) following the steps in [*Preparation to run the pipeline*](https://github.com/nioo-knaw/epiGBS2#preparation-to-run-the-pipeline).


To compare, the reports obtained after running the [epiGBS pipeline](https://github.com/nioo-knaw/epiGBS2) are also provided:
  - [*Csp*6-*Nsi*I](https://doi.org/10.5281/zenodo.6793166) and 
  - [*Ase*I-*Nsi*I](https://doi.org/10.5281/zenodo.6793166).

## Downstream analysis of methylation data

Below there are notebooks demonstrating how data was processing in the article using a subset of the sequenced cytosines.

|# |Script|Description| Notebook|
|:-:|----|:------:|:---:|
|1|[01_epiTree.R](https://github.com/VeronicaNoe/epiTree/blob/main/Rscripts/01_filterMethylation.R)| Generate a filtered methylation file |[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/VeronicaNoe/epiTree/blob/main/notebooks/01_filterMethylation.ipynb)|
|2|[02_epiTree.R](https://github.com/VeronicaNoe/epiTree/blob/main/Rscripts/02_characterizeOverallMethylation.R)| Characterize overall methylation levels|[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/VeronicaNoe/epiTree/blob/main/notebooks/02_characterizeOverallMethylation.ipynb)|
|3|[03_epiTree.R](https://github.com/VeronicaNoe/epiTree/blob/main/Rscripts/03_distances.R)| Generate distances matrices |[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/VeronicaNoe/epiTree/blob/main/notebooks/03_distances.ipynb)|
|4|[04_epiTree.R](https://github.com/VeronicaNoe/epiTree/blob/main/Rscripts/04_differentialCytosineMethylationWithDSS.R)| Obtain the differential Methylated Cytosines with DSS |[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/VeronicaNoe/epiTree/blob/main/notebooks/04_differentialCytosineMethylationWithDSS.ipynb)|
|5|[05_epiTree.R](https://github.com/VeronicaNoe/epiTree/blob/main/Rscripts/05_manhattanPlot.R)| Obtain Manhattan plots |[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/VeronicaNoe/epiTree/blob/main/notebooks/05_manhattanPlot.ipynb)|
|6|[06_epiTree.R](https://github.com/VeronicaNoe/epiTree/blob/main/Rscripts/06_dendrogram.R)| Obtain dendrograms |[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/VeronicaNoe/epiTree/blob/main/notebooks/06_dendrogram.ipynb)|
|7|[07_epiTree.R](https://github.com/VeronicaNoe/epiTree/blob/main/Rscripts/07_DMC_description.R)| Characterize DMC |[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/VeronicaNoe/epiTree/blob/main/notebooks/07_DMC_description.ipynb)|
|8|[08_epiTree.R](https://github.com/VeronicaNoe/epiTree/blob/main/Rscripts/08_mantelTest.R)| Mantel test |[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/VeronicaNoe/epiTree/blob/main/notebooks/08_mantelTest.ipynb)|
