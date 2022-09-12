# epiTree: Environmental and genealogical signals on DNA methylation in a widespread apomictic dandelion lineage

This repository contains the instructions for reproducing results presented in "Environmental and genealogical signals on DNA methylation in a widespread apomictic dandelion lineage" by V.N. Ibañez, M. van Antro, C. Peña Ponton, S. Ivanovic, C.A.M. Wagemaker, F. Gawehns, K.J.F. Verhoeven. [preprint version](link here).

DNA methylation that occurs in CG sequence context shows transgenerational stability and high epimutation rate, and can thus provide genealogical information at short time scales. 
Here we analysed DNA methylation variation between accessions from a geographically widespread, apomictic common dandelion *(Taraxacum officinale)* lineage, when grown experimentally under different light conditions. Using a [epiGBS2 protocol] (https://github.com/nioo-knaw/epiGBS2), we show that the light treatment induced  differentially methylated cytosines (DMCs) in all sequence contexts, with a bias toward transposable elements. 
Accession differences were associated with DMCs in CG context. Hierarchical clustering of samples based on total mCG profiles revealed a perfect clustering of samples by accession identity, irrespective of light conditions. 
Using microsatellite information as a benchmark of genetic divergence within the clonal lineage, we show that genetic divergence between accessions correlates strongly with overall mCG profiles. However, our results suggest that environmental effects that do occur in CG context might produce a heritable signal that partly dilutes the genealogical signal. 
Our study shows that methylation information can be used to reconstruct micro-evolutionary genealogy, providing a useful tool in systems that lack genetic variation such as clonal and vegetatively propagated plants.

## Processing of raw multiplexed read sequences

The 40 samples of this study were multiplexed together with 40 additional epiGBS samples in the same sequencing library. 
These additional samples consist of the same 40 experimental plants after digestion with the restriction enzymes AseI and NsiI. 
In this manuscript we present results from the *Csp6*I - *Nsi*I digested epiGBS samples, as these yielded higher sequencing output than the *Ase*I - *Nsi*I based samples. The results of the latter are presented as supplementary information.

In order to proceed, you will need the [raw multiplexed read sequences](link2ENA), [*Csp*6-*Nsi*I_barcode.tsv](link2zenodo) and [*Ase*6-*Nsi*I_barcode.tsv](link2zenodo) files. In addition, the [*Csp*6-*Nsi*I_config.yaml] and [*Ase*6-*Nsi*I_config.yaml](link2zenodo) configuration files are provided.
Raw read data should be processed following the steps in [*Preparation to run the pipeline*](https://github.com/nioo-knaw/epiGBS2#preparation-to-run-the-pipeline). 

The outputs used in the following section are [methylation.bed](link2zenodo), [consensus_cluster.renamed_csp6.fa](link2zenodo) and [consensus_cluster.renamed_aseI.fa](link2zenodo).

To compare, the reports for [*Csp*6-*Nsi*I](link2zenodo) and [*Ase*I-*Nsi*I](link2zenodo) obtained after run the epiGBS pipeline are also provided.

## Detecting differentially methylated cytosines (DMC)

We also provide demo notebooks with the code used in the article. Pleas, note that only a subset of sequenced cytosines are used due to limited resource from colab and to make a faster excution.

[SSR data.csv](link2zenodo)
