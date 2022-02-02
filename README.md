# Nextcast: a software suite to analyze and model toxicogenomics data

Nextcast is a collection of tools for the analysis of toxicogenomics and cheminformatics data developed at FHAIVE (https://www.fhaive.fi/).

## How to clone the repository

```git clone --recurse-submodules git@github.com:fhaive/nextcast.git```


## Tool descriptions

|**eUTOPIA**| |
| ------------- |--------------|
|**Input**| Microarray raw data |
|| Phenodata describing the samples in the raw data files|
|**Output**| Normalized an preprocessed transcriptomic matrix (e.g. gene expression matrix) |
|| List of differentially expressed genes, with log-fold-changes and p-values for each comparison performed during the analysis|
|**Summary**|Graphically accessible guided workflow for preprocessing and analysis of omics data. Supports Agilent 2-color, Agilent 1-color, Affymetrix, and Illumina methylation microarray platforms (Ongoing efforts to add support for RNA-Seq data). Discreetly separated steps in analysis designed in R Shiny, incorporates widely used microarray analysis practices and R packages. Reporting is and data interpretation is leveraged from dynamically generated plots.|
|**Details**|[Paper](https://doi.org/10.1186/s13029-019-0071-7), [Repository](https://github.com/Greco-Lab/eUTOPIA)|
|**Tags**|\#bioinformatics \#analytics \#preprocessing|
|**Language**|R|
|**GUI**|Shiny|

|**FunMappOne**| |
| ------------- |--------------|
|**Input**| Excel file containing lists of genes and (optionally) some scores (e.g. log-fold-change, BMD, etc) |
|**Output**| Heatmap with enriched pathways or gene ontology terms. Excel file containing the results of the enrichment|
|**Summary**|A user-friendly graphical interface that allows to visualize and summarize the functional annotations of one or multiple molecular biology experiments at once.|
|**Details**|[Paper](https://doi.org/10.1186/s12859-019-2639-2), [Github](https://github.com/Greco-Lab/FunMappOne)|
|**Tags**|\#bioinformatics \#analytics \#functional-annotation|
|**Language**|R|
|**GUI**|Shiny|

|**INfORM**| |
| ------------- |--------------|
|**Input**| A gene expression matrix |
||The result of a differential expression analysis as a table with genes, log-fold-change and p-values|
|**Output**| A gene network|
||The list of genes beloging to each community identified on the network, along with the enriched GO terms|
|**Summary**|A novel computational method and its R and web-based implementations, to perform inference of gene network from transcriptome data and prioritization of key genes with central functional and topological role in the network.|
|**Details**|[Paper](https://doi.org/10.1093/bioinformatics/bty063), [Github](https://github.com/Greco-Lab/INfORM)|
|**Tags**|\#bioinformatics \#analytics \#network-analysis|
|**Language**|R|
|**GUI**|Shiny|

|**VOLTA**| |
| ------------- |--------------|
|**Input**| One or multiple gene networks |
|**Output**| |
|**Summary**|VOLTA is a Python network analysis package, suited for different types of networks but with a focus on co-expression network analysis. The goal of VOLTA is to provide all functionalities needed for a comprehensive network analysis and comparison in a single package. Its aim is to expose all functionalities in order to allow users to build their own analysis pipelines and adjust functionalities as needed. Additional complete analysis pipelines are provided in the form of Jupyter Notebooks.|
|**Details**|[Paper](https://doi.org/10.1093/bioinformatics/btab642), [Github](https://github.com/fhaive/VOLTA)|
|**Tags**|\#bioinformatics \#analytics \#network-analysis|
|**Language**|Python|

|**BMDx**| |
| ------------- |--------------|
|**Input**| Gene expression data derived from dose-time exposure series. Multiple exposures can be analysed at the same time |
|**Output**| List of dose-dependent genes, with estimated effective doses, for each exposure time point.|
|**Summary**|BMDx is a R based software with a graphical interphace in Shiny for the dose-dependent analysis of toxicogenomic datasets.|
|**Details**|[Paper](https://doi.org/10.1093/bioinformatics/btaa030), [Github](https://github.com/Greco-Lab/BMDx)|
|**Tags**|\#bioinformatics \#analytics \#dose-response|
|**Language**|R|
|**GUI**|Shiny|

|**TinderMIX**| |
| ------------- |--------------|
|**Input**| Gene expression data derived from dose-time exposure series. Multiple exposures can be analysed at the same time |
|**Output**| List of dynamic dose-dependent genes, with estimated point of departure|
|**Summary**|TinderMIX is a new computational framework for dose- and time- dependent gene expression analysis which aims to identify groups of genes that show a time and dose-response behaviour.|
|**Details**|[Paper](https://doi.org/10.1093/gigascience/giaa055), [Github](https://github.com/grecolab/TinderMIX)|
|**Tags**|\#bioinformatics \#analytics \#dynamic dose-response|
|**Language**|R|

|**MVDA**| |
| ------------- |--------------|
|**Input**| Multiple omics data (views) for the same set of samples  |
|**Output**| A multi-view clustering of the samples|
|**Summary**|Multi-view clustering methodology in which the information from different data layers (views) is integrated at the levels of the results of each single view clustering iterations.|
|**Details**|[Paper](https://doi.org/10.1186/s12859-015-0680-3), [Github](https://github.com/Greco-Lab/MVDA_package)|
|**Tags**|\#bioinformatics \#modelling \#multi-view-clustering|
|**Language**|R|

|**MOSIM**| |
| ------------- |--------------|
|**Input**| Desired network topology parameters (e.g.: number of genes, number of miRNAs and desired degree distribution), as well as parameters about the expression data (e.g.: number of conditions, number of samples per condition, noise distribution) |
|**Output**| A gene regulatory network with the desired properties together with simulated expression data. |
|**Summary**| A genomic multi-view data simulator that generates expression data from a simulated regulatory network comprising genes and miRNAs.|
|**Details**|[Paper](https://doi.org/10.1186/s12859-015-0577-1), [Source code](https://static-content.springer.com/esm/art%3A10.1186%2Fs12859-015-0577-1/MediaObjects/12859_2015_577_MOESM1_ESM.zip)|
|**Tags**|\#bioinformatics \#modelling \#simulator|
|**Language**|R|

|**FPRF**| |
| ------------- |--------------|
|**Input**| Multiple omics data (views) for the same set of samples|
|**Output**| A predictive model of the interest variables based on a reduced list of prioritized features according to their relevance in the learning phase.|
|**Summary**| A feature selection algorithm for multi-omics data. The algorithm is based on the Random Forest algorithm and a novel feature selection strategy based on the data transformation process named fuzzy patterns. |
|**Details**|[Paper](https://doi.org/10.1371/journal.pone.0107801)||
|**Tags**|\#bioinformatics \#modelling \#feature-selection|

|**GARBO**| |
| ------------- |--------------|
|**Input**| A high-dimensional omics dataset and class labels|
|**Output**| A collection of the best performing classifiers built upon the maximally relevant and minimally redundant feature sets|
|**Summary**| A classification framework for high accuracy predictions and biomarker selection in high-dimensional omics data based on Random Forests and Genetic Algorithms|
|**Details**|[Paper](https://doi.org/10.1093/bioinformatics/btaa144), [Github](https://github.com/Greco-Lab/GARBO)|
|**Tags**|\#bioinformatics \#modelling \#feature-selection|
|**Language**|Python|

|**INSIdE NANO**| |
| ------------- |--------------|
|**Input**| Sets of phenotypic entities (nanomaterials, drugs, diseases and chemicals)|
|**Output**| List of cliques found in the network connecting the input phenotypic entities|
|**Summary**|INSIdE nano is a web-based tool that highlights connections between phenotypic entities based on their effects on genes. The database behind the INSIdE nano is a network whose nodes are grouped into four categories: Nanomaterial exposures Drug treatments Chemical exposures Diseases For each node in the network, information regarding their effects on the genes is collected. Weights of the edges in the network are correlated to the similarity of these effects between each couple of nodes. More information and usage tutorials are available at: http://inano.biobyte.de/help.cgi|
|**Details**|[Paper](https://doi.org/10.1038/s41598-018-37411-y), [Online Tool](http://inano.biobyte.de/)|
|**Tags**|\#read-across|

|**MANGA**| |
| ------------- |--------------|
|**Input**| Chemical compounds represented as molecular descriptors and their corresponding continuous endpoint activity|
|**Output**| Predictive models as well as reports of internal and external validation measures|
|**Summary**| A QSAR modelling framework based on Multi-Objective Genetic Algorithms for high robustness and stability of the predictive models|
|**Details**|[Paper](https://doi.org/10.1093/bioinformatics/btz521),[Github](https://github.com/Greco-Lab/MaNGA)|
|**Tags**|\#cheminformatics \#QSAR|
|**Language**|Python|

|**HyQSAR**| |
| ------------- |--------------|
|**Input**| Chemical compounds represented as molecular descriptors or substructure fingerprints and their corresponding gene expression data, as well as a numerical activity/property of interest|
|**Output**| Predictive QSAR models based on reduced set of features trained with regularized linear models|
|**Summary**|hyQSAR is a set of R scripts that can be used to perform hybrid QSAR modelling by integrating structural properties of chemical compounds and their molecular mechanism-of-action (MOA) information.|
|**Details**|[Paper](https://doi.org/10.1186/s13321-019-0359-2), [Source code](https://static-content.springer.com/esm/art%3A10.1186%2Fs13321-019-0359-2/MediaObjects/13321_2019_359_MOESM3_ESM.r)|
|**Tags**|\#cheminformatics \#QSAR|
|**Language**|R|



