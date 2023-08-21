# KEGG ML

Incorporating KEGG pathways into ML models for Phenotype prediction.

## Folders
- GSE: gene expression data.

- KEGG_pathways: Some kegg pathway files downloaded from the KEGG database. Feel free to add more

- pathway_new_topList_low: Some more KEGG pathways that were associated with the low dose classification problem

## files

- graph.py: underling implementation that extracts kegg graphs and saves them in Python-friendly format
    to be used, for example, with sparse transport maps code. 
- pltPathways.py: plot the pathways using graphViz
- geneSetUtils.py: utils for dealing with the genesets

## Downloading KEGG Pathway xml files

To download a new pathway, search it on google. For instance google search "hsa04137" and clicking the first link to https://www.genome.jp/entry/pathway+hsa04137. Click on the pathway map diagram to go to https://www.genome.jp/pathway/hsa04137. Click download>kgml and this will download the xml file. 


