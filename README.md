# Kaiko_publication
All the data and materials for the publication about the Kaiko, which is a deep learning-based *de novo* peptide sequencing tool for natural bacterial isolates. We have deposited here both the data and the code for various figures in the manuscript.

- Figure 1 - Bacteria represented in training and testing data.
- Figure 2 - Accuracy of spectrum annotation.
- Figure 3 - Sequence prediction accuracy by peptide length.
- Figure 4 - Longest correct substring.
- Figure 5 - Taxonomic classification of natural isolates via proteomics.

- Supplemental Figure 1. Improving deep neural networks with more training data.
- Supplemental Figure 2 - Training and validation error.
- Supplemental Figure 3 - Distribution of peptide lengths for the training and testing dataset.

## Figure 1
We generated a phylogenetic tree of the Kaiko benchmark datasets using graphlan, a python package. Please refer to [this link](https://bitbucket.org/nsegata/graphlan/) for more information. You can see the annotation file and the tree input xml file in the [phylo_vis](phylo_vis) folder.
To tie the annotation file to the input tree,
```bash
$ python python graphlan_annotate.py --annot ../rplB.Tree/annotation_1_tree.txt ../rplB.Tree/rplB_0.xml
```
To generate the output image in a png format with setting up a resolution with --dpi and a size with --size (in inches),
```bash
$ python graphlan.py ../rplB.Tree/rplB_0.xml test.png --dpi 250 --size 7
```

## Figure 2, 3, 4
Please refer to [analysis/figures.ipynb](analysis/figures.ipynb).

## Figure 5
To visualize the output results of taxonomy identification, we used the ETE Toolkit, which is another tool for analysis and visualization of any type of trees. Please refer to [this link](http://etetoolkit.org/) for more information. You can reproduce this from this notebook [taxon_identification/fig5.ipynb](taxon_identification/fig5.ipynb).


## Supplemental Figures
Please refer to [analysis/suppl_figures.ipynb](analysis/suppl_figures.ipynb).

## Contacts

Written by Joon-Yong Lee for the Department of Energy (PNNL, Richland, WA) \
E-mail: proteomics@pnnl.gov \
Website: https://panomics.pnl.gov/ or https://omics.pnl.gov

## License

Licensed under the Apache License, Version 2.0; 
you may not use this file except in compliance with the License.  You may obtain 
a copy of the License at https://opensource.org/licenses/Apache-2.0
