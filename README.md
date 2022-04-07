# CellDestiny (HadjAbed-et-al. 2022)
--------------

<p align="center" width="100%">
    <img width="40%" src="https://github.com/TeamPerie/HadjAbed-et-al._2022/blob/main/logos/logoCelldestiny.jpg">
</p>

This is the code used to analyse the lineage tracing datasets in the paper:


CellDestiny: A RShiny application for the visualization and analysis of single-cell lineage tracing data (L. Hadj-Abed et al. 2022). Full details of the approach are given in the methods section of the paper.


### Installing the package


devtools::install_github("Louisahadj/CellDestiny", quiet = TRUE)
note: We need to move CellDestiny from Louisa repository to perie repository



### Dependencies

shiny, shinydashboard, shinyWidgets, shinycustomloader, gplots, ggplot2, ggtern,
ggforce, ggpubr, gridExtra, plyr, dplyr, reshape2, vroom, stringr, tibble, vegan,
scales, tidyr, RColorBrewer, rlang, corrplot, stats, devtools


### What is included in the folders of this repository?

**LentiviralBarcodingData:** this folder contains all of the raw data and scripts necessary to reproduce figures from panels 2-5. The folder is subdivided into a QC visualisation script and a data analysis script, within which there are scripts to create metadata for the experiment and to perform analysis and visualisation of the data


**GeneTherapy:** this folder contains all of the raw data and scripts necessary to reproduce figures from panels 6-8. The folder contains a data analysis script, within which there are scripts to create metadata for the experiment and to perform analysis and visualisation of the data




### LICENSE

<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.
