# Multi-omic Signatures of Host Response Associated with Presence, Type, and Outcome of Enterococcal Bacteremia. 

Here all the data associated with the manuscript Multi-omic Signatures of Host Response Associated with Presence, Type, and Outcome of Enterococcal Bacteremia can be found.  

The data/ code contained within this repo is sufficient to conduct all of the analysis reported in the paper downstream of Fragpipe (proteomics), MZMine3/GNPS (metabolomics) processing.  


The directory structure is as follows:   

* 00_r_input_data = "raw data" (meaning outputs of Fragpipe, MZmine3, GNPS, etc) and metadata associated with the study that have not been further processed by me. Includes clinical metadata, various mapping files, etc.   
* 01_r_processed_data = data that was processed in R by me. Includes normalization, combining data and metadata, etc. Processing was performed in manuscript_figures.rmd.   
* 02_figures = .pdf figures for main text of manuscript as output by manuscript_figures.rmd. Post processing/ assembly performed in inkscape in the .svg file.  
* 03_supplemental = .pdf figures for the supplemental as output by manuscript_figures.rmd. Post processing/ assembly into performed in inkscape in the .svg file.   
* 04_shiny_app = code to build the R shiny companion web app hosted at: https://charliebayne.shinyapps.io/EcB_Multiomics/
* 05_cheatsheets = cheatsheets containing info about the genes/metabolites/GO terms. I manually compiled this and found it useful when writing the manuscript.    
* 06_manuscript_edits = edits to the manuscript after co-authors took a look.   
* 07_final_manuscript = The final product!   


manuscript_figures.rmd = rmarkdown file containing the code to generate figures.    
[manuscript_figures.html](manuscript_figures.html) = knitted .html files containing all the figures produced by the rmarkdown file.    
functions.R  = R functions I created and used to generate figures in manuscript.     



