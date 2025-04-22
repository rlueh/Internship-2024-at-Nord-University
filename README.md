Genome assemblies are  downloaded using downlaod_teleostei_annotations.sh.

nuc_complexity.R and its corresponding C code is not included in this repository, but can be downloaded from its own repository here:
https://github.com/lmjakt/nuc_complexity

Most important files and what they include: 
- functions.R contains all the functions used throughout the code. 
- assembly_info contains code to parse assembly data and print relevant info to .csv to be read by other scripts
- analyze_regions/analyze_regions.R contains main logic used to generate data for final analysis
- sp_lost_akrab.txt is used as a reference for species of Teleostei that have lost the aKrab domain

The other files include some quite chaotic code snippets used for various testing and analysis.
