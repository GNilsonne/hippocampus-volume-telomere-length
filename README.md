# hippocampus-volume-telomere-length
Analysis code for a meta-analysis of the relationship between leukocyte telomere length and hippocampus volume

This repository contains the following analysis scripts:
- Meta-analysis main script.R. This script performs meta-analysis and generates a forest plot, funnel plots, and Q-Q plots. Data are self-contained in the script.
- Meta-analysis updated 2017.R. This script performs an updated meta-analysis, adding data from the Stockholm Sleepy Brain study.
- Analyses of TA_TL ratio.R. This script analyses the ratio between telomerase activity and telomere length as predictor for hippocampus volume. Data files listed below are required. 

This repository contains the following data files, estimated from scatter plots in Jacobs et al. 2014, doi: 10.1001/jamaneurol.2014.870:
- JacobsL.csv
- JacobsR.csv

This repository contains the following data files, estimated from scatter plots in Wolkowitz et al. 2015, doi: 10.1016/j.pscychresns.2015.01.007:
- Wolkowitz_TA_Con.csv
- Wolkowitz_TA_MDD.csv
- Wolkowitz_TL_Con.csv
- Wolkowitz_TL_MDD.csv
- Wolkowitz_merge_Con.csv
- Wolkowitz_merge_MDD.csv

This repository contains the following pdf files, showing original scatterplots and pltos of estimated data side by side to facilitate inspection to validate data extraction accuracy:
- Jacobs2014 data extraction validation 1.pdf
- Jacobs2014 data extraction validation 2.pdf
- Jacobs2014 data extraction validation 3.pdf
- Jacobs2014 data extraction validation 4.pdf
- Jacobs2014 data extraction validation 5.pdf
- Wolkowitz2015 data extraction validation 1.pdf
- Wolkowitz2015 data extraction validation 2.pdf

This repository contains the list of PubMed search results, using the string "telomer*hippocampus", including titles and abstracts, from two different dates:
- PubMed search result 150714.txt
- PubMed search result 150922.txt
- PubMed search result 150922.txt (limited from 150921 to 170616)
