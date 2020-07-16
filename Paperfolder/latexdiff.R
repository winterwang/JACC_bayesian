##%######################################################%##
#                                                          #
####                  make traced pdf                   ####
#                                                          #
##%######################################################%##
# library(latexdiffr)
# setwd("/Paperfolder")
# latexdiff("Paperfolder_v0_0_0.Rmd", "Paperfolder.Rmd", 
#           output = "diff_v0_0_1", clean = FALSE)
system("latexdiff Paperfolder_v0_0_0.tex Paperfolder.tex > diff_v_0_0_1.tex")

# abstract diff: Paperfolder_v0_0_0.tex Paperfolder.tex 20200716
# \DIFaddbegin \DIFadd{In conclusion, %DIF >
#     data from the JACC study has provided strong evidence that daily milk %DIF >
#     intake among Japanese men was associated with delayed and lower hazard %DIF >
#     of mortality from stroke especially cerebral infarction. }\DIFaddend