# nest_arch_comp_analysis
This repository contains the code and data files needed to reproduce the analysis from "Foraging behavior affects nest architecture in a cross-species comparison of ant nests" by Sean O'Fallon, Kim Drager, Art Zhao, Andy Suarez, and Noa Pinter-Wollman.


Once the code and data are saved to your computer, you will need to change the file paths in the script to match your local directory. The lines that need editing should be obvious - they are the lines that call read.csv on a filepath through "C:/Users/seano/Desktop/Projects/...". Once those read.csv function calls are changed to reflect the location of the data files on your computer, the R script in this repo should run without further modification. You'll see that the script begins with loading in data files and cleaning data objects, then proceeds to the analyses performed in the study, and concludes with the outputs found in the main text and supplementary material.


The data files:

Chamber_Widths_per_Chamber - contains the widths (in cm) of all measured chambers in nests for which a network was constructed. Col A has nest ID, B has node type, C has Node ID, D has maximum width (cm).

e-chamber-widths_no-net - for some nests, we could not construct a network, but could extract an entrance chamber width measurement. The widths of the entrance chambers for those nests are contained here. Col A contains nest ID, B contains species, C is left blank, and D contains maximum width (cm).

foraging-by-species_upload - contains the foraging strategy designation for all species included in this study. Col A contains species, B contains the designation given by Lanan 2014, C contains the designation used in this study, D contains the taxonomic level for which the strategy is attributed based on evidence from the literature (either genus or species), E contains the primary literature source on which the designation is based, F contains the type of evidence the source provides, G contains links to the source, H and I contain information on corroborating reports, J contains notes (can be safely ignored), and K contains search terms used to find sources.

morpho_upload - contains the morphometric measurements used to scale nest features for each species used in this study. A contains the species, B contains the specimen ID, C contains a link to the image used in head measurements, D contains the width of the specimen's head at its eyes (cm), E contains the maximum width of the specimen's head (if different than at eyes) (cm), F contains the width of the specimen's mandibles, G contains a link to the image used to measure the specimen's head length and Weber's length, H contains the length of the specimen's head, I contains the specimen's Weber's length.

nest-summaries_upload - for each nest included in this study, indicates which features could be measured. Also includes nest depths in meters. Col A contains nest ID, B contains species, C and D contain information on the files from which the nest measurements were taken, E through H indicate whether the nest was networked, whether entrance chamber width was measured, how many entrances the nest has (not used in this study), and whether tunnel length was measured (not used in this study), respectively. Col I contains nest depth (m).

num-cham_no-net - for some nests used in this study, the number of chambers could be extracted despite the nest not being constructed as a network. The number of chambers in those nests are contained here. Col A contains nest ID, B contains species, and C contains the number of chambers that nest has.
