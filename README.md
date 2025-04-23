# missouri-method
Analysis scripts for Missouri method workflow with supervised admixture

April 23, 2025
Matt DeSaix
The majority of scripts relevant to the analyses I did are in the `./hpc/` directory, which is all the scripts and data pulled from the NBAF HPC where simulations were run. Files of the same name but different suffixes (i.e., .sh and .R) refer to bash scripts that were submitted on the server and they subsequently run the corresponding R script.
Scripts:
01a_gscramble-founders.* - simulate 1,000 individual founders for a specified population in gscramble
01b_get-mo-method-scrambled-founders – perform supervised admixture analysis on the g-scrambled founders.
02_gscramble-pedigrees-admix-founders – gscramble F1, F1B1, and F1B2 individuals from the scrambled founder dataset.
03_ghosts-gscramble-admix – gscramble ghost individuals
04_* - different pedigrees for specific requests Ben had on determining pedigree variance in estimates
get-mo-method-ghosts.sh – do “Missouri method”, i.e. supervised admixture, for ghost individuals
get-mo-method-loo.sh – do supervised admixture leave-one-out for reference individuals
get-mo-method-matrix.sh – do supervised admixture for matrix individuals
get-mo-method-nonspat10.sh – do supervised admixture for some 10 individuals Tim had identified that didn’t have spatial data

