#!/bin/bash
#set a job name
#SBATCH --job-name=03_ghost-gscramble-admixture
#SBATCH --output=./err-out/03_ghost-gscramble-admixture.%A_%a.out
#SBATCH --error=./err-out/03_ghost-gscramble-admixture.%A_%a.err
#SBATCH --partition=cpu_compute
################
#SBATCH --time=00:00:00
#SBATCH --mem=8G
#SBATCH --ntasks=4
#################
#SBATCH --array=1-91%50
#################

module load R
module load plink
module load admixture

# popA is the "main" pop, i.e. in backcrosses a simulated individual is crossed to popA
# 91 total
popA=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $1}' ./data/gsps/ghost-full-pop-gsp-iter-array.txt)

gsp=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $2}' ./data/gsps/ghost-full-pop-gsp-iter-array.txt)
start_iteration=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $3}' ./data/gsps/ghost-full-pop-gsp-iter-array.txt)
final_iteration=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $4}' ./data/gsps/ghost-full-pop-gsp-iter-array.txt)

# plink reference genotype prefix
ref_genos=/lustrefs/nwrc/projects/DeSaix/MO-method/data/genotypes/MOMethod-SpatialPops1672x27866
# plink ghost genotype prefix
ghost_genos=/lustrefs/nwrc/projects/DeSaix/MO-method/data/genotypes/MOMethod-Ghosts3382x27866

# final results out directory
summary_outdir=/lustrefs/nwrc/projects/DeSaix/MO-method/out/gscamble-workflow/summary/ghosts/full-scramble/subsets

################################
# Run gsp (scramble MO pop, then gsp with preserved individuals
################################

# create tmp outdirectory
outdir=/lustrefs/nwrc/projects/DeSaix/MO-method/out/gscamble-workflow/ghosts/full-scramble/tmp.ghost.${popA}.${gsp}.${start_iteration}.${final_iteration}
mkdir -p ${outdir}

Rscript 03_ghosts-gscramble-admix.R ${gsp} ${start_iteration} ${final_iteration} ${popA} ${outdir} ${ref_genos} ${ghost_genos} ${summary_outdir}

rm -rf ${outdir}
