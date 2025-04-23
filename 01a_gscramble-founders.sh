#!/bin/bash
#set a job name
#SBATCH --job-name=01a_gscramble-founders
#SBATCH --output=./err-out/01a_gscramble-founders.%A_%a.out
#SBATCH --error=./err-out/01a_gscramble-founders.%A_%a.err
#SBATCH --partition=cpu_compute
################
#SBATCH --time=00:00:00
#################
#SBATCH --array=1-13
#################

module load R
module load plink

# set iterations
iterations=1000

# popA is the "main" pop to scramble, through F1 of B and C
popA=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $1}' ./data/gsps/scrambled-pop-array.txt)
popB=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $2}' ./data/gsps/scrambled-pop-array.txt)
popC=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $3}' ./data/gsps/scrambled-pop-array.txt)

# out directory to save scrambled pops
outdir=/lustrefs/nwrc/projects/DeSaix/MO-method/out/gscamble-workflow/scrambled-founders

# tmp directory for intermediary files to be deleted
tmpdir=${outdir}/tmp.${popA}
mkdir -p ${tmpdir}

# plink reference genotype prefix
ref_genos=/lustrefs/nwrc/projects/DeSaix/MO-method/data/genotypes/MOMethod-SpatialPops1672x27866

################################
# Scramble popA
################################

# subset reference genotypes to pops of interest to speed up gscramble operations
prefix=${tmpdir}/${popA}.${popB}.${popC}
awk -v popA=$popA -v popB=$popB -v popC=$popC '$1 == popA || $1 == popB || $1 == popC {print $1, $2}' ${ref_genos}.fam > ${prefix}.plink_keep.txt
plink --noweb --bfile ${ref_genos} --keep ${prefix}.plink_keep.txt --recode --out ${prefix}

Rscript 01a_gscramble-founders.R ${popA} ${popB} ${popC} ${prefix} ${outdir} ${iterations}

rm -rf ${tmpdir}
