#!/bin/bash
#set a job name
#SBATCH --job-name=04c_5pops-complex-admixture
#SBATCH --output=./err-out/04c_5pops-complex-admixture.%A_%a.out
#SBATCH --error=./err-out/04c_5pops-complex-admixture.%A_%a.err
#SBATCH --partition=cpu_compute
################
#SBATCH --time=00:00:00
#SBATCH --mem=8G
#SBATCH --ntasks=8
#################
#SBATCH --array=1-10
#################
module load R
module load plink
module load admixture

start_iteration=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $1}' ./data/gsps/ghost-complex-iter-array.txt)
final_iteration=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $2}' ./data/gsps/ghost-complex-iter-array.txt)

# plink reference genotype prefix
ref_genos=/lustrefs/nwrc/projects/DeSaix/MO-method/data/genotypes/MOMethod-SpatialPops1672x27866

# final results out directory
summary_outdir=/lustrefs/nwrc/projects/DeSaix/MO-method/out/gscamble-workflow/summary/complex/5pops

################################
# Run gsp (scramble MO pop, then gsp with preserved individuals
################################

# create tmp outdirectory
outdir=/lustrefs/nwrc/projects/DeSaix/MO-method/out/gscamble-workflow/complex/5pops/tmp.5pops.${start_iteration}.${final_iteration}
mkdir -p ${outdir}

Rscript 04c_complex_5pops.R ${start_iteration} ${final_iteration} ${outdir} ${ref_genos} ${summary_outdir}

rm -rf ${outdir}
