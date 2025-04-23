#!/bin/bash
#set a job name
#SBATCH --job-name=02_gscramble-admixture
#SBATCH --output=./err-out/02_gscramble-admixture.%A_%a.out
#SBATCH --error=./err-out/02_gscramble-admixture.%A_%a.err
#SBATCH --partition=cpu_compute
################
#SBATCH --time=00:00:00
#SBATCH --ntasks=4
#SBATCH --mem=8G
#################
#SBATCH --array=1-156%50
#################

module load R
module load plink
module load admixture

# popA is the "main" pop, i.e. in backcrosses a simulated individual is crossed to popA
# 156 total
popA=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $1}' ./data/gsps/pop-array.txt)
popB=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $2}' ./data/gsps/pop-array.txt)

# create tmp outdirectory
outdir=/lustrefs/nwrc/projects/DeSaix/MO-method/out/gscamble-workflow/tmp.${popA}.${popB}
mkdir -p ${outdir}
# final results out directory
summary_outdir=/lustrefs/nwrc/projects/DeSaix/MO-method/out/gscamble-workflow/summary/subsets


# plink reference genotype prefix
ref_genos=/lustrefs/nwrc/projects/DeSaix/MO-method/data/genotypes/MOMethod-SpatialPops1672x27866

################################
# Run full gsp (F1, F1B1, F1B2)
################################

# subset reference genotypes to pops of interest to speed up gscramble operations
prefix=${outdir}/${popA}.${popB}
awk -v popA=$popA -v popB=$popB '$1 == popA || $1 == popB {print $1, $2}' ${ref_genos}.fam > ${prefix}.plink_keep.txt
plink --noweb --bfile ${ref_genos} --keep ${prefix}.plink_keep.txt --recode --out ${prefix}

# specify gsp and number of iterations to run gsp
gsp=full
start_iteration=1
final_iteration=500

Rscript 02_gscramble-pedigrees-admix-founders.R ${gsp} ${start_iteration} ${final_iteration} ${popA} ${popB} ${outdir} ${ref_genos} ${summary_outdir}

rm -rf ${outdir}

################################
# Run f1b gsp to supplement
################################

mkdir -p ${outdir}
# subset reference genotypes to pops of interest to speed up gscramble operations
prefix=${outdir}/${popA}.${popB}
awk -v popA=$popA -v popB=$popB '$1 == popA || $1 == popB {print $1, $2}' ${ref_genos}.fam > ${prefix}.plink_keep.txt
plink --noweb --bfile ${ref_genos} --keep ${prefix}.plink_keep.txt --recode --out ${prefix}

# specify gsp and number of iterations to run gsp
gsp=f1b
start_iteration=1
final_iteration=125

Rscript 02_gscramble-pedigrees-admix-founders.R ${gsp} ${start_iteration} ${final_iteration} ${popA} ${popB} ${outdir} ${ref_genos} ${summary_outdir}

rm -rf ${outdir}



