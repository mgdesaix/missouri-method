#!/bin/bash
#set a job name
#SBATCH --job-name=supervised-admixture-loo
#SBATCH --output=./err-out/supervised-admixture-loo.%A_%a.out
#SBATCH --error=./err-out/supervised-admixture-loo.%A_%a.err
#SBATCH --partition=cpu_compute
################
#SBATCH --time=00:00:00
#SBATCH --ntasks=8
#################
#SBATCH --array=1-17
#################

module load plink
module load admixture

# start and final are the indices of individuals to perform leave-one-out admixture on
start=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $1}' ./data/admixture/MO-method-start-end-array.txt)
final=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $2}' ./data/admixture/MO-method-start-end-array.txt)

genotypes=/lustrefs/nwrc/projects/DeSaix/MO-method/data/genotypes/MOMethod-SpatialPops1672x27866

# Final filename to store results
mo_final_loo=/lustrefs/nwrc/projects/DeSaix/MO-method/out/admixture/empirical/SpatialPops1672x27866-admixture-summary.txt

# Make tmp directory for intermediary files
tmpdir=/lustrefs/nwrc/projects/DeSaix/MO-method/out/admixture/tmp_array_${SLURM_ARRAY_TASK_ID}
mkdir -p ${tmpdir}
cd ${tmpdir}

# copy over data
cp ${genotypes}.bed ./Iteration.bed
cp ${genotypes}.bim ./Iteration.bim

k=$(cut -f1 -d' ' ${genotypes}.fam | sort | uniq | wc -l | sed 's/^[[:space:]]*//g')

echo ${start} ${final}
for i in $( seq $start $final )
do
    awk -v i=$i 'NR != i {print}; NR == i {print "-", $2, $3, $4, $5, $6}' ${genotypes}.fam > ./Iteration.fam
    cut -f1-2 -d' ' Iteration.fam > Iteration.pop
    admixture -j8 --supervised Iteration.bed $k

    #Find the row of the file the query individual is in
    IID=$(awk -v i=$i 'NR == i {print $2}' Iteration.fam)
    NR_query=$(awk -v IID=${IID} '$2==IID {print NR}' Iteration.fam)
    # Get the Q value
    Q_query=$(awk -v NR_query=${NR_query} 'NR == NR_query {print}' Iteration.${k}.Q)
    echo ${IID} ${Q_query} >> ${mo_final_loo}
done

rm -rf ${tmpdir}




