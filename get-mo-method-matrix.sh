#!/bin/bash
#set a job name
#SBATCH --job-name=supervised-admixture-matrix
#SBATCH --output=./err-out/supervised-admixture-matrix.%A_%a.out
#SBATCH --error=./err-out/supervised-admixture-matrix.%A_%a.err
#SBATCH --partition=cpu_compute
################
#SBATCH --time=00:00:00
#SBATCH --ntasks=8
#################
#SBATCH --array=1-7
#################

module load plink
module load admixture

# start and final are the indices of individuals to perform leave-one-out admixture on
start=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $1}' ./data/admixture/matrix-start-end-array.txt)
final=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $2}' ./data/admixture/matrix-start-end-array.txt)

genotypes=/lustrefs/nwrc/projects/DeSaix/MO-method/data/genotypes/MOMethod-SpatialPops1672x27866
matrix=/lustrefs/nwrc/projects/DeSaix/MO-method/data/genotypes/MOMethod-Matrix683x27866

# Final filename to store results
mo_final_matrix=/lustrefs/nwrc/projects/DeSaix/MO-method/out/admixture/empirical/Matrix683x27866-admixture-summary.txt

# Make tmp directory for intermediary files
tmpdir=/lustrefs/nwrc/projects/DeSaix/MO-method/out/admixture/tmp_Matrix_array_${SLURM_ARRAY_TASK_ID}
mkdir -p ${tmpdir}
cd ${tmpdir}

# copy over data
cp ${genotypes}.bed ./Iteration.bed
cp ${genotypes}.bim ./Iteration.bim

k=$(cut -f1 -d' ' ${genotypes}.fam | sort | uniq | wc -l | sed 's/^[[:space:]]*//g')

for i in $( seq $start $final )
do
    # split matrix individual i
    awk -v i=$i 'NR == i {print $1, $2}' ${matrix}.fam > Query_keep.txt
    plink --bfile ${matrix} --keep Query_keep.txt --make-bed --out Query
    plink --bfile ${genotypes} --bmerge Query --make-bed --out Iteration

    cut -f1-2 -d' ' Iteration.fam > Iteration.pop
    admixture -j8 --supervised Iteration.bed $k

    #Find the row of the file the query individual is in
    IID=$(awk '{print $2}' Query_keep.txt)
    NR_query=$(awk -v IID=${IID} '$2==IID {print NR}' Iteration.fam)
    # Get the Q value
    Q_query=$(awk -v NR_query=${NR_query} 'NR == NR_query {print}' Iteration.${k}.Q)
    echo ${IID} ${Q_query} >> ${mo_final_matrix}
    
    rm Query*
    rm Iteration*
done

rm -rf ${tmpdir}




