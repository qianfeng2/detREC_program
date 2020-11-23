/usr/local/module/spartan_old.sh
module load Singularity
module load snakemake


## Step1
cd /data/gpfs/projects/punim0609/qian_feng/simulation/6_AA_model/MTMAM
for jump in $(seq 0 1 99);do mkdir replicate_$jump;done
cd ../
for jump in $(seq 0 1 99);do cp -r data MTMAM/replicate_$jump;done
for jump in $(seq 0 1 99);do cp -r codes MTMAM/replicate_$jump;done
for jump in $(seq 0 1 99);do cp Snakefile MTMAM/replicate_$jump;done
#for jump in $(seq 0 1 99);do cp -r /data/gpfs/projects/punim0609/qian_feng/snake_pipeline/INDELibleV1.03/src /data/gpfs/projects/punim0609/qian_feng/simulation/2_average_number_recombination/replicate_$jump;done
#for jump in $(seq 0 1 99);do cp -r /data/gpfs/projects/punim0609/qian_feng/snake_pipeline/INDELibleV1.03/control.txt /data/gpfs/projects/punim0609/qian_feng/simulation/2_average_number_recombination/replicate_$jump;done


## Step2
#!/bin/bash
# module load GSL/1.16-spartan_intel-2017.u2
for jump in $(seq 0 1 99); do
cd /data/gpfs/projects/punim0609/qian_feng/simulation/6_AA_model/MTMAM/replicate_${jump} && snakemake msprime_tree -s Snakefile

done


## Step3
#!/bin/bash
cd ../../slurms
for jump in $(seq 0 1 99); do
echo sbatch -p mig --nodes=1 --job-name 6.MT.${jump} --account "punim0609" --ntasks=1 --cpus-per-task=1 --mem=8192 --mail-type=FAIL --mail-user=fengq2@student.unimelb.edu.au --time=8-0:0:00 -e "slurm-%A_%a.out" --wrap=\"cd /data/gpfs/projects/punim0609/qian_feng/simulation/6_AA_model/MTMAM/replicate_${jump}\"

sbatch -p mig --nodes=1 --job-name 6.MT.${jump} --account "punim0609" --ntasks=1 --cpus-per-task=1 --mem=8192 --mail-type=FAIL --mail-user=fengq2@student.unimelb.edu.au --time=8-0:0:00 -e \"slurm-%A_%a.out\" --wrap="cd /data/gpfs/projects/punim0609/qian_feng/simulation/6_AA_model/MTMAM/replicate_${jump} && nohup snakemake -d /data/gpfs/projects/punim0609/qian_feng/simulation/6_AA_model/MTMAM/replicate_${jump} -s Snakefile"

sleep 1

done

