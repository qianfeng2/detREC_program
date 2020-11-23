/usr/local/module/spartan_old.sh
module load Singularity
module load snakemake


## Step1
cd /data/gpfs/projects/punim0609/qian_feng/simulation/8_indel_size/size_3.7
for jump in $(seq 0 1 99);do mkdir replicate_$jump;done
cd ../
for jump in $(seq 0 1 99);do cp -r data size_3.7/replicate_$jump;done
for jump in $(seq 0 1 99);do cp -r codes size_3.7/replicate_$jump;done
for jump in $(seq 0 1 99);do cp Snakefile_indels size_3.7/replicate_$jump;done
for jump in $(seq 0 1 99);do cp -r INDELibleV1.03/src size_3.7/replicate_$jump;done
for jump in $(seq 0 1 99);do cp INDELibleV1.03/control.txt size_3.7/replicate_$jump;done


## Step2
#!/bin/bash
# module load GSL/1.16-spartan_intel-2017.u2
for jump in $(seq 0 1 99); do
cd /data/gpfs/projects/punim0609/qian_feng/simulation/8_indel_size/size_3.7/replicate_${jump} && snakemake msprime_tree -s Snakefile_indels

done


## Step3
#!/bin/bash
cd ../../slurms
for jump in $(seq 0 1 99); do
echo sbatch -p mig --nodes=1 --job-name 8_3.7_${jump} --account "punim0609" --ntasks=1 --cpus-per-task=1 --mem=8192 --mail-type=FAIL --mail-user=fengq2@student.unimelb.edu.au --time=4-0:0:00 -e "slurm-%A_%a.out" --wrap=\"cd /data/gpfs/projects/punim0609/qian_feng/simulation/8_indel_size/size_3.7/replicate_${jump}\"

sbatch -p mig --nodes=1 --job-name 8_3.7_${jump} --account "punim0609" --ntasks=1 --cpus-per-task=1 --mem=8192 --mail-type=FAIL --mail-user=fengq2@student.unimelb.edu.au --time=4-0:0:00 -e \"slurm-%A_%a.out\" --wrap="cd /data/gpfs/projects/punim0609/qian_feng/simulation/8_indel_size/size_3.7/replicate_${jump} && nohup snakemake -d /data/gpfs/projects/punim0609/qian_feng/simulation/8_indel_size/size_3.7/replicate_${jump} -s Snakefile_indels"

sleep 1

done

#all jobs with 8192 mem and 4days
