/usr/local/module/spartan_old.sh
module load Singularity
module load snakemake


## Step1
cd /data/gpfs/projects/punim0609/qian_feng/simulation/7_indel_rate/rate_0.1
for jump in $(seq 0 1 99);do mkdir replicate_$jump;done
cd ../
for jump in $(seq 0 1 99);do cp -r data rate_0.1/replicate_$jump;done
for jump in $(seq 0 1 99);do cp -r codes rate_0.1/replicate_$jump;done
for jump in $(seq 0 1 99);do cp Snakefile_indels rate_0.1/replicate_$jump;done
for jump in $(seq 0 1 99);do cp -r INDELibleV1.03/src rate_0.1/replicate_$jump;done
for jump in $(seq 0 1 99);do cp INDELibleV1.03/control.txt rate_0.1/replicate_$jump;done


## Step2
#!/bin/bash
# module load GSL/1.16-spartan_intel-2017.u2
for jump in $(seq 0 1 99); do
cd /data/gpfs/projects/punim0609/qian_feng/simulation/7_indel_rate/rate_0.1/replicate_${jump} && snakemake msprime_tree -s Snakefile_indels

done


## Step3
#!/bin/bash
cd ../../slurms
for jump in $(seq 0 1 99); do
echo sbatch -p mig --nodes=1 --job-name 7_0.1.${jump} --account "punim0609" --ntasks=1 --cpus-per-task=1 --mem=8192 --mail-type=FAIL --mail-user=fengq2@student.unimelb.edu.au --time=4-0:0:00 -e "slurm-%A_%a.out" --wrap=\"cd /data/gpfs/projects/punim0609/qian_feng/simulation/7_indel_rate/rate_0.1/replicate_${jump}\"

sbatch -p mig --nodes=1 --job-name 7_0.1.${jump} --account "punim0609" --ntasks=1 --cpus-per-task=1 --mem=8192 --mail-type=FAIL --mail-user=fengq2@student.unimelb.edu.au --time=4-0:0:00 -e \"slurm-%A_%a.out\" --wrap="cd /data/gpfs/projects/punim0609/qian_feng/simulation/7_indel_rate/rate_0.1/replicate_${jump} && nohup snakemake -d /data/gpfs/projects/punim0609/qian_feng/simulation/7_indel_rate/rate_0.1/replicate_${jump} -s Snakefile_indels"

sleep 1

done

#all jobs with 8192 mem and 4days
