#!/bin/sh

for FILE in ../*.fasta; do
echo sbatch -p sysgen --job-name qian.${FILE} --account SG0011 --ntasks=1 --cpus-per-task=1 --mem=76800 --mail-user=fengq2@student.unimelb.edu.au --time=0-6:0:00 -e "slurm-%A_%a.out" --wrap=\"/vlsci/SG0011/qian-feng/MZmosaic/mosaic -ma -seq ${FILE} -aa -tag ${FILE} -group 2 db target -target target -del 0.00806934718714 -eps 0.2283998284 -rec 0.015\"

sleep 1

sbatch -p sysgen --job-name qian.${FILE} --account SG0011 --ntasks=1 --cpus-per-task=1 --mem=76800 --mail-user=fengq2@student.unimelb.edu.au --time=0-6:0:00 -e \"slurm-%A_%a.out\" --wrap="/vlsci/SG0011/qian-feng/MZmosaic/mosaic -ma -seq ${FILE} -aa -tag ${FILE} -group 2 db target -target target -del 0.00806934718714 -eps 0.2283998284 -rec 0.015"

echo "Job submitted!\n"

sleep 1
done

echo "All jobs submitted!\n"
