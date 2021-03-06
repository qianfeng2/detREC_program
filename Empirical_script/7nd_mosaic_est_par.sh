for FILE in *.fasta; do
echo sbatch -p sysgen --job-name qian.${FILE} --account SG0011 --ntasks=1 --cpus-per-task=1 --mem=25600 --mail-user=fengq2@student.unimelb.edu.au --time=0-2:0:00 -e "slurm-%A_%a.out" --wrap=\"/vlsci/SG0011/qian-feng/MZmosaic/mosaic -seq ${FILE} -aa -psum -group 2 db target -target target -del 0.0081825316222 -eps 0.219652036378 -rec 0.0 > /vlsci/SG0011/qian-feng/UniMelb_shared-master/project/mosaic_processed_data/results_full_7/${FILE}_output.log\"

sleep 1

sbatch -p sysgen --job-name qian.${FILE} --account SG0011 --ntasks=1 --cpus-per-task=1 --mem=25600 --mail-user=fengq2@student.unimelb.edu.au --time=0-2:0:00 -e \"slurm-%A_%a.out\" --wrap="/vlsci/SG0011/qian-feng/MZmosaic/mosaic -seq ${FILE} -aa -psum -group 2 db target -target target -del 0.0081825316222 -eps 0.219652036378 -rec 0.0 > /vlsci/SG0011/qian-feng/UniMelb_shared-master/project/mosaic_processed_data/results_full_7/${FILE}_output.log"

echo "Job submitted!\n"

sleep 1
done

echo "All jobs submitted!\n"