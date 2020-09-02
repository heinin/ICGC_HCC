for SAMPLE in $(cat all_bams_for_rseqc.txt); 
do 
  sbatch --export=SAMPLE=$SAMPLE rseqc_inferexperiment.sh
  sleep 5
done;
