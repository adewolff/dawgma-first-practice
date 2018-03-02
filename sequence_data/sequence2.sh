for i in {1..24};
do
echo "Processing sample ${i}"
# Assuming salmon is one directory above the reads directory, ie
# in the directory where salmon resides there is a directory reads.
salmon quant -i athal_index -l A \
         -r ./reads/${i}_L001*/*_S*_L001_R1_001.fastq.gz \
            ./reads/${i}_L002*/*_S*_L002_R1_001.fastq.gz \
            ./reads/${i}_L003*/*_S*_L003_R1_001.fastq.gz \
            ./reads/${i}_L004*/*_S*_L004_R1_001.fastq.gz \
         -p 8 -o quants/${i}_quant
done 
