conda activate af
IDX_DIR=/reference/Mus_musculus/simpleaf

# Define filename pattern
reads1_pat="_R1_"
reads2_pat="_R2_"

for FASTQ_DIR in ../GEX/o3*/3*GEX*; do
s=`basename $FASTQ_DIR | cut -d"-" -f 2`
echo $s

# Obtain and sort filenames
reads1="$(find -L ${FASTQ_DIR} -name "*$reads1_pat*" -type f | sort | awk -v OFS=, '{$1=$1;print}' | paste -sd, -)"
reads2="$(find -L ${FASTQ_DIR} -name "*$reads2_pat*" -type f | sort | awk -v OFS=, '{$1=$1;print}' | paste -sd, -)"

mkdir -p $s

# simpleaf quantfication
simpleaf quant \
--reads1 $reads1 \
--reads2 $reads2 \
--threads 32 \
--index $IDX_DIR/index \
--chemistry 10xv3 --resolution cr-like \
--expected-ori fw \
-x /mnt/scratch/afhome/plist/737K-arc-v1.txt \
--t2g-map $IDX_DIR/index/t2g_3col.tsv \
--output $s

done
