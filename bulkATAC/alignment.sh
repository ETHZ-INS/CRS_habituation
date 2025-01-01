#!/bin/bash

ref=/reference/Mus_musculus/Ensembl/GRCm38/Sequence/BOWTIE2Index/genome
adapters=/reference/adapters/truseqPE.fa
mkdir -p trimmed
mkdir -p aligned
mkdir -p tracks

for f in X204SC23071098*/01.RawData/*; do

base=`basename $f`
bam=aligned/"$base".bam

echo $base
trimdir=trimmed

if [ -f "$bam" ]; then
    echo "$bam found; skipping"
else

if [ -f "$trimdir/"$base"_1.paired.fastq.gz" ]; then
  echo $trimdir/$base"_*fastq.gz found; skipping"
else
  trimmomatic PE -threads 8 -summary $trimdir/"$base".stats -phred33 $f/*"_1.fq.gz" $f/*"_2.fq.gz" $trimdir/"$base"_1.paired.fastq.gz $trimdir/"$base"_1.unpaired.fastq.gz $trimdir/"$base"_2.paired.fastq.gz $trimdir/"$base"_2.unpaired.fastq.gz ILLUMINACLIP:$adapters:2:15:4:4:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25
fi


(bowtie2 -p 12 --dovetail --no-mixed --no-discordant -I 15 -X 2000 -x $ref -1 $trimdir/"$base"_1.paired.fastq.gz -2 $trimdir/"$base"_2.paired.fastq.gz) 2> aligned/"$base".bowtie2 | samtools view -bS - | samtools sort -@4 -m 2G - > $bam

java -Xms6G -Xmx12G -jar /common/picard.jar MarkDuplicates I=$bam O=$base.bam.2 M=aligned/$base.picard.dupMetrics.txt && mv $base.bam.2 $bam

samtools index $bam && rm $trimdir/"$base"_*.fastq.gz
fi

done

