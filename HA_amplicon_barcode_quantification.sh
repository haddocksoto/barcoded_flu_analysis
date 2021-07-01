#!/bin/sh

#generate the appropiate directories

mkdir data_output 

#merging paired end reads

path_with_reads='fastq_files/' #provide the directory where all of the R1 and R2 FASTQ files are located

cd $path_with_reads 

for reads in *_R1.fastq.gz; 
do
sample=${reads%%_R1.fastq.gz}
mkdir data_output/${sample} #create one subdirectory for each sample
bbmerge.sh in1=${sample}_R1.fastq.gz in2=${sample}_R2.fastq.gz \
out=../data_output/${sample}/${sample}.fastq
done


#MAPPING TO REFERENCE FASTA FILE

path='../data_output/'
cd $path
#now we are working on the data_output directory, each containinig a 
#subdirectory for each sample

#specify path to reference and path to merged fastq_files
HA_ref="../references/bcCA07_HA_NNN_amplicon.fasta"
 

bbmap.sh ref=$HA_ref

for file in */*.fastq; 
do
sample=${file%%.fastq}
echo ${sample}
bbmap.sh in=${sample}.fastq outm=${sample}_mapped.sam -Xmx20g \
outu=${sample}_unmapped.sam pairedonly=t printunmappedcount=t
done


#TRIMMING READS TO AVERAGE INSERT SIZE OBTAINED WHEN PAIRING

for file in */*_mapped*; 
do
sample=${file%%._mapped.sam}
echo ${sample}
reformat.sh in=${sample} out=${sample}_trimmed.sam minlength=227 maxlength=227 overwrite=true
    for f in */*mapped.sam_trimmed.sam; 
    do 
    mv "$f" "${f%_m*trimmed.sam}_trimmed.sam" 
    done
done


#ALIGNING WITH SAMTOOLS

for file in */*_trimmed*; 
do
sample=${file%%._trimmed.sam}
echo ${sample}
samtools view -hSbo \
${sample}_aligned.bam \
${sample}
    for f in */*trimmed.sam_aligned.bam; 
    do 
    mv "$f" "${f%_t*.bam}_aligned.bam" 
    done
done


#SORTING WITH SAMTOOLS

for file in */*_aligned*; 
do
sample=${file%%._aligned.bam}
echo ${sample}
samtools sort -o \
${sample}_aligned_sorted.bam \
${sample}
    for f in */*aligned.bam_aligned_sorted.bam; 
    do 
    mv "$f" "${f%_aligned.bam_aligned_sorted.bam}_sorted.bam" 
    done
done

#INDEXING WITH SAMTOOLS

for file in */*_sorted*; 
do
sample=${file%%._sorted.bam}
echo ${sample}
samtools index ${sample}
done

########### GENERATING FASTQ FILES USING ALIGNED BAM FILES GENERATED ABOVE ##########

#CREATE FASTQ

for file in */*_sorted.bam; 
do
sample=${file%%._sorted.bam}
echo ${sample}
samtools fastq ${sample} > ${sample}_temp.fastq
    for f in */*sorted.bam_temp.fastq; 
    do 
    mv "$f" "${f%_sorted.bam_temp.fastq}_processed_temp.fastq" 
    done
done

#READS INTO 5-3 DIRECTION

for file in */*processed_temp.fastq; 
do
sample=${file%%.processed_temp.fastq}
echo ${sample}
reformat.sh in=${sample} out=${sample}_5_to_3_temp.fastq rcomp
    for f in */*_temp.fastq_5_to_3_temp.fastq; 
    do 
    mv "$f" "${f%_temp.fastq_5_to_3_temp.fastq}_5_to_3_temp.fastq" 
    done
done

#GENERATING FASTA FILE WITH ALL OF THE READS SO THAT WE CAN EXTRACT BARCODE REGION

for file in */*processed_5_to_3_temp.fastq; 
do
sample=${file%%.processed_5_to_3_temp.fastq}
echo ${sample}
seqtk seq -a ${sample} > ${sample}.fasta
    for f in */*processed_5_to_3_temp.fastq.fasta; 
    do 
    mv "$f" "${f%processed_5_to_3_temp.fastq.fasta}processed_5_to_3_temp.fasta" 
    done
done

#removing special characters

for file in */*processed_5_to_3_temp.fasta;
do
sample=${file%%.processed_5_to_3_temp.fasta}
echo ${sample}
sed -i '' '/>/d' ${sample}
sed -i '' '/+/d' ${sample}
done

#EXTRACTING BARCODES FROM HEMAGGLUTININ (HA) READS ONLY

for file in */*processed_5_to_3_temp.fasta;
do
sample=${file%%.processed_5_to_3_temp.fasta}
echo ${sample}
sed 's/^[A-Z]\{154\}//' ${sample} > ${sample}_left_154.fasta
sed 's/.\{63\}$//' ${sample}_left_154.fasta > ${sample}_barcodes_temp.fasta
    for f in */*_temp.fasta_barcodes_temp.fasta; 
    do 
    mv "$f" "${f%_temp.fasta_barcodes_temp.fasta}_barcodes_temp.fasta" 
    done
    for f in */*_temp.fasta_left_154.fasta; 
    do 
    mv "$f" "${f%_temp.fasta_left_154.fasta}_left_154.fasta" 
    done
done


#Clenaing reads to extract only the reads with 10 nucleotides (barcodes)

for file in */*_barcodes_temp.fasta;
do
sample=${file%%.barcodes_temp.fasta}
echo ${sample}
sed -i '' '/>/d' ${sample}
sed -i '' '/N/d' ${sample}
sed -i '' '/-/d' ${sample}
sed -i '' '/?/d' ${sample}
grep -x '.\{10,10\}' ${sample} > ${sample}_clean.txt
    for f in */*barcodes_temp.fasta_clean.txt; 
    do 
    mv "$f" "${f%barcodes_temp.fasta_clean.txt}_barcodes_clean.txt" 
    done
done
#output: txt file with all of our barcodes as a single column.
