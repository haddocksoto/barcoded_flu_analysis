#!/bin/sh


####################################################################################
#### FASTQ files for reads 1 and 2 should be in a directory called "fastq_files" ###
####################################################################################


#First we generate the appropiate directory where all the data output will be stored

mkdir data_output 

####################################################################################

# cd to a directory above your raw_fastq files and name it 'fastq_files'


path_with_reads='fastq_files/' #provide the directory where all of the R1 and R2 FASTQ files are located

cd $path_with_reads 

for reads in *_R1.fastq.gz; 
do
sample=${reads%%_R1.fastq.gz}
mkdir data_output/${sample}                                         #create one subdirectory for each sample
bbmerge.sh in1=${sample}_R1.fastq.gz in2=${sample}_R2.fastq.gz \    #use bbmerge to trim non overlapping regions
  out=../data_output/${sample}/${sample}_R1_merged.fastq \
    out2=../data_output/${sample}/${sample}_R2_merged.fastq \
      outu=../data_output/${sample}/${sample}_R1_unmerged.fastq \
        outu2=../data_output/${sample}/${sample}_R2_unmerged.fastq \
          minoverlap=30 ecco=t mix=f
done 


####################################################################################

path='../data_output/'
cd $path
#now we are working on the data_output directory, each containinig a subdirectory for each sample

####################################################################################

#specify path to reference 
Ref="PATH/TO/REFERENCE/fasta_file.fa" 

#Index reference file with BWA

bwa index $Ref

#Map to reference using BWA mem

for file in */*_R1_merged.fastq.gz; #still working on data_output
do
sample=${file%%_R1_merged.fastq.gz}
echo ${sample} # to make sure we are working with the correct samples
bwa mem -t 8 -B 10 $Nhe_ref \
  ${sample}_R1_merged.fastq.gz \ #reading R1
    ${sample}_R2 _merged.fastq.gz \ #reading R2
      > ${sample}_mapped.sam #creating SAM file
 done
 
 # We use Samtools to align the BAM file

for file in */*_mapped.sam; 
do
sample=${file%%_mapped.sam}
echo ${sample} # to make sure we are working with the correct samples
samtools sort -o \
  ${sample}_mapped.bam \
    ${sample}_mapped.sam
done

# We use LoFreq (v2.1.5) to call variants

for file in */*_mapped.bam; 
do
sample=${file%%_mapped.bam}
echo ${sample} # to make sure we are working with the correct samples
lofreq call -f $Ref -q 30 -Q 30 -C 500  --force-overwrite \ #make sure to add path to reference
  -o  ${sample}.vcf \
  ${sample}_mapped.bam 
done

# variants were annotated using SnpEff. 
for file in */*.vcf; 
do
sample=${file%%.vcf}
echo ${sample}
java -Xmx8g -jar /Users/PATH_TO_snpeff/snpEff/snpEff.jar \ #download SNPeff at https://pcingola.github.io/SnpEff/
  ${sample}.vcf > ${sample}.ann.vcf
done

#We used bgzip and tabix to compress the vcf file and make it easier to work with

for file in */*.ann.vcf; 
do
sample=${file%%.ann.vcf}
echo ${sample}
bgzip -f -c ${sample}.ann.vcf > ${sample}${sample}.ann.vcf.gz
tabix ${sample}.ann.vcf.gz
done

# Now we use bcftools to extract columns from the vcf file to make it easier to work with

for file in */*.ann.vcf.gz; 
do
sample=${file%%.ann.vcf.gz}
echo ${sample}
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/AF\t%INFO/ANN\n' \
  ${sample}.ann.vcf.gz > ${sample}_simplified.vcf
done

#From now on we used python code to clean up columns and extract variant frequencies, locations, ect.






