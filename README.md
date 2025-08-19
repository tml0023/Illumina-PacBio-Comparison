# Illumina-PacBio-Comparison
In this repository I summarize the tools used to call multiple variant types for Illumina short-read data and PacBio HiFi read data for identical samples. Outputs for both datasets and multiple different bioinformatics tools are compared to analyze each platforms ability to identify variants, and how this can affect the identification of genetic variants in breast cancer susceptibility genes. 
___



# PacBio HiFi data processing with SMRTtools

This document shows how SMRT grant PacBio HiFi whole genome sequencing samples were procecced for CNV calling and short variant calling. I am attaching a link to the [SMRT_tools manual](chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.pacb.com/wp-content/uploads/SMRT_Tools_Reference_Guide_v11.0.pdf) that is a great help. 
___
## SMRTtools Long Variant Calling 

### Step1 Trim fastq files
HiFi data is very high quality and does not need a whole bunch of trimming but it is always helpful to take a look at the data before hand to see what exactly you want to trim. The raw data in this study was analyzed using FastQC/MultiQC.

```
#!/bin/sh

module load trimmomatic/0.39

ls *merged*fastq.gz > FSamplesList.txt

FILELIST=`cat FSamplesList.txt`

for FILENAME in $FILELIST
do

SHORT=`echo $FILENAME | awk -F "_" '{print $1}'`

java -jar /tools/trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 48 -phred33 -trimlog "$SHORT".trim.log "$SHORT"_merged.hifi_reads.fastq.gz "$SHORT"_merged.hifi_reads.trim.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 CROP:33000 HEADCROP:20 SLIDINGWINDOW:4:10

done
```

### Step2 Align HiFi reads to reference genome with pbmm2
PacBio recommends using the alignment tool pbmm2 for its long-read data. This tools is available in SMRTtools provided by PacBio. 

```
#!/bin/sh

FILELIST=`cat FSamplesList.txt`

for FILE in $FILELIST; do 

SHORT=`echo $FILE | awk -F "_" '{print $1}'`

/tools/smrtlink/new/bundles/smrttools/smrtcmds/bin/pbmm2 align /genome.fa "$SHORT"_merged.hifi_reads.trim.fastq.gz "$SHORT".hifi_reads.bam --sort --preset CCS --sample "$SHORT" --rg '@RG\tID:"$SHORT"'

done
```

### Step3 Discover PacBio Structural Variants  (pbsv) 
pbsv identifies the following types of variants: Insertions, deletions, duplications, copy number variants, inversions, and translocations and outputs into a svsig.gz file.

```
#!/bin/sh

FILELIST=`cat FSamplesList.txt`

for FILE in $FILELIST; do 

SHORT=`echo $FILE | awk -F "_" '{print $1}'`

/tools/smrtlink/new/bundles/smrttools/smrtcmds/bin/pbsv discover "$SHORT".hifi_reads.bam "$SHORT".hifi_read.svsig.gz

done
```

### Step4 Call PacBio Structural Variants (pbsv)
pbsv is now used to call all structural variants that were discovered in the last step and organizes them into a vcf file. 

```
#!/bin/sh

FILELIST=`cat FSamplesList.txt`

for FILE in $FILELIST; do 

SHORT=`echo $FILE | awk -F "_" '{print $1}'`

/tools/smrtlink/new/bundles/smrttools/smrtcmds/bin/pbsv call --ccs /genome.fa "$SHORT".hifi_read.svsig.gz "$SHORT".hifi_read.vcf

done
```

### Step5 Annotate variants with svpack
This program is not included from with SMRT_tools and must be downloaded from the [Github page](https://github.com/PacificBiosciences/svpack). This tool is used to annotate structural variants. 

```
#!/bin/sh

for FILE in *vcf; do 

svpack-main/svpack filter --pass-only --min-svlen 50 $FILE |
svpack-main/svpack consequence --require-csq - ensembl.GRCh38.101.reformatted.gff3 > Annotated_"$FILE"

done
```

___
## SMRTtools short variant calling


### Step1 Trim fastq files
HiFi data is very high quality and does not need a whole bunch of trimming but it is always helpful to take a look at the data before hand to see what exactly you want to trim. 

```
#!/bin/sh

module load trimmomatic/0.39

ls *merged*fastq.gz > FSamplesList.txt

FILELIST=`cat FSamplesList.txt`

for FILENAME in $FILELIST
do

SHORT=`echo $FILENAME | awk -F "_" '{print $1}'`

java -jar /tools/trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 48 -phred33 -trimlog "$SHORT".trim.log "$SHORT"_merged.hifi_reads.fastq.gz "$SHORT"_merged.hifi_reads.trim.fastq.gz ILLUMINACLIP:/hosted/cvmpt/archive/WGS_Human/WGS3_Dec2022_TL/TruSeq3-PE.fa:2:30:10 CROP:33000 HEADCROP:20 SLIDINGWINDOW:4:10

done
```

### Step2 Align HiFi reads to reference genome with minimap2 (pbmm2)
PacBio recommends using the alignment tool minimap2 for its long-read data. This tools is available in SMRTtools provided by PacBio. 

```
#!/bin/sh

FILELIST=`cat FSamplesList.txt`

for FILE in $FILELIST; do 

SHORT=`echo $FILE | awk -F "_" '{print $1}'`

/tools/smrtlink/new/bundles/smrttools/smrtcmds/bin/pbmm2 align /genome.fa "$SHORT"_merged.hifi_reads.trim.fastq.gz "$SHORT".hifi_reads.bam --sort --preset CCS --sample "$SHORT" --rg '@RG\tID:"$SHORT"'

done
```

### Step3 Use DeepVariant to call short variants
Deepvariant is a program that is not included in SMRTtools. However this is the program recommended by Pacbio for short variant calling. This program can be downloaded from the [Githhub page](https://github.com/google/deepvariant). There is also an example of deepvariant located [here](https://github.com/google/deepvariant/blob/r1.0/docs/deepvariant-pacbio-model-case-study.md).

```
#!/bin/sh

module load singularity

BIN_VERSION="1.0.0"
singularity exec --bind ${PATH} \
docker://google/deepvariant:${BIN_VERSION} \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type PACBIO \
  --ref genome.fa \
  --reads BC-EAMC-209-1.hifi_reads.bam \
  --output_vcf BC-EAMC-209-1.hifi_reads_deepvariant.vcf.gz \
  --num_shards ${nproc}\
```
___

#  Illumina and PacBio data processing using the Genome Analysis ToolKit (GATK)
This pipeline was used for calling short variants in the datasets provided by both sequencing platforms.

You notice that the Illumina and PacBio trimming parameters are different. These trimming scripts follow recommended strategies for maximizing data quality for each of the sequencing platforms. Quality outputs for trimmed Illumina and PacBio data can be viewed in this repository.
Trimming parameters are the only different parameters between Illumina and PacBio data in this GATK pipeline in effort to keep the data processing as similar as possible. The only difference between all other scripts is the naming of the files. Note file names are indicative of the PacBio data. Illumina scripts are the exact same other than the naming of the files. 

### Step1A Trimming of Illumina fastq files 
```
#!/bin/sh

module load trimmomatic/0.39

FILELIST=`cat FSamplesList.txt`

for FILENAME in $FILELIST
do

SHORT=`echo $FILENAME | awk -F "_" '{print $1}'`

java -jar /tools/trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 48 -phred33 -trimlog $SHORTER.trim.log "$SHORT"_1_Apr2022.fastq.gz "$SHORT"_2_Apr2022.fastq.gz -baseout $SHORT.trim.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 HEADCROP:0 LEADING:0 TRAILING:0 SLIDINGWINDOW:4:10

done
```
### Step1B Trimming of PacBio fastq files
```
#!/bin/sh

module load trimmomatic/0.39

FILELIST=`cat FSamplesList.txt`

for FILENAME in $FILELIST
do

SHORT=`echo $FILENAME | awk -F "_" '{print $1}'`

java -jar /tools/trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 48 -phred33 -trimlog "$SHORT".trim.log "$SHORT"_merged.hifi_reads.fastq.gz "$SHORT"_merged.hifi_reads.trim.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 CROP:33000 HEADCROP:20 SLIDINGWINDOW:4:10

done
```
### Step2 Align reads to reference using Burrows-Wheeler Aligner 
```
#!/bin/sh

module load bwa/2.0
module load samtools/1.19

FILELIST=`cat FSamplesList.txt`
for FILE in $FILELIST; do 

SHORT=`echo $FILE | awk -F "_" '{print $1}'`

bwa mem -t 20 -M /hosted/cvmpt/archive/Human_Genome/genome "$SHORT".trim_1P.fastq \
        "$SHORT".trim_2P.fastq > "$SHORT".mem.sam

samtools view -Sb -@ 20 "$SHORT".mem.sam -o "$SHORT".mem.bam
samtools sort -@ 20 "$SHORT".mem.bam -o "$SHORT".memsorted.bam
samtools index -@ 20 "$SHORT".memsorted.bam

done
```
### Step3 Mark Duplicate reads
```
#!/bin/sh

module load picard/2.23.9 
module load samtools/1.19

FILELIST=`cat FSamplesList.txt`
for FILE in $FILELIST; do

SHORT=`echo $FILE | awk -F "_" '{print $1}'`     

java -jar /tools/picard-2.23.9/libs/picard.jar MarkDuplicates \
      I="$SHORT"_merged.hifi_reads.memsorted.bam \
      O="$SHORT"_merged.hifi_reads.markdup.bam \
      M="$SHORT".marked_dup_metrics.txt

samtools sort -@ 48 "$SHORT"_merged.hifi_reads.markdup.bam -o "$SHORT"_merged.hifi_reads.markdup.sorted.bam
samtools index -@ 48 "$SHORT"_merged.hifi_reads.markdup.sorted.bam
done
```
### Step4 Adding or replacing read groups 
Note the read groups for the Illumina script are changed to represent the Illumina platform used for sequencing.
```
#!/bin/sh

module load samtools/1.19
module load picard/2.23.9 

FILELIST=`cat FSamplesList.txt`
for FILE in $FILELIST; do

SHORT=`echo $FILE | awk -F "_" '{print $1}'`

 java -jar /tools/picard-2.23.9/libs/picard.jar AddOrReplaceReadGroups \
       I="$SHORT"_merged.hifi_reads.markdup.sorted.bam \
       O="$SHORT"_merged.hifi_reads.rg.bam \
       RGID="$SHORT"_PacBio \
       RGLB=Hfi_reads \
       RGPL=PacBio \
       RGPU=SeqII \
       RGSM="$SHORT"

samtools sort -@ 48 "$SHORT"_merged.hifi_reads.rg.bam -o "$SHORT"_merged.hifi_reads.rgsorted.bam
samtools index -@ 48 "$SHORT"_merged.hifi_reads.rgsorted.bam

done
```
### Step5 Base recalibrator 
```
#!/bin/sh

module load gatk/4.1.9.0
module load samtools/1.19
module load R/4.3.2

FILELIST=`cat FSamplesList.txt`  ##Can be used if a file list is needed
for FILE in $FILELIST; do

SHORT=`echo $FILE | awk -F "_" '{print $1}'`     

 gatk BaseRecalibrator --maximum-cycle-value 33000 \
   -I "$SHORT"_merged.hifi_reads.rgsorted.bam \
   -R /genome.fa \
   --known-sites /dbSNP150.hg38.vcf \
   --known-sites /Mills_and_1000G_gold_standard.indels.hg38.vcf \
   --known-sites /Homo_sapiens_assembly38.known_indels.vcf \
   --known-sites /hapmap_3.3.hg38.vcf \
   -O "$SHORT"_merged.hifi_reads.recal_data.table
 
 gatk ApplyBQSR \
   -R /hosted/cvmpt/archive/Human_Genome/genome.fa \
   -I "$SHORT"_merged.hifi_reads.rgsorted.bam \
   --bqsr-recal-file "$SHORT"_merged.hifi_reads.recal_data.table \
   -O "$SHORT"_merged.hifi_reads.recal.bam
 
 gatk AnalyzeCovariates \
   -bqsr "$SHORT"_merged.hifi_reads.recal_data.table \
   -plots "$SHORT"_merged.hifi_reads.AnalyzeCovariates.pdf

samtools sort -@ 48 "$SHORT"_merged.hifi_reads.recal.bam -o "$SHORT"_merged.hifi_reads.recal.sorted.bam
samtools index -@ 48 "$SHORT"_merged.hifi_reads.recal.sorted.bam

done
```
### Step6 Haplotype caller
```
#!/bin/sh

module load gatk/4.1.9.0

FILELIST=`cat FSamplesList.txt`  ##Can be used if a file list is needed
for FILE in $FILELIST; do

SHORT=`echo $FILE | awk -F "_" '{print $1}'`     

 gatk --java-options "-Xmx200g" HaplotypeCaller  \
   -R /hosted/cvmpt/archive/Human_Genome/genome.fa \
   -I "$SHORT"_merged.hifi_reads.recal.sorted.bam \
   -O "$SHORT"_Oct2023_merged.hifi_reads.vcf.gz \
   -A AlleleFraction \
   -A BaseQuality \
   -A MappingQuality \
   --native-pair-hmm-threads 48 

done
```
### Step7 Genomic database import
```
#!/bin/sh

module load picard/2.23.9
module load gatk/4.1.9.0

gatk --java-options "-Xmx550g -Xms550g" GenomicsDBImport \
-V BC-CR-50-1_merged.hifi_reads.vcf.gz \
-V BC-EAMC-130-1_merged.hifi_reads.vcf.gz \
-V BC-EAMC-209-1_merged.hifi_reads.vcf.gz \
-L wgs_calling_regions.hg38.interval_list \
--genomicsdb-workspace-path PacBio_GDBI_Workspace \
--tmp-dir PacBio_GDBI_tmp \
--reader-threads 48 \
```
### Step8 GenotypeGVCF
```
#!/bin/sh

module load gatk 

 gatk --java-options "-Xmx550g" GenotypeGVCFs \
   -R /hosted/cvmpt/archive/Human_Genome/genome.fa \
   -V gendb://PacBio_GDBI_Workspace \
   -O GenotypeGVCF_output_PacBio.g.vcf.gz \
   --tmp-dir GenotypeGVCF_PacBio_tmp
```
### Step9 Variant recalibrator 
```
#!/bin/sh

module load gatk/4.1.9.0
module load R/4.3.2

## recalibrate SNPs
 gatk VariantRecalibrator \
   -R /hosted/cvmpt/archive/Human_Genome/genome.fa \
   -V GenotypeGVCF_output_PacBio.g.vcf.gz \
   --resource:1000G,known=false,training=true,truth=false,prior=10.0 /hosted/cvmpt/Human_Research/KnownSites/Mills_and_1000G_gold_standard.indels.hg38.vcf \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /hosted/cvmpt/Human_Research/KnownSites/dbSNP150.hg38.vcf \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /hosted/cvmpt/Human_Research/KnownSites/hapmap_3.3.hg38.vcf \
   -an DP -an QD -an FS -an SOR -an MQ -an ReadPosRankSum -an MQRankSum \
   -mode SNP \
   -L /hosted/cvmpt/Human_Research/KnownSites/wgs_calling_regions.hg38.interval_list \
   -O VQSR_SNP_PacBio_Output.vcf.recal \
   --tranches-file VQSR_SNP_PacBio_Output.vcf.tranches \
   --rscript-file VQSR_SNP_PacBio_Output.vcf.plots.R

 gatk ApplyVQSR \
   -R /hosted/cvmpt/archive/Human_Genome/genome.fa \
   -V GenotypeGVCF_output_PacBio.g.vcf.gz \
   -O ApplyVQSR_SNPs_PacBio_Output.vcf.gz \
   --truth-sensitivity-filter-level 90.0 \
   --tranches-file VQSR_SNP_PacBio_Output.vcf.tranches \
   --recal-file VQSR_SNP_PacBio_Output.vcf.recal \
   -mode SNP

## reaclibrate for INDELs
 gatk VariantRecalibrator \
   -R /hosted/cvmpt/archive/Human_Genome/genome.fa \
   -V ApplyVQSR_SNPs_PacBio_Output.vcf.gz \
   --resource:1000G,known=false,training=true,truth=false,prior=10.0 /hosted/cvmpt/Human_Research/KnownSites/Mills_and_1000G_gold_standard.indels.hg38.vcf \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /hosted/cvmpt/Human_Research/KnownSites/dbSNP150.hg38.vcf \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /hosted/cvmpt/Human_Research/KnownSites/hapmap_3.3.hg38.vcf \
   -an DP -an QD -an FS -an SOR -an MQ -an ReadPosRankSum -an MQRankSum \
   -mode INDEL \
   -L /hosted/cvmpt/Human_Research/KnownSites/wgs_calling_regions.hg38.interval_list \
   -O ApplyVQSR_INDEL_SNP_PacBio_Output.vcf.recal \
   --tranches-file VQSR_INDEL_SNP_PacBio_Output.vcf.tranches \
   --rscript-file VQSR_INDEL_SNP_PacBio_Output.vcf.plots.R

 gatk ApplyVQSR \
   -R /hosted/cvmpt/archive/Human_Genome/genome.fa \
   -V ApplyVQSR_SNPs_PacBio_Output.vcf.gz \
   -O ApplyVQSR_INDEL_SNP_PacBio_Output.vcf.gz \
   --truth-sensitivity-filter-level 90.0 \
   --tranches-file VQSR_INDEL_SNP_PacBio_Output.vcf.tranches \
   --recal-file ApplyVQSR_INDEL_SNP_PacBio_Output.vcf.recal \
   -mode INDEL
```
# Minimap2 Alignment 
Minimap2 was used to align both PacBio and Illumina data. Comparisons on Minimap2 alignment were ran against BWA-mem and pbmm2 alignments using GATKs DepthOfCoverage and Samtools Flagstat.

### PacBio alignment 
```
#!/bin/sh

module load minimap2/2.26

minimap2 -ax map-hifi ref.fa Input_file.fastq > Output_file.sam
```

### Illumina alignment
```
#!/bin/sh

module load minimap2/2.26

minimap2 -ax sr ref.fa Input_file_1.fastq Input_file_2.fastq > Output_file.sam
```























