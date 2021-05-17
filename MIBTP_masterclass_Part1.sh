#################################################################################################################################
# DNaseI-Seq workshop part 1 - Read alignment and identification of open chromatin regions
#################################################################################################################################

#
# In this practical we will be processing two DNaseI-Seq samples obtained from mouse T-Cells
# Here we have two types of T-Cells (denoted as NAg and TAg) which are:
# NAg - Naive T-Cells which have been stimulted with antigen
# TAg - Tolerant T-Cells which have been stimulated by antigen
#
# By comparing the DHS profiles of these two samples we can show that the epigenetic response to stimulation with antigen of these two cell types is altered
#
# This data is published and is described in more detail in:
# Bevington et al. (2020). Chromatin Priming Renders T Cell Tolerance-Associated Genes Sensitive to Activation below the Signaling Threshold for Immune Response Genes. Cell Reports. 31(10), 107748
#
# In order to save time we will only be processing data from mouse chromosome 1
# In part 2 of this workshop you will work with processed data from these samples that contain data for the entire genome as well as additional bioloical replicates
#

#################################################################################################################################
# Step 0). Load the required software into your bluebear session
#################################################################################################################################

module load FastQC/0.11.9-Java-11
module load Trimmomatic/0.39-Java-11
module load Bowtie2/2.3.5.1-GCC-8.3.0
module load SAMtools/1.10-GCC-8.3.0
module load picard/2.21.1-Java-11
module load MACS2/2.2.7.1-foss-2019b-Python-3.7.4
module load BEDTools/2.29.2-GCC-8.3.0
module load Subread/2.0.1-GCC-8.3.0

#################################################################################################################################
# Step 1). Quality assessment and read trimming
#################################################################################################################################

# As with most next-generation sequencing projects we will begin with raw un-processed reads in the fastQ format
# These files can be found in:
# <PATH TO FASTQ FILES>

# To inspect the read quality we will use a program called fastqc
# Run this command for both samples

fastqc NAg_Rep1_sample.fastq.gz
fastqc TAg_Rep1_sample.fastq.gz

# This creates a html report that shows various read quality stastics - This will be shown in the workshop slides

# The quality of these reads is good - however we can still trim the reads to resolve any potential quality issues
# This trimming will carry out a number steps:
# 1). Removal of sequencing adaptors
# 2). Removal of low-quality bases from sequence reads (which might contain sequencing errors)
# 3). Any read that is too short after adaptor and low-quality base removal is filtered out of the fastq file
# 
# To do this we will use a tool called trimmomatic

# Set the path to the illumina adaptors as a variable
ADAPTORS=/rds/projects/k/knowletj-mibtp-masterclasses/MITBP_Masterclass_2021/data_for_workshop/Annotation/Trimmomatic_adapters/TruSeq3-SE.fa

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE NAg_Rep1_sample.fastq.gz NAg_Rep1_sample_trimmed.fastq.gz ILLUMINACLIP:$ADAPTORS:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE TAg_Rep1_sample.fastq.gz TAg_Rep1_sample_trimmed.fastq.gz ILLUMINACLIP:$ADAPTORS:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

# Let's have a closer look at what each of the paramters of this command means
# trimmomatic SE - This runs trimmomatic in single-end mode. If you have paired-end data you can run filter it with trimmomatic by specifying trimmomatic PE
# NAg_Rep1_sample.fastq.gz - The input file
# NAg_Rep1_sample_trimmed.fastq.gz - The output file
# ILLUMINACLIP:<path to adaptors>:2:30:10 - Here we provide a fastA file with the sequencing adaptor sequences. The numbers control how many base mis-mathces between your data and the fasta file are allowed
# LEADING:3 TRAILING:3 - Remove low-qualty bases (quality score below 3) from the beginning and end of the read
# SLIDINGWINDOW:4:20 - Trim reads from the end of the read until the average base quality in a 4bp window is above 20
# MINLEN:50 - Remove reads that are less than 50 bases long after trimming

#################################################################################################################################
# Step 2). Read alignment and PCR duplicate removal
#################################################################################################################################

# The next step is to align the reads to the genome. 
# Here we will only align data to mouse chromosome 1 (mouse genome version mm10)

# To do the alignment we will be using bowtie2
# bowtie2 requires an indexed genome to align the data to
# This has already been created for you so no need to do this here. But for reference a genome index can be created from the genome fastA file like so:
# bowtie2-build mm10_chr1.fa mm10
#
# Genome fasta files can be obtained from ensembl or the UCSC genome browser

# Next, align the data with bowtie2 - the results will be stores in the sequence alignment map (SAM) format

bowtie2 --very-sensitive-local -x /rds/projects/k/knowletj-mibtp-masterclasses/MITBP_Masterclass_2021/data_for_workshop/Genome/mm10_chr1 -U NAg_Rep1_sample_trimmed.fastq.gz -S NAg_Rep1_sample.sam
bowtie2 --very-sensitive-local -x /rds/projects/k/knowletj-mibtp-masterclasses/MITBP_Masterclass_2021/data_for_workshop/Genome/mm10_chr1 -U TAg_Rep1_sample_trimmed.fastq.gz -S TAg_Rep1_sample.sam

# Before we move to the next step, the aligned reads need to be sorted by chromosome and position
# We will also convert the files into the binary alignment map (BAM) format
# We can do this using samtools

samtools sort -o NAg_Rep1_sample.bam NAg_Rep1_sample.sam
samtools sort -o TAg_Rep1_sample.bam TAg_Rep1_sample.sam

# During the library prep step the fragments to be sequenced were amplified using PCR
# As a result, we might find that we have multiple identical reads that actually represent the same DNA fragment
# We may want to remove these to avoid any potential artefacts from over-amplified fragments
# Picard Tools is useful for this

java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=NAg_Rep1_sample.bam O=NAg_Rep1_sample_rmdup.bam M=NAg_Rep1_sample_picard_metrics.txt AS=true REMOVE_DUPLICATES=true
java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=TAg_Rep1_sample.bam O=TAg_Rep1_sample_rmdup.bam M=TAg_Rep1_sample_picard_metrics.txt AS=true REMOVE_DUPLICATES=true

# Here is what each of these options do:
# I=NAg_Rep1_sample.bam - The input file which is our sorted bam file
# O=NAg_Rep1_sample_rmdup.bam - The output file
# M=NAg_Rep1_sample_picard_metrics.txt - Ouput some metrics about the duplicate removal
# AS=true - Assume the reads are sorted 
# REMOVE_DUPLICATES=true - Remove PCR duplicated reads

# We also need to create an index for the bam file

samtools index NAg_Rep1_sample_rmdup.bam
samtools index TAg_Rep1_sample_rmdup.bam

#################################################################################################################################
# Step 3). Peak calling and filtering
#################################################################################################################################

# The next step will be to identify potential open chromatin regions - which we can also call peaks
# To do this we will use MACS2

macs2 callpeak -t NAg_Rep1_sample_rmdup.bam -n NAg_Rep1_sample --keep-dup all --nomodel -g mm -B --trackline
macs2 callpeak -t TAg_Rep1_sample_rmdup.bam -n TAg_Rep1_sample --keep-dup all --nomodel -g mm -B --trackline

# Here is what each of these options do:
# -t NAg_Rep1_sample_rmdup.bam - The aligned reads from the treatment file - in this case this is sorted, de-duplicated alignment that we have made
# -n NAg_Rep1_sample - What to call the output files
# --keep-dup all - Tells macs not to perform PCR duplicate removal - we have already done this
# --nomodel - Turns off MAC's peak shiting model. This model is useful for ChIP-Seq but for other types of data like  DNaseI of ATAC-Seq we do not need this
# -g mm - This specifies the effective genome size. MACS2 has some genome sizes for common model organisms built in built-in. Here we specific mouse (mus musculus) using the abbreviation mm
# -B --trackline - These options are used to produce a bedGraph file. This file can be uploaded to the UCSC genome browser

# This produces a number of output files - They are:
# NAg_Rep1_sample_peaks.narrowPeak - A BED file that contains the peak positions called by MACS2 along with p-value and q-value of those peaks
# NAg_Rep1_sample_summits.bed - A BED file with the summit position for every peak called
# NAg_Rep1_sample_peaks.xls - A tab-delimited file which contains the peak and summit positions as well as some other useful statistics (peak height, pvalue, qvalue)
# NAg_Rep1_sample_control_lambda.bdg - BedGraph for the input sample (used only with chip) 
# NAg_Rep1_sample_treat_pileup.bdg - BedGraph for the treatment sample (can be uploaded to UCSC genome browser)

# Count the number of peaks called by MACS2

wc -l NAg_Rep1_sample_summits.bed
wc -l TAg_Rep1_sample_summits.bed

# Next we will filter the peaks. We will do this in two steps:
# 1). Remove small background peaks
# 2). Remove blacklisted regions

# The below set of commands use the linux pipe to reads the peaks.xls file produced by macs and extract peaks with a height >= 10

cat NAg_Rep1_sample_peaks.xls | grep -v '#' | grep -v 'start' | awk '$6 >= 10' | awk '{print $1"\t"$5"\t"$5}' > NAg_Rep1_sample_filtered_summits.bed
cat TAg_Rep1_sample_peaks.xls | grep -v '#' | grep -v 'start' | awk '$6 >= 10' | awk '{print $1"\t"$5"\t"$5}' > TAg_Rep1_sample_filtered_summits.bed

# This is what each of these commands is doing
# cat NAg_Rep1_sample_peaks.xls - This reads the file
# grep -v '#' | grep -v 'start' - This removes the headers from the file (i.e. lines that begin with # or contain the word start)
# awk '$6 >= 10' - Keep only peaks with a value >= 10 in column 6
# awk '{print $1"\t"$5"\t$5}' - Keep only the 1st column (which contains the chromosome) and the 5th column (summit position) twice - this ensures we have BED formatted file
# > NAg_Rep1_sample_filtered_summits.bed - The output file - A BED file

# Count the number of peaks kept after filtering

wc -l NAg_Rep1_sample_filtered_summits.bed
wc -l TAg_Rep1_sample_filtered_summits.bed

# Next we will remove blacklisted regions
# These are regions which are known to give unusually high read coverage and are considered an artefact
# The mm10 blacklist can be found here as a bed file:
# It has already been downladed for you so you can do the filtering as so

# Set the path to the blacklist bed file
BLACKLIST=/rds/projects/k/knowletj-mibtp-masterclasses/MITBP_Masterclass_2021/data_for_workshop/Annotation/mm10-blacklist.v2.bed

bedtools intersect -v -a NAg_Rep1_sample_filtered_summits.bed -b $BLACKLIST > NAg_Rep1_sample_filtered_noBlackList_summits.bed
bedtools intersect -v -a TAg_Rep1_sample_filtered_summits.bed -b $BLACKLIST > TAg_Rep1_sample_filtered_noBlackList_summits.bed

# Count the number of peaks after removing blacklisted sites

wc -l NAg_Rep1_sample_filtered_noBlackList_summits.bed
wc -l TAg_Rep1_sample_filtered_noBlackList_summits.bed

#################################################################################################################################
# Step 4). Create a peak union and count reads
#################################################################################################################################

# Next, we will merge the peak from each of our samples to create a single bed file of all peaks
# To do this we will merge any peak which has a summit within 200bp of each other

# First, we need to extend our filtered summits by -/+ 100bp - this will create a 200bp peak
# We can use bedtools slop for this

CHROMSIZE=/rds/projects/k/knowletj-mibtp-masterclasses/MITBP_Masterclass_2021/data_for_workshop/Annotation/mm10.chrom.sizes

bedtools slop -b 100 -i NAg_Rep1_sample_filtered_noBlackList_summits.bed -g $CHROMSIZE > NAg_Rep1_sample_filtered_noBlackList_peaks.bed
bedtools slop -b 100 -i TAg_Rep1_sample_filtered_noBlackList_summits.bed -g $CHROMSIZE > TAg_Rep1_sample_filtered_noBlackList_peaks.bed

# Here is what these options do:
# -b 100 - extend summit by 100bp in both directions
# -i NAg_Rep1_sample_filtered_noBlackList_summits.bed - The input file
# -g mm10.chromSizes - A file that contains the size (in base pairs) of each of the chromosomes in mm10
# > NAg_Rep1_sample_filtered_noBlackList_peaks.bed - The output file

# Next, we can combine the two bed files in to a single file
# At this point we also need to sort the peaks and then merge peaks which are overlapping
# All of this can be done with bedtools

cat NAg_Rep1_sample_filtered_noBlackList_peaks.bed TAg_Rep1_sample_filtered_noBlackList_peaks.bed > Merged.bed
bedtools sort -i Merged.bed > Merged_sorted.bed
bedtools merge -i Merged_sorted.bed > Peak_Union.bed

# Count the number of peaks in the peak union

wc -l Peak_Union.bed

# Finally, we can count the number of reads per peak using featureCounts. 
# This is what we will use for the differential peak analysis

# featureCounts will not read BED files directly, so first we need to convert the peak union bed file into something featureCounts can read
# We will use the SAF file format here, which is extremely similar to a bed file with two extra columns

# For this we can use awk
cat Peak_Union.bed | awk '{print $1":"$2":"$3"\t"$1"\t"$2"\t"$3"\t+"}' > Peak_Union.saf

# This will create a file that looks something like:
# chr7:110166237:110166637        chr7    110166237       110166637       .
# chr17:79298879:79299279 chr17   79298879        79299279        .
# chr19:16789975:16790375 chr19   16789975        16790375        .

# The first column is just the peak coordinate separated by :
# The next 3 columns are the peak position in a BED like format
# The final column is the strand - since this is genomic data we do not actually have stranded data so a . will work fine here

# Now we can run featureCounts

featureCounts -F SAF -a Peak_Union.saf -o NAgTAg_DHS_rawCounts.tsv NAg_Rep1_sample_rmdup.bam TAg_Rep1_sample_rmdup.bam

# Here is what these options do:
# -F SAF - tells featureCounts our peaks are in the SAF format
# -a Peak_Union.saf - The input file
# -o NAgTAg-DHS_rawCounts.tsv - The output file
# NAg_Rep1_sample_rmdup.bam TAg_Rep1_sample_rmdup.bam - Here we provide the bam files we made earlier - use these to count the reads

# The resulting file can be used for downstream analysis

#################################################################################################################################
# End of part 1 - Move to workshop part 2
#################################################################################################################################
