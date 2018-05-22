#!/bin/bash

## Hard Code Params
DIR='/scratch/cgsb/gencore/out/Gresham/2015-10-23_HK5NHBGXX/lib1-26/'
REF='/scratch/work/cgsb/reference_genomes/Public/Fungi/Saccharomyces_cerevisiae/R64-1-1/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa'
SNPEFF_DB='Saccharomyces_cerevisiae'
PL='illumina'
PM='nextseq'
EMAIL=${USER}@nyu.edu

## Or Set Params via Command Line
#DIR=$1 #DIRECTORY WITH FASTQ FILES TO PROCESS
#REF=$2 #UNCOMPRESSED VERSION OF GENOMIC.FNA FILE
#SNPEFFDB=$3    #SNPEFF DB TO USE
#PL=$4  #EX: illumina
#PM=$5  #EX: nextseq
#EMAIL=$6 #EX: netID@nyu.edu


## Modules Used in this script
BEDTOOLS='bedtools/intel/2.26.0'
BWA='bwa/intel/0.7.15'
PICARD='picard/2.8.2'
GATK='gatk/3.7-0'
R='r/intel/3.4.2'
SAMTOOLS='samtools/intel/1.3.1'
SNPEFF='snpeff/4.3i'

## Paths to Modules
PICARD_JAR='/share/apps/picard/2.8.2/picard-2.8.2.jar'
GATK_JAR='/share/apps/gatk/3.7-0/GenomeAnalysisTK.jar'
SNPEFF_JAR='/share/apps/snpeff/4.3i/snpEff.jar'

## Get the current working directory
WD="$(pwd)"

## Steps of the Workflow
pre_process(){

com="cd $FWD && \
module load $BWA && \
bwa mem -M -R '@RG\tID:$file\tLB:$file\tPL:$PL\tPM:$PM\tSM:$file' \
$REF \
$INPUT_1 $INPUT_2 \
> ${ID}_aligned_reads.sam"
response=\
$(sbatch -J $ID.bwa -o $ID.bwa.out -e $ID.bwa.err --mail-user=$EMAIL --mail-type=FAIL --nodes=1 -t 4:00:00 --mem=8000 --wrap="$com")
stringarray=($response)
alignment=${stringarray[-1]}
echo $alignment > $ID.log
echo "ALIGNMENT: " $alignment
echo "Alignment Submitted"

com="cd $FWD && \ 
module load $PICARD && \
java -jar $PICARD_JAR \
SortSam \
INPUT=${ID}_aligned_reads.sam \
OUTPUT=${ID}_sorted_reads.bam \
SORT_ORDER=coordinate"
response=\
$(sbatch -J $ID.samToSortedBam -o $ID.samToSortedBam.out -e $ID.samToSortedBam.err --dependency=afterok:$alignment --kill-on-invalid-dep=yes --mail-user=$EMAIL --mail-type=FAIL --nodes=1 -t 4:00:00 --mem=8000 --wrap="$com")
stringarray=($response)
samToSortedBam=${stringarray[-1]}
echo $samToSortedBam >> $ID.log
echo "SAMTOSORTEDBAM: " $samToSortedBam
echo "samToSortedBam Submitted"

com="cd $FWD && \
module  load $PICARD && \
module load $R && \
module load $SAMTOOLS && \
java -jar $PICARD_JAR \
CollectAlignmentSummaryMetrics \
R=$REF \
I=${ID}_sorted_reads.bam \
O=${ID}_alignment_metrics.txt && \
java -jar $PICARD_JAR \
CollectInsertSizeMetrics \
INPUT=${ID}_sorted_reads.bam \
OUTPUT=${ID}_insert_metrics.txt \
HISTOGRAM_FILE=${ID}_insert_size_histogram.pdf && \
samtools depth -a ${ID}_sorted_reads.bam > ${ID}_depth_out.txt"
response=\
$(sbatch -J $ID.getMetrics -o $ID.getMetrics.out -e $ID.getMetrics.err --dependency=afterok:$samToSortedBam --kill-on-invalid-dep=yes --mail-user=$EMAIL --mail-type=FAIL --nodes=1 -t 4:00:00 --mem=8000 --wrap="$com")
stringarray=($response)
getMetrics=${stringarray[-1]}
echo $getMetrics >> $ID.log
echo "GETMETRICS: " $getMetrics
echo "getMetrics Submitted"


com="cd $FWD && \ 
rm ${ID}_aligned_reads.sam && \
module load $PICARD && \
java -jar $PICARD_JAR \
MarkDuplicates \
INPUT=${ID}_sorted_reads.bam \
OUTPUT=${ID}_dedup_reads.bam \
METRICS_FILE=${ID}_metrics.txt" 
response=\
$(sbatch -J $ID.markDuplicates -o $ID.markDuplicates.out -e $ID.markDuplicates.err --dependency=afterok:$samToSortedBam --kill-on-invalid-dep=yes --mail-user=$EMAIL --mail-type=FAIL --nodes=1 -t 4:00:00 --mem=60000 --wrap="$com")
stringarray=($response)
markDuplicates=${stringarray[-1]}
echo $markDuplicates >> $ID.log
echo "MARKDUPLICATES: " $markDuplicates
echo "Mark Duplicates Submitted"

com="cd $FWD && \ 
module load $PICARD && \
java -jar $PICARD_JAR \
BuildBamIndex \
INPUT=${ID}_dedup_reads.bam" 
response=\
$(sbatch -J $ID.buildBamIndex -o $ID.buildBamIndex.out -e $ID.buildBamIndex.err --dependency=afterok:$markDuplicates --kill-on-invalid-dep=yes --mail-user=$EMAIL --mail-type=FAIL --nodes=1 -t 4:00:00 --mem=60000 --wrap="$com")
stringarray=($response)
buildBamIndex=${stringarray[-1]}
echo $buildBamIndex >> $ID.log
echo "BUILDBAMINDEX: " $buildBamIndex
echo "BuildBamIndex Submitted"

com="cd $FWD && \ 
rm ${ID}_sorted_reads.bam && \
module load $GATK && \
java -jar $GATK_JAR \
-T RealignerTargetCreator \
-R $REF \
-I ${ID}_dedup_reads.bam \
-o ${ID}_realignment_targets.list"
response=\
$(sbatch -J $ID.realignTargetCreator -o $ID.realignTargetCreator.out -e $ID.realignTargetCreator.err --dependency=afterok:$buildBamIndex:$getMetrics --kill-on-invalid-dep=yes --mail-user=$EMAIL --mail-type=FAIL --nodes=1 -t 4:00:00 --mem=60000 --wrap="$com")
stringarray=($response)
realignTargetCreator=${stringarray[-1]}
echo $realignTargetCreator >> $ID.log
echo "REALIGNTARGETCREATOR: " $realignTargetCreator
echo "RealignTargetCreator Submitted"

com="cd $FWD && \
module load $GATK && \
java -jar $GATK_JAR \
-T IndelRealigner \
-R $REF \
-I ${ID}_dedup_reads.bam \
-targetIntervals ${ID}_realignment_targets.list \
-o ${ID}_realigned_reads.bam"
response=\
$(sbatch -J $ID.realignIndels -o $ID.realignIndels.out -e $ID.realignIndels.err --dependency=afterok:$realignTargetCreator --kill-on-invalid-dep=yes --mail-user=$EMAIL --mail-type=FAIL --nodes=1 -t 4:00:00 --mem=60000 --wrap="$com")
stringarray=($response)
realignIndels=${stringarray[-1]}
echo $realignIndels >> $ID.log
echo "REALIGNINDELS: " $realignIndels
echo "RealignIndels Submitted"	
}

call_variants(){
	ROUND=$1
	if [[ $ROUND -eq 1 ]];then
		INPUT=${ID}_realigned_reads.bam
		OUTPUT=${ID}_raw_variants.vcf
		AFTEROK=$realignIndels
	fi
	if [[ $ROUND -eq 2 ]];then
		INPUT=${ID}_recal_reads.bam
		OUTPUT=${ID}_raw_variants_recal.vcf
		AFTEROK=$applyBqsr
	fi

com="cd $FWD && \
module load $GATK && \
java -jar $GATK_JAR \
-T HaplotypeCaller \
-R $REF \
-I $INPUT \
-o $OUTPUT" 
response=\
$(sbatch -J $ID.callVariants$ROUND -o $ID.callVariants$ROUND.out -e $ID.callVariants$ROUND.err --dependency=afterok:$AFTEROK --kill-on-invalid-dep=yes --mail-user=$EMAIL --mail-type=FAIL --nodes=1 -t 4:00:00 --mem=60000 --wrap="$com")
stringarray=($response)
callVariants=${stringarray[-1]}
echo $callVariants >> $ID.log
echo "CALLVARIANTS: " $callVariants
echo "Variant Calling Round $ROUND Submitted"
	
	if [[ $ROUND -eq 1 ]];then
		callVariants_1=$callVariants
	fi
	if [[ $ROUND -eq 2 ]];then
		callVariants_2=$callVariants
	fi
}

## extracts snps AND indels into seperate vcf files
extract_snps(){
        ROUND=$1
        if [[ $ROUND -eq 1 ]];then
		V=${ID}_raw_variants.vcf
		OS=${ID}_raw_snps.vcf
		OI=${ID}_raw_indels.vcf
		AFTEROK=$callVariants_1
        fi
        if [[ $ROUND -eq 2 ]];then
		V=${ID}_raw_variants_recal.vcf
		OS=${ID}_raw_snps_recal.vcf
		OI=${ID}_raw_indels_recal.vcf
		AFTEROK=$callVariants_2
        fi
		
com="cd $FWD &&\
module load $GATK && \
java -jar $GATK_JAR \
-T SelectVariants \
-R $REF \
-V $V \
-selectType SNP \
-o $OS && \
java -jar $GATK_JAR \
-T SelectVariants \
-R $REF \
-V $V \
-selectType INDEL \
-o $OI"
response=\
$(sbatch -J $ID.extractSnps$ROUND -o $ID.extractSnps$ROUND.out -e $ID.extractSnps$ROUND.err --dependency=afterok:$AFTEROK --kill-on-invalid-dep=yes --mail-user=$EMAIL --mail-type=FAIL --nodes=1 -t 4:00:00 --mem=60000 --wrap="$com")
stringarray=($response)
extractSnps=${stringarray[-1]}
echo $extractSnps >> $ID.log
echo "EXTRACTSNPS: " $extractSnps
echo "Extract SNPs Round $ROUND Submitted"

        if [[ $ROUND -eq 1 ]];then
                extractSnps_1=$extractSnps
        fi
        if [[ $ROUND -eq 2 ]];then
                extractSnps_2=$extractSnps
        fi
}

filter_snps(){
        ROUND=$1
        if [[ $ROUND -eq 1 ]];then
                V=${ID}_raw_snps.vcf
		O=${ID}_filtered_snps.vcf
                AFTEROK=$extractSnps_1
        fi
        if [[ $ROUND -eq 2 ]];then
                V=${ID}_raw_snps_recal.vcf
                O=${ID}_filtered_snps_final.vcf
                AFTEROK=$extractSnps_2
        fi

com="cd $FWD && \
module load $GATK && \
java -jar $GATK_JAR \
-T VariantFiltration \
-R $REF \
-V $V \
-filterName \"QD_filter\" \
-filter \"QD < 2.0\" \
-filterName \"FS_filter\" \
-filter \"FS > 60.0\" \
-filterName \"MQ_filter\" \
-filter \"MQ < 40.0\" \
-filterName \"SOR_filter\" \
-filter \"SOR > 4.0\" \
-o $O"
response=\
$(sbatch -J $ID.filterSnps$ROUND -o $ID.filterSnps$ROUND.out -e $ID.filterSnps$ROUND.err --dependency=afterok:$AFTEROK --kill-on-invalid-dep=yes --mail-user=$EMAIL --mail-type=FAIL --nodes=1 -t 4:00:00 --mem=60000 --wrap="$com")
stringarray=($response)
filterSnps=${stringarray[-1]}
echo $filterSnps >> $ID.log
echo "FILTEREDSNPS: " $filterSnps
echo "Filter SNPs Round $ROUND Submitted"

        if [[ $ROUND -eq 1 ]];then
                filterSnps_1=$filterSnps
        fi
        if [[ $ROUND -eq 2 ]];then
                filterSnps_2=$filterSnps
        fi
}

filter_indels(){
        ROUND=$1
        if [[ $ROUND -eq 1 ]];then
                V=${ID}_raw_indels.vcf
                O=${ID}_filtered_indels.vcf
                AFTEROK=$extractSnps_1
        fi
        if [[ $ROUND -eq 2 ]];then
                V=${ID}_raw_indels_recal.vcf
                O=${ID}_filtered_indels_final.vcf
                AFTEROK=$extractSnps_2
        fi

com="cd $FWD && \
module load $GATK && \
java -jar $GATK_JAR \
-T VariantFiltration \
-R $REF \
-V $V \
-filterName \"QD_filter\" \
-filter \"QD < 2.0\" \
-filterName \"FS_filter\" \
-filter \"FS > 200.0\" \
-filterName \"SOR_filter\" \
-filter \"SOR > 10.0\" \
-o $O"
response=\
$(sbatch -J $ID.filterIndels$ROUND -o $ID.filterIndels$ROUND.out -e $ID.filterIndels$ROUND.err --dependency=afterok:$AFTEROK --kill-on-invalid-dep=yes --mail-user=$EMAIL --mail-type=FAIL --nodes=1 -t 4:00:00 --mem=60000 --wrap="$com")
stringarray=($response)
filterIndels=${stringarray[-1]}
echo $filterIndels >> $ID.log
echo "FILTEREDINDELS: " $filterIndels
echo "Filter Indels Round $ROUND Submitted"

}

do_bqsr(){
	#todo: knownSites input shouldnt be full raw_variants.vcf file but only the TOP variants!
	ROUND=$1
	if [[ $ROUND -eq 1 ]];then
		POST=''
		OUT=${ID}_recal_data.table
		AFTEROK="$filterSnps_1:$filterIndels"
	fi
	if [[ $ROUND -eq 2 ]];then
		POST='-BQSR '${ID}_recal_data.table
		OUT=${ID}_post_recal_data.table
		AFTEROK=$bqsr_1
	fi

com="cd $FWD && \
module load $GATK && \
java -jar $GATK_JAR \
-T BaseRecalibrator \
-R $REF \
-I ${ID}_realigned_reads.bam \
-knownSites ${ID}_filtered_snps.vcf \
-knownSites ${ID}_filtered_indels.vcf \
$POST \
-o $OUT"
response=\
$(sbatch -J $ID.bqsr$ROUND -o $ID.bqsr$ROUND.out -e $ID.bqsr$ROUND.err --dependency=afterok:$AFTEROK --kill-on-invalid-dep=yes --mail-user=$EMAIL --mail-type=FAIL --nodes=1 -t 4:00:00 --mem=60000 --wrap="$com")
stringarray=($response)
bqsr=${stringarray[-1]}
echo $bqsr >> $ID.log
echo "BQSR: " $bqsr
echo "BQSR Round $ROUND Submitted"

	if [[ $ROUND -eq 1 ]];then
		bqsr_1=$bqsr
	fi
	if [[ $ROUND -eq 2 ]];then
		bqsr_2=$bqsr
	fi
}

analyze_covariates(){
com="cd $FWD && \
module load $GATK && \
module load $R && \
java -jar $GATK_JAR \
-T AnalyzeCovariates \
-R $REF \
-before ${ID}_recal_data.table \
-after ${ID}_post_recal_data.table \
-plots ${ID}_recalibration_plots.pdf"
response=\
$(sbatch -J $ID.analyzeCovariates -o $ID.analyzeCovariates.out -e $ID.analyzeCovariates.err --dependency=afterok:$bqsr_2 --kill-on-invalid-dep=yes --mail-user=$EMAIL --mail-type=FAIL --nodes=1 -t 4:00:00 --mem=60000 --wrap="$com")
stringarray=($response)
analyzeCovariates=${stringarray[-1]}
echo $analyzeCovariates >> $ID.log
echo "ANALYZECOVARIATES: " $analyzeCovariates
echo "AnalyzeCovariates Submitted"
}

apply_bqsr(){
com="cd $FWD && \
rm ${ID}_dedup_reads.bam ${ID}_dedup_reads.bai && \
module load $GATK && \
java -jar $GATK_JAR \
-T PrintReads \
-R $REF \
-I ${ID}_realigned_reads.bam \
-BQSR ${ID}_recal_data.table \
-o ${ID}_recal_reads.bam"
response=\
$(sbatch -J $ID.applyBqsr -o $ID.applyBqsr.out -e $ID.applyBqsr.err --dependency=afterok:$bqsr_2 --kill-on-invalid-dep=yes --mail-user=$EMAIL --mail-type=FAIL --nodes=1 -t 4:00:00 --mem=60000 --wrap="$com")
stringarray=($response)
applyBqsr=${stringarray[-1]}
echo $applyBqsr >> $ID.log
echo "APPLYBQSR: " $applyBqsr
echo "Apply BQSR Submitted"
}

parse_metrics(){
com="cd $FWD && \
module load $BEDTOOLS && \
bedtools genomecov -bga -ibam ${ID}_recal_reads.bam > ${ID}_genomecov.bedgraph && \
sh /scratch/work/cgsb/scripts/variant_calling/parse_metrics.sh $ID"
response=\
$(sbatch -J $ID.parseMetrics -o $ID.parseMetrics.out -e $ID.parseMetrics.err --dependency=afterok:$filterSnps_2 --kill-on-invalid-dep=yes --mail-user=$EMAIL --mail-type=FAIL --nodes=1 -t 4:00:00 --mem=60000 --wrap="$com")
stringarray=($response)
parseMetrics=${stringarray[-1]}
echo $parseMetrics >> $ID.log
echo "PARSEMETRICS: " $parseMetrics
echo "ParseMetrics Submitted"
}

do_snpeff(){
com="cd $FWD && \
module load $SNPEFF && \
java -jar $SNPEFF_JAR \
-v $SNPEFF_DB \
${ID}_filtered_snps_final.vcf > ${ID}_filtered_snps_final.ann.vcf"
response=\
$(sbatch -J $ID.snpEff -o $ID.snpEff.out -e $ID.snpEff.err --dependency=afterok:$filterSnps_2 --kill-on-invalid-dep=yes --mail-user=$EMAIL --mail-type=FAIL --nodes=1 -t 4:00:00 --mem=60000 --wrap="$com")
stringarray=($response)
snpEff=${stringarray[-1]}
echo $snpEff >> $ID.log
echo "SNPEFF: " $snpEff
echo "SnpEFF submitted"
}

## Build Report Header and Create File
REPORT_HEADER="ID,# reads,aligned reads,% aligned,aligned bases,read length,% paired,mean insert size,# SNPs 1,# SNPs 1 filtered,# SNPs 2, # SNPs filtered 2,average coverage"
echo $REPORT_HEADER > report.csv

for file in $(ls $DIR/*fastq.gz | perl -pe 's/^.+n0\d_(.+)\.fastq\.gz/$1/g' | sort | uniq)
do
	INPUT_1=$(ls $DIR/*n01_$file.fastq.gz)
	INPUT_2=$(ls $DIR/*n02_$file.fastq.gz)
	ID=$file

	echo "Processing files: " $INPUT_1 ", " $INPUT_2
	FWD=$WD/$ID
	mkdir $FWD
	cd $FWD

	## CALL WORKFLOW
	pre_process # Only Done Once
	call_variants 1 # Call Variants Round 1
	extract_snps 1 # Round 1. Extracts snps AND indels, separately
	filter_snps 1 # Round 1
	filter_indels 1 # Round 1
	do_bqsr 1 # Do BQSR Round 1
	do_bqsr 2 # Do BQSR Round 2
	analyze_covariates # Only Done Once
	apply_bqsr # Only Done Once
	call_variants 2 # Call Variants Round 2
	extract_snps 2 # Round 2. Extracts snps AND indels, separately
	filter_snps 2 # Round 2
	filter_indels 2 # Round 2
	parse_metrics
	do_snpeff
done
