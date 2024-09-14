# analysis pipeline for ChIP-seq 

## 01. PREPARATION
# Firstly, download the fastq file and reference genome 

# Then, we need to build a genome index and clean the fastq file for mapping procedure.
# software like Bowtie/Hisat/BWA are alternative.

# we take BWA-MEM2 as an example
# Usage: bwa-mem2 index [-p prefix] <in.fasta>
# the genome index was build in 
# current path '/home/gujianhui/reference/Oryza_sativa/BWA_MEM_idx/' with prefix 'Os_BWA_idx'
/home/gujianhui/bio_tools/bwa/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index -p Os_BWA_idx   \
    /home/gujianhui/reference/Oryza_sativa/Oryza_sativa.fa

# We use fastp to make quality control for Rawdata
# Sometimes, we need to correct the pattern of input file name( -in1, -in2 )
sampleID=9-FLAG-CR
fastp --in1 ${sampleID}_1.fq.gz --in2 ${sampleID}_2.fq.gz \
    --out1 ${sampleID}_clean_R1.fq.gz --out2 ${sampleID}_clean_R2.fq.gz \
    --json ${sampleID}.json --html ${sampleID}.html \
    --trim_poly_g --poly_g_min_len 10 \
    --trim_poly_x --poly_x_min_len 10 \
    --cut_front --cut_tail \
    --cut_window_size 4 \
    --qualified_quality_phred 30 \
    --low_complexity_filter --complexity_threshold 30 \
    --length_required 30 \
    --thread 10 > ${sampleID}_fastp.log


## Mapping
# Use BWA-MEM2 to map all reads to indexed genome
genome_idx=/home/gujianhui/reference/Oryza_sativa/BWA_MEM_idx/Os_BWA_idx
data_path=/home/gujianhui/GRAS/02.Chip_seq/01.analysis/00.data
/home/gujianhui/bio_tools/bwa/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -t 20 \
    ${genome_idx}  \
    ${data_path}/${sampleID}_clean_R1.fq.gz ${data_path}/${sampleID}_clean_R2.fq.gz \
    -o ${sampleID}.sam &

# simple status of mapping result and sorted/filtered mapped bam result
samtools view -Sb -F 0x04 ${sampleID}.sam > ${sampleID}.bam
samtools sort --threads 20 ${sampleID}.bam -o sorted_${sampleID}.bam

samtools index sorted_${sampleID}.bam
bamCoverage -b sorted_${sampleID}.bam -o sorted_${sampleID}_raw.bw 

# ## Convert bam format into bed file format
# bedtools bamtobed -i sorted_${sampleID}.bam -bedpe > sorted_${sampleID}.bed
# ## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
# awk '$1==$4 && $6-$2 < 1000 {print $0}' sorted_${sampleID}.bed >sorted_${sampleID}.clean.bed
# ## Only extract the fragment related columns
# cut -f 1,2,6 sorted_${sampleID}.clean.bed | sort -k1,1 -k2,2n -k3,3n  >sorted_${sampleID}.fragments.bed

# ChIP-seq using exact sites of integration are affected by the accessibility of surrounding DNA. 
# For this reason fragments that share exact starting and ending positions are expected to be common, 
# and such ‘duplicates’ may not be due to duplication during PCR
# if deduplicated is required, run the following command:
# picard MarkDuplicates  REMOVE_DUPLICATES=true  \
#         MAX_FILE_HANDLES=800  VALIDATION_STRINGENCY=LENIENT  \
#         I=sorted_${sampleID}.bam  O=sorted_${sampleID}.rmdup.bam  \
#         M=${sampleID}.markdup.log.txt

# assess insert fragment size
samtools view -F 0x04 sorted_${sampleID}.bam | \
    awk -F '\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print $1"\t"abs($9)}'| \
    sort | uniq  | cut -f 2 | sort | uniq -c  | \
    awk -v OFS="\t" '{print $2, $1}' >${sampleID}_fragmentLen.txt  # extract insert size from bam file

# plot the fragment distribution
# conda activate Deeptools
# bamPEFragmentSize -hist ${sampleID}_fragmentSize.png -T "Fragment size of ${sampleID}" --maxFragmentLength 1000 \
#     -b sorted_${sampleID}.bam\
#     --samplesLabel ${sampleID}
python /home/gujianhui/GRAS/00.scripts/plot_frg.py ${sampleID} 1000



# skip the step of Spike-in calibration
spike_idx=/home/gujianhui/reference/E_coli/BWA_MEM_idx/Ec_spike_idx
data_path=/home/gujianhui/GRAS/02.Chip_seq/01.analysis/00.data

/home/gujianhui/bio_tools/bwa/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -t 20 \
    ${spike_idx}  \
    ${data_path}/${sampleID}_clean_R1.fq.gz ${data_path}/${sampleID}_clean_R2.fq.gz \
    -o ${sampleID}_spike.sam 

seqDepthDouble=`samtools view -F 0x04 ${sampleID}_spike.sam | wc -l`
seqDepth=$((seqDepthDouble/2))
echo $seqDepth >$projPath/alignment/sam/bowtie2_summary/${histName}_bowtie2_spikeIn.seqDepth


# Peak calling
# We should have control assays to remove the 
macs3 callpeak -t sorted_4-FLAG-CHIP.bam  -f BAMPE -n 4-ChIP -g 3.64e8 -B -q 0.01 --keep-dup all
macs3 callpeak -t sorted_4-FLAG-CR.bam  -f BAMPE -n 4-CR -g 3.64e8 -B -q 0.01 --keep-dup all
macs3 callpeak -t sorted_4-H3-CR.bam  -f BAMPE -n 4-H3 -g 3.64e8 -B -q 0.01 --keep-dup all
macs3 callpeak -t sorted_9-FLAG-CHIP.bam  -f BAMPE  -n 9-ChIP -g 3.64e8 -B -q 0.01 --keep-dup all
macs3 callpeak -t sorted_9-FLAG-CR.bam  -f BAMPE -n 9-CR -g 3.64e8 -B -q 0.01 --keep-dup all



# calculate the interaction matrix nearby TSS
computeMatrix scale-regions -S sorted_4-FLAG-CHIP_raw.bw \
                               sorted_4-FLAG-CR_raw.bw \
                               sorted_4-H3-CR_raw.bw \
                               sorted_9-FLAG-CHIP_raw.bw \
                               sorted_9-FLAG-CR_raw.bw \
                              -R /home/gujianhui/reference/Oryza_sativa/Oryza_sativa.gene.bed \
                              --beforeRegionStartLength 3000 \
                              --regionBodyLength 3000 \
                              --afterRegionStartLength 3000 \
                              --skipZeros -o raw_matrix_gene.mat.gz -p 20

plotHeatmap -m raw_matrix_gene.mat.gz -out raw_heatmap.png --sortUsing sum


# plot interaction matrix nearby PEAKs identified by MACS3
computeMatrix reference-point -S sorted_9-FLAG-CHIP_raw.bw \
                              -R 9-ChIP_peaks.narrowPeak \
                              -a 3000 \
                              --referencePoint center \
                              -b  3000 \
                              --skipZeros -o 9-ChIP_Peak.mat.gz -p 20

plotHeatmap -m 9-ChIP_Peak.mat.gz -out 9-CHIP_peak_heatmap.png --sortUsing sum --startLabel "Peak Start" -\
-endLabel "Peak End" --xAxisLabel "" --regionsLabel "Peaks" --samplesLabel "9-CHIP"

# plot peaks distribution compared to genes
# narrowPeak bed format convert to bw
cut -f 1,2 Oryza_sativa.fa.fai > Oryza_sativa.chrom.size
cut -f 1,2,3,5 4-CR_peaks.narrowPeak > 4-CR_narrowPeak.bedgraph

/home/gujianhui/bio_tools/bedGraphToBigWig 4-CR_narrowPeak.bedgraph \
    /home/gujianhui/reference/Oryza_sativa/Oryza_sativa.chrom.size \
    4-CR_narrowPeak.bw


# compute like raw.bw
computeMatrix scale-regions -S 4-ChIP_narrowPeak.bw \
                               4-CR_narrowPeak.bw \
                               4-H3_narrowPeak.bw\
                               9-ChIP_narrowPeak.bw \
                               9-CR_narrowPeak.bw \
                              -R /home/gujianhui/reference/Oryza_sativa/Oryza_sativa.gene.bed \
                              --beforeRegionStartLength 3000 \
                              --regionBodyLength 3000 \
                              --afterRegionStartLength 3000 \
                              --skipZeros -o all_peaks_matrix_gene.mat.gz -p 20

plotHeatmap -m all_peaks_matrix_gene.mat.gz -out all_peaks_heatmap.png --sortUsing sum


# extract the mapped reads on Peak region
# extraction
bedtools intersect -a  sorted_9-FLAG-CR.bam  -b 9-CR_peaks.narrowPeak  > 9-CR_peaks.bam
bedtools intersect -a  sorted_9-FLAG-CHIP.bam -b 9-ChIP_peaks.narrowPeak  > 9-ChIP_peaks.bam
bedtools intersect -a  sorted_4-H3-CR.bam  -b 4-H3_peaks.narrowPeak  > 4-H3_peaks.bam
bedtools intersect -a  sorted_4-FLAG-CR.bam  -b 4-CR_peaks.narrowPeak  > 4-CR_peaks.bam
bedtools intersect -a  sorted_4-FLAG-CHIP.bam  -b 4-ChIP_peaks.narrowPeak  > 4-ChIP_peaks.bam

# format conversion
samtools index sorted_${sampleID}.bam
bamCoverage -b 9-ChIP_peaks.bam -o 9-ChIP_peaks.bw

# Compute interaction matrix
computeMatrix scale-regions -S *.bw \
                              -R /home/gujianhui/reference/Oryza_sativa/Oryza_sativa.gene.bed \
                              --beforeRegionStartLength 3000 \
                              --regionBodyLength 3000 \
                              --afterRegionStartLength 3000 \
                              --skipZeros -o all_peaks_matrix_gene.mat.gz -p 20
                
plotHeatmap -m all_peaks_matrix_gene.mat.gz -out all_peaks_heatmap.png --sortUsing sum



# extract top0.01 peaks
a=$( wc -l 9-CR_peaks.narrowPeak | cut -d ' ' -f 1 ) ; 
sort -k 5nr 9-CR_peaks.narrowPeak | head -$(( ${a} / 100 )) | sort -k 1n -k 2n > peak_ananlysis/9-CR_top0.01_peaks.bed

# format conversion
for i in 4-ChIP 4-CR 4-H3 9-ChIP 9-CR;
    do bedtools intersect -a  ${i}_peaks.bam -b ${i}_top0.01_peaks.bed  > ${i}_top0.01.peaks.bam ; done

samtools index sorted_${sampleID}.bam
bamCoverage -b 4-ChIP_top0.01.peaks.bam -o 4-ChIP_top0.01.peaks.bw


computeMatrix scale-regions -S *top0.01.peaks.bw  \
                              -R /home/gujianhui/reference/Oryza_sativa/Oryza_sativa.gene.bed \
                              --beforeRegionStartLength 3000 \
                              --regionBodyLength 3000 \
                              --afterRegionStartLength 3000 \
                              --skipZeros -o top0.01_peaks_matrix_gene.mat.gz -p 20

plotHeatmap -m top0.01_peaks_matrix_gene.mat.gz -out top0.01_peaks_heatmap.png --sortUsing sum


# take gene as a point
computeMatrix reference-point -S *_peaks.bw \
                              -R /home/gujianhui/reference/Oryza_sativa/Oryza_sativa.gene.bed \
                              -b 3000 \
                              --referencePoint center \
                              -a 3000 \
                              --skipZeros -o all_peaks_matrix_gene_point.mat.gz -p 20

plotHeatmap -m all_peaks_matrix_gene_point.mat.gz -out all_peaks_matrix_gene_point.png --sortUsing sum