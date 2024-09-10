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

# simple status of mapping result and sorted bam result
samtools flagstat ${sampleID}.sam
samtools view -Sb ${sampleID}.sam > ${sampleID}.bam
samtools sort --threads 20 ${sampleID}.bam -o sorted_${sampleID}.bam