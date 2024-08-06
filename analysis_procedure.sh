# analysis_procedure
linker=TGGT,GAGA
barcodes=/home/gujianhui/mingxia_data/barcode_pools.txt

# cut R1 with P7 but trimmed reads < 100 bp were filtered
cutadapt -a AAGTCGGAGGCCAAGCGGTC \
        -j 20 -e 0.1 -O 5 -m 100 \
        --pair-filter=any \
        -o ${sampleID}_rmP_R1.fq.gz -p ${sampleID}_rmP_R2.fq.gz \
        ${sampleID}_R1.fq.gz  ${sampleID}_R2.fq.gz > ${sampleID}_rmP_log

# get reads with adapter
adapter=GTGAGTGATGGTTGAGGTAGTGTGGAG
cutadapt -g ${adapter} \
        -j 20 -e 0.1 -O 5 -m 50 \
        --discard-untrimmed \
        -o ${sampleID}_cutadapt_R1.fq.gz -p ${sampleID}_cutadapt_R2.fq.gz \
        ${sampleID}_rmP_R1.fq.gz  ${sampleID}_rmP_R2.fq.gz > ${sampleID}_cutadapt_log


# identify barcodes
python /home/gujianhui/mingxia_data/analysis_script/decompsing_barcode.py \
        -r1 ${sampleID}_cutadapt_R1.fq.gz  -r2 ${sampleID}_cutadapt_R2.fq.gz \
        -b  ${barcodes}\
        -l  ${linker}\
        -s 1 -t 20 \
        -o ${sampleID}_debar

# the second round of barcodes
adapter=CGGTGATGTGTGATAGGAGGAG
cutadapt -g ${adapter} \
        -j 20 -e 0.1 -O 5 -m 50 \
        --discard-untrimmed \
        -o ${sampleID}_sec_cutadapt_R1.fq.gz -p ${sampleID}_sec_cutadapt_R2.fq.gz \
        ${sampleID}_debar_trimmed_R1.fq.gz  ${sampleID}_debar_trimmed_R2.fq.gz > ${sampleID}_sec_cutadapt_log



# status of barcodes
total_read=$(wc -l ${sampleID}_debar_bar | cut -d ' ' -f 1)
invalid_read=$(grep -E 'NO|False' ${sampleID}_debar_bar | wc -l)
valid_read=$((${total_read}-${invalid_read}))
echo total_reads:${total_read}     invalid_read:${invalid_read}    valid_read:${valid_read}


# cut the second adapter
second_adpat=CGGTGATGTGTGATAGGAGGAG
cutadapt -g ${second_adpat} \
        -j 20 -e 0.1 -O 5 -m 50 \
        --discard-untrimmed \
        -o ${sampleID}_sec_cutadapt_R1.fq.gz -p ${sampleID}_sec_cutadapt_R2.fq.gz \
        ${sampleID}_debar_trimmed_R1.fq.gz  ${sampleID}_debar_trimmed_R2.fq.gz > ${sampleID}_sec_cutadapt_log




# quality control
fastp --in1 ${sampleID}_debar_trimmed_R1.fq.gz --in2 ${sampleID}_debar_trimmed_R2.fq.gz \
    --out1 ${sampleID}_clean_R1.fq.gz --out2 ${sampleID}_clean_R2.fq.gz \
    --json ${sampleID}.json --html ${sampleID}.html \
    --trim_poly_g --poly_g_min_len 10 \
    --trim_poly_x --poly_x_min_len 10 \
    --cut_front --cut_tail \
    --cut_window_size 4 \
    --qualified_quality_phred 30 \
    --low_complexity_filter --complexity_threshold 30 \
    --length_required 30 \
    --thread 10 


# bulid bwa-mem2 mapping index
# /home/gujianhui/bio_tools/bwa/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index \
#     -p mm10_bwa_mem_idx/mm10 Mus_musculus.GRCm38.dna_sm.toplevel.fa.gz

# map reads to genome
genome_idx=/home/gujianhui/reference/Mus_musculus/mm10_bwa_mem_idx/mm10
/home/gujianhui/bio_tools/bwa/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -t 20 \
    ${genome_idx}  \
    ${sampleID}_clean_R1.fq.gz ${sampleID}_clean_R2.fq.gz \
    -o ${sampleID}.sam &

samtools flagstat ${sampleID}.sam


samtools view -Sb ${sampleID}.sam > ${sampleID}.bam
samtools sort ${sampleID}.bam -o sorted_${sampleID}.bam
bedtools bamtobed -i sorted_${sampleID}.bam > ${sampleID}.bed



## build methylation Bismark index
# generate merged fasta
zact ../humo_spiens/human.fasta.gz | sed 's/^>/>hm/g' > hg_mm_genome.fasta
zcat ../Mus_musculus/mm10.fa.gz | sed 's/^>/>mm/g' >> hg_mm_genome.fasta

# build index
bismark_genome_preparation \
    --bowtie2 \
    --parallel 20 \
    --large-index  \
    ./ > bismark_genome_preparation.log 