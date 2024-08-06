# pipeline to analyze transcriptome

## Download reference genome fasta files and annotation gff files

## convert gff to gtf 
gffread ${gff_file} -T -o ${gtf_file}

## build STAR_idx
star_idx=Zm_star_idx
fa_file=/home/gujianhui/reference/Zea_mays/Zea_mays.fa
gtf_file=/home/gujianhui/reference/Zea_mays/Zea_mays.gtf

STAR --runMode genomeGenerate \
    --runThreadN 20 \
    --genomeDir ${star_idx} \
    --genomeFastaFiles ${fa_file}  \
    --sjdbGTFfile ${gtf_file} \
    --sjdbOverhang 149

## generate STAR manifest file  and manual adjust format
for i in $(ls -d /home/gujianhui/GRAS/01.transcriptome/01.data/MG*);
 do for j in $( ls ${i}| grep fq.gz) ;
    do echo -n -e ${i}/${j}"\t" ;
    done ;
   echo ${i##*/}  ;
done > maize_manifest

## mapping fq file by STAR
manifest=maize_manifest
out_prefix=MG_out
index_dir=/home/gujianhui/reference/Zea_mays/Zm_star_idx
STAR --twopassMode Basic \
    --quantMode TranscriptomeSAM GeneCounts \
    --runThreadN 20 \
    --genomeDir  ${index_dir} \
    --outSAMtype BAM SortedByCoordinate \
    --sjdbOverhang 149 \
    --outSAMattrRGline ID:sample SM:sample PL:ILLUMINA \
    --outFilterMismatchNmax 2 \
    --outSJfilterReads Unique \
    --outSAMmultNmax 1 \
    --outFileNamePrefix ${out_prefix} \
    --outSAMmapqUnique 60 \
    --readFilesCommand gunzip -c \
    --readFilesIn ${r1_files[0]},${r1_files[1]}   ${r2_files[0]},${r2_files[1]}

    --readFilesManifest ${manifest}
for i in MG-1 MG-2;
 do r1_files=($(grep ${i} MG_manifest | cut -f 1));
  r2_files=($(grep ${i} MG_manifest | cut -f 2)) ;
   out_prefix=${i}_out;
    STAR --twopassMode Basic   \
        --quantMode TranscriptomeSAM GeneCounts     \
        --runThreadN 20     \
        --genomeDir  ${index_dir}     \
        --outSAMtype BAM SortedByCoordinate     \
        --sjdbOverhang 149     \
        --outSAMattrRGline ID:sample SM:sample PL:ILLUMINA     \
        --outFilterMismatchNmax 2     \
        --outSJfilterReads Unique     \
        --outSAMmultNmax 1     \
        --outFileNamePrefix ${out_prefix}     \
        --outSAMmapqUnique 60     \
        --readFilesCommand gunzip -c     \
        --readFilesIn ${r1_files[0]},${r1_files[1]}   ${r2_files[0]},${r2_files[1]} ;
done

## status of mapping ratio of each sample
grep 'Uniquely mapped reads %' M*_outLog.final.out


## merge outReadsPerGene.out.tab files
for i in $(ls M*_outReadsPerGene.out.tab); do echo -n '"'${i}'"', ; done
python /home/gujianhui/pyscripts/merge_star_out.py  > maize_gene_count.tsv

## create self-made eggnog database
gffread Zea_mays.gff -g Zea_mays.fa -y Zea_mays.pep
emapper.py  -m diamond -i Zea_mays.pep  \
    --cpu 20 --data_dir /home/gujianhui/reference/eggnog_db/ \
    --dmnd_db /home/gujianhui/reference/eggnog_db/eggnog_proteins.dmnd  \
    --output_dir ./eggnog_idx  \
    -o Zm_eggnog

## remove the abundant annotation lines
sed -i '/^# /d' *.emapper.annotations 
sed -i 's/#//' *.emapper.annotations