#!/bin/bash
## hifiasm pipeline

sp=ldtc
python3.10 /home/gujianhui/bio_tools/juicer/misc/generate_site_positions.py   \
    DpnII ${sp} \
    ${sp}.hifi.fasta

awk 'BEGIN{OFS="\t"}{print $1, $NF}' ${sp}_DpnII.txt  > ${sp}.lens
bwa index ${sp}.hifi.fasta

# create a directory to instore fastq file

bash  /home/gujianhui/bio_tools/juicer/scripts/juicer.sh \
    -d /home/gujianhui/Poaceae_assembly/05.hic_test/${sp}_juicer \
    -g ${sp} \
    -D /home/gujianhui/bio_tools/juicer/ -z ${sp}.hifi.fasta \
    -y ${sp}_DpnII.txt  -p ${sp}.lens  \
    -s DpnII -t 20 \
    --assembly

# status of HIC data mapping results
for i in $(ls | grep juicer);
do echo ${i};
grep ':' ${i}/aligned/inter_30.txt | cut  -d ':' -f 2  ;
done 