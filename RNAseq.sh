#!/bin/bash
#SRR3401655 -> Strand specific RNA-seq on WT G27 strain nickel treated, replica B
#SRR3401654 -> Strand specific RNA-seq on WT G27 strain nickel treated, replica A
#SRR3401651 -> Strand specific RNA-seq on WT G27 strain untreated, replica B
#SRR3401650 -> Strand specific RNA-seq on WT G27 strain untreated, replica A

# SRR3401647	 -> Strand specific RNA-seq on deltaNikR strain nickel treated, replica B
#12. SRR3401643 -> Strand specific RNA-seq on deltaNikR strain nickel treated, replica A
#13. SRR3401641 -> Strand specific RNA-seq on deltaNikR strain untreated, replica B
#14. SRR3401620 -> Strand specific RNA-seq on deltaNikR strain untreated, replica A

#download
mkdir fastqRNA
mkdir tagDirsRNA
for SRA in SRR3401655 SRR3401654 SRR3401651 SRR3401650 SRR3401647 SRR3401643 SRR3401641 SRR3401620
do
    #download
    sratoolkit.3.0.0-ubuntu64/bin/fastq-dump -O fastqRNA ${SRA}
    #align
    bwa mem ncbi-genomes-2022-05-29/GCF_000021165.1_ASM2116v1_genomic.fna fastqRNA/${SRA}.fastq | samtools view -bS > ${SRA}.bam
    #make tag directory
    makeTagDirectory tagDirsRNA/${SRA} ${SRA}.bam
    #produce a bed graph using the directory
    makeUCSCfile tagDirsRNA/${SRA} -o auto
done
