#!/bin/bash

python3 Pysam_cnv_bed_creator.py "nano_genome.bam" "cnv_delete.bed" "cnv.csv"

samtools view nano_genome.bam -@ 8 -L cnv_delete.bed -U nano_genome_untarget.bam -o nano_genome_target.bam
samtools view nano_genome.bam -@ 8 -L cnv_delete.bed -M -o nano_genome_target.bam
samtools index nano_genome_target.bam

python3 Pysam_cnv_maker.py "nano_genome_target.bam" "nano_genome_mod.bam" "cnv.csv"

samtools cat -@ 8 nano_genome_untarget.bam nano_genome_mod.bam -o nano_genome_new_mod.bam
samtools sort -@ 8 nano_genome_new_mod.bam -o nano_genome_new_mod_sorted.bam
samtools index -@ 8 nano_genome_new_mod_sorted.bam
