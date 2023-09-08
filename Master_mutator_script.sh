#!/bin/bash

python3 Deleter_pysam.py "./genome38_192021.bam" "./snps.csv" "./delete.bed" 4

samtools view --threads 8 -L delete.bed -U genome38_192021_untarget.bam -o genome38_192021_target.bam genome38_192021.bam

samtools index genome38_192021_target.bam

python3 Mutator_pysam.py "./genome38_192021.bam" "./new_bam_order_mod.bam" "./snps.csv" "./genome38_192021_target.bam" "./genome38_192021_untarget.bam" "./delete.bed" 4 0.5

samtools index new_bam_order_mod.bam

samtools cat --threads 8 genome38_192021_untarget.bam new_bam_order_mod.bam -o genome38_192021_mod.bam

samtools sort --threads 8 genome38_192021_mod.bam -o genome38_192021_mod_sorted.bam

samtools index -@ 8 genome38_192021_mod_sorted.bam
