#!/usr/bin/env python

import pysam
import array
import argparse
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

parser = argparse.ArgumentParser()
parser.add_argument("input_bam", type=str)
parser.add_argument("snp_file", type=str)
parser.add_argument("delete_bed", type=str)
parser.add_argument("n_proc", type=int)
arguments = parser.parse_args()
input_bam = arguments.input_bam
snp_file = arguments.snp_file
delete_bed = arguments.delete_bed
n_proc = arguments.n_proc

# Читаем файл с мутациями
snips = pd.read_csv(snp_file)

# Создаем Bed-файл с интервалами
with pysam.AlignmentFile(input_bam, "rb") as samfile_input, open(delete_bed, "a") as bed:
    chroms = samfile_input.references
    for chrom in chroms:

        # Для каждой хромосомы, упомянутой и в файле с мутациями, и в баме, отбираем интервал между первым и последним снипом
        if chrom in set(snips["chromosome"]):
            snp_subset = snips.query("chromosome==@chrom")
            chrom_snps = snp_subset["position"]
            first_p, last_p = chrom_snps.min(), chrom_snps.max()
            bed.write(f"{chrom}\t{first_p}\t{last_p + 1}\n")
