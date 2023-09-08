#!/usr/bin/env python

import pysam
import array
import argparse
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

parser = argparse.ArgumentParser()
parser.add_argument("input_bam", type=str)
parser.add_argument("output_bam", type=str)
parser.add_argument("snp_file", type=str)
parser.add_argument("target_reads", type=str)
parser.add_argument("untarget_reads", type=str)
parser.add_argument("delete_bed", type=str)
parser.add_argument("n_proc", type=int)
parser.add_argument("probability", type=float)
arguments = parser.parse_args()
input_bam = arguments.input_bam
output_bam = arguments.output_bam
snp_file = arguments.snp_file
target_reads = arguments.target_reads
untarget_reads = arguments.untarget_reads
delete_bed = arguments.delete_bed
n_proc = arguments.n_proc
probability = arguments.probability

def prob_mutation(read, iters, number_list, snp_subset, samfile_output, probability):

    # По умолчанию создаем гетерозиготную мутацию (вероятность равна 0.05)
    if probability is None:
        probability = 0.5

    rs = read.query_sequence
    rq = read.query_qualities
    rc = read.cigartuples

    # Для каждого целевого нуклеотида случайным образом выбираем, будет ли он мутирован
    for pos in iters:
        hetero = np.random.choice([0,1], size=1, p=[1 - probability, probability])

        #print(read.query_name, pos, number_list, read.query_sequence,
              #len(read.query_sequence) - len(number_list), read.cigarstring, read.cigartuples)
        #print(read.query_alignment_sequence)
        #print()

        # Если нет -- переходим к следующей позиции
        if hetero == 0:
            #print("No")
            continue

        # Если да -- меняем нуклеотид
        else:
            #print("Yes")
            # Определяем позицию в риде и индекс в csv
            ind = number_list.index(pos)
            indq = ind
            length = 0
            #print(ind)
            for cigar_block in rc:
                if cigar_block[0] in [0, 7, 8]:
                    length += cigar_block[1]
                    if indq < length:
                        break
                elif cigar_block[0] in [1, 4]:
                    length += cigar_block[1]
                    indq += cigar_block[1]
                elif cigar_block[0] == 2:
                    length += cigar_block[1]
                    indq -= cigar_block[1]

            if indq < 0:
                print("Warning", read.query_name)

            number = snp_subset.index[snp_subset["position"] == pos][0]

            # Меняем данные
            rs = rs[:indq] + snp_subset["nucleotide"].loc[number] + rs[indq + 1:]
            rq[ind] = snp_subset["quality"].loc[number]

    # Переприсваиваем последовательность и качество
    read.query_sequence = rs
    read.query_qualities = rq
    samfile_output.write(read)

snips = pd.read_csv(snp_file)

with pysam.AlignmentFile(target_reads, "rb") as samfile_input, pysam.AlignmentFile(output_bam, "wb", template=samfile_input) as samfile_output:
    chroms1 = samfile_input.references # При расчленении бама хедер со всеми хромосомами (даже если их нет) уходит в дочерний таргетный файл
    chroms2 = snips["chromosome"]
    chroms = list(set(chroms1) & set(chroms2))
    #print(chroms)
    for chrom in chroms:
        #print(chrom)
        snp_subset = snips.query("chromosome==@chrom")
        chrom_snps = snp_subset["position"]
        #print(chrom_snps)
        first_p, last_p = chrom_snps.min(), chrom_snps.max()
        reads = samfile_input.fetch(chrom, first_p, last_p + 1)
        
        for read in reads:
            number_list = read.get_reference_positions()
            iters = list(set(chrom_snps) & set(number_list))
            if len(iters) == 0:
                samfile_output.write(read)
            else:
                prob_mutation(read, iters, number_list, snp_subset, samfile_output, probability)
