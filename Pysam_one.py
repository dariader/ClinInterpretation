#!/usr/bin/env python

# Importing packages
import pysam
import array
import argparse
import numpy as np
import pandas as pd


# Parsing arguments
parser = argparse.ArgumentParser()
parser.add_argument("input_bam", type=str)
parser.add_argument("output_bam", type=str)
parser.add_argument("snp_file", type=str)
arguments = parser.parse_args()
input_bam = arguments.input_bam
output_bam = arguments.output_bam
snp_file = arguments.snp_file


# Executing the code
snips = pd.read_csv(snp_file)

with pysam.AlignmentFile(input_bam, "rb") as samfile_input, pysam.AlignmentFile(output_bam, "wb", template=samfile_input) as samfile_output:
    
    # Инициализируем список хромосом
    chroms = samfile_input.references
    
    for chrom in chroms:
        
        # Если в хромосоме не планируются SNP, не трогаем ее
        if chrom not in set(snips["chromosome"]):
            reads = samfile_input.fetch(chrom)
            for read in reads:
                samfile_output.write(read)
                
        # А если планируются, трогаем
        else:
            
            # Отбираем SNP для данной хромосомы
            snp_subset = snips.query("chromosome==@chrom")
            chrom_snps = snp_subset["position"]
            first_p, last_p = chrom_snps.min(), chrom_snps.max()
            
            # Риды до SNP-позиций
            reads = samfile_input.fetch(chrom, 0, first_p)
            for read in reads:
                
                # Страховка для ситуаций, где картировался только один рид из пары
                if read.reference_end is None:
                    if (read.reference_start + read.query_length - 1) < first_p:
                        samfile_output.write(read)
                else:
                    if read.reference_end < first_p:
                        samfile_output.write(read)             
            
            # Риды, содержащие SNP-позиции
            reads = samfile_input.fetch(chrom, first_p, last_p + 1)
            for read in reads:
                number_list = read.get_reference_positions()
                
                # Сохраняем в переменные данные исходных ридов
                rs = read.query_sequence
                rq = read.query_qualities
                rc = read.cigarstring
                        
                for number in range(len(chrom_snps)):
                    if chrom_snps.iloc[number] in number_list:
                        
                        # Вычисляем индекс интересующего нуклеотида
                        ind = number_list.index(chrom_snps.iloc[number])
                                       
                        # Создаем данные модифицированных ридов
                        rs = rs[:ind] + snp_subset["nucleotide"].iloc[number] + rs[ind + 1:]
                        rq[ind] = snp_subset["quality"].iloc[number]
                
                # Изменяем парметры ридов, вносим в новый файл
                read.query_sequence = rs
                read.query_qualities = rq           
                samfile_output.write(read)
                
            # Риды после SNP-позиций
            reads = samfile_input.fetch(chrom, last_p + 1)
            for read in reads:
                if read.reference_start > last_p:
                    samfile_output.write(read)
        
        # Записываем невыровнявшиеся риды в конец
        reads = samfile_input.fetch("*")
        for read in reads:
            samfile_output.write(read)
