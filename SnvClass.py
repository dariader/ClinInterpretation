import pysam
import pandas as pd
import numpy as np

from BamClass import BamClass


class SnvClass(BamClass):
    """
    This class inserts mutations, and generates temporary file per chromosome it was run onto
    """
    def __init__(self, args, matched_bam):
        super().__init__(args)
        self.input_bam = self.bam_file
        self.output_bam = './TMP/temp_out.bam'  # constant name
        self.snp_file = self.mutation_file
        self.target_reads = matched_bam  # matched_reads, must be passed while class is initialised
        self.probability = 1  # constant value
        # Below adaptation of code written by @Bogdan
        # Variable names not changed, nothing changed except added `self` prefix where needed
        self.run_mutation()


    def prob_mutation(self, read, iters, number_list, snp_subset, samfile_output, probability):
        # По умолчанию создаем гетерозиготную мутацию (вероятность равна 0.05)
        if probability is None:
            probability = 0.5

        rs = read.query_sequence
        rq = read.query_qualities
        rc = read.cigartuples

        # Для каждого целевого нуклеотида случайным образом выбираем, будет ли он мутирован
        for pos in iters:
            hetero = np.random.choice([0, 1], size=1, p=[1 - probability, probability])
            # print(read.query_name, pos, number_list, read.query_sequence,
            # len(read.query_sequence) - len(number_list), read.cigarstring, read.cigartuples)
            # print(read.query_alignment_sequence)
            # print()

            # Если нет -- переходим к следующей позиции
            if hetero == 0:
                continue

            # Если да -- меняем нуклеотид
            else:
                # Определяем позицию в риде и индекс в csv
                ind = number_list.index(pos)
                indq = ind - 1   # Minus 1 is from SAM-BAM convertation (where SAM is 1-based and BAM is 0-based)
                length = 0
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
                rq[indq] = snp_subset["quality"].loc[number]

        # Переприсваиваем последовательность и качество
        read.query_sequence = rs
        read.query_qualities = rq
        samfile_output.write(read)

    def run_mutation(self):
        snips = self.snp_file
        with pysam.AlignmentFile(self.target_reads, "rb") as samfile_input, pysam.AlignmentFile(self.output_bam, "wb",
                                                                                                template=samfile_input) as samfile_output:
            chroms1 = samfile_input.references  # При расчленении бама хедер со всеми хромосомами (даже если их нет) уходит в дочерний таргетный файл
            chroms2 = snips["chromosome"].values.tolist()
            chroms = list(set(chroms1) & set(chroms2))
            # print(chroms)
            if len(chroms) == 0:
                raise ValueError('Check the data! Chrom names in file and in bam do not match!')
            for chrom in chroms:
                # print(chrom)
                chrom = chrom
                snp_subset = snips.query("chromosome==@chrom")
                chrom_snps = snp_subset["position"]
                # print(chrom_snps)
                first_p, last_p = chrom_snps.min(), chrom_snps.max()
                reads = samfile_input.fetch(chrom, first_p, last_p + 1)

                for read in reads:
                    number_list = read.get_reference_positions()
                    iters = list(set(chrom_snps) & set(number_list))
                    if len(iters) == 0:
                        samfile_output.write(read)
                    else:
                        self.prob_mutation(read, iters, number_list, snp_subset, samfile_output, self.probability)
