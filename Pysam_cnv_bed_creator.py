import pysam
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("input_bam", type=str)
parser.add_argument("delete_bed", type=str)
parser.add_argument("cnv_file", type=str)
arguments = parser.parse_args()
input_bam = arguments.input_bam
delete_bed = arguments.delete_bed
cnv_file = arguments.snp_file


cnvs = pd.read_csv(cnv_file)

with pysam.AlignmentFile(input_bam, "rb") as samfile_input, open(delete_bed, "w") as bed:
    chroms = list(set(samfile_input.references) & set(cnvs["chromosome"]))
    if len(chroms) == 0:
        # Поменять на raise
        print("Target chromosomes are absent in reference\nPlease, check correctness of your csv file or names of contigs")
    else:
        for chr_name in chroms:  # Проходимся по каждой хромосоме, которая есть и в выравнивании, и в csv
            chrom_subset = cnvs.query("chromosome==@chr_name")
            for cnv_num in chrom_subset.index:
                start, stop = chrom_subset.loc[cnv_num, "position_start"], chrom_subset.loc[cnv_num, "position_finish"]
                bed.write(f"{chr_name}\t{start}\t{stop}\n")