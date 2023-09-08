from Config import config
import subprocess
import pandas as pd
import pysam


class BamClass:
    def __init__(self, args):
        self.bam_file = args.bam_file
        self.mutation_file = pd.read_csv(args.mutation_file, converters={'position': int, 'chromosome': str})
        self.mutation_type = args.mutation_type
        self.sorted_bam_file_name = f"{self.bam_file}_sorted"
        self.threads = args.threads
        self.chrom_names = None
