import glob

from BamClass import BamClass
from CnvClass import CnvClass
from Config import config
import subprocess
import pandas as pd
import pysam

from SnvClass import SnvClass


def log_function_name(func):
    def wrapper(*args, **kwargs):
        class_name = args[0].__class__.__name__
        function_name = func.__name__
        print(f"Function {function_name} in class {class_name} is called.")
        return func(*args, **kwargs)

    return wrapper


class Runner(BamClass):
    def __init__(self, args):
        super().__init__(args)
        self.args = args
        # run steps of pipeline here
        self.run_pipeline()
        # self.bam_file.close()

    # todo: make assert mutationfile names == bam names

    def run_pipeline(self):
        self.get_chrom_names()
        if self.args.start < 1:  # if start == 1 -- start from bam splitting and skip these lines:
            self.sort_bam(self.bam_file, self.sorted_bam_file_name)
            self.index_bam(self.sorted_bam_file_name)
        if self.args.start < 2:  # if start == 2 -- start from bam select region of interest and skip these lines:
            self.split_bam()
        if self.args.start < 3:  # if start == 3 -- start from bam insert modification and skip these lines:
            self.select_roi()
        if self.args.start < 4:  # if start == 4 -- start from bam merging and skip these lines:
            self.insert_modification()
        self.merge_by_chromosome()
        #self.merge_bam()


    @log_function_name
    def sort_bam(self, bam, sorted_bam):
        command = f'samtools sort {bam} -o {sorted_bam} -@ {self.threads}'
        subprocess.run(command, shell=True)

    @log_function_name
    def index_bam(self, bam):
        command = f'samtools index {bam} -@ {self.threads}'
        subprocess.run(command, shell=True)

    @log_function_name
    def get_chrom_names(self):
        command = f"samtools idxstats {self.bam_file} | cut -f1"
        self.chrom_names = subprocess.run(command, capture_output=True, shell=True).stdout.decode().split('\n')

    @log_function_name
    def split_bam(self):
        for chr_name in self.chrom_names:
            if chr_name != '*':  # TODO: how to avoid parsing *.bam? add \*.bam?
                command = f'mkdir -p ./TMP; samtools view -b {self.sorted_bam_file_name} {chr_name} -@ {self.threads}> ./TMP/chr_{chr_name}.bam'
                subprocess.run(command, shell=True)

    @log_function_name
    def select_roi(self):
        """
        Walk through SNV coordinates, select chrom file, select region of interest, save into matched and unmatched bed file
        :return:
        """
        for chr_name in set(self.mutation_file['chromosome'].values.tolist()):
            chrom_subset = self.mutation_file.query('chromosome==@chr_name')
            if len(chrom_subset) == 0:
                raise NameError('Check chromosome names!!')
            if self.mutation_type == 'SNV':
                for snip_num in chrom_subset.index:
                    mut_site = snp_subset.loc[snip_num, "position"]
                    start, stop = mut_site, mut_site + 1
                    bed_string += [f'{chr_name}    {start}    {stop}']
                bed_string = "\n".join(bed_string)
            else:
                assert self.mutation_type == 'CNV', 'Mutation type MUST be one of: SNV or CNV!'
                positions = chrom_subset['position'].values.tolist() + chrom_subset['position_end'].values.tolist()
                start, stop = min(positions)-1, max(positions)-1
            command = (
                f'echo "{chr_name}    {start}    {stop}" > ./TMP/coords.bed;'
                f"samtools view -@ {self.threads} -L ./TMP/coords.bed -U ./TMP/chr_{chr_name}_unmatched.bam -o ./TMP/chr_{chr_name}_matched.bam ./TMP/chr_{chr_name}.bam;"
                f"rm ./TMP/chr_{chr_name}.bam")
            subprocess.run(command, shell=True)
        pass

    @log_function_name
    def insert_modification(self):
        for chr_name in set(self.mutation_file['chromosome'].values.tolist()):
            matched_bam = f"./TMP/chr_{chr_name}_matched.bam"
            self.index_bam(matched_bam)
            if self.mutation_type == 'SNV':
                SnvClass(self.args, matched_bam)
            else:
                CnvClass(self.args, matched_bam)
            self.sort_bam('./TMP/temp_out.bam', './TMP/temp_out_sorted.bam')
            command = f'mkdir -p ./OUT/; rm ./TMP/temp_out.bam; mv ./TMP/temp_out_sorted.bam ./OUT/chr_{chr_name}_matched.bam'
            subprocess.run(command, shell=True)

    @log_function_name
    def merge_bam(self):
        command = f'samtools merge -f -@ {self.threads} -o {self.bam_file.strip("bam")}.edited.bam ./TMP/chr_*ed.bam'
        subprocess.run(command, shell=True)

    @log_function_name
    def merge_by_chromosome(self):
        for chr_name in set(self.mutation_file['chromosome'].values.tolist()):
            matched_bam = f"./OUT/chr_{chr_name}_matched.bam"
            unmatched_bam = f"./TMP/chr_{chr_name}_unmatched.bam"
            self.index_bam(matched_bam)
            self.index_bam(unmatched_bam)
            command = f'samtools merge -f -@ {self.threads} -o ./OUT/chr_{chr_name}.edited.bam {matched_bam} {unmatched_bam}'
            subprocess.run(command, shell=True)
            self.index_bam(f'./OUT/chr_{chr_name}.edited.bam')

