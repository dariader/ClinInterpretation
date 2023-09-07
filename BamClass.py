from Config import config
import subprocess
import pandas as pd
import pysam

def log_function_name(func):
    def wrapper(*args, **kwargs):
        class_name = args[0].__class__.__name__
        function_name = func.__name__
        print(f"Function {function_name} in class {class_name} is called.")
        return func(*args, **kwargs)
    return wrapper


class BamClass:
    def __init__(self, args):
        self.bam_file = args.bam_file
        self.mutation_file = pd.read_csv(args.mutation_file, converters={'position': int})
        self.sorted_bam_file_name = f"{self.bam_file}_sorted"
        self.threads = args.threads
        self.chrom_names = None
        ## run steps of pipeline here
        self.sort_bam()
        self.index_bam()
        self.get_chrom_names()
        self.split_bam()
        self.select_roi()
        self.bam_file.close()

    @log_function_name
    def sort_bam(self):
        command = f'samtools sort {self.bam_file} -o {self.sorted_bam_file_name} -@ {self.threads}'
        subprocess.run(command, shell=True)

    @log_function_name
    def index_bam(self):
        command = f'samtools index {self.sorted_bam_file_name} -@ {self.threads}'
        subprocess.run(command, shell=True)

    @log_function_name
    def get_chrom_names(self):
        command = f"samtools idxstats {self.sorted_bam_file_name} | cut -f1"
        self.chrom_names = subprocess.run(command, capture_output=True, shell=True).stdout.decode().split('\n')

    @log_function_name
    def split_bam(self):
        for chr_name in self.chrom_names:
            if chr_name != '*':
                command = f'mkdir -p ./TMP; samtools view -b {self.sorted_bam_file_name} {chr_name} -@ {self.threads}> ./TMP/{chr_name}.bam'
                subprocess.run(command, shell=True)

    @log_function_name
    def select_roi(self):
        """
        Walk through SNV coordinates, select chrom file, select region of interest, save into matched and unmatched bed file
        :return:
        """
        for chr_name in self.mutation_file['chromosome']:
            if chr_name != '*':
                for chrom_subset in self.mutation_file.query('chromosome=@chr_name'):
                    start, stop = min(chrom_subset['position']), max(chrom_subset['position'])
                    command = (f"samtools view --threads {self.threads} -L - -U ./TMP/{chr_name}_unmatched.bam -o ./TMP/{chr_name}_matched.bam ./TMP/{chr_name}.bam < {chr_name}\t{start}\{stop}; rm ./TMP/{chr_name}.bam")
                    subprocess.run(command, shell=True)
        pass

    def insert_modification(self):
        pass

    def write_chm_tofile(self):
        pass

    def merge_bam(self):
        """
        samtools index new_bam_order_mod.bam

        samtools cat --threads 8 genome38_192021_untarget.bam new_bam_order_mod.bam -o genome38_192021_mod.bam

        samtools sort --threads 8 genome38_192021_mod.bam -o genome38_192021_mod_sorted.bam

        samtools index -@ 8 genome38_192021_mod_sorted.bam
        :return:
        """
        command = f'pass'
        subprocess.run(command)
        pass
