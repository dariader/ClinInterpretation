import os

import yaml
from dataclasses import dataclass


@dataclass
class Config:
    bam_file: str
    mutation_file: str
    mutation_type: str
    num_threads: int
    mem_size: int


def make_paths_absolute(base_dir, paths):
    return {key: os.path.join(base_dir, value) for key, value in paths.items()}


def load_config_from_yaml(yaml_file):
    with open(yaml_file, "r") as f:
        config_data = yaml.safe_load(f)
        config = Config(bam_file=config_data.get('bam_file'),
                        mutation_file=config_data.get('mutation_file'),
                        mutation_type=config_data.get('mutation_type'),
                        num_threads=config_data.get('num_threads'),
                        mem_size=config_data.get('mem_size'))
    return config


# Specify your base directory here
base_directory = "/home/daria/PycharmProjects/InSilico_test/"  # change it if necessary
config = load_config_from_yaml(f'{base_directory}/config.yaml')
