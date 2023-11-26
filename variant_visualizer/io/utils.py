import os, subprocess, dill
import pandas as pd
from .._config import config

def remove_chr(chromosome: str) -> str:
    """Return the chromosome string without 'chr' leading characters."""
    if chromosome[0:3] == 'chr':
                return chromosome[3:]
    else:
        return chromosome

def add_chr(chromosome: str) -> str:
    """Return the chromosome string with 'chr' as leading characters."""
    if chromosome[0:3] != 'chr':
        return 'chr'+chromosome
    else:
        return chromosome
    
def get_out_path(input_file, out_dir, file_ending=''):
    """Return a path where the file format from the input file is exchanged to file_ending, and the directory is out_dir."""
    name = os.path.splitext(os.path.basename(input_file))[0]
    outName = os.path.join(out_dir, name)
    outName = os.path.realpath(outName)
    return str(outName + file_ending)

def _load_bed(path, header=None):
    """Return bed as pandas.DataFrame"""
    return pd.read_csv(path, sep="\t",
                        header=header,
                        dtype={0: 'str'})

def load_converted_bed(path, header=None):
    """Load a bed-file to a pandas Dataframe and convert coordinates to base-1, start and end inclusive."""
    bed = _load_bed(path=path,
                    header=header)
    bed[0] = [remove_chr(c) for c in bed[0]]
    bed[1] = [start+1 for start in bed[1]]
    return bed

def liftover_GRCh37_GRCh38(input_bed_path, output_bed_path):
     # Liftover bed to GRCh38
    subprocess.run(f'{config["tcga_mc3"]["ucsc_liftover"]} {input_bed_path} {config["tcga_mc3"]["ucsc_liftover_chainfile"]} {output_bed_path} {output_bed_path}.unmap', 
                    shell=True, 
                    check=True, 
                    text=True)
    return

def dill_dump_object(out_path: str, object: object):
    """Save any object as binary file at the given path."""
    with open(out_path, 'wb') as f:
        f.write(dill.dumps(obj=object))
    return

def dill_load_object(object_path: str):
     """Load any object from a binary file, generated with dill_dump_object."""
     with open(object_path, 'rb') as f:
          out = dill.load(f)
     return out