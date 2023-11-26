import subprocess
from .._config import config
import os, dill, yaml, gzip, zipfile, shutil, sys
import pandas as pd

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

def liftover_GRCh37_GRCh38(input_bed_path, output_bed_path):
     # Liftover bed to GRCh38
    subprocess.run(f'{config["tcga_mc3"]["ucsc_liftover"]} {input_bed_path} {config["tcga_mc3"]["ucsc_liftover_chainfile"]} {output_bed_path} {output_bed_path}.unmap', 
                    shell=True, 
                    check=True, 
                    text=True)
    return

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
    
def update_config(updated_config, path=os.path.realpath('config.yml')):
    """Add new values to the config.yml"""
    with open(path, 'r+') as f:
        # reset file position
        f.seek(0)
        # write the updated YAML to the file
        yaml.dump(updated_config, f)
        f.truncate()

def decompress_gz(gz_file):
    """
    Decompress file and save it at the same location, 
    without .gz ending. 
    Returns path of decompressed file.
    """
    print(f'Decompressing {gz_file}')
    out_path = os.path.splitext(gz_file)[0]
    with gzip.open(gz_file, 'rb') as f_in:
        with open(out_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(gz_file)
    return out_path

def decompress_zip(zip_file, unzip_dir=None):
    """
    Unzips a .zip file and extracts contents to its directory.
    Returns the directory path.
    Uses the basename of the .zip file to name the directory at its original location. 
    Alternatively, define unzip_dir to name the output_directory.
    """
    print(f'Decompressing {zip_file}')
    if unzip_dir is None:
        unzip_dir = os.path.basename(zip_file)
    new_dir = os.path.join(os.path.dirname(zip_file),
                           os.path.splitext(unzip_dir)[0])
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)
    with zipfile.ZipFile(zip_file, 'r') as f:
        f.extractall(new_dir)
    os.remove(zip_file)
    return new_dir

def block_print():
    """Disable print statements"""
    sys.stdout = open(os.devnull, 'w')

# Restore
def enable_print():
    """Enable print statements"""
    sys.stdout = sys.__stdout__