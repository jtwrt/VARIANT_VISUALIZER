from ... import core
from ..._config import config
from ...io import utils as io_utils
from ._pas import PAS
import pandas as pd


def load_pas_atlas_db(pas_atlas_path=config['pas_atlas']['pas_atlas_path']) -> pd.DataFrame:
    return io_utils.load_converted_bed(path=pas_atlas_path,
                                     header=None
                                     )

def get_pas_atlas_pas(region: core.BioRegion, pas_atlas: pd.DataFrame) -> list:
    """Get a list of PAS BioRegions in the given region."""


    if not isinstance(region.reference, core.GenomicReference):
        raise ValueError(f'Region reference GenomicReference.')
        
    # subset pas_atlas for rows that can include pas hitting the region
    pas_atlas = pas_atlas.loc[
        (pas_atlas[0] == region.reference.chromosome)&
        (pas_atlas[1] >= region.start)&
        (pas_atlas[2] <= region.end)&
        (pd.notna(pas_atlas[10]))
        ]

    out = []
    for i,row in pas_atlas.iterrows():
        out.extend(_get_pas_in_row(i,row))        
    return out

def _get_pas_in_row(i: int, row) -> list:
    """
    Create pas BioRegion from i,row in pandas.DataFrame.iterrows()
    """
    if pd.isna(row[1]):
        raise ValueError('NA value in given row column 10')

    chromosome = str(row[0])
    start = int(row[1])
    end = int(row[2])
    strand = row[5]
    _ = core.check_strand(strand)
    signals = str(row[10])
    source = f'pas_atlas_row:{i}'

    reference = core.get_reference(
        reference_type='genomic',
        chromosome=chromosome,
        strand=strand)
    cleavage_site = core.BioRegion(
        start=start,
        end=end,
        reference=reference)
    
    out = []
    for signal in signals.split(';'):
        if strand == '+':
            pasStart = int(signal.split("@")[2])
            pasEnd = int(signal.split("@")[2]) + 5  # end inclusive
        elif strand == '-':
            pasStart = int(signal.split("@")[2]) - 5
            pasEnd = int(signal.split("@")[2])
                    
        out.append(PAS(
            start=pasStart,
            end=pasEnd,
            reference=reference,
            sequence=signal.split('@')[0],
            source=source,
            cleavage_site=cleavage_site,
            label=signal
            )) 
    return out