from ... import core
from ..._config import config
from ...io import utils as io_utils
from ._mirna_binding import MiRNABinding
import pandas as pd

def load_targetscan_db(path=config['targetscan']['target_locations_GRCh38_bed']) -> pd.DataFrame:
    """Load GRCh38 converted bed-file with human miRNA binding sites to a pandas Dataframe."""
    return io_utils.load_converted_bed(path=path)

def get_targetscan_mirna_binding(bio_region: core.BioRegion, targetscan_db: pd.DataFrame) -> list:
    """Get a list of BioRegions containing miRNA binding sites in the given region."""

    if not isinstance(bio_region.reference, core.GenomicReference):
        raise TypeError(f'BioRegion.reference must be GenomicReference.')

    db_slice = targetscan_db.loc[(targetscan_db[1] >= bio_region.start)&
                            (targetscan_db[2] <= bio_region.end)&
                            (targetscan_db[0].isin([bio_region.reference.chromosome, 'chr'+bio_region.reference.chromosome])), :]
    out = []
    for i, row in db_slice.iterrows():
        reference = core.get_reference(reference_type='genomic', 
                                   chromosome=io_utils.remove_chr(chromosome=row[0]),
                                   strand = row[5])
        out.append(MiRNABinding(
            start=row[1],
            end=row[2],
            reference=reference,
            source=f'targetscan_target_locations_row:{i}',
            binding=row[3].split(":")[-1],
            score=float(row[4])
            ))
    return out