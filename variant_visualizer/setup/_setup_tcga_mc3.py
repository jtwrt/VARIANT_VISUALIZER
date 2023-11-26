from .._config import config
import pandas as pd
from . import _setup_utils as s_utils
from .. import genome_annotation as ga
import os

def liftover_mc3() -> str:
    """
    Description
    ---
    Convert the mc3 from reference GRCh37 to GRCh38 using UCSC Liftover.
    Uses the mc3_GRCh37_maf file specified in config.yml and saves the converted file
    at the config data_dir per default.
    Returns the output path.
    """    
    out_dir=config['general']['data_dir']
    mc3_path=config['tcga_mc3']['mc3_GRCh37_maf']

    # Load original mc3_maf
    mc3 = pd.read_csv(mc3_path,
                       sep="\t",
                       dtype={'Chromosome': 'str'})
    
    # Convert to bed and save
    bed = pd.DataFrame()
    bed.insert(loc=0, column='Chromosome', value=[s_utils.add_chr(c) for c in mc3['Chromosome']])
    start = [int(s)-1 for s in mc3['Start_Position'].tolist()]
    end = [int(e) for e in mc3['End_Position'].tolist()]
    bed.insert(loc=1, column='Start', value=start)
    bed.insert(loc=2, column='End', value=end)
    index = [i for i in range(len(bed))]
    context = mc3['CONTEXT'].tolist()
    name = [f"{i}:{c}" for i, c in zip(index, context)]
    bed.insert(loc=3, column='Name', value=name)
    bed.insert(loc=4, column='Score', value=0)
    bed.insert(loc=5, column='Strand', value=mc3['Strand'].tolist())
    tmp_name = s_utils.get_out_path(
        input_file=mc3_path,
        out_dir=config['general']['tmp_dir'])
    mc3_GRCh37_bed = f'{tmp_name}.GRCh37.bed'
    bed.to_csv(mc3_GRCh37_bed,
               sep="\t",
               header=False,
               index=False)
    
    # Liftover bed to GRCh38
    mc3_GRCh38_bed = f'{tmp_name}.GRCh38.bed'
    s_utils.liftover_GRCh37_GRCh38(
        input_bed_path=mc3_GRCh37_bed,
        output_bed_path=mc3_GRCh38_bed
        )
        
    # Create new mc3_maf with successfully converted coordinates
    mc3_GRCh38_locs = (s_utils.load_converted_bed(mc3_GRCh38_bed)[lambda x: x[3].isin(name)])

    columns = dict(Chromosome=[],
                   Start_Position=[],
                   End_Position=[],
                   )
    mc3_index =[]
    for chr_GRCh38, start_GRCh38, end_GRCh38, mc3_id in mc3_GRCh38_locs[[0, 1, 2, 3]].values:
        columns['Chromosome'].append(chr_GRCh38)
        columns['Start_Position'].append(start_GRCh38)
        columns['End_Position'].append(end_GRCh38)
        mc3_index.append((int(mc3_id.split(':')[0])))

    columns['Chromosome'] = [s_utils.remove_chr(c) for c in columns['Chromosome']]
    columns['NCBI_Build'] = ['GRCh38' for _ in mc3_index]
    for key in mc3.columns:
        if key in ['Chromosome', 'Start_Position', 'End_Position', 'NCBI_Build']: 
            continue
        else:
            columns[key] = mc3.loc[mc3_index, key].tolist()

    out_path = s_utils.get_out_path(input_file=mc3_path,
                                out_dir=out_dir,
                                file_ending='.GRCh38.maf')
    
    columns['mc3_GRCh37_index'] = mc3_index
    pd.DataFrame(columns).to_csv(out_path,
                                  sep="\t",
                                  header=True,
                                  index=False)
    
    return out_path

def load_mc3():
    """Load the GRCh38 converted mc3 file specified in the config.yml. Convert mc3 with tcga_mc3.setup.liftover_mc3()."""
    mc3_path=config['tcga_mc3']['mc3_GRCh38_maf']
    mc3_GRCh38 = pd.read_csv(mc3_path,
                              sep="\t",
                              dtype={'Chromosome': 'str'})
    mc3_GRCh38.set_index('mc3_GRCh37_index',
                         inplace=True)
    
    return _clean_mc3(mc3_GRCh38)

def slice_mc3_cluster(cluster: int, mc3: pd.DataFrame, out_path=None) -> str:
    """
    Description
    ---
    Save a slice of the mc3 separately, which includes rows 
    that intersect the specified .gtf.GtfCluster.
    Returns the output Path of the file created.
    """
    
    cluster_region = ga.GtfCluster(cluster_id=cluster)
    mc3_slice = get_region_slice(region=cluster_region, 
                                 mc3=mc3)
    
    if out_path is None:
        out_dir = os.path.join(config['general']['data_dir'],
                            'cluster_inputs')
        os.makedirs(out_dir, exist_ok=True)
        out_path = s_utils.get_out_path(input_file=config['tcga_mc3']['mc3_GRCh38_maf'],
                                    out_dir=out_dir,
                                    file_ending=f'.{cluster}_cluster.maf')

    mc3_slice = _clean_mc3(mc3_slice)
    mc3_slice.to_csv(out_path, sep='\t')
    return out_path

def _clean_mc3(mc3):
    """Remove all rows from the mc3, that have the value 'nonpreferredpair' in the FILTER column."""
    return mc3.loc[mc3['FILTER']=='PASS']

def get_region_slice(region, mc3):
    mc3_slice = mc3.loc[(mc3['Chromosome'] == region.reference.chromosome)&
                        (mc3['Start_Position'] >= region.start)&
                        (mc3['End_Position'] <= region.end), :]

    return mc3_slice
