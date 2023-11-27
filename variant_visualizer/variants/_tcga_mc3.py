import pandas as pd
from .._config import config
from .. import core, setup
from collections import Counter
import os
from ._variants import Variant
from ..setup._setup_tcga_mc3 import load_mc3


tss_code_study_table = pd.read_csv(
    os.path.join(config['tcga_mc3']['tcga_code_tables'],'tissueSourceSite.tsv'),
    sep='\t')

study_disease_table = pd.read_csv(
    os.path.join(config['tcga_mc3']['tcga_code_tables'],'diseaseStudy.tsv'),
    sep='\t')

def _select_tcga_samples(tcga_sample_class):
    """
    Description
    ---
    Returns the chosen subset of tcga samples of the single base 
    substitution variantsal signatures specified in config.yml at
    'tcga_sbs_path'.

    Parameters
    ---
    tcga_sample_class : str
        'any' for all samples
        'mss' for micro satelite stable samples
        'msi' for micro satelite instable samples
        'pole' for Polymerase-epsilon muated samples
        'msi_pole' for msi or pole samples
        'mss_no_SKCM' for mss samples, excluding SKCM samples
        'mss_no_SKCM_DLBC' for mss samples, excluding SKCM and DLBC samples
    """
    t = _load_tcga_sbs()

    if tcga_sample_class == 'any':
        samples = t['Sample Names'].tolist()
    elif tcga_sample_class == 'mss':
        samples = t.loc[(t['MSI SBS'] == 0) & (
            t['PolE SBS'] == 0), 'Sample Names'].tolist()
    elif tcga_sample_class == 'msi':
        samples = t.loc[(t['MSI SBS'] > 0)&
                        (t['PolE SBS'] == 0), 'Sample Names'].tolist()
    elif tcga_sample_class == 'pole':
        samples = t.loc[(t['MSI SBS'] == 0)&
                        (t['PolE SBS'] > 0), 'Sample Names'].tolist()
    elif tcga_sample_class == 'msi_pole':
        samples = t.loc[(t['MSI SBS'] > 0)&
                        (t['PolE SBS'] > 0), 'Sample Names'].tolist()
    elif tcga_sample_class == 'mss_no_SKCM':
        samples = _load_tcga_sbs('mss')
        samples = [s for s in samples if get_tcga_cancer_type(s) != 'SKCM']
    elif tcga_sample_class == 'mss_no_SKCM_DLBC':
        samples = _load_tcga_sbs('mss_no_SKCM')
        samples = [s for s in samples if get_tcga_cancer_type(s) != 'DLBC']
    else:
        raise ValueError(
            f'Unknown tcga_sample_class {tcga_sample_class}. Valid options: any, mss, msi, pole, msi/pole')

    return samples


def _load_tcga_sbs():
    """Load the cosmic SBS signatures from the file specified in the config.yml under key: tcga_sbs_path."""
    sbs_df = pd.read_csv(config['tcga_mc3']['tcga_sbs_path'], sep=",")

    cosmic_mut_signatures = [c for c in sbs_df.columns if 'SBS' in c]
    cosmic_msi_mut_signatures = config['tcga_mc3']['cosmic_msi_sbs']
    cosmic_pole_mut_signatures = config['tcga_mc3']['cosmic_pole_sbs']

    sum_sbs = pd.DataFrame(sbs_df, columns=cosmic_mut_signatures).sum(axis=1).values
    sbs_df.insert(loc=0,
                   column='Total SBS',
                   value=sum_sbs)

    sbsMsi = pd.DataFrame(
        sbs_df, columns=cosmic_msi_mut_signatures).sum(axis=1).values
    sbs_df.insert(loc=1,
                   column='MSI SBS',
                   value=sbsMsi)

    sbsPole = pd.DataFrame(
        sbs_df, columns=cosmic_pole_mut_signatures).sum(axis=1).values
    sbs_df.insert(loc=2,
                   column='PolE SBS',
                   value=sbsPole)
    return sbs_df



def load_mc3_slice(cluster: int):
    """
    Load a slice of the converted mc3 file. Slices can be generated with slice_mc3_cluster().
    Default values load cluster slices from the directory where they are saved by default.
    """
    path = setup.get_out_path(
        input_file=config['tcga_mc3']['mc3_GRCh38_maf'],
        out_dir=os.path.join(config['general']['data_dir'],
                            'cluster_inputs'),
        file_ending=f'.{cluster}_cluster.maf')

    mc3_slice = pd.read_csv(path,
                            index_col=0,
                            sep="\t",
                            dtype={'Chromosome': 'str'})
    return mc3_slice  

def get_row_variant(i: int, row: pd.Series):
    reference = core.get_reference(
    reference_type='genomic',
    chromosome=str(row['Chromosome']),
    strand=row['Strand']
    )

    variant_type = row['Variant_Type']
    consequence = row['Consequence']
    if pd.isna(variant_type):
        variant_type = ''
    if pd.isna(consequence):
        consequence = ''
    return Variant(
        start=row['Start_Position'],
        end=row['End_Position'],
        reference=reference,
        variant_type=variant_type,
        consequence=consequence,
        sample_id=row['Tumor_Sample_Barcode'],
        normal_sample_id=row['Matched_Norm_Sample_Barcode'],
        disease=get_tcga_cancer_type(row['Tumor_Sample_Barcode']),
        source=f'mc3_row:{i}',
        ref_allele=row['Reference_Allele'],
        alt_allele_1=row['Tumor_Seq_Allele1'],
        alt_allele_2=row['Tumor_Seq_Allele2'],
        label = f'{row["Reference_Allele"]} > {row["Tumor_Seq_Allele1"]}, {row["Tumor_Seq_Allele2"]}; {consequence}'
    
        )

def get_mc3_variants(mc3_slice: pd.DataFrame) -> list:
    """Returns all variants in section of mc3"""

    out = []
    for i,row in mc3_slice.iterrows():
        out.append(get_row_variant(i, row))

    return out

def get_region_mc3_variants(region: core.BioRegion, variants: list):
    """Return all variants from the list that overlap with the region."""
    out = []
    variants = [m for m in variants if m.reference.chromosome == region.reference.chromosome]
    for m in sorted(variants, key=lambda x: x.start):
        if m.start > region.end:
            break
        if m.overlaps(region):
            out.append(m)
    return out

def get_region_variant_rates(region: core.BioRegion, region_variants: list) -> dict:
    """
    Description
    ---
    Screens a Genomic Region for intersecting variants. 
    Calculates variant rates for each variantal variant type represented in the region, and the overall variant rate.

    region : BioRegion
    region_variants : list
                       List of variants that are overlapping the region
                       Can be generated using get_region_variants()
    """

    out = {'n_variants': len(region_variants),
           'combined_variant_rate': len(region_variants)/region.get_length()}

    type_counter = Counter([m.group for m in region_variants])
    for m_type in type_counter.keys():
        out.update({f'n_{m_type}_variants': type_counter.get(m_type),
                    f'{m_type}_variant_rate': type_counter.get(m_type)/region.get_length()})
        
    return out

def get_tcga_cancer_type(tcga_sample_barcode):
    '''
    Returns the cancer type for a tcga_sample_barcode
    Returns str
    If tss_code is NA, returns ''
    '''
    id = tcga_sample_barcode

    tss_code = id.split('-')[1]
    if tss_code == 'NA':
        return ''
    else:
        study = tss_code_study_table.loc[tss_code_study_table['TSS Code']==tss_code, 'Study Name'].values[0]
        study_abbrev = study_disease_table.loc[study_disease_table['Study Name']==study, 'Study Abbreviation'].values[0]

    return study_abbrev

def get_region_slice(region, mc3):
    mc3_slice = mc3.loc[(mc3['Chromosome'] == region.reference.chromosome)&
                        (mc3['Start_Position'] >= region.start)&
                        (mc3['End_Position'] <= region.end), :]

    return mc3_slice