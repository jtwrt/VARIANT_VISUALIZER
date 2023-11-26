from .._config import config
import subprocess, warnings, os
import pandas as pd

def load_gtf(clustered=False, parse_attributes=True) -> pd.DataFrame:
    '''
    Load reference genome GTF
    Parses 'attribute' column to separate columns
    Returns GTF as data frame
    Designed to work with Ensemble Homo_sapiens.GRCh38.109.gtf
    '''
    if clustered == True:
        gtf_path = config['general']['clustered_gtf']
        header = None
        colnames = ['seqname', 'source', 'feature', 'start', 'end',
                    'score', 'strand', 'frame', 'attribute', 'cluster']
    else:
        gtf_path = config['general']['gtf']
        header=config['general']['gtf_header']
        colnames = ['seqname', 'source', 'feature', 'start',
                    'end', 'score', 'strand', 'frame', 'attribute']

    gtf = pd.read_csv(gtf_path,
                       sep='\t',
                       header=header,
                       names=colnames,
                       dtype={'seqname': 'str',
                              'attribute': 'str'})

    if parse_attributes != True:
        warnings.warn(
            'Loading gtf without splitting attribute column in separate columns. Column \'transcript_id\' is needed by functions using gtf as input.')
        return gtf

    attribute_type = ['gene_id', 'gene_version', 'gene_name', 'gene_source', 'gene_biotype',
                      'transcript_id', 'transcript_version', 'transcript_name', 'transcript_source', 'transcript_biotype',
                      'tag', 'ccds_id']

    parsed_attribute = {}
    for a in gtf['attribute']:
        split = a.split('; ')
        split[-1] = split[-1][:len(split[-1])-1]  # remove trailing ';'

        for t in attribute_type:
            tmp = ''
            for s in split:
                if s[:len(t)] == t:
                    # remove leading ' ' and quotation marks
                    parsed_s = s[len(t)+1:].replace("\"", "")
                    tmp = f'{tmp}{parsed_s};'

            if tmp != '':
                tmp = tmp[:len(tmp)-1]  # remove trailing ';'

            if parsed_attribute.get(t) == None:
                parsed_attribute[t] = [tmp]
            else:
                parsed_attribute[t].append(tmp)

    # replace attribute column from gtf with multiple parsed information columns
    del gtf['attribute']
    n = len(gtf.columns)
    for i, t in enumerate(attribute_type):
        gtf.insert(loc=n+i,
                   column=t,
                   value=parsed_attribute[t])
    return gtf

def cluster_gtf(min_gap=config['general']['gtf_cluster_min_gap']) -> str:
    """
    Description
    ---
    Cluster sequences in the gtf file together that are within min_gap nucleotides of each other using bedtools cluster -d min_gap.
    Uses the init_bedtools command and gtf path defined in the config.yml.
    Generates temporary and output files in the config tmp_dir and data_dir if not specified otherwise.
    """

    out_dir=config['general']['data_dir']
    tmp_dir=config['general']['tmp_dir']

    gtf_path = config['general']['gtf']
    init_bedtools = config['general']['init_bedtools']

    # Load gtf
    with warnings.catch_warnings():
        warnings.filterwarnings(
            'ignore', message='Loading gtf without splitting attribute column in separate columns. Column \'transcript_id\' is needed by functions using gtf as input.')
        gtf = load_gtf(clustered=False, parse_attributes=False)

    org_colnames = ['seqname', 'source', 'feature', 'start',
                    'end', 'score', 'strand', 'frame', 'attribute']

    # Save as bed
    gtf['start'] = [s-1 for s in gtf['start'].tolist()]
    colnames = ['seqname', 'start', 'end', 'feature',
                'score', 'strand', 'source', 'frame', 'attribute']
    gtf = gtf[colnames]  # reorder columns
    name = os.path.splitext(os.path.basename(gtf_path))[0]
    tmp_name = os.path.join(tmp_dir, name)
    gtf.to_csv(f'{tmp_name}.bed', sep='\t', index=False, header=False)

    # Calculate clusters
    if init_bedtools is not None:
        init_bedtools += ';'
    else:
        init_bedtools = ''

    subprocess.run(f'{init_bedtools} sort -k1,1 -k2,2n {tmp_name}.bed > {tmp_name}.sorted.bed; bedtools cluster -i {tmp_name}.sorted.bed -d {min_gap} > {tmp_name}.cluster.bed',
                    shell=True,
                    check=True,
                    text=True)

    # Convert clustered bed back to gtf
    colnames.append('cluster')
    gtf = pd.read_csv(f'{tmp_name}.cluster.bed', sep='\t',
                       names=colnames, dtype={'seqname': 'str'})
    gtf['start'] = [s+1 for s in gtf['start'].tolist()]
    # Restore original column order
    org_colnames.append('cluster')
    gtf = gtf[org_colnames]
    out_name = os.path.join(out_dir, name)
    out_path = os.path.realpath(f'{out_name}.d{min_gap}_clustered.gtf')
    gtf.to_csv(out_path, sep='\t', index=False, header=False)

    # Delete intermediates
    os.remove(f'{tmp_name}.bed')
    os.remove(f'{tmp_name}.sorted.bed')
    os.remove(f'{tmp_name}.cluster.bed')

    return out_path


