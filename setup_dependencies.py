import os, sys
from pathlib import Path
from urllib.request import urlretrieve
from variant_visualizer._config import config
import variant_visualizer.setup._setup_utils as s_utils
import variant_visualizer.setup._setup_tcga_mc3 as s_mc3
import variant_visualizer.setup._setup_gtf as s_gtf
from variant_visualizer.cre.rbp_binding import load_smeagol_pwms
from Bio import SeqIO

def main():

    print('VARIANT_VISUALIZER: Initiating setup ...')

    # Check required values in config file
    print('Preparing config.yml ...')

    ## define tmp_dir
    if config['general']['tmp_dir'] is None:
        config['general']['tmp_dir'] = os.path.realpath('.tmp/')
        s_utils.update_config(config)

    ## define data_dir
    if config['general']['data_dir'] is None:
        config['general']['data_dir'] = os.path.realpath('.data/')
        s_utils.update_config(config)


    # Download missing dependencies and update config file accordingly

    if config['general']['gtf'] is None:
        print('Required file not provided: Downloading gtf ...')
        config['general']['gtf_header'] = 4
        s_utils.update_config(config)
        url = 'https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz'
        filename = os.path.join(config['general']['data_dir'],
                                'Homo_sapiens.GRCh38.109.gtf.gz')
        Path(filename).touch()
        _ = urlretrieve(url=url,
                                         filename=filename)
        filename = s_utils.decompress_gz(filename)
        config['general']['gtf'] = filename
        print('Done')
    
    if config['general']['gtf_header'] is None:
        raise AttributeError('gtf_header value must be provided in config.yml')

    if config['pas_atlas']['pas_atlas_path'] is None:
        print('Required file not provided: Downloading PAS-Atlas ...')
        url = 'https://www.polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.bed.gz'
        filename = os.path.join(config['general']['data_dir'],
                                'atlas.clusters.2.0.GRCh38.96.bed.gz')
        Path(filename).touch()
        _ = urlretrieve(url=url, filename=filename)
        filename = s_utils.decompress_gz(filename)
        config['pas_atlas']['pas_atlas_path'] = filename
        s_utils.update_config(config)
        print('Done')

    if config['tcga_mc3']['mc3_GRCh37_maf'] is None:
        print('Required file not provided: Downloading TCGA MC3 ...')
        url = 'https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc'
        filename = os.path.join(config['general']['data_dir'],
                                'mc3.v0.2.8.PUBLIC.maf.gz')
        Path(filename).touch()
        _ = urlretrieve(url=url,
                                        filename=filename)
        filename = s_utils.decompress_gz(filename)
        config['tcga_mc3']['mc3_GRCh37_maf'] = filename
        s_utils.update_config(config)
        print('Done')

    if config['tcga_mc3']['ucsc_liftover_chainfile'] is None:
        print('Required file not provided: Downloading UCSC Liftover Chainfile for hg19 to hg38 ...')
        url = 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz'
        filename = os.path.join(config['general']['data_dir'],
                                'hg19ToHg38.over.chain.gz')
        Path(filename).touch()
        _ = urlretrieve(url=url,
                                        filename=filename)
        filename = s_utils.decompress_gz(filename)
        config['tcga_mc3']['ucsc_liftover_chainfile'] = filename
        s_utils.update_config(config)
        print('Done')

    if config['tcga_mc3']['tcga_code_tables'] is None:
        print('Required file not provided: Downloading TCGA code tables ...')
        url = 'https://gdc.cancer.gov/files/public/file/tcga_code_tables.zip'
        filename = os.path.join(config['general']['data_dir'],
                                'tcga_code_tables.zip')
        Path(filename).touch()
        _ = urlretrieve(url=url,
                                         filename=filename)
        dir = s_utils.decompress_zip(filename)
        config['tcga_mc3']['tcga_code_tables'] = dir
        s_utils.update_config(config)
        print('Done')

    if config['tcga_mc3']['tcga_sbs_path'] is None:
        print('Required file not provided: Downloading TCGA sbs signatures ...')
        url = 'https://dcc.icgc.org/api/v1/download?fn=/PCAWG/mutational_signatures/Signatures_in_Samples/SP_Signatures_in_Samples/TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv'
        filename = os.path.join(config['general']['data_dir'],
                                'TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv')
        Path(filename).touch()
        _ = urlretrieve(url=url,
                                         filename=filename)
        config['tcga_mc3']['tcga_sbs_path'] = filename
        s_utils.update_config(config)
        print('Done')

    if config['tcga_mc3']['tcga_ids_path'] is None:
        print('Required file not provided: Downloading TCGA id signatures ...')
        url = 'https://dcc.icgc.org/api/v1/download?fn=/PCAWG/mutational_signatures/Signatures_in_Samples/SP_Signatures_in_Samples/TCGA_WES_sigProfiler_ID_signatures_in_samples.csv'
        filename = os.path.join(config['general']['data_dir'],
                                'TCGA_WES_sigProfiler_ID_signatures_in_samples.csv')
        Path(filename).touch()
        _ = urlretrieve(url=url,
                                         filename=filename)
        config['tcga_mc3']['tcga_ids_path'] = filename
        s_utils.update_config(config)
        print('Done')

    if config['targetscan']['target_locations_hg19_bed'] is None:
        print('Required file not provided: Downloading miRNA binding site predictions from Targetscan database ...')
        url = 'https://www.targetscan.org/vert_80/vert_80_data_download/Predicted_Target_Locations.default_predictions.hg19.bed.zip'
        filename = os.path.join(config['general']['data_dir'],
                                'Predicted_Target_Locations.default_predictions.hg19.bed.zip')
        Path(filename).touch()
        _ = urlretrieve(url=url,
                                         filename=filename)
        dir = s_utils.decompress_zip(filename, 'Targetscan_Predicted_Target_Locations')
        config['targetscan']['target_locations_hg19_bed'] = os.path.join(dir, 'Predicted_Target_Locations.default_predictions.hg19.bed')
        s_utils.update_config(config)
        print('Done')

    if config['smeagol']['genome_fasta'] is None:
        print('Required file not provided: Downloading Human genome from Ensembl ...')
        url = 'https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz'
        filename = os.path.join(config['general']['data_dir'],
                                'Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz')
        Path(filename).touch()
        _ = urlretrieve(url=url,
                                         filename=filename)
        filename = s_utils.decompress_gz(filename)
        config['smeagol']['genome_fasta'] = filename
        s_utils.update_config(config)
        print('Done')

    if config['smeagol']['motif_dir'] is None:
        print('Required files not provided: Downloading PWMS for smeagol module ...')
        out_dir = os.path.join(config['general']['data_dir'], 'smeagol')
        os.makedirs(out_dir, exist_ok=True)
        _ = load_smeagol_pwms(data_dir=out_dir)
        config['smeagol']['motif_dir'] = out_dir
        s_utils.update_config(config)
        print('Done')

    if config['uniprotkb']['sprot_xml'] is None:
        print('Required file not provided: Downloading UniprotKB SwissProt human database ...')
        url = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_human.xml.gz'
        filename = os.path.join(config['general']['data_dir'],
                                'uniprot_sprot_human.xml.gz')
        Path(filename).touch()
        _ = urlretrieve(url=url,
                                         filename=filename)
        filename = s_utils.decompress_gz(filename)
        config['uniprotkb']['sprot_xml'] = filename
        s_utils.update_config(config)
        print('Done')

    if config['uniprotkb']['ensembl_sprot_map'] is None:
        print('Mapping ensembl transcript ids to uniprot entries ...')
        uniprotkb_record_dict = SeqIO.index(config['uniprotkb']['sprot_xml'], "uniprot-xml")
        # map to ensembl ids
        ensembl_map = dict()
        for k, v in uniprotkb_record_dict.items():
            ids = [ref.split(':')[1] for ref in v.dbxrefs if ref.split(':')[0]=='Ensembl']
            for id in ids:
                id = id.split('.')[0]
                if id not in ensembl_map:
                    ensembl_map[id] = v
                elif id in ensembl_map:
                    raise Exception(f'Multiple values with ensemble ensembl transcript id \'{k}\'')
        out_path = os.path.join(config['general']['data_dir'],
                                'ensembl_sprot_map.dill')
        s_utils.dill_dump_object(out_path=out_path,
                         object=ensembl_map)
        config['uniprotkb']['ensembl_sprot_map'] = out_path
        s_utils.update_config(config)
        print('Done')   
 
    ## Generate clustered gtf
    if config['general']['clustered_gtf'] is None:
        print('Generating clustered gtf ...')
        clustered_gtf_path = s_gtf.cluster_gtf()
        config['general']['clustered_gtf'] = clustered_gtf_path
        s_utils.update_config(config)
        print('Done')

    ## Liftover mc3 from GRCh37 to GRCh38
    if config['tcga_mc3']['mc3_GRCh38_maf'] is None:
        print('Initiating MC3 liftover to GRCh38 ...')
        mc3_GRCh38_path = s_mc3.liftover_mc3()
        config['tcga_mc3']['mc3_GRCh38_maf'] = mc3_GRCh38_path
        s_utils.update_config(config)
        print('Done')

    ## Liftover targetscan db to GRCh38
    if config['targetscan']['target_locations_GRCh38_bed'] is None:
        print('Initiating targetscan DB liftover to GRCh38 ...')
        input_bed_path = str(config['targetscan']['target_locations_hg19_bed'])
        output_bed_path = input_bed_path.replace('.hg19.', '.GRCh38.')
        s_utils.liftover_GRCh37_GRCh38(input_bed_path=input_bed_path,
                               output_bed_path=output_bed_path)
        config['targetscan']['target_locations_GRCh38_bed'] = output_bed_path
        s_utils.update_config(config)
        print('Done')

    print('VARIANT_VISUALIZER dependencies are complete!')

if __name__ == '__main__':
    main()
