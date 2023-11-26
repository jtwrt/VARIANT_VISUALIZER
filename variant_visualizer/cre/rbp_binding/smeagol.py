import smeagol.io, smeagol.models, smeagol.scan, smeagol.enrich, smeagol.visualize, smeagol.variant
from ..._config import config
import os
from ... import core
from ._rbp_binding import RBPBinding

def load_smeagol_pwms(data_dir=config['smeagol']['motif_dir']):
    wd = os.getcwd()
    os.chdir(data_dir)
    #pwms = smeagol.io.load_rbpdb() # this only load rbpdb data
    os.makedirs('motifs/', exist_ok=True)
    pwms = smeagol.io.load_smeagol_PWMset()
    os.chdir(wd)
    return pwms

def load_smeagol_model(smeagol_pwms, data_dir=config['smeagol']['motif_dir']):
    """Create smeagol model for the given pwms."""
    pwms = load_smeagol_pwms()
    return smeagol.models.PWMModel(pwms)

def load_smeagol_genome(genome_path=config['smeagol']['genome_fasta']) -> list:
    return smeagol.io.read_fasta(genome_path)

def get_region_sequence(bio_region: core.BioRegion, genome, as_string=False):
    """
    Returns the nucleotide sequence of the provided region as BioPython Seq object.
    Returns the sequence as String if as_string is True.
    """
    chr = bio_region.reference.chromosome
    
    # select correct sequence/chromosome from genome
    for seq in genome:
        if seq.description[0:len(chr)] == chr:
            chr_seq = seq
            break

    # select bio_region nucleotide sequence
    bio_region_seq = chr_seq[bio_region.start-1:bio_region.end]

    # convert nucleotide sequence to upper case letters only
    bio_region_seq.seq = bio_region_seq.seq.upper()

    if as_string is True:
        if bio_region.reference.strand == '+':
            sequence = str(bio_region_seq[:].seq)
        elif bio_region.reference.strand == '-':
            sequence = str(bio_region_seq.reverse_complement()[:].seq)
        else:
            raise ValueError
        return sequence
    else: 
        return bio_region_seq

def get_smeagol_rbp_binding(bio_region: core.BioRegion, pwms, model: smeagol.models.PWMModel, genome: list, bio_region_seq=None, threshold=config['smeagol']['threshold']) -> list:
    """
    Description
    ---
    Scans the genomic sequence of a given genomic BioRegion for RBP binding sites. Returns them in a list.
    If bio_region_seq is given, bypasses seq generation and uses the provided object instead.
    Generate bio_region_seq with get_region_sequence. 
    The bio_region.reference is still used to determine which strands
     sequence will be used (reverse complement for - strand).
    """
    # select bio_region nucleotide sequence
    bio_region_smeagol_seq = get_region_sequence(
                                bio_region=bio_region,
                                genome=genome)

    # allow overwrite of seq with provided seq
    if bio_region_seq is not None:
        bio_region_smeagol_seq.seq = bio_region_seq

    # convert nucleotide sequence to upper case letters only
    bio_region_smeagol_seq.seq = bio_region_smeagol_seq.seq.upper()

    # calculate binding sites
    if bio_region.reference.strand == '+':
        rcomp='none'
    elif bio_region.reference.strand == '-':
        rcomp='only'
    else:
        raise ValueError('Unsupported reference.strand value in bio_region.')
    sites = smeagol.scan.scan_sequences([bio_region_smeagol_seq], 
                                        model=model, 
                                        threshold=threshold, 
                                        rcomp=rcomp, 
                                        outputs=['sites'], 
                                        score=True)
    # convert sites to BioRegions
    binding_sites =[]
    for _,row in sites['sites'].iterrows():
        if row['sense'] == '+':
            sequence = str(bio_region_smeagol_seq[row['start']:row['end']].seq)
            start = row['start']+bio_region.start
            end = row['end']+bio_region.start-1
        elif row['sense'] == '-':
            sequence = str(bio_region_smeagol_seq.reverse_complement().seq[row['start']:row['end']])
            start = bio_region.end-row['end']+1# smeagol site locations on rcomp strand count from beginning of rcomp sequence, not of provided sequence! 
            end = bio_region.end-row['start']
        info = dict(
            raw_score=row['score'],
            max_score=row['max_score'],
            width=row['width'],
            sequence=sequence,
            matrix_id=str(row['Matrix_id'])        
            )
        binding_sites.append(
            RBPBinding(
                start=start,
                end=end,
                reference=core.get_reference(
                    'genomic',
                    chromosome=bio_region.reference.chromosome,
                    strand=row['sense']),
                binding=pwms[pwms.Matrix_id==str(row['Matrix_id'])]['Gene_name'].values[0],
                score=row['frac_score'],
                info=info,
                source='smeagol_prediction'
                ))
    return binding_sites