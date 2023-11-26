from __future__ import annotations
from variant_visualizer.core._bio_references import _BioReference
from .._config import config
from .. import core
import pandas as pd
import os
from copy import deepcopy
from ..setup._setup_gtf import load_gtf

def get_gtf_path(input_file, out_dir, file_ending=''):
    """Return a path where the file format from the input file is exchanged to file_ending, and the directory is out_dir."""
    name = os.path.splitext(os.path.basename(input_file))[0]
    outName = os.path.join(out_dir, name)
    outName = os.path.realpath(outName)
    return str(outName + file_ending)

def get_gtf_cluster_slice(cluster: int, clustered_gtf: pd.DataFrame):
    """Return the cluster features from the provided gtf."""
    return clustered_gtf.loc[(clustered_gtf['cluster'] == cluster)]

def get_gtf_slice_path(cluster_id: int):
    path = get_gtf_path(
        input_file=config['general']['clustered_gtf'],
        out_dir=os.path.join(
            config['general']['tmp_dir'],
            'cluster_inputs'),
        file_ending=f'.{cluster_id}_cluster.gtf')
    return path

def slice_clustered_gtf(cluster: int, clustered_gtf: pd.DataFrame) -> str:
    """
    Description
    ---
    Split the clustered gtf dataframe into separate dataframes, and save it.
    clustered_gtf must be loaded to memory core.io.load_gtf(clustered=True).
    To load the dataframe created to memory correctly, use load_gtf_cluster(path).
    """
    path = get_gtf_slice_path(cluster)
    get_gtf_cluster_slice(cluster, clustered_gtf).to_csv(path, sep='\t')
    return os.path.realpath(path)

def load_gtf_clustered_slice(cluster_id) -> pd.DataFrame:
    """Load a gtf slice generated with generateInput.sliceGtfCluster()"""
    path = get_gtf_slice_path(cluster_id)
        
    gtf = pd.read_csv(
        path,
        sep='\t',
        index_col=0,
        dtype={'seqname': 'str'})
    return gtf

class GtfFeature(core.BioRegion):
    """
    BioRegion including information of a gtf feature.
    """    
    def __init__(
        self, start: int, end: int, reference: core._BioReference, 
        feature: str, transcript_id: str, gene_id: str,
        transcript_biotype: str, gene_biotype: str):
        super().__init__(start, end, reference)
        gtf_attributes = [
            feature,
            transcript_id,
            gene_id,
            transcript_biotype,
            gene_biotype
        ]
        if all([isinstance(attribute, str) for attribute in gtf_attributes]):
            pass
        self.feature = feature
        self.transcript_id = transcript_id
        self.gene_id = gene_id
        self.transcript_biotype = transcript_biotype
        self.gene_biotype = gene_biotype

    def __key(self):
        return ('GtfFeature', self.start, self.end, self.reference, self.get_gtf_info())
    
    def __hash__(self) -> tuple:
        return hash(self.__key())

    def __repr__(self):
        return f'GtfFeature {self.start}-{self.end}, reference: {self.reference}, {self.get_gtf_info()}'
    
    def __eq__(self, other: GtfFeature):
            """Compares BioRegions. Must be implemented by subclasses."""
            this_class = GtfFeature
            if isinstance(other, this_class):
                return self.__key() == other.__key()
            elif isinstance(other, core.BioRegion) or other is None:
                return False
            else:
                raise TypeError(f'Testing equality not defined between given objects.')
        
    def get_gtf_info(self) -> str:
        """
        Description
        ---
        Returns information feature, id and biotype information of the GtfFeature in a tuple
        containing: 'feature','transcript_id','gene_id','transcript_biotype','gene_biotype'
        """
        out = tuple()
        for key in ['feature','transcript_id','gene_id','transcript_biotype','gene_biotype']:
            out += tuple([str(self.__dict__[key])])
        return out

class EmptyClusterException(Exception):
    pass

class GtfCluster(core.BioRegion):
    """
    BioRegion subclass. Contains information of gtf features that are clustered in the gtf file.
    Attributes:
    all_regions: list containing all GtfFeatures of the cluster
    cluster_segments: collapsed features for each strand
    transcript_ids: set of ensembl transcript ids for transcripts in the cluster
    gene_ids: set of ensembl gene ids for genes in the cluster.
    Raises EmptyClusterException if there are no features in the cluster.
    """
    def __init__(self, cluster_id: int, gtf_cluster_slice: pd.DataFrame) -> None:
        
        self.cluster_id = cluster_id

        # Load gtf slice of cluster
        cluster_slice = load_gtf_clustered_slice(cluster_id)
        self.all_regions = self._get_regions(gtf_slice=cluster_slice)

        self.cluster_segments = collapse_gtf_features(self.all_regions)

        if len(self.all_regions) == 0:
            raise EmptyClusterException('No features in cluster.')

        self.transcript_ids = set([r.transcript_id for r in self.all_regions])
        try:
            self.transcript_ids.remove('')
        except KeyError: pass
        self.gene_ids = set([r.gene_id for r in self.all_regions])
        try:
            self.gene_ids.remove('')
        except KeyError: pass


        super().__init__(start=min([r.start for r in self.all_regions]),
                         end=max([r.end for r in self.all_regions]),
                         reference=core.get_reference(reference_type='genomic', 
                                                  chromosome=self.all_regions[0].reference.chromosome
                                                  )
                         )
        
        self.strand_regions = {'+':[],
                               '-':[]}
        [self.strand_regions[r.reference.strand].append(r) for r in self.all_regions]

        # index features
        self.nt_regions = self._get_nt_regions()

    def __key(self):
        return ('GtfCluster', self.start, self.end, self.reference)
    
    def __hash__(self) -> tuple:
        return hash(self.__key())

    def __repr__(self):
        return f'GtfCluster {self.start}-{self.end}, cluster_id: {self.cluster_id}, reference: {self.reference}'
    
    def __eq__(self, other: GtfCluster):
            """Compares BioRegions. Must be implemented by subclasses."""
            this_class = GtfCluster
            # do not compare subclasses with this method
            if issubclass(self.__class__, this_class) or issubclass(other.__class__, this_class):
                raise NotImplementedError("Subclasses must implement __eq__ method to test equality.")
            if isinstance(other, this_class):
                return self.__key() == other.__key()
            elif other is None:
                return False
            else:
                raise TypeError(f'Testing equality not defined between given objects.')

    def _get_regions(self, gtf_slice: pd.DataFrame) -> list:
        """
        Parse rows from the gtf and return them as a list of GtfFeature objects.
        """

        chromosome = None
        out = []
        for i, row in gtf_slice.iterrows():
            if chromosome is None: chromosome = str(row['seqname'])
            reference = core.get_reference(reference_type='genomic',
                                    chromosome=chromosome,
                                    strand=row['strand'])
            
            info = {'feature':'',
                    'gene_id':'',
                    'transcript_id':'',
                    'gene_biotype':'',
                    'transcript_biotype':''}
            for key in info:
                if not pd.isna(row[key]):
                    info[key] = row[key]

            out.append(GtfFeature(start=row['start'],
                                end=row['end'],
                                reference=reference,
                                feature=info['feature'],
                                gene_id=info['gene_id'],
                                transcript_id=info['transcript_id'],
                                gene_biotype=info['gene_biotype'],
                                transcript_biotype=info['transcript_biotype']
                    ))
        
        # Update reference objects to reference transcript/coding regions
        #   One reference object per transcript_id to reduce memory usage

        reference_regions = {}
        transcripts = []
        empty_ref = {'transcript_region':None,
                    'coding_regions':[]
                    }
        for r in out:
            if r.transcript_id == '':
                continue
            if not reference_regions.get(r.transcript_id):
                reference_regions[r.transcript_id] = deepcopy(empty_ref)
            if r.feature == 'transcript':
                reference_regions[r.transcript_id]['transcript_region'] = r
                transcripts.append(r)
            elif r.feature == 'CDS':
                reference_regions[r.transcript_id]['coding_regions'].append(r)

        references = dict()
        for id in reference_regions:
            try:
                references[id] = core.get_reference(reference_type='genomic',
                                        chromosome=reference_regions[id]['transcript_region'].reference.chromosome,
                                        strand=reference_regions[id]['transcript_region'].reference.strand,
                                        transcript_region=core.Region(start=reference_regions[id]['transcript_region'].start,
                                                                end=reference_regions[id]['transcript_region'].end
                                                                ),
                                        coding_regions=[core.Region(c.start,c.end) for c in reference_regions[id]['coding_regions']]
                                        )
            except AttributeError:
                pass

        for r in out:
            if r.transcript_id != '' and references.get(r.transcript_id):
                r.update(reference=references[r.transcript_id])

        # Create regulatory sequence objects according to config.yml specifications
        for r in transcripts:
            if r.reference.strand == '+':
                start = r.start - config['general']['five_prime_regulatory_nts']
                end = r.end + config['general']['three_prime_regulatory_nts']
                out.append(GtfFeature(start=start,
                                    end=r.start-1,
                                    reference=r.reference,
                                    feature='five_prime_regulatory_sequence',
                                    gene_id=r.gene_id,
                                    transcript_id=r.transcript_id,
                                    gene_biotype=r.gene_biotype,
                                    transcript_biotype=r.transcript_biotype
                                    ))
                out.append(GtfFeature(start=r.end+1,
                                    end=end,
                                    reference=r.reference,
                                    feature='three_prime_regulatory_sequence',
                                    gene_id=r.gene_id,
                                    transcript_id=r.transcript_id,
                                    gene_biotype=r.gene_biotype,
                                    transcript_biotype=r.transcript_biotype
                                    ))
            elif r.reference.strand == '-':
                start = r.start - config['general']['three_prime_regulatory_nts']
                end = r.end + config['general']['five_prime_regulatory_nts']
                out.append(GtfFeature(start=start,
                                    end=r.start-1,
                                    reference=r.reference,
                                    feature='three_prime_regulatory_sequence',
                                    gene_id=r.gene_id,
                                    transcript_id=r.transcript_id,
                                    gene_biotype=r.gene_biotype,
                                    transcript_biotype=r.transcript_biotype
                                    ))
                out.append(GtfFeature(start=r.end+1,
                                    end=end,
                                    reference=r.reference,
                                    feature='five_prime_regulatory_sequence',
                                    gene_id=r.gene_id,
                                    transcript_id=r.transcript_id,
                                    gene_biotype=r.gene_biotype,
                                    transcript_biotype=r.transcript_biotype
                                    ))
            else:
                raise ValueError(f'Unsupported strand value: \'{r.reference.strand}\'')
        
        return out

    def _get_nt_regions(self) -> dict:

        
        nt_regions = {i: [] for i in range(self.start,self.end+1)}
        
        for r in self.all_regions:
            [nt_regions[i].append(r) for i in range(r.start,r.end+1)]

        return nt_regions
    
    def get_nt_info(self, nucleotide: int) -> set:
        if nucleotide not in self.nt_regions.keys():
            return set()
        out = set()    
        for r in self.nt_regions[nucleotide]:
            out.add(r.get_gtf_info())
        return out

    def get_overlap_gtf_features(self, region: core.BioRegion) -> set():
        """Return a set of gtf feature types which the given region is overlapping"""
        nt_features = set()
        for nt in range(region.start, region.end+1):
            [nt_features.add(gtf_info[0]) for gtf_info in self.get_nt_info(nt)]
        return nt_features

    def get_similar_nts(self, nucleotide: int):
        """Return list of nucleotides with identical gtf features as the one provided"""
        if nucleotide not in self.nt_regions.keys():
            return []
        out = []
        for i in self.nt_regions:
            if self.nt_regions[nucleotide] == self.nt_regions[i]:
                out.append(i)
        return out

    def _get_location_type(self, location: core.BioRegion) -> dict:
        """
        Description
        ---
        Checks which features are on either strand of the given location.
        Given location must be BioRegion with length 1.
        """
        if not isinstance(location.reference, core.GenomicReference):
            raise ValueError('location.reference must be of type GenomicReference.')
 
        out = {'+': 'none',
               '-': 'none'}
        for strand in self.strand_segments:
            break2 = False
            for feature_type in self.strand_segments[strand]:
                if break2 is True:
                    break
                for feature_region in self.strand_segments[strand][feature_type]:
                    if location.within(feature_region):
                        out[strand] = feature_type
                        break2 = True
                        break
        return out

def collapse_gtf_features(regions):
    """
    Description
    ---
    Combines BioRegions describing gtf features. CDS, start_codon and stop_codon are combined it they overlap.
    Regions that do not overlap with them are combined to five_prime_utr or three_prime_utr or five_and_three_prime_utr features.
    Regions within a transcript feature, which are not included in the previously generated features are defined as intron features
    if the transcript_biotype is protein_coding. Otherwise, they are defined as 'other_feature'.
    Additionaly, five_prime_regulatory_sequence, three_prime_regulatory_sequence and five_and_three_prime_regulatory_sequence
    features are generated 
    Features are only collapsed on identical strands.
    The output contains a dictionary with subdictionary for each strand, each of which contains a list of the previously generated features.
    """
    
    strand_regions = {'+': [],
                      '-': []}
    for r in regions:
        strand_regions[r.reference.strand].append(r)
    
    out = {'+': dict(),
           '-': dict()}
    for strand in strand_regions:
        out[strand]['CDS'] = core.combine_regions([r for r in strand_regions[strand] if r.feature=='CDS'])
        remove_regions = out[strand]['CDS']

        out[strand]['start_codon'] = core.combine_regions([r for r in strand_regions[strand] if r.feature=='start_codon'])
        remove_regions = core.combine_regions(remove_regions+
                                            out[strand]['start_codon'])

        out[strand]['stop_codon'] = core.combine_regions([r for r in strand_regions[strand] if r.feature=='stop_codon'])
        remove_regions = core.combine_regions(remove_regions+
                                            out[strand]['stop_codon'])

        out[strand]['five_prime_utr'] = core.combine_regions([r for r in strand_regions[strand] if r.feature=='five_prime_utr'])
        out[strand]['five_prime_utr'] = core.remove_regions(regions=out[strand]['five_prime_utr'],
                                                        remove_regions=remove_regions)
        out[strand]['three_prime_utr'] = core.combine_regions([r for r in strand_regions[strand] if r.feature=='three_prime_utr'])
        out[strand]['three_prime_utr'] = core.remove_regions(regions=out[strand]['three_prime_utr'],
                                                            remove_regions=remove_regions)
        out[strand]['five_and_three_prime_utr'] = core.get_all_overlaps(out[strand]['five_prime_utr'],
                                                                    out[strand]['three_prime_utr'])
        out[strand]['five_prime_utr'] = core.remove_regions(regions=out[strand]['five_prime_utr'],
                                                        remove_regions=out[strand]['five_and_three_prime_utr'])
        out[strand]['three_prime_utr'] = core.remove_regions(regions=out[strand]['three_prime_utr'],
                                                        remove_regions=out[strand]['five_and_three_prime_utr'])

        remove_regions = core.combine_regions(remove_regions+
                                            out[strand]['five_prime_utr']+
                                            out[strand]['three_prime_utr']+
                                            out[strand]['five_and_three_prime_utr']
                                            )
        
        out[strand]['intron'] = core.combine_regions([r for r in strand_regions[strand] if r.feature=='transcript' and r.transcript_biotype=='protein_coding'])
        out[strand]['intron'] = core.remove_regions(regions=out[strand]['intron'],
                                                remove_regions=remove_regions)
        remove_regions = core.combine_regions(remove_regions+
                                            out[strand]['intron'])

                                            

        out[strand]['five_prime_regulatory_sequence'] = core.combine_regions([r for r in strand_regions[strand] if r.feature=='five_prime_regulatory_sequence'])
        out[strand]['five_prime_regulatory_sequence'] = core.remove_regions(regions=out[strand]['five_prime_regulatory_sequence'],
                                                                        remove_regions=remove_regions)
        out[strand]['three_prime_regulatory_sequence'] = core.combine_regions([r for r in strand_regions[strand] if r.feature=='three_prime_regulatory_sequence'])
        out[strand]['three_prime_regulatory_sequence'] = core.remove_regions(regions=out[strand]['three_prime_regulatory_sequence'],
                                                                            remove_regions=remove_regions)
        out[strand]['five_and_three_prime_regulatory_sequence'] = core.get_all_overlaps(out[strand]['five_prime_regulatory_sequence'],
                                                                                    out[strand]['three_prime_regulatory_sequence'])
        out[strand]['five_prime_regulatory_sequence'] = core.remove_regions(regions=out[strand]['five_prime_regulatory_sequence'],
                                                                        remove_regions=out[strand]['five_and_three_prime_regulatory_sequence'])
        out[strand]['three_prime_regulatory_sequence'] = core.remove_regions(regions=out[strand]['three_prime_regulatory_sequence'],
                                                                            remove_regions=out[strand]['five_and_three_prime_regulatory_sequence'])
        
        out[strand]['other_feature'] = core.combine_regions([r for r in strand_regions[strand] if r.feature=='transcript'])
        out[strand]['other_feature'] = core.remove_regions(regions=out[strand]['other_feature'],
                                                                remove_regions=remove_regions)
        
    for strand in out:
        for feature in out[strand]:
            out[strand][feature] = [
                    core.BioRegion(
                        r.start,
                        r.end,
                        core.get_reference(
                            'genomic', 
                            chromosome=r.reference.chromosome,
                            strand=r.reference.strand)) 
                for r in out[strand][feature]]
    return out


def collapse_gtf_genes(gtf_regions):
    """
    Description
    ---
    Collapses gtf features with identical gene_id using the gtf.collapse_gtf_features function.
    The output contains a dict with gene_ids for keys and subdits containing the collapsed features.
    """
    gene_dict = dict()
    gene_strands = dict()
    for r in gtf_regions:
        if not gene_dict.get(r.gene_id):
            gene_dict[r.gene_id] = []
            gene_strands[r.gene_id] = r.reference.strand
        gene_dict[r.gene_id].append(r)

    for gene_id in gene_dict:
        collapsed_gene = collapse_gtf_features(gene_dict[gene_id])
        # select correct strand output
        gene_dict[gene_id] = collapsed_gene[gene_strands[gene_id]]
    return gene_dict
