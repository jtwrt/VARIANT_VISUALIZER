from copy import deepcopy as _deepcopy
from Bio import SeqIO
from ._protein_annotation import ProteinAnnotation
from ..setup import _setup_utils as s_utils
from .._config import config
from .. import core, clusters


class UniprotAnnotations():

    @classmethod
    def _load_sprot_ensembl_map(cls):
        return s_utils.dill_load_object(config['uniprotkb']['ensembl_sprot_map'])

    def __init__(self):
        self.ensembl_map = self.__class__._load_sprot_ensembl_map()
        # map from gene_names to assicated transcript_ids
        gene_name_map = dict()
        for transcript_id in self.ensembl_map.keys():
            gene_name = self.get_transcript_gene_name(transcript_id)
            if gene_name is None:
                continue
            if not gene_name_map.get(gene_name):
                gene_name_map[gene_name] = set()
            gene_name_map[gene_name].add(transcript_id)
        self.gene_name_map = gene_name_map
    
    def get_transcript_features(self, ensembl_id: str, cluster: clusters.Cluster, features:list=[]) -> list:
        """
        Description
        ---
        Return annotated regions as of the queried transcript from uniprotkb as BioRegions.
        If features is empty, returns all annotations. 
        Otherwise only returns annotations where type is in features.

        Valid feature types are:
            'DNA-binding region',
            'active site',
            'binding site',
            'chain',
            'coiled-coil region',
            'compositionally biased region',
            'cross-link',
            'disulfide bond',
            'domain',
            'glycosylation site',
            'helix',
            'initiator methionine',
            'intramembrane region',
            'lipid moiety-binding region',
            'modified residue',
            'mutagenesis site',
            'non-standard amino acid',
            'non-terminal residue',
            'peptide',
            'propeptide',
            'region of interest',
            'repeat',
            'sequence conflict',
            'sequence variant',
            'short sequence motif',
            'signal peptide',
            'site',
            'splice variant',
            'strand',
            'topological domain',
            'transit peptide',
            'transmembrane region',
            'turn',
            'zinc finger region'
        """
        gtf_cluster = cluster.gtf_cluster
        annotations = self.ensembl_map.get(ensembl_id).features
        if len(features) != 0:
            annotations = [a for a in annotations if a.type in features]
        
        coding_regions = set()
        transcript_region = None
        for r in gtf_cluster.all_regions:
            if r.transcript_id == ensembl_id and r.feature == 'CDS':
                coding_regions.add(r)
            if r.transcript_id == ensembl_id and r.feature == 'transcript':
                transcript_region = r
        if transcript_region is None:
            raise ValueError('Transcript is not in given gtf_cluster') 

        out = []
        for a in annotations:
            if a.qualifiers.get('description'):
                description = a.qualifiers['description']
            else:
                description = ''
            out.append(ProteinAnnotation(
                start=int(a.location._start)+1,
                end=int(a.location._end),
                reference=core.get_reference('protein', 
                                        transcript_region=transcript_region,
                                        coding_regions=coding_regions,
                                        chromosome=transcript_region.reference.chromosome,
                                        strand=transcript_region.reference.strand),
                annotation_type=a.type,
                description=description,
                source='UniprotKB_annotation',
                label=description
                ))
        return out

    def get_transcript_gene_name(self, ensembl_transcript_id) -> str:
        """
        Description
        ---
        Returns the annotated primary gene name if it exists for the given tanscript.
        Returns None if the ensembl_id is not part of self.ensembl_map or if no gene name is annotated.
        """
        if not self.ensembl_map.get(ensembl_transcript_id):
            return None
        if not self.ensembl_map[ensembl_transcript_id].annotations.get('gene_name_primary'):
            return None
        gene_name = self.ensembl_map[ensembl_transcript_id].annotations['gene_name_primary']
        if gene_name == '':
            return None
        else:
            return str(gene_name)

    def get_gene_name(self, ensembl_gene_id: str, cluster: clusters.Cluster):
        """
        Description
        ---
        Returns the annotated gene name. 
        If multiple gene_names are defined as primary in different transcripts, 
        returns a string containing all of them delimited by \'; \'.
        Returns the original gene_id if no gene name can be retrieved.
        """
        gtf_cluster = cluster.gtf_cluster
        gene_regions = [r for r in gtf_cluster.all_regions if r.gene_id == ensembl_gene_id]
        transcript_ids = set([r.transcript_id for r in gene_regions])
        
        out = set()
        for t in transcript_ids:
            gene_name = self.get_transcript_gene_name(t)
            if gene_name is not None:
                out.add(gene_name)
        if out != set():
            return '; '.join([s for s in sorted(out)])
        else:
            return ensembl_gene_id
        
    def get_ensembl_ids(self, gene_name: str) -> set:
        """
        Returns a set of associated ensembl transcript ids with the given gene.
        Raises KeyError for unkown gene_name (Not listed as uniprotkb primary gene name).
        """
        if gene_name not in self.gene_name_map.keys():
            raise KeyError('Unkown gene name provided.')
        return self.gene_name_map[gene_name]
    
    def get_transcript_gene_id(self, ensembl_transcript_id: str, cluster: clusters.Cluster):
        """
        Returns the gene_id associated with the transcript_id. 
        Raises ValueError if the no gene_id is associated with the given 
        transcript_id in the given Cluster.
        If a list/set of transcript_ids is given, searches associated gene_ids for each.
        The found gene_id must match for each.
        """
        if isinstance(ensembl_transcript_id, list):
            ensembl_transcript_id = set(ensembl_transcript_id)

        if isinstance(ensembl_transcript_id, set):
            out = set()
            [out.add(self.get_transcript_gene_id(id, cluster)) for id in ensembl_transcript_id]
        elif isinstance(ensembl_transcript_id, str):
            gtf_cluster = cluster.gtf_cluster
            transcript_regions = [r for r in gtf_cluster.all_regions if r.transcript_id == ensembl_transcript_id]
            out = set()
            [out.add(r.gene_id) for r in transcript_regions]
        else:
            raise TypeError('Provide ensembl transcript_id or a list thereof.')
        
        if len(out) == 0:
            raise ValueError('At least one transcript_id not associated with any gene_id in this cluster.')
        elif len(out) == 1:
            return list(out)[0]
        else: 
            raise IndexError('Mutliple gene_ids associated with the given transcript_ids.')
