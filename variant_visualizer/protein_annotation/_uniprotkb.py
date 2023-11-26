from copy import deepcopy as _deepcopy
from Bio import SeqIO
from ._protein_annotation import Annotation
from ..setup import _setup_utils as s_utils
from .._config import config
from .. import core


class UniprotAnnotations():

    @classmethod
    def _load_sprot_ensembl_map(cls):
        return s_utils.dill_load_object(config['uniprotkb']['ensembl_sprot_map'])

    def __init__(self):
        self.ensembl_map = self.__class__._load_sprot_ensembl_map()
    
    def get_transcript_features(self, ensembl_id: str, features: list, gtf_cluster) -> list:
        """
        Description
        ---
        Return annotated regions as of the queried transcript from uniprotkb as BioRegions.
        If features is empty, returns all annotations. 
        Otherwise only returns annotations where type is in features.
        """

        annotations = self.ensembl_map.get(ensembl_id).features
        if len(features) != 0:
            annotations = [a for a in annotations if a.type in features]
        
        coding_regions = []
        transcript_region = None
        for r in gtf_cluster.all_regions:
            if r.transcript_id == ensembl_id and r.feature == 'CDS':
                coding_regions.append(r)
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
            out.append(Annotation(
                start=int(a.location._start)+1,
                end=int(a.location._end),
                reference=core.get_reference('protein', 
                                        transcript_region=transcript_region,
                                        coding_regions=coding_regions,
                                        chromosome=transcript_region.reference.chromosome,
                                        strand=transcript_region.reference.strand),
                annotation_type=a.type,
                description=description
                ))
        return out

    def get_transcript_gene_name(self, ensembl_transcript_id):
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
            return gene_name

    def get_gene_name(self, ensembl_gene_id, gtf_cluster):
        """
        Description
        ---
        Returns the annotated gene name. 
        If multiple gene_names are defined as primary in different transcripts, 
        returns a string containing all of them delimited by \'; \'.
        Returns the original gene_id if no gene name can be retrieved.
        """
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
