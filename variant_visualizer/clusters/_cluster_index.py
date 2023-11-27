from .. import genome_annotation as ga
from .. import core
from .._config import config
import os, warnings
from ..setup._setup_utils import dill_dump_object, dill_load_object
from ._clusters import load_cluster, _get_cluster_path
from .. import protein_annotation as pa
from tqdm import tqdm

class IndexingError(Exception):
    """Raised, when errors occur during cluster indexing."""
    pass

class ClusterIndexGenerator():

    def __init__(self):
        print('Loading GTF ...')
        clustered_gtf = ga.load_gtf(clustered=True)
        all_gtf_cluster_ids = set(clustered_gtf['cluster'])

        print('Finding non-indexed Clusters ...')
        try:
            index = load_index()
            indexed_clusters = index._index['cluster_id'].keys()
            missing_clusters_ids = [i for i in all_gtf_cluster_ids if i not in indexed_clusters]
        except FileNotFoundError:
            index = ClusterIndex(dict()) # empty index
            missing_clusters_ids = all_gtf_cluster_ids

        uniprotkb = pa.UniprotAnnotations()

        index_in = dict()
        print('Reading clusters, gathering index data ...')
        
        missing_clusters = []
        to_be_indexed_clusters = []
        for cluster_id in missing_clusters_ids:
            if not os.path.isfile(_get_cluster_path(cluster_id)):
                missing_clusters.append(cluster_id)
            else:
                to_be_indexed_clusters.append(cluster_id)

        for cluster_id in tqdm(to_be_indexed_clusters):
            cluster = load_cluster(cluster_id)
            index_in[cluster_id] = dict(
                region=core.BioRegion(
                    cluster.gtf_cluster.start, cluster.gtf_cluster.end,
                    cluster.gtf_cluster.reference),
                ensembl_id=cluster.gtf_cluster.transcript_ids|cluster.gtf_cluster.gene_ids
        )
            index_in[cluster_id]['gene_name'] = []
            for gene_id in cluster.gtf_cluster.gene_ids:
                gene_names = uniprotkb.get_gene_name(gene_id, cluster.gtf_cluster)
                if gene_names == gene_id:
                    continue
                gene_names = gene_names.split(';')
                index_in[cluster_id]['gene_name'].extend(gene_names)

        index.update(index_in)
        dill_dump_object(_get_index_path(), index)
        if len(missing_clusters) > 0:
            warnings.warn('ClusterIndex not complete! {len(missing_clusters)} Clusters not indexed.')

class ClusterIndex():
    """
    Class that indexes genes and transcript and the genomic clusters
    they are located in.
    Used to find clusters which genes/transcript of interest.
    Generated when executing setup_index.py execution.
    """

    def __init__(self, input_data) -> None:
        self.config = config
        self._index = {
            'cluster_id':dict(),
            'ensembl_id':dict(),
            'gene_name':dict(),
            'bio_region':{'chromosome':dict()}
        }
        for cluster_id in input_data:
            self._update_index(cluster_id, input_data[cluster_id])

    def update(self, input_data):
        for cluster_id in input_data:
            self._update_index(cluster_id, input_data[cluster_id])
        return

    def _update_index(self, cluster_id, input_data):
        if cluster_id in self._index['cluster_id']:
            raise IndexingError('Input data includes information on previously indexed cluster.')

        # index by cluster_id
        self._index['cluster_id'][cluster_id] = input_data

        # index by ensembl ids
        for ensembl_id in input_data['ensembl_id']:
            if self._index['ensembl_id'].get(ensembl_id):
                raise IndexingError('Duplicate ensembl_id cannot be indexed.')
            self._index['ensembl_id'][ensembl_id] = cluster_id
        
        # index by region
        region = input_data['region']
        if not self._index['bio_region']['chromosome'].get(region.reference.chromosome):
            self._index['bio_region']['chromosome'][region.reference.chromosome] = []
        self._index['bio_region']['chromosome'][region.reference.chromosome].append((region, cluster_id))
        
        # index by gene names
        for gene_name in input_data['gene_name']:
            if gene_name in self._index['gene_name'].keys() and \
                    cluster_id != self._index['gene_name'][gene_name]:
                value = self._index['gene_name'][gene_name]
                if isinstance(value, int):
                    value = set([value])
                value.add(cluster_id)
                warnings.warn(f'Gene name {gene_name} indexed in multiple clusters: {value}')
                #raise IndexError(f'Gene can not be indexed in cluster {cluster_id} and {self._index["gene_name"][gene_name]}.')
            else:
                self._index['gene_name'][gene_name] = cluster_id

    def query_genomic_region(self, bio_region: core.BioRegion):
        """
        Returns the cluster_id of the cluster, in which the bio_region is located.
        Returns None if bio_region is not within any cluster.
        """
        if bio_region.get_reference_type != 'genomic':
            return ValueError('Can only select cluster from genomic regions.')
        chromosome = bio_region.reference.chromosome
        for cluster, cluster_id in self._index['bio_region'][chromosome]:
            bio_region.within(cluster)
            return cluster_id
        return None

    def query_ensembl_id(self, ensembl_id: str):
        """
        Returns the cluster_id of the cluster, in which the gene or transcript
        with the given id is located. 
        Returns None if ensembl_id is not indexed.
        """
        if self._index['ensembl_id'].get(ensembl_id):
            return self._index['ensembl_id'][ensembl_id]
        else:
            return None
        

    def query_gene_name(self, gene_name: str) -> int|set:
        """
        Returns the cluster_id of the cluster, in which the gene
        with the given name is located.
        !!! Returns a set of cluster_ids if gene_name was indexed multiple times.
        Raises None if it can not be found.
        """
        if gene_name in self._index['gene_name'].keys():
            return self._index['gene_name'][gene_name]
        else:
            return None

def _get_index_path() -> str:
    return os.path.join(config['general']['data_dir'], 'cluster_index.dill')

def _get_pre_generated_index_path() -> str:
    return os.path.join(config['general']['data_dir'], 'default_cluster_index.dill')

def load_index(index='own') -> ClusterIndex:
    """
    Load a ClusterIndex from a binary file. 
    Default option loads ClusterIndex generated with setup_index.py.
    Set index=\'pre-generated\' to load the pre-generated ClusterIndex.
    """
    if index == 'own':
        index_path = _get_index_path()
    elif index == 'pre-generated':
        index_path = _get_pre_generated_index_path()
    else:
        raise ValueError
    return dill_load_object(index_path)