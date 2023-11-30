from __future__ import annotations
from ... import core
from .._cre import CRE

class PAS(CRE):

    def __init__(self, start: int, end: int, reference: core._BioReference, 
                 source: str, sequence: str, cleavage_site: core.BioRegion, info=dict(),
                 label=None
                 ):
        """
        Polyadenylation and cleavage site.
        CRE object that requires additional input:
        source: str
            String allowing identification of the input source for this CRE
        sequence: str
            nucleotide sequence of the signal
        cleavage_site: core.BioRegion
            cleavage site associated with this PAS
        """
        super().__init__(start, end, reference, source, info, label=label)
        if not isinstance(sequence, str) or \
                not isinstance(cleavage_site, core.BioRegion) or \
                not isinstance(source, str):
            raise TypeError()
        self.sequence = sequence
        self.cleavage_site = cleavage_site
        self.source = source

    def __key(self):
        info = f'{" ".join([f"{key}:{self.info[key]}" for key in self.info])}'
        return ('PAS', self.start, self.end, self.reference, self.source, info, self.sequence, self.cleavage_site)
    
    def __hash__(self) -> tuple:
        return hash(self.__key())

    def __repr__(self):
        return f'''PAS {self.start}-{self.end}, reference: {self.reference}, 
            source: {self.source}, info: {self.info}, sequence: {self.sequence}, cleavage_site: {self.cleavage_site}'''
    
    def __eq__(self, other: PAS):
            """Compares BioRegions. Must be implemented by subclasses."""
            this_class = PAS
            if isinstance(other, this_class):
                return self.__key() == other.__key()
            elif isinstance(other, core.Region) or other is None:
                return False
            else:
                raise TypeError(f'Testing equality not defined between given objects.')