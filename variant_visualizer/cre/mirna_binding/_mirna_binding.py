from __future__ import annotations
from ... import core
from .._cre import CRE

class MiRNABinding(CRE):

    def __init__(self, start: int, end: int, reference: core._BioReference, 
                 source: str, binding: str, score: float, info=dict(), label=None
                 ):
        """
        CRE object that requires additional input:
        source: str
            String allowing identification of the input source for this CRE
        binding: str
            miRNA which is bound by region
        score: float
            attribute should be used for binding site prediction scores
        """
        super().__init__(start, end, reference, source, info, label=label)
        if not isinstance(binding, str) or not isinstance(score, float):
            raise TypeError()
        self.binding = binding
        self.score = score

    def __key(self):
        return ('BioRegion', self.start, self.end, self.reference, self.source, self.info, self.binding, self.score)
    
    def __hash__(self) -> tuple:
        return hash(self.__key())

    def __repr__(self):
        return f'''MiRNABinding {self.start}-{self.end}, reference: {self.reference}, 
            source: {self.source}, info: {self.info}, binding: {self.binding}, score: {self.score}'''
    
    def __eq__(self, other: MiRNABinding):
            """Compares BioRegions. Must be implemented by subclasses."""
            this_class = MiRNABinding
            if isinstance(other, this_class):
                return self.__key() == other.__key()
            elif isinstance(other, core.Region) or other is None:
                return False
            else:
                raise TypeError(f'Testing equality not defined between given objects.')