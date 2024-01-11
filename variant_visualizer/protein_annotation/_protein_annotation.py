from variant_visualizer.core._bio_references import _BioReference
from .. import core

class ProteinAnnotation(core.BioRegion):

    def __init__(self, start: int, end: int, reference: _BioReference, annotation_type: str, description: str, source: str, label=str):
        super().__init__(start, end, reference, label)
        self.annotation_type = annotation_type
        self.description = description
        self.source = source

    def __key(self):
        return ('ProteinAnnotation', self.start, self.end, self.reference, self.annotation_type,
                self.description, self.source)
    
    def __hash__(self) -> tuple:
        return hash(self.__key())

    def __repr__(self):
        return f'ProteinAnnotation {self.start}-{self.end}, reference: {self.reference}, annotation_type: {self.annotation_type}, description: {self.description}, source: {self.source}'
    
    def __eq__(self, other: core.BioRegion):
        """Compares BioRegions. Must be implemented by subclasses."""
        if isinstance(other, ProteinAnnotation):
            return self.__key() == other.__key()
        elif isinstance(other, core.BioRegion) or other is None:
            return False
        else:
            raise TypeError(f'Testing equality not defined between given objects.')
