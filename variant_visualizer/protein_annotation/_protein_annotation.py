from variant_visualizer.core._bio_references import _BioReference
from .. import core

class Annotation(core.BioRegion):

    def __init__(self, start: int, end: int, reference: _BioReference, annotation_type: str, description: str, source: str):
        super().__init__(start, end, reference)
        self.annotation_type = annotation_type
        self.description = description
        self.source = source

    def __key(self):
        return ('Annotation', self.start, self.end, self.reference, self.annotation_type,
                self.text, self.source)
    
    def __hash__(self) -> tuple:
        return hash(self.__key())

    def __repr__(self):
        return f'BioRegion {self.start}-{self.end}, reference: {self.reference}'
    
    def __eq__(self, other: core.BioRegion):
        """Compares BioRegions. Must be implemented by subclasses."""
        if isinstance(other, Annotation):
            return self.__key() == other.__key()
        elif isinstance(other, core.BioRegion) or other is None:
            return False
        else:
            raise TypeError(f'Testing equality not defined between given objects.')
