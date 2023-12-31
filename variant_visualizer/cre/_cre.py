from __future__ import annotations
from .. import core

class CRE(core.BioRegion):

    def __init__(self, start: int, end: int, reference: core._BioReference, source: str, info=dict(), label=None):
        """
        BioRegion object that requires additional input:
        source: str
            String allowing identification of the input source for this CRE
        info: str
            String attribute can be used to place additional information.
        """
        super().__init__(start, end, reference, label)
        if not isinstance(source, str) or not isinstance(info, dict):
            raise TypeError()
        self.source = source
        self.info = info

    def __key(self):
        info = f'{" ".join([f"{key}:{self.info[key]}" for key in self.info])}'
        return ('CRE', self.start, self.end, self.reference, self.source, info)
    
    def __hash__(self) -> tuple:
        return hash(self.__key())

    def __repr__(self):
        return f'BioRegion {self.start}-{self.end}, {self.reference.__repr__()}, source: {self.source}, info: {self.info}'
    
    def __eq__(self, other: CRE):
            """Compares BioRegions. Must be implemented by subclasses."""
            this_class = CRE
            if isinstance(other, this_class):
                return self.__key() == other.__key()
            elif isinstance(other, core.Region) or other is None:
                return False
            else:
                raise TypeError(f'Testing equality not defined between given objects.')