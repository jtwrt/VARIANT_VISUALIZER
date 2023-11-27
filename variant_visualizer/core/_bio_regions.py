from __future__ import annotations
from ._regions import Region
from ._bio_references import _BioReference, GenomicReference, TranscriptReference, ProteinReference, get_reference
from ._bio_conversion import convert_location
from copy import deepcopy

class BioRegion(Region):
    """
    Description
    ---
    Region within a protein, transcript or genome/chromosome.

    Parameters
    ---
    start : int
    end : int 
    reference : modules.core.bio_references._BioReference
                Assign reference object created using the get_reference function 
    """

    def __init__(self, start: int, end: int, reference: _BioReference, label=None):
        super().__init__(start, end)
        self.reference = reference
        self.label = label

    def __key(self):
        return ('BioRegion', self.start, self.end, self.reference)
    
    def __hash__(self) -> tuple:
        return hash(self.__key())

    def __repr__(self):
        return f'BioRegion {self.start}-{self.end}, reference: {self.reference}'
    
    def __eq__(self, other: BioRegion):
        """Compares BioRegions. Must be implemented by subclasses."""
        this_class = BioRegion
        if isinstance(other, this_class):
            return self.__key() == other.__key()
        elif isinstance(other, BioRegion) or other is None:
            return False
        else:
            raise TypeError(f'Testing equality not defined between given objects.')

    def get_reference_type(self):
        return self.reference.reference_type
    
    def convert(self, new_reference: _BioReference) -> BioRegion:
        """
        Description
        ---
        Converts this BioRegion from its BioReference to a new BioReference.
        Returns a new BioRegion with converted start and end boundaries,

        Raises
        ----------
        OutOfBoundsException
            If BioRegion is not encoded in Reference. 
        """
        missing_args = self.reference.list_missing_conversion_args(new_reference)
        if len(missing_args) > 0:
            raise ValueError(f'BioReference is missing required information to convert Region. Missing: \'{missing_args}\'.')  
        else:
            out = self.get_deepcopy()
            new_reference = deepcopy(new_reference)
            start = convert_location(location=self.start,
                                      input_reference=self.reference,
                                      output_reference=new_reference,
                                      encoding_base='first')
            end = convert_location(location=self.end,
                                    input_reference=self.reference,
                                    output_reference=new_reference,
                                    encoding_base='third')
            new_reference.update(self.reference)
            out.update(start=start,
                       end=end,
                       reference=new_reference)
            
        return out

    def get_opposite_regions(self, regions: list) -> list:
        """
        Description
        ---
        From a list of regions, selects regions 
        that are overlapping self on the opposite strand.
        self.reference and reference property of regions in provided list
        must be of type bio_references.GenomicReference.
        """

        if not isinstance(self.reference, GenomicReference):
            raise ValueError(f'self.reference must be of type \'GenomicReference\'.')
        opposite_strand = self.reference.get_opposite_strand()

        out = []
        for r in regions:
            if not isinstance(r.reference, GenomicReference):
                raise ValueError(f'region.reference must be of type \'GenomicReference\'.')
            if r.reference.strand != opposite_strand:
                continue
            if self.overlaps(r):
                out.append(r)
        return out
    
    def to_bed(self, name, score=0, strand=None, additional_columns='') -> str:
        """
        Description
        ---
        Returns the BioRegion as string in bed-format.
        self.reference must be genomic reference.
        If reference.strand is not defined, it may be given directly to this function.
        Returned string ends with '\n'

        Parameters
        ---
        name : str
               Value to be added to the name column (4th)
        score : int
                Value between 0 and 100 to be added to the score column (5th)
        strand : str
                 + or - if reference.strand is not definded
        additional_columns : str
                             columns to be appended to the bed-string, must start
                             with '\t' and use it as delimiter between columns 
        """
        if strand is None:
            strand = self.reference.strand
        if not isinstance(self.reference, GenomicReference):
            raise TypeError(f'Only genomic references can be converted to bed-style string. Given region: {self}')
        return f'{self.reference.chromosome}\t{self.start-1}\t{self.end}\t{name}\t{score}\t{strand}{additional_columns}\n'