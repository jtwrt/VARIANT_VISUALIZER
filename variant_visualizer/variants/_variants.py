from __future__ import annotations
from .. import core

class Variant(core.BioRegion):

    def __init__(self, start: int, end: int, reference: core._BioReference,
                 variant_type: str, consequence: str, sample_id: str, normal_sample_id: str,
                 disease: str, source: str, 
                 ref_allele: str, alt_allele_1: str, alt_allele_2: str, label=None
                 ):
        super().__init__(start, end, reference, label=label)
        variant_attributes = [
            variant_type, consequence, sample_id, normal_sample_id,
            disease, source, ref_allele, alt_allele_1, alt_allele_2
            ]
        if not all([isinstance(attribute, str) for attribute in variant_attributes]):
            raise TypeError()
        self.variant_type = variant_type
        self.consequence = consequence
        self.sample_id = sample_id
        self.normal_sample_id = normal_sample_id
        self.disease = disease
        self.source = source
        self.ref_allele = ref_allele
        self.alt_allele_1 = alt_allele_1
        self.alt_allele_2 = alt_allele_2

    def _get_variant_attributes(self):
        variant_attributes = [
            self.variant_type, self.consequence, self.sample_id, self.normal_sample_id,
            self.disease, self.source, self.ref_allele, self.alt_allele_1, self.alt_allele_2
            ]
        return variant_attributes

    def __key(self):
        out = ('Variant', self.start, self.end, self.reference)
        variant_attributes = self._get_variant_attributes()
        for attribute in variant_attributes:
            out += f', {attribute}'
        return out
    
    def __hash__(self) -> tuple:
        return hash(self.__key())

    def __repr__(self):
        out = f'Variant {self.start}-{self.end}, reference: {self.reference}'
        variant_attributes = self._get_variant_attributes()
        for attribute in variant_attributes:
            out += f', {attribute}'
        return out

    def __eq__(self, other: Variant):
        """Compares BioRegions. Must be implemented by subclasses."""
        this_class = Variant
        if isinstance(other, this_class):
            return self.__key() == other.__key()
        elif isinstance(other, core.BioRegion) or other is None:
            return False
        else:
            raise TypeError(f'Testing equality not defined between given objects.')

    @property
    def label(self):
        return f'{self.reference_allele} > {self.tumor_allele_1}, {self.tumor_allele_2}; {self.consequence}'
    
    def is_snv(self):
        if self.variant_type == 'SNP':
            return True
        else:
            return False
        
    def is_indel(self):
        if self.variant_type in ['INS','DEL']:
            return True
        else:
            return False