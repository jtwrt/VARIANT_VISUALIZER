from ._regions import Region
from ._conversion import map_location_to_regions, map_location_from_regions
from ._bio_references import _BioReference, GenomicReference, TranscriptReference, ProteinReference, get_reference
from math import ceil

def convert_location(location: int, input_reference: _BioReference, output_reference: _BioReference, encoding_base=None, _return_floats=False):
    """
    Convert a BioRegion between from input_refernce to output_reference.
    (Location in genomic region, location in transcript/genomic subregion or location in protein.)
    input_reference and output_refernce can be either GenomicReference, TranscripReference or ProteinReference.
    encoding_base can be either \'first\', \'second\' or \'third\' and allow to define which encoding_base will be returned,
    when converting from a ProteinReference (amino acid sequence) to a GenomicReference or TranscriptReference (nucleotide sequence).
    """
    if location <= 0:
        raise ValueError(f'Convertible locations can not be smaller than 1.')
    if not input_reference.convertible(output_reference):
        raise ValueError(f'Missing BioReference arguments nececarry for conversion.')

    if isinstance(input_reference, GenomicReference):
        # transfrom - strand regions to negative space
        if input_reference.strand == '+':
            coding_regions = input_reference.coding_regions
            transcript_region = input_reference.transcript_region
        elif input_reference.strand == '-':
            location = -location
            if input_reference.coding_regions is not None:
                coding_regions = [_invert_region(r) for r in input_reference.coding_regions]
            if input_reference.transcript_region is not None:
                transcript_region = _invert_region(input_reference.transcript_region) 
    

        if isinstance(output_reference, TranscriptReference):
            return map_location_to_regions(location=location,
                                            regions=set([transcript_region]))

        elif isinstance(output_reference, ProteinReference):
            base =  map_location_to_regions(location=location,
                                             regions=set(coding_regions))
            if _return_floats is True:
                return base / 3
            else:
                return ceil(base / 3)

    elif isinstance(input_reference, TranscriptReference):
        if isinstance(output_reference, GenomicReference):

            if input_reference.strand == '+':
                transcript_region = input_reference.transcript_region 
            elif input_reference.strand == '-' and input_reference.transcript_region is not None:
                transcript_region = _invert_region(input_reference.transcript_region) 
            else:
                raise ValueError

            out = map_location_from_regions(location=location,
                                              regions=set([transcript_region]))
            if input_reference.strand == '-': 
                return -out
            else: 
                return out

        elif isinstance(output_reference, ProteinReference):
            if _return_floats is True:
                return location/3
            else:
                return ceil(location / 3)

    elif isinstance(input_reference, ProteinReference):
        if  encoding_base is None or encoding_base not in ['first', 'second', 'third']:
            raise ValueError(f'Argument \'encoding_base\' must be specified when converting from protein location. Please specify \'first\', \'second\' or \'third\' base which encode the amino acid.')

        if encoding_base == 'third':
            i = 0
        elif encoding_base == 'second':
            i = -1
        elif encoding_base == 'first':
            i = -2

        if isinstance(output_reference, GenomicReference):
            
            if input_reference.strand == '+':
                coding_regions = input_reference.coding_regions
            elif input_reference.strand == '-' and input_reference.coding_regions is not None:
                coding_regions = [_invert_region(r) for r in input_reference.coding_regions]
            else:
                raise ValueError

            out = map_location_from_regions(location=location*3+i,
                                             regions=set(coding_regions))
            if input_reference.strand == '-': 
                return -out
            else: 
                return out

        elif isinstance(output_reference, TranscriptReference):
            # convert coding_regions from genomic reference to transcript reference
            transcript_coding_regions = []
            genomic_ref = get_reference('genomic', 
                    chromosome='1',  #arbitrary chromosome
                    strand=input_reference.strand, 
                    coding_regions=input_reference.coding_regions, 
                    transcript_region=input_reference.transcript_region) 
            for r in input_reference.coding_regions:
                new_reg = Region(start=convert_location(r.start, genomic_ref, output_reference),
                                  end=convert_location(r.end, genomic_ref, output_reference)
                                  )
                transcript_coding_regions.append(new_reg)
            # convert location to transcript_reference
            return map_location_from_regions(location=location*3+i,
                                              regions=set(transcript_coding_regions))
    
    raise NotImplementedError('Reference types are not supported for conversion.')
    
def _invert_region(region: Region):
    """Converts region with start a and end b to region with start -b and end -a"""
    return Region(start=-region.end,
                   end=-region.start)