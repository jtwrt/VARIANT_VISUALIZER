from ._regions import Region, remove_regions, combine_regions, ints_to_regions, get_all_overlaps, get_gap
from ._conversion import map_location_from_regions, map_location_to_regions, OutOfBoundsException
from ._utils import check_strand, get_other_strand
from ._bio_references import get_reference, _BioReference, GenomicReference, TranscriptReference, ProteinReference
from ._bio_conversion import convert_location
from ._bio_regions import BioRegion