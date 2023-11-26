import pytest
from ..core._regions import Region
from ..core._bio_regions import BioRegion
from ..core._conversion import OutOfBoundsException
from ..core._bio_references import get_reference

def test_convert_1():
    bio_region = BioRegion(start=4,
                           end=6,
                           reference=get_reference('genomic', transcript_region=Region(3,4), chromosome='10', strand='+')
                           )
    with pytest.raises(OutOfBoundsException):
        bio_region.convert(get_reference('transcript'))

def test_convert_2():
    bio_region = BioRegion(start=3,end=3,
                           reference=get_reference('genomic', transcript_region=Region(3,4), chromosome='10', strand='+')
                           )
    out = bio_region.convert(get_reference('transcript'))
    assert (out.start, out.end) == (1, 1)

def test_convert_3():
    bio_region = BioRegion(start=7,end=8,
                           reference=get_reference('genomic', coding_regions=[Region(3,4), Region(7,9)], chromosome='10', strand='+') 
                           )
    out = bio_region.convert(get_reference('protein'))
    assert (out.start, out.end) == (1, 2)

def test_convert_4():
    bio_region = BioRegion(start=2,end=2,
                           reference=get_reference('transcript', transcript_region=Region(3,4), chromosome='10', strand='+')
                           )
    out = bio_region.convert(get_reference('genomic', chromosome='10', strand='+'))
    assert (out.start, out.end) == (4, 4)

def test_convert_5():
    bio_region = BioRegion(start=3,
                           end=4,
                           reference=get_reference('transcript', strand='+'))
    out = bio_region.convert(get_reference('protein'))
    assert (out.start, out.end) == (1, 2)

def test_convert_6():
    bio_region = BioRegion(1,1,
                           reference=get_reference('protein', coding_regions=[Region(3,4), Region(7,9)], chromosome='10', strand='+'))
    out = bio_region.convert(get_reference('genomic', chromosome='10', strand='+'))
    assert (out.start, out.end) == (3, 7)

def test_convert_7():
    bio_region = BioRegion(2,4,
                           get_reference('protein', transcript_region=Region(1,12), coding_regions=[Region(1,12)], strand='+'))
    out = bio_region.convert(get_reference('transcript'))
    assert (out.start, out.end) == (4, 12)

def test_convert_8():
    bio_region = BioRegion(2,4,
                           get_reference('protein', transcript_region=Region(1,12), coding_regions=[Region(1,12)], strand='-'))
    out = bio_region.convert(get_reference('transcript'))
    assert (out.start, out.end) == (4, 12)

def test_convert_9():
    bio_region = BioRegion(2,4,
                           get_reference('protein', transcript_region=Region(1,15), coding_regions=[Region(1,2),Region(6,15)], strand='+'))
    out = bio_region.convert(get_reference('transcript'))
    assert (out.start, out.end) == (7, 15)

def test_convert_10():
    bio_region = BioRegion(2,4,
                           get_reference('protein', transcript_region=Region(1,15), coding_regions=[Region(1,2),Region(6,15)], strand='-'))
    out = bio_region.convert(get_reference('transcript'))
    assert (out.start, out.end) == (4, 15)  

def test_eq_1():
    a = BioRegion(1,2,get_reference('genomic',chromosome='1'))
    b = BioRegion(1,2,get_reference('genomic',chromosome='2'))
    assert a != b

def test_eq_2():
    a = BioRegion(1,2,get_reference('genomic',chromosome='1'))
    b = BioRegion(1,2,get_reference('transcript',chromosome='1'))
    assert a != b

def test_eq_3():
    a = BioRegion(1,2,get_reference('genomic',chromosome='1',strand='+'))
    b = BioRegion(1,2,get_reference('genomic',chromosome='1',strand='-'))
    assert a != b

def test_key_1():
    a = BioRegion(1,2,get_reference('genomic',chromosome='1',strand='+'))
    b = BioRegion(1,2,get_reference('genomic',chromosome='1',strand='-'))
    assert set([a,b]) == set([b,a])