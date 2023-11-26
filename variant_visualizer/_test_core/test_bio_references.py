from ..core._regions import Region
from ..core._bio_references import get_reference

def test_get_reference_1():
    tA = get_reference('transcript', 
                       transcript_region=Region(100, 2000))
    tB = get_reference('transcript')
    assert tA != tB

def test_get_reference_2():
    tA = get_reference('transcript', 
                       transcript_region=Region(100, 2000),
                       strand='+')
    pA = get_reference('protein')
    assert tA.convertible(pA) is True

def test_eq_1():
    gA = get_reference('genomic',
                       chromosome='Y')
    gB = get_reference('genomic',
                       chromosome='Y')
    assert gA == gB

def test_eq_2():
    gA = get_reference('genomic',
                       chromosome='Y')
    gB = get_reference('protein',
                       chromosome='Y')
    assert gA != gB
