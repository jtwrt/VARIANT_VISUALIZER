from copy import deepcopy
from ..core._regions import combine_regions, remove_regions, ints_to_regions
import pytest

from ..core._regions import Region

def test_eq_1():
    assert Region(2,4) == Region(4,2)

def test_eq_2():
    assert Region(2,4) != Region(3,5)

def test_lt_1():
    assert Region(2,4) < Region(5,6)

def test_lt_2():
    assert Region(4,4) < Region(5,6)

def test_update_1():
    regionA = deepcopy(Region(2,4))
    regionA.update(start=1, end=2)
    assert regionA == Region(1,2)

def test_get_length_1():
    assert Region(2,4).get_length() == 3

def test_overlaps_1():
    assert Region(2,4).overlaps(Region(3,5)) is True

def test_overlaps_2():
    assert Region(2,4).overlaps(Region(4,6)) is True

def test_overlaps_3():
    assert Region(2,4).overlaps(Region(5,6)) is False

def test_overlaps_4():
    assert Region(4,4).overlaps(Region(5,6)) is False

def test_overlaps_5():
    assert Region(4,4).overlaps(Region(2,4)) is True

def test_touches_1():
    assert Region(2,4).touches(Region(4,6)) is True

def test_touches_2():
    assert Region(4,4).touches(Region(4,6)) is True

def test_within_1():
    assert Region(5,6).within(Region(4,6)) is True

def test_within_2():
    assert Region(4,4).within(Region(2,4)) is True

def test_within_3():
    assert Region(2,4).within(Region(5,6)) is False

def test_within_4():
    assert Region(2,4).within(Region(4,2)) is True

def test_enlarge_with_1():
    assert Region(2,4).enlarge_with(Region(4,2)) is True

def test_enlarge_with_2():
    region = Region(2,4)
    out = region.enlarge_with(Region(4,4)) 
    assert (out, region._Region__key()) == (True, ('Region', 2, 4))

def test_enlarge_with_3():
    assert Region(5,6).enlarge_with(Region(4,4)) is True

def test_enlarge_with_4():
    assert Region(3,5).enlarge_with(Region(4,6)) is True

def test_get_overlap_1():
    assert Region(2,4).get_overlap(Region(3,5)) == Region(3,4)

def test_get_overlap_2():
    with pytest.raises(ValueError):
        Region(2,4).get_overlap(Region(5,6))

def test_remove_overlap_1():
    assert Region(4,6).remove_overlap(Region(2,4)) == [Region(5,6)]

def test_remove_overlap_2():
    assert Region(2,4).remove_overlap(Region(4,4)) == [Region(2,3)]

def test_remove_overlap_3():
    assert set(Region(3,5).remove_overlap(Region(4,4))) == set([Region(3,3), Region(5,5)])

def test_combine_regions_1():
    regions_list = [Region(2,4), Region(3,5)]
    assert combine_regions(regions_list) == [Region(2,5)]

def test_combine_regions_2():
    regions_list = [Region(1,3),
                    Region(3,6),
                    Region(8,9)]
    combined = combine_regions(regions_list)
    assert  combined == [Region(1,6),
                            Region(8,9)]

def test_remove_regions_1():
    regions_list = [Region(1,6),
                    Region(8,9)]
    out = remove_regions(regions=regions_list,
                                    remove_regions=[Region(4,6)])
    assert  set(out) == set([Region(1,3),
                                Region(8,9)])
    
def test_remove_regions_2():
    regions_list = [Region(1,6),
                    Region(8,9),
                    Region(11,15),
                    Region(12,16)]
    out = remove_regions(regions=regions_list,
                                    remove_regions=[Region(4,6),
                                                    Region(16,17)])
    assert  set(out) == set([Region(1,3),
                                Region(8,9),
                                Region(11,15)])
    
def test_ints_to_regions_1():
    ints = [1,2,3,4,6,7,8]
    out = ints_to_regions(ints=ints)
    assert set(out) == set([Region(1,4),
                           Region(6,8)])
    