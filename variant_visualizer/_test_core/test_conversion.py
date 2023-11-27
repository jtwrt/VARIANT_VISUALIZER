from ..core._regions import Region
from ..core._conversion import map_location_to_regions, map_location_from_regions

def test_map_location_to_regions_1():
    regions = set([Region(1,5),
               Region(10,15)])
    assert map_location_to_regions(location=3, regions=regions) == 3

def test_map_location_to_regions_2():
    regions = set([Region(1,5),
               Region(10,15)])
    assert map_location_to_regions(location=12, regions=regions) == 8

def test_map_location_from_regions_1():
    regions = set([Region(1,5),
               Region(10,15)])
    assert map_location_from_regions(location=3, regions=regions) == 3
    
def test_map_location_from_regions_2():
    regions = set([Region(1,5),
               Region(10,15)])
    assert map_location_from_regions(location=6, regions=regions) == 10
    
def test_map_location_from_regions_3():
    regions = set([Region(1,5),
               Region(10,15)])
    assert map_location_from_regions(location=8, regions=regions) == 12   
