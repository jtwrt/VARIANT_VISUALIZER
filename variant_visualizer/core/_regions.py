from __future__ import annotations
from copy import deepcopy

class Region(object):

    def __init__(self, start: int, end: int):
        """
        Description
        ---
        Create instance of Region class. Start and end values are inclusive.
        Child classes are to inplement the following functions at least:
        __init__, __key, __hash__, __repr__
        
        Parameters
        ---
        start : int
        end : int 
        """ 
        self.start, self.end = self._sort_locations(
            start=start,
            end=end)

    def _sort_locations(self, start, end) -> tuple:
        """Makes sure start is <= end. Returns start, end"""
        if not isinstance(start, int) or \
                not isinstance(end, int):
            raise TypeError("Boundary attributes \'start\' and \'end\' must be integers.")
        bounds = sorted([start, end])
        return bounds[0], bounds[1]

    def __key(self) -> tuple:
        return ('Region', self.start, self.end)
    
    def __hash__(self) -> tuple:
        """Allows creation of sets and dicts using Region as keys."""
        return hash(self.__key())
    
    def __eq__(self, other: Region):
            this_class = Region
            if isinstance(other, this_class):
                return self.__key() == other.__key()
            elif isinstance(other, Region) or other is None:
                return False
            else:
                raise TypeError(f'Testing equality not defined between given objects.')
    
    def __lt__(self, other: Region):
        """Allows sorting of regions to sort by start value."""
        if self.end < other.start:
            return True
        else:
            return False

    def __repr__(self) -> str:
        return f'Region {self.start}-{self.end}'
    
    def __str__(self) -> str:
        return repr(self)
    
    def get_deepcopy(self) -> Region:
        return deepcopy(self)
    
    def update(self, **kwargs):
        """Update properties of this instance."""
        self.__dict__.update(kwargs)
        if self.start > self.end:
            tmp = self.start
            self.start = self.end
            self.end = tmp

    def overlaps(self, other: Region):
        """Test if self overlaps other Region object."""
        if not self < other and not self > other:
            return True
        else:
            return False

    def touches(self, other: Region):
        """Test if self overlaps, or is right next to other Region object."""
        if self.end+1 == other.start or other.end+1 == self.start or self.overlaps(other):
            return True
        else:
            return False
        
    def within(self, other: Region):
        """Test if self start-end range is part of other start-end range."""
        if self.start >= other.start and self.end <= other.end:
            return True
        else:
            return False

    def enlarge_with(self, other:Region):
        """
        If self touches other, add other start-end range to own and return True
        Otherwise return False
        """
        if not self.touches(other):
            return False
        else:
            self.start = min((self.start, other.start))
            self.end = max((self.end, other.end))
            return True
    
    def get_length(self):
        return self.end - self.start + 1

    def get_overlap(self, other: Region):
        """
        Returns a new Region object with the start-end range 
        in which self and other overlap.
        Raises ValueError if self and other do not overlap.
        """
        if not self.overlaps(other):
            raise ValueError('Regions not overlapping.')
        out = self.get_deepcopy()
        if self.within(other):
            out.update(start=self.start, end=self.end)
        elif other.within(self):
            out.update(start=other.start, end=other.end)
        elif self.start < other.start:
            out.update(start=other.start, end=self.end)
        elif self.start > other.start:
            out.update(start=self.start, end=other.end)
        else: 
            raise BaseException
        return out
    
    def remove_overlap(self, other: Region) -> list:
        """Removes a overlapping Region start-end range from self start-end range. Returns a list of resulting Regions."""
        if not self.overlaps(other):
            raise ValueError('Regions not overlapping.')
        
        out = []

        if self == other:
            pass
        elif self.within(other):
            pass
        elif other.within(self) and self.start != other.start and self.end != other.end:
            out0 = self.get_deepcopy()
            out0.update(start=self.start, end=other.start-1)
            out1 = self.get_deepcopy()
            out1.update(start=other.end+1, end=self.end)
            out.extend([out0, out1])
        elif other.start <= self.start and other.end < self.end:
            out0 = self.get_deepcopy()
            out0.update(start=other.end+1, end=self.end)
            out.append(out0)
        elif other.start > self.start and other.end >= self.end:
            out0 = self.get_deepcopy()
            out0.update(start=self.start, end=other.start-1)
            out.append(out0)
        else:
            raise BaseException
        
        return out

    def split_to_locations(self) -> list:
        """
        Returns a list of copies of self with length 1. 
        The compination of the list is self.
        """
        out = []
        for location in range(self.start, self.end+1):
            subregion = self.get_deepcopy()
            subregion.update(start=location,
                             end=location)
            out.append(subregion)
        return out

def combine_regions(regions: list) -> list:
    """
    Description
    ---
    Tries to merge regions from the provided list.
    Regions with smallest start value are extended iteratively. 
    Attributes of the region with the smaller start value are kept, while others are lost.
    """
  
    regions = sorted(regions, key=lambda r: r.start)
    regions = deepcopy(regions)
    i=0
    while i < len(regions)-1 and len(regions) > 1:
        region0 = regions[i]
        for j in range(i+1, len(regions)):
            region1 = regions[j]
            if region0.enlarge_with(region1) is False:
                if j == len(regions)-1: 
                    i += 1 # continue while loop
                continue

            else:
                next_iter_regions = []
                for m,s in enumerate(regions):
                    if m == j:
                        continue
                    else:
                        next_iter_regions.append(s)
                regions = sorted(next_iter_regions, key=lambda r: r.start)
                i = 0 # start while loop again
                break
    return regions

def remove_regions(regions: list, remove_regions: list) -> list:
    """
    Description
    ---
    From a given list with Region objects, 
    removes or shortens regions that overlap with remove_regions.
    """
    regions = combine_regions(regions)
    remove_regions = combine_regions(remove_regions)
    new_regions = []

    for d in remove_regions:
        for r in regions:
            if r.overlaps(d):
                new_regions.extend(r.remove_overlap(d))
            else:
                new_regions.append(r)
        regions = new_regions
        new_regions = []
    return regions

def ints_to_regions(ints: list) -> list:
    """Creates a list of regions from integers. Combines length one regions whenever possible."""
    regions = [Region(i,i) for i in ints]
    return combine_regions(regions)

def get_all_overlaps(regionsA: list, regionsB: list) -> list:
    """Returns the overlaps between two lists of regions."""
    overlaps = []
    for a in regionsA:
        for b in regionsB:
            if a.overlaps(b):
                overlaps.append(a.get_overlap(b))
    return combine_regions(overlaps)

def get_gap(region_a: Region, region_b: Region):
    """Return a region that covers the gap between region_a and region_b"""
    if region_a.touches(region_b):
        raise ValueError('Given regions are touching.')
    
    regions = sorted([region_a, region_b], key=lambda r: r.start)
    return Region(start=regions[0].end+1,
                  end=regions[1].start-1)