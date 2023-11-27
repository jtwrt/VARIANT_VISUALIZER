class OutOfBoundsException(Exception):
    """
    Raised, when conversion to new reference is not possible,
    because the new location is not defined by the given reference.
    """
    pass

def clean_reference_regions(regions: set):
    """Sorts regions and makes sure that none are overlapping."""
    if not isinstance(regions, set):
        raise TypeError
    n_regions = len(regions)
    if n_regions == 0:
        raise ValueError(f'Empty list given.')
    elif n_regions == 1:
        return regions
    else:
        regions = sorted(regions, key=lambda r: r.start)
        for i, a in enumerate(regions[0:n_regions-1]):
            for b in regions[i+1:]:
                if a.overlaps(b):
                    raise ValueError(f'Cannot map to overlapping regions.')
        return regions

def map_location_to_regions(location: int, regions: set) -> int:
    """
    Description
    ----------
    Maps a location in one region 1-n to a 
    relative location in a group of region (1-a, b-c, d-e, ...)

    Parameters
    ----------
    location: int
    regions : list of regions
    """
    regions = clean_reference_regions(regions)

    n_from_start = 0
    for r in regions:
        if location >= r.start and location <= r.end:
            return location - r.start + n_from_start + 1
        n_from_start += r.end - r.start + 1  # Add one as end value is included in region

    raise OutOfBoundsException(f'Given location is not within provided regions.')

def map_location_from_regions(location: int, regions: set) -> int:
    """
    Description
    ----------
    Maps a location within a group of region (1-a, b-c, d-e, ...)
    to a absolute location 1-n

    Parameters
    ----------
    location: int
    regions : list of regions
    """
    
    regions = clean_reference_regions(regions)

    skip_n = location - 1
    skipped = 0  # Counter to keep track of skipped bases

    for r in regions:
        length = r.end - r.start + 1

        if skipped == skip_n:
            return r.start
        
        elif skipped < skip_n:
            if skipped + length <= skip_n:
                skipped += length
                continue
            else:  # Skip remaining part of m
                return r.start + skip_n - skipped
            
    raise OutOfBoundsException(f'Given location is not within provided regions.')