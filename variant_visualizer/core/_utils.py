def check_strand(strand: str):
    """
    Checks if the given strand is a str and its value is valid.
    Valid values are \'+\' and \'-\'.
    Returns True if the strand value is valid.
    Raises TypeError or ValueError otherwise.
    """
    if not isinstance(strand, str):
        raise TypeError("Strand values must be str.")
    if not strand in ['+','-']:
        raise ValueError("Strand value mus be either \'+\' or \'-\'.")
    return True

def get_other_strand(strand: str):
    """Returns opposite strand for a given \'+\' or \'-\'."""
    _ = check_strand(strand=strand)
    if strand == '+':
        return '-'
    if strand == '-':
        return '+'
    raise BaseException # never reached