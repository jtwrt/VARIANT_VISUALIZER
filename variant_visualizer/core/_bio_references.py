from __future__ import annotations
from ._regions import Region
from ._utils import check_strand, get_other_strand
from copy import deepcopy


def get_reference(reference_type: str, *args, **kwargs):
    """
    Description
    ----------
    Create a reference object used to initialize a BioRegion.
    This function returns objects of subclasses of the _BioReference class
    based on the reference_type chosen.
    Each subclass may require additional arguments.

    _BioReference instances allow references locations to be converted between
    subclasses, if all required parameters of the target subclass are set and 
    additional conversion arguments are set.
    
    Parameters
    ----------
    reference_type : str
                     Choose the locations that this object is referencing. 
                     Valid values are 'genomic', 'transcript' or 'protein'.
                     
        Additional parameters
        ---
        reference_type == 'genomic'
            chromosome : str

        Conversion parameters
        ---
        genomic <---> transcript
            transcript_region : list

        genomic <---> protein
            coding_regions : list
    """
    if reference_type == 'genomic':
        return GenomicReference(*args, **kwargs)
    elif reference_type == 'transcript':
        return TranscriptReference(*args, **kwargs)
    elif reference_type == 'protein':
        return ProteinReference(*args, **kwargs)
    else:
        raise NotImplementedError(f'Unknown reference_type \'{reference_type}\'.')
    
class _BioReference():
    """
    Parent class. Child classes are used as attributes of BioRegion and allow location conversion between references.
    Child classes must define the class attributes reference_type, mandatory_args and conversion_args.
    reference_type must be a str, and is used in bio_conversion module to identify the correct reference class.
    mandatory_args contains a list of tuples with (attribute_name, attribute_type).
    conversion_args contains a dict of lists with the same format as mandatory_args, the keys are reference_type of
    the different classes. conversion_args contains information of which attributes are needed for successfull conversion
    between references.
    """

    reference_type = None
    mandatory_args = []
    conversion_args = {}

    @classmethod
    def list_all_conversion_args(cls):
        all_conversion_args = set()
        for reference_type in cls.conversion_args:
            for arg in cls.conversion_args[reference_type]:
                all_conversion_args.add(arg)
        return list(all_conversion_args)

    @classmethod
    def list_all_args(cls):
        return cls.mandatory_args + cls.list_all_conversion_args()


    def __init__(self, *args, **kwargs):
        
        # Get list of mandatory arguments
        mandatory_args = self.__class__.mandatory_args

        # Add kwargs
        self.__dict__.update(kwargs)
        
        # Add args
        n_undefined = 0
        for kw, _ in mandatory_args:
            if kw not in self.__dict__:
                n_undefined += 1 
        if len(args) > n_undefined:
            raise ValueError(f'Too many arguments. Optional arguments must be defined with keyword.')
        elif len(args) < n_undefined:
            raise ValueError(f'Missing mandatory arguments.')
        else:
            for i, arg in enumerate(args):
                kw = mandatory_args[i][0]
                if kw not in self.__dict__:
                    self.__dict__[kw] = arg
                else: 
                    raise ValueError(f'Argument \'{kw}\' is given twice.')

        # Check if all mandatory arguments are defined and have the correct type
        for kw, arg_type in mandatory_args:
            if kw not in self.__dict__:
                raise ValueError(f'Missing mandatory argument \'{kw}\'.')
            if not isinstance(self.__dict__[kw], arg_type):
                raise TypeError(f'Type of argument \'{kw}\' must be \'{arg_type}\' not \'{type(self.__dict__[kw])}\'.')

        # Check type of defined conversion args and initialize missing values
        for kw, arg_type in self.list_all_conversion_args():
            if kw not in self.__dict__:
                self.__dict__[kw] = None
            elif not isinstance(self.__dict__[kw], arg_type) and self.__dict__[kw] is not None:
                raise TypeError(f'Type of argument \'{kw}\' must be \'{arg_type}\' not \'{type(self.__dict__[kw])}\'.')

        # Check strand value if provided, allow None value
        if "strand" in self.__dict__.keys() and self.strand is not None:
            _ = check_strand(strand=self.strand)

        # convert reference regions to basic Region objects
        self.convert_reference_regions()


    def convert_reference_regions(self):
        """Convert transcript_region and coding_regions to basic Region objects"""

        def convert_region(region: Region):
            if isinstance(region, Region):
                return Region(region.start, region.end)
            else:
                raise TypeError

        if 'transcript_region' in self.__dict__.keys() and isinstance(self.transcript_region, Region):
            self.transcript_region = convert_region(self.transcript_region)
        if 'coding_regions' in self.__dict__.keys() and \
                self.coding_regions is not None and \
                all([isinstance(r, Region) for r in self.coding_regions]):
            self.coding_regions = set([convert_region(r) for r in self.coding_regions])
        

    def __repr__(self) -> str:
        this_class = str(self.__class__).split('.')[-1][:-2]
        out = [str(this_class)]
        for kw, _ in self.list_all_args():
            info = self.__dict__[kw]
            if isinstance(info, list):
                info = f'{len(info)} objects'
            if isinstance(info, Region):
                info = f'{info.start}-{info.end}'
            out.append(f'{kw}: {info}')
        return ', '.join(out)
                
    def __str__(self) -> str:
        return repr(self)
    
    def __key(self):
        this_class = str(self.__class__).split('.')[-1][:-2]
        out = (this_class)
        for kw, _ in self.list_all_args():
            out += f'{kw}:{str(self.__dict__[kw])}'
        return out
    
    def __hash__(self):
        return hash(self.__key())
    
    def __eq__(self, other: _BioReference):
        """Two BioReferences (or childclasses) are equal, if all mandatory_args and conversion_args are equal."""
        if isinstance(other, type(self)):
            return self.__key() == other.__key()
        elif isinstance(other, _BioReference):
            return False
        else:
            raise NotImplementedError(f'Testing equality not defined with type \'{type(other)}\'.')
                
    def update(self, other: _BioReference, force=False):
        """
        Add new attributes from another reference to this reference.
        Throws ValueError if self.attribute and other.attribute are conflicting.
        """
        other_args = other.__class__.mandatory_args + other.__class__.list_all_conversion_args()
        for kw, _ in other_args:
            if other.__dict__[kw] is None:
                continue
            elif other.__dict__[kw] is not None:
                if kw in self.__dict__ and self.__dict__[kw] is None:
                    self.__dict__[kw] = other.__dict__[kw]
                elif kw in self.__dict__ and not self.__dict__[kw] == other.__dict__[kw]:
                    if not force is True:
                        raise ValueError(f'Properties \'{kw}\' are set in both references and differ.')
                    elif force is False:
                        self.__dict__[kw] = other.__dict__[kw]
        self.convert_reference_regions()
        
    def list_missing_conversion_args(self, other: _BioReference) -> list:
        """Returns attributes that are missing to allow conversion between _BioReference child classes."""
        if type(self) == type(other):
            raise ValueError(f'Cannot convert between BioReferences of the same type.')
        missing = []
        for kw, _ in other.__class__.mandatory_args:
            if kw not in self.__dict__ or self.__dict__[kw] is None:
                missing.append(kw)
        for kw, _ in self.__class__.conversion_args[other.__class__.reference_type]:
            if kw not in self.__dict__ or self.__dict__[kw] is None:
                missing.append(kw)
        return missing
    
    def convertible(self, other: _BioReference) -> bool:
        """
        Check if self has all required properties 
        to convert to other _BioReference child class.
        """
        missing_args = self.list_missing_conversion_args(other=other)
        if len(missing_args) == 0:
            return True
        else:
            return False
    
    def get_opposite_strand(self):
        """
        Returns the opposite strand. If this GenomicReference has no strand attribute,
        raises AttributeError.
        """
        if 'strand' not in self.__dict__.keys():
            raise AttributeError("No strand attribute defined in this GenomicReference.")
        return get_other_strand(self.strand)

class GenomicReference(_BioReference):
    """
    Contains reference information for a genomic region.
    """

    reference_type = 'genomic'
    mandatory_args = deepcopy(_BioReference.mandatory_args)
    mandatory_args.append(('chromosome', str))
    conversion_args = {'transcript': [('transcript_region', Region)],
                       'protein': [('coding_regions', set),
                                   ('strand', str)]
                       }
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

class TranscriptReference(_BioReference):
    
    reference_type = 'transcript'
    mandatory_args = deepcopy(_BioReference.mandatory_args)
    conversion_args = {'genomic': [('transcript_region', Region)],
                       'protein': [('strand', str)]
                       }
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
class ProteinReference(_BioReference):
    
    reference_type = 'protein'
    mandatory_args = deepcopy(_BioReference.mandatory_args)
    conversion_args = {'genomic': [('coding_regions', set)],
                       'transcript': [('coding_regions', set),
                                      ('transcript_region', Region),
                                      ('strand', str)]
                       }
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

