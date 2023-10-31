from __future__ import annotations
from collections import UserDict
import warnings


from aide.base import Variant

from typing import Iterable, Union, List, Dict

class IncompleteParentageWarning(Warning):
    """Warning for when the parent for a variant could not be constructed."""
    pass

class Library(UserDict):
    """A library or variants.
    
    
    """
    def __init__(self, variants: Iterable[Variant]):
        super().__init__()
        for v in variants:
            if not isinstance(v, Variant):
                raise TypeError(f"Library must be initialized with an iterable of Variants, not {type(v)}")

        self.data = {variant.id: variant for variant in variants}

    def __repr__(self):
        return f"Library({list(self.keys())})"

    def __contains__(self, key: object) -> bool:
        if isinstance(key, Variant):
            return key.id in self
        else:
            return super().__contains__(key)
        
    @property
    def single_parent(self) -> bool:
        """Check if all variants in the library have a single parent."""
        return len(set([variant.parent for variant in self.values()])) == 1
    
    @property
    def parent(self) -> Variant:
        """Get the parent of all variants in the library."""
        if not self.single_parent:
            raise ValueError("Library does not have a single parent")
        return list(self.values())[0].parent

    def get_unlabeled(self, names: Union[str, Iterable[str]]=None, round_idx: Union[int, Iterable[int]]=None) -> Library:
        """Get a library of variants that do not have a label with the given name and round_idx.
        
        Args:
            names: The name or names of the labels to check for.
            round_idx: The round_idx or round_idxs of the labels to check for.
        
        Returns:
            A library of variants that do not have a label with the given name and round_idx.
        """
        return Library([variant for variant in self.values() if not variant.labels.has_labels(names=names, round_idx=round_idx)])
    
    def get_labeled(self, names: Union[str, Iterable[str]]=None, round_idx: Union[int, Iterable[int]]=None) -> Library:
        """Get a library of variants that have a label with the given name and round_idx.
        
        Args:
            names: The name or names of the labels to check for.
            round_idx: The round_idx or round_idxs of the labels to check for.
        
        Returns:
            A library of variants that have a label with the given name and round_idx.
        """
        return Library([variant for variant in self.values() if variant.labels.has_labels(names=names, round_idx=round_idx)])
    
    def join(self, other: Library):
        """Join another library into this one.
        
        Args:
            other: The other library to join into this one.
        """
        for variant in other.values():
            if variant.id in self:
                assert self[variant.id] == variant
                self[variant.id].labels.add_labels(variant.labels)
            else:
                self.data[variant.id] = variant

    @classmethod
    def build_variants_from_lookup(
            cls,
            variant_ids: Iterable[str],
            lookup: Dict[str, Dict[str, str]]=None,
            existing_database: 'CampaignDatabase'=None,
        ) -> Library:

        """Build variants from a lookup dictlike.
        
        Params
        ------
        variant_ids: Iterable[str]
            The variant ids to build.
        lookup: DictLike = None
            A dictlike object that maps variant ids to variant data. The variant data
            should be a dict with schema: 'sequence', 'mutations', 'parent_id', 'labels',
            'round_putative', 'round_added', 'round_experiment'. 'labels' should be a dict
            of lists with keys 'name', 'value', 'round_idx'. Only 'base_sequence' is required. 
        existing_database: CampaignDatabase = None
            A database to use to fetch parent sequences from if not in the lookup.
        """
        all_variants = {}
        
        def initialize_variant(variant_id):
            # Check if variant already exists in all_variants
            if variant_id in all_variants:
                return all_variants[variant_id]

            # Get variant data from lookup if available
            variant_data = lookup.get(variant_id, None)
            
            parent_id = None
            if variant_data:
                parent_id = variant_data.get('parent_id', None)

            # If not in lookup, then check in existing_database
            if variant_data is None and existing_database is not None:
                variant = existing_database.get_variant(variant_id)
                all_variants[variant_id] = variant
                return variant
            elif variant_data is None:
                return None
            else:
                pass

            # Initialize parent if needed
            parent = None
            if parent_id:
                parent = initialize_variant(parent_id)
                if parent is None:
                    warnings.warn(f"Could not build parentage for {variant_id}, not in lookup or database. The descendent history is incomplete.", IncompleteParentageWarning)
                else:
                    if str(parent) != variant_data['sequence']:
                        raise ValueError(f"Parent sequence {parent} does not match sequence in lookup {variant_data['sequence']}. Something is wrong.")
                    else:
                        base_sequence=parent
            try:
                if parent is None:
                    base_sequence = variant_data['sequence']
            except KeyError:
                raise KeyError(f"No parent sequence for {variant_id} in lookup, and base sequence not specified")

            # Create Variant object
            labels = variant_data.get('labels', None)
            round_putative = variant_data.get('round_putative', None)
            round_added = variant_data.get('round_added', None)
            round_experiment = variant_data.get('round_experiment', None)
            mutations = variant_data.get('mutations', None)
            variant = Variant(sequence=base_sequence, mutations=mutations, labels=labels, round_putative=round_putative, round_added=round_added, round_experiment=round_experiment, id=variant_id)

            # Add to all_variants
            all_variants[variant_id] = variant

            return variant

        # Initialize all requested variants
        for variant_id in variant_ids:
            initialize_variant(variant_id)
            if variant_id not in all_variants:
                raise KeyError(f"Could not find variant {variant_id} in lookup or database")

        return cls([all_variants[variant_id] for variant_id in variant_ids])

