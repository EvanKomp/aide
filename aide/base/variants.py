"""Base classes for dealing with protein sequences.
"""
from __future__ import annotations
from typing import Union

class Variant:
    """A protein sequence with a unique identifier.
    
    Params
    ------
    sequence : Union[str, Variant, None]
        The sequence of the variant. If a Variant is passed, the sequence
        is interpreted as the parent sequence of the variant. If None is
        passed, the parent is attempted to be infered from the ids
    mutations : Union[str, MutationSet, Mutation, None]
        The mutations that define the variant. If a string is passed, it
        will be parsed into a MutationSet.
    id : Union[str, None]
        The unique identifier of the variant. If None is passed, the id

    Attributes
    ----------
    base_sequence : str
        The sequence of the variant without mutations.
    parent : Variant or None
        The parent variant of the variant, if present
    mutations : MutationSet
        The mutations that define the variant.
    """
    def __init__(
        self,
        sequence: Union[None, str, Variant] = None,
        mutations: Union[None, str, MutationSet, Mutation] = None,
        id: Union[None, str] = None,
    ):
        self._base_sequence = None
        self._parent = None
        self.mutations = None
        # any mutations passed need to be parsed to mutation sets
        # if they are already a MutationSet, then the parent is already
        # determined. Parse from string after we have a base sequence
        
        if isinstance(mutations, Mutation):
            self.mutations = MutationSet([mutations])
        elif mutations is None:
            pass
        elif isinstance(mutations, MutationSet):
            self.mutations = mutations
        elif type(mutations) == str:
            pass
        else:
            raise ValueError('mutations must be a MutationSet, Mutation, or str')
        
        # If we have just a sequence, we have our base sequence
        # if a parent was passed, we want to get it dynamically
        # if we have to infer it from mutations, mutations must be a set already
        if isinstance(sequence, str):
            self._base_sequence = sequence
        elif isinstance(sequence, Variant):
            self._parent = sequence
        else:
            if not isinstance(self.mutations, MutationSet):
                raise ValueError('If no sequence is passed, mutations must be a MutationSet so that we can parse the sequence.')
            self._parent = self.mutations.parent

        # if mutations were a string, we still need to parse it and assign the parent
        if isinstance(mutations, str):
            mutations = MutationSet.from_string(self, mutations)

        # check that any mutations have the same parent as here
        if self.mutations is not None:
            if self.mutations.parent != self:
                raise ValueError('The parent of the mutations must be the same as the parent of the variant.')
        else:
            self.mutations = MutationSet(self)

    @property
    def parent(self) -> Variant:
        """The parent variant of the variant, if present."""
        return self._parent

    @property
    def base_sequence(self) -> str:
        """The sequence of the variant without mutations."""
        if self.parent is not None:
            return self(self.parent)
        else:
            return self._base_sequence
        
    def __str__(self):
        self.mutations.get_variant_str()
    
    def add_mutations(self, mutations: Union[MutationSet, Mutation]):
        """Add mutations to the variant.
        
        Params
        ------
        mutations : Union[MutationSet, Mutation]
            The mutations to add to the variant.
        """
        # first convert to MutationSet
        if isinstance(mutations, Mutation):
            mutations = MutationSet([mutations])
        elif isinstance(mutations, MutationSet):
            pass
        else:
            raise ValueError('mutations must be a MutationSet or Mutation')
        
        self.mutations = self.mutations.union(mutations)

    def parse_mutations(self, other: Variant, expect_indels: bool=False, **blast_params):
        raise NotImplementedError('Need to write alignment method.')
    
    def __eq__(self, other: Variant):
        return str(self) == str(other)
        
