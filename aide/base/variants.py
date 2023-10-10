"""Base classes for dealing with protein sequences.
"""
from __future__ import annotations
import re
from typing import Union, List, Iterable
from dataclasses import dataclass
from collections.abc import MutableSet, Hashable

class Variant:
    """A protein sequence with a unique identifier.
    
    Params
    ------
    sequence : Union[str, Variant]
        The sequence of the variant. If a Variant is passed, the sequence
        is interpreted as the parent sequence of the variant.
    mutations : Union[str, MutationSet, Mutation, None]
        The mutations that define the variant. If a string is passed, it
        will be parsed into a MutationSet.
    id : Union[str, None]
        The unique identifier of the variant. If None is passed, the id
        is computed as a hash of the sequence and mutations.

    Attributes
    ----------
    base_sequence : str
        The sequence of the variant without mutations.
    parent : Variant or None
        The parent variant of the variant, if present
    mutations : MutationSet
        The mutations that define the variant.
    hash : int
        The hash of the variant. Determined by post mutation sequence.
    id : str
        The unique identifier of the variant. If none given, the id is
        computed as a hash of the sequence and mutations.
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
            self.mutations = MutationSet.from_string(self, mutations)
        else:
            raise ValueError('mutations must be a MutationSet, Mutation, or str')
        
        # If we have just a sequence, we have our base sequence
        # if a parent was passed, the base sequnce comes from that
        if isinstance(sequence, str):
            self._base_sequence = sequence
        elif isinstance(sequence, Variant):
            self._parent = sequence
        else:
            raise ValueError('sequence must be a str or Variant')

        # check that any mutations are valid for this sequnce
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
    
@dataclass(frozen=True, eq=True)
class Mutation:
    """A single mutation.
    
    Params
    ------
    position : int
        The position of the mutation.
    ref : str
        The reference amino acid.
    alt : str
        The alternate amino acid.
    """
    position: int
    ref: str
    alt: str

    def __repr__(self):
        return f'Mutation({str(self)})'

    @property
    def initial_width(self):
        """The width of the mutation before accounting for indels."""
        return len(self.ref)

    def _check_validity(self, variant: Variant):
        if variant.base_sequence[self.position:self.position + self.initial_width] != self.ref:
            raise ValueError('The reference amino acid does not match the parent sequence.')
    
    @staticmethod
    def _parse_mutation_string(string):
        """Parse a mutation string into a mutation.
        
        Params
        ------
        string : str
            The mutation string in the form of '<ref><pos: int><alt>'.
        
        Returns
        -------
        Mutation
            The mutation object.
        """
        # need regular expressions to parse the mutation
        # examples of ref, pos, and alt:
        # A2M -> A, 2, M
        # A2[MV] -> A, 2, MV
        # A2- -> A, 2, -
        # [AV]2[--] -> AV, 2, --
        pattern = r'(?P<ref>\[?[A-Z]+\]?)?(?P<pos>\d+)(?P<alt>\[?[A-Z-]+\]?)'
        match = re.match(pattern, string)
        
        if match:
            # Extracting matched groups
            ref = match.group('ref').replace('[', '').replace(']', '')
            pos = int(match.group('pos'))
            alt = match.group('alt').replace('[', '').replace(']', '')
        else:
            raise ValueError(f"Invalid mutation string: {string}")

        return ref, pos, alt

    @classmethod
    def from_string(cls, mutation: str, zero_indexed: bool=True):
        """Create a mutation from a string.
        
        Params
        ------
        mutation : str
            The mutation string in the form of '<ref><pos: int><alt>'.
        zero_indexed : bool
            Whether the position is zero indexed or not.
        
        Returns
        -------
        Mutation
            The mutation object.
        """
        ref, pos, alt = cls._parse_mutation_string(mutation)
        if not zero_indexed:
            pos -= 1
        return cls(pos, ref, alt)
    
    def __str__(self):
        if len(self.ref) == 1:
            ref = self.ref
        else:
            ref = f'[{self.ref}]'

        if len(self.alt) == 1:
            alt = self.alt
        else:
            alt = f'[{self.alt}]'
        return f'{ref}{self.position}{alt}'

    @staticmethod
    def _remove_gaps(sequence):
        """Remove gaps from a sequence.
        
        Params
        ------
        sequence : str
            The sequence to remove gaps from.
        
        Returns
        -------
        str
            The sequence without gaps.
        """
        return sequence.replace('-', '')

    def _get_variant_base_list(self, variant: Variant):
        """Convert the sequence of the variant into a list with the correct positions.
        
        This is a little bid hard when the reference in the mutation is multiple amino acids.
        Eg. [AM]1V. We need to keep AM together to replace it, but also we don't want to shift
        the other positions. So we need to add 'gaps' to the sequence to keep the positions
        """
        if self.initial_width == 1:
            return list(variant)
        else:
            # must keep the width together and append gaps
            # eg. for sequence AMV, and mutation [AM]0V -> ['AM', '-', 'V']
            output = []
            for i, aa in enumerate(variant.base_sequence):
                if i == self.position:
                    output.append(self.ref)
                    output.extend(['-'] * (self.initial_width - 1))
                else:
                    output.append(aa)
            return output  

    def _update_variant_list(self, aa_list: List[str]):
        """Update a list of amino acids by the mutation.
        """
        aa_list[self.position] = self.alt

    def get_variant_str(self, variant: Variant):
        """Parse the protein sequence after the mutation has been applied.

        Params
        ------
        variant : Variant
            The variant to parse.
        """
        # check the variant against the sequence
        self._check_validity(variant)

        aa_list = self._get_variant_base_list(variant)
        self._update_variant_list(aa_list)
        output_seq = ''.join(aa_list)
        output_seq = self._remove_gaps(output_seq)
        return output_seq
    
class MutationSet(MutableSet, Hashable):
    """A set of mutations.
    
    Params
    ------
    mutations : List[Mutation]
        The list of mutations.
    """
    __hash__ = MutableSet._hash

    def __init__(self, mutations: Iterable[Mutation]):
        mutations = list(mutations)
        # check the correct type
        for mutation in mutations:
            if not isinstance(mutation, Mutation):
                raise ValueError('mutations must be a list of Mutation objects.')
        mutation_set = set(mutations)
        # check than non of the positions overlap
        self._check_unique_mutation_positions(mutation_set)
        self.mutations = mutation_set

    def __repr__(self):
        return f'MutationSet({self.mutations})'

    @staticmethod
    def _check_unique_mutation_positions(mutations: Iterable[Mutation]):
        positions = [mutation.position for mutation in mutations]
        # some mutations are multiple wide, so we must add those positions
        for mutation in mutations:
            if mutation.initial_width > 1:
                positions.extend(range(mutation.position + 1, mutation.position + mutation.initial_width))

        if not len(set(positions)) == len(positions):
            raise ValueError('The mutation positions must be unique, but some are applied to the same position.')

    # MutableSet methods
    def __contains__(self, mutation: Mutation):
        return mutation in self.mutations
    
    def __iter__(self):
        return iter(self.mutations)
    
    def __len__(self):
        return len(self.mutations)
    
    def add(self, mutation: Mutation):
        putative_set = self.mutations | {mutation}
        self._check_unique_mutation_positions(putative_set)

        self.mutations.add(mutation)

    def discard(self, mutation: Mutation):
        self.mutations.discard(mutation)

    @property
    def position_map(self):
        """A map of positions in a sequence to mutations."""
        return {m.position: m for m in self.mutations}

    @classmethod
    def from_string(cls, mutations: Union[str, List[str]], zero_indexed: bool=False):
        """Create a mutation set from a string or list of strings.
        
        Params
        ------
        mutations : Union[str, List[str]]
            The mutation string or list of mutation strings.
            If a single string, the mutations are separated by semicolons.
        zero_indexed : bool
            Whether the position is zero indexed or not.
        """
        if isinstance(mutations, str):
            mutations = mutations.split(';')
        return cls([Mutation.from_string(mutation, zero_indexed) for mutation in mutations])
    
    def _check_validity(self, variant: Variant):
        """Check that the mutations are valid for the variant.
        
        Params
        ------
        variant : Variant
            The variant to check.
        """
        for mutation in self.mutations:
            mutation._check_validity(variant)

    def _get_variant_base_list(self, variant: Variant):
        """Convert the sequence of the variant into a list with the correct positions.
        
        This is a little bid hard when the reference in a mutation is multiple amino acids.
        Eg. [AM]1V. We need to keep AM together to replace it, but also we don't want to shift
        the other positions. So we need to add 'gaps' to the sequence to keep the positions
        """
        # we must loop through each mutation and keep multiple AA references together
        # eg. for sequence AMV, and mutation [AM]0V -> ['AM', '-', 'V']
        position_map = self.position_map
        output = []
        for i, aa in enumerate(variant.base_sequence):
            if i in position_map:
                mutation = position_map[i]
                output.append(mutation.ref)
                output.extend(['-'] * (mutation.initial_width - 1))
            else:
                output.append(aa)
        return output

    def get_variant_str(self, variant: Variant):
        """Parse the protein sequence after the mutation has been applied.

        Params
        ------
        variant : Variant
            The variant to parse.
        """
        self._check_validity(variant)
        
        aa_list = self._get_variant_base_list(variant)
        for mutation in self.mutations:
            mutation._update_variant_list(aa_list)
        output_seq = ''.join(aa_list)
        output_seq = Mutation._remove_gaps(output_seq)
        return output_seq
        

    



        

        
        
        
        
