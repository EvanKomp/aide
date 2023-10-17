"""Base classes for dealing with protein sequences.
"""
from __future__ import annotations
import re
from typing import Union, List, Iterable, Set
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
        The parent variant of the variant, if present. Immutable. 
    children : Set[Variant]
        The children variants of the variant. Immutable.
    mutations : MutationSet
        The mutations that define the variant.
    hash : int
        The hash of the variant. Determined by post mutation sequence.
    id : str
        The unique identifier of the variant. If none given, the id is
        computed as a hash of the sequence and mutations.
    mutatable : bool
        Whether the variant is mutatable. If False, the variant is
        immutable and cannot be mutated. Note that when a variant becomes
        a parent of another variant, it becomes immutaatable.
    """
    def __init__(
        self,
        sequence: Union[None, str, Variant] = None,
        mutations: Union[None, str, MutationSet, Mutation] = None,
        id: Union[None, str] = None,
        mutatable: bool = True,
        label: Union[None, float] = None,
    ):
        self._base_sequence = None
        self._parent = None
        self._id = id
        self._mutatable = mutatable
        self._children = {}
        self._label = label
        self.mutations = MutationSet()

        # we are going to keep track of a hidden variable that is a hash of the
        # mutations assigned to the variant. Since the string/output sequence of
        # a variant is defined dynamically based on the parent, a long line of
        # parantage could lead to suboptimal repeat computation. By keeping track
        # the mutation hash since the last time we called __str__, we can can check
        # if the mutations have changed and only recompute if they have. See __str__
        # below for more details.
        self._str = None
        self._mutation_hash_last_str_call = None

        # If we have just a sequence, we have our base sequence
        # if a parent was passed, the base sequnce comes from that
        if isinstance(sequence, str):
            self._base_sequence = sequence
        elif isinstance(sequence, Variant):
            self._parent = sequence
            sequence._add_child(self)
        else:
            raise ValueError('sequence must be a str or Variant')

        # any mutations passed need to be added
        if isinstance(mutations, Mutation):
            to_add = MutationSet([mutations])
        elif mutations is None:
            to_add = MutationSet()
        elif isinstance(mutations, MutationSet):
            to_add = mutations
        elif type(mutations) == str:
            to_add = MutationSet.from_string(mutations)
        else:
            raise ValueError('mutations must be a MutationSet, Mutation, or str')
        self.add_mutations(to_add)

    @property
    def parent(self) -> Variant:
        """The parent variant of the variant, if present."""
        return self._parent

    @property
    def base_sequence(self) -> str:
        """The sequence of the variant without mutations."""
        if self._base_sequence is not None:
            return self._base_sequence
        else:
            return str(self.parent)

    @property
    def id(self) -> str:
        """The unique identifier of the variant."""
        if self._id is None:
            return hash(self)
        else:
            return self._id

    @property
    def children(self) -> Set[Variant]:
        """The children variants of the variant."""
        return self._children
    
    def _add_child(self, child: Variant):
        """Add a child to the variant.
        
        Params
        ------
        child : Variant
            The child variant to add.
        """
        self._children[id(child)] = child

    def _remove_child(self, child: Variant):
        self._children.pop(id(child))

    def cut_children(self):
        """Remove the reference to all children.
        """
        for child in list(self.children.values()):
            child.cut_parent()

    def cut_parent(self):
        if self.parent is not None:
            parent_sequence = str(self.parent)
            self.parent._remove_child(self)
            self._parent = None
            self._base_sequence = parent_sequence
        else:
            raise ValueError('Cannot cut parent of a variant without a parent.')

    @property
    def mutatable(self) -> bool:
        return len(self.children) == 0 and self._mutatable
    
    @property
    def label(self) -> float:
        return self._label
    
    @label.setter
    def label(self, label: float):
        label = float(label)
        self._label = label
    
    def __repr__(self):
        return f'Variant(id={self.id})'

    def __str__(self):
        # if we have never called this before we need to bite the bullet
        if self._str is None:
            recompute = True
        # otherwise, we only need to recompute if mutations have changed
        else:
            if self._mutation_hash_last_str_call != hash(self.mutations):
                recompute = True
            else:
                recompute = False

        if recompute:
            self._str = self.mutations.get_variant_str(self)
            self._mutation_hash_last_str_call = hash(self.mutations)
        return self._str
    
    def __hash__(self):
        if self._id is not None:
            return hash(self._id)
        else:
            return hash((self.base_sequence, self.mutations))
    
    def __eq__(self, other: Variant):
        return str(self) == str(other)
    
    def __neq__(self, other: Variant):
        return str(self) != str(other)
    
    def __contains__(self, mutation: Mutation):
        return mutation in self.mutations
    
    def add_mutations(self, mutations: Union[MutationSet, Mutation]):
        """Add mutations to the variant.
        
        Params
        ------
        mutations : Union[MutationSet, Mutation]
            The mutations to add to the variant.
        """
        # check that the sequence is mutatable
        # if it has children, messing up with the sequence by adding mutations
        # will mess up the children. Does not make sense anyway, this sequence is set/
        if not self.mutatable:
            raise ValueError('Cannot mutate a variant with children.')

        # first convert to MutationSet
        if isinstance(mutations, Mutation):
            mutations = MutationSet([mutations])
        elif isinstance(mutations, MutationSet):
            pass
        else:
            raise ValueError('mutations must be a MutationSet or Mutation')
        
        new_mutations = self.mutations | mutations
        new_mutations._check_validity(self)
        self.mutations = new_mutations

    def parse_mutations(self, other: Variant, expect_indels: bool=False, **blast_params):
        raise NotImplementedError('Need to write alignment method.')
    
    def is_descendant_of(self, other: Variant, only_parent: bool=False):
        """Whether the variant is a descendant of another variant.
        
        Params
        ------
        other : Variant
            The other variant.
        only_parent : bool
            Whether to only check the parent variant.
        
        Returns
        -------
        bool
            Whether the variant is a descendant of the other variant.
        """
        if only_parent:
            return self.parent == other
        else:
            current = self
            while current.parent is not None:
                if current.parent == other:
                    return True
                current = current.parent
            return False
        
    def is_ancestor_of(self, other: Variant, only_parent: bool=False):
        """Whether the variant is an ancestor of another variant.
        
        Params
        ------
        other : Variant
            The other variant.
        only_parent : bool
            Whether to only check the parent variant.
        
        Returns
        -------
        bool
            Whether the variant is an ancestor of the other variant.
        """
        return other.is_descendant_of(self, only_parent=only_parent)
    
    def propegate(self, mutations: Union[MutationSet, Mutation, str, None] = None, id: str = None):
        """Create a new variant with the called variant as the parent.
        
        Params
        ------
        mutations : Union[MutationSet, Mutation]
            The mutations to propegate.
        """
        child = Variant(self, mutations=mutations, id=id)
        return child
    
    
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

    def __init__(self, mutations: Iterable[Mutation]=None):
        if mutations is None:
            mutations = []
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
        if variant.mutations == self:
            # in this case we ran the check wehn we assigned it already, no need to
            # do it again
            pass
        else:
            self._check_validity(variant)
        
        aa_list = self._get_variant_base_list(variant)
        for mutation in self.mutations:
            mutation._update_variant_list(aa_list)
        output_seq = ''.join(aa_list)
        output_seq = Mutation._remove_gaps(output_seq)
        return output_seq
        

    



        

        
        
        
        
