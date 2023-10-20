"""Base classes for dealing with protein sequences.
"""
from __future__ import annotations
import re
from typing import Union, List, Iterable, Set, Dict
from dataclasses import dataclass
from collections.abc import MutableSet, Hashable
from collections import UserDict

from aide.utils.alignment import biopython_align

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
    label : Union[float, None, Dict[float]]
        The label of the variant.
    round_putative : Union[int, None]
        The most recent round the variant was put up as a putative variant.
    round_experiment : Union[int, None]
        The round the variant selected for experiment.

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
    labels : float or dict of float
        The label of the variant.
    round_added : int
        The round the variant was added to the database.
    round_putative : int
        The most recent round the variant was put up as a putative variant.
    round_experiment : int
        The round the variant selected for experiment.

    Methods
    -------
    cut_children()
        Remove the reference to all children.
    cut_parent()
        Remove the reference to the parent.
    set_label(label: float)
        Set the label of the variant.
    parse_mutations(other: Variant, expect_indels: bool=False, aggregate_indels: bool=False, **blast_params)
        Parse the mutations self and another variants. Be default, it is assumed there are no indels and the sequences are the same length.
        If expect_indels is True, the sequences are aligned and indels are parsed. This is not a perfect process, so it is recommended to
        keep track of mutations manually instead of try to parse them.
    is_descendant_of(other: Variant, only_parent: bool=False)
        Whether the variant is a descendant of another variant.
    is_ancestor_of(other: Variant, only_parent: bool=False)
        Whether the variant is an ancestor of another variant.
    propegate(mutations: Union[MutationSet, Mutation, str, None] = None, id: str = None)
        Create a new variant with the called variant as the parent.

    """
    def __init__(
        self,
        sequence: Union[None, str, Variant] = None,
        mutations: Union[None, str, MutationSet, Mutation] = None,
        id: Union[None, str] = None,
        mutatable: bool = True,
        labels: VariantLabel = VariantLabel(),
        round_added: int=None,
        round_putative: Union[int, List[int]]=None,
        round_experiment: Union[int, List[int]]=None
    ):
        self._base_sequence = None
        self._parent = None
        self._id = id
        self._mutatable = mutatable
        self._children = {}
        if not isinstance(labels, VariantLabel):
            raise ValueError('label must be a VariantLabel')
        self._labels = labels
        self._round_added = None
        self._round_putative = []
        self._round_experiment = []

        self.round_added = round_added
        if round_putative is not None:
            if type(round_putative) == int:
                round_putative = [round_putative]
            self.round_putative = round_putative
        if round_experiment is not None:
            if type(round_experiment) == int:
                round_experiment = [round_experiment]
            self.round_experiment = round_experiment
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
    def labels(self) -> VariantLabel:
        return self._label

    @property
    def round_added(self) -> int:
        return self._round_added
    
    @round_added.setter
    def round_added(self, round_added: int):
        self._round_added = int(round_added)

    @property
    def round_putative(self) -> List[int]:
        return self._round_putative
    
    def add_round_putative(self, round_putative: int):
        if round_putative not in self.round_putative:
            self.round_putative.append(round_putative)

    @round_putative.setter
    def round_putative(self, round_putative: List[int]):
        self._round_putative = round_putative
        
    @property
    def round_experiment(self) -> List[int]:
        return self._round_experiment
    
    @round_experiment.setter
    def round_experiment(self, round_experiment: List[int]):
        self._round_experiment = round_experiment

    def add_round_experiment(self, round_experiment: int):
        if round_experiment not in self.round_experiment:
            self.round_experiment.append(round_experiment)

    def set_labels(self, labels: Union[dict, VariantLabel], enforce_signature: bool=True, round_idx: int=None):
        """Set the label of the variant.

        Params
        ------
        labels : Union[dict, VariantLabel]
            The labels to set.
        enforce_signature : bool
            Whether to enforce the signature of the label.
        round_idx : int
            The round the labels were measured.
        """
        if isinstance(labels, VariantLabel):
            labels = labels.data
        self._labels.add_labels(enforce_signature=enforce_signature, **labels)
    
    def __repr__(self):
        return f'Variant(id={self.id}, label={self.label})'

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
            return hash(str(self))
    
    def __eq__(self, other: Variant):
        return str(self) == str(other)
    
    def __neq__(self, other: Variant):
        return str(self) != str(other)
    
    def __contains__(self, mutation: Mutation):
        return mutation in self.mutations
    
    def __len__(self):
        return len(str(self))
    
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

    def parse_mutations(self, other: Variant, expect_indels: bool=False, aggregate_indels: bool=False, **blast_params):
        """Parse the mutations between two variants.
        
        Params
        ------
        other : Variant
            The other variant to parse.
        expect_indels : bool
            Whether to expect indels or not. If not, sequences must be same length and are compared position wise.
        aggregate_indels : bool
            Whether to aggregate any indels into a single event. Eg. if there is a deletion of 2 AAs, next to eachother
            in the sequence, the deletion will be aggregated into a single event.
        """
        if expect_indels:
            query, subject, score, begin, end = biopython_align(str(self), str(other), **blast_params)
        else:
            if len(self) != len(other):
                raise ValueError('The sequences must be the same length if expect_indels is False.')
            else:
                query = str(self)
                subject = str(other)

        mutations = []
        pos1 = 0  # Position counter for seq1
        order = 0 # Order counter for insertions

        print(query, subject)
        
        for s1, s2 in zip(query, subject):
            mutation_str = ""
            # always reset order if not an insertion
            if s1 != '-':
                order = 0

            if s1 != '-' and s2 != '-':
                if s1 != s2:
                    mutation_str = f"{s1}{pos1}{s2}"
                    pos1 += 1
                else:
                    pos1 += 1
            elif s1 == '-':
                mutation_str = f"{order}>{pos1}{s2}"
                order += 1
            elif s2 == '-':
                mutation_str = f"{s1}{pos1}-"
                pos1 += 1

            if mutation_str:
                mutation = Mutation.from_string(mutation_str, zero_indexed=True)
                mutations.append(mutation)

        mutationset = MutationSet(mutations, aggregate_indels=aggregate_indels)
        return mutationset
                
    
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


class VariantLabel(UserDict):
    """Stores labels for variants.
    
    Variants can have multiple measurements for multiple rounds. This class
    stores the labels for each round.

    Params
    ------
    initial_data : Dict[int, float]
        The initial labels.
    signature : Iterable[str]
        The signature of the labels. If None, the signature is inferred from the initial data.
    round_idx : int
        The round the labels were measured.

    """
    def __init__(self, initial_data: Dict[int, float]=None, signature: Iterable[str], round_idx: int=None):
        super().__init__()

        if signature is not None:
            self.signature = signature

        if initial_data is not None:
            self.add_labels(round_idx=round_idx, enforce_signature=True, **initial_data)

    @property
    def signature(self):
        return list(self.data.keys())
    
    @signature.setter
    def signature(self, signature: Iterable[str]):
        if len(self) > 0:
            raise ValueError('Cannot change the signature of a VariantLabel with data.')
        self.data = {label_name: [] for label_name in signature}

    def add_labels(self, round_idx: int=None, enforce_signature: bool=True, **kwargs):
        """Add labels to the variant.
        
        Params
        ------
        round_idx : int
            The round the labels were measured.
        enforce_signature : bool
            Whether to enforce the signature of the label.
        kwargs
            The labels to add.
        """
        if enforce_signature:
            if len(self.signature) == 0:
                # if we have no signature, no need to enfore, just set it
                self.signature = kwargs.keys()
            else:
                for label_name in kwargs:
                    if label_name not in self.signature:
                        raise ValueError(f'Label name {label_name} not in signature.')
        for label_name, label_value in kwargs.items():
            if label_name not in self.data:
                self.data[label_name] = []
            self.data[label_name].append((label_value, round_idx))

    def get_labels(self, label_name: Union[str, Iterable[str], None]=None, agg_func: callable=None):
        """Return label values.
        
        Params
        ------
        label_name : str
            The label name to return. If None, return all labels.
        agg_func : callable
            The function to aggregate the labels. If None, return the raw labels.
        """
        # need to extract only labels, not rounds
        def _data_to_labels_no_round(data):
            return [label[0] for label in data]

        if type(label_name) == str:
            label_names = [label_name]
        elif label_name is None:
            label_names = self.signature
        else:
            label_names = label_name

        if agg_func is None:
            agg_func = lambda x: x

        if len(label_names) == 1:
            return agg_func(_data_to_labels_no_round(self.data[label_names[0]]))
        else:
            return {label_name: agg_func(_data_to_labels_no_round(self.data[label_name])) for label_name in label_names}
        
    
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
    order : int, default None
        For insertions, specify the order of the insertions when other insertions occur at the same position.
    """
    position: int
    ref: str
    alt: str
    order: int=0

    def __repr__(self):
        return f'Mutation({str(self)})'

    @property
    def initial_width(self):
        """The width of the mutation before accounting for indels."""
        return len(self.ref)
    
    @property
    def is_insertion(self):
        """Whether the mutation is an insertion."""
        return self.ref == ''

    @property
    def is_deletion(self):
        """Whether the mutation is a deletion."""
        return self.alt == '-'

    def _check_validity(self, variant: Variant):
        if self.is_insertion and self.position > len(variant):
            raise ValueError('Insertions are applied before the position, but the position is greater than the length of the sequence.')

        elif variant.base_sequence[self.position:self.position + self.initial_width] != self.ref:
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
        # >2V -> '', 2, V (insertion)
        # 0>2V -> '', 2, V (insertion with order, eg. this insertion is to the left of 1>2M)
        pattern = r'(?P<order>\d+>)?(?P<ref>\[?[A-Z>]+\]?)?(?P<pos>\d+)(?P<alt>\[?[A-Z-]+\]?)'
        match = re.match(pattern, string)

        if match:
            # Extracting matched groups
            order = match.group('order')
            if order:
                order = int(order[:-1])  # Remove '>' and convert to int
            ref = match.group('ref')
            if ref:
                ref = ref.replace('[', '').replace(']', '')
            else:
                ref = ''
            pos = int(match.group('pos'))
            alt = match.group('alt').replace('[', '').replace(']', '')

            if ref == '>':
                ref = ''
            if order is None:
                order = 0

            return ref, pos, alt, order
        else:
            raise ValueError(f"Invalid mutation string: {string}")


    @classmethod
    def from_string(cls, mutation: str, zero_indexed: bool=False):
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
        ref, pos, alt, order = cls._parse_mutation_string(mutation)
        if not zero_indexed:
            pos -= 1

        if not pos >= 0:
            raise ValueError('The position must be greater than or equal to 0.')
        return cls(pos, ref, alt, order)
    
    def __str__(self):
        if len(self.ref) == 1:
            ref = self.ref
        elif self.ref == '':
            ref = f'{self.order}>'
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
            return list(str(variant))
        else:
            # must keep the width together and append gaps
            # eg. for sequence AMV, and mutation [AM]0V -> ['AM', '-', 'V']
            output = []
            last_mutation_excess = None
            for i, aa in enumerate(variant.base_sequence):
                if i == self.position:
                    if self.initial_width > 1:
                        output.append(self.ref)
                    else:
                        output.append(aa)
                    output.extend(['-'] * (self.initial_width - 1))
                    last_mutation_excess = self.initial_width - 1
                elif last_mutation_excess is not None and last_mutation_excess > 0:
                    # here, we appended a mutation that was bigger than 1, so we need to
                    # not double count AAs in the sequence
                    last_mutation_excess -= 1
                else:
                    output.append(aa)
            return output  

    def _update_variant_list(self, aa_list: List[str]):
        """Update a list of amino acids by the mutation.
        """
        if not self.is_insertion:
            aa_list[self.position] = self.alt
        else:
            # check if we are at the end of the sequence
            if self.position == len(aa_list):
                aa_list.append(self.alt)
            else:
                # insert before what is there
                aa_list[self.position] = self.alt + aa_list[self.position]

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

    def __init__(self, mutations: Iterable[Mutation]=None, aggregate_indels: bool=False):
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

        if aggregate_indels:
            self._aggregate_indels()  

    def __repr__(self):
        return f'MutationSet({self.mutations})'
    
    def __str__(self):
        mutations = sorted(list(self.mutations), key=lambda x: (x.position, x.order))
        return ';'.join([str(mutation) for mutation in mutations])

    @staticmethod
    def _check_unique_mutation_positions(mutations: Iterable[Mutation]):
        positions = [mutation.position for mutation in mutations if not mutation.is_insertion]
        insertion_positions = [mutation.position for mutation in mutations if mutation.is_insertion]
        # some mutations are multiple wide, so we must add those positions
        for mutation in mutations:
            if mutation.initial_width > 1:
                positions.extend(range(mutation.position + 1, mutation.position + mutation.initial_width))
                # also check if no insertion lives inside here
                for i in range(mutation.position + 1, mutation.position + mutation.initial_width):
                    if i in insertion_positions:
                        raise ValueError('Cannot have an insertion inside a multi AA mutation.')

        if not len(set(positions)) == len(positions):
            raise ValueError('The mutation positions must be unique, but some are applied to the same position.')
        
        # now look at insertions
        insertion_mutations = {}
        for mutation in mutations:
            if mutation.is_insertion:
                if mutation.position not in insertion_mutations:
                    insertion_mutations[mutation.position] = []
                insertion_mutations[mutation.position].append(mutation)
        # if any insertions occur at the same position, they must have an order
        for _, insertion_list in insertion_mutations.items():
            if len(insertion_list) > 1:
                orders = [m.order for m in insertion_list]
                if len(set(orders)) != len(orders):
                    raise ValueError('If multiple insertions occur at the same position, they must have a unique order.')

    def _aggregate_indels(self):
        # if insertions are next to insertions or deletions are next to deletions, we can aggregate them
        # into single larger events
        sorted_by_position = sorted(list(self.mutations), key=lambda x: (x.position, x.order))
        print(sorted_by_position)
        current_aggregate = None
        aggregated_mutations = []

        for mutation in sorted_by_position:
            if current_aggregate is None:
                current_aggregate = mutation
                continue

            # check if we are adjacent and of the same type
            if mutation.position == current_aggregate.position:
                # if this is true and we passed the init checks, then these are both insertions
                # we are on the correct order, so we can jsut append the insertion
                current_aggregate = Mutation(ref=current_aggregate.ref, alt=current_aggregate.alt + mutation.alt, position=current_aggregate.position, order=current_aggregate.order)
            elif mutation.position == current_aggregate.position + current_aggregate.initial_width:
                # these are adjacent, now check if they are both deletions
                if current_aggregate.alt.endswith('-') and mutation.is_deletion:
                    current_aggregate = Mutation(ref=current_aggregate.ref + mutation.ref, alt=current_aggregate.alt, position=current_aggregate.position, order=current_aggregate.order)
                    # this should update initial with too so that we catch the next deletion
                else:
                    # not a deletion, so we can't aggregate
                    aggregated_mutations.append(current_aggregate)
                    current_aggregate = mutation
            else:
                # not adjacent, so we can't aggregate
                aggregated_mutations.append(current_aggregate)
                current_aggregate = mutation
        if current_aggregate is not None:
            aggregated_mutations.append(current_aggregate)
        self.mutations = set(aggregated_mutations)

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
        last_mutation_excess = None
        for i, aa in enumerate(variant.base_sequence):
            if i in position_map:
                mutation = position_map[i]
                if mutation.initial_width > 1:
                    output.append(mutation.ref)
                else:
                    output.append(aa)
                output.extend(['-'] * (mutation.initial_width - 1))
                last_mutation_excess = mutation.initial_width - 1
            elif last_mutation_excess is not None and last_mutation_excess > 0:
                # here, we appended a mutation that was bigger than 1, so we need to
                # not double count AAs in the sequence
                last_mutation_excess -= 1
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
        # start with non insertions
        for mutation in self.mutations:
            if not mutation.is_insertion:
                mutation._update_variant_list(aa_list)

        # now do insertions
        # mutation update list will add the mutation by looking at what is currently in
        # the position and PREprending the mutation
        # thus if multiple insertions happen at the same position, we need to go in reverse order
        # we should have already checked that there are orders for these mutations
        insertion_mutations = {}
        for mutation in self.mutations:
            if mutation.is_insertion:
                if mutation.position not in insertion_mutations:
                    insertion_mutations[mutation.position] = []
                insertion_mutations[mutation.position].append(mutation)
        for _, mutations in insertion_mutations.items():
            mutations = sorted(mutations, key=lambda x: x.order)
            for mutation in mutations[::-1]:
                mutation._update_variant_list(aa_list)

        output_seq = ''.join(aa_list)
        output_seq = Mutation._remove_gaps(output_seq)
        return output_seq
        

    



        

        
        
        
        
