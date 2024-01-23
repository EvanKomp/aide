"""Base classes for representing a protein sequence.
"""
from __future__ import annotations
import hashlib
from typing import Union, List, Iterable, Set, Dict


from aide.utils.alignment import biopython_align
from aide.base.labels import VariantLabel, VariantLabels
from aide.base.mutations import Mutation, MutationSet

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
    labels : dict=None
        The labels of the variant. Must be of the form 
        {
            'name': Iterable of label names,
            'values': Iterable of label values,
            'round_idx': Iterable of round indices or int
        }
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
    labels : VariantLabels
        The labels of the variant.
    is_labelled : bool
        Whether the variant is labelled.
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
    add_label(self, name: str, value: float, round_idx: int=None)
        Add a label to the variant.
    add_labels(self, names: Iterable[str], values: Iterable[float], round_idx: int=None)
        Add multiple labels to the variant.
    remove_label(self, name: str, round_idx: int=None)
        Remove a label name from the variant.
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
        labels: Dict = None,
        round_added: int=None,
        round_putative: Union[int, List[int]]=None,
        round_experiment: Union[int, List[int]]=None
    ):
        self._base_sequence = None
        self._parent = None
        self._id = id
        self._mutatable = mutatable
        self._children = {}
        self._round_added = None
        self._round_putative = []
        self._round_experiment = []

        if round_added is not None:
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

        self._labels = VariantLabels()
        if labels is None:
            pass
        else:
            self.add_labels(**labels)

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
            return str(hash(self))
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
        return self._labels
    
    def has_labels(self, name: Union[str, Iterable[str]]=None, round_idx: int=None) -> bool:
        """Whether the variant has a label.
        
        Params
        ------
        name : Union[str, Iterable[str]]=None
            The name/s of the labels that the variant must have at least one of.
            if None, the variant must have at least one label of any name
        round_idx : int
            The round the label was measured.
        
        Returns
        -------
        bool
            Whether the variant has the label.
        """
        return self.labels.has_labels(name, round_idx=round_idx)

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

    def add_labels(self, names: Iterable[str], values: Iterable[float], round_idx: int=None):
        """Add multiple labels to the variant.
        
        Params
        ------
        names : Union[str, Iterable[str]]
            The names of the labels.
        values : Union[float, Iterable[float]]
            The values of the labels.
        round_idx : Union[int, List[int]]=None
            The round the label was measured.
        """
        if type(names) == str:
            names = [names]
        if not hasattr(values, '__iter__'):
            values = [values]
        if not hasattr(round_idx, '__iter__'):
            round_idx = [round_idx] * len(names)
        labels = [VariantLabel(variant_id=self.id, name=name, value=value, round_idx=rix) for name, value, rix in zip(names, values, round_idx)]
        for label in labels:
            if not label.variant_id == self.id:
                raise ValueError('Cannot add label with different variant_id.')
        self.labels.add_labels(labels)

    def get_label_df(self):
        """Get a dataframe of the labels.
        """
        return self.labels.df

    def remove_labels(self, names: Union[str, Iterable[str]], round_idx: int=None):
        """Remove multiple labels from the variant.
        
        Params
        ------
        names : Union[str, Iterable[str]]
            The names of the labels.
        round_idx : int
            The round the label was measured.
        """
        if type(names) == str:
            names = [names]
        # first we need to get the labels that match the query
        to_remove = self.labels.select(name=names, round_idx=round_idx)
        # then we need to remove them
        self.labels.remove_labels(to_remove)
                   
    def __repr__(self):
        return f'Variant(id={self.id}, label={self.labels})'

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
        m = hashlib.md5()
        if self._id is not None:
            string = self._id
        else:
            string = str(self)
        m.update(string.encode('utf-8'))
        return int(m.hexdigest(), 16)
    
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

    
        

    



        

        
        
        
        
