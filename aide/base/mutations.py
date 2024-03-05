"""Base classes for representing mutations"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Union, List, Iterable, Dict, Hashable
from collections.abc import MutableSet


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