from aide import Variant, Mutation, MutationSet, VariantLabel
import unittest

class TestVariant(unittest.TestCase):
    def setUp(self):
        self.parent_seq = Variant('MAGV')

    def test_add_mutations(self):
        mutation = Mutation.from_string('A2[TM]')
        self.parent_seq.add_mutations(mutation)
        self.assertEqual(str(self.parent_seq), 'MTMGV')

    def test_add_mutation_set(self):
        mutation_set = MutationSet.from_string('A2[TM];G3[ST]')
        self.parent_seq.add_mutations(mutation_set)
        self.assertEqual(str(self.parent_seq), 'MTMSTV')

    def test_parse_mutations(self):
        variant1 = Variant('MAVG')
        self.assertEqual(str(self.parent_seq.parse_mutations(variant1)), 'G2V;V3G')

    def test_parse_mutations_indels(self):
        variant1 = Variant('MMMMAAAAAAAMMMMAV')
        variant2 = Variant('MMMMMMMMAVG')
        self.assertEqual(
            str(variant1.parse_mutations(variant2, expect_indels=True, aggregate_indels=True)), '[AAAAAAA]4-;0>17G')

    def test_parse_mutations_identical(self):
        variant3 = Variant('MAGV')
        self.assertEqual(str(self.parent_seq.parse_mutations(variant3)), '')

    def test_eq(self):
        variant1 = Variant('MAGV')
        variant2 = Variant('MAGV')
        self.assertEqual(variant1, variant2)
        variant4 = Variant('MAMV', mutations='M3G')
        self.assertEqual(variant1, variant4)

    def test_ne(self):
        variant1 = Variant('MAGV')
        variant2 = Variant('MAVG')
        self.assertNotEqual(variant1, variant2)

        variant3 = Variant('MAGV', mutations='A2[TM]')
        self.assertNotEqual(variant1, variant3)

    def test_parent(self):
        variant = Variant(self.parent_seq, mutations='A2[TM]')
        self.assertTrue(variant.parent is self.parent_seq)

    def test_cut_parent(self):
        variant = Variant(self.parent_seq, mutations='A2[TM]')
        variant.cut_parent()
        self.assertTrue(variant.parent is None)
        self.assertTrue(str(variant) == 'MTMGV')

    def test_cut_children(self):
        variant = Variant(self.parent_seq, mutations='A2[TM]')
        self.parent_seq.cut_children()
        self.assertTrue(len(self.parent_seq.children) == 0)
        self.assertTrue(variant.parent is None)

    def test_propegate(self):
        variant = self.parent_seq.propegate(mutations = 'A2[TM]', id='test')
        self.assertEqual(str(variant), 'MTMGV')
        self.assertEqual(variant.id, 'test')
        self.assertTrue(variant.parent is self.parent_seq)
        self.assertTrue(id(variant) in self.parent_seq.children)
    
    def test_immutable_parent(self):
        variant = Variant(self.parent_seq, mutations='A2[TM]')
        with self.assertRaises(ValueError):
            self.parent_seq.add_mutations(Mutation.from_string('V3G'))

    def test_hash(self):
        # variants with different ids should always evaluate as different
        # but if no ids are given, and the variants have the same sequence
        # then they should be the same
        variant1 = Variant('MAGV', id='test')
        variant2 = Variant('MAGV', id='test2')
        variant3 = Variant('MAGV', mutations='A2[TM]')
        variant4 = Variant(variant1, mutations='A2[TM]')
        self.assertNotEqual(hash(variant1), hash(variant2))
        self.assertEqual(hash(variant4), hash(variant3))

    def test_labels(self):
        variant = Variant(self.parent_seq, mutations='A2[TM]')
        variant.assign_labels(labels={'test': 1})
        self.assertEqual(variant.labels.get_labels('test'), [1])

class TestVariantLabel(unittest.TestCase):

    def test_init(self):
        label = VariantLabel(initial_data={'a': [1,2], 'b': [3,4]})
        self.assertEqual(label.signature, ['a', 'b'])

        label = VariantLabel(signature = ['a', 'b'])
        self.assertEqual(label.signature, ['a', 'b'])

        label = VariantLabel(initial_data={'a': [1,2], 'b': [3,4]}, round_idx=2)
        # round idx is stored as the second element in the tuple
        # for each value added to the label
        for labels in label.values():
            for value in labels:
                self.assertEqual(value[1], 2)

    def test_get_labels(self):
        ## This method gets the actual data and leaves the round index
        label = VariantLabel(initial_data={'a': [1,2], 'b': [3,4]})

        outs = label.get_labels('a', agg_func=None)
        self.assertEqual(outs, [1,2])

        outs = label.get_labels(['a','b'], agg_func=None)
        outs2 = label.get_labels(None, agg_func=None)
        self.assertEqual(outs, outs2)
        self.assertEqual(outs, {'a': [1,2], 'b': [3,4]})

        def mean(x):
            return sum(x)/len(x)
        outs = label.get_labels('a', agg_func=mean)
        self.assertEqual(outs, 1.5)

        outs = label.get_labels(['a','b'], agg_func=mean)
        self.assertEqual(outs, {'a': 1.5, 'b': 3.5})

    def test_add_labels(self):
        label = VariantLabel(initial_data={'a': [1,2], 'b': [3,4]})
        label.add_labels(a=5, b=6)
        self.assertEqual(label.get_labels('a'), [1,2,5])
        self.assertEqual(label.get_labels('b'), [3,4,6])

        label.add_labels(a=7, b=8, round_idx=2)
        self.assertEqual(label.get_labels('a'), [1,2,5,7])
        self.assertEqual(label.get_labels('b'), [3,4,6,8])
        self.assertEqual(label['a'][-1][1], 2)

        with self.assertRaises(ValueError):
            label.add_labels(c=9)

        label.add_labels(c=9, enforce_signature=False)
        self.assertEqual(label.get_labels('c'), [9])

        
class TestMutation(unittest.TestCase):
    def setUp(self):
        self.parent_seq = Variant('MAGV')

    def test_order(self):
        mutation = Mutation.from_string('1>2T')
        self.assertEqual(mutation.order, 1)

    def test_subst(self):
        mutation = Mutation.from_string('A2M')
        variant_seq = mutation.get_variant_str(self.parent_seq)
        self.assertEqual(variant_seq, 'MMGV')

    def test_del(self):
        mutation = Mutation.from_string('A2[-]')
        variant_seq = mutation.get_variant_str(self.parent_seq)
        self.assertEqual(variant_seq, 'MGV')

    def test_multi_del(self):
        mutation = Mutation.from_string('[AG]2[--]')
        variant_seq = mutation.get_variant_str(self.parent_seq)
        self.assertEqual(variant_seq, 'MV')

    def test_ins(self):
        mutation = Mutation.from_string('A2[AT]')
        variant_seq = mutation.get_variant_str(self.parent_seq)
        self.assertEqual(variant_seq, 'MATGV')

        mutation = Mutation.from_string('>2T')
        variant_seq = mutation.get_variant_str(self.parent_seq)
        self.assertEqual(variant_seq, 'MTAGV')

    def test_eq(self):
        mutation1 = Mutation.from_string('A2M')
        mutation2 = Mutation.from_string('A2M')
        self.assertEqual(mutation1, mutation2)

    def test_ne(self):
        mutation1 = Mutation.from_string('A2M')
        mutation2 = Mutation.from_string('A2[-]')
        self.assertNotEqual(mutation1, mutation2)


class TestMutationSet(unittest.TestCase):
    def setUp(self):
        self.parent_seq = Variant('MAGV')
        self.mutation_set = MutationSet.from_string('A2[TM];G3[ST]')

    def test_same_positions(self):
        self.assertRaises(ValueError, MutationSet.from_string, 'A2[TM];A2M')
        self.assertRaises(ValueError, MutationSet.from_string, '>2T;>2M')

        mutation_set = MutationSet.from_string('1>2T;0>2M') # this should be fine since both are insertions

    def test_get_variant_str(self):
        variant_str = self.mutation_set.get_variant_str(self.parent_seq)
        self.assertEqual(str(variant_str), 'MTMSTV')

    def test_in(self):
        mutation = Mutation.from_string('A2[TM]')
        self.assertIn(mutation, self.mutation_set)

    def test_not_in(self):
        mutation = Mutation.from_string('A2M')
        self.assertNotIn(mutation, self.mutation_set)

    def test_eq(self):
        mutation_set1 = MutationSet.from_string('A2[TM];G3[ST]')
        mutation_set2 = MutationSet.from_string('A2[TM];G3[ST]')
        self.assertEqual(mutation_set1, mutation_set2)

    def test_ne(self):
        mutation_set1 = MutationSet.from_string('A2[TM];G3[ST]')
        mutation_set2 = MutationSet.from_string('A2[TM];G3[ST];V4G')
        self.assertNotEqual(mutation_set1, mutation_set2)

    def test_union(self):
        mutation_set1 = MutationSet.from_string('A2[TM];G3[ST]')
        mutation_set2 = MutationSet.from_string('A2[TM];G3[ST];V4G')
        mutation_set3 = mutation_set1 | mutation_set2
        self.assertEqual(mutation_set3, MutationSet.from_string('A2[TM];G3[ST];V4G'))

    def test_intersection(self):
        mutation_set1 = MutationSet.from_string('A2[TM];G3[ST]')
        mutation_set2 = MutationSet.from_string('A2[TM];G3[ST];V4G')
        mutation_set3 = mutation_set1 & mutation_set2
        self.assertEqual(mutation_set3, MutationSet.from_string('A2[TM];G3[ST]'))

    def test_difference(self):
        mutation_set1 = MutationSet.from_string('A2[TM];G3[ST]')
        mutation_set2 = MutationSet.from_string('A2[TM];G3[ST];V4G')
        mutation_set3 = mutation_set2 - mutation_set1
        self.assertEqual(mutation_set3, MutationSet.from_string('V4G'))

    def test_get_variant_str(self):
        variant_str = self.mutation_set.get_variant_str(self.parent_seq)
        self.assertEqual(str(variant_str), 'MTMSTV')

    def test_get_variant_str_hard(self):
        parent = Variant('MAVGGVMDEYLMGY')
        mutation_set = MutationSet.from_string(
            '[AVG]2-;G5M;D8[DYG];G13V;Y14-;>1G;1>13I;0>13J;>15I')
        variant_str = mutation_set.get_variant_str(parent)
        self.assertEqual(
            variant_str,
            'GMMVMDYGEYLMJIVI')


if __name__ == '__main__':
    unittest.main()
