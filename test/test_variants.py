from aide import Variant, Mutation, MutationSet
import unittest


class TestVariant(unittest.TestCase):
    def setUp(self):
        self.parent_seq = Variant('MAGV')

    def test_add_mutations(self):
        mutation = Mutation('A2[TM]')
        self.parent_seq.add_mutations(mutation)
        self.assertEqual(str(self.parent_seq), 'MTMGV')

    def test_add_mutation_set(self):
        mutation_set = MutationSet.from_string('A2[TM];G3[ST]')
        self.parent_seq.add_mutations(mutation_set)
        self.assertEqual(str(self.parent_seq), 'MTMSTV')

    def test_parse_mutations(self):
        variant1 = Variant('MAVG')
        self.assertEqual(str(self.parent_seq.parse_mutations(variant1)), 'G3V;V4G')

    def test_parse_mutations_indels(self):
        variant2 = Variant('MAGVV')
        self.assertEqual(str(self.parent_seq.parse_mutations(variant2, expect_indels=True)), 'G3[GV]')

    def test_parse_mutations_identical(self):
        variant3 = Variant('MAGV')
        self.assertEqual(str(self.parent_seq.parse_mutations(variant3)), '')

    def test_eq(self):
        variant1 = Variant('MAGV')
        variant2 = Variant('MAGV')
        self.assertEqual(variant1, variant2)
        variant4 = Variant('MAMV', mutation='M3G')
        self.assertEqual(variant1, variant4)

    def test_ne(self):
        variant1 = Variant('MAGV')
        variant2 = Variant('MAVG')
        self.assertNotEqual(variant1, variant2)

        variant3 = Variant('MAGV', mutation='A2[TM]')
        self.assertNotEqual(variant1, variant3)

    def test_parent(self):
        variant = Variant(self.parent_seq, mutation='A2[TM]')
        self.assertTrue(variant.parent is self.parent_seq)

    def test_cut_parent(self):
        variant = Variant(self.parent_seq, mutation='A2[TM]')
        variant.cut_parent()
        self.assertTrue(variant.parent is None)
        self.assertTrue(str(variant) == 'MTMGV')

    def test_cut_children(self):
        variant = Variant(self.parent_seq, mutation='A2[TM]')
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
        variant = Variant(self.parent_seq, mutation='A2[TM]')
        with self.assertRaises(ValueError):
            self.parent_seq.add_mutations(Mutation('V3G'))

        

class TestMutation(unittest.TestCase):
    def setUp(self):
        self.parent_seq = Variant('MAGV')

    def test_apply_subst(self):
        mutation = Mutation('A2M')
        variant_seq = mutation.get_variant_str(self.parent_seq)
        self.assertEqual(variant_seq, 'MMGV')

    def test_apply_del(self):
        mutation = Mutation('A2[-]')
        variant_seq = mutation.get_variant_str(self.parent_seq)
        self.assertEqual(variant_seq, 'AGV')

    def test_apply_multi_del(self):
        mutation = Mutation('[AG]2[--]')
        variant_seq = mutation.get_variant_str(self.parent_seq)
        self.assertEqual(variant_seq, 'MV')

    def test_apply_ins(self):
        mutation = Mutation('A2[AT]')
        variant_seq = mutation.get_variant_str(self.parent_seq)
        self.assertEqual(variant_seq, 'MATGV')

    def test_eq(self):
        mutation1 = Mutation('A2M')
        mutation2 = Mutation('A2M')
        self.assertEqual(mutation1, mutation2)

    def test_ne(self):
        mutation1 = Mutation('A2M')
        mutation2 = Mutation('A2[-]')
        self.assertNotEqual(mutation1, mutation2)


class TestMutationSet(unittest.TestCase):
    def setUp(self):
        self.parent_seq = Variant('MAGV')
        self.mutation_set = MutationSet.from_string('A2[TM];G3[ST]')

    def test_get_variant_str(self):
        variant_str = self.mutation_set.get_variant_str(self.parent_seq)
        self.assertEqual(str(variant_str), 'MTMSTV')

    def test_in(self):
        mutation = Mutation('A2[TM]')
        self.assertIn(mutation, self.mutation_set)

    def test_not_in(self):
        mutation = Mutation('A2M')
        self.assertNotIn(mutation, self.mutation_set)

    def test_eq(self):
        mutation_set1 = MutationSet.from_string('A2[TM];G3[ST]')
        mutation_set2 = MutationSet.from_string('A2[TM];G3[ST]')
        self.assertEqual(mutation_set1, mutation_set2)

    def test_ne(self):
        mutation_set1 = MutationSet.from_string('A2[TM];G3[ST]')
        mutation_set2 = MutationSet.from_string('A2[TM];G3[ST];V4G')
        self.assertNotEqual(mutation_set1, mutation_set2)

        other_parent_seq = Variant('MAVG')
        mutation_set3 = MutationSet.from_string(other_parent_seq, 'A2[TM];G3[ST]')
        self.assertNotEqual(mutation_set1, mutation_set3)

    def test_union(self):
        mutation_set1 = MutationSet.from_string('A2[TM];G3[ST]')
        mutation_set2 = MutationSet.from_string('A2[TM];G3[ST];V4G')
        mutation_set3 = mutation_set1 | mutation_set2
        self.assertEqual(mutation_set3), MutationSet.from_string('A2[TM];G3[ST];V4G')

    def test_intersection(self):
        mutation_set1 = MutationSet.from_string('A2[TM];G3[ST]')
        mutation_set2 = MutationSet.from_string('A2[TM];G3[ST];V4G')
        mutation_set3 = mutation_set1 & mutation_set2
        self.assertEqual(mutation_set3), MutationSet.from_string('A2[TM];G3[ST]')

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
            '[AVG]2-;G5M;D8[DYG];G13V;Y14-')
        variant_str = mutation_set.get_variant_str(parent)
        self.assertEqual(
            variant_str,
            'MMVMDYGEYLMV')


if __name__ == '__main__':
    unittest.main()
