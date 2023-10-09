from aide import Variant, Mutation, MutationSet
import unittest


class TestVariant(unittest.TestCase):
    def setUp(self):
        self.parent_seq = Variant('MAGV')

    def test_add_mutations(self):
        mutation = Mutation(self.parent_seq, 'A2[TM]')
        self.parent_seq.add_mutations(mutation)
        self.assertEqual(str(self.parent_seq), 'MTMGV')

    def test_add_mutation_set(self):
        mutation_set = MutationSet.from_string(self.parent_seq, 'A2[TM];G3[ST]')
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
        

class TestMutation(unittest.TestCase):
    def setUp(self):
        self.parent_seq = Variant('MAGV')

    def test_apply_subst(self):
        mutation = Mutation(self.parent_seq, 'A2M')
        variant = mutation.apply()
        self.assertEqual(variant, Variant('MMGV'))

    def test_apply_del(self):
        mutation = Mutation(self.parent_seq, 'A2[-]')
        variant = mutation.apply()
        self.assertEqual(variant, Variant('AGV'))

    def test_apply_multi_del(self):
        mutation = Mutation(self.parent_seq, '[AG]2[--]')
        variant = mutation.apply()
        self.assertEqual(variant, Variant('MV'))

    def test_apply_ins(self):
        mutation = Mutation(self.parent_seq, 'A2[AT]')
        variant = mutation.apply()
        self.assertEqual(variant, Variant('MATGV'))

    def test_eq(self):
        mutation1 = Mutation(self.parent_seq, 'A2M')
        mutation2 = Mutation(self.parent_seq, 'A2M')
        self.assertEqual(mutation1, mutation2)

    def test_ne(self):
        mutation1 = Mutation(self.parent_seq, 'A2M')
        mutation2 = Mutation(self.parent_seq, 'A2[-]')
        self.assertNotEqual(mutation1, mutation2)

        other_parent_seq = Variant('MAVG')
        mutation3 = Mutation(other_parent_seq, 'A2M')
        self.assertNotEqual(mutation1, mutation3)


class TestMutationSet(unittest.TestCase):
    def setUp(self):
        self.parent_seq = Variant('MAGV')
        self.mutation_set = MutationSet.from_string(self.parent_seq, 'A2[TM];G3[ST]')

    def test_apply(self):
        variant = self.mutation_set.apply()
        self.assertEqual(str(variant), 'MTMSTV')

    def test_in(self):
        mutation = Mutation(self.parent_seq, 'A2[TM]')
        self.assertIn(mutation, self.mutation_set)

    def test_not_in(self):
        mutation = Mutation(self.parent_seq, 'A2M')
        self.assertNotIn(mutation, self.mutation_set)

    def test_eq(self):
        mutation_set1 = MutationSet.from_string(self.parent_seq, 'A2[TM];G3[ST]')
        mutation_set2 = MutationSet.from_string(self.parent_seq, 'A2[TM];G3[ST]')
        self.assertEqual(mutation_set1, mutation_set2)

    def test_ne(self):
        mutation_set1 = MutationSet.from_string(self.parent_seq, 'A2[TM];G3[ST]')
        mutation_set2 = MutationSet.from_string(self.parent_seq, 'A2[TM];G3[ST];V4G')
        self.assertNotEqual(mutation_set1, mutation_set2)

        other_parent_seq = Variant('MAVG')
        mutation_set3 = MutationSet.from_string(other_parent_seq, 'A2[TM];G3[ST]')
        self.assertNotEqual(mutation_set1, mutation_set3)

    def test_union(self):
        mutation_set1 = MutationSet.from_string(self.parent_seq, 'A2[TM];G3[ST]')
        mutation_set2 = MutationSet.from_string(self.parent_seq, 'A2[TM];G3[ST];V4G')
        mutation_set3 = mutation_set1.union(mutation_set2)
        self.assertEqual(mutation_set3), MutationSet.from_string(self.parent_seq, 'A2[TM];G3[ST];V4G')

    def test_intersection(self):
        mutation_set1 = MutationSet.from_string(self.parent_seq, 'A2[TM];G3[ST]')
        mutation_set2 = MutationSet.from_string(self.parent_seq, 'A2[TM];G3[ST];V4G')
        mutation_set3 = mutation_set1.intersection(mutation_set2)
        self.assertEqual(mutation_set3), MutationSet.from_string(self.parent_seq, 'A2[TM];G3[ST]')

    def test_difference(self):
        mutation_set1 = MutationSet.from_string(self.parent_seq, 'A2[TM];G3[ST]')
        mutation_set2 = MutationSet.from_string(self.parent_seq, 'A2[TM];G3[ST];V4G')
        mutation_set3 = mutation_set2.difference(mutation_set1)
        self.assertEqual(mutation_set3, MutationSet.from_string(self.parent_seq, 'V4G'))

    def test_apply(self):
        mutation_set = MutationSet.from_string(self.parent_seq, 'A2M;G3S')
        variant = mutation_set.apply()
        self.assertEqual(variant, Variant('MMSG'))


if __name__ == '__main__':
    unittest.main()
