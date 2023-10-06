from aide import Variant, Mutation, MutationSet
import unittest


class TestVariant(unittest.TestCase):
    def setUp(self):
        self.parent_seq = Variant('MAGV')

    def test_add_mutations(self):
        mutation = Mutation(self.parent_seq, 'A2[TM]')
        self.parent_seq.add_mutations(mutation)
        self.assertEqual(str(self.parent_seq), 'MTMGV')

    def test_parse_mutations(self):
        variant1 = Variant('MAVG')
        self.assertEqual(str(self.parent_seq.parse_mutations(variant1)), 'G3V;V4G')


class TestMutation(unittest.TestCase):
    def setUp(self):
        self.parent_seq = Variant('MAGV')

    def test_apply_subst(self):
        mutation = Mutation(self.parent_seq, 'A2M')
        variant = mutation.apply()
        self.assertEqual(str(variant), 'MMGV')

    def test_apply_del(self):
        mutation = Mutation(self.parent_seq, 'A2[-]')
        variant = mutation.apply()
        self.assertEqual(str(variant), 'AGV')

    def test_apply_multi_del(self):
        mutation = Mutation(self.parent_seq, '[AG]2[--]')
        variant = mutation.apply()
        self.assertEqual(str(variant), 'MV')

    def test_apply_ins(self):
        mutation = Mutation(self.parent_seq, 'A2[AT]')
        variant = mutation.apply()
        self.assertEqual(str(variant), 'MATGV')


class TestMutationSet(unittest.TestCase):
    def setUp(self):
        self.parent_seq = Variant('MAGV')
        self.mutation_set = MutationSet.from_string(self.parent_seq, 'A2[TM];G3[ST]')

    def test_apply(self):
        variant = self.mutation_set.apply()
        self.assertEqual(str(variant), 'MTMSTV')


if __name__ == '__main__':
    unittest.main()
