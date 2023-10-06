from aide import Variant, Mutation, MutationSet, Library, CombinatorialLibrary

class TestLibrary(unittest.TestCase):
    def setUp(self):
        self.variant1 = Variant('MAGV', id='a')
        self.variant2 = Variant('MAGV', id='b')
        self.variants = [self.variant1, self.variant2]
        self.library = Library(self.variants)

    def test_get_statistics(self):
        stats = self.library.get_statistics()
        self.assertIsInstance(stats, dict)

    def test_set_labels(self):
        self.library.set_labels({'a': 1.0, 'b': 2.0})
        self.assertEqual(self.library.get_labeled(), [self.variant1, self.variant2])

    def test_get_unlabeled(self):
        self.assertIsInstance(self.library.get_unlabeled(), Library)

    def test_get_labeled(self):
        self.assertIsInstance(self.library.get_labeled(), Library)

    def test_join(self):
        variant3 = Variant('MAGV', id='c')
        other_library = Library([variant3])
        joined_library = self.library.join(other_library)
        self.assertIsInstance(joined_library, Library)
        self.assertEqual(len(joined_library.variants), 3)

    def test_combinatorial_library(self):
        mutation_set = MutationSet.from_string('MAGV', 'A2[TM];A3[ST]')
        comb_lib = CombinatorialLibrary(mutation_set)
        self.assertIsInstance(comb_lib, CombinatorialLibrary)
        self.assertEqual(len(comb_lib.variants), 4)

if __name__ == '__main__':
    unittest.main()
