import unittest
from aide import Variant, Library, VariantLabel

class TestLibrary:

    def setUp(self):
        parent = Variant('MAGV')
        variants = [
            Variant(parent, mutations='A2[TM]', id='myid'),
            Variant(parent, mutations='G3[ST]', label=VariantLabel({'a': [1]}, round_idx=1)),
        ]
        self.variants = variants
        self.library = Library(variants)

    def test_bad_input(self):
        not_variants = ['a', 1, 1.0]
        with self.assertRaises(TypeError):
            Library(not_variants)

    def test_in(self):
        self.assertIn(self.variants[0], self.library)
        other_variant = Variant(self.variants[0], mutations='A2M')
        self.assertNotIn(other_variant, self.library)

    def test_len(self):
        self.assertEqual(len(self.library), 2)

    def test_get(self):
        self.assertIs(self.library['myid'], self.variants[0])
        self.assertIs(self.library[self.variants[1].id], self.variants[1])

    def test_get_bad_id(self):
        with self.assertRaises(KeyError):
            self.library['bad_id']

    def test_get_unlabeled(self):
        unlabeled = self.library.get_unlabeled()
        self.assertIsInstance(unlabeled, Library)
        self.assertEqual(len(unlabeled), 1)
        self.assertIn(self.variants[0], unlabeled)

    def test_get_labeled(self):
        labeled = self.library.get_labeled()
        self.assertIsInstance(labeled, Library)
        self.assertEqual(len(labeled), 1)
        self.assertIn(self.variants[1], labeled)

    def test_join(self):
        # new variants should be added, old variants should have their labels
        # joined
        other_variants = [
            Variant(self.parent, mutations='V4T'),
            Variant(self.parent, mutations='G3[ST]', label=VariantLabel({'a': [2]}, round_idx=2)),
        ]
        other_library = Library(other_variants)

        self.library.join(other_library)

        self.assertEqual(len(self.library), 3)
        self.assertIn(self.variants[0], self.library)
        self.assertIn(self.variants[1], self.library)
        self.assertIn(other_variants[0], self.library)


