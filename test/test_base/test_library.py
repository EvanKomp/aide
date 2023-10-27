import unittest
from aide import Variant, Library, VariantLabel
from io import StringIO
import pandas as pd

class TestLibrary(unittest.TestCase):

    def setUp(self):
        parent = Variant('MAGV')
        self.parent = parent
        variants = [
            Variant(parent, mutations='A2[TM]', id='myid'),
            Variant(parent, mutations='G3[ST]', labels=dict(names=['a'], values=[1], round_idx=[1])),
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

        unlabeled = self.library.get_unlabeled(names=['b'])
        self.assertIsInstance(unlabeled, Library)
        self.assertEqual(len(unlabeled), 2)
        self.assertIn(self.variants[0], unlabeled)
        self.assertIn(self.variants[1], unlabeled)

        unlabeled = self.library.get_unlabeled(round_idx=2)
        self.assertIsInstance(unlabeled, Library)
        self.assertEqual(len(unlabeled), 2)
        self.assertIn(self.variants[0], unlabeled)
        self.assertIn(self.variants[1], unlabeled)

    def test_get_labeled(self):
        labeled = self.library.get_labeled()
        self.assertIsInstance(labeled, Library)
        self.assertEqual(len(labeled), 1)
        self.assertIn(self.variants[1], labeled)

        labeled = self.library.get_labeled(names=['a'])
        self.assertIsInstance(labeled, Library)
        self.assertEqual(len(labeled), 1)
        self.assertIn(self.variants[1], labeled)

        labeled = self.library.get_labeled(round_idx=2)
        self.assertIsInstance(labeled, Library)
        self.assertEqual(len(labeled), 0)


    def test_join(self):
        # new variants should be added, old variants should have their labels
        # joined
        other_variants = [
            Variant(self.parent, mutations='V4T'),
            Variant(self.parent, mutations='G3[ST]', labels=dict(names=['a'], values=[2], round_idx=[2])),
        ]
        other_library = Library(other_variants)

        self.library.join(other_library)

        self.assertEqual(len(self.library), 3)
        self.assertIn(self.variants[0], self.library)
        self.assertIn(self.variants[1], self.library)
        self.assertIn(other_variants[0], self.library)

        # also check label join
        self.assertEqual(len(self.library['-7013149887522224480'].labels), 2)

    def test_build_variants_from_lookup(self):

        lookup = {
            'id1': {
                'base_sequence': 'MAGV',
                'mutations': 'A2[TM]',
                'labels': {
                    'names': ['a'],
                    'values': [1],
                    'round_idx': [1], 
                },
                'parent_id': 'id0',
            },
            'id2': {
                'base_sequence': 'MAGV',
                'mutations': 'G3[ST]',
                'labels': {
                    'names': ['a'],
                    'values': [2],
                    'round_idx': [2],
                },
                'parent_id': 'id0',
            },
            'id0': {
                'base_sequence': 'MAGV',
                'mutations': None,
                'labels': {
                    'names': ['a'],
                    'values': [0],
                    'round_idx': [0],
                },
                'parent_id': None,
            },
            'id3': {
                'base_sequence': 'MAGV',
                'mutations': 'V4T',
                'labels': {
                    'names': ['a'],
                    'values': [3],
                    'round_idx': [3],
                },
                'parent_id': 'id_not_present',
            }

        }

        library = Library.build_variants_from_lookup(variant_ids=['id1', 'id2'], lookup=lookup)
        self.assertEqual(len(library), 2)
        self.assertIn('id1', library)
        self.assertIn('id2', library)
        self.assertEqual(library['id1'].parent.id, 'id0')
        self.assertEqual(library['id2'].parent.id, 'id0')

        # this one should throw a warning
        # because id3 has a parent not in the lookup.
        with self.assertWarns(UserWarning):
            library = Library.build_variants_from_lookup(variant_ids=['id1', 'id2', 'id3'], lookup=lookup)

        # TODO need to add a test with a database

    def test_save_to_file(self):

        outfile = StringIO()
        label_outfile = StringIO()
        self.library.save_to_file(
            variant_file=outfile,
            label_file=label_outfile,
        )

        outfile.seek(0)
        lines = outfile.readlines()
        self.assertEqual(len(lines), 3)
        self.assertEqual(
            lines[0],
            'id,sequence,mutations,parent_id,round_added,round_putative,round_experiment\n'
        )
        self.assertEqual(
            lines[1],
            'myid,MAGV,A2[TM],3640824701896916613,,,\n'
        )
        self.assertEqual(
            lines[2],
            '-7013149887522224480,MAGV,G3[ST],3640824701896916613,,,\n'
        )

        label_outfile.seek(0)
        lines = label_outfile.readlines()
        self.assertEqual(len(lines), 2)
        self.assertEqual(
            lines[0],
            'id,label_name,label_value,round_idx\n'
        )
        self.assertEqual(
            lines[1],
            'myid,a,1.0,1\n'
        )

    def test_load_from_file(self):
        infile = StringIO(
            'id,sequence,mutations,parent_id,round_added,round_putative,round_experiment\n'
            'myid,MAGV,A2[TM],3640824701896916613,,,\n'
            '-7013149887522224480,MAGV,G3[ST],3640824701896916613,,,\n'
        )
        label_infile = StringIO(
            'id,label_name,label_value,round_idx\n'
            'myid,a,1.0,1\n'
        )
        library = Library.load_from_file(
            variant_file=infile,
            label_file=label_infile,
        )

        self.assertEqual(len(library), 2)
        self.assertIn('myid', library)
        self.assertIn('-7013149887522224480', library)
        self.assertEqual(library['myid'].parent.id, '3640824701896916613')
        self.assertEqual(library['-7013149887522224480'].parent.id, '3640824701896916613')
        self.assertEqual(library['myid'].get_label_values(names='a', round_idx=1), [1.0])

    def test_single_parent(self):
        self.assertTrue(self.library.single_parent)

        other_library = Library([
            Variant('MATV'),
            Variant('MAYV'),
        ])
        self.library.join(other_library)
        self.assertFalse(self.library.single_parent)

    def test_parent(self):
        self.assertIs(self.library.parent, self.parent)

    def test_variable_residues(self):
        self.assertEqual(self.library.variable_residues, [1, 2])

    def test_add_labels_df(self):
        df = pd.DataFrame(
            {
                'id': ['myid', '-7013149887522224480'],
                'label_name': ['c', 'c'],
                'label_value': [1, 2],
                'round_idx': [2, 2],
            }
        )
        self.library.add_labels_df(df, label_name_col='label_name', label_value_col='label_value', round_idx_col='round_idx')
        self.assertTrue(len(self.library['myid'].get_label_values('c')) == 1)

    def test_add_labels(self):
        self.library.add_labels(
            variant_ids=['myid'],
            label_names=['c'],
            label_values=[1],
            round_idx=[2])
        self.assertTrue(len(self.library['myid'].get_label_values('c')) == 1)




        


