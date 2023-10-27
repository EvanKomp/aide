import unittest
from aide.base import VariantLabel, VariantLabels
from unittest.mock import MagicMock
import pandas as pd

class TestVariantLabel(unittest.TestCase):

    def test_init(self):
        label = VariantLabel(variant_id="v1", name="label1", value=1, round_idx=1)
        self.assertEqual(label.variant_id, "v1")
        self.assertEqual(label.name, "label1")
        self.assertEqual(label.value, 1)
        self.assertEqual(label.round_idx, 1)

class TestVariantLabels(unittest.TestCase):

    def setUp(self):
        self.mock_db = MagicMock()
        self.mock_db.get_labels_json.return_value = [
            {"variant_id": "v1", "name": "label1", "value": 1, "round_idx": 1},
            {"variant_id": "v1", "name": "label2", "value": 2, "round_idx": 1},
        ]
        self.mock_db.add_labels.return_value = None
        self.mock_db.remove_labels.return_value = None
        self.labels_list = [
            VariantLabel(variant_id="v1", name="label1", value=1, round_idx=1),
            VariantLabel(variant_id="v1", name="label2", value=2, round_idx=1),
        ]
        self.labels = VariantLabels(self.labels_list)

    def test_init(self):
        labels = [
            VariantLabel(variant_id="v1", name="label1", value=1, round_idx=1),
            VariantLabel(variant_id="v1", name="label2", value=2, round_idx=1),
        ]
        labels = VariantLabels(labels)
        self.assertIsInstance(labels, VariantLabels)

        labels = [
            VariantLabel(variant_id="v1", name="label1", value=1, round_idx=1),
            VariantLabel(variant_id="v1", name="label2", value=2, round_idx=1),
            VariantLabel(variant_id="v2", name="label1", value=1, round_idx=1),
        ]
        with self.assertRaises(ValueError):
            labels = VariantLabels(labels)
        
        labels = [
            {"variant_id": "v1", "name": "label1", "value": 1, "round_idx": 1},
            {"variant_id": "v1", "name": "label2", "value": 2, "round_idx": 1},
        ]
        labels = VariantLabels(labels)

        labels = [
            {"bad_key": "v1", "name": "label1", "value": 1, "round_idx": 1},
            {"variant_id": "v1", "name": "label2", "value": 2, "round_idx": 1},
        ]
        with self.assertRaises(TypeError):
            labels = VariantLabels(labels)

    def test_schema(self):
        self.assertEqual(self.labels.schema, set(["label2", "label1"]))

    def test_select(self):
        selected = self.labels.select(name="label1", round_idx=1)
        self.assertIsInstance(selected, VariantLabels)
        self.assertEqual(selected.schema, {"label1"})

    def test_get_values(self):
        values = self.labels.get_values("label1", round_idx=1)
        self.assertIsInstance(values, list)
        self.assertEqual(values[0], 1)

    def test_add_labels(self):
        labels = VariantLabels()
        labels.add_labels([VariantLabel(variant_id="v1", name="label1", value=1, round_idx=1)])
        labels.add_labels([{ "variant_id": "v1", "name": "label2", "value": 2, "round_idx": 1 }])
        self.assertTrue(labels.has_labels(["label1"]))
        self.assertTrue(labels.has_labels(["label2"]))

    def test_remove_labels(self):
        self.labels.remove_labels([dict(variant_id="v1", name="label2", value=2, round_idx=1)])
        self.assertFalse(self.labels.has_labels(["label2"]))

    def test_has_labels(self):
        self.assertTrue(self.labels.has_labels(["label1", "label2"], round_idx=1))
        self.assertFalse(self.labels.has_labels(["label1", "label2"], round_idx=2))
        self.assertFalse(self.labels.has_labels(["label1", "label3"], round_idx=1))

    def test_df(self):
        labels = VariantLabels()
        self.assertIsInstance(labels.df, pd.core.frame.DataFrame)

    def test_from_db(self):
        labels = VariantLabels.from_db(self.mock_db, variant_id="v1")
        self.assertIsInstance(labels, VariantLabels)
        self.mock_db.get_labels_json.assert_called_once_with(variant_id="v1")

    def test_to_db(self):
        self.labels.to_db(self.mock_db)
        self.mock_db.add_labels.assert_called_once_with(self.labels)

    def test_remove_from_db(self):
        self.labels.remove_from_db(self.mock_db)
        self.mock_db.remove_labels.assert_called_once_with(self.labels)