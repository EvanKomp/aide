from __future__ import annotations
from dataclasses import dataclass, field
from typing import Union, List, Iterable, Dict
from collections.abc import MutableSet
import pandas as pd
import json

@dataclass(frozen=True, eq=True)
class VariantLabel:
    variant_id: str
    name: str
    value: Union[float, int, str]
    round_idx: int = field(default=None, hash=True)

class VariantLabels(MutableSet):

    def __init__(self, labels: Iterable[Union[VariantLabel, Dict]] = None):
        self.labels = set()
        if labels:
            self.add_labels(labels)

    def add_labels(self, labels: Union[VariantLabels, List[Union[Dict, VariantLabel]]]):
        for label in labels:
            
            if isinstance(label, Dict):
                to_add = VariantLabel(**label)
            else:
                to_add = label

            self.add(to_add)

    def remove_labels(self, labels: Union[VariantLabels, List[Union[Dict, VariantLabel]]]):
        for label in labels:
            if isinstance(label, Dict):
                self.discard(VariantLabel(**label))
            else:
                self.discard(label)

    def has_labels(self, names: Union[str, Iterable[str]]=None, round_idx: Union[int, Iterable[int]] = None) -> bool:
        if names is None:
            return len([label for label in self.labels if label.round_idx == round_idx if round_idx is not None]) > 0
        elif type(names) == str:
            names = [names]
        return all([
            any(label.name == name and (label.round_idx == round_idx if round_idx is not None else True) for label in self.labels) for name in names
            ])

    @classmethod
    def from_jsons(cls, jsons: Iterable[str]) -> VariantLabels:
        dicts = [json.loads(j) for j in jsons]
        return cls(dicts)

    @property
    def df(self) -> pd.DataFrame:
        return pd.DataFrame([label.__dict__ for label in self.labels])
    
    @property
    def variant_id(self) -> List[str]:
        # all ids should be the same
        if self.labels:
            return list(set(label.variant_id for label in self.labels))[0]
        else:
            return None

    @classmethod
    def from_db(cls, db: 'CampaignDatabase', variant_id: str) -> VariantLabels:
        # Simulate database fetch, actual implementation may vary.
        labels_json = db.get_labels_json(variant_id=variant_id)
        return cls(labels_json)

    def to_db(self, db: 'CampaignDatabase'):
        db.add_labels(self.labels)

    def remove_from_db(self, db: 'CampaignDatabase'):
        db.remove_labels(self.labels)

    @property
    def schema(self) -> List[str]:
        return set(label.name for label in self.labels)

    def select(self, name: Union[str, List[str]], round_idx: Union[int, Iterable[int]] = None) -> 'VariantLabels':
        filtered_labels = {label for label in self.labels if label.name in name and (label.round_idx == round_idx if round_idx is not None else True)}
        return VariantLabels(filtered_labels)

    def get_values(self, name: str, round_idx: Union[int, Iterable[int]] = None) -> List[float]:
        return [label.value for label in self.labels if label.name == name and (label.round_idx == round_idx if round_idx is not None else True)]

    # Required by MutableSet
    def __contains__(self, x: object) -> bool:
        return x in self.labels

    def __iter__(self):
        return iter(self.labels)

    def __len__(self) -> int:
        return len(self.labels)

    def add(self, value: VariantLabel) -> None:
        if self.variant_id and value.variant_id != self.variant_id:
            raise ValueError("All labels must have the same variant_id")
        self.labels.add(value)

    def discard(self, value: VariantLabel) -> None:
        self.labels.discard(value)
