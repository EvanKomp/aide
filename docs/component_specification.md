# Component Specification for AIDE

## Table of Contents
1. Variant and Mutation Classes
2. Library Classes
3. Estimator and Transformation Classes
4. Acquisition Function
5. Library Generation
6. Lab Step
7. Runner
8. Database Interaction
9. Logging and Monitoring

---

## 1. Variant and Mutation Classes

### `Variant`
- `__init__(self, parent_sequence: str=None,  mutation: Union[Mutation, MutationSet], id: str=None)`: Initialize with a sequence.
    - Notes: if not given, id is hash of parent + mutated sequence.
- `add_mutations(self, mutation: Union[Mutation, MutationSet])`: Add a mutation to the variant.
- `__str__(self) -> str`: Return the current sequence. Applies each mutation to the parent sequence
- `parse_mutations(self, other: Variant, expect_indels: bool=False, **blast_params) -> MutationSet`: Return a list of mutations that differ from another variant. If all mutations are single point, this is easy, but if there are insertions or deletions, this is more complicated.

### `Mutation`
- `__init__(self, parent: Variant, mutation_string: str)`: Initialize with a mutation string (e.g., 'A132M').
    - Notes: 
      1. X can be used to indicate any amino acid, useful in combinatorial libraries.
      2. Brackets can be used to indicate indels, eg 'A132[AMVW]' inserts 3 AA after the 132nd position. '[AMWV]132[----]' deletes 4 AA starting at the 132nd position.
      3. Check mutation string is valid, eg. the parent sequence has the correct amino acid at the correct position.
- `apply(self) -> Variant`: Apply the mutation to the variant, returns a new Variant

### `MutationSet`
- `__init__(self, mutations: List[Mutation])`: Initialize with a list of mutations.
- `from_string(cls, parent: Union[Variant, str], mutation_string: str) -> MutationSet`: Initialize from a mutation string. Capable of handling semi colon seperated lists of mutations.
- `union(self, other: MutationSet) -> MutationSet`: Return the union of two mutation sets.
- `intersection(self, other: MutationSet) -> MutationSet`: Return the intersection of two mutation sets.
- `difference(self, other: MutationSet) -> MutationSet`: Return the difference of two mutation sets.
- `apply(self) -> Variant`: Apply the mutations to the variant, returns a new Variant. Careful with positioning, break into list first.

---

## 2. Library Classes

### `Library`
- `__init__(self, variants: List[Variant], labels: List[float]=None, round: Round=None)`: Initialize with a list of variants.
- `get_statistics(self) -> Dict`: Return descriptive statistics.
- `set_labels(self, mapping: Dict[str: float])`: Set supervised for each variant specified by its id.
- `get_unlabeled(self) -> Library`: Return a new library with only unlabeled variants.
- `get_labeled(self) -> Library`: Return a new library with only labeled variants.
- `join(self, other: Library) -> Library`: Return the union of two libraries. If labels are present, they are also joined.
- `save_to_file(self, filename: str)`: Save the library to a file.
- `load_from_file(cls, filename: str, parent: str=None, id_col: str=None, mutation_col: str=None, label_col: str=None)`: Load the library from a file.
- `db_save(self, db: CampaignDatabase)`: Save the library to the database.
- `db_load(cls, db: CampaignDatabase, idx: Union[int, None], only_labeled: bool=True) -> Library`: Load a library from the database.
- `_single_parent`: True if all variants have the same parent sequence.
- `_variable_residues`: Set of residues that are mutated in the library, only valid if `_single_parent` is True.

### `CombinatorialLibrary(Library)`
- `__init__(self, mutations: MutationSet)`: Initialize with a list of mutations.
- `generate(self)`: Generate the library.

---

## 3. Estimator and Transformation Classes

Use sklearn. Used for both encoding and models.

### `BaseEstimator`
- `fit(self, X, y)`: Fit the model.
- `predict(self, X) -> Array`: Make predictions.

### `BaseTransformation(BaseEstimator)`
- `transform(self, X) -> Array`: Apply the transformation.

### `SequenceEncoder(BaseTransformation)`: encodes full sequence, arbitrary number of residues
### `MutationEncoder(BaseTransformation)`: encodes variable residues of a fixed length library
  
---

## 4. Acquisition Function

### `AcquisitionFunction`
- `__init__(self, estimator: BaseEstimator, params: Dict)`: Initialize with an estimator and parameters.
- `acquire(self) -> Variant`: Perform the acquisition.

---

## 5. Library Generation

### `LibraryGeneration`
Abstract parent libary generation class.
- `__init__(self, estimator: BaseEstimator, strategy: str)`: Initialize with an estimator and a generation strategy.
- `generate(self) -> Library`: Generate a new library.

---

## 6. Round
### `Round`
Abstract parent round class. Used to define a step in the lab, eg. this might be generating a new library to test,
or  it might be filtering an existing library
- `__init__(self, params: Dict)`: Initialize with parameters. No library yet. Status is None
- `set_starting_library(self, library: Library)`: Set the library for this round. Save to database if present.
- `get_library_for_exp(self) -> Library`: Return the library for this round from the database if it exists, or generate it and save it to the database.
- `set_labels(self, notes: str=None, exp_time: str=None)`: Check that labels have been added to the libarary and save, update status to 'complete'.
- `to_db(self, db: CampaignDatabase)`: Save the round to the database.
- `from_db(cls, db: CampaignDatabase, idx: int) -> Round`: Load the round from the database.
- `_id`: Unique identifier for the round, None if not saved/from the database.
- `_status`: Status of the round, one of 
    1. 'unknown': Not connected to database at all
    2. 'not started': Connected to database, but it is a future round and we cannot generate the library yet
    3. 'ready': Connected to database, and we can generate the library
    4. 'started': Library has been generated, but not finished
    5. 'complete': We have labels saved back to the database.
_ `_metadata`: Dictionary of metadata for the round, eg. notes, start time, end time, etc.

### `RoundLibraryGeneration(Round)`, `RoundFiltering(Round)`, `RoundScreening(Round)`
These are subclassed further for specific library generation methods, filtering methods, and screening methods from the literature.

---

## 7. Runner
Defines multiple rounds, handles saving to database, logging, and monitoring. Should be subclassed to
### `Runner`
- `__init__(self, config: Dict)`: Initialize with a configuration dictionary. Sets steps
- `get_status(self) -> str`: Return the status of the campaign, eg. which round we are on and whether it is complete.
- `_steps`: List of Rounds.
- `step(self, data_path: str=None, id_col: str=None, seq_col: str=None, mutations_col: str=None, labels_col: str=None, exp_time: str=None, notes: str=None)`: Execute the next step in the pipeline. Will give library if that is the next step, or will take library and labels if that is the next step, depending on the round status.

### `RunnerConfig`
Dataclass that defines the configuration for a Runner. Load from json file.
Expected fields:
- `name: str`: Name of the campaign
- `description: str`: Description of the campaign
- `database: str`: Path to the database file.
- `default_columns`: Dictionary of default column names for loading libraries from file. Keys will be passed as kwargs to the library load function.
- `rounds: List[Dict]`: List of dictionaries defining the rounds. Each dictionary should have a 'type' field, which is the name of the Round subclass, and a 'params' field, which is a dictionary of parameters for the Round subclass.

---

## 8. Database Interaction
Libraries are loaded and saved to the database. Scheme for a campaign:

### `CampaignDatabase`
- __Table: Rounds__:
    1. idx: int
    2. notes: str
    3. start_time: datetime
    4. end_time: datetime nullable
    5. status: str, one of 'unknown', 'not started', 'ready', 'started', 'complete'. See Round class for details.
    6. generated: bool nullable
    7. size: int nullable
    8. labeled_size: int nullable
    9. type: str, name of Round class
    10. params: str, json string of parameters for round

- __Table: Variants__:
    1. idx: int
    2. id: str
    3. round: int
    4. parent: str, nullable
    5. mutations: str, nullable
    6. sequence: str
    7. labels: float, nullable

- `__init__(self, path: str)`: Initialize or load the database.
- `current_round_status(self) -> (int, status)`: Return the first the latest round that is at least ready, and its status.
- `get_round(self, idx: int) -> Round`: Return the round with the given index.
- `get_current_round(self) -> Round`: Return the latest round that is at least ready.



