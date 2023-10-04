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
- `__init__(self, parent_sequence: str,  mutation: Union[Mutation, List[Mutation]], id: str=None)`: Initialize with a sequence.
    - Notes: if not given, id is hash of parent + mutated sequence.
- `add_mutations(self, mutation: Union[Mutation, List[Mutation]])`: Add a mutation to the variant.
- `__repr__(self) -> str`: Return the current sequence. Applies each mutation to the parent sequence
- `parse_mutations(self, other: Variant, expect_indels: bool=False) -> List[Mutation]`: Return a list of mutations that differ from another variant. If all mutations are single point, this is easy, but if there are insertions or deletions, this is more complicated.

### `Mutation`
- `__init__(self, parent: Variant, mutation_string: str)`: Initialize with a mutation string (e.g., 'A132M').
    - Notes: 
      1. X can be used to indicate any amino acid, useful in combinatorial libraries.
      2. Brackets can be used to indicate indels, eg 'A132[AMVW]' inserts 3 AA after the 132nd position. '[AMWV]132[----]' deletes 4 AA starting at the 132nd position.
- `apply(self) -> Variant`: Apply the mutation to the variant, returns a new Variant

---

## 2. Library Classes

### `Library`
- `__init__(self, variants: List[Variant], labels: List[float]=None, round: Round=None)`: Initialize with a list of variants.
- `get_statistics(self) -> Dict`: Return descriptive statistics.
- `set_labels(self, labels: List[float])`: Set supervised for each variant.
- `save_to_file(self, filename: str)`: Save the library to a file.
- `load_from_file(cls, filename: str)`: Load the library from a file.
- `db_save(self, db: CampaignDatabase)`: Save the library to the database.

### `CombinatorialLibrary(Library)`
- `__init__(self, mutations: List[Mutation])`: Initialize with a list of mutations.
- `generate(self)`: Generate the library.

---

## 3. Estimator and Transformation Classes

Use sklearn.

### `BaseEstimator`
- `fit(self, X, y)`: Fit the model.
- `predict(self, X) -> Array`: Make predictions.

### `BaseTransformation(BaseEstimator)`
- `transform(self, X) -> Array`: Apply the transformation.
  
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
Abstract parent round class. Used to define a step in the labe, eg. this might be generating a new library to test,
or  it might be filtering an existing library
- `__init__(self, params: Dict)`: Initialize with parameters. No library yet
- `set_library(self, library: Library)`: Set the library for this round.
- `get_library(self, db: CampaignDatabase) -> Library`: Return the library for this round from the database if it exists, or generate it and save it to the database.
- `execute(self) -> Any`: Execute the lab step.
- `save(self, db: CampaignDatabase)`: Save the round to the database.
- `_id`: Unique identifier for the round, None if not saved/from the database.

### `RoundLibraryGeneration(Round)`, `RoundFiltering(Round)`, `RoundScreening(Round)`

---

## 7. Runner
Defines multiple rounds, handles saving to database, logging, and monitoring.
### `Runner`
- `__init__(self, config: Dict)`: Initialize with a configuration dictionary.
- `run(self)`: Execute the full pipeline.

---

## 8. Database Interaction
Libraries are loaded and saved to the database. Scheme for a campaign:

### `CampaignDatabase`
- __Table: Rounds__:
    1. idx: int
    2. notes: str
    3. start_time: datetime
    4. end_time: datetime nullable
    5. status: str, one of ['started', 'complete']. Eg. either we have generated/input a library but have no labels, or we have labels and are ready for the next round
    6. generated: bool nullable
    7. size: int nullable
    8. labeled_size: int nullable

- __Table: Variants__:
    1. idx: int
    2. id: str
    3. round: int
    4. parent: str, nullable
    5. mutations: str, nullable
    6. sequence: str
    7. labels: float, nullable

- `__init__(self, path: str)`: Initialize or load the database.



