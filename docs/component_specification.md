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
These classes define proteins and functionality to apply mutations to those sequence strings. Data classes are utilized here.

### `Variant`
- `__init__(self, sequence: Union[str, Variant],  mutation: Union[Mutation, MutationSet, str], id: str=None, mutatable: bool=True, labels: VariantLabel=VariantLabel(), round_produced: int=None, round_labeled: int=None)`: Initialize with a sequence.
    - Notes: if not given, id is hash of parent + mutated sequence.
    - Mutatable can be used to prevent new mutations being added. This is done automatically if the variant has children.
- `add_mutations(self, mutation: Union[Mutation, MutationSet])`: Add a mutation to the variant.
- `__str__(self) -> str`: Return the current sequence. Applies each mutation to the parent sequence
- `parse_mutations(self, other: Variant, expect_indels: bool=False, **blast_params) -> MutationSet`: Return a list of mutations that differ from another variant. If all mutations are single point, this is easy, but if there are insertions or deletions, this is more complicated.
- `parent -> Variant`: Return the parent variant if present
- `mutations -> MutationSet`: Return the mutations applied to the parent variant
- `id -> str`: Return the id of the variant if present, otherwise the hash of the variant
- `base_sequence`: Return the base sequence of the variant, without mutations applied
- `hash -> str`: Return the hash of the variant, determined by the sequence after mutations are applied. Or hash of ID if present.
- `mutatable -> bool`: Return whether the variant is mutatable
- `labels -> VariantLabels` Return the labels of the variant.
- `round_added -> int`: Return the round that the variant was created/added.
- `round_putative -> List[int]`: Return the rounds that the variant was putatively generated in.
- `round_experiment -> List[int]`: Return the round that the variant was put up for experimental testing.


### `DBVariant(Variant)`
Information is stored in a database and can be retrieved dynamically when called, but otherwise is not stored, such that we can track parents and children without storing all of the sequences in memory.

- `__init__(self, db: CampaignDatabase, id: str)` overloaded base variant to first check if the variant is in the database, then if so, construct.
- `base_sequence -> str`: Return the base sequence of the variant, without mutations applied, from the database. Does not store in memory.
- `__str__ -> str` Overloaded to retrieve the sequence from the database. 
- `mutations` Overloaded to retrieve the mutations from the database.
- `mutatable` always False
- `parent` overloaded to retrieve the parent from the database.
- `children` overloaded to retrieve the children from the database.
- `label` overloaded to set and retrieve from database and set. Only mutable quantity of a DBVariant.
- `from_variant(cls, variant: Variant, db: CampaignDatabase) -> DBVariant`: Construct a DBVariant from a Variant. Uses CampaignDatabase methods. If the variant is already in the database, return that, otherwise add it to the database and return that.

### `Mutation`
- `__init__(self, position: int, ref: str, alt: str)`
- `from_string(self, mutation_string: str)`: Initialize with a mutation string (e.g., 'A132M').
    - Notes: 
      1. X can be used to indicate any amino acid, useful in combinatorial libraries.
      2. Brackets can be used to indicate indels, eg 'A133[AMVW]' inserts 3 AA after the 132nd position. '[AMWV]132[----]' deletes 4 AA starting at the 132nd position.
      3. We can also represent insertions without the reference amino acid, eg. '>132M' inserts M __before__ the 132nd position. 
      4. Check mutation string is valid, eg. the parent sequence has the correct amino acid at the correct position.
- `get_variant_str(self, variant: Variant) -> str`: Apply the mutation to the variant and show its string.
- `is_insertion -> bool` True if the mutation is an insertion. Note that it only counts single amino acid insertions, eg. 'A132[AMVW]' does not return True
- `is_deletion -> bool` True if the mutation is a deletion. Note that it only counts single amino acid deletions, eg. 'AG132-' does not return True

### `MutationSet`
- `__init__(self, mutations: List[Mutation])`: Initialize with a list of mutations.
- `from_string(cls, parent: Union[Variant, str], mutation_string: str) -> MutationSet`: Initialize from a mutation string. Capable of handling semi colon seperated lists of mutations.
- `__str__(self) -> str`: Return a string of the mutations.

Note for set opertations, Mutations are hashed by position, ref, alt, order.
- `add(self, mutation: Mutation)`: Add a mutation to the set.
- `discard(self, mutation: Mutation)`: Remove a mutation from the set.
- `| -> MutationSet`: Return the union of two mutation sets.
- `& -> MutationSet`: Return the intersection of two mutation sets.
- `- -> MutationSet`: Return the difference of two mutation sets.

- `get_variant_str(self, variant: Variant) -> str`: Return the string of the variant. See apply



---

## 2. Library Classes

Defines a collection of variants and methods to split them, extract certain variants, save to file etc. The package revolvoes around working with these libraries eg. generating them, adding labels to them, filtering them, etc.

### `Library`
- `__init__(self, variants: Iterable[Variant])`: Initialize with a list of variants.
- `get_statistics(self) -> Dict`: Return descriptive statistics.
- `set_labels(self, mapping: Dict[str: float], round: int=None)`: Set supervised for each variant specified by its id.
- `get_unlabeled(self) -> Library`: Return a new library with only unlabeled variants.
- `get_labeled(self) -> Library`: Return a new library with only labeled variants.
- `join(self, other: Library) -> Library`: Return the union of two libraries. If labels are present, they are also joined.
- `save_to_file(self, filename: str)`: Save the library to a file.
- `load_from_file(cls, filename: str, parent_str: str=None, seq_col: str=None, parent_col: str=None, id_col: str=None, mutation_col: str=None, label_col: str=None)`: Load the library from a file.
- `db_save(self, db: CampaignDatabase)`: Save the library to the database.
- `db_load(cls, db: CampaignDatabase, round_idx: Union[int, None]) -> Library`: Load a library from the database.
- `single_parent`: True if all variants have the same parent sequence.
- `variable_residues`: Set of residues that are mutated in the library, only valid if `single_parent` is True.
- `parent`: Parent sequence of the library, only valid if `single_parent` is True.



### `CombinatorialLibrary(Library)`
- `__init__(self, mutations: MutationSet)`: Initialize with a list of mutations.
- `generate(self)`: Generate the library.

---

## 3. Estimator and Transformation Classes

These apply to Libraries. Some functionality include, making predictions, encoding proteins, generating proeins. Use sklearn estimators so that we don't have to rewrite.

### `BaseEstimator`
- `fit(self, X, y)`: Fit the model.
- `predict(self, X) -> Array`: Make predictions.

### `BaseTransformation(BaseEstimator)`
- `transform(self, X) -> Array`: Apply the transformation.

### `SequenceEncoder(BaseTransformation)`: encodes full sequence, arbitrary number of residues
### `MutationEncoder(BaseTransformation)`: encodes variable residues of a fixed length library
  

For AA physicochemical features: PyPEF
---

## 4. Acquisition Function

Define a math transform to select the next variants to test from the output of models. This is the core of the active learning process. This will be stacked into a pipeline after a predictor and a featurizer.

### `AcquisitionFunction`
- `__init__(self, **kwargs)`: Initialize with an estimator and parameters.
- `transform(self)`: Perform the acquisition.

---

## 5. Library Generation

### `LibraryGeneration`
Abstract parent libary generation class.
- `__init__(self, **kwargs)`: Initialize with an estimator and a generation strategy.
- `fit(self, Library)`: Fit the model to a library. __Need to think about inputs here, these are all NN based so nor sure if we can extract out the features from the library inside of the class or if we need to do it before hand with an encoder (prefered).__
- `predict(self) -> Library`: Generate a new library.
- `from_pretrained(cls, **kwargs) -> LibraryGeneration`: Load a pretrained model.

---

## 6. Round
### `BaseRound`
Abstract parent round class. Used to define a step in the lab, eg. this might be generating a new library to test, or  it might be filtering an existing library. These need to be modular and stackable. We need functionality to check if the round is ready to go
- `__init__(self, database: CampaignDatabase, params: Dict)`: Initialize with parameters. Status is None
- `commit_to_db(self)`: Save the round to the database on the top of the stack. Sets the round number.

- `_setup_library_generation(self) -> LibraryGeneration`: Do any setup for library creation step, eg. load or train a sequence generator. Return a LibraryGeneration instance. Should be done from instance attributes only. By default, return a AllUnlabeledGeneration which just returns all unlabeled variants in the database.
- `_library_generator` The setup library generation instance.
- `library_generation(self)`: Set the putative library for the round. Runs the library_generator which returns a Library. See `_setup_library_generation`. The library is saved to the database. If variants in the libary were not already in the database, their round_produced is set to the current round.
- `_setup_library_selection(self) -> AcquisitionFunction`: Do any setup for library selection step, eg. load or train a selection model. Return an AcquisitionFunction instance. Should be done from instance attributes only. By default, return a AllSelection which just returns all variants in the putative library.
- `library_selection(self) -> Library`: Get the library to be used in the experiment. See `_setup_library_selection`. Calls the AquisitionFunction on the putative library. 
- `set_labels(self, notes: str=None, exp_time: str=None)`: Check that labels have been added to the libarary and save, update status to 'complete'.
- `from_db(cls, db: CampaignDatabase, idx: int) -> Round`: Load the round from the database.
- `db`: CampaignDatabase object, None if not connected to a database.
- `round_idx`: Unique identifier for the round, None if not saved/from the database.
- `status`: Status of the round, one of 
    1. 'unknown': Not connected to database at all
    2. 'not ready': Connected to database, but it is a future round and we cannot generate the library yet
    3. 'ready': Connected to database, and we can generate the library
    4. 'generated': Library has been generated
    5. 'selected': library has been selected for exp
    6. 'labeled': library has been labeled and saved to database
    7. 'complete': round complete, made immutable.
- `notes` Notes about the round, eg. what was done, why, etc.
- `start_time` Time the round was started, None if not started.
- `end_time` Time the round was ended, None if not ended.

---

## 7. Runner
Defines multiple rounds, handles saving to database, logging, and monitoring. Should be subclassed to define a specific campaign eg. from a paper.
### `Runner`
- `__init__(self, config: Dict, overwrite: bool=False)`: Initialize with a configuration dictionary. Sets steps
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
    5. status: str, one of 'unknown', 'not ready', 'ready', 'generated', 'selected', 'labeled', 'complete'. See Round class for details.
    6. generated: bool nullable
    7. size: int nullable
    8. labeled_size: int nullable
    9. type: str, name of Round class
    10. params: str, json string of parameters for round

- __Table: Variants__:
    1. idx: int
    2. id: str
    2. round_added: int
    3. round_putative: int nullable
    4. round_experiment: int nullable
    5. parent: str, nullable
    6. mutations: str, nullable
    7. sequence: str
    8. labels: float, nullable

- `__init__(self, path: str, overwrite: bool=False)`: Initialize or load the database.
- `get_variant(self, id: str) -> Variant`: Return the variant with the given id.
- `save_variant(self, variant: Variant)`: Save the variant to the database. If it is already in the database, update it.
- `get_library(self, ids: List[str]=None, round_added_idc: int=None, round_putative_idx: int=None, round_experiment_idx: int=None, only_labeled: bool=True, only_unlabled: bool=False) -> Library`: Return a library.
- `save_library(self, library: Library)`: Save the library to the database. If it is already in the database, update it.
- `current_round_status(self) -> (int, status)`: Return the first the latest round that is at least ready, and its status.
- `get_round(self, idx: int) -> Round`: Return the round with the given index.
- `add_round(self, round: Round) -> int`: Add a round to the database and return its index.
- `save_round(self, round: Round)`: Update round data, eg. status, size, labeled_size, etc. Round must have an index, eg. it is assigned to the database.
- `_check_round(self, round: Round)`: Check that the round is in the database, and that the round has an index, and is the same type, etc.
- `get_rounds(self) -> List[Round]`: Return a list of all rounds in the database.
- `get_rounds_by_status(self, status: str) -> List[Round]`: Return a list of all rounds in the database with the given status.
- `get_current_round(self) -> Round`: Return the latest round that is at least ready.

## Package Structure

```
base/variants.py - Defines the Variant and mutations class, which is the basic unit of a library.
base/library.py - Defines the Library class, which is a collection of variants.
base/estimators.py - Defines abstract classes for estimators, which are used to estimate the fitness of variants or generate variants
base/library_generation.py - Defines abstract classes for library generation, which is used to generate a library from a parent library.
base/library_aquisition.py - Defines abstract classes for library aquisition, which is used to score a library
campaign/database.py - Defines the CampaignDatabase class, which is a database for a campaign.
campaign/round.py - Defines the Round class, which is a step in the campaign.
campaign/runner.py - Defines the Runner class, which is a collection of rounds.
zoo/estimators/* - defines base estimators
zoo/library_generation/* - defines library generation methods 
zoo/library_aquisition/* - defines library aquisition methods
zoo/rounds/* - defines rounds
zoo/runners/* - defines runners from the literature
```
