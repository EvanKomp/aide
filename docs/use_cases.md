# Use cases


## 1. Base interactions with the package

### 1.1 Defining a parent sequence or a standalone variant

```python
from aide import Variant

parent_seq = Variant('MAGV')
```

### 1.2 Defining a variant from a parent

```python
from aide import Variant, Mutation, MutationSet

parent_seq = Variant('MAGV')
mutation = Mutation(parent_seq, 'A2[TM]')
variant = mutation.apply()
# OR variant = Variant(mutation=mutation)
variant
>>> Variant('MAGV', mutation='A2[TM]')
str(variant)
>>> 'MTMGV'

# with mutliple mutations
mutations = MutationSet([Mutation(parent_seq, 'A2[TM]'), ...])
# OR mutations = MutationSet.from_string(parent_seq, 'A2[TM];A3[ST]')
variant = mutations.apply()
```

### 1.2 Labeling a variant and removing labels
    
```python
from aide import Variant

variant = Variant('MAGV')
variant.add_labels(names=['a', 'a', 'b'], values=[1, 2, 3], round=0)

# now remove the labels names 'a', removal probably a rare use case
variant.remove_labels(names=['a'], round=0)
```

### 1.3 Defining a library from a list of variants

```python
from aide import Library, Variant
variants = [Variant('MAGV', id='a'), ...]
library = Library(variants)
library[0].id
>>> 'a'
```

### 1.4 Defining a combinatorial library or mutations
    
```python
from aide import CombinatorialLibrary, MutationSet
mutations = MutationSet.from_string('MAGV', 'A2[TM];A3[ST]') # or [Mutation(...), ...
library = CombinatorialLibrary(mutations)
```

### 1.5 Loading a set of variants from eg. a CSV file from a saturation mutagenesis experiment

CSV contains 'id', and 'mutation' columns. 'id' is a unique identifier for the variant, and 'mutation' is a string of the form 'A123T' where the first character is the wildtype amino acid, the number is the position, and the last character is the mutant amino acid.

```python

from aide import Library

library = Library.from_file('mutations.csv', id_col='id', mutation_col='mutation', parent_seq='MAGV') # for mutagenesis experiment
library = Library.from_file('mutations.csv', id_col='id', seq_col='sequence') # for a library of sequences
```

### 1.6 Parsing mutations from two variants
    
```python
from aide import Variant

# easy when there is no indels
parent_seq = Variant('MAGV')
variant1 = Variant('MAVG')

parent_seq.parse_mutations(variant1)
>>> MutationSet('G3V;V4G')

# harder when there are indels
variant2 = Variant('MAVGV')
parent_seq.parse_mutations(variant2, expect_indels=True, **blast_params)
>>> MutationSet('A2[AV]')
```

### 1.7 Generating a library from a LibraryGenerator

```python
from aide import MSALibraryGenerator

# pretrained
generator = MSALibraryGenerator('path_to_model_save')

# OR train
generator = MSALibraryGenerator()
generator.train('path_to_msa_file', **hyperparams)

# generate
library = generator.generate(parent_seq, number_of_variants=100)
```

### 1.8 Using an aquisition function to filter a library

```python
from aide import Library
from aide.acquisition import GPAcquisitionFunction
from aide.features import ResiduePhysicalFeatures

training_library = Library.from_file('training_library.csv', id_col='id', seq_col='sequence', label_col='label')
library = Library.from_file('library.csv', id_col='id', seq_col='sequence')
featurizer = ResiduePhysicalFeatures() # figures out which positions are variable when fit

acquisition_function = GPAcquisitionFunction(featurizer=featurizer, strategy='UCB')
acquisition_function.fit(training_library)

filtered_library = acquisition_function.filter(library, number_of_variants=10)
```

## 2. Running DE campaigns

### 2.1 Running a single round of a selection DE campaign

```python
from aide import Library
from aide.rounds import GPSelectionRound

round = GPSelectionRound(features='ResiduePhysicalFeatures', strategy='UCB', number_of_variants=10)

library = Library.from_file('library.csv', id_col='id', seq_col='sequence', label_col='label') # some of these are labeled, some not
round.set_starting_library(library)

to_experiment = round.get_library_for_exp() # this is the library that should be tested in the lab
to_experiment.save_to_file('to_experiment.csv') # run it

# after lab
experimental_results = Library.from_file('experimental_results.csv', id_col='id', seq_col='sequence', label_col='label')
round.set_labels(experimental_results) 
```

### 2.2 Running a single round of a generation DE campaign

```python
from aide import Library
from aide.rounds import MSAGenerationRound

round = MSAGenerationRound('path_to_model_save', number_of_variants=100)
to_experiment = round.get_library_for_exp() # this is the library that should be tested in the lab
to_experiment.save_to_file('to_experiment.csv') # run it

# after lab
experimental_results = Library.from_file('experimental_results.csv', id_col='id', seq_col='sequence', label_col='label')
round.set_labels(experimental_results) 
```

### 2.3 Running a full predefined DE campaign

```python
from aide.zoo import Saito2018Runner, Saito2018Config

config = Saito2018Config(
    name='my_campaign',
    description='my campaign description',
    database='./my_campaign.db'
    number_of_rounds=2,
    variants_per_round=10,
    ) # initial library determined by saturation mutagenesis

runner = Saito2018Runner(config)
runner.step()
>>> CampaignError: Cannot take step of round 0, status is 'not started'
runner.step(data_path='initial_library.csv', parent_seq='MVKMG', id_col='id', mutations_col='mut') # the library from mutagenesis, etc
>>> Round 0 of type RandomGenerationRound() status is 'ready'
runner.step()
>>> Round 0 of type RandomGenerationRound() status is 'started', library has 10 variants at file path 'library_for_exp.csv'
runner.step(data_path='experimental_results.csv', id_col='id', mutations_col='mut', label_col='label')
>>> Round 0 of type RandomGenerationRound() status is 'complete'
runner.step()
>>> Round 1 of type GPSelectionRound(features='ResiduePhysicalFeatures', strategy='InfoMax', number_of_variants=10) status is 'started', library has 10 variants at file path 'library_for_exp.csv'
runner.step(data_path='experimental_results.csv', id_col='id', mutations_col='mut', label_col='label')
>>> Round 1 of type GPSelectionRound(features='ResiduePhysicalFeatures', strategy='InfoMax', number_of_variants=10) status is 'complete'
runner.step()
>>> Round 2 of type GPSelectionRound(features='ResiduePhysicalFeatures', strategy='greedy', number_of_variants=10) status is 'started', library has 10 variants at file path 'library_for_exp.csv'

# the final 10 are selected greedily
```

### 2.4 Running a full custom DE campaign

Note that the above example would is just a wrapper of code below to manually define a campaign

```python
from aide import Runner, RunnerConfig

config = RunnerConfig(
    name='my_campaign',
    description='my campaign description',
    database='./my_campaign.db'
    rounds = [
        {
            'type': 'RandomGenerationRound',
            'params': {
                'number_of_variants': 10,
            }
        },
        {
            'type': 'GPSelectionRound',
            'params': {
                'number_of_variants': 10,
                'features': 'ResiduePhysicalFeatures',
                'strategy': 'InfoMax'
            }
        },
        {
            'type': 'GPSelectionRound',
            'params': {
                'number_of_variants': 10,
                'features': 'ResiduePhysicalFeatures',
                'strategy': 'greedy'
            }
        }     

    ]
)

runner = Runner(config)
runner.step()
>>> CampaignError: Cannot take step of round 0, status is 'not started'
runner.step(data_path='initial_library.csv', parent_seq='MVKMG', id_col='id', mutations_col='mut') # the library from mutagenesis, etc
>>> Round 0 of type RandomGenerationRound() status is 'ready'
runner.step()
>>> Round 0 of type RandomGenerationRound() status is 'started', library has 10 variants at file path 'library_for_exp.csv'
runner.step(data_path='experimental_results.csv', id_col='id', mutations_col='mut', label_col='label')
>>> Round 0 of type RandomGenerationRound() status is 'complete'
runner.step()
>>> Round 1 of type GPSelectionRound(features='ResiduePhysicalFeatures', strategy='InfoMax', number_of_variants=10) status is 'started', library has 10 variants at file path 'library_for_exp.csv'
runner.step(data_path='experimental_results.csv', id_col='id', mutations_col='mut', label_col='label')
>>> Round 1 of type GPSelectionRound(features='ResiduePhysicalFeatures', strategy='InfoMax', number_of_variants=10) status is 'complete'
runner.step()
>>> Round 2 of type GPSelectionRound(features='ResiduePhysicalFeatures', strategy='greedy', number_of_variants=10) status is 'started', library has 10 variants at file path 'library_for_exp.csv'

# the final 10 are selected greedily
```






