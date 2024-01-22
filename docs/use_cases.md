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
mutation = Mutation(position=2, ref='A', alt='TM')
variant = Variant(parent_seq, mutation=mutation)

variant
>>> Variant('MAGV', mutation='A2[TM]')
str(variant)
>>> 'MTMGV'

# with mutliple mutations
mutations = MutationSet([Mutation.from_string('A2[TM]'), ...])
# OR mutations = MutationSet.from_string(parent_seq, 'A2[TM];A3[ST]')
variant = Variant(parent_seq, mutations=mutations)
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

### 1.4 Loading a set of variants from eg. a CSV file from a saturation mutagenesis experiment or a list of sequences

Labels are tracked seperately.

```python

from aide import VariantLibrary, MutationLibrary

# for a list of sequences
library = VaraintLibrary.load_from_file('sequences.csv', schema={
    'id_col': 'id',
    'seq_col': 'sequence'
    })
# also include labels
# here there is a second file with mutiple tyes of labels for each sequence
library = VariantLibrary.load_from_file(
    'sequences.csv',
    schema={
        'id_col': 'id',
        'seq_col': 'sequence',
    },
    label_file='labels.csv',
    label_schema={
        'id_col': 'id',
        'name_col': 'label_type'
        'value_col': 'label_value'
    }
)

# from a saturation mutagenesis experiment that has been phenotyped
# lets say there is just one columns, 'mutation'
library = MutationLibrary.load_from_file(
    'saturation_mutagenesis.csv',
    schema={
        'id_col': 'mutation',
        'mutation_col': 'mutation',
    },
)
# convert to a library of variants with a parent sequence
variant_library = library.apply_parent(parent='MAGV')

# can also create a combinatorial library, though I am not sure this is
# relevant experimentally, as on the becnh we do not have control over
# the genotype? So this may be an unused feature also of type MutationLibrary
# this is a generator of mutation combinations
comb_library = library.combinations()
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
library
>>> VariantLibrary(size=100, parent=Variant('MAGV'), fixed_size=False)
```

### 1.8 Using an aquisition function to filter a library

```python
from aide import VariantLibrary
from aide.acquisition import GPAcquisitionFunction
from aide.features import ResiduePhysicalFeatures

# this is a library of variants that have been tested in the lab
# assume it has multiple labels, eg. 'a', 'b', 'c'
# but we only care about 'a' and 'b'
training_library = VariantLibrary.load_from_file(
    'training_library.csv',
    schema={
        'id_col': 'id',
        'seq_col': 'sequence',
    },
    label_file='labels.csv',
    label_schema={
        'id_col': 'id',
        'name_col': 'label_type'
        'value_col': 'label_value'
    })

acquisition_function = GPAcquisitionFunction(
    features=ResiduePhysicalFeatures(), # this does not work if the library is not a fixed sequence length
    labels=['a', 'b'],
    number_of_variants=10,
    strategy='InfoMax',
    multiobjective_aggregation='sum',
)
acquisition_function.fit(training_library)

# this is a library of variants that have not been tested in the lab
library = VariantLibrary.load_from_file(
    'library.csv',
    schema={
        'id_col': 'id',
        'seq_col': 'sequence',
    },
)
to_test = acquisition_function.filter(library)
to_test
>>> VariantLibrary(size=10, parent=Variant('MAGV'), fixed_size=True)
```

## 2. Running DE campaigns

### 2.1 Running a single round of a selection step
```python
from aide import VariantLibrary
from aide.rounds import GPSelectionRound
from aide.campaign.database import CampaignDatabase

library = VariantLibrary.load_from_file(
    'library.csv',
    schema={
        'id_col': 'id',
        'seq_col': 'sequence',
    },
)

database = CampaignDatabase('path_to_db')
database.save_library(library)

# assumes you already have a saved GP model of base type AcquisitionFunction
# from sklearn API
# otherwise pass 
round = GPSelectionRound(
    database=database,
    pretrained_acquisition_model='path_to_model_save',
    number_of_variants=10)
round.library_generation()
round.library_selection()
>>> VariantLibrary(size=10, parent=Variant('MAGV'), fixed_size=True)
```
PICKUP HERE

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
from aide import MutationLibrary
from aide.zoo import Saito2018Runner, Saito2018Config

config = Saito2018Config(
    name='my_campaign',
    description='my campaign description',
    database='./my_campaign.db'
    user_input_every_step=False, # if False, status updates that do not require user input are automatically executed
    number_of_rounds=2,
    variants_per_round=10,
) # initial library determined by saturation mutagenesis

runner = Saito2018Runner(config)
# Here 4 residues are modified. This includes more than just single point mutations
# The authors provide a csv file with the 4 residues as a string eg 'GVAY'
# it is assumed that the file has been updated to mark mutations as full specification
# eg. 'S65G;S72V;H77A;T203Y'
# let's say size of library is 10000
initial_library = MutationLibrary.load_from_file(
    'initial_library.csv',
    schema={
        'id_col': 'id',
        'mutation_col': 'mutation',
    }
)
initial_library = initial_library.apply_parent(parent='...') # whatever it is

runner.step()
>>> CampaignError: Cannot take step of round 0, status is 'not ready'
runner.step(library=initial_library)
>>> Round 0 of type RandomSelectionRound() status is 'ready'
>>> Round 0 of type RandomSelectionRound() status is 'generated', putative library contains 10000 variants.
>>> Round 0 of type RandomSelectionRound() status is 'selected', library for experiment has 10 variants at file path 'library_for_exp.csv' # this includes the ids of the variants to test
# after testing them and creating a csv file with the resulting experimental labels
runner.step(label_path='experimental_results.csv', id_col='id', name_col='label_name', value_col='label_value')
>>> Round 0 of type RandomSelectionRound() status is 'labeled'
>>> Round 0 of type RandomSelectionRound() status is 'complete' and is immutable in database.
>>> Round 1 of type GPSelectionRound(features='ResiduePhysicalFeatures', strategy='InfoMax', number_of_variants=10) status is 'ready'

# start next round
runner.step()
>>> Round 1 of type GPSelectionRound(features='ResiduePhysicalFeatures', strategy='InfoMax', number_of_variants=10) status is 'generated', putative library contains 9990 variants.
>>> Round 1 of type GPSelectionRound(features='ResiduePhysicalFeatures', strategy='InfoMax', number_of_variants=10) status is 'selected', library for experiment has 10 variants at file path 'library_for_exp.csv'
# after testing them and creating a csv file with the resulting experimental labels
runner.step(label_path='experimental_results.csv', id_col='id', name_col='label_name', value_col='label_value')
>>> Round 1 of type GPSelectionRound(features='ResiduePhysicalFeatures', strategy='InfoMax', number_of_variants=10) status is 'labeled'
>>> Round 1 of type GPSelectionRound(features='ResiduePhysicalFeatures', strategy='InfoMax', number_of_variants=10) status is 'complete' and is immutable in database.
>>> Round 2 of type GPSelectionRound(features='ResiduePhysicalFeatures', strategy='greedy', number_of_variants=10) status is 'ready'

# start next round
runner.step()
>>> Round 2 of type GPSelectionRound(features='ResiduePhysicalFeatures', strategy='greedy', number_of_variants=10) status is 'generated', putative library contains 9980 variants.
>>> Round 2 of type GPSelectionRound(features='ResiduePhysicalFeatures', strategy='greedy', number_of_variants=10) status is 'selected', library has 10 variants at file path 'library_for_exp.csv'
# the final 10 are selected greedily
runner.step(label_path='experimental_results.csv', id_col='id', name_col='label_name', value_col='label_value')
>>> Round 2 of type GPSelectionRound(features='ResiduePhysicalFeatures', strategy='greedy', number_of_variants=10) status is 'labeled'
>>> Round 2 of type GPSelectionRound(features='ResiduePhysicalFeatures', strategy='greedy', number_of_variants=10) status is 'complete' and is immutable in database.

# you can then select the best variant
runner.get_best_variant()
```

Here is another one, eg. Fox et. all which operates on a per mutation basis instead of a persequence basis

```python
from aide import Library
from aide.zoo import Fox2018Runner, Fox2018Config

config = Fox2018Config(
    name='my_campaign',
    description='my campaign description',
    database='./my_campaign.db'
    user_input_every_step=False, # if False, status updates that do not require user input are automatically executed
    number_of_rounds=2,
    n_libraries_per_round=1
    n_mutations_per_library=15,
    max_prosar_mutations_per_library=10, # the rest are randomly included
    n_variants_per_library=None,
    p_critereon=1e-4,
    max_parent_adjustments=5,
    rule_out_percentile=0.25,
) 
# this library is of a bunch of mutations
# Fox et al consider the process on a per mutation basis
# the model is used to determine which mutations to either apply to the backbone,
# keep in a pool for future rounds, or rule out
# they did mutliple "campaigns" in paralle but this is out of scope for now.
# instead, do one linear run where the starting pool of mutations is all mutations
# and we wittle it down.
runner = Fox2018Runner(config)

initial_library = MutationLibrary.load_from_file(
    'initial_library.csv',
    schema={
        'id_col': 'id',
        'mutation_col': 'mutation',
    }
)
initial_library = initial_library.apply_parent(parent='...') # whatever it is

# start by loading the initial library
runner.step(library=initial_library)
>>> Round 0 of type RandomCombinationsRound(n_libraries_per_round=1, n_mutations_per_library=15, n_variants_per_library=None) status is 'ready', 
>>> Round 0 of type RandomCombinationsRound(n_libraries_per_round=1, n_mutations_per_library=15, n_variants_per_library=None) status is 'generated', putative combinatorial library contains 15! variants. Mutations for combinations and parent sequence are at file path 'mutations_for_combinations.txt'
runner.step()
>>> CampaignError: Cannot take step of round 0, status is 'not started'
# we have to load the phenotypes for the sample of the combinatorial library that were created
combinatorial_library = MutationLibrary.load_from_file(
    'combinatorial_library.csv', # assuming we created 1000 variants
    schema={
        'id_col': 'id',
        'mutation_col': 'mutation',
    }
)
runner.step(library=combinatorial_library) # it knows the parent for the round so we can convert to VariantLibrary internally
>>> Round 0 of type RandomCombinationsRound(n_libraries_per_round=1, n_mutations_per_library=15, n_variants_per_library=None) status is 'selected', library for experiment has 1000 variants at file path 'library_for_exp.csv' # note that since n_variants_per_library=None, the entire library is tested
# after testing them and creating a csv file with the resulting experimental labels
runner.step(label_path='experimental_results.csv', id_col='id', name_col='label_name', value_col='label_value')
>>> Round 0 of type RandomCombinationsRound(n_libraries_per_round=1, n_mutations_per_library=15, n_variants_per_library=None) status is 'labeled'
>>> Round 0 of type RandomCombinationsRound(n_libraries_per_round=1, n_mutations_per_library=15, n_variants_per_library=None) status is 'complete' and is immutable in database.
>>> Round 1 of type ProSARCombinationsRound(n_libraries_per_round=1, n_mutations_per_library=15, n_prosar_mutations_per_library=10, n_variants_per_library=10, p_critereon=1e-4, max_parent_adjustments=5, rule_out_percentile=0.25) status is 'ready'

# start next round
runner.step()
>>> Round 1 of type ProSARCombinationsRound(n_libraries_per_round=1, n_mutations_per_library=15, n_prosar_mutations_per_library=10, n_variants_per_library=10, p_critereon=1e-4, max_parent_adjustments=5, rule_out_percentile=0.25) status is 'generated', putative combinatorial library contains 15! variants. Mutations for combinations and parent sequence are at file path 'mutations_for_combinations.txt'
# For this round, the parent is the best performer from the previous round.
# The mutations tested in the previous round are considered. All labeled data are included, and featurized
# on those mutations. From the analysis, the mutations from the previous round are either fixed if already present
# in the new parent, included in the next library, or ruled out. Random mutations are also added
combinatorial_library = MutationLibrary.load_from_file(
    'combinatorial_library.csv', # assuming we created 1000 variants
    schema={
        'id_col': 'id',
        'mutation_col': 'mutation',
    }
)
runner.step(library=combinatorial_library) # it knows the parent for the round so we can convert to VariantLibrary internally
>>> Round 1 of type ProSARCombinationsRound(n_libraries_per_round=1, n_mutations_per_library=15, n_prosar_mutations_per_library=10, n_variants_per_library=10, p_critereon=1e-4, max_parent_adjustments=5, rule_out_percentile=0.25) status is 'selected', library for experiment has 10 variants at file path 'library_for_exp.csv' # note that since n_variants_per_library=10, only the top 10 are tested
# after testing them and creating a csv file with the resulting experimental labels
runner.step(label_path='experimental_results.csv', id_col='id', name_col='label_name', value_col='label_value')
>>> Round 1 of type ProSARCombinationsRound(n_libraries_per_round=1, n_mutations_per_library=15, n_prosar_mutations_per_library=10, n_variants_per_library=10, p_critereon=1e-4, max_parent_adjustments=5, rule_out_percentile=0.25) status is 'labeled'
>>> Round 1 of type ProSARCombinationsRound(n_libraries_per_round=1, n_mutations_per_library=15, n_prosar_mutations_per_library=10, n_variants_per_library=10, p_critereon=1e-4, max_parent_adjustments=5, rule_out_percentile=0.25) status is 'complete' and is immutable in database.

# you can then select the best variant
runner.get_best_variant()
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






