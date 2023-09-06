# Catalog of recent ML + DE publications

Design considerations from the methods unto the framework are noted.

## Publications
### 1. Fox, R. J. et al. Improving catalytic function by ProSAR-driven enzyme evolution. Nat Biotechnol 25, 338–344 (2007).

- __Summary__: The authors utilize a number of techniques to generate variation (saturation, random, inclusion of close homologs, and gene shuffling) and from the generated variation the indentify each possible mutation to design a combinatorial library. The use semi synthetic oligonucleotides to sample from permutations of the library to screen. A linear PLAS model is trained with a single coefficient for each possible binary mutation eg. A123T has a different coefficient than A123G. After each round, large coefficients from the model that are sufficiently confident (over multiple variants with that mutation) are used to filter out bad or neutral mutations. Beneficial mutations are fixed and cycled back into the library, along with some fresh variation.
- __Design considerations__: 
    1. This method does not filter variants to test, it saturates the library with good mutations. This suggests that "estimators" may need to be used in the context of variant generation, not only aquisition functions. 
    2. This method would require one to be able to join libraries. Eg. the user inputs a library of variants with trained scores, A library ot mutations is produced by the PLAS estimator. This library needs to be joined with new varitation like random mutagenesis


### Saito, Y. et al. Machine-Learning-Guided Mutagenesis for Directed Evolution of Fluorescent Proteins. ACS Synth. Biol. 7, 2014–2022 (2018).
- __Summary__: An initial random library is devleloped via point saturation mutagenesis and random mutagensis at four residue positions, which are used to train a GP. Used physicochemical residue features (residue wise features) concatenated to train a GP (COMBO) and used probability-of-improvement as an aquisition function. Eg. assuming guassian, the cumulative mass function of performance greater than the current step.
- __Design considerations__: 
    1. We should have an abstract class to represent aquisition functions, and take at a minimum, predictions as input, and maybe other parameters. They use these parameters to output a score. These should likely wrap around "estimators" and should check that the estimator can output the necessary parameters. For example, a probability-of-improvement aquisition function relies on the internal estimator outputing standard deviation as well as prediction
    2. We need to represent an "encoder" class that converts proteins or mutations sets into features, as well as the "estimators" class.
    3. Since some methods rely on a fixed set of residues being explored, while others suggest full proteins, we must differentiate between mutation libraries and variant libraries. Further, we must differentiate between mutation set encoders and variant encoders.
    4. We may not have control over the libraries for mutagensis methods, eg we know which residues can change, but the changes are done randomly and we later sequence. This suggests that we need a method that produces libraries by data input after sequencing.

### 1. Wu, Z., Kan, S. B. J., Lewis, R. D., Wittmann, B. J. & Arnold, F. H. Machine learning-assisted directed protein evolution with combinatorial libraries. Proc Natl Acad Sci U S A 116, 8852–8858 (2019).
- __Summary__: Used ML on one-hot encoded mutations over a combinatorial library. Randomly screened some variants from the library. With trained model, find top 1000 performing variants by prediction over combinatorial space, and compute the AA distributions over each position from thos to develop degenerate codon scheme to sample from new variants. Note this differs from directly selection top N predicted variants because those variants would have to be fully synthesized. This possible but more expensive than random site directed mutagensis with degenerate oligonucleotides.
__Design considerations__:
    1. While the model was able to filter possible variants from the combinatorial library, the authors went a cheaper route and did some analysis over the distribution. This suggests that We may need some sort of analysis component that acts on a library. Eg. to fit within the framework, the aquisition function in this case is simply a top N greedy model filter, outputing a filtered Library of variants, but maybe we want an optional downstream component to analyze a Library and extract some information. In this case, distribution over mutations. May even take it one step further and allow Mutation class to contain some sort of statistics, so the framework could output a Combinatorial library with statistics. For this use case: Aquisiton function outputs a VariantLibrary, and the Analysis component can convert it to a CombinatorialLibrary with each Mutation having statistics. 




## Design notes

- The framework will require that ML estimators be used in the context of variant generation as well as aquisition functions.
- We should center the framework around a "Library" class that records "Variants" and "Mutations". A variant maps a full sequence. Variants can be directly defined or defined by combining Mutations. Variant generation steps will produce Libraries. Aquisition function methods act on Libraries to produce smaller Libraries.
    - Mutation class must track wild type sequence (a Variant object), position and change of amino acids
    - A Variant can be generated from the Mutation
    - Variants can be generated combinatorially from mutations
    - Variants should track (if any) Mutations that resulted in them as well as the original Variant
- There will have to be a few subclasses of Library, eg. VariantLibrary that is just a list of variants, CombinatorialLibrary which is defined by a set of mutations
- Variant generation should probably have a few types, dictated by the experimental methods. Some methods do not have full control over the specific set of mutations and are instead random. For the random case, the variation is manually input to the framework producing a manual Library for downstream aquisition functions. For the experimental setups that allow us to express a specific protein, we can have Libraries defined in silico by the framework, eg. by random mutagenesis, recombination, deep learning generation, or model filtered mutations.
- We should be writing things to a duck db database and libraries should be able to be saved and loaded to the database. So should variant scores. We must somehow also track which round we are in.