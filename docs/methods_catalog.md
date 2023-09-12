# Catalog of recent ML + DE publications

Works that help define the contours of the space are retained. Work that is either not seminal or is captured by the other work is not given. Design considerations from the methods unto the framework are noted.

## Learning the fitness landscape publications

### Romero, P. A., Krause, A. & Arnold, F. H. Navigating the protein fitness landscape with Gaussian processes. Proceedings of the National Academy of Sciences 110, E193–E201 (2013).
- __Summary__: Authors show that GPs can predict arbitrary protein properties from sequences, and the uncertainty from the model can be used to explore holes in the fitness landscape.

### Hopf, T. A. et al. Mutation effects predicted from sequence co-variation. Nat Biotechnol 35, 128–135 (2017).
- __Summary__: MSA model (trained per each MSA) predicts liklihood of mutations incorporating pairwise interactions. "EVMutation"
- __Design Considerations__: 
    1. See the DeepSequence notes below

### Riesselman, A. J., Ingraham, J. B. & Marks, D. S. Deep generative models of genetic variation capture the effects of mutations. Nat Methods 15, 816–822 (2018).
- __Summary__: MSA model (trained per each MSA) predicts liklihood of mutations incorporating high order interaction effects by relating the mutations occuring to latent distribution. Incorporates some bio structure into the architecture/training. "DeepSequence"
- __Design Considerations__: 
    1. Can be used to evaluate mutations zero shot, as was done to design training set in ftMLDE. Could be wrapped into an Estimator
    2. The model can only provide information for specific positions within the WT protein, variants are described by mutation strings and breaks for positions outside of the focus of the MSA. In our framework, to use a tool like this, variants must be able to be described as a set of positonal mutation without indels. We should have a check for this for some models. Eg. have properties of the Variant class that indicate if it has indels.

### Riesselman, A. et al. Accelerating Protein Design Using Autoregressive Generative Models. 757252 Preprint at https://doi.org/10.1101/757252 (2019).
- __Summary__: Family model (trained per each Family) predicts liklihood of variants incorporating high order interaction effects by modeling the causal probability of AA strings in the family. Eg. Train causal language model (Not a transformer, comparitavely low parameter) to a protein family. probability of a variant is the product probability of all causal token probs. Differs from their DeepSequence work in that it is not a position/mutation based model and can handle indels.
- __Design Considerations__: 
    1. Can be used to evaluate mutations zero shot, as was done to design training set in ftMLDE. Could be wrapped into an Estimator



### HSU combine models

## DE Publications
### Fox, R. J. et al. Improving catalytic function by ProSAR-driven enzyme evolution. Nat Biotechnol 25, 338–344 (2007).
- __Summary__: The authors utilize a number of techniques to generate variation (saturation, random, inclusion of close homologs, and gene shuffling) and from the generated variation the indentify each possible mutation to design a combinatorial library. The use semi synthetic oligonucleotides to sample from permutations of the library to screen. A linear PLAS model is trained with a single coefficient for each possible binary mutation eg. A123T has a different coefficient than A123G. After each round, large coefficients from the model that are sufficiently confident (over multiple variants with that mutation) are used to filter out bad or neutral mutations. Beneficial mutations are fixed and cycled back into the library, along with some fresh variation.
- __Design considerations__: 
    1. This method does not filter variants to test, it saturates the library with good mutations. This suggests that "estimators" may need to be used in the context of variant generation, not only aquisition functions. 
    2. This method would require one to be able to join libraries. Eg. the user inputs a library of variants with tested scores, A library of mutations is produced by the PLAS estimator. This library needs to be joined with other mutations from variation like random mutagenesis

### Saito, Y. et al. Machine-Learning-Guided Mutagenesis for Directed Evolution of Fluorescent Proteins. ACS Synth. Biol. 7, 2014–2022 (2018).
- __Summary__: An initial random library is developed via point saturation mutagenesis and random mutagensis at four residue positions, which are used to train a GP. Used physicochemical residue features (residue wise features) concatenated to train a GP (COMBO) and used probability-of-improvement as an aquisition function. Eg. assuming guassian, the cumulative mass function of performance greater than the current step.
- __Design considerations__: 
    1. We should have an abstract class to represent aquisition functions, and take at a minimum, predictions as input, and maybe other parameters. They use these parameters to output a score. These should likely wrap around "estimators" and should check that the estimator can output the necessary parameters. For example, a probability-of-improvement aquisition function relies on the internal estimator outputing standard deviation as well as mean prediction.
    2. We need to represent an "encoder" class that converts proteins or mutations sets into features, as well as the "estimators" class.
    3. Since some methods rely on a fixed set of residues being explored, while others suggest full proteins, we must differentiate between mutation libraries and variant libraries. Further, we must differentiate between mutation set encoders and variant encoders.
    4. We will not have control over the libraries most of the time for eg. mutagensis methods, eg we know which residues can change, but the changes are done randomly and we later sequence. This suggests that we need a method that produces libraries by data input after sequencing.

### Wu, Z., Kan, S. B. J., Lewis, R. D., Wittmann, B. J. & Arnold, F. H. Machine learning-assisted directed protein evolution with combinatorial libraries. Proc Natl Acad Sci U S A 116, 8852–8858 (2019).
- __Summary__: Used ML on one-hot encoded mutations over a combinatorial library. Randomly screened some variants from the library. With trained model, find top 1000 performing variants by prediction over combinatorial space, and compute the AA distributions over each position from this to develop degenerate codon scheme to sample from new variants. Note this differs from directly selection top N predicted variants because those variants would have to be fully synthesized. This possible but more expensive than random site directed mutagensis with degenerate oligonucleotides. "MLDE" Authors use a number of model architectures with cross validation, picking the top ensemble of models based on CV score for actual protein filtering.
- __Design considerations__:
    1. While the model was able to filter possible variants from the combinatorial library, the authors went a cheaper route and did some analysis over the distribution. This suggests that We may need some sort of analysis component that acts on a library. Eg. to fit within the framework, the aquisition function in this case is simply a top N greedy model filter, outputing a filtered Library of variants, but maybe we want an optional downstream component to analyze a Library and extract some information. In this case, distribution over mutations. May even take it one step further and allow Mutation class to contain some sort of statistics, so the framework could output a Combinatorial library with statistics. For this use case: Aquisiton function outputs a VariantLibrary, and the Analysis component can convert it to a CombinatorialLibrary with each Mutation having statistics. 
    2. Ensembles of models should be considered, eg. a protocol to test out a bunch of models in CV and make an ensemble of the best should be an "Estimator" in the library.
- __Datasets used__: GB1

### Wittmann, B. J., Yue, Y. & Arnold, F. H. Informed training set design enables efficient machine learning-assisted directed protein evolution. Cell Syst 12, 1026-1045.e7 (2021).
- __Summary__: Extends strategy of MLDE from Wu. et al. above by using better-than-random initial data sampling for experimenting. Zero shot models are used to identify likely functional proteins for experiment, avoiding taking data on broken proteins. A number of zero shot models are included. They also explore the effect of encodings on the supervised model. Primary takeaway is that using Zero shot models to design the training set is critical, regardless of downstream SML architecture/embeddings.
- __Design Considerations__:
    1. Code included. Not sure how clean it is.
    2. This does not intriduce anything new to our framework. It is just a novel aquisition function for round 0 instead of random.

### Qiu, Y., Hu, J. & Wei, G.-W. Cluster learning-assisted directed evolution. Nat Comput Sci 1, 809–818 (2021).
- __Summary__: The authors incorporate unsupervised clustering into sample selection as follows. Variants (encoded) are clustered. Optionally cluster further within cluster. An aquisition function is used over clusters, then over sequences, to select a round for experiments. Initialliy this is random. After first round, clusters are selected with probability proportional to average performance of tested variants. Then, sequences are selected from cluster using eps-greedy or UCB etc. After some number of rounds, a supervised model is trained and the last round of experiment is greedily selected from the global library. They compared this method to global strategies eg. MLDE and GP. They can also us Zero shot predictors to aid with initial random sampling like ftMLDE, which is global but starts with a non-random training set.
- __Design Considerations__: 
    1. They have code, hopefully most of the functionality can be imported.
    2. If not, we need to have the Aquisition function be highly modular b/c it involves multiple clustering models. 

### Emami, P., Perreault, A., Law, J., Biagioni, D. & St. John, P. Plug & play directed evolution of proteins with gradient-based discrete MCMC. Mach. Learn.: Sci. Technol. 4, 025014 (2023).
- __Summary__: Capitalizing on recent work (Hsu et al.) showing that unsupervised evolutionary probability models can be combined with supervised models to increase predictive performance, the authors propose a DE mutation sampler that can explore greater than single point change paths. The sampler relies on an arbitrary product of expert models (POE) contingent that the models are differentiable. For some random max distance R (eg 5 mutations away) compute POE derivative, use to sample a single AA (where AA pointing more steeply to POE optima more likely), repeat to R. Keep variant with liklihood propto product of increase in POE liklihood compared to parent.
- __Design Considerations__: 
    1. We need an operation that can combine estimators, check they are differentiable, then run the sampler to output a library.
    2. This method is general in that it could work to sample from random positions, but also a constrained set. Thus the method should at minimum take the parent seq and the expert models as input, but could also take a mutation set or list of positions to constrain the sampling.
    3. Inputs to all models must be one-hot encoded to ensure the categorical softmax of derivatives directly describes an AA at a position.


### LINDER



## Design notes
- Some strategies are pure active learning strategies that use the same algorithm for each round, but most have multiple stages: exploration the exploitation. Thus suggests that when we rebuild the above strategies in our framework, we need to build multiple rounds each pulling from the same pool of components but maybe differing between rounds. eg. we might have an ftMLDEPipeline component that defines the estimators and aquisition functions for a single round based on inpuit parameters. Then we have A runner run multiple.
- The framework will require that ML estimators be used in the context of variant generation as well as aquisition functions.
- We should center the framework around a "Library" class that records "Variants" and "Mutations". A variant maps a full sequence. Variants can be directly defined or defined by combining Mutations. Variant generation steps will produce Libraries. Aquisition function methods act on Libraries to produce smaller Libraries.
    - Mutation class must track wild type sequence (a Variant object), position and change of amino acids
    - A Variant can be generated from the Mutation or set of mutations
    - Variants libraries can be generated combinatorially from mutations
    - Variants should track (if any) Mutations that resulted in them as well as the original Variant
    - Variants can be of type "Mutational" where it is exactly described by a set of mutations from another variant or "Abstract" where it is just a sequence. We can try to intuit a Mutational Variant from abstratc by eg. alignment, but it will be inexact.
- There will have to be a few subclasses of Library, eg. VariantLibrary that is just a list of variants, CombinatorialLibrary which is defined by a set of mutations
- Variant generation should probably have a few types, dictated by the experimental methods. Most methods do not have full control over the specific set of mutations and are instead random. For the random case, the variation is manually input to the framework producing a manual Library for downstream aquisition functions. For the experimental setups that allow us to express a specific protein, we can have Libraries defined in silico by the framework, eg. by random mutagenesis, recombination, deep learning generation, or model filtered mutations.
- We should be writing things to a duck db database and libraries should be able to be saved and loaded to the database. So should variant scores. We must somehow also track which directoed evolution round we are in. This way it is easy to load libraries from each round, as well as the combination of all evaluated variants for model training from all rounds.
- It would be nice to leverage existing frameworks like Sklearn for estimators, transformations, pipelines
- We should be able to define a Pipeline using the modular pieces if the library that represent one directed evolution round. Fully verbose parameterization
- Some high level "Runner" class should handle saving and loading to file over mutliple evolution rounds, as well as calling the pipelines with varying parameters as the experiments procede. Eg. not all methods are exactly the same parameters for each round. Many start off more random, then become more refined as the models learns.


## Benchmark Datasets

- GB1 dataset: 4 site combinatorial lbrary. About 90% empiracally tested. https://www.ncbi.nlm.nih.gov/bioproject/PRJNA278685/
- PhoQ
- Hsu et al.