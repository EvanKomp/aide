
# AI-assisted Directed Evolution Unifying Framework (AIDE) Summary
![ci_pass](https://github.com/EvanKomp/aide/blob/main/.github/workflows/ci.yaml/badge.svg) [![codecov](https://codecov.io/gh/EvanKomp/aide/graph/badge.svg?token=JNCIVORJJ6)](https://codecov.io/gh/EvanKomp/aide)

AIDE is a modular software project designed to integrate current and future AI techniques to expedite directed evolution campaigns. These campaigns aim to optimize proteins for specific properties using iterative mutation and evaluation, a process that can be resource-intensive and time-consuming. By incorporating AI methods, AIDE aims to reduce the number of evaluations needed.

The framework integrates several key strategies:

1. **In Silico Estimators**: These tools predict properties to filter out low-performing variants, improving through iterative training.
2. **Uncertainty Estimators**: These act as acquisition functions, balancing the exploration and exploitation of the fitness landscape.
3. **Library Generation Tools**: These tools produce manageable library sizes from vast fitness landscapes.

__While these recent applications are all somewhat unique, I believe that they all breakdown as above__. This library attempts to define a framework into which these modern tools can be aggregated.


# Design Considerations

AIDE is centered around a "Library" class that records "Variants" and "Mutations". A variant can be defined directly or by combining mutations. Variant generation steps produce Libraries, and acquisition function methods act on Libraries to produce smaller Libraries. The project also requires that machine learning estimators be used in the context of variant generation and acquisition functions.

The framework also accommodates different subclasses of Library, such as VariantLibrary (a list of variants) and CombinatorialLibrary (defined by a set of mutations).

Variant generation should be flexible, with a few types based on the experimental methods. Most methods do not have full control over the specific set of mutations and are instead random. For these cases, the variation is manually input to the framework, producing a manual Library for downstream acquisition functions.

AIDE uses a database, to save and load Libraries and track the directed evolution round. Existing frameworks like Sklearn could be leveraged for estimators, transformations, and pipelines. A high-level "Runner" class handles saving and loading to file over multiple evolution rounds, as well as calling the pipelines with varying parameters as the experiments proceed.

## Getting Started
TODO

## Summary of the Framework
TODO

## References
TODO

## License
