# Unifying frameworkd for AI assisted Directed Evolution (AIDE)
An attempt to create a modular, unifying code base for incorporating recent (and future) AI techniques to accelerate directed evolution campaigns. 

## Motivation

Directed evolution campaigns leverage "survival of the fittest" in an artifical setting in order to drive a protein to target properties such as higher stability, activity, etc. This process involves multiple rounds of mutation and evaluation. The variants must be evaluated and screened. A number of techniques exist to express and screen the variants or otherwise select for the target property, descibed elsewhere[]. Regardless, screening in the lab is expensive and time consuming, and it is desired to limit the number of evaluations as much as possible while still achieving an increase in fitness. __To this end, a number of different strategies have been proposed recently to accelerate this process with machine learning__. 

These strategies broadly fall into one or more of the following categories:

1. In silico estimators of the desired property/ies, used to filter (predicted) poor performing variants and reduce the number of evaluations. They are sometimes "improved" by further training after each round.
2. In silico estimators of the desired property/ies with uncertainty estimation act as aquisition functions that balance exploration of fitness space with exploitation
3. Library generation tools. The fitness landscape is infinately large and cannot be explored in full, even with in silico estimators, and these tools produce a tractable library size to explore.

__While these recent applications are all somewhat unique, I believe that they all fit within the above framework__. This library attempts to define a framework into which these mmodern tools can be aggregated.

## Getting started

TODO

## Summary of the framework

TODO

## References

TODO

## License