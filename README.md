<div align="center">    
 
# Quality-Diversity Optimisation for Molecular Design
   ![Logo](/data/figures/logo.png "Logo")

</div>

[![GitHub issues](https://img.shields.io/github/issues/Jonas-Verhellen/argenomic)](https://github.com/Jonas-Verhellen/argenomic/issues)
![GitHub](https://img.shields.io/github/license/Jonas-Verhellen/Argenomic)
[![DOI](https://img.shields.io/badge/DOI-10.1039/D0SC03544K-blue)](https://pubs.rsc.org/en/content/articlehtml/2020/sc/d0sc03544k)


## Description

Argenomic is an open-source implementation of an illumination algorithm for optimization of small organic molecules. Argenomic provides a holistic overview of how high-performing molecules are distributed throughout a search space. This novel approach produces potent but qualitatively different molecules, illuminates the distribution of optimal solutions, and improves search efficiency compared to both machine learning and traditional genetic algorithm approaches. This implementation is based on an open-source, [graph-based genetic algorithm](https://github.com/jensengroup/GB-GA) for molecular optimisation, and influenced by state-of-the-art concepts from [soft robot design](https://github.com/resibots/pymap_elites). For more information, see the accompanying [blog post](https://jonas-verhellen.github.io/posts/2020/07/argenomic/).

<p align="center">
  <img src="https://github.com/Jonas-Verhellen/jonas-verhellen.github.io/blob/main/images/video.gif"/>
</p>

*Example of iterative illumination in a 2D representation of chemical space. In this case the fitness is determined as the molecular similarity to Troglitazone.*

## Getting Started

After installing the software and running the tests, a basic usage example of argenomic (i.e. the rediscovery of Troglitazone) can be called upon in the following manner:
```
python3 illuminate.py generations=100
```

### Installing

Download the source code from Github to your local machine and create the environment from the environment.yml file:
```
conda env create -f environment.yml
```
Activate the new environment:
```
conda activate argenomic-stable
```
Verify that the new environment was installed correctly:
```
conda env list
```

### Running the tests

To run the unit tests:

```
pytest ./tests
```

## Authors
Based on the paper *[Illuminating elite patches of chemical space.](https://pubs.rsc.org/en/content/articlehtml/2020/sc/d0sc03544k) Chemical science 11.42 (2020): 11485-11491.*
* **Jonas Verhellen** - *Concept, implementation, and development*
* **Jeriek Van den Abeele** - *implementation, and development*

## Dependencies

Important dependencies of the Argenomic software environment and where to find the source.

* [Python](https://www.python.org/) - Python is a widely used scientific and numeric programming language.
* [RDKit](https://github.com/rdkit/rdkit) - Cheminformatics and machine-learning software toolkit.
* [Scikit-learn](https://github.com/scikit-learn/scikit-learn) - Data science and deep learning toolset in Python.
* [Omegaconf](https://github.com/omry/omegaconf) - Configuration system for multiple sources, providing a consistent API.

## Acknowledgments

* Jan Jensen for his work in developing and open-sourcing a graph-based genetic algorithm for molecular optimisation, which served as impetus for this project.

* Jean-Baptiste Mouret and Jeff Clune for their breakthrough invention of illumination algorithms, providing a holistic view of high-performing solutions throughout a search space.

* Pat Walters for his scripts indicating how to run structural alerts using the RDKit and ChEMBL, and for his many enlightening medicinal chemistry blog posts.

## Copyright License

This project is licensed under the GNU Affero General Public License v3.0 (AGPLv3).
