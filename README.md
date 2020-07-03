# Argenomic
Argenomic is  an open-source implementation of an illumination algorithm for optimization of small organic molecules.

## Getting Started

To run the the argenomic software (rediscovery of Thiotixene):
```
python3 illuminate.py configuration_file=./configuration/config.yaml generations=50
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

## Authors

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

This project is licensed under the GNU General Public License v3 (GPLv3).
