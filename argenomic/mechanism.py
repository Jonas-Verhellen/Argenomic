import sys
import numpy as np
import pandas as pd
from typing import List, Tuple

from rdkit import Chem
from rdkit import rdBase
from rdkit import RDConfig
rdBase.DisableLog('rdApp.error')

from rdkit.Chem import AllChem
from rdkit.Chem import Crippen
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity

class Descriptor:
    """
    A strategy class for calculating the descriptor vector of a molecule.
    """
    def __init__(self, config_descriptor) -> None:
        self.properties = []
        self.ranges = config_descriptor.ranges   
        self.property_names = config_descriptor.properties
        for name in self.property_names:
            module, fuction = name.split(".")
            module = getattr(sys.modules[__name__], module)
            self.properties.append(getattr(module, fuction))
        return None

    def __call__(self, molecule) -> None:
        """
        Updates the descriptor vector of a molecule.
        """
        descriptor = []
        molecular_graph = Chem.MolFromSmiles(molecule.smiles)
        for property, range in zip(self.properties, self.ranges):
            descriptor.append(self.rescale(property(molecular_graph), range))
        molecule.descriptor = descriptor
        return molecule

    @staticmethod
    def rescale(feature: List[float], range: List[float]) -> List[float]:
        """
        Rescales the feature to the unit range.
        """
        rescaled_feature = (feature - range[0])/(range[1] - range[0])
        return rescaled_feature

class Fitness:
    """
    A strategy class for calculating the fitness of a molecule.
    """
    def __init__(self, config_fitness) -> None:
        self.fingerprint_type = config_fitness.type
        self.target = Chem.MolFromSmiles(config_fitness.target)
        self.target_fingerprint = self.get_fingerprint(self.target, self.fingerprint_type)
        return None

    def __call__(self, molecule) -> None:
        """
        Updates the fitness value of a molecule.
        """
        molecular_graph = Chem.MolFromSmiles(Chem.CanonSmiles(molecule.smiles))
        molecule_fingerprint = self.get_fingerprint(molecular_graph, self.fingerprint_type)
        fitness = TanimotoSimilarity(self.target_fingerprint, molecule_fingerprint)
        molecule.fitness = fitness
        return molecule

    def get_fingerprint(self, molecular_graph: Chem.Mol, fingerprint_type: str):
        method_name = 'get_' + fingerprint_type
        method = getattr(self, method_name)
        if method is None:
            raise Exception('{} is not a supported fingerprint type.'.format(fingerprint_type))
        return method(molecular_graph)

    def get_ECFP4(self, molecular_graph: Chem.Mol):
        return AllChem.GetMorganFingerprint(molecular_graph, 2)

    def get_ECFP6(self, molecular_graph: Chem.Mol):
        return AllChem.GetMorganFingerprint(molecular_graph, 3)

    def get_FCFP4(self, molecular_graph: Chem.Mol):
        return AllChem.GetMorganFingerprint(molecular_graph, 2, useFeatures=True)

    def get_FCFP6(self, molecular_graph: Chem.Mol):
        return AllChem.GetMorganFingerprint(molecular_graph, 3, useFeatures=True)

