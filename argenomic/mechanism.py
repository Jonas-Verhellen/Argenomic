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

class descriptor:
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

    def __call__(self, molecule: Chem.Mol) -> List[float]:
        """
        Calculating the descriptor vector of a molecule.
        """
        descriptor = []
        for property, range in zip(self.properties, self.ranges):
            descriptor.append(self.rescale(property(molecule), range))
        return descriptor

    @staticmethod
    def rescale(feature: List[float], range: List[float]) -> List[float]:
        """
        Rescaling the feature to the unit range.
        """
        rescaled_feature = (feature - range[0])/(range[1] - range[0])
        return rescaled_feature

class fitness:
    """
    A strategy class for calculating the fitness of a molecule.
    """
    def __init__(self, config_fitness) -> None:
        self.memoized_cache = dict()
        self.fingerprint_type = config_fitness.type
        self.target = Chem.MolFromSmiles(config_fitness.target)
        self.target_fingerprint = self.get_fingerprint(self.target, self.fingerprint_type)
        return None

    def __call__(self, molecule: Chem.Mol) -> float:
        smiles = Chem.MolToSmiles(molecule)
        if smiles in self.memoized_cache:
            fitness = self.memoized_cache[smiles]
        else:
            molecule_fingerprint = self.get_fingerprint(molecule, self.fingerprint_type)
            fitness = TanimotoSimilarity(self.target_fingerprint, molecule_fingerprint)
            self.memoized_cache[smiles] = fitness
        return fitness

    def get_fingerprint(self, molecule: Chem.Mol, fingerprint_type: str):
        method_name = 'get_' + fingerprint_type
        method = getattr(self, method_name)
        if method is None:
            raise Exception('{} is not a supported fingerprint type.'.format(fingerprint_type))
        return method(molecule)

    def get_ECFP4(self, molecule: Chem.Mol):
        return AllChem.GetMorganFingerprint(molecule, 2)

    def get_ECFP6(self, molecule: Chem.Mol):
        return AllChem.GetMorganFingerprint(molecule, 3)

    def get_FCFP4(self, molecule: Chem.Mol):
        return AllChem.GetMorganFingerprint(molecule, 2, useFeatures=True)

    def get_FCFP6(self, molecule: Chem.Mol):
        return AllChem.GetMorganFingerprint(molecule, 3, useFeatures=True)
