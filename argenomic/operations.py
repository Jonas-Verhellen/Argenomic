import hydra
import random
import logging
import numpy as np
import pandas as pd
from typing import List, Tuple

from rdkit import Chem
from rdkit import rdBase
rdBase.DisableLog('rdApp.error')

from rdkit.Chem import AllChem
from rdkit.Chem import rdMMPA

class mutator:
    """
    A catalog class containing and implementing mutations to small molecules
    according to the principles of positional analogue scanning.
    """
    def __init__(self) -> None:
        self.mutation_data = pd.read_csv(hydra.utils.to_absolute_path("data/smarts/mutation_collection.tsv"), sep='\t')

    def __call__(self, molecule:Chem.Mol) -> List[Chem.Mol]:
        sampled_mutation = self.mutation_data.sample(n=1, weights='probability').iloc[0]
        reaction = AllChem.ReactionFromSmarts(sampled_mutation['smarts'])
        try:
            molecules = [products[0] for products in reaction.RunReactants([molecule])]
        except:
            molecules = []
        return molecules

class crossover:
    """
    A strategy class implementing a parent-centric crossover of small molecules.
    """
    def __init__(self):
        pass

    def __call__(self, molecule_pair:Tuple[Chem.Mol, Chem.Mol]) -> List[Chem.Mol]:
        molecule_cores, molecule_sidechains = self.fragmentate(molecule_pair)
        molecules = self.merge(molecule_cores, molecule_sidechains)
        return molecules

    def merge(self, molecule_cores:List[Chem.Mol], molecule_sidechains:List[Chem.Mol]) -> List[Chem.Mol]:
        molecules = []
        random.shuffle(molecule_sidechains)
        reaction = AllChem.ReactionFromSmarts('[*:1]-[1*].[1*]-[*:2]>>[*:1]-[*:2]')
        for core, sidechain in zip(molecule_cores, molecule_sidechains):
            molecules.append(reaction.RunReactants((core, sidechain))[0][0])
        return molecules

    def fragmentate(self, molecule_pair:Tuple[Chem.Mol, Chem.Mol]) -> Tuple[List[Chem.Mol], List[Chem.Mol]]:
        molecule_cores = []
        molecule_sidechains = []
        for molecule in molecule_pair:
            molecule_frags = rdMMPA.FragmentMol(molecule, maxCuts=1, resultsAsMols=False)
            if len(molecule_frags) > 0:
                _, molecule_frags = map(list, zip(*molecule_frags))
                for molecule_pair in molecule_frags:
                    core, sidechain = molecule_pair.split(".")
                    molecule_cores.append(Chem.MolFromSmiles(core.replace("[*:1]", "[1*]")))
                    molecule_sidechains.append(Chem.MolFromSmiles(sidechain.replace("[*:1]", "[1*]")))
        return molecule_cores, molecule_sidechains
