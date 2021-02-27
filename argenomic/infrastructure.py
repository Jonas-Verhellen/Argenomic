import os
import csv
import hydra
import random
import itertools

import numpy as np
import pandas as pd
from typing import List, Tuple
from dataclasses import dataclass

from datetime import datetime
from sklearn.cluster import KMeans
from sklearn.neighbors import KDTree

from rdkit import Chem
from rdkit import rdBase
from rdkit.Chem import AllChem
rdBase.DisableLog('rdApp.error')
from rdkit.Chem import Lipinski

from argenomic.base import Molecule, Elite

class Archive:
    """
    A composite class containing the current elite molecules in a CVT tree structure. Allows for processing of 
    new molecules, sampling of the existing elite molecules, and disk storage of the current state of the archive. 
    The CVT centers are either loaded from or deposited to cache disk storage. 
    """
    def __init__(self, archive_config, descriptor_config) -> None:
        self.archive_size = archive_config.size
        self.archive_accuracy = archive_config.accuracy
        self.archive_dimensions = len(descriptor_config.properties)
        self.cache_string = "cache_{}_{}.csv".format(self.archive_dimensions, self.archive_accuracy)
        self.cvt_location = hydra.utils.to_absolute_path("data/cvt/" + self.cache_string)
        if os.path.isfile(self.cvt_location):
            self.cvt_centers = np.loadtxt(self.cvt_location)
        else:
            kmeans = KMeans(n_clusters=self.archive_size)
            kmeans = kmeans.fit(np.random.rand(archive_config.accuracy, self.archive_dimensions))
            self.cvt_centers = kmeans.cluster_centers_
            np.savetxt(self.cvt_location, self.cvt_centers)        
        self.cvt = KDTree(self.cvt_centers, metric='euclidean')
        self.elites = [Elite(index) for index, _ in enumerate(self.cvt_centers, start=0)]
        return None

    def cvt_index(self, descriptor: List[float]) -> int:
        """
        Returns CVT index for the niche nearest to the given discriptor. 
        """
        return self.cvt.query([descriptor], k=1)[1][0][0]

    def add_to_archive(self, molecules) -> None:
        """
        Takes in a list of molecules and adds them to the archive as prescribed by the MAP-Elites algorithm, 
        i.e. each niche only contains the most fit molecule. Other molecules are discarded. 
        """
        for molecule in molecules:
            self.elites[self.cvt_index(molecule.descriptor)].update(molecule)
        return None

    def sample(self, size: int) -> List[Chem.Mol]:
        """
        Returns a list of elite molecules of the requisted length. 
        The elite molecules are randomly drawn, weighted by their fitness. 
        """
        pairs = [(elite.molecule, elite.molecule.fitness) for elite in self.elites if elite.molecule]
        molecules, weights = map(list, zip(*pairs))
        return random.choices(molecules, k=size, weights=weights)

    def sample_pairs(self, size: int) -> List[Tuple[Chem.Mol, Chem.Mol]]:
        """
        Returns a list of pairs of elite molecules of the requisted length. 
        The elite molecules are randomly drawn, weighted by their fitness. 
        """
        pairs = [(elite.molecule, elite.molecule.fitness) for elite in self.elites if elite.molecule]
        molecules, weights = map(list, zip(*pairs))
        sample_molecules = random.choices(molecules, k=size, weights=weights)
        sample_pairs = np.random.choice(list(filter(None, sample_molecules)), size=(size, 2), replace=True)
        sample_pairs = [tuple(sample_pair) for sample_pair in sample_pairs]       
        return sample_pairs

    def store_data(self, generation: float) -> None:
        """
        Creates a dataframe representing the archive and writes it to disk. In addtion, basic statistics about 
        the state of the archive are saved to disk and printed to the IO stream.           
        """
        archive_data = self.get_archive_data()
        fractional_size = len(archive_data["smiles"])/self.archive_size
        max_fitness, mean_fitness = np.max(archive_data["fitnesses"]), np.mean(archive_data["fitnesses"])
        if os.path.isfile('statistics.csv'):
            with open('statistics.csv', 'a') as file:
                csv.writer(file).writerow([generation, max_fitness, mean_fitness, fractional_size])
                file.close()
        else:
            with open('statistics.csv', 'w') as file:
                file.close()
        pd.DataFrame(data=archive_data).to_csv("archive_{}.csv".format(generation), index=False)
        print('Generation: {}, Size: {:.2f}'.format(generation, fractional_size))
        print('Fitness Max: {:.5f}, Fitness Mean: {:.5f}'.format(max_fitness, mean_fitness))
        return None
    
    def get_archive_data(self) -> None:
        elite_indices = [elite.index for elite in self.elites if elite.molecule]
        elite_molecules = [elite.molecule for elite in self.elites if elite.molecule]
        elites_smiles = [molecule.smiles for molecule in elite_molecules]
        elites_pedigree = [molecule.pedigree for molecule in elite_molecules]
        elites_descriptors = [molecule.descriptor for molecule in elite_molecules]
        elites_fitnesses = [molecule.fitness for molecule in elite_molecules]
        archive_data = {'index': elite_indices, 'smiles': elites_smiles, 'pedigree': elites_pedigree, 'descriptors': elites_descriptors, 'fitnesses': elites_fitnesses}
        return archive_data

class Arbiter:
    """
    A catalog class containing different druglike filters for small molecules.
    Includes the option to run the structural filters from ChEMBL.
    """
    def __init__(self, arbiter_config) -> None:
        self.cache_smiles = []
        self.rules_dict = pd.read_csv(hydra.utils.to_absolute_path("data/smarts/alert_collection.csv"))
        self.rules_dict= self.rules_dict[self.rules_dict.rule_set_name.isin(arbiter_config.rules)]
        self.rules_list = self.rules_dict["smarts"].values.tolist()
        self.tolerance_list = pd.to_numeric(self.rules_dict["max"]).values.tolist()
        self.pattern_list = [Chem.MolFromSmarts(smarts) for smarts in self.rules_list]

    def __call__(self, molecules):
        """
        Applies the chosen filters (hologenicity, veber_infractions,
        ChEMBL structural alerts, ...) to a list of molecules and removes duplicates.
        """
        filtered_molecules = []
        molecules = self.unique_molecules(molecules)
        for molecule in molecules:
            molecular_graph = Chem.MolFromSmiles(molecule.smiles)
            if self.molecule_filter(molecular_graph):
                filtered_molecules.append(molecule)
        return filtered_molecules

    def unique_molecules(self, molecules: List[Molecule]) -> List[Molecule]:
        """
        Checks if a molecule in a lost of molcules is duplicated, either in this batch or before.
        """
        unique_molecules = []
        for molecule in molecules:
            if molecule.smiles not in self.cache_smiles:
                unique_molecules.append(molecule)
                self.cache_smiles.append(molecule.smiles)
        return unique_molecules

    def molecule_filter(self, molecular_graph: Chem.Mol) -> bool:
        """
        Checks if a given molecular structure passes through the chosen filters (hologenicity,
        veber_infractions, ChEMBL structural alerts, ...).
        """
        toxicity = self.toxicity(molecular_graph)
        hologenicity = self.hologenicity(molecular_graph)
        veber_infraction = self.veber_infraction(molecular_graph)
        validity = not (toxicity or hologenicity or veber_infraction)
        if molecular_graph.HasSubstructMatch(Chem.MolFromSmarts('[R]')):
            ring_infraction = self.ring_infraction(molecular_graph)
            validity = validity and not (ring_infraction)
        return validity

    def toxicity(self, molecular_graph: Chem.Mol) -> bool:
        """
        Checks if a given molecule fails the structural filters.
        """
        for (pattern, tolerance) in zip(self.pattern_list, self.tolerance_list):
            if len(molecular_graph.GetSubstructMatches(pattern)) > tolerance:
                return True
        return False

    @staticmethod
    def hologenicity(molecular_graph: Chem.Mol) -> bool:
        """
        Checks if a given molecule fails the hologenicity filters.
        """
        fluorine_saturation = len(molecular_graph.GetSubstructMatches(Chem.MolFromSmarts('[F]'))) > 6
        bromide_saturation = len(molecular_graph.GetSubstructMatches(Chem.MolFromSmarts('[Br]'))) > 3
        chlorine_saturation = len(molecular_graph.GetSubstructMatches(Chem.MolFromSmarts('[Cl]'))) > 3
        return chlorine_saturation or bromide_saturation or fluorine_saturation

    @staticmethod
    def ring_infraction(molecular_graph: Chem.Mol) -> bool:
        """
        Checks if a given molecule fails the ring infraction filters.
        """
        ring_allene = molecular_graph.HasSubstructMatch(Chem.MolFromSmarts('[R]=[R]=[R]'))
        macro_cycle = max([len(j) for j in molecular_graph.GetRingInfo().AtomRings()]) > 6
        double_bond_in_small_ring = molecular_graph.HasSubstructMatch(Chem.MolFromSmarts('[r3,r4]=[r3,r4]'))
        return ring_allene or macro_cycle or double_bond_in_small_ring

    @staticmethod
    def veber_infraction(molecular_graph: Chem.Mol) -> bool:
        """
        Checks if a given molecule fails the veber infraction filters.
        """
        rotatable_bond_saturation = Lipinski.NumRotatableBonds(molecular_graph) > 10
        hydrogen_bond_saturation = Lipinski.NumHAcceptors(molecular_graph) + Lipinski.NumHDonors(molecular_graph) > 10
        return rotatable_bond_saturation or hydrogen_bond_saturation
