import csv
import random
import itertools

import numpy as np
import pandas as pd
from typing import List, Tuple

from datetime import datetime
from sklearn.cluster import KMeans
from sklearn.neighbors import KDTree

from rdkit import Chem
from rdkit import rdBase
from rdkit.Chem import AllChem
rdBase.DisableLog('rdApp.error')
from rdkit.Chem import Lipinski

class elite():
    def __init__(self, index, descriptor):
        self.index = index
        self.fitness = 0.0
        self.molecule = None
        self.descriptor = descriptor

    def update(self, fitness, molecule, descriptor):
        if self.fitness < fitness:
            self.fitness = fitness
            self.molecule = molecule
            self.descriptor = descriptor
        return None

class archive:
    def __init__(self, archive_config, descriptor_config) -> None:
        self.archive_name = archive_config.name
        self.archive_size = archive_config.size
        kmeans = KMeans(n_clusters=self.archive_size, n_jobs=-1)
        kmeans = kmeans.fit(np.random.rand(archive_config.accuracy, len(descriptor_config.properties)))
        self.cvt_centers = kmeans.cluster_centers_
        self.cvt = KDTree(self.cvt_centers, metric='euclidean')
        self.elites = [elite(index, cvt_center) for index, cvt_center in enumerate(self.cvt_centers, start=0)]
        with open('{}/statistics.csv'.format(self.archive_name), 'w') as file:
            file.write("## Argenomic Statistics File: {}".format(datetime.now()))
            file.close()
        return None

    def cvt_index(self, descriptor: List[float]) -> int:
        return self.cvt.query([descriptor], k=1)[1][0][0]

    def add_to_archive(self, molecules: List[Chem.Mol], descriptors: List[List[float]], fitnesses: List[float]) -> None:
        for molecule, descriptor, fitness in zip(molecules, descriptors, fitnesses):
            self.elites[self.cvt_index(descriptor)].update(fitness, molecule, descriptor)
        return None

    def sample(self, size: int) -> List[Chem.Mol]:
        pairs = [(elite.molecule, elite.fitness) for elite in self.elites if elite.fitness > 0.0]
        molecules, weights = map(list, zip(*pairs))
        return random.choices(molecules, k=size, weights=weights)

    def sample_pairs(self, size: int) -> List[Tuple[Chem.Mol, Chem.Mol]]:
        pairs = [(elite.molecule, elite.fitness) for elite in self.elites if elite.fitness > 0.0]
        molecules, weights = map(list, zip(*pairs))
        sample_molecules = random.choices(molecules, k=size, weights=weights)
        sample_pairs = np.random.choice(list(filter(None, sample_molecules)), size=(size, 2), replace=True)
        sample_pairs = [tuple(sample_pair) for sample_pair in sample_pairs]
        return sample_pairs

    def store_archive(self, generation: float) -> None:
        elites_smiles, elites_descriptors, elites_fitnesses = self.elites_data()
        data = {'elites': elites_smiles, 'descriptors': elites_descriptors, 'fitnesses': elites_fitnesses}
        pd.DataFrame(data=data).to_csv("{}/archive_{}.csv".format(self.archive_name, generation), index=False)
        return None

    def store_statistics(self, generation: float) -> None:
        elites_smiles, elites_descriptors, elites_fitnesses = self.elites_data()
        fractional_size = len(elites_smiles)/self.archive_size
        statistics = [generation, np.max(elites_fitnesses), np.mean(elites_fitnesses), np.std(elites_fitnesses), fractional_size]
        with open('{}/statistics.csv'.format(self.archive_name), 'a') as file:
            csv.writer(file).writerow(statistics)
        print('Generation: {}, Size: {:.2f}'.format(statistics[0], statistics[4]))
        print('Fitness Max: {:.7f}, Mean: {:.7f}, Std: {:.7f}'.format(statistics[1], statistics[2], statistics[3]))
        return None

    def elites_data(self) -> Tuple[List[str], List[float], List[float]]:
        elites_list = [elite for elite in self.elites if elite.molecule]
        elites_smiles = [Chem.MolToSmiles(elite.molecule) for elite in elites_list]
        elites_descriptors = [elite.descriptor for elite in elites_list]
        elites_fitnesses = [elite.fitness for elite in elites_list]
        return elites_smiles, elites_descriptors, elites_fitnesses


class arbiter:
    """
    A catalog class containing different druglike filters for small molecules.
    Includes the option to run the structural filters from ChEMBL.
    """
    def __init__(self, arbiter_config) -> None:
      self.rules_dict = pd.read_csv("../data/smarts/alert_collection.csv")
      self.rules_dict= self.rules_dict[self.rules_dict.rule_set_name.isin(arbiter_config.rules)]
      self.rules_list = self.rules_dict["smarts"].values.tolist()
      self.tolerance_list = pd.to_numeric(self.rules_dict["max"]).values.tolist()
      self.pattern_list = [Chem.MolFromSmarts(smarts) for smarts in self.rules_list]

    def __call__(self, molecules:List[Chem.Mol]) -> List[Chem.Mol]:
      """
      Applies the chosen filters (hologenicity, veber_infractions,
      ChEMBL structural alerts, ...) to a list of molecules.
      """
      filtered_molecules = []
      for molecule in molecules:
        if self.molecule_validity(molecule):
          filtered_molecules.append(molecule)
      return filtered_molecules

    def molecule_validity(self, molecule: Chem.Mol) -> bool:
      """
      Checks if a given molecule passes through the chosen filters (hologenicity,
      veber_infractions, ChEMBL structural alerts, ...).
      """
      toxicity = self.toxicity(molecule)
      hologenicity = self.hologenicity(molecule)
      veber_infraction = self.veber_infraction(molecule)
      validity = not (toxicity or hologenicity or veber_infraction)
      if molecule.HasSubstructMatch(Chem.MolFromSmarts('[R]')):
        ring_infraction = self.ring_infraction(molecule)
        validity = validity and not (ring_infraction)
      return validity

    def toxicity(self, molecule: Chem.Mol) -> bool:
      """
      Checks if a given molecule fails the structural filters.
      """
      for (pattern, tolerance) in zip(self.pattern_list, self.tolerance_list):
            if len(molecule.GetSubstructMatches(pattern)) > tolerance:
              return True
      return False

    @staticmethod
    def hologenicity(molecule: Chem.Mol) -> bool:
      """
      Checks if a given molecule fails the hologenicity filters.
      """
      fluorine_saturation = len(molecule.GetSubstructMatches(Chem.MolFromSmarts('[F]'))) > 6
      bromide_saturation = len(molecule.GetSubstructMatches(Chem.MolFromSmarts('[Br]'))) > 3
      chlorine_saturation = len(molecule.GetSubstructMatches(Chem.MolFromSmarts('[Cl]'))) > 3
      return chlorine_saturation or bromide_saturation or fluorine_saturation

    @staticmethod
    def ring_infraction(molecule: Chem.Mol) -> bool:
      """
      Checks if a given molecule fails the ring infraction filters.
      """
      ring_allene = molecule.HasSubstructMatch(Chem.MolFromSmarts('[R]=[R]=[R]'))
      macro_cycle = max([len(j) for j in molecule.GetRingInfo().AtomRings()]) > 6
      double_bond_in_small_ring = molecule.HasSubstructMatch(Chem.MolFromSmarts('[r3,r4]=[r3,r4]'))
      return ring_allene or macro_cycle or double_bond_in_small_ring

    @staticmethod
    def veber_infraction(molecule: Chem.Mol) -> bool:
      """
      Checks if a given molecule fails the veber infraction filters.
      """
      rotatable_bond_saturation = Lipinski.NumRotatableBonds(molecule) > 10
      hydrogen_bond_saturation = Lipinski.NumHAcceptors(molecule) + Lipinski.NumHDonors(molecule) > 10
      return rotatable_bond_saturation or hydrogen_bond_saturation
