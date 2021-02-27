import hydra
import numpy as np
import pandas as pd
from typing import List, Tuple, Type

from rdkit import Chem
from rdkit.Chem import PandasTools as pdtl

from dask import bag
from dask.distributed import Client

from argenomic.base import Molecule
from argenomic.operations import Mutator, Crossover
from argenomic.infrastructure import Arbiter, Archive
from argenomic.mechanism import Fitness, Descriptor

class Illuminate:
    def __init__(self, config) -> None:
        self.data_file = config.data_file
        self.generations = config.generations
        self.batch_size = config.batch_size
        self.initial_size = config.initial_size

        self.arbiter = Arbiter(config.arbiter)
        self.fitness = Fitness(config.fitness)
        self.mutator = Mutator(config.mutator)
        self.crossover = Crossover()
        self.descriptor = Descriptor(config.descriptor)
        self.archive = Archive(config.archive, config.descriptor)

        self.client = Client(n_workers=config.workers, threads_per_worker=config.threads)
        return None

    def __call__(self) -> None:
        self.initial_population()
        for generation in range(1, self.generations):
            molecules = self.generate_molecules(generation)
            molecules = self.process_molecules(molecules)
            self.archive.add_to_archive(molecules)
            self.archive.store_data(generation)
        return None

    def initial_population(self) -> None:
        molecules = self.arbiter(self.load_from_database())
        molecules = self.calculate_descriptors(molecules)
        molecules = self.calculate_fitnesses(molecules)
        self.archive.add_to_archive(molecules)
        self.archive.store_data(0)
        return None

    def load_from_database(self) -> List[Molecule]:
        dataframe = pd.read_csv(hydra.utils.to_absolute_path(self.data_file))
        smiles_list = dataframe['smiles'].sample(n=self.initial_size).tolist()
        pedigree = ("database", "no reaction", "no parent")   
        molecules = [Molecule(Chem.CanonSmiles(smiles), pedigree) for smiles in smiles_list]
        return molecules

    def generate_molecules(self, generation) -> List[Molecule]:
        molecules = []
        molecule_samples = self.archive.sample(self.batch_size)
        molecule_sample_pairs = self.archive.sample_pairs(self.batch_size, generation)
        for molecule in molecule_samples:
            molecules.extend(self.mutator(molecule)) 
        for molecule_pair in molecule_sample_pairs:
            molecules.extend(self.crossover(molecule_pair)) 
        return molecules

    def process_molecules(self, molecules: List[Molecule]) -> List[Molecule]:
        molecules = self.arbiter(molecules)
        molecules = self.calculate_descriptors(molecules)
        molecules = self.calculate_fitnesses(molecules)
        return molecules

    def calculate_fitnesses(self, molecules: List[Molecule]) -> List[Molecule]:
        molecules = bag.map(self.fitness, bag.from_sequence(molecules)).compute()
        return molecules

    def calculate_descriptors(self, molecules: List[Molecule]) -> List[Molecule]:
        molecules = bag.map(self.descriptor, bag.from_sequence(molecules)).compute()
        molecules = [molecule for molecule in molecules if all(1.0 > property > 0.0 for property in molecule.descriptor)]
        return molecules

@hydra.main(config_path="configuration", config_name="config.yaml")
def launch(config) -> None:
    print(config.pretty())
    current_instance = Illuminate(config)
    current_instance()
    current_instance.client.close()

if __name__ == "__main__":
    launch()
