import omegaconf
import pandas as pd
from typing import List, Tuple

from rdkit import Chem
from rdkit.Chem import PandasTools as pdtl

from argenomic.operations import crossover, mutator
from argenomic.mechanism import descriptor, fitness
from argenomic.infrastructure import archive, arbiter

class illumination:
    def __init__(self, config: omegaconf.DictConfig) -> None:
        self.data_file = config.data_file
        self.batch_size = config.batch_size
        self.initial_size = config.initial_size

        self.mutator = mutator()
        self.crossover = crossover()
        self.arbiter = arbiter(config.arbiter)
        self.descriptor = descriptor(config.descriptor)
        self.archive = archive(config.archive, config.descriptor)
        self.fitness = fitness(config.fitness)
        return None

    def __call__(self, generations: int) -> None:
        self.initial_population()
        for generation in range(generations):
            molecules = self.generate_molecules()
            molecules, descriptors, fitnesses = self.process_molecules(molecules)
            self.archive.add_to_archive(molecules, descriptors, fitnesses)
            self.archive.store_statistics(generation)
            self.archive.store_archive(generation)
        return None

    def initial_population(self) -> None:
        dataframe = pd.read_csv(self.data_file)
        pdtl.AddMoleculeColumnToFrame(dataframe, 'smiles', 'molecule')
        molecules = dataframe['molecule'].sample(n=self.initial_size).tolist()
        molecules = self.arbiter(self.unique_molecules(molecules))
        molecules, descriptors, fitnesses = self.process_molecules(molecules)
        self.archive.add_to_archive(molecules, descriptors, fitnesses)
        return None

    def generate_molecules(self) -> None:
        molecules = []
        sample_molecules = self.archive.sample(self.batch_size)
        sample_molecule_pairs = self.archive.sample_pairs(self.batch_size)
        for molecule in sample_molecules:
            molecules.extend(self.mutator(molecule))
        for molecule_pair in sample_molecule_pairs:
            molecules.extend(self.crossover(molecule_pair))
        molecules = self.arbiter(self.unique_molecules(molecules))
        return molecules

    def process_molecules(self, molecules: List[Chem.Mol]) -> Tuple[List[List[float]],List[float]]:
        descriptors = list(map(self.descriptor, molecules))
        molecules, descriptors = zip(*[(molecule, descriptor) for molecule, descriptor in zip(molecules, descriptors)\
                if all(1.0 > property > 0.0 for property in descriptor)])
        molecules, descriptors = list(molecules), list(descriptors)
        fitnesses = list(map(self.fitness, molecules))
        return molecules, descriptors, fitnesses

    @staticmethod
    def unique_molecules(molecules: List[Chem.Mol]) -> List[Chem.Mol]:
        molecules = [Chem.MolFromSmiles(Chem.MolToSmiles(molecule)) for molecule in molecules if molecule is not None]
        molecule_records = [(molecule, Chem.MolToSmiles(molecule)) for molecule in molecules if molecule is not None]
        molecule_dataframe = pd.DataFrame(molecule_records, columns = ['molecules', 'smiles'])
        molecule_dataframe.drop_duplicates('smiles', inplace = True)
        return molecule_dataframe['molecules']

if __name__ == "__main__":
    configuration = omegaconf.OmegaConf.from_cli()
    configuration_file = omegaconf.OmegaConf.load(configuration.configuration_file)
    illumination(omegaconf.OmegaConf.load(configuration.configuration_file))(configuration.generations)
