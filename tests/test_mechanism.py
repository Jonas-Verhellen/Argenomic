import pytest
import omegaconf
from rdkit import Chem
from argenomic.mechanism import descriptor, fitness

@pytest.fixture
def default_descriptor():
    '''
    Returns a descriptor instance, set-up to calculate normalised ExactMolWt and MolLogP.
    '''
    configuration_file = omegaconf.OmegaConf.load("./tests/test_config.yaml")
    return descriptor(configuration_file.descriptor)

@pytest.fixture
def default_fitness():
    configuration_file = omegaconf.OmegaConf.load("./tests/test_config.yaml")
    return fitness(configuration_file.fitness)

@pytest.fixture
def default_molecules():
    '''
    Returns a list of two molecules.
    '''
    smiles = ["Clc1ccc(cc1)C(c2ccccc2)N3CCN(CC3)CCOCC(=O)O", "CC1=CC(Cl)=CC(C(=O)N[C@@H]2C[C@@H]3CCCC[C@@H]32)=C1C"]
    molecules = [Chem.MolFromSmiles(individual_smiles) for individual_smiles in smiles]
    return molecules

def test_default_descriptor(default_descriptor, default_molecules):
    descriptors = default_descriptor(default_molecules)
    for descriptor in descriptors:
        assert 0.00 <= descriptor
        assert descriptor <= 1.00

def test_default_descriptor(default_fitness, default_molecules):
    for molecule in default_molecules:
        fitness = default_fitness(molecule)
        assert 0.00 <= fitness
        assert fitness <= 1.00
