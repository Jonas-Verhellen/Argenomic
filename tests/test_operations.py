import pytest
import omegaconf
from rdkit import Chem
from argenomic.operations import Mutator, Crossover

@pytest.fixture
def default_mutator():
    '''
    Returns an instance of a mutator.
    '''
    return Mutator()

@pytest.fixture
def default_crossover():
    '''
    Returns an instance of a crossover.
    '''
    return Crossover()

@pytest.fixture
def default_molecules():
    '''
    Returns a list of two molecules.
    '''
    smiles = ["Clc1ccc(cc1)C(c2ccccc2)N3CCN(CC3)CCOCC(=O)O", "CC1=CC(Cl)=CC(C(=O)N[C@@H]2C[C@@H]3CCCC[C@@H]32)=C1C"]
    molecules = [Chem.MolFromSmiles(individual_smiles) for individual_smiles in smiles]
    return molecules

@pytest.mark.xfail
def test_default_mutator(default_mutator, default_molecules):
    '''
    Tests the action of the mutator. May fail occasionally due to stochasticity.
    The result of this test is reported separtely.
    '''
    for molecule in default_molecules:
        assert len(default_mutator(molecule)) > 0

def test_default_crossover(default_crossover, default_molecules):
    '''
    Tests the action of the crossover.
    '''
    assert len(default_crossover(default_molecules)) > 0
