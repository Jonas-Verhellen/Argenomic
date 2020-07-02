import pytest
import omegaconf
from rdkit import Chem
from argenomic.infrastructure import archive, arbiter

@pytest.fixture
def default_archive():
    configuration_file = omegaconf.OmegaConf.load("./test_config.json")
    '''
    Returns an archive instance of a 150 niches, spanned by ExactMolWt and MolLogP.
    '''
    return archive(configuration_file.archive, configuration_file.descriptor)

@pytest.fixture
def default_arbiter():
    configuration_file = omegaconf.OmegaConf.load("./test_config.json")
    '''
    Returns an arbiter instance, initialised with GSK structural alerts.
    '''
    return arbiter(configuration_file.arbiter)

@pytest.fixture
def default_molecules():
    smiles = ["Clc1ccc(cc1)C(c2ccccc2)N3CCN(CC3)CCOCC(=O)O", "CC1=CC(Cl)=CC(C(=O)N[C@@H]2C[C@@H]3CCCC[C@@H]32)=C1C"]
    molecules = [Chem.MolFromSmiles(individual_smiles) for individual_smiles in smiles]
    return molecules

def test_default_archive(default_archive, default_molecules):
    default_archive.add_to_archive(default_molecules, [[0.1, 0.1], [0.9, 0.9]], [0.0, 1.0])
    assert len(default_archive.sample(2)) == 2
    assert len(default_archive.sample_pairs(5)) == 5

def test_default_arbiter(default_arbiter, default_molecules):
    assert len(default_arbiter.pattern_list) == 55
    assert len(default_arbiter(default_molecules)) == 2
