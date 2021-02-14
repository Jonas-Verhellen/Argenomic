from typing import List, Tuple
from dataclasses import dataclass

class Elite:
    def __init__(self, index):
        self.index = index
        self.molecule = None

    def update(self, molecule):
        if self.molecule is None or (molecule.fitness - self.molecule.fitness) > 0.0:
            self.molecule = molecule
        return None

@dataclass
class Molecule:
    smiles: str
    pedigree: Tuple[str, str ,str] 
    fitness: float = None
    descriptor: List[float] = None
