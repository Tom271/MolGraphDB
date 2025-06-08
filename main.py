"""Molecular graph database using SQLite and RDKit."""

from __future__ import annotations

import logging
import sqlite3
from collections import defaultdict
from contextlib import contextmanager
from pathlib import Path

from rdkit.Chem import FindAllSubgraphsOfLengthN, GetPeriodicTable, PathToSubmol
from rdkit.Chem.MolStandardize.rdMolStandardize import LargestFragmentChooser
from rdkit.Chem.rdchem import Mol, RWMol
from rdkit.Chem.rdmolfiles import MolFromSmiles, MolToSmiles
from rdkit.Chem.rdmolops import SanitizeMol
from rdkit.RDLogger import DisableLog, EnableLog
from rich.logging import RichHandler


@contextmanager
def suppress_rdkit_logs():
    """Context manager to suppress RDKit logs."""
    DisableLog("rdApp.*")
    try:
        yield
    finally:
        EnableLog("rdApp.*")


periodic_table = GetPeriodicTable()
chooser = LargestFragmentChooser()
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    handlers=[RichHandler()],
    force=True,
)
logger = logging.getLogger(__name__)


class Molecule:
    """Class representing a molecule with SMILES representation."""

    def __init__(self, smiles: str):
        self.smiles = smiles
        self.mol: Mol = MolFromSmiles(smiles)
        self.n_atoms = self.mol.GetNumAtoms()
        self.n_bonds = self.mol.GetNumBonds()

    def __repr__(self):
        return self.smiles

    def get_subgraphs(self, min_atoms: int = 2):
        """Generate all subgraphs of the molecule."""

        for subgraph in FindAllSubgraphsOfLengthN(self.mol, length=min_atoms - 1):
            mol = PathToSubmol(self.mol, subgraph)
            yield from self._check_and_return(mol, min_atoms=min_atoms)

    def get_bond_deletions(self, min_atoms: int = 2):
        """Generate all possible molecules by deleting one bond."""
        for bond in self.mol.GetBonds():
            rw_mol = RWMol(self.mol)
            rw_mol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            with suppress_rdkit_logs():
                chooser.chooseInPlace(rw_mol)
                yield from self._check_and_return(rw_mol, min_atoms=min_atoms)

    def get_atom_deletions(self, min_atoms: int = 2):
        """Generate all possible molecules by deleting one atom."""
        for atom in self.mol.GetAtoms():
            rw_mol = RWMol(self.mol)
            rw_mol.RemoveAtom(atom.GetIdx())
            yield from self._check_and_return(rw_mol, min_atoms=min_atoms)

    def get_atom_substitutions(self, elem_from="C", elem_to="N", min_atoms: int = 2):
        """Generate all possible molecules by substituting one atom."""
        num_from = periodic_table.GetAtomicNumber(elem_from)
        num_to = periodic_table.GetAtomicNumber(elem_to)
        for atom in self.mol.GetAtoms():
            rw_mol = RWMol(self.mol)
            if atom.GetAtomicNum() == num_from:
                rw_mol.GetAtomWithIdx(atom.GetIdx()).SetAtomicNum(num_to)
                yield from self._check_and_return(rw_mol, min_atoms=min_atoms, subs=1)

    def _check_and_return(
        self, next_mol: Mol, min_atoms: int = 2, min_bonds: int = 1, subs: int = 0
    ):
        """Check the next generated molecule and return if valid."""

        for next_mol in [next_mol]:
            with suppress_rdkit_logs():
                try:
                    SanitizeMol(next_mol, catchErrors=False)
                    next_smiles = MolToSmiles(next_mol)
                    next_n_atoms = next_mol.GetNumAtoms()
                    next_n_bonds = next_mol.GetNumBonds()
                except Exception as exception:
                    continue

            if "." in next_smiles:
                continue

            if next_n_atoms < min_atoms or next_n_bonds < min_bonds:
                continue

            next_mol = Molecule(next_smiles)

            diff_atoms = self.n_atoms - next_n_atoms
            diff_bonds = self.n_bonds - next_n_bonds
            substitution = subs

            yield next_mol, self, diff_atoms, diff_bonds, substitution
            yield self, next_mol, -diff_atoms, -diff_bonds, substitution


class PersistentDataBase:
    """Persistent molecular database using SQLite."""

    def __init__(self, db_path: str = "molecular_graph.db"):
        self.db_path = Path(db_path)

        # Initialize SQLite connection
        self.conn = sqlite3.connect(str(self.db_path))
        self.conn.execute("PRAGMA journal_mode=WAL")  # Better concurrency
        self.conn.execute("PRAGMA synchronous=NORMAL")  # Better performance
        self.conn.execute("PRAGMA cache_size=10000")  # Larger cache

        # Create tables
        self._init_database()

        # In-memory cache for performance
        self._cache = defaultdict(set)
        self._load_cache()

        logger.info(f"Connected to SQLite database at {self.db_path}")

    def _init_database(self):
        """Initialize the SQLite database schema."""
        cursor = self.conn.cursor()

        # Create main table for molecular relations
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS molecular_relations (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                smiles TEXT NOT NULL,
                target_smiles TEXT NOT NULL,
                diff_atom INTEGER NOT NULL,
                diff_bond INTEGER NOT NULL,
                subs INTEGER NOT NULL DEFAULT 0,
                UNIQUE(smiles, target_smiles)
            )
        """)

        # Create indexes for better performance
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_smiles ON molecular_relations(smiles)
        """)
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_target ON molecular_relations(target_smiles)
        """)

        # Create metadata table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS metadata (
                key TEXT PRIMARY KEY,
                value TEXT
            )
        """)

        self.conn.commit()

    def _load_cache(self):
        """Load data from SQLite into memory cache."""
        logger.info("Loading database cache...")
        cursor = self.conn.cursor()

        cursor.execute("""
            SELECT smiles, target_smiles, diff_atom, diff_bond, subs
            FROM molecular_relations
        """)

        count = 0
        for smiles, target_smiles, diff_atom, diff_bond, subs in cursor:
            self._cache[smiles].add((target_smiles, diff_atom, diff_bond, subs))
            count += 1

        # logger.info("Loaded %d molecular relations from database", count)

    def _save_batch_to_db(self, batch_data):
        """Save a batch of relations to SQLite."""
        cursor = self.conn.cursor()

        # Prepare data for batch insert
        insert_data = []
        for smiles in batch_data:
            for results in batch_data[smiles]:
                insert_data.append((smiles, *results))

        if insert_data:
            cursor.executemany(
                """
                INSERT OR IGNORE INTO molecular_relations
                (smiles, target_smiles, diff_atom, diff_bond, subs)
                VALUES (?, ?, ?, ?, ?)
            """,
                insert_data,
            )

            self.conn.commit()

    def add_relations(
        self, molecules: list[Molecule], min_atoms: int = 2, depth=100
    ) -> int:
        """Add all subgraphs of given molecules to the database."""
        wavefront = set(mol for mol in molecules if mol.smiles not in self._cache)
        if len(wavefront) == 0:
            return 0

        logger.info(f"Adding {len(wavefront)} molecules to the database...")
        wave = 0
        batch_data = defaultdict(set)

        while wavefront and wave < depth:
            wave += 1
            logger.info(f"Wave {wave}: {len(wavefront)} molecules")
            wavefront = self.add_next_wavefront(wavefront, min_atoms=min_atoms)

            # Collect new data for batch saving
            for smiles in self._cache:
                for results in self._cache[smiles]:
                    batch_data[smiles].add(results)

            # Batch save to SQLite every waves
            self._save_batch_to_db(batch_data)
            batch_data.clear()

        # Final save
        if batch_data:
            self._save_batch_to_db(batch_data)

        return len(self._cache)

    def add_next_wavefront(
        self, mols: set[Molecule], min_atoms: int = 2
    ) -> set[Molecule]:
        """Add next wavefront of molecules to the graph."""
        added_mols = set()
        for mol in mols:
            # for mol_src, rel, mol_tgt in mol.get_subgraphs(length=min_atoms - 1):
            #     added_mols.update(self.add_and_return_new(mol_src,   mol_tgt))
            for mod_results in mol.get_atom_substitutions(
                "C", "N", min_atoms=min_atoms
            ):
                added_mols.update(self.add_and_return_new(*mod_results))
            for mod_results in mol.get_bond_deletions(min_atoms=min_atoms):
                added_mols.update(self.add_and_return_new(*mod_results))
            for mod_results in mol.get_atom_deletions(min_atoms=min_atoms):
                added_mols.update(self.add_and_return_new(*mod_results))
        return added_mols

    def add_and_return_new(
        self, mol_src: Molecule, mol_tgt: Molecule, d_atoms, d_bonds, subs
    ):
        """Add a new molecular relation if it doesn't exist."""
        results = (mol_tgt.smiles, d_atoms, d_bonds, subs)
        if results not in self._cache.get(mol_src.smiles, set()):
            self._cache[mol_src.smiles].add(results)
            return {mol_tgt}
        return set()

    def query(self, molecule_from: Molecule, molecule_to: Molecule) -> int:
        """Find the shortest path between two molecules in the database."""
        # Ensure that related molecules are in the database
        self.add_relations([molecule_from, molecule_to])

        # Search for the shortest path breadth-first
        logger.info(
            "Searching for path from %s to %s...",
            molecule_from.smiles,
            molecule_to.smiles,
        )
        previous_wavefront = set()
        current_wavefront = {molecule_from.smiles}
        next_wavefront = set()
        wave = 0

        while current_wavefront:
            wave += 1
            logger.info("Wave %d: %d molecules", wave, len(current_wavefront))
            for smiles in current_wavefront:
                if smiles == molecule_to.smiles:
                    return wave
                for results in self._cache[smiles]:
                    related_smiles = results[0]
                    if related_smiles in previous_wavefront:
                        continue
                    if related_smiles in current_wavefront:
                        continue
                    next_wavefront.add(related_smiles)
            previous_wavefront.update(current_wavefront)
            current_wavefront = next_wavefront
            next_wavefront = set()

        return -1

    def get_stats(self) -> dict:
        """Get database statistics."""
        cursor = self.conn.cursor()

        # # Total relations
        # cursor.execute("SELECT COUNT(*) FROM molecular_relations")
        # total_relations = cursor.fetchone()[0]

        # Unique molecules
        cursor.execute("SELECT COUNT(DISTINCT smiles) FROM molecular_relations")
        unique_molecules = cursor.fetchone()[0]

        # # Relations by type
        # cursor.execute("""
        #     SELECT relation, COUNT(*)
        #     FROM molecular_relations
        #     GROUP BY relation
        # """)
        # relations_by_type = dict(cursor.fetchall())

        return {
            # "total_relations": total_relations,
            "unique_molecules": unique_molecules,
            # "relations_by_type": relations_by_type,
            "cache_size": len(self._cache),
        }

    @property
    def number_entries(self) -> int:
        """Return the number of molecules in the database."""
        return len(self._cache)

    def close(self):
        """Close the database connection."""
        if hasattr(self, "conn"):
            # Final save
            batch_data = defaultdict(set)
            for smiles in self._cache:
                for results in self._cache[smiles]:
                    batch_data[smiles].add(results)

            if batch_data:
                self._save_batch_to_db(batch_data)

            self.conn.close()
            logger.info("Database connection closed")

    def __repr__(self):
        return f"PersistentDataBase(num_entries={self.number_entries})"

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


# Backward compatibility alias
DataBase = PersistentDataBase


if __name__ == "__main__":
    # Initialize persistent database
    with PersistentDataBase("molecular_graph.db") as database:
        # Show database statistics
        stats = database.get_stats()
        logger.info(f"Database statistics: {stats}")

        # Test with simple molecules
        logger.info("Testing with simple molecules:")
        database.add_relations(
            [Molecule("CCC"), Molecule("CC"), Molecule("C(C)C")], min_atoms=1
        )
        result = database.query(Molecule("CCCCCC"), Molecule("CC(CC)C"))
        logger.info(f"Path length: {result}")

        # Test with pharmaceutical molecules
        logger.info("Testing with pharmaceutical molecules:")
        smiles_dict = {
            "sildenafil": r"CCCC1=NN(C2=C1N=C(NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C",
            "vardenafil": r"O=C2\N=C(/Nn1c(nc(c12)C)CCC)c3cc(ccc3OCC)S(=O)(=O)N4CCN(CC)CC4",
            "tadalafil": r"CN1CC(=O)N2[C@H](Cc3c([nH]c4ccccc34)[C@H]2c2ccc3c(c2)OCO3)C1=O",
            "aspirin": r"O=C(C)Oc1ccccc1C(=O)O",
            "ibuprofen": r"CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O",
        }

        aspirin = Molecule(smiles_dict["aspirin"])
        ibuprofen = Molecule(smiles_dict["ibuprofen"])
        logger.info(f"Aspirin atoms: {aspirin.n_atoms}")
        logger.info(f"Ibuprofen atoms: {ibuprofen.n_atoms}")

        # Add molecules to database with high minimum atom count
        database.add_relations([aspirin])
        database.add_relations([ibuprofen])

        # Query for path between sildenafil and vardenafil
        result = database.query(aspirin, ibuprofen)
        logger.info(f"Path length between aspirin and ibuprofen: {result}")

        sildenafil = Molecule(smiles_dict["sildenafil"])
        vardenafil = Molecule(smiles_dict["vardenafil"])
        tadalafil = Molecule(smiles_dict["tadalafil"])
        logger.info(f"Sildenafil atoms: {sildenafil.n_atoms}")
        logger.info(f"Vardenafil atoms: {vardenafil.n_atoms}")
        logger.info(f"Tadalafil atoms: {tadalafil.n_atoms}")

        database.add_relations([sildenafil], depth=5)
        database.add_relations([vardenafil], depth=5)
        database.add_relations([tadalafil], depth=5)

        result = database.query(sildenafil, vardenafil)
        logger.info(f"Path length between sildenafil and vardenafil: {result}")
